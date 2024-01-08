#!/usr/bin/env python3
#
# PmUtils/pmmerge
#

"""
Created on Mar 17, 2021

@author: kovi
"""

from getopt import getopt, GetoptError
from sys import argv
from astropy.table import Table, join, setdiff
from astropy.io import fits
from os import remove
from os.path import exists
from math import log10

import pmbase as pm
from pmviz import VizUCAC4


class CatalogMatcher:
    folder = './'
    baseName = None
    colors = ['Ri', 'Gi', 'Bi']
    showGraphs = False
    saveGraphs = False

    # TODO: ezt a kettot konfigba tenni
    HMG_MAX_ERR = 0.15  # mg

    LIMIT_THRE = 0.5  # mg

    IMAGE_BORDER_SIZE = 10  # pixels

    MG_LIMIT = 17.0  # TODO: nem lehet-e ezt szamolni?

    COORD_MATCH_THRE = 2.0 / 60.0  # degrees

    viz = None

    def __init__(self, baseName, folder, colors, show, save, loggingMode=pm.Logger.LOG_MODE_INFO):
        self.baseName = baseName
        if folder is not None:
            self.folder = folder
        if colors is not None:
            self.colors = colors
        self.viz = VizUCAC4(self.MG_LIMIT)
        self.showGraphs = show
        self.saveGraphs = save
        self.log = pm.Logger(mode=loggingMode)

    def hmg(self, pmTable, nr):
        mgIsoCorrArr = pmTable['MAG_ISOCOR']
        mgIsoCorrErrArr = pmTable['MAGERR_ISOCOR']
        mgBestArr = pmTable['MAG_BEST']
        mgBestErrArr = pmTable['MAGERR_BEST']
        logm1 = []
        loge1 = []
        logm2 = []
        loge2 = []
        for j in range(len(pmTable)):
            if mgIsoCorrArr[j] <= 0.0:
                logm1.append(log10(-mgIsoCorrArr[j]))
                loge1.append(log10(mgIsoCorrErrArr[j]))
            if mgBestArr[j] <= 0.0:
                logm2.append(log10(-mgBestArr[j]))
                loge2.append(log10(mgBestErrArr[j]))

        m1, b1 = pm.linfit(logm1, loge1)
        hmg1 = -1.0 * 10.0 ** ((log10(self.HMG_MAX_ERR) - b1) / m1)

        m2, b2 = pm.linfit(logm2, loge2)
        hmg2 = -1.0 * 10.0 ** ((log10(self.HMG_MAX_ERR) - b2) / m2)

        plotcolor = 'gbr'[nr]
        self.hmgPlot.add(logm1, loge1, [m1, b1], f'log mi ({plotcolor.upper()}i)', f'log ei ({plotcolor.upper()}i)',
                         plotcolor, pm.Plot.INV_X, title="Mg limit")

        return [hmg1, hmg2]

    def calcMgLimit(self, fileName, nr):
        t = Table.read(fileName, format='ascii.sextractor')
        mgs = self.hmg(t, nr)
        return mgs[0]

    def readFrameSize(self, fitsFileName):
        h = pm.getFitsHeaders(fitsFileName, ['NAXIS1', 'NAXIS2'])
        xf = int(h['NAXIS1'])
        yf = int(h['NAXIS2'])
        return xf, yf

    def determineOnFrameStatus(self, x_obj, y_obj, im_w, im_h):
        if x_obj < 0 or x_obj >= im_w or y_obj < 0 or y_obj >= im_h:
            return 'O'
        if x_obj < self.IMAGE_BORDER_SIZE or x_obj >= im_w - self.IMAGE_BORDER_SIZE or y_obj < self.IMAGE_BORDER_SIZE \
                or y_obj >= im_h - self.IMAGE_BORDER_SIZE:
            return 'B'
        return '-'

    def mergeCatalogs(self, refFile, matchFile, outFile):
        # match ref and input
        r = self.invokeGrMatch(refFile, matchFile, outFile)

        if not r.startswith('ERROR'):
            # merge exclusions
            self.mergeExclusions(outFile)

        return r

    def mergeExclusions(self, outFile):
        # TODO
        xrfExclusionFile = outFile + '.xrf'
        inpExclusionFile = outFile + '.xin'

        cat = Table.read(outFile, format='ascii.no_header', delimiter=' ')
        refEx = Table.read(xrfExclusionFile, format='ascii.no_header', delimiter=' ')
        refExCount = len(refEx)

        pass

    def invokeGrMatch(self, refFile, matchFile, outFile):
        # TODO: add non-matched objects from xrf and xin file
        matchOptions = '--match-coord --col-ref 14,15 --col-inp 14,15'
        return pm.invoke(
                f"grmatch {matchOptions} -r {refFile} -i {matchFile} -o {outFile} --output-excluded-reference {outFile + '.xrf'} --output-excluded-input {outFile + '.xin'}")

    def calcVisFlag(self, mg, limit):
        """
        flags I - INVISIBLE, B - BELOW_LIMIT, N - NEAR_LIMIT, S - SATURATED
        """
        # TODO 1.3: saturated objects
        mgf = float(mg)
        if mgf > 0.0:
            return 'I'
        elif mgf > limit:
            return 'B'
        elif mgf > limit - self.LIMIT_THRE:
            return 'N'
        return '-'

    def calcPosFlag(self, x_obj, y_obj, im_w, im_h):
        """
        x_obj # object x position
        y_obj # object y position
        im_w  # image width
        im_h  # image height

        flags B - object is near of image border
        """

        return self.determineOnFrameStatus(x_obj, y_obj, im_w, im_h)

    def matchCatalogs(self):
        catFileBase = f"{self.folder}/{pm.setup['PHOT_FOLDER_NAME']}/{self.baseName}-"
        tempFileBase = f"{self.folder}/Temp/{self.baseName}-"
        pm.assureFolder(f"{self.folder}/Temp")

        # query image size
        img_w, img_h = self.readFrameSize(catFileBase + 'Gi.ast.fits')
        self.log.debug(f'Image size: {img_w} x {img_h}')

        # calculate mg limits
        self.hmgPlot = pm.Plot(len(self.colors), self.showGraphs, self.saveGraphs)

        #        if len(self.opt['color']) == 3:
        #        hmgGi = self.calcMgLimit(catFileBase + 'Gi.cat', 0)
        #        hmgBi = self.calcMgLimit(catFileBase + 'Bi.cat', 1)
        #        hmgRi = self.calcMgLimit(catFileBase + 'Ri.cat', 2)
        #        self.mgLimits = { 'Gi': hmgGi, 'Bi': hmgBi, 'Ri': hmgRi }
        #        self.log.info(f'Limit mg Gi: {hmgGi}, Bi: {hmgBi}, Ri: {hmgRi}')
        self.mgLimits = {}
        s = ''
        for j, clr in enumerate(self.colors):
            self.mgLimits[clr] = self.calcMgLimit(catFileBase + clr + '.cat', j)
            s += f'{clr}: {self.mgLimits[clr]}, '
        self.log.info(f'Limit mg {s[:-2]}')

        self.hmgPlot.showOrSave(self.folder + '/magnitude_limit.png')

        # TODO: Azok a csillagok, amik nem latszanak Gi-ben, csak Bi-ben vagy Ri-ben, azok is keruljenek bele

        # match Bi cat to Gi
        self.mergeCatalogs(catFileBase + 'Gi.cat', catFileBase + 'Bi.cat', tempFileBase + 'GiBi.cat')

        # match Ri cat to combined Gi + Bi
        self.mergeCatalogs(tempFileBase + 'GiBi.cat', catFileBase + 'Ri.cat', tempFileBase + 'GiBiRi.cat')

        # load combined cat

        table = Table.read(tempFileBase + 'GiBiRi.cat', format='ascii.no_header', delimiter=' ')
        count = len(table)
        self.log.debug(f'{count} matched objects found')

        # create cmb catalog
        cmb = self.createCombinedCat(count)

        # fill cmb data from combined cat
        cmb['RA_DEG'] = table['col14']
        cmb['DEC_DEG'] = table['col15']
        cmb['MAG_GI'] = table['col4']
        cmb['ERR_GI'] = table['col5']
        cmb['MAG_BI'] = table['col27']
        cmb['ERR_BI'] = table['col28']
        cmb['MAG_RI'] = table['col50']
        cmb['ERR_RI'] = table['col51']

        ras = []
        des = []
        auid = []
        n_auid = 1
        visFlags = []
        posFlags = []
        for r in cmb:
            auid.append('%03d-FFF-%03d' % (n_auid // 1000, n_auid % 1000))
            n_auid += 1
            ra_s = pm.deg2hexa(float(r['RA_DEG']) / 15.0)
            de_s = pm.deg2hexa(float(r['DEC_DEG']))
            ras.append(ra_s)
            des.append(de_s)
            posFlags.append(
                    self.calcPosFlag(float(table[r.index]['col16']), float(table[r.index]['col17']), img_w, img_h))
            visFlags.append((self.calcVisFlag(r['MAG_GI'], self.mgLimits['Gi']) if 'Gi' in self.colors else " ") +
                            (self.calcVisFlag(r['MAG_BI'], self.mgLimits['Bi']) if 'Bi' in self.colors else " ") +
                            (self.calcVisFlag(r['MAG_RI'], self.mgLimits['Ri']) if 'Ri' in self.colors else " "))

        cmb['AUID'] = auid
        cmb['RA'] = ras
        cmb['DEC'] = des
        cmb['POS_FLAG'] = posFlags
        cmb['VIZ_FLAG'] = visFlags

        return cmb

    def createCombinedCat(self, count=0):
        fnames = ['AUID', 'VIZ_ID', 'ROLE', 'LABEL',
                  'VIZ_FLAG', 'POS_FLAG', 'MATCH_FLAG',
                  'RA', 'RA_DEG', 'DEC', 'DEC_DEG',
                  'MAG_GI', 'ERR_GI', 'MAG_BI', 'ERR_BI', 'MAG_RI', 'ERR_RI',
                  'MAG_B', 'ERR_B', 'MAG_V', 'ERR_V', 'MAG_R', 'ERR_R']
        fdtype = ['U11', 'U16', 'U2', 'U50',
                  'U1', 'U1', 'U1',
                  'U12', 'U12', 'U12', 'U12',
                  'U6', 'U5', 'U6', 'U5', 'U6', 'U5',
                  'U6', 'U5', 'U6', 'U5', 'U6', 'U5']
        if count > 0:
            rows = [['AUID'] * count, ['-'] * count, ['F'] * count, ['-'] * count,
                    ['-'] * count, ['-'] * count, ['T'] * count,
                    [''] * count, [0.0] * count, [''] * count, [0.0] * count,
                    ['99.0'] * count, ['-'] * count, ['99.0'] * count, ['-'] * count, ['99.0'] * count, ['-'] * count,
                    ['99.0'] * count, ['-'] * count, ['99.0'] * count, ['-'] * count, ['99.0'] * count, ['-'] * count]
            return Table(rows, names=fnames, dtype=fdtype)
        else:
            return Table(names=fnames, dtype=fdtype)

    def xmatchViz(self, combined):
        # xmatch cmb to UCAC4 catalog
        # TODO: err values from xmatch
        xmt = self.viz.xmatch(combined, 'RA_DEG', 'DEC_DEG')
        if not xmt:
            self.log.error("Accessing Vizier service if failed.")
            return
        self.log.debug(f'Xmatch table contains {len(xmt)} records')

        combined.add_index('AUID')

        for xr in xmt:
            auid = xr['AUID']
            r = combined.loc[auid]
            r['VIZ_ID'] = 'UCAC4-' + xr['UCAC4']
            r['MATCH_FLAG'] = '-'
            r['MAG_B'] = xr['Bmag']
            r['MAG_V'] = xr['Vmag']
            r['MAG_R'] = xr['Rmag']

    # TODO: query errors from vizier
    #            r['ERR_B'] = xr['e_Bmag']
    #            r['ERR_V'] = xr['e_Vmag']
    #            r['ERR_R'] = xr['e_Rmag']

    def loadRefcat(self):
        refcatFileName = f"{self.folder}/ref.cat"
        if exists(refcatFileName):
            return Table.read(refcatFileName, format='ascii')
        else:
            self.log.warning(f'No ref.cat in the folder {self.folder}')
            return None

    def addFrameCoords(self, cat, fitsFileName):
        # 1. convert cat to fits format
        # baseFolder = fitsFileName.partition(pm.setup['PHOT_FOLDER_NAME'])[0]
        tempFile = self.folder + '/Temp/tempCat.fits'
        if not exists(tempFile):
            cat.write(tempFile)

        # 2. calculate ref objects' frame xy points
        wcsFile = fitsFileName.replace('.ast.fits', '.wcs')
        axyFile = fitsFileName.replace('.ast.fits', '.ref.axy')
        pm.invoke("wcs-rd2xy -w %s -i %s -o %s -R RA_DEG -D DEC_DEG" % (
            wcsFile, tempFile, axyFile))  # -f option need argument

        remove(tempFile)

        # 3. merge frame xy point to refCat
        tlen = len(cat)
        cat['X'] = [0.0] * tlen
        cat['Y'] = [0.0] * tlen

        f = fits.open(axyFile)
        d = f[1].data

        for j in range(tlen):
            x, y = d[j]
            cat[j]['X'] = x
            cat[j]['Y'] = y

    def filterOutFrameStars(self, cat, fitsFileName):
        # add frame coords to the catalog
        self.addFrameCoords(cat, fitsFileName)

        # read frame size from fits file
        fx, fy = self.readFrameSize(fitsFileName)

        cat['FRAME_STATUS'] = ['-'] * len(cat)
        for row in cat:
            row['FRAME_STATUS'] = self.determineOnFrameStatus(row['X'], row['Y'], fx, fy)

        foMask = cat['FRAME_STATUS'] != 'O'
        foCat = cat[foMask]
        #        print('Out of frame objects:')
        #        print(foCat)

        # foMask = cat['FRAME_STATUS'] != 'O'
        # foCat = cat[foMask]
        return foCat

    def matchCatalogByCoords(self, cat, ra, dec):
        dmin = 99.0
        auid = None
        for row in cat:
            if row['VIZ_ID'] == '-':
                cra = float(row['RA_DEG'])
                cdec = float(row['DEC_DEG'])
                d = pm.quad(cra - ra, cdec - dec)
                if d < dmin:
                    dmin = d
                    auid = row['AUID']
        #        print(f"Coord matched - auid: {auid}, d: {sqrt(dmin)*60}'")
        return auid, dmin

    def updateCmb(self, cmb, r, auid=None, vizid=None):
        if vizid:
            cr = cmb.loc['VIZ_ID', vizid]
        else:
            cr = cmb.loc['AUID', auid]
        cr['ROLE'] = r['ROLE']
        cr['AUID'] = r['AUID']
        cr['LABEL'] = r['LABEL']
        cr['MAG_B'] = r['MAG_B']
        cr['MAG_V'] = r['MAG_V']
        cr['MAG_R'] = r['MAG_R']
        cr['ERR_B'] = r['ERR_B']
        cr['ERR_V'] = r['ERR_V']
        cr['ERR_R'] = r['ERR_R']

    #        r['MATCH_FLAG'] = xr['MATCH_FLAG']

    def insertCmb(self, cmb, r):
        n = [r['AUID'], '-', r['ROLE'], r['LABEL'], 'III', '-', 'N', r['RA'], r['RA_DEG'], r['DEC'], r['DEC_DEG'],
             '99.0', '99.0', '99.0', '99.0', '99.0', '99.0', r['MAG_B'], r['ERR_B'], r['MAG_V'], r['ERR_V'], r['MAG_R'],
             r['ERR_R']]
        cmb.add_row(n)

    def addRefcatData(self, cmb, refcat):

        del refcat.meta['comments']

        fitsFileName = f"{self.folder}/{pm.setup['PHOT_FOLDER_NAME']}/{self.baseName}-Gi.ast.fits"
        foRefcat = self.filterOutFrameStars(refcat, fitsFileName)

        # xmatch refcat to UCAC4
        # TODO: drop out records, where angDist comlumn value is grater then a threshold (bad match)
        xmt = self.viz.xmatch(foRefcat, 'RA_DEG', 'DEC_DEG')
        self.log.debug(f'Refcat xmatch table contains {len(xmt)} records')

        missingRefcat = setdiff(foRefcat, xmt, 'AUID')

        cmb.add_index('VIZ_ID')

        for xr in xmt:
            vizId = 'UCAC4-' + xr['UCAC4']
            if vizId in cmb['VIZ_ID']:
                try:
                    self.updateCmb(cmb, xr, vizid=vizId)

                except KeyError:
                    self.log.debug(xr)

            # else:
            # self.log.debug('xmt record not matched in cmb:')
            # print(xr)  # RBL

        mergedTable = join(xmt, cmb, keys='AUID', join_type='left')
        mask = mergedTable['VIZ_ID'] == '-'  # or mergedTable['VIZ_ID'] == '' or mergedTable['VIZ_ID'] == None
        nt = mergedTable[mask]

        self.log.debug(f'Negative transients: {len(nt)} found.')

        # refCat: AUID         ROLE RA             RA_DEG         DEC            DEC_DEG        MAG_B     ERR_B     MAG_V     ERR_V     MAG_R     ERR_R     LABEL
        # cmb:    'AUID','VIZ_ID','ROLE','LABEL','VIZ_FLAG','POS_FLAG','MATCH_FLAG','RA','RA_DEG','DEC','DEC_DEG',
        #        'MAG_GI','ERR_GI','MAG_BI','ERR_BI','MAG_RI','ERR_RI','MAG_B','ERR_B','MAG_V','ERR_V','MAG_R','ERR_R']

        for r in nt:
            n = [r['AUID'], '-', r['ROLE_1'], r['LABEL_1'], r['VIZ_FLAG'], r['POS_FLAG'], r['MATCH_FLAG'], r['RA_1'],
                 r['RA_DEG_1'], r['DEC_1'], r['DEC_DEG_1'],
                 r['MAG_GI'], r['ERR_GI'], r['MAG_BI'], r['ERR_BI'], r['MAG_RI'], r['ERR_RI'],
                 r['MAG_B_1'], r['ERR_B_1'], r['MAG_V_1'], r['ERR_V_1'], r['MAG_R_1'], r['ERR_R_1']]
            cmb.add_row(n)

        self.log.debug(f'Stars missing from refcat: {len(missingRefcat)} found.')
        for r in missingRefcat:
            auid, d = self.matchCatalogByCoords(cmb, float(r['RA_DEG']), float(r['DEC_DEG']))
            self.log.debug(f"Match missingRefcat object {r['AUID']} {r['LABEL']} for {auid} within {d * 60}'")
            # print(r)  # RBL
            if d < 2.0 / 60.0:
                self.updateCmb(cmb, r, auid=auid)
            else:
                self.insertCmb(cmb, r)

    def addMgLimits(self, cmb):
        for cc in self.mgLimits.keys():
            pm.addTableComment(cmb, 'MgLimitInst' + cc, '%7.3f' % (self.mgLimits[cc]))

    def save(self, combined):
        combinedFileName = f"{self.folder}/{pm.setup['PHOT_FOLDER_NAME']}/{self.baseName}.cmb"
        self.log.debug(f"Combined file name: {combinedFileName}")
        combined.write(combinedFileName, format='ascii.fixed_width', delimiter=' ', overwrite=True)

    def process(self):

        # match cat files together
        self.log.print('Combine measurement catalogs all together')
        combined = self.matchCatalogs()

        # xmatch combined catalog to UCAC4, and store UCAC4 id
        self.log.print('Match combined catalog with UCAC4')
        self.xmatchViz(combined)

        # load refcat
        refcat = self.loadRefcat()

        # add refcat data to combined catalog
        if refcat is not None:
            self.log.print('Add ref.cat data to combined catalog')
            self.addRefcatData(combined, refcat)

        self.addMgLimits(combined)

        # save combined catalog
        self.log.print('Save combined catalog')
        self.save(combined)


if __name__ == '__main__':

    class MainApp:

        color = None
        folder = None
        showGraphs = False
        saveGraphs = False

        def __init__(self, argv):
            self.argv = argv
            pm.printInfo('pmmerge, v1.2')

        def processCommands(self):
            try:
                optlist, args = getopt(argv[1:], "c:", ['color=', 'show-graph', 'save-graph'])
            except GetoptError:
                pm.printError('Invalid command line options')
                return

            for o, a in optlist:
                if a[:1] == ':':
                    a = a[1:]
                elif o == '-c' or o == '--color':
                    self.color = a
                elif o == '--show-graph':
                    self.showGraphs = True
                elif o == '--save-graph':
                    self.saveGraphs = True

            self.folder = args[0]
            if not self.folder.endswith('/'):
                self.folder += '/'

        def run(self):
            self.processCommands()

            colors = ['Ri', 'Gi', 'Bi'] if self.color is None else [self.color]

            matcher = CatalogMatcher('Combined', self.folder, colors, self.showGraphs, self.saveGraphs,
                                     pm.Logger.LOG_MODE_DEBUG)
            matcher.process()


    app = MainApp(argv)
    app.run()

# end main.
