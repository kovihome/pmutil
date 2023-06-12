#!/usr/bin/env python3
#
# PmUtils/pmmerge
#

'''
Created on Mar 17, 2021

@author: kovi
'''

from getopt import getopt, GetoptError
from sys import argv
from astropy.table import Table, join, setdiff
from astropy.io import fits
from os import remove
from os.path import exists
from math import sqrt, log10
from numpy import arange

import pmbase as pm
from pmviz import VizUCAC4


class CatalogMatcher:

    folder = './'
    baseName = None
    colors = ['Ri','Gi','Bi']
    showGraphs = False
    saveGraphs = False

    # TODO: ezt a kettot konfigba tenni
    HMG_MAX_ERR = 0.15 # mg

    LIMIT_THRE = 0.5 # mg

    IMAGE_BORDER_SIZE = 10  # pixels

    MG_LIMIT = 17.0  # TODO: nem lehet-e ezt szamolni?

    COORD_MATCH_THRE = 2.0/60.0 # degrees

    viz = None

    def __init__(self, baseName, folder, colors, show, save):
        self.baseName = baseName
        if folder != None:
            self.folder = folder
        if colors != None:
            self.colors = colors
        self.viz =  VizUCAC4(self.MG_LIMIT)
        self.showGraphs = show
        self.saveGraphs = save

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

        m1,b1 = pm.linfit(logm1, loge1)
        hmg1 = -1.0 * 10.0 ** ((log10(self.HMG_MAX_ERR) - b1) / m1);

        m2,b2 = pm.linfit(logm2, loge2)
        hmg2 = -1.0 * 10.0 ** ((log10(self.HMG_MAX_ERR) - b2) / m2);

        plotcolor = 'gbr'[nr]
        self.hmgPlot.add(logm1, loge1, [m1, b1], f'log mi ({plotcolor.upper()}i)', f'log ei ({plotcolor.upper()}i)', plotcolor, pm.Plot.INV_X)

        return [hmg1, hmg2]

    def calcMgLimit(self, fileName, nr):
        t = Table.read(fileName, format='ascii.sextractor')
        mgs = self.hmg(t, nr)
        return mgs[0]

    def readFrameSize(self, fitsFileName):
        h = pm.getFitsHeaders(fitsFileName, ['NAXIS1', 'NAXIS2'])
        xf = int(h['NAXIS1'])
        yf = int(h['NAXIS2'])
        return (xf, yf)

    def determineOnFrameStatus(self, x_obj, y_obj, im_w, im_h):
        if x_obj < 0 or x_obj >= im_w or y_obj < 0 or y_obj >= im_h:
            return 'O'
        if x_obj < self.IMAGE_BORDER_SIZE or x_obj >= im_w - self.IMAGE_BORDER_SIZE or y_obj < self.IMAGE_BORDER_SIZE or y_obj >= im_h - self.IMAGE_BORDER_SIZE:
            return 'B'
        return '-'

    def invokeGrMatch(self, refFile, matchFile, outFile):
        matchOptions = '--match-coord --col-ref 14,15 --col-inp 14,15'
#        matchOptions = '--match-points --col-ref 16,17 --col-inp 16,17'
        return pm.invoke("grmatch %s -r %s -i %s -o %s --output-excluded-reference %s --output-excluded-input %s" % \
                      (matchOptions, refFile, matchFile, outFile, outFile + '.xrf', outFile + '.xin'))

    def calcVisFlag(self, mg, limit):
        '''
        flags I - INVISIBLE, B - BELOW_LIMIT, N - EAR_LIMIT, S - SATURATED
        '''
        # TODO: saturated objects

        if mg > 0.0:
            return 'I'
        elif mg > limit:
            return 'B'
        elif mg > limit - self.LIMIT_THRE:
            return 'N'
        return '-'

    def calcPosFlag(self, x_obj, y_obj, im_w, im_h):
        '''
        x_obj # object x position
        y_obj # object y position        
        im_w  # image width
        im_h  # image height

        flags B - object is near of image border
        '''

        return self.determineOnFrameStatus(x_obj, y_obj, im_w, im_h)


    def matchCatalogs(self):
        catFileBase = f"{self.folder}/{pm.setup['PHOT_FOLDER_NAME']}/{self.baseName}-"
        tempFileBase = f"{self.folder}/Temp/{self.baseName}-"
        pm.assureFolder(f"{self.folder}/Temp")

        # query image size
        img_w, img_h = self.readFrameSize(catFileBase + 'Gi.ast.fits')
        print(f'Image size: {img_w} x {img_h}')

        # calculate mg limits
        self.hmgPlot = pm.Plot(3, self.showGraphs, self.saveGraphs)

        hmgGi = self.calcMgLimit(catFileBase + 'Gi.cat', 0)
        hmgBi = self.calcMgLimit(catFileBase + 'Bi.cat', 1)
        hmgRi = self.calcMgLimit(catFileBase + 'Ri.cat', 2)
        self.mgLimits = { 'Gi': hmgGi, 'Bi': hmgBi, 'Ri': hmgRi }
        print(f'Limit mg Gi: {hmgGi}, Bi: {hmgBi}, Ri: {hmgRi}')
  
        
        self.hmgPlot.showOrSave(self.folder + '/magnitude_limit.png')

        # TODO: Azok a csillagok, amik nem latszanak Gi-ben, csak Bi-ben vagy Ri-ben, azok is keruljenek bele

        # match Bi cat to Gi
        self.invokeGrMatch(catFileBase + 'Gi.cat', catFileBase + 'Bi.cat', tempFileBase + 'GiBi.cat')
        
        # match Ri cat to combined Gi + Bi
        self.invokeGrMatch(tempFileBase + 'GiBi.cat', catFileBase + 'Ri.cat', tempFileBase + 'GiBiRi.cat')

        # load combined cat

        table = Table.read(tempFileBase + 'GiBiRi.cat', format='ascii.no_header', delimiter=' ')
        count = len(table)
        print(f'{count} matched objects found')

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
            posFlags.append(self.calcPosFlag(float(table[r.index]['col16']), float(table[r.index]['col17']), img_w, img_h))
            visFlags.append(self.calcVisFlag(r['MAG_GI'], hmgGi) + self.calcVisFlag(r['MAG_BI'], hmgBi) + self.calcVisFlag(r['MAG_RI'], hmgRi))

        cmb['AUID'] = auid
        cmb['RA'] = ras
        cmb['DEC'] = des
        cmb['POS_FLAG'] = posFlags
        cmb['VIZ_FLAG'] = visFlags

        return cmb


    def createCombinedCat(self, count=0):
        fnames=['AUID','VIZ_ID','ROLE','LABEL',
               'VIZ_FLAG','POS_FLAG','MATCH_FLAG',
               'RA','RA_DEG','DEC','DEC_DEG',
               'MAG_GI','ERR_GI','MAG_BI','ERR_BI','MAG_RI','ERR_RI',
               'MAG_B','ERR_B','MAG_V','ERR_V','MAG_R','ERR_R']
        fdtype=['U11','U16','U2','U50',
               'U1','U1','U1',
               'U12','U12','U12','U12',
               'U6','U5','U6','U5','U6','U5',
               'U6','U5','U6','U5','U6','U5']
        if count > 0:
            rows = [['AUID']*count,['-']*count,['F']*count,['-']*count,
               ['-']*count,['-']*count,['T']*count,
               ['']*count,[0.0]*count,['']*count,[0.0]*count,
               ['99.0']*count,['-']*count,['99.0']*count,['-']*count,['99.0']*count,['-']*count,
               ['99.0']*count,['-']*count,['99.0']*count,['-']*count,['99.0']*count,['-']*count]
            return Table(rows, names=fnames, dtype=fdtype)
        else:
            return Table(names=fnames, dtype=fdtype)


    def xmatchViz(self, combined):
        # xmatch cmb to UCAC4 catalog
        # TODO: err values from xmatch
        xmt = self.viz.xmatch(combined, 'RA_DEG', 'DEC_DEG')
        if not xmt:
            printError("Accessing Vizier service if failed.")
            return
        print(f'Xmatch table contains {len(xmt)} records')

        combined.add_index('AUID')

        for xr in xmt:
            auid = xr['AUID']
            r = combined.loc[auid]
            r['VIZ_ID'] = 'UCAC4-' + xr['UCAC4']
            r['MATCH_FLAG'] = '-'
            r['MAG_B'] = xr['Bmag']
            r['MAG_V'] = xr['Vmag']
            r['MAG_R'] = xr['Rmag']
#            r['ERR_B'] = xr['e_Bmag']
#            r['ERR_V'] = xr['e_Vmag']
#            r['ERR_R'] = xr['e_Rmag']

    def loadRefcat(self):
        refcatFileName = f"{self.folder}/ref.cat"
        if exists(refcatFileName):
            return Table.read(refcatFileName, format='ascii')
        else:
            pm.printWarning(f'No ref.cat in the folder {self.folder}')
            return None

    def addFrameCoords(self, cat, fitsFileName):
        # 1. convert cat to fits format
        baseFolder = fitsFileName.partition(pm.setup['PHOT_FOLDER_NAME'])[0]
        tempFile = self.folder + '/Temp/tempCat.fits'
        if not exists(tempFile):
            cat.write(tempFile)

        # 2. calculate ref objects' frame xy points
        wcsFile = fitsFileName.replace('.ast.fits', '.wcs')
        axyFile = fitsFileName.replace('.ast.fits', '.ref.axy')
        pm.invoke("wcs-rd2xy -w %s -i %s -o %s -R RA_DEG -D DEC_DEG" % (wcsFile, tempFile, axyFile))  # -f option need argument

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

        foMask = cat['FRAME_STATUS'] == 'O'
        foCat = cat[foMask]
#        print('Out of frame objects:')
#        print(foCat)

        foMask = cat['FRAME_STATUS'] != 'O'
        foCat = cat[foMask]
        return foCat

    def matchCatalogByCoords(self, cat, ra, dec):
        dmin = 99.0
        auid = None
        for row in cat:
            if row['VIZ_ID'] == '-':
                cra = float(row['RA_DEG'])
                cdec = float(row['DEC_DEG'])
                d = pm.quad(cra-ra, cdec-dec)
                if d < dmin:
                    dmin = d
                    auid = row['AUID']
#        print(f"Coord matched - auid: {auid}, d: {sqrt(dmin)*60}'")
        return auid, dmin

    def updateCmb(self, cmb, r, auid = None, vizid = None):
        if vizid:
            cr = cmb.loc['VIZ_ID', vizid]
        else:
            cr = cmb.loc['AUID', auid]
        cr['ROLE']  = r['ROLE']
        cr['AUID']  = r['AUID']
        cr['LABEL'] = r['LABEL']
        cr['MAG_B'] = r['MAG_B']
        cr['MAG_V'] = r['MAG_V']
        cr['MAG_R'] = r['MAG_R']
        cr['ERR_B'] = r['ERR_B']
        cr['ERR_V'] = r['ERR_V']
        cr['ERR_R'] = r['ERR_R']
#        r['MATCH_FLAG'] = xr['MATCH_FLAG']

    def insertCmb(self, cmb, r):
        n = [  r['AUID'], '-', r['ROLE'], r['LABEL'], 'NNN', '-', 'N', r['RA'], r['RA_DEG'], r['DEC'], r['DEC_DEG'], \
                '-', '-', '-', '-', '-', '-', r['MAG_B'], r['ERR_B'], r['MAG_V'], r['ERR_V'], r['MAG_R'], r['ERR_R'] ]
        cmb.add_row(n)
       

    def addRefcatData(self, cmb, refcat):

        del refcat.meta['comments']

        fitsFileName = f"{self.folder}/{pm.setup['PHOT_FOLDER_NAME']}/{self.baseName}-Gi.ast.fits"
        foRefcat = self.filterOutFrameStars(refcat, fitsFileName)

        # xmatch refcat to UCAC4
        # TODO: drop out records, where angDist comlumn value is grater then a threshold (bad match)
        xmt = self.viz.xmatch(foRefcat, 'RA_DEG', 'DEC_DEG')
        print(f'Refcat xmatch table contains {len(xmt)} records')

#        print(xmt)

        missingRefcat = setdiff(foRefcat, xmt, 'AUID')
#        print(missingRefcat)

        cmb.add_index('VIZ_ID')

        for xr in xmt:
            vizId = 'UCAC4-' + xr['UCAC4']
            if vizId in cmb['VIZ_ID']:
                try:
                    self.updateCmb(cmb, xr, vizid=vizId)

                except KeyError:
                    print(r)
                    print(xr)

            else:
                print('xmt record not matched in cmb:')
                print(xr)

        mergedTable = join(xmt, cmb, keys='AUID', join_type='left')
        mask = mergedTable['VIZ_ID'] == '-' # or mergedTable['VIZ_ID'] == '' or mergedTable['VIZ_ID'] == None
        nt = mergedTable[mask]

        print(f'Negative transients: {len(nt)} found.')

#refCat: AUID         ROLE RA             RA_DEG         DEC            DEC_DEG        MAG_B     ERR_B     MAG_V     ERR_V     MAG_R     ERR_R     LABEL
#cmb:    'AUID','VIZ_ID','ROLE','LABEL','VIZ_FLAG','POS_FLAG','MATCH_FLAG','RA','RA_DEG','DEC','DEC_DEG',
#        'MAG_GI','ERR_GI','MAG_BI','ERR_BI','MAG_RI','ERR_RI','MAG_B','ERR_B','MAG_V','ERR_V','MAG_R','ERR_R']

        for r in nt:
            self.insertCmb(cmb, r)

        for r in missingRefcat:
             auid, d = self.matchCatalogByCoords(cmb, float(r['RA_DEG']), float(r['DEC_DEG']))
             print(f"Match missingRefcat object {r['AUID']} {r['LABEL']} for {auid} within {d*60}'")
             if d < 2.0/60.0:
                 self.updateCmb(cmb, r, auid=auid)
             else:
                 self.insertCmb(cmb, r)

    def addMgLimits(self, cmb):
        for cc in self.mgLimits.keys():
            pm.addTableComment(cmb, 'MgLimitInst' + cc, '%7.3f' % (self.mgLimits[cc]))

    def save(self, combined):
        combinedFileName = f"{self.folder}/{pm.setup['PHOT_FOLDER_NAME']}/{self.baseName}.cmb"
        combined.write(combinedFileName, format = 'ascii.fixed_width', delimiter = ' ', overwrite=True)

    def process(self):

        # match cat files together
        print('Combine measurement catalogs all together')
        combined = self.matchCatalogs()
        
        # xmatch combined catalog to UCAC4, and store UCAC4 id
        print('Match combined catalog with UCAC4')
        self.xmatchViz(combined)

        # load refcat
        refcat = self.loadRefcat()

        # add refcat data to combined catalog
        if refcat != None:
            print('Add ref.cat data to combined catalog')
            self.addRefcatData(combined, refcat)

        self.addMgLimits(combined)

        # save combined catalog
        print('Save combined catalog')
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
                optlist, args = getopt (argv[1:], "c:", ['color=', 'show-graph', 'save-graph'])
            except GetoptError:
                print ('Invalid command line options')
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

            colors = ['Ri','Gi','Bi'] if self.color == None else [ self.color ]

            matcher = CatalogMatcher('Combined', self.folder, colors, self.showGraphs, self.saveGraphs)
            matcher.process()


    app = MainApp(argv)
    app.run()

# end main.
