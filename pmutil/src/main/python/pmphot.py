#
# PmUtils/pmphot
#
'''
Created on Dec 28, 2019

@author: kovi
'''
from numpy import median, sqrt, arange  # , polyfit
from getopt import getopt, GetoptError
from sys import argv  # , stdout
# from glob import glob
from astropy.table import Table, Column
from os.path import exists

import matplotlib.pyplot as plt

from pmbase import jd, getFitsHeaders, printError, loadPplSetup, quad


def loadCatalog(refFileName):
    '''
        Load a catalog with one line header
    '''

    refCat = { 'header': [], 'cat': [] }

    ref = open(refFileName, 'r')
    if ref:
        headerLine = ref.readline()
        refCat['header'] = headerLine.split()

        for refLine in ref:
            refCatLine = refLine.split()
            refCat['cat'].append(refCatLine)
        ref.close()

        print (("Photometric result catalog: %s") % refFileName)
        return refCat

    else:

        print (("Reference catalog %s not found") % refFileName)
        return None


def getHeaderPos(cat, headerName):
    index = 0
    for header in cat['header']:
        if header == headerName:
            return index
        index = index + 1
    return None


def addHeader(cat, headerName):
    pos = getHeaderPos(cat, headerName)
    if pos != None:
        return pos

    pos = len(cat['header'])
    cat['header'].append(headerName)
    return pos


fieldMgInstrumental = 'MAG_BEST'
fieldMgTrue = 'MAG_V'
fieldMgErrInstrumental = 'MAGERR_BEST'
fieldMgErrTrue = 'ERR_V'
fieldAuid = 'AUID'
fieldRole = 'ROLE'

MAG_ERR_DEFAULT = 0.2


#def getDateObs(fileName):
#    '''
#    fileName: path/[Seq_nnn|Combined]_c.fits.cat.cat
#    '''
#    astFileName = fileName.split('.')[0]
#    astFileName = astFileName + '.ast.fits'
#
#    return getFitsHeader(astFileName, 'DATE-OBS')


class StdCoeffs:

    coeffs = None  # astropy.table.Table of std coeffs
    fileName = None  # coeff fiel name
    observer = None  # observer name
    date = None  # date of the current image
    camera = None  # camera using to crate current image
    telescope = None  # telescope using to create current image

    def __init__(self, coeffFileName, observer, date, camera, telescope):
        self.fileName = coeffFileName
        self.observer = observer
        self.date = date
        self.camera = camera
        self.telescope = telescope

    def open(self):
        if not exists(self.fileName):
            self.coeffs = Table()
            self.coeffs.add_column(Column(name = 'OBSERVER', dtype = 'U5', format = '%-10s'))
            self.coeffs.add_column(Column(name = 'DATE', dtype = 'U10', format = '%-14s'))
            self.coeffs.add_column(Column(name = 'TV', dtype = 'f4', format = '%11.4f'))
            self.coeffs.add_column(Column(name = 'ERR_TV', dtype = 'f4', format = '%11.4f'))
            self.coeffs.add_column(Column(name = 'TVR', dtype = 'f4', format = '%11.4f'))
            self.coeffs.add_column(Column(name = 'ERR_TVR', dtype = 'f4', format = '%11.4f'))
            self.coeffs.add_column(Column(name = 'TBV', dtype = 'f4', format = '%11.4f'))
            self.coeffs.add_column(Column(name = 'ERR_TBV', dtype = 'f4', format = '%11.4f'))
            self.coeffs.add_column(Column(name = 'CAMERA', dtype = 'U24', format = '%-24s'))
            self.coeffs.add_column(Column(name = 'TELESCOPE', dtype = 'U24', format = '%-24s'))
            self.coeffs.add_column(Column(name = 'FIELD', dtype = 'U24', format = '%-24s'))
        else:
            self.coeffs = Table.read(self.fileName, format = 'ascii')
            self.coeffs['OBSERVER'].format = '%-10s'
            self.coeffs['DATE'].format = '%-14s'
            self.coeffs['TV'].format = '%11.4f'
            self.coeffs['ERR_TV'].format = '%11.4f'
            self.coeffs['TVR'].format = '%11.4f'
            self.coeffs['ERR_TVR'].format = '%11.4f'
            self.coeffs['TBV'].format = '%11.4f'
            self.coeffs['ERR_TBV'].format = '%11.4f'
            self.coeffs['CAMERA'].format = '%-24s'
            self.coeffs['TELESCOPE'].format = '%-24s'
            self.coeffs['FIELD'].format = '%-24s'

    def addCoeffs(self, Tv, Tvr, Tbv, errTv = None, errTvr = None, errTbv = None, field = None):
        lastCoeffs = self.getCoeffs()
        if not lastCoeffs:
            self.coeffs.add_row((self.observer, self.date, Tv, errTv, Tvr, errTvr, Tbv, errTbv, self.camera, self.telescope, field))
        else:
            lastCoeffs['TV'] = Tv
            lastCoeffs['ERR_TV'] = errTv
            lastCoeffs['TVR'] = Tvr
            lastCoeffs['ERR_TVR'] = errTvr
            lastCoeffs['TBV'] = Tbv
            lastCoeffs['ERR_TBV'] = errTbv
            lastCoeffs['FIELD'] = field

    def getCoeffs(self):
        flt = (self.coeffs['OBSERVER'] == self.observer) & (self.coeffs['DATE'] == self.date) & (self.coeffs['CAMERA'] == self.camera) & (self.coeffs['TELESCOPE'] == self.telescope)
        rows = self.coeffs[flt]
        if len(rows) == 0:
            return None
        else:
            return rows[0]

    def getBestCoeffs(self):
        # TODO:
        flt = (self.coeffs['OBSERVER'] == self.observer) & (self.coeffs['CAMERA'] == self.camera) & (self.coeffs['TELESCOPE'] == self.telescope)
        rows = self.coeffs[flt]
        if len(rows) == 0:
            return None
        elif len(rows) == 1:
            return rows[0]
        else:
            jdnow = jd(self.date)
            delta = abs(jdnow - jd(rows[0]['DATE']))
            ix = 0
            for j in range(len(rows) - 1):
                d = abs(jdnow - jd(rows[j + 1]['DATE']))
                if d < delta:
                    delta = d
                    ix = j + 1
            return rows[ix]

    def getAvgCoeffs(self):
        # TODO:
        return None

    def save(self):
#        self.coeffs.write(stdout, format='ascii.fixed_width_two_line', delimiter=' ')
        self.coeffs.write(self.fileName, format = 'ascii.fixed_width_two_line', delimiter = ' ')


class Photometry:

    opt = {}  # command line options
    ppl = {}  # PPL setup from ppl-setup config file
    pos = None  # catalog header positions
    imageProps = {}  # image properties from FITS header
    coeffTable = None
    fits = {}

    def __init__(self, opt, ppl):
        self.opt = opt
        self.ppl = ppl

    def loadFitsHeaders(self, fileName):
        '''
        fileName: path/[Seq_nnn|Combined]_c.fits.cat.cat
        '''
        astFileName = fileName.split('.')[0]
        astFileName = astFileName + '.ast.fits'
        self.fits = getFitsHeaders(astFileName, ['DATE-OBS', 'INSTRUME', 'TELESCOP'])

    def getCamera(self):
        if 'INSTRUME' in self.fits.keys() and self.fits['INSTRUME'] and len(self.fits['INSTRUME']) > 0:
            return self.fits['INSTRUME']
        if 'camera' in self.opt and self.opt['camera'] and len(self.opt['camera']) > 0:
            return self.opt['camera']
        if 'DEF_CAMERA' in self.ppl.keys() and self.ppl['DEF_CAMERA'] and len(self.ppl['DEF_CAMERA']) > 0:
            return self.ppl['DEF_CAMERA']
        return 'Generic Camera'

    def getTelescope(self):
        if 'TELESCOP' in self.fits.keys() and self.fits['TELESCOP'] and len(self.fits['TELESCOP']) > 0:
            return self.fits['TELESCOP']
        if 'telescope' in self.opt and self.opt['telescope'] and len(self.opt['telescope']) > 0:
            return self.opt['telescope']
        if 'DEF_TELESCOPE' in self.ppl.keys() and self.ppl['DEF_TELESCOPE'] and len(self.ppl['DEF_TELESCOPE']) > 0:
            return self.ppl['DEF_TELESCOPE']
        return 'Generic Telescope'

    def transformMgs(self, refCat, stdcolor, p, err):
        mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
        miPos = getHeaderPos(refCat, fieldMgInstrumental)
        evPos = getHeaderPos(refCat, 'ERR_' + stdcolor)
        eiPos = getHeaderPos(refCat, fieldMgErrInstrumental)
        rolePos = getHeaderPos(refCat, fieldRole)

        result = []
        for pm in refCat['cat']:
            if pm[rolePos] == 'V' or pm[rolePos] == 'VF' or pm[rolePos] == 'K':
                mv = p[0] * float(pm[miPos]) + p[1]
                pm[mvPos] = mv
                if err and pm[rolePos] != 'VF':
                    pm[evPos] = quad(float(pm[eiPos]), err)
                else:
                    pm[evPos] = 0.0
                result.append(pm)
        return result

    def stdColor(self, color):
        stdcolor = color[:1].upper()
        if stdcolor == 'G':
            stdcolor = 'V'
        return stdcolor

    def calculateMgsRobustAveraging(self, refCat, color, bestCompId):
        '''
        Calculate unit slope linear fit regression with robust averaging
        see: http://gcx.sourceforge.net/html/node11.html
        inputs:    I[k] - instrumental mgs              MAG_ISOCORR | MAG_BEST
                   S[k] - standard mgs                  MAG_R \ MAG_V | MAG_B
                   ei[k] - error of instrumental mgs    MAGERR_ISOCORR | MAGERR_BEST
                   es[k] - error of standard mgs        ERR_R | ERR_V | ERR_B
        '''

        stdcolor = self.stdColor(color)

        mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
        miPos = getHeaderPos(refCat, fieldMgInstrumental)
        evPos = getHeaderPos(refCat, 'ERR_' + stdcolor)
        eiPos = getHeaderPos(refCat, fieldMgErrInstrumental)
        rolePos = getHeaderPos(refCat, fieldRole)

        bestComp = self.findRow(refCat['cat'], bestCompId)
        bestComp[rolePos] = 'K'

        y = []
        e2 = []
        w = []
        mis = []
        mvs = []
        for pm in refCat['cat']:
            if pm[rolePos] == 'C' and pm[mvPos] != '-':
                mi = float(pm[miPos])
                mv = float(pm[mvPos])
                mis.append(mi)
                mvs.append(mv)
                ei = float(pm[eiPos])
                if pm[evPos][0].isdigit():
                    ev = float(pm[evPos])
                else:
                    ev = MAG_ERR_DEFAULT  # TODO: what to do, if no mg error value?
                y.append(mv - mi)
                ek2 = ei * ei + ev * ev
                e2.append(ek2)
                w.append(1.0 / ek2)

        z = median(y)
        r = [0.0] * len(y)
        r_ = [0.0] * len(y)
        w_ = [1.0] * len(y)

        # iteration starts
        for it in range(10):

            zsum = 0.0
            wsum = 0.0
            beta = 1.0
            alpha = 1.0
            for k in range(len(y)):
                r[k] = y[k] - z
                w_k = w[k] / (1 + (r[k] / alpha) ** beta)
                w_[k] = w_k
                zsum = zsum + y[k] * w_k
                wsum = wsum + w_k
            z = zsum / wsum
        # iteration ends

        su = 0.0
        sl = 0.0
        se2 = 0.0
        for k in range(len(y)):
            su = su + r[k] * r[k] * w_[k]
            sl = sl + w_[k]
            se2 = se2 + 1.0 / (w_[k] * w_[k])
        ez_2 = su / sl
        mel_2 = su / float(len(y) - 1)
        err = sqrt(se2 / len(y))
        print ("(gcx) color: %s, zp: %7.4f, mel^2: %7.4f, ez^2: %7.4f, ez: %7.4f, N: %d, ev: %7.4f" % (color, z, mel_2, ez_2, sqrt(ez_2), len(y), err))

        if self.opt['showGraphs']:
            plotcolor = color[0].lower()
            plt.subplot(311 + self.opt['color'].index(color))
            xr = arange(min(mis), max(mis), 0.01)
            plt.plot(mis, mvs, plotcolor + 'o', xr, xr + z, 'k')
            plt.xlabel(color)
            plt.ylabel(stdcolor)

        # apply result to variables
        result = self.transformMgs(refCat, stdcolor, [1.0, z], err)
        return result

    def calculateMgsLinearFit(self, refCat, color, bestCompId):
        '''
        Calculate magnitudes with naive linear fit
        '''
        # TODO: calculate comp mag error

        stdcolor = self.stdColor(color)

        mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
        miPos = getHeaderPos(refCat, fieldMgInstrumental)
        evPos = getHeaderPos(refCat, 'ERR_' + stdcolor)
        eiPos = getHeaderPos(refCat, fieldMgErrInstrumental)
        rolePos = getHeaderPos(refCat, fieldRole)

        bestComp = self.findRow(refCat['cat'], bestCompId)
        bestComp[rolePos] = 'K'

        # TODO
        mi = []
        mv = []
        er = []
#    p = [1.0, 0.0]
        ep = 0.0
        N = 0
        for pm in refCat['cat']:
            if pm[rolePos] == 'C' and pm[mvPos] != '-':
                mi.append(float(pm[miPos]))
                mv.append(float(pm[mvPos]))
                ei = float(pm[eiPos])
                ev = float(pm[evPos]) if pm[evPos] != '-' else MAG_ERR_DEFAULT
                e2 = ei * ei + ev * ev
                ep = ep + e2
                N = N + 1
                er.append(1.0 / e2)

#        p = Polynomial.fit(mi, mv, 1, w = er)
        coef = self.linfit(mi, mv, er)
        ep = sqrt(ep / float(N))

        print ('polyfit result:', coef, 'error:', ep)

        if self.opt['showGraphs']:
            plotcolor = color[0].lower()
            plt.subplot(311 + self.opt['color'].index(color))
            xr = arange(min(mi), max(mi), 0.01)
            plt.plot(mi, mv, plotcolor + 'o', xr, coef[0] * xr + coef[1], 'k')
            plt.xlabel(color)
            plt.ylabel(stdcolor)

        # apply result to variables
        result = self.transformMgs(refCat, stdcolor, coef, ep)
        return result

    def findBestCompStar(self, allCatalogs, color):

        refCat = allCatalogs[color]
        stdcolor = self.stdColor(color)

        idPos = getHeaderPos(refCat, fieldAuid)
        eiPos = getHeaderPos(refCat, fieldMgErrInstrumental)
        rolePos = getHeaderPos(refCat, fieldRole)

        mvPosList = {}
        evPosList = {}
        for c in self.opt['color']:
            mvPosList[c] = getHeaderPos(refCat, 'MAG_' + self.stdColor(c))
            evPosList[c] = getHeaderPos(refCat, 'ERR_' + self.stdColor(c))

        emin = 1.0
        bestComp = None
        for pm in refCat['cat']:
            mvs = pm[mvPosList[color]]
            evs = pm[evPosList[color]]
            if pm[rolePos] == 'C' and mvs != '-' and evs != '-' and evs != '0.0':
                ei = float(pm[eiPos]) if pm[eiPos] != '-' else MAG_ERR_DEFAULT
                ev = float(evs)
                e = ei * ei + ev * ev
                if e < emin:
                    haveAllMags = True
                    for c in self.opt['color']:
                       if c != color:
                           if pm[mvPosList[c]] == '-' or pm[evPosList[c]] == '-':
                               haveAllMags = False
                               break
                    if not haveAllMags:
                        continue
                    emin = e
                    bestComp = pm

        if bestComp:
            ei = float(bestComp[eiPos]) if bestComp[eiPos] != '-' else MAG_ERR_DEFAULT
            ev = float(bestComp[evPosList[color]]) if bestComp[evPosList[color]] != '-' else MAG_ERR_DEFAULT
            e = quad(ei, ev)
            bestComp[rolePos] = 'K'
            print('Best comp star: %s, mag: %s, err: %4.3f' % (bestComp[idPos], bestComp[mvPosList[color]], e))
            return bestComp[idPos]
        else:
            printError('No usable comp star found ; check the comp stars if they exist in all colors you need')
            return None

    def calculateMgsComparision(self, refCat, color, bestCompId):
        stdcolor = self.stdColor(color)

#        idPos = getHeaderPos(refCat, fieldAuid)
        mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
        miPos = getHeaderPos(refCat, fieldMgInstrumental)
        evPos = getHeaderPos(refCat, 'ERR_' + stdcolor)
        eiPos = getHeaderPos(refCat, fieldMgErrInstrumental)
        rolePos = getHeaderPos(refCat, fieldRole)

        # find best comp star
        bestComp = self.findRow(refCat['cat'], bestCompId)
        bestComp[rolePos] = 'K'

        zp = float(bestComp[mvPos]) - float(bestComp[miPos])
        p = [ 1.0, zp ]
        ei = float(bestComp[eiPos]) if bestComp[eiPos] != '-' else MAG_ERR_DEFAULT
        ev = float(bestComp[evPos]) if bestComp[evPos] != '-' else MAG_ERR_DEFAULT
        ep = quad(ei, ev)

        if self.opt['showGraphs']:
            mi = []
            mv = []
            for pm in refCat['cat']:
                if pm[rolePos] == 'C' and pm[mvPos] != '-':
                    mi.append(float(pm[miPos]))
                    mv.append(float(pm[mvPos]))

            plotcolor = color[0].lower()
            plt.subplot(311 + self.opt['color'].index(color))
            xr = arange(min(mi), max(mi), 0.01)
            plt.plot(mi, mv, plotcolor + 'o', xr, xr + zp, 'k')
            plt.xlabel(color)
            plt.ylabel(stdcolor)

        # apply result to variables
        result = self.transformMgs(refCat, stdcolor, p, ep)
        return result

    def getHeaderPositions(self, refCat):
        pos = {}
        pos[fieldAuid] = getHeaderPos(refCat, fieldAuid)
        pos['LABEL'] = getHeaderPos(refCat, 'LABEL')
        pos['NUMBER'] = getHeaderPos(refCat, 'NUMBER')
        pos['MAG_V'] = getHeaderPos(refCat, 'MAG_V')
        pos['MAG_R'] = getHeaderPos(refCat, 'MAG_R')
        pos['MAG_B'] = getHeaderPos(refCat, 'MAG_B')
        pos['ERR_V'] = getHeaderPos(refCat, 'ERR_V')
        pos['ERR_R'] = getHeaderPos(refCat, 'ERR_R')
        pos['ERR_B'] = getHeaderPos(refCat, 'ERR_B')
        pos[fieldMgInstrumental] = getHeaderPos(refCat, fieldMgInstrumental)
        pos[fieldMgErrInstrumental] = getHeaderPos(refCat, fieldMgErrInstrumental)
        pos['RA'] = getHeaderPos(refCat, 'RA')
        pos['DEC'] = getHeaderPos(refCat, 'DEC')
        pos[fieldRole] = getHeaderPos(refCat, fieldRole)
        return pos

    def reportResult(self, result, refCat, outFileName, color, dateObs):

        stdcolor = self.stdColor(color)

        idPos = getHeaderPos(refCat, fieldAuid)
        labelPos = getHeaderPos(refCat, 'LABEL')
        astIdPos = getHeaderPos(refCat, 'NUMBER')
        mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
        evPos = getHeaderPos(refCat, 'ERR_' + stdcolor)
        raPos = getHeaderPos(refCat, 'RA')
        decPos = getHeaderPos(refCat, 'DEC')
        rolePos = getHeaderPos(refCat, fieldRole)

        if outFileName == None:
            outFileName = 'result.pm'
        r = open(outFileName, 'w')

        # write headers
        r.write('AUID'.ljust(15))
        r.write('LABEL'.ljust(30))
        r.write('ASTID'.ljust(10))
        # ref cat: RA_DEG, DEC_DEG
        # pm cat: ALPHA_J2000, DELTA_J2000
        r.write('RA'.ljust(20))
        r.write('DEC'.ljust(20))
        r.write('DATE-OBS'.ljust(20))
        r.write('JD'.ljust(20))
        r.write(('FLAG').ljust(5))
        r.write(('MAG').ljust(20))
        r.write(('COL').ljust(5))
        r.write(('ERR').ljust(20))
        r.write('\n')

        # write result data
        jdObs = jd(dateObs)
        jds = "%12.4f" % jdObs
        for res in result:
            r.write(res[idPos].ljust(15))
            r.write(res[labelPos].ljust(30))
            r.write(res[astIdPos].ljust(10))
            r.write(res[raPos].ljust(20))
            r.write(res[decPos].ljust(20))
            r.write(dateObs.ljust(20))
            r.write(jds.ljust(20))
            if res[rolePos] == 'VF':
                flag = 'F'
            elif res[rolePos] == 'K':
                flag = 'K'
            else:
                flag = '-'
            r.write(flag.ljust(5))
            r.write(("%6.3f" % (float(res[mvPos]))).ljust(20))
            c = stdcolor if self.opt['useCoeffs'] else color
            r.write(c.ljust(5))
            r.write(("%5.3f" % (float(res[evPos])) if res[evPos] != '-' and flag != 'F' else '-').ljust(20))
            r.write('\n')

        r.close()

    # standardization functions

    def findRow(self, cdata, rid):
        for row in cdata:
            if row[self.pos[fieldAuid]] == rid:
                return row
        return None

    def mergeCatalogs(self, allCatalogs):
        auidPos = self.pos[fieldAuid]
        rolePos = self.pos['ROLE']
        mgiPos = self.pos[fieldMgInstrumental]

        # collect ids
        allIds = set()
        for j in range(len(self.opt['color'])):
            color = self.opt['color'][j]
            stdcolor = self.stdColor(color)

            mgsPos = self.pos['MAG_' + stdcolor]

            result = allCatalogs[color]['cat']
            ids = set()
            for row in result:
                role = row[rolePos]
                mgi = row[mgiPos]
                mgs = row[mgsPos]
                auid = row[auidPos]
                if role == 'C' and mgs and mgs != '-' and mgi and mgi != '-':
                    ids.add(auid)

            if len(allIds) == 0:
                allIds = ids
            else:
                allIds = allIds.intersection(ids)

        # find records for all ids
        # merged catalog fields: ID, B, V, R, bi, gi, ri
        merged = []
        for rid in allIds:
            grow = self.findRow(allCatalogs['Gi']['cat'], rid)
            brow = self.findRow(allCatalogs['Bi']['cat'], rid)
            rrow = self.findRow(allCatalogs['Ri']['cat'], rid)

            mV = grow[self.pos['MAG_V']]
            mR = rrow[self.pos['MAG_R']]
            mB = brow[self.pos['MAG_B']]
            mgi = grow[self.pos[fieldMgInstrumental]]
            mbi = brow[self.pos[fieldMgInstrumental]]
            mri = rrow[self.pos[fieldMgInstrumental]]

            eB = brow[self.pos['ERR_B']]
            eV = grow[self.pos['ERR_V']]
            eR = rrow[self.pos['ERR_R']]
            egi = grow[self.pos[fieldMgErrInstrumental]]
            ebi = brow[self.pos[fieldMgErrInstrumental]]
            eri = rrow[self.pos[fieldMgErrInstrumental]]

            merged.append([ rid, mB, mV, mR, mbi, mgi, mri, eB, eV, eR, ebi, egi, eri ])

        return merged

    def linfit(self, x, y, w = None):
        N = len(x)
        X = 0.0
        Y = 0.0
        XY = 0.0
        X2 = 0.0
        M = 0.0
        for j in range(N):
            wj = w[j] if w else 1.0
            X += x[j] * wj
            Y += y[j] * wj
            XY += x[j] * y[j] * wj
            X2 += x[j] * x[j] * wj
            M += wj
        m = (Y * X - M * XY) / (X * X - M * X2)
        b = (Y - m * X) / M
        return m, b

    def calculateCoeffs(self, mergedCat):
        # Tv:  slope of (V-v) -> (V-R)
        # Tvr: 1/slope of (v-r) -> (V-R)
        # Tbv: 1/slope of (b-v) -> (B-V)
        V_vi = []
        V_R = []
        vi_ri = []
        bi_vi = []
        B_V = []
        w_Tv = []
        w_Tvr = []
        w_Tbv = []
        for row in mergedCat:
            V_vi.append(float(row[2]) - float(row[5]))
            V_R.append(float(row[2]) - float(row[3]))
            vi_ri.append(float(row[5]) - float(row[6]))
            bi_vi.append(float(row[4]) - float(row[5]))
            B_V.append(float(row[1]) - float(row[2]))

            for p in range(7, 12):
                if row[p] == '-':
                    row[p] = "0.1"
            e_V_vi = float(row[8]) + float(row[11])
            e_V_R = float(row[8]) + float(row[9])
            e_vi_ri = float(row[11]) + float(row[12])
            e_bi_vi = float(row[10]) + float(row[11])
            e_B_V = float(row[7]) + float(row[8])
            w_Tv.append(1.0 / (e_V_vi * e_V_vi + e_V_R * e_V_R))
            w_Tvr.append(1.0 / (e_vi_ri * e_vi_ri + e_V_R * e_V_R))
            w_Tbv.append(1.0 / (e_bi_vi * e_bi_vi + e_B_V * e_B_V))

        coef = self.linfit(V_R, V_vi, w_Tv)
        Tv = coef[0]
        Bv = coef[1]  # ##
        coef = self.linfit(V_R, vi_ri, w_Tvr)
        Tvr = 1.0 / coef[0]
        Bvr = coef[1]  # ##
        coef = self.linfit(B_V, bi_vi, w_Tbv)
        Tbv = 1.0 / coef[0]
        Bbv = coef[1]  # ##

        if self.opt['showGraphs']:
            plt.figure()
            plt.subplot(311)
            xr = arange(min(V_R), max(V_R), 0.01)
            plt.plot(V_R, V_vi, 'go', xr, Tv * xr + Bv, 'k')
            plt.xlabel('V-R')
            plt.ylabel('V-v')

            plt.subplot(312)
            plt.plot(V_R, vi_ri, 'ro', xr, xr / Tvr + Bvr, 'k')
            plt.xlabel('V-R')
            plt.ylabel('v-r')

            plt.subplot(313)
            xr = arange(min(B_V), max(B_V), 0.01)
            plt.plot(B_V, bi_vi, 'bo', xr, xr / Tbv + Bbv, 'k')
            plt.xlabel('B-V')
            plt.ylabel('b-v')
            plt.show()

        print ('standard coeffs: Tv =', Tv, 'Tvr =', Tvr, 'Tbv = ', Tbv)
        return [Tv, Tvr, Tbv]

    def openCoeffs(self, obsDate):
        if not self.coeffTable:
            configFolder = self.ppl['CONFIG_FOLDER'].rstrip('/')
            coeffFile = configFolder + '/stdcoeffs.cat'
            self.coeffTable = StdCoeffs(coeffFile, self.opt['observer'], obsDate, self.getCamera().replace(' ', '_'), self.getTelescope().replace(' ', '_'))
            self.coeffTable.open()

    def saveCoeffs(self, coeffs, obsDate):
        self.openCoeffs(obsDate)
        self.coeffTable.addCoeffs(coeffs[0], coeffs[1], coeffs[2], 0.0, 0.0, 0.0, 'SA104'.replace(' ', '_'))
#        self.coeffTable.save(overwrite=True)
        self.coeffTable.save()

    def loadCoeffs(self, obsDate):
        self.openCoeffs(obsDate)
        bestCoeffs = self.coeffTable.getBestCoeffs()
        if not bestCoeffs:
            printWarning ('standard coeffs not found')
            return None
        print ('standard coeffs: Tv =', bestCoeffs['TV'], 'Tvr =', bestCoeffs['TVR'], 'Tbv = ', bestCoeffs['TBV'])
        return [bestCoeffs['TV'], bestCoeffs['TVR'], bestCoeffs['TBV']]

    def calculateStdMgs(self, coeffs, allCatalogs, bestComp):
        '''
        allResults: dict of results, key is th color
        '''
        Tv = coeffs[0]
        Tvr = coeffs[1]
        Tbv = coeffs[2]
        cg_row = self.findRow(allCatalogs['Gi']['cat'], bestComp)
        cg_row[self.pos['ROLE']] = 'K'
        cb_row = self.findRow(allCatalogs['Bi']['cat'], bestComp)
        cb_row[self.pos['ROLE']] = 'K'
        cr_row = self.findRow(allCatalogs['Ri']['cat'], bestComp)
        cr_row[self.pos['ROLE']] = 'K'
        vc = float(cg_row[self.pos[fieldMgInstrumental]])
        bc = float(cb_row[self.pos[fieldMgInstrumental]])
        rc = float(cr_row[self.pos[fieldMgInstrumental]])
        Vc = float(cg_row[self.pos['MAG_V']])
        Bc = float(cb_row[self.pos['MAG_B']])
        Rc = float(cr_row[self.pos['MAG_R']])
        Vc_Rc = Vc - Rc

        err_vc = float(cg_row[self.pos[fieldMgErrInstrumental]])
        err_bc = float(cb_row[self.pos[fieldMgErrInstrumental]])
        err_rc = float(cr_row[self.pos[fieldMgErrInstrumental]])
        err_Vc = float(cg_row[self.pos['ERR_V']])
        err_Bc = float(cb_row[self.pos['ERR_B']])
        err_Rc = float(cr_row[self.pos['ERR_R']])

        allResults = { 'Gi': [], 'Bi': [], 'Ri': [] }
        for g_row in allCatalogs['Gi']['cat']:
            role = g_row[self.pos['ROLE']]
            if role == 'V' or role == 'VF' or role == 'K':
                auid = g_row[self.pos['AUID']]
                b_row = self.findRow(allCatalogs['Bi']['cat'], auid)
                r_row = self.findRow(allCatalogs['Ri']['cat'], auid)

                v0 = float(g_row[self.pos[fieldMgInstrumental]])
                b0 = float(b_row[self.pos[fieldMgInstrumental]])
                r0 = float(r_row[self.pos[fieldMgInstrumental]])

                V_R = (Vc_Rc) + Tvr * ((v0 - r0) - (vc - rc))
                V = v0 + (Vc - vc) + Tv * ((V_R) - (Vc_Rc))
                R = V - (V_R)
                B = V + (Bc - Vc) + Tbv * ((b0 - v0) - (bc - vc))

                g_row[self.pos['MAG_V']] = V
                b_row[self.pos['MAG_B']] = B
                r_row[self.pos['MAG_R']] = R

                if role != 'VF':
                    err_v0 = float(g_row[self.pos[fieldMgErrInstrumental]])
                    err_b0 = float(b_row[self.pos[fieldMgErrInstrumental]])
                    err_r0 = float(r_row[self.pos[fieldMgErrInstrumental]])

                    errV = sqrt(err_v0 * err_v0 + err_Vc * err_Vc + err_vc * err_vc)
                    errB = sqrt(err_b0 * err_b0 + err_Bc * err_Bc + err_bc * err_bc)
                    errR = sqrt(err_r0 * err_r0 + err_Rc * err_Rc + err_rc * err_rc)

                    g_row[self.pos['ERR_V']] = errV
                    b_row[self.pos['ERR_B']] = errB
                    r_row[self.pos['ERR_R']] = errR

                allResults['Gi'].append(g_row)
                allResults['Bi'].append(b_row)
                allResults['Ri'].append(r_row)

        return allResults

    def process(self):
        if not self.opt['files']:
            printError('No input files are given.')
            return

        allCatalogs = {}

#        dateObs = None

        # load refcats
        for j in range(len(self.opt['color'])):

            fileName = self.opt['files'][j]
            color = self.opt['color'][j]

            if not exists(fileName):
                printError("File %s is not exists." % (fileName))
                return

            refCat = loadCatalog(fileName)
            if not self.pos:
                self.pos = self.getHeaderPositions(refCat)

            allCatalogs[color] = refCat

        # find best comp star in Gi frame
        color = 'Gi'  # TODO: what if no Gi color, just G, g, gi, or V ?
        indexGi = self.opt['color'].index(color)
        #dateObs = getDateObs(self.opt['files'][indexGi])
        self.loadFitsHeaders(self.opt['files'][indexGi])
        dateObs = self.fits['DATE-OBS']
        dateObsDate = dateObs.split('T')[0]
        bestComp = self.findBestCompStar(allCatalogs, color)  # refcat record of best comp in Gi
        if not bestComp:
            return False

        # calculate instrumental mgs
        allResults = {}
        if self.opt['useCoeffs'] or self.opt['makeCoeffs']:
            # standardization
            coeffs = None

            if self.opt['makeCoeffs']:
                mergedCat = self.mergeCatalogs(allCatalogs)

                coeffs = self.calculateCoeffs(mergedCat)

                if self.opt['saveCoeffs']:
                    self.saveCoeffs(coeffs, dateObsDate)

            elif self.opt['loadCoeffs']:

                coeffs = self.loadCoeffs(dateObsDate)

            if not coeffs:
                printError('No std coefficients for transformation ; use Tv = 0, Tvr = -1, Tbv = 1 for comp star method.')
                coeffs = [ 0.0, -1.0, 1.0 ]  # this is for comp star method

            allResults = self.calculateStdMgs(coeffs, allCatalogs, bestComp)

        else:

            if self.opt['showGraphs']:
                plt.figure()

            for color in self.opt['color']:

                if self.opt['method'] == 'comp':
                    allResults[color] = self.calculateMgsComparision(refCat, color, bestComp)

                elif self.opt['method'] == 'lfit':
                    allResults[color] = self.calculateMgsLinearFit(refCat, color, bestComp)

                else:
                    allResults[color] = self.calculateMgsRobustAveraging(refCat, color, bestComp)

            if self.opt['showGraphs']:
                plt.show()

        # save results
        for j in range(len(self.opt['color'])):

            fileName = self.opt['files'][j]
            color = self.opt['color'][j]

            self.reportResult(allResults[color], allCatalogs[color], fileName + '.pm', color, dateObs)

        return True


class MainApp:

    opt = {
        'out' : None,  # output photometry file name, if None, '*.pm' will be created
        'method': 'gcx',  # mg calculation method: comp - use best comp star, gcx - m=1 linear fit ensemble, lfit - general linear fit ensemble
        'color' : ['Gi'],  # photometry bands
        'loadCoeffs' : False,  # load std coeffs
        'useCoeffs': False,  # use std coeffs
        'makeCoeffs': False,  # create std coeffs
        'saveCoeffs': False,  # save std coeffs
        'showGraphs': False,  # show standard coefficient graphs
        'observer': 'XXX',  # observer code
        'camera'   : None,   # camera name
        'telescope': None,   # telescope name
        'files': None,
        }

    def __init__(self, argv):
        self.argv = argv
        pass

    def processCommands(self):
        try:
            optlist, args = getopt (argv[1:], "c:smo:p:", ['color=', 'std', 'make-std', 'out=', 'comp='])
        except GetoptError:
            printError('Invalid command line options')
            return

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-c' or o == '--color':
                self.opt['color'] = [a]
            elif o == '-s' or o == '--std':
                self.opt['std'] = True
            elif o == '-m' or o == '--make-std':
                self.opt['makestd'] = True
            elif o == '-p' or o == '--comp':
                self.opt['comp'] = a
            elif o == '-o' or o == '--':
                self.opt['out'] = a

        self.opt['files'] = args

    def run(self):
        self.processCommands()

        ppl = loadPplSetup()

        phot = Photometry(self.opt, ppl)
        phot.process()


if __name__ == '__main__':

    app = MainApp(argv)
    app.run()

# end __main__
