#
# PmUtils/pmphot
#
'''
Created on Dec 28, 2019

@author: kovi
'''
from numpy import median, sqrt  # , polyfit
from numpy.polynomial.polynomial import Polynomial
from getopt import getopt, GetoptError
from sys import argv
from glob import glob
from astropy.table import Table

from pmbase import jd, getFitsHeader, printError, loadPplSetup


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
        # print (("%d reference object loaded." % len(refCat['cat'])))
        # print (refCat)
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


def getDateObs(fileName):
    '''
    fileName: path/[Seq_nnn|Combined]_c.fits.cat.cat
    '''
    astFileName = fileName.split('.')[0]
    astFileName = astFileName + '.ast.fits'

    return getFitsHeader(astFileName, 'DATE-OBS')


class Photometry:

    opt = {}  # command line options
    ppl = {}  # PPL setup from ppl-setup config file
    pos = None  # catalog header positions

    def __init__(self, opt, ppl):
        self.opt = opt
        self.ppl = ppl

    def transformMgs(self, refCat, stdcolor, p, err):
        mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
        miPos = getHeaderPos(refCat, fieldMgInstrumental)
        evPos = getHeaderPos(refCat, 'ERR_' + stdcolor)
        eiPos = getHeaderPos(refCat, fieldMgErrInstrumental)
        rolePos = getHeaderPos(refCat, fieldRole)

        result = []
        for pm in refCat['cat']:
            if pm[rolePos] == 'V' or pm[rolePos] == 'VF':
#            if pm[rolePos] == 'V':
                mv = p[0] * float(pm[miPos]) + p[1]
                pm[mvPos] = mv
                if err:
                    pm[evPos] = max(float(pm[eiPos]), err)
                result.append(pm)
        # print (result)
        return result

    def calculateMgsRobustAveraging(self, refCat, color):
        '''
        Calculate unit slope linear fit regression with robust averaging
        see: http://gcx.sourceforge.net/html/node11.html
        inputs:    I[k] - instrumental mgs              MAG_ISOCORR | MAG_BEST
                   S[k] - standard mgs                  MAG_R \ MAG_V | MAG_B
                   ei[k] - error of instrumental mgs    MAGERR_ISOCORR | MAGERR_BEST
                   es[k] - error of standard mgs        ERR_R | ERR_V | ERR_B
        '''
        stdcolor = color[:1].upper()
        if stdcolor == 'G':
            stdcolor = 'V'

        mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
        miPos = getHeaderPos(refCat, fieldMgInstrumental)
        evPos = getHeaderPos(refCat, 'ERR_' + stdcolor)
        eiPos = getHeaderPos(refCat, fieldMgErrInstrumental)
        rolePos = getHeaderPos(refCat, fieldRole)

        y = []
        e2 = []
        w = []
        for pm in refCat['cat']:
            if pm[rolePos] == 'C' and pm[mvPos] != '-':
                mi = float(pm[miPos])
                mv = float(pm[mvPos])
                ei = float(pm[eiPos])
#            if not pm[evPos].startswith('-'):
                if pm[evPos][0].isdigit():
#                print("ev:[" + pm[evPos] + "]")
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

#    print("y:", y)
#    print("w:", w)

        # iteration starts
        for it in range(10):

#        print("r:", r)
#        print("r_:", r_)
#        print("w_:", w_)
#        print("z:", z)

            zsum = 0.0
            wsum = 0.0
            beta = 1.0
            alpha = 1.0
            for k in range(len(y)):
                r[k] = y[k] - z
                r_k = (y[k] - z) * sqrt(w[k])
                w_k = w[k] / (1 + (r[k] / alpha) ** beta)
                w_[k] = w_k
                r_[k] = r_k
                zsum = zsum + (y[k]) * w_k
                wsum = wsum + w_k
            z = zsum / wsum
        # iteration ends

#    print("r:", r)
#    print("r_:", r_)
#    print("w_:", w_)
#    print("zp:", z)

        su = 0.0
        sl = 0.0
        for k in range(len(y)):
            su = su + r[k] * r[k] * w_[k]
            sl = sl + w_[k]
        ez_2 = su / sl
        mel_2 = su / float(len(y) - 1)
        print ("zp:", z, "mel^2:", mel_2, "ez^2:", ez_2, "ez:", sqrt(ez_2), "len y:", len(y))

        # apply result to variables
        result = self.transformMgs(refCat, stdcolor, [1.0, z], sqrt(ez_2))
#        result = []
#        for pm in refCat['cat']:
#            if pm[rolePos] == 'V' or pm[rolePos] == 'VF':
#                mv = float(pm[miPos]) + z
#                pm[mvPos] = mv
#                pm[evPos] = sqrt(ez_2)
#                result.append(pm)
        # print (result)
        return result

    def calculateMgsLinearFit(self, refCat, color):
        '''
        Calculate magnitudes with naive linear fit
        '''
        # TODO: calculate comp mag error

        stdcolor = color[:1].upper()
        if stdcolor == 'G':
            stdcolor = 'V'

        mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
        miPos = getHeaderPos(refCat, fieldMgInstrumental)
        evPos = getHeaderPos(refCat, 'ERR_' + stdcolor)
        eiPos = getHeaderPos(refCat, fieldMgErrInstrumental)
        rolePos = getHeaderPos(refCat, fieldRole)

        mi = []
        mv = []
        er = []
#    p = [1.0, 0.0]
        for pm in refCat['cat']:
            if pm[rolePos] == 'C' and pm[mvPos] != '-':
                mi.append(float(pm[miPos]))
                mv.append(float(pm[mvPos]))
                ei = float(pm[eiPos])
                ev = float(pm[evPos])
                er.append(1.0 / sqrt(ei * ei + ev * ev))

        p = Polynomial.fit(mi, mv, 1, w = er)
        # print ('polyfit result:', p)

        # apply result to variables
        result = self.transformMgs(refCat, stdcolor, [p.coef[1], p.coef[0]])
#        result = []
#        for pm in refCat['cat']:
#            if pm[rolePos] == 'V':
#                mv = p.coef[1] * float(pm[miPos]) + p.coef[0]
#                pm[mvPos] = mv
#                result.append(pm)
        # print (result)
        return result

    def calculateMgsComparision(self, refCat, color):
        stdcolor = color[:1].upper()
        if stdcolor == 'G':
            stdcolor = 'V'

        idPos = getHeaderPos(refCat, fieldAuid)
        mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
        miPos = getHeaderPos(refCat, fieldMgInstrumental)
        evPos = getHeaderPos(refCat, 'ERR_' + stdcolor)
        eiPos = getHeaderPos(refCat, fieldMgErrInstrumental)
        rolePos = getHeaderPos(refCat, fieldRole)

        # find best comp star
        p = [1.0, 0.0]
        emin = 1.0
        ep = 0.0
        compId = None
        for pm in refCat['cat']:
            if pm[rolePos] == 'C' and pm[mvPos] != '-':
#               if comp != None and pm[idPos] == comp:
#                   p = float(pm[mvPos]) - float(pm[miPos])
#                   ep = max(float(pm[evPos]), float(pm[eiPos])
                ei = float(pm[eiPos]) if pm[eiPos] != '-' else MAG_ERR_DEFAULT
                ev = float(pm[evPos]) if pm[evPos] != '-' else MAG_ERR_DEFAULT
                e = ei * ei + ev * ev
                if e < emin:
                    emin = e
                    p[1] = float(pm[mvPos]) - float(pm[miPos])
                    ep = max(ev, ei)
                    compId = pm[idPos]
        print('Best comp star: %s, poly: [ %4.3f, %4.3f ], err: %4.3f' % (compId, p[0], p[1], ep))

        # apply result to variables
        result = self.transformMgs(refCat, stdcolor, p, ep)
#        result = []
#        for pm in refCat['cat']:
#            if pm[rolePos] == 'V':
#                mv = p[0] * float(pm[miPos]) + p[1]
#                pm[mvPos] = mv
#                pm[evPos] = max(float(pm[eiPos]), ep)
#                result.append(pm)
        # print (result)
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

        # print(allResults)
        stdcolor = color[:1].upper()
        if stdcolor == 'G':
            stdcolor = 'V'

        idPos = getHeaderPos(refCat, fieldAuid)
        labelPos = getHeaderPos(refCat, 'LABEL')
        astIdPos = getHeaderPos(refCat, 'NUMBER')
        mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
        evPos = getHeaderPos(refCat, fieldMgErrInstrumental)
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
#        for dateObs in allResults.keys():
        jdObs = jd(dateObs)
        jds = "%12.4f" % jdObs
#            result = allResults[dateObs]
        for res in result:
            r.write(res[idPos].ljust(15))
            r.write(res[labelPos].ljust(30))
            r.write(res[astIdPos].ljust(10))
            r.write(res[raPos].ljust(20))
            r.write(res[decPos].ljust(20))
            r.write(dateObs.ljust(20))
            r.write(jds.ljust(20))
            flag = "F" if res[rolePos] == "VF" else "-"
            r.write(flag.ljust(5))
            r.write(("%6.3f" % float(res[mvPos])).ljust(20))
            r.write(color.ljust(5))
            r.write(("%5.3f" % float(res[evPos])).ljust(20))
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
            stdcolor = color[:1].upper()
            if stdcolor == 'G':
                stdcolor = 'V'

            mgsPos = self.pos['MAG_' + stdcolor]

            result = allCatalogs[color]['cat']
#            print(self.opt['color'][j],'result',result)
            ids = set()
            for row in result:
#                print('mgipos',mgiPos,'row',row)
                role = row[rolePos]
                mgi = row[mgiPos]
                mgs = row[mgsPos]
                auid = row[auidPos]
                if role == 'C' and mgs and mgs != '-' and mgi and mgi != '-':
                    ids.add(auid)

#           print('ids',self.opt['color'][j],ids)
            if len(allIds) == 0:
                allIds = ids
            else:
                allIds = allIds.intersection(ids)

        # find records for all ids
        # merged catalog fields: ID, B, V, R, bi, gi, ri
#        print('allIds:',allIds)
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
#        print(mergedCat)
        for row in mergedCat:
            V_vi.append(float(row[2]) - float(row[5]))
            V_R.append(float(row[2]) - float(row[3]))
            vi_ri.append(float(row[5]) - float(row[6]))
            bi_vi.append(float(row[4]) - float(row[5]))
            B_V.append(float(row[1]) - float(row[2]))

            for p in range(7, 12):
                if row[p] == '-':
                    row[p] = "0.1"
            e_V_vi = max(float(row[8]), float(row[11]))
            e_V_R = max(float(row[8]), float(row[9]))
            e_vi_ri = max(float(row[11]), float(row[12]))
            e_bi_vi = max(float(row[10]), float(row[11]))
            e_B_V = max(float(row[7]), float(row[8]))
            w_Tv.append(1.0 / sqrt(e_V_vi * e_V_vi + e_V_R * e_V_R))
            w_Tvr.append(1.0 / sqrt(e_vi_ri * e_vi_ri + e_V_R * e_V_R))
            w_Tbv.append(1.0 / sqrt(e_bi_vi * e_bi_vi + e_B_V * e_B_V))

        print('V-vi:', V_vi)
        print('V-R:', V_R)
        print('vi-ri:', vi_ri)
        print('bi-vi:', bi_vi)
        print('B-V:', B_V)

        print('w(Tv):' , w_Tv)
        print('w(Tvr):' , w_Tvr)
        print('w(Tbv):' , w_Tbv)

        p = Polynomial.fit(V_vi, V_R, 1, w = w_Tv)
        print(p.coef)
        Tv = p.coef[1]
        p = Polynomial.fit(vi_ri, V_R, 1, w = w_Tvr)
        Tvr = 1.0 / p.coef[1]
        p = Polynomial.fit(bi_vi, B_V, 1, w = w_Tbv)
        Tbv = p.coef[1]

        print ('standard coeffs: Tv =', Tv, 'Tvr =', Tvr, 'Tbv = ', Tbv)
        return [Tv, Tvr, Tbv]

    def saveCoeffs(self, coeffs):
        configFolder = self.ppl['CONFIG_FOLDER'].strip('/')
        coeffFile = configFolder + '/stdcoeffs.cat'
        if not exists(coeffFile):
            f = open(coeffFile, "w")
            f.write('DATE          TV         ERR_TV    TVR        ERR_TVR   TBV        ERR_TBV   CAMERA               TELESCOPE                 FIELD\n')
        else:
            f = open(coeffFIle, "a+")
        f.write("%-14s" % ('2020-02-02')) # TODO: date of image
        f.write("%-11.4f" % (coeff[0]))   # Tv
        f.write("%-11.4f" % (0.0))        # TODO: error of Tv
        f.write("%-11.4f" % (coeff[1]))   # Tvr
        f.write("%-11.4f" % (0.0))        # TODO: error of Tvr
        f.write("%-11.4f" % (coeff[2]))   # Tbv
        f.write("%-11.4f" % (0.0))        # TODO: error of Tbv
        f.write("%-21s" % ('Canon EOS 1100D')) # TODO: camera type from image
        f.write("%-21s" % ('250/1200 T')) # TODO: telescope
        f.write("%-21s" % ('SA104'))      # TODO: standard field or target object
        f.write("\n")
        f.close()

    def loadCoeffs(self):
        # TODO: search for appropriate coeffs (i.e. matching camera and telescope, near of image date, etc.), or average coeffs of same camera/telescope
        configFolder = self.ppl['CONFIG_FOLDER'].strip('/')
        coeffFile = configFolder + '/stdcoeffs.cat'
        if not exists(coeffFile):
            return None
#        f = open(coeffFile,)
        coeffTable = Table.read(coeffFile, format='ascii')
        lastRow = coeffTable[len(coeffTable)-1]
        return [lastRow['TV'], lastRow['TVR'], lastRow['TBV']]

    def calculateStdMgs(self, coeffs, mergedCat):
        pass

    def process(self):
        if not self.opt['files']:
            f = glob('Phot/Seq_*.fits.cat.cat')
            f.sort()
            self.opt['files'] = f

        allResults = {}
        allCatalogs = {}

        for j in range(len(self.opt['color'])):

            fileName = self.opt['files'][j]
            color = self.opt['color'][j]

            refCat = loadCatalog(fileName)
            if not self.pos:
                self.pos = self.getHeaderPositions(refCat)

            dateObs = getDateObs(fileName)

#            if self.opt['comp']:
            if self.opt['method'] == 'comp':
                result = self.calculateMgsComparision(refCat, color)

            elif self.opt['method'] == 'lfit':
                result = self.calculateMgsLinearFit(refCat, color)

            else:
                result = self.calculateMgsRobustAveraging(refCat, color)

            allResults[color] = result
            allCatalogs[color] = refCat

            self.reportResult(result, refCat, fileName + '.pm', color, dateObs)

        coeffs = None

        # if not option == load_coeffs
        if self.opt['makeCoeffs']:
            #   merge catalogs -> id, B, V, R, bi, vi, ri
            mergedCat = self.mergeCatalogs(allCatalogs)

            #   calculate coeffs -> Tv, Tvr, Tbv
            coeffs = self.calculateCoeffs(mergedCat)

            #   if option == make_coeffs
            if self.opt['saveCoeffs']:

                #     save coeffs
                self.saveCoeffs(coeffs)

        elif self.opt['loadCoeffs']:

            #   load coeffs
            coeffs = self.loadCoeffs()

        # if option == use_coeffs
        if coeffs and self.opt['useCoeffs']:

            #   calculate std mgs
            self.calculateStdMgs(coeffs, mergedCat)

            # save result


class MainApp:

    opt = {
        'out' : None,  # output photometry file name, if None, '*.pm' will be created
        'method': 'gcx',  # mg calculation method: comp - use best comp star, gcx - m=1 linear fit ensemble, lfit - general linear fit ensemble
        'color' : ['Gi'],  # photometry bands
        'loadCoeffs' : False,  # load std coeffs
        'useCoeffs': False,  # use std coeffs
        'makeCoeffs': False,  # create std coeffs
        'saveCoeffs': False,  # save std coeffs
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
