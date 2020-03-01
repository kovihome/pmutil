#!/usr/bin/env python3
#
# PmUtils/pmfilter
#

'''
Created on Apr 28, 2019

@author: kovi
'''

from getopt import getopt, GetoptError
from sys import argv
from subprocess import Popen, PIPE

from pmbase import invoke


def loadRefCatalog(refFileName):

    refCat = { 'header': [], 'cat': [] }

    ref = open(refFileName, 'r')
    if ref:
        headerLine = ref.readline()
        refCat['header'] = headerLine.split()

        for refLine in ref:
            if not refLine.startswith('#') and refLine.strip() != '':
                refCatLine = refLine.split()
                refCat['cat'].append(refCatLine)
        ref.close()

        print (("Reference catalog: %s") % refFileName)
        print (("%d reference object loaded." % len(refCat['cat'])))
        # print (refCat)
        return refCat

    else:

        print (("Reference catalog %s not found") % refFileName)
        return None


def loadPmCatalog(pmFileName):

    pmCat = { 'header': [], 'cat': [] }
    pm = open(pmFileName, 'r')
    if pm:

        headerLine = pm.readline()
        while headerLine.startswith('#'):
            header = headerLine.split()
            pmCat['header'].append(header[2])
            headerLine = pm.readline()

        pmCatLine = headerLine.split()
        pmCat['cat'].append(pmCatLine)

        for pmLine in pm:
            pmCatLine = pmLine.split()
            pmCat['cat'].append(pmCatLine)
        pm.close()

        print (("Photometry catalog: %s") % pmFileName)
        print (("%d measured object loaded." % len(pmCat['cat'])))
        # print (pmCat)
        return pmCat

    else:

        print (("Reference catalog %s not found") % pmFileName)
        return None

def loadIdCatalog(idFileName):

    idCat = []

    ids = open(idFileName, 'r')
    if ids:
        for idLine in ids:
            if idLine.strip() != '':
                idCatLine = idLine.split()
                idCat.append(idCatLine)
        ids.close()

        print ("%d matched object loaded." % (len(idCat)))
        return idCat

    else:

        print ("Id catalog %s not found" % (idFileName))
        return None




def getHeaderPos(cat, headerName):
    index = 0
    for header in cat['header']:
        if header == headerName:
            return index
        index = index + 1
    return None

def findId(idCat, refId):
    for idr in idCat:
        if refId == idr[0]:
            return idr[1]
    return None

def findPmRec(pmCat, pmId, pmcIdPos):
    for pmr in pmCat:
        if pmId == pmr[pmcIdPos]:
            return pmr
    return None

class CatalogMatcher:

    opt = {}                 # command line options
#    pplSetup = {}            # PPL setup from ppl-setup config file
    pos = None               # catalog header positions

    refHeaders = ['AUID', 'ROLE', 'RA', 'RA_DEG', 'DEC', 'DEC_DEG', 'MAG_V', 'ERR_V', 'MAG_B', 'ERR_B', 'MAG_R', 'ERR_R']

    pmHeaders = ['NUMBER', 'MAG_ISOCOR', 'MAGERR_ISOCOR', 'MAG_BEST', 'MAGERR_BEST', 'ALPHA_J2000', 'DELTA_J2000']

    refTrailerHeaders = [ 'LABEL']

    HMG_MAX_ERR = 0.4
    HMG_SIGMA = 3.0

    def __init__(self, opt):
        self.opt = opt

    def hmg(self, pmCat):
        mg1Pos = getHeaderPos(pmCat, 'MAG_ISOCOR')
        mg1ErrPos = getHeaderPos(pmCat, 'MAGERR_ISOCOR')
        mg2Pos = getHeaderPos(pmCat, 'MAG_BEST')
        mg2ErrPos = getHeaderPos(pmCat, 'MAGERR_BEST')

        hmg1 = -99.99
        hmg2 = -99.99
        hmg1Err = 99.99
        hmg2Err = 99.99
        for pm in pmCat['cat']:
            mg1 = float(pm[mg1Pos])
            mg1Err = float(pm[mg1ErrPos])
            mg2 = float(pm[mg2Pos])
            mg2Err = float(pm[mg2ErrPos])
            if mg1Err < self.HMG_MAX_ERR and mg2 < self.HMG_MAX_ERR:
                if mg1 - self.HMG_SIGMA * mg1Err > hmg1:
                    hmg1 = mg1 - self.HMG_SIGMA * mg1Err
                    hmg1Err = mg1Err
                if mg2 - self.HMG_SIGMA * mg2Err > hmg2:
                    hmg2 = mg2 - self.HMG_SIGMA * mg2Err
                    hmg2Err = mg2Err

#        return [hmg1 - self.HMG_SIGMA * hmg1Err, hmg2 - self.HMG_SIGMA * hmg2Err]
        return [hmg1, hmg2]

    def matchCatalogsByGrmatch(self, refCatFile, pmCatFile, outFile):

        idFile = pmCatFile + ".idmatch"
        cmd = "grmatch -r %s -i %s --match-coord --col-ref 4,6 --col-inp 14,15 --output-id %s --col-ref-id 1 --col-inp-id 1" % (refCatFile, pmCatFile, idFile)
        result = invoke(cmd)
#        p = Popen(cmd, stdout = PIPE, shell = True)
#        (output, errno) = p.communicate()
#        p.wait()

        #output.decode('ascii').strip()

        refCat = loadRefCatalog(self.opt['ref'])
        pmCat = loadPmCatalog(self.opt['files'][0])
        idCat = loadIdCatalog(idFile)

        refIdPos = getHeaderPos(refCat, 'AUID')
        labelPos = getHeaderPos(refCat, 'LABEL')
        pmcIdPos = getHeaderPos(pmCat, 'NUMBER')
        rolePos = getHeaderPos(refCat, 'ROLE')
        mg1Pos = getHeaderPos(pmCat, 'MAG_ISOCOR')
        mg2Pos = getHeaderPos(pmCat, 'MAG_BEST')

        rlen = len(refCat['header'])
        pmlen = len(pmCat['cat'])
        hmgs = self.hmg(pmCat)
        for ref in refCat['cat']:
            refId = ref[refIdPos]
            pmId = findId(idCat, refId)
            pm = findPmRec(pmCat['cat'], pmId, pmcIdPos)
            ref.append(99.99)
            ref.append(-1)
            if pmId != None:
                if ref[rolePos] == 'V' and float(pm[mg2Pos]) > hmgs[1]:
                    ref[rlen + 1] = pmlen
                    ref[rolePos] = 'VF'
                else:             
                   ref[rlen + 1] = int(pmId) - 1
            else:
                if ref[rolePos] == 'V':
                    ref[rlen + 1] = pmlen
                    ref[rolePos] = 'VF'
                print("Ref object not matched:", refId, ref[labelPos])

#   1 NUMBER                 Running object number                                     
#   2 FLUX_ISOCOR            Corrected isophotal flux                                   [count]
#   3 FLUXERR_ISOCOR         RMS error for corrected isophotal flux                     [count]
#   4 MAG_ISOCOR             Corrected isophotal magnitude                              [mag]
#   5 MAGERR_ISOCOR          RMS error for corrected isophotal magnitude                [mag]
#   6 FLUX_AUTO              Flux within a Kron-like elliptical aperture                [count]
#   7 FLUXERR_AUTO           RMS error for AUTO flux                                    [count]
#   8 FLUX_BEST              Best of FLUX_AUTO and FLUX_ISOCOR                          [count]
#   9 FLUXERR_BEST           RMS error for BEST flux                                    [count]
#  10 MAG_BEST               Best of MAG_AUTO and MAG_ISOCOR                            [mag]
#  11 MAGERR_BEST            RMS error for MAG_BEST                                     [mag]
#  12 X_WORLD                Barycenter position along world x axis                     [deg]
#  13 Y_WORLD                Barycenter position along world y axis                     [deg]
#  14 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
#  15 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
#  16 XWIN_IMAGE             Windowed position estimate along x                         [pixel]
#  17 YWIN_IMAGE             Windowed position estimate along y                         [pixel]
#  18 ERRAWIN_IMAGE          RMS windowed pos error along major axis                    [pixel]
#  19 ERRBWIN_IMAGE          RMS windowed pos error along minor axis                    [pixel]
#  20 ERRTHETAWIN_IMAGE      Windowed error ellipse pos angle (CCW/x)                   [deg]
#  21 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
#  22 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
        pmCat['cat'].append([str(pmlen), "0.0", "0.0", "%7.4f" % (hmgs[0]), "0.0", "0.0", "0.0", "0.0", "0.0", "%7.4f" % (hmgs[1]), "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0"])
        print("Instumental HMG: HMG_ISOCORR = %8.4f, HMG_BEST = %8.4f" % (hmgs[0], hmgs[1]))

        self.dumpResult(refCat, pmCat, outFile)


    def matchCatalogsNaiveApproach(self, refCatFile, pmCatFile, outFile):

        refCat = loadRefCatalog(self.opt['ref'])
        pmCat = loadPmCatalog(self.opt['files'][0])

        refcRAFieldPos = getHeaderPos(refCat, 'RA_DEG')
        refcDECFieldPos = getHeaderPos(refCat, 'DEC_DEG')
        pmcRAFieldPos = getHeaderPos(pmCat, 'ALPHA_J2000')
        pmcDECFieldPos = getHeaderPos(pmCat, 'DELTA_J2000')

        rlen = len(refCat['header'])
        for ref in refCat['cat']:
            ref.append(99.99)
            ref.append(-1)

        pmi = 0
        for pm in pmCat['cat']:
            pmRA = float(pm[pmcRAFieldPos])
            pmDEC = float(pm[pmcDECFieldPos])
            for ref in refCat['cat']:
                refRA = float(ref[refcRAFieldPos])
                refDEC = float(ref[refcDECFieldPos])
                dRA = refRA - pmRA
                dDEC = refDEC - pmDEC
                dist2 = dRA * dRA + dDEC * dDEC
                if dist2 < ref[len]:
                    ref[rlen] = dist2
                    ref[rlen + 1] = pmi
            pmi = pmi + 1

        self.dumpResult(refCat, pmCat, outFile)


    def dumpResult(self, refCat, pmCat, outFileName):

        outf = open(outFileName, 'w')
        for h in self.refHeaders:
            outf.write(h.ljust(15))
        for h in self.pmHeaders:
            outf.write(h.ljust(15))
        for h in self.refTrailerHeaders:
            outf.write(h.ljust(30))
        outf.write('\n')

        rlen = len(refCat['header'])
        for ref in refCat['cat']:
            if ref[rlen + 1] != -1:
                pm = pmCat['cat'][ref[rlen + 1]]
    
                for h in self.refHeaders:
                    pos = getHeaderPos(refCat, h)
                    outf.write(ref[pos].ljust(15))
                for h in self.pmHeaders:
                    pos = getHeaderPos(pmCat, h)
                    outf.write(pm[pos].ljust(15))
                for h in self.refTrailerHeaders:
                    pos = getHeaderPos(refCat, h)
                    outf.write(ref[pos].ljust(30))
                outf.write('\n')

        outf.close()
        print ('outpuf file: ', outFileName)

    def process(self):
        # if self.opt['ref']:

        self.matchCatalogsByGrmatch(self.opt['ref'], self.opt['files'][0], self.opt['out'])

        # self.matchCatalogsNaiveApproach(self.opt['ref'], self.opt['files'][0], self.opt['out'])


class MainApp:
    
    opt = {
        'ref' : None,
        'out' : None,
        'color': 'Gi',
        'files': None,
        }


    def __init__(self, argv):
        self.argv = argv
        pass


    def processCommands(self):
        try:
            optlist, args = getopt (argv[1:], "r:o:c:", ['--ref', '--out', '--color'])
        except GetoptError:
            print ('Invalid command line options')
            return

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-r':
                self.opt['ref'] = a
            elif o == '-o':
                self.opt['out'] = a
            elif o == '-c':
                self.opt['color'] = a

        self.opt['files'] = args
        if not self.opt['out']:
            self.opt['out'] = args[0] + '.refout'


    def run(self):
        self.processCommands()

        matcher = CatalogMatcher(self.opt)
        matcher.process()


if __name__ == '__main__':

    app = MainApp(argv)
    app.run()

# end main.
