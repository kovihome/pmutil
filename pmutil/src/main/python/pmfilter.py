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
from astropy.table import Table
from astropy.io import fits
from os import remove
from os.path import exists
from math import sqrt, log10
from numpy import arange

import matplotlib.pyplot as plt

from pmbase import invoke, getFitsHeaders, linfit


class CatalogMatcher:

    opt = {}  # command line options
#    pplSetup = {}            # PPL setup from ppl-setup config file
    pos = None  # catalog header positions

    refHeaders = ['AUID', 'ROLE_', 'RA', 'RA_DEG', 'DEC', 'DEC_DEG', 'MAG_V', 'ERR_V', 'MAG_B', 'ERR_B', 'MAG_R', 'ERR_R']

    pmHeaders = ['NUMBER', 'MAG_ISOCOR', 'MAGERR_ISOCOR', 'MAG_BEST', 'MAGERR_BEST', 'ALPHA_J2000', 'DELTA_J2000']

    refTrailerHeaders = [ 'LABEL']

    HMG_MAX_ERR = 0.15

    IMAGE_BORDER_SIZE = 10  # pixels

    def __init__(self, opt):
        self.opt = opt

    def hmg(self, pmTable):
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

        m1,b1 = linfit(logm1, loge1)
        hmg1 = -1.0 * 10.0 ** ((log10(self.HMG_MAX_ERR) - b1) / m1);

        m2,b2 = linfit(logm2, loge2)
        hmg2 = -1.0 * 10.0 ** ((log10(self.HMG_MAX_ERR) - b2) / m2);

#        plotcolor = 'b'
#        xr = arange(min(logm), max(logm), 0.01)
#        plt.plot(logm, loge, plotcolor + 'o', xr, m * xr + b, 'k')
#        plt.xlabel('log mi')
#        plt.ylabel('log ei')
#        plt.show()

        return [hmg1, hmg2]

    def readFrameSize(self, fitsFileName):
        h = getFitsHeaders(fitsFileName, ['NAXIS1', 'NAXIS2'])
        xf = int(h['NAXIS1'])
        yf = int(h['NAXIS2'])
        return (xf, yf)

    def determineOnFrameStatus(self, x_obj, y_obj, im_w, im_h):
        if x_obj < 0 or x_obj >= im_w or y_obj < 0 or y_obj >= im_h:
            return 'O'
        if x_obj < self.IMAGE_BORDER_SIZE or x_obj >= im_w - self.IMAGE_BORDER_SIZE or y_obj < self.IMAGE_BORDER_SIZE or y_obj >= im_h - self.IMAGE_BORDER_SIZE:
            return 'B'
        return ''

    def matchCatalogsByGrmatchPoints(self, refCatFile, pmCatFile, outFile):
        '''
        Match reference catalog with sextractor's .cat file by frame xy points
        '''

        # 1. convert refCat to fits format
        refFitsFile = refCatFile + '.fits'
        table = Table.read(refCatFile, format = 'ascii')
        if not exists(refFitsFile):
            table.write(refFitsFile)

        # 2. calculate ref objects' frame xy points
        wcsFile = pmCatFile.replace('.cat', '.wcs')
        axyFile = pmCatFile.replace('.cat', '.ref.axy')
        invoke("wcs-rd2xy -w %s -i %s -o %s -R RA_DEG -D DEC_DEG" % (wcsFile, refFitsFile, axyFile))  # -f option need argument

        remove(refFitsFile)

        # 3. merge frame xy point to refCat
        tlen = len(table)
        table['X'] = [0.0] * tlen
        table['Y'] = [0.0] * tlen
        table['ID'] = [0] * tlen
        table['ROLE_'] = ['AA'] * tlen
        table['HMG_'] = [0.0] * tlen

        f = fits.open(axyFile)
        d = f[1].data

        for j in range(tlen):
            x, y = d[j]
            table[j]['X'] = x
            table[j]['Y'] = y

        # 4. match refCat with sextractor's cat by frame xy points
        idFile = pmCatFile + ".idmatch"
        result = invoke("grmatch -r %s -i %s --match-coord --col-ref 4,6 --col-inp 14,15 --output-id %s --col-ref-id 1 --col-inp-id 1" % (refCatFile, pmCatFile, idFile))

        # get image size
        astFile = pmCatFile.replace('.cat', '.ast.fits')
        im_x, im_y = self.readFrameSize(astFile)

        idTable = Table.read(idFile, format = 'ascii')
        pmTable = Table.read(pmCatFile, format = 'ascii')

        hmgs = self.hmg(pmTable)
        print("Instumental HMG: HMG_ISOCORR = %8.4f, HMG_BEST = %8.4f" % (hmgs[0], hmgs[1]))

        for j in range(tlen):

            onFrameStatus = self.determineOnFrameStatus(table[j]['X'], table[j]['Y'], im_x, im_y)

            matched = False
            for k in range(len(idTable)):
                if table[j]['AUID'] == idTable[k][0]:
                    table[j]['ID'] = int(idTable[k][1])
                    matched = True
                    break

            pmid = int(table[j]['ID'])
            d = -1.0
            if id != 0:
                pmrow = pmTable[pmid - 1]
                dx = pmrow['XWIN_IMAGE'] - table[j]['X']
                dy = pmrow['YWIN_IMAGE'] - table[j]['Y']
                d = sqrt(dx * dx + dy * dy)

            if onFrameStatus == '' and matched:
                # print('AUID: %s, ID:%d, d:%f, OK, on-frame and mathed' % (table[j]['AUID'], pmid, d))
                if d > 2.0:
                    print('AUID: %s, ID:%d, d:%f, OK, on-frame and mathed' % (table[j]['AUID'], pmid, d))
                    print('   and too large distance')
                if table[j]['ROLE'] == 'V' and float(pmrow['MAG_BEST']) > hmgs[1]:
                    print('   and under limit, mi:%7.3f, hmg:%7.3f' % (float(pmrow['MAG_BEST']), hmgs[1]))
                    table[j]['ROLE_'] = 'VF'
                    table[j]['HMG_'] = hmgs[1]
                else:
                    table[j]['ROLE_'] = table[j]['ROLE']
            elif onFrameStatus == 'B' and matched:
                table[j]['ROLE_'] = table[j]['ROLE'] + 'B'
                print('AUID: %s, ID:%d, d:%f, OK, on-frame-border and mathed' % (table[j]['AUID'], pmid, d))
                if d > 2.0:
                    print('   and too large distance')
                if table[j]['ROLE'] == 'VB' and float(pmrow['MAG_BEST']) > hmgs[1]:
                    print('   and under limit, mi:%7.3f, hmg:%7.3f' % (float(pmrow['MAG_BEST']), hmgs[1]))
                    table[j]['ROLE_'] = 'VF'
            elif onFrameStatus == 'O' and matched:
                table[j]['ROLE_'] = table[j]['ROLE'] + 'O'
                print('AUID: %s, ID:%d, d:%f, BAD-MATCH, out-of-frame and mathed' % (table[j]['AUID'], pmid, d))
            elif onFrameStatus == '' and not matched:
                table[j]['ROLE_'] = table[j]['ROLE'] + 'F'
                table[j]['HMG_'] = hmgs[1]
                print('AUID: %s, FAINTER, on-frame and not mathed' % (table[j]['AUID']))
                print('   and under limit, hmg:%7.3f' % (hmgs[1]))
            elif onFrameStatus == 'B' and not matched:
                table[j]['ROLE_'] = table[j]['ROLE'] + 'O'
                print('AUID: %s, FAINTER-ON-BORDER, on-frame-border and not mathed' % (table[j]['AUID']))
            elif onFrameStatus == 'O' and not matched:
                table[j]['ROLE_'] = table[j]['ROLE'] + 'O'
                print('AUID: %s, OK, out-of-frame and not mathed' % (table[j]['AUID']))

        # save result
        self.dumpResult(table, pmTable, outFile)

    def dumpResult(self, refTable, pmTable, outFileName):

        outf = open(outFileName, 'w')
        for h in self.refHeaders:
            if h.endswith('_'):
                h = h[:-1]
            outf.write(h.ljust(15))
        for h in self.pmHeaders:
            outf.write(h.ljust(15))
        for h in self.refTrailerHeaders:
            outf.write(h.ljust(30))
        outf.write('\n')

        goodRoles = [ 'C', 'V', 'VF' ]
        for ref in refTable:
            if ref['ROLE_'] not in goodRoles:
                continue
            pm = None
            if ref['ID'] != 0:
                pm = pmTable[ref['ID'] - 1]
            else:
                pm = { 'NUMBER':'-', 'MAG_ISOCOR': ref['HMG_'], 'MAGERR_ISOCOR':'-', 'MAG_BEST': ref['HMG_'], 'MAGERR_BEST':'-', 'ALPHA_J2000':'-', 'DELTA_J2000':'-' }

            for h in self.refHeaders:
                v = str(ref[h])
                if v == '-1.0':
                    v = '-'
                outf.write(v.ljust(15))
            for h in self.pmHeaders:
                outf.write(str(pm[h]).ljust(15))
            for h in self.refTrailerHeaders:
                outf.write(str(ref[h]).ljust(30))
            outf.write('\n')

        outf.close()
        print ('outpuf file: ', outFileName)

    def process(self):

        self.matchCatalogsByGrmatchPoints(self.opt['ref'], self.opt['files'][0], self.opt['out'])


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
            optlist, args = getopt (argv[1:], "r:o:c:", ['ref=', 'out=', 'color='])
        except GetoptError:
            print ('Invalid command line options')
            return

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-r' or o == '--ref':
                self.opt['ref'] = a
            elif o == '-o' or o == '--out':
                self.opt['out'] = a
            elif o == '-c' or o == '--color':
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
