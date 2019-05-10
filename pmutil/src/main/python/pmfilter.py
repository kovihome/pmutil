#!/usr/bin/env python3
#
# PmUtils/pmfilter
#

'''
Created on Apr 28, 2019

@author: kovi
'''
import getopt
import sys


def loadRefCatalog(refFileName):

    refCat = { 'header': [], 'cat': [] }

    ref = open(refFileName, 'r')
    if ref:
        headerLine = ref.readline()
        refCat['header'] = headerLine.split()

        for refLine in ref:
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


def getHeaderPos(cat, headerName):
    index = 0
    for header in cat['header']:
        if header == headerName:
            return index
        index = index + 1
    return None


def matchCatalogs(refCat, pmCat):

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
            if dist2 < ref[rlen]:
                ref[rlen] = dist2
                ref[rlen + 1] = pmi
        pmi = pmi + 1

    # print (refCat)


refHeaders = ['AUID', 'RA', 'RA_DEG', 'DEC', 'DEC_DEG', 'LABEL', 'MAG_V', 'ERR_V', 'B-V', 'ERR_BV']
pmHeaders = ['NUMBER', 'MAG_ISOCOR', 'MAGERR_ISOCOR', 'MAG_BEST', 'MAGERR_BEST', 'ALPHA_J2000', 'DELTA_J2000']


def dumpResult(refCat, pmCat, outFileName):

    outf = open(outFileName, 'w')
    for h in refHeaders:
        outf.write(h.ljust(15))
    for h in pmHeaders:
        outf.write(h.ljust(15))
    outf.write('\n')

    rlen = len(refCat['header'])
    for ref in refCat['cat']:
        pm = pmCat['cat'][ref[rlen + 1]]

        for h in refHeaders:
            pos = getHeaderPos(refCat, h)
            outf.write(ref[pos].ljust(15))
        for h in pmHeaders:
            pos = getHeaderPos(pmCat, h)
            outf.write(pm[pos].ljust(15))
        outf.write('\n')

    outf.close()
    print ('outpuf file: ', outFileName)


commandLineOptions = {
    'ref' : None,
    'out' : None
    }


def processCommands():
    try:
        optlist, args = getopt.getopt (sys.argv[1:], "r:o:", ['--ref', '--out'])
    except getopt.GetoptError:
        print ('Invalid command line options')
        return

    for o, a in optlist:
        if a[:1] == ':':
            a = a[1:]
        elif o == '-r':
            commandLineOptions['ref'] = a
        elif o == '-o':
            commandLineOptions['out'] = a

    commandLineOptions['files'] = args
    if not commandLineOptions['out']:
        commandLineOptions['out'] = args[0] + '.refout'


if __name__ == '__main__':

    processCommands()

    # if commandLineOptions['ref']:
    refCat = loadRefCatalog(commandLineOptions['ref'])

    pmCat = loadPmCatalog(commandLineOptions['files'][0])

    matchCatalogs(refCat, pmCat)

    dumpResult(refCat, pmCat, commandLineOptions['out'])

