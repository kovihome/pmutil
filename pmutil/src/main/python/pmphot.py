#!/usr/bin/env python3
#
# PmUtils/pmphot
#
'''
Created on Dec 28, 2019

@author: kovi
'''
from numpy import median, sqrt, polyfit
from subprocess import Popen, PIPE
from getopt import getopt, GetoptError
from sys import argv
from datetime import datetime
from glob import glob


def jd(dateObs):
    # d = datetime.fromisoformat(dateObs) # python 3.7+
    dt = dateObs.split('T')
    dd = dt[0].split('-')
    tt = dt[1].split(':')
    d = datetime(int(dd[0]), int(dd[1]), int(dd[2]), int(tt[0]), int(tt[1]), int(tt[2]))

    if d.month == 1 or d.month == 2:
        d.year -= 1
        d.month += 12
    if d.month > 2:
        year = d.year
        month = d.month
    if year >= 1582:
        a = int(year / 100)
        b = 2 - a + int(a / 4)
    else:
        b = 0

    day = d.day + d.hour / 24 + d.minute / (24 * 60) + d.second / (24 * 3600)
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + b - 1524.5
    return jd


def loadCatalog(refFileName):

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


fieldMgInstrumental = 'MAG_BEST'
fieldMgTrue = 'MAG_V'
fieldMgErrInstrumental = 'MAGERR_BEST'
fieldMgErrTrue = 'ERR_V'
fieldAuid = 'AUID'
fieldRole = 'ROLE'


def calculateMgsRobustAveraging(refCat, color):
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
                ev = 0.2    # TODO: what to do, if no mg error value?
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

    result = []
    for pm in refCat['cat']:
        if pm[rolePos] == 'V' or pm[rolePos] == 'VF':
            mv = float(pm[miPos]) + z
            pm[mvPos] = mv
            pm[evPos] = sqrt(ez_2)
            result.append(pm)
    # print (result)
    return result


def calculateMgsLinearFit(refCat, color):
    '''
    Calculate magnitudes with naive linear fit
    '''

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
            er.append(sqrt(ei * ei + ev * ev))

    p = polyfit(mi, mv, 1, w=er)
    # print ('polyfit result:', p)

    result = []
    for pm in refCat['cat']:
        if pm[rolePos] == 'V':
            mv = p[0] * float(pm[miPos]) + p[1]
            pm[mvPos] = mv
            result.append(pm)
    # print (result)
    return result


def calculateMgsComparision(refCat, color, comp):
    stdcolor = color[:1].upper()
    if stdcolor == 'G':
        stdcolor = 'V'

    idPos = getHeaderPos(refCat, fieldAuid)
    mvPos = getHeaderPos(refCat, 'MAG_' + stdcolor)
    miPos = getHeaderPos(refCat, fieldMgInstrumental)
    evPos = getHeaderPos(refCat, 'ERR_' + stdcolor)
    eiPos = getHeaderPos(refCat, fieldMgErrInstrumental)
    rolePos = getHeaderPos(refCat, fieldRole)

    p = 0.0
    ep = 0.0
    for pm in refCat['cat']:
        if pm[rolePos] == 'C' and pm[mvPos] != '-':
            if comp != None and pm[idPos] == comp:
                p = float(pm[mvPos]) - float(pm[miPos])
                ep = max(float(pm[evPos]), float(pm[eiPos]))

    result = []
    for pm in refCat['cat']:
        if pm[rolePos] == 'V':
            mv = float(pm[miPos]) + p
            pm[mvPos] = mv
            pm[evPos] = max(float(pm[eiPos]), ep)
            result.append(pm)
    # print (result)
    return result


def getDateObs(fileName):
    '''
    fileName: path/[Seq_nnn|Combined]_c.fits.cat.cat
    '''
    astFileName = fileName.split('.')[0]
#    if seqFileName.find('/'):
#        seqFileName = seqFileName.rpartition('/')[2]
    astFileName = astFileName + '.ast.fits'

    p = Popen('fiheader --read DATE-OBS ' + astFileName, stdout = PIPE, shell = True)
    (output, errno) = p.communicate()
    p.wait()

    return output.decode('ascii').strip()


def reportResult(allResults, refCat, outFileName, color):

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
    rolePos = getHeaderPos(refCat, 'ROLE')

    if outFileName == None:
        outFileName = 'result.txt'
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
    for dateObs in allResults.keys():
        jdObs = jd(dateObs)
        jds = "%12.4f" % jdObs
        result = allResults[dateObs]
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


commandLineOptions = {
    'out' : None,  # output photometry file name, if None, '*.pm' will be created
    'comp': None,  # comparision star id, if None, ensemble will proceed
    'color' : 'Gi',  # photometry band
    'std': False  # use std coeffs
    }


def processCommands():
    try:
        optlist, args = getopt (argv[1:], "c:so:p:", ['--color', '--std', '--out', '--comp'])
    except GetoptError:
        print ('Invalid command line options')
        return

    for o, a in optlist:
        if a[:1] == ':':
            a = a[1:]
        elif o == '-c':
            commandLineOptions['color'] = a
        elif o == '-s':
            commandLineOptions['std'] = True
        elif o == '-m':
            commandLineOptions['makestd'] = True
        elif o == '-c':
            commandLineOptions['comp'] = a
        elif o == '-o':
            commandLineOptions['out'] = a

    commandLineOptions['files'] = args


if __name__ == '__main__':

    processCommands()

    if not commandLineOptions['files']:
        f = glob('Phot/Seq_*.fits.cat.pm')
        f.sort()
        commandLineOptions['files'] = f

    allResults = {}

    for fileName in commandLineOptions['files']:

        refCat = loadCatalog(fileName)

        dateObs = getDateObs(fileName)

        if commandLineOptions['comp']:
            result = calculateMgsComparision(refCat, commandLineOptions['color'], commandLineOptions['comp'])

        else:
            result = calculateMgsRobustAveraging(refCat, commandLineOptions['color'])
#           result = calculateMgsLinearFit(refCat, commandLineOptions['color'])

        allResults[dateObs] = result

    reportResult(allResults, refCat, commandLineOptions['out'], commandLineOptions['color'])

# end __main__
