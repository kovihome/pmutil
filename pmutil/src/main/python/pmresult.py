#!/usr/bin/env python3
#
# PmUtils/pmresult
#
import numpy
import subprocess
import glob
import datetime

'''
Created on Apr 28, 2019

@author: kovi
'''
import getopt
import sys


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
fieldAuid = 'AUID'


def calculateMgs(refCat, comp, alsoInteresting):
    idPos = getHeaderPos(refCat, fieldAuid)
    miPos = getHeaderPos(refCat, fieldMgInstrumental)
    mvPos = getHeaderPos(refCat, fieldMgTrue)
    mi = []
    mv = []
    p = [1.0, 0.0]
    for pm in refCat['cat']:
        if pm[mvPos] != '-':
            mi.append(float(pm[miPos]))
            mv.append(float(pm[mvPos]))
            if comp != None and pm[idPos] == comp:
                p[1] = float(pm[mvPos]) - float(pm[miPos])

    if comp == None:
        p = numpy.polyfit(mi, mv, 1)
        # print ('polyfit result:', p)

    result = []
    for pm in refCat['cat']:
        if pm[mvPos] == '-' or pm[idPos] in alsoInteresting:
            mv = p[0] * float(pm[miPos]) + p[1]
            pm[mvPos] = mv
            result.append(pm)
    # print (result)
    return result


def getDateObs(fileName):
    seqFileName = fileName.split('.')[0]
    if seqFileName.find('/'):
        seqFileName = seqFileName.rpartition('/')[2]
    seqFileName = 'Sequence/' + seqFileName + '.fits'

    p = subprocess.Popen('fiheader --read DATE-OBS ' + seqFileName, stdout = subprocess.PIPE, shell = True)
    (output, errno) = p.communicate()
    p.wait()

    return output.decode('ascii').strip()


def reportResult(allResults, refCat, outFileName):

    # print(allResults)

    idPos = getHeaderPos(refCat, 'AUID')
    labelPos = getHeaderPos(refCat, 'LABEL')
    mvPos = getHeaderPos(refCat, 'MAG_V')
    errPos = getHeaderPos(refCat, 'MAGERR_BEST')

    if outFileName == None:
        outFileName = 'result.txt'
    r = open(outFileName, 'w')

    r.write('DATE-OBS'.ljust(20))
    result = allResults[list(allResults.keys())[0]]
    for res in result:
        auid = res[idPos] + " " + res[labelPos]
        r.write(auid.ljust(40))
    r.write('\n')

    r.write(' ' * 20)
    for j in range(len(result)):
        r.write('MAG_V'.ljust(20))
        r.write('MAGERR_V'.ljust(20))
    r.write('\n')

    for dateObs in allResults.keys():
        r.write(dateObs.ljust(20))
        result = allResults[dateObs]
        for res in result:
            r.write(str(res[mvPos]).ljust(20))
            r.write(str(res[errPos]).ljust(20))
        r.write('\n')
    r.close()


commandLineOptions = {
    'out' : None,  # output result file name, if None, 'result.txt' will be created
    'comp': None,  # comparision star id, if None, ensemble will proceed
    'also': [],  # id list of also interesting stars
    'rpt' : False
    }


def processCommands():
    try:
        optlist, args = getopt.getopt (sys.argv[1:], "o:a:c:r", ['--out', '--also', '--comp', '--report'])
    except getopt.GetoptError:
        print ('Invalid command line options')
        return

    for o, a in optlist:
        if a[:1] == ':':
            a = a[1:]
        elif o == '-o':
            commandLineOptions['out'] = a
        elif o == '-c':
            commandLineOptions['comp'] = a
        elif o == '-a':
            commandLineOptions['also'] = a.split(',')
        elif o == '-r':
            commandLineOptions['rpt'] = True

    commandLineOptions['files'] = args
#    if not commandLineOptions['out']:
#        commandLineOptions['out'] = args[0] + '.refout'


def jd(dateObs):
    # d = datetime.fromisoformat(dateObs) # python 3.7+
    dt = dateObs.split('T')
    dd = dt[0].split('-')
    tt = dt[1].split(':')
    d = datetime.datetime(int(dd[0]), int(dd[1]), int(dd[2]), int(tt[0]), int(tt[1]), int(tt[2]))

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


def reportForAAVSO(allResults, refCat, comp):
    '''
    header:
    #TYPE=Extended
    #OBSCODE=?
    #SOFTWARE=
    #DELIM=,
    #DATE=JD
    #OBSTYPE=DSLR
    
    data:
    STARID: The star's identifier. It can be the AAVSO Designation, the AAVSO Name, or the AAVSO Unique Identifier, but NOT more than one of these. Limit: 30 characters.
    DATE: The date of the observation, in the format specified in the DATE parameter. Limit: 16 characters.
    MAGNITUDE: The magnitude of the observation. Prepend with < if a fainter-than.  A dot is required (e.g. "9.0" rather than "9"). Limit: 8 characters.
    MAGERR: Photometric uncertainty associated with the variable star magnitude. If not available put "na". Limit: 6 characters.
    FILTER: The filter used for the observation. This can be one of the following letters (in bold):
        B: Johnson B
        V: Johnson V
        R: Cousins R
        TG: Green Filter (or Tri-color green). This is commonly the "green-channel" in a DSLR or color CCD camera. These observations use V-band comp star magnitudes.
        TB: Blue Filter (or Tri-color blue). This is commonly the "blue-channel" in a DSLR or color CCD camera. These observations use B-band comp star magnitudes.
        TR: Red Filter (or Tri-color red). This is commonly the "red-channel" in a DSLR or color CCD camera. These observations use R-band comp star magnitudes.
    TRANS: YES if transformed using the Landolt Standards or those fields that contain secondary standards, or NO if not. Document the method used to transform in the "NOTES" section.
    MTYPE: Magnitude type. STD if standardized (Click here for definition of standardized) or DIF if differential (very rare). If you are currently using ABS for 'absolute' we will still accept it. Differential requires the use of CNAME.
    CNAME: Comparison star name or label such as the AUID (much preferred) or the chart label for the comparison star used. If not present, use "na". Limit: 20 characters.
    CMAG: Instrumental magnitude of the comparison star. If "ensemble" see example below. If not present, use "na". Limit: 8 characters.
    KNAME: Check star name or label such as the AUID (much preferred) or the chart label for the check star. If not present, use "na". Limit: 20 characters.
    KMAG: Instrumental magnitude of the check star. If "ensemble" see example below. If not present, use "na".Limit: 8 characters.
    AIRMASS: Airmass of observation Limit 7 characters - entry will be truncated if longer than that. If not present, use "na".
    GROUP: Grouping identifier (maximum 5 characters). It is used for grouping multiple observations together, usually an observation set that was taken through multiple filters. It makes it easier to retrieve all magnitudes from a given set in the database, say, if someone wanted to form color indices such as (B-V) with them. If you are just doing time series, or using the same filter for multiple stars, etc., just set GROUP to "na." For cases where you want to group observations, GROUP should be an integer, identical for all observations in a group, and unique for a given observer for a given star on a given Julian Date. Limit: 5 characters.
    CHART: Please use the sequence ID you will find written in Bold print near the top of the photometry table in a sentence that reads "Report this sequence as [ID] in the chart field of your observation report." If a non–AAVSO sequence was used, please describe it as clearly as possible. Limit: 20 characters.
    NOTES: Comments or notes about the observation. If no comments, use "na". This field has a maximum length of several thousand characters, so you can be as descriptive as necessary. One convention for including a lot of information as concisely as possible is to use subfields in the format |A=B; the '|' character is the separator, A is a keyvalue name and B is its value. If you need an alternative delimiter from '|', use it but preceed the first instance with "DELIM=". Using this mechanismnyou can document your transform process in more detail. Here is an example as used by TransformApplier:
        5 records aggregated|VMAGINS=-7.244|VERR=0.006|CREFMAG=13.793|CREFERR=0.026|
        KREFMAG=14.448|KREFERR=0.021|VX=1.1501|CX=1.1505|KX=1.1500|Tv_bv=0.0090|Tv_bvErr=0.0100|
        TAver=2.47
    '''

    # idPos = getHeaderPos(refCat, 'AUID')
    labelPos = getHeaderPos(refCat, 'LABEL')
    mvPos = getHeaderPos(refCat, 'MAG_V')
    errPos = getHeaderPos(refCat, 'MAGERR_BEST')

    outFileName = 'result.extended.aavso'
    r = open(outFileName, 'w')

    r.write('#TYPE=Extended\n')
    r.write('#OBSCODE=?\n')
    r.write('#SOFTWARE=pmresult.py\n')
    r.write('#DELIM=,\n')
    r.write('#DATE=JD\n')
    r.write('#OBSTYPE=DSLR\n')

    for dateObs in allResults.keys():
        result = allResults[dateObs]
        for res in result:

            starid = res[labelPos].replace('_', ' ')
            date = jd(dateObs)
            magnitude = res[mvPos]
            magerr = res[errPos]
            flt = 'TG'  # TODO: select appropriate filter for multicolor photometry
            trans = 'NO'
            mtype = 'STD'
            cname = 'ENSEMBLE' if comp == None else comp
            cmag = 'na'  # TODO: set comp star mg if not ensemble
            kname = 'na'  # TODO: specify ensemble reference comp star when ensemble, or set a check star
            kmag = 'na'  # TODO: set measured mg of check star
            airmass = 'na'  # TODO: calculate airmass
            group = 'na'
            chart = 'na'
            notes = 'na'
            r.write(('%s,%7.4f,%4.3f,%4.3f,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n') % (starid, date, float(magnitude), float(magerr), flt, trans, mtype, cname, cmag, kname, kmag, airmass, group, chart, notes))

    r.close()


if __name__ == '__main__':

    processCommands()

    if not commandLineOptions['files']:
        f = glob.glob('Phot/Seq_*.fits.cat.pm')
        f.sort()
        commandLineOptions['files'] = f

    allResults = {}

    for fileName in commandLineOptions['files']:

        refCat = loadCatalog(fileName)

        dateObs = getDateObs(fileName)

        result = calculateMgs(refCat, commandLineOptions['comp'], commandLineOptions['also'])

        allResults[dateObs] = result

    reportResult(allResults, refCat, commandLineOptions['out'])

    if commandLineOptions['rpt']:
        reportForAAVSO(allResults, refCat, commandLineOptions['comp'])

