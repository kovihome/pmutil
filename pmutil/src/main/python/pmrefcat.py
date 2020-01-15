#!/usr/bin/env python3
#
# PmUtils/pmrefcat
#
'''
Created on Jan 1, 2020

@author: kovi

1. load photometry for comp stars from AAVSO VSP fotometry table. 
'''

from sys import argv
from getopt import getopt, GetoptError
from urllib import request
import json
import xmltodict
import re

defMgLimit = 18.0
defFov = 60

auidOnly = True

Color_Off = '\033[0m'  # Text Reset
BRed = "\033[1;31m"  # Red
BGreen = '\033[1;32m'  # Green
Blue = '\033[0;34m'  # Blue
BCyan = '\033[1;36m'  # Cyan


def printError(s):
    print(BRed + "Error: " + s + Color_Off)


def printWarning(s):
    print(BGreen + "Warning: " + s + Color_Off)


def printInfo(s):
    print(BCyan + s + Color_Off)


def printTitle():
    print()
    print(BGreen + "ppl-refcat, version 0.1.0 " + Color_Off)
    print(Blue + "Create reference catalog for photometry." + Color_Off)


def hexa2deg(s):
    r = s.split(':')
    d = float(r[0])
    m = float(r[1])
    s = float(r[2])
    return d + m / 60.0 + s / 3600.0


def deg2hexa(d):
    dd = int(d)
    r = (d - float(dd)) * 60.0
    mm = int(r)
    ss = (r - float(mm)) * 60.0
    return "%02d:%02d:%04.1f" % (dd, mm, ss)


def loadVspPhotometryData(outFile, objectName, ra = None, dec = None, fov = defFov):
    '''
    star*: name of the star to plot. You must provide EITHER the star parameter OR ra and dec parameters (see below)
    ra*: right ascension
    dec*: declination
    fov*: field of view, in arcminutes.
    maglimit*: magnitude limit for the chart; stars fainter than the maglimit will not be plotted.
    other: other variable stars to mark on the chart. Omit or leave blank for none, "gcvs" for GCVS variables, "all" for all variables
    '''
    vspUrl = "https://www.aavso.org/apps/vsp/api/chart/?format=json&star=%s&fov=%d&maglimit=%4.1f&other=all" % (objectName.replace(' ', '+'), fov, defMgLimit)
    f = request.urlopen(vspUrl)
    respJson = f.read().decode('UTF-8')

    resp = json.loads(respJson)

    # print("Star:", resp['star'])
    # print("AUID:", resp['auid'])
    # print("Coord:", resp['ra'], resp['dec'])
    # print("ChartID:", resp['chartid'])
    # print("Photometry table:")
    for pm in resp['photometry']:

        # 'V', 'B', 'Rc', 'Ic', 'J', 'H', 'K'
        bands = {}
        for b in pm['bands']:
            bands[b['band']] = b

        auid = pm['auid']
        role = "C"
        ra = pm['ra']
        dec = pm['dec']
        raDeg = "%10.8f" % (hexa2deg(ra) * 15.0)
        if not dec.startswith('-'):
            dec = "+" + dec
        decDeg = "%+10.8f" % (hexa2deg(dec))
        if 'B' in bands:
            mgB = str(bands['B']['mag'])
            mgerrB = str(bands['B']['error'])
            if mgerrB == 'None':
                mgerrB = "-"
        else:
            mgB = "-"
            mgerrB = "-"
        if 'V' in bands:
            mgV = str(bands['V']['mag'])
            mgerrV = str(bands['V']['error'])
            if mgerrV == 'None':
                mgerrV = "-"
        else:
            mgV = "-"
            mgerrV = "-"
        if 'Rc' in bands:
            mgR = str(bands['Rc']['mag'])
            mgerrR = str(bands['Rc']['error'])
            if mgerrR == 'None':
                mgerrR = "-"
        else:
            mgR = "-"
            mgerrR = "-"
        label = str(pm['label']).replace(' ', '_')

        s = auid.ljust(13) + role.ljust(5) + ra.ljust(12) + raDeg.ljust(15) + dec.ljust(12) + decDeg.ljust(14) + mgB.ljust(7) + mgerrB.ljust(8) + mgV.ljust(7) + mgerrV.ljust(8) + mgR.ljust(7) + mgerrR.ljust(8) + label
        if outFile != None:
            outFile.write(s + "\n")
        else:
            print(s)


def loadVsxCatalogData(outFile, objectName, ra = None, dec = None, fov = defFov):
    '''
    &coords   The central J2000 RA/DEC coordinates for a radius search, expressed sexagesimally (by default), or in decimal degrees if format is set to d. Northern hemisphere coordinates must have the plus () sign URL-encoded as %2B. Space characters between all other figures must be replaced with the URL whitespace character (). The order argument (which see) must also be included in the query string with its value set to 9 in order to prompt VSX to display distances from the central coordinates in the results listing. Default is empty string (no radius search).
    &ident    Object identification for name searches. Space characters must be replaced with the URL whitespace character (+). Other special characters may also need to be URL-encoded. Default is empty string (no name search).
    &format   Explicit specification for format of coords. For sexagesimal, this value should be s. For decimal degrees, this value should be d. Default is s (sexagesimal).
    &geom     The geometry for central coordinate-based searches. For radius searches, this value should be r. For box searches, this value should be b. Default is r (radius search).
    &size     For box searches (geom=b), the width of the box. For radius searches (geom=r), the radius of the circle. Expressed in the units specified by unit (see next). Default is 10.0.
    &unit     The unit of measurement used for the value given in size (see above). For arc degrees, this value should be 1. For arc minutes, this value should be 2. For arc seconds, this value should be 3. Default is 2 (arc minutes).
    '''
    # https://www.aavso.org/vsx/index.php?view=results.csv&ident=R+Car
    vsxUrl = "https://www.aavso.org/vsx/index.php?view=api.object&format=json&ident=%s" % (objectName.replace(' ', '+'))
    # print(vsxUrl)
    f = request.urlopen(vsxUrl)
    respJson = f.read().decode('UTF-8')
    resp = json.loads(respJson)
    # print(resp)
    # print("Object:", resp['VSXObject']['Name'], "RA:", resp['VSXObject']['RA2000'], "Dec:", resp['VSXObject']['Declination2000'])

    ra = deg2hexa(float(resp['VSXObject']['RA2000']) / 15.0)
    dec = deg2hexa(float(resp['VSXObject']['Declination2000']))
    coords = ra + "+" + dec
    vsxUrl = "https://www.aavso.org/vsx/index.php?view=query.votable&format=d&coords=%s&size=%d&unit=2&geom=b" % (coords, fov)
    # https://www.aavso.org/vsx/index.php?view=query.votable&coords=19:54:17+32:13:08&format=d&size=60&unit=2
    # print(vsxUrl)
    f = request.urlopen(vsxUrl)
    respVOTable = f.read().decode('UTF-8')
    # resp = json.loads(respJson)
    # print(respVOTable)
    votable = xmltodict.parse(respVOTable)
    # print(votable)
    # print(votable['VOTABLE']['RESOURCE']['TABLE'].keys())
    # print(votable['VOTABLE']['RESOURCE']['TABLE']['FIELD'])
    # print(votable['VOTABLE']['RESOURCE']['TABLE']['DATA']['TABLEDATA'])
    trs = votable['VOTABLE']['RESOURCE']['TABLE']['DATA']['TABLEDATA']['TR']
    nrUnknown = 1
    for tr in trs:
        # print(tr['TD'])
        # [None, 'PS1-3PI J185203.12+325154.5', 'Lyr', '283.01304000,32.86516000', 'RRAB', '16.810', 'r', '0.400', 'r', '57000.85600', None, '0.839122', None, None, None]
        # print("AUID:", tr['TD'][0], "Coord:", tr['TD'][3], "LABEL", tr['TD'][1])

        # #0.00    Variable    AY Lyr                000-BCD-108    18 44 26.69 +37 59 51.9    Lyr    UGSU    0.0733    12.5 - 18.4 B
        # distance, varstate, label(1), auid(0), coords(3), constell(2), type(4), period(11), max(5,6) - min(7,8)
        # print("   Constell:", tr['TD'][2], "Type:", tr['TD'][4], "Period:", tr['TD'][11], "Mg:", tr['TD'][5] + tr['TD'][6], "-", tr['TD'][7] + tr['TD'][8])

        auid = tr['TD'][0]
        if auid == None:
            if auidOnly:
                continue
            auid = '999-VAR-%03d' % (nrUnknown)
            nrUnknown = nrUnknown + 1
            tr['TD'][0] = auid
        role = 'V'
        raDeg, decDeg = tr['TD'][3].split(',')
        ra = deg2hexa(float(raDeg) / 15.0)
        dec = deg2hexa(float(decDeg))
        if not dec.startswith('-'):
            dec = "+" + dec
            decDeg = "+" + decDeg
        label = tr['TD'][1].replace(' ', '_')

        s = auid.ljust(13) + role.ljust(5) + ra.ljust(12) + raDeg.ljust(15) + dec.ljust(12) + decDeg.ljust(14) + "-".ljust(7) + "-".ljust(8) + "-".ljust(7) + "-".ljust(8) + "-".ljust(7) + "-".ljust(8) + label
        if outFile != None:
            outFile.write(s + "\n")
        else:
            print(s)

    nrUnknown = 1
    for tr in trs:

        auid = tr['TD'][0]
        if auid == None:
            if auidOnly:
                continue
            auid = '999-VAR-%03d' % (nrUnknown)
            nrUnknown = nrUnknown + 1
        role = 'V'
        raDeg, decDeg = tr['TD'][3].split(',')
        ra = deg2hexa(float(raDeg) / 15.0)
        dec = deg2hexa(float(decDeg))
        label = tr['TD'][1]
        constell = tr['TD'][2]
        vartype = tr['TD'][4]
        if vartype == None:
            vartype = "-"
        period = tr['TD'][11]
        maxmg = tr['TD'][5] + tr['TD'][6]
        minmg = tr['TD'][7] + tr['TD'][8]

        s = "#" + label.ljust(29) + auid.ljust(14) + ra.ljust(15) + dec.ljust(14) + constell.ljust(4) + vartype.ljust(11) + period.ljust(11) + maxmg + " - " + minmg
        if outFile != None:
            outFile.write(s + "\n")
        else:
            print(s)

def usage():
    print()
    print(BGreen + "ppl-refcat, version 0.1.0 " + Color_Off)
    print()
    print("Usage: ppl-refcat [OPTIONS]... CATALOG_FILE_NAME")
    print("Create reference catalog for photometry.")
    print()
    print("Mandatory arguments to long options are mandatory for short options too.")
    print("  -c,  --coords ra,decl       coordinates of the center of reference frame. Valid format is 12:34:56.7,-12:34:56.7")
    print("  -s,  --source catalog       source catalog for field stars")
    print("  -o,  --object object_name   object (variable star) name")
    print("  -f,  --field size           field size in arcmin")
    print("  -h,  --help                 print this page")


commandLineOptions = {
    'coords' : None,  # coordinates of reference field
    'source': None,  # source catalog of field stars
    'object': None,  # object (variable star) name
    'field' : defFov,  # reference field size in arcmins
    'file'  : None  # reference catalog file name
    }


def processCommands():
    try:
        optlist, args = getopt (argv[1:], "c:s:o:f:h", ['--coord', '--source', '--object', '--field', '--help'])
    except GetoptError:
        print ('Invalid command line options')
        exit(1)

    for o, a in optlist:
        if a[:1] == ':':
            a = a[1:]
        elif o == '-c':
            commandLineOptions['coords'] = a
        elif o == '-s':
            commandLineOptions['source'] = a
        elif o == '-o':
            commandLineOptions['object'] = a.replace('_', ' ')
        elif o == '-f':
            if a.isdigit():
                commandLineOptions['field'] = int(a)
            else:
                print("Invalid field size parameter: %s. Use default 60 arcmin instead." % (a))
        elif o == '-h':
            usage()
            exit(0)

    if len(args) > 0:
        commandLineOptions['file'] = args[0]

#     if commandLineOptions['file'] == None:
#         printError("Catalog file name is mandatory.")
#         exit(1)

    if commandLineOptions['object'] == None and commandLineOptions['coords'] == None:
        printError("Either object name (-o) or coordinates (-c) must be given.")
        exit(1)

    if commandLineOptions['coords'] != None:
        c = commandLineOptions['coords'].split(',', maxsplit = 2)
        ok = True
        if len(c) != 2:
            ok = False
        else:
            if re.fullmatch("(\d){2}:(\d){2}:(\d){2}(\.\d)*", c[0]) == None:
                ok = False
            if re.fullmatch("[+-]*(\d){2}:(\d){2}:(\d){2}(\.\d)*", c[1]) == None:
                ok = False
        if not ok:
            printError("Invalid coordinate format")
            exit(1)

        commandLineOptions['ra'] = c[0]
        commandLineOptions['dec'] = c[1]


catalogHeader = "AUID         ROLE RA          RA_DEG         DEC         DEC_DEG       MAG_B  ERR_B   MAG_V  ERR_V   MAG_R  ERR_R   LABEL"

if __name__ == '__main__':

    processCommands()

    printTitle()
    print()

    outFile = None
    if commandLineOptions['file'] != None:
        outFile = open(commandLineOptions['file'], 'w')
        outFile.write(catalogHeader + "\n")
    else:
        print(catalogHeader)

    if commandLineOptions['object'] != None:

        loadVspPhotometryData(outFile, commandLineOptions['object'], fov = commandLineOptions['field'])

        loadVsxCatalogData(outFile, commandLineOptions['object'], fov = commandLineOptions['field'])

    else:

        loadVspPhotometryData(outFile, ra = commandLineOptions['ra'], dec = commandLineOptions['dec'], fov = commandLineOptions['field'])

        loadVsxCatalogData(outFile, ra = commandLineOptions['ra'], dec = commandLineOptions['dec'], fov = commandLineOptions['field'])

    if outFile != None:
        outFile.close()

    printInfo("Reference catalog file %s created." % (commandLineOptions['file']))

# end main.
