#!/usr/bin/env python3
#
# PmUtils/pmbase
#
'''
Created on Jan 1, 2020

@author: kovi
'''

from sys import argv
from os import getenv, makedirs
from datetime import datetime
from glob import glob
from astropy.io import fits
import subprocess

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


def loadPplSetup():
    pplSetup = {}
    userhome = getenv("HOME")
    f = open(userhome + "/bin/ppl-setup")
    pplSetup['HOME'] = userhome
    for line in f:
        if not line.startswith('#') and line.strip() != "":
            r = line.split('=')
            value = r[1].rstrip()[1:-1]
            if '$' in value:
                for k in pplSetup.keys():
                    if '$' + k in value:
                        value = value.replace('$' + k, pplSetup[k])
            pplSetup[r[0].strip()] = value
    return pplSetup


def invoke(cmd):
    r = cmd.split()
    a = subprocess.check_output(r)
    return a.decode('ascii')[:-1]


def hexa2deg(s):
    r = s.split(':')
    sign = 1.0
    if r[0].startswith('-'):
        sign = -1.0
        r[0] = r[0][1:]
    d = float(r[0])
    m = float(r[1])
    s = float(r[2])
    return sign * (d + m / 60.0 + s / 3600.0)


def deg2hexa(d):
    sign = ''
    if d < 0.0:
        d = -d
        sign = '-'
    dd = int(d)
    r = (d - float(dd)) * 60.0
    mm = int(r)
    ss = (r - float(mm)) * 60.0
    return "%s%02d:%02d:%04.1f" % (sign, dd, mm, ss)


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

def getPair(s, delim = ','):
    ss = s.split(delim)
    return [ss[0].strip(), ss[1].strip()]


def determineCoordsFromImage(imageFileName):
    configFolder = pplSetup["PMLIB"]
    if configFolder == None:
        configFolder = "$HOME/.pmlib"
    SOLVE_ARGS = "-O --config " + configFolder + "/astrometry.cfg --use-sextractor --sextractor-path sextractor -r -y -p"

    # /usr/local/astrometry/bin/solve-field $SOLVE_ARGS -D $2 -N $AST_FILE $f

    makedirs("temp", exist_ok = True)

    invoke("cp " + imageFileName + " temp/src.fits")

    astFolder = pplSetup["AST_BIN_FOLDER"]
    if astFolder == None:
        astFolder = "/usr/local/astrometry/bin"

    astResult = invoke(astFolder + "/solve-field " + SOLVE_ARGS + " -D temp -N temp/ast.fits temp/src.fits")

    r = astResult.split('\n')
    sCoords = ["0.0", "0.0"]
    for line in r:
        if line.startswith("Field center"):
            if line.find("RA,Dec") > -1:
                sCoords = getPair(line.split('(')[2].split(')')[0])
                printInfo("Image center: %s %s" % (sCoords[0], sCoords[1]))
            else:
                sCoordsSexa = getPair(line.split('(')[2].split(')')[0])
                printInfo("Image center: %s %s" % (sCoordsSexa[0], sCoordsSexa[1]))
        elif line.startswith("Field size"):
            sFieldSize = getPair(line.split(':')[1].split('a')[0], 'x')
            printInfo("Image size: %s' x %s'" % (sFieldSize[0], sFieldSize[1]))

    return sCoords


def saveCommand(basePath, argv, cmdName):
    for j in range(len(argv)):
        if " " in argv[j]:
            argv[j] = "\"" + argv[j] + "\""
    cmd = " ".join(argv)
    cmd = cmd[cmd.find("ppl-"):]

    path = ""
    j = basePath.rfind("/")
    if j != -1:
        path = basePath[:j+1]
	
    f = open(path + cmdName, "w+")
    f.write(cmd + "\n")
    f.close()


def discoverFolders(baseFolder, seqFolderName):
    seqFolders = []
    if baseFolder != None:
        seqFolders = glob('*' + baseFolder + '*/' + seqFolderName + '*')
    else:
        seqFolders = glob(seqFolderName + '*')
    return seqFolders


def getFitsHeaders(fitsFileName, headers):
    try:
        hdul = fits.open(fitsFileName)
    except Exception as e:
        print('Exception:',e)
        return None

    fitsHeaders = {}
    for h in headers:
        fitsHeaders[h] = hdul[0].header[h]

    hdul.close()
    return fitsHeaders

def getFitsHeader(fitsFileName, header):
    hdrs = getFitsHeaders(fitsFileName, [header])
    if hdrs:
        return hdrs[header]
    else:
        return {}




# end pmbase.
