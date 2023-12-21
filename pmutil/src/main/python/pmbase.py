#!/usr/bin/env python3
#
# PmUtils/pmbase
#
"""
Created on Jan 1, 2020

@author: kovi
"""

from os import getenv, makedirs, chmod, getcwd
from os.path import isfile, exists
from glob import glob
from math import sqrt
import subprocess
import stat

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground
import matplotlib.pyplot as plt

Color_Off = '\033[0m'  # Text Reset
BRed = "\033[1;31m"  # Red
BGreen = '\033[1;32m'  # Green
Blue = '\033[0;34m'  # Blue
BCyan = '\033[1;36m'  # Cyan
Color_Yellow = '\033[0;93m'  # Light yellow

pmlog = None

setup = None


class Logger:
    LOG_MODE_SILENT = 0
    LOG_MODE_ERROR = 1
    LOG_MODE_WARNING = 2
    LOG_MODE_INFO = 3
    LOG_MODE_PRINT = 4
    LOG_MODE_DEBUG = 5

    logFileName = None
    logMode = LOG_MODE_INFO

    def __init__(self, logFileName=None, mode=LOG_MODE_INFO):
        global pmlog
        self.logFileName = logFileName
        self.logMode = mode
        pmlog = self

    def write(self, text):
        if self.logFileName:
            f = open(self.logFileName, "a+")
            f.write(text + '\n')
            f.close()

    def log(self, prefix, text, mode, color):
        s = (color if color else '') + (prefix if prefix else '') + text + Color_Off
        if self.logMode >= mode:
            print(s)
        self.write(s)

    def error(self, text):
        self.log("Error: ", text, self.LOG_MODE_ERROR, BRed)

    #        s = BRed + "Error: " + text + Color_Off
    #        if self.logMode >= self.LOG_MODE_ERROR:
    #            print(s)
    #        self.write(s)

    def warning(self, text):
        self.log("Warning: ", text, self.LOG_MODE_WARNING, BGreen)

    #        s = BGreen + "Warning: " + text + Color_Off
    #        if self.logMode >= self.LOG_MODE_WARNING:
    #            print(s)
    #        self.write(s)

    def info(self, text):
        self.log(None, text, self.LOG_MODE_INFO, BCyan)

    #        s = BCyan + text + Color_Off
    #        if self.logMode >= self.LOG_MODE_INFO:
    #            print(s)
    #        self.write(s)

    def debug(self, text):
        self.log(None, text, self.LOG_MODE_DEBUG, Color_Yellow)

    #        s = Color_Yellow + text + Color_Off
    #        if self.logMode >= self.LOG_MODE_DEBUG:
    #            print(s)
    #        self.write(s)

    def print(self, text):
        self.log(None, text, self.LOG_MODE_PRINT, None)


#        if self.logMode >= self.LOG_MODE_PRINT:
#            print(text)
#        self.write(text)

def printError(s):
    print(BRed + "Error: " + s + Color_Off)


def printWarning(s):
    print(BGreen + "Warning: " + s + Color_Off)


def printInfo(s):
    print(BCyan + s + Color_Off)


def printDebug(s):
    print(Color_Yellow + s + Color_Off)


# RBL: it is moving to pmconventions as loadConfig
def loadPplSetup():
    pplSetup = {}
    userhome = getenv("HOME")
    #    f = open(userhome + "/bin/ppl-setup")
    f = open(userhome + "/.pmlib/ppl.cfg")
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
    if pmlog:
        pmlog.write('invoke: ' + cmd)
    r = cmd.split()
    try:
        a = subprocess.check_output(r)
        return a.decode('ascii')[:-1]
    except subprocess.CalledProcessError as e:
        return 'ERROR: exit code: %d, answer: %s' % (e.returncode, e.output)


def invokep(cmds):
    r = cmds[0].split()
    p = []
    p.append(subprocess.Popen(r, stdout=subprocess.PIPE))
    j = 0
    for cmd in cmds[1:]:
        r = cmd.split()
        p.append(subprocess.Popen(r, stdin=p[j].stdout, stdout=subprocess.PIPE))
        j += 1
    out = p[j].communicate()
    return out[0].decode('ascii')[:-1]


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
    t = Time(dateObs, format='isot')
    return t.jd


def getPair(s, delim=','):
    ss = s.split(delim)
    return [ss[0].strip(), ss[1].strip()]


# TODO: move to pmfits
def determineCoordsFromImage(imageFileName, pplSetup=None):
    configFolder = pplSetup["PMLIB"]
    if configFolder is None:
        configFolder = "$HOME/.pmlib"
    SOLVE_ARGS = "-O --config " + configFolder + "/astrometry.cfg --use-sextractor --sextractor-path sextractor -r -y -p"

    makedirs("temp", exist_ok=True)

    invoke("cp " + imageFileName + " temp/src.fits")

    astFolder = pplSetup["AST_BIN_FOLDER"]
    if astFolder is None:
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

    if not basePath:
        FOLDERS = ['./']
    elif isfile(basePath):
        path = ""
        j = basePath.rfind("/")
        if j != -1:
            path = basePath[:j + 1]
        FOLDERS = [path]
    else:
        FOLDERS = glob('*' + basePath.rstrip('/') + '*/')

    for folder in FOLDERS:
        if not folder.endswith('/'):
            folder += '/'

        f = open(folder + cmdName, "w+")
        f.write(f"#!/bin/bash\n#\n# ppl command: {cmdName}\n#\n\n")
        f.write(f"cd {getcwd()}\n")
        f.write(cmd + "\n")
        f.close()

        chmod(folder + cmdName, stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH)


def assureFolder(folder):
    if not exists(folder):
        makedirs(folder)
    return folder


def discoverFolders(baseFolder, folderName):
    folders = []
    if baseFolder is not None:
        folders = glob('*' + baseFolder + '*/' + folderName + '*')
    else:
        folders = glob(folderName + '*')
    return folders


# TODO: move tp pmfits
def getFitsHeaders(fitsFileName, headers):
    try:
        hdul = fits.open(fitsFileName)
    except Exception as e:
        print('Exception:', e)
        return None

    fitsHeaders = {}
    for h in headers:
        if h in hdul[0].header:
            fitsHeaders[h] = hdul[0].header[h]

    hdul.close()
    return fitsHeaders


# TODO: move tp pmfits
def getFitsHeader(fitsFileName, header):
    hdrs = getFitsHeaders(fitsFileName, [header])
    if hdrs and header in hdrs:
        return hdrs[header]
    else:
        return None


# TODO: move tp pmfits
def setFitsHeader(fitsFileName, headerName, headerValue, comment=None):
    if comment is not None:
        headers = {headerName: (headerValue, comment)}
    else:
        headers = {headerName: headerValue}
    setFitsHeaders(fitsFileName, headers)


# TODO: move tp pmfits
def setFitsHeaders(fitsFileName, headers):
    try:
        hdul = fits.open(fitsFileName, mode='update')
    except Exception as e:
        print(fitsFileName)
        print('Exception:', e)
        return

    h = hdul[0].header
    for key in headers.keys():
        h[key] = headers[key]
    hdul.flush()
    hdul.close()


# TODO: move tp pmfits
def subtractFitsBackground(fitsFileName):
    if exists(fitsFileName):
        hdul = fits.open(fitsFileName, mode='update')
        data = hdul[0].data
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        bkg = Background2D(data, (50, 50), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        # print("Bkg median: %7.4f, RMS median: %7.4f" % (bkg.background_median, bkg.background_rms_median))
        data -= bkg.background
        hdul.flush()
        hdul.close()


def findInFile(fileName, s):
    f = open(fileName)
    for line in f:
        if s in line:
            f.close()
            return line
    f.close()
    return None


def quad(a, b):
    return sqrt(a * a + b * b)


def guess(path):
    """
    try to identify some data from folder or file name: target object, observation date, observer, filter
    """

    # guess target
    r = path.strip('/').split('/')
    f = None
    for s in r:
        if '_' in s:
            f = s
    target_s = f.split('_', 1)[1]
    cats = glob('/home/kovi/.pmlib/cat/' + target_s + '*')
    if len(cats) > 0:
        cat = min(cats, key=len) if len(cats) > 1 else cats[0]
        target_s = cat.rsplit('/', 1)[1].split('.')[0]
    target = target_s.replace('_', ' ')
    return {'target': target}


def linfit(x, y, w=None):
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


class Plot:
    INV_X = 1
    INV_Y = 2

    def __init__(self, count, show, save):
        self.show = show
        self.save = save
        if show or save:
            fig = plt.figure(figsize=[6.4, 4.8 * count])
            #            fig.tight_layout(h_pad=1.7)
            self.plot_index = 100 * count + 11

    def add(self, xdata, ydata, coef, xlabel, ylabel, dotcolor, invaxis=None, title=None):
        if self.show or self.save:
            ax = plt.subplot(self.plot_index)
            if invaxis is not None:
                if invaxis & 1 == 1:
                    ax.invert_xaxis()
                if invaxis & 2 == 2:
                    ax.invert_yaxis()
            # ax[self.plot_index % 10].set_title('Plot title')
            xr = np.arange(min(xdata), max(xdata), 0.01)
            plt.plot(xdata, ydata, dotcolor + 'o', xr, coef[0] * xr + coef[1], 'k')
            if title:
                plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            self.plot_index += 1

    def showOrSave(self, fileName):
        if self.save or self.show:
            plt.subplots_adjust(hspace=0.3)
        if self.save:
            plt.savefig(fileName, format='png', dpi=72.0)
        if self.show:
            plt.show()
        elif self.save:
            plt.close()


def addTableComment(tbl, key, value):
    cmt = key + ': ' + value
    print(f'add table comment: {cmt}')
    if 'comments' not in tbl.meta:
        tbl.meta['comments'] = []
    tbl.meta['comments'].append(cmt)


def getTableComment(tbl, key):
    if 'comments' in tbl.meta:
        cmt = next(filter(lambda x: x.startswith(key + ':'), tbl.meta['comments']))
        if cmt:
            return cmt.partition(':')[2].strip()
    return None


if setup is None:
    setup = loadPplSetup()

# end pmbase.
