#!/usr/bin/env python3
#
# PmUtils/pmclean
#
'''
Created on Mar 22, 2020

@author: kovi
'''

from sys import argv
from getopt import getopt, GetoptError
from os import remove
from os.path import isdir, exists
from glob import glob
from datetime import datetime

from pmbase import printError, printInfo, loadPplSetup, Blue, Color_Off, BGreen, invoke

class Cleaner:

    opt = {}  # command line options
    ppl = {}  # PPL setup from ppl-setup config file

    TEMPDIR="temp"

    def __init__(self, opt):
        self.opt = opt

    def cleanFolder(self, subfolder, filePattern = '*'):
        folderPattern = '*' + self.opt['baseFolder'] + '*/' + subfolder
        removePattern = folderPattern + '/' + filePattern
        fs = glob(removePattern)
        if len(fs) > 0:
            printInfo("Cleaning folders: %s" % (folderPattern))
            for f in fs:
                invoke("rm -f %s" % (f))
        else:
            print("No %s files in folders %s" % (filePattern, folderPattern))

    def removeFolder(self, subfolder):
        folderPattern = '*' + self.opt['baseFolder'] + '*/' + subfolder
        fs = glob(folderPattern)
        printInfo("Removing folders: %s" % (folderPattern))
        for f in fs:
            invoke("rmdir %s" % (f))
        
    def cleanCalibrationFolders(self):

        # clean and remove Temp folder
        self.cleanFolder(self.TEMPDIR)
        self.removeFolder(self.TEMPDIR)

        self.cleanFolder(self.ppl['BIAS_FOLDER_NAME'], self.ppl['BIAS_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(self.ppl['DARK_FOLDER_NAME'], self.ppl['DARK_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(self.ppl['FLAT_FOLDER_NAME'], self.ppl['FLAT_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(self.ppl['FLAT_BIAS_FOLDER_NAME'], self.ppl['BIAS_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(self.ppl['FLAT_DARK_FOLDER_NAME'], self.ppl['DARK_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(self.ppl['LIGHT_FOLDER_NAME'], self.ppl['LIGHT_FILE_PREFIX'] + '*.fits')

        # clean all files in Calibrated folder
        self.cleanFolder(self.ppl['CALIB_FOLDER_NAME'])
        self.removeFolder(self.ppl['CALIB_FOLDER_NAME'])

    def cleanPhotometryFolders(self):
        print(glob('*' + self.opt['baseFolder'] + '*/' + self.ppl['PHOT_FOLDER_NAME']))		# 
        pass

    def execute(self):
        self.ppl = loadPplSetup()

        self.cleanCalibrationFolders()

        self.cleanPhotometryFolders()


class MainApp:

    appName = 'ppl-clean'
    appVersion = '1.1.0'
    appDescription = 'Clean all generated FITS and other files.'

    opt = {
        'lightsAlso' : False,  # remove light FITS too - ignored
        'baseFolder' : '.',    # base folder
        }

    def __init__(self, argv):
        self.argv = argv

    def printVersion(self):
        print(BGreen + self.appName + ", version " + self.appVersion + Color_Off)

    def printTitle(self):
        print()
        self.printVersion()
        print(Blue + self.appDescription + Color_Off)
        print()

    def usage(self):
        print("Usage: " + self.appName + " [OPTIONS]... FOLDER_NAME")
        print(self.appDescription)
        print()
        #print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -l,  --lights     remove FITS files in Light folder too")
        print("  -h,  --help       print this page")
        print()

    def processCommands(self):

        try:
            optlist, args = getopt (argv[1:], "lh", ['lights', 'help'])
        except GetoptError:
            print ('Invalid command line options')
            exit(1)

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-l' or o == '--lights':
                self.opt['lightsAlso'] = True
            elif o == '-h' or o == '--help':
                self.usage()
                exit(0)

        if len(args) > 0:
#            if not isdir(args[0]):
#                printError('%s is not a folder.' % (args[0]))
#                exit(1)
            self.opt['baseFolder'] = args[0]

#        if not self.opt['baseFolder'].endswith('/'):
#            self.opt['baseFolder'] += '/'

    def run(self):
        self.printTitle()
        self.processCommands()

        start = datetime.now()

        cleaner = Cleaner(self.opt)
        cleaner.execute()

        exectime = (datetime.now() - start).total_seconds()
        print("%sexecution time was %d seconds.%s" % (Blue, exectime, Color_Off))


if __name__ == '__main__':

    app = MainApp(argv)
    app.run()

# end main.
