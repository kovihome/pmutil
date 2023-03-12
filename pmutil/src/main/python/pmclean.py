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
from os import remove, removedirs, symlink
from os.path import isdir, exists, abspath, basename
from glob import glob
from datetime import datetime
from zipfile import ZipFile

from pmbase import printError, printInfo, Blue, Color_Off, BGreen, assureFolder
from pmconventions import findBaseFolders, loadConfig


class Cleaner:

    opt = {}  # command line options
    ppl = {}  # PPL setup from ppl-setup config file

    TEMPDIR="temp"

    def __init__(self, opt):
        self.opt = opt

    def removeFiles(self, filePattern, exceptions = None):
        fs = glob(filePattern)
        if exceptions != None:
            fs = [f for f in fs if not f.endswith(exceptions)]
        for f in fs:
            remove(f)

    def cleanFolder(self, baseFolder, subfolder, filePattern = '*'):
        folderPattern = baseFolder + '/' + subfolder
        if not exists(folderPattern):
            return
        print(f"Cleaning folder: {folderPattern}")

        filePatterns = [filePattern] if type(filePattern) is str else filePattern
        for fp in filePatterns:
            self.removeFiles(folderPattern + '/'+ fp)

    def removeFolder(self, baseFolder, subfolder):
        folderPattern = baseFolder + '/' + subfolder
        if not exists(folderPattern):
            return
        print(f"Removing folder: {folderPattern}")

        self.removeFiles(folderPattern + "/*")

        removedirs(folderPattern)
        
    def clean(self, baseFolder):

        # clean and remove Temp folder
        self.removeFolder(baseFolder, self.TEMPDIR)

        self.cleanFolder(baseFolder, self.ppl['BIAS_FOLDER_NAME'], self.ppl['BIAS_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(baseFolder, self.ppl['DARK_FOLDER_NAME'], self.ppl['DARK_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(baseFolder, self.ppl['FLAT_FOLDER_NAME'], self.ppl['FLAT_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(baseFolder, self.ppl['FLAT_BIAS_FOLDER_NAME'], self.ppl['BIAS_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(baseFolder, self.ppl['FLAT_DARK_FOLDER_NAME'], self.ppl['DARK_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(baseFolder, self.ppl['LIGHT_FOLDER_NAME'], self.ppl['LIGHT_FILE_PREFIX'] + '*.fits')

        # clean all files in Calibrated folder
        self.removeFolder(baseFolder, self.ppl['CALIB_FOLDER_NAME'])

        # clean photometry folder
        self.cleanFolder(baseFolder, self.ppl['PHOT_FOLDER_NAME'], ['*.axy', '*.corr', '*.idmatch', '*.match', '*.rdls', '*.solved', '*.wcs', '*.xyls'])


    def removeAll(self, baseFolder):
        print(f"Remove all folders in {baseFolder}")

        # remove all subfolders
        self.removeFolder(baseFolder, self.TEMPDIR)
        self.removeFolder(baseFolder, self.ppl['BIAS_FOLDER_NAME'])
        self.removeFolder(baseFolder, self.ppl['DARK_FOLDER_NAME'])
        self.removeFolder(baseFolder, self.ppl['FLAT_FOLDER_NAME'])
        self.removeFolder(baseFolder, self.ppl['FLAT_BIAS_FOLDER_NAME'])
        self.removeFolder(baseFolder, self.ppl['FLAT_DARK_FOLDER_NAME'])
        self.removeFolder(baseFolder, self.ppl['LIGHT_FOLDER_NAME'])
        self.removeFolder(baseFolder, self.ppl['CALIB_FOLDER_NAME'])
        self.removeFolder(baseFolder, self.ppl['SEQ_FOLDER_NAME'])
        self.removeFolder(baseFolder, self.ppl['PHOT_FOLDER_NAME'])

        # clean up base folder, except archive files
        self.removeFiles(baseFolder + "/*", exceptions = ".zip")

    def archiveFolder(self, archiveFile, baseFolder, subfolder, filePattern):
        files = f"{baseFolder}/{subfolder}/{filePattern}" if subfolder else f"{baseFolder}/{filePattern}"
        fs = glob(files)
        if len(fs) == 0:
            return
        print(f"Adding to archive: {files}")

        for f in fs:
            archiveFile.write(f)


    def archive(self, baseFolder):

        # calculate archive file name
        archiveFolder = self.opt['archiveFolder']
        if not archiveFolder:
            archiveFolder = self.ppl['ARCHLIB'] if 'ARCHLIB' in self.ppl and self.ppl['ARCHLIB'] != '' else baseFolder
        assureFolder(archiveFolder)
        archiveFileName = archiveFolder + '/' + basename(abspath(baseFolder)) + '.zip'
        print(f"Archive file name: {archiveFileName}")

        archiveFile = ZipFile(archiveFileName, "w")

        self.archiveFolder(archiveFile, baseFolder, self.ppl['BIAS_FOLDER_NAME'], self.ppl['BIAS_FILE_PREFIX'] + '*.cr2')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['BIAS_FOLDER_NAME'], 'master-bias-*.fits')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['DARK_FOLDER_NAME'], self.ppl['DARK_FILE_PREFIX'] + '*.cr2')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['DARK_FOLDER_NAME'], 'master-dark.*.fits')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['FLAT_FOLDER_NAME'], self.ppl['FLAT_FILE_PREFIX'] + '*.cr2')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['FLAT_FOLDER_NAME'], 'master-flat-*.fits')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['FLAT_BIAS_FOLDER_NAME'], self.ppl['BIAS_FILE_PREFIX'] + '*.cr2')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['FLAT_BIAS_FOLDER_NAME'], 'master-bias-*.fits')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['FLAT_DARK_FOLDER_NAME'], self.ppl['DARK_FILE_PREFIX'] + '*.cr2')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['FLAT_DARK_FOLDER_NAME'], 'master-dark.*.fits')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['LIGHT_FOLDER_NAME'], self.ppl['LIGHT_FILE_PREFIX'] + '*.cr2')

        self.archiveFolder(archiveFile, baseFolder, self.ppl['SEQ_FOLDER_NAME'], 'Seq-*.fits')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['SEQ_FOLDER_NAME'], 'Combined-*.fits')

        self.archiveFolder(archiveFile, baseFolder, self.ppl['PHOT_FOLDER_NAME'], '*.ast.fits')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['PHOT_FOLDER_NAME'], '*.pm')
        self.archiveFolder(archiveFile, baseFolder, self.ppl['PHOT_FOLDER_NAME'], '*[BGR]i.cat')

        self.archiveFolder(archiveFile, baseFolder, None, 'calibration')
        self.archiveFolder(archiveFile, baseFolder, None, 'photometry')
        self.archiveFolder(archiveFile, baseFolder, None, 'refcat')
        self.archiveFolder(archiveFile, baseFolder, None, 'ref.cat')
        self.archiveFolder(archiveFile, baseFolder, None, '*.aavso')
        self.archiveFolder(archiveFile, baseFolder, None, '*.jpg')

        archiveFile.close()

        return archiveFileName

    def createLink(self, baseFolder, archiveFile):
        if exists(baseFolder + "/archive.zip"):
            return
        symlink(archiveFile, baseFolder + "/archive.zip")

    def execute(self):
        self.ppl = loadConfig()
        baseFolders = findBaseFolders(self.opt['baseFolder'])
        for bf in baseFolders:
            if self.opt['archive']:
                printInfo(f"Archive base folder: {bf}")
                archiveFileName = self.archive(bf)
                self.removeAll(bf)
                if basename(archiveFileName) != basename(bf):
                    self.createLink (bf, archiveFileName)
            else:
                printInfo(f"Clean base folder: {bf}")
                self.clean(bf)


class MainApp:

    appName = 'ppl-clean'
    appVersion = '1.2.0'
    appDescription = 'Clean/archive all generated FITS and other files.'

    opt = {
        'lightsAlso'    : False,  # remove light FITS too - ignored
        'archive'       : False,  # save result files into archive folder
        'archiveFolder' : None,   # archive folder
        'debug'         : False,  # debug print
        'baseFolder'    : '.',    # base folder
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
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -l,  --lights                 remove FITS files in Light folder too")
        print("  -a,  --archive                save result files into the archive folder")
        print("       --archive-folder folder  set archive folder instead of the configured one")
        # print("  -w,  --overwrite              overwrite archive file")
        # print("       --debug                  print debug info too")
        print("  -h,  --help                   print this page")
        print()

    def processCommands(self):

        try:
            optlist, args = getopt (argv[1:], "lah", ['lights', 'archive', 'archive-folder=', 'help'])
        except GetoptError:
            print ('Invalid command line options')
            exit(1)

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-l' or o == '--lights':
                self.opt['lightsAlso'] = True
            elif o == '-a' or o == '--archive':
                self.opt['archive'] = True
            elif o == '--archive-folder':
                self.opt['archiveFolder'] = a
            elif o == '-h' or o == '--help':
                self.usage()
                exit(0)

        if len(args) > 0:
            self.opt['baseFolder'] = args[0]

        if self.opt['baseFolder'].endswith('/'):
            self.opt['baseFolder'] = self.opt['baseFolder'][:-1]
        if self.opt['archiveFolder'] and self.opt['archiveFolder'].endswith('/'):
            self.opt['archiveFolder'] = self.opt['archiveFolder'][:-1]

        if self.opt['baseFolder']:
            print(f"Base folder: {self.opt['baseFolder']}")
        if self.opt['archive']:
            print(f"Archive valuable files")
        if self.opt['archiveFolder']:
            print(f"Archive folder: {self.opt['archiveFolder']}")

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
