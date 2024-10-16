#!/usr/bin/env python3
#
# PmUtils/pmclean
#
"""
Created on Mar 22, 2020

@author: kovi
"""

from datetime import datetime
from getopt import getopt, GetoptError
from glob import glob
from os import remove, removedirs, symlink
from os.path import exists, abspath, basename
from sys import argv
from zipfile import ZipFile

import pmbase as pm
from pmconventions import findBaseFolders, RAW_FILE_EXTENSIONS, PMUTIL_VERSION


class Cleaner:
    opt = {}  # command line options

    TEMPDIR = "temp"

    def __init__(self, opt):
        self.opt = opt

    def removeFiles(self, filePattern, exceptions=None):
        fs = glob(filePattern)
        if exceptions is not None:
            fs = [f for f in fs if not f.endswith(exceptions)]
        for f in fs:
            remove(f)

    def cleanFolder(self, baseFolder, subfolder, filePattern='*'):
        folderPattern = baseFolder + '/' + subfolder
        if not exists(folderPattern):
            return
        print(f"Cleaning folder: {folderPattern}")

        filePatterns = [filePattern] if type(filePattern) is str else filePattern
        for fp in filePatterns:
            self.removeFiles(folderPattern + '/' + fp)

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

        self.cleanFolder(baseFolder, pm.setup['BIAS_FOLDER_NAME'], pm.setup['BIAS_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(baseFolder, pm.setup['DARK_FOLDER_NAME'], pm.setup['DARK_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(baseFolder, pm.setup['FLAT_FOLDER_NAME'], pm.setup['FLAT_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(baseFolder, pm.setup['FLAT_BIAS_FOLDER_NAME'], pm.setup['BIAS_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(baseFolder, pm.setup['FLAT_DARK_FOLDER_NAME'], pm.setup['DARK_FILE_PREFIX'] + '*.fits')
        self.cleanFolder(baseFolder, pm.setup['LIGHT_FOLDER_NAME'], pm.setup['LIGHT_FILE_PREFIX'] + '*.fits')

        # clean all files in Calibrated folder
        self.removeFolder(baseFolder, pm.setup['CALIB_FOLDER_NAME'])

        # clean photometry folder
        for file_pattern in ['*.axy', '*.corr', '*.idmatch', '*.match', '*.rdls', '*.solved', '*.wcs', '*.xyls']:
            self.cleanFolder(baseFolder, pm.setup['PHOT_FOLDER_NAME'], file_pattern)

    def removeAll(self, baseFolder):
        print(f"Remove all folders in {baseFolder}")

        # remove all subfolders
        self.removeFolder(baseFolder, self.TEMPDIR)
        self.removeFolder(baseFolder, pm.setup['BIAS_FOLDER_NAME'])
        self.removeFolder(baseFolder, pm.setup['DARK_FOLDER_NAME'])
        self.removeFolder(baseFolder, pm.setup['FLAT_FOLDER_NAME'])
        self.removeFolder(baseFolder, pm.setup['FLAT_BIAS_FOLDER_NAME'])
        self.removeFolder(baseFolder, pm.setup['FLAT_DARK_FOLDER_NAME'])
        self.removeFolder(baseFolder, pm.setup['LIGHT_FOLDER_NAME'])
        self.removeFolder(baseFolder, pm.setup['CALIB_FOLDER_NAME'])
        self.removeFolder(baseFolder, pm.setup['SEQ_FOLDER_NAME'])
        self.removeFolder(baseFolder, pm.setup['PHOT_FOLDER_NAME'])

        # clean up base folder, except archive files
        self.removeFiles(baseFolder + "/*", exceptions=".zip")

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
            archiveFolder = pm.setup['ARCHLIB'] if 'ARCHLIB' in pm.setup and pm.setup['ARCHLIB'] != '' else baseFolder
        pm.assureFolder(archiveFolder)
        archiveFileName = archiveFolder + '/' + basename(abspath(baseFolder)) + '.zip'
        print(f"Archive file name: {archiveFileName}")

        archiveFile = ZipFile(archiveFileName, "w")

        for ext in RAW_FILE_EXTENSIONS:
            self.archiveFolder(archiveFile, baseFolder, pm.setup['BIAS_FOLDER_NAME'],
                               pm.setup['BIAS_FILE_PREFIX'] + f'*.{ext}')
            self.archiveFolder(archiveFile, baseFolder, pm.setup['DARK_FOLDER_NAME'],
                               pm.setup['DARK_FILE_PREFIX'] + f'*.{ext}')
            self.archiveFolder(archiveFile, baseFolder, pm.setup['FLAT_FOLDER_NAME'],
                               pm.setup['FLAT_FILE_PREFIX'] + f'*.{ext}')
            self.archiveFolder(archiveFile, baseFolder, pm.setup['FLAT_BIAS_FOLDER_NAME'],
                               pm.setup['BIAS_FILE_PREFIX'] + f'*.{ext}')
            self.archiveFolder(archiveFile, baseFolder, pm.setup['FLAT_DARK_FOLDER_NAME'],
                               pm.setup['DARK_FILE_PREFIX'] + f'*.{ext}')
            self.archiveFolder(archiveFile, baseFolder, pm.setup['LIGHT_FOLDER_NAME'],
                               pm.setup['LIGHT_FILE_PREFIX'] + f'*.{ext}')

        self.archiveFolder(archiveFile, baseFolder, pm.setup['BIAS_FOLDER_NAME'], 'master-bias-*.fits')
        self.archiveFolder(archiveFile, baseFolder, pm.setup['DARK_FOLDER_NAME'], 'master-dark.*.fits')
        self.archiveFolder(archiveFile, baseFolder, pm.setup['FLAT_FOLDER_NAME'], 'master-flat-*.fits')
        self.archiveFolder(archiveFile, baseFolder, pm.setup['FLAT_BIAS_FOLDER_NAME'], 'master-bias-*.fits')
        self.archiveFolder(archiveFile, baseFolder, pm.setup['FLAT_DARK_FOLDER_NAME'], 'master-dark.*.fits')

        self.archiveFolder(archiveFile, baseFolder, pm.setup['SEQ_FOLDER_NAME'], 'Seq-*.fits')
        self.archiveFolder(archiveFile, baseFolder, pm.setup['SEQ_FOLDER_NAME'], 'Combined-*.fits')

        self.archiveFolder(archiveFile, baseFolder, pm.setup['PHOT_FOLDER_NAME'], '*.ast.fits')
        self.archiveFolder(archiveFile, baseFolder, pm.setup['PHOT_FOLDER_NAME'], '*.pm')
        self.archiveFolder(archiveFile, baseFolder, pm.setup['PHOT_FOLDER_NAME'], '*[BGR]i.cat')

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
        baseFolders = findBaseFolders(self.opt['baseFolder'])
        for bf in baseFolders:
            if self.opt['archive']:
                pm.printInfo(f"Archive base folder: {bf}")
                archiveFileName = self.archive(bf)
                self.removeAll(bf)
                if basename(archiveFileName) != basename(bf):
                    self.createLink(bf, archiveFileName)
            else:
                pm.printInfo(f"Clean base folder: {bf}")
                self.clean(bf)


class MainApp:
    appName = 'ppl-clean'
    appVersion = PMUTIL_VERSION
    appDescription = 'Clean/archive all generated FITS and other files.'

    opt = {
        'lightsAlso': False,  # remove light FITS too - ignored
        'archive': False,  # save result files into archive folder
        'archiveFolder': None,  # archive folder
        'debug': False,  # debug print
        'baseFolder': '.',  # base folder
    }

    def __init__(self, argv):
        self.argv = argv

    def printVersion(self):
        print(pm.BGreen + self.appName + ", version " + self.appVersion + pm.Color_Off)

    def printTitle(self):
        print()
        self.printVersion()
        print(pm.Blue + self.appDescription + pm.Color_Off)
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
            optlist, args = getopt(argv[1:], "lah", ['lights', 'archive', 'archive-folder=', 'help'])
        except GetoptError:
            print('Invalid command line options')
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
        if self.opt['archiveFolder'] is not None and self.opt['archiveFolder'].endswith('/'):
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
        print("%sexecution time was %d seconds.%s" % (pm.Blue, exectime, pm.Color_Off))


if __name__ == '__main__':
    app = MainApp(argv)
    app.run()

# end main.
