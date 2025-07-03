#!/usr/bin/env python3
#
# PmUtils/pmconventions (renamed from pmdisco)
#
"""
Created on Mar 2, 2020

@author: kovi
"""
from os import getenv
from os.path import isdir, exists
from glob import glob

from astropy.table import Table

PMUTIL_VERSION = "1.2.0"
PMUTIL_VERSION_SHORT = "1.2"

RAW_FILE_EXTENSIONS = ["cr2", "cr3"]

FITS_FILE_EXTENSIONS = ["fits", "fit", "fts"]

class Discovery:
    BIAS_FOLDER = None
    DARK_FOLDER = None
    FLAT_FOLDER = None
    FLAT_BIAS_FOLDER = None
    FLAT_DARK_FOLDER = None
    LIGHT_FOLDERS = None

    def __init__(self, opt, pplSetup):
        self.flatOnly = opt['flatOnly']
        self.nonFlatFolderOnly = opt['useMasterFlat']
        self.folderPattern = opt['baseFolder']
        self.calibFolder = opt['calibFolder']
        self.ppl = pplSetup

    @staticmethod
    def findFolders(baseFolder, folderName):
        return glob('*' + baseFolder + '*/' + folderName) if baseFolder is not None else glob(folderName)

    def discoverFolder(self, baseFolder, folderName, title):
        FOLDERS = self.findFolders(baseFolder, folderName)

        if len(FOLDERS) == 0:
            print(f"Error: no {title} folder found; add one, and rerun this script.")
            exit(1)

        if len(FOLDERS) != 1:
            print(f"Error: more than one {title} folder found; remove on of them, and rerun this script.")
            print("  --> " + ' '.join(FOLDERS))
            exit(1)

        if not isdir(FOLDERS[0]):
            print(f"Error: {title} folder {FOLDERS[0]} is not exist ot not a directory.")
            exit(1)

        print(f"{title} folder discovered: {FOLDERS[0]}")
        return FOLDERS[0]

    # discover Bias folder
    def discoverBiasFolder(self, baseFolder):
        self.BIAS_FOLDER = self.discoverFolder(baseFolder if not self.calibFolder else self.calibFolder,
                                               self.ppl['BIAS_FOLDER_NAME'], 'Bias')

    # discover Dark folder
    def discoverDarkFolder(self, baseFolder):
        self.DARK_FOLDER = self.discoverFolder(baseFolder if not self.calibFolder else self.calibFolder,
                                               self.ppl['DARK_FOLDER_NAME'], 'Dark')

    def discover(self):
        if not self.flatOnly:
            self.discoverBiasFolder(self.folderPattern)
            self.discoverDarkFolder(self.folderPattern)

        if not self.nonFlatFolderOnly:

            # discover flat bias folder
            FLAT_BIAS_FOLDERS = self.findFolders(self.folderPattern if not self.calibFolder else self.calibFolder,
                                                 self.ppl['FLAT_BIAS_FOLDER_NAME'])
            if len(FLAT_BIAS_FOLDERS) > 1:
                print(f"Error: more than one {'Flat Bias'} folder found; remove on of them, and rerun this script.")
                print("  --> " + ' '.join(FLAT_BIAS_FOLDERS))
                exit(1)

            if len(FLAT_BIAS_FOLDERS) == 0 or not isdir(FLAT_BIAS_FOLDERS[0]):
                if self.flatOnly:
                    self.discoverBiasFolder(self.folderPattern)
                self.FLAT_BIAS_FOLDER = self.BIAS_FOLDER
            else:
                self.FLAT_BIAS_FOLDER = FLAT_BIAS_FOLDERS[0]
            print(f"Flat Bias folder discovered: {self.FLAT_BIAS_FOLDER}")

            # discover flat dark folder
            FLAT_DARK_FOLDERS = self.findFolders(self.folderPattern if not self.calibFolder else self.calibFolder,
                                                 self.ppl['FLAT_DARK_FOLDER_NAME'])
            if len(FLAT_DARK_FOLDERS) > 1:
                print(f"Error: more than one {'Flat Dark'} folder found; remove on of them, and rerun this script.")
                print("  --> " + ' '.join(FLAT_DARK_FOLDERS))
                exit(1)

            if len(FLAT_DARK_FOLDERS) == 0 or not isdir(FLAT_DARK_FOLDERS[0]):
                if self.flatOnly:
                    self.discoverDarkFolder(self.folderPattern)
                self.FLAT_DARK_FOLDER = self.DARK_FOLDER
            else:
                self.FLAT_DARK_FOLDER = FLAT_DARK_FOLDERS[0]
            print(f"Flat Dark folder discovered: {self.FLAT_DARK_FOLDER}")

            # discover flat folder
            self.FLAT_FOLDER = self.discoverFolder(self.folderPattern if not self.calibFolder else self.calibFolder,
                                                   self.ppl['FLAT_FOLDER_NAME'], 'Flat')

        if not self.flatOnly:
            # discovery Light folders
            self.LIGHT_FOLDERS = self.findFolders(self.folderPattern, self.ppl['LIGHT_FOLDER_NAME'])

            print("Light folders discovered:" + ' '.join(self.LIGHT_FOLDERS))


def findBaseFolders(baseFolderPattern):
    if baseFolderPattern == '' or baseFolderPattern == '.' or baseFolderPattern == './':
        return ['.']
    baseFolders = glob('*' + baseFolderPattern + '*')
    return baseFolders


def loadConfig():
    pplSetup = {}
    userhome = getenv("HOME")
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


# Refcat functions

# TODO: pmmerge, pplphotometry, pmresult, pmphot
def loadRefcat(baseFolder):
    refcatFileName = baseFolder.rstrip('/') + '/ref.cat'
    return Table.read(refcatFileName, format='ascii') if exists(refcatFileName) else None
    
