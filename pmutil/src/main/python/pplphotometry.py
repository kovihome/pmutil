#!/usr/bin/env python3
#
# PmUtils/pplphotometry
#
'''
Created on Feb 22, 2020

@author: kovi
'''
from getopt import getopt, GetoptError
from sys import argv
from datetime import datetime
from glob import glob
from os import getcwd, mkdir
from os.path import isdir, exists, basename

from pmbase import printError, printWarning, printInfo, saveCommand, loadPplSetup, invoke, Blue, Color_Off, BGreen, discoverFolders
from pmphot import Photometry



class Pipeline:

    opt = {}                 # command line options
    pplSetup = {}            # PPL setup from ppl-setup config file

    AST_CFG = ''             # astrometry.net config file
    SEX_CFG = ''             # sextractor config file
    SOLVE_ARGS = ''          # astrometry.net command line arguments

    def __init__(self, opt):
        self.opt = opt

    def photometry(self, seqFolder, photFolder, refcatFileName, color):
    
        # Names of the sequence/combined files:
        SEQLIST = glob(seqFolder + '/' + self.pplSetup['SEQ_FILE_PREFIX'] + '*-' + color + '.fits')
        if len(SEQLIST) == 0:
            SEQLIST = glob(seqFolder + '/Combined-' + color + '.fits')
            if len(SEQLIST) == 0:
              printWarning("No files for photometry in folder %s" % (seqFolder))
              return False

        for f in SEQLIST:

            print("photometry of %s" % (f))

            AST_FILE = photFolder + '/' + basename(f).replace('.fits', '.ast.fits')
            PMCAT_FILE = photFolder + '/' + basename(f).replace('.fits', '.cat')

            printInfo("Make astrometry for %s" % (f))
            invoke("%s/solve-field %s -D %s -N %s %s" % (self.pplSetup['AST_BIN_FOLDER'], self.SOLVE_ARGS, photFolder, AST_FILE, f))

            printInfo("Make photometry for %s" % (f))
            invoke("sextractor %s -c %s -CATALOG_NAME %s -CATALOG_TYPE ASCII_HEAD" % (AST_FILE, self.SEX_CFG, PMCAT_FILE))

            PMCAT_FILE_FLT = PMCAT_FILE + ".cat"
            printInfo("Filtering result catalog to %s" % (PMCAT_FILE_FLT))
            # TODO: call directly
            invoke("pmfilter -c %s -r %s -o %s %s" % (color, refcatFileName, PMCAT_FILE_FLT, PMCAT_FILE))

            self.calculateMags(color, PMCAT_FILE_FLT)

            #PM_FILE = PMCAT_FILE_FLT + '.pm'
            #printInfo("Calculate real magnitudes to %s" % (PM_FILE))
            #invoke("pmphot -c %s -o %s %s" % (color, PM_FILE, PMCAT_FILE_FLT))    

        return True

    def calculateMags(self, photFolder, colors):

        # Names of the sequence/combined files:
        color = colors[0]
        PHOTLIST = glob(photFolder + '/' + self.pplSetup['SEQ_FILE_PREFIX'] + '*-' + color + '.cat.cat')
        if len(PHOTLIST) == 0:
            PHOTLIST = glob(photFolder + '/Combined-' + color + '.cat.cat')
            if len(PHOTLIST) == 0:
              printWarning("No files for calculate magnitudes in folder %s" % (photFolder))
              return False

        for f in PHOTLIST:
            inputFiles = []
            for c in colors:
                inputFiles.append(f.replace('-' + color, '-' + c))

            #PM_FILE = PMCAT_FILE_FLT + '.pm'
            printInfo("Calculate real magnitudes to %s" % (f.replace('-' + color, '-*') + '.pm'))
            opt = {
                'out' : None,
                'comp': None,
                'color' : colors,
                'std': False,
                'makestd': False,
                'files': inputFiles,
            }
            phot = Photometry(opt)
            phot.process()


    def process_photometry(self, seqFolder, photFolder, title):

        printInfo("%s: Make photometry on sequence file(s)." % (title))

        # Create the photometry dir, if not exists
        if not exists(photFolder):
            mkdir(photFolder)

        # Photometry reference file:
        sf = seqFolder
        basedir = sf.replace(self.pplSetup['SEQ_FOLDER_NAME'], '')
        PMREFS = glob(basedir + '/*.cat')
        if len(PMREFS) > 0 and exists(PMREFS[0]):
            PMREF = PMREFS[0]
        else:
            printError("No reference catalog file (.cat) in folder %s" % (basedir))
            printInfo("Use ppl-refcat command to create reference catalog for the object.")
            return False    

        for color in self.opt['color']:
            self.photometry(seqFolder, photFolder, PMREF, color)

        self.calculateMags(photFolder, self.opt['color'])

        # TODO: cleanup - delete light FITS files

        return True


    def execute(self):

        ##########################
        # step 0. setup photometry
        ##########################
        self.pplSetup = loadPplSetup()

        seqFolders = discoverFolders(self.opt['baseFolder'], self.pplSetup['SEQ_FOLDER_NAME'])
        if len(seqFolders) == 0:
            printError('No sequence folder found in base folder %s' % (self.opt['baseFolder']))
            exit(1)
        print("Sequence folders discovered:", seqFolders)
        
        ### Common arguments (saturation level, image section & trimming, etc.):
        self.AST_CFG = self.pplSetup['CONFIG_FOLDER'] + "/astrometry.cfg"
        self.SEX_CFG = self.pplSetup['CONFIG_FOLDER'] + "/sex.cfg"
        #SOLVE_ARGS="-O --config $AST_CFG --use-sextractor --sextractor-path sextractor -i ${PHOTDIR}/scamp.cat -n ${PHOTDIR}/scamp.cfg -j ${PHOTDIR}/scamp.ref -r -y -p"
        self.SOLVE_ARGS = "-O --config %s --use-sextractor --sextractor-path sextractor -r -y -p" % (self.AST_CFG)


        ####################################
        # step 1. do photometry for all file
        ####################################

        PMERROR = False
        for sf in seqFolders:

            pf = sf.replace(self.pplSetup['SEQ_FOLDER_NAME'], self.pplSetup['PHOT_FOLDER_NAME'])
    
            success = self.process_photometry(sf, pf, "PHOTOMETRY")
            if not success:
                PMERROR = True
                break	
  

        ########################################
        # step 2. create report from all results
        ########################################

        if not PMERROR:
            if self.opt['baseFolder'] != None:
                PM_FILES = glob('*' + self.opt['baseFolder'] + '*/' + self.pplSetup['PHOT_FOLDER_NAME'] + '/*.pm')
                BASE_FOLDERS = glob('*' + self.opt['baseFolder'] + '*')
                REPORT_FOLDER=BASE_FOLDERS[0]
            else:
                PM_FILES = glob(self.pplSetup['PHOT_FOLDER_NAME'] + '/*.pm')
                REPORT_FOLDER = getcwd()

            printInfo("Create report into %s folder." % (REPORT_FOLDER))
            # TODO: call directly
            invoke("pmresult -o %s -n %s %s" % (REPORT_FOLDER, self.opt['nameCode'], ' '.join(PM_FILES)))

        print()
        print(Blue + "Photometry is ready." + Color_Off)


class MainApp:
    
    opt = {
        'color' : ['Gi'],    # photometry band, mandatory
        'nameCode' : None,   # observer code for the AAVSO report, mandatory
        'overwrite': False,  # force to overwrite existing results, optional
        'files': None,       
        'baseFolder': None,  # base folder, optional
        }

    availableBands = ['gi', 'g', 'bi', 'b', 'ri', 'r', 'all']


    def __init__(self, argv):
        self.argv = argv
        pass

    def printTitle(self):
        print()
        print(BGreen + "pplphotometry, version 1.1.0 " + Color_Off)
        print(Blue + "Make photometry on calibrated FITS images." + Color_Off)
        print()


    def usage(self):
        print("Usage: pplphotometry [OPTIONS]... [BASE_FOLDER]")
        print("Make photometry on calibrated FITS images.")
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -c,  --color arg        set filter(s), arg is the color code, default color is 'Gi', for available color codes see below")
        print("  -n,  --name nameCode    set observer code for the AAVSO report")
        print("  -h,  --help             print this page")
        print("  -w,  --overwrite        force to overwrite existing results")
        print()
        print("Available filter color codes are:")
        print("  Gi | G | gi | g         green channel")
        print("  Bi | B | bi | b         blue channel")
        print("  Ri | R | ri | r         red channel")
        print("  all | ALL | All         all channels, results 3 separate frame")
        print()


    def processCommands(self):
        try:
            optlist, args = getopt (self.argv[1:], "c:n:wh", ['--color', '--name', '--help', '--overwrite'])
        except GetoptError:
            printError('Invalid command line options.')
            return

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-c':
                color = a.lower()
                if not color in self.availableBands:
                    printError('Invalid color: %s, use on of these: Gi, g, Bi, b, Ri, r, all' % (a))
                    exit(1)
                if color == 'all':
                    self.opt['color'] = ['Ri', 'Gi', 'Bi']
                else:
                    self.opt['color'] = [a]
            elif o == '-n':
                self.opt['nameCode'] = a.upper()
            elif o == '-w':
                self.opt['makestd'] = True
            elif o == '-h':
                self.usage()
                exit(0)

        if len(args) > 0:
            self.opt['baseFolder'] = args[0]
            if args[0].endswith('/'):
                self.opt['baseFolder'] = args[0][:-1]

        if self.opt['nameCode'] == None:
            printWarning('No observer code was given. Use \'XYZ\' instead.')
            self.opt['nameCode'] = 'XYZ'
    

    def run(self):
        self.printTitle()
        self.processCommands()
        saveCommand(self.opt['baseFolder'], self.argv, 'photometry')

        start = datetime.now()

        ppl = Pipeline(self.opt)
        ppl.execute()

        exectime = (datetime.now() - start).total_seconds()
        print("%sexecution time was %d seconds.%s" % (Blue, exectime, Color_Off))



if __name__ == '__main__':

    app = MainApp(argv)
    app.run()

# end __main__
