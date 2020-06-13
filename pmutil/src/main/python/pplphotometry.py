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
from os.path import exists, basename
from astropy.table import Table

from pmbase import printError, printWarning, printInfo, saveCommand, loadPplSetup, invoke, Blue, Color_Off, BGreen, discoverFolders
from pmphot import Photometry
from pmresult import ReportProcessor
from pmfilter import CatalogMatcher


class Pipeline:

    opt = {}  # command line options
    pplSetup = {}  # PPL setup from ppl-setup config file

    AST_CFG = ''  # astrometry.net config file
    SEX_CFG = ''  # sextractor config file
    SOLVE_ARGS = ''  # astrometry.net command line arguments

    def __init__(self, opt):
        self.opt = opt

    def inspectRefCat(self, baseFolder):

        colors = ''
        if 'Bi' in self.opt['color']:
            colors += 'b'
        if 'Gi' in self.opt['color'] or 'Vi' in self.opt['color']:
            colors += 'v'
        if 'Ri' in self.opt['color']:
            colors += 'r'

        refcatFileName = baseFolder.rstrip('/') + '/ref.cat'
        print('Inspect refcat %s' % (refcatFileName))
        if not exists(refcatFileName):
            printWarning("Reference catalog %s is missing." % (refcatFileName))
            return False
        rc = Table.read(refcatFileName, format='ascii')
        #print(rc)
        varCount = 0
        compCount = 0
        compColorCount = { 'bvr': 0, 'bv': 0, 'vr': 0, 'b': 0, 'v': 0, 'r': 0 }
        for r in rc:
            if r['ROLE'] == 'C':
                compCount += 1
                has_b = r['MAG_B'] != -1.0 and r['MAG_B'] != '-'
                has_v = r['MAG_V'] != -1.0 and r['MAG_V'] != '-'
                has_r = r['MAG_R'] != -1.0 and r['MAG_R'] != '-'
                if has_b:
                    compColorCount['b'] += 1
                if has_v:
                    compColorCount['v'] += 1
                if has_r:
                    compColorCount['r'] += 1
                if has_b and has_v:
                    compColorCount['bv'] += 1
                if has_v and has_r:
                    compColorCount['vr'] += 1
                if has_b and has_v and has_r:
                    compColorCount['bvr'] += 1
            elif r['ROLE'] == 'V':
                varCount += 1
        
        if varCount == 0:
            printWarning('No variable star in reference catalog %s' % (refcatFileName))
            return False
        if compCount == 0:
            printWarning('No comp star in reference catalog %s' % (refcatFileName))
            return False
        if compColorCount[colors] == 0:
            if colors == 'bvr':
                if compColorCount['bv'] > 0:
                    self.opt['color'].remove('Ri')
                    self.opt['useStd'] = False
                    printWarning('No R comp stars ; downgrade colors to BV and disable std photometry.')
                    return True
                elif compColorCount['vr'] > 0:
                    self.opt['color'].remove('Bi')
                    self.opt['useStd'] = False
                    printWarning('No B comp stars ; downgrade colors to VR and disable std photometry.')
                    return True
                elif compColorCount['v'] > 0:
                    self.opt['color'].remove('Bi')
                    self.opt['color'].remove('Ri')
                    self.opt['useStd'] = False
                    printWarning('No B and R comp stars ; downgrade colors to V and disable std photometry.')
                    return True
            printWarning('No comp stars to achieve %s photometry with reference catalog %s' % (colors.upper(), refcatFileName))
            # TODO: more color downgrade rule
            return False
        elif compColorCount[colors] < 3:
            printWarning('Only 1 comp star found to achieve %s photometry with reference catalog %s ; only the comp star method is possible' % (colors.upper(), refcatFileName))
            self.opt['method'] = 'comp'
            printWarning('Downgrade method to best comp star.')
            return True
        elif compColorCount[colors] < 5 and (self.opt['method'] != 'comp' or self.opt['makeStd'] or self.opt['adhocStd']):
            printWarning('Only %d comp star found to achieve %s photometry with reference catalog %s; do ensemble or ad-hoc standardization carefully' % (compColorCount[colors], colors.upper(), refcatFileName))
            return False

        return True

    def photometry(self, seqFolder, photFolder, refcatFileName, color):

        # Names of the sequence/combined files:
        SEQLIST = glob(seqFolder + '/' + self.pplSetup['SEQ_FILE_PREFIX'] + '*-' + color + '.fits')
        SEQLIST.sort()
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
            if not exists(AST_FILE) or self.opt['overwrite']:
                print("%s/solve-field %s -D %s -N %s %s" % (self.pplSetup['AST_BIN_FOLDER'], self.SOLVE_ARGS, photFolder, AST_FILE, f))
                invoke("%s/solve-field %s -D %s -N %s %s" % (self.pplSetup['AST_BIN_FOLDER'], self.SOLVE_ARGS, photFolder, AST_FILE, f))
                fsolved = PMCAT_FILE.replace('.cat', '.solved')
                if not exists(fsolved):
                    printError("Astrometry of %s failed." % (fsolved))
                    break
            else:
                print("astrometry file %s already exists." % (AST_FILE))

            printInfo("Make photometry for %s" % (f))
            if not exists(PMCAT_FILE) or self.opt['overwrite']:
                invoke("sextractor %s -c %s -CATALOG_NAME %s -CATALOG_TYPE ASCII_HEAD" % (AST_FILE, self.SEX_CFG, PMCAT_FILE))
            else:
                print("photometry file %s already exists." % (PMCAT_FILE))

            PMCAT_FILE_FLT = PMCAT_FILE + ".cat"
            printInfo("Filtering result catalog to %s" % (PMCAT_FILE_FLT))
            if not exists(PMCAT_FILE_FLT) or self.opt['overwrite']:
                opt = {
                    'ref' : refcatFileName,
                    'out' : PMCAT_FILE_FLT,
                    'color': color,
                    'files': [ PMCAT_FILE ],
                    }
                matcher = CatalogMatcher(opt)
                matcher.process()

            else:
                print("filtered photometry file %s already exists." % (PMCAT_FILE_FLT))

        return True

    def calculateMags(self, photFolder, colors):

        # Names of the sequence/combined files:
        color = colors[0]
        PHOTLIST = glob(photFolder + '/' + self.pplSetup['SEQ_FILE_PREFIX'] + '*-' + color + '.cat.cat')
        PHOTLIST.sort()
        if len(PHOTLIST) == 0:
            PHOTLIST = glob(photFolder + '/Combined-' + color + '.cat.cat')
#            print(PHOTLIST)
            if len(PHOTLIST) == 0:
                printWarning("No files for calculate magnitudes in folder %s" % (photFolder))
                return False

        for f in PHOTLIST:
            inputFiles = []
            for c in colors:
                inputFiles.append(f.replace('-' + color, '-' + c))

            printInfo("Calculate real magnitudes to %s" % (f.replace('-' + color, '-*') + '.pm'))
            pmopt = {
                'out' : None,
                'comp': None,
                'color' : colors,
                'method' : self.opt['method'],
                'loadCoeffs' : self.opt['useStd'],
                'useCoeffs': self.opt['useStd'] or self.opt['adhocStd'],
                'makeCoeffs': self.opt['makeStd'] or self.opt['adhocStd'],
                'saveCoeffs': self.opt['makeStd'],
                'showGraphs': self.opt['showGraphs'],
                'observer': self.opt['nameCode'],
                'files': inputFiles,
            }
            phot = Photometry(pmopt, self.pplSetup)
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

        # ## Common arguments (saturation level, image section & trimming, etc.):
        self.AST_CFG = self.pplSetup['CONFIG_FOLDER'] + "/astrometry.cfg"
        self.SEX_CFG = self.pplSetup['CONFIG_FOLDER'] + "/sex.cfg"
        # SOLVE_ARGS="-O --config $AST_CFG --use-sextractor --sextractor-path sextractor -i ${PHOTDIR}/scamp.cat -n ${PHOTDIR}/scamp.cfg -j ${PHOTDIR}/scamp.ref -r -y -p"
        self.SOLVE_ARGS = "-O --config %s --use-sextractor --sextractor-path sextractor -r -y -p" % (self.AST_CFG)

        ####################################
        # step 1. do photometry for all file
        ####################################

        PMERROR = False
        requestedColors = self.opt['color']
        requestedStd = self.opt['useStd']
        for sf in seqFolders:

            baseFolder = sf.replace('/' + self.pplSetup['SEQ_FOLDER_NAME'], '')
            self.opt['color'] = requestedColors
            self.opt['useStd'] = requestedStd
            refcatAvailable = self.inspectRefCat(baseFolder)
            if not refcatAvailable:
                printError('Reference catalog is not usable for photometry ; skip folder %s' % (sf))
                continue

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
                REPORT_FOLDER = BASE_FOLDERS[0]
            else:
                PM_FILES = glob(self.pplSetup['PHOT_FOLDER_NAME'] + '/*.pm')
                REPORT_FOLDER = getcwd()

            printInfo("Create report into %s folder." % (REPORT_FOLDER))
            opt = {
                'out' : REPORT_FOLDER,  # output folder
                'rpt' : 'aavso',  # report format, default: aavso extended
                'name': self.opt['nameCode'],  # observer name code
                'method' : self.opt['method'], # mg calculation method, comp - single com star, gcx/lfit - ensemble
                'files': PM_FILES,
                }
            proc = ReportProcessor(opt)
            proc.process()

        print()
        print(Blue + "Photometry is ready." + Color_Off)


class MainApp:

    opt = {
        'color' : ['Gi'],  # photometry band, mandatory
        'nameCode' : None,  # observer code for the AAVSO report, mandatory
        'makeStd'  : False,  # make std coeffs
        'useStd'   : False,  # use std coeffs
        'adhocStd' : False,  # make and use std coeffs for this image only
        'showGraphs': False, # show standard coefficient graphs
        'method'   : 'gcx',  # mg calculation method: comp, gcx, lfit
        'camera'   : None,   # camera name
        'telescope': None,   # telescope name
        'overwrite': False,  # force to overwrite existing results, optional
        'files': None,
        'baseFolder': None,  # base folder, optional
        }

    availableBands = ['gi', 'g', 'bi', 'b', 'ri', 'r', 'all']

    mgCalcMethods = {
        'comp': 'Best comparision star',
        'gcx' : 'GCX''s robust averaging ensemble',
        'lfit': 'Linear fit ensemble',
        }

    def __init__(self, argv):
        self.argv = argv
        pass

    def printTitle(self):
        print()
        print(BGreen + "ppl-photometry, version 1.1.0 " + Color_Off)
        print(Blue + "Make photometry on calibrated FITS images." + Color_Off)
        print()

    def usage(self):
        print("Usage: ppl-photometry [OPTIONS]... [BASE_FOLDER]")
        print("Make photometry on calibrated FITS images.")
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -c,  --color arg        set filter(s), arg is the color code, default color is 'Gi', for available color codes see below")
        print("  -n,  --name nameCode    set observer code for the AAVSO report")
        print("  -t,  --method method    magnitude calculation method ; values are: comp, gcx, lfit")
        print("  -h,  --help             print this page")
        print("  -w,  --overwrite        force to overwrite existing results")
        print("standardization:")
        print("  -m,  --make-std         create standard coefficients from a Standard Area and save them (for all color photometry)")
        print("  -s,  --use-std          use standard coefficients ; calculate standard magnitudes (for all color photometry)")
        print("  -a,  --adhoc-std        create standard coefficients and use them for calculate standard magnitudes (for all color photometry)")
        print("       --show-graph       show ensemble or standard coefficient graphs for diagnostic or illustration purpose")
        print("       --camera           set camera name ; this overrides DEF_CAMERA settings in ppl.cfg, but does not override the INSTRUME FITS header value")
        print("       --telescope        set telescope name ; this overrides DEF_TELESCOPE settings in ppl.cfg, but does not override the TELESCOP FITS header value")
        print()
        print("Available filter color codes are:")
        print("  Gi | G | gi | g         green channel")
        print("  Bi | B | bi | b         blue channel")
        print("  Ri | R | ri | r         red channel")
        print("  all | ALL | All         all channels, results 3 separate frame")
        print()

    def processCommands(self):
        try:
            optlist, args = getopt (self.argv[1:], "c:n:msat:wh", ['color=', 'name=', 'make-std', 'use-std', 'adhoc-std', 'show-graph', 'method=', 'overwrite', 'help'])
        except GetoptError:
            printError('Invalid command line options.')
            return

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-c' or o == '--color':
                color = a.lower()
                if not color in self.availableBands:
                    printError('Invalid color: %s, use on of these: Gi, g, Bi, b, Ri, r, all' % (a))
                    exit(1)
                if color == 'all':
                    self.opt['color'] = ['Ri', 'Gi', 'Bi']
                else:
                    self.opt['color'] = [a]
            elif o == '-n' or o == '--name':
                self.opt['nameCode'] = a.upper()
            elif o == '-m' or o == '--make-std':
                self.opt['makeStd'] = True
            elif o == '-s' or o == '--use-std':
                self.opt['useStd'] = True
            elif o == '-a' or o == '--adhoc-std':
                self.opt['adhocStd'] = True
            elif o == '--show-graph':
                self.opt['showGraphs'] = True
            elif o == '--camera':
                self.opt['camera'] = a
            elif o == '--telescope':
                self.opt['telescope'] = a
            elif o == '-t' or o == '--method':
                if not a in self.mgCalcMethods.keys():
                    printWarning('Invalid mg calculation method %s ; use gcx instead.')
                else:
                    self.opt['method'] = a

            elif o == '-w' or o == '--overwrite':
                self.opt['overwrite'] = True
            elif o == '-h' or o == '--help':
                self.usage()
                exit(0)

        if self.opt['adhocStd'] and (self.opt['makeStd'] or self.opt['useStd']):
            printWarning('Both -a and either -m or -s option cannot be used at once.')
            exit(0)

        if len(args) > 0:
            self.opt['baseFolder'] = args[0]
            if args[0].endswith('/'):
                self.opt['baseFolder'] = args[0][:-1]

        if self.opt['nameCode'] == None:
            printWarning('No observer code was given. Use \'XYZ\' instead.')
            self.opt['nameCode'] = 'XYZ'

        print('Mg calculation method: ' + self.mgCalcMethods[self.opt['method']])

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
