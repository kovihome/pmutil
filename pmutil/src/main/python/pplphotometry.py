#!/usr/bin/env python3
#
# PmUtils/pplphotometry
#
"""
Created on Feb 22, 2020

@author: kovi
"""
from getopt import getopt, GetoptError
from sys import argv
from datetime import datetime
from glob import glob
from os import getcwd, mkdir
from os.path import exists, basename
from astropy.table import Table

import pmbase as pm
from pmphot import Photometry
from pmresult import ReportProcessor
from pmmerge import CatalogMatcher


class Pipeline:
    opt = {}  # command line options

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
        print(f'Inspect refcat {refcatFileName}')
        if not exists(refcatFileName):
            pm.printWarning(f"Reference catalog {refcatFileName} is missing.")
            return False
        rc = Table.read(refcatFileName, format='ascii')
        # print(rc)
        varCount = 0
        compCount = 0
        compColorCount = {'bvr': 0, 'bv': 0, 'vr': 0, 'b': 0, 'v': 0, 'r': 0}
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
            pm.printWarning(f'No variable star in reference catalog {refcatFileName}')
            return False
        if compCount == 0:
            pm.printWarning(f'No comp star in reference catalog {refcatFileName}')
            return False
        if compColorCount[colors] == 0:
            if colors == 'bvr':
                if compColorCount['bv'] > 0:
                    self.opt['color'].remove('Ri')
                    self.opt['useStd'] = False
                    pm.printWarning('No R comp stars ; downgrade colors to BV and disable std photometry.')
                    return True
                elif compColorCount['vr'] > 0:
                    self.opt['color'].remove('Bi')
                    self.opt['useStd'] = False
                    pm.printWarning('No B comp stars ; downgrade colors to VR and disable std photometry.')
                    return True
                elif compColorCount['v'] > 0:
                    self.opt['color'].remove('Bi')
                    self.opt['color'].remove('Ri')
                    self.opt['useStd'] = False
                    pm.printWarning('No B and R comp stars ; downgrade colors to V and disable std photometry.')
                    return True
            pm.printWarning(
                f'No comp stars to achieve {colors.upper()} photometry with reference catalog {refcatFileName}')
            # TODO: more color downgrade rule
            return False
        elif compColorCount[colors] < 3:
            pm.printWarning(
                    f'Only 1 comp star found to achieve {colors.upper()} photometry with reference catalog {refcatFileName} ; only the comp star method is possible')
            self.opt['method'] = 'comp'
            pm.printWarning('Downgrade method to best comp star.')
            return True
        elif compColorCount[colors] < 5 and (
                self.opt['method'] != 'comp' or self.opt['makeStd'] or self.opt['adhocStd']):
            pm.printWarning(
                    f'Only {compColorCount[colors]:d} comp star found to achieve {colors.upper()} photometry with reference catalog {refcatFileName}; do ensemble or ad-hoc standardization carefully')
            return False

        return True

    def do_astrometry(self, seqFile, astFile, photFolder):
        pm.printInfo(f"Make astrometry for {seqFile}")
        if not exists(astFile) or self.opt['overwrite']:
            # print("%s/solve-field %s -D %s -N %s %s" % (pm.setup['AST_BIN_FOLDER'], self.SOLVE_ARGS, photFolder, AST_FILE, f))
            # pm.invoke("%s/solve-field %s -D %s -N %s %s" % (pm.setup['AST_BIN_FOLDER'], self.SOLVE_ARGS, photFolder, AST_FILE, f))
            astCommand = f"solve-field {self.SOLVE_ARGS} -D {photFolder} -N {astFile} {seqFile}"
            print(astCommand)
            pm.invoke(astCommand)
            fsolved = photFolder + '/' + basename(seqFile).replace('.fits', '.solved')
            if not exists(fsolved):
                pm.printError(f"Astrometry of {fsolved} failed.")
                return False
        else:
            print(f"astrometry file {astFile} already exists.")
        return True

    def do_photometry(self, astFile, pmcatFile):
        pm.printInfo(f"Make photometry for {astFile}")
        if not exists(pmcatFile) or self.opt['overwrite']:
            pmCommand = f"sextractor {astFile} -c {self.SEX_CFG} -CATALOG_NAME {pmcatFile} -CATALOG_TYPE ASCII_HEAD"
            print(pmCommand)
            pm.invoke(pmCommand)
        else:
            print(f"photometry file {pmcatFile} already exists.")

    def do_merge(self, baseFile, photFolder, colors):
        combinedCatalog = f"{photFolder}/{baseFile}.cmb"
        pm.printInfo(f"Combine result catalogs to {combinedCatalog}")
        if not exists(combinedCatalog) or self.opt['overwrite']:
            matcher = CatalogMatcher(baseFile, photFolder, colors, self.opt['showGraphs'], self.opt['saveGraphs'],
                                     self.opt['logMode'])
            matcher.process()
        else:
            print(f"filtered photometry file {combinedCatalog} already exists.")

    #    def photometry(self, seqFolder, photFolder, refcatFileName, color):

    # Names of the sequence/combined files:
    #        SEQLIST = glob(seqFolder + '/' + pm.setup['SEQ_FILE_PREFIX'] + '*-' + color + '.fits')
    #        SEQLIST.sort()
    #        if len(SEQLIST) == 0:
    #            SEQLIST = glob(seqFolder + '/Combined-' + color + '.fits')
    #            if len(SEQLIST) == 0:
    #                pm.printWarning("No files for photometry in folder %s" % (seqFolder))
    #                return False

    #        for f in SEQLIST:

    #            print(f"photometry of {f}")

    #            AST_FILE = photFolder + '/' + basename(f).replace('.fits', '.ast.fits')
    #            PMCAT_FILE = photFolder + '/' + basename(f).replace('.fits', '.cat')

    # make astrometry
    #            self.do_astrometry(f, AST_FILE, photFolder)

    # make photometry

    #            self.do_photometry(AST_FILE, PMCAT_FILE)

    #            PMCAT_FILE_FLT = PMCAT_FILE + ".cat"
    #            pm.printInfo("Filtering result catalog to %s" % (PMCAT_FILE_FLT))
    #            if not exists(PMCAT_FILE_FLT) or self.opt['overwrite']:
    #                opt = {
    #                    'ref' : refcatFileName,
    #                    'out' : PMCAT_FILE_FLT,
    #                    'color': color,
    #                    'files': [ PMCAT_FILE ],
    #                    }
    #                matcher = CatalogMatcher(opt)
    #                matcher.process()

    #            else:
    #                print("filtered photometry file %s already exists." % (PMCAT_FILE_FLT))

    #        return True

    def calculateMags(self, baseName, photFolder, colors):

        # Names of the sequence/combined files:
        #        color = colors[0]
        #        PHOTLIST = glob(photFolder + '/' + pm.setup['SEQ_FILE_PREFIX'] + '*-' + color + '.cmb')
        #        PHOTLIST.sort()
        #        if len(PHOTLIST) == 0:
        #            PHOTLIST = glob(f"{photFolder}/{baseName}-{color}.cmb")
        #            if len(PHOTLIST) == 0:
        #                pm.printWarning("No files for calculate magnitudes in folder %s" % (photFolder))
        #                return False

        #        for f in PHOTLIST:
        #            inputFiles = []
        #            for c in colors:
        #                inputFiles.append(f.replace('-' + color, '-' + c))

        cmbFileName = f"{photFolder}/{baseName}.cmb"

        #        pm.printInfo("Calculate real magnitudes to %s" % (f.replace('-' + color, '-*') + '.pm'))
        pm.printInfo(f"Calculate real magnitudes to {cmbFileName}.pm")
        pmopt = {
            'out': None,
            'comp': None,
            'color': colors,
            'method': self.opt['method'],
            'loadCoeffs': self.opt['useStd'],
            'useCoeffs': self.opt['useStd'] or self.opt['adhocStd'],
            'makeCoeffs': self.opt['makeStd'] or self.opt['adhocStd'],
            'saveCoeffs': self.opt['makeStd'],
            'showGraphs': self.opt['showGraphs'],
            'saveGraphs': self.opt['saveGraphs'],
            'observer': self.opt['nameCode'],
            'files': [cmbFileName],
        }
        phot = Photometry(pmopt)
        phot.process()

    def process_photometry(self, seqFolder, photFolder, title):

        pm.printInfo(f"{title}: Make photometry on sequence file(s).")

        # Create the photometry dir, if not exists
        if not exists(photFolder):
            mkdir(photFolder)

        # Photometry reference file:
        sf = seqFolder
        basedir = sf.replace(pm.setup['SEQ_FOLDER_NAME'], '')
        PMREFS = glob(basedir + '/*.cat')
        if len(PMREFS) <= 0 or not exists(PMREFS[0]):
            pm.printError(f"No reference catalog file (.cat) in folder {basedir}")
            pm.printInfo("Use ppl-refcat command to create reference catalog for the object.")
            return False
        # else:
        #     PMREF = PMREFS[0]

        # collect base sequence file names
        baseSeqFileNames = self.collectBaseFileNames(sf)

        for baseFile in baseSeqFileNames:
            print(f"Do photometry on {baseFile}-*.fits")

            for color in self.opt['color']:
                seqFile = f"{seqFolder}/{baseFile}-{color}.fits"
                astFile = f"{photFolder}/{baseFile}-{color}.ast.fits"
                # pmcatFile = f"{photFolder}/{baseFile}-{color}.cat"

                self.do_astrometry(seqFile, astFile, photFolder)

                pmcatFile = f"{photFolder}/{baseFile}-{color}.cat"

                self.do_photometry(astFile, pmcatFile)

            baseFolder = seqFolder.replace('/' + pm.setup['SEQ_FOLDER_NAME'], '')
            self.do_merge(baseFile, baseFolder, self.opt['color'])

            self.calculateMags(baseFile, photFolder, self.opt['color'])

        #            self.photometry(seqFolder, photFolder, PMREF, color)

        #        self.calculateMags(photFolder, self.opt['color'])

        # TODO: cleanup

        return True

    def collectBaseFileNames(self, seqFolder):
        s = glob(f"{seqFolder}/*.fits")
        b = [basename(x).split('-')[0] for x in s]
        return list(set(b))

    def execute(self):

        ##########################
        # step 0. setup photometry
        ##########################

        seqFolders = pm.discoverFolders(self.opt['baseFolder'], pm.setup['SEQ_FOLDER_NAME'])
        if len(seqFolders) == 0:
            pm.printError('No sequence folder found in base folder %s' % (self.opt['baseFolder']))
            exit(1)
        print("Sequence folders discovered:", seqFolders)

        # ## Common arguments (saturation level, image section & trimming, etc.):
        self.AST_CFG = pm.setup['CONFIG_FOLDER'] + "/astrometry.cfg"
        self.SEX_CFG = pm.setup['CONFIG_FOLDER'] + "/sex.cfg"
        # SOLVE_ARGS="-O --config $AST_CFG --use-sextractor --sextractor-path sextractor -i ${PHOTDIR}/scamp.cat -n ${PHOTDIR}/scamp.cfg -j ${PHOTDIR}/scamp.ref -r -y -p"
        self.SOLVE_ARGS = f"-O --config {self.AST_CFG} --use-source-extractor --source-extractor-path sextractor -r -y -p"

        ####################################
        # step 1. do photometry for all file
        ####################################

        PMERROR = False
        requestedColors = self.opt['color']
        requestedStd = self.opt['useStd']
        for sf in seqFolders:

            baseFolder = sf.replace('/' + pm.setup['SEQ_FOLDER_NAME'], '')
            self.opt['color'] = requestedColors
            self.opt['useStd'] = requestedStd
            refcatAvailable = self.inspectRefCat(baseFolder)
            if not refcatAvailable:
                pm.printError(f'Reference catalog is not usable for photometry ; skip folder {sf}')
                continue

            pf = sf.replace(pm.setup['SEQ_FOLDER_NAME'], pm.setup['PHOT_FOLDER_NAME'])

            success = self.process_photometry(sf, pf, "PHOTOMETRY")
            if not success:
                PMERROR = True
                break

        ########################################
        # step 2. create report from all results
        ########################################

        if not PMERROR:
            if self.opt['baseFolder'] is not None:
                PM_FILES = glob('*' + self.opt['baseFolder'] + '*/' + pm.setup['PHOT_FOLDER_NAME'] + '/*.cmb.pm')
                BASE_FOLDERS = glob('*' + self.opt['baseFolder'] + '*')
                REPORT_FOLDER = BASE_FOLDERS[0]
            else:
                PM_FILES = glob(pm.setup['PHOT_FOLDER_NAME'] + '/*.cmb.pm')
                # PM_FILES = [ f"{pm.setup['PHOT_FOLDER_NAME']}/{baseName}.cmb.pm" for baseName in self.collectBaseFileNames(pm.setup['PHOT_FOLDER_NAME'])  ]
                REPORT_FOLDER = getcwd()

            pm.printInfo(f"Create report into {REPORT_FOLDER} folder.")
            opt = {
                'out': REPORT_FOLDER,  # output folder
                'rpt': 'aavso',  # report format, default: aavso extended
                'name': self.opt['nameCode'],  # observer name code
                'method': self.opt['method'],  # mg calculation method, comp - single com star, gcx/lfit - ensemble
                'showGraphs': self.opt['showGraphs'],
                'saveGraphs': self.opt['saveGraphs'],
                'files': PM_FILES,
            }
            proc = ReportProcessor(opt)
            proc.process()

        print()
        print(pm.Blue + "Photometry is ready." + pm.Color_Off)


class MainApp:
    opt = {
        'color': ['Gi'],  # photometry band, mandatory
        'nameCode': None,  # observer code for the AAVSO report, mandatory
        'makeStd': False,  # make std coeffs
        'useStd': False,  # use std coeffs
        'adhocStd': False,  # make and use std coeffs for this image only
        'showGraphs': False,  # show standard coefficient graphs
        'saveGraphs': False,  # save graphs
        'method': 'gcx',  # mg calculation method: comp, gcx, lfit
        'camera': None,  # camera name
        'telescope': None,  # telescope name
        'overwrite': False,  # force to overwrite existing results, optional
        'logMode': pm.Logger.LOG_MODE_INFO,  # logging mode
        'files': None,
        'baseFolder': None,  # base folder, optional
    }

    availableBands = ['gi', 'g', 'bi', 'b', 'ri', 'r', 'all']

    mgCalcMethods = {
        'comp': 'Best comparision star',
        'gcx': 'GCX''s robust averaging ensemble',
        'lfit': 'Linear fit ensemble',
    }

    def __init__(self, argv):
        self.argv = argv

    def printTitle(self):
        print()
        print(pm.BGreen + "ppl-photometry, version 1.2.0 " + pm.Color_Off)
        print(pm.Blue + "Make photometry on calibrated FITS images." + pm.Color_Off)
        print()

    def usage(self):
        print("Usage: ppl-photometry [OPTIONS]... [BASE_FOLDER]")
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print(
            "  -c,  --color arg        set filter(s), arg is the color code, default color is 'Gi', for available color codes see below")
        print("  -n,  --name nameCode    set observer code for the AAVSO report")
        print("  -t,  --method method    magnitude calculation method ; values are: comp, gcx, lfit")
        print("  -h,  --help             print this page")
        print("  -w,  --overwrite        force to overwrite existing results")
        print("       --debug            print debug info")
        print("standardization:")
        print(
            "  -m,  --make-std         create standard coefficients from a Standard Area and save them (for all color photometry)")
        print(
            "  -s,  --use-std          use standard coefficients ; calculate standard magnitudes (for all color photometry)")
        print(
            "  -a,  --adhoc-std        create standard coefficients and use them for calculate standard magnitudes (for all color photometry)")
        print("       --show-graph       show graphs or plots for diagnostic or illustration purpose")
        print("       --save-graph       save graphs or plots for diagnostic or illustration purpose")
        print(
            "       --camera           set camera name ; this overrides DEF_CAMERA settings in ppl.cfg, but does not override the INSTRUME FITS header value")
        print(
            "       --telescope        set telescope name ; this overrides DEF_TELESCOPE settings in ppl.cfg, but does not override the TELESCOP FITS header value")
        print()
        print("Available filter color codes are:")
        print("  Gi | G | gi | g         green channel")
        print("  Bi | B | bi | b         pm.Blue channel")
        print("  Ri | R | ri | r         red channel")
        print("  all | ALL | All         all channels, results 3 separate frame")
        print()

    def processCommands(self):
        try:
            optlist, args = getopt(self.argv[1:], "c:n:msat:wdh",
                                   ['color=', 'name=', 'make-std', 'use-std', 'adhoc-std', 'show-graph', 'save-graph',
                                    'method=', 'overwrite', 'debug', 'help'])
        except GetoptError:
            pm.printError('Invalid command line options.')
            return

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-c' or o == '--color':
                color = a.lower()
                if color not in self.availableBands:
                    pm.printError(f'Invalid color: {a}, use on of these: Gi, g, Bi, b, Ri, r, all')
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
            elif o == '--save-graph':
                self.opt['saveGraphs'] = True
            elif o == '--camera':
                self.opt['camera'] = a
            elif o == '--telescope':
                self.opt['telescope'] = a
            elif o == '-t' or o == '--method':
                if not a in self.mgCalcMethods.keys():
                    pm.printWarning('Invalid mg calculation method %s ; use gcx instead.')
                else:
                    self.opt['method'] = a

            elif o == '-w' or o == '--overwrite':
                self.opt['overwrite'] = True
            elif o == '--debug':
                self.opt['logMode'] = pm.Logger.LOG_MODE_DEBUG
            elif o == '-h' or o == '--help':
                self.usage()
                exit(0)

        if self.opt['adhocStd'] and (self.opt['makeStd'] or self.opt['useStd']):
            pm.printWarning('Both -a and either -m or -s option cannot be used at once.')
            exit(0)

        if len(args) > 0:
            self.opt['baseFolder'] = args[0]
            if args[0].endswith('/'):
                self.opt['baseFolder'] = args[0][:-1]

        if self.opt['nameCode'] is None:
            pm.printWarning('No observer code was given. Use \'XYZ\' instead.')
            self.opt['nameCode'] = 'XYZ'

        print('Mg calculation method: ' + self.mgCalcMethods[self.opt['method']])

    def run(self):
        self.printTitle()
        self.processCommands()
        pm.saveCommand(self.opt['baseFolder'], self.argv, 'photometry')

        start = datetime.now()

        ppl = Pipeline(self.opt)
        ppl.execute()

        exectime = (datetime.now() - start).total_seconds()
        print("%sexecution time was %d seconds.%s" % (pm.Blue, exectime, pm.Color_Off))


if __name__ == '__main__':
    app = MainApp(argv)
    app.run()

# end __main__
