#!/usr/bin/env python3
#
# PmUtils/pplcalibration
#
import os
'''
Created on Mar 1, 2020

@author: kovi
'''

from getopt import getopt, GetoptError
from sys import argv
from datetime import datetime, timedelta
import glob
from os import makedirs
from os.path import isdir, exists, basename
from shutil import copyfile
import time

from pmbase import printError, printWarning, printInfo, saveCommand, loadPplSetup, invoke, Blue, Color_Off, BGreen, getFitsHeaders, getFitsHeader, invokep, setFitsHeaders, findInFile
from pmdisco import Discovery


class Pipeline:

    opt = {}  # command line options
    pplSetup = {}  # PPL setup from ppl-setup config file

    # Common arguments (saturation level, image section & trimming, etc.):
    COMMON_ARGS = "--saturation 16000 --trim"
    FISTAR_ARGS = "--algorithm uplink --prominence 0.0 --model elliptic --format id,x,y,s,d,k,amp,flux"
    GRMATCH_ARGS = "--col-ref 2,3 --col-ref-ordering +8 --col-inp 2,3 --col-inp-ordering +8 --weight reference,column=4,magnitude,power=2 --triangulation maxinp=100,maxref=100,conformable,auto,unitarity=0.002 --order 2 --max-distance 3 --comment"
    FICOMBINE_ARGS = "-m sum -n"

    TEMPDIR = "temp"

    BIAS_FOLDER = None
    DARK_FOLDER = None
    FLAT_FOLDER = None
    FLAT_BIAS_FOLDER = None
    FLAT_DARK_FOLDER = None
    LIGHT_FOLDERS = None

    def __init__(self, opt):
        self.opt = opt

    def discoverFolders(self):
        disco = Discovery(self.opt['flatOnly'], self.opt['useMasterFlat'], self.opt['baseFolder'], self.pplSetup)
        disco.discover()
        self.BIAS_FOLDER = disco.BIAS_FOLDER
        self.DARK_FOLDER = disco.DARK_FOLDER
        self.FLAT_FOLDER = disco.FLAT_FOLDER
        self.FLAT_BIAS_FOLDER = disco.FLAT_BIAS_FOLDER
        self.FLAT_DARK_FOLDER = disco.FLAT_DARK_FOLDER
        self.LIGHT_FOLDERS = disco.LIGHT_FOLDERS

    def imSizeArg(self, fitsFileName):
        h = getFitsHeaders(fitsFileName, ['NAXIS1', 'NAXIS2'])
        return "--image 0:0:%d:%d" % (int(h['NAXIS1']) - 1, int(h['NAXIS2']) - 1)

    def lt2ut(self, utcDate):
        if utcDate:
            s = utcDate
            t = time.strptime(s, '%Y-%m-%dT%H:%M:%S')
            d = datetime(*time.gmtime(time.mktime(t))[:6])
        else:
            d = datetime.utcnow()
        if d.microsecond > 500000:
            td = timedelta(microseconds = 1000000 - d.microsecond)
            d = d + td
        else:
            td = timedelta(microseconds = d.microsecond)
            d = d - td
        return d.isoformat()

    def sumDateObs(self, targetFile, sourceFiles):
        e_sum = 0
        m_sum = 0
        for sourceFile in sourceFiles:
            h = getFitsHeaders(sourceFile, ['DATE-OBS', 'EXPTIME'])
            t = time.mktime(time.strptime(h['DATE-OBS'], '%Y-%m-%dT%H:%M:%S'))
            e = int(h['EXPTIME'])
            e_sum = e_sum + e
            m_sum = m_sum + (t + e / 2) * e

        m_avg = m_sum / e_sum
        d = datetime(*time.gmtime(m_avg)[:6])
        setFitsHeaders(targetFile, { 'EXPTIME': e_sum, 'DATE-MID': (d.isoformat(), 'Effective center of time of observation')  })
        print('observation center time:', d.isoformat(), 'cumulated exptime:', str(e_sum), 'sec')

    def raw2fitsFile(self, rawfile, color):
        FITS_NAME = "%s-%s.fits" % (rawfile[:rawfile.rfind('.')], color)
        if self.opt['overwrite'] or not exists(FITS_NAME):
            print("%s -> %s" % (rawfile, FITS_NAME))

            # convert raw image to fits
            os.system("rawtran -c %s -o %s -B 16 -C '-4 -D' -X '-q 3 -w' %s" % (color, FITS_NAME, rawfile))

            # read image temperature from raw file
            exif = invoke("exiftool -s -g %s" % (rawfile))
            IMAGETEMP = None
            for line in exif.split('\n'):
                if line.startswith('CameraTemperature'):
                    IMAGETEMP = line.split(':')[1].strip()
                    break

            if IMAGETEMP:
                # write image temperature to fits file
                setFitsHeaders(FITS_NAME, { 'CCD-TEMP': (IMAGETEMP, 'CCD Temperature (Celsius)') })

            # convert observation time from local time to UT
            if self.opt['imageTime'] == 'LT':
                IMAGE_DATE = getFitsHeader(FITS_NAME, 'DATE-OBS')
                IMAGE_DATE_UTC = self.lt2ut(IMAGE_DATE)
                setFitsHeaders(FITS_NAME, { 'DATE-OBS': IMAGE_DATE_UTC, 'DATE-IMG': (IMAGE_DATE, 'Original image date (in local time)') })

        else:
            print("%s file already exists." % (FITS_NAME))

    def raw2fits(self, folder):
        # list of raw CR2 file in the input directory
        rawFiles = glob.glob(folder + "/*.cr2")
        rawFiles.sort()

        # for all raw file do the conversion
        for rawfile in rawFiles:
            for color in self.opt['color']:
                self.raw2fitsFile(rawfile, color)

    def makeMasterBias(self, biasFolder, color):
        BIAS_PATTERN = "%s/%s*-%s.fits" % (biasFolder, self.pplSetup['BIAS_FILE_PREFIX'], color)
        BIASLIST = glob.glob(BIAS_PATTERN)
        BIASLIST.sort()

        masterBiasFile = "%s/%s-%s.fits" % (biasFolder, self.pplSetup['MASTER_BIAS_FILE'], color)

        if len(BIASLIST) == 0:
            return 1

        print("%s -> %s" % (BIAS_PATTERN, masterBiasFile))

        # Calibrated images: all the images have the same name but put into a separate directory ($TARGET):
        R_BIASLIST = list(map(lambda x: self.TEMPDIR + '/' + basename(x), BIASLIST))

        # The calibration of the individual bias frames, followed by their combination into a single master image:
        invoke("ficalib -i %s %s %s -o %s" % (' '.join(BIASLIST), self.COMMON_ARGS, self.imSizeArg(BIASLIST[0]), ' '.join(R_BIASLIST)))
        invoke("ficombine %s --mode median -o %s" % (' '.join(R_BIASLIST), masterBiasFile))

        # cleanup: remove temp files
        for f in R_BIASLIST:
            os.remove(f)

        return 0

    def makeMasterDark(self, darkFolder, biasFolder, color):

        # Names of the individual files storing the raw bias, dark, flat and object frames:
        DARK_PATTERN = "%s/%s*-%s.fits" % (darkFolder, self.pplSetup['DARK_FILE_PREFIX'], color)
        DARKLIST = glob.glob(DARK_PATTERN)

        masterDarkFile = "%s/%s-%s.fits" % (darkFolder, self.pplSetup['MASTER_DARK_FILE'], color)
        masterBiasFile = "%s/%s-%s.fits" % (biasFolder, self.pplSetup['MASTER_BIAS_FILE'], color)

        if len(DARKLIST) == 0:
            return 1

        print("%s -> %s" % (DARK_PATTERN, masterDarkFile))

        # Calibrated images: all the images have the same name but put into a separate directory ($TARGET):
        R_DARKLIST = list(map(lambda x: self.TEMPDIR + '/' + basename(x), DARKLIST))

        # The calibration of the individual bias frames, followed by their combination into a single master image:
        invoke("ficalib -i %s %s %s -o %s --input-master-bias %s" % (' '.join(DARKLIST), self.COMMON_ARGS, self.imSizeArg(DARKLIST[0]), ' '.join(R_DARKLIST), masterBiasFile))
        invoke("ficombine %s --mode median -o %s" % (' '.join(R_DARKLIST), masterDarkFile))

        # Calculate average ccd temperature from .cr2 files, and set it into the master dark
        tsum = 0
        count = 0
        for f in DARKLIST:
            hdr = getFitsHeader(f, 'CCD-TEMP')
            if hdr:
                ccdtemp = int()
                tsum += ccdtemp
                count += 1
        if count > 0:
            AVGTEMP = (tsum + (count / 2)) / count
            print("average dark temperature: %d C" % (AVGTEMP))
            setFitsHeaders(masterDarkFile, { 'CCD-TEMP': ("%d." % (AVGTEMP), "CCD Temperature (Celsius)") })

        # cleanup: remove temp files
        for f in R_DARKLIST:
            os.remove(f)

        return 0

    def makeMasterFlat(self, flatfolder, biasFolder, darkFolder, color):
        # Names of the individual files storing the raw bias, dark, flat and object frames:
        FLAT_PATTERN = "%s/%s*-%s.fits" % (flatfolder, self.pplSetup['FLAT_FILE_PREFIX'], color)
        FLATLIST = glob.glob(FLAT_PATTERN)

        masterFlatFile = "%s/%s-%s.fits" % (flatfolder, self.pplSetup['MASTER_FLAT_FILE'], color)
        masterDarkFile = "%s/%s-%s.fits" % (darkFolder, self.pplSetup['MASTER_DARK_FILE'], color)
        masterBiasFile = "%s/%s-%s.fits" % (biasFolder, self.pplSetup['MASTER_BIAS_FILE'], color)

        if len(FLATLIST) == 0:
            return 1

        print("%s -> %s" % (FLAT_PATTERN, masterFlatFile))

        # Calibrated images: all the images have the same name but put into a separate directory ($TARGET):
        R_FLATLIST = list(map(lambda x: self.TEMPDIR + '/' + basename(x), FLATLIST))

        # The calibration of the individual flat frames, followed by their combination into a single master image:
        invoke("ficalib -i %s %s %s --post-scale 20000 -o %s --input-master-bias %s --input-master-dark %s" % (' '.join(FLATLIST), self.COMMON_ARGS, self.imSizeArg(FLATLIST[0]), ' '.join(R_FLATLIST), masterBiasFile, masterDarkFile))
        invoke("ficombine %s --mode median -o %s" % (' '.join(R_FLATLIST), masterFlatFile))

        # remove temp files
        for f in R_FLATLIST:
            os.remove(f)

        return 0

    def calibrate(self, lightFolder, calibFolder, color):
        MB = "%s/%s-%s.fits" % (self.BIAS_FOLDER, self.pplSetup['MASTER_BIAS_FILE'], color)
        MD = "%s/%s-%s.fits" % (self.DARK_FOLDER, self.pplSetup['MASTER_DARK_FILE'], color)
        MF = "%s/%s-%s.fits" % (self.FLAT_FOLDER, self.pplSetup['MASTER_FLAT_FILE'], color)

        # Names of the individual files storing the raw bias, dark, flat and object frames:
        LIGHT_PATTERN = "%s/%s*-%s.fits" % (lightFolder, self.pplSetup['LIGHT_FILE_PREFIX'], color)
        IOBJLIST = glob.glob(LIGHT_PATTERN)
        if len(IOBJLIST) == 0:
            return 1

        # echo "$1/${LIGHT_FILE_PREFIX}*-$3.fits -> $2/${LIGHT_FILE_PREFIX}*-$3.fits"
        CALIB_PATTERN = "%s/%s*-%s.fits" % (calibFolder, self.pplSetup['LIGHT_FILE_PREFIX'], color)
        print("%s -> %s" % (LIGHT_PATTERN, CALIB_PATTERN))

        # Calibrated images: all the images have the same name but put into a separate directory ($TARGET):
        R_IOBJLIST = list(map(lambda x: calibFolder + '/' + basename(x), IOBJLIST))

        # The calibration of the object images:
        invoke("ficalib -i %s %s %s -o %s --input-master-bias %s --input-master-dark %s --input-master-flat %s" % (' '.join(IOBJLIST), self.COMMON_ARGS, self.imSizeArg(IOBJLIST[0]), ' '.join(R_IOBJLIST), MB, MD, MF))

        return 0

    def registrate(self, calibFolder, seqFolder, color):
        # Names of the individual files storing the raw bias, dark, flat and object frames:
        CALIB_PATTERN = "%s/%s*-%s.fits" % (calibFolder, self.pplSetup['LIGHT_FILE_PREFIX'], color)
        COBJLIST = glob.glob(CALIB_PATTERN)
        if len(COBJLIST) == 0:
            return 1

        COMBINED_FILE = "%s/Combined-%s.fits" % (seqFolder, color)
        print("%s -> %s" % (CALIB_PATTERN, COMBINED_FILE))

        # ------------------------------------

        # make reference image from the first one
        print("set reference image: " + COBJLIST[0])
        REF = "%s/ref-%s.stars" % (self.TEMPDIR, color)
        invoke("fistar %s %s -o %s" % (COBJLIST[0], self.FISTAR_ARGS, REF))
        copyfile(COBJLIST[0], "%s/%s" % (self.TEMPDIR, basename(COBJLIST[0])))

        # Registration of the source images:
        for f in COBJLIST[1:]:

            bn = basename(f)
            pbn = self.TEMPDIR + '/' + bn

            print("transform %s -> %s" % (f, pbn))
            fstars = pbn + ".stars"
            ftrans = pbn + ".trans"

            invoke("fistar %s %s -o %s" % (f, self.FISTAR_ARGS, fstars))

            invoke("grmatch %s --input %s --reference %s --output-transformation %s" % (self.GRMATCH_ARGS, fstars, REF, ftrans))

            ln = findInFile(ftrans, 'Match failed')
            print(ln)
            matchFailed = 1 if ln != None else 0
            if matchFailed == 1:
                printError("Match failed on %s." % (f))
                if self.opt['onError'] == "stop":
                    self.exitOnError()

            if matchFailed == 0  or self.opt['onError'] == "noop":
                invoke("fitrans %s --input-transformation %s --reverse -k -o %s" % (f, ftrans, pbn))

        # list of transformed frames
        TR_PATTERN = "%s/%s*-%s.fits" % (self.TEMPDIR, self.pplSetup['LIGHT_FILE_PREFIX'], color)
        TROBJLIST = glob.glob(TR_PATTERN)

        # create a sequence of combined frames
        cc = self.opt['countCombine']
        if cc != 0:
            for a in range(0, len(TROBJLIST), cc):
                fi = "%03d" % (a / cc)
                combined = "%s/%s%s-%s.fits" % (seqFolder, self.pplSetup['SEQ_FILE_PREFIX'], fi, color)
                print("combine " + combined)
                invoke("ficombine %s %s -o %s" % (' '.join(TROBJLIST[a:a + cc]), self.FICOMBINE_ARGS, combined))
                self.sumDateObs(combined, TROBJLIST[a:a + cc])

        # create a combined frame of all images
        print("combine %s -> %s" % (TR_PATTERN, COMBINED_FILE))
        invoke("ficombine %s %s -o %s" % (' '.join(TROBJLIST), self.FICOMBINE_ARGS, COMBINED_FILE))
        self.sumDateObs(COMBINED_FILE, TROBJLIST)

        # cleanup - remove temp files
        for f in TROBJLIST:
            os.remove(f)
        os.remove(REF)
        for f in glob.glob(self.TEMPDIR + '/*.stars'):
            os.remove(f)
        for f in glob.glob(self.TEMPDIR + '/*.trans'):
            os.remove(f)

    def mastersExist(self, folder, prefix):
        for color in self.opt['color']:
            MB = "%s/%s-%s.fits" % (folder, prefix, color)
            if not exists(MB):
                return False
        return True

    def processBias(self, biasFolder, title):
        ex = self.mastersExist(biasFolder, self.pplSetup['MASTER_BIAS_FILE'])
        if self.opt['overwrite'] or not ex:
            printInfo(title + ": Convert bias RAW files to FITS.")
            self.raw2fits(biasFolder)

            printInfo(title + ": Create master bias file.")
            for color in self.opt['color']:
                self.makeMasterBias(biasFolder, color)
        else:
            printInfo(title + ": Master bias file(s) are already created.")
        # TODO: cleanup - delete bias FITS files

    def processDark(self, darkFolder, biasFolder, title):
        ex = self.mastersExist(darkFolder, self.pplSetup['MASTER_DARK_FILE'])
        if self.opt['overwrite'] or not ex:
            printInfo(title + ": Convert dark RAW files to FITS.")
            self.raw2fits(darkFolder)

            printInfo(title + ": Create master dark file(s).")
            for color in self.opt['color']:
                self.makeMasterDark(darkFolder, biasFolder, color)
        else:
            printInfo(title + ": Master dark file(s) are already created.")
        # TODO: cleanup - delete dark FITS files

    def processFlat(self, flatFolder, biasFolder, darkFolder, title):
        printInfo(title + ": Convert flat RAW files to FITS.")
        self.raw2fits(flatFolder)

        printInfo(title + ": Create master dark file(s).")
        for color in self.opt['color']:
            self.makeMasterFlat(flatFolder, biasFolder, darkFolder, color)

    def processCalibration(self, lightFolder, calibFolder, title):
        printInfo(title + ": Convert light RAW files to FITS.")
        self.raw2fits(lightFolder)

        printInfo(title + ": Create calibrated light file(s).")

        # Create dir for the calibrated images, if not exists
        if not exists(calibFolder):
            makedirs(calibFolder)

        for color in self.opt['color']:
            self.calibrate(lightFolder, calibFolder, color)
        # TODO: cleanup - delete light FITS files

    def processRegistration(self, calibFolder, seqFolder, title):

        printInfo(title + ": Create calibrated light file(s).")

        # Create the sequence dir, if not exists
        if not exists(seqFolder):
            makedirs(seqFolder)

        for color in self.opt['color']:
            self.registrate(calibFolder, seqFolder, color)
        # TODO: cleanup - delete light FITS files

    def execute(self):

        ##########################
        # step 0. setup photometry
        ##########################
        self.pplSetup = loadPplSetup()

        self.discoverFolders()

        if self.opt['masterFlat']:
            if not isdir(self.opt['masterFlat']):
                printError("Master-flat folder %s not exists." % (self.opt['masterFlat']))
                exit(1)
            else:
                hasMaster = True
                for c in self.opt['color']:
                    mfc = "%s/%s-%s.fits" % (self.opt['masterFlat'], self.pplSetup['MASTER_FLAT_FILE'], c)
                    print(mfc)
                    if not exists(mfc):
                        printError("No %s master-flat file in the directory %s" % (c, self.opt['masterFlat']))
                        hasMaster = False
                if hasMaster:
                    self.FLAT_FOLDER = self.opt['masterFlat']
                    print("master-flat: " + self.opt['masterFlat'])
                else:
                    exit(1)

        # create the temp dir, if not exists
        if not exists(self.TEMPDIR):
            makedirs(self.TEMPDIR)

        ####################################
        # step 1. create master bias frame
        ####################################

        if not self.opt['flatOnly']:

            self.processBias(self.BIAS_FOLDER, "BIAS")

        ####################################
        # step 2. create master dark frame
        ####################################

            self.processDark(self.DARK_FOLDER, self.BIAS_FOLDER, "DARK")

        if not self.opt['useMasterFlat']:

            # process flat bias, flat dark and flat, only if flat master is not exist
            ex = self.mastersExist(self.FLAT_FOLDER, self.pplSetup['MASTER_FLAT_FILE'])
            if self.opt['overwrite'] or not ex:

        ############################################################################
        # step 3. create master flat bias frame, if it is differs from master bias
        ############################################################################
                self.processBias(self.FLAT_BIAS_FOLDER, "FLAT BIAS")

        ############################################################################
        # step 4. create master flat dark frame, if it is differs from mastre dark
        ############################################################################
                self.processDark(self.FLAT_DARK_FOLDER, self.FLAT_BIAS_FOLDER, "FLAT DARK")

        ##############################
        # step 5. create master flat
        ##############################
                self.processFlat(self.FLAT_FOLDER, self.FLAT_BIAS_FOLDER, self.FLAT_DARK_FOLDER, "FLAT")

            else:
                printInfo("FLAT: Master flat file(s) are already created.")

        ##################################
        # step 6. calibrate light frames
        ##################################
        if not self.opt['flatOnly']:

            for lf in self.LIGHT_FOLDERS:

                cf = lf.replace(self.pplSetup['LIGHT_FOLDER_NAME'], self.pplSetup['CALIB_FOLDER_NAME'])

                self.processCalibration(lf, cf, "CALIBRATE")

        ###############################################
        # step 7. registrate and combine light frames
        ###############################################

                sf = lf.replace(self.pplSetup['LIGHT_FOLDER_NAME'], self.pplSetup['SEQ_FOLDER_NAME'])

                self.processRegistration(cf, sf, "REGISTRATE")

        print()
        print(Blue + "Calibration is ready." + Color_Off)


class MainApp:

    opt = {
        'color' : ['Gi'],  # photometry band, mandatory
        'countCombine' : 0,  #
        'flatOnly' : False,
        'useMasterFlat'  : False,  #
        'imageTime' : 'LT',  #
        'masterFlat' : None,  # make and use std coeffs for this image only
        'onError'   : 'noop',  # mg calculation method: comp, gcx, lfit
        'overwrite': False,  # force to overwrite existing results, optional
        'baseFolder': None,  # base folder, optional
        }

    availableBands = ['gi', 'g', 'bi', 'b', 'ri', 'r', 'all']

    def __init__(self, argv):
        self.argv = argv
        pass

    def printTitle(self):
        print()
        print(BGreen + "ppl-calibration, version 1.1.0" + Color_Off)
        print(Blue + "Calibrate a set of RAW or FITS images." + Color_Off)
        print()

    def usage(self):
        print("Usage: ppl-calibration [OPTIONS]... [BASE_FOLDER]")
        print("Make calibration process for raw or fits images.")
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -c,  --color arg        set filter(s), arg is the color code, default color is 'Gi', for available color codes see below")
        print("  -n, --count-combine n          set number of frames to combine in the sequence, 0 means all frames, default is 0")
        print("  -f,  --flat                     make master flat frame only")
        print("  -m, --master-flat folder       use the given master-flat folder")
        print("  -t,  --image-time LT|UT         specify orignal image time zone, LT=local time, UT=universal time")
        print("  -w,  --overwrite        force to overwrite existing results")
        print("  -e,  --on-error noop|skip|stop  specify what to do on error: noop=nothing to do; skip=remove the file on processing; stop=stop processing at all")
        print("  -h,  --help             print this page")
        print()
        print("Available filter color codes are:")
        print("  Gi | G | gi | g         green channel")
        print("  Bi | B | bi | b         blue channel")
        print("  Ri | R | ri | r         red channel")
        print("  all | ALL | All         all channels, results 3 separate frame")
        print()

    def processCommands(self):
        try:
            optlist, args = getopt (self.argv[1:], "c:n:fm:t:we:h", ['--color', '--count-combine', '--flat', '--master-flat', '--image-time', '--overwrite', '--on-error', '--help'])
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
            elif o == '-f':
                self.opt['flatOnly'] = True
            elif o == '-m':
                self.opt['useMasterFlat'] = True
                self.opt['masterFlat'] = a
            elif o == '-n':
                self.opt['countCombine'] = int(a)
            elif o == '-t':
                if a in ['LT', 'UT']:
                    self.opt['imageTime'] = a
                else:
                    printWarning("Bad image time zone value: %s, using default %s instead." % (a, self.opt['imageTime']))
            elif o == '-r':
                if not a in ['noop', 'skip', 'stop']:
                    printWarning("Bad on-error instruction; available values are: noop, stop, skip.")
                else:
                    self.opt['onError'] = a
            elif o == '-w':
                self.opt['overwrite'] = True
            elif o == '-h':
                self.usage()
                exit(0)

        if len(args) > 0:
            self.opt['baseFolder'] = args[0]
            if args[0].endswith('/'):
                self.opt['baseFolder'] = args[0][:-1]

            if not isdir(self.opt['baseFolder']):
                printError("Base folder %s not exists or not a directory." % (self.opt['baseFolder']))
                exit(1)
            else:
                print("base folder: %s" % (self.opt['baseFolder']))

        print("colors: " + ' '.join(self.opt['color']))
        print("overwrite: " + str(self.opt['overwrite']))
        print("count combine: " + str(self.opt['countCombine']))
        print("image timezone: " + str(self.opt['imageTime']))

    def run(self):
        self.printTitle()
        self.processCommands()
        saveCommand(self.opt['baseFolder'], self.argv, 'calibration')

        start = datetime.now()

        ppl = Pipeline(self.opt)
        ppl.execute()

        exectime = (datetime.now() - start).total_seconds()
        print("%sexecution time was %d seconds.%s" % (Blue, exectime, Color_Off))


if __name__ == '__main__':

    app = MainApp(argv)
    app.run()

# end __main__
