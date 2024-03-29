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
from os.path import isdir, exists, basename, expanduser
from shutil import copyfile
import time

import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.stats import SigmaClip
import astroalign as aa
from photutils import Background2D, MedianBackground

from pmbase import printError, printWarning, printInfo, printDebug, saveCommand, loadPplSetup, invoke, Blue, Color_Off, BGreen, getFitsHeaders, getFitsHeader, setFitsHeaders, findInFile, subtractFitsBackground, assureFolder
from pmdisco import Discovery
from pmhotpix import BadPixelDetector, BadPixelEliminator


class Pipeline:

    opt = {}  # command line options
    pplSetup = {}  # PPL setup from ppl-setup config file
    badPixels = {}

    # Common arguments (saturation level, image section & trimming, etc.):
    COMMON_ARGS = "--saturation 16000 --trim"

    FLATLIB_FOLDER = expanduser("~/.pmlib/flat")

    TEMPDIR = "temp"

    BIAS_FOLDER = None
    DARK_FOLDER = None
    FLAT_FOLDER = None
    FLAT_BIAS_FOLDER = None
    FLAT_DARK_FOLDER = None
    LIGHT_FOLDERS = None

    REF = None

    def __init__(self, opt):
        self.opt = opt

    def discoverFolders(self):
        disco = Discovery(self.opt, self.pplSetup)
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
        t_min = None
        count = 0
        for sourceFile in sourceFiles:
            h = getFitsHeaders(sourceFile, ['DATE-OBS', 'EXPTIME'])
            t = time.mktime(time.strptime(h['DATE-OBS'], '%Y-%m-%dT%H:%M:%S'))
            e = int(h['EXPTIME'])
            e_sum = e_sum + e
            m_sum = m_sum + (t + e / 2) * e
            if t_min == None or t < t_min:
                t_min = t
            count += 1

        m_avg = m_sum / e_sum
        d = datetime(*time.gmtime(m_avg)[:6])
        t_obs = datetime(*time.gmtime(t_min)[:6])
        headers = { 
            'DATE-OBS' : t_obs.isoformat(), 
            'EXPTIME'  : e_sum, 
            'DATE-MID' : (d.isoformat(), 'Effective center of time of observation'),
            'NCOMBINE' : (count, 'number of images used for stacking'),
            'MCOMBINE' : ('sum', 'combination mode')
        }
        setFitsHeaders(targetFile, headers)
        print(f'observation start time: {t_obs.isoformat()}, observation center time: {d.isoformat()}, cumulated exptime: {str(e_sum)} sec')

    def invoked(self, cmd):
        if self.opt['debug']:
            printDebug(cmd)
        result = invoke(cmd)
        if result.startswith('ERROR:'):
            printError(result[len('ERROR: '):])
        elif self.opt['debug']:
            printDebug(result)
        return result

    def raw2fitsFile(self, rawfile, color):
        FITS_NAME = "%s-%s.fits" % (rawfile[:rawfile.rfind('.')], color)
        if self.opt['overwrite'] and exists(FITS_NAME):
            os.remove(FITS_NAME)

        if not exists(FITS_NAME):
            print("%s -> %s" % (rawfile, FITS_NAME))

            # convert raw image to fits
            # TODO: use invoke, parse options with ' or " correctly, or use preparsed options
            os.system("rawtran -c %s -o %s -B 16 -C '-4 -D -t 0' -X '-q 3 -w' %s" % (color, FITS_NAME, rawfile))

            # read image temperature from raw file
            exif = invoke("exiftool -s -g %s" % (rawfile))
            IMAGETEMP = None
            for line in exif.split('\n'):
                if line.startswith('CameraTemperature'):
                    IMAGETEMP = line.split(':')[1].strip()
                    break

            # add extra headers for fits file
            h = getFitsHeaders(FITS_NAME, ['INSTRUME', 'TELESCOP', 'DATE-OBS'])
            hx = {}

            if IMAGETEMP:
                # write image temperature to fits file
                hx['CCD-TEMP'] = (IMAGETEMP, 'CCD Temperature (Celsius)')

            if self.opt['imageTime'] == 'LT':
                # convert observation time from local time to UT
                IMAGE_DATE = h['DATE-OBS']
                IMAGE_DATE_UTC = self.lt2ut(IMAGE_DATE)
                hx['DATE-OBS'] = IMAGE_DATE_UTC
                hx['DATE-IMG'] = (IMAGE_DATE, 'Original image date (in local time)')

            if 'INSTRUME' not in h:
                # write camera name to fits file
                hx['INSTRUME'] = (self.pplSetup['DEF_CAMERA'] if 'DEF_CAMERA' in self.pplSetup else 'Generic Camera', 'Camera manufacturer and model')

            if 'TELESCOP' not in h:
                # write telescope name to fits file
                hx['TELESCOP'] = (self.pplSetup['DEF_TELESCOPE'] if 'DEF_TELESCOPE' in self.pplSetup else 'Generic Telescope', 'Telescope manufacturer and model')

            setFitsHeaders(FITS_NAME, hx)

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
        DARKLIST.sort()

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
        FLATLIST.sort()

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

    def getFitsHeadersForFlat(self, fitsFile):
        h = getFitsHeaders(fitsFile, ["INSTRUME", "TELESCOP", "DATE-OBS"])
        instrument = (h["INSTRUME"] if "INSTRUME" in h else self.pplSetup["DEF_CAMERA"]).translate({ord(c): '_' for c in " /."})
        telescope = (h["TELESCOP"] if "TELESCOP" in h else self.pplSetup["DEF_TELESCOPE"]).translate({ord(c): '_' for c in " /."})
        date = h["DATE-OBS"].split('T')[0].translate({ord(c): None for c in ":-"})
        return instrument, telescope, date

    def saveMasterFlat(self, flatFolder, color):
        flatFileName = "%s/%s-%s.fits" % (flatFolder, self.pplSetup['MASTER_FLAT_FILE'], color)
        instrument, telescope, date = self.getFitsHeadersForFlat(flatFileName)
        flatLibFolder = assureFolder(self.FLATLIB_FOLDER)
        libFlatFileName = "%s/%s-%s-%s-%s-%s.fits" % (flatLibFolder, self.pplSetup['MASTER_FLAT_FILE'], instrument, telescope, date, color)
        copyfile(flatFileName, libFlatFileName)
        print(f'Master flat {flatFileName} was saved into flat library as {libFlatFileName}')

    def findBestLibraryFlat (self, color, lightFileName):
        instrument, telescope, date = self.getFitsHeadersForFlat(lightFileName)
        masterFlatList = glob.glob("%s/%s-%s-%s-*-%s.fits" % (self.FLATLIB_FOLDER, self.pplSetup['MASTER_FLAT_FILE'], instrument, telescope, color))
        idate = int(date)
        date = None
        delta = 999
        for flat in masterFlatList:
            d = flat.split('/')[-1].split('-')[-2]
            dd = idate - int(d)
            if dd >= 0 and dd < delta:
                date = d
                delta = dd
        if date is not None:
            return "%s/%s-%s-%s-%s-%s.fits" % (self.FLATLIB_FOLDER, self.pplSetup['MASTER_FLAT_FILE'], instrument, telescope, date, color)
        else:
            return None

    def locateMasterFlat(self, color, lightFileName):
        if self.opt['useMasterFlat']:
            if self.opt['masterFlat'] == 'flatlib':
                # search in flat library (-M)
                libraryFlat = self.findBestLibraryFlat(color, lightFileName)
                if libraryFlat is not None:
                    return libraryFlat
                # fallback to local master flat
                localFlat = "%s/%s-%s.fits" % (self.FLAT_FOLDER, self.pplSetup['MASTER_FLAT_FILE'], color)
                return localFlat if exists(localFlat) else None

            else:
                # user given flat folder (-m)
                folderFlat = "%s/%s-%s.fits" % (self.opt['masterFlat'], self.pplSetup['MASTER_FLAT_FILE'], color)
                if exists(folderFlat):
                    return folderFlat
                # fallback to local flat
                localFlat = "%s/%s-%s.fits" % (self.FLAT_FOLDER, self.pplSetup['MASTER_FLAT_FILE'], color)
                if exists(localFlat):
                    return localFlat
                # fallback to library flat
                return self.findBestLibraryFlat(color, lightFileName)

        else:
            # local master flat
            localFlat = "%s/%s-%s.fits" % (self.FLAT_FOLDER, self.pplSetup['MASTER_FLAT_FILE'], color)
            if exists(localFlat):
                return localFlat
            # fallback to library flat
            return self.findBestLibraryFlat(color, lightFileName)

    def calibrate(self, lightFolder, calibFolder, color):
        # Names of the individual files storing the raw bias, dark, flat and object frames:
        LIGHT_PATTERN = "%s/%s*-%s.fits" % (lightFolder, self.pplSetup['LIGHT_FILE_PREFIX'], color)
        IOBJLIST = glob.glob(LIGHT_PATTERN)
        IOBJLIST.sort()
        if len(IOBJLIST) == 0:
            return 1

        # get calibration master images
        MB = "%s/%s-%s.fits" % (self.BIAS_FOLDER, self.pplSetup['MASTER_BIAS_FILE'], color)
        MD = "%s/%s-%s.fits" % (self.DARK_FOLDER, self.pplSetup['MASTER_DARK_FILE'], color)
#        MF = "%s/%s-%s.fits" % (self.FLAT_FOLDER, self.pplSetup['MASTER_FLAT_FILE'], color)
        MF = self.locateMasterFlat(color, IOBJLIST[0])
        print(f'Master flat: {MF}')

        # echo "$1/${LIGHT_FILE_PREFIX}*-$3.fits -> $2/${LIGHT_FILE_PREFIX}*-$3.fits"
        CALIB_PATTERN = "%s/%s*-%s.fits" % (calibFolder, self.pplSetup['LIGHT_FILE_PREFIX'], color)
        print("%s -> %s" % (LIGHT_PATTERN, CALIB_PATTERN))

        # Calibrated images: all the images have the same name but put into a separate directory ($TARGET):
        R_IOBJLIST = list(map(lambda x: calibFolder + '/' + basename(x), IOBJLIST))

        # The calibration of the object images:
        self.invoked("ficalib -i %s %s %s -o %s --input-master-bias %s --input-master-dark %s --input-master-flat %s" % (' '.join(IOBJLIST), self.COMMON_ARGS, self.imSizeArg(IOBJLIST[0]), ' '.join(R_IOBJLIST), MB, MD, MF))

        # post process calibrated images
        print("Subtract background for %s" % (CALIB_PATTERN))
        bpe = BadPixelEliminator()
        bpe.loadBadPixelsForDark(MD, color)

        for fitsFileName in R_IOBJLIST:

            # subtract background
            subtractFitsBackground(fitsFileName)

            # remove bad pixels
            bpe.process(fitsFileName)

        return 0

    def alignImages(self, imageList, refIndex):
        alignedImageList = []
        for n, image in enumerate(imageList):
            if n != refIndex:
                alignedImage, footprint = aa.register(image, imageList[refIndex], detection_sigma=12)
                alignedImageList.append(alignedImage)
            else:
                alignedImageList.append(image)
        return alignedImageList

    def combineImages(self, alignedImageList, refIndex):
        refImage = alignedImageList[refIndex]
        for n, image in enumerate(alignedImageList):
            if n != refIndex:
                refImage = refImage + image
        return self.subtractBackground(refImage)

    def subtractBackground(self, image):
        sigma_clip = SigmaClip(sigma = 3.0)
        bkg_estimator = MedianBackground()
        bkg = Background2D(image, (50, 50), filter_size = (3, 3), sigma_clip = sigma_clip, bkg_estimator = bkg_estimator)
        print("Bkg median: %7.4f, RMS median: %7.4f" % (bkg.background_median, bkg.background_rms_median))
        image -= bkg.background
        return image


    def registrate(self, calibFolder, seqFolder, color):

        # names of the images files to registrate
        CALIB_PATTERN = "%s/%s*-%s.fits" % (calibFolder, self.pplSetup['LIGHT_FILE_PREFIX'], color)
        COBJLIST = glob.glob(CALIB_PATTERN)
        COBJLIST.sort()
        if len(COBJLIST) == 0:
            return 1

        COMBINED_FILE = "%s/Combined-%s.fits" % (seqFolder, color)

        # load images
        fitsList = []
        for f in COBJLIST:
            hdul = fits.open(f)
            fitsList.append(hdul[0].data.byteswap().newbyteorder())

        # search for reference image by the best average FWHM
        refIndex = int(len(fitsList)/2) # TODO

	# align
        alignedImages = self.alignImages(fitsList, refIndex)

        # create a sequence of combined frames
        cc = self.opt['countCombine']
        if cc != 0:
            for a in range(0, len(alignedImages), cc):
                seqRefIndex = int(cc/2)
                combinedImage = self.combineImages(alignedImages[a:a+cc], seqRefIndex)

                fi = "%03d" % (a / cc)
                combined = "%s/%s%s-%s.fits" % (seqFolder, self.pplSetup['SEQ_FILE_PREFIX'], fi, color)
                print("combine " + combined)
                hdul = fits.open(COBJLIST[seqRefIndex], mode='update')
                hdul[0].data = combinedImage
                hdul.writeto(combined, overwrite=True)
                hdul.close()

                self.sumDateObs(combined, COBJLIST[a:a + cc])

        # create a combined frame of all images
        print("combine %s -> %s" % (CALIB_PATTERN, COMBINED_FILE))
        combinedImage = self.combineImages(alignedImages, refIndex)
        hdul = fits.open(COBJLIST[refIndex], mode='update')
        hdul[0].data = combinedImage

        # TODO add new FITS headers
        #   source FITS file names
        #   clipping area
        hdul.writeto(COMBINED_FILE, overwrite=True)
        hdul.close()

        self.sumDateObs(COMBINED_FILE, COBJLIST)

        # cleanup - remove temp files
        #  nothing to delete


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

        printInfo(title + ": Create master flat file(s).")
        for color in self.opt['color']:
            self.makeMasterFlat(flatFolder, biasFolder, darkFolder, color)
            if self.opt['saveFlat']:
                self.saveMasterFlat(flatFolder, color)

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

        printInfo(title + ": Register and stack calibrated light file(s).")

        # Create the sequence dir, if not exists
        if not exists(seqFolder):
            makedirs(seqFolder)

        self.REF = None

        for color in self.opt['color']:
            self.registrate(calibFolder, seqFolder, color)

        # TODO: cleanup - delete light FITS files

    def execute(self):

        ##########################
        # step 0. setup photometry
        ##########################
        self.pplSetup = loadPplSetup()

        self.discoverFolders()

        if self.opt['masterFlat'] is not None and self.opt['masterFlat'] != 'flatlib':
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
        # step 4. create master flat dark frame, if it is differs from master dark
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
        'color' : ['Gi'],          # photometry band, mandatory
        'countCombine' : 0,        # number of images to combine in a sequence
        'flatOnly' : False,        # create master flat only
        'saveFlat' : False,        # save master flat into flat library
        'useMasterFlat'  : False,  # 
        'imageTime' : 'LT',        #
        'masterFlat' : None,       # 'flatlib' or path for master flat
        'onError'   : 'noop',      # mg calculation method: comp, gcx, lfit
        'overwrite': False,        # force to overwrite existing results, optional
        'baseFolder': None,        # base folder, optional
        'calibFolder': None,       # optional folder for calibration frames (bias, dark, flat)
        'debug': False,            # debug mode
        }

    availableBands = ['gi', 'g', 'bi', 'b', 'ri', 'r', 'all']

    def __init__(self, argv):
        self.argv = argv
        pass

    def printTitle(self):
        print()
        print(BGreen + "ppl-calibration, version 1.2.0" + Color_Off)
        print(Blue + "Calibrate a set of RAW or FITS images." + Color_Off)
        print()

    def usage(self):
        print("Usage: ppl-calibration [OPTIONS]... [BASE_FOLDER]")
        print("Make calibration process for raw or fits images.")
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -c,  --color arg               set filter(s), arg is the color code, default color is 'Gi', for available color codes see below")
        print("  -n,  --count-combine n         set number of frames to combine in the sequence, 0 means all frames, default is 0")
        print("  -f,  --flat                    make master flat frame only")
        print("  -F,  --save-flat               save master flat into flat library")
        print("  -m,  --master-flat folder      use the given master-flat folder")
        print("  -M,  --use-flat                use master flat from flat library")
        print("  -t,  --image-time LT|UT        specify orignal image time zone, LT=local time, UT=universal time")
        print("       --calib-folder folder     alternative folder for calibration frames (bias, dark, flat)")
        print("  -w,  --overwrite               force to overwrite existing results")
        print("  -e,  --on-error noop|skip|stop specify what to do on error: noop=nothing to do; skip=remove the file on processing; stop=stop processing at all")
        print("       --debug                   print useful debugging informations ; reserve temp folder content")
        print("  -h,  --help                    print this page")
        print()
        print("Available filter color codes are:")
        print("  Gi | G | gi | g         green channel")
        print("  Bi | B | bi | b         blue channel")
        print("  Ri | R | ri | r         red channel")
        print("  all | ALL | All         all channels, results 3 separate frame")
        print()

    def processCommands(self):
        try:
            optlist, args = getopt (self.argv[1:], "c:n:fFm:Mt:we:h", ['color=', 'count-combine=', 'flat', 'save-flat', 'master-flat=', 'use-flat', 'image-time=', 'overwrite', 'on-error=', 'help', 'calib-folder=', 'alt-stack',  'debug'])
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
                    self.opt['color'] = ['Gi', 'Ri', 'Bi']
                else:
                    self.opt['color'] = [a]
            elif o == '-f' or o == '--flat':
                self.opt['flatOnly'] = True
            elif o == '-F' or o == '--save-flat':
                self.opt['saveFlat'] = True
            elif o == '-m' or o == '--master-flat':
                self.opt['useMasterFlat'] = True
                if self.opt['masterFlat'] != None:
                    printWarning ('Cannot use both flat folder and flat library in the same time; flat folder will be used.')
                self.opt['masterFlat'] = a
            elif o == '-M' or o == '--use-flat':
                self.opt['useMasterFlat'] = True
                if self.opt['masterFlat'] != None:
                    printWarning ('Cannot use both flat folder and flat library in the same time; flat folder will be used.')
                else:
                    self.opt['masterFlat'] = 'flatlib'
            elif o == '-n' or o == '--count-combine':
                self.opt['countCombine'] = int(a)
            elif o == '-t' or o == '--image-time':
                if a in ['LT', 'UT']:
                    self.opt['imageTime'] = a
                else:
                    printWarning("Bad image time zone value: %s, using default %s instead." % (a, self.opt['imageTime']))
            elif o == '-e' or o == '--on-error':
                if not a in ['noop', 'skip', 'stop']:
                    printWarning("Bad on-error instruction; available values are: noop, stop, skip.")
                else:
                    self.opt['onError'] = a
            elif o == '--calib-folder':
                if not isdir(a):
                    printError("%s is not a folder. Calibration folder option is ignored." % (a))
                else:
                    self.opt['calibFolder'] = a
            elif o == '-w' or o == '--overwrite':
                self.opt['overwrite'] = True
            elif o == '--debug':
                self.opt['debug'] = True
            elif o == '-h' or o == '--help':
                self.usage()
                exit(0)


        # checking invalid options
        if self.opt['flatOnly'] and self.opt['useMasterFlat']:
            self.opt['useMasterFlat'] = False
            self.opt['masterFlat'] = None
            printWarning('Options -m or -M are useless when -f option is set.')
        if self.opt['saveFlat'] and self.opt['useMasterFlat']:
            self.opt['saveFlat'] = False
            printWarning('Option -F is inconsistent with options -m or -M; master flat will not be saved into flat library.')

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
