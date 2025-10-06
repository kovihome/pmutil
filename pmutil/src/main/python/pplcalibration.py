#!/usr/bin/env python3
#
# PmUtils/pplcalibration
#
"""
Created on Mar 1, 2020

@author: kovi
"""

import glob
import os
import os.path
import shutil
import time
from datetime import datetime
from getopt import getopt, GetoptError
from shutil import copyfile
from sys import argv

import astroalign as aa
from astropy.io import fits
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground

import pmbase as pm
import pmconventions as pmc
from pmhotpix import BadPixelEliminator
from pmraw import RawConverter
from pmfits import FitsSeparator

# FITS header names
FITS_HEADER_NAXIS = "NAXIS"
FITS_HEADER_DATE_OBS = "DATE-OBS"
FITS_HEADER_EXPTIME = "EXPTIME"
FITS_HEADER_STACKCNT = "STACKCNT"
FITS_HEADER_BAYER = "BAYERPAT"
FITS_HEADER_OBSERVER = "OBSERVER"
FITS_HEADER_INSTRUMENT = "INSTRUME"
FITS_HEADER_TELESCOPE = "TELESCOP"
# 'IMAGETYP'
# 'NAXIS1'
# 'NAXIS2'
# 'CREATOR'
# 'DATE-MID'
# 'NCOMBINE'
# 'MCOMBINE'
# 'CCD-TEMP'


class Pipeline:
    opt = {}  # command line options

    badPixels = {}

    # Common arguments (saturation level, image section & trimming, etc.):
    COMMON_ARGS = "--saturation 16000 --trim"

    FLATLIB_FOLDER = os.path.expanduser("~/.pmlib/flat")

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
        disco = pmc.Discovery(self.opt, pm.setup)
        disco.discover()
        self.BIAS_FOLDER = disco.BIAS_FOLDER
        self.DARK_FOLDER = disco.DARK_FOLDER
        self.FLAT_FOLDER = disco.FLAT_FOLDER
        self.FLAT_BIAS_FOLDER = disco.FLAT_BIAS_FOLDER
        self.FLAT_DARK_FOLDER = disco.FLAT_DARK_FOLDER
        self.LIGHT_FOLDERS = disco.LIGHT_FOLDERS

    @staticmethod
    def imSizeArg(fitsFileName):
        h = pm.getFitsHeaders(fitsFileName, ['NAXIS1', 'NAXIS2'])
        return f"--image 0:0:{int(h['NAXIS1']) - 1:d}:{int(h['NAXIS2']) - 1:d}"

    @staticmethod
    def sumDateObs(targetFile, sourceFiles):
        e_sum = 0
        m_sum = 0
        t_min = None
        count = 0
        for sourceFile in sourceFiles:
            h = pm.getFitsHeaders(sourceFile, [FITS_HEADER_DATE_OBS, FITS_HEADER_EXPTIME])
            if h is not None:
                fmt = '%Y-%m-%dT%H:%M:%S.%f' if h[FITS_HEADER_DATE_OBS].find(".") != -1 else '%Y-%m-%dT%H:%M:%S'
                t = time.mktime(time.strptime(h[FITS_HEADER_DATE_OBS], fmt))
                e = int(h[FITS_HEADER_EXPTIME])
                e_sum = e_sum + e
                m_sum = m_sum + (t + e / 2) * e
                if t_min is None or t < t_min:
                    t_min = t
                count += 1

        m_avg = m_sum / e_sum
        d = datetime(*time.gmtime(m_avg)[:6])
        t_obs = datetime(*time.gmtime(t_min)[:6])
        headers = {
            FITS_HEADER_DATE_OBS: t_obs.isoformat(),
            FITS_HEADER_EXPTIME: e_sum,
            'DATE-MID': (d.isoformat(), 'Effective center of time of observation'),
            'NCOMBINE': (count, 'number of images used for stacking'),
            FITS_HEADER_STACKCNT: (count, 'number of images used for stacking'),
            'MCOMBINE': ('sum', 'combination mode')
        }
        pm.setFitsHeaders(targetFile, headers)
        print(
                f'observation start time: {t_obs.isoformat()}, observation center time: {d.isoformat()}, cumulated '
                f'exptime: {str(e_sum)} sec')

    def invoked(self, cmd):
        if self.opt['debug']:
            pm.printDebug(cmd)
        result = pm.invoke(cmd)
        if result.startswith('ERROR:'):
            pm.printError(result[len('ERROR: '):])
        elif self.opt['debug']:
            pm.printDebug(result)
        return result

    @staticmethod
    def def_fits_headers(imtype):
        hx = {'CREATOR': f"pmutil {pmc.PMUTIL_VERSION_SHORT}", 'IMAGETYP': imtype}
        if 'DEF_NAMECODE' in pm.setup:
            hx[FITS_HEADER_OBSERVER] = pm.setup['DEF_NAMECODE']
        if 'DEF_CAMERA' in pm.setup and pm.setup['DEF_CAMERA'] != 'Generic Camera':
            hx[FITS_HEADER_INSTRUMENT] = pm.setup['DEF_CAMERA']
        if 'DEF_TELESCOPE' in pm.setup and pm.setup['DEF_TELESCOPE'] != 'Generic Telescope':
            hx[FITS_HEADER_TELESCOPE] = pm.setup['DEF_TELESCOPE']
        return hx

    def raw2fits(self, folder:str) -> int:
        """
        Convert all RAW images in a folder into FITS format.
        FITS images is creating in the same folder.
        Args:
            folder (str): The folder name containing RAW images.
        Returns:
            int: number of RAW file converted, 0 if there are no RAW files in the folder at all.
        """
        # list of raw files in the input directory
        rawFiles = []
        for ext in pmc.RAW_FILE_EXTENSIONS:
            rawFiles.extend(glob.glob(folder + f"/*.{ext}"))
        rawFiles.sort()

        if len(rawFiles) == 0:
            return 0

        # determine image type from filename
        fp = os.path.basename(rawFiles[0])
        fn: str = os.path.splitext(fp)[0]
        imtype = fn.split('_')[0].upper()
        if imtype not in ["LIGHT", "DARK", "BIAS", "FLAT"]:
            imtype = "LIGHT"

        hx = self.def_fits_headers(imtype)

        # for all raw file do the conversion
        raw_converter = RawConverter(self.opt, hx)
        for n, rawfile in enumerate(rawFiles):
            print(f"Convert {n+1} of {len(rawFiles)}\r", end="")
            # for color in self.opt['color']:
            #     self.raw2fitsFile(rawfile, color)
            raw_converter.convert(rawfile)

        return len(rawFiles)

    def makeMasterBias(self, biasFolder:str, color:str) -> int:
        BIAS_PATTERN = f"{biasFolder}/{pm.setup['BIAS_FILE_PREFIX']}*-{color}.fits"
        BIASLIST = glob.glob(BIAS_PATTERN)
        BIASLIST.sort()

        masterBiasFile = f"{biasFolder}/{pm.setup['MASTER_BIAS_FILE']}-{color}.fits"

        if len(BIASLIST) == 0:
            return 1

        print(f"{BIAS_PATTERN} -> {masterBiasFile}")

        # Calibrated images: all the images have the same name but put into a separate directory ($TARGET):
        R_BIASLIST = list(map(lambda x: self.TEMPDIR + '/' + os.path.basename(x), BIASLIST))

        # The calibration of the individual bias frames, followed by their combination into a single master image:
        pm.invoke(
                f"ficalib -i {' '.join(BIASLIST)} {self.COMMON_ARGS} {self.imSizeArg(BIASLIST[0])} -o {' '.join(R_BIASLIST)}")
        pm.invoke(f"ficombine {' '.join(R_BIASLIST)} --mode median -o {masterBiasFile}")

        # cleanup: remove temp files
        for f in R_BIASLIST:
            os.remove(f)

        return 0

    def makeMasterDark(self, darkFolder, biasFolder, color):

        # Names of the individual files storing the raw bias, dark, flat and object frames:
        DARK_PATTERN = f"{darkFolder}/{pm.setup['DARK_FILE_PREFIX']}*-{color}.fits"
        DARKLIST = glob.glob(DARK_PATTERN)
        DARKLIST.sort()

        masterDarkFile = f"{darkFolder}/{pm.setup['MASTER_DARK_FILE']}-{color}.fits"
        masterBiasFile = f"{biasFolder}/{pm.setup['MASTER_BIAS_FILE']}-{color}.fits"

        if len(DARKLIST) == 0:
            return 1

        print(f"{DARK_PATTERN} -> {masterDarkFile}")

        # Calibrated images: all the images have the same name but put into a separate directory ($TARGET):
        R_DARKLIST = list(map(lambda x: self.TEMPDIR + '/' + os.path.basename(x), DARKLIST))

        # The calibration of the individual bias frames, followed by their combination into a single master image:
        pm.invoke(
                f"ficalib -i {' '.join(DARKLIST)} {self.COMMON_ARGS} {self.imSizeArg(DARKLIST[0])} -o {' '.join(R_DARKLIST)} --input-master-bias {masterBiasFile}")
        pm.invoke(f"ficombine {' '.join(R_DARKLIST)} --mode median -o {masterDarkFile}")

        # Calculate average ccd temperature from .cr2 files, and set it into the master dark
        tsum = 0
        count = 0
        for f in DARKLIST:
            hdr = pm.getFitsHeader(f, 'CCD-TEMP')
            if hdr:
                if ' ' in hdr:
                    hdr = hdr.split()[0]
                ccdtemp = int(hdr)
                tsum += ccdtemp
                count += 1
        if count > 0:
            AVGTEMP = int((tsum + (count / 2)) / count)
            print(f"average dark temperature: {AVGTEMP:d} C")
            pm.setFitsHeaders(masterDarkFile, {'CCD-TEMP': (f"{AVGTEMP:d}.", "CCD Temperature (Celsius)")})

        # cleanup: remove temp files
        for f in R_DARKLIST:
            os.remove(f)

        return 0

    def makeMasterFlat(self, flatFolder, biasFolder, darkFolder, color):
        # Names of the individual files storing the raw bias, dark, flat and object frames:
        FLAT_PATTERN = f"{flatFolder}/{pm.setup['FLAT_FILE_PREFIX']}*-{color}.fits"
        FLATLIST = glob.glob(FLAT_PATTERN)
        FLATLIST.sort()

        masterFlatFile = f"{flatFolder}/{pm.setup['MASTER_FLAT_FILE']}-{color}.fits"
        masterDarkFile = f"{darkFolder}/{pm.setup['MASTER_DARK_FILE']}-{color}.fits"
        masterBiasFile = f"{biasFolder}/{pm.setup['MASTER_BIAS_FILE']}-{color}.fits"

        if len(FLATLIST) == 0:
            return 1

        print(f"{FLAT_PATTERN} -> {masterFlatFile}")

        # Calibrated images: all the images have the same name but put into a separate directory ($TARGET):
        R_FLATLIST = list(map(lambda x: self.TEMPDIR + '/' + os.path.basename(x), FLATLIST))

        # The calibration of the individual flat frames, followed by their combination into a single master image:
        pm.invoke(
                f"ficalib -i {' '.join(FLATLIST)} {self.COMMON_ARGS} {self.imSizeArg(FLATLIST[0])} --post-scale 20000 -o {' '.join(R_FLATLIST)} --input-master-bias {masterBiasFile} --input-master-dark {masterDarkFile}")
        pm.invoke(f"ficombine {' '.join(R_FLATLIST)} --mode median -o {masterFlatFile}")

        # remove temp files
        for f in R_FLATLIST:
            os.remove(f)

        return 0

    @staticmethod
    def getFitsHeadersForFlat(fitsFile):
        h = pm.getFitsHeaders(fitsFile, [FITS_HEADER_INSTRUMENT, FITS_HEADER_TELESCOPE, FITS_HEADER_DATE_OBS])
        instrument = (h[FITS_HEADER_INSTRUMENT] if FITS_HEADER_INSTRUMENT in h else pm.setup["DEF_CAMERA"]).translate(
                {ord(c): '_' for c in " /."})
        telescope = (h[FITS_HEADER_TELESCOPE] if FITS_HEADER_TELESCOPE in h else pm.setup["DEF_TELESCOPE"]).translate(
                {ord(c): '_' for c in " /."})
        date = h[FITS_HEADER_DATE_OBS].split('T')[0].translate({ord(c): None for c in ":-"})
        return instrument, telescope, date

    def saveMasterFlat(self, flatFolder, color):
        flatFileName = f"{flatFolder}/{pm.setup['MASTER_FLAT_FILE']}-{color}.fits"
        instrument, telescope, date = self.getFitsHeadersForFlat(flatFileName)
        flatLibFolder = pm.assureFolder(self.FLATLIB_FOLDER)
        libFlatFileName = f"{flatLibFolder}/{pm.setup['MASTER_FLAT_FILE']}-{instrument}-{telescope}-{date}-{color}.fits"
        copyfile(flatFileName, libFlatFileName)
        pm.printInfo(f'Master flat {flatFileName} was saved into flat library as {libFlatFileName}')

    def findBestLibraryFlat(self, color, lightFileName):
        instrument, telescope, date = self.getFitsHeadersForFlat(lightFileName)
        masterFlatList = glob.glob(
                f"{self.FLATLIB_FOLDER}/{pm.setup['MASTER_FLAT_FILE']}-{instrument}-{telescope}-*-{color}.fits")
        idate = int(date)
        date = None
        delta = 999
        for flat in masterFlatList:
            d = flat.split('/')[-1].split('-')[-2]
            dd = idate - int(d)
            if 0 <= dd < delta:
                date = d
                delta = dd
        if date is not None:
            return f"{self.FLATLIB_FOLDER}/{pm.setup['MASTER_FLAT_FILE']}-{instrument}-{telescope}-{date}-{color}.fits"
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
                localFlat = f"{self.FLAT_FOLDER}/{pm.setup['MASTER_FLAT_FILE']}-{color}.fits"
                return localFlat if os.path.exists(localFlat) else None

            else:
                # user given flat folder (-m)
                folderFlat = f"{self.opt['masterFlat']}/{pm.setup['MASTER_FLAT_FILE']}-{color}.fits"
                if os.path.exists(folderFlat):
                    return folderFlat
                # fallback to local flat
                localFlat = f"{self.FLAT_FOLDER}/{pm.setup['MASTER_FLAT_FILE']}-{color}.fits"
                if os.path.exists(localFlat):
                    return localFlat
                # fallback to library flat
                return self.findBestLibraryFlat(color, lightFileName)

        else:
            # local master flat
            localFlat = f"{self.FLAT_FOLDER}/{pm.setup['MASTER_FLAT_FILE']}-{color}.fits"
            if os.path.exists(localFlat):
                return localFlat
            # fallback to library flat
            return self.findBestLibraryFlat(color, lightFileName)

    def calibrate(self, lightFolder, calibFolder, color):
        # Names of the individual files storing the raw bias, dark, flat and object frames:
        LIGHT_PATTERN = f"{lightFolder}/{pm.setup['LIGHT_FILE_PREFIX']}*-{color}.fits"
        IOBJLIST = glob.glob(LIGHT_PATTERN)
        IOBJLIST.sort()
        if len(IOBJLIST) == 0:
            return 1

        # get calibration master images
        MB = f"{self.BIAS_FOLDER}/{pm.setup['MASTER_BIAS_FILE']}-{color}.fits"
        MD = f"{self.DARK_FOLDER}/{pm.setup['MASTER_DARK_FILE']}-{color}.fits"
        MF = self.locateMasterFlat(color, IOBJLIST[0])
        print(f'Master flat: {MF}')

        # echo "$1/${LIGHT_FILE_PREFIX}*-$3.fits -> $2/${LIGHT_FILE_PREFIX}*-$3.fits"
        CALIB_PATTERN = f"{calibFolder}/{pm.setup['LIGHT_FILE_PREFIX']}*-{color}.fits"
        print(f"{LIGHT_PATTERN} -> {CALIB_PATTERN}")

        # Calibrated images: all the images have the same name but put into a separate directory ($TARGET):
        R_IOBJLIST = list(map(lambda x: calibFolder + '/' + os.path.basename(x), IOBJLIST))

        # The calibration of the object images:
        self.invoked(
                f"ficalib -i {' '.join(IOBJLIST)} {self.COMMON_ARGS} {self.imSizeArg(IOBJLIST[0])} -o {' '.join(R_IOBJLIST)} --input-master-bias {MB} --input-master-dark {MD} --input-master-flat {MF}")

        # post process calibrated images
        print(f"Subtract background for {CALIB_PATTERN}")
        h = pm.getFitsHeaders(R_IOBJLIST[0], ['NAXIS1', 'NAXIS2'])
        bp_max = [int(h['NAXIS1']), int(h['NAXIS2'])]
        bpe = BadPixelEliminator()
        bpe.loadBadPixelsForDark(MD, color, bp_max)

        for fitsFileName in R_IOBJLIST:
            # subtract background
            pm.subtractFitsBackground(fitsFileName)

            # remove bad pixels
            bpe.process(fitsFileName)

        return 0

    @staticmethod
    def alignImages(imageList, refIndex):
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

    @staticmethod
    def subtractBackground(image):
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        bkg = Background2D(image, (50, 50), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        pm.printDebug(f"Bkg median: {bkg.background_median:7.4f}, RMS median: {bkg.background_rms_median:7.4f}")
        image -= bkg.background
        return image

    def alignAndCombine(self, refImage, imageFileNameList):
        stackedImage = None
        for n, imageFileName in enumerate(imageFileNameList):
            #            print(f'Register image {imageFileNameList[n]}')
            if "_bad" in imageFileName:
                pm.printWarning(f"Image {imageFileName} is assigned as bad; exclude from stacking")
                continue
            try:
                image = self.loadImageForAlign(imageFileName)
                alignedImage, footprint = aa.register(image, refImage, detection_sigma=12)
                if stackedImage is None:
                    stackedImage = alignedImage
                else:
                    stackedImage = stackedImage + alignedImage
            except Exception as e:
                pm.printWarning(f"Registration failed on image {imageFileName}; exclude from stacking")
                print(f"Exception type={type(e)}, {str(e)}")
                # assign as bad image
                badImageFileName = imageFileName.replace(".fits", "_bad.fits")
                os.rename(imageFileName, badImageFileName)

        return self.subtractBackground(stackedImage)

    @staticmethod
    def loadImageForAlign(fitsName):
        hdul = fits.open(fitsName)
        return hdul[0].data.byteswap().newbyteorder()

    def registrate(self, calibFolder, seqFolder, color):

        # names of the images files to registrate
        CALIB_PATTERN = f"{calibFolder}/{pm.setup['LIGHT_FILE_PREFIX']}*-{color}.fits"
        COBJLIST = glob.glob(CALIB_PATTERN)
        COBJLIST.sort()
        if len(COBJLIST) == 0:
            return 1

        COMBINED_FILE = f"{seqFolder}/Combined-{color}.fits"

        # select reference image
        # TODO: search for reference image by the best average FWHM
        refIndex = int(len(COBJLIST) / 2)
        refImage = self.loadImageForAlign(COBJLIST[refIndex])  # TODO: select ref image more sophisticated

        # create a sequence of combined frames
        cc = self.opt['countCombine']
        if cc != 0:
            for a in range(0, len(COBJLIST), cc):
                print(f'Stack sequence {a} to {a + cc}')
                combinedImage = self.alignAndCombine(refImage, COBJLIST[a:a + cc])

                fi = f"{int(a/cc):03d}"
                combined = f"{seqFolder}/{pm.setup['SEQ_FILE_PREFIX']}{fi}-{color}.fits"
                print("Combine " + combined)
                hdul = fits.open(COBJLIST[refIndex], mode='update')
                hdul[0].data = combinedImage
                hdul.writeto(combined, overwrite=True)
                hdul.close()

                self.sumDateObs(combined, COBJLIST[a:a + cc])

        # search for reference image by the best average FWHM
        #        refImage = fitsList[int(len(fitsList)/2)] # TODO

        # create a combined frame of all images
        print(f"Combine {CALIB_PATTERN} -> {COMBINED_FILE}")
        combinedImage = self.alignAndCombine(refImage, COBJLIST)

        # TODO add new FITS headers
        #   source FITS file names
        #   clipping area
        hdul = fits.open(COBJLIST[refIndex], mode='update')
        hdul[0].data = combinedImage
        hdul.writeto(COMBINED_FILE, overwrite=True)
        hdul.close()

        self.sumDateObs(COMBINED_FILE, COBJLIST)

        # cleanup - remove temp files
        #  nothing to delete

        return 0

    def mastersExist(self, folder, prefix):
        for color in self.opt['color']:
            MB = f"{folder}/{prefix}-{color}.fits"
            if not os.path.exists(MB):
                return False
        return True

    def processBias(self, biasFolder, title):
        ex = self.mastersExist(biasFolder, pm.setup['MASTER_BIAS_FILE'])
        if self.opt['overwrite'] or not ex:
            pm.printInfo(title + ": Convert bias RAW files to FITS.")
            ident = self.preprocessRawFrames(biasFolder)
            raw_count = ident["count"]
            if raw_count == 0 and not ex:
                pm.printError(f"No Bias RAW images was found in folder {biasFolder}")
            elif raw_count == 0 and self.opt['overwrite'] and ex:
                pm.printWarning(f"No Bias RAW images was found in folder {biasFolder}, but Master Bias images exist.")
            else:
                pm.printInfo(title + ": Create master bias file.")
                for color in self.opt['color']:
                    self.makeMasterBias(biasFolder, color)
        else:
            pm.printInfo(title + ": Master bias file(s) are already created.")
        # TODO: cleanup - delete bias FITS files

    def processDark(self, darkFolder, biasFolder, title):
        ex = self.mastersExist(darkFolder, pm.setup['MASTER_DARK_FILE'])
        if self.opt['overwrite'] or not ex:
            pm.printInfo(title + ": Convert dark RAW files to FITS.")
            ident = self.preprocessRawFrames(darkFolder)
            raw_count = ident["count"]
            if raw_count == 0 and not ex:
                pm.printError(f"No Dark RAW images was found in folder {darkFolder}")
            elif raw_count == 0 and self.opt['overwrite'] and ex:
                pm.printWarning(f"No Dark RAW images was found in folder {darkFolder}, but Master Dark images exist.")
            else:
                pm.printInfo(title + ": Create master dark file(s).")
                for color in self.opt['color']:
                    self.makeMasterDark(darkFolder, biasFolder, color)
        else:
            pm.printInfo(title + ": Master Dark file(s) are already created.")
        # TODO: cleanup - delete dark FITS files

    def processFlat(self, flatFolder, biasFolder, darkFolder, title):
        pm.printInfo(title + ": Convert flat RAW files to FITS.")
        ident = self.preprocessRawFrames(flatFolder)
        raw_count = ident["count"]
        if raw_count == 0:
            pm.printError(f"No Flat RAW images was found in folder {darkFolder}")
        else:
            pm.printInfo(title + ": Create master flat file(s).")
            for color in self.opt['color']:
                self.makeMasterFlat(flatFolder, biasFolder, darkFolder, color)
                if self.opt['saveFlat']:
                    self.saveMasterFlat(flatFolder, color)

    def processCalibration(self, lightFolder, calibFolder, title):
        pm.printInfo(title + ": Convert light RAW files to FITS.")
        self.preprocessRawFrames(lightFolder)

        pm.printInfo(title + ": Create calibrated light file(s).")

        # Create dir for the calibrated images, if not exists
        if not os.path.exists(calibFolder):
            os.makedirs(calibFolder)

        for color in self.opt['color']:
            self.calibrate(lightFolder, calibFolder, color)
        # TODO: cleanup - delete light FITS files

    def processRegistration(self, calibFolder, seqFolder, title):

        pm.printInfo(title + ": Register and stack calibrated light file(s).")

        # Create the sequence dir, if not exists
        if not os.path.exists(seqFolder):
            os.makedirs(seqFolder)

        self.REF = None

        for color in self.opt['color']:
            self.registrate(calibFolder, seqFolder, color)

        # TODO: cleanup - delete light FITS files

    def preprocessRawFrames(self, folder:str) -> dict:
        pm.printDebug(f"Preprocess RAW Frames in folder {folder}")
        ident = self.identifyRawFrames(folder)
        pm.printDebug(f"Identify RA Frames as {ident}")
        if ident["type"] == "raw":
            pm.printDebug(f"Convert RAW to FITS")
            file_count = self.raw2fits(folder)
            ident["count"] = file_count

        elif ident["type"] == "fits":
            pm.printDebug(f"Separate FITS to color channels")
            fs = FitsSeparator(self.opt)
            fitsFiles = glob.glob(f"{folder}/*.{ident["ext"]}")
            for fitsFile in fitsFiles:
                fs.separate(fitsFile, folder)
            ident["count"] = len(fitsFiles)

        else:
            ident["count"] = 0

        return ident

    def identifyRawFrames(self, folder:str) -> dict:
        identification = {
            "type": "",
            "ext": "",
            "pre": False,
            "flat": False,
            "stack": False
        }
        # identify as RAW images
        for ext in pmc.RAW_FILE_EXTENSIONS:
            rawFiles = glob.glob(f"{folder}/*.{ext}")
            if len(rawFiles) > 0:
                identification["type"] = "raw"
                identification["ext"] = ext
                return identification

        # identify as FITS images
        fitsFiles = []
        for ext in pmc.FITS_FILE_EXTENSIONS:
            fitsFiles = glob.glob(f"{folder}/*.{ext}")
            if len(fitsFiles) > 0:
                identification["type"] = "fits"
                identification["ext"] = ext
                break

        if identification["type"] == "fits":

            # get FITS headers
            headers = pm.getFitsHeaders(fitsFiles[0], [FITS_HEADER_STACKCNT])

            # check flat existence
            if (self.FLAT_FOLDER is None or self.FLAT_FOLDER == "") and not self.opt["useMasterFlat"] and self.opt[
                "masterFlat"] is None:
                identification["flat"] = True

            # check precalibration
            if FITS_HEADER_STACKCNT in headers:
                identification["stack"] = True
                identification["pre"] = True  # not 100% sure, but for Seestar is OK

        return identification

    def execute(self):

        ##########################
        # step 0. setup photometry
        ##########################

        self.discoverFolders()

        ident = self.identifyRawFrames(self.LIGHT_FOLDERS[0])
        print(f"Light identification: {ident}")

        if not ident["flat"] and self.opt['masterFlat'] is not None and self.opt['masterFlat'] != 'flatlib':
            if not os.path.isdir(self.opt['masterFlat']):
                pm.printError("Master-flat folder %s not exists." % (self.opt['masterFlat']))
                exit(1)
            else:
                hasMaster = True
                for c in self.opt['color']:
                    mfc = "%s/%s-%s.fits" % (self.opt['masterFlat'], pm.setup['MASTER_FLAT_FILE'], c)
                    print(mfc)
                    if not os.path.exists(mfc):
                        pm.printError("No %s master-flat file in the directory %s" % (c, self.opt['masterFlat']))
                        hasMaster = False
                if hasMaster:
                    self.FLAT_FOLDER = self.opt['masterFlat']
                    print("master-flat: " + self.opt['masterFlat'])
                else:
                    exit(1)

        # create the temp dir, if not exists
        if not os.path.exists(self.TEMPDIR):
            os.makedirs(self.TEMPDIR)

        if not self.opt['flatOnly'] and not ident["pre"]:
            ####################################
            # step 1. create master bias frame
            ####################################

            self.processBias(self.BIAS_FOLDER, "BIAS")

            ####################################
            # step 2. create master dark frame
            ####################################

            self.processDark(self.DARK_FOLDER, self.BIAS_FOLDER, "DARK")

        if not self.opt['useMasterFlat'] and not ident["flat"]:

            # process flat bias, flat dark and flat, only if flat master is not exist
            ex = self.mastersExist(self.FLAT_FOLDER, pm.setup['MASTER_FLAT_FILE'])
            if self.opt['overwrite'] or not ex:

                ############################################################################
                # step 3. create master flat bias frame, if it differs from master bias
                ############################################################################
                if self.opt['flatOnly'] or self.FLAT_BIAS_FOLDER != self.BIAS_FOLDER:
                    self.processBias(self.FLAT_BIAS_FOLDER, "FLAT BIAS")

                ############################################################################
                # step 4. create master flat dark frame, if it differs from master dark
                ############################################################################
                if self.opt['flatOnly'] or self.FLAT_DARK_FOLDER != self.DARK_FOLDER:
                    self.processDark(self.FLAT_DARK_FOLDER, self.FLAT_BIAS_FOLDER, "FLAT DARK")

                ##############################
                # step 5. create master flat
                ##############################
                self.processFlat(self.FLAT_FOLDER, self.FLAT_BIAS_FOLDER, self.FLAT_DARK_FOLDER, "FLAT")

            else:
                pm.printInfo("FLAT: Master flat file(s) are already created.")

        if not self.opt['flatOnly']:

            for lf in self.LIGHT_FOLDERS:
                ##################################
                # step 6. calibrate light frames
                ##################################
                cf = lf.replace(pm.setup['LIGHT_FOLDER_NAME'], pm.setup['CALIB_FOLDER_NAME'])
                if not ident["pre"]:
                    self.processCalibration(lf, cf, "CALIBRATE")
                else:
                    shutil.copy(lf, cf)

                ###############################################
                # step 7. registrate and combine light frames
                ###############################################
                if not ident["stack"]:
                    sf = lf.replace(pm.setup['LIGHT_FOLDER_NAME'], pm.setup['SEQ_FOLDER_NAME'])
                    self.processRegistration(cf, sf, "REGISTRATE")

        print()
        print(pm.Blue + "Calibration is ready." + pm.Color_Off)


class MainApp:
    opt = {
        'color': ['Gi'],  # photometry band, mandatory
        'countCombine': 0,  # number of images to combine in a sequence
        'flatOnly': False,  # create master flat only
        'saveFlat': False,  # save master flat into flat library
        'useMasterFlat': False,  #
        'imageTime': 'LT',  #
        'masterFlat': None,  # 'flatlib' or path for master flat
        'onError': 'noop',  # action on error: noop=nothing to do; skip=remove the file on processing; stop=stop processing at all
        'overwrite': False,  # force to overwrite existing results, optional
        'baseFolder': None,  # base folder, optional
        'calibFolder': None,  # optional folder for calibration frames (bias, dark, flat)
        'debug': False,  # debug mode
    }

    availableBands = ['gi', 'g', 'bi', 'b', 'ri', 'r', 'all']

    def __init__(self, argv):
        self.argv = argv
        pass

    @staticmethod
    def printTitle():
        print()
        print(f"{pm.BGreen}ppl-calibration, version {pmc.PMUTIL_VERSION}{pm.Color_Off}" )
        print(pm.Blue + "Calibrate a set of RAW or FITS images." + pm.Color_Off)
        print()

    @staticmethod
    def usage():
        print("Usage: ppl-calibration [OPTIONS]... [BASE_FOLDER]")
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print(
                "  -c,  --color arg               set filter(s), arg is the color code, default color is 'Gi', for available color codes see below")
        print(
                "  -n,  --count-combine n         set number of frames to combine in the sequence, 0 means all frames, default is 0")
        print("  -f,  --flat                    make master flat frame only")
        print("  -F,  --save-flat               save master flat into flat library")
        print("  -m,  --master-flat folder      use the given master-flat folder")
        print("  -M,  --use-flat                use master flat from flat library")
        print("  -t,  --image-time LT|UT        specify original image time zone, LT=local time, UT=universal time")
        print("       --calib-folder folder     alternative folder for calibration frames (bias, dark, flat)")
        print("  -w,  --overwrite               force to overwrite existing results")
        print(
                "  -e,  --on-error noop|skip|stop specify what to do on error: noop=nothing to do; skip=remove the file on processing; stop=stop processing at all")
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
            optlist, args = getopt(self.argv[1:], "c:n:fFm:Mt:we:h",
                                   ['color=', 'count-combine=', 'flat', 'save-flat', 'master-flat=', 'use-flat',
                                    'image-time=', 'overwrite', 'on-error=', 'help', 'calib-folder=', 'alt-stack',
                                    'debug'])
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
                    self.opt['color'] = ['Gi', 'Ri', 'Bi']
                else:
                    self.opt['color'] = [a]
            elif o == '-f' or o == '--flat':
                self.opt['flatOnly'] = True
            elif o == '-F' or o == '--save-flat':
                self.opt['saveFlat'] = True
            elif o == '-m' or o == '--master-flat':
                self.opt['useMasterFlat'] = True
                if self.opt['masterFlat'] is not None:
                    pm.printWarning(
                            'Cannot use both flat folder and flat library in the same time; flat folder will be used.')
                self.opt['masterFlat'] = a
            elif o == '-M' or o == '--use-flat':
                self.opt['useMasterFlat'] = True
                if self.opt['masterFlat'] is not None:
                    pm.printWarning(
                            'Cannot use both flat folder and flat library in the same time; flat folder will be used.')
                else:
                    self.opt['masterFlat'] = 'flatlib'
            elif o == '-n' or o == '--count-combine':
                self.opt['countCombine'] = int(a)
            elif o == '-t' or o == '--image-time':
                if a in ['LT', 'UT']:
                    self.opt['imageTime'] = a
                else:
                    pm.printWarning(
                            "Bad image time zone value: %s, using default %s instead." % (a, self.opt['imageTime']))
            elif o == '-e' or o == '--on-error':
                if a not in ['noop', 'skip', 'stop']:
                    pm.printWarning("Bad on-error instruction; available values are: noop, stop, skip.")
                else:
                    self.opt['onError'] = a
            elif o == '--calib-folder':
                if not os.path.isdir(a):
                    pm.printError(f"{a} is not a folder. Calibration folder option is ignored.")
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
            pm.printWarning('Options -m or -M are useless when -f option is set.')
        if self.opt['saveFlat'] and self.opt['useMasterFlat']:
            self.opt['saveFlat'] = False
            pm.printWarning(
                    'Option -F is inconsistent with options -m or -M; master flat will not be saved into flat library.')

        if len(args) > 0:
            self.opt['baseFolder'] = args[0]
            if args[0].endswith('/'):
                self.opt['baseFolder'] = args[0][:-1]

            if not os.path.isdir(self.opt['baseFolder']):
                pm.printError("Base folder %s not exists or not a directory." % (self.opt['baseFolder']))
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
        pm.saveCommand(self.opt['baseFolder'], self.argv, 'calibration')

        start = datetime.now()

        ppl = Pipeline(self.opt)
        ppl.execute()

        exectime = (datetime.now() - start).total_seconds()
        print("%sexecution time was %d seconds.%s" % (pm.Blue, exectime, pm.Color_Off))


if __name__ == '__main__':
    app = MainApp(argv)
    app.run()

# end __main__
