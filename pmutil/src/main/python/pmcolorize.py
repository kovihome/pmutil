#!/usr/bin/env python3
#
# PmUtils/pmcolorize
#
'''
Created on Feb 22, 2020

@author: kovi
'''
from sys import argv
from getopt import getopt, GetoptError
from glob import glob
from os import remove
from os.path import isdir, exists

import numpy as np
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.visualization import make_lupton_rgb, simple_norm
from matplotlib.image import imsave
import matplotlib.pyplot as plt
from photutils import Background2D, MedianBackground
from PIL import Image, ImageDraw

from pmbase import printError, printInfo, loadPplSetup, Blue, Color_Off, BGreen, invoke
import pmconventions as pmc
import img_scale


class Colorize:

    opt = {}  # command line options
    ppl = {}  # PPL setup from ppl-setup config file

    scaleMethods = ['linear', 'sqrt', 'log', 'asinh', 'hist']

    wcs = None

#    wb = [1.15, 1.0, 0.79]   # for Canon EOS 1100D
    wb = [1.25, 1.0, 0.96]  # for Canon EOS 1100D sRGB
#   wb = [1.12, 1.0, 0.77]   # for Canon EOS 350D
#   wb = [1.04, 1.0, 0.72]   # for Canon EOS 450D
#   wb = [0.98, 1.0, 0.82]   # for Canon EOS 1300D
#   wb = [1.01, 1.0, 0.87]   # for Canon EOS 5D Mark III
#   wb = [1.20, 1.0, 0.72]   # for Canon EOS 60D

    SKY_SIGMA = 3.0
    SKY_CONVERGENCE = 0.01
    SCALE_RANGE = 256.0
    SCALE_METHOD = 'linear'

    role2color = {
        'V' : 'crimson',
        'C' : 'gold',
        'F' : 'grey',
        'D' : 'lightskyblue',
        'P' : 'grey',
        'T' : 'forestgreen',
        'M' : 'sandybrown'
        }

    def __init__(self, opt):
        self.opt = opt
        if 'scaleMethod' in self.opt:
            self.SCALE_METHOD = self.opt['scaleMethod']

    def doSubtractBackground(self, rgb):
        sigma_clip = SigmaClip(sigma = 3.0)
        bkg_estimator = MedianBackground()
        for j in range(3):
            bkg = Background2D(rgb[j], (50, 50), filter_size = (3, 3), sigma_clip = sigma_clip, bkg_estimator = bkg_estimator)
#            bkg = Background2D(rgb[j], (100, 100), filter_size = (5, 5), sigma_clip = sigma_clip, bkg_estimator = bkg_estimator)
            #print("Bkg median: %7.4f, RMS median: %7.4f" % (bkg.background_median, bkg.background_rms_median))
            rgb[j] -= bkg.background
#            plt.imshow(bkg.background, origin='lower', cmap='Greys_r')
#            plt.show()

    def doHistScaling(self, img):
        # Contrast stretching
        p2, p98 = np.percentile(img, (2, 98))
        img_rescale = exposure.rescale_intensity(img, in_range=(p2, p98))
        return img_rescale

        # Equalization
        #return exposure.equalize_hist(img)

        # Adaptive Equalization
        #return exposure.equalize_adapthist(img, clip_limit=0.03)

    def doScaling(self, img, rgb):
        for j in range(3):
            smin, it = img_scale.sky_mean_sig_clip(rgb[j], self.SKY_SIGMA, self.SKY_CONVERGENCE)
#            pxmax = max(rgb[j].flatten())

            smax = int(self.SCALE_RANGE / (self.opt['scaling'] * self.wb[j]) + smin)
            smin = max(smin, 0)
#            smax = 512

            p2, p98 = np.percentile(rgb[j], (2, 98))

            print(f"Scaling sky mean min={smin}, max={smax}, percentile min={p2}, max={p98}")

#            hist, bin_edges = np.histogram(rgb[j], bins=256, range=(p2, p98))

            # configure and draw the histogram figure
#            plt.figure()
#            plt.title("Grayscale Histogram")
#            plt.xlabel("grayscale value")
#            plt.ylabel("pixel count")
#            plt.xlim([0.0, 1.0])  # <- named arguments do not work here

#            plt.plot(hist)  # <- or here
#            plt.show()

            if self.SCALE_METHOD == self.scaleMethods[0]:
                img[:, :, j] = img_scale.linear(rgb[j], scale_min = smin, scale_max = smax)
            elif self.SCALE_METHOD == self.scaleMethods[1]:
                img[:, :, j] = img_scale.sqrt(rgb[j], scale_min = smin, scale_max = smax)
            elif self.SCALE_METHOD == self.scaleMethods[2]:
                img[:, :, j] = img_scale.log(rgb[j], scale_min = smin, scale_max = smax)
            elif self.SCALE_METHOD == self.scaleMethods[3]:
                img[:, :, j] = img_scale.asinh(rgb[j], scale_min = smin, scale_max = smax)
            elif self.SCALE_METHOD == self.scaleMethods[4]:
                img[:, :, j] = self.doHistScaling(rgb[j])
           

    def doScalingLupton(self, rgb, imgFileName):
#        return make_lupton_rgb(rgb[0], rgb[1], rgb[2], stretch=(1/self.opt['scaling']), filename = imgFileName)
        return make_lupton_rgb(rgb[0], rgb[1], rgb[2], filename = imgFileName)

    def calcXY(self, rc, raName, decName):
        bf = self.opt['baseFolder'].rstrip('/')

        # 1. convert refCat to fits format
        refFitsFile = bf + '/ref.cat.fits'
        rc.write(refFitsFile, overwrite=True)

        # 2. calculate ref objects' frame xy points
        pmCatFileBase = bf + '/' + self.ppl['PHOT_FOLDER_NAME'] + '/Combined-Gi'
        wcsFile = pmCatFileBase + '.wcs'      
        axyFile = pmCatFileBase + '.ref.axy'
        if not exists(wcsFile):
            return None

        invoke(f"wcs-rd2xy -w {wcsFile} -i {refFitsFile} -o {axyFile} -R {raName} -D {decName}")  # -f option need argument

        remove(refFitsFile)

        # 3. merge frame xy point to refCat
        tlen = len(rc)
        rc['X'] = [0.0] * tlen
        rc['Y'] = [0.0] * tlen

        f = fits.open(axyFile)
        d = f[1].data

        for j in range(tlen):
            x, y = d[j]
            rc[j]['X'] = x
            rc[j]['Y'] = y

        return rc

    def plotStars(self, img, imgFileName):
        # TODO: plot photometry result stars (p option)
        # TODO: plot transient result stars and moving objects (t and m options)

        annImgFileName = imgFileName.replace('.jpg', '-ann.jpg')
        print(f"Create annotated image {annImgFileName}")
        
        # load refcat file
        rc = pmc.loadRefcat(self.opt['baseFolder'])
        if rc is None:
            printError(f'No ref.cat in the folder {self.folder} ; use ppl-refcat first.')
            return

        # convert world coords to pixel coords
        rcxy = self.calcXY(rc, 'RA_DEG', 'DEC_DEG')
        if rcxy is None:
            printError(f'No .wcs file in the photometry folder ; use ppl-photometry first.')
            return

        # open original image
        image = Image.open(imgFileName)
        draw = ImageDraw.Draw(image)

        # for all objects:
        for j in range(len(rcxy)):
            c = rcxy[j]
            if c['ROLE'] in self.opt['plot']:
                x = c['X']
                y = c['Y']
                r = 15
                color = self.role2color[c['ROLE']]
                draw.ellipse((x-r, y-r, x+r, y+r), outline=color, width=1)

        # save annotated image
        image.convert(mode='RGB').save(annImgFileName)

    def saveImage(self, img, imgFileName):
        imsave(imgFileName, img)
        if self.opt['color'] == 'all':
            printInfo("Color image %s created." % (imgFileName))
        else:
            printInfo("Monochrome %s image %s created." % (self.opt['color'], imgFileName))

    def createImage(self, baseFolder, subFolder, ext, colors):
        # TODO: before addition, images should be reregistered with Gi as reference, for more exact alignment

        rgb = []

        imgFileName = "%s/%s.jpg" % (baseFolder, baseFolder.split('/')[-1])

        for c in colors:
            fitsFileName = "%s/%s/Combined-%s.%s" % (baseFolder, subFolder, c, ext)
            print("%s color source image: %s" % (c, fitsFileName))
            if not exists(fitsFileName):
                printError("image file %s not exists." % (fitsFileName))
                return
            hdu = fits.open(fitsFileName)
            rgb.append(hdu[0].data)

        img = np.zeros((rgb[1].shape[0], rgb[1].shape[1], 3), dtype = float)

        # TODO: with log method, save is failed: Floating point image RGB values must be in the 0..1 range.
        # need some normalization

        try:
            self.doSubtractBackground(rgb)
            if self.SCALE_METHOD == 'lupton':
                img = self.doScalingLupton(rgb, imgFileName)
            else:
                self.doScaling(img, rgb)
                plt.imshow(img)
                self.saveImage(img, imgFileName)
        except ValueError as e:
            printError("creating color image with %s scaling mathod is failed: %s" % (self.SCALE_METHOD, str(e)))
            if self.SCALE_METHOD != self.scaleMethods[0]:
                self.SCALE_METHOD = self.scaleMethods[0]
                try:
                    self.doScaling(img, rgb)
                    self.saveImage(img, imgFileName)
                except ValueError as ex:
                    printError("creating color image with alternative %s scaling mathod is failed: %s" % (self.SCALE_METHOD, str(ex)))
                    return

        if self.opt['plot'] != None:
            self.plotStars(img, imgFileName)

    def execute(self):
        print('Colorize.execute starts.')
        self.ppl = loadPplSetup()
        pmFolder = self.ppl['PHOT_FOLDER_NAME']
        seqFolder = self.ppl['SEQ_FOLDER_NAME']

        baseFolderPattern = self.opt['baseFolder']
        print(f'baseFolderPattern: {baseFolderPattern}')

        baseFolders = glob('*' + baseFolderPattern + '*')
        print(f'baseFolders: {baseFolders}')

        colors = ['Ri', 'Gi', 'Bi']
        if self.opt['color'] != 'all':
            colors = [self.opt['color']] * 3
        print(f'colors: {colors}')

        for folder in baseFolders:
            if isdir(folder + '/' + pmFolder):
                self.createImage(folder, pmFolder, 'ast.fits', colors)
            elif isdir(folder + '/' + seqFolder):
                self.createImage(folder, seqFolder, 'fits', colors)


class MainApp:

    # TODO: wb should be command line option somehow
    opt = {
        'scaleMethod': 'linear',
        'scaling': 1.0,
        'color': 'all',
        'baseFolder': None,
        'plot': None,
        }

    def __init__(self, argv):
        self.argv = argv
        pass

    def printTitle(self):
        print()
        print(BGreen + "ppl-colorize, version 1.2.0 " + Color_Off)
        print(Blue + "Make color jpeg image from calibrated FITS images." + Color_Off)
        print()

    def usage(self):
        print("Usage: ppl-colorize [OPTIONS]... [BASE_FOLDER]")
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -m,  --method arg   scaling method, available values are: linear, sqrt, log, asinh")
        print("  -c,  --color arg    use selected color for all channel ; it results monochrome image")
        print("       --scale value  scaling constant")
        print("       --plot arg     plot objects on image ; v - variables, c - comp stars, t - tranzients, f - field stars")
        print("  -h,  --help         print this page")
        print()

    def processCommands(self):
        try:
            optlist, args = getopt (self.argv[1:], "c:m:h", ['color=', 'method=', 'scale=', 'plot=', 'help'])
        except GetoptError:
            printError ('Invalid command line options')
            return

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-m' or o == '--method':
                if not a in Colorize.scaleMethods:
                    printError("Invalid scaling method value: %s" % (a))
                self.opt['scaleMethod'] = a
            elif o == '--scale':
                self.opt['scaling'] = float(a)
            elif o == '-c' or o == '--color':
                if a in ['Gi', 'Bi', 'Ri']:
                   self.opt['color'] = a
                else:
                   print("Invalid color code %s. Available color codes are: Gi, Bi, Ri")
            elif o == '--plot':
                self.opt['plot'] = a.upper()
            elif o == '-h' or o == '--help':
                self.usage()
                exit(0)

        self.opt['baseFolder'] = args[0]
        if self.opt['baseFolder'].endswith('/'):
            self.opt['baseFolder'] = args[0][:-1]

        print("scaling method: %s" % (self.opt['scaleMethod']))

    def run(self):
        self.printTitle()
        self.processCommands()

        col = Colorize(self.opt)
        col.execute()

# main program


if __name__ == '__main__':

    app = MainApp(argv)
    app.run()

# end main.

