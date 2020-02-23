#!/usr/bin/env python3
#
# PmUtils/pmcolorize
#
'''
Created on Feb 22, 2020

@author: kovi
'''
from getopt import getopt, GetoptError
from sys import argv
from matplotlib.image import imsave
from astropy.io import fits
from numpy import zeros

from pmbase import printError, printWarning, printInfo, loadPplSetup, Blue, Color_Off, BGreen
import img_scale

class Colorize:

    opt = {}                 # command line options
    pplSetup = {}            # PPL setup from ppl-setup config file

    scaleMethods = ['linear', 'sqrt', 'log', 'asinh']

    wb = [1.15, 1.0, 0.79]   # for Canon EOS 1100D
#   wb = [1.12, 1.0, 0.77]   # for Canon EOS 350D
#   wb = [1.04, 1.0, 0.72]   # for Canon EOS 450D
#   wb = [0.98, 1.0, 0.82]   # for Canon EOS 1300D
#   wb = [1.01, 1.0, 0.87]   # for Canon EOS 5D Mark III
#   wb = [1.20, 1.0, 0.72]   # for Canon EOS 60D

    SKY_SIGMA = 3.0
    SKY_CONVERGENCE = 0.01
    SCALE_RANGE = 10000.0
    SCALE_METHOD = 'linear'

    def __init__(self, opt):
        self.opt = opt
        if 'scaleMethod' in self.opt:
            self.SCALE_METHOD = self.opt['scaleMethod']

    def execute(self):
        self.pplSetup = loadPplSetup()

        baseFolder = self.opt['baseFolder']

        colors = ['Ri', 'Gi', 'Bi']
        rgb = []
        seqFolder = self.pplSetup['SEQ_FOLDER_NAME']
        # TODO: before addition, images should be reregistered with Gi as reference, for more exact alignment

        for c in colors:
            fitsFileName = "%s%s/Combined-%s.fits" % (baseFolder, seqFolder, c)
            print("%s color source image: %s" % (c, fitsFileName))
            rgb.append(fits.open(fitsFileName)[0].data)

        img = zeros((rgb[1].shape[0], rgb[1].shape[1], 3), dtype = float)
        for j in range(3):
            smin, it = img_scale.sky_mean_sig_clip(rgb[j], self.SKY_SIGMA, self.SKY_CONVERGENCE)
            smax = int(self.wb[j] * (self.SCALE_RANGE + smin))
            if self.SCALE_METHOD == self.scaleMethods[0]:
                img[:,:,j] = img_scale.linear(rgb[j], scale_min=smin, scale_max=smax)
            elif self.SCALE_METHOD == self.scaleMethods[1]:
                img[:,:,j] = img_scale.sqrt(rgb[j], scale_min=smin, scale_max=smax)
            elif self.SCALE_METHOD == self.scaleMethods[2]:
                img[:,:,j] = img_scale.log(rgb[j], scale_min=smin, scale_max=smax)
            elif self.SCALE_METHOD == self.scaleMethods[3]:
                img[:,:,j] = img_scale.asinh(rgb[j], scale_min=smin, scale_max=smax)

        # TODO: with log method, save is failed: Floating point image RGB values must be in the 0..1 range.
        # need som normalization

        imgFileName = "%s%s/Combined-color.jpg" % (baseFolder, seqFolder)
        imsave(imgFileName, img)
        printInfo("Color image %s created." % (imgFileName))



class MainApp:
    
    # TODO: wb should be command line option somehow
    opt = {
        'scaleMethod': 'linear',
        'baseFolder': None,
        }


    def __init__(self, argv):
        self.argv = argv
        pass

    def printTitle(self):
        print()
        print(BGreen + "pmcolorize, version 1.1.0 " + Color_Off)
        print(Blue + "Make color jpeg image from calibrated FITS images." + Color_Off)
        print()


    def usage(self):
        print("Usage: pmcolorize [OPTIONS]... [BASE_FOLDER]")
        print("Make color jpeg image from calibrated FITS images.")
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -m,  --scale-method     scaling method, available values are: linear, sqrt, log, asinh")
        print("  -h,  --help             print this page")
        print()

    def processCommands(self):
        try:
            optlist, args = getopt (self.argv[1:], "m:h", ['--scale-method', '--help'])
        except GetoptError:
            printError ('Invalid command line options')
            return

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-m':
                if not a in Colorize.scaleMethods:
                    printError("Invalid scaling method value: %s" % (a))
                self.opt['scaleMethod'] = a
            elif o == '-h':
                self.usage()
                exit(0)

        self.opt['baseFolder'] = args[0]
        if not self.opt['baseFolder'].endswith('/'):
            self.opt['baseFolder'] += '/'

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

