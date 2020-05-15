#!/usr/bin/env python3
#
# PmUtils/pmhotpix
#
'''
Created on Feb 22, 2020

@author: kovi
'''
from getopt import getopt, GetoptError
from sys import argv
from os.path import exists
from astropy.io import fits
from astropy.table import Table
import numpy as np
from scipy.ndimage import gaussian_filter, median_filter
from astropy.stats import sigma_clipped_stats

from pmbase import printError, printInfo, loadPplSetup, Blue, Color_Off, BGreen


class BadPixelDetector:

    SKY_SIGMA = 3.0

    def __init__(self):
        pass

    def detectBadPixels(self, fileName):
        hdu = fits.open(fileName)[0]
        print("image size: %d x %d" % (hdu.header['NAXIS1'], hdu.header['NAXIS2']))
        
        image = hdu.data

        #TODO: ez nem lesz jo flat-nel, oda komolyabb bg estimation kell
        mean,median,std = sigma_clipped_stats(image, sigma=self.SKY_SIGMA)
        print("mean: %4.1f, median: %4.1f, std: %4.1f" % (mean, median, std))

        blurredImage = gaussian_filter(image, sigma=2)
        diff = image - blurredImage
        pthre = mean + 16 * std
        nthre = mean - 16 * std
        print("thre [16 * std + mean]: %4.1f / %4.1f" % (pthre, nthre))
        maxpx = 2**14-1
        hotPixels = np.nonzero((diff>pthre) | (diff<nthre))
        #print('hot pixels:')
        #print('    X      Y      ADU      SNR   Level  Type')
        badPixels = []
        countHot = 0
        countDead = 0
        for j in range(len(hotPixels[0])):
            v = image[hotPixels[0][j]][hotPixels[1][j]]
            if v > pthre:
                #print("%5d  %5d  %7.1f  %7.1f  %5.1f%%   hot" % (hotPixels[1][j], hotPixels[0][j], v, (v-mean)/std, 100*v/maxpx))
                badPixels.append({'x': hotPixels[1][j], 'y': hotPixels[0][j], 'adu': v, 'snr': (v-mean)/std, 'level': 100*v/maxpx, 'type': 'hot'})
                countHot += 1
            elif v < nthre:
                #print("%5d  %5d  %7.1f  %7.1f  %5.1f%%  dead" % (hotPixels[1][j], hotPixels[0][j], v, (mean-v)/std, -100*v/maxpx))
                badPixels.append({'x': hotPixels[1][j], 'y': hotPixels[0][j], 'adu': v, 'snr': (mean-v)/std, 'level': -100*v/maxpx, 'type': 'dead'})
                countDead += 1
        print()
        print("Hot pxs: %d, dead pixels: %d" % (countHot, countDead))
        return badPixels

    def saveBadPixels(self, badPixels, cosmeticFileName):
        cf = open(cosmeticFileName, 'w')
        cf.write("    X      Y      ADU      SNR   LEVEL  TYPE\n")
        for bp in badPixels:
            cf.write("%5d  %5d  %7.1f  %7.1f  %5.1f%%  %s\n" % (bp['x'], bp['y'], bp['adu'], bp['snr'], bp['level'], bp['type']))
        cf.close()

    def execute(self, fileName, color):
        bpxs = self.detectBadPixels(fileName)
        cosmeticFileName = fileName[:fileName.rindex('/')] + '/bad_pixels-' + color
        self.saveBadPixels(bpxs, cosmeticFileName)
        return bpxs
        


class BadPixelEliminator:

    badPixels = []

    def __init__(self, badPixelFileName = None):
        if badPixelFileName:
            self.loadBadPixels(badPixelFileName)

    def loadBadPixels(self, badPixelFileName):
        table = Table.read(badPixelFileName, format = 'ascii')
        for j in range(len(table)):
            self.badPixels.append({ 'x': int(table[j]['X']), 'y': int(table[j]['Y']) })

    def process(self, fitsFileName):
        if exists(fitsFileName):
            hdul = fits.open(fitsFileName, mode='update')
            data = hdul[0].data
            print('data size',data.shape)
            median = median_filter(data, 3)
            print('median size',median.shape)
            for px in self.badPixels:
                #print(px['x'],px['y'])
                data[px['y'], px['x']] = median[px['y'], px['x']]
            hdul.flush()    
            hdul.close()    

            

class MainApp:

    opt = {
        'use'     : None,
        'fileName': None,
        }

    def __init__(self, argv):
        self.argv = argv
        pass

    def printTitle(self):
        print()
        print(BGreen + "pmhotpix, version 1.1.0 " + Color_Off)
        print(Blue + "Detect and eliminate bad pixels." + Color_Off)
        print()

    def usage(self):
        print("Usage: pmhotpix [OPTIONS]... [FILE_NAME]")
        print("Detect and eliminate bad pixels in a FITS image.")
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -h,  --help        print this page")
        print()

    def processCommands(self):
        try:
            optlist, args = getopt (self.argv[1:], "u:h", ['use=', 'help'])
        except GetoptError:
            printError ('Invalid command line options')
            return

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            if o == '-u' or '--use':
                self.opt['use'] = a
            elif o == '-h' or o == '--help':
                self.usage()
                exit(0)

        self.opt['fileName'] = args[0]

    def run(self):
        self.printTitle()
        self.processCommands()
        
        colors = ['Gi','Bi','Ri']
        color = 'Gi'
        for c in colors:
            if c in self.opt['fileName']:
                color = c

        if self.opt['use'] == None:

            bpd = BadPixelDetector()
            bpd.execute(self.opt['fileName'], color)

        else:

            bpe = BadPixelEliminator(self.opt['use'])
            bpe.process(self.opt['fileName'])


# main program

if __name__ == '__main__':

    app = MainApp(argv)
    app.run()

# end main.

