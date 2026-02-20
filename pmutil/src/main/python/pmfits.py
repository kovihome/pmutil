#!/usr/bin/env python3
#
# PmUtils/pmfits
#
"""
Created on Jun 24, 2024
@author: kovi
@description: FITS manipulation
"""
from getopt import GetoptError, getopt
from sys import argv
from astropy.io import fits
from astropy.time import Time

FITS_HEADER_BZERO = "BZERO"
FITS_HEADER_BSCALE = "BSCALE"
FITS_HEADER_NAXIS = "NAXIS"
FITS_HEADER_NAXIS1 = "NAXIS1"
FITS_HEADER_NAXIS2 = "NAXIS2"
FITS_HEADER_HISTORY = "HISTORY"
FITS_HEADER_FILTER = "FILTER"
FITS_HEADER_BAYER = "BAYERPAT"
FITS_HEADER_CAMERA_TEMPERATURE = "CCD-TEMP"
FITS_HEADER_DATE_OBS = "DATE-OBS"
FITS_HEADER_DATE_END = "DATE-END"
FITS_HEADER_DATE_AVG = "DATE-AVG"
FITS_HEADER_MJD_OBS = "MJD-OBS"
FITS_HEADER_MJD_END = "MJD-END"
FITS_HEADER_MJD_AVG = "MJD-AVG"
FITS_HEADER_EXPTIME = "EXPTIME"
FITS_HEADER_INSTRUMENT = "INSTRUME"
FITS_HEADER_OBSERVER = "OBSERVER"
FITS_HEADER_TELESCOPE = "TELESCOP"

FITS_HEADER_CREATOR = "CREATOR"
FITS_HEADER_IMAGE_TYPE = "IMAGETYP"
FITS_HEADER_PHOTOMETRY_SYSTEM = "PHOTSYS"
FITS_HEADER_XBIN = "XBINNING"
FITS_HEADER_YBIN = "YBINNING"
FITS_HEADER_DATE_LOCAL = "DATE-LOC"
FITS_HEADER_DATE_IMAGE = "DATE-IMG"
FITS_HEADER_EXPOSURE = "EXPSURE"
FITS_HEADER_GAIN = "GAIN"
FITS_HEADER_ISO = "ISO"
FITS_HEADER_XPIX_SIZE = "XPIXSZ"
FITS_HEADER_YPIX_SIZE = "YPIXSZ"
FITS_HEADER_WHITE_BALANCE = "WB"

FITS_HEADER_DATE_MID = "DATE-MID" # --> DATE-AVG
FITS_HEADER_STACK_COUNT = "STACKCNT"
FITS_HEADER_COMBINATION_COUNT = "NCOMBINE"
FITS_HEADER_COMBINATION_MODE = "MCOMBINE"
FITS_HEADER_CAMERA_TEMPERATURE_STD = "STD-TEMP"

# WCS generated headers
FITS_HEADER_WCS_CRVAL1 = "CRVAL1"
FITS_HEADER_WCS_CRVAL2 = "CRVAL2"
FITS_HEADER_WCS_CTYPE1 = "CTYPE1"
FITS_HEADER_WCS_CTYPE2 = "CTYPE2"

# SeeStar specific headers
FITS_HEADER_TOTALEXP = "TOTALEXP"


class FitsSeparator:
    def __init__(self, opt):
        self.colors = opt["color"]
        self.overwrite = opt["overwrite"]

    def saveChannelFits(self, channel_data, original_header, channel, base_filename):
        if channel in self.colors:
            new_header = original_header.copy()
            new_header[FITS_HEADER_FILTER] = channel
            if FITS_HEADER_BAYER in new_header:
                bayer_pattern = new_header[FITS_HEADER_BAYER]
                del new_header[FITS_HEADER_BAYER]
                new_header[FITS_HEADER_HISTORY] = f"Debayering with Bayer pattern {bayer_pattern}"
            hdu = fits.PrimaryHDU(data=channel_data, header=new_header)
            hdu.writeto(f'{base_filename}-{channel}.fits', overwrite=True)

    def separate(self, fits_filename:str, output_path:str):
        # open FITS file
        with fits.open(fits_filename) as hdul:
            new_base_filename = f"{output_path}/{fits_filename[:fits_filename.rfind('.')]}"
            # Get header and image data
            data = hdul[0].data  # vagy hdul[ext].data, ha nem az első HDU-ban van
            original_header = self.adjust_headers(hdul[0].header)
            naxis = int(original_header[FITS_HEADER_NAXIS])
            bayer_pattern = original_header[FITS_HEADER_BAYER] if FITS_HEADER_BAYER in original_header else None
            if naxis == 3:
                if bayer_pattern is not None:
                    self.separateBayer3D(data, original_header, new_base_filename)
                else:
                    self.separateFlat3D(data, original_header, new_base_filename)
            else:
                if bayer_pattern is not None:
                    self.separateBayer2D(data, original_header, new_base_filename)
                else:
                    self.copyWithColor(data, original_header, new_base_filename)

    @staticmethod
    def adjust_headers(header):
        # adjust EXPTIME, when STACKCNT and TOTALEXP exist
        if FITS_HEADER_STACK_COUNT in header:
            if FITS_HEADER_TOTALEXP in header:
                totalexp = header[FITS_HEADER_TOTALEXP]
                exp = header[FITS_HEADER_EXPTIME]
                if int(totalexp) > int(exp):
                    header[FITS_HEADER_EXPTIME] = totalexp
                    header[FITS_HEADER_EXPOSURE] = totalexp
            if FITS_HEADER_COMBINATION_COUNT not in header:
                header[FITS_HEADER_COMBINATION_COUNT] = header[FITS_HEADER_STACK_COUNT]
                header[FITS_HEADER_COMBINATION_MODE] = "sum"

        # calculate MJD-OBS
        if FITS_HEADER_MJD_OBS not in header:
            dateobs = header[FITS_HEADER_DATE_OBS]
            utc = Time(dateobs, format='isot', scale='utc')
            mjd = utc.mjd
            header[FITS_HEADER_MJD_OBS] = (mjd, "MJD of observation start")

        return header

    @staticmethod
    def separate_channel(data, bayer_pattern, color, second = False):
        if not bayer_pattern:
            return data
        ci = bayer_pattern.index(color) if not second else bayer_pattern.rindex(color)
        xpi = ci // 2
        ypi = ci % 2
        return data[xpi::2, ypi::2]

    def separateBayer3D(self, data, original_header, new_base_filename:str):
        bayer_pattern = original_header[FITS_HEADER_BAYER] if FITS_HEADER_BAYER in original_header else None

        if 'Bi' in self.colors:
            B = self.separate_channel(data[0, :, :], bayer_pattern, 'B')
            self.saveChannelFits(B, original_header, 'Bi', new_base_filename)

        if 'Gi' in self.colors:
            if bayer_pattern:
                G1 = self.separate_channel(data[1, :, :], bayer_pattern, 'G')
                G2 = self.separate_channel(data[1, :, :], bayer_pattern, 'G', True)
                G = (G1 + G2) // 2
            else:
                G =  data[1, :, :]
            self.saveChannelFits(G, original_header, 'Gi', new_base_filename)

        if 'Ri' in self.colors:
            R = self.separate_channel(data[2, :, :], bayer_pattern, 'R')
            self.saveChannelFits(R, original_header, 'Ri', new_base_filename)


    def separateBayer2D(self, data, original_header, new_base_filename:str):

        # Separate color channels
        # ch = [data[0::2, 0::2], data[0::2, 1::2], data[1::2, 0::2], data[1::2, 1::2]]

        # Apply Bayer pattern to resolve channel colors
        bayer_pattern = original_header[FITS_HEADER_BAYER] if FITS_HEADER_BAYER in original_header else None
        # B = ch[bayer_pattern.index('B')]
        # G1 = ch[bayer_pattern.index('G')]
        # G2 = ch[bayer_pattern.rindex('G')]
        # R = ch[bayer_pattern.index('R')]
        if 'Bi' in self.colors:
            B = self.separate_channel(data, bayer_pattern, 'B')
            self.saveChannelFits(B, original_header, 'Bi', new_base_filename)

        if 'Gi' in self.colors:
            G1 = self.separate_channel(data, bayer_pattern, 'G')
            G2 = self.separate_channel(data, bayer_pattern, 'G', True)
            G = (G1 + G2) // 2
            self.saveChannelFits(G, original_header, 'Gi', new_base_filename)

        if 'Ri' in self.colors:
            R = self.separate_channel(data, bayer_pattern, 'R')
            self.saveChannelFits(R, original_header, 'Ri', new_base_filename)

    def separateFlat3D(self, data, original_header, new_base_filename):
        for i, color in enumerate(['Bi', 'Gi', 'Ri']):
            if color in self.colors:
                self.saveChannelFits(data[i, :, :], original_header, color, new_base_filename)

    def copyWithColor(self, data, original_header, new_base_filename):
        color = original_header[FITS_HEADER_FILTER] if FITS_HEADER_FILTER in original_header else 'Gi'
        if color in self.colors:
            self.saveChannelFits(data, original_header, color, new_base_filename)

if __name__ == '__main__':
    print("pmfits - separate 3D or 2D Bayered FITS into 3 color channels")
    print()

    opt = {
        'color': ['Ri', 'Gi', 'Bi'],
        'overwrite': True
    }

    try:
        optlist, args = getopt(argv[1:], "c:h", ['color=', 'help'])
    except GetoptError:
        print('Invalid command line options')
        exit(0)

    for o, a in optlist:
        if a[:1] == ':':
            a = a[1:]
        elif o == '-c' or o == '--color':
            if a in ['Gi', 'Bi', 'Ri']:
                opt['color'] = [a]
            else:
                print(f"Invalid color code {a}. Available color codes are: Gi, Bi, Ri")
        elif o == '-h' or o == '--help':
            print("Usage: pmfits [OPTIONS]... [FITS_FILE] [DESTINATION_FOLDER]")
            print()
            print("Mandatory arguments to long options are mandatory for short options too.")
            print("  -c,  --color arg    use selected color for all channel ; it results monochrome image")
            print("  -h,  --help         print this page")
            print()
            exit(0)

    print(args)
    filename = args[0]
    folder = args[1] if len(args) > 1 else "."
    if folder.endswith('/'):
        folder = folder.rstrip('/')

    separator = FitsSeparator(opt)
    separator.separate(filename, folder)

# end main.
