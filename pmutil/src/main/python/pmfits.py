#!/usr/bin/env python3
#
# PmUtils/pmfits
#
"""
Created on Jun 24, 2024
@author: kovi
@description: FITS manipulation
"""
from sys import argv
from astropy.io import fits

class FitsSeparator:
    def __init__(self, opt):
        self.colors = opt["color"]
        self.overwrite = opt["overwrite"]

    def saveChannelFits(self, channel_data, original_header, channel, base_filename):
        new_header = original_header.copy()
        new_header['FILTER'] = channel
        hdu = fits.PrimaryHDU(data=channel_data, header=new_header)
        hdu.writeto(f'{base_filename}-{channel}.fits', overwrite=True)
        pass

    def separate(self, fits_filename:str, output_path:str):

        # FITS file megnyitása
        hdul = fits.open(fits_filename)
        original_header = hdul[0].header
        data = hdul[0].data  # vagy hdul[ext].data, ha nem az első HDU-ban van

        # Feltételezve, hogy a data shape-je (3, Y, X) - azaz 3 színcsatorna
        # Ellenőrizzük:
        print(data.shape)  # Például: (3, 1024, 1024)
        new_base_filename = f"{output_path}/{fits_filename[:fits_filename.rfind('.')]}"

        for n, channel in enumerate(['Bi', 'Gi', 'Ri']):
            self.saveChannelFits(data[n, :, :], original_header, channel, new_base_filename)

        hdul.close()


if __name__ == '__main__':
    print("pmfits - separate 3D FITS into 3 color channels")
    print()

    filename = argv[1]
    folder = argv[2] if len(argv) > 2 else "."
    if folder.endswith('/'):
        folder = folder.rstrip('/')
    opt = {
        'color': ['Ri', 'Gi', 'Bi'],
        'overwrite': True
    }
    separator = FitsSeparator(opt)
    separator.separate(filename, folder)

# end main.
