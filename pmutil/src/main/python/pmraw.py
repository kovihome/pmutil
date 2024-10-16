#!/usr/bin/env python3
#
# PmUtils/pmraw
#
"""
Created on Jun 24, 2024
@author: kovi
"""
import os
from datetime import datetime, timezone
from os.path import exists
from sys import argv

import numpy as np
import rawpy
from astropy.io import fits

import pmbase as pm

class RawConverter:
    def __init__(self, opt, defhdrs=None):
        """
        default headers: CREATOR, IMAGETYPE, OBSERVER, TELESCOP, INSTRUME
        """
        self.colors = opt["color"]
        self.overwrite = opt["overwrite"]
        self.def_headers = defhdrs if defhdrs is not None else {}

    def add_headers(self, exif, header, color):
        def hdr(key, data, comment):
            header[key] = data
            header.comments[key] = comment

        # FITS description comment
        header.add_comment(
            "FITS (Flexible Image Transport System) format is defined in 'Astronomy and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H ")

        # BZERO, BSCALE
        hdr("BZERO", float(0), "offset data range to that of unsigned short")
        hdr("BSCALE", float(1), "default scaling factor")

        # CREATOR
        # TODO: app name and version from parameters
        cr = self.def_headers["CREATOR"] if "CREATOR" in self.def_headers else "pmutils 1.2"
        hdr("CREATOR", cr, f"Created by {cr}")

        # IMAGETYP
        # TODO: image type from file name
        imt = self.def_headers["IMAGETYP"] if "IMAGETYP" in self.def_headers else "LIGHT"
        hdr("IMAGETYP", imt, "Type of exposure")

        # PHOTSYS
        hdr("PHOTSYS", "Instrumental", "Photometry filter system")

        # FILTER
        hdr("FILTER", color, "Spectral filter or colorspace component")

        # XBINNING, YBINNING
        hdr("XBINNING", 1, "X axis binning factor")
        hdr("YBINNING", 1, "Y axis binning factor")

        # DATE-OBS
        # DATE-IMG, DATE-LOC
        dt_str = None
        if "TimeStamp" in exif:
            dt_str = exif["TimeStamp"] + "0000"  # datetime string
            if "OffsetTime" in exif:  # timezone
                tz_str = exif["OffsetTime"].replace(":", "")  # timezone string
                dt_str += tz_str
                dt = datetime.strptime(dt_str, "%Y:%m:%d %H:%M:%S.%f%z")  # datetime object
            else:
                dt = datetime.strptime(dt_str, "%Y:%m:%d %H:%M:%S.%f")  # datetime object
        elif "SubSecDateTimeOriginal" in exif:
            dt_str = exif["SubSecDateTimeOriginal"]
            dt = datetime.strptime(dt_str, "%Y:%m:%d %H:%M:%S.%f")  # datetime object

        if dt_str is not None:
            dt_s = dt.replace(tzinfo=None).isoformat().rstrip('0')
            dt_utc_s = dt.astimezone(timezone.utc).replace(tzinfo=None).isoformat().rstrip('0')
            hdr("DATE-LOC", dt_s, "Time of observation (local)")
            hdr("DATE-IMG", dt_s, "Time of observation (local)")
            hdr("DATE-OBS", dt_utc_s, "Time of observation (UTC)")

        # EXPTIME, EXPSURE
        if "ExposureTime" in exif:
            et = str(exif["ExposureTime"])
            if "/" in et:
                exposure_time = float(et.split("/")[0]) / float(et.split("/")[1])
            else:
                exposure_time = float(et)
            hdr("EXPSURE", exposure_time, "[s] Exposure time")
            hdr("EXPTIME", exposure_time, "[s] Exposure time")

        # INSTRUME
        inst = None
        if "INSTRUME" in self.def_headers:
            inst = self.def_headers["INSTRUME"]
        elif "Model" in exif:
            inst = exif["Model"]
        if inst is not None:
            hdr("INSTRUME", inst, "Imaging instrument name")

        # TELESCOP
        scope = None
        if "TELESCOP" in self.def_headers:
            scope = self.def_headers["TELESCOP"]
        elif "LensModel" in exif and exif["LensModel"] != "":
            scope = exif["LensModel"]
        if scope is not None:
            hdr("TELESCOP", scope, "Name of telescope")

        # ISO, GAIN
        if "ISO" in exif:
            hdr("GAIN", int(exif["ISO"]), "Sensor gain (ISO)")
            hdr("ISO", int(exif["ISO"]), "ISO speed")

        # CCD-TEMP
        if "CameraTemperature" in exif:
            hdr("CCD-TEMP", exif["CameraTemperature"], "CCD Temperature (Celsius)")

        # XPIXSZ, YPIXSZ
        if "FocalPlaneXResolution" in exif and "FocalPlaneYResolution" in exif:
            hdr("XPIXSZ", 25400 / float(exif["FocalPlaneXResolution"]), "[um] Pixel X axis size")
            hdr("YPIXSZ", 25400 / float(exif["FocalPlaneYResolution"]), "[um] Pixel Y axis size")

        # APERTURE
        # TODO: aperture

        # FOCUS, FOCALLEN
        # if "FocalLength" in exif:
        #     fl = exif["FocalLength"]
        #     if "/" in fl:
        #         focal_length = float(fl.split("/")[0]) / float(fl.split("/")[1])
        #     else:
        #         focal_length = float(fl)
        #     hdr("FOCALLEN", focal_length, "[mm] Focal length")
        #     hdr("FOCUS", focal_length, "[mm] Focal length")

        # OBSERVER
        ob = None
        if "OBSERVER" in self.def_headers:
            ob = self.def_headers["OBSERVER"]
        elif "Artist" in exif:
            ob = exif["Artist"]
        if ob is not None:
            hdr("OBSERVER", ob, "Observer name")

        # WB from BlueBalance and RedBalance
        wb = None
        if color == "Bi" and "BlueBalance" in exif:
            wb = float(exif["BlueBalance"])
        elif color == "Ri" and "RedBalance" in exif:
            wb = float(exif["RedBalance"])
        elif color == "Gi":
            wb = 1.0
        if wb is not None:
            hdr("WB", wb, f"White balance for {color} filter")

    def get_exif_data(self, raw_filename: str):
        exif_data = {}
        exif = pm.invoke(f"exiftool -s -g {raw_filename}")
        for line in exif.split('\n'):
            r = line.split(':', 1)
            # print(r)
            if len(r) > 1:
                tag = r[0].strip()
                data = r[1].strip()
                exif_data[tag] = data
        return exif_data

    def convert_channel(self, ch, color, exif_data, fileBaseName):
        FITS_NAME = f"{fileBaseName}-{color}.fits"
        if self.overwrite and exists(FITS_NAME):
            os.remove(FITS_NAME)

        if not exists(FITS_NAME):
            # channel = ch.astype(np.float32)
            channel = ch.astype(np.int16)
            hdu = fits.PrimaryHDU(channel)
            self.add_headers(exif_data, hdu.header, color)
            hdu.writeto(FITS_NAME, overwrite=True)
        else:
            print(f"{FITS_NAME} file already exists.")

    def convert(self, rawFilename: str):
        fileBaseName = rawFilename[:rawFilename.rindex('.')]

        # Load the CR3 image using rawpy
        print(f"Converting {rawFilename} to {fileBaseName}-*.fits")
        with rawpy.imread(rawFilename) as raw:
            # Postprocess the raw image to get a numpy array, with half size and no demosaic
            rgb = raw.postprocess(half_size=True, output_bps=16, output_color=rawpy.ColorSpace.XYZ, user_flip=2,
                                  no_auto_bright=True, no_auto_scale=True, gamma=(1,1), use_camera_wb=True)

        # Get EXIF data from the CR3 file
        exif_data = self.get_exif_data(rawFilename)

        for index, color in enumerate(["Ri", "Gi", "Bi"]):
            if color in self.colors:
                self.convert_channel(rgb[:, :, index], color, exif_data, fileBaseName)


if __name__ == '__main__':
    print("pmraw - convert raw to FITS format, using rawpy/LibRaw")
    print()

    filename = argv[1]
    opt = {
        'color': ['Ri', 'Gi', 'Bi'],
        'overwrite': True
    }
    converter = RawConverter(opt)
    converter.convert(filename)

# end main.
