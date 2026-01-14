#!/usr/bin/env python3
#
# PmUtils/pmmerge
#
"""
Created on Mar 17, 2021

@author: kovi
"""

from getopt import getopt, GetoptError
from sys import argv
from os import remove
from os.path import exists
from math import log10, cos

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table, join, setdiff
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree

import pmbase as pm
from pmviz import VizierQuery, UCAC4, APASS

SNR = 7.2                       # SNR to determine limiting magnitude
PX_PREC = 2.0                   # Pixel coordinate matching precision
SEP_PREC = 2.0 / 3600.0         # Coordinate matching precisin in degrees

# sextractor catalog field names
CAT_ID_FIELD = "NUMBER"
MAG_FIELD = "MAG_ISOCOR"
MAGERR_FIELD = "MAGERR_ISOCOR"
FLUX_FIELD = "FLUX_ISOCOR"
FLUXERR_FIELD = "FLUXERR_ISOCOR"
XPOS_FIELD = "XWIN_IMAGE"
YPOS_FIELD = "YWIN_IMAGE"
RA_FIELD = "ALPHA_J2000"
DEC_FIELD = "DELTA_J2000"

NO_MAG = "99.0"

class CatalogMatcher:
    folder = './'
    baseName = None
    colors = ['Ri', 'Gi', 'Bi']
    showGraphs = False
    saveGraphs = False

    # TODO: ezt a kettot konfigba tenni
    HMG_MAX_ERR = 0.15  # mg

    LIMIT_THRE = 0.5  # mg

    IMAGE_BORDER_SIZE = 10  # pixels

    MG_LIMIT = 17.0  # TODO: nem lehet-e ezt szamolni?

    COORD_MATCH_THRE = 2.0 / 60.0  # degrees

    catUCAC4 = None

    def __init__(self, baseName, folder, colors, show, save, loggingMode=pm.Logger.LOG_MODE_INFO):
        self.baseName = baseName
        if folder is not None:
            self.folder = folder
        if colors is not None:
            self.colors = colors
        # self.viz = VizUCAC4(self.MG_LIMIT)
        self.catUCAC4 = VizierQuery(cat=UCAC4, limit=self.MG_LIMIT)
        self.showGraphs = show
        self.saveGraphs = save
        self.log = pm.Logger(mode=loggingMode)
        self.mgLimits = None
        self.hmgPlot = None

    def hmg(self, pmTable, nr):
        mgIsoCorrArr = pmTable['MAG_ISOCOR']
        mgIsoCorrErrArr = pmTable['MAGERR_ISOCOR']
        mgBestArr = pmTable['MAG_BEST']
        mgBestErrArr = pmTable['MAGERR_BEST']
        logm1 = []
        loge1 = []
        logm2 = []
        loge2 = []
        for j in range(len(pmTable)):
            if mgIsoCorrArr[j] <= 0.0:
                logm1.append(log10(-mgIsoCorrArr[j]))
                loge1.append(log10(mgIsoCorrErrArr[j]))
            if mgBestArr[j] <= 0.0:
                logm2.append(log10(-mgBestArr[j]))
                loge2.append(log10(mgBestErrArr[j]))

        m1, b1 = pm.linfit(logm1, loge1)
        hmg1 = -1.0 * 10.0 ** ((log10(self.HMG_MAX_ERR) - b1) / m1)

        m2, b2 = pm.linfit(logm2, loge2)
        hmg2 = -1.0 * 10.0 ** ((log10(self.HMG_MAX_ERR) - b2) / m2)

        plotcolor = 'gbr'[nr]
        self.hmgPlot.add(logm1, loge1, [m1, b1], f'log mi ({plotcolor.upper()}i)', f'log ei ({plotcolor.upper()}i)',
                         plotcolor, title="Mg limit")

        return [hmg1, hmg2]

    def calcMgLimit(self, fileName, nr):
        t = Table.read(fileName, format='ascii.sextractor')
        mgs = self.hmg(t, nr)
        return mgs[0]

    def readFrameSize(self, fitsFileName):
        h = pm.getFitsHeaders(fitsFileName, ['NAXIS1', 'NAXIS2'])
        xf = int(h['NAXIS1'])
        yf = int(h['NAXIS2'])
        return xf, yf

    def determineOnFrameStatus(self, x_obj, y_obj, im_w, im_h):
        if x_obj < 0 or x_obj >= im_w or y_obj < 0 or y_obj >= im_h:
            return 'O'
        if x_obj < self.IMAGE_BORDER_SIZE or x_obj >= im_w - self.IMAGE_BORDER_SIZE or y_obj < self.IMAGE_BORDER_SIZE \
                or y_obj >= im_h - self.IMAGE_BORDER_SIZE:
            return 'B'
        return '-'

    def calcVisFlag(self, mg, limit):
        """
        flags I - INVISIBLE, B - BELOW_LIMIT, N - NEAR_LIMIT, S - SATURATED
        """
        # TODO 1.3: saturated objects
        mgf = float(mg)
        if mgf > 0.0:
            return 'I'
        elif mgf > limit:
            return 'B'
        elif mgf > limit - self.LIMIT_THRE:
            return 'N'
        return '-'

    def calcPosFlag(self, x_obj, y_obj, im_w, im_h):
        """
        x_obj # object x position
        y_obj # object y position
        im_w  # image width
        im_h  # image height

        flags B - object is near of image border
        """
        return self.determineOnFrameStatus(x_obj, y_obj, im_w, im_h)

    def calcMgLimits(self, catFileBase: str):
        # calculate mg limits
        self.hmgPlot = pm.Plot(len(self.colors), self.showGraphs, self.saveGraphs)

        self.mgLimits = {}
        s = ''
        for j, clr in enumerate(self.colors):
            self.mgLimits[clr] = self.calcMgLimit(catFileBase + clr + '.cat', j)
            s += f'{clr}: {self.mgLimits[clr]}, '
        self.log.info(f'Limit mg {s[:-2]}')

        self.hmgPlot.showOrSave(self.folder + '/magnitude_limit.png')

    def loadAndFiltercatalogFile(self, catFileBase: str, color: str) -> Table:
        catFn = f"{catFileBase}{color}.cat"
        cat = Table.read(catFn, format="ascii.sextractor")
        self.log.debug(f"Original catalog file size for {color}: {len(cat)}")

        # filter catalog by SNR
        cat_snr = cat[cat[FLUX_FIELD] / cat[FLUXERR_FIELD] > SNR]
        self.log.debug(f"SNR filtered catalog size for Gi: {len(cat_snr)}")

        return cat_snr

    @staticmethod
    def matchCatalogsByPixelCoords(cats: dict[str,Table]) -> Table:
        gCat_snr = cats["Gi"]
        rCat_snr = cats["Ri"]
        bCat_snr = cats["Bi"]
        # Match SNR filtered Gi and Ri tables
        # find nearest Ri object for all Gi ones
        g_coords = np.vstack([gCat_snr[XPOS_FIELD], gCat_snr[YPOS_FIELD]]).T
        r_coords = np.vstack([rCat_snr[XPOS_FIELD], rCat_snr[YPOS_FIELD]]).T
        tree = cKDTree(r_coords)
        dist, idx2 = tree.query(g_coords, distance_upper_bound=PX_PREC)
        di = np.vstack([dist, idx2]).T
        # print(f"match size: {len(di)}")

        # filter out large distances
        d_mask = di[:, 0] < PX_PREC
        # di_masked = di[d_mask]

        # Get matched records in Gi
        g_matched_to_r = gCat_snr[d_mask]

        # Get non-matched records in Gi
        g_not_matched_to_r = gCat_snr[~d_mask]

        # Get matched records in Ri
        # r_matched_to_g = rCat_snr[idx2[d_mask]]

        # Get Ri records not matched to Gi
        r_ix = idx2[d_mask]
        fi2 = list(set(range(0, len(rCat_snr))) - set(r_ix))
        r_not_matched_to_g = rCat_snr[fi2]

        # merge G and R table indices
        from astropy.table import vstack
        rg_merged = Table({
            "G_NUMBER": g_matched_to_r[CAT_ID_FIELD],
            "R_NUMBER": rCat_snr[r_ix][CAT_ID_FIELD],
            "X": g_matched_to_r[XPOS_FIELD],
            "Y": g_matched_to_r[YPOS_FIELD]})
        g_only = Table({
            "G_NUMBER": g_not_matched_to_r[CAT_ID_FIELD],
            "R_NUMBER": [0] * len(g_not_matched_to_r),
            "X": g_not_matched_to_r[XPOS_FIELD],
            "Y": g_not_matched_to_r[YPOS_FIELD]
        })
        r_only = Table({
            "G_NUMBER": [0] * len(r_not_matched_to_g),
            "R_NUMBER": r_not_matched_to_g[CAT_ID_FIELD],
            "X": r_not_matched_to_g[XPOS_FIELD],
            "Y": r_not_matched_to_g[YPOS_FIELD]
        })
        rg_merged_xref = vstack([rg_merged, g_only, r_only])

        # Find nearest Bi object for all Gi/Ri ones
        gr_coords = np.vstack([rg_merged_xref['X'], rg_merged_xref['Y']]).T
        b_coords = np.vstack([bCat_snr[XPOS_FIELD], bCat_snr[YPOS_FIELD]]).T
        tree = cKDTree(b_coords)
        dist3, idx23 = tree.query(gr_coords)
        di3 = np.vstack([dist3, idx23]).T

        # Filter out large distances
        d3_mask = di3[:, 0] < PX_PREC
        # di3_masked = di3[d3_mask]

        # Get matched records in Gi/Ri
        rg_matched_to_b = rg_merged_xref[d3_mask]

        # Get non-matched records in Gi/Ri
        rg_not_matched_to_b = rg_merged_xref[~d3_mask]

        # Get matched records in Ri
        # b_matched_to_gr = rCat_snr[idx23[d3_mask]]

        # Get Ri records not matched to Gi
        r3_ix = idx23[d3_mask]
        fi2 = list(set(range(0, len(bCat_snr))) - set(r3_ix))
        b_not_matched_to_gr = bCat_snr[fi2]

        # merge Gi/Ri and Bi table indices
        rgb_merged = Table({
            "G_NUMBER": rg_matched_to_b["G_NUMBER"],
            "R_NUMBER": rg_matched_to_b["R_NUMBER"],
            "B_NUMBER": bCat_snr[r3_ix][CAT_ID_FIELD],
        })
        gr_only = Table({
            "G_NUMBER": rg_not_matched_to_b["G_NUMBER"],
            "R_NUMBER": rg_not_matched_to_b["R_NUMBER"],
            "B_NUMBER": [0] * len(rg_not_matched_to_b),
        })
        b_only = Table({
            "G_NUMBER": [0] * len(b_not_matched_to_gr),
            "R_NUMBER": [0] * len(b_not_matched_to_gr),
            "B_NUMBER": b_not_matched_to_gr[CAT_ID_FIELD],
        })
        rgb_merged_xref = vstack([rgb_merged, gr_only, b_only])

        return rgb_merged_xref

    def buildCombinedCatalog(self, cats: dict[str,Table], indices: Table, imageSize: tuple[int, int]) -> Table:
        gCat = cats["Gi"]
        rCat = cats["Ri"]
        bCat = cats["Bi"]

        cmb = self.createCombinedCat()

        # combine catalogs
        n_auid = 1
        for ix in indices:
            auid = '%03d-FFF-%03d' % (n_auid // 1000, n_auid % 1000)
            gRec = gCat[gCat[CAT_ID_FIELD] == ix["G_NUMBER"]][0] if ix["G_NUMBER"] > 0 else None
            rRec = rCat[rCat[CAT_ID_FIELD] == ix["R_NUMBER"]][0] if ix["R_NUMBER"] > 0 else None
            bRec = bCat[bCat[CAT_ID_FIELD] == ix["B_NUMBER"]][0] if ix["B_NUMBER"] > 0 else None
            ra = gRec[RA_FIELD] if gRec else rRec[RA_FIELD] if rRec else bRec[RA_FIELD]
            decl = gRec[DEC_FIELD] if gRec else rRec[DEC_FIELD] if rRec else bRec[DEC_FIELD]
            ra_s = pm.deg2hexa(float(ra) / 15.0)
            de_s = pm.deg2hexa(float(decl))
            magG = gRec[MAG_FIELD] if gRec else NO_MAG
            magR = rRec[MAG_FIELD] if rRec else NO_MAG
            magB = bRec[MAG_FIELD] if bRec else NO_MAG
            xpos = gRec[XPOS_FIELD] if gRec else rRec[XPOS_FIELD] if rRec else bRec[XPOS_FIELD]
            ypos = gRec[YPOS_FIELD] if gRec else rRec[YPOS_FIELD] if rRec else bRec[YPOS_FIELD]
            # Calculate flags
            # FLAGS & 4 -> VIZ_FLAG Saturated object
            # TODO: FLAGS & 3 -> VIZ_FLAG??? Blended object
            # FLAGS & 24 -> POS_FLAG Close to boundary
            gFlag = int(gRec["FLAGS"]) if gRec else 0
            rFlag = int(rRec["FLAGS"]) if rRec else 0
            bFlag = int(bRec["FLAGS"]) if bRec else 0
            closeToBoundary = ((gFlag | rFlag | bFlag) & 24) > 0
            posFlags = 'B' if closeToBoundary else self.calcPosFlag(float(xpos), float(ypos), imageSize[0], imageSize[1])
            gVisFlag = 'S' if gFlag & 4 > 0 else self.calcVisFlag(magG, self.mgLimits['Gi'])
            bVisFlag = 'S' if bFlag & 4 > 0 else self.calcVisFlag(magB, self.mgLimits['Bi'])
            rVisFlag = 'S' if rFlag & 4 > 0 else self.calcVisFlag(magR, self.mgLimits['Ri'])
            visFlags = gVisFlag + bVisFlag + rVisFlag
            #####
            rec = [auid, '-', 'F', '-', visFlags, posFlags, 'T', ra_s, str(ra), de_s, str(decl),
                   str(magG) if gRec else NO_MAG, str(gRec[MAGERR_FIELD]) if gRec else '-',
                   str(magB) if bRec else NO_MAG, str(bRec[MAGERR_FIELD]) if bRec else '-',
                   str(magR) if rRec else NO_MAG, str(rRec[MAGERR_FIELD]) if rRec else '-',
                   NO_MAG, '-',
                   NO_MAG, '-',
                   NO_MAG, '-']
            cmb.add_row(rec)
            n_auid += 1
        return cmb

    def matchCatalogsAllColors(self, catFileBase: str, imageSize: tuple[int, int]) -> Table:
        # open catalog files
        pmCatalogs = {}
        for cc in ["Gi", "Ri", "Bi"]:
            pmCatalogs[cc] = self.loadAndFiltercatalogFile(catFileBase, cc)

        # match catalogs by pixel coords, return indices
        indices = self.matchCatalogsByPixelCoords(pmCatalogs)

        # merge catalog by indices
        return self.buildCombinedCatalog(pmCatalogs, indices, imageSize)

    def getIndicesForOneColor(self, cat: Table) -> Table:
        cc = self.colors[0]
        indexColumn = cat[CAT_ID_FIELD]
        emptyColumn = [0] * len(cc)
        return Table({
            "G_NUMBER": indexColumn if cc == "Gi" else emptyColumn,
            "R_NUMBER": indexColumn if cc == "Ri" else emptyColumn,
            "B_NUMBER": indexColumn if cc == "Bi" else emptyColumn,
        })

    def matchCatalogsOneColor(self, catFileBase, imageSize: tuple[int, int]):
        # filter catalog by SNR
        pmCatalog = self.loadAndFiltercatalogFile(catFileBase, self.colors[0])

        # create indices for this one catalog
        indices = self.getIndicesForOneColor(pmCatalog)

        # merge catalog by indices
        pmCatalogs = {}
        for cc in ["Gi", "Ri", "Bi"]:
            pmCatalogs[cc] = pmCatalog if cc == self.colors[0] else None
        return self.buildCombinedCatalog(pmCatalogs, indices, imageSize)

    def matchCatalogs(self) -> Table :
        # create photometry files base name
        catFileBase = f"{self.folder}/{pm.setup['PHOT_FOLDER_NAME']}/{self.baseName}-"

        # query image size
        imgSize = self.readFrameSize(f"{catFileBase}{self.colors[0]}.ast.fits")
        self.log.debug(f'Image size: {imgSize[0]} x {imgSize[1]}')

        # calculate mg limits
        self.calcMgLimits(catFileBase)

        # match photometry catalogs by pixel coords
        if len(self.colors) == 3:
            cmb = self.matchCatalogsAllColors(catFileBase, imgSize)
        elif len(self.colors) == 1:
            cmb = self.matchCatalogsOneColor(catFileBase, imgSize)
        else:
            self.log.error(f"No merge strategy for colors {self.colors}")
            return self.createCombinedCat()

        return cmb

    @staticmethod
    def createCombinedCat():
        fnames = ['AUID', 'VIZ_ID', 'ROLE', 'LABEL',
                  'VIZ_FLAG', 'POS_FLAG', 'MATCH_FLAG',
                  'RA', 'RA_DEG', 'DEC', 'DEC_DEG',
                  'MAG_GI', 'ERR_GI', 'MAG_BI', 'ERR_BI', 'MAG_RI', 'ERR_RI',
                  'MAG_B', 'ERR_B', 'MAG_V', 'ERR_V', 'MAG_R', 'ERR_R']
        fdtype = ['U11', 'U16', 'U2', 'U50',
                  'U1', 'U1', 'U1',
                  'U12', 'U12', 'U12', 'U12',
                  'U6', 'U5', 'U6', 'U5', 'U6', 'U5',
                  'U6', 'U5', 'U6', 'U5', 'U6', 'U5']
        return Table(names=fnames, dtype=fdtype)

    def xmatchViz(self, combined):
        # xmatch cmb to UCAC4 catalog
        xmt = self.catUCAC4.xmatch(combined, 'RA_DEG', 'DEC_DEG')
        if xmt is None:
            self.log.error("Accessing Vizier xmatch service for UCAC4 catalog was failed.")
            return False
        self.log.debug(f'UCAC4 xmatch results {len(xmt)} objects')

        combined.add_index('AUID')

        for xr in xmt:
            auid = xr['AUID']
            r = combined.loc[auid]
            r['VIZ_ID'] = 'UCAC4-' + xr['UCAC4']
            r['MATCH_FLAG'] = '-'
            r['MAG_B'] = xr['Bmag']
            r['MAG_V'] = xr['Vmag']
            r['MAG_R'] = xr['Rmag']

    # TODO: query errors from vizier
    #            r['ERR_B'] = xr['e_Bmag']
    #            r['ERR_V'] = xr['e_Vmag']
    #            r['ERR_R'] = xr['e_Rmag']

        # get valid mags and errors from APASS
        catAPASS = VizierQuery(APASS, self.MG_LIMIT)
        apassMatchTable = catAPASS.xmatch(combined, 'RA_DEG', 'DEC_DEG')
        if apassMatchTable is None:
            self.log.error("Accessing Vizier xmatch service for APASS catalog was failed.")
            return False
        self.log.debug(f'APASS xmatch results {len(xmt)} objects')

        fmapApass = {
            "MAG_V": "Vmag",
            "ERR_V": "e_Vmag",
            "MAG_B": "Bmag",
            "ERR_B": "e_Bmag",
            "MAG_R": "rpmag",
            "ERR_R": "e_rpmag"
        }
        ff = "{:.3f}.format"
        for r in apassMatchTable:
            auid = r['AUID']
            br = combined.loc[auid]
            for f in fmapApass.keys():
                if type(r[fmapApass[f]]) != np.ma.core.MaskedConstant:
                    br[f] = np.round(np.float64(r[fmapApass[f]]), 3)

        return True

    def loadRefcat(self):
        refcatFileName = f"{self.folder}/ref.cat"
        if exists(refcatFileName):
            return Table.read(refcatFileName, format='ascii')
        else:
            self.log.warning(f'No ref.cat in the folder {self.folder}')
            return None

    def addFrameCoords(self, cat, fitsFileName):
        # 1. convert cat to fits format
        # baseFolder = fitsFileName.partition(pm.setup['PHOT_FOLDER_NAME'])[0]
        tempFile = self.folder + '/temp/tempCat.fits'
        if not exists(tempFile):
            cat.write(tempFile)

        # 2. calculate ref objects' frame xy points
        wcsFile = fitsFileName.replace('.ast.fits', '.wcs')
        axyFile = fitsFileName.replace('.ast.fits', '.ref.axy')
        pm.invoke("wcs-rd2xy -w %s -i %s -o %s -R RA_DEG -D DEC_DEG" % (
            wcsFile, tempFile, axyFile))  # -f option need argument

        remove(tempFile)

        # 3. merge frame xy point to refCat
        tlen = len(cat)
        cat['X'] = [0.0] * tlen
        cat['Y'] = [0.0] * tlen

        # f = fits.open(axyFile)
        # d = f[1].data
        with fits.open(axyFile) as f:
            d = f[1].data

        for j in range(tlen):
            x, y = d[j]
            cat[j]['X'] = x
            cat[j]['Y'] = y

    @staticmethod
    def getWCS(fitsFileName):
        with fits.open(fitsFileName) as f:
            h = f[0].header
            return WCS(h)

    def filterOutFrameStars(self, cat, fitsFileName):
        # get WCS from the image
        wcs = self.getWCS(fitsFileName)

        # add frame coords to the catalog
        # self.addFrameCoords(cat, fitsFileName)

        # read frame size from fits file
        # fx, fy = self.readFrameSize(fitsFileName)

        # cat['FRAME_STATUS'] = ['-'] * len(cat)
        # for row in cat:
        #     row['FRAME_STATUS'] = self.determineOnFrameStatus(row['X'], row['Y'], fx, fy)
        #
        # foMask = cat['FRAME_STATUS'] != 'O'
        # foCat = cat[foMask]
        # return foCat

        sky = SkyCoord(ra=cat["RA_DEG"], dec=cat["DEC_DEG"], unit="deg")
        pxc = Table(wcs.world_to_pixel(sky))
        # pxc.add_column(list(range(len(sky))), name="INDEX")
        # xmax,ymax = wcs.pixel_shape
        inframe_mask = (pxc['col0'] > 0) & (pxc['col0'] < wcs.pixel_shape[0]) & (pxc['col1'] > 0) & (pxc['col1'] < wcs.pixel_shape[1])
        # pxcout = pxc[~inframe_mask]
        # pxcin = pxc[inframe_mask]
        return cat[inframe_mask]

    @staticmethod
    def matchCatalogByCoords(cat, ra, dec):
        dmin = SEP_PREC
        auid = None
        for row in cat:
            if row['VIZ_ID'] == '-':
                cra = float(row['RA_DEG'])
                cdec = float(row['DEC_DEG'])
                d = pm.quad((cra - ra) * cos(dec), cdec - dec)
                if d < dmin:
                    dmin = d
                    auid = row['AUID']
        #        print(f"Coord matched - auid: {auid}, d: {sqrt(dmin)*60}'")
        return auid, dmin

    @staticmethod
    def updateCmb(cmb, r, auid=None, vizid=None):
        if vizid:
            cr = cmb.loc['VIZ_ID', vizid]
        else:
            cr = cmb.loc['AUID', auid]
        cr['ROLE'] = r['ROLE']
        cr['AUID'] = r['AUID']
        cr['LABEL'] = r['LABEL']
        cr['MAG_B'] = r['MAG_B']
        cr['MAG_V'] = r['MAG_V']
        cr['MAG_R'] = r['MAG_R']
        cr['ERR_B'] = r['ERR_B']
        cr['ERR_V'] = r['ERR_V']
        cr['ERR_R'] = r['ERR_R']

    @staticmethod
    def insertCmb(cmb: Table, r):
        n = [r['AUID'], '-', r['ROLE'], r['LABEL'], 'III', '-', 'N', r['RA'], r['RA_DEG'], r['DEC'], r['DEC_DEG'],
             NO_MAG, NO_MAG, NO_MAG, NO_MAG, NO_MAG, NO_MAG, r['MAG_B'], r['ERR_B'], r['MAG_V'], r['ERR_V'], r['MAG_R'],
             r['ERR_R']]
        cmb.add_row(n)

    def addRefcatData(self, cmb: Table, refcat: Table) -> None:

        del refcat.meta['comments']

        # filter out refcat objects that are not on the frame
        fitsFileName = f"{self.folder}/{pm.setup['PHOT_FOLDER_NAME']}/{self.baseName}-{self.colors[0]}.ast.fits"
        foRefcat = self.filterOutFrameStars(refcat, fitsFileName)

        # match refcat LABEL to cmb VIZ_ID
        refcat_x_cmb = join(foRefcat, cmb, join_type='left', keys_left='LABEL', keys_right='VIZ_ID')
        miss_mask = refcat_x_cmb["AUID_2"].mask == True
        missingRefcat = refcat_x_cmb[miss_mask]

        # match missing records to cmb by coords
        cmb_coords = SkyCoord(ra=cmb["RA_DEG"], dec=cmb["DEC_DEG"], unit="deg")
        refcat_coords = SkyCoord(ra=missingRefcat["RA_DEG_1"], dec=missingRefcat["DEC_DEG_1"], unit="deg")
        match_idx, sep, dist = refcat_coords.match_to_catalog_sky(cmb_coords)
        index = list(range(len(match_idx)))
        # TODO: col2 (dist) is unused
        xref = Table([match_idx, sep, dist, index])
        xm = xref["col1"] < SEP_PREC
        xref_good = xref[xm]
        # missing_to_upgrade = missingRefcat[xref_good["col3"]]

        # xmatch refcat to UCAC4
        # xmt = self.catUCAC4.xmatch(foRefcat, 'RA_DEG', 'DEC_DEG')
        # if xmt is None :
        #     self.log.error("Adding refcat data to combined catalog failed.")
        #     return
        # self.log.debug(f'Refcat xmatch table contains {len(xmt)} records')
        #
        # missingRefcat = setdiff(foRefcat, xmt, 'AUID')
        #
        # cmb.add_index('VIZ_ID')
        #
        # for xr in xmt:
        #     vizId = 'UCAC4-' + xr['UCAC4']
        #     if vizId in cmb['VIZ_ID']:
        #         try:
        #             self.updateCmb(cmb, xr, vizid=vizId)
        #
        #         except KeyError:
        #             self.log.debug(xr)
        #

        # update cmb from refcat data found by coord matching
        # TODO: a mezőcsillagokat meg kell nézni, hogy olyat talált-e meg, aminek nincs VIZ_ID-je, és ha van, de eltérő, akkor Warning, és nem update-elni
        for j in range(len(xref_good)):
            br = cmb[xref_good[j]["col0"]]
            ur = missingRefcat[xref_good[j]["col3"]]
            for f in ["AUID", "ROLE", "LABEL", "MAG_V", "ERR_V", "MAG_B", "ERR_B", "MAG_R", "ERR_R"]:
                f_ix = f"{f}_1"
                if ur[f_ix] != '-':
                    br[f] = ur[f_ix]
        # mergedTable = join(xmt, cmb, keys='AUID', join_type='left')
        # mask = mergedTable['VIZ_ID'] == '-'  # or mergedTable['VIZ_ID'] == '' or mergedTable['VIZ_ID'] == None
        # nt = mergedTable[mask]

        # refcat records not found in cmb at all
        xref_bad = xref[~xm]
        # missing_missing = missingRefcat[xref_bad["col3"]]
        nt = missingRefcat[xref_bad['col3']]
        fields = ["AUID", "ROLE", "RA", "RA_DEG", "DEC", "DEC_DEG", "MAG_V", "ERR_V", "MAG_B", "ERR_B", "MAG_R",
                  "ERR_R", "LABEL"]
        cols = [nt[f"{f}_1"] for f in fields]
        negative_transients = Table(cols, names=fields)

        # save negative transients
        ntFileName = f"{self.folder}/{pm.setup['PHOT_FOLDER_NAME']}/{self.baseName}-{self.colors[0]}.nt"
        negative_transients.write(ntFileName, format='ascii.fixed_width', overwrite=True)

        # self.log.debug(f'Negative transients: {len(nt)} found.')

        # refCat: AUID         ROLE RA             RA_DEG         DEC            DEC_DEG        MAG_B     ERR_B     MAG_V     ERR_V     MAG_R     ERR_R     LABEL
        # cmb:    'AUID','VIZ_ID','ROLE','LABEL','VIZ_FLAG','POS_FLAG','MATCH_FLAG','RA','RA_DEG','DEC','DEC_DEG',
        #        'MAG_GI','ERR_GI','MAG_BI','ERR_BI','MAG_RI','ERR_RI','MAG_B','ERR_B','MAG_V','ERR_V','MAG_R','ERR_R']

        # for r in nt:
        #     n = [r['AUID'], '-', r['ROLE_1'], r['LABEL_1'], r['VIZ_FLAG'], r['POS_FLAG'], r['MATCH_FLAG'], r['RA_1'],
        #          r['RA_DEG_1'], r['DEC_1'], r['DEC_DEG_1'],
        #          r['MAG_GI'], r['ERR_GI'], r['MAG_BI'], r['ERR_BI'], r['MAG_RI'], r['ERR_RI'],
        #          r['MAG_B_1'], r['ERR_B_1'], r['MAG_V_1'], r['ERR_V_1'], r['MAG_R_1'], r['ERR_R_1']]
        #     cmb.add_row(n)
        #
        # self.log.debug(f'Stars missing from refcat: {len(missingRefcat)} found.')
        # for r in missingRefcat:
        #     auid, d = self.matchCatalogByCoords(cmb, float(r['RA_DEG']), float(r['DEC_DEG']))
        #     self.log.debug(f"Match missingRefcat object {r['AUID']} {r['LABEL']} for {auid} within {d * 60}'")
        #     # print(r)  # RBL
        #     if d < 2.0 / 60.0:
        #         self.updateCmb(cmb, r, auid=auid)
        #     else:
        #         self.insertCmb(cmb, r)

    def addMgLimits(self, cmb):
        self.log.debug(f"Instrumental mg limits: Gi = {self.mgLimits['Gi']}, Bi = {self.mgLimits['Bi']}, Ri = {self.mgLimits['Ri']}")
        for cc in self.mgLimits.keys():
            pm.addTableComment(cmb, 'MgLimitInst' + cc, '%7.3f' % (self.mgLimits[cc]))

    def save(self, combined):
        combinedFileName = f"{self.folder}/{pm.setup['PHOT_FOLDER_NAME']}/{self.baseName}.cmb"
        self.log.debug(f"Combined file name: {combinedFileName}")
        combined.write(combinedFileName, format='ascii.fixed_width', delimiter=' ', overwrite=True)

    def process(self):

        # match cat files together
        self.log.print('Combine measurement catalogs all together')
        combined = self.matchCatalogs()

        # match combined catalog to UCAC4, APASS catalogs
        self.log.print('Match combined catalog with UCAC4')
        result = self.xmatchViz(combined)
        if not result:
            print("Error: something wrong with the xmatchViz")
            return False

        # load refcat
        refcat = self.loadRefcat()

        # add refcat data to combined catalog
        if refcat is not None:
            self.log.print('Add ref.cat data to combined catalog')
            self.addRefcatData(combined, refcat)

        self.addMgLimits(combined)

        # save combined catalog
        self.log.print('Save combined catalog')
        self.save(combined)
        
        return True


if __name__ == '__main__':

    class MainApp:

        color = None
        folder = None
        showGraphs = False
        saveGraphs = False

        def __init__(self, argv):
            self.argv = argv
            pm.printInfo('pmmerge, v1.2')

        def processCommands(self):
            try:
                optlist, args = getopt(argv[1:], "c:", ['color=', 'show-graph', 'save-graph'])
            except GetoptError:
                pm.printError('Invalid command line options')
                return

            for o, a in optlist:
                if a[:1] == ':':
                    a = a[1:]
                if o == '-c' or o == '--color':
                    self.color = a
                elif o == '--show-graph':
                    self.showGraphs = True
                elif o == '--save-graph':
                    self.saveGraphs = True

            if args:
                self.folder = args[0]
                if not self.folder.endswith('/'):
                    self.folder += '/'
            else:
                pm.printError('No folder specified in the command-line arguments.')
                return

        def run(self):
            self.processCommands()

            colors = ['Ri', 'Gi', 'Bi'] if self.color is None else [self.color]

            matcher = CatalogMatcher('Combined', self.folder, colors, self.showGraphs, self.saveGraphs,
                                     pm.Logger.LOG_MODE_DEBUG)
            matcher.process()


    app = MainApp(argv)
    app.run()

# end main.
