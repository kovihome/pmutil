#!/usr/bin/env python3
#
# PmUtils/pmphot
#
"""
Created on Dec 28, 2019

@author: kovi
"""
from getopt import getopt, GetoptError
from sys import argv
from os.path import exists
from copy import deepcopy

import numpy as np
from astropy.table import Table, Column

import pmbase as pm

fieldMgInstrumental = 'MAG_BEST'
fieldMgTrue = 'MAG_V'
fieldMgErrInstrumental = 'MAGERR_BEST'
fieldMgErrTrue = 'ERR_V'
fieldAuid = 'AUID'
fieldRole = 'ROLE'

MAG_ERR_DEFAULT = 0.15


class StdCoeffs:
    coeffs = None  # astropy.table.Table of std coeffs
    fileName = None  # coeff fiel name
    observer = None  # observer name
    date = None  # date of the current image
    camera = None  # camera using to crate current image
    telescope = None  # telescope using to create current image
    std_area = None  # standard area or other field name

    fields = [
        {'name': 'OBSERVER', 'format': '%-10s', 'dtype': 'U5'},
        {'name': 'DATE', 'format': '%-14s', 'dtype': 'U10'},
        {'name': 'TV', 'format': '%11.4f', 'dtype': 'f4'},
        {'name': 'ERR_TV', 'format': '%11.4f', 'dtype': 'f4'},
        {'name': 'TVR', 'format': '%11.4f', 'dtype': 'f4'},
        {'name': 'ERR_TVR', 'format': '%11.4f', 'dtype': 'f4'},
        {'name': 'TBV', 'format': '%11.4f', 'dtype': 'f4'},
        {'name': 'ERR_TBV', 'format': '%11.4f', 'dtype': 'f4'},
        {'name': 'CAMERA', 'format': '%-24s', 'dtype': 'U24'},
        {'name': 'TELESCOPE', 'format': '%-24s', 'dtype': 'U24'},
        {'name': 'FIELD', 'format': '%-24s', 'dtype': 'U24'},
    ]

    def __init__(self, coeffFileName, observer, date, camera, telescope, std_area):
        self.fileName = coeffFileName
        self.observer = observer
        self.date = date
        self.camera = camera
        self.telescope = telescope
        self.std_area = std_area

    def open(self):
        if not exists(self.fileName):
            self.coeffs = Table()
            for field in self.fields:
                self.coeffs.add_column(Column(name=field['name'], dtype=field['dtype'], format=field['format']))

        else:
            self.coeffs = Table.read(self.fileName, format='ascii')
            for field in self.fields:
                self.coeffs[field['name']].format = field['format']

    def addCoeffs(self, Tv, Tvr, Tbv, errTv=None, errTvr=None, errTbv=None, field=None):
        lastCoeffs = self.getCoeffs()
        if not lastCoeffs:
            self.coeffs.add_row((self.observer, self.date, Tv, errTv, Tvr, errTvr, Tbv, errTbv, self.camera,
                                 self.telescope, self.std_area))
        else:
            lastCoeffs['TV'] = Tv
            lastCoeffs['ERR_TV'] = errTv
            lastCoeffs['TVR'] = Tvr
            lastCoeffs['ERR_TVR'] = errTvr
            lastCoeffs['TBV'] = Tbv
            lastCoeffs['ERR_TBV'] = errTbv

    #            lastCoeffs['FIELD'] = field

    def getCoeffs(self):
        flt = (self.coeffs['OBSERVER'] == self.observer) & (self.coeffs['DATE'] == self.date) & (
                self.coeffs['CAMERA'] == self.camera) & (self.coeffs['TELESCOPE'] == self.telescope)
        rows = self.coeffs[flt]
        if len(rows) == 0:
            return None
        else:
            return rows[0]

    def exists(self):
        return self.getCoeffs() is not None

    def getBestCoeffs(self):
        flt = (self.coeffs['OBSERVER'] == self.observer) & (self.coeffs['CAMERA'] == self.camera) & (
                self.coeffs['TELESCOPE'] == self.telescope)
        rows = self.coeffs[flt]
        if len(rows) == 0:
            return None
        elif len(rows) == 1:
            return rows[0]
        else:
            jdnow = pm.jd(self.date)
            delta = abs(jdnow - pm.jd(rows[0]['DATE']))
            ix = 0
            for j in range(len(rows) - 1):
                d = abs(jdnow - pm.jd(rows[j + 1]['DATE']))
                if d < delta:
                    delta = d
                    ix = j + 1
            return rows[ix]

    def getAvgCoeffs(self):
        flt = (self.coeffs['OBSERVER'] == self.observer) & (self.coeffs['CAMERA'] == self.camera) & (
                self.coeffs['TELESCOPE'] == self.telescope)
        rows = self.coeffs[flt]
        if len(rows) == 0:
            return None
        avgc = deepcopy(rows[0])
        avgc['TV'] = np.average(rows['TV'])
        avgc['TVR'] = np.average(rows['TVR'])
        avgc['TBV'] = np.average(rows['TBV'])

        # TODO: calculate average errors

        return avgc

    def save(self):
        self.coeffs.write(self.fileName, format='ascii.fixed_width_two_line', delimiter=' ')


class Photometry:
    opt = {}  # command line options
    pos = None  # catalog header positions
    imageProps = {}  # image properties from FITS header
    coeffTable = None
    fits = {}

    mif = {'Gi': 'MAG_GI', 'Bi': 'MAG_BI', 'Ri': 'MAG_RI'}
    eif = {'Gi': 'ERR_GI', 'Bi': 'ERR_BI', 'Ri': 'ERR_RI'}
    mvf = {'Gi': 'MAG_V', 'Bi': 'MAG_B', 'Ri': 'MAG_R'}
    evf = {'Gi': 'ERR_V', 'Bi': 'ERR_B', 'Ri': 'ERR_R'}

    def __init__(self, opt):
        self.opt = opt

    def loadFitsHeaders(self, fileName):
        """
        fileName: path/[Seq_nnn|Combined].cmb.pm
        """
        self.fits = pm.getFitsHeaders(fileName, ['DATE-OBS', 'INSTRUME', 'TELESCOP'])

    def getCamera(self):
        if 'INSTRUME' in self.fits.keys() and self.fits['INSTRUME'] and len(self.fits['INSTRUME']) > 0:
            return self.fits['INSTRUME']
        if 'camera' in self.opt and self.opt['camera'] and len(self.opt['camera']) > 0:
            return self.opt['camera']
        if 'DEF_CAMERA' in pm.setup.keys() and pm.setup['DEF_CAMERA'] and len(pm.setup['DEF_CAMERA']) > 0:
            return pm.setup['DEF_CAMERA']
        return 'Generic Camera'

    def getTelescope(self):
        if 'TELESCOP' in self.fits.keys() and self.fits['TELESCOP'] and len(self.fits['TELESCOP']) > 0:
            return self.fits['TELESCOP']
        if 'telescope' in self.opt and self.opt['telescope'] and len(self.opt['telescope']) > 0:
            return self.opt['telescope']
        if 'DEF_TELESCOPE' in pm.setup.keys() and pm.setup['DEF_TELESCOPE'] and len(pm.setup['DEF_TELESCOPE']) > 0:
            return pm.setup['DEF_TELESCOPE']
        return 'Generic Telescope'

    def stdColor(self, color):
        stdcolor = color[:1].upper()
        if stdcolor == 'G':
            stdcolor = 'V'
        return stdcolor

    def calculateMgsRobustAveraging(self, comps, color, bestCompId):
        """
        Calculate unit slope linear fit regression with robust averaging
        see: http://gcx.sourceforge.net/html/node11.html
        inputs:    I[k] - instrumental mgs              MAG_ISOCORR | MAG_BEST
                   S[k] - standard mgs                  MAG_R \ MAG_V | MAG_B
                   ei[k] - error of instrumental mgs    MAGERR_ISOCORR | MAGERR_BEST
                   es[k] - error of standard mgs        ERR_R | ERR_V | ERR_B
        """

        stdcolor = self.stdColor(color)

        # bestComp = comps.loc['AUID', bestCompId]

        y = []
        e2 = []
        w = []
        mis = []
        mvs = []
        for r in comps:
            if r[self.mvf[color]] != '-' and r['VIZ_FLAG'] == '---':
                mi = float(r[self.mif[color]])
                mv = float(r[self.mvf[color]])
                mis.append(mi)
                mvs.append(mv)
                ei = float(r[self.eif[color]])
                ev = float(r[self.evf[color]]) if r[self.evf[color]][0].isdigit() and r[
                    self.evf[color]] != '0.0' else ei
                y.append(mv - mi)
                ek2 = ei * ei + ev * ev
                e2.append(ek2)
                w.append(1.0 / ek2)

        z = np.median(y)
        r = [0.0] * len(y)
        # r_ = [0.0] * len(y)
        w_ = [1.0] * len(y)

        # iteration starts
        for it in range(10):

            zsum = 0.0
            wsum = 0.0
            beta = 1.0
            alpha = 1.0
            for k in range(len(y)):
                r[k] = y[k] - z
                w_k = w[k] / (1 + (r[k] / alpha) ** beta)
                w_[k] = w_k
                zsum = zsum + y[k] * w_k
                wsum = wsum + w_k
            z = zsum / wsum
        # iteration ends

        su = 0.0
        sl = 0.0
        se2 = 0.0
        for k in range(len(y)):
            su = su + r[k] * r[k] * w_[k]
            sl = sl + w_[k]
            se2 = se2 + 1.0 / (w_[k] * w_[k])
        ez_2 = su / sl
        mel_2 = su / float(len(y) - 1)
        err = np.sqrt(se2 / len(y))
        pm.printDebug(
                f"(gcx) color: {color}, zp: {z:7.4f}, mel^2: {mel_2:7.4f}, ez^2: {ez_2:7.4f}, ez: {np.sqrt(ez_2):7.4f}, N: {len(y):d}, ev: {err:7.4f}")

        self.ePlot.add(mis, mvs, [1, z], color, stdcolor, color[0].lower(), pm.Plot.INV_X + pm.Plot.INV_Y,
                       title=f"Inst {color} vs. Cat {stdcolor}")

        return [1.0, z], err

    def calculateMgsLinearFit(self, comps, color, bestCompId):
        """
        Calculate magnitudes with naive linear fit
        """
        # TODO: calculate comp mag error

        stdcolor = self.stdColor(color)

        # bestComp = cmb.loc['AUID', bestCompId]

        # TODO
        mi = []
        mv = []
        er = []
        #    p = [1.0, 0.0]
        ep = 0.0
        N = 0
        for pmr in comps:
            if pmr[self.mvf[color]] != '-' and pmr['VIZ_FLAG'] == '---':
                mi.append(float(pmr[self.mif[color]]))
                mv.append(float(pmr[self.mvf[color]]))
                ei = float(pmr[self.eif[color]])
                ev = float(pmr[self.evf[color]]) if pmr[self.evf[color]] != '-' and pmr[
                    self.evf[color]] != '0.0' else ei
                e2 = ei * ei + ev * ev
                ep = ep + e2
                N = N + 1
                er.append(1.0 / e2)

        if len(mi) < 2:
            pm.printError("Not enough comp stars for linear fit ensemble")
            return []

        #        p = Polynomial.fit(mi, mv, 1, w = er)
        coef = pm.linfit(mi, mv, er)
        ep = np.sqrt(ep / float(N))

        pm.printDebug(f"polyfit result: {coef} error: {ep}")

        self.ePlot.add(mi, mv, coef, color, stdcolor, color[0].lower(), title=f"Inst {color} vs. Cat {stdcolor}")

        return coef, ep

    def findBestCompStar(self, comps, color):
        # filter out com stars not measured
        emin = 1.0
        bestComp = None
        for pms in comps:
            mvs = pms[self.mvf[color]]
            evs = pms[self.evf[color]]
            if mvs != '-' and mvs != '99.0':
                ei = float(pms[self.eif[color]]) if pms[self.eif[color]] != '-' else MAG_ERR_DEFAULT
                ev = float(evs) if evs != '-' and evs != '0.0' else ei
                e = ei * ei + ev * ev
                if e < emin:
                    haveAllMags = True
                    for c in self.opt['color']:
                        if c != color:
                            if pms[self.mvf[c]] == '-' or pms[self.evf[c]] == '-':
                                haveAllMags = False
                                break
                    if not haveAllMags:
                        continue
                    emin = e
                    bestComp = pms

        if bestComp:
            ei = float(bestComp[self.eif[color]]) if bestComp[self.eif[color]] != '-' else MAG_ERR_DEFAULT
            ev = float(bestComp[self.evf[color]]) if bestComp[self.evf[color]] != '-' else MAG_ERR_DEFAULT
            e = pm.quad(ei, ev)
            bestComp['ROLE'] = 'K'
            pm.printDebug(f'Best comp star: {bestComp["AUID"]}, mag: {bestComp[self.mvf[color]]}, err: {e:4.3f}')
            return bestComp['AUID']
        else:
            pm.printError('No usable comp star found ; check the comp stars if they exist in all colors you need')
            return None

    def calculateMgsComparision(self, comps, color, bestCompId):
        stdcolor = self.stdColor(color)

        # find best comp star
        bestComp = comps.loc['AUID', bestCompId]

        zp = float(bestComp[self.mvf[color]]) - float(bestComp[self.mif[color]])
        p = [1.0, zp]
        ei = float(bestComp[self.eif[color]]) if bestComp[self.eif[color]] != '-' else MAG_ERR_DEFAULT
        ev = float(bestComp[self.evf[color]]) if bestComp[self.evf[color]] != '-' and bestComp[
            self.evf[color]] != '0.0' else ei
        ep = pm.quad(ei, ev)

        if self.opt['showGraphs'] or self.opt['saveGraphs']:
            mi = []
            mv = []
            for pmr in comps:
                if pmr[self.mvf[color]] != '-' and pmr['VIZ_FLAG'] == '---':
                    mi.append(float(pmr[self.mif[color]]))
                    mv.append(float(pmr[self.mvf[color]]))

            self.ePlot.add(mi, mv, [1, zp], color, stdcolor, color[0].lower(), title=f"Inst {color} vs. Cat {stdcolor}")

        return p, ep

    def reportResult(self, cmb, outFileName, dateObs, coeffs):

        jdObs = "%12.4f" % pm.jd(dateObs)

        dateObsCol = [dateObs] * len(cmb)
        jdObsCol = [jdObs] * len(cmb)
        cmb['DATE-OBS'] = dateObsCol
        cmb['JD'] = jdObsCol

        # put std coeffs into table comments
        pm.addTableComment(cmb, 'StdTv', str(coeffs[0]))
        pm.addTableComment(cmb, 'StdTbv', str(coeffs[1]))
        pm.addTableComment(cmb, 'StdTvr', str(coeffs[2]))

        cmb.write(outFileName, format='ascii.fixed_width', delimiter=' ', overwrite=True)

        return

    # standardization functions

    def mgd(self, s):
        return float(s) if s != '-' and s != '0.0' else MAG_ERR_DEFAULT

    def calculateCoeffs(self, mergedCat):
        # Tv:  slope of (V-v) -> (V-R)
        # Tvr: 1/slope of (v-r) -> (V-R)
        # Tbv: 1/slope of (b-v) -> (B-V)
        V_vi = []
        V_R = []
        vi_ri = []
        bi_vi = []
        B_V = []
        w_Tv = []
        w_Tvr = []
        w_Tbv = []
        for row in mergedCat:
            if row['MAG_V'] == '-' or row['MAG_B'] == '-' or row['MAG_R'] == '-' or \
                    row['MAG_GI'] == '-' or row['MAG_BI'] == '-' or row['MAG_RI'] == '-':
                continue

            V_vi.append(float(row['MAG_V']) - float(row['MAG_GI']))
            V_R.append(float(row['MAG_V']) - float(row['MAG_R']))
            vi_ri.append(float(row['MAG_GI']) - float(row['MAG_RI']))
            bi_vi.append(float(row['MAG_BI']) - float(row['MAG_GI']))
            B_V.append(float(row['MAG_B']) - float(row['MAG_V']))

            e_V_vi = self.mgd(row['ERR_V']) + self.mgd(row['ERR_GI'])
            e_V_R = self.mgd(row['ERR_V']) + self.mgd(row['ERR_R'])
            e_vi_ri = self.mgd(row['ERR_GI']) + self.mgd(row['ERR_RI'])
            e_bi_vi = self.mgd(row['ERR_BI']) + self.mgd(row['ERR_GI'])
            e_B_V = self.mgd(row['ERR_B']) + self.mgd(row['ERR_V'])
            w_Tv.append(1.0 / (e_V_vi * e_V_vi + e_V_R * e_V_R))
            w_Tvr.append(1.0 / (e_vi_ri * e_vi_ri + e_V_R * e_V_R))
            w_Tbv.append(1.0 / (e_bi_vi * e_bi_vi + e_B_V * e_B_V))

        coef = pm.linfit(V_R, V_vi, w_Tv)
        Tv = coef[0]
        Bv = coef[1]  # ##
        coef = pm.linfit(V_R, vi_ri, w_Tvr)
        Tvr = 1.0 / coef[0]
        Bvr = coef[1]  # ##
        coef = pm.linfit(B_V, bi_vi, w_Tbv)
        Tbv = 1.0 / coef[0]
        Bbv = coef[1]  # ##

        # show/save graph of std coeffs
        stdPlot = pm.Plot(3, self.opt['showGraphs'], self.opt['saveGraphs'])
        stdPlot.add(V_R, V_vi, [Tv, Bv], 'V-R', 'V-v', 'g', title="V-R vs. V-v")
        stdPlot.add(V_R, vi_ri, [1.0 / Tvr, Bvr], 'V-R', 'v-r', 'r', title="V-R vs. v-r")
        stdPlot.add(B_V, bi_vi, [1.0 / Tbv, Bbv], 'B-V', 'b-v', 'b', title="B-V vs. b-v")
        folder, a, b = self.opt['files'][0].partition('Photometry')
        stdPlot.showOrSave(folder + 'std_coefficients.png')

        pm.printDebug(f"Standard coeffs: Tv = {Tv} Tvr = {Tvr} Tbv = {Tbv}")
        return [Tv, Tvr, Tbv]

    def openCoeffs(self, obsDate, target):
        if not self.coeffTable:
            configFolder = pm.setup['CONFIG_FOLDER'].rstrip('/')
            coeffFile = configFolder + '/stdcoeffs.cat'
            self.coeffTable = StdCoeffs(coeffFile, self.opt['observer'], obsDate, self.getCamera().replace(' ', '_'),
                                        self.getTelescope().replace(' ', '_'), target.replace(' ', '_').upper())
            self.coeffTable.open()

    def saveCoeffs(self, coeffs, obsDate, target):
        self.openCoeffs(obsDate, target)
        if not self.opt['overwrite'] and self.coeffTable.exists():
            pm.printError(
                f'Standard coefficients for {self.coeffTable.observer}/{self.coeffTable.date}/{self.coeffTable.camera}/{self.coeffTable.telescope} is exists ; to overwrite it use -w option')
            return
        self.coeffTable.addCoeffs(coeffs[0], coeffs[1], coeffs[2], 0.0, 0.0, 0.0, target.replace(' ', '_').upper())
        #        self.coeffTable.save(overwrite=True)
        self.coeffTable.save()

    def loadCoeffs(self, obsDate, target):
        self.openCoeffs(obsDate, target)
        bestCoeffs = self.coeffTable.getBestCoeffs()
        if not bestCoeffs:
            pm.printWarning('standard coeffs not found')
            return None
        pm.printDebug(f"standard coeffs: Tv = {bestCoeffs['TV']} Tvr = {bestCoeffs['TVR']} Tbv = {bestCoeffs['TBV']}")
        return [bestCoeffs['TV'], bestCoeffs['TVR'], bestCoeffs['TBV']]

    def calcVComp(self, colors, cmb, bestComp, p_a, err_a):
        vcomp = {}
        c_row = cmb.loc['AUID', bestComp]
        for c in colors:
            vcomp[c] = {}
            vcomp[c]['mi'] = float(c_row[self.mif[c]])
            vcomp[c]['ei'] = float(c_row[self.eif[c]])
            vcomp[c]['mc'] = p_a[c][0] * vcomp[c]['mi'] + p_a[c][1]
            vcomp[c]['ec'] = pm.quad(vcomp[c]['ei'], err_a[c]) if err_a[c] else vcomp[c]['ei'] * np.sqrt(2.0)

        return vcomp

    def setEmptyStd(self, cmb):
        cmb['MAG_STDB'] = ['-'] * len(cmb)
        cmb['ERR_STDB'] = ['-'] * len(cmb)
        cmb['MAG_STDV'] = ['-'] * len(cmb)
        cmb['ERR_STDV'] = ['-'] * len(cmb)
        cmb['MAG_STDR'] = ['-'] * len(cmb)
        cmb['ERR_STDR'] = ['-'] * len(cmb)
        return cmb

    def calculateSimpleMgs(self, cmb, vcomp):
        ma = []
        ea = []
        for c in ['Ri', 'Gi', 'Bi']:
            cc = c[:1]
            for row in cmb:
                isFainter = row['VIZ_FLAG'] == 'III'  # TODO: I es B flagek kezelese is

                m0 = float(row['MAG_' + c])
                M = m0 + vcomp[c]['mc'] - vcomp[c]['mi']
                ma.append(M)

                err_m0 = float(row['ERR_' + c])
                errM = np.sqrt(err_m0 * err_m0 + vcomp[c]['ec'] * vcomp[c]['ec'] + vcomp[c]['ei'] * vcomp[c]['ei'])
                ea.append(errM)

            cmb['MAG_T' + cc] = ma
            cmb['ERR_T' + cc] = ea

        return self.setEmptyStd(cmb)

    def calculateStdMgs(self, coeffs, cmb, vcomp, std):
        """
        allResults: dict of results, key is th color
        """
        vc = vcomp['Gi']['mi']
        bc = vcomp['Bi']['mi']
        rc = vcomp['Ri']['mi']
        Vc = vcomp['Gi']['mc']
        Bc = vcomp['Bi']['mc']
        Rc = vcomp['Ri']['mc']
        Vc_Rc = Vc - Rc

        err_vc = vcomp['Gi']['ei']
        err_bc = vcomp['Bi']['ei']
        err_rc = vcomp['Ri']['ei']
        err_Vc = vcomp['Gi']['ec']
        err_Bc = vcomp['Bi']['ec']
        err_Rc = vcomp['Ri']['ec']

        v_result = []
        b_result = []
        r_result = []
        ev_result = []
        eb_result = []
        er_result = []

        for g_row in cmb:
            role = g_row['ROLE']
            if 1 == 1:  # role == 'V' or role == 'VF' or role == 'K':
                # auid = g_row['AUID']

                isFainter = g_row['VIZ_FLAG'] == 'III'  # TODO: I es B flagek kezelese is

                Tv = coeffs[0] if not isFainter else 0.0
                Tvr = coeffs[1] if not isFainter else 1.0
                Tbv = coeffs[2] if not isFainter else 1.0

                v0 = float(g_row['MAG_GI'])
                b0 = float(g_row['MAG_BI'])
                r0 = float(g_row['MAG_RI'])

                V_R = (Vc_Rc) + Tvr * ((v0 - r0) - (vc - rc))
                V = v0 + (Vc - vc) + Tv * ((V_R) - (Vc_Rc))
                R = V - (V_R)
                B = V + (Bc - Vc) + Tbv * ((b0 - v0) - (bc - vc))

                v_result.append(V)
                b_result.append(B)
                r_result.append(R)

                if not isFainter:

                    err_v0 = float(g_row['ERR_GI']) \
                        if g_row['ERR_GI'] != '-' and g_row['ERR_GI'] != '0.0' \
                        else err_vc
                    err_b0 = float(g_row['ERR_BI']) \
                        if g_row['ERR_BI'] != '-' and g_row['ERR_BI'] != '0.0' \
                        else err_bc
                    err_r0 = float(g_row['ERR_RI']) \
                        if g_row['ERR_RI'] != '-' and g_row['ERR_RI'] != '0.0' \
                        else err_rc

                    errV = np.sqrt(err_v0 * err_v0 + err_Vc * err_Vc + err_vc * err_vc)
                    errB = np.sqrt(err_b0 * err_b0 + err_Bc * err_Bc + err_bc * err_bc)
                    errR = np.sqrt(err_r0 * err_r0 + err_Rc * err_Rc + err_rc * err_rc)

                    ev_result.append(errV)
                    eb_result.append(errB)
                    er_result.append(errR)

                else:

                    ev_result.append(-1.0)
                    eb_result.append(-1.0)
                    er_result.append(-1.0)

            else:

                v_result.append(99.0)
                b_result.append(99.0)
                r_result.append(99.0)
                ev_result.append(-1.0)
                eb_result.append(-1.0)
                er_result.append(-1.0)

        if std:

            cmb['MAG_STDB'] = b_result
            cmb['ERR_STDB'] = eb_result
            cmb['MAG_STDV'] = v_result
            cmb['ERR_STDV'] = ev_result
            cmb['MAG_STDR'] = r_result
            cmb['ERR_STDR'] = er_result

        else:

            cmb['MAG_TB'] = b_result
            cmb['ERR_TB'] = eb_result
            cmb['MAG_TG'] = v_result
            cmb['ERR_TG'] = ev_result
            cmb['MAG_TR'] = r_result
            cmb['ERR_TR'] = er_result

        return cmb

    def calcHmg(self, cmb, vcomp):
        # print(f'vcomp: {vcomp}')
        for cc in ['Gi', 'Bi', 'Ri']:
            hmgInst = pm.getTableComment(cmb, 'MgLimitInst' + cc)
            hmg = float(hmgInst) + vcomp[cc]['mc'] - vcomp[cc]['mi']
            pm.addTableComment(cmb, 'MgLimit' + cc, '%7.3f' % hmg)

    def calculateEnsembleParams(self, cmb, bestComp):

        self.ePlot = pm.Plot(3, self.opt['showGraphs'], self.opt['saveGraphs'])

        p_a = {}
        err_a = {}
        for color in self.opt['color']:

            if self.opt['method'] == 'comp':
                p, ep = self.calculateMgsComparision(cmb, color, bestComp)
            elif self.opt['method'] == 'lfit':
                p, ep = self.calculateMgsLinearFit(cmb, color, bestComp)
            else:
                p, ep = self.calculateMgsRobustAveraging(cmb, color, bestComp)

            p_a[color] = p
            err_a[color] = ep

        folder, a, b = self.opt['files'][0].partition('Photometry')
        self.ePlot.showOrSave(folder + 'ensemble_parameters.png')

        return p_a, err_a

    def process(self):
        if not self.opt['files']:
            pm.printError('No input files are given.')
            return

        # load cmb catalog
        cmbFileName = self.opt['files'][0]
        pm.printDebug(f' Input files: {cmbFileName}')
        cmb = Table.read(cmbFileName, format='ascii')
        cmb.add_index('AUID')

        # filter comp stars
        compStarMask = cmb['ROLE'] == 'C'
        compStars = cmb[compStarMask]

        # find best comp star in Gi frame
        color = 'Gi' if 'Gi' in self.opt['color'] else self.opt['color'][0]
        fitsFileName = cmbFileName.replace('.cmb', f'-{color}.ast.fits')
        pm.printDebug(f'.ast.fits file name: {fitsFileName}')
        self.loadFitsHeaders(fitsFileName)
        pm.printDebug(f' Headers: {self.fits}')
        dateObs = self.fits['DATE-OBS']
        dateObsDate = dateObs.split('T')[0]

        bestComp = self.findBestCompStar(compStars, color)  # refcat record of best comp in Gi
        if not bestComp:
            return False

        # calculate virtual comp star
        p_a, err_a = self.calculateEnsembleParams(compStars, bestComp)
        vcomp = self.calcVComp(self.opt['color'], compStars, bestComp, p_a, err_a)

        # calculate or load std coeffs
        std = True
        if self.opt['useCoeffs'] or self.opt['makeCoeffs']:
            # standardization
            coeffs = None

            print("************** files")
            print(self.opt['files'])
            target = pm.guess(self.opt['files'][0])['target']

            if self.opt['makeCoeffs']:

                coeffs = self.calculateCoeffs(compStars)

                if self.opt['saveCoeffs']:
                    self.saveCoeffs(coeffs, dateObsDate, target)

            elif self.opt['loadCoeffs']:

                coeffs = self.loadCoeffs(dateObsDate, target)

            if not coeffs:
                pm.printError(
                        'No std coefficients for transformation ; use Tv = 0, Tvr = 1, Tbv = 1 for comp star method.')
                coeffs = [0.0, 1.0, 1.0]  # this is for non-std transformation
                std = False

        else:

            coeffs = [0.0, 1.0, 1.0]
            std = False

        # apply virtual comp and std coeffs to calculate magnitudes
        if len(self.opt['color']) == 3:
            cmb = self.calculateStdMgs([0.0, 1.0, 1.0], cmb, vcomp, False)
            if std:
                cmb = self.calculateStdMgs(coeffs, cmb, vcomp, True)
            else:
                cmb = self.setEmptyStd(cmb)
        else:
            cmb = self.calculateSimpleMgs(cmb, vcomp)

        # calculate real mg limit
        self.calcHmg(cmb, vcomp)

        # save compStar auid
        if bestComp:
            pm.addTableComment(cmb, 'CompStar', bestComp)

        # save results
        fileName = self.opt['files'][0]
        self.reportResult(cmb, fileName + '.pm', dateObs, coeffs)

        return True


if __name__ == '__main__':

    class MainApp:

        opt = {
            'out': None,  # output photometry file name, if None, '*.pm' will be created
            'method': 'gcx',
            # mg calculation method: comp - use best comp star, gcx - m=1 linear fit ensemble, lfit - general linear fit ensemble
            'color': ['Gi'],  # photometry bands
            'loadCoeffs': False,  # load std coeffs
            'useCoeffs': False,  # use std coeffs
            'makeCoeffs': False,  # create std coeffs
            'saveCoeffs': False,  # save std coeffs
            'showGraphs': False,  # show standard coefficient graphs
            'saveGraphs': False,  # save standard coefficient graphs
            'camera': None,  # camera name
            'telescope': None,  # telescope name
            'files': None,

            # from pplphotomety
            'makeStd': False,  # make std coeffs
            'useStd': False,  # use std coeffs
            'adhocStd': False,  # make and use std coeffs for this image only
            'overwrite': False,  # force to overwrite existing results, optional
            'baseFolder': None,  # base folder, optional

        }

        availableBands = ['gi', 'g', 'bi', 'b', 'ri', 'r', 'all']

        mgCalcMethods = {
            'comp': 'Best comparision star',
            'gcx': 'GCX''s robust averaging ensemble',
            'lfit': 'Linear fit ensemble',
        }

        def __init__(self, argv):
            self.argv = argv
            pass

        def processCommands(self):
            try:
                optlist, args = getopt(self.argv[1:], "c:msat:w",
                                       ['color=', 'make-std', 'use-std', 'adhoc-std', 'show-graph', 'save-graph',
                                        'method=', 'overwrite'])
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
                        self.opt['color'] = ['Ri', 'Gi', 'Bi']
                    else:
                        self.opt['color'] = [a]
                elif o == '-m' or o == '--make-std':
                    self.opt['makeStd'] = True
                elif o == '-s' or o == '--use-std':
                    self.opt['useStd'] = True
                elif o == '-a' or o == '--adhoc-std':
                    self.opt['adhocStd'] = True
                elif o == '--show-graph':
                    self.opt['showGraphs'] = True
                elif o == '--save-graph':
                    self.opt['saveGraphs'] = True
                elif o == '--camera':
                    self.opt['camera'] = a
                elif o == '--telescope':
                    self.opt['telescope'] = a
                elif o == '-t' or o == '--method':
                    if a not in self.mgCalcMethods.keys():
                        pm.printWarning(f'Invalid mg calculation method {a} ; use gcx instead.')
                    else:
                        self.opt['method'] = a

                elif o == '-w' or o == '--overwrite':
                    self.opt['overwrite'] = True

            if self.opt['adhocStd'] and (self.opt['makeStd'] or self.opt['useStd']):
                pm.printWarning('Both -a and either -m or -s option cannot be used at once.')
                exit(0)

            if len(args) > 0:
                self.opt['baseFolder'] = args[0]
                if args[0].endswith('/'):
                    self.opt['baseFolder'] = args[0][:-1]

            pm.printInfo('Mg calculation method: ' + self.mgCalcMethods[self.opt['method']])

            self.opt['loadCoeffs'] = self.opt['useStd']
            self.opt['useCoeffs'] = self.opt['useStd'] or self.opt['adhocStd']
            self.opt['makeCoeffs'] = self.opt['makeStd'] or self.opt['adhocStd']
            self.opt['saveCoeffs'] = self.opt['makeStd']

        def run(self):
            self.processCommands()

            self.opt['files'] = [self.opt['baseFolder'] + '/Photometry/Combined.cmb']

            phot = Photometry(self.opt)
            phot.process()


    app = MainApp(argv)
    app.run()

# end __main__
