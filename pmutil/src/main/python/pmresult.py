#!/usr/bin/env python3
#
# PmUtils/pmresult
#
"""
Created on Apr 28, 2019
@author: kovi
"""

import getopt
import glob
import sys
from os.path import exists

import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
from astropy.table import Table
from astropy.time import Time
from matplotlib import pyplot as plt

import pmbase as pm

REPORT_TYPE_AAVSO = "aavso"
REPORT_TYPE_VSNET = "vsnet"

# manage catalog file
def loadCatalog(refFileName):
    table = Table.read(refFileName, format='ascii')
    if not table:
        print(f"Reference catalog {refFileName} not found")
        return None
    return table


# create report
class ReportProcessor:
    NOT_AVAILABLE = 'na'

    opt = {}  # command line options
    pos = None  # catalog header positions
    airmasses = {}  # airmass cache
    loc = None

    lcdata = {}

    chartId = NOT_AVAILABLE

    aavso_colors = {
        'Gi': 'TG',
        'Ri': 'TR',
        'Bi': 'TB'
    }
    aavso_std_colors = {
        'Gi': 'V',
        'Ri': 'R',
        'Bi': 'B'
    }

    def __init__(self, opt):
        self.mgLimits = None
        self.airmasses = {}
        self.opt = opt

    def findChartId(self, baseFolder):
        catfile = baseFolder + 'ref.cat'
        if not exists(catfile):
            catfiles = glob.glob(baseFolder + '*.cat')
            if len(catfiles) < 1:
                return self.NOT_AVAILABLE
            catfile = catfiles[0]

        cf = open(catfile)
        chartId = self.NOT_AVAILABLE
        for line in cf:
            if line.find('### chartid:') > -1:
                r = line.split()
                chartId = r[2].strip()
                break
        cf.close()
        return chartId

    def getAirmass(self, ra, dec, obsDate):
        key = f"{ra}{dec}{obsDate}"
        if key in self.airmasses:
            return self.airmasses[key]
        else:
            airmass = self.calcAirmass(ra, dec, obsDate)
            self.airmasses[key] = airmass
            return airmass

    def calcAirmass(self, ra, dec, obsDate):
        co = SkyCoord(ra, dec, frame=ICRS, unit=(u.hourangle, u.deg))
        t = Time(obsDate.replace('T', ' '))
        # TODO: store observer's location in config
        if self.loc is None:
            self.loc = EarthLocation.of_address('Budapest')
        hc = co.transform_to(AltAz(obstime=t, location=self.loc))
        return float(hc.secz)

    def convertToAAVSOData(self, res, color, compStar):

        # MAG_Tc, ERR_Tc, MAG_STDc, ERR_STDc

        cc = color[:1]
        sc = self.aavso_std_colors[color]

        starid = res['LABEL'].replace('_', ' ')
        date = res['JD']

        visFlag = res['VIZ_FLAG']['GBR'.index(cc)]
        fainter = "<" if visFlag in "IB" else ""

        if fainter == '<':
            magnitude = self.mgLimits[color]
            magerr = self.NOT_AVAILABLE
            aavsoColor = self.aavso_colors[color]
        else:
            # std = True
            m = res['MAG_STD' + sc]
            e = res['ERR_STD' + sc]
            aavsoColor = self.aavso_std_colors[color]
            if str(m) == '-' or str(m) == '99.0' or m == 99.0:
                m = res['MAG_T' + cc]
                e = res['ERR_T' + cc]
                aavsoColor = self.aavso_colors[color]
                # std = False
            if str(m) == '-':
                return None

            magnitude = f"{float(m):4.3f}"
            magerr = f"{float(e):4.3f}" if e != '-' else self.NOT_AVAILABLE

        flt = aavsoColor
        trans = 'NO'
        mtype = 'STD'

        if self.opt['method'] == 'comp':
            cname = compStar['AUID'] if compStar else self.NOT_AVAILABLE
            cmag = compStar['MAG_' + sc] if compStar else self.NOT_AVAILABLE
            kname = self.NOT_AVAILABLE  # TODO: set a check star
            kmag = self.NOT_AVAILABLE  # TODO: set measured mg of check star
        else:
            cname = 'ENSEMBLE'
            cmag = self.NOT_AVAILABLE
            kname = compStar['AUID'] if compStar else self.NOT_AVAILABLE
            kmag = compStar['MAG_' + sc] if compStar else self.NOT_AVAILABLE
        airmass = f"{self.getAirmass(res['RA'], res['DEC'], res['DATE-OBS']):6.4f}"
        group = self.NOT_AVAILABLE
        chart = self.chartId
        notes = self.NOT_AVAILABLE

        return (starid.upper(), date, fainter, magnitude, magerr, flt, trans, mtype, cname, cmag,
                kname, kmag, airmass, group, chart, notes)

    @staticmethod
    def convertToAAVSORecord(data):
        # (starid.upper(), date, fainter, magnitude, magerr, flt, trans, mtype,
        #  cname, cmag, kname, kmag, airmass, group, chart, notes)
        return f'{data[0]},{data[1]},{data[2]}{data[3]},{data[4]},{data[5]},{data[6]},{data[7]},{data[8]},{data[9]},{data[10]},{data[11]},{data[12]},{data[13]},{data[14]},{data[15]}\n'

    def convertToAAVSORecords(self, pmt, colors):
        # self.tick()
        compStar = self.findCompStar(pmt)

        varMask = pmt['ROLE'] == 'V'
        vt = pmt[varMask]
        aavsoResults = []
        for v in vt:
            for c in colors:
                aavsoData = self.convertToAAVSOData(v, c, compStar)
                if aavsoData:
                    aavso = self.convertToAAVSORecord(aavsoData)
                    aavsoResults.append(aavso)
                    self.addToLCData(aavsoData)
        return aavsoResults

    @staticmethod
    def saveAAVSOReport(aavso, outFolder, obsName, postfix=""):
        r = outFolder.split('/')
        fileName = r[-1]
        outFileName = outFolder + '/' + fileName + postfix + '.extended.aavso'
        pm.printInfo(f"Save AAVSO extended report to the file {outFileName}")

        with open(outFileName, 'w') as f:
            f.write('#TYPE=Extended\n')
            f.write(f'#OBSCODE={obsName}\n')
            f.write('#SOFTWARE=pmutil v1.2\n')
            f.write('#DELIM=,\n')
            f.write('#DATE=JD\n')
            f.write('#OBSTYPE=DSLR\n')
            f.write(aavso)

    def convertToVSNETData(self, res, color: str):
        """
            VSNET formátum:
                - mező szeparátor: space, lehet több is, úgy táblázatot lehet csinálni belőle
                - 4 mező: NAME, OBSDATE, MAG, OBSERVER
                - NAME: constellation + id, greek letters and 'nova' lowercase, otherwise uppercase (CYGSS)
                - DATE: in UT, numbers only (20000720.123)
                    - 1 mp pontosság 5 tizedes lenne
                    - 4 tizedes 9 mp pontosság
                - MAG: nincs megkötés, lehet 14.999 formátumban is akár, utána szűrőkód
                - OBSERVER: nincs megkötés, lehet kód és név is
        """
        cc = color[:1]
        sc = self.aavso_std_colors[color]

        CONSTELLATIONS = [
            "AND","ANT","APS","AQR","AQL","ARA","ARI","AUR","BOO","CAE",
            "CAM","CNC","CVN","CMA","CMI","CAP","CAR","CAS","CEN","CEP",
            "CET","CHA","CIR","COL","COM","CRA","CRB","CRV","CRT","CRU",
            "CYG","DEL","DOR","DRA","EQU","ERI","FOR","GEM","GRU","HER",
            "HOR","HYA","HYI","IND","LAC","LEO","LMI","LEP","LIB","LUP",
            "LYN","LYR","MEN","MIC","MON","MUS","NOR","OCT","OPH","ORI",
            "PAV","PEG","PER","PHE","PIC","PSC","PSA","PUP","PYX","RET",
            "SGE","SGR","SCO","SCL","SCT","SER","SEX","TAU","TEL","TRI",
            "TRA","TUC","UMA","UMI","VEL","VIR","VOL","VUL"
        ]

        GREEK_LETTERS = [
            "ALF", "BET", "GAM", "DEL", "EPS", "ZET", "ETA", "TET", "IOT", "KAP", "LAM", "MIU",
            "NIU", "KSI", "OMI", "PI", "RHO", "SIG", "TAU", "UPS", "PHI", "KHI", "PSI", "OME"
        ]

        # format star id into vsnet style
        star_name = res['LABEL'].replace('_', ' ')
        sn = star_name.rsplit(' ', 1)
        if len(sn) == 2 and len(sn[1]) == 3 and sn[1].upper() in CONSTELLATIONS:
            sn_id = sn[0].lower() if (len(sn[0]) > 2 and sn[0].upper() in GREEK_LETTERS) or "nova" in sn[0].lower() else sn[0].upper()
            star_name = sn[1].upper() + sn_id
        star_name = star_name.replace(' ', '')

        # format observation date
        obs_date = Time(res['JD'], format='jd', scale='utc')
        date_s = f"{obs_date.strftime("%Y%m%d")}{round(obs_date.mjd % 1, 5)}"

        visFlag = res['VIZ_FLAG']['GBR'.index(cc)]
        fainter = "<" if visFlag in "IB" else ""

        if fainter == '<':
            magnitude = self.mgLimits[color]
            aavsoColor = self.aavso_colors[color]
        else:
            # std = True
            m = res['MAG_STD' + sc]
            aavsoColor = self.aavso_std_colors[color]
            if str(m) == '-' or str(m) == '99.0' or m == 99.0:
                m = res['MAG_T' + cc]
                aavsoColor = self.aavso_colors[color]
                # std = False
            if str(m) == '-':
                return None

            magnitude = f"{float(m):4.3f}"

        return star_name, date_s, fainter, magnitude, aavsoColor, self.opt['name']

    @staticmethod
    def convertToVSNETRecord(data):
        # (starid, date, fainter, magnitude, flt, observer)
        return f'{data[0]:<20} {data[1]} {data[2]}{data[3]}{data[4]:<2} {data[5]}\n'

    def convertToVSNETRecords(self, pmt: Table, colors: list[str]) -> list:
        varMask = pmt['ROLE'] == 'V'
        vt = pmt[varMask]
        vsnetResults = []
        for v in vt:
            for c in colors:
                vsnetData = self.convertToVSNETData(v, c)
                if vsnetData:
                    vsnet = self.convertToVSNETRecord(vsnetData)
                    vsnetResults.append(vsnet)
        return vsnetResults

    @staticmethod
    def saveVSNETReport(vsnet, outFolder, postfix=""):
        r = outFolder.split('/')
        fileName = r[-1]
        outFileName = outFolder + '/' + fileName + postfix + '.vsnet'
        pm.printInfo(f"Save VSNET report to the file {outFileName}")

        with open(outFileName, 'w') as f:
            f.write(vsnet)

    def addToLCData(self, aavsoData):
        name = aavsoData[0]
        #        jd = aavsoData[1]
        #        fainter = aavsoData[2]
        #        mg = aavsoData[3]
        #        err = aavsoData[4]
        #        color = aavsoData[5]
        lcr = aavsoData[1:6]
        if name not in self.lcdata:
            self.lcdata[name] = [lcr]
        else:
            self.lcdata[name].append(lcr)

    def drawLCs(self, save_folder):
        for name in self.lcdata:
            if len(self.lcdata[name]) > 3:
                self.drawLC(name, save_folder)

    def drawLC(self, name, save_folder):
        if self.opt["showGraphs"] or self.opt['saveGraphs']:
            lc = {
                'r': [[], [], []],
                'g': [[], [], []],
                'b': [[], [], []],
                'fr': [[], []],
                'fg': [[], []],
                'fb': [[], []]
            }
            for data in self.lcdata[name]:
                color = data[4].lower()
                if color[0] == 't':  # TG, TR, TB
                    color = color[1:]
                elif color == 'v':
                    color = 'g'
                if data[1] != '<':
                    lc[color][0].append(data[0])
                    lc[color][1].append(float(data[2]))
                    err = float(data[3]) if data[3] != "na" else 0.0
                    if err > 1.0:
                        err = 1.0
                    lc[color][2].append(err)
                else:
                    lc['f' + color][0].append(data[0])
                    lc['f' + color][1].append(float(data[2]))

            plot = pm.Plot(1, self.opt["showGraphs"], self.opt['saveGraphs'])
            ax = plt.subplot(111)
            ax.invert_yaxis()
            # plt.locator_params(axis='y', nbins=10)
            plt.xlabel("JD")
            plt.ylabel("mg")
            plt.title(name)
            for color in ['r', 'g', 'b']:
                if len(lc[color][0]) > 0:
                    # plt.plot(lc[color][0], lc[color][1], color + 'o')
                    plt.errorbar(lc[color][0], lc[color][1], lc[color][2], None, color + 'o')
                if len(lc['f' + color][0]) > 0:
                    plt.plot(lc['f' + color][0], lc['f' + color][1], color + 'v')

            lcfile = save_folder + ("/" if save_folder[-1] != '/' else "") + name + ".png"
            pm.printDebug(f"Save LC of {name} into file {lcfile}")
            plot.showOrSave(lcfile)

    @staticmethod
    def findCompStar(results):
        auid = pm.getTableComment(results, 'CompStar')
        return results.loc['AUID', auid] if auid else None

    def getMgLimits(self, cmb):
        hmgGi = pm.getTableComment(cmb, 'MgLimitGi')
        hmgBi = pm.getTableComment(cmb, 'MgLimitBi')
        hmgRi = pm.getTableComment(cmb, 'MgLimitRi')
        self.mgLimits = {'Gi': hmgGi, 'Bi': hmgBi, 'Ri': hmgRi}

    def processFiles(self, fileNames, postfix=""):
        report_types = self.opt['rpt'] if isinstance(self.opt['rpt'], list) else [self.opt['rpt']]
        aavso = []
        vsnet = []
        for fileName in fileNames:
            pmResult = loadCatalog(fileName)
            pmResult.add_index('AUID')

            baseFolder = fileName[:fileName.find(pm.setup['PHOT_FOLDER_NAME'])]
            chartId = self.findChartId(baseFolder)
            self.chartId = chartId if chartId else self.NOT_AVAILABLE

            self.getMgLimits(pmResult)

            if REPORT_TYPE_AAVSO in report_types:
                aavso += self.convertToAAVSORecords(pmResult, ['Ri', 'Gi', 'Bi'])
            if REPORT_TYPE_VSNET in report_types:
                vsnet += self.convertToVSNETRecords(pmResult, ['Ri', 'Gi', 'Bi'])

        aavso.sort()
        vsnet.sort()

        if REPORT_TYPE_AAVSO in report_types:
            self.saveAAVSOReport(''.join(aavso), self.opt['out'], self.opt['name'], postfix)
        if REPORT_TYPE_VSNET in report_types:
            self.saveVSNETReport(''.join(vsnet), self.opt['out'], postfix)

    def process(self):

        if not self.opt['files']:
            f = glob.glob(pm.setup['PHOT_FOLDER_NAME'] + '/' + pm.setup['SEQ_FILE_PREFIX'] + '*.cmb.pm')
            f.sort()
            self.opt['files'] = f

        seqFilenames = list(filter(lambda fn: "Seq_" in fn, self.opt['files']))
        cmbFilenames = list(filter(lambda fn: "Combined" in fn, self.opt['files']))

        if len(seqFilenames) > 0:
            self.processFiles(seqFilenames)

            baseFolder = self.opt['files'][0][:self.opt['files'][0].find(pm.setup['PHOT_FOLDER_NAME'])]
            self.drawLCs(baseFolder)

        if len(cmbFilenames) > 0:
            postfix = "_cmb" if len(seqFilenames) > 0 else ""
            self.processFiles(cmbFilenames, postfix)


if __name__ == '__main__':

    class MainApp:

        opt: dict[str, str|list[str]|None] = {
            'out': None,  # output folder
            'rpt': REPORT_TYPE_AAVSO,  # report format, default: aavso extended
            'name': '?',  # observer name code
            'method': 'gcx',
            # mg calculation method: comp - use best comp star, gcx - m=1 linear fit ensemble, lfit - general linear fit ensemble
            'files': None,
        }

        def __init__(self, argv):
            self.argv = argv
            pass

        def processCommands(self):
            try:
                optlist, args = getopt.getopt(sys.argv[1:], "o:n:", ['out=', 'name=', 'vsn'])
            except getopt.GetoptError:
                pm.printError('Invalid command line options')
                return

            for o, a in optlist:
                if a[:1] == ':':
                    a = a[1:]
                if o == '-o' or o == '--out':
                    self.opt['out'] = a
                elif o == '--vsn':
                    self.opt['rpt'] = [REPORT_TYPE_AAVSO, REPORT_TYPE_VSNET]
                elif o == '-n' or o == '--name':
                    self.opt['name'] = a

            baseFolder = args[0].replace("\\", "/")
            if baseFolder.endswith("/"):
                baseFolder = baseFolder[:-1]
            if self.opt['out'] is None:
                self.opt['out'] = baseFolder

            self.opt['files'] = [f"{baseFolder}/Photometry/Combined.cmb.pm"]

        def run(self):
            self.processCommands()

            proc = ReportProcessor(self.opt)
            proc.process()


    app = MainApp(sys.argv)
    app.run()

# end main.
