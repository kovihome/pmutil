#!/usr/bin/env python3
#
# PmUtils/pmresult
#
'''
Created on Apr 28, 2019

@author: kovi

AAVSO extended format spec:

        header:
        #TYPE=Extended
        #OBSCODE=?
        #SOFTWARE=
        #DELIM=,
        #DATE=JD
        #OBSTYPE=DSLR
    
        data:
        STARID: The star's identifier. It can be the AAVSO Designation, the AAVSO Name, or the AAVSO Unique Identifier, but NOT more than one of these. 
		Limit: 30 characters.
        DATE: The date of the observation, in the format specified in the DATE parameter. Limit: 16 characters.
        MAGNITUDE: The magnitude of the observation. Prepend with < if a fainter-than.  A dot is required (e.g. "9.0" rather than "9"). Limit: 8 characters.
        MAGERR: Photometric uncertainty associated with the variable star magnitude. If not available put "na". Limit: 6 characters.
        FILTER: The filter used for the observation. This can be one of the following letters (in bold):
            B: Johnson B
            V: Johnson V
            R: Cousins R
            TG: Green Filter (or Tri-color green). This is commonly the "green-channel" in a DSLR or color CCD camera. These observations use V-band comp star magnitudes.
            TB: Blue Filter (or Tri-color blue). This is commonly the "blue-channel" in a DSLR or color CCD camera. These observations use B-band comp star magnitudes.
            TR: Red Filter (or Tri-color red). This is commonly the "red-channel" in a DSLR or color CCD camera. These observations use R-band comp star magnitudes.
        TRANS: YES if transformed using the Landolt Standards or those fields that contain secondary standards, or NO if not. 
		Document the method used to transform in the "NOTES" section.
        MTYPE: Magnitude type. STD if standardized (Click here for definition of standardized) or DIF if differential (very rare). 
		If you are currently using ABS for 'absolute' we will still accept it.
            Differential requires the use of CNAME.
        CNAME: Comparison star name or label such as the AUID (much preferred) or the chart label for the comparison star used. 
		If not present, use "na". Limit: 20 characters.
        CMAG: Instrumental magnitude of the comparison star. If "ensemble" see example below. If not present, use "na". Limit: 8 characters.
        KNAME: Check star name or label such as the AUID (much preferred) or the chart label for the check star. If not present, use "na". Limit: 20 characters.
        KMAG: Instrumental magnitude of the check star. If "ensemble" see example below. If not present, use "na".Limit: 8 characters.
        AIRMASS: Airmass of observation Limit 7 characters - entry will be truncated if longer than that. If not present, use "na".
        GROUP: Grouping identifier (maximum 5 characters). It is used for grouping multiple observations together, usually an observation set 
		that was taken through multiple filters. 
            It makes it easier to retrieve all magnitudes from a given set in the database, say, if someone wanted to form color indices such as (B-V) with them. 
		If you are just doing time series, 
            or using the same filter for multiple stars, etc., just set GROUP to "na." For cases where you want to group observations, GROUP should be an integer, 
		identical for all observations in a group,
            and unique for a given observer for a given star on a given Julian Date. Limit: 5 characters.
        CHART: Please use the sequence ID you will find written in Bold print near the top of the photometry table in a sentence that reads 
		"Report this sequence as [ID] in the chart field of your observation report." If a non–AAVSO sequence was used, please describe it as clearly as possible. 			Limit: 20 characters.
        NOTES: Comments or notes about the observation. If no comments, use "na". This field has a maximum length of several thousand characters, 
		so you can be as descriptive as necessary. 
            One convention for including a lot of information as concisely as possible is to use subfields in the format |A=B; the '|' character is the separator, 
		A is a keyvalue name and B is its value. 
            If you need an alternative delimiter from '|', use it but preceed the first instance with "DELIM=". 
		Using this mechanismnyou can document your transform process in more detail. 
            Here is an example as used by TransformApplier:
            5 records aggregated|VMAGINS=-7.244|VERR=0.006|CREFMAG=13.793|CREFERR=0.026|
            KREFMAG=14.448|KREFERR=0.021|VX=1.1501|CX=1.1505|KX=1.1500|Tv_bv=0.0090|Tv_bvErr=0.0100|
            TAver=2.47
'''

import getopt
import sys
import glob
from os.path import exists
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
from astropy.time import Time
import astropy.units as u

import pmbase as pm

# manage catalog file
def loadCatalog(refFileName):
    table = Table.read(refFileName, format='ascii')
    if not table:
        print (("Reference catalog %s not found") % refFileName)
        return None
    return table


# create report
class ReportProcessor:

    NOT_AVAILABLE = 'na'

    opt = {}  # command line options
    pos = None  # catalog header positions

    chartId = NOT_AVAILABLE

    aavso_colors = {
        'Gi': 'TG',
        'Ri': 'TR',
        'Bi': 'TB' 
    }
    aavso_std_colors = {
        'Gi' : 'V',
        'Ri' : 'R',
        'Bi' : 'B'
    }

    def __init__(self, opt):
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

    def calcAirmass(self, ra, dec, obsDate):
        co = SkyCoord(ra, dec, frame = ICRS, unit = (u.hourangle, u.deg))
        t = Time(obsDate.replace('T', ' '))
        # TODO: store observer's location in config
        loc = EarthLocation.of_address('Budapest')
        hc = co.transform_to(AltAz(obstime = t, location = loc))
        return float(hc.secz)

    def convertToAAVSORecord(self, res, color, compStar):

        # MAG_Tc, ERR_Tc, MAG_STDc, ERR_STDc

        cc = color[:1]
        sc = self.aavso_std_colors[color]

        starid = res['LABEL'].replace('_', ' ')
        date = res['JD']

        visFlag = res['VIZ_FLAG']['BGR'.index(cc)]
        fainter = "<" if visFlag in "IB" else ""

        if fainter == '<':
            magnitude = self.mgLimits[color]
            magerr = self.NOT_AVAILABLE
            aavsoColor = self.aavso_colors[color]
        else:
            std = True
            m = res['MAG_STD' + sc]
            e = res['ERR_STD' + sc]
            aavsoColor = self.aavso_std_colors[color]
            if str(m) == '-':
               m = res['MAG_T' + cc]
               e = res['ERR_T' + cc]
               aavsoColor = self.aavso_colors[color]
               std = False
            if str(m) == '-':
               return None
       
            magnitude = "%4.3f" % (float(m))
            magerr = "%4.3f" % (float(e)) if e != '-' else self.NOT_AVAILABLE

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
        airmass = "%6.4f" % (self.calcAirmass(res['RA'], res['DEC'], res['DATE-OBS']))
        group = self.NOT_AVAILABLE
        chart = self.chartId
        notes = self.NOT_AVAILABLE
        return '%s,%s,%s%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % \
            (starid.upper(), date, fainter, magnitude, magerr, flt, trans, mtype, cname, cmag, kname, kmag, airmass, group, chart, notes)


    def convertToAAVSORecords(self, pmt, colors):

        compStar = self.findCompStar(pmt)

        varMask = pmt['ROLE'] == 'V'
        vt = pmt[varMask]
        aavsoResults = ''
        for v in vt:
            for c in colors:
                aavso = self.convertToAAVSORecord(v, c, compStar)
                if aavso:
                    aavsoResults += aavso
        return aavsoResults


    def saveAAVSOReport(self, aavso, outFolder, obsName):
        r = outFolder.split('/')
        fileName = r[-1]

        outFileName = outFolder + '/' + fileName + '.extended.aavso'
        print("Print AAVSO extended report to the file", outFileName)
        r = open(outFileName, 'w')

        r.write('#TYPE=Extended\n')
        r.write('#OBSCODE=%s\n' % (obsName))
        r.write('#SOFTWARE=pmutil v1.2\n')
        r.write('#DELIM=,\n')
        r.write('#DATE=JD\n')
        r.write('#OBSTYPE=DSLR\n')

        r.write(aavso)

        r.close()


    def findCompStar(self, results):
        auid = pm.getTableComment(results, 'CompStar')
        return results.loc['AUID', auid] if auid else None


    def getMgLimits(self, cmb):
        hmgGi = pm.getTableComment(cmb, 'MgLimitGi')
        hmgBi = pm.getTableComment(cmb, 'MgLimitBi')
        hmgRi = pm.getTableComment(cmb, 'MgLimitRi')
        self.mgLimits = { 'Gi' : hmgGi, 'Bi': hmgBi, 'Ri': hmgRi }


    def process(self):

        if not self.opt['files']:
            f = glob.glob(pm.setup['PHOT_FOLDER_NAME'] + '/' + pm.setup['SEQ_FILE_PREFIX'] + '*.cmb.pm')
            f.sort()
            self.opt['files'] = f

        headers = []
        aavso = ''

        print(self.opt['files'])

        for fileName in self.opt['files']:

            pmResult = loadCatalog(fileName)
            pmResult.add_index('AUID')

            baseFolder = fileName[:fileName.find(pm.setup['PHOT_FOLDER_NAME'])]
            chartId = self.findChartId(baseFolder)
            self.chartId = chartId if chartId else self.NOT_AVAILABLE

            self.getMgLimits(pmResult)

            aavso += self.convertToAAVSORecords(pmResult, ['Ri', 'Gi', 'Bi'])

        if self.opt['rpt'] == 'aavso':
            self.saveAAVSOReport(aavso, self.opt['out'], self.opt['name'])


if __name__ == '__main__':

    class MainApp:

        opt = {
            'out' : None,     # output folder
            'rpt' : 'aavso',  # report format, default: aavso extended
            'name': '?',      # observer name code
            'method': 'gcx',  # mg calculation method: comp - use best comp star, gcx - m=1 linear fit ensemble, lfit - general linear fit ensemble
            'files': None,
            }

        def __init__(self, argv):
            self.argv = argv
            pass

        def processCommands(self):
            try:
                optlist, args = getopt.getopt (sys.argv[1:], "o:r:n:", ['--out', '--report', '--name'])
            except getopt.GetoptError:
                print ('Invalid command line options')
                return

            for o, a in optlist:
                if a[:1] == ':':
                    a = a[1:]
                elif o == '-o':
                    self.opt['out'] = a
                elif o == '-r':
                    if a == 'aavso':
                        self.opt['rpt'] = a
                    else:
                        print("Invalid report type: " + a + ". Create default aavso extended report instead")
                elif o == '-n':
                    self.opt['name'] = a

            self.opt['files'] = args

        def run(self):
            self.processCommands()

            proc = ReportProcessor(self.opt)
            proc.process()


    app = MainApp(sys.argv)
    app.run()

# end main.

