#!/usr/bin/env python3
#
# PmUtils/pmresult
#
'''
Created on Apr 28, 2019

@author: kovi
'''
import getopt
import sys
import glob
from os.path import exists

from pmbase import loadPplSetup

# manage catalog file


def loadCatalog(refFileName):

    refCat = { 'header': [], 'cat': [] }

    ref = open(refFileName, 'r')
    if ref:
        headerLine = ref.readline()
        refCat['header'] = headerLine.split()

        for refLine in ref:
            refCatLine = refLine.split()
            refCat['cat'].append(refCatLine)
        ref.close()

        print (("Photometric result catalog: %s") % refFileName)
        # print (("%d reference object loaded." % len(refCat['cat'])))
        # print (refCat)
        return refCat

    else:

        print (("Reference catalog %s not found") % refFileName)
        return None


def getHeaderPos(hdr, headerName):
    index = 0
    for header in hdr:
        if header == headerName:
            return index
        index = index + 1
    return None

# create report
class ReportProcessor:

    NOT_AVAILABLE = 'na'

    opt = {}                 # command line options
    pplSetup = {}            # PPL setup from ppl-setup config file
    pos = None               # catalog header positions

    chartId = NOT_AVAILABLE

    aavso_colors = {
        'Gi': 'TG',
        'Ri': 'TR',
        'Bi': 'TB',
        'V' : 'V',
        'R' : 'R',
        'B' : 'B'
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

    def reportForAAVSO(self, allResults, header, outFolder, obsName):
        '''
        header:
        #TYPE=Extended
        #OBSCODE=?
        #SOFTWARE=
        #DELIM=,
        #DATE=JD
        #OBSTYPE=DSLR
    
        data:
        STARID: The star's identifier. It can be the AAVSO Designation, the AAVSO Name, or the AAVSO Unique Identifier, but NOT more than one of these. Limit: 30 characters.
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
        TRANS: YES if transformed using the Landolt Standards or those fields that contain secondary standards, or NO if not. Document the method used to transform in the "NOTES" section.
        MTYPE: Magnitude type. STD if standardized (Click here for definition of standardized) or DIF if differential (very rare). If you are currently using ABS for 'absolute' we will still accept it.
            Differential requires the use of CNAME.
        CNAME: Comparison star name or label such as the AUID (much preferred) or the chart label for the comparison star used. If not present, use "na". Limit: 20 characters.
        CMAG: Instrumental magnitude of the comparison star. If "ensemble" see example below. If not present, use "na". Limit: 8 characters.
        KNAME: Check star name or label such as the AUID (much preferred) or the chart label for the check star. If not present, use "na". Limit: 20 characters.
        KMAG: Instrumental magnitude of the check star. If "ensemble" see example below. If not present, use "na".Limit: 8 characters.
        AIRMASS: Airmass of observation Limit 7 characters - entry will be truncated if longer than that. If not present, use "na".
        GROUP: Grouping identifier (maximum 5 characters). It is used for grouping multiple observations together, usually an observation set that was taken through multiple filters. 
            It makes it easier to retrieve all magnitudes from a given set in the database, say, if someone wanted to form color indices such as (B-V) with them. If you are just doing time series, 
            or using the same filter for multiple stars, etc., just set GROUP to "na." For cases where you want to group observations, GROUP should be an integer, identical for all observations in a group,
            and unique for a given observer for a given star on a given Julian Date. Limit: 5 characters.
        CHART: Please use the sequence ID you will find written in Bold print near the top of the photometry table in a sentence that reads "Report this sequence as [ID] in the chart field 
            of your observation report." If a non–AAVSO sequence was used, please describe it as clearly as possible. Limit: 20 characters.
        NOTES: Comments or notes about the observation. If no comments, use "na". This field has a maximum length of several thousand characters, so you can be as descriptive as necessary. 
            One convention for including a lot of information as concisely as possible is to use subfields in the format |A=B; the '|' character is the separator, A is a keyvalue name and B is its value. 
            If you need an alternative delimiter from '|', use it but preceed the first instance with "DELIM=". Using this mechanismnyou can document your transform process in more detail. 
            Here is an example as used by TransformApplier:
            5 records aggregated|VMAGINS=-7.244|VERR=0.006|CREFMAG=13.793|CREFERR=0.026|
            KREFMAG=14.448|KREFERR=0.021|VX=1.1501|CX=1.1505|KX=1.1500|Tv_bv=0.0090|Tv_bvErr=0.0100|
            TAver=2.47
        '''

        # idPos = getHeaderPos(header, 'AUID')
        labelPos = getHeaderPos(header, 'LABEL')
        mvPos = getHeaderPos(header, 'MAG')
        colorPos = getHeaderPos(header, 'COL')
        errPos = getHeaderPos(header, 'ERR')
        jdPos = getHeaderPos(header, 'JD')
        flagPos = getHeaderPos(header, 'FLAG')

        r = outFolder.split('/')
        fileName = r[-1]

        outFileName = outFolder + '/' + fileName + '.extended.aavso'
        print("Print AAVSO extended report to the file", outFileName)
        r = open(outFileName, 'w')

        r.write('#TYPE=Extended\n')
        r.write('#OBSCODE=%s\n' % (obsName))
        r.write('#SOFTWARE=pmutil v1.0\n')
        r.write('#DELIM=,\n')
        r.write('#DATE=JD\n')
        r.write('#OBSTYPE=DSLR\n')

        for res in allResults:

            color = self.aavso_colors[res[colorPos]]

            starid = res[labelPos].replace('_', ' ')
            date = res[jdPos]
            magnitude = res[mvPos]
            magerr = res[errPos]
            fainter = "<" if res[flagPos] == "F" else ""
            flt = color
            trans = 'NO'
            mtype = 'STD'
            cname = 'ENSEMBLE'
            cmag = self.NOT_AVAILABLE  # TODO: set comp star mg if not ensemble
            kname = self.NOT_AVAILABLE  # TODO: specify ensemble reference comp star when ensemble, or set a check star
            kmag = self.NOT_AVAILABLE  # TODO: set measured mg of check star
            airmass = self.NOT_AVAILABLE  # TODO: calculate airmass
            group = self.NOT_AVAILABLE
            chart = self.chartId
            notes = self.NOT_AVAILABLE
            r.write(('%s,%s,%s%4.3f,%4.3f,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n') % (starid.upper(), date, fainter, float(magnitude), float(magerr), flt, trans, mtype, cname, cmag, kname, kmag, airmass, group, chart, notes))

        r.close()

    def process(self):

        self.pplSetup = loadPplSetup()

        if not self.opt['files']:
            f = glob.glob(self.pplSetup['PHOT_FOLDER_NAME'] + '/' + self.pplSetup['SEQ_FILE_PREFIX'] + '*.fits.cat.pm')
            f.sort()
            self.opt['files'] = f

        allResults = []
        headers = []

        for fileName in self.opt['files']:

            pmResult = loadCatalog(fileName)

            headers = pmResult['header']

            baseFolder = fileName[:fileName.find(self.pplSetup['PHOT_FOLDER_NAME'])]
            chartId = self.findChartId(baseFolder)
            self.chartId = chartId if chartId else self.NOT_AVAILABLE

            allResults.extend(pmResult['cat'])

        if self.opt['rpt'] == 'aavso':
            self.reportForAAVSO(allResults, headers, self.opt['out'], self.opt['name'])


class MainApp:
    
    opt = {
        'out' : None,  # output folder
        'rpt' : 'aavso',  # report format, default: aavso extended
        'name': '?',  # observer name code
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
#        if not self.opt['out']:
#            self.opt['out'] = args[0] + '.refout'

    def run(self):
        self.processCommands()

        proc = ReportProcessor(self.opt)
        proc.process()


if __name__ == '__main__':

    app = MainApp(argv)
    app.run()

# end main.

