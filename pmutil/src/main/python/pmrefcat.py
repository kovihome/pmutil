#!/usr/bin/env python3
#
# PmUtils/pmrefcat
#
'''
Created on Jan 1, 2020

@author: kovi
'''

from sys import argv
from getopt import getopt, GetoptError
from urllib import request
from datetime import datetime
from os import makedirs, symlink, remove
from os.path import isdir, exists
from glob import glob
import json
import xmltodict
import re

from pmbase import printError, printWarning, printInfo, saveCommand, loadPplSetup, invoke, Blue, Color_Off, BGreen, hexa2deg, deg2hexa

class RefCat:

    opt = {}                 # command line options
    pplSetup = {}            # PPL setup from ppl-setup config file

    defMgLimit = 18.0
    defFov = 60

    chartId = None

    catalogHeader = "AUID         ROLE RA          RA_DEG         DEC         DEC_DEG       MAG_B  ERR_B   MAG_V  ERR_V   MAG_R  ERR_R   LABEL"

    def __init__(self, opt):
        self.opt = opt

    def loadVspPhotometryData(self, outFile, objectName, ra = None, dec = None, fov = defFov):
        '''
        star*: name of the star to plot. You must provide EITHER the star parameter OR ra and dec parameters (see below)
        ra*: right ascension
        dec*: declination
        fov*: field of view, in arcminutes.
        maglimit*: magnitude limit for the chart; stars fainter than the maglimit will not be plotted.
        other: other variable stars to mark on the chart. Omit or leave blank for none, "gcvs" for GCVS variables, "all" for all variables
        '''
        if objectName != None:
            ids = "star=%s" % (objectName.replace(' ', '+'))
        else:
            ids = "ra=%s&dec=%s" % (ra, dec)
        vspUrl = "https://www.aavso.org/apps/vsp/api/chart/?format=json&%s&fov=%d&maglimit=%4.1f&other=all" % (ids, fov, self.defMgLimit)
        f = request.urlopen(vspUrl)
        respJson = f.read().decode('UTF-8')

        resp = json.loads(respJson)

        # print("Star:", resp['star'])
        # print("AUID:", resp['auid'])
        # print("Coord:", resp['ra'], resp['dec'])
        # print("ChartID:", resp['chartid'])
        # print("Photometry table:")

        self.chartId = resp['chartid']

        for pm in resp['photometry']:

            # 'V', 'B', 'Rc', 'Ic', 'J', 'H', 'K'
            bands = {}
            for b in pm['bands']:
                bands[b['band']] = b

            auid = pm['auid']
            role = "C"
            ra = pm['ra']
            dec = pm['dec']
            raDeg = "%10.8f" % (hexa2deg(ra) * 15.0)
            if not dec.startswith('-'):
                dec = "+" + dec
            decDeg = "%+10.8f" % (hexa2deg(dec))
            if 'B' in bands:
                mgB = str(bands['B']['mag'])
                mgerrB = str(bands['B']['error'])
                if mgerrB == 'None':
                    mgerrB = "-"
            else:
                mgB = "-"
                mgerrB = "-"
            if 'V' in bands:
                mgV = str(bands['V']['mag'])
                mgerrV = str(bands['V']['error'])
                if mgerrV == 'None':
                    mgerrV = "-"
            else:
                mgV = "-"
                mgerrV = "-"
            if 'Rc' in bands:
                mgR = str(bands['Rc']['mag'])
                mgerrR = str(bands['Rc']['error'])
                if mgerrR == 'None':
                    mgerrR = "-"
            else:
                mgR = "-"
                mgerrR = "-"
            label = str(pm['label']).replace(' ', '_')

            s = auid.ljust(13) + role.ljust(5) + ra.ljust(12) + raDeg.ljust(15) + dec.ljust(12) + decDeg.ljust(14) + mgB.ljust(7) + mgerrB.ljust(8) + mgV.ljust(7) + mgerrV.ljust(8) + mgR.ljust(7) + mgerrR.ljust(8) + label
            if outFile != None:
                outFile.write(s + "\n")
            else:
                print(s)


    def loadVsxCatalogData(self, outFile, objectName, ra = None, dec = None, fov = defFov, auidOnly = True):
        '''
        &coords   The central J2000 RA/DEC coordinates for a radius search, expressed sexagesimally (by default), or in decimal degrees if format is set to d. Northern hemisphere coordinates must have the plus () sign URL-encoded as %2B. Space characters between all other figures must be replaced with the URL whitespace character (). The order argument (which see) must also be included in the query string with its value set to 9 in order to prompt VSX to display distances from the central coordinates in the results listing. Default is empty string (no radius search).
        &ident    Object identification for name searches. Space characters must be replaced with the URL whitespace character (+). Other special characters may also need to be URL-encoded. Default is empty string (no name search).
        &format   Explicit specification for format of coords. For sexagesimal, this value should be s. For decimal degrees, this value should be d. Default is s (sexagesimal).
        &geom     The geometry for central coordinate-based searches. For radius searches, this value should be r. For box searches, this value should be b. Default is r (radius search).
        &size     For box searches (geom=b), the width of the box. For radius searches (geom=r), the radius of the circle. Expressed in the units specified by unit (see next). Default is 10.0.
        &unit     The unit of measurement used for the value given in size (see above). For arc degrees, this value should be 1. For arc minutes, this value should be 2. For arc seconds, this value should be 3. Default is 2 (arc minutes).
        '''
        if objectName != None:
            # https://www.aavso.org/vsx/index.php?view=results.csv&ident=R+Car
            vsxUrl = "https://www.aavso.org/vsx/index.php?view=api.object&format=json&ident=%s" % (objectName.replace(' ', '+'))
            # print(vsxUrl)
            f = request.urlopen(vsxUrl)
            respJson = f.read().decode('UTF-8')
            resp = json.loads(respJson)
            # print(resp)
            # print("Object:", resp['VSXObject']['Name'], "RA:", resp['VSXObject']['RA2000'], "Dec:", resp['VSXObject']['Declination2000'])

            sra = deg2hexa(float(resp['VSXObject']['RA2000']) / 15.0)
            sdec = deg2hexa(float(resp['VSXObject']['Declination2000']))
        else:
            sra = deg2hexa(float(ra) / 15.0)
            sdec = deg2hexa(float(dec))

        coords = sra + "+" + sdec

        vsxUrl = "https://www.aavso.org/vsx/index.php?view=query.votable&format=d&coords=%s&size=%d&unit=2&geom=b" % (coords, fov)
        # https://www.aavso.org/vsx/index.php?view=query.votable&coords=19:54:17+32:13:08&format=d&size=60&unit=2
        # print(vsxUrl)
        f = request.urlopen(vsxUrl)
        respVOTable = f.read().decode('UTF-8')
        # resp = json.loads(respJson)
        # print(respVOTable)
        votable = xmltodict.parse(respVOTable)
        # print(votable)
        # print(votable['VOTABLE']['RESOURCE']['TABLE'].keys())
        # print(votable['VOTABLE']['RESOURCE']['TABLE']['FIELD'])
        # print(votable['VOTABLE']['RESOURCE']['TABLE']['DATA']['TABLEDATA'])
        trs = votable['VOTABLE']['RESOURCE']['TABLE']['DATA']['TABLEDATA']['TR']
        nrUnknown = 1
        for tr in trs:
            # [None, 'PS1-3PI J185203.12+325154.5', 'Lyr', '283.01304000,32.86516000', 'RRAB', '16.810', 'r', '0.400', 'r', '57000.85600', None, '0.839122', None, None, None]
            # print("AUID:", tr['TD'][0], "Coord:", tr['TD'][3], "LABEL", tr['TD'][1])
    
            # #0.00    Variable    AY Lyr                000-BCD-108    18 44 26.69 +37 59 51.9    Lyr    UGSU    0.0733    12.5 - 18.4 B
            # distance, varstate, label(1), auid(0), coords(3), constell(2), type(4), period(11), max(5,6) - min(7,8)
            # print("   Constell:", tr['TD'][2], "Type:", tr['TD'][4], "Period:", tr['TD'][11], "Mg:", tr['TD'][5] + tr['TD'][6], "-", tr['TD'][7] + tr['TD'][8])
    
            auid = tr['TD'][0]
            if auid == None:
                if auidOnly:
                    continue
                auid = '999-VAR-%03d' % (nrUnknown)
                nrUnknown = nrUnknown + 1
                tr['TD'][0] = auid
            role = 'V'
            raDeg, decDeg = tr['TD'][3].split(',')
            ra = deg2hexa(float(raDeg) / 15.0)
            dec = deg2hexa(float(decDeg))
            if not dec.startswith('-'):
                dec = "+" + dec
                decDeg = "+" + decDeg
            label = tr['TD'][1].replace(' ', '_')

            s = auid.ljust(13) + role.ljust(5) + ra.ljust(12) + raDeg.ljust(15) + dec.ljust(12) + decDeg.ljust(14) + "-".ljust(7) + "-".ljust(8) + "-".ljust(7) + "-".ljust(8) + "-".ljust(7) + "-".ljust(8) + label
            if outFile != None:
                outFile.write(s + "\n")
            else:
                print(s)

        nrUnknown = 1
        for tr in trs:

            auid = tr['TD'][0]
            if auid == None:
                if auidOnly:
                    continue
                auid = '999-VAR-%03d' % (nrUnknown)
                nrUnknown = nrUnknown + 1
            role = 'V'
            raDeg, decDeg = tr['TD'][3].split(',')
            ra = deg2hexa(float(raDeg) / 15.0)
            dec = deg2hexa(float(decDeg))
            label = tr['TD'][1]
            constell = tr['TD'][2]
            vartype = tr['TD'][4]
            if vartype == None:
                vartype = "-"
            period = tr['TD'][11]
            maxmg = tr['TD'][5]
            if tr['TD'][6] != None:
                maxmg += tr['TD'][6]
            minmg = tr['TD'][7]
            if tr['TD'][8] != None:
                minmg += tr['TD'][8]

            s = "#" + label.ljust(29) + auid.ljust(14) + ra.ljust(15) + dec.ljust(14) + constell.ljust(4) + vartype.ljust(11) + period.ljust(11) + maxmg + " - " + minmg
            if outFile != None:
                outFile.write(s + "\n")
            else:
                print(s)

    def getPair(self, s, delim = ','):
        ss = s.split(delim)
        return [ss[0].strip(), ss[1].strip()]


    def determineCoordsFromImage(self, imageFileName):
        configFolder = self.pplSetup["PMLIB"]
        if configFolder == None:
            configFolder = "$HOME/.pmlib"
        SOLVE_ARGS = "-O --config " + configFolder + "/astrometry.cfg --use-sextractor --sextractor-path sextractor -r -y -p"

        # /usr/local/astrometry/bin/solve-field $SOLVE_ARGS -D $2 -N $AST_FILE $f

        makedirs("temp", exist_ok = True)

        invoke("cp " + imageFileName + " temp/src.fits")

        astFolder = self.pplSetup["AST_BIN_FOLDER"]
        if astFolder == None:
            astFolder = "/usr/local/astrometry/bin"

        astResult = invoke(astFolder + "/solve-field " + SOLVE_ARGS + " -D temp -N temp/ast.fits temp/src.fits")

        r = astResult.split('\n')
        sCoords = ["0.0", "0.0"]
        for line in r:
            if line.startswith("Field center"):
                if line.find("RA,Dec") > -1:
                    sCoords = self.getPair(line.split('(')[2].split(')')[0])
                    printInfo("Image center: %s %s" % (sCoords[0], sCoords[1]))
                else:
                    sCoordsSexa = self.getPair(line.split('(')[2].split(')')[0])
                    printInfo("Image center: %s %s" % (sCoordsSexa[0], sCoordsSexa[1]))
            elif line.startswith("Field size"):
                sFieldSize = self.getPair(line.split(':')[1].split('a')[0], 'x')
                printInfo("Image size: %s' x %s'" % (sFieldSize[0], sFieldSize[1]))

        return sCoords

    def execute(self):
        self.pplSetup = loadPplSetup()

        if self.opt['image'] != None:

            coords = self.determineCoordsFromImage (self.opt['image'])
            self.opt['ra'] = coords[0]
            self.opt['dec'] = coords[1]

#        print("coords:", coords)

        outFileName = self.pplSetup['CONFIG_FOLDER']
        if not outFileName.endswith('/'):
            outFileName += '/'
        outFileName += 'cat/'
        if not exists(outFileName):
            makedirs(outFileName)

        if self.opt['object'] != None:
            outFileName += self.opt['object'].lower().replace(' ', '_') + '.cat'
        else:
            if not self.opt['dec'].startswith('+') and not self.opt['dec'].startswith('-'):
                self.opt['dec'] = '+' + self.opt['dec']
            outFileName += self.opt['ra'].replace(':', '')[:4] + self.opt['dec'].replace(':', '')[:5] + '.cat'

        if not exists(outFileName) or self.opt['overwrite']:

            outFile = open(outFileName, 'w')
            outFile.write(self.catalogHeader + "\n")

            if self.opt['object'] != None:

                self.loadVspPhotometryData(outFile, self.opt['object'], fov = self.opt['field'])

                self.loadVsxCatalogData(outFile, self.opt['object'], fov = self.opt['field'], auidOnly = self.opt['auidOnly'])

            else:

                self.loadVspPhotometryData(outFile, None, ra = self.opt['ra'], dec = self.opt['dec'], fov = self.opt['field'])

                self.loadVsxCatalogData(outFile, None, ra = self.opt['ra'], dec = self.opt['dec'], fov = self.opt['field'], auidOnly = self.opt['auidOnly'])

            outFile.write('### chartid: ' + self.chartId + '\n')
            outFile.write('### fov: ' + str(self.opt['field']) + ' arcmin\n')

            outFile.close()

            printInfo("Reference catalog file %s created." % (outFileName))

        else:

            printInfo("Reference catalog file %s is exists." % (outFileName))

        if self.opt['folder']:
        
            for f in glob(self.opt['folder'] + '*.cat'):
                remove(f)
            symlink(outFileName, self.opt['folder'] + 'ref.cat')


class MainApp:

    appName        = 'ppl-refcat'
    appVersion     = '1.1.0'
    appDescription = 'Create reference catalog for photometry.'

    opt = {
        'coords'   : None,   # coordinates of reference field
        'source'   : None,   # source catalog of field stars
        'object'   : None,   # object (variable star) name
        'image'    : None,   # image file
        'field'    : RefCat.defFov, # reference field size in arcmins
        'auidOnly' : True,   # variables having auid only
        'overwrite': False,  # overwrite catalog
        'folder'   : None,   # reference catalog file name
        }

    def __init__(self, argv):
        self.argv = argv

    def printVersion(self):
        print(BGreen + self.appName + ", version " + self.appVersion + Color_Off)

    def printTitle(self):
        print()
        self.printVersion()
        print(Blue + self.appDescription + Color_Off)

    def usage(self):
        print()
        self.printVersion()
        print()
        print("Usage: " + self.appName + " [OPTIONS]... FOLDER_NAME")
        print(self.appDescription)
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -o,  --object object_name   object (variable star) name")
        print("  -c,  --coords ra,decl       coordinates of the center of reference frame, valid format is 12:34:56.7,-12:34:56.7")
        print("  -i,  --image filename       image file name")
        print("  -s,  --source catalog       source catalog for field stars")
        print("  -f,  --field size           field size in arcmin, default is 60 arcmin")
        print("  -a,  --all                  collect all variables ; if not set, collect variables having AUID only")
        print("  -w,  --overwrite            overwrite catalog file, if exists")
        print("  -h,  --help                 print this page")

    def processCommands(self):

        try:
            optlist, args = getopt (argv[1:], "ac:s:o:i:f:wh", ['--all', '--coord', '--source', '--object', '--image', '--field', '--overwrite', '--help'])
        except GetoptError:
            print ('Invalid command line options')
            exit(1)

        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-c':
                self.opt['coords'] = a
            elif o == '-s':
                self.opt['source'] = a
            elif o == '-o':
                self.opt['object'] = a.replace('_', ' ')
            elif o == '-i':
                self.opt['image'] = a
            elif o == '-f':
                if a.isdigit():
                    self.opt['field'] = int(a)
                else:
                    print("Invalid field size parameter: %s. Use default 60 arcmin instead." % (a))
            elif o == '-a':
                self.opt['auidOnly'] = False
            elif o == '-w':
                self.opt['overwrite'] = True
            elif o == '-h':
                self.usage()
                exit(0)

        if len(args) > 0:
            if not isdir(args[0]):
                printError('%s is not a folder.' % (args[0]))
                exit(1)
            self.opt['folder'] = args[0]
            if not self.opt['folder'].endswith('/'):
                self.opt['folder'] += '/'

        if self.opt['object'] == None and self.opt['coords'] == None and self.opt['image'] == None:
            printError("Either object name (-o), coordinates (-c) or image file (-i) must be given.")
            exit(1)

        if self.opt['coords'] != None:
            c = self.opt['coords'].split(',', maxsplit = 2)
            ok = True
            if len(c) != 2:
                ok = False
            else:
                if re.fullmatch("(\d){2}:(\d){2}:(\d){2}(\.\d)*", c[0]) == None:
                    ok = False
                if re.fullmatch("[+-]*(\d){2}:(\d){2}:(\d){2}(\.\d)*", c[1]) == None:
                    ok = False
            if not ok:
                printError("Invalid coordinate format")
                exit(1)

            self.opt['ra'] = c[0]
            self.opt['dec'] = c[1]

    def run(self):
        self.printTitle()
        self.processCommands()
        saveCommand(self.opt['folder'], self.argv, 'refcat')

        start = datetime.now()

        cat = RefCat(self.opt)
        cat.execute()

        exectime = (datetime.now() - start).total_seconds()
        print("%sexecution time was %d seconds.%s" % (Blue, exectime, Color_Off))


if __name__ == '__main__':

    app = MainApp(argv)
    app.run()

# end main.
