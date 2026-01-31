#!/usr/bin/env python3
#
# PmUtils/pmrefcat
#
"""
Created on Jan 1, 2020

@author: kovi
"""

from sys import argv
from getopt import getopt, GetoptError
from urllib import request
from datetime import datetime
from os import makedirs, symlink, remove
from os.path import isdir, exists
from glob import glob
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import json
import xmltodict
import re
import urllib

import pmbase as pm
from pmviz import VizierQuery, UCAC4

# Old URL                  New URL                      Location
# www.aavso.org/vsx/*      vsx.aavso.org/*              pmrefcat
# www.aavso.org/apps/*     apps.aavso.org/*             pmrefcat
# www.aavso.org/vsp/*      apps.aavso.org/vsp/*
# www.aavso.org/cgi-bin/*  archive.aavso.org/cgi-bin/*

AAVSO_VSX_URL = "https://vsx.aavso.org"
AAVSO_VSP_URL = "https://apps.aavso.org/vsp"


class RefCat:

    opt = {}  # command line options

    defMgLimit = 18.0
    defFov = 60

    chartId = None

    cache = []

    xmatchTable = Table(names=['AUID', 'RA_DEG', 'DEC_DEG', 'LABEL'], dtype=['U16','U16','U16','str'])

    origRefcat = None

    def __init__(self, opt):
        self.opt = opt

    def loadVspPhotometryData(self, objectName, ra = None, dec = None, fov = defFov):
        """
        star*: name of the star to plot. You must provide EITHER the star parameter OR ra and dec parameters (see below)
        ra*: right ascension
        dec*: declination
        fov*: field of view, in arcminutes.
        maglimit*: magnitude limit for the chart; stars fainter than the maglimit will not be plotted.
        other: other variable stars to mark on the chart. Omit or leave blank for none, "gcvs" for GCVS variables, "all" for all variables
        """
        if objectName is not None:
            ids = "star=%s" % (urllib.parse.quote_plus(objectName))
        else:
            ids = "ra=%s&dec=%s" % (pm.deg2hexa(ra/15.0), pm.deg2hexa(dec))
        vspUrl = f"{AAVSO_VSP_URL}/api/chart/?format=json&%s&fov=%d&maglimit=%4.1f&other=all" % (ids, fov, self.defMgLimit)
        f = request.urlopen(vspUrl)
        respJson = f.read().decode('UTF-8')

        resp = json.loads(respJson)

        self.chartId = resp['chartid']

        for pmr in resp['photometry']:

            # 'V', 'B', 'Rc', 'Ic', 'J', 'H', 'K'
            bands = {}
            for b in pmr['bands']:
                bands[b['band']] = b

            auid = pmr['auid']
            role = "C"
            ra = pmr['ra']
            dec = pmr['dec']
            raDeg = "%10.8f" % (pm.hexa2deg(ra) * 15.0)
            if not dec.startswith('-'):
                dec = "+" + dec
            decDeg = "%+10.8f" % (pm.hexa2deg(dec))
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
            label = str(pmr['label']).replace(' ', '_')


            self.writeRecord(auid, role, ra, raDeg, dec, decDeg, mgB, mgerrB, mgV, mgerrV, mgR, mgerrR, label)

            self.xmatchTable.add_row([auid, raDeg, decDeg, label])

        print("%d comparision stars were found in VSP database" % (len(resp['photometry'])))

    def dismissSameAUID(self, auid):
        for j in range(len(self.cache)):
            s = self.cache[j]
            if s.startswith(auid):
                self.cache[j] = '#' + s
                break

    def loadVsxCatalogData(self, objectName, ra = None, dec = None, fov = defFov, auidOnly = True):
        """
        &coords   The central J2000 RA/DEC coordinates for a radius search, expressed sexagesimally (by default), or in decimal degrees if format is set to d. Northern hemisphere coordinates must have the plus () sign URL-encoded as %2B. Space characters between all other figures must be replaced with the URL whitespace character (). The order argument (which see) must also be included in the query string with its value set to 9 in order to prompt VSX to display distances from the central coordinates in the results listing. Default is empty string (no radius search).
        &ident    Object identification for name searches. Space characters must be replaced with the URL whitespace character (+). Other special characters may also need to be URL-encoded. Default is empty string (no name search).
        &format   Explicit specification for format of coords. For sexagesimal, this value should be s. For decimal degrees, this value should be d. Default is s (sexagesimal).
        &geom     The geometry for central coordinate-based searches. For radius searches, this value should be r. For box searches, this value should be b. Default is r (radius search).
        &size     For box searches (geom=b), the width of the box. For radius searches (geom=r), the radius of the circle. Expressed in the units specified by unit (see next). Default is 10.0.
        &unit     The unit of measurement used for the value given in size (see above). For arc degrees, this value should be 1. For arc minutes, this value should be 2. For arc seconds, this value should be 3. Default is 2 (arc minutes).
        """
        if objectName is not None:
            # https://www.aavso.org/vsx/index.php?view=results.csv&ident=R+Car
            vsxUrl = f"{AAVSO_VSX_URL}/index.php?view=api.object&format=json&ident={urllib.parse.quote_plus(objectName)}"
            # ids = "star=%s" % (urllib.parse.quote_plus(objectName))
            f = request.urlopen(vsxUrl)
            respJson = f.read().decode('UTF-8')
            resp = json.loads(respJson)

            sra = pm.deg2hexa(float(resp['VSXObject']['RA2000']) / 15.0)
            sdec = pm.deg2hexa(float(resp['VSXObject']['Declination2000']))
        else:
            sra = pm.deg2hexa(ra / 15.0)
            sdec = pm.deg2hexa(dec)

        if not sdec.startswith('-') and not sdec.startswith('+'):
            sdec = "+" + sdec
        coords = sra + sdec

        vsxUrl = f"{AAVSO_VSX_URL}/index.php?view=query.votable&format=d&coords=%s&size=%d&unit=2&geom=b" % (coords, fov)
        # https://www.aavso.org/vsx/index.php?view=query.votable&coords=19:54:17+32:13:08&format=d&size=60&unit=2
        f = request.urlopen(vsxUrl)
        respVOTable = f.read().decode('UTF-8')
        votable = xmltodict.parse(respVOTable)
        tableData = votable['VOTABLE']['RESOURCE']['TABLE']['DATA']['TABLEDATA']
        if tableData is None:
            pm.printWarning('No VSX data for this field')
            return
        trs = tableData['TR']
        nrUnknown = 1
        if not isinstance(trs, list):
            trs = [trs]
        for tr in trs:
            # [None, 'PS1-3PI J185203.12+325154.5', 'Lyr', '283.01304000,32.86516000', 'RRAB', '16.810', 'r', '0.400', 'r', '57000.85600', None, '0.839122', None, None, None]
            auid = tr['TD'][0]
            if auid is None:
                if auidOnly:
                    continue
                auid = '999-VAR-%03d' % nrUnknown
                nrUnknown = nrUnknown + 1
                tr['TD'][0] = auid
            role = 'V'
            raDeg, decDeg = tr['TD'][3].split(',')
            ra = pm.deg2hexa(float(raDeg) / 15.0)
            dec = pm.deg2hexa(float(decDeg))
            if not dec.startswith('-'):
                dec = "+" + dec
                decDeg = "+" + decDeg
            label = tr['TD'][1].replace(' ', '_')

            self.dismissSameAUID(auid)

            self.writeRecord(auid, role, ra, raDeg, dec, decDeg, '-', '-', '-', '-', '-', '-', label)

            self.xmatchTable.add_row([auid, raDeg, decDeg, label])

        nrUnknown = 1
        count = 0
        for tr in trs:

            auid = tr['TD'][0]
            if auid is None:
                if auidOnly:
                    continue
                auid = f'999-VAR-{nrUnknown:03d}'
                nrUnknown = nrUnknown + 1
            # role = 'V'
            raDeg, decDeg = tr['TD'][3].split(',')
            ra = pm.deg2hexa(float(raDeg) / 15.0)
            dec = pm.deg2hexa(float(decDeg))
            label = tr['TD'][1]
            constell = tr['TD'][2]
            vartype = tr['TD'][4]
            if vartype is None:
                vartype = "-"
            period = tr['TD'][11]
            maxmg = tr['TD'][5]
            if tr['TD'][6] is not None:
                maxmg += tr['TD'][6]
            minmg = tr['TD'][7]
            if tr['TD'][8] is not None:
                minmg += tr['TD'][8]

            s = "#" + label.ljust(29) + auid.ljust(14) + ra.ljust(15) + dec.ljust(14) + constell.ljust(4) + vartype.ljust(11) + period.ljust(11) + maxmg + " - " + minmg
            self.cache.append(s)
            count = count + 1
        print(f"{count} variable stars were found in VSX database")

    def loadStdFieldData(self):
        saName = self.opt['stdFieldName']
        configFolder = pm.setup["CONFIG_FOLDER"]
        if configFolder is None:
            configFolder = "$HOME/.pmlib"

        # load stdArea file
        saFieldsFile = configFolder + "/landolt_fields.txt"
        if not exists(saFieldsFile):
            pm.printError(f"Landolt standard field catalog file {saFieldsFile} not found")
            return 'NoFile'
        saFields = Table.read(saFieldsFile, format = 'ascii')

        # find std area
        rows = list(filter(lambda r: r['FIELD_NAME'] == saName, saFields))
        if len(rows) == 0:
            pm.printError(f"No field name {saName} found in standard area file {saFieldsFile}.")
            return 'NoField'
        sa = rows[0]
        pm.printInfo(f"StdArea: {saName}, Coords: {sa['RA_2000']} {sa['DEC_2000']}, NumStars: {sa['NS']:d}, Comment: {sa['COMMENT']}")

        self.opt['ra'] = pm.hexa2deg(sa['RA_2000']) * 15.0
        self.opt['dec'] = pm.hexa2deg(sa['DEC_2000'])
        self.chartId = saName

        # load stdStars file
        saStarsFile = configFolder + "/landolt_stars.txt"
        if not exists(saStarsFile):
            pm.printError(f"Landolt standard stars catalog file {saStarsFile} not found")
            return 'NoFile'
        saStars = Table.read(saStarsFile, format = 'ascii')

        # find std stars
        stars = list(filter(lambda r: r['FIELD_NAME'] == saName, saStars))
        if len(stars) == 0:
            pm.printError("No stars for standard field %s in standard area file %s." % (saName, saStarsFile))
            return 'NoStars'

        # write stars to catalog file
        nrUnknown = 1
        for star in stars:
            auid = f'999-STD-{nrUnknown:03d}'
            nrUnknown = nrUnknown + 1
            role = "C"
            ra = star['RA_2000']
            raDeg = "%10.8f" % (pm.hexa2deg(ra) * 15.0)
            dec = star['DEC_2000']
            decDeg = f"{pm.hexa2deg(dec):+10.8f}"
            mv = float(star['MAG_V'])
            ev = float(star['ERR_V'])
            magB = "%+6.4f" % (float(star['MAG_BV']) + mv)
            errB = "%+7.4f" % (max(ev, float(star['ERR_BV'])))
            magV = "%+6.4f" % mv
            errV = "%+7.4f" % ev
            magR = "%+6.4f" % (mv - float(star['MAG_VR']))
            errR = "%+7.4f" % (max(ev, float(star['ERR_VR'])))
            label = saName + ":" + star['STAR']

            self.writeRecord(auid, role, ra, raDeg, dec, decDeg, magB, errB, magV, errV, magR, errR, label)

            self.xmatchTable.add_row([auid, raDeg, decDeg, label])

        print("%d standard star found for standard area %s" % (len(stars), saName))
        return None


    def loadFieldStars(self, target, fieldSize, mgLimit):
        v = VizierQuery(UCAC4, mgLimit)
        t = v.query(target, fieldSize)
        xt = v.xmatch(self.xmatchTable, 'RA_DEG', 'DEC_DEG')
        if t is None or xt is None:
            pm.printError("Accessing Vizier service if failed.")
            return

        matchCount = 0
        for row in t:
            auid = row['AUID']
            label = row['LABEL']
            mask = xt[UCAC4.cols['LABEL']] == row['LABEL'][6:]
            f = xt[mask]
            if len(f) > 0:
                auid = '#' + f[0]['AUID']
                label = label + '(' + f[0]['LABEL'] + ')'
                matchCount += 1

            self.writeRecord(auid, row['ROLE'], row['RA'], row['RA_DEG'], row['DEC'], row['DEC_DEG'],
                             row['MAG_B'], row['ERR_B'], row['MAG_V'], row['ERR_V'], row['MAG_R'], row['ERR_R'], label)

        print("%d field stars found in the %s catalog ; %d matched with comps and vars" % (len(t), v.catName, matchCount))

    def writeRecord(self, auid, role, ra, raDeg, dec, decDeg, magB, errB, magV, errV, magR, errR, label):
        s = auid.ljust(13) + role.ljust(5) + ra.ljust(15) + raDeg.ljust(15) + dec.ljust(15) + decDeg.ljust(15) + \
            magB.ljust(10) + errB.ljust(10) + magV.ljust(10) + errV.ljust(10) + magR.ljust(10) + errR.ljust(10) + label
        self.cache.append(s)

    def getPair(self, s, delim = ','):
        ss = s.split(delim)
        return [ss[0].strip(), ss[1].strip()]

    def determineCoordsFromImage(self, imageFileName):
        configFolder = pm.setup["CONFIG_FOLDER"]
        if configFolder is None:
            configFolder = "$HOME/.pmlib"
        SOLVE_ARGS = "-O --config " + configFolder + "/astrometry.cfg --use-sextractor --sextractor-path sextractor -r -y -p"

        makedirs("temp", exist_ok = True)

        pm.invoke("cp " + imageFileName + " temp/src.fits")

        astFolder = pm.setup["AST_BIN_FOLDER"]
        if astFolder is None:
            astFolder = "/usr/local/astrometry/bin"

        astResult = pm.invoke(astFolder + "/solve-field " + SOLVE_ARGS + " -D temp -N temp/ast.fits temp/src.fits")

        r = astResult.split('\n')
        coords = [0.0, 0.0]
        # sCoordsSexa = ["00:00:00", "+00:00:00"]
        for line in r:
            if line.startswith("Field center"):
                if line.find("RA,Dec") > -1:
                    sCoords = self.getPair(line.split('(')[2].split(')')[0])
                    coords = [ float(sCoords[0]), float(sCoords[1]) ]
                    pm.printInfo("Image center: %s %s" % (sCoords[0], sCoords[1]))
                else:
                    sCoordsSexa = self.getPair(line.split('(')[2].split(')')[0])
                    coords = [ pm.hexa2deg(sCoordsSexa[0]) * 15.0, pm.hexa2deg(sCoordsSexa[1]) ]
                    pm.printInfo("Image center: %s %s" % (sCoordsSexa[0], sCoordsSexa[1]))
            elif line.startswith("Field size"):
                sFieldSize = self.getPair(line.split(':')[1].split('a')[0], 'x')
                pm.printInfo("Image size: %s' x %s'" % (sFieldSize[0], sFieldSize[1]))

        return coords

    def loadOrigRefcat(self, fileName):
        self.origRefcat = Table.read(fileName, format='ascii')

    def updateUserModifications(self):
        # commented records
        for rec in [ x for x in self.origRefcat.meta['comments'] if x.startswith('000-') and x[4:7] != 'FFF' ]:
            r = rec.split()
            auid = r[0]
            label = r[-1]
            for j in range(len(self.cache)):
                if self.cache[j].startswith(auid):
                    self.cache[j] = '#' + self.cache[j]
                    print(f"Update: preserve user disabled object {auid} {label}")
                    break

        # user defined variables
        var_mask = self.origRefcat['ROLE'] == 'V'
        var_list = self.origRefcat[var_mask]
        for rec in var_list:
            auid = rec['AUID']
            label = rec['LABEL']
            missing = True
            for j in range(len(self.cache)):
                if self.cache[j][13] == 'V' and self.cache[j].endswith(label):
                    missing = False
                    break
            if missing:
                self.writeRecord(auid, rec['ROLE'], rec['RA'], "%10.8f" % (rec['RA_DEG']), rec['DEC'], "%+10.8f" % (rec['DEC_DEG']),
                                 "-", "-", "-", "-", "-", "-", label)
                for cmt in [ x for x in self.origRefcat.meta['comments'] if x.startswith(label) ]:
                    self.cache.append(f"#{cmt}")
                print(f"Update: preserve user added variable {auid} {label}")
        

    def execute(self):
        if self.opt['mgLimit'] is None:
            ml = float(pm.setup['DEF_FIELD_STAR_MG_LIMIT'])
            self.opt['mgLimit'] = ml if 0.0 < ml < 25.0 else 17.0

        if self.opt['image'] is not None:
            coords = self.determineCoordsFromImage (self.opt['image'])
            self.opt['ra'] = coords[0]
            self.opt['dec'] = coords[1]

        outFileName = pm.setup['CONFIG_FOLDER']
        if not outFileName.endswith('/'):
            outFileName += '/'
        outFileName += 'cat/'
        if not exists(outFileName):
            makedirs(outFileName)

        if self.opt['object'] is not None:
            outFileName += self.opt['object'].lower().replace(' ', '_') + '.cat'
        elif self.opt['stdFieldName'] is not None:
            outFileName += self.opt['stdFieldName'].lower().replace(' ', '_') + '.cat'
        else:
            s_ra = pm.deg2hexa(self.opt['ra'] / 15.0)
            s_dec = pm.deg2hexa(float(self.opt['dec']))
            if not s_dec.startswith('+') and not s_dec.startswith('-'):
                s_dec = '+' + s_dec
            outFileName += s_ra.replace(':', '')[:4] + s_dec.replace(':', '')[:5] + '.cat'

        # save original Refcat table for update
        if exists(outFileName) and self.opt['update']:
            self.loadOrigRefcat(outFileName)

        if not exists(outFileName) or self.opt['overwrite']:

            self.writeRecord('AUID', 'ROLE', 'RA', 'RA_DEG', 'DEC', 'DEC_DEG', 'MAG_B', 'ERR_B', 'MAG_V', 'ERR_V', 'MAG_R', 'ERR_R', 'LABEL')

            if self.opt['object'] is not None:

                self.loadVspPhotometryData(self.opt['object'], fov = self.opt['field'])

                self.loadVsxCatalogData(self.opt['object'], fov = self.opt['field'], auidOnly = self.opt['auidOnly'])

            elif self.opt['stdFieldName'] is not None:

                error = self.loadStdFieldData()

                if not error:
                    self.loadVsxCatalogData(None, ra = self.opt['ra'], dec = self.opt['dec'], fov = self.opt['field'], auidOnly = self.opt['auidOnly'])

            else:

                self.loadVspPhotometryData(None, ra = self.opt['ra'], dec = self.opt['dec'], fov = self.opt['field'])

                self.loadVsxCatalogData(None, ra = self.opt['ra'], dec = self.opt['dec'], fov = self.opt['field'], auidOnly = self.opt['auidOnly'])

            # collect field stars
            if self.opt['fieldStars']:

                if self.opt['object'] is not None:
                    target = self.opt['object']
                else:
                    target = SkyCoord(ra=self.opt['ra'], dec=self.opt['dec'], unit=(u.deg, u.deg), frame='icrs')

                self.loadFieldStars(target, self.opt['field'], self.opt['mgLimit'])

            # update user modifications from original Refcat table
            if self.origRefcat:
                self.updateUserModifications()

            # write cataglog to file
            outFile = open(outFileName, 'w')
            for s in self.cache:
                outFile.write(s + '\n')
            outFile.write('### chartid: ' + self.chartId + '\n')
            outFile.write('### fov: ' + str(self.opt['field']) + ' arcmin\n')
            outFile.close()

            pm.printInfo(f"Reference catalog file {outFileName} created.")

        else:

            pm.printInfo(f"Reference catalog file {outFileName} is exists.")

        if self.opt['folder']:

            # create link to refcat in given folder
            for f in glob(self.opt['folder'] + '*.cat'):
                remove(f)
            symlink(outFileName, self.opt['folder'] + 'ref.cat')


class MainApp:

    appName = 'ppl-refcat'
    appVersion = '1.2.0'
    appDescription = 'Create reference catalog for photometry.'

    opt = {
        'coords'      : None,           # coordinates of reference field
        'source'      : None,           # source catalog of field stars
        'object'      : None,           # object (variable star) name
        'stdFieldName': None,           # standard area name
        'image'       : None,           # image file
        'field'       : RefCat.defFov,  # reference field size in arcmins
        'auidOnly'    : True,           # variables having auid only
        'overwrite'   : False,          # overwrite catalog
        'folder'      : None,           # reference catalog file name
        'fieldStars'  : False,          # collect field stars into the catalog
        'mgLimit'     : None,           # mg limit for field star selection
        'update'      : False,          # update catalog
        }

    def __init__(self, argv):
        self.argv = argv

    def printVersion(self):
        print(pm.BGreen + self.appName + ", version " + self.appVersion + pm.Color_Off)

    def printTitle(self):
        print()
        self.printVersion()
        print(pm.Blue + self.appDescription + pm.Color_Off)
        print()

    def usage(self):
        print("Usage: " + self.appName + " [OPTIONS]... FOLDER_NAME")
        print()
        print("Mandatory arguments to long options are mandatory for short options too.")
        print("  -o,  --object object_name     object (variable star) name")
        print("  -c,  --coords ra,decl         coordinates of the center of reference frame, valid format is 12:34:56.7,-12:34:56.7")
        print("  -n,  --field-name field_name  standard field name")
        print("  -i,  --image filename         image file name")
        print("  -f,  --field size             field size in arcmin, default is 60 arcmin")
        print("  -a,  --all                    collect all variables ; if not set, collect variables having AUID only")
        print("  -r,  --field-stars            collect field stars")
        print("  -s,  --source catalog         source catalog for field stars ; default catalog is UCAC-4")
        print("  -l,  --limit mg               magnitude limit for field star selection")
        print("  -u,  --update                 update catalog file, if exists")
        print("  -w,  --overwrite              overwrite catalog file, if exists")
        print("  -h,  --help                   print this page")

    def processCommands(self):

        try:
            optlist, args = getopt (argv[1:], "ac:s:o:n:i:f:rl:uwh", ['all', 'coord=', 'source=', 'object=', 'field-name=', 'image=', 'field=', 'field-stars', 'limit=', 'update', 'overwrite', 'help'])
        except GetoptError:
            print ('Invalid command line options')
            exit(1)

        fieldStarsIndicated = False
        for o, a in optlist:
            if a[:1] == ':':
                a = a[1:]
            elif o == '-c' or o == '--coord':
                self.opt['coords'] = a
            elif o == '-s' or o == '--source':
                self.opt['source'] = a
                fieldStarsIndicated = True
            elif o == '-o' or o == '--object':
                self.opt['object'] = a.replace('_', ' ')
            elif o == '-n' or o == '--field-name':
                self.opt['stdFieldName'] = a
            elif o == '-i' or o == '--image':
                self.opt['image'] = a
            elif o == '-f' or o == '--field':
                if a.isdigit():
                    self.opt['field'] = int(a)
                else:
                    pm.printError(f"Invalid field size parameter: {a}. Use default 60 arcmin instead.")
            elif o == '-a' or o == '--all':
                self.opt['auidOnly'] = False
            elif o == '-r' or o == '--field-stars':
                self.opt['fieldStars'] = True
            elif o == '-l' or o == '--limit':
                err = False
                try:
                    ml = float(a)
                    if ml > 0.0 or ml <= 25.0:
                        self.opt['mgLimit'] = float(a)
                    else:
                        err = True
                except ValueError:
                    err = True
                if err:
                   pm.printError(f"Invalid mg limit: {a}. Use default instead.")
                fieldStarsIndicated = True
            elif o == '-u' or o == '--update':
                self.opt['update'] = True
                self.opt['overwrite'] = True
            elif o == '-w' or o == '--overwrite':
                self.opt['overwrite'] = True
            elif o == '-h' or o == '--help':
                self.usage()
                exit(0)

        if len(args) > 0:
            if not isdir(args[0]):
                pm.printError('%s is not a folder.' % (args[0]))
                exit(1)
            self.opt['folder'] = args[0]
            if not self.opt['folder'].endswith('/'):
                self.opt['folder'] += '/'

        if self.opt['object'] is None and self.opt['coords'] is None and self.opt['stdFieldName'] is None and self.opt['image'] is None:
            pm.printError("Either object name (-o), coordinates (-c), standard field name (-n) or image file (-i) must be given.")
            exit(1)

        if self.opt['coords'] is not None:
            c = self.opt['coords'].split(',', maxsplit = 2)
            ok = True
            if len(c) != 2:
                ok = False
            else:
                if re.fullmatch("(\\d){2}:(\\d){2}:(\\d){2}(\\.\\d)*", c[0]) is None:
                    ok = False
                if re.fullmatch("[+-]*(\\d){2}:(\\d){2}:(\\d){2}(\\.\\d)*", c[1]) is None:
                    ok = False
            if not ok:
                pm.printError("Invalid coordinate format")
                exit(1)

            self.opt['ra'] = pm.hexa2deg(c[0]) * 15.0
            self.opt['dec'] = pm.hexa2deg(c[1])

        if not self.opt['fieldStars'] and fieldStarsIndicated:
            pm.printError("Options -s or -l indicates collecting field stars, but -r option was not given. These options will be ignored.")


    def run(self):
        self.printTitle()
        self.processCommands()
        pm.saveCommand(self.opt['folder'], self.argv, 'refcat')

        start = datetime.now()

        cat = RefCat(self.opt)
        cat.execute()

        exectime = (datetime.now() - start).total_seconds()
        print("%sexecution time was %d seconds.%s" % (pm.Blue, exectime, pm.Color_Off))


if __name__ == '__main__':

    app = MainApp(argv)
    app.run()

# end main.
