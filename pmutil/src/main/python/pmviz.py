#!/usr/bin/env python3
#
# PmUtils/viz
#
'''
Created on Mar 03, 2021

@author: kovi
'''

#from sys import argv
#from getopt import getopt, GetoptError
#from os import remove
#from os.path import isdir, exists
#from glob import glob
#from datetime import datetime


from requests.exceptions import ConnectionError
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.coordinates import SkyCoord

from pmbase import hexa2deg, deg2hexa


class Refcat(Table):
    names=['AUID','ROLE','RA','RA_DEG','DEC','DEC_DEG','MAG_B','ERR_B','MAG_V','ERR_V','MAG_R','ERR_R','LABEL']
    dtype=['U11','U2','U12','U12','U12','U12','U6','U5','U6','U5','U6','U5','U50']
    def __init__(self, masked=None):
        Table.__init__(self, names=self.names, dtype=self.dtype, masked=masked)
    def fmg(self, mg, er):
        fmg = str("%-6.3f" % (float(mg))) if mg != '--' else '-'
        fer = ("%-5.2f" % (float(er) / 100.0)) if mg != '--' and er != '' else '-'
        return fmg,fer
    def load(self, table, abbrev, mapping, auid_start=1):
        n_auid = auid_start
        for row in table:
            v,ev = self.fmg(row[mapping['MAG_V']], row[mapping['ERR_V']])
            b,eb = self.fmg(row[mapping['MAG_B']], row[mapping['ERR_B']])
            r,er = self.fmg(row[mapping['MAG_R']], row[mapping['ERR_R']])
            fra = float(row[mapping['RA_DEG']])
            fdec = float(row[mapping['DEC_DEG']])
            ra = "%12.8f" % fra
            dec = "%+12.8f" % fdec
            sra = deg2hexa(fra / 15.0)
            sdec = deg2hexa(fdec)
            if not sdec.startswith('-') and not sdec.startswith('+'):
                sdec = "+" + sdec
            auid = '%03d-FFF-%03d' % (n_auid // 1000, n_auid % 1000)
            self.add_row([auid, 'F', sra, ra, sdec, dec, b, eb, v, ev, r, er, abbrev + '-' + row[mapping['LABEL']]])
            n_auid += 1

class VizUCAC4:
    cat = 'I/322A' # UCAC-4
    catFull = 'I/322A/out'
    catName = 'UCAC4'
    cs = {'LABEL':'UCAC4','RA_DEG':'RAJ2000','DEC_DEG':'DEJ2000','MAG_B':'Bmag','ERR_B':'e_Bmag','MAG_V':'Vmag','ERR_V':'e_Vmag','MAG_R':'rmag','ERR_R':'e_rmag'}

    def __init__(self, limit):
        slimit = "<" + ("%.1f" % limit)
        print(f'slimit: {slimit}')
        self.viz = Vizier(catalog=self.cat, columns=list(self.cs.values()), column_filters={'Vmag': slimit}, row_limit=-1)

    def query(self, target, size):
        t = Refcat()
        try:
            result = self.viz.query_region(target, width=size*u.arcmin)
        except ConnectionError:
            return None
        t.load(result[0], self.catName, self.cs)
        return t

    def xmatch(self, srcTable, raColName, decColName):
        try:
            return XMatch.query(cat1=srcTable, cat2='vizier:' + self.catFull, max_distance=5 * u.arcsec, \
                colRA1=raColName, colDec1=decColName, colRA2='RAJ2000', colDec2='DEJ2000')
        except ConnectionError:
            return None

if __name__ == '__main__':

    # M31

    sz  = 30

    v = VizUCAC4(17.0)

    # Test 1: query for object
    obj = "R Crb"
    print(f'Test 1: query for object {obj}')

    t = v.query(obj, sz)

    print(f'{len(t)} objects found')
    print(t[0:10])

    # Test 2: query for coordinates
    r = 10.6846
    d = 41.2694
    obj = SkyCoord(ra=r, dec=d, unit=(u.deg, u.deg), frame='icrs')
    print(f'Test 2: query for coordinates {obj.ra}, {obj.dec}')

    t = v.query(obj, sz)

    print(f'{len(t)} objects found')
    print(t[0:10])


    # Test 3: xmatch for two objects
    srcTable = Table(names=['AUID', 'RA_DEG', 'DEC_DEG'], dtype=['U16','U16','U16'])
    srcTable.add_row(['R Boo', '219.29825000', '+26.73658000'])
    srcTable.add_row(['ASASSN', '219.45444000', '+27.11915000'])
    print('Test 3: xmatch for two objects')
    print(srcTable)

    t = v.xmatch(srcTable, 'RA_DEG', 'DEC_DEG')

    print(f'{len(t)} objects found')
    print(t[0:10])

# end main.




