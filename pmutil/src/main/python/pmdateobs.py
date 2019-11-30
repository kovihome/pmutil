#!/usr/bin/env python3
#
# PmUtils/pmdateobs
#
'''
Calculate effective middle observation time, full exposure time

Created on Dec 1, 2019

@author: kovi
'''

import sys
import time
import datetime
import subprocess

def invoke(cmd):
    r = cmd.split()
    a = subprocess.check_output(r)
    return a.decode('ascii')[:-1]
    

if __name__ == '__main__':

    if len(sys.argv) > 3:
        target_fn = sys.argv[1]

        e_sum = 0
        m_sum = 0
        for source_fn in sys.argv[2:]:

            date_obs = invoke('fiheader --get DATE-OBS ' + source_fn)
            exptime = invoke('fiheader --get EXPTIME ' + source_fn)

            t = time.mktime(time.strptime(date_obs, '%Y-%m-%dT%H:%M:%S'))
            e = int(exptime)

            e_sum = e_sum + e
            m_sum = m_sum + (t + e / 2) * e

        m_avg = m_sum / e_sum
        d = datetime.datetime(*time.gmtime(m_avg)[:6])
        invoke('sethead -s _ ' + target_fn + ' EXPTIME=' + str(e_sum) + ' DATE-MID=' + d.isoformat() + ' / Effective_center_of_time_of_observation')
        print('observation center time:', d.isoformat(), 'cumulated exptime:', str(e_sum))

# end __main__

