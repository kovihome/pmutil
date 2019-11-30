#!/usr/bin/env python3
#
# PmUtils/pmlt2ut
#
'''
Convert Local Time to UTC, format is both of them is ISO 8601
If no argument passed, the current UTC time will returns.

Created on Dec 1, 2019

@author: kovi
'''

import sys
import time
import datetime
from _datetime import timedelta

if __name__ == '__main__':

    if len(sys.argv) > 1:
        s = sys.argv[1]
        t = time.strptime(s, '%Y-%m-%dT%H:%M:%S')
        d = datetime.datetime(*time.gmtime(time.mktime(t))[:6])
    else:
        d = datetime.datetime.utcnow()
    if d.microsecond > 500000:
        td = timedelta(microseconds = 1000000 - d.microsecond)
        d = d + td
    else:
        td = timedelta(microseconds = d.microsecond)
        d = d - td
    print(d.isoformat())

# end __main__
