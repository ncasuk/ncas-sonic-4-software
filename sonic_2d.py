#!/usr/bin/env python
# vim: set fileencoding=utf-8 expandtab ts=4:
'''A set of functions to convert a Gill 2D Sonic datafile to NetCDF'''
from pandas.io.parsers import read_csv, read_fwf
from pandas.tseries.offsets import DateOffset

import pandas as pd
import numpy as np
import glob
from datetime import datetime
from os.path import join, getmtime, basename

from io import StringIO


def polar(x, y, deg=True): # radian if deg=False; degree if deg=True
    """ Convert from rectangular (x,y) to polar (r,w)
    r = sqrt(x^2 + y^2)
    w = arctan(y/x) = [-\pi,\pi] = [-180,180] 
    """ 
    if deg:
        r = np.hypot(x,y)
        #theta = np.degrees(np.arctan2(y, x))
        theta = 180.0 * np.arctan2(y, x) / np.pi
        normtheta = theta.where(theta > 0, theta + 360)
        return r, normtheta
    else:
        return np.hypot(x, y), np.arctan2(y, x)


def get_sonic_data(infiles):
    '''Takes Gill Windsonic 2D datafile(s), prefixed by ISO-format system date, timezone, and timestamp of reading, and returns a Pandas dataframe containing T, U, V, r, θ.
    
    Note that Gill Windsonics use unconventional defns for U & V. 
    (see http://novalynx.com/products/download/WindSonicWebManual.pdf )

    Example good dataline:

    ::

        2017-08-30 01:17:51,UTC,2017-08-30T01:17:52.906838 ^BQ,+002.03,+000.64,M,00,^C

    ``^B`` and ``^C`` here are ASCII "Start of Text" (STX) and "End of Text" (ETX)
    
    Incoming datafiles are stripped of datalines that:

        #. Are not 77 characters long
        #. Do not end in the Gill "working" statuscode
        #. Contain "odd" characters (i.e. are corrupted in some way)

        :param infiles: list(-like) of data filenames
    '''

    sonic = pd.DataFrame()
    for infile in infiles:
        with open(infile,'rb') as f:
            '''pre-process to remove dodgy lines. Good lines are 77 chrs long and end in M,00,^C '''
            okchars = dict.fromkeys("0123456789: \nUTCQM+,.-"+chr(2)+chr(3))

            memory_file = StringIO()
            for line in f:
                line_dec = line.decode('ISO-8859-1')
                if len(line) == 77 and line[-7:-1] == b'M,00,\x03' and all(c in okchars for c in line_dec): #00 is "no error" status
                    memory_file.write(line.decode('ISO-8859-1').replace(' \x02Q',''));

            #back to top of memory file
            memory_file.seek(0)

            sonic_in = read_csv(memory_file, header=None, usecols=[2,3,4], names=['DateTime','UGILL','VGILL'], parse_dates=['DateTime'], index_col='DateTime', dtype={'VGILL': np.float64, 'UGILL': np.float64}, error_bad_lines=False, encoding="ISO-8859-1", na_values={'UGILL':['M','Ó'],'VGILL':['M','Ó'],'DateTime':['M','Ó']})

            sonic_in['U'] = -sonic_in.UGILL
            sonic_in['V'] = sonic_in.VGILL
            #del sonic_in['VGILL']
            #del sonic_in['UGILL']

            (sonic_in['r'], sonic_in['theta']) = polar(sonic_in.U,sonic_in.V)

            #if successful, add to output
            if isinstance(sonic_in, pd.DataFrame):
                sonic = pd.concat([sonic,sonic_in])

    return sonic


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser()
    parser.add_argument('--dir', dest="DATADIR", help="base data directory (e.g. /data/2d-sonic/ )", default='.')
    parser.add_argument('infiles',nargs='+')

    args = parser.parse_args()
    DATADIR=args.DATADIR

    print(get_sonic_data(args.infiles).theta.mean())
