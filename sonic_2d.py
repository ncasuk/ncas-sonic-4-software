#!/usr/bin/env python
# vim: set fileencoding=utf-8 expandtab ts=4:
"""A set of functions to convert a Gill 2D Sonic datafile to NetCDF"""
import time
import csv
from io import StringIO
import pandas as pd
import numpy as np
import glob
from datetime import datetime, timedelta
from os.path import join, getmtime, basename

from netCDF4 import Dataset
from pandas.io.parsers import read_csv, read_fwf
from pandas.tseries.offsets import DateOffset



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
    """Takes Gill Windsonic 2D datafile(s), prefixed by ISO-format system date, timezone, and timestamp of reading, and returns a Pandas dataframe containing T, U, V, r, θ.
    
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
        :return: a Pandas DataFrame of the sonic data
    """

    sonic = pd.DataFrame()
    for infile in infiles:
        with open(infile,'rb') as f:
            """pre-process to remove dodgy lines. Good lines are 77 chrs long and end in M,00,^C """
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

def read_amf_variables(csv_var_file):
    """
    Reads an AMF data project CSV-format variable list into a structure.
    """
    out = {}
    with open(csv_var_file,'r') as f:
        varfile = csv.DictReader(f)
        for line in varfile:
            if len(line['Variable']) >0:
                out[line['Variable']] = {}
                current_var = line['Variable']
            else:
                out[current_var][line['Attribute']] = line['Value']

    return out


def sonic_netcdf(sonic, output_file = "sonic_2d_data.nc"):
    """
    Takes a DataFrame with 2D sonic data and outputs a well-formed NetCDF
    using appropriate conventions.

    :param sonic: DataFrame with a Pandas DateTime index, U / ms¯¹, V / ms¯¹, wind speed / ms¯¹, and wind direction / degrees

    Test plot with (e.g.) ``cis plot wind_speed:sonic_2d_data.nc``
    """

    #instantiate NetCDF output
    dataset = Dataset(output_file, "w", format="NETCDF4_CLASSIC")

    # Create the time dimension - with unlimited length
    time_dim = dataset.createDimension("time", None)

    # Create the time variable
    base_time = sonic.index[0]
    sonic['timeoffsets'] = (sonic.index - base_time).total_seconds()

    time_units = "seconds since " + base_time.strftime('%Y-%m-%d %H:%M:%S')
    time_var = dataset.createVariable("time", np.float64, ("time",))
    time_var.units = time_units
    time_var.standard_name = "time"
    time_var.calendar = "standard"
    time_var[:] = sonic.timeoffsets.values

    #get variable descriptions
    amfvars = read_amf_variables("mean-winds.xlsx - Variables - Specific.csv")

    tempvar = {}
    #Create wind speed and wind direction vars
    for each in ['wind_speed','wind_from_direction','eastward_wind','northward_wind']:  
        tempvar[each] = dataset.createVariable(amfvars[each]['name'], amfvars[each]['type'], (amfvars[each]['dimension'],))
        tempvar[each].long_name = amfvars[each]['long_name']
        tempvar[each].units = amfvars[each]['units']
        tempvar[each].standard_name = amfvars[each]['standard_name']

    tempvar['wind_speed'][:] = sonic.r.values
    tempvar['wind_from_direction'][:] = sonic.theta.values
    tempvar['eastward_wind'][:] = sonic.U.values
    tempvar['northward_wind'][:] = sonic.V.values

    #  Set   the   global   attributes
    dataset.Conventions  =  "CF-1.6" 
    dataset.institution  =  "NCAS"   
    dataset.title  =  "2D Sonic NetCDF file" 
    dataset.history = "%s:  Written  with  script:  sonic_2d.py" % (datetime.now().strftime("%x  %X"))

    dataset.close()

def arguments():
    """
    Processes command-line arguments, returns parser.
    """
    from argparse import ArgumentParser
    parser=ArgumentParser()
    parser.add_argument('--outfile', dest="output_file", help="NetCDF output filename", default='sonic_2d_data.nc')
    parser.add_argument('infiles',nargs='+', help="Gill 2D Windsonic data files" )

    return parser

if __name__ == '__main__':
    args = arguments.parse_args()
    sonic_netcdf(get_sonic_data(args.infiles), args.output_file)
