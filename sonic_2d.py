#!/usr/bin/env python
# vim: set fileencoding=utf-8 expandtab ts=4:
import time
import csv
import os
import subprocess 
from io import StringIO
import pandas as pd
import numpy as np
import glob
from datetime import datetime, timedelta

from os.path import join, getmtime, basename
from amfutils.instrument import AMFInstrument
from netCDF4 import Dataset
from pandas.io.parsers import read_csv, read_fwf
from pandas.tseries.offsets import DateOffset

class GillWindSonic(AMFInstrument):
    """
    class to convert a Gill 2D Sonic datafile to NetCDF
    """
    
    amf_variables_file = "mean-winds.xlsx - Variables - Specific.csv"
    progname = __file__

    def polar(self, x, y, deg=True): # radian if deg=False; degree if deg=True
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
    
    
    def get_sonic_data(self, infiles):
        """Takes Gill Windsonic 2D datafile(s), prefixed by ISO-format system date, timezone, and timestamp of reading, and returns a Pandas dataframe containing time, U, V, r, θ.
        
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
    
                (sonic_in['r'], sonic_in['theta']) = self.polar(sonic_in.U,sonic_in.V)
    
                #if successful, add to output
                if isinstance(sonic_in, pd.DataFrame):
                    sonic = pd.concat([sonic,sonic_in])

        #set start and end times
        self.time_coverage_start = sonic.index[0].strftime(self.timeformat)
        self.time_coverage_end = sonic.index[-1].strftime(self.timeformat)

        return sonic
    
    def read_dataset_attributes(self, comvarfile):
        """
        Reads a csv file of the form 
           attrib1name,attrib1value
           attrib2name,attrib2value
           ...
           attribNname,attribNvalue
    
        and returns a dict of the results
            :param comvarfile: CSV file of attribute/value pairs
            :return: dictionary of attribte/value pairs
    
        """
        out= {}
        with open(comvarfile, 'r') as g:
            comvarfile_reader = csv.reader(g)
            for line in comvarfile_reader:
                if len(line) == 2:
                    out[line[0]] = line[1]
    
        return out
    
    def sonic_netcdf(self, sonic, output_dir="./netcdf/", metadata="2d-sonic-metadata"):
    
        """
        Takes a DataFrame with 2D sonic data and outputs a well-formed NetCDF
        using appropriate conventions.
    
        :param sonic: DataFrame with a Pandas DateTime index, U / ms¯¹, V / ms¯¹, wind speed / ms¯¹, and wind direction / degrees
    
        Test plot with (e.g.) ``cis plot wind_speed:sonic_2d_data.nc``
        """
    
        #instantiate NetCDF output
        self.dataset = Dataset(os.path.join(output_dir, self.filename("mean-winds","1")), "w", format="NETCDF4_CLASSIC")
    
        # Create the time dimension - with unlimited length
        time_dim = self.dataset.createDimension("time", None)
    
        # Create the time variable
        base_time = datetime(1970,1,1,0,0,0)
        sonic['timeoffsets'] = (sonic.index - base_time).total_seconds()
    
        #create the location dimensions - length 1 for stationary devices
        lat  = self.dataset.createDimension('latitude', 1)
        lon  = self.dataset.createDimension('longitude', 1)
    
        #create the location variables
        latitudes = self.dataset.createVariable('latitude', np .float32,  ('latitude',))
        latitudes.units = 'degrees_north'
        latitudes.standard_name = 'latitude'
        latitudes.long_name = 'Latitude'
    
        longitudes = self.dataset.createVariable('longitude', np .float32,  ('longitude',))
        longitudes.units = 'degrees_east'
        longitudes.standard_name = 'longitude'
        longitudes.long_name = 'Longitude'
    
        time_units = "seconds since " + base_time.strftime('%Y-%m-%d %H:%M:%S')
        time_var = self.dataset.createVariable("time", np.float64, dimensions=("time",))
        time_var.units = time_units
        time_var.axis = 'T'
        time_var.standard_name = "time"
        time_var.long_name = "Time (%s)" % time_units
        time_var.calendar = "standard"
        time_var.type = "float64"
        time_var.dimension = "time"
        time_var[:] = sonic.timeoffsets.values
    
        
    
        longitudes[:] = [self.raw_metadata['platform_longitude']]
        latitudes[:] = [self.raw_metadata['platform_latitude']]
    
        #remove lat/long
        self.raw_metadata.pop('platform_longitude',None)
        self.raw_metadata.pop('platform_latitude',None)
    
        tempvar = {}
        #Create wind speed and wind direction vars
        if 'r' in sonic:
            #r/θ polar windspeed/dir
            for each in ['wind_speed','wind_from_direction']:
                tempvar[each] = self.amf_var_to_netcdf_var(each)
            tempvar['wind_speed'][:] = sonic.r.values
            tempvar['wind_from_direction'][:] = sonic.theta.values
        if 'U' in sonic:
            #U/V Cartesian windspeed/dir
            for each in ['eastward_wind','northward_wind']:  
                tempvar[each] = self.amf_var_to_netcdf_var(each)
            tempvar['eastward_wind'][:] = sonic.U.values
            tempvar['northward_wind'][:] = sonic.V.values
    
        #  Set   the   global   attributes
        self.dataset.institution  =  "NCAS"   
        self.dataset.title  =  "2D Sonic NetCDF file" 
        self.dataset.history = "%s:  Written  with  script:  %s" % (datetime.now().strftime("%x  %X"),self.progname)
        self.dataset.processing_software_url = subprocess.check_output(["git", "remote", "-v"]).split()[1]#
        self.dataset.processing_software_url = self.dataset.processing_software_url.replace('git@github.com:','https://github.com/') # get the git repository URL
        self.dataset.processing_software_version = subprocess.check_output(['git','rev-parse', '--short', 'HEAD']).strip() #record the Git revision
        self.dataset.time_coverage_start = sonic.index[0].strftime('%Y-%m-%dT%H:%M:%S')
        self.dataset.time_coverage_end = sonic.index[-1].strftime('%Y-%m-%dT%H:%M:%S')
        #if latitudes.shape == (1,):
        #    dataset.geospatial_bounds = '('+latitudes[0].min().astype('str')+'N ' + longitudes[0].min().astype('str')+'E)'
        #else: #for future proofing, handles moving platform
        #    dataset.geospatial_bounds = '('+latitudes[:].min().astype('str')+'N ' + longitudes[:].min().astype('str')+'E, '+latitudes[:].max().astype('str')+'N ' + longitudes[:].max().astype('str')+'E)'
    
        #add all remaining attribs
        self.dataset.setncatts(self.raw_metadata)
    
        self.dataset.close()
    
if __name__ == '__main__':
    args = GillWindSonic.arguments().parse_args()
    sn = GillWindSonic(args.metadata)
   
    try:
        os.makedirs(args.outdir,mode=0o755)
    except OSError:
         #Dir already exists, probably
         pass
    else:
        print ("Successfully create dirctory %s" % args.outdir)
    sn.sonic_netcdf(sn.get_sonic_data(args.infiles), args.outdir, args.metadata)
