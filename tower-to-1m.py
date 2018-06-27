#!/usr/bin/env python
# vim: set fileencoding=utf-8 expandtab ts=3:
'''Produce a 1m rolling average for the tower RHT/2D sonic pairs and 3d sonic
Outputs HTML.'''
from pandas.io.parsers import read_csv, read_fwf
from pandas.tseries.offsets import DateOffset
#from pandas.stats.moments import rolling_mean
import pandas as pd
import numpy as np
import glob
from datetime import datetime
from os.path import join, getmtime, basename

from optparse import OptionParser
from io import StringIO
from file_get_contents import file_get_contents

yuns=file_get_contents('/opt/scripts/tower-processing/yuns').strip().split("\n")

print(yuns)

#DATADIR='/data/shareddata/tower/'
#UNIT='000'
#DATE='2014-08-20'
parser=OptionParser()
parser.add_option('--dir', dest="DATADIR", help="base data directory (e.g. /data/shareddata/tower/ )", default='/data/shareddata/tower/')
parser.add_option('--outfile', dest="OUTFILE", help="Output HTML file", default='/var/www/html/vmast/index.html')
parser.add_option('--template', dest="TEMPLATE", help="Input HTML file", default='/opt/scripts/tower-processing/template')
#parser.add_option('-u','--unit', dest="UNIT", help="Unit number of the sensor device, e.g. 000, 002", default="000")
#parser.add_option('--day', dest="DATE", help="date in 10-character ISO style, e.g. 2014-08-20", default=date.today().isoformat())

(options, args) = parser.parse_args()

DATADIR=options.DATADIR
OUTFILE=options.OUTFILE
TEMPLATE=options.TEMPLATE
#DATE=options.DATE

#print DATADIR, UNIT, DATE

""" Convert from rectangular (x,y) to polar (r,w)
 r = sqrt(x^2 + y^2)
 w = arctan(y/x) = [-\pi,\pi] = [-180,180] 
""" 
def polar(x, y, deg=True): # radian if deg=0; degree if deg=1
   from math import hypot, atan2, pi
   if deg:
      r = hypot(x,y)
      theta = 180.0 * atan2(y, x) / pi
      if theta < 0:
         theta = theta +360
      return r, theta
   else:
      return hypot(x, y), atan2(y, x)


out = {}

#get T, RH, U, V from 2d sonics and RHT sensors
for UNIT in yuns:
   newest_RHT = max(glob.iglob(join(DATADIR,UNIT+'w','data','RHT_'+UNIT+'_*.txt')),key =getmtime)
   print(newest_RHT, UNIT)
   
   rht_in= read_csv(newest_RHT, header=None, names=['DateTime','Nonsense','TempSH','RH'], parse_dates=True, index_col=0)

   newest_sonic = max(glob.glob(join(DATADIR,UNIT+'w','data','2D-sonic-'+UNIT+'-*.tsv')))
   #print newest_sonic
   #names=['DateTime','node','U','V','Units','Status','Checksum']
   print(newest_sonic)
   sonic_in = read_csv(StringIO(file_get_contents(newest_sonic).replace("\t",",")), header=None,  names=['DateTime','node','UGILL','VGILL','Units','Status','Checksum'], parse_dates=True, index_col=0)
   #sonic_in = read_csv(newest_sonic, header=None,  names=['DateTime','node','U','V','Units','Status','Checksum'], parse_dates=True, index_col=0)
   #Gill Windsonics use...unconventional..defns for U & V... 
   #(see http://novalynx.com/products/download/WindSonicWebManual.pdf)
   sonic_in['U'] = -sonic_in.UGILL
   sonic_in['V'] = sonic_in.VGILL
   del sonic_in['VGILL']
   del sonic_in['UGILL']

   #out[UNIT] = rolling_mean(rht_in,60).last('5S')
   #out[UNIT] = rht_in.resample('60S', how='mean').last('60S')
 
   #combines previous data; will put NaN in if data are missing 
   out_df = pd.concat([rht_in.resample('60S').mean(), sonic_in.resample('60S').mean()], axis=1).last('60S')
   (r, theta) = polar(out_df.U[-1],out_df.V[-1])
   out[UNIT] = {
                  'Time':out_df.index[-1].isoformat(), 
                  'U':'%.1f' % out_df.U[-1], 
                  'V':'%.1f' % out_df.V[-1], 
                  'WindSpeed':'%.1f' % r,
                  'WindAngle':'%.1f' % theta,
                  #'TempTH':'%.1f' % out_df.TempTH[-1], 
                  'TempSH':'%.1f' % out_df.TempSH[-1], 
                  'RH':'%.1f' % out_df.RH[-1]
               }

'''#get 3D sonic data
newest_3d= max(glob.glob(join(DATADIR,'sonw','data','*.metek')))
#fwf stands for "fixed width fields" colspec is the defn
print newest_3d
colspec = [(5,11),(15,21),(25,31),(35,41)]
threeD_in=read_fwf(newest_3d, colspecs=colspec, header=None, names=['U','V','W','T'])#, converters = {0: int, 1: int})
basedate = datetime.strptime(basename(newest_3d)[0:13],'%y%m%d_%H%M%S')

rnf = pd.date_range(basedate, periods=len(threeD_in.index), freq=DateOffset(seconds=0.1))
#prepend timeseries
threeD_in.insert(0, 'Time', rnf)

threeD_in.U.astype(np.float64)
threeD_in.V.astype(np.float64)

#reset index
threeD_in.set_index('Time', inplace=True)

print threeD_in.resample('60S', how='mean').last('60S')'''

template = "".join(open(TEMPLATE,'r').readlines())
output = open(OUTFILE,'wb')
output.write(bytes(template.format(
   z5rh=out['005']['RH'], z5ws=out['005']['WindSpeed'], z5wd=out['005']['WindAngle'],z5temp2=out['005']['TempSH'],
   z4rh=out['004']['RH'], z4ws=out['004']['WindSpeed'], z4wd=out['004']['WindAngle'],z4temp2=out['004']['TempSH'],
   z3rh=out['003']['RH'], z3ws=out['003']['WindSpeed'], z3wd=out['003']['WindAngle'],z3temp2=out['003']['TempSH'],
   z2rh=out['002']['RH'], z2ws=out['002']['WindSpeed'], z2wd=out['002']['WindAngle'],z2temp2=out['002']['TempSH'],
   z1rh=out['001']['RH'], z1ws=out['001']['WindSpeed'], z1wd=out['001']['WindAngle'],z1temp2=out['001']['TempSH'],
   z0rh=out['000']['RH'], z0ws=out['000']['WindSpeed'], z0wd=out['000']['WindAngle'],z0temp2=out['000']['TempSH'],
),'UTF-8')

)
output.close()

