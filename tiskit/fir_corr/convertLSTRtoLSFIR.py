#!/usr/bin/env python3
""" run fir2caus.py on selected data
"""
import os
import re
import sys
# Gives access to my fir2caus module
sys.path.append('/Users/crawford/bin')

import obspy.core as obc
from obspy.core import UTCDateTime
import fir2caus

# How to reload fir2caus if debugging/modifying
# import imp; imp.reload(fir2caus)

basedir="/Volumes/Wayne_Data/SEISAN/LSTR/WAV/"
database='LSTR_'
firCorrFileName='/Users/crawford/soft/FIR_CORR/corr_coeffs/lc2000_fir3_0.json'
firDecimation=2
startdate=UTCDateTime('2015-04-21T00:00:00')
enddate=UTCDateTime('2016-05-28T00:00:00')

for ydir in os.listdir(os.path.join(basedir,database)) :
  ypath=os.path.join(basedir,database,ydir)
  if os.path.isdir(ypath):
    for mdir in os.listdir(ypath) :
      mpath=os.path.join(ypath,mdir)
      if os.path.isdir(mpath):
        for fn in os.listdir(os.path.join(basedir,database,ydir,mdir)) :
          if re.match('\d{4}-\d{2}-\d{2}-\d{4}-\d{2}M\..........',fn):
            filedate=UTCDateTime(fn[:18])
            if filedate>=startdate and filedate<=enddate:
                fpath=os.path.join(mpath,fn)
                print(fpath)  
                st = obc.read(fpath)
                st_corr=st.copy()
                # Set location code to 01, to distinguish from original
                for tr in st_corr:
                  tr.stats.location='01'
                st_corr = fir2caus.fir2caus(st_corr,firCorrFileName,firDecimation)
                newfn=fn.replace('LSTR_','LSFIR')
                st_corr.write(newfn, format='MSEED')    # Uses same encoding as input
                # print('Running date {}, to {}'.format(str(filedate),newfn))                

