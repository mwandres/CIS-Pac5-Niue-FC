# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 13:12:34 2023

@author: antonioh
"""

import os
import datetime
import numpy as np
import pandas as pd
import datetime as dt
import netCDF4 as nc
import matplotlib.dates as mdates
# import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# import tide programs
import pyTMD.io
import pyTMD.time
import pyTMD.predict
import pyTMD.tools
import pyTMD.utilities

###############################################################################
def gen_tide(now):
    
    print('Generating tide from TPOX9 Global Tidal Model https://www.tpxo.net/global')
        
    folder_tmp ='../tmp/' 
    # tide_dir = '../extras/'
    
  
    
    # default coordinates to use
    LAT,LON = (-21.2100,200.225)
    
    # define time
    now_ = now - dt.timedelta(2)# NCEP needs at least 3 hours to upload their forecast from UTC00
    then_ = now_ + dt.timedelta(days=9.5,minutes=1)
    time_tide = mdates.drange(now_,then_,dt.timedelta(minutes=1))
    time_tide_hourly = mdates.drange(now_,then_,dt.timedelta(hours=1))
    date_time = pd.to_datetime(mdates.num2date(time_tide))
    date_time_hourly = pd.to_datetime(mdates.num2date(time_tide_hourly))
    time_tmd=time_tide-mdates.date2num(np.datetime64('1992-01-01T00:00:00')) #-- convert time from MJD to days relative to Jan 1, 1992 (48622 MJD)
    
   
    
    # OTIS_path = "../extras/TPXO9/"
    grid_file = "../extras/TPXO9/grid_tpxo9v2"
    h_file = "../extras/TPXO9/h_tpxo9.v2"
  
    
    # model =  pyTMD.io.OTIS.read_otis_grid(grid_file)
    


    amp,ph,D,c = pyTMD.io.extract_constants(np.atleast_1d(LON), np.atleast_1d(LAT),
                               grid_file,h_file,EPSG='4326',TYPE='z',METHOD='spline',GRID='OTIS',extrapolate = True,cutoff=50)
   
    DELTAT = np.zeros_like(time_tmd)
   
    # calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    # calculate constituent oscillation
    hc = amp*np.exp(cph)
    # predict tidal elevations at time 1 and infer minor corrections
    TIDE = pyTMD.predict.time_series(time_tmd, hc, c,
    deltat=DELTAT)
    MINOR = pyTMD.predict.infer_minor(time_tmd, hc, c,
    deltat=DELTAT)
    TIDE.data[:] += MINOR.data[:]
    # convert to centimeters
    TIDE.data[:] *= 100.0
    
    f = interp1d(time_tide,TIDE)
    TIDE_hourly = f(time_tide_hourly)
    

    #-- Save netcdf with the 10 days tide forecast hourly 
    fn = folder_tmp  + 'tide_hourly.nc'
    try:
        os.remove(fn)
    except:
        print('The system cannot find the file specified')
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    time = ds.createDimension('time', None)
    times = ds.createVariable('time', 'f8', ('time',))
    times.units='hours since 1950-01-01 00:00:00'
    times.calendar='gregorian'
    # times_min = ds.createVariable('time_min', 'f4', ('time',))
    # times_min.units='hours since 1950-01-01 00:00:00'
    # times_min.calendar='gregorian'
    tide= ds.createVariable('tide', 'f4', ('time',))
    tide.units = 'cm'
    # tide_min= ds.createVariable('tide_min', 'f4', ('time',))
    # tide_min.units = 'cm'
    times[:]=[nc.date2num(x,units=times.units,calendar=times.calendar) for x in date_time_hourly]
    # time_or = nc.num2date(times,units=times.units,calendar=times.calendar)
    tide[:]=TIDE_hourly
    # times_min[:]=[nc.date2num(x,units=times.units,calendar=times.calendar) for x in date_time]
    # tide_min[:]=TIDE.data
    ds.close()
    print('1 hour tides stored as ' + fn)
    
    #-- Save  tide forecast in netcdf in 1 minute and one hour resolution 
    fn = folder_tmp  + 'tide_minute.nc'
    try:
        os.remove(fn)
    except:
        print('The system cannot find the file specified')
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    time = ds.createDimension('time', None)
    times = ds.createVariable('time', 'f8', ('time',))
    times.units='hours since 1950-01-01 00:00:00'
    times.calendar='gregorian'
    tide= ds.createVariable('tide', 'f4', ('time',))
    tide.units = 'cm'
    times[:]=[nc.date2num(x,units=times.units,calendar=times.calendar) for x in date_time]
    tide[:]=TIDE.data
    ds.close()
    
    print('1 minute tides stored as ' + fn)
    
    
    
    # #-- differentiate to calculate high and low tides
    # diff = np.zeros_like(time_tide, dtype=np.float64)
    # #-- forward differentiation for starting point
    # diff[0] = TIDE.data[1] - TIDE.data[0]
    # #-- backward differentiation for end point
    # diff[-1] = TIDE.data[-1] - TIDE.data[-2]
    # #-- centered differentiation for all others
    # diff[1:-1] = (TIDE.data[2:] - TIDE.data[0:-2])/2.0
    # htindex, = np.nonzero((np.sign(diff[0:-1]) >= 0) & (np.sign(diff[1:]) < 0))
    # ltindex, = np.nonzero((np.sign(diff[0:-1]) <= 0) & (np.sign(diff[1:]) > 0))
    
    # #-- create plot with tidal displacements, high and low tides and dates
    # plt.figure(figsize=(15,5)) 
    # fig,ax1 = plt.subplots(num=1)
    # ax1.plot(date_time,TIDE.data,'k')
    # ax1.plot(date_time_hourly,TIDE_hourly,'k*')
    # ax1.plot(date_time[htindex],TIDE.data[htindex],'r*')
    # ax1.plot(date_time[ltindex],TIDE.data[ltindex],'b*')
    # ax1.plot(np.matlib.repmat(now,5,1),np.linspace(-100,100,5),'r')
    # plt.ylim((-100,100))
    # for h in range(0,len(date_time),60*24):
    #     ax1.axvline(date_time[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
    # ylbl = ax1.axes.get_yticks()
    # for h in range(0,len(ylbl)):
    #     ax1.axhline(ylbl[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
    # # ax1.set_xlim(0,7*24)
    # ax1.set_ylabel('{0} Tidal Displacement [cm]'.format('tpxo8'))
    # ax1.set_title(u'{0:0.6f}\u00b0N {1:0.6f}\u00b0W'.format(LAT,LON))
    # fig.subplots_adjust(left=0.10,right=0.98,bottom=0.30,top=0.7)
    # plt.show()
    