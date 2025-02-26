# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:06:28 2024

@author: antonioh
"""

import os
import netCDF4
import numpy as np
import matplotlib.dates as mdates
from scipy.interpolate import interp1d
import netCDF4 as nc
import pandas as pd
import datetime as dt

# Functions
def read_netcdf(nc_fname,dtt,length_dt):
# Function to read Sea Level Pressure  contained in the netCDF files from step1_download_NCEP.py


    nc_fnameX = nc_fname + '_' + "{0:0>3}".format(0) +'.nc'
    nc = netCDF4.Dataset(nc_fnameX)
    # kk = np.array(nc['ordered_sequence_of_data'])
    ln = np.array(nc['lon'])
    lt = np.array(nc['lat'])
   
    
    # lon,lat =np.meshgrid(ln,lt)
    
    # from matplotlib import pyplot as plt
    # plt.pcolor(lon,lat,slp[0,:,:], cmap='viridis')
    # plt.colorbar() 
    # plt.plot(lon[5,4],lat[5,4],'*')
    # plt.show()
    
    
   
    
    tt =  dt.datetime.strptime(nc['time'].units[11:-1],"%Y-%m-%dT%H:%M:%S") +  dt.timedelta(hours=int(nc['time'][:].data[0]))
    slp = np.nan_to_num(np.array(nc['Pressure_surface'])[:,5,4])# Position 5,4 is the closest point to the west of Rarotonga

   
   # read SLP from hour=3 to hour=180 every 3 hours
    for i in range(dtt,length_dt,dtt):
        nc_fnameX = nc_fname + '_' +  "{0:0>3}".format(i) +'.nc'
        nc = netCDF4.Dataset(nc_fnameX)
        tt_ = dt.datetime.strptime(nc['time'].units[11:-1],"%Y-%m-%dT%H:%M:%S") +  dt.timedelta(hours=int(nc['time'][:].data[0]))
        tt = np.append(tt,tt_)
        slp = np.nan_to_num(np.append(slp,np.array(nc['Pressure_surface'])[:,5,4],axis=0))
        print(nc_fnameX)
        
  

    return(tt,slp)



###############################################################################
def gen_IB(now):
    
    print('Proccessing Inverted Barometer Effect')  
    
    # define the output folder were SWAN will be run
    out_name = '../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H") 

    folder_tmp ='../tmp/' 
    # directory of the nc files
    slp_nc = '../tmp/slp_tmp'
    h_slp_nc = '../tmp/h_slp_tmp'
    hh_slp_nc = '../tmp/hh_slp_tmp'
     
    slp_dt = 3
    time_length = 181
    h_time_length = 22
     
    # read slp nc files
    tt,slp = read_netcdf(slp_nc,slp_dt,time_length) 
    h_tt,h_slp = read_netcdf(h_slp_nc,slp_dt,h_time_length) 
    hh_tt,hh_slp = read_netcdf(hh_slp_nc,slp_dt,h_time_length) 
    
    tt=np.concatenate((hh_tt,h_tt,tt), axis=0)
    slp=np.concatenate((hh_slp,h_slp,slp), axis=0)
    ss = (1013 - slp/100)# 0.01 m per Mb
        
     
        
    tt_h = np.array([pd.Timestamp(d).timestamp() for d in tt])
    dates = pd.date_range(start=tt[0], end=tt[-1], freq='h')
    tt_hour = np.array([pd.Timestamp(d).timestamp() for d in dates])
    f = interp1d(tt_h,ss)
    ss_hour = f(tt_hour)
        
     
    
    #-- Save netcdf with the 10 days IB forecast hourly 
    fn = folder_tmp  + 'IB_hourly.nc'
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
    IB= ds.createVariable('IB', 'f4', ('time',))
    IB.units = 'cm'
    # tide_min= ds.createVariable('tide_min', 'f4', ('time',))
    # tide_min.units = 'cm'
    times[:]=[nc.date2num(x,units=times.units,calendar=times.calendar) for x in dates]
    # time_or = nc.num2date(times,units=times.units,calendar=times.calendar)
    IB[:]=ss_hour
    # times_min[:]=[nc.date2num(x,units=times.units,calendar=times.calendar) for x in date_time]
    # tide_min[:]=TIDE.data
    ds.close()
    print('1 hour IB stored as ' + fn)
    
    
    
    dates = pd.date_range(start=tt[0], end=tt[-1], freq='min')
    tt_min = np.array([pd.Timestamp(d).timestamp() for d in dates])
    f = interp1d(tt_h,ss)
    ss_min = f(tt_min)
    
     
    #-- Save  IB forecast in netcdf in 1 minute
    fn = folder_tmp  + 'IB_minute.nc'
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
    times[:]=[nc.date2num(x,units=times.units,calendar=times.calendar) for x in dates]
    tide[:]=ss_min
    ds.close() 
    print('1 minute IB stored as ' + fn)
 
 