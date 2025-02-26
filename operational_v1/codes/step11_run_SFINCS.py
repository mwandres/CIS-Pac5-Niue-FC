# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:22:16 2024

@author: antonioh
"""
import os
import numpy as np
import xarray as xr
from scipy.signal import correlate, find_peaks
from scipy.interpolate import interp1d
import pandas as pd
import subprocess
import shutil

# import matplotlib.pyplot as plt


def load_netcdf(fpath):
       with xr.open_dataset(fpath) as fncdf:
           data = { key: item.data.squeeze() for key, item in fncdf.items() }
           coords = { key: item.data.squeeze() for key, item in fncdf.coords.items() }
       return data, coords


def generate_sfincs_forcing(Xb_path,index):

    # Storing the time serie of water level at SFINCS boundary, loop through all profiles
    tr_water_level = []
    for ii in range(len(index)):
        name = 'Prof_'+ str(ii+1) +'/xboutput.nc'
        fpath = os.path.join(Xb_path, name)
        data, coords = load_netcdf(fpath)
        tr_water_level.append(data['zs'][:,index[ii]-1] .copy())# indices calculated with matlab
    
    # Pre-processing the data for SFINCS
    tr_water_level = np.array(tr_water_level)
    fbzs = np.nanmean(tr_water_level, axis=1)
    fbzi = tr_water_level - fbzs[:,None]
    fbzi = fbzi.T
    

    ## Syncronize free surface elevation time series
    series = fbzi  # Example data for 5 series of length 3601
    Lags = np.zeros(series.shape[1])
    for i in range(1, series.shape[1]):
        A = series[:, i-1]
        B = series[:, i]
        thresholdA = np.percentile(A, 50, axis=None)
        thresholdB = np.percentile(B, 50, axis=None)
        peaks_A, _ = find_peaks(A, height=thresholdA)
        peaks_B, _ = find_peaks(B, height=thresholdB)
        A_weighted = np.zeros_like(A)
        B_weighted = np.zeros_like(B)
        A_weighted[peaks_A] = A[peaks_A]  # Assign the original value to detected peaks
        B_weighted[peaks_B] = B[peaks_B]
        cross_corr_weighted = correlate(A_weighted, B_weighted, mode='full')
        lags = np.arange(-len(A) + 1, len(A))
        lag = lags[np.argmax(cross_corr_weighted)]
        Lags[i] =  lag
        Lags[(Lags>100) | (Lags<-100)] = 0
    
    
    timexb = np.arange(1, 3602)  
    time_sf = np.arange(1, 3002) + 350  
    
    series_new = np.zeros((len(time_sf), series.shape[1]))
    lag = 0
    for i in range(series.shape[1]):
        if i == 0:
            lag = Lags[0]
        else:
            lag += Lags[i]    # Shift the time base by the lag
        time2 = timexb + lag    # Perform interpolation using interp1d
        interpolator = interp1d(time2, series[:, i], bounds_error=False, fill_value=0)
        series_new[:, i] = interpolator(time_sf)
    
    interpolator = interp1d(time2, series[:, i], bounds_error=False, fill_value=0)
    series_new[:, i] = interpolator(time_sf)
    
    # Repeating the array, same water level at every timestep
    fbzs = np.repeat(fbzs[None,:], len(time_sf), axis=0)


    # Adding the time
    fbzs = np.hstack((time_sf.reshape(-1, 1), fbzs))
    fbzi = np.hstack((time_sf.reshape(-1, 1), series_new))

    return  fbzs, fbzi


###############################################################################

def run_SFINCS(now):

#sfincs_exe = '../models/SFINCS/runsfincsmpi.bat'
    
    out_name ='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  +'/'
    
    
    
    templ_path  ='../extras/rarotonga/sfincs'
    run_path = out_name + 'SFINCS'
    run_pathc = os.path.abspath(run_path)
    Xb_path = out_name + 'XBeach'
    
    
    print('Running SFINCS in: ' + run_path)
    
    try:
        shutil.copytree(templ_path, run_path)
    except OSError as error:
         print(error)    
 
    
    # load intersection points XBeach prifiles and SFINCS open boundary
    intersect_file = '../extras/rarotonga/gis/XB_SFINCS_intersection.csv'
    df = pd.read_csv(intersect_file)
    
    fbnd = np.array([df['Lon'],df['Lat']]).T
    index  = df['Index'][:]
    
    fbzs, fbzi = generate_sfincs_forcing(Xb_path,index)
        
    np.savetxt(os.path.join(run_path, 'sfincs.bnd'), fbnd, fmt='%4f')     # Boundary points (x,y) in UTM
    np.savetxt(os.path.join(run_path, 'sfincs.bzs'), fbzs, fmt='%4f')     # (slowly varying) Water level at boundary
    np.savetxt(os.path.join(run_path, 'sfincs.bzi'), fbzi, fmt='%4f')     # (quickly varying) Wave heights at boundary
    
           
    command = f'cd {run_pathc} && runsfincsmpi.bat'
    result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)





