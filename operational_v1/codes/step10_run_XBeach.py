# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 06:56:11 2024

@author: antonioh
"""
import numpy as np
import shutil
import os
import pandas as pd
import subprocess
from time import sleep


def wirite_Xb_wave_boundary(wparams,fpath):
    
    lines = [
    "Hm0        =  %(hm0)s\n" % wparams,
    "fp         =  %(fp)s\n" % wparams,
    "mainang    =  %(mainang)s\n" % wparams,
    "gammajsp   =  %(gammajsp)s\n" % wparams,
    "fnyq       =  %(fnyq)s\n" % wparams]
    

    with open(fpath + '/jonswap_1.txt', 'w') as file:
        file.writelines(lines)
        
    with open(fpath + '/jonswap_2.txt', 'w') as file:
        file.writelines(lines)
        
    return()


def write_Xb_water_level(level,fpath):
    
  inpfile = fpath + '/params.txt'
  
  with open(inpfile, 'r') as file:
      lines = file.readlines()
      
  with open(inpfile, 'w') as file:
      for line in lines:
          if line.startswith('%z0%'):
              new_line = 'zs0           = ' + level + '\n'
              file.write(new_line)
          else:
              file.write(line)  
  return()





###############################################################################
def run_XBeach_profiles(now):
    
    out_name ='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  +'/'
    
    xbeach_exe = '../models/XBeach/xbeach.exe'
    xbeach_exec = os.path.abspath(xbeach_exe)
    
    ## read XBeach forcings
    df = pd.read_csv(out_name + 'results/Points_results.csv')

    templ_path  ='../extras/rarotonga/xbeach'
    run_path = out_name + 'XBeach'
    
    print('Running XBeach in: ' + run_path)

    # Generate XBeach forcings
    try:
        shutil.copytree(templ_path, run_path)
    except OSError as error:
         print(error)    
    
    xb_run_path = []
    for ii in range(len(df['Hs_max[m]'])):
        name = 'Prof_'+ str(ii+1)
        fpath = os.path.join(run_path, name)
        xb_run_path.append(fpath)
        xb_run_path.append(fpath)
        

        wparams = {"hm0": str(round(df['Hs_max[m]'][ii],3)),
        "fp": str(round(1/df['Tp_max[s]'][ii],3)),
        "mainang": str(round(df['Dp_max[deg north]'][ii],3)),
        "gammajsp": str(round(1/df['Fsp_max[]'][ii],3)),
        "s": str(round((2/(np.pi/180*df['Dsp_max[deg]'][ii])**2)-1,3)),
        "fnyq": str(round(1,3))}

        wirite_Xb_wave_boundary(wparams,fpath)
        write_Xb_water_level(str(round(df['Zeta_max[m]'][ii],3)),fpath)
        
        
        
        # xb_run_path = [f'..\\runs\\024101400\\XBeach\\Prof_{i}' for i in range(1, 9)]
        # xbeach_exec = '..//models//XBeach//xbeach.exe'
        
    
    # Running XBeach
    processes = []
    for f in xb_run_path:

        #p = subprocess.Popen(xbeach_exec, cwd=f, shell=True)
        p = subprocess.Popen(xbeach_exec, cwd=f, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
        processes.append(p)
        
    while any([p.poll() is None for p in processes]):
        sleep(0.1)

