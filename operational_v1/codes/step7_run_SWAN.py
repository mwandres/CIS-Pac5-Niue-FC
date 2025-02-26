# -*- coding: utf-8 -*-
"""
Created on Thu May 20 17:51:04 2021

@author: antonioh
"""
import datetime as dt
import os, shutil
from adcircpy import AdcircMesh
# from AdcircPy import read_mesh
import numpy as np
from subprocess import Popen, PIPE
import subprocess

def execute(cmd):
    popen = Popen(cmd, stdout=PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    return()

def startenddate(start_str,end_str):
    fmt = '%Y%m%d%H'
    dt1 = dt.datetime.strptime(start_str, fmt)
    dt2 = dt.datetime.strptime(end_str, fmt)
    dt1f = dt1.strftime('%Y%m%d.%H%M%S')
    dt2f = dt2.strftime('%Y%m%d.%H%M%S')
    return(dt1f,dt2f)

def parallelization_fort26(folder_name,t1,t2):
    
    
    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(15,5)) 
    # fig,ax1 = plt.subplots(num=1)
    # ax1.plot(    pbnd_x,    pbnd_y,'* k')
    # ax1.plot(    subbnd_x,    subbnd_y,'* r')
    # plt.show()
    
    
    # directory = os.getcwd()
    # fpath = os.path.join(directory,folder_name)
    flist = os.listdir(folder_name)
    
    pmesh = AdcircMesh.open(folder_name+"/fort.14",crs=4326)
    #pmesh = read_mesh(folder_name+"/fort.14")
    px = pmesh.x
    py = pmesh.y
    pbnd = pmesh.ocean_boundaries
    pbnd_val = pbnd.indexes[0]
    pbnd_val = np.array(pbnd_val)+1
    pbnd_x = px[pbnd_val]
    pbnd_y = py[pbnd_val]

    for fl in flist:
        if fl.startswith("PE"):
            subpath = os.path.join(folder_name,fl)
            # sublist = os.listdir(subpath)
            # for fsub in sublist:
            submesh = AdcircMesh.open(subpath+"/fort.14",crs=4326)
            subbnd = submesh.ocean_boundaries
            
            if  len(subbnd.indexes)>0:
                subx = submesh.x
                suby = submesh.y
                subbnd_val = subbnd.indexes[0]
                subbnd_val = np.array(subbnd_val)+1
                subbnd_val2 = list(set(subbnd_val))        
                subbnd_x = subx[subbnd_val2]
                subbnd_y = suby[subbnd_val2]
                
                pid = np.empty(0, int)
                for sbx,sby in zip(subbnd_x,subbnd_y):
                    pos = np.sqrt((pbnd_x-sbx)**2+(pbnd_y-sby)**2).argmin()
                    # pos = np.sqrt((pbnd_x-subbnd_x[0])**2+(pbnd_y-subbnd_y[0])**2).argmin()
                    pid_val = pbnd_val[pos]
                    pid = np.append(pid,pid_val)
                    
                pidf = pid
                subidf = np.array(subbnd_val)
                bnd_txt = []
                for pf14id,subid in zip(pidf,subidf):
                    s = ("BOUN SEGMENT IJ %s VAR FILE 0 '%s.sp2' 1, &\n" % (subid,pf14id))
                    bnd_txt.append(s)
                bnd_txt[-1] = bnd_txt[-1][:-4].strip()
                bnd_str = ''.join(bnd_txt)
        

                # Read in the file
                with open(subpath+"/fort.26", 'r') as f:
                    filedata = f.read()

                # Replace the target string
                filedata = filedata.replace('$%%BOUND_COMMAND%%', bnd_str)
                filedata = filedata.replace('%%BOUND_STARTDATE%%', t1)
                filedata = filedata.replace('%%BOUND_ENDDATE%%', t2)

                # Write the file out again
                with open(subpath+"/fort.26", 'w') as fout:
                    fout.write(filedata)
                    
                   
            else:
                # Read in the file
                with open(subpath+"/fort.26", 'r') as f:
                    filedata = f.read()

                # Replace the target string
                filedata = filedata.replace('%%BOUND_STARTDATE%%', t1)
                filedata = filedata.replace('%%BOUND_ENDDATE%%', t2)

                # Write the file out again
                with open(subpath+"/fort.26", 'w') as fout:
                    fout.write(filedata)
    # return(print("fort26 edited for parallelization"))

def main_fort26(folder_name,t1,t2):
    # directory = os.getcwd()
    # fpath = os.path.join(directory,folder_name)
    # flist = os.listdir(folder_name)
    
    pmesh = AdcircMesh.open(folder_name+"/fort.14",crs=4326)
    #pmesh = read_mesh(folder_name+"/fort.14")
    # px = pmesh.x
    # py = pmesh.y
    pbnd = pmesh.ocean_boundaries
    pbnd_val = pbnd.indexes[0]
    pbnd_val = np.array(pbnd_val)+1
    # pbnd_x = px[pbnd_val]
    # pbnd_y = py[pbnd_val]
    
    
    # bnd_txt = []
    # for pf14id in pbnd_val:
    #     s = ("BOUN SEGMENT IJ %s VAR FILE 0 'Pto_%s.sp2' 1, &\n" % (pf14id ,pf14id))
    #     bnd_txt.append(s)
        
    # bnd_txt[-1] = bnd_txt[-1][:-4].strip()
    # bnd_str = ''.join(bnd_txt)
    
    
    
    # Read in the file
    with open(folder_name+"/swanconf.swn", 'r') as f:
        filedata = f.read()

    # Replace the target string
    #filedata = filedata.replace('$%%BOUND_COMMAND%%', bnd_str)
    filedata = filedata.replace('%%BOUND_STARTDATE%%', t1)
    filedata = filedata.replace('%%BOUND_ENDDATE%%', t2)

    # Write the file out again
    with open(folder_name+"/swanconf.swn", 'w') as fout:
        fout.write(filedata)


############################################################################################################################################
def run_SWAN(now):
  
    out_name='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H")  +'/SWAN'
    # cygfolder = '/cygdrive/d/Projects_SPC/Kiribati/CREWS-KI-Forecast-System-main/operational_v1/runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  +'/'
    ini = now - dt.timedelta(days=2)
    start_str = ini.strftime("%Y") + ini.strftime("%m") + ini.strftime("%d") + ini.strftime("%H") 
    end = ini + dt.timedelta(days=9.5)
    end_str = end.strftime("%Y") + end.strftime("%m") + end.strftime("%d") + end.strftime("%H") 
    
    
    t1,t2 = startenddate(start_str,end_str)
    
    source_dir ='../extras/niue/swan'

    
        
    # List of specific files to copy (with full names, not paths)
    files_to_copy = ['fort.14', 'swaninit', 'swanconf.swn']
    

    # Loop through the list of files
    for filename in files_to_copy:
        full_file_name = os.path.join(source_dir, filename)
        if os.path.isfile(full_file_name):
            shutil.copy(full_file_name, out_name)


    main_fort26(out_name,t1,t2)
   
    out_name_w= os.path.join('D:\\CISPac-5\\CIS-PAC5-Niue\\operational_v1','runs', now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H"),'SWAN')
    command = f"cd {out_name_w} && swanrun swanconf"

    dll_path = r"C:\Program Files (x86)\Common Files\Intel\Shared Libraries\intel64"
    
    # Copy current environment
    env = os.environ.copy()
    env["PATH"] = dll_path + os.pathsep + env["PATH"]  # Prepend correct path
    env["OMP_NUM_THREADS"] = "24"
    # Run your command
    result = subprocess.run(command,shell=True,text=True,capture_output=True,env=env)
   
    # Print the output
    print(result.stdout)
##############################################################################
 
