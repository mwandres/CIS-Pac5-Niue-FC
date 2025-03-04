# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 08:51:33 2021

@author: judithg
"""
# import math
from adcircpy import AdcircMesh
import numpy as np
import datetime as dt
from matplotlib import pyplot as plt
import matplotlib.tri as tri
import scipy.io as sio
import os
# os.environ['PROJ_LIB'] = 'C:/Users/antonioh/Anaconda3/Library/share/proj/'
from matplotlib.colors import ListedColormap#, LinearSegmentedColormap
import netCDF4 as nc

import geopandas as gpd





def load_SWAN_data(fl_name,var_name):
     if var_name == 'time':
        time= True
        avname = 'Hs'
     else :
        time = False
        avname = var_name
        
     infos = sio.whosmat(fl_name)
     names, sizes, types  = zip(*infos)
     vnames = []
     for n in names:
        if avname in n[0:len(avname):1]:
            vnames.append(n)
     vnames.sort()
     if time :
        data = []
        for i in range(len(vnames)): 
            aux = vnames[i][len(vnames[0])-15:]
            data.append(aux)
           
     else:    
        MAT = sio.loadmat(fl_name)         
        data = np.zeros([len(vnames),len(MAT[vnames[0]].T)])
        for i in range(len(vnames)):            
            data[i,:] = MAT[vnames[i]]   
            
     return(data)
            
def inter2points2nc(folder_tmp,file_name,time_waves,points_all,point_type,triang,Hs,Tm,Tp,Fspr,Dir,Dspr,Windv_x,Windv_y):
    
    
    Hsp = np.zeros((len(points_all),len(time_waves)))
    Tmp = np.zeros((len(points_all),len(time_waves)))
    Tpp = np.zeros((len(points_all),len(time_waves)))
    Fsprp = np.zeros((len(points_all),len(time_waves)))
    Dirp = np.zeros((len(points_all),len(time_waves)))
    Dsprp = np.zeros((len(points_all),len(time_waves)))
    Windv_xp= np.zeros((len(points_all),len(time_waves)))
    Windv_yp = np.zeros((len(points_all),len(time_waves)))
    
    for ts in range(0,len(Hs)):
        
        #interpolate to the points surroounding the land but still too deep for wave breaking
        intsurf = tri.LinearTriInterpolator(triang,Hs[ts])
        Hsp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Tp[ts])
        Tpp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Tm[ts])
        Tmp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Fspr[ts])
        Fsprp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Dir[ts])
        Dirp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Dspr[ts])
        Dsprp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Windv_x[ts])
        Windv_xp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Windv_y[ts])
        Windv_yp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
    
    fn = folder_tmp  + file_name
    try:
        os.remove(fn)
        print('removing wind_and_waves.nc of previous run')
    except:
        print('saving time serie at buoy location')
        
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    time = ds.createDimension('time', None)
    index = ds.createDimension('index', None)
    times = ds.createVariable('time', 'f8', ('time',))
    times.units='hours since 1950-01-01 00:00:00'
    times.calendar='gregorian'
    times[:] = [nc.date2num(x,units=times.units,calendar=times.calendar) for x in time_waves]
    lonnc = ds.createVariable('lon', 'f4', ('index',))
    lonnc.units ='degrees_east'
    lonnc[:] = points_all[:,0]
    latnc = ds.createVariable('lat', 'f4', ('index',))
    latnc.units ='degrees_north'
    latnc[:] = points_all[:,1]
    Hsnc= ds.createVariable('Hs', 'f4', ('index','time'))
    Hsnc.units = 'm'
    Hsnc[:,:] =  Hsp
    Tpnc= ds.createVariable('Tp', 'f4', ('index','time'))
    Tpnc.units = 's'
    Tpnc[:,:] =  Tpp
    Tmnc= ds.createVariable('Tm', 'f4', ('index','time'))
    Tmnc.units = 's'
    Tmnc[:,:] =  Tmp
    Fsprnc= ds.createVariable('Fspr', 'f4', ('index','time'))
    Fsprnc.units = 'frequency spreading swan'
    Fsprnc[:,:] =  Fsprp
    Dirnc= ds.createVariable('Dir', 'f4', ('index','time'))
    Dirnc.units = 'degress from north (north=0, east=90)'
    Dirnc[:,:] =  Dirp
    Dsprnc= ds.createVariable('Dspr', 'f4', ('index','time'))
    Dsprnc.units = 'one-dided directional width swan degress'
    Dsprnc[:,:] =  Dsprp
    Windxnc= ds.createVariable('Windx', 'f4', ('index','time'))
    Windxnc.units = 'm/s'
    Windxnc[:,:] =  Windv_xp
    Windync= ds.createVariable('Windy', 'f4', ('index','time'))
    Windync.units = 'm/s'
    Windync[:,:] =  Windv_yp
    ds.setncattr_string('point_type', point_type)
    ds.close()

def plot_Hs_Dir_maps(nat,tar,dtt,x,y,triang,Time,Hs,Dir,result_folder):


    f = open('../extras/rarotonga/swan/Hs_colormap.dat', 'r') # 'r' = read
    colors = np.genfromtxt(f, delimiter='  ')
    f.close()
    colors= np.hstack((colors,np.ones((len(colors[:,1]),1))))
    newcmp = ListedColormap(colors)

    ## Mesh grid values for interpolation
    # mesh grid national scale

    # linear interpolation - Dir
    Dir_cor = (270-Dir)*np.pi/180
    xDir = np.cos(Dir_cor)
    yDir = np.sin(Dir_cor)
    
    
    levels = np.arange(0.,6, 0.1)
    
    
    # Gilbert grid
    # xllg = 172.06
    # yllg = -2.15
    # xurg = 176.6
    # yurg = 3.55   
    # lgrid_x, lgrid_y = np.mgrid[xllg:xurg:20j,yllg:yurg:20j]
    
    
    # Tarawa grid

 

    # xllt = -160.2
    # yllt = -21.6
    # xurt = -159.3
    # yurt = -20.8


    xllt = -159.85
    yllt = -21.29
    xurt = -159.70
    yurt = -21.18
    # xllt = -160
    # yllt = -21
    # xurt = -159.6
    # yurt = -21.4
    
    
    sgrid_x, sgrid_y = np.mgrid[xllt:xurt:30j,yllt:yurt:30j]


    
    plt.ioff()
    for ts in range(24,len(Hs),dtt):  
        
        fxDir = tri.LinearTriInterpolator(triang,xDir[ts])
        fyDir = tri.LinearTriInterpolator(triang,yDir[ts])
        
    #     if nat==1:
    #         lxdir = fxDir(lgrid_x,lgrid_y)
    #         lydir = fyDir(lgrid_x,lgrid_y)

    #         fig, ax = plt.subplots(figsize=(10,10))
    #         ax.set_aspect('equal')
    #         tcf = ax.tricontourf(triang, Hs[ts],levels=levels,cmap=newcmp,zorder=1)
    #         cbar = fig.colorbar(tcf, ax=ax)
    #         cbar.solids.set_edgecolor("face")
    #         ax.quiver(lgrid_x,lgrid_y,lxdir,lydir,scale_units='xy',angles='xy',zorder=2) #Dir x and 
    #         ax.set_title('')
    #         ax.set_xlabel('Longitude (degrees)')
    #         ax.set_ylabel('Latitude (degrees)')
    #         ax.set_xlim(xllg,xurg)
    #         ax.set_ylim(yllg,yurg)
    #         plt.title(Time[ts])
    #         fileprint=result_folder  + Time[ts]
    #         plt.savefig(fileprint)
    #         #plt.show()
    #         plt.close(fig)
    
        if tar==1:
            

            sxdir = fxDir(sgrid_x,sgrid_y)
            sydir = fyDir(sgrid_x,sgrid_y)

            
            fig2, ax2 = plt.subplots(figsize=(10,10))
            ax2.set_aspect('equal')
            tcf2 = ax2.tricontourf(triang, Hs[ts],levels=levels,cmap=newcmp,zorder=1)
            cbar2 = fig2.colorbar(tcf2, ax=ax2)
            cbar2.solids.set_edgecolor("face")
            ax2.quiver(sgrid_x,sgrid_y,sxdir,sydir,scale_units='xy',angles='xy',zorder=2) #Dir x and y
            ax2.set_xlabel('Longitude (degrees)')
            ax2.set_ylabel('Latitude (degrees)')
            ax2.set_xlim(xllt,xurt)
            ax2.set_ylim(yllt,yurt)
            plt.title(Time[ts])
            fileprint=result_folder + 'hr_' + Time[ts]
            plt.savefig(fileprint)
            plt.close(fig2)


def SWAN2nc(x,y,triang,time_waves,Hs,Tm,Tp,Dir,Windv_x,Windv_y,xll,yll,xur,yur,inc,namenc):

    
    x = np.arange(xll,xur,inc) 
    y = np.arange(yll,yur,inc)
    lgrid_x, lgrid_y = np.meshgrid(x,y)
    
    
    
    

    Hsp = np.ones((len(time_waves)-23,lgrid_x.shape[0],lgrid_x.shape[1]), dtype=np.float16)*-999
    Tmp = np.ones((len(time_waves)-23,lgrid_x.shape[0],lgrid_x.shape[1]), dtype=np.float16)*-999
    Tpp = np.ones((len(time_waves)-23,lgrid_x.shape[0],lgrid_x.shape[1]), dtype=np.float16)*-999
    Dirp = np.ones((len(time_waves)-23,lgrid_x.shape[0],lgrid_x.shape[1]), dtype=np.float16)*-999
    # Dirp_u = np.ones((len(time_waves)-23,lgrid_x.shape[0],lgrid_x.shape[1]), dtype=np.float16)*-999 
    # Dirp_v = np.ones((len(time_waves)-23,lgrid_x.shape[0],lgrid_x.shape[1]), dtype=np.float16)*-999 
    Windp = np.ones((len(time_waves)-23,lgrid_x.shape[0],lgrid_x.shape[1]), dtype=np.float16)*-999
    WindDirp = np.ones((len(time_waves)-23,lgrid_x.shape[0],lgrid_x.shape[1]), dtype=np.float16)*-999

    for ts in range(23,len(time_waves)):# change the starting hour into range to the lenght of the spin-up 
        
        intsurf = tri.LinearTriInterpolator(triang,Hs[ts])
        Hsp[ts-23,:,:] =intsurf(lgrid_x,lgrid_y)
        intsurf = tri.LinearTriInterpolator(triang,Tp[ts])
        Tpp[ts-23,:,:] =intsurf(lgrid_x,lgrid_y)
        intsurf = tri.LinearTriInterpolator(triang,Tm[ts])
        Tmp[ts-23,:,:] =intsurf(lgrid_x,lgrid_y)
        intsurf = tri.LinearTriInterpolator(triang,Dir[ts])
        Dirp[ts-23,:,:] = intsurf(lgrid_x,lgrid_y)
        # dir_ = intsurf(lgrid_x,lgrid_y)
        # nan_pos = np.isnan(dir_)
        # Dir_cor = (270-dir_)*np.pi/180
        # xDir = np.cos(Dir_cor)
        # xDir[nan_pos] = np.nan
        # yDir = np.sin(Dir_cor)
        # yDir[nan_pos] = np.nan
        
        # Dirp_u[ts-23,:,:] = xDir
        # Dirp_v[ts-23,:,:] = yDir
     
        # plt.figure()
        # plt.pcolor(dir_, cmap='turbo')
        # plt.show()
        # plt.colorbar()
        # Dirp[ts-23,:,:] =  dir_
        
        intsurf = tri.LinearTriInterpolator(triang,Windv_x[ts])
        wx=intsurf(lgrid_x,lgrid_y)
        intsurf = tri.LinearTriInterpolator(triang,Windv_y[ts])
        wy =intsurf(lgrid_x,lgrid_y)
        Windp[ts-23,:,:]= np.sqrt(wx**2 + wy**2)
        WindDirp[ts-23,:,:] = 270 - np.arctan2(wy,wx)*180/np.pi
  
    
    ds = nc.Dataset(namenc, 'w', format='NETCDF4')
    time = ds.createDimension('time', None)
    lon = ds.createDimension('lon', None)
    lat = ds.createDimension('lat', None)
    times = ds.createVariable('time',  'f8', ('time',))
    times.units='hours since 1950-01-01 00:00:00'
    times.calendar='gregorian'
    times[:] = [nc.date2num(ti,units=times.units,calendar=times.calendar) for ti in time_waves]
    times = times[23:]
    lonnc =ds.createVariable('lon',  'f4', ('lon',))
    lonnc.units ='degrees_east'
    lonnc[:] = x
    latnc = ds.createVariable('lat','f4', ('lat',))
    latnc.units ='degrees_north'
    latnc[:] = y
    Hsnc= ds.createVariable('Hs', 'f4', ('time','lat','lon'))
    Hsnc.units = 'm'
    Hsnc[:] =  Hsp
    Tpnc= ds.createVariable('Tp', 'f4', ('time','lat','lon'))
    Tpnc.units = 's'
    Tpnc[:] =  Tpp
    Tmnc= ds.createVariable('Tm', 'f4', ('time','lat','lon'))
    Tmnc.units = 's'
    Tmnc[:] =  Tmp
    Dirnc= ds.createVariable('sea_surface_wave_from_direction', 'f4', ('time','lat','lon'))
    Dirnc.units = 'degress from north (north=0, east=90)'
    Dirnc[:] =  Dirp
    # Dirncu= ds.createVariable('Dir_u', 'f4', ('time','lat','lon'))
    # Dirncu.units = 'u component'#degress from north (north=0, east=90)'
    # Dirncu[:] =  Dirp_u
    # Dirncv= ds.createVariable('Dir_v', 'f4', ('time','lat','lon'))
    # Dirncv.units = 'v component'#degress from north (north=0, east=90)'
    # Dirncv[:] =  Dirp_v
    Windc= ds.createVariable('Wind','f4', ('time','lat','lon'))## fix wind direction compliant with nvWMS
    Windc.units = 'm/s'
    Windc[:] =  Windp
    DirWindc= ds.createVariable('DirWind', 'f4', ('time','lat','lon'))
    DirWindc.units = 'degress from north (north=0, east=90)'
    DirWindc[:] =  WindDirp
    ds.close()
    
    
   

def readvars_nc(nc_fname,varname):
    ncc = nc.Dataset(nc_fname,'r') 

    if varname == 'time':
        time_var = ncc.variables[varname] 
        time_or = nc.num2date(ncc[varname],time_var.units, time_var.calendar)
        time_str=[time_or[i].strftime('%Y-%m-%d %H:%M') for i in range(0,len(time_or))]
        var = [dt.datetime.strptime(time_str[i],'%Y-%m-%d %H:%M') for i in range(0,len(time_str))]
    else:
        var = np.array(ncc[varname])
        
    ncc.close()
    return(var)

def buoyhindcastforecast(now):
    folder_name='../archives/'
    flist = os.listdir(folder_name)
    dates = []
    for file in flist:
        if file.endswith("_buoy.nc"):
            nb = file.split("_")[0]
            dates.append(dt.datetime.strptime(nb+'0000',"%Y%m%d%H%M%S"))       
    dates.sort()
    

    output =open('../tmp/buoy_hindcast_forecast.txt','w')
    for r in dates:
         
        nc_fname = folder_name +  r.strftime("%Y") + r.strftime("%m") + r.strftime("%d") + r.strftime("%H") + '_buoy.nc'
        time = readvars_nc(nc_fname,'time')
        year = np.array([t.strftime('%Y') for t in time])
        month = np.array([t.strftime('%m') for t in time])
        day = np.array([t.strftime('%d') for t in time])
        hour = np.array([t.strftime('%H') for t in time])
        Hs = readvars_nc(nc_fname,'Hs')
        Tm = readvars_nc(nc_fname,'Tm')
        Tp = readvars_nc(nc_fname,'Tp')
        Dir = readvars_nc(nc_fname,'Dir')
        Wx = readvars_nc(nc_fname,'Windx')
        Wy = readvars_nc(nc_fname,'Windy')
        Wind = (Wx**2+Wy**2)**(1/2)
        WindDir= 270 - (np.arctan2(Wy,Wx)*180/np.pi)
        
        if r != dates[-1]:
            for i in range(47,53):
                s="{}".format(year[i])+" {}".format(month[i])+" {}".format(day[i])+" {}".format(hour[i])+ \
                " {:.3f}".format(Hs[0,i]) +" {:.3f}".format(Tm[0,i]) + " {:.3f}".format(Tp[0,i])+" {:.3f}".format(Dir[0,i])  +\
                " {:.3f}".format(Wind[0,i]) +" {:.3f}".format(WindDir[0,i]) + '\n'
                output.write(s)
                          
        elif r==dates[-1]:
            
            for i in range(47,len(time)):
                s="{}".format(year[i])+" {}".format(month[i])+" {}".format(day[i])+" {}".format(hour[i])+ \
                " {:.3f}".format(Hs[0,i]) +" {:.3f}".format(Tm[0,i]) + " {:.3f}".format(Tp[0,i])+" {:.3f}".format(Dir[0,i])  +\
                " {:.3f}".format(Wind[0,i]) +" {:.3f}".format(WindDir[0,i]) + '\n'
                output.write(s)
    output.close()
        

###################################################################################################################
def postprocess_SWAN(now,plot):

    folder_tmp ='../tmp/'
    out_name='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  +'/SWAN/'
    result_folder = '../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  + '/results/'
    try :
        os.mkdir(result_folder)
    except:
        print("results already exists")
        
        
    fl_name = out_name +'output.mat'

      
    shapefile_path = '../extras/rarotonga/gis/500m_transects.shp'
    gdf = gpd.read_file(shapefile_path)
    gdf = gdf.to_crs(epsg=4326)
      
    points_all = []
    for i in range(0, len(gdf['geometry'])):
        ptos = np.array(gdf['geometry'][i].coords[:])
        points_all.append( ptos[0,:])
    
    
    points_all = np.array(points_all)
    point_type = ['forereef' for i in range(len(points_all))]

    
    f14 = out_name + 'fort.14'
    # open mesh file
    mesh = AdcircMesh.open(f14,crs=4326)
    
    x=np.array(mesh.x)
    y=np.array(mesh.y)
    
    # fig, ax = plt.subplots(figsize=(10,10))
    # ax.plot(x,y,'.')
    
    Time = load_SWAN_data(fl_name,"time")
    time_waves = [dt.datetime.strptime(Time[i],'%Y%m%d_%H%M%S') for i in range(0,len(Time))]
    
    
    Hs = load_SWAN_data(fl_name,"Hsig")#*2###################################### Watch out!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Tp = load_SWAN_data(fl_name,"RTpeak")#+3####################################
    Tm = load_SWAN_data(fl_name,"Tm02")
    Fspr = load_SWAN_data(fl_name,"FSpr_")
    Dir = load_SWAN_data(fl_name,"Dir_")
    Dspr = load_SWAN_data(fl_name,"Dspr_")
    
    Windv_x = load_SWAN_data(fl_name,"Windv_x")
    Windv_y = load_SWAN_data(fl_name,"Windv_y")
    
    
    Hs = np.nan_to_num(Hs)
    Tp = np.nan_to_num(Tp)
    Tm = np.nan_to_num(Tm)
    Fspr = np.nan_to_num(Fspr)
    Dir = np.nan_to_num(Dir)
    Dspr = np.nan_to_num(Dspr)
    Windv_x = np.nan_to_num(Windv_x)
    Windv_y = np.nan_to_num(Windv_y)
   

    
     ## Trimesh triangulation used for plots and interpolation
    triangulation = mesh.triangles[:]-1
    triang = tri.Triangulation(x,y,triangulation)
    
 

     
    ## Mesh grid values for interpolation island scale 
    
    

    # xll = -160.2
    # yll = -21.6
    # xur = -159.3
    # yur = -20.8

    xll = -159.85
    yll = -21.29
    xur = -159.70
    yur = -21.18
    inc = 0.00025# about 10 m 
    
    
    
    namenc = folder_tmp  + '/Rarotonga.nc'  
    SWAN2nc(x,y,triang,time_waves,Hs,Tm,Tp,Dir,Windv_x,Windv_y,xll,yll,xur,yur,inc,namenc)
   
    file_name='wind_and_waves.nc'
    inter2points2nc(folder_tmp,file_name,time_waves,points_all,point_type,triang,Hs,Tm,Tp,Fspr,Dir,Dspr,Windv_x,Windv_y)
    
    
    # buoy =np.array([[172.911383,   1.42725 ]])
    # archfolder='../archives/' 
    # file_name = now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H") + '_buoy.nc'
    # inter2points2nc(archfolder,file_name,time_waves,buoy,'buoy',triang,Hs,Tm,Tp,Dir,Windv_x,Windv_y)
    
    
    # buoyhindcastforecast(now)

    if plot==1:
        
        nat=0 #plot national scale
        tar=1 #plot Tarawa scale
        dtt = 1 #time interval in hours
      
        plot_Hs_Dir_maps(nat,tar,dtt,x,y,triang,Time,Hs,Dir,result_folder)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    