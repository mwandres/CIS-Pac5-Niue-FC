# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 08:58:04 2021

@author: antonioh
"""

import pandas as pd
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import datetime as dt
import netCDF4 as nc
import os
from scipy.interpolate import interp1d
import geopandas as gpd
import contextily as ctx
from shapely.geometry import Point
import geojson


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


def rbfphi_gaussian(r, const):
    return np.exp(-0.5*r*r/(const*const))

def RBF_Interpolation(rbf_constant, rbf_coeff, nodes, x):

    phi = rbfphi_gaussian   # gaussian RBFs
    rbf_coeff = rbf_coeff.flatten()

    dim, n = nodes.shape
    dim_p, n_p = x.shape

    f = np.zeros(n_p)
    r = np.zeros(n)

    for i in range(n_p):
        r = np.linalg.norm(
            np.repeat([x[:,i]], n, axis=0)-nodes.T,
            axis=1
        )
        s = rbf_coeff[n] + np.sum(rbf_coeff[:n] * phi(r, rbf_constant))

        # linear part
        for k in range(dim):
            s = s + rbf_coeff[k+n+1] * x[k,i]

        f[i] = s

    return(f)

def plot_risk_maps(coordinates,time_w,Risk_Cat,thresholds,fig_folder):
    
  
    corners =np.array([[-159.859, -21.308],[-159.691, -21.169]])
    gdf = gpd.GeoDataFrame(geometry=[Point(xy) for xy in corners], crs="EPSG:4326")
    
    
    riskcolors=np.array([[0.0745, 0.6235, 1.0000], [0.8706, 0.8706, 0.3725],[1, 0.4118, 0.1608]])
    newcmp = ListedColormap(riskcolors)   
    plt.ioff()
    # for ts in range(24,len(time_w)):
       
    risk = np.amax(Risk_Cat,0)

    
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_aspect('equal')
    tcf=ax.scatter(coordinates[:,0],coordinates[:,1],40,risk,cmap=newcmp,zorder=1,vmin=-0.5, vmax=2.5)
    ctx.add_basemap(ax, crs="EPSG:4326", source=ctx.providers.OpenStreetMap.Mapnik)
    ax.set_xlim(-159.842, -159.713)
    ax.set_ylim(-21.284, -21.185)
    xlbl = ax.axes.get_xticks()
    for h in range(0,len(xlbl)):
        ax.axvline(xlbl[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
    ylbl = ax.axes.get_yticks()
    for h in range(0,len(ylbl)):
        ax.axhline(ylbl[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
    cbar = fig.colorbar(tcf, ax=ax)
    cbar.solids.set_edgecolor("face")
    cbar. set_ticks([0,1,2])
    cbar. set_ticklabels(["No flood", "Minor flood", "Moderate flood"])
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_title(' Max flood risk from '+time_w[48].strftime('%Y%m%d%H%M')+' to '+time_w[-1].strftime('%Y%m%d%H%M'))
    # fileprint=fig_folder + 'Tarawa_fld_' + time_w[ts].strftime('%Y%m%d%H%M') 
    fileprint=fig_folder + 'Rarotonga_Max_fld_risk'
    plt.savefig(fileprint)
    # plt.show()
    plt.close(fig)

    return()

def risk_time_series(now,TWL,time_w,time_tide_min,sla,tide_min,thresholds,fig_folder):
    
    plt.ioff()
    for npo in range(len(TWL[1,:])):
        
        xor = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_w]
        xnew = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_tide_min]
        # xor = np.linspace(0, len(time_w)-1, len(time_w), endpoint=True).T
        # xnew = np.linspace(0, len(time_w)-1, len(time_w)*60, endpoint=True).T
        f = interp1d(xor,TWL[:,npo],kind='cubic')
        twl_min=f(xnew)
        f = interp1d(xor,sla,kind='cubic')
        sla_min=f(xnew)
        # plot time series
        
        #-- create plot with tidal displacements, high and low tides and dates
        
 
        #-- differentiate to calculate high and low tides
        diff = np.zeros_like(time_tide_min, dtype=np.float64)
        #-- forward differentiation for starting point
        diff[0] = tide_min[1] - tide_min[0]
        #-- backward differentiation for end point
        diff[-1] = tide_min[-1] - tide_min[-2]
        #-- centered differentiation for all others
        diff[1:-1] = (tide_min[2:] - tide_min[0:-2])/2.0
        htindex, = np.nonzero((np.sign(diff[0:-1]) >= 0) & (np.sign(diff[1:]) < 0))
        # ltindex, = np.nonzero((np.sign(diff[0:-1]) <= 0) & (np.sign(diff[1:]) > 0))
        
        
        
        
        plt.figure(figsize=(25,10)) 
        # fig,ax1 = plt.subplots(num=1)
        plt.fill_between(time_tide_min,twl_min,0,color='paleturquoise',label='Total water level')
        plt.fill_between(time_tide_min,tide_min+sla_min,0,color='c',label='Ocean water level')
        plt.plot(time_tide_min,twl_min,color='paleturquoise',linewidth=0.5)
        plt.plot(time_tide_min,tide_min,'--k',linewidth=0.5,label='Astronomical tide')
        plt.plot(np.matlib.repmat(now,5,1),np.linspace(-1,thresholds[npo,4]+0.6,5),color= 'gray',linewidth=1)
        plt.plot(time_w,np.matlib.repmat(thresholds[npo,2],len(time_w),1),color=[1, 0.4118, 0.1608])
        plt.text(time_tide_min[-4000], thresholds[npo,2]+0.02, 'Moderate flood threshold',color=[1, 0.4118, 0.1608],size=15)
        plt.plot(time_w,np.matlib.repmat(thresholds[npo,0],len(time_w),1),color=[0.8706, 0.8706, 0.3725])
        plt.text(time_tide_min[-3800], thresholds[npo,0]+0.02, 'Minor flood threshold',color=[0.8706, 0.8706, 0.3725],size=15)
        #ax1.plot(time_tide_min[htindex],twl_min[htindex],"v",color= 'gray',markersize=5)
        for h in range(2,len(htindex)-1):
            text=time_tide_min[htindex[h]].strftime("%H:%M")
            plt.plot(time_tide_min[htindex[h]],twl_min[htindex[h]],"v",color= 'gray',markersize=5)
            plt.text(time_tide_min[htindex[h]]-dt.timedelta(hours=3),twl_min[htindex[h]]+0.08,text,color= 'gray',size=15) 
        plt.plot(now,thresholds[npo,4]+0.6,"v",color= 'gray',markersize=25)
        plt.text(now-dt.timedelta(hours=3),thresholds[npo,4]+0.6+0.09, 'Now',color='gray',size=15)
        plt.ylim((-1,thresholds[npo,4]+0.6))       
        plt.xlim((time_tide_min[60*24],time_tide_min[-1]))
        plt.xlabel('UTC Time',fontsize=15)
        plt.ylabel('Water level [m]',fontsize=15)
        plt.yticks(fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        for h in range(0,len(time_w),24):
            plt.axvline(time_w[h].replace(hour=0),color='gray',lw=0.5,ls='dashed',dashes=(11,5))
        ylbl = plt.yticks()[0]
        for h in range(0,len(ylbl)):
            plt.axhline(ylbl[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
        plt.legend(loc='upper center', ncol = 3,fontsize=20)
        # plt.show()
        fileprint=fig_folder + 'TSeries_' + str(npo) 
        plt.savefig(fileprint)
        plt.close()
    return()

  
def TWL_points2nc(result_folder,Risk_Cat,TWL,time_w,time_tide_min,sla,tide_min,thresholds):

    x=thresholds[:,0]
    y=thresholds[:,1]
    thresholdss = thresholds[:,[2,4]]
    time_10min = [time_tide_min[i] for i in range(1440, len(time_tide_min), 10)]# vhange the starting time into range to the lenght of the spin-up  in minutes
        
    xorhour = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_w]
    xormin = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_tide_min]
    xor10min = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_10min]
    
    
    f = interp1d(xormin,tide_min)
    tide_10min=f(xor10min)
    f = interp1d(xorhour,sla)
    sla_10min=f(xor10min)   
    
    riskmax=np.max(Risk_Cat,axis=0)

    TWL_10min = np.zeros((len(time_10min),len(x)))       
    for npo in range(len(x)):          
        f = interp1d(xorhour,TWL[:,npo],kind='cubic')
        TWL_10min[:,npo] =f(xor10min)

    
    
    fn = result_folder  + 'Risk_results.nc'          
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    time = ds.createDimension('time', None)
    index = ds.createDimension('index', len(riskmax))
    thres = ds.createDimension('thres', None)
    index = np.arange(0, len(riskmax))
    times = ds.createVariable('time', 'f8', ('time',))
    times.units='hours since 1950-01-01 00:00:00'
    times.calendar='gregorian'
    times[:] = [nc.date2num(x,units=times.units,calendar=times.calendar) for x in time_10min]
    lonnc = ds.createVariable('lon', 'f4', ('index',))
    lonnc.units ='degrees_east'
    lonnc[:] = x
    latnc = ds.createVariable('lat', 'f4', ('index',))
    latnc.units ='degrees_north'
    latnc[:] = y
    TWLnc= ds.createVariable('TWL', 'f8', ('index','time'))
    TWLnc.units = 'm'
    TWLnc[:,:] =  TWL_10min.T
    SLAnc= ds.createVariable('SLA', 'f8', ('time'))
    SLAnc.units = 'm'
    SLAnc[:] =  sla_10min
    Tidenc= ds.createVariable('Tide', 'f8', ('time'))
    Tidenc.units = 'm'
    Tidenc[:] =  tide_10min
    Riskmaxnc= ds.createVariable('Riskmax', 'f4', ('index'))
    Riskmaxnc.units = 'No risk=0, Minor=1, Moderate=2'
    Riskmaxnc[:] =  riskmax
    Tresholdsnc= ds.createVariable('Thresholds', 'f8', ('index','thres'))
    Tresholdsnc.units = 'm'
    Tresholdsnc[:,:] =  thresholdss
    ds.close()



def estimate_wave_runup(Hs_,Tp_,Dir_, Dn_):
        '''
        Wave runup parametric model as proposed in Merrifield et al., 2014
        It has been orrected based on the empirical rellation with XBeach profiles in Rarotonga
        
        2024/10/22: Fixed an issue related to not solving the full circle when looking
        at waves directions, now removing angles out of the [-90;90] range
        
        Parameters
        ----------
        Dn : np.ndarray (N,)
            Transects orientation with regard to true north [rad].
         
        Returns
        -------
        eta : np.ndarray (M,N)
            Maximum wave setup [m].
            
            
            
        ''' 
        gamma = 1.2
        g = 9.81
        
        
        angle  = Dir_ - Dn_
        
        
        term1 = Hs_**2 * Tp_ / (4*np.pi)
        term2 = abs(np.cos(np.radians(angle)))* np.sqrt(gamma * g)
        term2[(angle < -90) | (angle > 90)] = 0
       
        hb = (term1 * term2) ** (2/5)
        hb[hb<0] = 0
        runup = 0.33 * hb - 0.11
        ## add empirical correction for Rarotonga, considered in the definition of the thresholds
        runup = 2*runup-0.54  
        runup[runup<0] = 0
        
        return runup  
    
    
def riskCat2geojson(df,result_folder):
        
    feature_collection = {"type": "FeatureCollection", "features": []}
 
    for ind in df.index:
        lat = df['lat'][ind]
        lon = df['lon'][ind]
        risk = df['Risk[0=no risk,1=minor,2=moderate]'][ind]
        # name = df['name'][ind]
        # id = df['id'][ind]
        feature = {"type": "Feature", "geometry": {"type": "Point", "coordinates": [lon, lat]},
                   "properties": {"Risk": risk,"link":f"root-path/TSeries_{ind}.png"}}
        feature_collection['features'].append(feature)
       
    with open(result_folder+'Points_results.geojson', 'w') as fp:
        geojson.dump(feature_collection, fp)   
    
    return()

###############################################################################
def make_flood_points(now,fmap,fseries):
    
    out_name='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  +'/'
    result_folder = out_name + 'results/'
    try :
        os.mkdir(result_folder)
    except:
        print("results already exists")
    
    folder_tmp ='../tmp/'
    
    ##- load TWL thresholds from the analysis of the 44 years hindcast data
    f = open('../extras/rarotonga/gis/operational_thresholds_1_2_5_10_25_50_100.dat', 'r') 
    thresholds = np.genfromtxt(f, delimiter='  ')
    f.close()
       
    ##- load transects and select the ocean points to calculate TWL proxy
    shapefile_path = '../extras/rarotonga/gis/500m_transects.shp'
    gdf = gpd.read_file(shapefile_path)
    gdf = gdf.to_crs(epsg=4326)
    
  
  
    points = np.empty((len(gdf['geometry']),3))
   
    for i in range(0, len(gdf['geometry'])):

        ptos = np.array(gdf['geometry'][i].coords[:])
        angle_rad = np.arctan2(ptos[1,1]  - ptos[0,1] , ptos[1,0]  - ptos[0,0])
        angle_deg = np.degrees(angle_rad)  
        angle_deg = 270 - angle_deg  % 360    
        if angle_deg < 0:
            angle_deg += 360
        points[i] = (ptos[0,:][0], ptos[0,:][1], angle_deg)
        
    npoints = len(points)   
        
    ##- load evething else
    nc_fname = folder_tmp  + 'sla_hourly.nc'
    time_sla = readvars_nc(nc_fname,'time')
    sla = readvars_nc(nc_fname,'SLA')/100
    
    nc_fname = folder_tmp  + 'tide_hourly.nc'
    time_tide = readvars_nc(nc_fname,'time')
    tide = readvars_nc(nc_fname,'tide')/100
    
    nc_fname = folder_tmp  + 'tide_minute.nc'
    time_tide_min = readvars_nc(nc_fname,'time')
    tide_min = readvars_nc(nc_fname,'tide')/100
    
    nc_fname = folder_tmp  + 'IB_hourly.nc'
    time_ib = readvars_nc(nc_fname,'time')
    ib = readvars_nc(nc_fname,'IB')/100
    
    nc_fname = folder_tmp  + 'IB_hourly.nc'
    time_ib = readvars_nc(nc_fname,'time')
    ib_min = readvars_nc(nc_fname,'IB')/100
    

    nc_fname = folder_tmp  + 'wind_and_waves.nc'
    time_w = readvars_nc(nc_fname,'time')
    Hs = readvars_nc(nc_fname,'Hs')
    #Tm = readvars_nc(nc_fname,'Tm')
    Fspr = readvars_nc(nc_fname,'Fspr')
    Tp = readvars_nc(nc_fname,'Tp')
    Dir = readvars_nc(nc_fname,'Dir')
    Dspr = readvars_nc(nc_fname,'Dspr')
        
    
    Ru_ = np.empty((len(points),len(Hs[0,:])))  # Calculate wave runup all points
    for i in range(0, len(points)):
        Ru_[i,:] = estimate_wave_runup(Hs[i,:],Tp[i,:],Dir[i,:],points[i,2])
       
        
    
    Sla_ = sla + ib # add inverted baromenter to sea level anomaly 
    Zetaoff_ = Sla_ + tide # add tide 
    
    TWL_ = np.tile(Sla_,(len(points),1)) + np.tile(tide,(len(points),1)) + Ru_ # add wave runup
    
    
    #max_TWL_ = np.nanmax(TWL_,axis=1) 
    index_max = np.nanargmax(TWL_,axis=1) # find indices of maximum TWL at each point

     
    # store wave and sea levels values leading to the maximun TWL at each point 
    max_TWL_ = np.zeros(npoints)
    max_hs = np.zeros(npoints)
    max_tp = np.zeros(npoints)
    max_fsp = np.zeros(npoints)
    max_dp = np.zeros(npoints)
    max_dsp = np.zeros(npoints)
    max_zetaoff = np.zeros(npoints)
    
    for i in range(npoints):
        max_hs[i] = Hs[i,index_max[i]]
        max_tp[i] = Tp[i,index_max[i]]
        max_fsp[i] = Fspr[i,index_max[i]]
        max_dp[i] = Dir[i,index_max[i]]
        max_dsp[i] = Dspr[i,index_max[i]]
        max_zetaoff[i] = Zetaoff_[index_max[i]]
        max_TWL_[i] = TWL_[i,index_max[i]]
       

    
    Risk_Cat=np.zeros(TWL_.shape)
    for ts in range(len(TWL_[:,0])):
        Risk_Cat[ts,(TWL_[ts,:]>=thresholds[ts,0])]=2
        Risk_Cat[ts,(TWL_[ts,:]>=thresholds[ts,0]) & (TWL_[ts,:]<=thresholds[ts,4])]=1      

    Risk_Cat_max =  np.max(Risk_Cat,axis=1)

    # write results into csv
    df = pd.DataFrame({"lon": points[:,0], "lat": points[:,1], "Zeta_max[m]": max_zetaoff,
                       "Hs_max[m]": max_hs, "Tp_max[s]":max_tp,"Fsp_max[]":max_fsp,
                       "Dp_max[deg north]":max_dp,"Dsp_max[deg]":max_dsp,
                       "TWL_nearshore_max[m]": max_TWL_, "Risk[0=no risk,1=minor,2=moderate]": Risk_Cat_max})
    df.to_csv(result_folder+'Points_results.csv',index=False)
            
    riskCat2geojson(df,result_folder)
  
    if fmap==1:

        plot_risk_maps(points[:,0:2],time_w,Risk_Cat.T,thresholds,result_folder)
    
    if fseries==1:
        risk_time_series(now,TWL_.T,time_w,time_tide_min,Sla_,tide_min,thresholds,result_folder)
      
 
    
 
    
 
    
 
    
 
    
 