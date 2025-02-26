# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 16:26:56 2024

@author: antonioh
"""


# import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import os
import pyproj
import xarray as xr
import rasterio  
from rasterio.transform import  from_origin

def geotiff(bckground_file):
    # Open the GeoTIFF file
    with rasterio.open(bckground_file) as dataset:
        # Read the RGB bands (assume the first three bands are Red, Green, Blue)
        red = dataset.read(1)   # Red band
        green = dataset.read(2) # Green band
        blue = dataset.read(3)  # Blue band
    
        # Stack the RGB bands into a single array with 3 channels
        rgb = np.stack([red, green, blue], axis=-1)  # Shape will be (4961, 2082, 3)
    
        # Normalize the values to the range 0-255 for display
        rgb_normalized = (rgb / rgb.max()) * 255
        mapa = rgb_normalized.astype(np.uint8)
    
        # Get the extent (coordinates) to plot it correctly
        extent = [dataset.bounds.left, dataset.bounds.right, dataset.bounds.bottom, dataset.bounds.top]

    return(mapa,extent)


def grid2tiff(lgrid_x,lgrid_y,var,fullfilepath,epsg):
    

    fullfilepath = fullfilepath + '.tif'
    x0 = lgrid_x.min()
    y0 = lgrid_y.min()
    deltax =  np.abs(lgrid_x[0,0]-lgrid_x[0,1])
    deltay = -np.abs(lgrid_y[0,0]-lgrid_y[1,0])
    x0 = lgrid_x.min()-deltax#/2
    y0 = lgrid_y.min()+deltay#/2
    
    transform = from_origin(x0,y0,deltax,deltay)
    
    ds_tif = rasterio.open(fullfilepath,
                       'w+',
                        driver = 'GTiff',
                        height = lgrid_x[:,0].size,
                        width =  lgrid_x[0,:].size,
                        dtype =  var.dtype ,
                        count = 1,
                        crs = epsg,
                        nodata=-9999,
                        transform = transform)
    ds_tif.write(var, 1)
    ds_tif.close()
    
    return()






def reproject(netcdf_in,netcdf_out,epgs_in,epgs_out):
  
    ds = xr.open_dataset(netcdf_in)
    
    # Define the source and target coordinate reference systems (CRS)
    src_crs = pyproj.CRS('EPSG:32704')
    tgt_crs = pyproj.CRS('EPSG:4326')
    
    # Transform the coordinates
    transformer = pyproj.Transformer.from_crs(src_crs, tgt_crs, always_xy=True)
    
    
    # Transform the coordinates
    transformer = pyproj.Transformer.from_crs(src_crs, tgt_crs, always_xy=True)
    lon, lat = transformer.transform(ds['x'].values, ds['y'].values)
    
    # Assign the transformed coordinates to new variables
    ds.coords['lon'] = (('y', 'x'), lon)
    ds.coords['lat'] = (('y', 'x'), lat)
    
    ds.to_netcdf(netcdf_out)
    
    
    # # Plot the transformed coordinates
    # plt.figure(figsize=(10, 6))
    # plt.scatter(ds['lon'], ds['lat'], c=ds['hmax'], cmap='viridis')
    # plt.colorbar(label='Variable')
    # plt.xlabel('Longitude')
    # plt.ylabel('Latitude')
    # plt.title('Reprojected NetCDF Data')
    
    # plt.show()
    return()

###############################################################################

def plot_SFINCS(now):

    
    out_name ='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  +'/'
    result_folder = '../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  + '/results/'
    #templ_path  ='../extras/rarotonga/sfincs'
    run_path = out_name + 'SFINCS'
    XBrun_path = out_name + 'XBeach'
    
    # Background map
    bckground = '../extras/rarotonga/gis/osm_rarotonga.tif'
    Map,extent = geotiff(bckground)

    
    # Loading xbeach outputs
    keylist = ['name', 'zsmax',  'x', 'y']
    data = { key: [] for key in keylist }
    
    
    xbdata = []
    xblist = os.listdir(XBrun_path)
    for ii, name in enumerate(xblist):
        fpath = os.path.join(XBrun_path, name, 'xboutput.nc')
        tmp = data.copy()
        
        with xr.open_dataset(fpath) as fncdf:
            tmp['name'] = name
            tmp['zsmax'] = np.nanmax(fncdf['zs'].data.squeeze(), axis=0)
            tmp['x'] = fncdf.coords['globalx'].data.squeeze()
            tmp['y'] = fncdf.coords['globaly'].data.squeeze()
            
        xbdata.append(tmp)
        
    netcdf_in = run_path + '/sfincs_map.nc'
    netcdf_out =  result_folder + '/sfincs_map_reproj.nc'
        
    reproject(netcdf_in,netcdf_out,'EPSG:32704','EPSG:4326') 
        
    # SFINCS inundation map

    with xr.open_dataset(netcdf_in) as fncdf:
        data = { key: item.data.squeeze() for key, item in fncdf.items() }
        coords = { key: item.data.squeeze() for key, item in fncdf.coords.items() }
    
    # Creating a mask to remove the data beneath MSL
    hmax = data['zsmax'].copy()
    # Setting up the levels used to plot the data
    levels = np.arange(0, 4.1, 0.1)
    msize = 10
    sc_kwargs = {'cmap': 'jet',
                 'vmin': levels.min(),
                 'vmax': levels.max(),
                 'marker': 'o'}
    
    
    ## Plot maximum water elevation (above MSL) of SFCINS and XBeach for verification
    fig, axs = plt.subplots(1,1, layout='compressed',figsize=(15 ,12))
    # Ploting the background map
    axs.imshow(Map, extent=extent)
    axs.set_xlim(412645.42284697585, 425898.8303943549)
    axs.set_ylim(7646583.127815168, 7656352.306208582)
    # Ploting maximun water level
    im = axs.contourf(coords['x'], coords['y'], hmax, levels, cmap='jet')
    # Ploting the transects
    for ii, v in enumerate(xbdata):
        axs.scatter(v['x'], v['y'], s=msize, c=v['zsmax'], **sc_kwargs)
    
    
    axs.set(xlabel='x-UTM [m]', ylabel='y-UTM [m]', title='Maximun water elevation above MSL [m]')
    axs.legend()
    fig.colorbar(im, ax=axs, label='Maximum water level [m]', location='right')
    # fig.suptitle(run)
    # plt.show()
    fileprint = result_folder + 'Rarotonga_XB_SFINCS'
    plt.savefig(fileprint)
    plt.close('all')
    
    # Creating a mask to remove the data beneath MSL
    msk_b = data['zb'] < 0.2
    hmax = data['hmax'].copy()
    hmax[msk_b] = np.NaN
    # Setting up the levels used to plot the data
    levels = np.arange(0, 4.1, 0.1)
    msize = 10
    sc_kwargs = {'cmap': 'jet',
                 'vmin': levels.min(),
                 'vmax': levels.max(),
                 'marker': 'o'}
    
    
    ## Plot inundation extent
    fig, axs = plt.subplots(1,1, layout='compressed',figsize=(15 ,12))
    # Ploting the background map
    axs.imshow(Map, extent=extent)
    axs.set_xlim(412645.42284697585, 425898.8303943549)
    axs.set_ylim(7646583.127815168, 7656352.306208582)
    # Ploting the inundation extent
    im = axs.contourf(coords['x'], coords['y'], hmax, levels, cmap='jet')
    axs.set(xlabel='x-UTM [m]', ylabel='y-UTM [m]', title='Inundation depth [m]')
    axs.legend()
    fig.colorbar(im, ax=axs, label='Maximum water level [m]', location='right')
    # fig.suptitle(run)
    #plt.show()
    fileprint = result_folder + 'Rarotonga_inundation_depth'
    plt.savefig(fileprint)
    plt.close('all')
    

    grid2tiff(coords['x'],coords['y'],hmax,fileprint,'EPSG:32704')








