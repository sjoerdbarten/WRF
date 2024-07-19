import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import csv
import netCDF4 as nc
from os.path import dirname, join as pjoin
import scipy.io as sio
import pandas as pd
from matplotlib.lines import Line2D
import tkinter

def MakeDomainMap():
    #Define basemap
    geopath = '/archive/ESG/barte035/MOSAiC/WarmAirIntrusion/input_61levels/geo_em.d01.nc'
    geo = nc.Dataset(geopath,'r')
    wrflat = geo.variables['XLAT_M'][0,:,:]
    wrflon = geo.variables['XLONG_M'][0,:,:]
    we  = geo.variables['XLAT_M'].shape[2]
    sn  = geo.variables['XLAT_M'].shape[1]
    lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
    lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
    dx_dom = 27000
    wrffile = '/archive/ESG/barte035/MOSAiC/WarmAirIntrusion/input_61levels/wrflowinp_d01'
    ds = nc.Dataset(wrffile,'r')
    
    shippos_mosaic = pd.read_csv("./df_for_trj_file_SCM.csv",
                                sep='\t', skiprows=0, header=None, skipinitialspace=True)
    shippos_mosaic = shippos_mosaic.loc[(shippos_mosaic[3] == 4)]
    shippos_mosaic = shippos_mosaic.loc[4056+(24*4):4056+(24*4)+(20*24)]
    print(shippos_mosaic)
    lon_path = shippos_mosaic[10].values
    lat_path = shippos_mosaic[9].values
    #print(lon_path)
    #print(lat_path)
    
    wrfout = nc.Dataset('/archive/ESG/barte035/MOSAiC/WarmAirIntrusion/input_61levels/wrfout/wrfout_warmairintrusion_chem_v9_d01_2020-04-05_00:00:00','r')

       
    fig, ax = plt.subplots(figsize=(24,24))    
    m = Basemap(width=7750000,height=7750000,          #6750000
            resolution='l',projection='stere',\
            lat_ts=84.3,lat_0=84.3,lon_0=0,
	    #llcrnrlon=wrfout.variables['XLONG'][0,0,0], llcrnrlat=wrfout.variables['XLAT'][0,0,0], urcrnrlon=wrfout.variables['XLONG'][0,-1,-1], urcrnrlat=wrfout.variables['XLAT'][0,-1,-1]
	    )
    
    x1, y1 = m(wrflon[0,:], wrflat[0,:])
    x2, y2 = m(wrflon[:,0], wrflat[:,0])
    x3, y3 = m(wrflon[:,-1], wrflat[:,-1])
    x4, y4 = m(wrflon[-1,:], wrflat[-1,:])
    x_path,y_path= m(lon_path,lat_path)

    geo2 = nc.Dataset('/archive/ESG/barte035/MOSAiC/WarmAirIntrusion/input_61levels/geo_em.d02.nc','r')
    geo3 = nc.Dataset('/archive/ESG/barte035/MOSAiC/WarmAirIntrusion/input_61levels/geo_em.d03.nc','r')   
    x1d2, y1d2 = m(geo2.variables['XLONG_M'][0,0,:], geo2.variables['XLAT_M'][0,0,:])
    x2d2, y2d2 = m(geo2.variables['XLONG_M'][0,:,0], geo2.variables['XLAT_M'][0,:,0])
    x3d2, y3d2 = m(geo2.variables['XLONG_M'][0,:,-1], geo2.variables['XLAT_M'][0,:,-1])
    x4d2, y4d2 = m(geo2.variables['XLONG_M'][0,-1,:], geo2.variables['XLAT_M'][0,-1,:])    
    x1d3, y1d3 = m(geo3.variables['XLONG_M'][0,0,:], geo3.variables['XLAT_M'][0,0,:])
    x2d3, y2d3 = m(geo3.variables['XLONG_M'][0,:,0], geo3.variables['XLAT_M'][0,:,0])
    x3d3, y3d3 = m(geo3.variables['XLONG_M'][0,:,-1], geo3.variables['XLAT_M'][0,:,-1])
    x4d3, y4d3 = m(geo3.variables['XLONG_M'][0,-1,:], geo3.variables['XLAT_M'][0,-1,:])

    lonswrf,latswrf = m(geo.variables['XLONG_M'][0,:,:],geo.variables['XLAT_M'][0,:,:],inverse=False)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.1,1.5],colors='w',zorder=2,alpha=0.1)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.2,1.5],colors='w',zorder=2,alpha=0.2)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.3,1.5],colors='w',zorder=2,alpha=0.3)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.4,1.5],colors='w',zorder=2,alpha=0.4)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=2,alpha=0.5)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.6,1.5],colors='w',zorder=2,alpha=0.6)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.7,1.5],colors='w',zorder=2,alpha=0.7)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.8,1.5],colors='w',zorder=2,alpha=0.8)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.9,1.5],colors='w',zorder=2,alpha=0.9)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[1.0,1.5],colors='w',zorder=2,alpha=1.0)   
    m.contourf(lonswrf,latswrf,wrfout.variables['SNOWC'][0,:,:],[0.5,1.5],colors='w',zorder=2,alpha=0.9)
    m.plot(x1,y1, color='black', linewidth=2) 
    m.plot(x2,y2, color='black', linewidth=2) 
    m.plot(x3,y3, color='black', linewidth=2) 
    m.plot(x4,y4, color='black', linewidth=2) 
    m.plot(x1d2,y1d2, color='black', linewidth=2) 
    m.plot(x2d2,y2d2, color='black', linewidth=2) 
    m.plot(x3d2,y3d2, color='black', linewidth=2) 
    m.plot(x4d2,y4d2, color='black', linewidth=2) 
    m.plot(x1d3,y1d3, color='black', linewidth=2) 
    m.plot(x2d3,y2d3, color='black', linewidth=2) 
    m.plot(x3d3,y3d3, color='black', linewidth=2) 
    m.plot(x4d3,y4d3, color='black', linewidth=2)
    m.plot(x_path,y_path, linewidth=3, color='black')
    m.drawcoastlines(linewidth=1.0)
    m.drawcountries(linewidth=1.0)
    m.fillcontinents(color='peru',lake_color='cornflowerblue')
    m.drawparallels(np.arange(-80.,81.,10.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    plt.annotate('d01', xy=(0.92, 0.91), xycoords='axes fraction', fontsize=30)
    plt.annotate('d02', xy=(0.62, 0.62), xycoords='axes fraction', fontsize=30)
    plt.annotate('d03', xy=(0.53, 0.53), xycoords='axes fraction', fontsize=30)
    #m.etopo()
    m.drawmapboundary(fill_color='cornflowerblue')
    fig.tight_layout()
    plt.savefig('2020_domainmap.png', format='png', dpi=100)
    plt.show()
    
MakeDomainMap()
