""" OMI_read.py

This script reads in QA4ECV OMI data on a 0.125x0.125-degree regular 
grid and transforms it to the 20x20km2 equidistant WRF-Chem grid.

Author: Auke Visser
date: 18 Sep 2017

"""

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm

def OMI_read(yr):
	no2path = '/archive/ESG/barte035/Colombia/OMI/gridded/'+str(yr)+'/'
	geopath = '/archive/ESG/barte035/WPS/geo_em.d01.nc'
	no2files = [os.path.join(no2path,filename) for filename in os.listdir(no2path) if filename.endswith('.asc')]
	no2files_nc = [os.path.join(no2path,filename) for filename in os.listdir(no2path) if filename.endswith('.nc')]
   
	#Extract WRF-Chem lat/lon and regrid data
	geo = nc.Dataset(geopath,'r')
	we  = geo.variables['XLAT_M'].shape[2]
	sn  = geo.variables['XLAT_M'].shape[1]
	lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
	lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
	dx_dom = 27000

    	#m = Basemap(width=(dx_dom*we)-1,height=(dx_dom*sn)-1,resolution='l',area_thresh=1000.,projection='lcc',lat_0=lat,lon_0=lon)
    	m = Basemap(width=(we*dx_dom)-1,height=(sn*dx_dom)-1,resolution='l',area_thresh=100.,projection='lcc',lat_0=4.891087,lon_0=-71.069204,lat_1=40.,lat_2=70.,lon_1=-70.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120)
  
	#nx = int((m.xmax-m.xmin)/dx_dom)+1
	#ny = int((m.ymax-m.ymin)/dx_dom)+1
	nx = 99
	ny = 99
	x2_trans = np.linspace(-179.9375,179.9375,2880)
	y2_trans = np.linspace(-89.9375,89.9375,1440)

	no2stack_nc = np.zeros((len(no2files_nc),we,sn))
                           
	for i in range(0,len(no2files_nc)):
	   ncfile = nc.Dataset(no2files_nc[i],'r')
	   if yr == 2014:
	    	no2_nc = ncfile.variables['no2_col'][:]
	   elif yr == '2014_mod':
 	   	no2_nc = ncfile.variables['mod_no2_col'][:]
 	   no2_nc.data[no2_nc.data == 9.96920997e+36] = np.nan
                          
 	   #Transform data to WRF-Chem grid
 	   no2stack_nc[i,:,:] = m.transform_scalar(no2_nc.data,lons=x2_trans,lats=y2_trans,nx=nx,ny=ny,order=0)
 	   no2stack_nc[i,:,:][no2stack_nc[i,:,:] < 0.] = np.nan
    
	nonans = (~np.isnan(no2stack_nc)).sum(axis=0)
        
	#Generate mean of the data
	no2_nc_trans = np.nanmean(no2stack_nc,axis=0)
	no2_nc_trans[nonans < 2] = np.nan
    
	print 'Plotting mean NO2 column density from .nc'
	plt.figure(figsize=(15,9))
	#im = m.imshow(no2_nc_trans,norm=LogNorm(),vmin=0.01,vmax=30,interpolation='none')
	im = m.imshow(no2_nc_trans,vmin=0.1,vmax=8,norm=LogNorm(),interpolation='none')
	m.drawcoastlines()
	m.drawcountries()
	m.drawmapboundary(fill_color='lightgray')
	#c = plt.colorbar()
	c = plt.colorbar(extend='max')
	c.set_ticks([0.1,0.5,1,2,3,4,6,8])
	c.set_ticklabels([0.1,0.5,1,2,3,4,6,8])
	c.set_label(u'Mean NO_2 column density [$10^{15}$ molec. cm$^2$]',size=15)
	if yr == 2014:
		plt.title(u'Mean NO_2 column OMI',size=18,weight='bold')
	elif yr == '2014_mod':
		plt.title(u'Mean recalculated NO_2 column OMI',size=18,weight='bold')
	#plt.title('NO2 VCD from .nc',size=18,weight='bold')
	plt.show()
        
	#return regridded QA4ECV NO2 columns
	return no2stack_nc,no2_nc_trans
