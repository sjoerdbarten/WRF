import netCDF4 as nc
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ncfile = nc.Dataset('/home/WUR/barte035/WRFChem/WRFV3/run/output/wrfout_d01_2014-01-01_00:00:00','r')

daysim = 30 #30
timesim = daysim*24
timearray = np.arange(0,timesim+1)

ISOP = np.zeros(shape=(timesim+1,99,99))
for i in range(0,timesim):
	ISOP[i,:,:] = ncfile.variables['EBIO_ISO'][i,:,:]
	
ISOP_Total = np.sum(ISOP,axis=0)

lat = ncfile.variables['XLAT'][timesim,:,:]
lon = ncfile.variables['XLONG'][timesim,:,:]
we = len(ncfile.dimensions['west_east'])
sn = len(ncfile.dimensions['south_north'])
dx = 27000.
width = (we*dx)-1
height = (sn*dx)-1
lat_0 = 4.891087
lon_0 = -71.069204

plt.figure(figsize=(15,9))
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
m.imshow(ISOP_Total,cmap='plasma_r',interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(-20,30,10),labels=[True,False,True])
m.drawmeridians(np.arange(-90,-50,20),labels=[True,False,True])
c = plt.colorbar(extend='max')
c.set_label(r'Isoprene flux [mol km$^{-2}$ month$^{-1}$]',size=14)
plt.show()
