"""WRF_spatialplot_avg.py

Averages wrfout files at one time instance
and plots results on a map
"""
import os
import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#def plot_spatial_avg(path,plotvar,plottime,figsave):
path = '/home/WUR/barte035/WRFChem/WRFV3/run/'
plotvar = 'iso'
plottime = 24.
filelist = [os.path.join(path,filename) for filename in os.listdir(path) if filename.startswith('wrfout')]
filelist = filelist[0:-1]
filelist.sort()

#Define/extract dimension variables needed for plotting
ncfile = nc.Dataset(filelist[0],'r')
lat = ncfile.variables['XLAT'][plottime,:,:]
lon = ncfile.variables['XLONG'][plottime,:,:]
we = len(ncfile.dimensions['west_east'])
sn = len(ncfile.dimensions['south_north'])
dx = 20000
width = (we*dx)-1
height = (sn*dx)-1
lat_0 = 4.891087
lon_0 = -71.069204
ncfile.close()

if plotvar == 'DEP_FLX_O3':
    axlabel = r'$F_{O_{3},tot} \: [\mu mol]$'
elif plotvar == 'no2':
    axlabel = r'[NO$_2$] [ppb]'
elif plotvar == 'VD_O3':
    axlabel = r'$V_d(O_3) \: [cm \: s^{-1}]$'
elif plotvar == 'o3':
    axlabel = r'$[O_3] \: [ppb]$'
elif plotvar == 'T':
    axlabel = r'$T \ [K]$'
elif plotvar == 'T2':
    axlabel = r'$T2m \: [K]$'
elif plotvar == 'U10':
    axlabel = r'$u \: wind \: component \: [m \: s^{-1}]$'
elif plotvar == 'PBLH':
    axlabel = r'$PBL \: height \: [m]$'
elif plotvar == 'MEBIO_NO':
    axlabel = r'$MEGAN biogenic NO emission [mol \: km^{-2} \: hr^{-1}]$'
elif plotvar == 'EBIO_ISO':
    axlabel = r'$MEGAN biogenic NO emission [mol \: km^{-2} \: hr^{-1}]$'
elif plotvar == 'iso':
    axlabel = r'$MEGAN biogenic NO emission [mol \: km^{-2} \: hr^{-1}]$'

var_stack = np.zeros((len(filelist),we,sn))

#Loop over different files and extract the variable of interest
for i in range(0,len(filelist)):
    ncfile = nc.Dataset(filelist[i],'r')
    if plotvar in ['no2','o3','T','iso']:
       #var_stack[i,:,:] = np.nanmean(ncfile.variables[plotvar][plottime,0:5,:,:],axis=0)
        var_stack[i,:,:] = ncfile.variables[plotvar][plottime,0,:,:]
    else:
        var_stack[i,:,:] = ncfile.variables[plotvar][plottime,:,:]

if plotvar == 'DEP_FLX_O3':
    var_stack *= 1.e6
elif plotvar in ['no2','o3']:
    var_stack *= 1.e3

#Average the values
var_avg = np.nanmean(var_stack,axis=0)

#Generate a plot
plt.figure(figsize=(15,9))
m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,projection='lcc',lat_0=lat_0,lon_0=lon_0)
if plotvar == 'DEP_FLX_O3':
    m.imshow(var_avg,cmap='jet',vmin=0,vmax=3,interpolation='none')
elif plotvar == 'no2':
    m.imshow(var_avg,cmap='jet',vmin=0,vmax=0.05,interpolation='none')
elif plotvar == 'VD_O3':
    m.imshow(var_avg,cmap='jet',vmin=0,vmax=1.2,interpolation='none')
elif plotvar == 'o3':
    m.imshow(var_avg,cmap='jet',vmin=0,vmax=50,interpolation='none')
elif plotvar == 'T':
    m.imshow(var_avg,cmap='jet',vmin=273,vmax=323,interpolation='none')
elif plotvar == 'T2':
    m.imshow(var_avg,cmap='jet',vmin=273,vmax=305,interpolation='none')
elif plotvar == 'U10':
    m.imshow(var_avg,cmap='jet',vmin=0,vmax=10,interpolation='none')
elif plotvar == 'PBLH':
    m.imshow(var_avg,cmap='jet',vmin=0,vmax=3000,interpolation='none')
elif plotvar == 'MEBIO_NO':
    m.imshow(var_avg,cmap='jet',interpolation='none')
elif plotvar == 'EBIO_ISO':
    m.imshow(var_avg,cmap='jet',vmin=-0.00000000000000000000000000001,vmax=0.00000000000000000000000000001,interpolation='none')
elif plotvar == 'iso':
    m.imshow(var_avg,cmap='jet',vmin=0,vmax=100,interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(-20,30,10),labels=[True,False,True])
m.drawmeridians(np.arange(-90,-50,20),labels=[True,False,True])
if plotvar == 'DEP_FLX_O3':
    c = plt.colorbar(ticks=[0,1,2,3,4,5])
elif plotvar == 'no2':
    c = plt.colorbar(ticks=np.arange(0,0.051,0.01))
elif plotvar == 'VD_O3':
    c = plt.colorbar(ticks=np.arange(0,1.4,0.2))
elif plotvar == 'o3':
    c = plt.colorbar(ticks=np.arange(0,51,10))
elif plotvar == 'T':
    c = plt.colorbar(ticks=np.arange(280,321,10))
elif plotvar == 'T2':
    c = plt.colorbar(ticks=np.arange(270,310,5))
elif plotvar == 'U10':
    c = plt.colorbar(ticks=np.arange(0,11,2))
elif plotvar == 'PBLH':
    c = plt.colorbar(ticks=np.arange(0,3001,1000))
elif plotvar == 'MEBIO_NO':
    c = plt.colorbar(ticks=np.arange())
elif plotvar == 'EBIO_ISO':
    c = plt.colorbar(ticks=np.arange(-0.00000000000000000000000000001,0.00000000000000000000000000001,0.000000000000000000000000000001))
elif plotvar == 'iso':
    c = plt.colorbar(ticks=np.arange(0,100,10))
c.set_label(axlabel,size=18)

plt.savefig('/home/WUR/barte035/WRFChem/Python/Figures/Example_'+ str(plotvar) +'_time'+ str(int(plottime)) +'.png',bbox_inches='tight')

plt.show()
