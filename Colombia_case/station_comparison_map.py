""" station_comparison_map.py

This script takes directory with WRF-Chem output and a 
directory with station observations, and plots mean
12:00 h UTC values on the map to facilitate model-
observation comparison.

"""

import os
import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import genfromtxt

def station_comparison_map(path,plotvar,plottime,obsdir,figsave):
    #Declare constants
    R = 8.314
    speclist = ['no2','no','o3']
    #Observations for NO and NO2 are in ug N/m3!
    molarmass = [14.01,14.01,48.]
    
    filelist = [os.path.join(path,filename) for filename in os.listdir(path) if filename.startswith('wrfout')]
    filelist.sort()
    filelist = filelist[0:-1]
    
    #Define/extract dimension variables needed for plotting
    ncfile = nc.Dataset(filelist[0],'r')
    lat = ncfile.variables['XLAT'][plottime,:,:]
    lon = ncfile.variables['XLONG'][plottime,:,:]
    we = len(ncfile.dimensions['west_east'])
    sn = len(ncfile.dimensions['south_north'])
    dx = 20000
    width = (we*dx)-1
    height = (sn*dx)-1
    lat_0 = 51.98769
    lon_0 = 5.66746
    ncfile.close()
    
    if plotvar == 'o3':
        axlabel = str(plottime+1) + r':00 h UTC [O$_3$] [$\mu$g m$^{-3}$]'
    elif plotvar == 'no2':
        axlabel = str(plottime+1) + r':00 h UTC [NO$_2$] [$\mu$g N m$^{-3}$]'
    elif plotvar == 'no':
        axlabel = str(plottime+1) + r':00 h UTC [NO] [$\mu$g N m$^{-3}$]'
    
    var_stack = np.zeros((len(filelist),we,sn))
    #var_stack = np.zeros((len(filelist),5,we,sn))
    
    #Loop over different files and extract the variable of interest
    for i in range(0,len(filelist)):
        ncfile = nc.Dataset(filelist[i],'r')
        
        var_stack[i,:,:] = np.nanmean(ncfile.variables[plotvar][plottime,0:5,:,:],axis=0)
        #var_stack[i,:,:,:] = ncfile.variables[plotvar][plottime,0:5,:,:]

        #Conversion to ug/m3
        #Pres = ncfile.variables['P'][plottime,0:5,:,:] + ncfile.variables['PB'][plottime,0:5,:,:]
        #T = (ncfile.variables['T'][plottime,0:5,:,:] + 300.) * ((Pres / 100000.) ** 0.2854)
        T2 = ncfile.variables['T2'][plottime,:,:]
        PSFC = ncfile.variables['PSFC'][plottime,:,:]
        molmass = molarmass[speclist.index(plotvar)]

        var_stack[i,:,:] *= molmass * PSFC / (R * T2)
        #var_stack[i,:,:,:] *= molmass * Pres / (R * T)
        
        
    #Average the values
    var_avg = np.nanmean(var_stack,axis=0) #* 1000.
    #var_avg = np.nanmean(np.nanmean(var_stack,axis=1),axis=0)
    
    #Obtain observations
    station_list = ['AT0030','AT0040','AT0041','BG0053','CH0002','CH0003','CZ0001','CZ0003','CZ0005','DE0002','DE0003','DE0007','DE0008','DE0009','DE0044','DK0012','DK0031','EE0009','ES0012','ES0013','ES0014','FI0037','FR0008','FR0009','FR0013','FR0018','FR0023','FR0030','GB0002','GB0006','GB0015','GB0031','GB0036','GB0045','GR0001','HU0002','IE0031','IT0001','IT0004','LV0010','LV0016','NL0007','NL0009','NL0010','NL0091','NL0644','NO0002','NO0039','NO0056','PL0002','PL0003','PL0005','SE0012','SE0018','SE0019','SE0039','SI0008','SK0002','SK0006','SK0007']
    lat_list = [48.721,47.348,47.973,41.696,46.813,47.480,49.733,49.573,49.067,52.802,47.915,53.167,50.650,54.433,51.530,55.694,56.290,59.500,39.086,41.283,41.400,62.583,48.500,49.900,49.617,48.633,44.569,45.767,55.313,54.443,57.734,52.504,51.573,52.298,38.367,46.967,53.326,42.100,45.800,56.162,57.135,52.083,53.334,51.541,52.300,51.974,58.389,62.783,60.372,51.817,50.733,54.150,58.800,57.165,57.953,59.782,45.567,48.833,49.050,47.960]
    
    #lat_list = [48.721,47.348,47.973,41.696,46.813,47.480,49.733,49.573,49.067,52.802,47.915,53.167,50.650,54.433,51.530,55.694,56.290,59.500,39.086,41.283,41.400,62.583,48.500,49.900,49.617,48.633,44.569,45.767,55.313,54.443,57.734,52.504,51.573,52.298,38.367,53.326,42.100,56.162,45.800,56.162,57.135,52.083,53.334,51.541,52.300,51.974,58.389,62.783,60.327,61.817,50.733,54.150,58.800,57.165,57.953,59.782,45.567,48.933,49.050,47.960] 
    lon_list = [15.942,15.882,13.016,24.739,6.945,8.905,16.050,15.080,13.600,10.759,7.909,13.033,10.767,12.733,12.934,12.086,8.427,25.900,-1.102,-5.867,0.717,24.183,7.133,4.633,0.183,-0.450,5.279,2.950,-3.204,-7.870,-4.774,-3.033,-1.317,-0.293,23.083,19.583,-9.899,12.633,8.633,21.173,25.906,6.567,6.277,5.854,4.500,4.924,8.252,8.883,11.078,21.983,15.733,22.067,17.383,14.782,12.403,15.472,14.867,19.583,22.267,17.861]
    
    obsfiles = [os.path.join(obsdir,filename) for filename in os.listdir(obsdir) if filename.endswith('_201507.csv')]
    obsfiles.sort()
    
    out_arr = np.zeros((len(station_list)))
    
    for i in range(len(obsfiles)):
        
        #Load data
        csv_arr = genfromtxt(obsfiles[i],delimiter=',')
        if plotvar == 'no2':
            obsdata = csv_arr[1:,1]
        elif plotvar == 'no':
            obsdata = csv_arr[1:,2]
        elif plotvar == 'o3':
            obsdata = csv_arr[1:,3]
        
        obsdata = obsdata[0:len(filelist)*24]
        out_arr[i] = np.nanmean(obsdata[plottime::24])
    
    #Generate a plot
    plt.figure(figsize=(15,9))
    m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,projection='lcc',lat_0=lat_0,lon_0=lon_0)
    if plotvar == 'o3':
        m.imshow(var_avg,cmap='jet',vmin=0,vmax=180,interpolation='none')
    elif plotvar in ['no','no2']:
        m.imshow(var_avg,cmap='jet',vmin=0.1,vmax=10,norm=LogNorm(),interpolation='none')
    m.drawcoastlines()
    m.drawcountries()
    m.drawparallels(np.arange(40,71,10),labels=[True,False,True])
    m.drawmeridians(np.arange(-20,51,20),labels=[True,False,True])
    if plotvar == 'o3':
        c = plt.colorbar(ticks=np.arange(0,181,30))
    elif plotvar in ['no','no2']:
        c = plt.colorbar(ticks=[0.01,0.05,0.1,0.5,1.,5.,10.])
        c.set_ticklabels([0.01,0.05,0.1,0.5,1.,5.,10.])
    c.set_label(axlabel,size=18)
    
    lons,lats = m(lon_list,lat_list)
    if plotvar == 'o3':
        m.scatter(lons,lats,marker='o',s=50,lw=2,c=out_arr,cmap='jet',vmin=0,vmax=180,edgecolors='k',zorder=5)
    elif plotvar in ['no','no2']:
        m.scatter(lons,lats,marker='o',s=50,lw=2,c=out_arr,cmap='jet',vmin=0.1,vmax=10,norm=LogNorm(),edgecolors='w',zorder=5)
    
    if figsave:
        plt.savefig('/home/WUR/visse218/python/figures/WRF/ModEmis/F_O3_1-14Jul.png',bbox_inches='tight')
    
    plt.show()
    return m,var_avg,lons,lats,out_arr
    