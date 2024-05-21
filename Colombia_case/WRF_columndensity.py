""" WRF_columndensity.py

Author: Auke Visser
Last modified: 29-9-2017

This script calculates the total column density for a 
chemical species from WRF-Chem output to facilitate
comparison with satellite-retrieved columns.

This script applies the ideal gas law to find the number
of NO2 molecules per unit volume in a given tropospheric 
column, and then multiplies by the cell height (approxi-
mated using the scale height for pressure) to obtain the
vertical column density. 

"""

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import netCDF4 as nc
import numpy as np
import os
from mpl_toolkits.basemap import Basemap
from numpy import ma
from matplotlib.colors import LogNorm


def wrf_columns(wrf_path,max_lev):
    
    outfiles = [os.path.join(wrf_path,filename) for filename in os.listdir(wrf_path) if filename.startswith('wrfout_d01')]
    outfiles.sort()
    R = 8.31
    Ma = 0.02884
    N_A = 6.022e23
    g = 9.81
    t = 11
    
    #wrf_column_stack = np.zeros((len(outfiles),170,170)) #SHOULD BE CHANGED DEPENDING ON THE SELECTED TIME RANGE!
    wrf_column_stack = np.zeros((len(outfiles),99,99)) #SHOULD BE CHANGED DEPENDING ON THE SELECTED TIME RANGE!
   
    #Retrieve monthly average WRF NO2 columns
    for i in range(0,len(outfiles)):
        print('Currently processing: '+ str(outfiles[i]))
        ds = nc.Dataset(outfiles[i],'r')
        
        #Access variables
        wrf_no2 = ds.variables['no2'][t,0:max_lev,:,:]
        wrf_Pres = ds.variables['P'][t,0:max_lev,:,:] + ds.variables['PB'][t,0:max_lev,:,:]
        wrf_T = (ds.variables['T'][t,0:max_lev,:,:] + 300.) * ((wrf_Pres / 100000.) ** 0.2854)
        wrf_Psfc = ds.variables['PSFC'][t,:,:]
        wrf_zs = ds.variables['HGT'][0,:,:]
        
        #convert to molec/cm2
        n_air = wrf_Pres / (R * wrf_T)
        no2_molec = 1.e-6 * wrf_no2 * N_A * n_air
        
        #wrf_dz = -7290. * np.log(wrf_Pres[1:,:,:] / wrf_Pres[0:-1,:,:])
        
        #Approximate z at mass coordinates using scale height for Pressure
        #P = P_0 * exp(-z/H_p)
        #Check if this works: scale height depends on T!!!
        #wrf_z_m = -R / (Ma * g) * wrf_T * np.log(wrf_Pres / wrf_Psfc)
        #wrf_z_s = np.zeros((wrf_z_m.shape[0]+1,wrf_z_m.shape[1],wrf_z_m.shape[2]))
        #wrf_z_s[0,:,:] = ds.variables['HGT'][0,:,:]
        
        #stagger in the vertical
        #for j in range(1,wrf_z_s.shape[0]):
        #    wrf_z_s[j,:,:] = 2 * wrf_z_m[j-1,:,:] - wrf_z_s[j-1,:,:]
            
        #wrf_dz = wrf_z_s[1:,:,:] - wrf_z_s[0:-1,:,:]
        
        wrf_dzs = -R / (Ma * g) * wrf_T * np.log(wrf_Pres / wrf_Psfc)
        wrf_zg = wrf_dzs + wrf_zs
        
        wrf_z = np.zeros((wrf_dzs.shape[0]+1,wrf_dzs.shape[1],wrf_dzs.shape[2]))
        #wrf_z[0,:,:] = wrf_zs + 2 * wrf_dzs[0,:,:]
        wrf_z[0,:,:] = wrf_zs
        #print np.count_nonzero(wrf_zs < 0)
        for j in range(1,wrf_z.shape[0]):
            wrf_z[j,:,:] = wrf_z[j-1,:,:] + 2 * (wrf_zg[j-1,:,:] - wrf_z[j-1,:,:])
            
        wrf_dz = wrf_z[1:,:,:] - wrf_z[0:-1,:,:]
        #print np.count_nonzero(wrf_dz < 0)
        
        wrf_column_stack[i,:,:] = np.sum(no2_molec * wrf_dz, axis=0) / 1.e4 #conversion from molec/m2 to molec/cm2
        
    return wrf_zs,wrf_dz, wrf_Pres, wrf_column_stack
