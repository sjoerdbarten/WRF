""" WRF_columndensity.py

Author: Auke Visser
Last modified: 29-9-2017

This script calculates the total column density for a 
chemical species from WRF-Chem output to facilitate
comparison with satellite-retrieved columns.

This script calculates cell mass from differences between
staggered pressure levels, which, combined with the molar
mass, results in the mole number of air in model cells.
This can then be multiplied by the NO2 mixing ratio and
Avogadro's number to obtain the number of cells per unit
area.

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
    #outfiles = outfiles[0:-1]
    Ma = 0.02884 #kg/mol
    N_A = 6.022e23 #molec./mol
    A = 27000.*27000. #m2/WRF cell
    g = 9.81
    t = 18
    days = 31
    
    #wrf_column_stack = np.zeros((len(outfiles),170,170)) #SHOULD BE CHANGED DEPENDING ON THE SELECTED TIME RANGE!
    wrf_column_stack = np.zeros((days,99,99)) #SHOULD BE CHANGED DEPENDING ON THE SELECTED TIME RANGE!
   
    #Retrieve monthly average WRF NO2 columns
    for i in range(0,days):
        print('Currently processing: day '+ str(i))
        ds = nc.Dataset(outfiles[0],'r')
        
        #Access variables
        wrf_no2 = ds.variables['no2'][i*24+t,0:max_lev,:,:]
        wrf_Pres = ds.variables['P'][i*24+t,0:max_lev,:,:] + ds.variables['PB'][i*24+t,0:max_lev,:,:]
        wrf_T = (ds.variables['T'][i*24+t,0:max_lev,:,:] + 300.) * ((wrf_Pres / 100000.) ** 0.2854)
        wrf_Psfc = ds.variables['PSFC'][i*24+t,:,:]
        wrf_zs = ds.variables['HGT'][0,:,:]
        
        #Calculate staggered pressure, i.e. not at the center of the cell
        #but at the top and bottom of the cells.
        Ps = np.zeros((wrf_Pres.shape[0]+1,wrf_Pres.shape[1],wrf_Pres.shape[2]))
        Ps[0,:,:] = wrf_Psfc
        for j in range(1,Ps.shape[0]):
            Ps[j,:,:] = wrf_Pres[j-1,:,:] - Ps[j-1,:,:] + wrf_Pres[j-1,:,:]
        
        dP = Ps[0:-1,:,:] - Ps[1:,:,:]
        
        #Calculate cell masses
        m_c = dP * A / g #Divide by area of grid cell
        n_a = m_c / Ma
        
        no2_molec = n_a * 1.e-6 * wrf_no2 * N_A / ((A * 1.e4)*(27000.*27000./(20000.*20000.)))
        
        wrf_column_stack[i,:,:] = np.sum(no2_molec, axis=0)
        
    return wrf_Pres, wrf_column_stack
