from pylab import *
from numpy import *
import netCDF4 as nc
import sys
import subprocess
import os
from scipy import stats
import itertools
from collections import OrderedDict
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker as mtick
from mpl_toolkits.basemap import Basemap

execfile('gridlist_clean.py')

def OMI_regrid(WRFpath,OMIpath):
    WRF = nc.Dataset(WRFpath,'r')
    OMI = nc.Dataset(OMIpath,'r')
    prod = OMI.groups['PRODUCT']
    
    #extract dimensions
    im0 = 180.*2/0.25
    jm0 = 90.*2/0.25
    
    #Extract parameters from QA4ECV OMI file
    param_list = ['amf_trop']
    
    nparam = len(param_list)
    
    scanl  = prod.variables['amf_trop'].shape[1]
    gr_pix = prod.variables['amf_trop'].shape[2]
    nobs = scanl * gr_pix
    
    field = np.empty([im0,jm0,nparam])
    nrfield = np.empty([im0,jm0])
    icnt_o = 0
    
    #loop over pixels
    for iobs in range(scanl):
        for jobs in range(gr_pix):
            
            icnt_o += 1
            
            #copy corner latitudes
            tmplat4 = prod.groups['SUPPORT_DATA'].groups['GEOLOCATIONS'].variables['latitude_bounds'][0,iobs,jobs,:]
            tmplon4 = prod.groups['SUPPORT_DATA'].groups['GEOLOCATIONS'].variables['longitude_bounds'][0,iobs,jobs,:]
            
            nind, ind, jind = gridlist(tmplon4, tmplat4, im0, jm0)
            
            for i in range(nind):
                igx = iind[i]
                jg1 = jind[i]
                
                nrfield[igx,jg1] += 1.
            
            for iparam in range(nparam):
                val = prod.variables[param_list[iparam]][0,iobs,jobs]
                print val
                
            for i in range(nind):
                igx = iind[i]
                jg1 = jind[i]
                
                field[igx,jg1,iparam] += val
                
    print 'Main: Number of pixels found in this sequence = ', icnt_o

    if(icnt_o > 0):

        for ig in range(im0):
            for jg in range(jm0):

                if (nrfield[ig,jg] > min_nr):

                    for iparam in range(nparam):
                        field[ig,jg,iparam] /= nrfield[ig,jg]

                if (nrfield[ig,jg] <= min_nr):

                    for iparam in range(nparam):
                        field[ig,jg,iparam] = np.nan
                    nrfield[ig,jg] = 0.
    
    return field, nrfield