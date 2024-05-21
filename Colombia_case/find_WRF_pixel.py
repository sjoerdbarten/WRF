""" find_WRF_pixel.py

takes:
- a lat/lon pair from OMI
- a WRF netCDF dataset

calculates:
- the index corresponding to the OMI latitude 
  at a certain search distance (dist)
  
and returns:
- wrf_lat
- wrf_lon
- dist

"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy import interpolate

def find_pixel(wds,ods,omi_lat,omi_lon):
    
    #Read in variables
    #ods = nc.Dataset(OMI_path,'r')
    #wds = nc.Dataset(WRF_path,'r')
    prod = ods.groups['PRODUCT']
    detres = ods.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS']
    WRF_lat = wds.variables['XLAT'][0,:,:]
    WRF_lon = wds.variables['XLONG'][0,:,:]
    
    #Find OMI pixel that corresponds to WRF lat/lon
    #using an automatically modified search distance
    latlon_tup = (np.array([]),np.array([]))
    dist = 0.05
    
    while(latlon_tup[0].shape[0] != 1 and dist <= 0.25):
        latlon_tup = np.where((WRF_lat > (omi_lat-dist)) & (WRF_lat < (omi_lat+dist)) & \
                              (WRF_lon > (omi_lon-dist/2.)) & (WRF_lon < (omi_lon+dist/2.)) )
        
        if latlon_tup[0].shape[0] > 1:
            dist -= 0.005
        elif latlon_tup[0].shape[0] < 1:
            dist += 0.005
    
    if dist > 0.25:
        print 'No WRF-Chem pixel found nearby OMI pixel'
        pix_lat = np.nan
        pix_lon = np.nan
        return pix_lat, pix_lon, dist
    
    i_lat = latlon_tup[0][0]
    i_lon = latlon_tup[1][0]
    
    pix_lat = WRF_lat[i_lat,i_lon]
    pix_lon = WRF_lon[i_lat,i_lon]
    
    print 'pix_lat = ' + str(pix_lat) + ', pix_lon = ' + str(pix_lon) + ', dist = ' + str(dist)
    return scanl, gr_pix, dist