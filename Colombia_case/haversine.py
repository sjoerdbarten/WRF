""" haversine.py

Takes:
- pixel/point lat/lon information
- WRF-Chem lat/lon arrays

Calculates:
- a distance matrix between all WRF-Chem
  pixels and this point using the
  Haversine function
  
And returns:
- the distance matrix

source: https://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
"""

import netCDF4 as nc
import numpy as np
from math import radians, cos, sin, asin, sqrt

def haversine(lon1, lat1, lon2, lat2):
    
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2. * asin(sqrt(a))
    r = 6371.
    return c * r

def haversine_array(lon, lat, ar_lon, ar_lat):
    
    dist = np.zeros((ar_lon.shape[0],ar_lon.shape[1]))
    
    for i in range(ar_lon.shape[0]):
        for j in range(ar_lon.shape[1]):
            dist[i,j] = haversine(lon, lat, ar_lon[i,j], ar_lat[i,j])
        
    i,j = np.unravel_index(dist.argmin(),dist.shape)
    
    return i, j, dist[i,j], ar_lon[i,j], ar_lat[i,j]