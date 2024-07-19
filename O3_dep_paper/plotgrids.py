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

def plotgrids():
    #Define basemap
    geopath = '/home/WUR/barte035/WRFChem/WPS/geo_em.d01.nc'
    geo = nc.Dataset(geopath,'r')
    wrflat = geo.variables['XLAT_M'][0,:,:]
    wrflon = geo.variables['XLONG_M'][0,:,:]
    we  = geo.variables['XLAT_M'].shape[2]
    sn  = geo.variables['XLAT_M'].shape[1]
    lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
    lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
    dx_dom = 27000
    
    m = Basemap(projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0,resolution='l')
    lonswrf,latswrf = m(geo.variables['XLAT_M'][0,:,:],geo.variables['XLONG_M'][0,:,:],inverse=False)
    m.drawcoastlines(linewidth=0.3)
    m.drawcountries(linewidth=0.3)
    m.fillcontinents(color='peru',lake_color='cornflowerblue')
    m.drawparallels(np.arange(-80.,81.,10.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='cornflowerblue')
    ax = plt.gca()
    plt.show()
    
plotgrids()
