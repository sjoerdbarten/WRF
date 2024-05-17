import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import wrf
np.bool = np.bool_

#--- User settings---#


# for which point (latitude,longitude) the pbl is computed
lat_1 = 51.97
lon_1 = 4.92

# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,00)

datestring = str(yy)+'-'+str(mm).zfill(2)+'-'+str(dd).zfill(2)+'_'+str(HH).zfill(2)+':00:00'

# which domain is selected

domain        = "d02"

#--- User settings---#

home          = os.path.expanduser("~") 
infilename = home+'/atmmod/WRF/run/wrfout_d01_2018-06-05_00:00:00'
print(type(infilename))
print('opening %s ...'%infilename)

infile = nc.Dataset(infilename,'r')

xloc,yloc = wrf.to_np(wrf.ll_to_xy(infile,lat_1,lon_1))

print('')
print('**********************************************')
print('You have entered the following coordinates:')
print('Latitude = ',lat_1)
print('Longitude = ',lon_1)
print('**********************************************')
print('This corresponds to the follwing grid cells in the WRF output:')
print('yloc =', yloc)
print('xloc = ',xloc)
print('**********************************************')
print('This corresponds to the follwing coordinates in the WRF output:')
print('Latitude:', wrf.getvar(infile,'lat',timeidx=0).values[yloc,xloc])
print('Longitude:', wrf.getvar(infile,'lon',timeidx=0).values[yloc,xloc])
print('**********************************************')
print('')
print('***Successfully finished determine_ylocxloc.py script***')
