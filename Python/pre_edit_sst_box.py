#
# in this python script, the sst is changed
#
#

import os
import netCDF4 as nc
import numpy as np
import pwd

home = os.path.expanduser("~") # set home directory
datapath = "%s/atmmod/WRF/run"%home # set full datapath

domain = 'd02' # select domain

sst_file = nc.Dataset(datapath+"/wrfinput_"+domain,"r+")
min_lat = 51.
max_lat = 52.
min_lon = 5.
max_lon = 6.

LU = sst_file["LU_INDEX"][:]
lat = sst_file['XLAT'][:]
lon = sst_file['XLONG'][:]
SST = sst_file["SST"]
TSK = sst_file["TSK"]
SST_update = np.copy(SST)
TSK_update = np.copy(TSK)


for iy in range(0,SST.shape[1]): # loop in y-direction
   for ix in range(0,SST.shape[2]): # loop in x-direction
      
#
# select grid cells where the dominant land
# use is water
#
      if LU[0,iy,ix] == 17: 

# change for the sea points
# the variables SST and TSK
#

         if (lat[0,iy,ix] > min_lat) & (lat[0,iy,ix] < max_lat) & (lon[0,iy,ix] > min_lon) & (lon[0,iy,ix] < max_lon):
            SST_update[:,iy,ix] = SST[:,iy,ix] + 2.
            TSK_update[:,iy,ix] = TSK[:,iy,ix] + 2.
         else:
            SST_update[:,iy,ix] = SST[:,iy,ix] + 1.
            TSK_update[:,iy,ix] = TSK[:,iy,ix] + 1.

sst_file["SST"][:] = SST_update[:]
sst_file["TSK"][:] = TSK_update[:]

sst_file.close()
