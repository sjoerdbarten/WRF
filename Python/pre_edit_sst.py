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

LU = sst_file["LU_INDEX"][:]
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

         SST_update[:,iy,ix] = SST[:,iy,ix] + 2.
         TSK_update[:,iy,ix] = TSK[:,iy,ix] + 2.

sst_file["SST"][:] = SST_update[:]
sst_file["TSK"][:] = TSK_update[:]

sst_file.close()
