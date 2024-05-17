#
# in this python script, the land use is changed
# in the geo_em files
#
#
import os
import netCDF4 as nc
import numpy as np
import pwd

home = os.path.expanduser("~") # set home directory
datapath = "%s/atmmod/WPS"%home # set full datapath

domain = 'd01' # select domain

geo = nc.Dataset(datapath+"/geo_em."+domain+".nc","r+") # open the geo_em file

LU = geo["LU_INDEX"][:]
LU_update = np.copy(LU)
LEFF = geo["LANDUSEF"][:]
LEFF_update = np.copy(LEFF)

tar_lu = 16 # select the land use to which it needs to be adapted

for iy in range(0,LU.shape[1]): # loop in y-direction
   for ix in range(0,LU.shape[2]): # loop in x-direction
      
#
# select grid cells where the dominant land
# use is water
#
      if LU[0,iy,ix] == 17:
         LEFF_update[:,:,iy,ix] =0.

# change LANDUSEF and LU_INDEX
# that water points are replaced
# by land use type defined in tar_lu
#
      else:  
         LEFF_update[:,:,iy,ix] = 0.
         LEFF_update[:,tar_lu-1,iy,ix] = 1.
         LU_update[:,iy,ix] = tar_lu
	  
geo["LU_INDEX"][:] = LU_update[:]
geo["LANDUSEF"][:] = LEFF_update[:]
geo.close()
