import os
import netCDF4 as nc
import numpy as np
import pwd

home = os.path.expanduser("~") # set home directory
datapath = "%s/atmmod/WPS"%home # set full datapath

domain = 'd02' # select domain

geo = nc.Dataset(datapath+"/geo_em."+domain+".nc","r+") # open the geo_em file

height = geo["HGT_M"][:]
height_up = np.copy(height)

#
# change in geo files the terrain height
# be careful: the geo_em file is overwritten
#

for iy in range(0,height.shape[1]): # loop over y-direction
  for ix in range(0,height.shape[2]): # loop over x-direction
  
      height_up[:,iy,ix] =  1.1*height[:,iy,ix] # height 10 % higher

geo["HGT_M"][:] = height_up[:]

geo.close()
