import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import wrf
np.bool = np.bool_

#--- User settings---#

home = os.path.expanduser('~') # set home directory
datapath = '%s/atmmod/WRF/run'%home # set full datapath

# for which point (latitude,longitude) the pbl is computed

# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,00)

# which domain is selected

domain        = "d02"

#--- User settings---#

infilename = '%s/wrfout_%s_%4d-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)

infile = nc.Dataset(infilename,'r')

qvapor = wrf.getvar(infile,"QVAPOR",-1)

print(qvapor)
