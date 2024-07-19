import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import datetime as dt
import time

filename = './wrfinput_d01_test'

wrfinput = nc.Dataset(filename,'r+',format='NETCDF4')

print(wrfinput)
print(wrfinput.variables)

co2_add = wrfinput.createVariable('co2','f4',('Time','bottom_top','south_north','west_east'))
co2_add.units = wrfinput.variables['co'].units
ch4_add = wrfinput.createVariable('ch4','f4',('Time','bottom_top','south_north','west_east'))
ch4_add.units = wrfinput.variables['co'].units

print(co2_add)
print(ch4_add)

co2_add[:,:,:] = 380.
ch4_add[:,:,:] = 1.7

print(wrfinput)
print(wrfinput.variables)

wrfinput.close()
