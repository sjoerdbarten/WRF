import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import datetime as dt
import time

'''
filename = '../WRFV411-Polar_ERA5/run/wrfinput_d01'

wrfinput = nc.Dataset(filename,'r+',format='NETCDF4')

print(wrfinput)
print(wrfinput.variables)

u10wind = wrfinput.variables['U10']
uwind = wrfinput.variables['U']

#u10wind[u10wind > 1e10] = 0.
#uwind[uwind > 1e10] = 0.

#j = 124; i = 7,17,22,46,203,227,232,242
u10wind[:,124,7] = 0.
u10wind[:,124,17] = 0.
u10wind[:,124,22] = 0.
u10wind[:,124,46] = 0.
u10wind[:,124,203] = 0.
u10wind[:,124,227] = 0.
u10wind[:,124,232] = 0.
u10wind[:,124,242] = 0.

uwind[:,:,124,7] = 0.
uwind[:,:,124,17] = 0.
uwind[:,:,124,22] = 0.
uwind[:,:,124,46] = 0.
uwind[:,:,124,203] = 0.
uwind[:,:,124,227] = 0.
uwind[:,:,124,232] = 0.
uwind[:,:,124,242] = 0.

print(u10wind)
print(uwind)

print(u10wind.shape)
print(uwind.shape)

wrfinput.close()
'''


#SAME FOR WRFFDDA FILE!
filenamefdda = '../WRFV411-unpolar/run/wrffdda_d01'

wrffdda = nc.Dataset(filenamefdda,'r+',format='NETCDF4')

print(wrffdda)
print(wrffdda.variables)

uoldwind = wrffdda.variables['U_NDG_OLD']
unewwind = wrffdda.variables['U_NDG_NEW']

#j = 124; i = 7,17,22,46,203,227,232,242
uoldwind[:,:,124,7] = uoldwind[:,:,125,7]
uoldwind[:,:,124,17] = uoldwind[:,:,125,17]
uoldwind[:,:,124,22] = uoldwind[:,:,125,22]
uoldwind[:,:,124,46] = uoldwind[:,:,125,46]
uoldwind[:,:,124,203] = uoldwind[:,:,125,203]
uoldwind[:,:,124,227] = uoldwind[:,:,125,227]
uoldwind[:,:,124,232] = uoldwind[:,:,125,232]
uoldwind[:,:,124,242] = uoldwind[:,:,125,242]

unewwind[:,:,124,7] = unewwind[:,:,125,7]
unewwind[:,:,124,17] = unewwind[:,:,125,17]
unewwind[:,:,124,22] = unewwind[:,:,125,22]
unewwind[:,:,124,46] = unewwind[:,:,125,46]
unewwind[:,:,124,203] = unewwind[:,:,125,203]
unewwind[:,:,124,227] = unewwind[:,:,125,227]
unewwind[:,:,124,232] = unewwind[:,:,125,232]
unewwind[:,:,124,242] = unewwind[:,:,125,242]

print(uoldwind)
print(unewwind)

print(uoldwind.shape)
print(unewwind.shape)

wrffdda.close()

