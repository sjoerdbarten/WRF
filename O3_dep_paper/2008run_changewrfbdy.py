import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import datetime as dt
import time

filename = './wrfbdy_d01_test'

wrfinput = nc.Dataset(filename,'r+',format='NETCDF4')

varnames = wrfinput.variables.keys()[3:]
print(wrfinput.variables.keys())
varnames.remove('o3_BTXE')
varnames.remove('o3_BTXS')
varnames.remove('o3_BTYE')
varnames.remove('o3_BTYS')
varnames.remove('o3_BXE')
varnames.remove('o3_BXS')
varnames.remove('o3_BYE')
varnames.remove('o3_BYS')
print(varnames)

print(len(wrfinput.variables.keys()))
print(len(varnames))

for i in range(len(varnames)):
	print('Processing variable '+str(i)+' out of '+str(len(varnames))+': '+str(varnames[i]))
	#Change following timesteps
	wrfinput.variables[varnames[i]][20,:,:] = wrfinput.variables[varnames[i]][19,:,:]
	wrfinput.variables[varnames[i]][37,:,:] = wrfinput.variables[varnames[i]][38,:,:]
	wrfinput.variables[varnames[i]][46,:,:] = wrfinput.variables[varnames[i]][47,:,:]
	wrfinput.variables[varnames[i]][87,:,:] = wrfinput.variables[varnames[i]][86,:,:]
	wrfinput.variables[varnames[i]][88,:,:] = wrfinput.variables[varnames[i]][89,:,:]

print('Processing ozone variables:')
wrfinput.variables['o3_BTXE'][26,:,:] = wrfinput.variables['o3_BTXE'][25,:,:]
wrfinput.variables['o3_BTXS'][26,:,:] = wrfinput.variables['o3_BTXS'][25,:,:]
wrfinput.variables['o3_BTYE'][26,:,:] = wrfinput.variables['o3_BTYE'][25,:,:]
wrfinput.variables['o3_BTYS'][26,:,:] = wrfinput.variables['o3_BTYS'][25,:,:]
wrfinput.variables['o3_BXE'][26,:,:] = wrfinput.variables['o3_BXE'][25,:,:]
wrfinput.variables['o3_BXS'][26,:,:] = wrfinput.variables['o3_BXS'][25,:,:]
wrfinput.variables['o3_BYE'][26,:,:] = wrfinput.variables['o3_BYE'][25,:,:]
wrfinput.variables['o3_BYS'][26,:,:] = wrfinput.variables['o3_BYS'][25,:,:]
	
wrfinput.close()
