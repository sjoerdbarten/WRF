import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


ncfile = nc.Dataset('/home/WUR/barte035/WRFChem/WRFV3/run/wrfout_d01_2014-01-01_00:00:00','r')
daysim = 10
timesim = daysim*24
IC_LNOx = ncfile.variables['IC_FLASHCOUNT'][:,:,:] #flashes
CG_LNOx = ncfile.variables['CG_FLASHCOUNT'][:,:,:] #flashes
N_IC = 500. #moles emitted per lightning strike IC
N_CG = 500. #moles emitted per lightning strike CG


times = np.zeros(timesim)
ICstrikes =  np.zeros(timesim,dtype=object)
CGstrikes =  np.zeros(timesim,dtype=object)
ICstrikesdomain =  np.zeros(timesim,dtype=object)
CGstrikesdomain =  np.zeros(timesim,dtype=object)

for i in range (0,timesim):
	times[i] = i
	ICstrikes[i] = sum(IC_LNOx[i,:,:])
	CGstrikes[i] = sum(CG_LNOx[i,:,:])
	ICstrikesdomain[i] = sum(ICstrikes[i])*N_IC*14.0067
	CGstrikesdomain[i] = sum(CGstrikes[i])*N_IC*14.0067
	


plt.plot(times,ICstrikesdomain)
plt.plot(times,CGstrikesdomain)
plt.show()
