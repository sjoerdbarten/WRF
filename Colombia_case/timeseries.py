import netCDF4 as nc
from pylab import *
import matplotlib.pyplot as plt
import numpy as np


print 'Use timeseries(variable,location). Variable = no2, o3. Location: 1=Bogota, 2=Orinoco. Or use timeseriesboth(variable) for both locations in 1 plot.'

def timeseries(variable,location):
	ncfile = nc.Dataset('/home/WUR/barte035/WRFChem/WRFV3/run/wrfout_d01_2014-01-01_00:00:00','r')
	if location == 1:
		gridy=48
		gridx=34
	if location == 2:
		gridy=53
		gridx=50
	plotvar=ncfile.variables[variable][:,0,gridy,gridx]
	plot(plotvar)
	show()

def timeseriesboth(variable):
	ncfile = nc.Dataset('/home/WUR/barte035/WRFChem/WRFV3/run/wrfout_d01_2014-01-01_00:00:00','r')
	plotvarBog=ncfile.variables[variable][:,0,48,34]*1000
	plotvarOri=ncfile.variables[variable][:,0,53,50]*1000
	plt.plot(plotvarBog)
	plt.plot(plotvarOri)
	plt.xlabel('time since start of simulation [h]', fontsize=14)
	plt.ylabel(''+str(variable)+' concentration [ppbv]', fontsize=14)
	plt.legend(['Bogota','Orinoco'], loc='best')
	plt.savefig('Figures/timeseries'+str(variable)+'',dpi=300)
	plt.show()
