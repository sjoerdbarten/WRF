import netCDF4 as nc
import scipy.interpolate
from pylab import *
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def verticalprofile(plotvar,location):
	ncfile = nc.Dataset('/home/WUR/barte035/WRFChem/WRFV3/run/wrfout_d01_2014-01-01_00:00:00','r')
	timesim = 240.
	
	if location == 'Bogota':
		var = ncfile.variables[plotvar][:,:,48,34]*1000.
		p = ncfile.variables['P'][:,:,48,34] + ncfile.variables['PB'][:,:,48,34] 
	if location == 'Orinoco':
		var = ncfile.variables[plotvar][:,:,53,50]*1000.
		p = ncfile.variables['P'][:,:,53,50] + ncfile.variables['PB'][:,:,53,50]
	if location == 'Amazone':
		var = ncfile.variables[plotvar][:,:,10,85]*1000.
		p = ncfile.variables['P'][:,:,10,85] + ncfile.variables['PB'][:,:,10,85]
		
	Times = p[:,0]
	
	#array = np.zeros(shape=(240.,59.,59.))
	#
	#print array.shape
	#array[:,:,:] = np.arange(timesim)
	#print array[0]
	
	array = np.arange(timesim)
	array = array[:,np.newaxis,np.newaxis]
	np.reshape(array[:,:,:],[timesim,60,60])
	print array.shape
	print array[:,0,0]
	
	for i in range(0,240):
		array[i,:,0] = ncfile.variables['P'][:,:,53,50] + ncfile.variables['PB'][:,:,53,50]
		
	print array[0,:,0]
	
	'''
	plt.figure()
	plt.contourf(Times[0,:],POrinoco[:,0]/100.,NO2Orinoco[:,:])
	plt.ylim(50,1050)
	plt.gca().invert_yaxis()
	plt.xlabel(u'[NO$_2$] [ppb]')
	plt.ylabel('Pressure [hPa]')
	plt.title(r'[NO$_2$] profile')
	plt.legend(loc='upper right')
	plt.grid(True)
	plt.show()
	
	xlist = PBOrinoco[:,0]
	print xlist
	ylist = PBOrinoco[0,:] 
	print ylist
	#X, Y = np.meshgrid(xlist, ylist)
	Z = PBOrinoco[:,:]
	plt.figure()
	cp = plt.contourf(xlist, ylist, Z)
	plt.gca().invert_xaxis()
	plt.colorbar(cp)
	plt.title('Filled Contours Plot')
	plt.xlabel('Time (hours)')
	plt.ylabel('Pressure (hPa)')
	plt.show()
	
	for i in range(0,240,10):
		A=i
		plt.plot(NO2Bogota[i,:],list(range(60)),label='Line 2')
		plt.plot(NO2Orinoco[i,:],list(range(60)))
		plt.plot(NO2Amazone[i,:],list(range(60)))
		plt.legend(['Bogota','Orinoco','Amazone'], loc='best')
		plt.text(14,17,'timestep=')
		plt.text(17,17,A)
		plt.text(17,15,plotvar)
		plt.text(14,15,'plotvar=')
		plt.axis([0, 50, 0, 60])
		plt.show()
	'''
