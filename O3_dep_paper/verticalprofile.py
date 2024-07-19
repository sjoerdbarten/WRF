import netCDF4 as nc
import scipy.interpolate
from pylab import *
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.axes as ax
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm

print 'This script makes a contourplot of vertical profiles over the whole output timeseries (or selected times). To use this use verticalprofile(wrffile,plotvar,locx,locy,spinup,figsave) where plotvar is any WRF-Chem variable (or NOx). Figsave = 1 saves figure.'

def verticalprofile(wrffile,plotvar,locx,locy,spinup,figsave):
	ncfile = nc.Dataset(wrffile,'r')
	print(ncfile)
	print(ncfile[plotvar].shape)
	timesim = ncfile[plotvar].shape[0]
	nlev = ncfile[plotvar].shape[1]
	
	if plotvar == 'o3':
		plotvarlegend = 'O$_3$'
	if plotvar == 'co':
		plotvarlegend = 'CO'
	if plotvar == 'no':
		plotvarlegend = 'NO'
	if plotvar == 'no2':
		plotvarlegend = 'NO$_2$'
	if plotvar == 'nox':
		plotvarlegend = 'NO$_x$'
	if plotvar == 'iso':
		plotvarlegend = 'Isoprene'
	if plotvar == 'lnox_total':
		plotvarlegend = 'LNO$_x$ additional'
	if plotvar == 'CLDFRA':
		plotvarlegend = 'Cloudcover'
	if plotvar == 'ho':
		plotvarlegend = 'OH'
	if plotvar == 'PHOTR2':
		plotvarlegend = 'PHOTR2'
	if plotvar == 'PHOTR4':
		plotvarlegend = 'PHOTR4'
	if plotvar == 'W':
		plotvarlegend = 'Vertical wind speed'

	Times = np.zeros((timesim-spinup,nlev))
	for i in range(0,int(nlev),1):
		Times[:,i] = np.arange(spinup,timesim,1)/24.
	pressure = np.zeros((timesim-spinup,nlev))
	concentration = np.zeros((timesim-spinup,nlev))
	cldfra = np.zeros((timesim-spinup,nlev))
	
	yloc = locy
	xloc = locx
	
	for i in range(int(spinup),int(timesim),1):
		print 'Processing time',i
		for j in range(0,nlev):
			pressure[i,j] = ncfile.variables['P'][i,j,yloc,xloc]/100. + ncfile.variables['PB'][i,j,yloc,xloc]/100.
			cldfra[i,j] = ncfile.variables['CLDFRA'][i,j,yloc,xloc]*100.
			if plotvar == 'nox':
				concentration[i,j] = ncfile.variables['no'][i,j,yloc,xloc]*1000.+ncfile.variables['no2'][i,j,yloc,xloc]*1000.
			if plotvar != 'nox':
				concentration[i,j] = ncfile.variables[plotvar][i,j,yloc,xloc]*1000.
	
	#This is a quick fix for colorbar labels and will not do well for any case. Sorry.
	lvls = np.logspace(1,2,25)
	lvlslabel=['','','','',0.1,'','','','','','','',1.,'','','','','','','',10.,'','','','']
	
	plt.figure()
	ax = plt.axes()
	plt.contourf(Times,pressure,concentration,lvls,norm=LogNorm(),cmap=plt.cm.Oranges,origin='lower')
	#plt.title(''+str(plotvarlegend)+' concentration vertical profile')
	plt.gca().invert_yaxis()
	plt.xlabel('Time since start of simulation [days]')
	plt.ylabel('Pressure [hPa]')
	cbar = plt.colorbar(ticks=lvls)
	cbar.set_ticklabels(lvlslabel)
	plt.xlim([spinup/24, timesim/24])
	plt.xticks(np.arange(spinup/24, timesim/24, 3))
	plt.yscale('log')
	plt.ylim([1000,50])
	plt.yticks([1000,900,800,700,600,500,400,300,200,100,50])
	ax.set_yticks([1000,900,800,700,600,500,400,300,200,100,50])
	ax.set_yticklabels(["1000","900","800","700","600","500","400","300","200","100","50"])
	cbar.ax.set_ylabel(''+str(plotvarlegend)+' mixing ratio [ppb]')
	#plt.contour(Times,pressure,cldfra,levels=[90,91],colors='r')
	if figsave == 1:
		plt.savefig('/home/WUR/barte035/WRFChem/Python/Figures/Verticalprofile_'+str(plotvar)+'.png',bbox_inches='tight')
	plt.show()


	
	#return pressure,concentration,temp #CHANGE THIS AND COMMENT PLOT SECTION ABOVE IF YOU DO LOOP

	#fig = plt.figure()
	#ax = Axes3D(fig)
	#ax.plot_surface(Times,pressure,concentration,rstride=1,cstride=1,cmap='hot')
	#ax.contourf(Times,pressure,concentration, zdir='z', offset=-2, cmap=plt.cm.hot)
	#plt.show()
	
	#plt.figure()
	#plt.scatter(concentration,cldfra,c='r')
	#plt.show()
