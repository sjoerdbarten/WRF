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


print 'This script makes a contourplot of vertical profiles over the whole output timeseries (or selected times). To use this use verticalprofilev2(plotvar,location,area,figsave) where plotvar is any WRF-Chem variable (or NOx), location = Bogota, Orinoco or Amazone and area = (0 (1 pixel),1 (area averaged)). Figsave = 1 saves figure.'

def verticalprofilev2(plotvar,location,area,figsave):
	ncfile = nc.Dataset('/home/WUR/barte035/WRFChem/WRFV3/run/output/wrfout_d01_2014-01-01_00:00:00','r')
	spinup = 24
	timesim = 31*24-6
	nlev = 60
	
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

	if area == 0:
		Times = np.zeros((timesim-spinup,nlev))
		for i in range(0,int(nlev),1):
			Times[:,i] = np.arange(spinup,timesim,1)/24.
		pressure = np.zeros((timesim-spinup+1,nlev))
		concentration = np.zeros((timesim-spinup+1,nlev))
		cldfra = np.zeros((timesim-spinup+1,nlev))
		
		if location == 'Bogota':
			yloc = 48
			xloc = 34
		if location == 'Orinoco':
			yloc = 53
			xloc = 50
		if location == 'Amazone':
			yloc = 10
			xloc = 85
	
		for i in range(int(spinup),int(timesim),1):
			print 'Processing time',i
			for j in range(0,nlev):
				pressure[i,j] = ncfile.variables['P'][i,j,yloc,xloc]/100. + ncfile.variables['PB'][i,j,yloc,xloc]/100.
				cldfra[i,j] = ncfile.variables['CLDFRA'][i,j,yloc,xloc]*100.
				if plotvar == 'nox':
					concentration[i,j] = ncfile.variables['no'][i,j,yloc,xloc]*1000.+ncfile.variables['no2'][i,j,yloc,xloc]*1000.
				if plotvar != 'nox':
					concentration[i,j] = ncfile.variables[plotvar][i,j,yloc,xloc]*1000.
	
	if area == 1:
		Times = np.zeros((timesim-spinup,nlev))
		for i in range(0,int(nlev),1):
			Times[:,i] = np.arange(spinup,timesim,1)/24.
		pressure = np.zeros((timesim-spinup,nlev))
		pressure1 = np.zeros((timesim-spinup,nlev,4))
		cldfra = np.zeros((timesim-spinup,nlev))
		cldfra1 = np.zeros((timesim-spinup,nlev,4))
		concentration = np.zeros((timesim-spinup,nlev))
		concentration1 = np.zeros((timesim-spinup,nlev,4))
		temp = np.zeros((timesim-spinup,nlev))
		
		
		if location == 'Bogota':
			for i in range(spinup,int(timesim),1):
				print 'Processing time',i
				for j in range(0,nlev):
					pressure1[i-spinup,j,0] = ncfile.variables['P'][i,j,48,34]/100.+ncfile.variables['PB'][i,j,48,34]/100.
					pressure1[i-spinup,j,1] = ncfile.variables['P'][i,j,47,34]/100.+ncfile.variables['PB'][i,j,47,34]/100.
					pressure1[i-spinup,j,2] = ncfile.variables['P'][i,j,47,33]/100.+ncfile.variables['PB'][i,j,47,33]/100.
					pressure1[i-spinup,j,3] = ncfile.variables['P'][i,j,48,33]/100.+ncfile.variables['PB'][i,j,48,33]/100.
					#temp[i-spinup,j] =  ncfile.variables['T'][i,j,48,34]
					#cldfra1[i,j,0] = ncfile.variables['CLDFRA'][i,j,48,34]*100.
					#cldfra1[i,j,1] = ncfile.variables['CLDFRA'][i,j,47,34]*100.
					#cldfra1[i,j,2] = ncfile.variables['CLDFRA'][i,j,47,33]*100.
					#cldfra1[i,j,3] = ncfile.variables['CLDFRA'][i,j,48,33]*100.
					if plotvar == 'nox':
						concentration1[i-spinup,j,0] = ncfile.variables['no'][i,j,48,34]*1000.+ncfile.variables['no2'][i,j,48,34]*1000.
						concentration1[i-spinup,j,1] = ncfile.variables['no'][i,j,47,34]*1000.+ncfile.variables['no2'][i,j,47,34]*1000.
						concentration1[i-spinup,j,2] = ncfile.variables['no'][i,j,47,33]*1000.+ncfile.variables['no2'][i,j,47,33]*1000.
						concentration1[i-spinup,j,3] = ncfile.variables['no'][i,j,48,33]*1000.+ncfile.variables['no2'][i,j,48,33]*1000.
					if plotvar != 'nox':
						concentration1[i-spinup,j,0] = ncfile.variables[plotvar][i,j,48,34]*1000.
						concentration1[i-spinup,j,1] = ncfile.variables[plotvar][i,j,47,34]*1000.
						concentration1[i-spinup,j,2] = ncfile.variables[plotvar][i,j,47,33]*1000.
						concentration1[i-spinup,j,3] = ncfile.variables[plotvar][i,j,48,33]*1000.
				
		if location == 'Orinoco':
			for i in range(spinup,int(timesim),1):
				print 'Processing time',i
				for j in range(0,nlev):
					pressure1[i-spinup,j,0] = ncfile.variables['P'][i,j,51,46]/100.+ncfile.variables['PB'][i,j,51,46]/100.
					pressure1[i-spinup,j,1] = ncfile.variables['P'][i,j,51,47]/100.+ncfile.variables['PB'][i,j,51,47]/100.
					pressure1[i-spinup,j,2] = ncfile.variables['P'][i,j,52,46]/100.+ncfile.variables['PB'][i,j,52,46]/100.
					pressure1[i-spinup,j,3] = ncfile.variables['P'][i,j,52,47]/100.+ncfile.variables['PB'][i,j,52,47]/100.
					#temp[i-spinup,j] =  ncfile.variables['T'][i,j,51,46]
					#cldfra1[i,j,0] = ncfile.variables['CLDFRA'][i,j,51,46]*100.
					#cldfra1[i,j,1] = ncfile.variables['CLDFRA'][i,j,51,47]*100.
					#cldfra1[i,j,2] = ncfile.variables['CLDFRA'][i,j,52,46]*100.
					#cldfra1[i,j,3] = ncfile.variables['CLDFRA'][i,j,52,47]*100.
					if plotvar == 'nox':
						concentration1[i-spinup,j,0] = ncfile.variables['no'][i,j,51,46]*1000.+ncfile.variables['no2'][i,j,51,46]*1000.
						concentration1[i-spinup,j,1] = ncfile.variables['no'][i,j,51,47]*1000.+ncfile.variables['no2'][i,j,51,47]*1000.
						concentration1[i-spinup,j,2] = ncfile.variables['no'][i,j,47,33]*1000.+ncfile.variables['no2'][i,j,52,46]*1000.
						concentration1[i-spinup,j,3] = ncfile.variables['no'][i,j,48,33]*1000.+ncfile.variables['no2'][i,j,52,47]*1000.
					if plotvar != 'nox':
						concentration1[i-spinup,j,0] = ncfile.variables[plotvar][i,j,51,46]*1000.
						concentration1[i-spinup,j,1] = ncfile.variables[plotvar][i,j,51,47]*1000.
						concentration1[i-spinup,j,2] = ncfile.variables[plotvar][i,j,52,46]*1000.
						concentration1[i-spinup,j,3] = ncfile.variables[plotvar][i,j,52,47]*1000.
				
		if location == 'Amazone':
			for i in range(spinup,int(timesim),1):
				print 'Processing time',i
				for j in range(0,nlev):
					pressure1[i-spinup,j,0] = ncfile.variables['P'][i,j,11,89]/100.+ncfile.variables['PB'][i,j,11,89]/100.
					pressure1[i-spinup,j,1] = ncfile.variables['P'][i,j,11,90]/100.+ncfile.variables['PB'][i,j,11,90]/100.
					pressure1[i-spinup,j,2] = ncfile.variables['P'][i,j,12,89]/100.+ncfile.variables['PB'][i,j,12,89]/100.
					pressure1[i-spinup,j,3] = ncfile.variables['P'][i,j,12,90]/100.+ncfile.variables['PB'][i,j,12,90]/100.
					#temp[i-spinup,j] =  ncfile.variables['T'][i,j,11,89]
					#cldfra1[i,j,0] = ncfile.variables['CLDFRA'][i,j,11,89]*100.
					#cldfra1[i,j,1] = ncfile.variables['CLDFRA'][i,j,11,90]*100.
					#cldfra1[i,j,2] = ncfile.variables['CLDFRA'][i,j,12,89]*100.
					#cldfra1[i,j,3] = ncfile.variables['CLDFRA'][i,j,12,90]*100.
					if plotvar == 'nox':
						concentration1[i-spinup,j,0] = ncfile.variables['no'][i,j,11,89]*1000.+ncfile.variables['no2'][i,j,11,89]*1000.
						concentration1[i-spinup,j,1] = ncfile.variables['no'][i,j,11,90]*1000.+ncfile.variables['no2'][i,j,11,90]*1000.
						concentration1[i-spinup,j,2] = ncfile.variables['no'][i,j,12,89]*1000.+ncfile.variables['no2'][i,j,12,89]*1000.
						concentration1[i-spinup,j,3] = ncfile.variables['no'][i,j,12,90]*1000.+ncfile.variables['no2'][i,j,12,90]*1000.
					if plotvar != 'nox':
						concentration1[i-spinup,j,0] = ncfile.variables[plotvar][i,j,11,89]*1000.
						concentration1[i-spinup,j,1] = ncfile.variables[plotvar][i,j,11,90]*1000.
						concentration1[i-spinup,j,2] = ncfile.variables[plotvar][i,j,12,89]*1000.
						concentration1[i-spinup,j,3] = ncfile.variables[plotvar][i,j,12,90]*1000.
			
		pressure = np.mean(pressure1,axis=2)
		#cldfra = np.mean(cldfra1,axis=2)
		concentration = np.mean(concentration1,axis=2)
	

	#This is a quick fix for colorbar labels and will not do well for any case. Sorry.
	lvls = np.logspace(-1.5,1.5,25)
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
	if location == 'Bogota':
		plt.text(2,60,'(a)',fontsize=9)
	if location == 'Orinoco':
		plt.text(2,60,'(b)',fontsize=9)
	if location == 'Amazone':
		plt.text(2,60,'(c)',fontsize=9)
	if figsave == 1:
		plt.savefig('/home/WUR/barte035/WRFChem/Python/Figures/Vertical'+str(location)+''+str(plotvar)+'.png',bbox_inches='tight')
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
