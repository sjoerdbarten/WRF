from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

print('This script makes a timeseries of maps of wrf output (surface variable) using maptimeseries(wrffile,plotvar,onlyocean,savefig). wrffile = path to wrffile; plotvar = wrf variable to plot; onlyocean = [0,1] (0 is plot over entire domain, 1 is plot over ocean only); savefig = [0,1] (0 to not save, 1 to save)')
def maptimeseries(wrffile,plotvar,onlyocean,savefig):
	#DEFINE VARIABLE AND WRF OUTPUT
	ds = nc.Dataset(wrffile,'r')
	if plotvar in ['DEPO3COAREG','DEPCO2COAREG','DEPDMSCOAREG','DEPCH4COAREG']:
		multfactor = 1.e12
		labelunit = ' [10$^{-12}$ kg m$^{-2}$ s$^{-1}$]'
	if plotvar in ['co2']:
		multfactor = 1.
		labelunit = ' [ppm]'
	if plotvar in  ['o3','ch4','so2']:
		multfactor = 1000.
		labelunit = ' [ppb]'
	if plotvar in ['dms']:
		multfactor = 1000000.
		labelunit = ' [ppt]'
	if plotvar in ['ho']:
		multfactor = 1000000000.
		labelunit = ' [10$^{-3}$ ppt]'
	if plotvar in ['VTC','VTCCO2','VTCDMS','VTCCH4','DEP_VEL']:
		multfactor = 100000.
		labelunit = ' [10$^{-3}$ cm s$^{-1}$]'
	print(ds.variables[plotvar].ndim)
	if ds.variables[plotvar].ndim == 3:
		plotvar2 = ds.variables[plotvar][:,:,:]*multfactor
	if ds.variables[plotvar].ndim == 4:
		plotvar2 = ds.variables[plotvar][:,0,:,:]*multfactor	#Get bot level
	minplotvar = np.round(np.min((plotvar2)),1)
	maxplotvar = np.round(np.max((plotvar2)),1)
	#if plotvar in ['co2']:
	#	minplotvar = 360.
	#	maxplotvar = 381.
	#if plotvar in ['ch4']:
	#	minplotvar = 1600.
	#	maxplotvar = 1710.
	if plotvar in ['so2']:
		minplotvar = 0.
		maxplotvar = 2.
	if plotvar in ['ho']:
		minplotvar = 0.
		maxplotvar = 10.


	print(minplotvar)
	print(maxplotvar)
	diff = maxplotvar-minplotvar
	ticks = [minplotvar,minplotvar+diff*0.1,minplotvar+diff*0.2,minplotvar+diff*0.3,minplotvar+diff*0.4,minplotvar+diff*0.5,minplotvar+diff*0.6,minplotvar+diff*0.7,minplotvar+diff*0.8,minplotvar+diff*0.9,maxplotvar]
	#if (maxplotvar+minplotvar)/2 > 0:
	#	ticks = [minplotvar,maxplotvar*0.1,maxplotvar*0.2,maxplotvar*0.3,maxplotvar*0.4,maxplotvar*0.5,maxplotvar*0.6,maxplotvar*0.7,maxplotvar*0.8,maxplotvar*0.9,maxplotvar]
	#if (minplotvar+maxplotvar)/2 < 0:
	#	ticks = [minplotvar,minplotvar*0.9,minplotvar*0.8,minplotvar*0.7,minplotvar*0.6,minplotvar*0.5,minplotvar*0.4,minplotvar*0.3,minplotvar*0.2,minplotvar*0.1,maxplotvar]
	
	for i in range(0,plotvar2.shape[0],1):
		print('Doing timestep '+str(i))
		##MAKE FIGURE
		plt.figure(figsize=(15,10))
		##MAKE BASEMAP
		m = Basemap(projection='npstere',boundinglat=65.8,lon_0=0.,resolution='l')
		#m.drawcoastlines()
		##DRAW COASTLINES, SEA ICE (BOUNDARY), AND CONTINENTS
		m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m(ds.variables['XLONG'][i,:,:],ds.variables['XLAT'][i,:,:],inverse=False)
		clevs = [0,1]
		if onlyocean == 0:
			m.contour(lons,lats,ds.variables['SEAICE'][i,:,:],clevs,linewidths=1.,colors='k',linestyles='dashed',zorder=1002)	#For sea ice boundary (black dashed line)
		if onlyocean == 1:
			m.contourf(lons,lats,ds.variables['SEAICE'][i,:,:],[0.5,1.5],colors='w',zorder=1002)					#For filled sea ice (white)
			m.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)							#For filled continents
		##DRAW PARALLELS AND MERIDIANS
		m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		#m.drawmapboundary(fill_color='aqua')
		##PLOT WRF DATA
		m.imshow(plotvar2[i,:,:],cmap='jet',interpolation='none',zorder=1000,vmin=minplotvar,vmax=maxplotvar)
		m.colorbar(label=plotvar+labelunit,ticks=ticks)
		if savefig == 1:
			if i < 10:
				plt.savefig('./Figures/Maptimeseries/'+plotvar+'/'+plotvar+'time00'+str(i)+'.png',dpi=300)
			if i < 100 and i >= 10:
				plt.savefig('./Figures/Maptimeseries/'+plotvar+'/'+plotvar+'time0'+str(i)+'.png',dpi=300)
			if i < 1000 and i >= 100:
				plt.savefig('./Figures/Maptimeseries/'+plotvar+'/'+plotvar+'time'+str(i)+'.png',dpi=300)
		plt.close()				
		#plt.show()

def maptimeseriesdiff(wrffile1,wrffile2,plotvar,onlyocean,savefig):
	#DEFINE VARIABLE AND WRF OUTPUT
	ds1 = nc.Dataset(wrffile1,'r')
	ds2 = nc.Dataset(wrffile2,'r')
	if plotvar in ['DEPO3COAREG','DEPCO2COAREG','DEPDMSCOAREG','DEPCH4COAREG']:
		multfactor = 1.e12
		labelunit = ' [10$^{-12}$ kg m$^{-2}$ s$^{-1}$]'
	if plotvar in ['co2']:
		multfactor = 1.
		labelunit = ' [ppm]'
	if plotvar in  ['o3','ch4','so2']:
		multfactor = 1000.
		labelunit = ' [ppb]'
	if plotvar in ['dms']:
		multfactor = 1000000.
		labelunit = ' [ppt]'
	if plotvar in ['ho']:
		multfactor = 1000000000.
		labelunit = ' [10$^{-3}$ ppt]'
	if plotvar in ['VTC','VTCCO2','VTCDMS','VTCCH4','DEP_VEL']:
		multfactor = 100000.
		labelunit = ' [10$^{-3}$ cm s$^{-1}$]'
	print(ds1.variables[plotvar].ndim)
	if ds1.variables[plotvar].ndim == 3:
		plotvar2 = (ds1.variables[plotvar][:,:,:]-ds2.variables[plotvar][:,:,:])*multfactor
	if ds1.variables[plotvar].ndim == 4:
		plotvar2 = (ds1.variables[plotvar][:,0,:,:]-ds2.variables[plotvar][:,0,:,:])*multfactor	#Get bot level
	minplotvar = np.round(np.min((plotvar2)),1)
	maxplotvar = np.round(np.max((plotvar2)),1)
	if plotvar in ['o3']:
		minplotvar = -8.
		maxplotvar = 8.
	
	#Make sure 0 is the center of data.	
	absmax = max(abs(minplotvar),abs(maxplotvar))	
	minplotvar = -absmax
	maxplotvar = absmax

	print(minplotvar)
	print(maxplotvar)
	diff = maxplotvar-minplotvar
	ticks = [minplotvar,minplotvar+diff*0.1,minplotvar+diff*0.2,minplotvar+diff*0.3,minplotvar+diff*0.4,minplotvar+diff*0.5,minplotvar+diff*0.6,minplotvar+diff*0.7,minplotvar+diff*0.8,minplotvar+diff*0.9,maxplotvar]
	#if (maxplotvar+minplotvar)/2 > 0:
	#	ticks = [minplotvar,maxplotvar*0.1,maxplotvar*0.2,maxplotvar*0.3,maxplotvar*0.4,maxplotvar*0.5,maxplotvar*0.6,maxplotvar*0.7,maxplotvar*0.8,maxplotvar*0.9,maxplotvar]
	#if (minplotvar+maxplotvar)/2 < 0:
	#	ticks = [minplotvar,minplotvar*0.9,minplotvar*0.8,minplotvar*0.7,minplotvar*0.6,minplotvar*0.5,minplotvar*0.4,minplotvar*0.3,minplotvar*0.2,minplotvar*0.1,maxplotvar]
	
	for i in range(0,plotvar2.shape[0],1):
		print('Doing timestep '+str(i))
		##MAKE FIGURE
		plt.figure(figsize=(15,10))
		##MAKE BASEMAP
		m = Basemap(projection='npstere',boundinglat=65.8,lon_0=0.,resolution='l')
		#m.drawcoastlines()
		##DRAW COASTLINES, SEA ICE (BOUNDARY), AND CONTINENTS
		m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m(ds1.variables['XLONG'][i,:,:],ds1.variables['XLAT'][i,:,:],inverse=False)
		clevs = [0,1]
		if onlyocean == 0:
			m.contour(lons,lats,ds1.variables['SEAICE'][i,:,:],clevs,linewidths=1.,colors='k',linestyles='dashed',zorder=1002)	#For sea ice boundary (black dashed line)
		if onlyocean == 1:
			m.contourf(lons,lats,ds1.variables['SEAICE'][i,:,:],[0.5,1.5],colors='w',zorder=1002)					#For filled sea ice (white)
			m.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)							#For filled continents
		##DRAW PARALLELS AND MERIDIANS
		m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		#m.drawmapboundary(fill_color='aqua')
		##PLOT WRF DATA
		m.imshow(plotvar2[i,:,:],cmap='RdBu',interpolation='none',zorder=1000,vmin=minplotvar,vmax=maxplotvar)
		m.colorbar(label=plotvar+labelunit,ticks=ticks)
		if savefig == 1:
			if i < 10:
				plt.savefig('./Figures/Maptimeseries/diff'+plotvar+'/diff'+plotvar+'time00'+str(i)+'.png',dpi=300)
			if i < 100 and i >= 10:
				plt.savefig('./Figures/Maptimeseries/diff'+plotvar+'/diff'+plotvar+'time0'+str(i)+'.png',dpi=300)
			if i < 1000 and i >= 100:
				plt.savefig('./Figures/Maptimeseries/diff'+plotvar+'/diff'+plotvar+'time'+str(i)+'.png',dpi=300)
		plt.close()				
		#plt.show()
