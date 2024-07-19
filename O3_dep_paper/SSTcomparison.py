from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib import pyplot
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn.metrics import r2_score

print('')
def WRF_SST(wrffile):
	ds = nc.Dataset(wrffile,'r')
	sst_wrf = ds.variables['SST'][:,:,:]
	sst_wrf[sst_wrf < 200] = np.nan
	plt.figure(figsize=(15,10))
	##MAKE BASEMAP
	m = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l')
	#m.drawcoastlines()
	##DRAW COASTLINES, SEA ICE (BOUNDARY), AND CONTINENTS
	m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m.contourf(lons,lats,ds.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002)
	m.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	##DRAW PARALLELS AND MERIDIANS
	m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	##PLOT WRF DATA
	m.imshow(sst_wrf[0,:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
	m.colorbar(label='Sea surface temperature WRF [K]')
	plt.show()
	
def AMSRE_SST(amsrefile):
	es = nc.Dataset(amsrefile,'r')
	sst_amsre = es.variables['Low_res_sst'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
	sst_amsre = np.delete(sst_amsre,np.s_[0:540],axis=0).astype(float)
	sst_amsre[sst_amsre == -9999] = np.nan
	sst_amsre = sst_amsre/100.+273.15
	plt.figure(figsize=(15,10))
	##MAKE BASEMAP
	m = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l')
	sst_amsre = m.transform_scalar(sst_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	#m.drawcoastlines()
	##DRAW COASTLINES, SEA ICE (BOUNDARY), AND CONTINENTS
	m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	##DRAW PARALLELS AND MERIDIANS
	m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	##PLOT WRF DATA
	m.imshow(sst_amsre[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
	m.colorbar(label='Sea surface temperature AMSR-E [K]')
	plt.show()
		
def meandailysst(wrffile,amsremap):
	ds = nc.Dataset(wrffile,'r')
	for day in range(0,int(np.trunc(ds.variables['SST'][:,:,:].shape[0]/24.)),1):
		if day <= 21:
			amsrefile = amsremap+'/AMSR_E_L3_DailyOcean_V04_200808'+str(10+day)+'.nc'
		else:
			amsrefile = amsremap+'/AMSR_E_L3_DailyOcean_V04_2008090'+str(day-21)+'.nc'
		es = nc.Dataset(amsrefile,'r')
		sst_wrf = ds.variables['SST'][(day*24):(day*24)+23,:,:]
		sst_wrf[sst_wrf < 200] = np.nan
		sst_wrf = np.nanmean(sst_wrf,axis=0)	#get mean of day
		seaice_wrf = np.mean(ds.variables['SEAICE'][(day*24):(day*24)+23,:,:],axis=0) #get mean of day
		sst_wrf[seaice_wrf >= 0.5] = np.nan
		sst_amsre = es.variables['Low_res_sst'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
		sst_amsre = np.delete(sst_amsre,np.s_[0:540],axis=0).astype(float)
		sst_amsre[sst_amsre == -9999] = np.nan
		sst_amsre = sst_amsre/100.+273.15
		sst_amsre_vlr = es.variables['Very_low_res_sst'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
		sst_amsre_vlr = np.delete(sst_amsre_vlr,np.s_[0:540],axis=0).astype(float)
		sst_amsre_vlr[sst_amsre_vlr == -9999] = np.nan
		sst_amsre_vlr = sst_amsre_vlr/100.+273.15
		
		fig,axes = plt.subplots(2,3,figsize=(20,10))
		
		m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
		m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m1(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
		m1.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
		m1.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
		m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m1.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		m1.imshow(sst_wrf[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
		
		m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
		sst_amsre_m2 = m2.transform_scalar(sst_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m2.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
		m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m2.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		m2.imshow(sst_amsre_m2[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
		
		m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
		sst_amsre_vlr_m3 = m3.transform_scalar(sst_amsre_vlr,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m3.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
		m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m3.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		m3.imshow(sst_amsre_vlr_m3[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
		
		m4 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,1])
		sst_amsre_m4 = m4.transform_scalar(sst_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m4.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m4(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
		m4.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
		m4.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
		m4.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m4.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		m4.imshow(sst_wrf[:,:]-sst_amsre_m4[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
		
		m5 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,2])
		sst_amsre_vlr_m5 = m5.transform_scalar(sst_amsre_vlr,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m5.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m5(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
		m5.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
		m5.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
		m5.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m5.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		m5.imshow(sst_wrf[:,:]-sst_amsre_vlr_m5[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
		
		#axes[1,2].hist([(sst_wrf[:,:]-sst_amsre_m4[:,:])[~np.isnan(sst_wrf[:,:]-sst_amsre_m4[:,:])],(sst_wrf[:,:]-sst_amsre_vlr_m5[:,:])[~np.isnan(sst_wrf[:,:]-sst_amsre_vlr_m5[:,:])]],bins=50,label=['WRF minus AMSR-E low res.','WRF minus AMSR-E very low res.'])
		#axes[1,2].legend(loc='upper right')
		#axes[1,2].set_xlabel('Sea surface temperature difference [K]')
		#axes[1,2].set_ylabel('Frequency')
		
		ax=axes[1,2]
		pyplot.hist((sst_wrf[:,:]-sst_amsre_m4[:,:])[~np.isnan(sst_wrf[:,:]-sst_amsre_m4[:,:])],bins=50,alpha=0.5,label='WRF minus AMSR-E low res.')
		pyplot.hist((sst_wrf[:,:]-sst_amsre_vlr_m5[:,:])[~np.isnan(sst_wrf[:,:]-sst_amsre_vlr_m5[:,:])],bins=50,alpha=0.5,label='WRF minus AMSR-E very low res.')
		pyplot.legend(loc='upper right')
		pyplot.xlabel('Sea surface temperature difference [K]')
		pyplot.ylabel('Frequency')
		
		plt.savefig('Figures/SSTplots/Junk/SSTdiff_metonly'+str(day))
		#plt.show()
		plt.close()

		diffm4 = (sst_wrf[:,:]-sst_amsre_m4[:,:])[~np.isnan(sst_wrf[:,:]-sst_amsre_m4[:,:])]
		diffm5 = (sst_wrf[:,:]-sst_amsre_vlr_m5[:,:])[~np.isnan(sst_wrf[:,:]-sst_amsre_vlr_m5[:,:])]
		lensim = int(np.trunc(ds.variables['SST'][:,:,:].shape[0]/24.))
		if day == 0:
			mean4 = np.zeros(lensim)
			mean5 = np.zeros(lensim)
			mae4 = np.zeros(lensim)
			mae5 = np.zeros(lensim)
			dayarr = np.zeros(lensim)
		dayarr[day] = day
		mean4[day] = np.round(np.nanmean(diffm4),2)
		mean5[day] = np.round(np.nanmean(diffm5),2)
		mae4[day] = np.round(np.nanmean(np.abs(diffm4)),2)
		mae5[day] = np.round(np.nanmean(np.abs(diffm5)),2)

 		print('Processing timestep '+str(day)+' out of '+str(lensim))
		if day == lensim-1:
			fig,ax = plt.subplots()
			ax.plot(dayarr,mean4,color='black',marker='.')
			ax.plot(dayarr,mean5,color='black',marker='s')
			ax.set_xlabel('Time since start of simulation [days]')
			ax.set_ylabel("Bias [K]",color="black")
			ax.set_xlim([0,27])
			ax2=ax.twinx()
			ax2.plot(dayarr,mae4,color='red',marker='.')
			ax2.plot(dayarr,mae5,color='red',marker='s')
			ax2.set_ylabel("Mean bias error [K]",color='red')
			plt.savefig('Figures/SSTplots/SSTerrorpropagation.png',dpi=300)
			plt.show()
		
def monthlyaverage(wrffile,amsremap):
	ds = nc.Dataset(wrffile,'r')
	amsrefile = amsremap+'/AMSR_E_L3_MonthlyOcean_V04_200808.nc'
	es = nc.Dataset(amsrefile,'r')
	sst_wrf = ds.variables['SST'][:504,:,:]
	sst_wrf[sst_wrf < 200] = np.nan
	sst_wrf = np.nanmean(sst_wrf,axis=0)	#get mean of month
	seaice_wrf = np.mean(ds.variables['SEAICE'][:504,:,:],axis=0) #get mean of day
	sst_wrf[seaice_wrf >= 0.5] = np.nan
	sst_amsre = es.variables['Low_res_sst'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
	sst_amsre = np.delete(sst_amsre,np.s_[0:540],axis=0).astype(float)
	sst_amsre[sst_amsre == -9999] = np.nan
	sst_amsre = sst_amsre/100.+273.15
	sst_amsre_vlr = es.variables['Very_low_res_sst'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
	sst_amsre_vlr = np.delete(sst_amsre_vlr,np.s_[0:540],axis=0).astype(float)
	sst_amsre_vlr[sst_amsre_vlr == -9999] = np.nan
	sst_amsre_vlr = sst_amsre_vlr/100.+273.15
		
	fig,axes = plt.subplots(2,3,figsize=(20,10))
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m1.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
	m1.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	m1.imshow(sst_wrf[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
		
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
	sst_amsre_m2 = m2.transform_scalar(sst_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	m2.imshow(sst_amsre_m2[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
		
	m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
	sst_amsre_vlr_m3 = m3.transform_scalar(sst_amsre_vlr,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m3.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m3.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	m3.imshow(sst_amsre_vlr_m3[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
		
	m4 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,1])
	sst_amsre_m4 = m4.transform_scalar(sst_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	m4.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m4(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m4.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
	m4.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m4.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m4.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	m4.imshow(sst_wrf[:,:]-sst_amsre_m4[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
		
	m5 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,2])
	sst_amsre_vlr_m5 = m5.transform_scalar(sst_amsre_vlr,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	m5.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m5(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m5.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
	m5.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m5.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m5.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	m5.imshow(sst_wrf[:,:]-sst_amsre_vlr_m5[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
	
	ax=axes[1,2]
	pyplot.hist((sst_wrf[:,:]-sst_amsre_m4[:,:])[~np.isnan(sst_wrf[:,:]-sst_amsre_m4[:,:])],bins=50,alpha=0.5,label='WRF minus AMSR-E low res.')
	pyplot.hist((sst_wrf[:,:]-sst_amsre_vlr_m5[:,:])[~np.isnan(sst_wrf[:,:]-sst_amsre_vlr_m5[:,:])],bins=50,alpha=0.5,label='WRF minus AMSR-E very low res.')
	pyplot.legend(loc='upper right')
	pyplot.xlabel('Sea surface temperature difference [K]')
	pyplot.ylabel('Frequency')
	
	plt.savefig('Figures/SSTplots/Junk/MonthlySSTdiffAug_metonly',dpi=300)
	#plt.show()
	plt.close()
	
	
	#NOW FOR SEPTEMBER
	amsrefile = amsremap+'/AMSR_E_L3_MonthlyOcean_V04_200809.nc'
	es = nc.Dataset(amsrefile,'r')
	sst_wrf = ds.variables['SST'][504:,:,:]
	sst_wrf[sst_wrf < 200] = np.nan
	sst_wrf = np.nanmean(sst_wrf,axis=0)	#get mean of month
	seaice_wrf = np.mean(ds.variables['SEAICE'][504:,:,:],axis=0) #get mean of day
	sst_wrf[seaice_wrf >= 0.5] = np.nan
	sst_amsre = es.variables['Low_res_sst'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
	sst_amsre = np.delete(sst_amsre,np.s_[0:540],axis=0).astype(float)
	sst_amsre[sst_amsre == -9999] = np.nan
	sst_amsre = sst_amsre/100.+273.15
	sst_amsre_vlr = es.variables['Very_low_res_sst'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
	sst_amsre_vlr = np.delete(sst_amsre_vlr,np.s_[0:540],axis=0).astype(float)
	sst_amsre_vlr[sst_amsre_vlr == -9999] = np.nan
	sst_amsre_vlr = sst_amsre_vlr/100.+273.15
		
	fig,axes = plt.subplots(2,3,figsize=(20,10))
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m1.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
	m1.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	m1.imshow(sst_wrf[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
		
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
	sst_amsre_m2 = m2.transform_scalar(sst_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	m2.imshow(sst_amsre_m2[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
		
	m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
	sst_amsre_vlr_m3 = m3.transform_scalar(sst_amsre_vlr,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m3.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m3.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	m3.imshow(sst_amsre_vlr_m3[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
		
	m4 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,1])
	sst_amsre_m4 = m4.transform_scalar(sst_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	m4.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m4(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m4.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
	m4.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m4.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m4.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	m4.imshow(sst_wrf[:,:]-sst_amsre_m4[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
		
	m5 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,2])
	sst_amsre_vlr_m5 = m5.transform_scalar(sst_amsre_vlr,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	m5.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m5(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m5.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
	m5.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m5.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m5.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	m5.imshow(sst_wrf[:,:]-sst_amsre_vlr_m5[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
	
	ax=axes[1,2]
	pyplot.hist((sst_wrf[:,:]-sst_amsre_m4[:,:])[~np.isnan(sst_wrf[:,:]-sst_amsre_m4[:,:])],bins=50,alpha=0.5,label='WRF minus AMSR-E low res.')
	pyplot.hist((sst_wrf[:,:]-sst_amsre_vlr_m5[:,:])[~np.isnan(sst_wrf[:,:]-sst_amsre_vlr_m5[:,:])],bins=50,alpha=0.5,label='WRF minus AMSR-E very low res.')
	pyplot.legend(loc='upper right')
	pyplot.xlabel('Sea surface temperature difference [K]')
	pyplot.ylabel('Frequency')
	
	plt.savefig('Figures/SSTplots/Junk/MonthlySSTdiffSep_metonly',dpi=300)
	#plt.show()
	plt.close()

def monthlyaverage_final_vlr(wrffile,amsremap):
	ds = nc.Dataset(wrffile,'r')
	amsrefile = amsremap+'/AMSR_E_L3_MonthlyOcean_V04_200808.nc'
	es = nc.Dataset(amsrefile,'r')
	sst_wrf = ds.variables['SST'][:504,:,:]
	sst_wrf[sst_wrf < 200] = np.nan
	sst_wrf = np.nanmean(sst_wrf,axis=0)	#get mean of month
	seaice_wrf = np.mean(ds.variables['SEAICE'][:504,:,:],axis=0) #get mean of day
	sst_wrf[seaice_wrf >= 0.5] = np.nan
	sst_amsre = es.variables['Very_low_res_sst'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
	sst_amsre = np.delete(sst_amsre,np.s_[0:540],axis=0).astype(float)
	sst_amsre[sst_amsre == -9999] = np.nan
	sst_amsre = sst_amsre/100.+273.15
		
	fig,axes = plt.subplots(2,2,figsize=(17,15), gridspec_kw = {'wspace':0.24, 'hspace':0.02})
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m1.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
	m1.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im = m1.imshow(sst_wrf[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)

	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
	sst_amsre_trans = m2.transform_scalar(sst_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m2.imshow(sst_amsre_trans[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
	
	m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
	m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m3(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m3.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
	m3.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
	m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im3 = m3.imshow(sst_wrf[:,:]-sst_amsre_trans[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
				
	ax=axes[1,1]
	x = sst_wrf[:,:][~np.isnan(sst_wrf[:,:]-sst_amsre_trans[:,:])]
	y = sst_amsre_trans[:,:][~np.isnan(sst_wrf[:,:]-sst_amsre_trans[:,:])]
	r2wrf,r2amsre=x,y
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)
	idx = z.argsort()
	x, y, z = x[idx], y[idx], z[idx]
	ax.scatter(x, y, c=z, cmap='viridis' ,s=5, edgecolor='')
	ax.plot([0,1000],[0,1000],'r-',label='1:1 line')
	ax.yaxis.set_label_position("right")
	ax.yaxis.tick_right()
	ax.set_xlim([272,294])
	ax.set_ylim([272,294])
	ax.set_xlabel('WRF SST [K]',fontsize=13)
	ax.set_ylabel('AMSR-E SST [K]',fontsize=13)
	ax.text(0.01,0.95,'R$^{2}$: '+str(np.round(r2_score(r2amsre,r2wrf),2)),fontsize=12,transform = ax.transAxes)
	
	cax = fig.add_axes([axes[0,0].get_position().x1+0.01,axes[0,0].get_position().y0+0.01,0.02,axes[0,0].get_position().y1-axes[0,0].get_position().y0-0.02])
	cax2 = fig.add_axes([axes[1,0].get_position().x1+0.01,axes[1,0].get_position().y0+0.01,0.02,axes[1,0].get_position().y1-axes[1,0].get_position().y0-0.02])
	cbar = fig.colorbar(im,cax=cax,ticks=[271,275,279,283,287,291])
	cbar.set_label("SST [K]",fontsize=13)
	cbar2 = fig.colorbar(im3,cax=cax2,extend='both',ticks=[-5,-4,-3,-2,-1,0,1,2,3,4,5])
	cbar2.set_label("SST difference (WRF minus AMSR-E) [K]",fontsize=13)
	plt.savefig('Figures/SSTplots/Final_MonthlySSTdiffAug_metonly',dpi=300)
	plt.show()
	#plt.close()
	
	
	#NOW FOR SEPTEMBER
	amsrefile = amsremap+'/AMSR_E_L3_MonthlyOcean_V04_200809.nc'
	es = nc.Dataset(amsrefile,'r')
	sst_wrf = ds.variables['SST'][504:,:,:]
	sst_wrf[sst_wrf < 200] = np.nan
	sst_wrf = np.nanmean(sst_wrf,axis=0)	#get mean of month
	seaice_wrf = np.mean(ds.variables['SEAICE'][504:,:,:],axis=0) #get mean of day
	sst_wrf[seaice_wrf >= 0.5] = np.nan
	sst_amsre = es.variables['Very_low_res_sst'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
	sst_amsre = np.delete(sst_amsre,np.s_[0:540],axis=0).astype(float)
	sst_amsre[sst_amsre == -9999] = np.nan
	sst_amsre = sst_amsre/100.+273.15
		
	fig,axes = plt.subplots(2,2,figsize=(17,15), gridspec_kw = {'wspace':0.24, 'hspace':0.05})
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m1.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
	m1.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im = m1.imshow(sst_wrf[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)

	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
	sst_amsre_trans = m2.transform_scalar(sst_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m2.imshow(sst_amsre_trans[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
	
	m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
	m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m3(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m3.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
	m3.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
	m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im3 = m3.imshow(sst_wrf[:,:]-sst_amsre_trans[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
				
	ax=axes[1,1]
	x = sst_wrf[:,:][~np.isnan(sst_wrf[:,:]-sst_amsre_trans[:,:])]
	y = sst_amsre_trans[:,:][~np.isnan(sst_wrf[:,:]-sst_amsre_trans[:,:])]
	r2wrf,r2amsre=x,y
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)
	idx = z.argsort()
	x, y, z = x[idx], y[idx], z[idx]
	ax.scatter(x, y, c=z, cmap='viridis' ,s=5, edgecolor='')
	ax.plot([0,1000],[0,1000],'r-',label='1:1 line')
	ax.yaxis.set_label_position("right")
	ax.yaxis.tick_right()
	ax.set_xlim([272,294])
	ax.set_ylim([272,294])
	ax.set_xlabel('WRF SST [K]',fontsize=13)
	ax.set_ylabel('AMSR-E SST [K]',fontsize=13)
	ax.text(0.01,0.95,'R$^{2}$: '+str(np.round(r2_score(r2amsre,r2wrf),2)),fontsize=12,transform = ax.transAxes)
	
	cax = fig.add_axes([axes[0,0].get_position().x1+0.01,axes[0,0].get_position().y0+0.01,0.02,axes[0,0].get_position().y1-axes[0,0].get_position().y0-0.02])
	cax2 = fig.add_axes([axes[1,0].get_position().x1+0.01,axes[1,0].get_position().y0+0.01,0.02,axes[1,0].get_position().y1-axes[1,0].get_position().y0-0.02])
	cbar = fig.colorbar(im,cax=cax,ticks=[271,275,279,283,287,291])
	cbar.set_label("SST [K]",fontsize=13)
	cbar2 = fig.colorbar(im3,cax=cax2,extend='both',ticks=[-5,-4,-3,-2,-1,0,1,2,3,4,5])
	cbar2.set_label("SST difference (WRF minus AMSR-E) [K]",fontsize=13)
	plt.savefig('Figures/SSTplots/Final_MonthlySSTdiffSep',dpi=300)
	plt.show()
	#plt.close()

def meandailysst_final(wrffile,amsremap):
	ds = nc.Dataset(wrffile,'r')
	for day in range(0,int(np.trunc(ds.variables['SST'][:,:,:].shape[0]/24.)),1):
		if day <= 21:
			amsrefile = amsremap+'/AMSR_E_L3_DailyOcean_V04_200808'+str(10+day)+'.nc'
		else:
			amsrefile = amsremap+'/AMSR_E_L3_DailyOcean_V04_2008090'+str(day-21)+'.nc'
		es = nc.Dataset(amsrefile,'r')
		sst_wrf = ds.variables['SST'][(day*24):(day*24)+23,:,:]
		sst_wrf[sst_wrf < 200] = np.nan
		sst_wrf = np.nanmean(sst_wrf,axis=0)	#get mean of day
		seaice_wrf = np.mean(ds.variables['SEAICE'][(day*24):(day*24)+23,:,:],axis=0) #get mean of day
		sst_wrf[seaice_wrf >= 0.5] = np.nan
		sst_amsre = es.variables['Very_low_res_sst'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
		sst_amsre = np.delete(sst_amsre,np.s_[0:540],axis=0).astype(float)
		sst_amsre[sst_amsre == -9999] = np.nan
		sst_amsre = sst_amsre/100.+273.15
		
		fig,axes = plt.subplots(2,2,figsize=(17,15), gridspec_kw = {'wspace':0.24, 'hspace':0.02})
		
		m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
		m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m1(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
		m1.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
		m1.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		im = m1.imshow(sst_wrf[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)

		m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
		sst_amsre_trans = m2.transform_scalar(sst_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m2.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		m2.imshow(sst_amsre_trans[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=271,vmax=293)
	
		m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
		m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m3(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
		m3.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
		m3.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		im3 = m3.imshow(sst_wrf[:,:]-sst_amsre_trans[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
				
		ax=axes[1,1]
		x = sst_wrf[:,:][~np.isnan(sst_wrf[:,:]-sst_amsre_trans[:,:])]
		y = sst_amsre_trans[:,:][~np.isnan(sst_wrf[:,:]-sst_amsre_trans[:,:])]
		r2wrf,r2amsre=x,y
		xy = np.vstack([x,y])
		z = gaussian_kde(xy)(xy)
		idx = z.argsort()
		x, y, z = x[idx], y[idx], z[idx]
		ax.scatter(x, y, c=z, cmap='viridis' ,s=5, edgecolor='')
		ax.plot([0,1000],[0,1000],'r-',label='1:1 line')
		ax.yaxis.set_label_position("right")
		ax.yaxis.tick_right()
		ax.set_xlim([272,294])
		ax.set_ylim([272,294])
		ax.set_xlabel('WRF SST [K]',fontsize=13)
		ax.set_ylabel('AMSR-E SST [K]',fontsize=13)
		ax.text(0.01,0.95,'R$^{2}$: '+str(np.round(r2_score(r2amsre,r2wrf),2)),fontsize=12,transform = ax.transAxes)
	
		cax = fig.add_axes([axes[0,0].get_position().x1+0.01,axes[0,0].get_position().y0+0.01,0.02,axes[0,0].get_position().y1-axes[0,0].get_position().y0-0.02])
		cax2 = fig.add_axes([axes[1,0].get_position().x1+0.01,axes[1,0].get_position().y0+0.01,0.02,axes[1,0].get_position().y1-axes[1,0].get_position().y0-0.02])
		cbar = fig.colorbar(im,cax=cax,ticks=[271,275,279,283,287,291])
		cbar.set_label("SST [K]",fontsize=13)
		cbar2 = fig.colorbar(im3,cax=cax2,extend='both',ticks=[-5,-4,-3,-2,-1,0,1,2,3,4,5])
		cbar2.set_label("SST difference (WRF minus AMSR-E) [K]",fontsize=13)
		plt.savefig('Figures/SSTplots/Final_SSTdiff'+str(day),dpi=300)
		#plt.show()
		plt.close()
				
				
#WRF_SST('/home/WUR/barte035/WRFChem/WRFV411-Polar/run/wrfout_metonlyv3_d01_2008-08-10_00:00:00')
#AMSRE_SST('/home/WUR/barte035/WRFChem/AMSR-E/AMSR_E_L3_DailyOcean_V04_20080810.nc')
#meandailysst('/lustre/nobackup/WUR/ESG/barte035/wrfout_metonlyv3_d01_2008-08-10_00:00:00','/home/WUR/barte035/WRFChem/AMSR-E/')
#monthlyaverage('/lustre/nobackup/WUR/ESG/barte035/wrfout_metonlyv3_d01_2008-08-10_00:00:00','/home/WUR/barte035/WRFChem/AMSR-E/MonthlyAverage')
monthlyaverage_final_vlr('/lustre/nobackup/WUR/ESG/barte035/wrfout_metonlyv3_d01_2008-08-10_00:00:00','/home/WUR/barte035/WRFChem/AMSR-E/MonthlyAverage')
meandailysst_final('/lustre/nobackup/WUR/ESG/barte035/wrfout_metonlyv3_d01_2008-08-10_00:00:00','/home/WUR/barte035/WRFChem/AMSR-E/')
