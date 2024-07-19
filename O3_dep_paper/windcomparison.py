from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib import pyplot
from sklearn.metrics import mean_squared_error
from math import sqrt
from scipy.stats import gaussian_kde
from sklearn.metrics import r2_score
from matplotlib.dates import DateFormatter
import datetime


def WRF_wind(wrffile):
	ds = nc.Dataset(wrffile,'r')
	wind_wrf = (ds.variables['U10'][:,:,:]**2+ds.variables['V10'][:,:,:]**2)**0.5
	plt.figure(figsize=(15,10))
	##MAKE BASEMAP
	m = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l')
	#m.drawcoastlines()
	##DRAW COASTLINES, SEA ICE (BOUNDARY), AND CONTINENTS
	m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
	m.contour(lons,lats,ds.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	##DRAW PARALLELS AND MERIDIANS
	m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	##PLOT WRF DATA
	m.imshow(wind_wrf[0,:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=20.)
	m.colorbar(label='10m wind WRF [m s$^{-1}$]',extend='max')
	plt.show()
	
def AMSRE_wind(amsrefile):
	es = nc.Dataset(amsrefile,'r')
	wind_amsre = es.variables['Low_res_wind'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
	wind_amsre = np.delete(wind_amsre,np.s_[0:540],axis=0).astype(float)
	wind_amsre[wind_amsre == -9999] = np.nan
	wind_amsre = wind_amsre/100.
	plt.figure(figsize=(15,10))
	##MAKE BASEMAP
	m = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l')
	wind_amsre = m.transform_scalar(wind_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
	#m.drawcoastlines()
	##DRAW COASTLINES, SEA ICE (BOUNDARY), AND CONTINENTS
	m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	##DRAW PARALLELS AND MERIDIANS
	m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
	##PLOT WRF DATA
	m.imshow(wind_amsre[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=20.)
	m.colorbar(label='Surface wind AMSR-E low resolution [m s$^{-1}$]',extend='max')
	plt.show()
		
def meandailywind(wrffile,amsremap):
	ds = nc.Dataset(wrffile,'r')
	ds_nudged = nc.Dataset('/lustre/nobackup/WUR/ESG/barte035/wrfout_metnudgedBL_d01_2008-08-10_00:00:00','r')
	for day in range(0,int(np.trunc(ds.variables['U10'][:,:,:].shape[0]/24.)),1):
		if day <= 21:
			amsrefile = amsremap+'/AMSR_E_L3_DailyOcean_V04_200808'+str(10+day)+'.nc'
		else:
			amsrefile = amsremap+'/AMSR_E_L3_DailyOcean_V04_2008090'+str(day-21)+'.nc'
		es = nc.Dataset(amsrefile,'r')
		wind_wrf = (ds.variables['U10'][(day*24):(day*24)+23,:,:]**2+ds.variables['V10'][(day*24):(day*24)+23,:,:]**2)**0.5
		wind_wrf = np.nanmean(wind_wrf,axis=0)	#get mean of day
		seaice_wrf = np.mean(ds.variables['SEAICE'][(day*24):(day*24)+23,:,:],axis=0) #get mean of day
		wind_wrf[seaice_wrf >= 0.5] = np.nan
		wind_amsre = es.variables['Low_res_wind'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
		wind_amsre = np.delete(wind_amsre,np.s_[0:540],axis=0).astype(float)
		wind_amsre[wind_amsre == -9999] = np.nan
		wind_amsre = wind_amsre/100.
		wind_amsre_mr = es.variables['Med_res_wind'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
		wind_amsre_mr = np.delete(wind_amsre_mr,np.s_[0:540],axis=0).astype(float)
		wind_amsre_mr[wind_amsre_mr == -9999] = np.nan
		wind_amsre_mr = wind_amsre_mr/100.
		wind_wrf_nudged = (ds_nudged.variables['U10'][(day*24):(day*24)+23,:,:]**2+ds_nudged.variables['V10'][(day*24):(day*24)+23,:,:]**2)**0.5
		wind_wrf_nudged = np.nanmean(wind_wrf_nudged,axis=0)	#get mean of day
		wind_wrf_nudged[seaice_wrf >= 0.5] = np.nan
		
		fig,axes = plt.subplots(2,3,figsize=(20,10))
		
		m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
		m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m1(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
		m1.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
		m1.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
		m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m1.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		m1.imshow(wind_wrf[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=20.)
		
		m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
		wind_amsre_m2 = m2.transform_scalar(wind_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m2.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
		m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m2.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		m2.imshow(wind_amsre_m2[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=20.)
		
		m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
		wind_amsre_mr_m3 = m3.transform_scalar(wind_amsre_mr,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m3.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
		m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m3.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		m3.imshow(wind_amsre_mr_m3[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=20.)
		
		m4 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,1])
		wind_amsre_m4 = m4.transform_scalar(wind_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m4.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m4(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
		m4.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
		m4.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
		m4.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m4.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		m4.imshow(wind_wrf[:,:]-wind_amsre_m4[:,:],cmap='RdBu',interpolation='none',zorder=1000,vmin=-3.,vmax=3.)
		
		m5 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,2])
		wind_amsre_mr_m5 = m5.transform_scalar(wind_amsre_mr,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m5.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m5(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
		m5.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
		m5.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
		m5.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m5.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
		m5.imshow(wind_wrf[:,:]-wind_amsre_mr_m5[:,:],cmap='RdBu',interpolation='none',zorder=1000,vmin=-3.,vmax=3.)
				
		ax=axes[1,2]
		diffm4 = (wind_wrf[:,:]-wind_amsre_m4[:,:])[~np.isnan(wind_wrf[:,:]-wind_amsre_m4[:,:])]
		diffm5 = (wind_wrf[:,:]-wind_amsre_mr_m5[:,:])[~np.isnan(wind_wrf[:,:]-wind_amsre_mr_m5[:,:])]
		diffm4_nudged = (wind_wrf_nudged[:,:]-wind_amsre_m4[:,:])[~np.isnan(wind_wrf_nudged[:,:]-wind_amsre_m4[:,:])]
		diffm5_nudged = (wind_wrf_nudged[:,:]-wind_amsre_mr_m5[:,:])[~np.isnan(wind_wrf_nudged[:,:]-wind_amsre_mr_m5[:,:])]
		pyplot.hist(diffm4,bins=50,alpha=0.5,label='WRF minus AMSR-E low res.')
		pyplot.hist(diffm5,bins=50,alpha=0.5,label='WRF minus AMSR-E medium res.')
		pyplot.legend(loc='upper right')
		pyplot.xlabel('Wind difference [m s$^{-1}$]')
		pyplot.ylabel('Frequency')
		lensim = int(np.trunc(ds.variables['U10'][:,:,:].shape[0]/24.))
		if day == 0:
			mean4 = np.zeros(lensim)
			mean5 = np.zeros(lensim)
			median4 = np.zeros(lensim)
			median5 = np.zeros(lensim)
			mae4 = np.zeros(lensim)
			mae5 = np.zeros(lensim)
			rmse4 = np.zeros(lensim)
			rmse5 = np.zeros(lensim)
			dayarr = np.zeros(lensim)
			mean4_nudged = np.zeros(lensim)
			mean5_nudged = np.zeros(lensim)
			mae4_nudged = np.zeros(lensim)
			mae5_nudged = np.zeros(lensim)
		dayarr[day] = day
		mean4[day] = np.round(np.nanmean(diffm4),2)
		mean5[day] = np.round(np.nanmean(diffm5),2)
		median4[day] = np.round(np.nanmedian(diffm4),2)
		median5[day] = np.round(np.nanmedian(diffm5),2)
		mae4[day] = np.round(np.nanmean(np.abs(diffm4)),2)
		mae5[day] = np.round(np.nanmean(np.abs(diffm5)),2)
		rmse4[day] = np.round(sqrt(mean_squared_error(wind_wrf[:,:][~np.isnan(wind_wrf[:,:]-wind_amsre_m4[:,:])],wind_amsre_m4[:,:][~np.isnan(wind_wrf[:,:]-wind_amsre_m4[:,:])])),2)
		rmse5[day] = np.round(sqrt(mean_squared_error(wind_wrf[:,:][~np.isnan(wind_wrf[:,:]-wind_amsre_mr_m5[:,:])],wind_amsre_mr_m5[:,:][~np.isnan(wind_wrf[:,:]-wind_amsre_mr_m5[:,:])])),2)
		mean4_nudged[day] = np.round(np.nanmean(diffm4_nudged),2)
		mean5_nudged[day] = np.round(np.nanmean(diffm5_nudged),2)
		mae4_nudged[day] = np.round(np.nanmean(np.abs(diffm4_nudged)),2)
		mae5_nudged[day] = np.round(np.nanmean(np.abs(diffm5_nudged)),2)
		plt.text(0.01,0.5,'Mean: '+str(mean4[day])+'/'+str(mean5[day]),fontsize=12,transform = ax.transAxes)
		plt.text(0.01,0.45,'Median: '+str(median4[day])+'/'+str(median5[day]),fontsize=12,transform = ax.transAxes)
		plt.text(0.01,0.4,'MAE: '+str(mae4[day])+'/'+str(mae5[day]),fontsize=12,transform = ax.transAxes)
		plt.text(0.01,0.35,'RMSE: '+str(rmse4[day])+'/'+str(rmse5[day]),fontsize=12,transform = ax.transAxes)
		
		plt.savefig('Figures/SSTplots/Junk/WINDdiff_metonly_'+str(day))
		#plt.show()
		plt.close()
		
		'''
		x = wind_wrf[:,:][~np.isnan(wind_wrf[:,:]-wind_amsre_m4[:,:])]
		y = wind_amsre_m4[:,:][~np.isnan(wind_wrf[:,:]-wind_amsre_m4[:,:])]
		xy = np.vstack([x,y])
		z = gaussian_kde(xy)(xy)
		idx = z.argsort()
		x, y, z = x[idx], y[idx], z[idx]
		fig, ax = plt.subplots()
		ax.scatter(x, y, c=z, s=50, edgecolor='')
		ax.plot([0,100],[0,100],'r-',label='1:1 line')
		ax.set_xlim([0,30])
		ax.set_ylim([0,30])
		plt.show()
		'''
		
		print('Processing timestep '+str(day)+' out of '+str(lensim))
		if day == lensim-1:
			#plt.plot(dayarr,mean4,color='black',marker='o')
			#plt.plot(dayarr,mean5,color='blue',marker='o')
			#plt.plot(dayarr,mae4,color='black',marker='p')
			#plt.plot(dayarr,mae5,color='blue',marker='p')
			#plt.xlabel('Time since start of simulation [days]')
			#plt.xlim([0,27])
			#plt.savefig('Figures/SSTplots/WINDerrorpropagation.png',dpi=300)
			#plt.show()
			fig,ax = plt.subplots()
			ax.plot(dayarr,mean4,color='black',marker='.')
			ax.plot(dayarr,mean5,color='black',marker='s')
			ax.plot(dayarr,mean4_nudged,color='black',marker='.',linestyle='dashed')
			ax.plot(dayarr,mean5_nudged,color='black',marker='s',linestyle='dashed')
			ax.plot([0,30],[0,0],color='black',alpha=0.5)
			ax.set_xlabel('Time since start of simulation [days]')
			ax.set_ylabel("Bias [m s$^{-1}$]",color="black")
			ax.set_xlim([0,30])
			ax2=ax.twinx()
			ax2.plot(dayarr,mae4,color='red',marker='.')
			ax2.plot(dayarr,mae5,color='red',marker='s')
			ax2.plot(dayarr,mae4_nudged,color='red',marker='.',linestyle='dashed')
			ax2.plot(dayarr,mae5_nudged,color='red',marker='s',linestyle='dashed')
			ax2.set_ylabel("Mean absolute error [m s$^{-1}$]",color='red')
			ax2.tick_params(axis='y', labelcolor='red')
			plt.savefig('Figures/SSTplots/WINDerrorpropagation_nudgedBL.png',dpi=300)
			plt.show()

def meandailywind_final(wrffile,amsremap):
	ds = nc.Dataset(wrffile,'r')
	for day in range(0,int(np.trunc(ds.variables['U10'][:,:,:].shape[0]/24.)),1):
		if day <= 21:
			amsrefile = amsremap+'/AMSR_E_L3_DailyOcean_V04_200808'+str(10+day)+'.nc'
		else:
			amsrefile = amsremap+'/AMSR_E_L3_DailyOcean_V04_2008090'+str(day-21)+'.nc'
		es = nc.Dataset(amsrefile,'r')
		wind_wrf = (ds.variables['U10'][(day*24):(day*24)+23,:,:]**2+ds.variables['V10'][(day*24):(day*24)+23,:,:]**2)**0.5
		#wind_wrf = (ds.variables['U'][(day*24):(day*24)+23,0,:,:]**2+ds.variables['V'][(day*24):(day*24)+23,0,:,:]**2)**0.5
		wind_wrf = np.nanmean(wind_wrf,axis=0)	#get mean of day
		seaice_wrf = np.mean(ds.variables['SEAICE'][(day*24):(day*24)+23,:,:],axis=0) #get mean of day
		wind_wrf[seaice_wrf >= 0.5] = np.nan
		wind_amsre = es.variables['Med_res_wind'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
		wind_amsre = np.delete(wind_amsre,np.s_[0:540],axis=0).astype(float)
		wind_amsre[wind_amsre == -9999] = np.nan
		wind_amsre = wind_amsre/100.
		
		fig,axes = plt.subplots(2,2,figsize=(17,15), gridspec_kw = {'wspace':0.24, 'hspace':0.02})
		
		m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
		m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m1(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
		m1.contourf(lons,lats,seaice_wrf,[0.5,1.5],colors='w',zorder=1002)
		m1.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		im = m1.imshow(wind_wrf[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=20.)
		
		m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
		wind_amsre_trans = m2.transform_scalar(wind_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m2.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		m2.imshow(wind_amsre_trans[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=20.)
		
		m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
		m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m3.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		im3 = m3.imshow(wind_wrf[:,:]-wind_amsre_trans[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
		
		ax=axes[1,1]
		x = wind_wrf[:,:][~np.isnan(wind_wrf[:,:]-wind_amsre_trans[:,:])]
		y = wind_amsre_trans[:,:][~np.isnan(wind_wrf[:,:]-wind_amsre_trans[:,:])]
		r2wrf,r2amsre=x,y
		xy = np.vstack([x,y])
		z = gaussian_kde(xy)(xy)
		idx = z.argsort()
		x, y, z = x[idx], y[idx], z[idx]
		ax.scatter(x, y, c=z, cmap='viridis' ,s=5, edgecolor='')
		ax.plot([0,1000],[0,1000],'r-',label='1:1 line')
		ax.yaxis.set_label_position("right")
		ax.yaxis.tick_right()
		ax.set_xlim([0,25])
		ax.set_ylim([0,25])
		ax.set_xlabel('WRF Wind speed [m s$^{-1}$]',fontsize=13)
		ax.set_ylabel('AMSR-E Wind speed [m s$^{-1}$]',fontsize=13)
		ax.text(0.01,0.95,'Bias: '+str(np.round(np.nanmean(r2wrf-r2amsre),2)),fontsize=12,transform = ax.transAxes)
		ax.text(0.01,0.90,'MAE: '+str(np.round(np.nanmean(np.abs(r2wrf-r2amsre)),2)),fontsize=12,transform = ax.transAxes)
				
		cax = fig.add_axes([axes[0,0].get_position().x1+0.01,axes[0,0].get_position().y0+0.01,0.02,axes[0,0].get_position().y1-axes[0,0].get_position().y0-0.02])
		cax2 = fig.add_axes([axes[1,0].get_position().x1+0.01,axes[1,0].get_position().y0+0.01,0.02,axes[1,0].get_position().y1-axes[1,0].get_position().y0-0.02])
		cbar = fig.colorbar(im,cax=cax,ticks=[0,5,10,15,20])
		cbar.set_label("Wind speed [m s$^{-1}$]",fontsize=13)
		cbar2 = fig.colorbar(im3,cax=cax2,extend='both',ticks=[-5,-4,-3,-2,-1,0,1,2,3,4,5])
		cbar2.set_label("Wind speed difference (WRF minus AMSR-E) [m s$^{-1}$]",fontsize=13)		
		#plt.savefig('Figures/SSTplots/Final_WINDdiff_nudgedBL'+str(day),dpi=300)
		#plt.show()
		plt.close()
		
		diffm5 = (wind_wrf[:,:]-wind_amsre_trans[:,:])[~np.isnan(wind_wrf[:,:]-wind_amsre_trans[:,:])]
		
		lensim = int(np.trunc(ds.variables['U10'][:,:,:].shape[0]/24.))
		if day == 0:
			bias = np.zeros(lensim)
			mae = np.zeros(lensim)
			dayarr = np.zeros(lensim)
		dayarr[day] = day
		bias[day] = np.round(np.nanmean(diffm5),2)
		mae[day] = np.round(np.nanmean(np.abs(diffm5)),2)
		
		print('Processing timestep '+str(day)+' out of '+str(lensim))
		if day == lensim-1:
			dayarr=dayarr[1:]
			bias=bias[1:]
			mae=mae[1:]
			dayarrdt=np.arange('2008-08-11', '2008-09-07', dtype='datetime64[D]')
			
			myFmt = DateFormatter("%d-%b")

			fig,ax = plt.subplots()
			ax.plot(dayarrdt,bias,color='black',marker='s')
			ax.plot([min(dayarrdt),max(dayarrdt)],[0,0],color='black',alpha=0.5)
			#ax.set_xlabel('Time since start of simulation [days]')
			ax.set_ylabel("Bias [m s$^{-1}$]",color="black")
			ax.set_xlim([datetime.date(2008, 8, 11), datetime.date(2008, 9, 7)])
			print(ax.get_xticklabels())
			#ax.set_xticks(np.arange('2008-08-11', '2008-09-07', np.timedelta64(4,'D'), dtype='datetime64[D]'))
			#ax.set_xticklabels([datetime.date(2008, 8, 11),datetime.date(2008, 8, 15),datetime.date(2008, 8, 19),datetime.date(2008, 8, 23),datetime.date(2008, 8, 27),datetime.date(2008, 8, 31),datetime.date(2008, 9, 3), datetime.date(2008, 9, 7)])
			#ax.set_xticklabels(np.arange('2008-08-11', '2008-09-07', np.timedelta64(4,'D'), dtype='datetime64[D]'))
			ax.xaxis.set_major_formatter(myFmt)
			ax2=ax.twinx()
			ax2.plot(dayarrdt,mae,color='red',marker='s')
			ax2.set_ylabel("Mean Absolute Error [m s$^{-1}$]",color='red')
			ax2.tick_params(axis='y', labelcolor='red')
			plt.savefig('Figures/SSTplots/WINDerrorpropagation_v_paper.png',dpi=300)
			plt.show()
		
def ecmwf_windcomparison(ecmwfmap,amsremap):
	for day in range(0,28,1):
		if day <= 21:
			amsrefile = amsremap+'/AMSR_E_L3_DailyOcean_V04_200808'+str(10+day)+'.nc'
		else:
			amsrefile = amsremap+'/AMSR_E_L3_DailyOcean_V04_2008090'+str(day-21)+'.nc'
		es = nc.Dataset(amsrefile,'r')
		wind_amsre = es.variables['Med_res_wind'][::-1,:]  #[720:1440] (0.25x0.25deg) + flip axis
		wind_amsre = np.delete(wind_amsre,np.s_[0:540],axis=0).astype(float)
		wind_amsre[wind_amsre == -9999] = np.nan
		wind_amsre = wind_amsre/100.
		#if day <= 21:
		#	ds1 = nc.Dataset(ecmwfmap+'met_em.d01.2008-08-'+str(10+day)+'_00:00:00.nc','r')
		#	ds2 = nc.Dataset(ecmwfmap+'met_em.d01.2008-08-'+str(10+day)+'_06:00:00.nc','r')
		#	ds3 = nc.Dataset(ecmwfmap+'met_em.d01.2008-08-'+str(10+day)+'_12:00:00.nc','r')
		#	ds4 = nc.Dataset(ecmwfmap+'met_em.d01.2008-08-'+str(10+day)+'_18:00:00.nc','r')
		#else:
		#	ds1 = nc.Dataset(ecmwfmap+'met_em.d01.2008-09-'+str(day-21)+'_00:00:00.nc','r')
		#	ds2 = nc.Dataset(ecmwfmap+'met_em.d01.2008-09-'+str(day-21)+'_06:00:00.nc','r')
		#	ds3 = nc.Dataset(ecmwfmap+'met_em.d01.2008-09-'+str(day-21)+'_12:00:00.nc','r')
		#	ds4 = nc.Dataset(ecmwfmap+'met_em.d01.2008-09-'+str(day-21)+'_18:00:00.nc','r')
		#wind_ecmwf = np.nanmean([(ds1.variables['UU'][0,:,:]**2+ds1.variables['VV'][0,:,:]**2)**0.5,(ds2.variables['UU'][0,:,:]**2+ds2.variables['VV'][0,:,:]**2)**0.5,(ds3.variables['UU'][0,:,:]**2+ds3.variables['VV'][0,:,:]**2)**0.5,(ds4.variables['UU'][0,:,:]**2+ds4.variables['VV'][0,:,:]**2)**0.5],axis=0)
		#seaice_ecmwf = np.nanmean([ds1.variables['SEAICE'][:,:],ds2.variables['SEAICE'][:,:],ds3.variables['SEAICE'][:,:],ds4.variables['SEAICE'][:,:]],axis=0)
		ds = nc.Dataset(ecmwfmap+'wrffdda_d01','r')
		dswrf = nc.Dataset('/lustre/nobackup/WUR/ESG/barte035/wrfout_metnudgedBL_d01_2008-08-10_00:00:00','r')
		wind_ecmwf = np.mean((ds.variables['U_NDG_OLD'][day*4:(day*4)+3,0,:,:]**2+ds.variables['V_NDG_OLD'][day*4:(day*4)+3,0,:,:]**2)**0.5,axis=0)	
		seaice_ecmwf = np.mean(dswrf.variables['SEAICE'][(day*24):(day*24)+23,:,:],axis=0) #get mean of day
		wind_ecmwf[seaice_ecmwf >= 0.5] = np.nan
		
		fig,axes = plt.subplots(2,2,figsize=(17,15), gridspec_kw = {'wspace':0.24, 'hspace':0.02})
		
		m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
		m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m1(dswrf.variables['XLONG'][0,:,:],dswrf.variables['XLAT'][0,:,:],inverse=False)
		m1.contourf(lons,lats,seaice_ecmwf,[0.5,1.5],colors='w',zorder=1002)
		m1.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		im = m1.imshow(wind_ecmwf[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=20.)
		
		m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
		wind_amsre_trans = m2.transform_scalar(wind_amsre,lons=np.arange(-180.,180.,0.25),lats=np.arange(45.,90.,0.25),nx=249,ny=249,order=0)
		m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m2.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		m2.imshow(wind_amsre_trans[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=20.)
		
		m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
		m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m3.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		im3 = m3.imshow(wind_ecmwf[:,:]-wind_amsre_trans[:,:],cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-5.,vmax=5.)
		
		ax=axes[1,1]
		x = wind_ecmwf[:,:][~np.isnan(wind_ecmwf[:,:]-wind_amsre_trans[:,:])]
		y = wind_amsre_trans[:,:][~np.isnan(wind_ecmwf[:,:]-wind_amsre_trans[:,:])]
		r2wrf,r2amsre=x,y
		xy = np.vstack([x,y])
		z = gaussian_kde(xy)(xy)
		idx = z.argsort()
		x, y, z = x[idx], y[idx], z[idx]
		ax.scatter(x, y, c=z, cmap='viridis' ,s=5, edgecolor='')
		ax.plot([0,1000],[0,1000],'r-',label='1:1 line')
		ax.yaxis.set_label_position("right")
		ax.yaxis.tick_right()
		ax.set_xlim([0,25])
		ax.set_ylim([0,25])
		ax.set_xlabel('ECMWF Wind speed [m s$^{-1}$]',fontsize=13)
		ax.set_ylabel('AMSR-E Wind speed [m s$^{-1}$]',fontsize=13)
		ax.text(0.01,0.95,'Bias: '+str(np.round(np.nanmean(r2wrf-r2amsre),2)),fontsize=12,transform = ax.transAxes)
		ax.text(0.01,0.90,'MAE: '+str(np.round(np.nanmean(np.abs(r2wrf-r2amsre)),2)),fontsize=12,transform = ax.transAxes)
				
		cax = fig.add_axes([axes[0,0].get_position().x1+0.01,axes[0,0].get_position().y0+0.01,0.02,axes[0,0].get_position().y1-axes[0,0].get_position().y0-0.02])
		cax2 = fig.add_axes([axes[1,0].get_position().x1+0.01,axes[1,0].get_position().y0+0.01,0.02,axes[1,0].get_position().y1-axes[1,0].get_position().y0-0.02])
		cbar = fig.colorbar(im,cax=cax,ticks=[0,5,10,15,20])
		cbar.set_label("Wind speed [m s$^{-1}$]",fontsize=13)
		cbar2 = fig.colorbar(im3,cax=cax2,extend='both',ticks=[-5,-4,-3,-2,-1,0,1,2,3,4,5])
		cbar2.set_label("Wind speed difference (ECMWF minus AMSR-E) [m s$^{-1}$]",fontsize=13)		
		plt.savefig('Figures/SSTplots/Final_WINDdiff_ECMWF'+str(day),dpi=300)
		#plt.show()
		plt.close()
		
						
#WRF_wind('/lustre/nobackup/WUR/ESG/barte035/wrfout_metonly_d01_2008-08-10_00:00:00')
#AMSRE_wind('/home/WUR/barte035/WRFChem/AMSR-E/AMSR_E_L3_DailyOcean_V04_20080810.nc')
#meandailywind('/lustre/nobackup/WUR/ESG/barte035/wrfout_metonlyv3_d01_2008-08-10_00:00:00','/home/WUR/barte035/WRFChem/AMSR-E/')
meandailywind_final('/lustre/backup/WUR/ESG/barte035/wrfout_metonlyv3_d01_2008-08-10_00:00:00','/home/WUR/barte035/WRFChem/AMSR-E/')
#meandailywind_final('/lustre/backup/WUR/ESG/barte035/wrfout_metnudgedBL_d01_2008-08-10_00:00:00','/home/WUR/barte035/WRFChem/AMSR-E/')
#ecmwf_windcomparison('/home/WUR/barte035/WRFChem/WRFV411-Polar/run/','/home/WUR/barte035/WRFChem/AMSR-E/')
