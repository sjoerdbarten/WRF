import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import datetime as dt
from sklearn.metrics import mean_squared_error
from math import sqrt
from scipy.stats import gaussian_kde

def coareg_analysis():
	bf1 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-10_00:00:00','r')
	bf2 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-13_01:00:00','r')
	bf3 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-21_01:00:00','r')
	bf4 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-31_01:00:00','r')
	cf1 = nc.Dataset('/lustre/nobackup/WUR/ESG/barte035/wrfout_coaregv2_d01_2008-08-10_00:00:00','r')
	cf2 = nc.Dataset('/lustre/nobackup/WUR/ESG/barte035/wrfout_coaregv2_d01_2008-08-15_01:00:00','r')	
	wrfdata_sst_base = np.concatenate((bf1['SST'][1:73,:,:],bf2['SST'][0:192,:,:],bf3['SST'][0:240,:,:],bf4['SST'][:,:,:]),axis=0)
	wrfdata_sst_coareg = np.concatenate((cf1['SST'][1:121,:,:],cf2['SST'][:,:,:]),axis=0)
	wrfdata_dep_base = np.concatenate((bf1['DEP_VEL'][1:73,:,:],bf2['DEP_VEL'][0:192,:,:],bf3['DEP_VEL'][0:240,:,:],bf4['DEP_VEL'][:,:,:]),axis=0)
	wrfdata_dep_coareg = np.concatenate((cf1['VTC'][1:121,:,:],cf2['VTC'][:,:,:]),axis=0)
	wrfdata_lu_base = np.concatenate((bf1['LU_INDEX'][1:73,:,:],bf2['LU_INDEX'][0:192,:,:],bf3['LU_INDEX'][0:240,:,:],bf4['LU_INDEX'][:,:,:]),axis=0)
	wrfdata_lu_coareg = np.concatenate((cf1['LU_INDEX'][1:121,:,:],cf2['LU_INDEX'][:,:,:]),axis=0)
	wrfdata_u10_base = np.concatenate(((bf1['U10'][1:73,:,:]**2+bf1['V10'][1:73,:,:]**2)**0.5,(bf2['U10'][0:192,:,:]**2+bf2['V10'][0:192,:,:]**2)**0.5,(bf3['U10'][0:240,:,:]**2+bf3['V10'][0:240,:,:]**2)**0.5,(bf4['U10'][:,:,:]**2+bf4['V10'][:,:,:]**2)**0.5),axis=0)
	wrfdata_u10_coareg = np.concatenate(((cf1['U10'][1:121,:,:]**2+cf1['V10'][1:121,:,:]**2)**0.5,(cf2['U10'][:,:,:]**2+cf2['V10'][:,:,:]**2)**0.5),axis=0)
		
	print('WE ARE HERE')
	for i in range(wrfdata_dep_base.shape[0]):
		wrfdata_sst_base[i,:,:] = np.where(wrfdata_sst_base[i,:,:] >= 271.,np.where(wrfdata_lu_base[i,:,:] == 16, wrfdata_sst_base[i,:,:], np.nan),np.nan)
		wrfdata_sst_coareg[i,:,:] = np.where(wrfdata_sst_coareg[i,:,:] >= 271.,np.where(wrfdata_lu_coareg[i,:,:] == 16, wrfdata_sst_coareg[i,:,:], np.nan),np.nan)
		wrfdata_dep_base[i,:,:] = np.where(wrfdata_sst_base[i,:,:] >= 271.,np.where(wrfdata_lu_base[i,:,:] == 16, wrfdata_dep_base[i,:,:], np.nan),np.nan)
		wrfdata_dep_coareg[i,:,:] = np.where(wrfdata_sst_coareg[i,:,:] >= 271.,np.where(wrfdata_lu_coareg[i,:,:] == 16, wrfdata_dep_coareg[i,:,:], np.nan),np.nan)
		wrfdata_u10_base[i,:,:] = np.where(wrfdata_sst_base[i,:,:] >= 271.,np.where(wrfdata_lu_base[i,:,:] == 16, wrfdata_u10_base[i,:,:], np.nan),np.nan)
		wrfdata_u10_coareg[i,:,:] = np.where(wrfdata_sst_coareg[i,:,:] >= 271.,np.where(wrfdata_lu_coareg[i,:,:] == 16, wrfdata_u10_coareg[i,:,:], np.nan),np.nan)
	print('WE ARE HERE')
	
	wrfdata_sst_base = wrfdata_sst_base[~np.isnan(wrfdata_sst_base)]
	wrfdata_sst_coareg = wrfdata_sst_coareg[~np.isnan(wrfdata_sst_coareg)]
	wrfdata_dep_base = wrfdata_dep_base[~np.isnan(wrfdata_dep_base)]
	wrfdata_dep_coareg = wrfdata_dep_coareg[~np.isnan(wrfdata_dep_coareg)]
	wrfdata_u10_base = wrfdata_u10_base[~np.isnan(wrfdata_u10_base)]
	wrfdata_u10_coareg = wrfdata_u10_coareg[~np.isnan(wrfdata_u10_coareg)]
	
	schmidtnum_base = (44./48.)**0.5*np.exp(-0.055*wrfdata_sst_base+22.63)
	schmidtnum_coareg = (44./48.)**0.5*np.exp(-0.055*wrfdata_sst_coareg+22.63)
		
	print('WE ARE HERE!!!')
	fig, axes = plt.subplots(2, 2)
	ax1 = axes[0,0]
	ax1.scatter(wrfdata_sst_base,wrfdata_dep_base*100.,color='red')
	ax1.scatter(wrfdata_sst_coareg,wrfdata_dep_coareg*100.,color='blue')
	ax1.set_xlabel('Sea surface temperature [K]')
	ax1.set_ylabel('Deposition velocity [cm s$^{-1}$]')
	ax2 = axes[0,1]
	ax2.scatter(wrfdata_u10_base,wrfdata_dep_base*100.,color='red')
	ax2.scatter(wrfdata_u10_coareg,wrfdata_dep_coareg*100.,color='blue')
	ax2.set_xlabel('10m wind speed [m s$^{-1}$]')
	ax2.set_ylabel('Deposition velocity [cm s$^{-1}$]')
	ax3 = axes[1,0]
	ax3.scatter(wrfdata_sst_base,schmidtnum_base,color='red')
	ax3.scatter(wrfdata_sst_coareg,schmidtnum_coareg,color='blue')
	ax3.set_xlabel('Sea surface temperature [K]')
	ax3.set_ylabel('Schmidt number [-]')
	ax4 = axes[1,1]
	ax4.scatter(wrfdata_u10_base,(wrfdata_dep_base*100.)*(schmidtnum_base/660.)**0.5,color='red')
	ax4.scatter(wrfdata_u10_coareg,(wrfdata_dep_coareg*100.)*(schmidtnum_coareg/660.)**0.5,color='blue')
	ax4.set_xlabel('10m wind speed [m s$^{-1}$]')
	ax4.set_ylabel('Deposition velocity normalized for Sc$_{660}$ [-]')
	plt.savefig('Figures/COAREGplots/sstwinddependency.png')
	plt.show()

def iodide_analysis():
	cf1 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_coaregv2_d01_2008-08-10_00:00:00','r')
	cf2 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_coaregv2_d01_2008-08-15_01:00:00','r')	
	wrfdata_sst_coareg = np.concatenate((cf1['SST'][1:121,:,:],cf2['SST'][:,:,:]),axis=0)
	wrfdata_dep_coareg = np.concatenate((cf1['VTC'][1:121,:,:],cf2['VTC'][:,:,:]),axis=0)
	wrfdata_lu_coareg = np.concatenate((cf1['LU_INDEX'][1:121,:,:],cf2['LU_INDEX'][:,:,:]),axis=0)
	
	for i in range(wrfdata_dep_coareg.shape[0]):
		wrfdata_sst_coareg[i,:,:] = np.where(wrfdata_sst_coareg[i,:,:] >= 271.,np.where(wrfdata_lu_coareg[i,:,:] == 16, wrfdata_sst_coareg[i,:,:], np.nan),np.nan)
		wrfdata_dep_coareg[i,:,:] = np.where(wrfdata_sst_coareg[i,:,:] >= 271.,np.where(wrfdata_lu_coareg[i,:,:] == 16, wrfdata_dep_coareg[i,:,:], np.nan),np.nan)
	
	wrfdata_sst_aug = wrfdata_sst_coareg[:503,:,:]
	wrfdata_dep_aug = wrfdata_dep_coareg[:503,:,:]
	wrfdata_sst_sep = wrfdata_sst_coareg[503:,:,:]
	wrfdata_dep_sep = wrfdata_dep_coareg[503:,:,:]

	wrfdata_sst_aug_mean = np.nanmean(wrfdata_sst_aug,axis=0)
	wrfdata_dep_aug_mean = np.nanmean(wrfdata_dep_aug,axis=0)
	wrfdata_sst_sep_mean = np.nanmean(wrfdata_sst_sep,axis=0)
	wrfdata_dep_sep_mean = np.nanmean(wrfdata_dep_sep,axis=0)
	wrf_iodide_aug = 1.46e15*np.exp(-9134./wrfdata_sst_aug_mean)
	wrf_iodide_sep = 1.46e15*np.exp(-9134./wrfdata_sst_sep_mean)

	#Load data
	iodidepath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/Iodide/predicted_iodide_0.125x0.125_Ns_All_Ensemble_members.nc'
	Iodide = nc.Dataset(iodidepath,'r')
	iodide_data_aug = Iodide.variables['Ensemble_Monthly_mean'][7,:,:] #BE SURE TO SELECT THE RIGHT DATA (Data is monthly, from January (=0) to December (=11))
	iodide_data_sep = Iodide.variables['Ensemble_Monthly_mean'][8,:,:] #BE SURE TO SELECT THE RIGHT DATA (Data is monthly, from January (=0) to December (=11))

	#Define basemap
	geopath = '/home/WUR/barte035/WRFChem/WPS/geo_em.d01.nc'
	geo = nc.Dataset(geopath,'r')
	wrflat = geo.variables['XLAT_M'][0,:,:]
	wrflon = geo.variables['XLONG_M'][0,:,:]
	we  = geo.variables['XLAT_M'].shape[2]
	sn  = geo.variables['XLAT_M'].shape[1]
	lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
	lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
	dx_dom = 30000

	'''
	fig,axes = plt.subplots(2,3,figsize=(20,10))
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(cf1.variables['XLONG'][0,:,:],cf1.variables['XLAT'][0,:,:],inverse=False)
	m1.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m1.imshow(wrf_iodide_aug[:,:],vmin=0,vmax=130,cmap='viridis_r',zorder=1000)
		
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
	#Select correct datapoints
	iodide_data_aug = np.delete(iodide_data_aug,np.s_[0:1081],axis=0)
	Iodide_trans_aug = m2.transform_scalar(iodide_data_aug,lons=np.arange(-180.,180.,0.125),lats=np.arange(45.,90.,0.125),nx=249,ny=249,order=1)
	Iodide_trans_aug = np.where(geo.variables['LANDMASK'][0,:,:] != 1, Iodide_trans_aug, np.nan)
	m2.imshow(Iodide_trans_aug,vmin=0,vmax=130,cmap='viridis_r')
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	
	m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,2])
	m3.imshow(wrf_iodide_aug[:,:]-Iodide_trans_aug,vmin=-50,vmax=50,cmap='PiYG')
	m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m3.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		
	m4 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
	m4.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m4.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m4.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m4.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m4.imshow(wrf_iodide_sep[:,:],vmin=0,vmax=130,cmap='viridis_r',zorder=1000)
		
	m5 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,1])
	#Select correct datapoints
	iodide_data_sep = np.delete(iodide_data_sep,np.s_[0:1081],axis=0)
	Iodide_trans_sep = m4.transform_scalar(iodide_data_sep,lons=np.arange(-180.,180.,0.125),lats=np.arange(45.,90.,0.125),nx=249,ny=249,order=1)
	Iodide_trans_sep = np.where(geo.variables['LANDMASK'][0,:,:] != 1, Iodide_trans_sep, np.nan)
	m5.imshow(Iodide_trans_aug,vmin=0,vmax=130,cmap='viridis_r')
	m5.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m5.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m5.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m5.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)

	m6 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,2])
	m6.imshow(wrf_iodide_sep[:,:]-Iodide_trans_sep,vmin=-50,vmax=50,cmap='PiYG')
	m6.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m6.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m6.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m6.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#plt.show()
	plt.close()
	'''


	fig,axes = plt.subplots(1,2,figsize=(20,10))
	print(axes)
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(cf1.variables['XLONG'][0,:,:],cf1.variables['XLAT'][0,:,:],inverse=False)
	m1.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im = m1.imshow(wrf_iodide_aug[:,:],vmin=0,vmax=120,cmap='viridis_r',zorder=1000)
		
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1])
	iodide_data_aug = np.delete(iodide_data_aug,np.s_[0:1081],axis=0)
	Iodide_trans_aug = m2.transform_scalar(iodide_data_aug,lons=np.arange(-180.,180.,0.125),lats=np.arange(45.,90.,0.125),nx=249,ny=249,order=1)
	Iodide_trans_aug = np.where(geo.variables['LANDMASK'][0,:,:] != 1, Iodide_trans_aug, np.nan)
	m2.imshow(Iodide_trans_aug,vmin=0,vmax=120,cmap='viridis_r')
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.fillcontinents(color='coral',lake_color='aqua',zorder=1001,alpha=1)
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	
	cax = fig.add_axes([axes[0].get_position().x1+0.01,axes[0].get_position().y0+0.05,0.02,axes[0].get_position().y1-axes[0].get_position().y0-0.1])
	cbar = fig.colorbar(im,cax=cax,ticks=[0,30,60,90,120],extend='max')
	cbar.set_label("[I$^{-}$]$_{sw}$ [nM]",fontsize=13)
	plt.show()


#coareg_analysis()
iodide_analysis()
