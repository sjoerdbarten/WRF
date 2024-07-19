from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib import pyplot
from sklearn.metrics import mean_squared_error
from math import sqrt
from scipy.stats import gaussian_kde
from sklearn.metrics import r2_score
import matplotlib.colors as colors
import pandas as pd
import datetime as dt


def dep_vel():
	bf1 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-10_00:00:00','r')
	bf2 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-13_01:00:00','r')
	bf3 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-21_01:00:00','r')
	bf4 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-31_01:00:00','r')
	wrfdata_base = np.concatenate((bf1['DEP_VEL'][0:73,:,:]*100.,bf2['DEP_VEL'][0:192,:,:]*100.,bf3['DEP_VEL'][0:240,:,:]*100.,bf4['DEP_VEL'][:,:,:]*100.),axis=0)
	cf1 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chem_nudgedBL_coareg_sherweniodide_d01_2008-08-10_00:00:00','r')
	cf2 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chem_nudgedBL_coareg_sherweniodide_d01_2008-08-11_01:00:00','r')
	wrfdata_fixeddep = np.concatenate((cf1['DEP_VEL'][0:25,:,:]*100.+cf1['VTC'][0:25,:,:]*100.,cf2['DEP_VEL'][:,:,:]*100.+cf2['VTC'][:,:,:]*100.),axis=0)
	wrfdata_fullfile = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_fixeddep_d01_2008-08-10_00:00:00','r')
	
	for day in range(0,int(np.trunc(wrfdata_fullfile.variables['DEP_VEL'][:,:,:].shape[0]/24.)),1):
		seaice_wrf = np.mean(wrfdata_fullfile.variables['SEAICE'][(day*24):(day*24)+23,:,:],axis=0) #get mean of day
		depvel_base = wrfdata_base[(day*24):(day*24)+23,:,:]
		depvel_base = np.nanmean(depvel_base,axis=0)	#get mean of day
		depvel_fixeddep = wrfdata_fixeddep[(day*24):(day*24)+23,:,:]
		depvel_fixeddep = np.nanmean(depvel_fixeddep,axis=0)	#get mean of day
	
		fig,axes = plt.subplots(2,2,figsize=(17,15), gridspec_kw = {'wspace':0.24, 'hspace':0.02})
		
		m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
		m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m1(wrfdata_fullfile.variables['XLONG'][0,:,:],wrfdata_fullfile.variables['XLAT'][0,:,:],inverse=False)
		m1.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
		#m1.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		im = m1.imshow(depvel_base[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=0.01)
		
		m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
		m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m2.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
		#m2.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		m2.imshow(depvel_fixeddep[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=0.01)
		
		m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
		m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m3.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
		#m3.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		im3 = m3.imshow((depvel_base[:,:]-depvel_fixeddep[:,:])/depvel_base[:,:]*100.,cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-100.,vmax=100.)
		
		ax=axes[1,1]
		ax.hist(depvel_fixeddep[:,:].ravel(),bins=50,range=(0.,1.),normed=True,cumulative=True,alpha=0.5,color='green')
		ax.hist(depvel_base[:,:].ravel(),bins=50,range=(0.,1.),normed=True,cumulative=True,alpha=0.5,color='black')
		ax.yaxis.set_label_position("right")
		ax.yaxis.tick_right()
		ax.set_xlim([0,1])
		ax.set_ylim([0,1])
		ax.set_xlabel('Deposition velocity [cm s$^{-1}$]',fontsize=13)
		ax.set_ylabel('Cumulative frequency [-]',fontsize=13)		
		
		cax = fig.add_axes([axes[0,0].get_position().x1+0.01,axes[0,0].get_position().y0+0.01,0.02,axes[0,0].get_position().y1-axes[0,0].get_position().y0-0.02])
		cax2 = fig.add_axes([axes[1,0].get_position().x1+0.01,axes[1,0].get_position().y0+0.01,0.02,axes[1,0].get_position().y1-axes[1,0].get_position().y0-0.02])
		cbar = fig.colorbar(im,cax=cax,ticks=[0,0.25,0.5,0.75,1.0])
		cbar.set_label("Deposition velocity [cm s$^{-1}$]",fontsize=13)
		cbar2 = fig.colorbar(im3,cax=cax2,extend='both',ticks=[-100,-75,-50,-25,0,25,50,75,100])
		cbar2.set_label("Deposition velocity difference (Base minus Adjusted) [%]",fontsize=13)		
		plt.savefig('Figures/DEPVELplots/DEP_VEL'+str(day),dpi=150)
		#plt.show()
		plt.close()

def ozone():
	#BASE
	bf1 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-10_00:00:00','r')
	bf2 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-13_01:00:00','r')
	bf3 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-21_01:00:00','r')
	bf4 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-31_01:00:00','r')
	wrfdata_base = np.concatenate((bf1['o3'][0:73,0,:,:]*1000.,bf2['o3'][0:192,0,:,:]*1000.,bf3['o3'][0:240,0,:,:]*1000.,bf4['o3'][:,0,:,:]*1000.),axis=0)
	#Nudged
	wrfdata_nudgedBL_nofixeddep = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chem_nudgedBL_nofixeddep_d01_2008-08-10_00:00:00','r')
	wrfdata_nudged = wrfdata_nudgedBL_nofixeddep.variables['o3'][:,0,:,:]*1000.
	#NudgedCOAREG
	df1 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chem_nudgedBL_coareg_sherweniodide_d01_2008-08-10_00:00:00','r')
	df2 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chem_nudgedBL_coareg_sherweniodide_d01_2008-08-11_01:00:00','r')
	wrfo3data_nudgedBL_coareg = np.concatenate((df1['o3'][0:25,0,:,:]*1000.,df2['o3'][:,0,:,:]*1000.),axis=0)
	#Fullfile
	wrfdata_fullfile = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_fixeddep_d01_2008-08-10_00:00:00','r')
	#CAMS
	camsdata = np.load('/home/WUR/barte035/WRFChem/o3_analysis_DATA/CAMS-MACC/MACCnpo3array.npy')
	camsdatanew = np.load('/home/WUR/barte035/WRFChem/o3_analysis_DATA/CAMS-MACC/MACCnpo3arraynew.npy')
	
	run1_o3mean = np.nanmean(wrfdata_base[:,:,:],axis=0)
	run2_o3mean = np.nanmean(wrfdata_nudged[:,:,:],axis=0)
	run3_o3mean = np.nanmean(wrfo3data_nudgedBL_coareg[:,:,:],axis=0)
	cams_o3mean = np.nanmean(camsdata[:,:,:],axis=0)
	camsnew_o3mean = np.nanmean(camsdatanew[:,:,:],axis=0)
	
	fig,axes = plt.subplots(2,3,figsize=(17,11.33), gridspec_kw = {'wspace':0.02, 'hspace':0.022})
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(wrfdata_fullfile.variables['XLONG'][0,:,:],wrfdata_fullfile.variables['XLAT'][0,:,:],inverse=False)
	m1.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im = m1.imshow(run1_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	axes[0,0].set_title("DEFAULT")
		
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m2.imshow(run2_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	axes[0,1].set_title("NUDGED")
	
	m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,2])
	m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m3.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m3.imshow(run3_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	axes[0,2].set_title("COAREG")
	
	m4 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
	m4.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m4.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m4.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m4.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m4.imshow(cams_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	axes[1,0].set_title("CAMS")
	
	#load obs
	stationdata = pd.read_csv('/home/WUR/barte035/WRFChem/o3_analysis_DATA/o3_station_data_new.csv')
	stationdata.index = stationdata['Unnamed: 0']
	stationdata.index = pd.to_datetime(stationdata.index)
	stationdata.drop(['Unnamed: 0'],axis=1)
	wrfstarttime = dt.datetime(2008,8,10,00,00,00)
	maxtimestep = wrfdata_base.shape[0]
	timearr = np.array([wrfstarttime + dt.timedelta(hours=i) for i in xrange(maxtimestep)])
	stationdata = stationdata[timearr[0]:timearr[-1]]
	loclist = ['barrow','storhofdi','summit','ahtari','bredkalen','esrange','karasjok','karvatn','lerwick','oulanka','pallas','tustervatn','villum','vindeln','virolahti','zeppelin','ascos','whitehorse','yellowknife','normanwells','fortliard','inuvik','denalinp','alert','hurdal']
	coordlist = [[188,96],[30,89],[74,84],[30,166],[29,150],[46,154],[54,157],[22,140],[11,122],[46,168],[48,158],[36,146],[93,115],[32,157],[25,176],[85,132],[114,123],[202,46],[167,29],[180,49],[186,30],[180,65],[209,73],[111,99],[14,146]]
	coordlistlatslons = [[71.3230,-156.6114],[63.400,-20.288],[72.5800018311,-38.4799995422],[62.583333,24.183333],[63.85,15.333333],[67.883333,21.066667],[69.466667,25.216667],[62.783333,8.883333],[60.13922,-1.185319],[66.320278,29.401667],[67.973333333,24.116111111],[65.833333,13.916667],[81.6,-16.67],[64.25,19.766667],[60.526667,27.686111],[78.90715,11.88668],[87.5178,-7.882],[60.718609,-135.049193],[62.45207,-114.364],[65.27926,-126.813],[60.23583,-123.467],[68.36005,-133.727],[63.72,-148.97],[82.4991455078,-62.3415260315],[60.372386,11.078142]]
	
	m5 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,1])
	m5.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	#m5.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m5.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m5.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m5.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=999)
	m5.drawmapboundary(fill_color='cornflowerblue',zorder=999)
	lonswrf,latswrf = m5(wrfdata_fullfile.variables['XLONG'][0,:,:],wrfdata_fullfile.variables['XLAT'][0,:,:],inverse=False)
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.1,1.5],colors='w',zorder=2,alpha=0.1)
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.2,1.5],colors='w',zorder=2,alpha=0.2)
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.3,1.5],colors='w',zorder=2,alpha=0.3)
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.4,1.5],colors='w',zorder=2,alpha=0.4)
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=2,alpha=0.5)
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.6,1.5],colors='w',zorder=2,alpha=0.6)
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.7,1.5],colors='w',zorder=2,alpha=0.7)
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.8,1.5],colors='w',zorder=2,alpha=0.8)
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.9,1.5],colors='w',zorder=2,alpha=0.9)
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SEAICE'][0,:,:],[1.0,1.5],colors='w',zorder=2,alpha=1.0)   
    	m5.contourf(lonswrf,latswrf,wrfdata_fullfile.variables['SNOWC'][0,:,:],[0.5,1.5],colors='w',zorder=2,alpha=0.9)
	for i in range(len(loclist)):
		print(coordlistlatslons[i][1],coordlistlatslons[i][0],stationdata[loclist[i]+'_o3'].mean())
		m5.scatter(coordlistlatslons[i][1],coordlistlatslons[i][0],c=stationdata[loclist[i]+'_o3'].mean(),s=65,marker='o',cmap='jet',zorder=1003,latlon=True,vmin=0,vmax=50.)
		#m5.scatter(coordlistlats[i][1],coordlistlats[i][0],100,marker='o',color='Red',zorder=1003,latlon=True)
	axes[1,1].set_title("Observations")
	
	ax=axes[1,2]
	#density = gaussian_kde(wrfdata_base[:,:,:].ravel())
	#density.covariance_factor = lambda : .25
	#density._compute_covariance()
	#density2 = gaussian_kde(wrfdata_nudged[:,:,:].ravel())
	#density2.covariance_factor = lambda : .25
	#density2._compute_covariance()
	#density3 = gaussian_kde(wrfo3data_nudgedBL_coareg[:,:,:].ravel())
	#density3.covariance_factor = lambda : .25
	#density3._compute_covariance()
	#densitycams = gaussian_kde(camsdata[:,:,:].ravel())
	#densitycams.covariance_factor = lambda : .25
	#densitycams._compute_covariance()
	#densitycamsnew = gaussian_kde(camsdatanew[:,:,:].ravel())
	#densitycamsnew.covariance_factor = lambda : .25
	#densitycamsnew._compute_covariance()
	#xs = np.linspace(0,60,120)
	#np.save('/home/WUR/barte035/WRFChem/Python/Arrays/xso3dens.npy', xs)
	#np.save('/home/WUR/barte035/WRFChem/Python/Arrays/densityxs.npy', density(xs))
	#np.save('/home/WUR/barte035/WRFChem/Python/Arrays/density2xs.npy', density2(xs))
	#np.save('/home/WUR/barte035/WRFChem/Python/Arrays/density3xs.npy', density3(xs))
	#np.save('/home/WUR/barte035/WRFChem/Python/Arrays/densitycamsxs.npy', densitycams(xs))
	#np.save('/home/WUR/barte035/WRFChem/Python/Arrays/densitycamsnewxs.npy', densitycamsnew(xs))
	xs = np.load('/home/WUR/barte035/WRFChem/Python/Arrays/xso3dens.npy')
	density = np.load('/home/WUR/barte035/WRFChem/Python/Arrays/densityxs.npy')
	density2 = np.load('/home/WUR/barte035/WRFChem/Python/Arrays/density2xs.npy')
	density3 = np.load('/home/WUR/barte035/WRFChem/Python/Arrays/density3xs.npy')
	densitycams = np.load('/home/WUR/barte035/WRFChem/Python/Arrays/densitycamsxs.npy')
	densitycamsnew = np.load('/home/WUR/barte035/WRFChem/Python/Arrays/densitycamsnewxs.npy')
	#ax.plot(xs,density(xs),color='red',label='DEFAULT')
	#ax.plot(xs,density2(xs),color='orange',label='NUDGED')
	#ax.plot(xs,density3(xs),color='cyan',label='COAREG')
	#ax.plot(xs,densitycams(xs),color='blue',label='CAMS')
	ax.plot(xs,density,color='red',label='DEFAULT')
	ax.plot(xs,density2,color='orange',label='NUDGED')
	ax.plot(xs,density3,color='cyan',label='COAREG')
	ax.plot(xs,densitycams,color='blue',label='CAMS')
	ax.legend(loc=1)
	#ax.yaxis.set_label_position("left")
	#ax.yaxis.tick_right()
	ax.set_xlim([0,60])
	ax.set_ylim([0,0.08])
	ax.set_xlabel('Mean surface ozone mixing ratio [ppb]',fontsize=13)
	ax.yaxis.tick_right()
	ax.yaxis.set_ticks_position('both')
	ax.yaxis.set_label_position("right")
	ax.set_ylabel('Frequency [-]',fontsize=13)
	
	cax = fig.add_axes([axes[0,2].get_position().x1+0.008,axes[0,2].get_position().y0+0.005,0.02,axes[0,2].get_position().y1-axes[0,2].get_position().y0-0.015])
	cbar = fig.colorbar(im,cax=cax,ticks=[0,5,10,15,20,25,30,35,40,45,50],extend='max')
	cbar.set_label("Mean surface ozone mixing ratio [ppb]",fontsize=13)
	#fig.tight_layout()
	xtext,ytext = m1(-133.,48.)
	axes[0,0].text(xtext,ytext,s='(a)',fontsize=14,zorder=1010)
	axes[0,1].text(xtext,ytext,s='(b)',fontsize=14,zorder=1010)
	axes[0,2].text(xtext,ytext,s='(c)',fontsize=14,zorder=1010)
	axes[1,0].text(xtext,ytext,s='(d)',fontsize=14,zorder=1010)
	axes[1,1].text(xtext,ytext,s='(e)',fontsize=14,zorder=1010)
	axes[1,2].text(2.5,0.075,s='(f)',fontsize=14,zorder=1010)
	#plt.savefig('Figures/DEPVELplots/MEAN_OZONE_final2',dpi=150)
	plt.show()
	#plt.close()
	
	#NEW PLOT
	fig,axes = plt.subplots(2,2,figsize=(17,17), gridspec_kw = {'wspace':0.02, 'hspace':0.022})
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(wrfdata_fullfile.variables['XLONG'][0,:,:],wrfdata_fullfile.variables['XLAT'][0,:,:],inverse=False)
	m1.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im = m1.imshow(run1_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	for i in range(len(loclist)):
		print(coordlistlatslons[i][1],coordlistlatslons[i][0],stationdata[loclist[i]+'_o3'].mean())
		m1.scatter(coordlistlatslons[i][1],coordlistlatslons[i][0],c=stationdata[loclist[i]+'_o3'].mean(),s=85,marker='o',cmap='jet',zorder=1003,latlon=True,vmin=0,vmax=50.,linewidth=1.5)
	axes[0,0].set_title("DEFAULT",fontsize=20)		
	
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m2.imshow(run3_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	for i in range(len(loclist)):
		print(coordlistlatslons[i][1],coordlistlatslons[i][0],stationdata[loclist[i]+'_o3'].mean())
		m2.scatter(coordlistlatslons[i][1],coordlistlatslons[i][0],c=stationdata[loclist[i]+'_o3'].mean(),s=85,marker='o',cmap='jet',zorder=1003,latlon=True,vmin=0,vmax=50.,linewidth=1.5)
	axes[0,1].set_title("COAREG",fontsize=20)
	
	m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
	m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m3.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#m3.imshow(cams_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	m3.imshow(camsnew_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	for i in range(len(loclist)):
		print(coordlistlatslons[i][1],coordlistlatslons[i][0],stationdata[loclist[i]+'_o3'].mean())
		m3.scatter(coordlistlatslons[i][1],coordlistlatslons[i][0],c=stationdata[loclist[i]+'_o3'].mean(),s=85,marker='o',cmap='jet',zorder=1003,latlon=True,vmin=0,vmax=50.,linewidth=1.5)
	axes[1,0].set_title("CAMS",fontsize=20)
		
	ax=axes[1,1]
	ax.plot(xs,density,color='red',label='DEFAULT',linewidth=3)
	ax.plot(xs,density2,color='orange',label='NUDGED',linewidth=3)
	ax.plot(xs,density3,color='green',label='COAREG',linewidth=3)
	#ax.plot(xs,densitycams,color='blue',label='CAMS',linewidth=3)
	ax.plot(xs,densitycamsnew,color='blue',label='CAMS',linewidth=3)
	ax.legend(loc=1,fontsize=18)
	ax.set_xlim([0,60])
	ax.set_ylim([0,0.07])
	ax.set_xlabel('Surface ozone mixing ratio [ppb]',fontsize=18)
	ax.set_ylabel('Frequency [-]',fontsize=18)
	ax.yaxis.tick_right()
	ax.yaxis.set_ticks_position('both')
	ax.yaxis.set_label_position("right")
	ax.set_xticks([0,5,10,15,20,25,30,35,40,45,50,55,60])
	ax.set_xticklabels([0,5,10,15,20,25,30,35,40,45,50,55,60],size=16)
	ax.set_yticks([0.01,0.02,0.03,0.04,0.05,0.06,0.07])
	ax.set_yticklabels([0.01,0.02,0.03,0.04,0.05,0.06,0.07],size=16)
	cax = fig.add_axes([axes[0,1].get_position().x1+0.008,axes[0,1].get_position().y0+0.005,0.02,axes[0,1].get_position().y1-axes[0,1].get_position().y0-0.015])
	cbar = fig.colorbar(im,cax=cax,ticks=[0,5,10,15,20,25,30,35,40,45,50],extend='max')
	cbar.set_label("Mean surface ozone mixing ratio [ppb]",fontsize=18)
	cbar.ax.tick_params(labelsize=16)
	#fig.tight_layout()
	xtext,ytext = m1(-133.,48.)
	axes[0,0].text(xtext,ytext,s='(a)',fontsize=15,zorder=1010)
	axes[0,1].text(xtext,ytext,s='(b)',fontsize=15,zorder=1010)
	axes[1,0].text(xtext,ytext,s='(c)',fontsize=15,zorder=1010)
	axes[1,1].text(2.5,0.066,s='(d)',fontsize=15,zorder=1010)
	#plt.savefig('Figures/DEPVELplots/MEAN_OZONE_final_4panel_new',dpi=500)
	#plt.savefig('Figures/DEPVELplots/MEAN_OZONE_final_4panel_newcamsdata_new',dpi=500)
	plt.show()
	#plt.close()
	
	#NEW PLOT
	fig,axes = plt.subplots(1,2,figsize=(17,8), gridspec_kw = {'wspace':0.02, 'hspace':0.022})
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(wrfdata_fullfile.variables['XLONG'][0,:,:],wrfdata_fullfile.variables['XLAT'][0,:,:],inverse=False)
	m1.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im = m1.imshow(run2_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	for i in range(len(loclist)):
		print(coordlistlatslons[i][1],coordlistlatslons[i][0],stationdata[loclist[i]+'_o3'].mean())
		m1.scatter(coordlistlatslons[i][1],coordlistlatslons[i][0],c=stationdata[loclist[i]+'_o3'].mean(),s=85,marker='o',cmap='jet',zorder=1003,latlon=True,vmin=0,vmax=50.,linewidth=1.5)
	axes[0].set_title("NUDGED",fontsize=20)		
	
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1])
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m2.imshow(run3_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	for i in range(len(loclist)):
		print(coordlistlatslons[i][1],coordlistlatslons[i][0],stationdata[loclist[i]+'_o3'].mean())
		m2.scatter(coordlistlatslons[i][1],coordlistlatslons[i][0],c=stationdata[loclist[i]+'_o3'].mean(),s=85,marker='o',cmap='jet',zorder=1003,latlon=True,vmin=0,vmax=50.,linewidth=1.5)
	axes[1].set_title("COAREG",fontsize=20)

	cax = fig.add_axes([axes[1].get_position().x1+0.008,axes[1].get_position().y0+0.005,0.02,axes[1].get_position().y1-axes[1].get_position().y0-0.015])
	cbar = fig.colorbar(im,cax=cax,ticks=[0,5,10,15,20,25,30,35,40,45,50],extend='max')
	cbar.set_label("Mean surface ozone mixing ratio [ppb]",fontsize=18)
	cbar.ax.tick_params(labelsize=16)
	#fig.tight_layout()
	xtext,ytext = m1(-133.,48.)
	axes[0].text(xtext,ytext,s='(a)',fontsize=15,zorder=1010)
	axes[1].text(xtext,ytext,s='(b)',fontsize=15,zorder=1010)
	#plt.savefig('Figures/DEPVELplots/MEAN_OZONE_final_4panel_new',dpi=500)
	#plt.savefig('Figures/DEPVELplots/MEAN_OZONE_final_4panel_newcamsdata_new_2panel',dpi=500)
	plt.show()
	#plt.close()
	
	#NEW PLOT TO GET DIFF (PANEL 2)
	fig,axes = plt.subplots(1,2,figsize=(17,8), gridspec_kw = {'wspace':0.02, 'hspace':0.022})
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(wrfdata_fullfile.variables['XLONG'][0,:,:],wrfdata_fullfile.variables['XLAT'][0,:,:],inverse=False)
	m1.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im = m1.imshow(run2_o3mean[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
	for i in range(len(loclist)):
		print(coordlistlatslons[i][1],coordlistlatslons[i][0],stationdata[loclist[i]+'_o3'].mean())
		m1.scatter(coordlistlatslons[i][1],coordlistlatslons[i][0],c=stationdata[loclist[i]+'_o3'].mean(),s=85,marker='o',cmap='jet',zorder=1003,latlon=True,vmin=0,vmax=50.,linewidth=1.5)
	axes[0].set_title("NUDGED",fontsize=20)		
	
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1])
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im2 = m2.imshow((run3_o3mean[:,:]-run2_o3mean[:,:])/run3_o3mean[:,]*100.,cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-35.,vmax=35.)
	axes[1].set_title("DIFF",fontsize=20)

	cax = fig.add_axes([axes[1].get_position().x1+0.008,axes[1].get_position().y0+0.005,0.02,axes[1].get_position().y1-axes[1].get_position().y0-0.015])
	cbar = fig.colorbar(im2,cax=cax,ticks=[-35.,-17.5,0.,17.5,35.],extend='both')
	cbar.set_label("Surface ozone increase [%]",fontsize=18)
	cbar.ax.tick_params(labelsize=16)
	plt.savefig('Figures/DEPVELplots/MEAN_OZONE_final_4panel_newcamsdata_new_2panel_DIFF',dpi=500)
	plt.show()
	
	
	'''
	for day in range(0,int(np.trunc(wrfdata_fullfile.variables['o3'][:,:,:,:].shape[0]/24.)),1):
		seaice_wrf = np.mean(wrfdata_fullfile.variables['SEAICE'][(day*24):(day*24)+23,:,:],axis=0) #get mean of day
		run1 = wrfdata_base[(day*24):(day*24)+23,:,:]
		histrun1 = run1
		run1 = np.nanmean(run1,axis=0)	#get mean of day
		run2 = wrfdata_nudged[(day*24):(day*24)+23,:,:]
		histrun2 = run2
		run2 = np.nanmean(run2,axis=0)	#get mean of day
		run3 = wrfo3data_nudgedBL_coareg[(day*24):(day*24)+23,:,:]
		histrun3 = run3
		run3 = np.nanmean(run3,axis=0)	#get mean of day
		CAMSo3 = camsdata[(day*6):(day*6)+5,:,:]
	
		fig,axes = plt.subplots(2,2,figsize=(17,15), gridspec_kw = {'wspace':0.24, 'hspace':0.02})
		
		m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
		m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		lons,lats = m1(wrfdata_fullfile.variables['XLONG'][0,:,:],wrfdata_fullfile.variables['XLAT'][0,:,:],inverse=False)
		m1.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
		#m1.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		im = m1.imshow(run1[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
		
		m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
		m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m2.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
		#m2.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		m2.imshow(run3[:,:],cmap='jet',interpolation='none',zorder=1000,vmin=0,vmax=50.)
		
		m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
		m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
		m3.contour(lons,lats,wrfdata_fullfile.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
		#m3.fillcontinents(color='peru',lake_color='cornflowerblue',zorder=1001,alpha=1)
		m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
		m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
		im3 = m3.imshow((run3[:,:]-run1[:,:])/run1[:,:]*100.,cmap='RdBu_r',interpolation='none',zorder=1000,vmin=-100.,vmax=100.)
		
		ax=axes[1,1]
		ax.hist(histrun3[:,:,:].ravel(),bins=50,range=(0.,50.),normed=True,cumulative=False,alpha=0.5,color='green')
		ax.hist(histrun1[:,:,:].ravel(),bins=50,range=(0.,50.),normed=True,cumulative=False,alpha=0.5,color='black')
		ax.yaxis.set_label_position("right")
		ax.yaxis.tick_right()
		ax.set_xlim([0,60])
		ax.set_ylim([0,0.1])
		ax.set_xlabel('Surface ozone mixing ratio [ppb]',fontsize=13)
		ax.set_ylabel('Frequency [-]',fontsize=13)		
		
		cax = fig.add_axes([axes[0,0].get_position().x1+0.01,axes[0,0].get_position().y0+0.01,0.02,axes[0,0].get_position().y1-axes[0,0].get_position().y0-0.02])
		cax2 = fig.add_axes([axes[1,0].get_position().x1+0.01,axes[1,0].get_position().y0+0.01,0.02,axes[1,0].get_position().y1-axes[1,0].get_position().y0-0.02])
		cbar = fig.colorbar(im,cax=cax,ticks=[0,10,20,30,40,50],extend='max')
		cbar.set_label("Daily averaged surface ozone mixing ratio [ppb]",fontsize=13)
		cbar2 = fig.colorbar(im3,cax=cax2,extend='max',ticks=[-100,-75,-50,-25,0,25,50,75,100])
		cbar2.set_label("Daily averaged surface ozone mixing ratio difference [%]",fontsize=13)		
		plt.savefig('Figures/DEPVELplots/OZONE_coaregbase'+str(day),dpi=150)
		#plt.show()
		plt.close()
	
	print(np.mean(histrun1[:,:,:]),np.mean(histrun2[:,:,:]),np.mean(histrun3[:,:,:]))
	fig,ax=plt.subplots(1,1,figsize=(10,10))
	density = gaussian_kde(histrun1[:,:,:].ravel())
	density.covariance_factor = lambda : .25
	density._compute_covariance()
	density2 = gaussian_kde(histrun2[:,:,:].ravel())
	density2.covariance_factor = lambda : .25
	density2._compute_covariance()
	density3 = gaussian_kde(histrun3[:,:,:].ravel())
	density3.covariance_factor = lambda : .25
	density3._compute_covariance()
	densitycams = gaussian_kde(camsdata[:,:,:].ravel())
	densitycams.covariance_factor = lambda : .25
	densitycams._compute_covariance()
	xs = np.linspace(0,50,100)
	ax.plot(xs,density(xs),color='red',label='Base')
	ax.plot(xs,density2(xs),color='orange',label='Nudged')
	ax.plot(xs,density3(xs),color='cyan',label='NudgedCOAREG')
	ax.plot(xs,densitycams(xs),color='blue',label='MACC')
	ax.yaxis.set_label_position("right")
	ax.yaxis.tick_right()
	ax.set_xlim([0,50])
	ax.set_ylim([0,0.1])
	ax.set_xlabel('Surface ozone mixing ratio [ppb]',fontsize=13)
	ax.set_ylabel('Frequency [-]',fontsize=13)
	ax.legend(loc='upper right')
	plt.savefig('Figures/DEPVELplots/OZONEtotalhist_final.png',dpi=150)
	plt.show()
	'''
	
def budgets():
	#BASE
	bf1 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-10_00:00:00','r')
	bf2 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-13_01:00:00','r')
	bf3 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-21_01:00:00','r')
	bf4 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chemdt15_d01_2008-08-31_01:00:00','r')
	wrfdata_base = np.concatenate((bf1['o3'][0:73,0,:,:]*1000.,bf2['o3'][0:192,0,:,:]*1000.,bf3['o3'][0:240,0,:,:]*1000.,bf4['o3'][:,0,:,:]*1000.),axis=0)
	wrfdata_depvel_base = np.concatenate((bf1['DEP_VEL'][0:73,:,:],bf2['DEP_VEL'][0:192,:,:],bf3['DEP_VEL'][0:240,:,:],bf4['DEP_VEL'][:,:,:]),axis=0)
	#Nudged
	wrfdata_nudgedBL_nofixeddep = nc.Dataset('//archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chem_nudgedBL_nofixeddep_d01_2008-08-10_00:00:00','r')
	wrfdata_nudged = wrfdata_nudgedBL_nofixeddep.variables['o3'][:,0,:,:]*1000.
	wrfdata_depvel_nudged = wrfdata_nudgedBL_nofixeddep.variables['DEP_VEL'][:,:,:]
	#NudgedCOAREG
	df1 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chem_nudgedBL_coareg_sherweniodide_d01_2008-08-10_00:00:00','r')
	df2 = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chem_nudgedBL_coareg_sherweniodide_d01_2008-08-11_01:00:00','r')
	wrfo3data_nudgedBL_coareg = np.concatenate((df1['o3'][0:25,0,:,:]*1000.,df2['o3'][:,0,:,:]*1000.),axis=0)
	wrfdata_depvel_coareg = np.concatenate((df1['DEP_VEL'][0:25,:,:],df2['DEP_VEL'][:,:,:]),axis=0)
	wrfdata_vtc_coareg = np.concatenate((df1['VTC'][0:25,:,:],df2['VTC'][:,:,:]),axis=0)
	wrfdata_seaice_coareg = np.concatenate((df1['SEAICE'][0:25,:,:],df2['SEAICE'][:,:,:]),axis=0)
	
	wrfdata_fixeddep = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_fixeddep_d01_2008-08-10_00:00:00','r')
	wrfdata_depvel_fixeddep = wrfdata_fixeddep.variables['DEP_VEL'][:,:,:]
	wrfo3data_fixeddep = wrfdata_fixeddep.variables['o3'][:,0,:,:]*1000.	
	
	#seaicezerotoonematrix = wrfdata_seaice_coareg
	#print(seaicezerotoonematrix)
	#print(seaicezerotoonematrix[1,:,125])
	#print(seaicezerotoonematrix[1,125,:])
	#seaicezerotoonematrix[seaicezerotoonematrix == 0.] = 1.
	
	#print(seaicezerotoonematrix)
	#print(seaicezerotoonematrix[1,:,125])
	#print(seaicezerotoonematrix[1,125,:])
		
	#Deposition in m/s
	run1_depo3 = wrfdata_depvel_base[:,:,:]
	run2_depo3 = wrfdata_depvel_nudged[:,:,:]
	run3_depo3 = wrfdata_depvel_coareg[:,:,:]+wrfdata_vtc_coareg[:,:,:]
	fixeddep_depo3 = wrfdata_depvel_fixeddep[:,:,:]
	
	run1_depo3[run1_depo3 < 3.e-4] = 3.e-4
	run2_depo3[run2_depo3 < 3.e-4] = 3.e-4
	run3_depo3[run3_depo3 < 9.95e-5] = 9.95e-5
	fixeddep_depo3[fixeddep_depo3 < 3.e-4] = 3.e-4

	#run1_depo3 = wrfdata_depvel_base[:,:,:]*(1./seaicezerotoonematrix)
	#run2_depo3 = wrfdata_depvel_nudged[:,:,:]*(1./seaicezerotoonematrix)
	#run3_depo3 = wrfdata_depvel_coareg[:,:,:]*(1./seaicezerotoonematrix)+wrfdata_vtc_coareg[:,:,:]


	#print(1./seaicezerotoonematrix)
	#print(1./seaicezerotoonematrix[1,:,125])
	#print(1./seaicezerotoonematrix[1,125,:])


	#Ozone in ppb
	run1_o3 = wrfdata_base[:,:,:]
	run2_o3 = wrfdata_nudged[:,:,:]
	run3_o3 = wrfo3data_nudgedBL_coareg[:,:,:]
	fixeddep_o3 = wrfo3data_fixeddep[:,:,:]
		
	#Temperatures in K
	run1_t = np.concatenate((bf1['T2'][0:73,:,:],bf2['T2'][0:192,:,:],bf3['T2'][0:240,:,:],bf4['T2'][:,:,:]),axis=0)
	run2_t = wrfdata_nudgedBL_nofixeddep.variables['T2'][:,:,:]
	run3_t = np.concatenate((df1['T2'][0:25,:,:],df2['T2'][:,:,:]),axis=0)
	fixeddep_t = wrfdata_fixeddep.variables['T2'][:,:,:]
	
	
	#Pressure in Pa
	run1_p = np.concatenate((bf1['PSFC'][0:73,:,:],bf2['PSFC'][0:192,:,:],bf3['PSFC'][0:240,:,:],bf4['PSFC'][:,:,:]),axis=0)
	run2_p = wrfdata_nudgedBL_nofixeddep.variables['PSFC'][:,:,:]
	run3_p = np.concatenate((df1['PSFC'][0:25,:,:],df2['PSFC'][:,:,:]),axis=0)
	fixeddep_p = wrfdata_fixeddep.variables['PSFC'][:,:,:]
	
	#Calculate deposition fluxes from piston velocity, mixing ratio, temp and pressure
	dep_run1,deppersec_run1 = calcdepbudget(run1_depo3,run1_o3,run1_t,run1_p)
	dep_run2,deppersec_run2 = calcdepbudget(run2_depo3,run2_o3,run2_t,run2_p)
	dep_run3,deppersec_run3 = calcdepbudget(run3_depo3,run3_o3,run3_t,run3_p)
	dep_fixeddep,deppersec_fixeddep = calcdepbudget(fixeddep_depo3,fixeddep_o3,fixeddep_t,fixeddep_p)
	
	#Recalculate and mean deposition fluxes
	meanflux_run1 = np.mean(dep_run1,axis=0)*1e9 #kg O3 m-2 s-1 to ug O3 m-2 s-1
	meanflux_run2 = np.mean(dep_run2,axis=0)*1e9 #kg O3 m-2 s-1 to ug O3 m-2 s-1
	meanflux_run3 = np.mean(dep_run3,axis=0)*1e9 #kg O3 m-2 s-1 to ug O3 m-2 s-1
	meanflux_fixeddep = np.mean(fixeddep_depo3,axis=0)*1e9 #kg O3 m-2 s-1 to ug O3 m-2 s-1
			
	#Plot mean deposition flux
	fig,axes = plt.subplots(2,2,figsize=(17,15), gridspec_kw = {'wspace':0.24, 'hspace':0.02})
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(wrfdata_nudgedBL_nofixeddep.variables['XLONG'][0,:,:],wrfdata_nudgedBL_nofixeddep.variables['XLAT'][0,:,:],inverse=False)
	m1.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#im = m1.imshow(meanflux_run1,cmap='jet',interpolation='none',zorder=1000, vmin=0, vmax=0.5)
	im = m1.imshow(meanflux_run1,cmap='jet',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.001,vmax=0.5))
		
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#m2.imshow(meanflux_run2,cmap='jet',interpolation='none',zorder=1000, vmin=0, vmax=0.5)
	m2.imshow(meanflux_run2,cmap='jet',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.001,vmax=0.5))
	
	m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
	m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m3.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#m3.imshow(meanflux_run3,cmap='jet',interpolation='none',zorder=1000, vmin=0, vmax=0.5)
	m3.imshow(meanflux_run3,cmap='jet',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.001,vmax=0.5))
	
	ax=axes[1,1]
	ax.hist(meanflux_run1[:,:].ravel(),bins=50,range=(0.,0.5),normed=True,cumulative=False,alpha=0.5,color='Red',label='Base')
	ax.hist(meanflux_run2[:,:].ravel(),bins=50,range=(0.,0.5),normed=True,cumulative=False,alpha=0.5,color='Orange',label='Nudged')
	ax.hist(meanflux_run3[:,:].ravel(),bins=50,range=(0.,0.5),normed=True,cumulative=False,alpha=0.5,color='Cyan',label='NudgedCOAREG')
	ax.yaxis.set_label_position("right")
	ax.yaxis.tick_right()
	ax.set_xlim([0,0.5])
	ax.set_ylim([0,100.])
	ax.legend(loc=1)
	ax.set_xlabel('Mean O$_{3}$ deposition flux [ug m$^{-2}$ s$^{-1}$]',fontsize=13)
	ax.set_ylabel('Frequency [-]',fontsize=13)	
		
	cax = fig.add_axes([axes[0,0].get_position().x1+0.01,axes[0,0].get_position().y0+0.01,0.02,axes[0,0].get_position().y1-axes[0,0].get_position().y0-0.02])
	cbar = fig.colorbar(im,cax=cax,extend='max')
	cbar.set_label('Mean O$_{3}$ deposition flux [ug m$^{-2}$ s$^{-1}$]',fontsize=13)
	plt.savefig('Figures/DEPVELplots/DEPFLUX',dpi=150)
	plt.show()
	
	
	print('WE ARE HERE')
	
	#Recalculate and mean deposition fluxes
	meandepvel_run1 = np.mean(run1_depo3,axis=0)*100. #mean dep_vel in cm/s
	meandepvel_run2 = np.mean(run2_depo3,axis=0)*100. #mean dep_vel in cm/s
	meandepvel_run3 = np.mean(run3_depo3,axis=0)*100. #mean dep_vel in cm/s
	stddepvel_run1 = np.std(run1_depo3,axis=0)*100. #std dep_vel in cm/s
	stddepvel_run2 = np.std(run2_depo3,axis=0)*100. #std dep_vel in cm/s
	stddepvel_run3 = np.std(run3_depo3,axis=0)*100. #std dep_vel in cm/s
	meandepvel_fixeddep = np.mean(fixeddep_depo3,axis=0)*100.
	stddepvel_fixeddep = np.std(fixeddep_depo3,axis=0)*100.
	
	np.save('/home/WUR/barte035/WRFChem/Python/Arrays/longitudes.npy',wrfdata_nudgedBL_nofixeddep.variables['XLONG'][0,:,:])
	np.save('/home/WUR/barte035/WRFChem/Python/Arrays/latitudes.npy',wrfdata_nudgedBL_nofixeddep.variables['XLAT'][0,:,:])
	np.save('/home/WUR/barte035/WRFChem/Python/Arrays/seaice.npy',wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:])
	np.save('/home/WUR/barte035/WRFChem/Python/Arrays/depvel_run1_new.npy',run1_depo3)
	np.save('/home/WUR/barte035/WRFChem/Python/Arrays/depvel_run2_new.npy',run2_depo3)
	np.save('/home/WUR/barte035/WRFChem/Python/Arrays/depvel_run3_new.npy',run3_depo3)
	np.save('/home/WUR/barte035/WRFChem/Python/Arrays/depvel_fixeddep_new.npy',fixeddep_depo3)

	print('WE ARE HERE')

	#Plot mean deposition velocities and std
	fig,axes = plt.subplots(1,2,figsize=(17,7), gridspec_kw = {'wspace':0.24, 'hspace':0.02})

	print('WE ARE HERE')
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(wrfdata_nudgedBL_nofixeddep.variables['XLONG'][0,:,:],wrfdata_nudgedBL_nofixeddep.variables['XLAT'][0,:,:],inverse=False)
	m1.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#im = m1.imshow(meanflux_run1,cmap='jet',interpolation='none',zorder=1000, vmin=0, vmax=0.5)
	im = m1.imshow(meandepvel_run1,cmap='jet',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.008,vmax=0.7))
		
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1])
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#m2.imshow(meanflux_run2,cmap='jet',interpolation='none',zorder=1000, vmin=0, vmax=0.5)
	m2.imshow(meandepvel_run3,cmap='jet',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.008,vmax=0.7))
		
	print('WE ARE HERE')

	cax = fig.add_axes([axes[0].get_position().x1+0.01,axes[0].get_position().y0+0.01,0.02,axes[0].get_position().y1-axes[0].get_position().y0-0.02])
	cbar = fig.colorbar(im,cax=cax,extend='max')
	cbar.set_label('Mean deposition velocity [cm s$^{-1}$]',fontsize=13)
	plt.savefig('Figures/DEPVELplots/DEPVELMEAN',dpi=150)
	plt.show()
	
	#Plot mean deposition velocities and total flux
	fig,axes = plt.subplots(2,2,figsize=(15,13), gridspec_kw = {'wspace':0.25, 'hspace':0.008})

	print('WE ARE HERE')
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(wrfdata_nudgedBL_nofixeddep.variables['XLONG'][0,:,:],wrfdata_nudgedBL_nofixeddep.variables['XLAT'][0,:,:],inverse=False)
	m1.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#im = m1.imshow(meanflux_run1,cmap='jet',interpolation='none',zorder=1000, vmin=0, vmax=0.5)
	im = m1.imshow(meandepvel_run1,cmap='viridis',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.009,vmax=0.5))
		
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0,1])
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#m2.imshow(meanflux_run2,cmap='jet',interpolation='none',zorder=1000, vmin=0, vmax=0.5)
	m2.imshow(meandepvel_run3,cmap='viridis',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.009,vmax=0.5))
	
	m3 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,0])
	m3.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(wrfdata_nudgedBL_nofixeddep.variables['XLONG'][0,:,:],wrfdata_nudgedBL_nofixeddep.variables['XLAT'][0,:,:],inverse=False)
	m3.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m3.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m3.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im3 = m3.imshow(meanflux_run1,cmap='cool',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.004,vmax=0.3))
		
	m4 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1,1])
	m4.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m4.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m4.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m4.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	m4.imshow(meanflux_run3,cmap='cool',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.004,vmax=0.3))
		
	print('WE ARE HERE')

	cax = fig.add_axes([axes[0,0].get_position().x1+0.01,axes[0,0].get_position().y0,0.02,axes[0,0].get_position().y1-axes[0,0].get_position().y0])
	cbar = fig.colorbar(im,cax=cax,extend='max',ticks=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5])
	cbar.ax.tick_params(labelsize=15,width=2)
	cbar.set_label('Mean deposition velocity [cm s$^{-1}$]',fontsize=16)
	cax3 = fig.add_axes([axes[1,0].get_position().x1+0.01,axes[1,0].get_position().y0,0.02,axes[1,0].get_position().y1-axes[1,0].get_position().y0])
	cbar3 = fig.colorbar(im3,cax=cax3,extend='max',ticks=[0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3])
	cbar3.ax.tick_params(labelsize=15,width=2)
	cbar3.set_label('Mean O$_{3}$ deposition flux [$\mu$g m$^{-2}$ s$^{-1}$]',fontsize=16)
	xtext,ytext = m1(-133.,48.)
	axes[0,0].set_title("DEFAULT")
	axes[0,1].set_title("COAREG")
	axes[0,0].text(xtext,ytext,s='(a)',fontsize=14,zorder=1010)
	axes[0,1].text(xtext,ytext,s='(b)',fontsize=14,zorder=1010)
	axes[1,0].text(xtext,ytext,s='(c)',fontsize=14,zorder=1010)
	axes[1,1].text(xtext,ytext,s='(d)',fontsize=14,zorder=1010)
	plt.savefig('Figures/DEPVELplots/DEPVELMEANANDFLUX',dpi=150)
	plt.show()
	
	
	fig,axes = plt.subplots(1,2,figsize=(17,7), gridspec_kw = {'wspace':0.24, 'hspace':0.02})

	print('WE ARE HERE')
		
	m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0])
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(wrfdata_nudgedBL_nofixeddep.variables['XLONG'][0,:,:],wrfdata_nudgedBL_nofixeddep.variables['XLAT'][0,:,:],inverse=False)
	m1.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#im = m1.imshow(meanflux_run1,cmap='jet',interpolation='none',zorder=1000, vmin=0, vmax=0.5)
	im = m1.imshow(stddepvel_run1,cmap='jet',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.005,vmax=0.1))
		
	m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1])
	m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m2.contour(lons,lats,wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	#m2.imshow(meanflux_run2,cmap='jet',interpolation='none',zorder=1000, vmin=0, vmax=0.5)
	m2.imshow(stddepvel_run3,cmap='jet',interpolation='none',zorder=1000, norm=colors.LogNorm(vmin=0.0001,vmax=0.1))
		
	print('WE ARE HERE')

	cax = fig.add_axes([axes[0].get_position().x1+0.01,axes[0].get_position().y0+0.01,0.02,axes[0].get_position().y1-axes[0].get_position().y0-0.02])
	cbar = fig.colorbar(im,cax=cax,extend='max')
	cbar.set_label('Standard deviation of deposition velocity [cm s$^{-1}$]',fontsize=13)
	plt.savefig('Figures/DEPVELplots/DEPVELSTD',dpi=150)
	plt.show()

	
	#Landuse masks
	lutyp = wrfdata_nudgedBL_nofixeddep.variables['ISLTYP'][:,:,:]
	seaice = wrfdata_nudgedBL_nofixeddep.variables['SEAICE'][:,:,:]
	truearr = np.ones((lutyp.shape[0], lutyp.shape[1], lutyp.shape[2]), dtype=bool)
	falsearr = np.zeros((lutyp.shape[0], lutyp.shape[1], lutyp.shape[2]), dtype=bool)
	ocean_mask = np.where(lutyp == 14, True, False)
	land_mask = np.where(lutyp != 14, np.where(lutyp != 16,True,False), False)
	landice_mask = np.where(lutyp == 16, np.where(seaice == 0,True,False), False)
	seaice_mask = np.where(lutyp == 16, np.where(seaice > 0, True,False),False)
	#Calculate percentage of gridcells
	print('Ocean',(np.count_nonzero(ocean_mask)/672.)/(249*249)*100,'%')
	print('Land',(np.count_nonzero(land_mask)/672.)/(249*249)*100,'%')
	print('Land ice',(np.count_nonzero(landice_mask)/672.)/(249*249)*100,'%')
	print('Sea ice',(np.count_nonzero(seaice_mask)/672.)/(249*249)*100,'%')
	
	a0,a1,a2 = maskarray(dep_run1,ocean_mask,'run1','ocean',run1_depo3)
	a0,a1,a2 = maskarray(dep_run1,land_mask,'run1','land',run1_depo3)
        a0,a1,a2 = maskarray(dep_run1,landice_mask,'run1','landice',run1_depo3)
	a0,a1,a2 = maskarray(dep_run1,seaice_mask,'run1','seaice',run1_depo3)
	a0,a1,a2 = maskarray(dep_run2,ocean_mask,'run2','ocean',run2_depo3)
	a0,a1,a2 = maskarray(dep_run2,land_mask,'run2','land',run2_depo3)
	a0,a1,a2 = maskarray(dep_run2,landice_mask,'run2','landice',run2_depo3)
	a0,a1,a2 = maskarray(dep_run2,seaice_mask,'run2','seaice',run2_depo3)
	a0,a1,a2 = maskarray(dep_run3,ocean_mask,'run3','ocean',run3_depo3)
	a0,a1,a2 = maskarray(dep_run3,land_mask,'run3','land',run3_depo3)
	a0,a1,a2 = maskarray(dep_run3,landice_mask,'run3','landice',run3_depo3)
	a0,a1,a2 = maskarray(dep_run3,seaice_mask,'run3','seaice',run3_depo3)
	a0,a1,a2 = maskarray(dep_fixeddep,ocean_mask,'fixeddep','ocean',fixeddep_depo3)
	a0,a1,a2 = maskarray(dep_fixeddep,land_mask,'fixeddep','land',fixeddep_depo3)
	a0,a1,a2 = maskarray(dep_fixeddep,landice_mask,'fixeddep','landice',fixeddep_depo3)
	a0,a1,a2 = maskarray(dep_fixeddep,seaice_mask,'fixeddep','seaice',fixeddep_depo3)	
		
def maskarray(depvalues,mask,runname,surfacename,pistvel):
	print('RUN: ',runname,', Surface type: ',surfacename)
	maskedarr = depvalues[1:,:,:][mask[1:,:,:]==True]
	maskedarrpistvel = pistvel[1:,:,:][mask[1:,:,:]==True]
	meandep = np.nanmean(maskedarr,axis=0) #kg O3 m-2 s-1
	stddep = np.nanstd(maskedarr,axis=0)
	print('Mean deposition in ug O3 m-2 s-1: ',meandep*1e9)
	print('Std deposition in ug O3 m-2 s-1: ',stddep*1e9)
	print('Mean deposition in nmol m-2 s-1: ',meandep*(1./48.e-12)) #kg O3 m-2 s-1 to nmol m-2 s-1 (* nmol/kg) 48 g/mol = 48e-3 kg/mol = 48e-12 kg/nmol
	sumarr = np.nansum(maskedarr*(30000.*30000.)*(60.*60.),axis=0) #kg O3 m-2 s-1 to kg O3 simulation-1
	print('Deposition in kg O3 over simulation: ',sumarr)
	sumarr2 = sumarr*1.e-9/671.*24.*365.
	print('Deposition in Tg O3 yr-1: ',sumarr2)
	meanpistvel = np.nanmean(maskedarrpistvel,axis=0)
	stdpistvel = np.nanstd(maskedarrpistvel,axis=0)
	print('Deposition velocity mean + std [cm s-1]: ',meanpistvel*100.,stdpistvel*100.)
	return(meandep,sumarr,sumarr2)
	
def calcdepbudget(pistonvelocity,surfaceozone,temperature,pressure):
	#Calculate deposition
	mwdry = 28.9647
	mwozo = 48.
	dx = 30000.
	sectohour = 60.*60.

	#rho = P/(Rgas*T)
	#kg/m3 = Pa/(J kg-1 K-1 * K)
	rho = pressure[:,:,:]/(287.058*temperature[:,:,:])
	#fac = 1e-6 * rho * 1/mwdry (mwdry = g/mol)
	#(mol/m3)/ppmv = kg/m3 * mol/kg
	fac = 1.e-6 * rho * (1./(mwdry/1.e-3))
	#dep = piston velocity * mixing ratio * factor * o3_mw
	#kg O3 m-2 s-1 = m s-1 * mol/mol * mol/m3 * kg/mol
	dep = pistonvelocity[:,:,:] * (surfaceozone[:,:,:]/1000. * fac) * mwozo/1.e-3 #kg O3 m-2 s-1
	#deptotal = dep * dx * dx * hourtosec
	#kg O3 s-1 = kg O3 m-2 s-1 * m * m
	deppersec = dep * (dx * dx)
	'''
	#airmas=-1.0*(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/grav
	#kg = Pa * m * m / m s-2
	airmas = -1.*pressure[:,:,:]*30000.*30000./9.81
	#depppbhour = dep*dx*dx*sectohour / airmas * mwdry / o3_mw)*1e6
	#ppb hr-1 or (mmol/mol) hr-1 = kg O3 m-2 s-1 * m * m * s h-1 / kg * g mol-1 / g mol-1
	depppbhour = (dep*30000.*30000.*sectohour) / airmas * mwdry / mwozo)*1e6
	'''	
	return(dep,deppersec)
	
#dep_vel()
ozone()
#budgets()
