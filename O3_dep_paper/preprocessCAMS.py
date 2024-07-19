import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
from sklearn.metrics import mean_squared_error
from math import sqrt
from scipy.stats import gaussian_kde
import time
from datetime import datetime
from dateutil.relativedelta import relativedelta
import dateutil.parser

data = nc.Dataset('/home/WUR/barte035/WRFChem/o3_analysis_DATA/CAMS-MACC/MACC_08-092008.nc','r')
camso3all = data.variables['go3'][:,:,:,:]
camstimearr2 = data.variables['time']
time_converted = nc.num2date(camstimearr2[:], camstimearr2.units, camstimearr2.calendar)

datanew = nc.Dataset('/home/WUR/barte035/WRFChem/o3_analysis_DATA/CAMS-MACC/levtype_ml.nc','r')
camso3allnew = datanew.variables['go3'][:-159,0:61,:]
camstimearr2new = datanew.variables['time']
time_convertednew = nc.num2date(camstimearr2new[:-159], camstimearr2new.units, camstimearr2new.calendar)
temp = np.zeros((camso3allnew.shape[0],camso3allnew.shape[1],camso3allnew.shape[2]))
temp[:,:,0:240] = camso3allnew[:,:,240:480]
temp[:,:,240:480] = camso3allnew[:,:,0:240]    
camso3allnew[:,:,:]=temp[:,:,:]

print(time_converted)
print(time_convertednew)

print(camso3all.shape,camso3allnew.shape)

plottime = 35

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
plt.figure(figsize=(15,10))

m = Basemap(width=(dx_dom*(we-10))-1,height=(dx_dom*(sn-10))-1,resolution='l',area_thresh=1000.,projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0.)
cams_trans = m.transform_scalar(camso3all[plottime,-1,::-1,:],lons=np.arange(-180.,180.,0.75),lats=np.arange(45.,90.,0.75),nx=249,ny=249,order=0)
m.imshow(cams_trans,vmin=0.,vmax=1E-7,cmap='jet')
m.colorbar(label='Surface ozone [kg kg-1]',extend='max')
m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
plt.savefig('TEMPO3PLOT')
plt.show()
m = Basemap(width=(dx_dom*(we-10))-1,height=(dx_dom*(sn-10))-1,resolution='l',area_thresh=1000.,projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0.)
cams_trans_ppb = m.transform_scalar((camso3all[plottime,-1,::-1,:]*(29.91e-3/48e-3)*1e9),lons=np.arange(-180.,180.,0.75),lats=np.arange(45.,90.,0.75),nx=249,ny=249,order=0)
m.imshow(cams_trans_ppb,vmin=0.,vmax=60,cmap='jet')
m.colorbar(label='Surface ozone [ppb]',extend='max')
m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
plt.savefig('TEMPO3PLOT_ppb')
plt.show()
m = Basemap(width=(dx_dom*(we-10))-1,height=(dx_dom*(sn-10))-1,resolution='l',area_thresh=1000.,projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0.)
cams_trans_ppb_new = m.transform_scalar((camso3allnew[plottime,::-1,:]*(28.97e-3/48e-3)*1e9),lons=np.arange(-180.,180.,0.75),lats=np.arange(45.,90.,0.75),nx=249,ny=249,order=0)
m.imshow(cams_trans_ppb_new,vmin=0.,vmax=60,cmap='jet')
m.colorbar(label='Surface ozone [ppb]',extend='max')
m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
plt.savefig('TEMPO3PLOTnew_ppb')
plt.show()

m = Basemap(width=(dx_dom*(we-10))-1,height=(dx_dom*(sn-10))-1,resolution='l',area_thresh=1000.,projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0.)
cams_trans_ppb = m.transform_scalar(np.mean((camso3all[:,-1,::-1,:]*(29.91e-3/48e-3)*1e9),axis=0),lons=np.arange(-180.,180.,0.75),lats=np.arange(45.,90.,0.75),nx=249,ny=249,order=0)
m.imshow(cams_trans_ppb,vmin=0.,vmax=50,cmap='jet')
m.colorbar(label='Surface ozone [ppb]',extend='max')
m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
plt.savefig('CAMSPREPROCESSMEANO3_ppb')
plt.show()

m = Basemap(width=(dx_dom*(we-10))-1,height=(dx_dom*(sn-10))-1,resolution='l',area_thresh=1000.,projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0.)
cams_trans_ppb_new = m.transform_scalar(np.mean((camso3allnew[:,::-1,:]*(28.97e-3/48e-3)*1e9),axis=0),lons=np.arange(-180.,180.,0.75),lats=np.arange(45.,90.,0.75),nx=249,ny=249,order=0)
m.imshow(cams_trans_ppb_new,vmin=0.,vmax=50,cmap='jet')
m.colorbar(label='Surface ozone [ppb]',extend='max')
m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
plt.savefig('CAMSPREPROCESSMEANO3new_ppb')
plt.show()



cams_trans_series = []
cams_trans_series_new = []


for i in range(camso3all.shape[0]):
	m = Basemap(width=(dx_dom*(we-10))-1,height=(dx_dom*(sn-10))-1,resolution='l',area_thresh=1000.,projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0.)
	cams_trans = m.transform_scalar(camso3all[i,-1,::-1,:],lons=np.arange(-180.,180.,0.75),lats=np.arange(45.,90.,0.75),nx=249,ny=249,order=0)
	camsppb = cams_trans*(29.91e-3/48e-3)*1e9 #kg/kg to ppb
	cams_trans_series.append(camsppb)
	print(i,camso3all.shape[0])
	
for i in range(camso3allnew.shape[0]):
	m = Basemap(width=(dx_dom*(we-10))-1,height=(dx_dom*(sn-10))-1,resolution='l',area_thresh=1000.,projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0.)
	cams_trans_new = m.transform_scalar(camso3allnew[i,::-1,:],lons=np.arange(-180.,180.,0.75),lats=np.arange(45.,90.,0.75),nx=249,ny=249,order=0)
	camsppbnew = cams_trans_new*(28.97e-3/48e-3)*1e9 #kg/kg to ppb
	cams_trans_series_new.append(camsppbnew)
	print(i,camso3allnew.shape[0])


cams_trans_series = np.array(cams_trans_series)[35:,:,:][:-16,:,:]
time_converted = time_converted[35:][:-16]
cams_trans_series_new = np.array(cams_trans_series_new)[72:,:,:][:-32,:,:]
time_convertednew = time_convertednew[72:][:-32]

print(time_converted)
print(time_convertednew)


#np.save('/home/WUR/barte035/WRFChem/o3_analysis_DATA/CAMS-MACC/MACCnpo3array',cams_trans_series)
#np.save('/home/WUR/barte035/WRFChem/o3_analysis_DATA/CAMS-MACC/MACCnptimearray',time_converted)
np.save('/home/WUR/barte035/WRFChem/o3_analysis_DATA/CAMS-MACC/MACCnpo3arraynew',cams_trans_series_new)
np.save('/home/WUR/barte035/WRFChem/o3_analysis_DATA/CAMS-MACC/MACCnptimearraynew',time_convertednew)

print(cams_trans_series.shape,cams_trans_series_new.shape,time_converted.shape,time_convertednew.shape)


#print(cams_trans_series)
#print(time_converted)
