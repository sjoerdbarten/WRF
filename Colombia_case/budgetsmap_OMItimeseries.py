from __future__ import unicode_literals

import netCDF4 as nc
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ncfile = nc.Dataset('/archive/ESG/barte035/WRFV3/wrfout_d01_2014-01-01_00:00:00','r')
GFEDpath = '/archive/ESG/barte035/Colombia/GFED/CBMZdaily/'
GFEDfiles = [os.path.join(GFEDpath,filename) for filename in os.listdir(GFEDpath) if filename.startswith('wrffire')]
landmask = ncfile.variables['LANDMASK'][0,:,:]
for j in range(50,99):
	for i in range(0,25):
		landmask[i,j]=1

daysim = 30 #30
timesim = daysim*24
timearray = np.arange(0,timesim+1)
timearrayfire = np.arange(0,timesim,24)
areacell = (27000.*27000.)*1e-6 #km2
N_IC = 250. #moles emitted per lightning strike IC
N_CG = 250. #moles emitted per lightning strike CG

NOAnthro = np.zeros(shape=(timesim+1,1,99,99))
NO2Anthro = np.zeros(shape=(timesim+1,1,99,99))
NOxBio = np.zeros(shape=(timesim+1,1,99,99))
LNOxTimeseriesDomain = np.zeros(shape=(timesim+1,1,99,99))
NOxFire = np.zeros(shape=(daysim,1,99,99))
for i in range(0,timesim):
	NOAnthro[i,:,:,:] = ncfile.variables['E_NO'][i,0,:,:]
	NO2Anthro[i,:,:,:] = ncfile.variables['E_NO2'][i,0,:,:]
	NOxAnthro = NOAnthro+NO2Anthro
	NOxBio[i,:,:,:] = ncfile.variables['EBIO_NO'][i,:,:]
	if i >= 1:
		LNOxTimeseriesDomain[i] = ((ncfile.variables['IC_FLASHCOUNT'][i,:,:]*N_IC/areacell)+(ncfile.variables['CG_FLASHCOUNT'][i,:,:]*N_CG/areacell))-((ncfile.variables['IC_FLASHCOUNT'][i-1,:,:]*N_IC/areacell)+(ncfile.variables['CG_FLASHCOUNT'][i-1,:,:]*N_CG/areacell))
for i in range(0,daysim):
	GFEDfile = nc.Dataset(GFEDfiles[i],'r')
	NOxFire[i,:,:,:] = GFEDfile.variables['ebu_in_no'][0,0,:,:]*24
IC_LNOx = ncfile.variables['IC_FLASHCOUNT'][timesim,:,:]*N_IC/areacell
CG_LNOx = ncfile.variables['CG_FLASHCOUNT'][timesim,:,:]*N_CG/areacell

NOxAnthroTimes = np.sum(NOxAnthro,axis=1)
NOxAnthroTotal = np.sum(NOxAnthroTimes,axis=0)
NOxBioTimes = np.sum(NOxBio,axis=1)
NOxBioTotal = np.sum(NOxBioTimes,axis=0)
NOxFireTimes = np.sum(NOxFire,axis=1)
NOxFireTotal = np.sum(NOxFireTimes,axis=0)
NOxLightningTotal = IC_LNOx + CG_LNOx

NOxALL = NOxAnthroTotal+NOxBioTotal+NOxFireTotal+NOxLightningTotal
NOxALLConvert = ((NOxAnthroTotal+NOxBioTotal+NOxFireTotal+NOxLightningTotal)*14.0067*areacell*31*1e-6)/daysim
NOxALL[NOxALL==0.]=np.nan

NOxAnthroFrac = NOxAnthroTotal/NOxALL*100
NOxBioFrac = NOxBioTotal/NOxALL*100
NOxFireFrac = NOxFireTotal/NOxALL*100
NOxLightningFrac = NOxLightningTotal/NOxALL*100

print 'NOxAnthroFrac=',np.nanmean(NOxAnthroFrac)
print 'NOxBioFrac=',np.nanmean(NOxBioFrac)
print 'NOxFireFrac=',np.nanmean(NOxFireFrac)
print 'NOxLightningFrac=',np.nanmean(NOxLightningFrac)

NOxMax = np.zeros(shape=(99,99))
NOxMaxAnthro = np.zeros(shape=(99,99))
NOxMaxBio = np.zeros(shape=(99,99))
NOxMaxFire = np.zeros(shape=(99,99))
NOxMaxLightning = np.zeros(shape=(99,99))
NOxEmpty = np.zeros(shape=(99,99))
NOxEmpty[NOxEmpty==0.]=np.nan


NOxAnthroTimeseries = np.sum(np.sum(np.sum(NOxAnthro,axis=3),axis=2),axis=1)
NOxBioTimeseries = np.sum(np.sum(np.sum(NOxBio,axis=3),axis=2),axis=1)
NOxFireTimeseries = np.sum(np.sum(np.sum(NOxFire,axis=3),axis=2),axis=1)
LNOxTimeseries = np.sum(np.sum(np.sum(LNOxTimeseriesDomain,axis=3),axis=2),axis=1)

plt.figure(figsize=(15,9))
plt.plot(timearray/24.,NOxAnthroTimeseries*14.0067*1e-6)
plt.plot(timearray/24.,NOxBioTimeseries*14.0067*1e-6)
plt.plot(timearrayfire/24.,NOxFireTimeseries*14.0067*1e-6/24)
plt.plot(timearray/24.,LNOxTimeseries*14.0067*1e-6)
plt.ylabel('N-flux [Mg N km-2 hr-1]',size=12)
plt.xlabel('Time since start of simulation [days]',size=12)
plt.show()

for i in range(0,99):
	for j in range(0,99):
		NOxMax[i,j] = np.maximum.reduce([NOxAnthroFrac[i,j],NOxBioFrac[i,j],NOxFireFrac[i,j],NOxLightningFrac[i,j]])
		if NOxMax[i,j] == NOxAnthroFrac[i,j]:
			NOxMaxAnthro[i,j] = NOxMax[i,j]
			NOxMaxBio[i,j] = np.nan
			NOxMaxFire[i,j] = np.nan
			NOxMaxLightning[i,j] = np.nan
			#NOxMaxIncl[i,j,0] = 444
		if NOxMax[i,j] == NOxBioFrac[i,j]:
			NOxMaxBio[i,j] = NOxMax[i,j]
			NOxMaxAnthro[i,j] = np.nan
			NOxMaxFire[i,j] = np.nan
			NOxMaxLightning[i,j] = np.nan
			#NOxMaxIncl[i,j,0] = 555
		if NOxMax[i,j] == NOxFireFrac[i,j]:
			NOxMaxFire[i,j] = NOxMax[i,j]
			NOxMaxAnthro[i,j] = np.nan
			NOxMaxBio[i,j] = np.nan
			NOxMaxLightning[i,j] = np.nan
			#NOxMaxIncl[i,j,0] = 666
		if NOxMax[i,j] == NOxLightningFrac[i,j]:
			NOxMaxLightning[i,j] = NOxMax[i,j]
			NOxMaxAnthro[i,j] = np.nan
			NOxMaxBio[i,j] = np.nan
			NOxMaxFire[i,j] = np.nan
			#NOxMaxIncl[i,j,0] = 777
#		if NOxALLConvert[i,j] <= 5:
#			NOxMaxAnthro[i,j] = np.nan
#			NOxMaxBio[i,j] = np.nan
#			NOxMaxFire[i,j] = np.nan
#			NOxMaxLightning[i,j] = np.nan
			
			
for i in range(0,4):
	NOxMax[i,:] = np.nan
	NOxMaxAnthro[i,:] = np.nan
	NOxMaxBio[i,:] = np.nan
	NOxMaxFire[i,:] = np.nan
	NOxMaxLightning[i,:] = np.nan
for j in range(0,4):
	NOxMax[:,j] = np.nan
	NOxMaxAnthro[:,j] = np.nan
	NOxMaxBio[:,j] = np.nan
	NOxMaxFire[:,j] = np.nan
	NOxMaxLightning[:,j] = np.nan
for i in range(95,99):
	NOxMax[i,:] = np.nan
	NOxMaxAnthro[i,:] = np.nan
	NOxMaxBio[i,:] = np.nan
	NOxMaxFire[i,:] = np.nan
	NOxMaxLightning[i,:] = np.nan
for j in range(95,99):
	NOxMax[:,j] = np.nan
	NOxMaxAnthro[:,j] = np.nan
	NOxMaxBio[:,j] = np.nan
	NOxMaxFire[:,j] = np.nan
	NOxMaxLightning[:,j] = np.nan

lat = ncfile.variables['XLAT'][timesim,:,:]
lon = ncfile.variables['XLONG'][timesim,:,:]
we = len(ncfile.dimensions['west_east'])
sn = len(ncfile.dimensions['south_north'])
dx = 27000.
width = (we*dx)-1
height = (sn*dx)-1
lat_0 = 4.891087
lon_0 = -71.069204


plt.figure(figsize=(15,9))
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
m.imshow(ncfile.variables['HGT'][0,:,:],cmap='terrain',vmin=0,vmax=4000,interpolation='none')
m.drawcoastlines()
m.drawcountries()
#m.drawparallels(np.arange(-20,30,10),labels=[True,False,True])
#m.drawmeridians(np.arange(-90,-50,20),labels=[True,False,True])
c = plt.colorbar(ticks=[0,500,1000,1500,2000,2500,3000,3500,4000],extend='max')
c.set_label('Terrain height [m]',size=15)
plt.show()

plt.figure(figsize=(15,9))
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
m.imshow(NOxALLConvert,cmap='plasma_r',vmin=0,vmax=100,interpolation='none')
m.drawcoastlines()
m.drawcountries()
#m.drawparallels(np.arange(-20,30,10),labels=[True,False,True])
#m.drawmeridians(np.arange(-90,-50,20),labels=[True,False,True])
c = plt.colorbar(ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='max')
c.set_label(r'NO$_{x}$ Flux [Mg N month$^{-1}$]',size=14)
#plt.text(100000,2.8e6,'(a)',fontsize=9)
plt.show()


plt.figure(figsize=(15,9))
#m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,projection='lcc',lat_0=lat_0,lon_0=lon_0)
#m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,llcrnrlon=-81.684021,llcrnrlat=-7.9757309,urcrnrlon=-60.025818,urcrnrlat=15.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0)
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
m.imshow(np.where((landmask==1),NOxMaxBio,np.nan),cmap='Greens',vmin=0,vmax=100,interpolation='none')
m.imshow(np.where((landmask==1),NOxMaxFire,np.nan),cmap='Oranges',vmin=0,vmax=100,interpolation='none')
m.imshow(np.where((landmask==1),NOxMaxLightning,np.nan),cmap='Blues',vmin=0,vmax=100,interpolation='none')
m.imshow(np.where((landmask==1),NOxMaxAnthro,np.nan),cmap='Reds',vmin=0,vmax=100,interpolation='none')
#c = m.colorbar(ticks=[0,10,20,30,40,50,60,70,80,90,100],location="left",pad='0.1%',size='2%')
#c.set_label('FRACTIONAL CONTRIBUTION',size=18)
#m.imshow(NOxEmpty,cmap='Greys',vmin=0,vmax=100,interpolation='none')
m.drawcoastlines()
m.drawcountries()
#m.drawparallels(np.arange(-20,30,10),labels=[True,False,True])
#m.drawmeridians(np.arange(-90,-50,20),labels=[True,False,True])
#c = plt.colorbar(ticks=[0,10,20,30,40,50,60,70,80,90,100])
#c.set_label('Maximum contribution [%]',size=14)
custom_legend = [Patch(facecolor='red',label='Anthropogenic'),Patch(facecolor='green',label='Biogenic'),Patch(facecolor='orange',label='Biomass burning'),Patch(facecolor='blue',label='Lightning')]
plt.legend(custom_legend,['Anthropogenic','Biogenic','Biomass burning','Lightning'],loc='upper left',prop={'size': 9})
#plt.text(2.5e6,2.8e6,'(b)',fontsize=9)
plt.show()


fig=plt.figure(figsize=(15,9))
ax=fig.add_subplot(211)
#m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,projection='lcc',lat_0=lat_0,lon_0=lon_0)
#m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,llcrnrlon=-81.684021,llcrnrlat=-7.9757309,urcrnrlon=-60.025818,urcrnrlat=15.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0)
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
m.imshow(np.where((landmask==1),NOxMaxBio,np.nan),cmap='Greens',vmin=0,vmax=100,interpolation='none')
m.imshow(np.where((landmask==1),NOxMaxFire,np.nan),cmap='Oranges',vmin=0,vmax=100,interpolation='none')
m.imshow(np.where((landmask==1),NOxMaxLightning,np.nan),cmap='Blues',vmin=0,vmax=100,interpolation='none')
m.imshow(np.where((landmask==1),NOxMaxAnthro,np.nan),cmap='Reds',vmin=0,vmax=100,interpolation='none')
#c = m.colorbar(ticks=[0,10,20,30,40,50,60,70,80,90,100],location="left",pad='0.1%',size='2%')
#c.set_label('FRACTIONAL CONTRIBUTION',size=18)
#m.imshow(NOxEmpty,cmap='Greys',vmin=0,vmax=100,interpolation='none')
m.drawcoastlines()
m.drawcountries()
#m.drawparallels(np.arange(-20,30,10),labels=[True,False,True])
#m.drawmeridians(np.arange(-90,-50,20),labels=[True,False,True])
#c = plt.colorbar(ticks=[0,10,20,30,40,50,60,70,80,90,100])
#c.set_label('Maximum contribution [%]',size=14)
custom_legend = [Patch(facecolor='red',label='Anthropogenic'),Patch(facecolor='green',label='Biogenic'),Patch(facecolor='orange',label='Biomass burning'),Patch(facecolor='blue',label='Lightning')]
plt.legend(custom_legend,['Anthropogenic','Biogenic','Biomass burning','Lightning'],loc='upper left',prop={'size': 7})
plt.text(2.5e6,2.8e6,'(b)',fontsize=9)


ax=fig.add_subplot(221)
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
m.imshow(NOxALLConvert,cmap='plasma_r',vmin=0,vmax=100,interpolation='none')
m.drawcoastlines()
m.drawcountries()
#m.drawparallels(np.arange(-20,30,10),labels=[True,False,True])
#m.drawmeridians(np.arange(-90,-50,20),labels=[True,False,True])
c = plt.colorbar(ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='max')
c.set_label(r'NO$_{x}$ Flux [Mg N month$^{-1}$]',size=14)
plt.text(100000,2.8e6,'(a)',fontsize=9)
plt.show()

print np.count_nonzero(~np.isnan(np.where((landmask==1),NOxMaxBio,np.nan)))
print np.count_nonzero(~np.isnan(np.where((landmask==1),NOxMaxFire,np.nan)))
print np.count_nonzero(~np.isnan(np.where((landmask==1),NOxMaxAnthro,np.nan)))
print np.count_nonzero(~np.isnan(np.where((landmask==1),NOxMaxLightning,np.nan)))

#return(NOxMaxLightning,NOxMaxAnthro,NOxMaxBio,NOxMaxFire)




import numpy as np
from os import listdir
from os.path import isfile, join
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

def read_grd(filename):
    with open(filename) as infile:
        ncols = int(infile.readline().split()[1])
        nrows = int(infile.readline().split()[1])
        xllcorner = float(infile.readline().split()[1])
        yllcorner = float(infile.readline().split()[1])
        cellsize = float(infile.readline().split()[1])
        nodata_value = int(infile.readline().split()[1])
        version = float(infile.readline().split()[1])
    longitude = xllcorner + cellsize * np.arange(ncols)
    latitude = yllcorner + cellsize * np.arange(nrows)
    value = np.loadtxt(filename, skiprows=7)
    print(value.shape)
    #print(longitude[97.625*(1/0.125)],longitude[120*(1/0.125)])
    #print(latitude[85*(1/0.125)],latitude[104.625*(1/0.125)])
    #print(value[85*(1/0.125):104.625*(1/0.125),97.625*(1/0.125):120*(1/0.125)])
    #print(value[85*(1/0.125):104.625*(1/0.125),97.625*(1/0.125):120*(1/0.125)].shape)
    return longitude[97.625*(1./0.125):120*(1./0.125)], latitude[85*(1./0.125):104.625*(1./0.125)], value
    #[85*(1/0.125):104.625*(1/0.125),97.625*(1/0.125):120*(1/0.125)]

#'''
path = '/archive/ESG/barte035/Colombia/OMI/MonthlyMeans/'
files = sorted([f for f in listdir(path) if isfile(join(path,f))])
print(files)
yeararray = ['2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019']
meanarray = np.zeros(len(yeararray))
meanarray_filtered = np.zeros(len(yeararray))
meanarray_filtered_land = np.zeros(len(yeararray))
meanarray_filtered_land_anthro = np.zeros(len(yeararray))
meanarray_filtered_land_bio = np.zeros(len(yeararray))
meanarray_filtered_land_fire = np.zeros(len(yeararray))
meanarray_filtered_land_lightning = np.zeros(len(yeararray))
array_bogota = np.zeros(len(yeararray))
array_caracas = np.zeros(len(yeararray))
geopath = '/archive/ESG/barte035/WPS/geo_em.d01.nc'
geo = nc.Dataset(geopath,'r')
for i in range(len(files)):
	print(files[i])
	longitude,latitude,value = read_grd(path+files[i])
	
	array_bogota[i] = max(value[682,846],value[683,846],value[683,847],value[682,847],value[681,847],value[681,846],value[681,845],value[682,845],value[683,845])
	array_caracas[i] = max(value[635,903],value[636,903],value[636,904],value[635,904],value[634,904],value[634,903],value[634,902],value[635,902],value[636,902])
	
	#Extract WRF-Chem lat/lon and regrid data
	geo = nc.Dataset(geopath,'r')
	we  = geo.variables['XLAT_M'].shape[2]
	sn  = geo.variables['XLAT_M'].shape[1]
	lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
	lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
	dx_dom = 27000
	
	#plt.figure(figsize=(15,10))
    	m = Basemap(width=(we*dx_dom)-1,height=(sn*dx_dom)-1,resolution='l',area_thresh=100.,projection='lcc',lat_0=4.891087,lon_0=-71.069204,lat_1=40.,lat_2=70.,lon_1=-70.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120)
	nx = 99
	ny = 99
	x2_trans = np.linspace(-179.9375,179.9375,2880)
	y2_trans = np.linspace(-89.9375,89.9375,1440)
	
	value = np.where(value < 0, np.nan, value)

	#value_trans = m.transform_scalar(value[::-1,:],lons=np.arange(-82.375,-60.,0.125),lats=np.arange(-5.,14.625,0.125),nx=100.,ny=100.,order=1)
	value_trans = m.transform_scalar(value[::-1],lons=x2_trans,lats=y2_trans,nx=nx,ny=ny,order=0)
	value_trans_filtered = np.where(value_trans < 0, np.nan, value_trans)
	value_trans_filtered_land = np.where(landmask==1, value_trans_filtered, np.nan)
	value_trans_filtered_land_anthro = np.where(~np.isnan(NOxMaxAnthro) > 0., value_trans_filtered_land, np.nan)
	value_trans_filtered_land_bio = np.where(~np.isnan(NOxMaxBio) > 0., value_trans_filtered_land, np.nan)
	value_trans_filtered_land_fire = np.where(~np.isnan(NOxMaxFire) > 0., value_trans_filtered_land, np.nan)
	value_trans_filtered_land_lightning = np.where(~np.isnan(NOxMaxLightning) > 0., value_trans_filtered_land, np.nan)
	#array_bogota[i] = value_trans_filtered[48,34]
	#array_caracas[i] = value_trans_filtered[76,69]
	
	print(value_trans.shape)
	
	#lons=np.arange(-180,180.1,0.1),lats=np.arange(60,90.1,0.1)
	#m.imshow(value_trans)
	#m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	#m.drawparallels(np.arange(-90.,-50,10.),labels=[False,False,False,False],zorder=1001)
	#m.drawmeridians(np.arange(-10,20,5.),labels=[True,False,False,True],zorder=1001)
	#plt.show()
	#plt.close()
	
	meanarray[i] = np.nanmean(value_trans)
	meanarray_filtered[i] = np.nanmean(value_trans_filtered)
	meanarray_filtered_land[i] = np.nanmean(value_trans_filtered_land)
	meanarray_filtered_land_anthro[i] = np.nanmean(value_trans_filtered_land_anthro)
	meanarray_filtered_land_bio[i] = np.nanmean(value_trans_filtered_land_bio)
	meanarray_filtered_land_fire[i] = np.nanmean(value_trans_filtered_land_fire)
	meanarray_filtered_land_lightning[i] = np.nanmean(value_trans_filtered_land_lightning)

print(array_bogota)
print(array_caracas)
mean_bogota = np.mean(array_bogota)/100.
mean_caracas = np.mean(array_caracas)/100.

print((meanarray_filtered_land_anthro/100.)/(meanarray_filtered/100.))
print((meanarray_filtered_land_anthro/100.)/(meanarray_filtered_land/100.))

plt.rcParams['xtick.labelsize']=20
plt.rcParams['ytick.labelsize']=20

plt.figure(figsize=(15,9))
#plt.plot(yeararray,meanarray/100.,label='Non-Filtered')
#plt.plot(yeararray,meanarray_filtered/100.,label='Whole domain',color='black',marker='o')
plt.plot(yeararray,meanarray_filtered_land/100.,label='Domain',color='black',marker='o',linewidth=3)
plt.plot(yeararray,meanarray_filtered_land_anthro/100.,label='Anthropogenic',color='red',marker='o',linewidth=3)
plt.plot(yeararray,meanarray_filtered_land_bio/100.,label='Biogenic',color='green',marker='o',linewidth=3)
plt.plot(yeararray,meanarray_filtered_land_fire/100.,label='Biomass burning',color='orange',marker='o',linewidth=3)
plt.plot(yeararray,meanarray_filtered_land_lightning/100.,label='Lightning',color='blue',marker='o',linewidth=3)
plt.plot(yeararray,array_bogota/100.,label='Bogot\xe1',color='brown',marker='o',linewidth=3)
plt.plot(yeararray,array_caracas/100.,label='Caracas',color='cyan',marker='o',linewidth=3)
plt.plot(yeararray,[mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota],label='Bogot\xe1 mean',color='brown',linestyle='dashed')
plt.plot(yeararray,[mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas],label='Caracas mean',color='cyan',linestyle='dashed')
plt.legend(loc='best',fontsize=14)
plt.xlim([2005,2019])
plt.axvspan(2009.8,2010.2,color='black',alpha=0.2)
plt.axvspan(2013.8,2014.2,color='black',alpha=0.2)
plt.axvspan(2004.9,2005.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.axvspan(2006.9,2007.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.axvspan(2009.9,2010.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.axvspan(2014.9,2015.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.axvspan(2015.9,2016.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.axvspan(2018.9,2019.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.xlabel('Year',size=22)
plt.ylabel('January monthly mean NO$_2$ column [10$^{15}$ molec. cm$^{-2}$]',size=22)
plt.show()

plt.figure(figsize=(15,9))
#plt.plot(yeararray,meanarray/100.,label='Non-Filtered')
#plt.plot(yeararray,meanarray_filtered/100.,label='Whole domain',color='black',marker='o')
plt.plot(yeararray,meanarray_filtered_land/100.,label='Domain',color='black',marker='o',linewidth=3)
plt.plot(yeararray,meanarray_filtered_land_anthro/100.,label='Anthropogenic',color='red',marker='o',linewidth=3)
plt.plot(yeararray,meanarray_filtered_land_bio/100.,label='Biogenic',color='green',marker='o',linewidth=3)
plt.plot(yeararray,meanarray_filtered_land_fire/100.,label='Biomass burning',color='purple',marker='o',linewidth=3)
plt.plot(yeararray,meanarray_filtered_land_lightning/100.,label='Lightning',color='blue',marker='o',linewidth=3)
#plt.plot(yeararray,array_bogota/100.,label='Bogot\xe1',color='brown',marker='o',linewidth=3)
#plt.plot(yeararray,array_caracas/100.,label='Caracas',color='cyan',marker='o',linewidth=3)
#plt.plot(yeararray,[mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota],label='Bogot\xe1 mean',color='brown',linestyle='dashed')
#plt.plot(yeararray,[mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas],label='Caracas mean',color='cyan',linestyle='dashed')
plt.legend(loc='best',fontsize=14)
plt.xlim([2005,2019])
#plt.axvspan(2009.8,2010.2,color='black',alpha=0.2)
plt.axvspan(2013.8,2014.2,color='black',alpha=0.2)
plt.axvspan(2004.9,2005.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.axvspan(2006.9,2007.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.axvspan(2009.9,2010.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.axvspan(2014.9,2015.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.axvspan(2015.9,2016.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.axvspan(2018.9,2019.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.xlabel('Year',size=22)
plt.ylabel('January monthly mean NO$_2$ column [10$^{15}$ molec. cm$^{-2}$]',size=22)
plt.show()

plt.figure(figsize=(15,9))
#plt.plot(yeararray,meanarray/100.,label='Non-Filtered')
#plt.plot(yeararray,meanarray_filtered/100.,label='Whole domain',color='black',marker='o')
#plt.plot(yeararray,meanarray_filtered_land/100.,label='Domain',color='black',marker='o',linewidth=3)
#plt.plot(yeararray,meanarray_filtered_land_anthro/100.,label='Anthropogenic',color='red',marker='o',linewidth=3)
#plt.plot(yeararray,meanarray_filtered_land_bio/100.,label='Biogenic',color='green',marker='o',linewidth=3)
#plt.plot(yeararray,meanarray_filtered_land_fire/100.,label='Biomass burning',color='orange',marker='o',linewidth=3)
#plt.plot(yeararray,meanarray_filtered_land_lightning/100.,label='Lightning',color='blue',marker='o',linewidth=3)
plt.plot(yeararray,array_bogota/100.,label='Bogot\xe1',color='brown',marker='o',linewidth=3)
plt.plot(yeararray,array_caracas/100.,label='Caracas',color='cyan',marker='o',linewidth=3)
plt.plot(yeararray,[mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota,mean_bogota],label='Bogot\xe1 mean',color='brown',linestyle='dashed')
plt.plot(yeararray,[mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas,mean_caracas],label='Caracas mean',color='cyan',linestyle='dashed')
plt.legend(loc='best',fontsize=14)
plt.xlim([2005,2019])
plt.axvspan(2009.8,2010.2,color='black',alpha=0.2)
plt.axvspan(2013.8,2014.2,color='black',alpha=0.2)
#plt.axvspan(2004.9,2005.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
#plt.axvspan(2006.9,2007.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
#plt.axvspan(2009.9,2010.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
#plt.axvspan(2014.9,2015.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
#plt.axvspan(2015.9,2016.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
#plt.axvspan(2018.9,2019.1,ymin=0,ymax=0.05,color='red',alpha=0.2)
plt.xlabel('Year',size=22)
plt.ylabel('January monthly mean NO$_2$ column [10$^{15}$ molec. cm$^{-2}$]',size=22)
plt.show()
#'''


#!SEASONAL ANALYSIS
plt.rcParams['xtick.labelsize']=20
plt.rcParams['ytick.labelsize']=20

'''
path = '/archive/ESG/barte035/Colombia/OMI/MonthlyMeans_ForSeasonalAnalysis/'
files = sorted([f for f in listdir(path) if isfile(join(path,f))])
print(files)
yeararray = ['2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019']
#yeararray = ['2005','2006']
montharray = ['01','02','03','04','05', '06','07','08','09','10','11','12']

fullarray = np.zeros(len(montharray)*len(yeararray))
fullarray_filtered = np.zeros(len(montharray)*len(yeararray))
fullarray_filtered_land = np.zeros(len(montharray)*len(yeararray))
fullarray_filtered_land_anthro = np.zeros(len(montharray)*len(yeararray))
fullarray_filtered_land_bio = np.zeros(len(montharray)*len(yeararray))
fullarray_filtered_land_fire = np.zeros(len(montharray)*len(yeararray))
fullarray_filtered_land_lightning = np.zeros(len(montharray)*len(yeararray))

meanarray = np.zeros(len(montharray))
meanarray_filtered = np.zeros(len(montharray))
meanarray_filtered_land = np.zeros(len(montharray))
meanarray_filtered_land_anthro = np.zeros(len(montharray))
meanarray_filtered_land_bio = np.zeros(len(montharray))
meanarray_filtered_land_fire = np.zeros(len(montharray))
meanarray_filtered_land_lightning = np.zeros(len(montharray))

stdarray = np.zeros(len(montharray))
stdarray_filtered = np.zeros(len(montharray))
stdarray_filtered_land = np.zeros(len(montharray))
stdarray_filtered_land_anthro = np.zeros(len(montharray))
stdarray_filtered_land_bio = np.zeros(len(montharray))
stdarray_filtered_land_fire = np.zeros(len(montharray))
stdarray_filtered_land_lightning = np.zeros(len(montharray))

fullarrayformap_lightning = np.zeros((99,99,len(files)))
arrayformap_lightning = np.zeros((99,99,len(montharray)))
arrayformap_lightning_std = np.zeros((99,99,len(montharray)))

geopath = '/archive/ESG/barte035/WPS/geo_em.d01.nc'
geo = nc.Dataset(geopath,'r')
for i in range(len(files)):
	print(files[i])
	longitude,latitude,value = read_grd(path+files[i])
		
	#Extract WRF-Chem lat/lon and regrid data
	geo = nc.Dataset(geopath,'r')
	we  = geo.variables['XLAT_M'].shape[2]
	sn  = geo.variables['XLAT_M'].shape[1]
	lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
	lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
	dx_dom = 27000
	
	#plt.figure(figsize=(15,10))
    	m = Basemap(width=(we*dx_dom)-1,height=(sn*dx_dom)-1,resolution='l',area_thresh=100.,projection='lcc',lat_0=4.891087,lon_0=-71.069204,lat_1=40.,lat_2=70.,lon_1=-70.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120)
	nx = 99
	ny = 99
	x2_trans = np.linspace(-179.9375,179.9375,2880)
	y2_trans = np.linspace(-89.9375,89.9375,1440)
	
	value = np.where(value < 0, np.nan, value)

	#value_trans = m.transform_scalar(value[::-1,:],lons=np.arange(-82.375,-60.,0.125),lats=np.arange(-5.,14.625,0.125),nx=100.,ny=100.,order=1)
	value_trans = m.transform_scalar(value[::-1],lons=x2_trans,lats=y2_trans,nx=nx,ny=ny,order=0)
	value_trans_filtered = np.where(value_trans < 0, np.nan, value_trans)
	value_trans_filtered_land = np.where(landmask==1, value_trans_filtered, np.nan)
	value_trans_filtered_land_anthro = np.where(~np.isnan(NOxMaxAnthro) > 0., value_trans_filtered_land, np.nan)
	value_trans_filtered_land_bio = np.where(~np.isnan(NOxMaxBio) > 0., value_trans_filtered_land, np.nan)
	value_trans_filtered_land_fire = np.where(~np.isnan(NOxMaxFire) > 0., value_trans_filtered_land, np.nan)
	value_trans_filtered_land_lightning = np.where(~np.isnan(NOxMaxLightning) > 0., value_trans_filtered_land, np.nan)
	fullarrayformap_lightning[:,:,i] = value_trans_filtered_land_lightning
	
	print(value_trans.shape)
	
	fullarray[i] = np.nanmean(value_trans)
	fullarray_filtered[i] = np.nanmean(value_trans_filtered)
	fullarray_filtered_land[i] = np.nanmean(value_trans_filtered_land)
	fullarray_filtered_land_anthro[i] = np.nanmean(value_trans_filtered_land_anthro)
	fullarray_filtered_land_bio[i] = np.nanmean(value_trans_filtered_land_bio)
	fullarray_filtered_land_fire[i] = np.nanmean(value_trans_filtered_land_fire)
	fullarray_filtered_land_lightning[i] = np.nanmean(value_trans_filtered_land_lightning)
	
	print(fullarray_filtered_land_lightning.shape)
	
	
for j in range(len(montharray)):
	meanarray[j] = np.mean(fullarray[j::12])
	meanarray_filtered[j] = np.mean(fullarray_filtered[j::12])
	meanarray_filtered_land[j] = np.mean(fullarray_filtered_land[j::12])
	meanarray_filtered_land_anthro[j] = np.mean(fullarray_filtered_land_anthro[j::12])
	meanarray_filtered_land_bio[j] = np.mean(fullarray_filtered_land_bio[j::12])
	meanarray_filtered_land_fire[j] = np.mean(fullarray_filtered_land_fire[j::12])
	meanarray_filtered_land_lightning[j] = np.mean(fullarray_filtered_land_lightning[j::12])	

	stdarray[j] = np.std(fullarray[j::12])
	stdarray_filtered[j] = np.std(fullarray_filtered[j::12])
	stdarray_filtered_land[j] = np.std(fullarray_filtered_land[j::12])
	stdarray_filtered_land_anthro[j] = np.std(fullarray_filtered_land_anthro[j::12])
	stdarray_filtered_land_bio[j] = np.std(fullarray_filtered_land_bio[j::12])
	stdarray_filtered_land_fire[j] = np.std(fullarray_filtered_land_fire[j::12])
	stdarray_filtered_land_lightning[j] = np.std(fullarray_filtered_land_lightning[j::12])
	
	arrayformap_lightning[:,:,j] = np.mean(fullarrayformap_lightning[:,:,j::12],axis=2)
	arrayformap_lightning_std[:,:,j] = np.std(fullarrayformap_lightning[:,:,j::12],axis=2)
'''

'''
np.save('/archive/ESG/barte035/Colombia/OMI/arrayformap_lightning.pkl',arrayformap_lightning)
np.save('/archive/ESG/barte035/Colombia/OMI/arrayformap_lightning_std.pkl',arrayformap_lightning_std)
'''

montharrint = [1,2,3,4,5,6,7,8,9,10,11,12]
'''
np.save('/archive/ESG/barte035/Colombia/OMI/meanarrland.pkl',meanarray_filtered_land)
np.save('/archive/ESG/barte035/Colombia/OMI/stdarrland.pkl',stdarray_filtered_land)
np.save('/archive/ESG/barte035/Colombia/OMI/meanarrlandanthro.pkl',meanarray_filtered_land_anthro)
np.save('/archive/ESG/barte035/Colombia/OMI/stdarrlandanthro.pkl',stdarray_filtered_land_anthro)
np.save('/archive/ESG/barte035/Colombia/OMI/meanarrlandbio.pkl',meanarray_filtered_land_bio)
np.save('/archive/ESG/barte035/Colombia/OMI/stdarrlandbio.pkl',stdarray_filtered_land_bio)
np.save('/archive/ESG/barte035/Colombia/OMI/meanarrlandfire.pkl',meanarray_filtered_land_fire)
np.save('/archive/ESG/barte035/Colombia/OMI/stdarrlandfire.pkl',stdarray_filtered_land_fire)
np.save('/archive/ESG/barte035/Colombia/OMI/meanarrlandlightning.pkl',meanarray_filtered_land_lightning)
np.save('/archive/ESG/barte035/Colombia/OMI/stdarrlandlightning.pkl',stdarray_filtered_land_lightning)
'''

meanarray_filtered_land = np.load('/archive/ESG/barte035/Colombia/OMI/meanarrland.pkl.npy')
stdarray_filtered_land = np.load('/archive/ESG/barte035/Colombia/OMI/stdarrland.pkl.npy')
meanarray_filtered_land_anthro = np.load('/archive/ESG/barte035/Colombia/OMI/meanarrlandanthro.pkl.npy')
stdarray_filtered_land_anthro = np.load('/archive/ESG/barte035/Colombia/OMI/stdarrlandanthro.pkl.npy')
meanarray_filtered_land_bio = np.load('/archive/ESG/barte035/Colombia/OMI/meanarrlandbio.pkl.npy')
stdarray_filtered_land_bio = np.load('/archive/ESG/barte035/Colombia/OMI/stdarrlandbio.pkl.npy')
meanarray_filtered_land_fire = np.load('/archive/ESG/barte035/Colombia/OMI/meanarrlandfire.pkl.npy')
stdarray_filtered_land_fire = np.load('/archive/ESG/barte035/Colombia/OMI/stdarrlandfire.pkl.npy')
meanarray_filtered_land_lightning = np.load('/archive/ESG/barte035/Colombia/OMI/meanarrlandlightning.pkl.npy')
stdarray_filtered_land_lightning = np.load('/archive/ESG/barte035/Colombia/OMI/stdarrlandlightning.pkl.npy')

plt.figure(figsize=(15,9))
plt.plot(montharrint,meanarray_filtered_land/100.,label='Domain',color='black',marker='o',linewidth=3)
plt.plot(montharrint,meanarray_filtered_land_anthro/100.,label='Anthropogenic',color='red',marker='o',linewidth=3)
plt.plot(montharrint,meanarray_filtered_land_bio/100.,label='Biogenic',color='green',marker='o',linewidth=3)
plt.plot(montharrint,meanarray_filtered_land_fire/100.,label='Biomass burning',color='purple',marker='o',linewidth=3)
plt.plot(montharrint,meanarray_filtered_land_lightning/100.,label='Lightning',color='blue',marker='o',linewidth=3)
plt.fill_between(montharrint,(meanarray_filtered_land-stdarray_filtered_land)/100.,(meanarray_filtered_land+stdarray_filtered_land)/100.,alpha=0.2,color='black')
plt.fill_between(montharrint,(meanarray_filtered_land_anthro-stdarray_filtered_land_anthro)/100.,(meanarray_filtered_land_anthro+stdarray_filtered_land_anthro)/100.,alpha=0.2,color='red')
plt.fill_between(montharrint,(meanarray_filtered_land_bio-stdarray_filtered_land_bio)/100.,(meanarray_filtered_land_bio+stdarray_filtered_land_bio)/100.,alpha=0.2,color='green')
plt.fill_between(montharrint,(meanarray_filtered_land_fire-stdarray_filtered_land_fire)/100.,(meanarray_filtered_land_fire+stdarray_filtered_land_fire)/100.,alpha=0.2,color='purple')
plt.fill_between(montharrint,(meanarray_filtered_land_lightning-stdarray_filtered_land_lightning)/100.,(meanarray_filtered_land_lightning+stdarray_filtered_land_lightning)/100.,alpha=0.2,color='blue')
plt.legend(loc='best',fontsize=14)
plt.xlim([1,12])
plt.xticks(montharrint,('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Okt','Nov','Dec'))
plt.xlabel('Month',size=22)
plt.ylabel('Monthly mean NO$_2$ column [10$^{15}$ molec. cm$^{-2}$]',size=22)
plt.savefig('/home/WUR/barte035/WRFChem/Python/Colombia/FiguresColombia/fig8b_afterproofread.png',dpi=300)
plt.show()

arrayformap_lightning = np.load('/archive/ESG/barte035/Colombia/OMI/arrayformap_lightning.pkl.npy')
arrayformap_lightning_std = np.load('/archive/ESG/barte035/Colombia/OMI/arrayformap_lightning_std.pkl.npy')

for i in range(len(montharrint)):
	plt.figure(figsize=(15,9))
	m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
	m.imshow(arrayformap_lightning[:,:,i]/100.,cmap='viridis',vmin=0.,vmax=1.0,interpolation='none')
	m.drawcoastlines()
	m.drawcountries()
	c = plt.colorbar(ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9,1.0])
	c.set_label('OMI NO$_2$ column [10$^{15}$ molec. cm$^{-2}$]',size=14)
	plt.savefig('/home/WUR/barte035/lightningcolumn_'+str(montharrint[i])+'.png',dpi=300)
	plt.show()

plt.figure(figsize=(15,9))
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
m.imshow(np.nanmean(arrayformap_lightning[:,:,:],axis=2)/100.,cmap='viridis',vmin=0.,vmax=1.0,interpolation='none')
m.drawcoastlines()
m.drawcountries()
c = plt.colorbar(ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9,1.0])
c.set_label('OMI NO$_2$ column [10$^{15}$ molec. cm$^{-2}$]',size=14)
plt.savefig('/home/WUR/barte035/lightningcolumn_mean.png',dpi=300)
plt.show()
