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

daysim = 30 #30100
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
#m.imshow(np.where((landmask==1),NOxMaxFire,np.nan),cmap='Oranges',vmin=0,vmax=100,interpolation='none')
m.imshow(np.where((landmask==1),NOxMaxFire,np.nan),cmap='Purples',vmin=0,vmax=100,interpolation='none')
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
custom_legend = [Patch(facecolor='red',label='Anthropogenic'),Patch(facecolor='green',label='Biogenic'),Patch(facecolor='purple',label='Biomass burning'),Patch(facecolor='blue',label='Lightning')]
plt.legend(custom_legend,['Anthropogenic','Biogenic','Biomass burning','Lightning'],loc='upper left',prop={'size': 9})
#plt.text(2.5e6,2.8e6,'(b)',fontsize=9)
plt.show()


fig=plt.figure(figsize=(15,9))
ax=fig.add_subplot(211)
#m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,projection='lcc',lat_0=lat_0,lon_0=lon_0)
#m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,llcrnrlon=-81.684021,llcrnrlat=-7.9757309,urcrnrlon=-60.025818,urcrnrlat=15.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0)
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
m.imshow(np.where((landmask==1),NOxMaxBio,np.nan),cmap='Greens',vmin=0,vmax=100,interpolation='none')
m.imshow(np.where((landmask==1),NOxMaxFire,np.nan),cmap='Purples',vmin=0,vmax=100,interpolation='none')
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
custom_legend = [Patch(facecolor='red',label='Anthropogenic'),Patch(facecolor='green',label='Biogenic'),Patch(facecolor='purple',label='Biomass burning'),Patch(facecolor='blue',label='Lightning')]
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

#_____________________________________________

lat = ncfile.variables['XLAT'][timesim,:,:]
lon = ncfile.variables['XLONG'][timesim,:,:]
we = len(ncfile.dimensions['west_east'])
sn = len(ncfile.dimensions['south_north'])
dx = 27000.
width = ((we-11.)*dx)-1
height = ((sn-11.)*dx)-1
lat_0 = 4.891087
lon_0 = -71.069204

NOxMaxBio = NOxMaxBio[5:len(NOxMaxBio)-5,5:len(NOxMaxBio)-5]
NOxMaxFire = NOxMaxFire[5:len(NOxMaxFire)-5,5:len(NOxMaxFire)-5]
NOxMaxLightning = NOxMaxLightning[5:len(NOxMaxLightning)-5,5:len(NOxMaxLightning)-5]
NOxMaxAnthro = NOxMaxAnthro[5:len(NOxMaxAnthro)-5,5:len(NOxMaxAnthro)-5]
NOxALLConvert = NOxALLConvert[5:len(NOxALLConvert)-5,5:len(NOxALLConvert)-5]
landmask = landmask[5:len(landmask)-5,5:len(landmask)-5]

fig=plt.figure(figsize=(15,9))
ax=fig.add_subplot(211)
#m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,projection='lcc',lat_0=lat_0,lon_0=lon_0)
#m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,llcrnrlon=-81.684021,llcrnrlat=-7.9757309,urcrnrlon=-60.025818,urcrnrlat=15.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0)
#m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
#m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-78.464021,llcrnrlat=-3.7557309,urcrnrlon=-62.245818,urcrnrlat=13.422120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=40.,lat_2=70.,lon_1=-70.)
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.164021,llcrnrlat=-4.1557309,urcrnrlon=-61.445818,urcrnrlat=13.852120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=-10.,lat_2=80.,lon_1=-70.)
m.imshow(np.where((landmask==1),NOxMaxBio,np.nan),cmap='Greens',vmin=0,vmax=100,interpolation='none')
m.imshow(np.where((landmask==1),NOxMaxFire,np.nan),cmap='Purples',vmin=0,vmax=100,interpolation='none')
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
custom_legend = [Patch(facecolor='red',label='Anthropogenic'),Patch(facecolor='green',label='Biogenic'),Patch(facecolor='purple',label='Biomass burning'),Patch(facecolor='blue',label='Lightning')]
plt.legend(custom_legend,['Anthropogenic','Biogenic','Biomass burning','Lightning'],loc='upper left',prop={'size': 7})
plt.text(1.5e6,1.58e6,'(b)',fontsize=9)


ax=fig.add_subplot(221)
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.164021,llcrnrlat=-4.1557309,urcrnrlon=-61.445818,urcrnrlat=13.852120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=-10.,lat_2=80.,lon_1=-70.)
m.imshow(NOxALLConvert/(27.*27.)*1000.,cmap='plasma_r',vmin=0,vmax=140,interpolation='none')
m.drawcoastlines()
m.drawcountries()
#m.drawparallels(np.arange(-20,30,10),labels=[True,False,True])
#m.drawmeridians(np.arange(-90,-50,20),labels=[True,False,True])
c = plt.colorbar(ticks=[0,20,40,60,80,100,120,140],extend='max')
c.set_label(r'NO$_{x}$ Flux [kg N km$^{-2}$ month$^{-1}$]',size=14)
plt.text(90000,1.58e6,'(a)',fontsize=9)
plt.show()

'''
ax=fig.add_subplot(221)
m = Basemap(width=width,height=height,resolution='l',area_thresh=100.,llcrnrlon=-79.164021,llcrnrlat=-4.1557309,urcrnrlon=-61.445818,urcrnrlat=13.852120,projection='lcc',lat_0=lat_0,lon_0=lon_0,lat_1=-10.,lat_2=80.,lon_1=-70.)
m.imshow(NOxALLConvert,cmap='plasma_r',vmin=0,vmax=100,interpolation='none')
m.drawcoastlines()
m.drawcountries()
#m.drawparallels(np.arange(-20,30,10),labels=[True,False,True])
#m.drawmeridians(np.arange(-90,-50,20),labels=[True,False,True])
c = plt.colorbar(ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='max')
c.set_label(r'NO$_{x}$ Flux [Mg N month$^{-1}$]',size=14)
plt.text(90000,1.58e6,'(a)',fontsize=9)
plt.show()
'''

print np.count_nonzero(~np.isnan(np.where((landmask==1),NOxMaxBio,np.nan)))
print np.count_nonzero(~np.isnan(np.where((landmask==1),NOxMaxFire,np.nan)))
print np.count_nonzero(~np.isnan(np.where((landmask==1),NOxMaxAnthro,np.nan)))
print np.count_nonzero(~np.isnan(np.where((landmask==1),NOxMaxLightning,np.nan)))



#return(NOxMaxLightning,NOxMaxAnthro,NOxMaxBio,NOxMaxFire)
