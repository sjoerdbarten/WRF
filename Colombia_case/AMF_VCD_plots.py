"AMF recalculation, cosampling, plotting"

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from matplotlib import cm
from scipy.stats import linregress
import matplotlib.style
import matplotlib as mpl
from matplotlib import pyplot
from matplotlib.offsetbox import AnchoredText
from pylab import *
from sklearn.metrics import r2_score
from matplotlib.gridspec import GridSpec


mpl.style.use('classic')
mpl.style.available
mpl.rcParams['figure.facecolor'] = '1.0'
mpl.rcParams['mathtext.default'] = 'regular'

'''
ncfile = '/home/WUR/barte035/WRFChem/OMI/gridded/qa4ecv_amf_201401.nc'
ds = nc.Dataset(ncfile,'r')
AMF_QA = ds.variables['amf_retreival'][:][680:840,776:960] #[960:1296,1167:1808] = EUROPE

plt.figure(figsize=(15,9))
m = Basemap(projection='cyl',llcrnrlat=-5,urcrnrlat=15,llcrnrlon=-83,urcrnrlon=-60,resolution='l')
m.drawcoastlines(zorder=4)
m.drawcountries(zorder=4)
#plot_update_range(m)
m.drawmapboundary(fill_color='lightgray')
#m.fillcontinents(color='white',lake_color='lightgray') 
m.imshow(AMF_QA,vmin=0,vmax=2,interpolation='none')
c = plt.colorbar()
c.set_label('AMF [-]')
plt.title('QA4ECV Air Mass Factors')
plt.show()

ncfile = '/home/WUR/barte035/WRFChem/OMI/gridded/qa4ecv_amf_201401_mod.nc'
ds = nc.Dataset(ncfile,'r')
AMF_mod = ds.variables['mod_amf'][:][680:840,776:960]

plt.figure(figsize=(15,9))
#m = Basemap(projection='cyl',llcrnrlat=30,urcrnrlat=70,llcrnrlon=-15,urcrnrlon=40,resolution='l')
m = Basemap(projection='cyl',llcrnrlat=-5,urcrnrlat=15,llcrnrlon=-83,urcrnrlon=-60,resolution='l')
#m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='l')
m.drawcoastlines(zorder=4)
m.drawcountries(zorder=4)
#plot_update_range(m)
m.drawmapboundary(fill_color='lightgray')
#m.fillcontinents(color='white',lake_color='lightgray')
m.imshow(AMF_mod,vmin=0,vmax=2,interpolation='none')
c = plt.colorbar()
c.set_label('AMF [-]')
plt.title('Re-calculated Air Mass Factors')
plt.show()

print "np.nanmin(AMF_mod - AMF_QA)", np.nanmin(AMF_mod - AMF_QA), "np.nanmax(AMF_mod - AMF_QA)", np.nanmax(AMF_mod - AMF_QA)

plt.figure(figsize=(15,9))
m = Basemap(projection='cyl',llcrnrlat=-5,urcrnrlat=15,llcrnrlon=-83,urcrnrlon=-60,resolution='l')
m.drawcoastlines(zorder=4)
m.drawcountries(zorder=4)
#plot_update_range(m)
m.drawmapboundary(fill_color='lightgray')
#m.fillcontinents(color='white',lake_color='lightgray') 
m.imshow(AMF_mod - AMF_QA,vmin=-1,vmax=1,cmap='bwr',interpolation='none')
c = plt.colorbar(extend='both')
c.set_label(r'$\Delta$AMF [-]')
plt.title('Air Mass Factor difference (Recalculated - Retreival)')
plt.show()

print np.nanmean(AMF_mod - AMF_QA)
print np.nanstd(AMF_mod - AMF_QA)


plt.figure(figsize=(15,9))
m.drawcoastlines(zorder=4)
m.drawcountries(zorder=4)
#plot_update_range(m)
m.drawmapboundary(fill_color='lightgray')
#m.fillcontinents(color='white',lake_color='lightgray') 
m.imshow((AMF_mod - AMF_QA)/AMF_QA,vmin=-0.5,vmax=0.5,cmap='bwr',interpolation='none',zorder=3)
c = plt.colorbar(extend='both')
c.set_label(r'$\Delta$AMF [-]')
plt.title('Air Mass Factor relative difference (Recalculated - Retreival)/Retreival')
plt.show()

print "np.nanmin(AMF_mod - AMF_QA)", np.nanmin(AMF_mod - AMF_QA), "np.nanmax(AMF_mod - AMF_QA)", np.nanmax(AMF_mod - AMF_QA)
'''

execfile("OMI_read.py")

no2stack_mod,no2trans_mod = OMI_read('2014_mod')
no2stack,no2trans = OMI_read(2014)
diff = no2trans_mod - no2trans

geopath = '/archive/ESG/barte035/WPS/geo_em.d01.nc'

'''
#Extract WRF-Chem lat/lon and regrid data
geo = nc.Dataset(geopath,'r')
we  = geo.variables['XLAT_M'].shape[2]
sn  = geo.variables['XLAT_M'].shape[1]
lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
dx_dom = 27000

m = Basemap(width=(we*dx_dom)-1,height=(sn*dx_dom)-1,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=4.891087,lon_0=-71.069204,lat_1=40.,lat_2=70.,lon_1=-70.)
plt.figure(figsize=(15,9))
#im = m.imshow(no2_nc_trans,norm=LogNorm(),vmin=0.01,vmax=30,interpolation='none')
im = m.imshow(diff,vmin=-0.5,vmax=0.5,cmap='bwr',interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
#c = plt.colorbar()
c = plt.colorbar(extend='both')
c.set_ticks(np.arange(-0.5,0.51,0.25))
c.set_ticklabels(np.arange(-0.5,0.51,0.25))
c.set_label(u'$\Delta$C$_{NO_2}$ [$10^{15}$ molec. cm$^{-2}$]')
plt.title('NO_2 VCD difference (recalculated - QA4ECV)',weight='bold')
plt.show()

print np.nanmean(diff)
print np.nanstd(diff)




fig=plt.figure(figsize=(15,9))
m = Basemap(width=(we*dx_dom)-1,height=(sn*dx_dom)-1,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=4.891087,lon_0=-71.069204,lat_1=40.,lat_2=70.,lon_1=-70.)
ax=fig.add_subplot(211)
im = m.imshow(diff,vmin=-0.5,vmax=0.5,cmap='bwr',interpolation='none')
m.drawcoastlines(zorder=4)
m.drawcountries(zorder=4)
m.drawmapboundary(fill_color='lightgray')
c = plt.colorbar(extend='both')
c.set_ticks(np.arange(-0.5,0.51,0.25))
c.set_ticklabels(np.arange(-0.5,0.51,0.25))
c.set_label(u'$\Delta$C$_{NO_2}$ [$10^{15}$ molec. cm$^{-2}$]')

ax=fig.add_subplot(221)
m = Basemap(projection='cyl',llcrnrlat=-5,urcrnrlat=15,llcrnrlon=-83,urcrnrlon=-60,resolution='l')
m.drawcoastlines(zorder=4)
m.drawcountries(zorder=4)
m.drawmapboundary(fill_color='lightgray')
m.imshow(AMF_mod - AMF_QA,vmin=-1,vmax=1,cmap='bwr',interpolation='none')
c = plt.colorbar(extend='both')
c.set_label(r'$\Delta$AMF [-]')
plt.text(-82.5,13.5,'(a)',fontsize=9)
plt.text(-50,13.5,'(b)',fontsize=9)
plt.show()







print "np.nanmin(diff)", np.nanmin(diff), "np.nanmax(diff)", np.nanmax(diff)

plt.figure(figsize=(15,9))
#im = m.imshow(no2_nc_trans,norm=LogNorm(),vmin=0.01,vmax=30,interpolation='none')
im = m.imshow(diff/no2trans,vmin=-0.5,vmax=0.5,cmap='bwr',interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
#c = plt.colorbar()
c = plt.colorbar(extend='both')
c.set_ticks(np.arange(-0.5,1.,0.5))
c.set_ticklabels(np.arange(-0.5,1.,0.5))
c.set_label(u'$\\Delta$C$_{NO_2}$ [$10^{15}$ molec. cm$^{-2}$]',size=15)
plt.title('NO_2 relative VCD difference (recalculated - QA4ECV)',size=18,weight='bold')
plt.show()

print "np.nanmin(diff/no2trans)", np.nanmin(diff/no2trans), "np.nanmean(diff/no2trans)", np.nanmean(diff/no2trans), "np.nanmax(diff/no2trans)", np.nanmax(diff/no2trans)
'''

## WRF-Chem-OMI comparison
execfile("WRF_columndensity_alt.py")

wrfdir = '/archive/ESG/barte035/WRFV3/'  #THIS
P,colstack = wrf_columns(wrfdir,50) #wrfdir, leveltropopause
#Co-sample WRF-Chem with OMI
print colstack.shape
print no2stack_mod.shape
colstack_cosample = colstack
colstack_cosample[np.isnan(no2stack_mod)] = np.nan
nonans = (~np.isnan(colstack_cosample)).sum(axis=0)
colmean_cosample = np.nanmean(colstack_cosample,axis=0)
colmean_cosample[nonans < 2] = np.nan #Analogous to OMI_read.py
#colmean = np.nanmean(colstack,axis=0)
#print 'Mean NO2 column density from WRF', np.nansum(colmean * 7.29e12 / 6.022e23 * 14 * 1.e-12) #changed based on area WRF cell
#Extract WRF-Chem lat/lon and regrid data
geo = nc.Dataset(geopath,'r')
we  = geo.variables['XLAT_M'].shape[2]
sn  = geo.variables['XLAT_M'].shape[1]
lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
dx_dom = 27000


'''
m = Basemap(width=(we*dx_dom)-1,height=(sn*dx_dom)-1,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=4.891087,lon_0=-71.069204,lat_1=40.,lat_2=70.,lon_1=-70.)
plt.figure(figsize=(15,9))
m.imshow(colmean_cosample/1.e15,vmin=0.1,vmax=8,norm=LogNorm(),interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
#m.drawparallels(np.arange(40,61,10),labels=[True,False,True])
#m.drawmeridians(np.arange(-10,51,20),labels=[True,False,True])
c = plt.colorbar(extend='max')
c.set_ticks([0.1,0.5,1,2,3,4,6,8])
c.set_ticklabels([0.1,0.5,1,2,3,4,6,8])
c.set_label(r'NO_2 column density [$10^{15}$ molec. cm$^{-2}$]',size=15)
plt.title('WRF NO_2 column',weight='bold')
plt.show()

plt.figure(figsize=(15,9))
m.imshow(colmean_cosample/1.e15 - no2trans_mod,cmap='bwr',vmin=-5,vmax=5,interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
#m.drawparallels(np.arange(40,61,10),labels=[True,False,True],size=10)
#m.drawmeridians(np.arange(-10,51,20),labels=[True,False,True],size=10)
c = plt.colorbar(ticks=np.arange(-5,5.1,0.5),extend='both')
c.set_label(r'$\Delta \hat{x}_{NO_2}$ [$10^{15}$ molec. cm$^{-2}$]',size=15)
plt.title('WRF - OMI',weight='bold')
plt.show()

print "np.nanmin((colmean_cosample/1.e15 - no2trans_mod)/no2trans_mod)", np.nanmin((colmean_cosample/1.e15 - no2trans_mod)/no2trans_mod), "np.nanmax((colmean_cosample/1.e15 - no2trans_mod)/no2trans_mod)", np.nanmax((colmean_cosample/1.e15 - no2trans_mod)/no2trans_mod)

plt.figure(figsize=(15,9))
m.imshow((colmean_cosample/1.e15 - no2trans_mod)/no2trans_mod,cmap='bwr',vmin=-5,vmax=5,interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
#m.drawparallels(np.arange(40,61,10),labels=[True,False,True])
#m.drawmeridians(np.arange(-10,51,20),labels=[True,False,True])
c = plt.colorbar(ticks=np.arange(-5,5.1,0.5),extend='both')
c.set_label(r'(WRF-Chem - OMI)/OMI [$10^{15}$ molec cm$^{-2}$]',size=15)
plt.title('(WRF - OMI) / OMI',weight='bold')
plt.show()
'''

var = colmean_cosample/1.e15 - no2trans_mod

'''
print "np.nanmin(var)",np.nanmin(var),"np.nanpercentile(var,1)",np.nanpercentile(var,1),"np.nanpercentile(var,10)",np.nanpercentile(var,10),"np.nanmean(var)",np.nanmean(var),"np.nanmedian(var)",np.nanmedian(var),"np.nanpercentile(var,90)",np.nanpercentile(var,90),"np.nanpercentile(var,99)",np.nanpercentile(var,99),"np.nanmax(var)",np.nanmax(var)

plt.boxplot([(colmean_cosample/1.e15)[~np.isnan(var)],no2trans_mod[~np.isnan(var)]])
#plt.boxplot((var)[~np.isnan(var)])
plt.title('Boxplot WRF and OMI')
plt.xlabel(['WRF-Chem','OMI'])
plt.ylabel('NO_2 Column Density [$10^{15}$ molec cm$^{-2}$]')
plt.show()

plt.hist([(colmean_cosample/1.e15)[~np.isnan(var)],no2trans_mod[~np.isnan(var)]], bins=50,label=['WRF','OMI'])
plt.legend(loc='upper right')
plt.xlabel('NO_2 Column Density [$10^{15}$ molec cm$^{-2}$]')
plt.ylabel('Frequency')
plt.show()

pyplot.hist((colmean_cosample/1.e15)[~np.isnan(var)], bins=60, range=(0,3), alpha=0.5, label='WRF-Chem')
pyplot.hist(no2trans_mod[~np.isnan(var)], bins=60, range=(0,3), alpha=0.5, label='OMI')
pyplot.legend(loc='upper right')
pyplot.xlabel('NO_2 Column Density [$10^{15}$ molec cm$^{-2}$]')
pyplot.ylabel('Frequency')
pyplot.show()

plt.hist((var)[~np.isnan(var)],bins=100,range=(-5,5),color='c',edgecolor='k',alpha=0.7)
plt.title('Histogram WRF-OMI + mean + median, 95% percentiles')
plt.axvline((var)[~np.isnan(var)].mean(), color='k', linestyle='dashed', linewidth=1)
plt.axvline(np.median((var)[~np.isnan(var)]), color='k', linestyle='dotted', linewidth=1)
plt.axvline(np.nanpercentile(var,5), color='k', linestyle='dashdot', linewidth=1)
plt.axvline(np.nanpercentile(var,95), color='k', linestyle='dashdot', linewidth=1)
plt.axvline(0, color='k', linewidth=2)
plt.xlim((-2,2))
plt.xlabel('VCD$_{WRF-Chem}$-VCD$_{OMI}$ [$10^{15}$ molec cm$^{-2}$]')
plt.ylabel('Frequency')
plt.show()

plt.scatter((colmean_cosample/1.e15)[~np.isnan(var)],no2trans_mod[~np.isnan(var)],s=2,linewidths=0)
plt.plot([0,100],[0,100],'r-')
#plt.xlim([0,2])
#plt.ylim([0,2])
plt.xlabel('WRF-Chem VCD [$10^{15}$ molec cm$^{-2}$]')
plt.ylabel('OMI VCD [$10^{15}$ molec cm$^{-2}$]')
#plt.yscale('log')
#plt.xscale('log')
plt.xlim([0,12])
plt.ylim([0,12])
plt.legend(['1:1 line','Data'],fontsize=7)
plt.show()
'''

m = Basemap(width=(we*dx_dom)-1,height=(sn*dx_dom)-1,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=4.891087,lon_0=-71.069204,lat_1=40.,lat_2=70.,lon_1=-70.)

fig=plt.figure(figsize=(11,15))
subplots_adjust(wspace=0.2,hspace=0.2)
ax=fig.add_subplot(321)
m.imshow(colmean_cosample/1.e15,vmin=0.1,vmax=8,norm=LogNorm(),interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
c = plt.colorbar(extend='max')
c.set_ticks([0.1,0.5,1,2,3,4,6,8])
c.set_ticklabels([0.1,0.5,1,2,3,4,6,8])
c.set_label(r'VCD$_{WRF-Chem}$ [$10^{15}$ molec. cm$^{-2}$]',size=11)
plt.text(100000,2.8e6,'(a)',fontsize=9)

ax=fig.add_subplot(322)
im = m.imshow(no2trans_mod,vmin=0.1,vmax=8,norm=LogNorm(),interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
c = plt.colorbar(extend='max')
c.set_ticks([0.1,0.5,1,2,3,4,6,8])
c.set_ticklabels([0.1,0.5,1,2,3,4,6,8])
c.set_label(r'VCD$_{OMI}$ [$10^{15}$ molec. cm$^{-2}$]',size=11)
plt.text(100000,2.8e6,'(b)',fontsize=9)


ax=fig.add_subplot(323)
m.imshow(colmean_cosample/1.e15 - no2trans_mod,cmap='bwr',vmin=-3,vmax=3,interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
c = plt.colorbar(ticks=np.arange(-3,3.1,1.),extend='both')
c.set_label(r'VCD$_{WRF-Chem}$-VCD$_{OMI}$ [$10^{15}$ molec. cm$^{-2}$]',size=11)
plt.text(100000,2.8e6,'(c)',fontsize=9)

'''
ax=fig.add_subplot(323)
m.imshow(((colmean_cosample/1.e15 - no2trans_mod)/no2trans_mod)*100.,cmap='bwr',vmin=-100.,vmax=100.,interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
c = plt.colorbar(ticks=np.arange(-100.,100.,25.),extend='both')
c.set_label(r'VCD$_{WRF-Chem}$-VCD$_{OMI}$ / VCD$_{OMI}$ [%]',size=11)
#plt.text(100000,2.8e6,'(c)',fontsize=9)
'''

ax=fig.add_subplot(324)
plt.scatter((colmean_cosample/1.e15)[~np.isnan(var)],no2trans_mod[~np.isnan(var)],s=2,linewidths=0)
plt.plot([0,100],[0,100],'r-')
plt.xlabel('VCD$_{WRF-Chem}$ [$10^{15}$ molec. cm$^{-2}$]',size=11)
plt.ylabel('VCD$_{OMI}$ [$10^{15}$ molec. cm$^{-2}$]',size=11)
plt.xlim([0,9.9])
plt.ylim([0,9.9])
plt.legend(['1:1 line','Data'],fontsize=9)
plt.text(0.3,9.2,'(d)',fontsize=9)

ax=fig.add_subplot(325)
pyplot.hist((colmean_cosample/1.e15)[~np.isnan(var)], bins=60, range=(0,3), alpha=0.5, label='WRF-Chem')
pyplot.hist(no2trans_mod[~np.isnan(var)], bins=60, range=(0,3), alpha=0.5, label='OMI')
pyplot.legend(loc='upper right',fontsize=9)
pyplot.xlabel(r'VCD [$10^{15}$ molec. cm$^{-2}$]',size=11)
pyplot.ylabel('Frequency',size=11)
plt.text(0.1,1110,'(e)',fontsize=9)

ax=fig.add_subplot(326)
plt.hist((var)[~np.isnan(var)],bins=100,range=(-5,5),color='c',edgecolor='k',alpha=0.7)
plt.axvline((var)[~np.isnan(var)].mean(), color='k', linestyle='dashed', linewidth=1)
plt.axvline(np.median((var)[~np.isnan(var)]), color='k', linestyle='dotted', linewidth=1)
plt.axvline(np.nanpercentile(var,5), color='k', linestyle='dashdot', linewidth=1)
plt.axvline(np.nanpercentile(var,95), color='k', linestyle='dashdot', linewidth=1)
plt.axvline(0, color='k', linewidth=2)
plt.xlim((-2,2))
plt.xlabel(r'VCD$_{WRF-Chem}$-VCD$_{OMI}$ [$10^{15}$ molec. cm$^{-2}$]',size=11)
plt.ylabel(' ',size=11)
plt.legend(['Mean','Median','90% conf.'],fontsize=9,loc=2)
plt.text(1.72,1500,'(f)',fontsize=9)
#plt.savefig('/home/WUR/barte035/WRFChem/Python/Figures/WRFOMIcomparison.png',bbox_inches='tight')
	
'''
ax=fig.add_subplot(326)
plt.hist((var/no2trans_mod*100)[~np.isnsafean(var)],bins=100,range=(-200.,200.),color='c',edgecolor='k',alpha=0.7)
plt.axvline((var/no2trans_mod*100)[~np.isnan(var)].mean(), color='k', linestyle='dashed', linewidth=1)
plt.axvline(np.median((var/no2trans_mod*100)[~np.isnan(var)]), color='k', linestyle='dotted', linewidth=1)
plt.axvline(np.nanpercentile(var/no2trans_mod*100,5), color='k', linestyle='dashdot', linewidth=1)
plt.axvline(np.nanpercentile(var/no2trans_mod*100,95), color='k', linestyle='dashdot', linewidth=1)
plt.axvline(0, color='k', linewidth=2)
plt.xlim((-200.,200.))
plt.xlabel(r'VCD$_{WRF-Chem}$-VCD$_{OMI}$/VCD$_{OMI}$ [%]',size=11)
plt.ylabel(' ',size=11)
plt.legend(['Mean','Median','90% conf.'],fontsize=9,loc=2)
'''
plt.show()



m = Basemap(width=(we*dx_dom)-1,height=(sn*dx_dom)-1,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=4.891087,lon_0=-71.069204,lat_1=40.,lat_2=70.,lon_1=-70.)

#CHANGE THIS ONE!!!!
fig = plt.figure()
gs = GridSpec(2, 4)
ax = fig.add_subplot(gs[0, 0])
m.imshow(colmean_cosample/1.e15,vmin=0.1,vmax=8,norm=LogNorm(),interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
plt.text(100000,2.8e6,'(a)',fontsize=9)
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0+0.01,0.02,ax.get_position().y1-ax.get_position().y0-0.02])
cbar = fig.colorbar(im,cax=cax,ticks=[0.1,0.5,1,2,3,4,6,8],extend='max')
cbar.ax.set_yticklabels([0.1,0.5,1,2,3,4,6,8])
cbar.set_label(r'VCD [$10^{15}$ molec. cm$^{-2}$]',fontsize=13)
ax = fig.add_subplot(gs[0, 1])
im = m.imshow(no2trans_mod,vmin=0.1,vmax=8,norm=LogNorm(),interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
plt.text(100000,2.8e6,'(b)',fontsize=9)
ax = fig.add_subplot(gs[0:1, :])
m.imshow(colmean_cosample/1.e15 - no2trans_mod,cmap='bwr',vmin=-3,vmax=3,interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
c = plt.colorbar(ticks=np.arange(-3,3.1,1.),extend='both')
c.set_label(r'VCD$_{WRF-Chem}$-VCD$_{OMI}$ [$10^{15}$ molec. cm$^{-2}$]',size=11)
plt.text(100000,2.8e6,'(c)',fontsize=9)
ax = fig.add_subplot(gs[1, 0])
pyplot.hist((colmean_cosample/1.e15)[~np.isnan(var)], bins=60, range=(0,3), alpha=0.5, label='WRF-Chem')
pyplot.hist(no2trans_mod[~np.isnan(var)], bins=60, range=(0,3), alpha=0.5, label='OMI')
pyplot.legend(loc='upper right',fontsize=9)
pyplot.xlabel(r'VCD [$10^{15}$ molec. cm$^{-2}$]',size=11)
pyplot.ylabel('Frequency',size=11)
plt.text(0.1,1110,'(d)',fontsize=9)
ax = fig.add_subplot(gs[1, 1])
plt.hist((var)[~np.isnan(var)],bins=100,range=(-5,5),color='c',edgecolor='k',alpha=0.7)
plt.axvline((var)[~np.isnan(var)].mean(), color='k', linestyle='dashed', linewidth=1)
plt.axvline(np.median((var)[~np.isnan(var)]), color='k', linestyle='dotted', linewidth=1)
plt.axvline(np.nanpercentile(var,5), color='k', linestyle='dashdot', linewidth=1)
plt.axvline(np.nanpercentile(var,95), color='k', linestyle='dashdot', linewidth=1)
plt.axvline(0, color='k', linewidth=2)
plt.xlim((-2,2))
plt.xlabel(r'VCD$_{WRF-Chem}$-VCD$_{OMI}$ [$10^{15}$ molec. cm$^{-2}$]',size=11)
plt.ylabel(' ',size=11)
plt.legend(['Mean','Median','90% conf.'],fontsize=9,loc=2)
plt.text(1.72,1500,'(e)',fontsize=9)
plt.show()

plt.plot()
m.imshow(colmean_cosample/1.e15 - no2trans_mod,cmap='bwr',vmin=-3,vmax=3,interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
c = plt.colorbar(ticks=np.arange(-3,3.1,1.),extend='both')
c.ax.tick_params(labelsize=15) 
c.set_label(r'VCD$_{WRF-Chem}$-VCD$_{OMI}$ [$10^{15}$ molec. cm$^{-2}$]',size=15)
plt.text(100000,2.8e6,'(c)',fontsize=10)
plt.savefig('/home/WUR/barte035/WRFChem/Python/Colombia/FiguresColombia/fig5c_proofread.png',dpi=1000)
plt.show()

print np.nanpercentile(colmean_cosample/1.e15,5)
print np.nanpercentile(colmean_cosample/1.e15,95)
print np.nanpercentile(no2trans_mod,5)
print np.nanpercentile(no2trans_mod,95)
print np.median((colmean_cosample/1.e15)[~np.isnan(var)])
print np.median((no2trans_mod)[~np.isnan(var)])
print r2_score(no2trans_mod[~np.isnan(var)],colmean_cosample[~np.isnan(var)]/1.e15)

print np.nanpercentile(var,5)
print np.nanpercentile(var,95)
print np.median((var)[~np.isnan(var)])
print (var)[~np.isnan(var)].mean()
print np.nanstd(var)


OMIcount = np.zeros((99,99))

for i in range(0,99):
	for j in range(0,99):
		OMIcount[i,j] = np.count_nonzero(~np.isnan(no2stack_mod[:,i,j]),axis=0)

print OMIcount
print np.mean(OMIcount)
print np.max(OMIcount)
print (np.mean(OMIcount)/np.max(OMIcount))*100


plt.figure(figsize=(15,9))
m.imshow(OMIcount/29.*100.,vmin=5,vmax=80,interpolation='none',cmap='viridis')
m.drawcoastlines(zorder=4,color='white')
m.drawcountries(zorder=4,color='white')
m.drawmapboundary(fill_color='lightgray')
m.drawparallels(np.arange(40,61,10),labels=[True,False,True])
m.drawmeridians(np.arange(-10,51,20),labels=[True,False,True])
c = plt.colorbar(extend='both')
#c.set_ticks([2,6,10,14,18,22])
#c.set_ticklabels([2,6,10,14,18,22])
c.set_ticks([5,20,40,60,80])
c.set_ticklabels([5,20,40,60,80])
c.cmap.set_under('red')
c.set_label(r'Available OMI measurements [% of month]',size=14)
#plt.title('Count of filtered OMI-Data',weight='bold')
plt.show()

print(len(OMIcount))
OMIcount2 = OMIcount
OMIcount4 = OMIcount
OMIcount6 = OMIcount
print(len(OMIcount2[OMIcount2 < 1.1]))
print(len(OMIcount2[OMIcount2 < 3.1]))
print(len(OMIcount2[OMIcount2 < 5.1]))


wrf_col_rat = colmean_cosample / 1.e15
omi_col_rat = no2trans_mod


outfile = '/home/WUR/barte035/WRFChem/WRFV3/run/wrfout_d01_2014-01-01_00:00:00'
ds = nc.Dataset(outfile,'r')
landmask = ds.variables['LANDMASK'][0,:,:]
mask = ~np.isnan(wrf_col_rat) & ~np.isnan(omi_col_rat) & (landmask == 1)
wrf_ratio_samp = np.where((landmask == 1), wrf_col_rat, np.nan)
omi_ratio_samp = np.where((landmask == 1), omi_col_rat, np.nan)
geopath = '/home/WUR/barte035/WRFChem/WPS/geo_em.d01.nc'
#Extract WRF-Chem lat/lon and regrid data
geo = nc.Dataset(geopath,'r')
we  = geo.variables['XLAT_M'].shape[2]
sn  = geo.variables['XLAT_M'].shape[1]
lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
dx_dom = 27000

m = Basemap(width=(we*dx_dom)-1,height=(sn*dx_dom)-1,resolution='l',area_thresh=100.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120,projection='lcc',lat_0=4.891087,lon_0=-71.069204,lat_1=40.,lat_2=70.,lon_1=-70.)
plt.figure(figsize=(15,9))
#m = Basemap(width=width,height=height,resolution='l',area_thresh=1000.,projection='lcc',lat_0=lat_0,lon_0=lon_0)
m.imshow(np.where((landmask==1),wrf_ratio_samp/omi_ratio_samp,np.nan),cmap='jet',vmin=0.1,vmax=10,norm = LogNorm(),interpolation='none')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='lightgray')
m.drawparallels(np.arange(40,61,10),labels=[True,False,True])
m.drawmeridians(np.arange(-10,51,20),labels=[True,False,True])
c = plt.colorbar()
c.set_ticks([0.1,0.2,0.25,0.33,0.5,0.75,1,1.33,2,3,4,5,10])
c.set_ticklabels([0.1,0.2,0.25,0.33,0.5,0.75,1,1.33,2,3,4,5,10])
c.set_label(r'NO2 column density [$10^{15}$ molec cm$^{-2}$]')
plt.show()

'''
execfile("budgetsmap.py")

varisnan = var[~np.isnan(var)]

plt.plot()
plt.hist((varisnan)[~np.isnan(NOxMaxLightning)],bins=100,range=(-5,5),color='c',edgecolor='k',alpha=0.7)
plt.axvline((varisnan)[~np.isnan(NOxMaxLightning)].mean(), color='k', linestyle='dashed', linewidth=1)
plt.axvline(np.median((varisnan)[~np.isnan(NOxMaxLightning)]), color='k', linestyle='dotted', linewidth=1)
plt.axvline(np.nanpercentile(varisnan[~np.isnan(NOxMaxLightning)],5), color='k', linestyle='dashdot', linewidth=1)
plt.axvline(np.nanpercentile(varisnan[~np.isnan(NOxMaxLightning)],95), color='k', linestyle='dashdot', linewidth=1)
plt.axvline(0, color='k', linewidth=2)
plt.xlim((-2,2))
plt.xlabel(r'VCD$_{WRF-Chem}$-VCD$_{OMI}$ lightning [$10^{15}$ molec. cm$^{-2}$]',size=11)
plt.ylabel('Frequency',size=11)
plt.legend(['Mean','Median','90% conf.'],fontsize=9,loc=2)
#plt.text(1.72,1500,'(f)',fontsize=9)
plt.show()
'''
