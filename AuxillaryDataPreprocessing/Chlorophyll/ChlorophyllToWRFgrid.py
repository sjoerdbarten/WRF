import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time

#Load data
chlorophyllpath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/Chlorophyll/aug_A20022132019243.L3m_MC_CHL_chlor_a_9km.nc'
Chlorophyll = nc.Dataset(chlorophyllpath,'r')
chlorophyll_data = Chlorophyll.variables['chlor_a'][:,:]

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

#Select correct datapoints  ([2160:4320])
chlorophyll_data = np.delete(chlorophyll_data,np.s_[540:],axis=0)[::-1,:]
chlorophyll_data = np.where(chlorophyll_data != -32767, chlorophyll_data, np.nan)

print(chlorophyll_data.shape)
print(np.arange(-180.,180.,1./12.).shape)
print(np.arange(45.,90.,1./12.).shape)

chlorophyll_trans = m.transform_scalar(chlorophyll_data,lons=np.arange(-180.,180.,1./12.),lats=np.arange(45.,90.,1./12.),nx=249,ny=249,order=1)

print(chlorophyll_trans[:,:].shape)

im = m.imshow(chlorophyll_trans,cmap='brg',vmin=0,vmax=30)
cbar = m.colorbar(im,extend='max')
cbar.set_label("Seawater chlorophyll concentration [mg m$^{-3}$]",fontsize=18)
cbar.ax.tick_params(labelsize=16)
m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
plt.savefig('ChlorophyllMapAugust.png',dpi=300)
time.sleep(5)
plt.show()

chlorophyll_trans = np.where(np.isnan(chlorophyll_trans), 0.1, chlorophyll_trans)

'''
#Add to wrflowinput
wrflowinputpath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/wrflowinp_d01_test'
wrflowinput = nc.Dataset(wrflowinputpath,'r+',format='NETCDF4')
print(wrflowinput.dimensions.values())
chlorophyll = wrflowinput.createVariable('CHL_OCEAN','f4',('Time','south_north','west_east'))
chlorophyll.units = 'mg m^-3'
for i in range(chlorophyll.shape[0]):
	chlorophyll[i,:,:] = chlorophyll_trans[:,:]
	print(i)
print(chlorophyll.shape)
print(chlorophyll_trans.shape)
wrflowinput.close()
'''

'''
#Add to wrfinput
wrfinputpath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/wrfinput_d01_test'
wrfinput = nc.Dataset(wrfinputpath,'r+',format='NETCDF4')
print(wrfinput.dimensions.values())
chlorophyll = wrfinput.createVariable('CHL_OCEAN','f4',('Time','south_north','west_east'))
chlorophyll.FieldType = int(104)
chlorophyll.MemoryOrder = "XY "
chlorophyll.description = "[Chl a] oceanic concentration"
chlorophyll.units = 'mg m^-3'
chlorophyll.stagger = ""
chlorophyll.coordinates = "XLONG XLAT XTIME"
chlorophyll[0,:,:] = chlorophyll_trans[:,:]
print(chlorophyll.shape)
print(chlorophyll_trans.shape)
wrfinput.close()
'''
