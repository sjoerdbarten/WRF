import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time

#Load data
nitratepath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/Nitrate/woa18_all_n08_01.nc'
Nitrate = nc.Dataset(nitratepath,'r')
nitrate_data = Nitrate.variables['n_an'][0,0,:,:]

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
m = Basemap(resolution='l',area_thresh=1000.,projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0.)

#Select correct datapoints  ([2160:4320])
nitrate_data = np.delete(nitrate_data,np.s_[0:135],axis=0)

print('NITRATE DATA SHAPE',nitrate_data.shape)
print(np.arange(-180.,180.,1.).shape)
print(np.arange(45.,90.,1.).shape)

nitrate_data = np.where(nitrate_data == 9.96921e+36, np.nan, nitrate_data)

nitrate_trans = m.transform_scalar(nitrate_data,lons=np.arange(-180.,180.,1.),lats=np.arange(45.,90.,1.),nx=249,ny=249,order=1)

print(nitrate_trans[:,:].shape)

m.imshow(nitrate_trans,vmin=0,vmax=10,cmap='viridis_r')
m.colorbar(label='Seawater NO$_3$$^-$ concentration [mmol kg$^{-1}$]',extend='max')
m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
plt.savefig('NitrateMapAugust.png',dpi=300)
time.sleep(5)
plt.show()

nitrate_trans = np.where(np.isnan(nitrate_trans), 0.001, nitrate_trans)

#Add to wrfinput
wrfinputpath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/wrfinput_d01_test'
wrfinput = nc.Dataset(wrfinputpath,'r+',format='NETCDF4')
print(wrfinput.dimensions.values())
nitrate = wrfinput.createVariable('Nitrate_bilinear','f4',('Time','south_north','west_east'))
nitrate.units = 'mmol kg^-1'
nitrate[0,:,:] = nitrate_trans[:,:]
print(nitrate.shape)
print(nitrate_trans.shape)
wrfinput.close()
