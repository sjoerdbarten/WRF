import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time
import pandas as pd

#Load data
dmspath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/DMS/DMSclim_AUG.csv'
DMS = pd.read_csv(dmspath,header=None)
dms_data = DMS.values

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
dms_data = np.delete(dms_data,np.s_[45:],axis=0)[::-1,:]

print(dms_data.shape)
print(np.arange(-180.,180.,1.).shape)
print(np.arange(45.,90.,1.).shape)

dms_trans = m.transform_scalar(dms_data,lons=np.arange(-180.,180.,1.),lats=np.arange(45.,90.,1.),nx=249,ny=249,order=1)

print(dms_trans[:,:].shape)

im = m.imshow(dms_trans,cmap='brg',vmin=0,vmax=10)
cbar = m.colorbar(im,extend='max')
cbar.set_label("Seawater DMS concentration [nM]",fontsize=18)
cbar.ax.tick_params(labelsize=16)
m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
plt.savefig('DMSMapAugust.png',dpi=300)
time.sleep(5)
plt.show()

dms_trans = np.where(np.isnan(dms_trans), 0.1, dms_trans)

'''
#Add to wrflowinput
wrflowinputpath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/wrflowinp_d01_test'
wrflowinput = nc.Dataset(wrflowinputpath,'r+',format='NETCDF4')
print(wrflowinput.dimensions.values())
dms = wrflowinput.createVariable('DMS_OCEAN','f4',('Time','south_north','west_east'))
dms.units = 'nM'
for i in range(dms.shape[0]):
	dms[i,:,:] = dms_trans[:,:]
	print(i)
print(dms.shape)
print(dms_trans.shape)
wrflowinput.close()
'''

#Add to wrfinput
wrfinputpath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/wrfinput_d01_test'
wrfinput = nc.Dataset(wrfinputpath,'r+',format='NETCDF4')
print(wrfinput.dimensions.values())
dms = wrfinput.createVariable('DMS_OCEAN','f4',('Time','south_north','west_east'))
dms.FieldType = int(104)
dms.MemoryOrder = "XY "
dms.description = "[DMS] oceanic concentration"
dms.units = 'nM'
dms.stagger = ""
dms.coordinates = "XLONG XLAT XTIME"
dms[0,:,:] = dms_trans[:,:]
print(dms.shape)
print(dms_trans.shape)
wrfinput.close()
