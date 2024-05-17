import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

#Load data
pco2path = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/spco2_MPI_SOM-FFN_v2018.nc'
pco2 = nc.Dataset(pco2path,'r')
print(pco2.variables['spco2_raw'].shape)
pco2_data = pco2.variables['spco2_raw'][-1,:,:] #BE SURE TO SELECT THE RIGHT DATE (Data is monthly, from January 1982 to December 2017)
print(pco2_data.shape)

#Define basemap
geopath = '/home/WUR/barte035/WRFChem/WPS/geo_em.d01.nc'
geo = nc.Dataset(geopath,'r')
wrflat = geo.variables['XLAT_M'][0,:,:]
wrflon = geo.variables['XLONG_M'][0,:,:]
we  = geo.variables['XLAT_M'].shape[2]
sn  = geo.variables['XLAT_M'].shape[1]
lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
dx_dom = 25000

plt.figure(figsize=(15,10))
m = Basemap(width=(dx_dom*(we-10))-1,height=(dx_dom*(sn-10))-1,resolution='l',area_thresh=1000.,projection='npstere',boundinglat=65.8,lat_0=90.,lon_0=0.)

#Select correct datapoints
pco2_data = np.delete(pco2_data,np.s_[0:150],axis=0)

print('aaaaa',pco2_data.shape)
print(np.arange(-180.,180.,1.).shape)
print(np.arange(60.,90.,1.).shape)

#pco2_trans = m.transform_scalar(pco2_data[::-1,:][(60+90)*(10*1.25):(90+90)*(10*1.25),(-180+180)*(10*1.25):(180+180)*(10*1.25)],lons=np.arange(-180,180,0.125),lats=np.arange(60,90,0.125),nx=220.,ny=220.,order=0)
pco2_trans = m.transform_scalar(pco2_data,lons=np.arange(-180.,180.,1.),lats=np.arange(60.,90.,1.),nx=219,ny=219,order=1)

print(pco2_trans[:,0])

#lons=np.arange(-180,180.1,0.1),lats=np.arange(60,90.1,0.1)
m.imshow(pco2_trans)
m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
plt.show()


#Add to wrfinput
wrfinputpath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/wrfinput_d01_test'
wrfinput = nc.Dataset(wrfinputpath,'r+',format='NETCDF4')
print(wrfinput.dimensions.values())
pco2ocean = wrfinput.createVariable('PCO2_bilinear','f4',('Time','south_north','west_east'))
pco2ocean.units = 'muatm'
pco2ocean[0,:,:] = pco2_trans[:,:]
print(pco2ocean.shape)
print(pco2_trans.shape)
#pco2ocean[0,:,:] = 380.
wrfinput.close()

#print(wrfinput)
