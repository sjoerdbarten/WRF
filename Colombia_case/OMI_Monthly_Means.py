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

path = '/archive/ESG/barte035/Colombia/OMI/MonthlyMeans/'
files = sorted([f for f in listdir(path) if isfile(join(path,f))])
print(files)
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
	
	plt.figure(figsize=(15,10))
    	m = Basemap(width=(we*dx_dom)-1,height=(sn*dx_dom)-1,resolution='l',area_thresh=100.,projection='lcc',lat_0=4.891087,lon_0=-71.069204,lat_1=40.,lat_2=70.,lon_1=-70.,llcrnrlon=-79.684021,llcrnrlat=-4.9757309,urcrnrlon=-60.025818,urcrnrlat=14.642120)
	nx = 99
	ny = 99
	x2_trans = np.linspace(-179.9375,179.9375,2880)
	y2_trans = np.linspace(-89.9375,89.9375,1440)
	
	value = np.where(value < 0, np.nan, value)

	#value_trans = m.transform_scalar(value[::-1,:],lons=np.arange(-82.375,-60.,0.125),lats=np.arange(-5.,14.625,0.125),nx=100.,ny=100.,order=1)
	value_trans = m.transform_scalar(value[::-1],lons=x2_trans,lats=y2_trans,nx=nx,ny=ny,order=0)
	np.where(value_trans < 0, np.nan, value_trans)
	
	print(value_trans.shape)
	
	#lons=np.arange(-180,180.1,0.1),lats=np.arange(60,90.1,0.1)
	m.imshow(value_trans)
	m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	m.drawparallels(np.arange(-90.,-50,10.),labels=[False,False,False,False],zorder=1001)
	m.drawmeridians(np.arange(-10,20,5.),labels=[True,False,False,True],zorder=1001)
	plt.show()
