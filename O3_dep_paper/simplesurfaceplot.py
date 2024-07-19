from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

def simplesurfaceplot():
	wrfnetcdf = nc.Dataset('/archive/ESG/barte035/nov19era5_d01_2019-11-01_00:00:00','r')
	print(wrfnetcdf.variables)
	sfcp = wrfnetcdf.variables['PSFC'][0,:,:] #Surface pressure timestep 0
	
	m1 = Basemap(projection='npstere',boundinglat=62.5,lon_0=0.,resolution='l')
	m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
	lons,lats = m1(wrfnetcdf.variables['XLONG'][0,:,:],wrfnetcdf.variables['XLAT'][0,:,:],inverse=False)
	m1.contour(lons,lats,wrfnetcdf.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
	m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
	m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
	im = m1.imshow(sfcp[:,:],cmap='jet',interpolation='none',zorder=1000)
	plt.show()

simplesurfaceplot()
