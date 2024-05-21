import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm 

ncfile = nc.Dataset('/home/WUR/barte035/v432_CH4_2010.0.1x0.1.nc','r')

CH4 = ncfile.variables['emi_ch4'][:,:]
CH4year = CH4*60.*60.*24.*365. #convert kg m-2 s-1 to kg m-2 yr-1

print(CH4)
print(CH4.shape)

latmin = 850
latmax = 1046
lonmin = 2800
lonmax = 3000

plt.figure(figsize=(12,8))
m = Basemap(projection='cyl',llcrnrlon=-80,llcrnrlat=-5,urcrnrlon=-60,urcrnrlat=14.6)
#m.imshow(CH4year[latmin:latmax,lonmin:lonmax],norm=LogNorm(vmin=0.000001,vmax=0.01),interpolation=None,cmap='Greens')
m.imshow(CH4year[850:1046,2800:3000],vmin=0.,vmax=0.01,interpolation=None,cmap='Greens')
m.drawcoastlines()
m.drawcountries()
c = plt.colorbar(extend='max')
c.set_label(r'CH$_{4}$ Flux [kg m$^{-2}$ year$^{-1}$]',size=14)
m.drawparallels(np.arange(-20,30,5),labels=[True,False,True])
m.drawmeridians(np.arange(-90,-50,5),labels=[True,False,True])
plt.show()

surfacecell = 11000.*11000. #m2
CH4flux = CH4year[latmin:latmax,lonmin:lonmax]*surfacecell #kg cell-1 yr-1
CH4flux = np.sum(CH4flux)*1e-9 #kg yr-1 whole domain
print(CH4flux)
