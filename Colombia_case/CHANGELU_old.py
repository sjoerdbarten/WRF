import netCDF4 as nc
import numpy as np

geopath = 'D:'
geofilename = '%s/geo_em.d02.nc'%geopath
geofile = nc.Dataset(geofilename,'r+')
lu_index_wrf = geofile.variables['LU_INDEX'][:]
I = np.where(lu_index_wrf == 10)
lu_index_wrf[I] = 2
geofile.variables['LU_INDEX'][:] = lu_index_wrf
geofile.close()
print 'conversion complete'
