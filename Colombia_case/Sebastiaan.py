import netCDF4 as nc
import matplotlib.pyplot as plt

def timeseries3D(varname,xcoor,ycoor,zcoor,datapath):
	ncfile = nc.Dataset(datapath+'wrfout_d01_2016-12-20_00:00:00','r')
	
	print(ncfile.variables['T'].shape)
	
	plotvar = ncfile.variables[varname][:,zcoor,xcoor,ycoor]
	
	print(plotvar.shape)
	
	plt.plot(plotvar)
	plt.show()
