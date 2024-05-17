import os
import numpy as np
import wrf
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.backends.backend_pdf import PdfPages
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
np.bool = np.bool_

##############
#User settings
##############
home = os.path.expanduser("~") # set home directory
datapath = "%s/atmmod/WRF/run"%home # set full datapath
print(datapath)
# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,00)

# which domain is selected

domain     = 'd02'

#Uncomment the variable to visualize below
#isfilevar = 'QVAPOR'
isfilevar = 'QCLOUD'
#isfilevar = 'QRAIN'
#isfilevar = 'QICE'     #Not in the sea breeze output file
#isfilevar = 'QGRAUPEL' #Not in the sea breeze output file

#Select the time instances and model level to visualize
#Select the time step for plotting

# step time in timesteps of wrfout 
# (i.e . if value is 6 every sixth time stamp in wrfout is plotted)

step_times = 12 # step time in timesteps of wrfout

#Select the model level to visualize. Lowest model level is default

plot_level = 0  #Select the model level to visualize

########################################################
#Processing (no need to change anything below this line)
########################################################

infilename = '%s/wrfout_%s_%i-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)

infile = nc.Dataset(infilename, "r")

name_project = infile.getncattr('TITLE')
times = wrf.getvar(infile,"Times",timeidx=wrf.ALL_TIMES)
plot_times = range(0,len(times),step_times)

XLONG = infile.variables['XLONG'][0,:,:]  # Longitude, deg
XLAT  = infile.variables['XLAT' ][0,:,:]  # Latitude, deg

inittime = times[0]
init_str = np.datetime_as_string(inittime.values,unit='s',timezone='UTC')

#Extract data and Plot

with PdfPages('plt_cloud_lev%i.pdf'%plot_level) as pdf:
    for t in plot_times:
        
        curtime = times[t]
        cur_str = np.datetime_as_string(curtime.values,unit='s',timezone='UTC')
	
        print('Processing time: ',cur_str)
        
        var = wrf.getvar(infile, isfilevar, timeidx=t)
        cart_proj = wrf.get_cartopy(var)
        my_projection = cart_proj
	
        var_np = wrf.to_np(var)[plot_level,:,:]
        var_np *= 1000.
        
        fig = plt.figure(figsize=(10,8))
        
        ax1 = plt.subplot(111, projection=cart_proj)
        
        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linewidth=0.2)
	
        ax1.set_xlim(wrf.cartopy_xlim(var))
        ax1.set_ylim(wrf.cartopy_ylim(var))
	
        gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
        gl.top_labels = gl.right_labels = False
	
        desc = var.attrs['description']
        units = 'g/kg'
        
        color_levels = np.arange(0., 1., 0.1)
 
        cb = ax1.contourf(XLONG, XLAT, var_np, color_levels, transform=ccrs.PlateCarree(), cmap=plt.cm.Blues)

        fig.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')', ticks=color_levels)
        
        ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
 
        ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)

        pdf.savefig(fig)
        plt.close()
	
infile.close()
