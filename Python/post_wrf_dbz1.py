import os
import numpy as np
import wrf
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
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
datapath = "%s/atmmod/WRF/run/"%home # set full datapath


# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,00)

# which domain is selected

domain     = 'd02'

#Select the time instances and model level to visualize

# step time in timesteps of wrfout 
# (i.e . if value is 6 every sixth time stamp in wrfout is plotted)

step_times = 1 

#Select the model level to visualize. For reflectivity the lowest level (0) is default.

plot_level = 0  #Select the model level to visualize. For reflectivity the lowest level (0) is default.

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

#Extract data and Plot

inittime = times[0]
init_str = np.datetime_as_string(inittime.values,unit='s',timezone='UTC')

with PdfPages('plt_dbz_lev%i.pdf'%plot_level) as pdf:
    
    for t in plot_times:

        curtime = times[t]
        cur_str = np.datetime_as_string(curtime.values,unit='s',timezone='UTC')
	
        print('Processing time: ',cur_str)
       
        dbz = wrf.getvar(infile, 'dbz', timeidx=t)
        dbz_np = wrf.to_np(dbz)[plot_level,:,:]
  	
        mdbz = wrf.getvar(infile, 'mdbz', timeidx=t)
        mdbz_np = wrf.to_np(mdbz)
	
        #FIGURE 1: Reflectivity
        fig = plt.figure(figsize=(10,8))
        cart_proj = wrf.get_cartopy(dbz)
	
        ax1 = plt.subplot(111, projection=cart_proj)
        
        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linewidth=0.2)
        
        ax1.set_xlim(wrf.cartopy_xlim(dbz))
        ax1.set_ylim(wrf.cartopy_ylim(dbz))
	
        gl = ax1.gridlines(draw_labels=True,
                   linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
        gl.top_labels = gl.right_labels = False
                
        cb = ax1.contourf(XLONG, XLAT, dbz_np, [5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75.], transform=ccrs.PlateCarree(), cmap=plt.cm.viridis, extend='max')
        
        desc = dbz.attrs['description']
        units = dbz.attrs['units']

        fig.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')
        
        ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
 
 # add title to plot and info on init time and the time of validaity of the forecast
  
        
            
        ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)

 
 # Save the figure
 
        pdf.savefig(fig)
        plt.close()
        
        #FIGURE 2: Maximum reflectivity	
        fig = plt.figure(figsize=(10,8))
        cart_proj = wrf.get_cartopy(mdbz)

        ax1 = plt.subplot(111, projection=cart_proj)
        
        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linewidth=0.2)

        ax1.set_xlim(wrf.cartopy_xlim(mdbz))
        ax1.set_ylim(wrf.cartopy_ylim(mdbz))
	
        gl = ax1.gridlines(draw_labels=True,
                   linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
        gl.top_labels = gl.right_labels = False
    
        cb = ax1.contourf(XLONG, XLAT, mdbz_np, [5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75.], transform=ccrs.PlateCarree(), cmap=plt.cm.viridis, extend='max')
        desc = mdbz.attrs['description']
        units = mdbz.attrs['units']

        fig.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')
        
        ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
         #fig.tight_layout() # With this function we make the figure fit as good as possible.
        ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)

        pdf.savefig(fig)
        plt.close()
        
infile.close()
