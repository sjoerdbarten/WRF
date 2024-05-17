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

# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,00)

# which domain is selected

domain     = 'd02'


#Select the time instances and pressure levels to visualize

# step time in timesteps of wrfout 
# (i.e . if value is 6 every sixth time stamp in wrfout is plotted)

step_times = 6 # step time in timesteps of wrfout

#Select the variable and the pressure level

var_pres = {'avo':500.0,'pvo':300.0}

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

with PdfPages('plt_pvo_avo.pdf') as pdf:
    for t in plot_times:


       curtime = times[t]
       cur_str = np.datetime_as_string(curtime.values,unit='s',timezone='UTC')
	
       print('Processing time: ',cur_str)
    
       for isfilevar in var_pres: 
         
        levp = var_pres[isfilevar]
        print('Processing: '+isfilevar+' on level: '+str(levp))
        
           
        var_get = wrf.getvar(infile, isfilevar, timeidx=t)
        
        pres = wrf.getvar(infile,'pressure', timeidx=t)
        var_int = wrf.interplevel(var_get, pres, levp)
	
	#var = wrf.to_np(wrf.getvar(infile, isfilevar, timeidx=t)[plot_level,:,:])
        var = wrf.to_np(var_int)
	
        cart_proj = wrf.get_cartopy(var_get)
        
        fig = plt.figure(figsize=(10,8))
        ax1 = plt.subplot(111, projection=cart_proj)
        
        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linewidth=0.2)
	

        ax1.set_xlim(wrf.cartopy_xlim(var_get))
        ax1.set_ylim(wrf.cartopy_ylim(var_get))
 
        gl = ax1.gridlines(draw_labels=True,
                   linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
        gl.top_labels = gl.right_labels = False
  
        cb = ax1.contourf(XLONG, XLAT, var, 12, transform=ccrs.PlateCarree(), cmap=plt.cm.Blues, extend='max')
 
        desc = var_get.attrs['description']
        units = var_get.attrs['units']

        fig.colorbar(cb, orientation='horizontal', pad=0.07, label='%s (%s)'%(isfilevar,units))
  
        ax1.set_title(name_project[-21:]+': %s (%s) on pressure level %s (hPa) \n' %(isfilevar,units,levp)) # We add a title to the plot 
        ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)
      
        pdf.savefig(fig)
        plt.close()
	
infile.close()
