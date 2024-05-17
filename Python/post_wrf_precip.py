import os
import numpy as np
import wrf
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy
from precip_paris import paris_precip_data
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
datapath = "%s/atmmod/Python"%home # set full datapath
 
 # for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2022,5,6,00)

# which domain is selected

domain     = 'd03'

# step time in timesteps of wrfout 
# (i.e . if value is 6 every sixth time stamp in wrfout is plotted)

step_times = 12 

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

i = 0

mask_ar, precip_ar, lon_ar, lat_ar = paris_precip_data()

with PdfPages('plt_precip.pdf') as pdf:
    for t in plot_times:

        curtime = times[t]
        cur_str = np.datetime_as_string(curtime.values,unit='s',timezone='UTC')
	
        print('Processing time: ',cur_str)
        
        slp = wrf.getvar(infile, "slp", timeidx=t)
        
        slp_np = wrf.to_np(slp)
        slp_np = wrf.smooth2d(slp_np, 3)	#smooth slp
        rain_exp = wrf.to_np(wrf.getvar(infile, "RAINNC", timeidx=t))
        rain_con = wrf.to_np(wrf.getvar(infile, "RAINC", timeidx=t))
        rain_tot = rain_exp + rain_con
	
        cart_proj = wrf.get_cartopy(slp)
        
        if t == plot_times[0]:
            rain_exp_save = rain_exp
            rain_con_save = rain_con
            rain_tot_save = rain_tot
        else:
            rain_exp_save = wrf.to_np(wrf.getvar(infile, "RAINNC", timeidx=t-1))
            rain_con_save = wrf.to_np(wrf.getvar(infile, "RAINC", timeidx=t-1))
            rain_tot_save = rain_exp_save + rain_con_save
            times_sav = plot_times[i-1]
        
        i = i + 1
	
        rain_exp_tend = rain_exp - rain_exp_save
        rain_con_tend = rain_con - rain_con_save
        rain_tot_tend = rain_tot - rain_tot_save
	
        # Bookkeeping, just to allow the tendency at the next time step
        rain_exp_save = rain_exp
        rain_con_save = rain_con
        rain_tot_save = rain_tot
	
        color_levels = np.arange(1.,15.,1.)
        cmap = plt.cm.Oranges
        norm = mpl.colors.BoundaryNorm(color_levels, cmap.N)
        slp_levels = np.arange(900.,1100.,2.)
	
        if t != plot_times[0]: #We will skip the first time    
            #FIGURE 1: Total precipitation   
            fig = plt.figure(figsize=(10,8))
            ax1 = plt.subplot(111, projection=cart_proj)
            
            ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
            ax1.add_feature(cartopy.feature.BORDERS, linewidth=0.2)
            
            ax1.set_xlim(wrf.cartopy_xlim(slp))
            ax1.set_ylim(wrf.cartopy_ylim(slp))
	
            gl = ax1.gridlines(draw_labels=True,
                   linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
            gl.top_labels = gl.right_labels = False
           
            cb = ax1.contourf(XLONG, XLAT, rain_tot, color_levels, transform=ccrs.PlateCarree(), cmap=plt.cm.Oranges, extend='max')
            for p in range(precip_ar.shape[0]):
                        if mask_ar[p] == False:
                                    ax1.scatter(lon_ar[p], lat_ar[p], 30, c=precip_ar[p], transform=ccrs.PlateCarree(), cmap=cmap, norm=norm, marker='o', linewidth=0.5, edgecolor='black',zorder=10+int(precip_ar[p]))
            fig.colorbar(cb, orientation='horizontal', pad=0.07, label='Total precipitation (mm)')
            ax1.set_title(name_project[-21:]+': Total precipitation (mm) \n') # We add a title to the plot 
            #fig.tight_layout() # With this function we make the figure fit as good as possible.
            ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
            ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)

            pdf.savefig(fig)
            plt.close()
	    
            #FIGURE 2: Precipitation tendency
            fig = plt.figure(figsize=(10,8))
            ax1 = plt.subplot(111, projection=cart_proj)
            
            ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
            ax1.add_feature(cartopy.feature.BORDERS, linewidth=0.2)
	    
            ax1.set_xlim(wrf.cartopy_xlim(slp))
            ax1.set_ylim(wrf.cartopy_ylim(slp))
	    
            gl = ax1.gridlines(draw_labels=True,
                   linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
            gl.top_labels = gl.right_labels = False
            
            
            cb = ax1.contourf(XLONG, XLAT, rain_tot_tend, color_levels, transform=ccrs.PlateCarree(), cmap=plt.cm.Oranges, extend='max')
            fig.colorbar(cb, orientation='horizontal', pad=0.07, label='Precipitation tendency (mm)')
            cs = ax1.contour(XLONG, XLAT, slp_np, slp_levels, transform=ccrs.PlateCarree(), linewidths=3, colors='blue')
            ax1.clabel(cs, inline=1, fontsize=10)
            ax1.set_title(name_project[-21:]+': Precipitation tendency (mm) \n') # We add a title to the plot 
	    
            ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
            ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)

            pdf.savefig(fig)
            plt.close()
            
            #FIGURE 3: Explicit precipitation tendency   
            fig = plt.figure(figsize=(10,8))
            ax1 = plt.subplot(111, projection=cart_proj)
            
            ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
            ax1.add_feature(cartopy.feature.BORDERS, linewidth=0.2)

            ax1.set_xlim(wrf.cartopy_xlim(slp))
            ax1.set_ylim(wrf.cartopy_ylim(slp))
	
            gl = ax1.gridlines(draw_labels=True,
                   linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
            gl.top_labels = gl.right_labels = False
            
            
            cb = ax1.contourf(XLONG, XLAT, rain_exp_tend, color_levels, transform=ccrs.PlateCarree(), cmap=plt.cm.Oranges, extend='max')
            fig.colorbar(cb, orientation='horizontal', pad=0.07, label='Explicit precipitation tendency (mm)')
            cs = ax1.contour(XLONG, XLAT, rain_con_tend, color_levels, transform=ccrs.PlateCarree(), linewidths=3, colors='red')
            ax1.clabel(cs, inline=1, fontsize=10)
            ax1.set_title(name_project[-21:]+': Explicit tendency (mm) \n') # We add a title to the plot 
            ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
            ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)

            pdf.savefig(fig)
            plt.close()

infile.close()
