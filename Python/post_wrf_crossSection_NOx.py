import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import wrf
import cartopy.crs as ccrs
import cartopy
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
np.bool = np.bool_

##############
#User settings
##############

home       = os.path.expanduser('~')    # set home directory 
datapath   = '%s/atmmod/WRF/run'%home # set full datapath

# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,0)

# which domain is selected
domain = "d01"

#Select the time instances and model level to visualize

# step time in timesteps of wrfout 
# (i.e . if value is 6 every sixth time stamp in wrfout is plotted)
 
step_times = 1 

# Select the coordinates of the starting and end point of the cross
# section

cross_start = wrf.CoordPair(lat=52, lon=3) # start coordinates
cross_end = wrf.CoordPair(lat=52, lon=9.5) # end coordinates

# Select the height for which vertical profile is determined

bot_height = 0. # lowest height
top_height = 5000. # top height
step_size = 10. # height step

########################################################
#Processing (no need to change anything below this line)
########################################################

infilename = '%s/wrfout_%s_%4d-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)

print('opening %s ...'%infilename)
infile     = nc.Dataset(infilename,'r')
name_project = infile.getncattr('TITLE')
times = wrf.getvar(infile,"Times",timeidx=wrf.ALL_TIMES)  # get times in the file


XLONG = infile.variables['XLONG'][0,:,:]  # Longitude, deg
XLAT  = infile.variables['XLAT' ][0,:,:]  # Latitude, deg

plot_times = range(0,len(times),step_times)

inittime = times[0]
init_str = np.datetime_as_string(inittime.values,unit='s',timezone='UTC')
ter = wrf.getvar(infile, "ter", timeidx=0)
levels = np.arange(bot_height,top_height,step_size)

with PdfPages('plt_crossSection_NOx.pdf') as pdf:
   
   fig = plt.figure(figsize=(10,8))
   cart_proj = wrf.get_cartopy(ter)
	
   ax1 = plt.subplot(111, projection=cart_proj)
        
   ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
   ax1.add_feature(cartopy.feature.BORDERS, linewidth=0.2)
        
   ax1.set_xlim(wrf.cartopy_xlim(ter))
   ax1.set_ylim(wrf.cartopy_ylim(ter))
   
   gl = ax1.gridlines(draw_labels=True,
                   linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
   gl.top_labels = gl.right_labels = False
                
   cb = ax1.contourf(XLONG, XLAT, ter,transform=ccrs.PlateCarree(), cmap=plt.cm.viridis, extend='max')
   
   var = wrf.getvar(infile,'NOx',timeidx=0)
   ht = wrf.getvar(infile, "z", timeidx=0)    
   z_cross = wrf.vertcross(var, ht, levels=levels,wrfin=infile,
                    start_point=cross_start,
                    end_point=cross_end,
                    latlon=True, meta=True)
   
   for elem,elemp1 in zip(z_cross.coords['xy_loc'][0:-1],z_cross.coords['xy_loc'][1:]):
      lat = wrf.to_np(elem).tolist().lat
      lon = wrf.to_np(elem).tolist().lon
      latp1 = wrf.to_np(elemp1).tolist().lat
      lonp1 = wrf.to_np(elemp1).tolist().lon
         
      ax1.plot([lon,lonp1], [lat, latp1],
          color='red', linestyle='-',lw=3,
         transform=ccrs.PlateCarree())
      ax1.set_title(name_project[-21:]+': crossSection transect')
   pdf.savefig()    
   
   for tracer in ['NOx']: 
      print('Processing variable: '+tracer)
      for t  in plot_times[1:]: # range(0,len(times),6): #,time_str in enumerate(times):

        curtime = times[t]
        cur_str = np.datetime_as_string(curtime.values,unit='s',timezone='UTC')
	
        print('Processing time: ',cur_str)
      
         
        ht = wrf.getvar(infile, "z", timeidx=t)
        var = wrf.getvar(infile,tracer,timeidx=t)
        desc = var.attrs['description']
        units = 'ppbv'
     
        z_cross = wrf.vertcross(var*1000., ht, levels=levels,wrfin=infile,
                    start_point=cross_start,
                    end_point=cross_end,
                    latlon=True, meta=True)

        z_cross_filled = np.ma.copy(wrf.to_np(z_cross))

# For each cross section column, find the first index with non-missing
# values and copy these to the missing elements below.


        for i in range(z_cross_filled.shape[-1]):
           column_vals = z_cross_filled[:,i]
           first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
           z_cross_filled[:first_idx, i] = z_cross_filled[first_idx, i]


        ter_line = wrf.interpline(ter, wrfin=infile, start_point=cross_start,
                      end_point=cross_end)

        fig = plt.figure(figsize=(10,8))
        ax_cross = plt.axes()

        xs = np.arange(0, z_cross.shape[-1], 1)
        ys = wrf.to_np(z_cross.coords["vertical"])
        z_contours = ax_cross.contourf(xs,
                                 ys,
                                 wrf.to_np(z_cross_filled))

        cb_z = fig.colorbar(z_contours,orientation="horizontal",pad=0.1)
        cb_z.set_label(label=desc+' ('+units+')')

        ht_fill = ax_cross.fill_between(xs, 0, wrf.to_np(ter_line),
                                facecolor="saddlebrown")
				
        coord_pairs = wrf.to_np(z_cross.coords["xy_loc"])
        x_ticks = np.arange(coord_pairs.shape[0])
        x_labels = [pair.latlon_str() for pair in wrf.to_np(coord_pairs)]
        num_ticks = 5
        thin = int((len(x_ticks) / num_ticks) + .5)
        ax_cross.set_xticks(x_ticks[::thin])
        ax_cross.set_xticklabels(x_labels[::thin], rotation=0)

# Set the x-axis and  y-axis labels
        ax_cross.set_ylim([ys[0],ys[-1]])
        ax_cross.set_xlabel("Latitude, Longitude")
        ax_cross.set_ylabel("Height (m)")
        ax_cross.set_title(name_project[-21:]+': '+desc+' ('+units+') \n')
        ax_cross.text(0,1,'Init: '+init_str,va='bottom',transform=ax_cross.transAxes)
        ax_cross.text(1,1,'valid: '+cur_str,ha='right',va='bottom',transform=ax_cross.transAxes)
 
        pdf.savefig()
        plt.close()

infile.close()
