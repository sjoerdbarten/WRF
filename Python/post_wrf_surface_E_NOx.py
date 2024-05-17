import os
import numpy as np
import matplotlib.pyplot as plt
import wrf
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
np.bool = np.bool_

# User settings
home          = os.path.expanduser("~") # set home directory
datapath      = "%s/atmmod/WRF/run"%home # set full datapath

# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,00)

# which domain is selected

domain        = "d02"

#Select the time instances

# step time in timesteps of wrfout 
# (i.e . if value is 6 every sixth time stamp in wrfout is plotted)

step_times = 1 

# User settings

# Create custom color map similar to the NCAR NCL WhiteBlueGreenYellowRed
cdict = {'red':   ((   0/253., 255./255., 255./255.), (  36/253., 157./255., 157./255.),
                   (  72/253.,  72./255.,  72./255.), ( 108/253.,  73./255.,  73./255.),
                   ( 145/253., 250./255., 250./255.), ( 181/253., 245./255., 245./255.),
                   ( 217/253., 211./255., 211./255.), ( 253/253., 146./255., 146./255.)),
         'green': ((   0/253., 255./255., 255./255.), (  36/253., 218./255., 218./255.),
                   (  72/253., 142./255., 142./255.), ( 108/253., 181./255., 181./255.),
                   ( 145/253., 232./255., 232./255.), ( 181/253., 106./255., 106./255.),
                   ( 217/253.,  31./255.,  31./255.), ( 253/253.,  21./255.,  21./255.)),
         'blue':  ((   0/253., 255./255., 255./255.), (  36/253., 247./255., 247./255.),
                   (  72/253., 202./255., 202./255.), ( 108/253.,  70./255.,  70./255.),
                   ( 145/253.,  92./255.,  92./255.), ( 181/253.,  45./255.,  45./255.),
                   ( 217/253.,  40./255.,  40./255.), ( 253/253.,  25./255.,  25./255.))}

my_cmap = LinearSegmentedColormap('my_colormap', cdict, 256)

infilename = '%s/wrfout_%s_%4d-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)

infile = nc.Dataset(infilename, 'r')

name_project = infile.getncattr('TITLE')
times = wrf.getvar(infile,"Times",timeidx=wrf.ALL_TIMES)
plot_times = range(1,len(times),step_times)


XLONG = infile.variables['XLONG'][0,:,:]  # Longitude, deg
XLAT  = infile.variables['XLAT' ][0,:,:]  # Latitude, deg
#infile.close()



#Extract data and Plot

inittime = times[0]
init_str = np.datetime_as_string(inittime.values,unit='s',timezone='UTC')
uneven_levels = [1.00e-01, 1.78e-01, 3.16e-01, 5.62e-01, 1.00e0, 1.78e0, 3.16e0, 5.62e0, 1.00e1, 1.78e1, 3.16e1, 5.62e1, 1.00e2, 1.78e2, 3.16e2, 5.62e2, 1.00e3, 1.78e3, 3.16e3, 5.62e3]


with PdfPages('plt_surface_E_NOx.pdf') as pdf:
    for it in plot_times: # range(1, nmodeltimesteps):

        curtime = times[it]
        cur_str = np.datetime_as_string(curtime.values,unit='s',timezone='UTC')
	
        print('Processing time: ',cur_str)
 
        u     = wrf.getvar(infile,'ua',    it)  
        v     = wrf.getvar(infile,'va',    it)
        E_NOx = wrf.getvar(infile,'E_NOx', it)
        NOx   = wrf.getvar(infile,'NOx',   it)
        
        NOx_mugm3 = wrf.to_np(NOx) *1000. # *1.35
	
        # NOx emissions at surface

        fig1 = plt.figure(figsize=(10,8))
        cart_proj = wrf.get_cartopy(E_NOx)
	
        ax1  = plt.subplot(111, projection=cart_proj)
       
        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

        ax1.set_xlim(wrf.cartopy_xlim(E_NOx))
        ax1.set_ylim(wrf.cartopy_ylim(E_NOx))
       
        colors = my_cmap(np.linspace(0, 1, len(uneven_levels)-1))
        my_cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)
        cb = ax1.pcolormesh(XLONG, XLAT, E_NOx[0,:,:], cmap=my_cmap, norm=norm,transform=ccrs.PlateCarree()) 
       
        desc = E_NOx.attrs['description']
        units = E_NOx.attrs['units']

        colbar = fig1.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')
        colbar.ax.set_xticklabels(["{:.3}".format(i) for i in colbar.get_ticks()]) # set ticks of your format

        ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+') at lowest model layer')) # We add a title to the plot 
      
        ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)

        pdf.savefig()
        plt.close()
	
        # NOx mixing ratio at surface
        fig2 = plt.figure(figsize=(10,8))
        cart_proj = wrf.get_cartopy(NOx)
	
        ax1  = plt.subplot(111, projection=cart_proj)
	
        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

        ax1.set_xlim(wrf.cartopy_xlim(NOx))
        ax1.set_ylim(wrf.cartopy_ylim(NOx))
      #  levels = np.linspace(0,100,10) 
        cb = ax1.contourf(XLONG, XLAT, NOx_mugm3[0,:,:],cmap=plt.cm.viridis,transform=ccrs.PlateCarree(),extend='max') 
 
        desc = NOx.attrs['description']
        units = 'ppb'

        fig2.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')

        nq = 7
        plt.barbs(XLONG[::nq, ::nq], XLAT[::nq, ::nq], u.values[0,::nq,::nq], v.values[0,::nq,::nq], \
        pivot='tip', sizes=dict(emptybarb=0.01, spacing=0.2, height=0.3), linewidth=0.3, \
	transform=ccrs.PlateCarree())

        ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+') at surface')) # We add a title to the plot 
         

        ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)

        pdf.savefig()
        plt.close()


    
    
    
