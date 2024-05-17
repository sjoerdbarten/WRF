import os
import numpy as np
import matplotlib.pyplot as plt
import wrf
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import metpy.calc as mpcalc
from metpy.units import units
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
np.bool = np.bool_

##############
#User settings
##############
home          = os.path.expanduser("~") # set home directory
datapath      = "%s/atmmod/WRF/run"%home # set full datapath

# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,00)

# which domain is selected

domain        = "d01"

#Select the time instances

# step time in timesteps of wrfout 
# (i.e . if value is 6 every sixth time stamp in wrfout is plotted)

step_times = 1

eta_level = 3 						#Select for which eta level you want to plot the divergence
color_levels_divergence = np.arange(-0.5,0.5,0.01)	#Color levels divergence Red White Blue

########################################################
#Processing (no need to change anything below this line)
########################################################

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

infilename = '%s/wrfout_%s_%i-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)

infile = nc.Dataset(infilename, 'r')

name_project = infile.getncattr('TITLE')
times = wrf.getvar(infile,"Times",timeidx=wrf.ALL_TIMES)
plot_times = range(0,len(times),step_times)


XLONG = infile.variables['XLONG'][0,:,:]  # Longitude, deg
XLAT  = infile.variables['XLAT' ][0,:,:]  # Latitude, deg

dx_wrf = infile.getncattr('DX')/1000.		#Your model grid spacing in km
dy_wrf = infile.getncattr('DY')/1000.		#Your model grid spacing in km

#Extract data and Plot

inittime = times[0]
init_str = np.datetime_as_string(inittime.values,unit='s',timezone='UTC')

slp_levels = np.arange(900.,1100.,2.)

#ESG_SB_20230411+ Load parameters needed for calculation divergence
uvmet = wrf.getvar(infile,"uvmet",timeidx=wrf.ALL_TIMES)
U = uvmet.values[0,:,eta_level,:,:]
V = uvmet.values[1,:,eta_level,:,:]
divergence = mpcalc.divergence(U, V, dx=dx_wrf*units.km, dy=dy_wrf*units.km)
#ESG_SB_20230411- Load parameters needed for calculation divergence

with PdfPages('plt_divergence.pdf') as pdf:
    # Figure 1: Divergence

    for it in plot_times:
        
        curtime = times[it]
        cur_str = np.datetime_as_string(curtime.values,unit='s',timezone='UTC')
	
        print('Processing time: ',cur_str)
 
        SLP = wrf.getvar(infile,'slp', it)  
        SLP = wrf.smooth2d(SLP, 3)	#smooth slp

        T2  = wrf.getvar(infile,'T2',  it)
        u10 = wrf.getvar(infile,'U10', it)
        v10 = wrf.getvar(infile,'V10', it)
        rh2 = wrf.getvar(infile,'rh2', it)
        
        fig5 = plt.figure(figsize=(10,8))
        cart_proj = wrf.get_cartopy(SLP)
	
        ax1  = plt.subplot(111, projection=cart_proj)
	
        ax1.add_feature(cartopy.feature.COASTLINE.with_scale('10m'), linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS.with_scale('10m'), linestyle='-', linewidth=.2)

        ax1.set_xlim(wrf.cartopy_xlim(SLP))
        ax1.set_ylim(wrf.cartopy_ylim(SLP))

        gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
        gl.top_labels = gl.right_labels = False
        
        #ESG_SB_20230411+ Changed plotting of divergence
        cb = ax1.contourf(XLONG, XLAT, wrf.to_np(divergence[it,:,:]), color_levels_divergence,
				transform=ccrs.PlateCarree(), cmap=plt.cm.seismic, extend='both')
        
        desc = 'WRF Divergence, eta_level = '+str(eta_level)
        units = '$10^{-5} s^{-1}$'
        #ESG_SB_20230411- Changed plotting of divergence
	
        fig5.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')
	
        cz = ax1.contour(XLONG, XLAT, wrf.to_np(SLP), slp_levels, colors='k', linewidths=3,transform=ccrs.PlateCarree())
        ax1.clabel(cz, fmt='%1.0f', fontsize=10.)
        nq = 7
        plt.barbs(XLONG[::nq, ::nq], XLAT[::nq, ::nq], u10.values[::nq,::nq], v10.values[::nq,::nq], \
	pivot='tip', sizes=dict(emptybarb=0.01, spacing=0.2, height=0.3), linewidth=0.3,transform=ccrs.PlateCarree())
 
        ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
        
        ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)

        pdf.savefig()
        plt.close()
infile.close()
