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
datapath = "%s/atmmod/WRF/run"%home # set full datapath

# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,0)

# which domain is selected

domain     = 'd02'

# select timestep for which cape needs to be calculated (-1 is last timestep)

timestep = -1

########################################################
#Processing (no need to change anything below this line)
########################################################

infilename = '%s/wrfout_%s_%i-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)
infile = nc.Dataset(infilename, "r")
name_project = infile.getncattr('TITLE')

XLONG = infile.variables['XLONG'][0,:,:]  # Longitude, deg
XLAT  = infile.variables['XLAT' ][0,:,:]  # Latitude, deg

#Extract data

cape_3d = wrf.getvar(infile, 'cape_3d', timeidx=timestep)
cape  = wrf.to_np(cape_3d)[0,:,:,:]
cin   = wrf.to_np(cape_3d)[1,:,:,:]

cape_2d = wrf.getvar(infile, 'cape_2d', timeidx=timestep)
mcape = wrf.to_np(cape_2d)[0,:,:]
mcin  = wrf.to_np(cape_2d)[1,:,:]
lcl   = wrf.to_np(cape_2d)[2,:,:]
lfc   = wrf.to_np(cape_2d)[3,:,:]

inittime = wrf.getvar(infile,"Times",timeidx=0)
curtime = wrf.getvar(infile,"Times",timeidx=timestep)
init_str = np.datetime_as_string(inittime.values,unit='s',timezone='UTC')
cur_str = np.datetime_as_string(curtime.values,unit='s',timezone='UTC')    
 
#######
#Output
#######

with PdfPages('plt_cape.pdf') as pdf:
    
    #Figure 1: MCAPE
    fig1 = plt.figure(figsize=(10,8))
    
    cart_proj = wrf.get_cartopy(cape_2d)
    
    ax1 = plt.subplot(111, projection=cart_proj)
    ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)
    
    ax1.set_xlim(wrf.cartopy_xlim(cape_2d))
    ax1.set_ylim(wrf.cartopy_ylim(cape_2d))
    
    gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
    gl.top_labels = gl.right_labels = False

    cb = ax1.contourf(XLONG, XLAT, mcape, 10, vmin=0, vmax=3000, transform=ccrs.PlateCarree(), cmap=plt.cm.RdYlGn)

    desc = cape_2d.attrs['description'].split(';')[0]
    units = cape_2d.attrs['units'].split(';')[0]
    
    fig1.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')
        
    ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
    
    ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
    ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)
    
    pdf.savefig()
    plt.close()
    
    #Figure 2: MCIN
    fig2 = plt.figure(figsize=(10,8))

    cart_proj = wrf.get_cartopy(cape_2d)

    ax1 = plt.subplot(111, projection=cart_proj)
    ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)
    
    ax1.set_xlim(wrf.cartopy_xlim(cape_2d))
    ax1.set_ylim(wrf.cartopy_ylim(cape_2d))
    
    gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
    gl.top_labels = gl.right_labels = False
    
    cb = ax1.contourf(XLONG, XLAT, mcin, 10, vmin=0, vmax=300, transform=ccrs.PlateCarree(), cmap=plt.cm.RdYlGn)

    desc = cape_2d.attrs['description'].split(';')[1]
    units = cape_2d.attrs['units'].split(';')[1]
    
    fig2.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')

    ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
 
    #fig2.tight_layout() # With this function we make the figure fit as good as possible.
    ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
    ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)
 
    pdf.savefig(fig2)
    plt.close()
    
    #Figure 3: Level of Free Convection
    fig3 = plt.figure(figsize=(10,8))
 
    cart_proj = wrf.get_cartopy(cape_2d)
    ax1 = plt.subplot(111, projection=cart_proj)

    ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

    ax1.set_xlim(wrf.cartopy_xlim(cape_2d))
    ax1.set_ylim(wrf.cartopy_ylim(cape_2d))
    gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
    gl.top_labels = gl.right_labels = False
    
    cb = ax1.contourf(XLONG, XLAT, lfc, levels=np.arange(200.,4000.,400.), transform=ccrs.PlateCarree(), cmap=plt.cm.RdYlGn, extend='max')

    desc = cape_2d.attrs['description'].split(';')[2]
    units = cape_2d.attrs['units'].split(';')[2]
    
    fig3.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')

    ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 

    ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
    ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)
 
    pdf.savefig(fig3)
    plt.close()  
    
    #Figure 4: Lifting-Condensation level
    fig4 = plt.figure(figsize=(10,8))

    cart_proj = wrf.get_cartopy(cape_2d)
    ax1 = plt.subplot(111, projection=cart_proj)

    ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

    ax1.set_xlim(wrf.cartopy_xlim(cape_2d))
    ax1.set_ylim(wrf.cartopy_ylim(cape_2d))
    
    gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
    gl.top_labels = gl.right_labels = False
    
    cb = ax1.contourf(XLONG, XLAT, lcl, levels=np.arange(200.,4000.,400.), transform=ccrs.PlateCarree(), cmap=plt.cm.RdYlGn, extend='max')

    desc = cape_2d.attrs['description'].split(';')[3]
    units = cape_2d.attrs['units'].split(';')[3]
    
    fig4.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')

    ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
    
    ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
    ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)
 
    pdf.savefig(fig4)
    plt.close()  
    
    #Figs. 5-14: CAPE in lowermost 10 layers
    for lev in range(10):
        fig = plt.figure(figsize=(10,8))

        cart_proj = wrf.get_cartopy(cape_3d)
        ax1 = plt.subplot(111, projection=cart_proj)

        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)
        
        ax1.set_xlim(wrf.cartopy_xlim(cape_3d))
        ax1.set_ylim(wrf.cartopy_ylim(cape_3d))
        gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
        gl.top_labels = gl.right_labels = False
   
        cb = ax1.contourf(XLONG, XLAT, cape[lev,:,:], levels=np.arange(0.,3300.,300.), transform=ccrs.PlateCarree(), cmap=plt.cm.RdYlGn, extend='max')

        desc = cape_3d.attrs['description'].split(';')[0]
        units = cape_3d.attrs['units'].split(';')[0]
        
        fig.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')

        ax1.set_title(name_project[-21:]+': %s on level %s \n'%(desc+' ('+units+')',str(lev))) # We add a title to the plot 
        
        ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)
 
        pdf.savefig(fig)
        plt.close()
	
    #Figs. 15-24: CAPE in lowermost 10 layers
    for lev in range(10):
        fig = plt.figure(figsize=(10,8))
  
        cart_proj = wrf.get_cartopy(cape_3d)
        ax1 = plt.subplot(111, projection=cart_proj)
	
        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)
        
        ax1.set_xlim(wrf.cartopy_xlim(cape_3d))
        ax1.set_ylim(wrf.cartopy_ylim(cape_3d))
        gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
        gl.top_labels = gl.right_labels = False

        cb = ax1.contourf(XLONG, XLAT, cin[lev,:,:], levels=np.arange(0.,330.,30.), transform=ccrs.PlateCarree(), cmap=plt.cm.RdYlGn, extend='max')
	
        desc = cape_3d.attrs['description'].split(';')[1]
        units = cape_3d.attrs['units'].split(';')[1]
        
        fig.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')

        ax1.set_title(name_project[-21:]+': %s on level %s \n'%(desc+' ('+units+')',str(lev))) # We add a title to the plot 

        ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)
 
        pdf.savefig(fig)
        plt.close()
    
infile.close()

