import os
import numpy as np              # Numpy is the fundamental package for
import netCDF4 as nc            # NetCDF is the data format of the 
import matplotlib.pyplot as plt # Matplotlib is a scientific plotting package.
import datetime                 # Datetime is a package to deal with dates.
import cartopy.crs as ccrs
import cartopy
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import wrf
import pwd
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

domain        = "d01"

#Select the time instances

# step time in timesteps of wrfout 
# (i.e . if value is 6 every sixth time stamp in wrfout is plotted)

step_times = 6


thr = 750.

########################################################
#Processing (no need to change anything below this line)
########################################################

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

inittime = times[0]
init_str = np.datetime_as_string(inittime.values,unit='s',timezone='UTC')

u_cart = wrf.getvar(infile,'ua',0)

# load variables from ERA5    
with PdfPages('plt_IVT.pdf') as pdf:
    # Figure 1: 2m air temperature

    for it in plot_times[:]:
        
        curtime = times[it]
        cur_str = np.datetime_as_string(curtime.values,unit='s',timezone='UTC')

        print('Processing time: ',cur_str)

        v     = wrf.to_np(wrf.getvar(infile,"va",it))
        u     = wrf.to_np(wrf.getvar(infile,"ua",it))
        q     = wrf.to_np(wrf.getvar(infile,"QVAPOR",it))
        p     = wrf.to_np(wrf.getvar(infile,'pres',it))
        
        bar = np.empty(u.shape[0]-1)
        bar_u = np.empty(u.shape[0]-1)
        bar_v = np.empty(v.shape[0]-1)

        VT = np.empty(u.shape[0])
        VT_u = np.empty(u.shape[0])
        VT_v = np.empty(v.shape[0])

        IVT = np.empty((u.shape[1],u.shape[2]))
        IVT_u = np.empty((u.shape[1],u.shape[2]))
        IVT_v = np.empty((v.shape[1],v.shape[2]))
	
        for x in range(u.shape[2]):
         for y in range(u.shape[1]):
           for z in range(u.shape[0]):
                    VT[z] = np.sqrt(np.square(q[z,y,x]*v[z,y,x])+np.square(q[z,y,x]*u[z,y,x]))
                    VT_u[z] = q[z,y,x]*u[z,y,x] 
                    VT_v[z] = q[z,y,x]*v[z,y,x]
		    
           for z in range(u.shape[0]-1):
                    bar[z] = ((VT[z]+VT[z+1])/2)*(p[z+1,y,x]-p[z,y,x])
                    bar_u[z] = ((VT_u[z]+VT_u[z+1])/2)*(p[z+1,y,x]-p[z,y,x])
                    bar_v[z] = ((VT_v[z]+VT_v[z+1])/2)*(p[z+1,y,x]-p[z,y,x])
           IVT[y,x]   = -sum(bar)/9.81
           IVT_u[y,x] = -sum(bar_u)/9.81
           IVT_v[y,x] = -sum(bar_v)/9.81
	   
         #  print(IVT[y,x])
        IVT_mask = (IVT >= thr)*1.
        fig5 = plt.figure(figsize=(10,8))
        cart_proj = wrf.get_cartopy(u_cart)
	
        ax1  = plt.subplot(111, projection=cart_proj)
	
        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

        ax1.set_xlim(wrf.cartopy_xlim(u_cart))
        ax1.set_ylim(wrf.cartopy_ylim(u_cart))

        gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
        gl.top_labels = gl.right_labels = False
       
        cb = ax1.pcolormesh(XLONG, XLAT, IVT, transform=ccrs.PlateCarree(),cmap=my_cmap) 
        ax1.quiver(XLONG,XLAT,IVT_u,IVT_v,transform=ccrs.PlateCarree())
        ax1.contour(XLONG,XLAT,IVT_mask,transform=ccrs.PlateCarree())
        desc = 'IVT'
        units = 'kg m-1 s-1'
        fig5.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')
        ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
        ax1.text(0,1,'Init: '+init_str,va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(1,1,'valid: '+cur_str,ha='right',va='bottom',fontsize=10,transform=ax1.transAxes)
        ax1.text(0,1,'Contouring: IVT => '+str(thr)+' kg m-1 s-1',va='top',fontsize=10,transform=ax1.transAxes)
	
        pdf.savefig()
        plt.close()
    
infile.close()
        
        
