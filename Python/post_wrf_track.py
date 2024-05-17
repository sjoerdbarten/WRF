#********************************************************
# Plot storm stracks from wrfout files.
#********************************************************

import os
import numpy as np
import matplotlib.pyplot as plt
import wrf
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
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

########################################################
#Processing (no need to change anything below this line)
########################################################

# Create custom color map similar to the NCAR NCL WhiteBlueGreenYellowRed
infilename = '%s/wrfout_%s_%i-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)

infile = nc.Dataset(infilename, 'r')

name_project = infile.getncattr('TITLE')
times = wrf.getvar(infile,"Times",timeidx=wrf.ALL_TIMES)
plot_times = range(0,len(times),step_times)


XLONG = infile.variables['XLONG'][0,:,:]  # Longitude, deg
XLAT  = infile.variables['XLAT' ][0,:,:]  # Latitude, deg

#Extract data and Plot

inittime = times[0]
init_str = np.datetime_as_string(inittime.values,unit='s',timezone='UTC')

lines = []
slp = wrf.getvar(infile,'slp', 0) 

fig = plt.figure(figsize=(10,8))    
cart_proj = wrf.get_cartopy(slp)
ax1  = plt.subplot(111, projection=cart_proj)

ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)
    
ax1.set_xlim(wrf.cartopy_xlim(slp))
ax1.set_ylim(wrf.cartopy_ylim(slp))

gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
gl.top_labels = gl.right_labels = False
   
for it in plot_times:
        
    curtime = times[it]
    cur_str = np.datetime_as_string(curtime.values,unit='s')
    hour = np.datetime64(cur_str).astype(object).hour
    
    print('Processing time: ',cur_str)
 
    slp = wrf.getvar(infile,'slp', it) 
    slp_oned = wrf.to_np(slp).flatten()
    imin = np.argmin(slp_oned)
    iy,ix = np.unravel_index((imin), slp.shape)
    
    
    lines.append(','.join((cur_str,'{:8.2f}'.format(wrf.to_np(slp)[iy,ix]),  \
                            '{:8d}'.format(ix),'{:8d}'.format(iy), \
			    '{:8.2f}'.format(XLONG[iy,ix]), \
			    '{:8.2f}'.format(XLAT[iy,ix]))))

    if it > plot_times[0]:
      if hour in [0,6,12,18]: 
         ax1.plot([xlongm1,XLONG[iy,ix]],[xlatm1,XLAT[iy,ix]],"k-",transform = ccrs.PlateCarree())
         print(hour,XLONG[iy,ix],XLAT[iy,ix])
         ax1.text(XLONG[iy,ix], XLAT[iy,ix], cur_str+': {:8.2f}'.format(wrf.to_np(slp)[iy,ix]),transform=ccrs.PlateCarree()) 
         xlongm1 = XLONG[iy,ix]
         xlatm1 = XLAT[iy,ix]    
    else:
       xlongm1 = XLONG[iy,ix]
       xlatm1 = XLAT[iy,ix]
       
plt.savefig("plt_track.pdf",format="pdf")    

outF = open("minslp.txt", "w")
outF.write("date/time,pressure (hPa),ix,iy,lon,lat")
for line in lines:
  # write line to output file
  outF.write(line)
  outF.write("\n")
outF.close()
     
infile.close()
