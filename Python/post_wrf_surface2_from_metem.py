import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import wrf
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy
from matplotlib.backends.backend_pdf import PdfPages
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
np.bool = np.bool_

# User settings
home          = os.path.expanduser("~")
datapath      = "%s/atmmod/WPS"%home

(yy,mm,dd,HH) = (2018,6,5,00)

domain        = "d01"
infilename = '%s/met_em.d02.2018-06-05_00:00:00.nc'%datapath
#infilename = '%s/wrfout_%s_%4d-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)

infile = nc.Dataset(infilename, 'r')
name_project = infile.getncattr('TITLE')

XLONG = infile.variables['XLONG_M'][0,:,:]  # Longitude, deg
XLAT  = infile.variables['XLAT_M' ][0,:,:]  # Latitude, deg

ter    = wrf.getvar(infile,'ter',timeidx=0)  
HGT = wrf.to_np(ter)  # Terrain height, m

LAND   = wrf.getvar(infile, 'LU_INDEX', timeidx=0)  # Land use cathegory
SEAICE = wrf.getvar(infile, 'SEAICE',   timeidx=0)  # Sea ice flag

#infile.close()

my_projection = ccrs.PlateCarree(central_longitude=0)

# Plot the global yearly net flux.
with PdfPages('plt_surface2_metem.pdf') as pdf:
    # Figure 1: Terrain height
    fig1 = plt.figure(figsize=(10,8))

    cart_proj = wrf.get_cartopy(ter)
    ax1  = plt.subplot(111, projection=cart_proj)
    ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

    ax1.set_xlim(wrf.cartopy_xlim(ter))
    ax1.set_ylim(wrf.cartopy_ylim(ter))

    gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
    gl.top_labels = gl.right_labels = False

    desc = ter.attrs['description']
    units = ter.attrs['units']
 
    cb = ax1.pcolormesh(XLONG, XLAT, HGT, cmap=plt.cm.RdYlGn,transform=ccrs.PlateCarree()) 
    fig1.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')

    ax1.set_title(name_project[-17:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
    
    pdf.savefig()
    plt.close()
    
    # Figure 2: Land use category

    luts = ['Evergreen needleleaf','Evergreen broadleaf','Deciduous needleleaf', \
             'Deciduous broadleaf','Mixed forest','Closed shrubland','Open shrubland','Woody savanna', \
	     'Savanna','Grassland','Permanent wetland','Cropland','Urban','Cropland/Natural Mosaic', \
	     'Snow and ice','Barren/sparsely veg.','Water','Wooded tundra','Mixed tundra','Barren tundra']
    
    cmap = colors.ListedColormap(["darkgreen", "chartreuse", "olivedrab","palegreen",\
		"darkseagreen","mediumpurple","brown","burlywood","gold",\
		"lightgreen","royalblue","palegoldenrod","red","goldenrod","white",\
		"lightsteelblue","aqua","darkgoldenrod","lightslategray","darkgray"])
 
    bounds = range(1,len(luts)+2)
    norm = colors.BoundaryNorm(bounds,cmap.N)
  

    fig2 = plt.figure(figsize=(10,8))

    cart_proj = wrf.get_cartopy(LAND)
    ax1  = plt.subplot(111, projection=cart_proj) 
    ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

    ax1.set_xlim(wrf.cartopy_xlim(LAND))
    ax1.set_ylim(wrf.cartopy_ylim(LAND))
    
    gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
    gl.top_labels = gl.right_labels = False


    cb = ax1.pcolormesh(XLONG, XLAT, LAND, cmap=cmap,norm=norm,transform=ccrs.PlateCarree()) 
    cbar = fig2.colorbar(cb)
    cbar.ax.get_yaxis().set_ticks([])
    for j,lab in enumerate(luts):
         cbar.ax.text(1.1,(j+.5)/(len(luts)),lab,ha='left',va='center',transform=cbar.ax.transAxes)

 #   fig2.colorbar(cb, orientation='horizontal', pad=0.07, label='Land use')
    desc = LAND.attrs['description']
    units = '-'

    ax1.set_title(name_project[-17:]+': %s \n'%(desc)) # We add a title to the plot 
    pdf.savefig(fig2)
    plt.close()
    
    # Figure 3: Sea ice flag

    seasice = ['0','1']
    
    cmap_seaice = colors.ListedColormap(["white", "black"])
 
    bounds = range(len(seasice)+1)
    norm = colors.BoundaryNorm(bounds,cmap.N)
     
 
    fig3 = plt.figure(figsize=(10,8))

    cart_proj = wrf.get_cartopy(SEAICE)
    ax1  = plt.subplot(111, projection=cart_proj) 
    ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

    ax1.set_xlim(wrf.cartopy_xlim(SEAICE))
    ax1.set_ylim(wrf.cartopy_ylim(SEAICE))
    
    gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
    gl.top_labels = gl.right_labels = False
 
    desc = SEAICE.attrs['description']
    units = SEAICE.attrs['units']

    cb = ax1.pcolormesh(XLONG, XLAT, SEAICE, cmap=cmap_seaice,norm=norm,transform=ccrs.PlateCarree()) 
    cbar = fig3.colorbar(cb)
    cbar.ax.get_yaxis().set_ticks([])
    for j,lab in enumerate(seasice):
         cbar.ax.text(1.1,(j+.5)/(len(seasice)),lab,ha='left',va='center')
    ax1.set_title(name_project[-17:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
    
    pdf.savefig(fig3)
    plt.close()   



	


    
    
    
