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
datapath      = "%s/atmmod/WRF/run"%home

(yy,mm,dd,HH) = (2018,6,5,00)

domain        = "d02"

#infilename = '/lustre/shared/atmos2020/tarfile/wrfout_d02_2018-06-05_00:00:00'
infilename = '%s/wrfout_%s_%4d-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)

infile = nc.Dataset(infilename, 'r')
name_project = infile.getncattr('TITLE')

XLONG = infile.variables['XLONG'][0,:,:]  # Longitude, deg
XLAT  = infile.variables['XLAT' ][0,:,:]  # Latitude, deg

ter    = wrf.getvar(infile,'ter',timeidx=0)  
HGT = wrf.to_np(ter)  # Terrain height, m
LU = wrf.getvar(infile,'LU_INDEX',timeidx=0)
LAND   = wrf.to_np(LU)  # Land use cathegory

soiltype = wrf.getvar(infile,'ISLTYP',timeidx=0)
ISOIL  = wrf.to_np(soiltype)  # Dominant soil cathegory

frac_veg = wrf.getvar(infile,'VEGFRA',timeidx=0)
VEGFRA = wrf.to_np(frac_veg)  # Vegetation fraction

tsb = wrf.getvar(infile, 'TSLB',timeidx=0)
TSOIL = wrf.to_np(tsb)  # Soil temperature, K
sms = wrf.getvar(infile,'SMOIS',timeidx=0)
SMOIS = wrf.to_np(sms)  # Soil moisture, m^3/m^3

#infile.close()


# Plot the global yearly net flux.
with PdfPages('plt_surface2.pdf') as pdf:
    # Figure 1: Terrain height
   
    print('Processing variable: '+ter.attrs['description'])
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

    ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
    
    pdf.savefig()
    plt.close()

    # Figure 2: Land use category
    print('Processing variable: '+LU.attrs['description'])
    
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

    cart_proj = wrf.get_cartopy(LU)
    ax1  = plt.subplot(111, projection=cart_proj) 
    ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

    ax1.set_xlim(wrf.cartopy_xlim(LU))
    ax1.set_ylim(wrf.cartopy_ylim(LU))
    
    gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
    gl.top_labels = gl.right_labels = False

 
    cb = ax1.pcolormesh(XLONG, XLAT, LAND, cmap=cmap,norm=norm,transform=ccrs.PlateCarree()) 
    cbar = fig2.colorbar(cb)
    cbar.ax.get_yaxis().set_ticks([])
    for j,lab in enumerate(luts):
        cbar.ax.text(1.1,(j+.5)/(len(luts)),lab,ha='left',va='center',transform=cbar.ax.transAxes)

 #   fig2.colorbar(cb, orientation='horizontal', pad=0.07, label='Land use')
    desc = LU.attrs['description']

    ax1.set_title(name_project[-21:]+': %s \n'%(desc)) # We add a title to the plot 
    pdf.savefig(fig2)
    plt.close()

    # Figure 3: Dominant soil category
    print('Processing variable: '+soiltype.attrs['description'])
    soilts = ['Sand','Loamy sand','Sandy loam', \
             'Silt loam','Silt','Loam','Sandy clay loam','Silty clay loam', \
	     'Clay loam','Sandy clay','Silty clay','Clay','Organic matter','Water', \
	     'Bedrock','Other','Playa','Lava','White sand']
    
    cmap_soil = colors.ListedColormap(["yellow", "chartreuse", "olivedrab","palegreen",\
		"darkseagreen","mediumpurple","brown","burlywood","gold",\
		"lightgreen","royalblue","palegoldenrod","goldenrod","aqua","lightslategray",\
		"lightsteelblue","red","darkgoldenrod","white"])
    bounds = range(1,len(soilts)+2)
    norm = colors.BoundaryNorm(bounds,cmap_soil.N)
  

    fig3 = plt.figure(figsize=(10,8))
    cart_proj = wrf.get_cartopy(soiltype)
    
    ax1  = plt.subplot(111, projection=cart_proj)
    
    ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

    ax1.set_xlim(wrf.cartopy_xlim(soiltype))
    ax1.set_ylim(wrf.cartopy_ylim(soiltype))
    
    gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
    gl.top_labels = gl.right_labels = False
 
    cb = ax1.pcolormesh(XLONG, XLAT, ISOIL, cmap=cmap_soil,norm=norm,transform=ccrs.PlateCarree()) 
   
    cbar = fig3.colorbar(cb)
    cbar.ax.get_yaxis().set_ticks([])
    for j,lab in enumerate(soilts):
         cbar.ax.text(1.1,(j+.5)/(len(soilts)),lab,ha='left',va='center',transform=cbar.ax.transAxes)

 #   fig2.colorbar(cb, orientation='horizontal', pad=0.07, label='Land use')
    desc = soiltype.attrs['description']

    ax1.set_title(name_project[-21:]+': '+desc) # We add a title to the plot 
    pdf.savefig(fig3)
    plt.close()

    # Figure 4: Vegetation fraction
    print('Processing variable: '+frac_veg.attrs['description'])

    fig4 = plt.figure(figsize=(10,8))
    cart_proj = wrf.get_cartopy(frac_veg)
    
    ax1  = plt.subplot(111, projection=cart_proj)

    ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

    ax1.set_xlim(wrf.cartopy_xlim(frac_veg))
    ax1.set_ylim(wrf.cartopy_ylim(frac_veg))
    
    gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
    gl.top_labels = gl.right_labels = False

    desc = frac_veg.attrs['description']
    units = '%'
 
    cb = ax1.pcolormesh(XLONG, XLAT, VEGFRA, cmap=plt.cm.RdYlGn,transform=ccrs.PlateCarree()) 
 
  
    fig4.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')

    ax1.set_title(name_project[-21:]+': %s \n'%(desc+' ('+units+')')) # We add a title to the plot 
    
    pdf.savefig(fig4)
    plt.close()

    # Figure 5: Soil temperature
    for it in range(4):
    
        print('Processing variable: '+tsb.attrs['description']+' at level',str(it))
  
        fig5 = plt.figure(figsize=(10,8))
        cart_proj = wrf.get_cartopy(tsb)
	
        ax1  = plt.subplot(111, projection=cart_proj)
        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

        ax1.set_xlim(wrf.cartopy_xlim(tsb))
        ax1.set_ylim(wrf.cartopy_ylim(tsb))
    
        gl = ax1.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline=False)
        gl.top_labels = gl.right_labels = False

        desc = tsb.attrs['description']
        units = tsb.attrs['units']
 
        cb = ax1.pcolormesh(XLONG, XLAT, TSOIL[it,:,:], vmin=285, vmax=295, cmap=plt.cm.RdYlGn,transform=ccrs.PlateCarree()) 
        
        fig5.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')

        ax1.set_title(name_project[-21:]+': %s at soil level %s \n'%(desc+' ('+units+')',str(it+1))) # We add a title to the plot 

        pdf.savefig(fig5)
        plt.close()
	
    # Figure 6: Soil moisture
    for it in range(4):
       
        print('Processing variable: '+sms.attrs['description']+' at level',str(it))
 
        fig6 = plt.figure(figsize=(10,8))
        cart_proj = wrf.get_cartopy(sms)
	
        ax1  = plt.subplot(111, projection=cart_proj)
        ax1.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

        ax1.set_xlim(wrf.cartopy_xlim(sms))
        ax1.set_ylim(wrf.cartopy_ylim(sms))	
     

        cb = ax1.pcolormesh(XLONG, XLAT, SMOIS[it,:,:], vmin=0, vmax=0.4, cmap=plt.cm.coolwarm,transform=ccrs.PlateCarree()) 

        desc = sms.attrs['description']
        units = sms.attrs['units']

        fig6.colorbar(cb, orientation='horizontal', pad=0.07, label=desc+' ('+units+')')

        ax1.set_title(name_project[-21:]+': %s at soil level %s \n'%(desc+' ('+units+')',str(it+1))) # We add a title to the plot 
              
        pdf.savefig(fig6)
        plt.close()


    
    
    
