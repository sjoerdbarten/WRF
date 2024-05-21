#!/usr/bin/env python

######################    initializing    ########################
from matplotlib.pyplot import *
from StringIO import StringIO
from StringIO import StringIO
from matplotlib.ticker import OldScalarFormatter
from matplotlib.mathtext import *
import pickle 
import datetime
import netCDF4 as nc
from pylab import *
from shutil import copyfile

######################    reading    ###############################
# WRF file
datapathwrf      = '/home/WUR/barte035/WRFChem/WPS'
#filelist = [os.path.join(datapathwrf,filename) for filename in os.listdir(datapathwrf) if filename.startswith('geo_em')]
#filelist = filelist[0:-1]
#filelist.sort()
#print filelist
wrfoldfilename   = '%s/geo_em.d02.nc'%datapathwrf #do this for all your domains CHANGE SAVE FILE IN BOTTOM
wrfoldfile       = nc.Dataset(wrfoldfilename,'r')
wrflat           = wrfoldfile.variables['XLAT_M'][0,:,:]
wrflon           = wrfoldfile.variables['XLONG_M'][0,:,:]
lu_index_wrf     = wrfoldfile.variables['LU_INDEX'][0,:,:]
lu_fraction_wrf  = wrfoldfile.variables['LANDUSEF'][0,:,:]
lu_index_wrf_new = lu_index_wrf
lu_fraction_wrf_new = lu_fraction_wrf
wrfoldfile.close()

nwrflat,nwrflon=wrflat.shape
#nwrflon=size(wrflon)


datapathimage    = '/home/WUR/barte035/WRFChem/IMAGE'
#filelistIMAGE = [os.path.join(datapathimage,filename) for filename in os.listdir(datapathimage) if filename.startswith('GLCT')]
#filelistIMAGE = filelistIMAGE[0:-1]
#filelistIMAGE.sort()
#print filelistIMAGE
imagefilename    = '%s/GLCT.nc'%datapathimage
imagefile        = nc.Dataset(imagefilename,'r')
print imagefile
imglat           = imagefile.variables['latitude'][:]
imglon           = imagefile.variables['longitude'][:]
imgtime          = imagefile.variables['time'][:]
print(imgtime)
# 0 for present-day data (1970), , 6=2000, 8=2010, 16=2050, 26=2100
lu_index_img     = imagefile.variables['GLCT'][8,:,:]
dlat_img         = 1./12.
dlon_img         = 1./12.
imagefile.close()

nimglat=size(imglat)
nimglon=size(imglon)

for iwrflat in range(nwrflat):
    for iwrflon in range(nwrflon):
        ilat=int(nimglat+1-((wrflat[iwrflat,iwrflon]+ 90.)/180.)*nimglat)
        ilon=int(((wrflon[iwrflat,iwrflon]+180.)/360.)*nimglon)
        print('WRF (lon,lat)= (%7.2f,%7.2f); IMG (lon,lat)=(%7.2f,%7.2f)'%(wrflon[iwrflat,iwrflon],wrflat[iwrflat,iwrflon],imglon[ilon],imglat[ilat]))
        lu_img = lu_index_img[ilat,ilon]

        # resetting the land use fraction
        lu_fraction_wrf_new[:,iwrflat,iwrflon]=0.

        # IMAGE class #1: agricultural land -> USGS class 2: dryland cropland and pasture
        if lu_img == 1:
           lu_index_wrf_new[iwrflat,iwrflon] = 2
        # IMAGE class #2: extensive grassland -> USGS class 7: grassland
        elif lu_img == 2:
           lu_index_wrf_new[iwrflat,iwrflon] = 7
        # IMAGE class #3: carbon plantation -> USGS class 6: cropland/woodland mosaic
        elif lu_img == 3:
           lu_index_wrf_new[iwrflat,iwrflon] = 6    
        # IMAGE class #4: regrowth forest abandoning -> USGS class 6: cropland/woodland mosaic
        elif lu_img == 4:
           lu_index_wrf_new[iwrflat,iwrflon] = 6
        # IMAGE class #5: regrowth forest timber -> USGS class 6: cropland/woodland mosaic
        elif lu_img == 5:
           lu_index_wrf_new[iwrflat,iwrflon] = 6
        # IMAGE class #6: biofuels -> USGS class 6: cropland/woodland mosaic
        elif lu_img == 6:
           lu_index_wrf_new[iwrflat,iwrflon] = 6    
        # IMAGE class #7: ice -> USGS class 24: snow or ice
        elif lu_img == 7:
           lu_index_wrf_new[iwrflat,iwrflon] = 24
        # IMAGE class #8: tundra -> USGS class 22: bare ground tundra
        elif lu_img == 8:
           lu_index_wrf_new[iwrflat,iwrflon] = 23    
        # IMAGE class #9: wooded tundra -> USGS class 21: wooded tundra 
        elif lu_img == 9:
           lu_index_wrf_new[iwrflat,iwrflon] = 21    
        # IMAGE class #10: boreal forest -> USGS class 15: mixed forest (alternatively evergreen needleleaf??)
        elif lu_img == 10:
           lu_index_wrf_new[iwrflat,iwrflon] = 15
        # IMAGE class #11: cool coniferous forest -> USGS class 14: evergreen needleleaf
        elif lu_img == 14:
           lu_index_wrf_new[iwrflat,iwrflon] = 15
        # IMAGE class #12: temperate mixed foresy -> USGS class 15: mixed forest (alternatively evergreen needleleaf??)
        elif lu_img == 12:
           lu_index_wrf_new[iwrflat,iwrflon] = 15
        # IMAGE class #13: temperate decideous forest -> USGS class 15: mixed forest (alternatively evergreen needleleaf??)
        elif lu_img == 13:
           lu_index_wrf_new[iwrflat,iwrflon] = 15
        # IMAGE class #14: warm mixed forest -> USGS class 13: Evergreen broadleaf
        elif lu_img == 14:
           lu_index_wrf_new[iwrflat,iwrflon] = 13
        # IMAGE class #15: grassland/steppe -> USGS class 7: grassland
        elif lu_img == 15:
           lu_index_wrf_new[iwrflat,iwrflon] = 7
        # IMAGE class #16: hot desert -> USGS class 23: bare ground tundra (?)
        elif lu_img == 16:
           lu_index_wrf_new[iwrflat,iwrflon] = 23
        # IMAGE class #17: scrubland -> USGS class 8: shrubland
        elif lu_img == 17:
           lu_index_wrf_new[iwrflat,iwrflon] = 8
        # IMAGE class #18: savanna -> USGS class 10: savanna
        elif lu_img == 18:
           lu_index_wrf_new[iwrflat,iwrflon] = 10
        # IMAGE class #19: tropical woodland -> USGS class 10: savanna
        elif lu_img == 19:
           lu_index_wrf_new[iwrflat,iwrflon] = 10
        # IMAGE class #20: tropical forest -> USGS class 13: evergreen broadleaf (or dedideous broadleaf forest?)
        elif lu_img == 20:
           lu_index_wrf_new[iwrflat,iwrflon] = 13
        # IMAGE class #21: water -> USGS class 16: water bodies
        elif lu_img == 21:
           lu_index_wrf_new[iwrflat,iwrflon] = 16

        # resetting the land use fraction
        lu_fraction_wrf_new[lu_index_wrf_new[iwrflat,iwrflon],iwrflat,iwrflon]=1.
        #print(lu_fraction_wrf_new[:,iwrflat,iwrflon])
wrfnewfilename = '%s/geo_em.d02.IMAGE2010.nc'%datapathwrf
copyfile(wrfoldfilename,wrfnewfilename)


wrfnewfile     = nc.Dataset(wrfnewfilename,'r+')
print(wrfnewfilename)

wrfnewfile.variables['LU_INDEX'][0,:,:]   = lu_index_wrf_new
wrfnewfile.variables['LANDUSEF'][0,:,:,:] = lu_fraction_wrf_new

print 'Conversion complete'

wrfnewfile.close()

