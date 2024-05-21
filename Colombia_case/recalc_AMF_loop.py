""" recalc_AMF_loop.py

This script takes:
- a WRF-Chem output file
- a directory with OMI files

calculates:
- the lat/lon range of the WRF-Chem output file
- which OMI files (corresponding to orbits) have lat/lons in this range
- per relevant OMI orbit, re-calculate the OMI AMFs for the pixels that
  correspond to a WRF-Chem cell by calling OMI_recalc_AMF.py. This re-
  quires netCDF datasets from WRF-Chem and OMI and WRF-Chem lat/lon in-
  formation

and adds the following variables to files in the OMI directory:
- an array with AMFs, some re-calculated
- an array with value 1 at indices where the AMF and VCD have been 
  re-calculated
- an array with dimensions of the WRF file with search distances in 
  degrees

"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import datetime as dt
import time

#execfile('find_WRF_pixel.py')
execfile('haversine.py')
execfile('OMI_recalc_AMF.py')

def select_time(OMIpath):
    t = dt.datetime(2014,1,1,int(OMIpath[-24:-22]),int(OMIpath[-22:-20]))
    bs = dt.datetime(2014,1,1)
    arr = np.array([bs + dt.timedelta(hours=i) for i in np.arange(0,24)])
    t_sel = min(arr, key=lambda x: abs(x - t))
    return t_sel.hour


#def recalc_AMF_loop(WRFpath,WRFfile,OMIpath):
def recalc_AMF_loop(WRFpath,OMIpath):
    
    #wds = nc.Dataset(WRFpath + WRFfile,'r')
    #wrf_lat = wds.variables['XLAT'][0,:,:]
    #wrf_lon = wds.variables['XLONG'][0,:,:]
    
    #Extract OMI files (on day of interest)
    #OMIfiles = [os.path.join(OMIpath,filename) for filename in os.listdir(OMIpath) if                 filename.startswith('QA4ECV_L2_NO2_OMI_' + WRFfile[11:15]+WRFfile[16:18]+WRFfile[19:21])]
    OMIfiles = [os.path.join(OMIpath,filename) for filename in os.listdir(OMIpath) if                 filename.startswith('QA4ECV_L2_NO2_OMI_')]
    OMIfiles.sort()
    
    #Loop over all OMI files
    for file in OMIfiles:
        tstart = time.time()
        
        ods    = nc.Dataset(file,'r+',format='NETCDF4')
        prod   = ods.groups['PRODUCT']
        detres = ods.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS']
        inp    = ods.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['INPUT_DATA']
        geo    = ods.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['GEOLOCATIONS']
        
        OMI_lat  = prod.variables['latitude'][0,:,:]
        OMI_lon  = prod.variables['longitude'][0,:,:]
        err_flag = prod.variables['processing_error_flag'][0,:,:]
        amf_trop = prod.variables['amf_trop'][0,:,:]
        amf_geo  = detres.variables['amf_geo'][0,:,:]
        siflag   = inp.variables['snow_ice_flag'][0,:,:]
        crf_no2  = detres.variables['cloud_radiance_fraction_no2'][0,:,:]
        solar_za = geo.variables['solar_zenith_angle'][0,:,:]
        cstripe = detres.variables['scd_no2_stripe_correction_omi'][0,:]
        
        #Implement conditions from Product Specification Document by 
        #Boersma et al. (2017) here, as well as the lat/lon range
        inds = np.where((OMI_lat > -5.) & (OMI_lat < 15.) & (OMI_lon > -83.) & (OMI_lon < -60.) & (err_flag == 0) & \
                        ((siflag < 10) | (siflag == 255)) & (amf_trop / amf_geo > 0.2) & (crf_no2 <= 0.5) & (solar_za < 80.))
        print 'Orbit: ' + file + ', valid pixels: ' + str(inds[0].shape[0])
        
        ###Create modified variables###
        
        #Create variables, add units and long names
        amf_trop_mod = prod.createVariable('AMF_trop_mod_AV_vz','f4',('time','scanline','ground_pixel'))
        amf_trop_mod.units = prod.variables['amf_trop'].units
        amf_trop_mod.long_name = u'modified tropospheric air mass factor (Auke Visser, WUR)'
        
        vcd_trop_mod = prod.createVariable('VCD_trop_mod_AV_v2','f4',('time','scanline','ground_pixel'))
        vcd_trop_mod.units = prod.variables['tropospheric_no2_vertical_column'].units
        vcd_trop_mod.long_name = u'modified tropospheric vertical column of nitrogen dioxide (Auke Visser, WUR)'
        
        amf_rec_flag = prod.createVariable('amf_rec_flag_AV_v2','i1',('time','scanline','ground_pixel'))
        amf_rec_flag.units = u'1'
        amf_rec_flag.long_name = u'flag indicating where AMF is recalculated'
        
        #Fill new variables with data
        amf_trop_mod[:,:,:] = ods.groups['PRODUCT'].variables['amf_trop']
        vcd_trop_mod[:,:,:] = ods.groups['PRODUCT'].variables['tropospheric_no2_vertical_column']
        amf_rec_flag[:,:,:] = np.zeros((1,amf_trop_mod.shape[1],amf_trop_mod.shape[2]))
        
        #Open corresponding WRF-Chem file
        m = file[-29:-27]; d = file[-27:-25]
        #print WRFpath + 'wrfout_d01_2014-' + m + '-' + d + '_00:00:00'
        #wds = nc.Dataset(WRFpath + 'wrfout_d01_2014-' + m + '-' + d + '_00:00:00','r')
	wds = nc.Dataset(WRFpath + 'wrfout_d01_2014-01-01_00:00:00','r')
        wrf_lat = wds.variables['XLAT'][0,:,:]
        wrf_lon = wds.variables['XLONG'][0,:,:]
        
        ###Recalculate AMFs/VCDs for the valid pixels###
        
        #Find time for comparison
        comptime = select_time(file)
        
        #loop over indices
        #loop over indices in lat/lon range
        for i in np.arange(inds[0].shape[0]):
            if np.mod(i,500) == 0:
                print 'Processed ' + str(i) + ' OMI pixels'
                
            #ilat = OMI_inds[0][i]; ilon = OMI_inds[1][i]
            ilat = inds[0][i]; ilon = inds[1][i]
            
            lon_pix = OMI_lon[ilat,ilon]
            lat_pix = OMI_lat[ilat,ilon]

            xind,yind,d,lon,lat = haversine_array(lon_pix,lat_pix,wrf_lon,wrf_lat)
            
            #Only update AMF and VCD for particular scanlines and if
            #the pixel is sufficiently close to a WRF-Chem cell.
            if ( ilon >= 5 and ilon <= 25 and d < 20. ):
                #update AMF recalculation flag
                amf_rec_flag[0,ilat,ilon] += 1
                
                amf_trop_mod[0,ilat,ilon],vcd_trop_mod[0,ilat,ilon] = recalc_AMF(ods,wds,ilat,ilon,lat,lon,comptime,51) #>50 is mass balance check, >100 is figures
                
                #Apply stripe correction to adjusted and re-calculated pixel
                vcd_trop_mod[0,ilat,ilon] -= cstripe[ilon] * 1.e15
            #If the pixel is outside the scanline range but within the WRF-Chem
            #domain, raise an error to ensure these pixels do not end up in the
            #monthly average produced by the gridding routine.
            elif (d < 20. ):
                err_flag[ilat,ilon] += 1
                
        ods.close()
        wds.close()
        
        tend = time.time()
        print 'Elapsed time for orbit: ' + str(tend - tstart) + ' s'
    
    """
    #Generate list of files with relevant lat/lon data
    for file in OMIfiles:
        ods    = nc.Dataset(file,'r+')
        prod   = ods.groups['PRODUCT']
        detres = ods.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS']
        inp    = ods.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['INPUT_DATA']
        geo    = ods.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['GEOLOCATIONS']
        
        OMI_lat  = prod.variables['latitude'][0,:,:]
        OMI_lon  = prod.variables['longitude'][0,:,:]
        err_flag = prod.variables['processing_error_flag'][0,:,:]
        amf_trop = prod.variables['amf_trop'][0,:,:]
        amf_geo  = detres.variables['amf_geo'][0,:,:]
        siflag   = inp.variables['snow_ice_flag'][0,:,:]
        crf_no2  = detres.variables['cloud_radiance_fraction_no2'][0,:,:]
        solar_za = geo.variables['solar_zenith_angle'][0,:,:]
        
        #Implement conditions from Product Specification Document by 
        #Boersma et al. (2017) here, as well as the lat/lon range
        #inds = np.where((OMI_lat[:,5:25] > 40.) & (OMI_lat[:,5:25] < 60.) & (OMI_lon[:,5:25] > -10.) & (OMI_lon[:,5:25] < 25.))
        inds = np.where((OMI_lat > 35.) & (OMI_lat < 65.) & (OMI_lon > -10.) & (OMI_lon < 35.) & (err_flag == 0) & \
                        ((siflag < 10) | (siflag == 255)) & (amf_trop / amf_geo > 0.2) & (crf_no2 <= 0.5) & (solar_za < 80.))
        
        if inds[0].shape[0] > 0:
            OMIfiles_EUR.append(file)
            OMI_inds.append(inds)
            
    #Current solution to omit backward orbits
    print OMIfiles_EUR
    #OMIfiles_EUR = [OMIfiles_EUR[1]]
    #OMI_inds = [OMI_inds[1]]
    
    #Initialize output arrays
    #AMF_trop     = np.zeros((len(OMIfiles_EUR),1644,60))
    #AMF_trop_adj = np.zeros((len(OMIfiles_EUR),1644,60))
    #VCD_trop     = np.zeros((len(OMIfiles_EUR),1644,60))
    #VCD_trop_adj = np.zeros((len(OMIfiles_EUR),1644,60))
    #VCD_trop_rec = np.zeros((len(OMIfiles_EUR),1644,60))
    #distance     = np.zeros((len(OMIfiles_EUR),1644,60))
    #clats        = np.zeros((len(OMIfiles_EUR),1644,60,4))
    #clons        = np.zeros((len(OMIfiles_EUR),1644,60,4))
    AMF_trop     = []
    AMF_trop_adj = []
    VCD_trop     = []
    VCD_trop_adj = []
    VCD_trop_rec = []
    distance     = []
    clats        = []
    clons        = []
    cstripe      = np.zeros((len(OMIfiles_EUR),60))
    
    #Loop over OMI swaths
    for k in range(len(OMIfiles_EUR)):
        
        print OMI_inds[k][0].shape
        print 'Currently processing: ' + OMIfiles_EUR[k]
        ods     = nc.Dataset(OMIfiles_EUR[k],'r')
        prod    = ods.groups['PRODUCT']
        geo     = ods.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['GEOLOCATIONS']
        detres  = ods.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS']
        OMI_lat = prod.variables['latitude'][0,:,:]
        OMI_lon = prod.variables['longitude'][0,:,:]
        cstripe[k,:] = detres.variables['scd_no2_stripe_correction_omi'][0,:]
        
        AMF_trop.append(np.zeros((OMI_lat.shape[0],60)))    
        AMF_trop_adj.append(np.zeros((OMI_lat.shape[0],60)))
        VCD_trop.append(np.zeros((OMI_lat.shape[0],60)))    
        VCD_trop_adj.append(np.zeros((OMI_lat.shape[0],60)))
        VCD_trop_rec.append(np.zeros((OMI_lat.shape[0],60)))
        distance.append(np.zeros((OMI_lat.shape[0],60)))    
        clats.append(np.zeros((OMI_lat.shape[0],60)))       
        clons.append(np.zeros((OMI_lat.shape[0],60)))       
        
        #Extract variables for selection criteria
        #clats[k,:,:,:] = geo.variables['latitude_bounds'][0,:,:,:]
        #clons[k,:,:,:] = geo.variables['longitude_bounds'][0,:,:,:]
        clats[k] = geo.variables['latitude_bounds'][0,:,:,:]
        clons[k] = geo.variables['longitude_bounds'][0,:,:,:]
        
        comptime = select_time(OMIfiles_EUR[k])
        print comptime
        
        #loop over indices in lat/lon range
        for i in np.arange(OMI_inds[k][0].shape[0]):
            if np.mod(i,100) == 0:
                print 'Processed ' + str(i) + ' OMI pixels'
                
            ilat = OMI_inds[k][0][i]; ilon = OMI_inds[k][1][i]
            
            
            if ( ilon >= 5 and ilon <= 25 ):
                lon_pix = OMI_lon[ilat,ilon]
                lat_pix = OMI_lat[ilat,ilon]

                xind,yind,d,lon,lat = haversine_array(lon_pix,lat_pix,wrf_lon,wrf_lat)
                distance[k][ilat,ilon] = d
                
                AMF_trop[k][ilat,ilon], AMF_trop_adj[k][ilat,ilon],VCD_trop[k][ilat,ilon],VCD_trop_adj[k][ilat,ilon],VCD_trop_rec[k][ilat,ilon] = recalc_AMF(ods,wds,ilat,ilon,lat,lon,comptime,0)
                
                #Apply stripe correction to adjusted and re-calculated pixel
                VCD_trop_adj[k][ilat,ilon] -= cstripe[k,ilon] * 1.e15
                VCD_trop_rec[k][ilat,ilon] -= cstripe[k,ilon] * 1.e15

        ods.close()
        
    wds.close()
    """        
    #return AMF_trop, AMF_trop_adj, VCD_trop, VCD_trop_adj, VCD_trop_rec, distance, clats, clons
    return 
