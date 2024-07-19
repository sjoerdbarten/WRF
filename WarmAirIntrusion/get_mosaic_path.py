import numpy as np
import netCDF4 as nc
import datetime as dt
import pandas as pd
from math import cos, asin, sqrt
from itertools import repeat
from tempfile import TemporaryFile
import os

print('This script takes a WRFout file and calculates from a preprocessed dataframe the grid cells corresponding to the MOSAiC drift')

#------------------------------------------------------------------------------------------------------------- 
#The following functions are defined and called later in the script
#------------------------------------------------------------------------------------------------------------- 

#This is a function to find the distance between two given lat/lon combinations
def distance(lat1, lon1, lat2, lon2):
    p = 0.017453292519943295
    a = 0.5 - cos((lat2-lat1)*p)/2 + cos(lat1*p)*cos(lat2*p) * (1-cos((lon2-lon1)*p)) / 2
    return 12742 * asin(sqrt(a))

#This is a function that returns the closest lat/lon combination
def closest(data, v):
    return min(data, key=lambda p: distance(v['lat'],v['lon'],p['lat'],p['lon']))

#This function puts the WRF and MOSAiC lat/lon data in the right format (list of dictionaries) to ust in the functions above
def getclosest(latmosaic,lonmosaic,wrflats,wrflons):
    dlist=[]                                         #Create empty list
    for i in range(len(wrflats)):                    #Loop over all latitudes
        for j in range(len(wrflons)):                #Loop over all longitudes
            d={}                                     #Create an empty dictionary
            d['lat']=wrflats[i,j]                      #Fill the dictionary with WRF latitudes
            d['lon']=wrflons[i,j]                      #Fill the dictionary with WRF longitudes
            dlist.append(d)                          #Fill the list with the dictionary
        
    v = {'lat': latmosaic, 'lon': lonmosaic}         #Create a dictionary of the MOSAiC latitude/longitudes
    #print(v)
    closestlatlon = closest(dlist,v)                 #Use the closest (and distance) functions to find the closest lat/lon combination
    latindex = np.where(wrflats == closestlatlon['lat'])  #Now we find the index of the CAMS latarray that corresponds to that pair
    lonindex = np.where(wrflons == closestlatlon['lon'])  #Now we find the index of the CAMS lonarray that corresponds to that pair
    return(latindex,lonindex)                       #We return the indexes because we need them for later

def make_locs_array(filename,appendixname):
	#wrffile = '/archive/ESG/barte035/MOSAiC/WarmAirIntrusion/input_61levels/wrfout/wrfout_warmairintrusion_chem_d03_2020-04-05_00:00:00'
	wrffile = filename
	ds = nc.Dataset(wrffile,'r')

	wrflats = ds.variables['XLAT'][0,:,:]
	wrflons = ds.variables['XLONG'][0,:,:]
	timeswrf = ds.variables['XTIME'][:]
	print(timeswrf.shape)
	timeswrfdatetime = np.arange(dt.datetime(2020,4,5,0,0,0), dt.datetime(2020,4,22,1,0,0,0), dt.timedelta(hours=1)).astype(dt.datetime)
	print(timeswrfdatetime.shape)
	
	#Here we read in the MOSAiC data: we need the times, latituded and longitudes
	mosaiclatlondata = pd.read_csv('df_lat_lon_for_WRF_postprocessing.csv')		#Read the .csv file to a pd dataframe with lat lon and datetime
	#mosaiclatlondata['datetime'] = pd.to_datetime(mosaiclatlondata['datetime'])
	#mosaictimelist = mosaiclatlondata.datetime                                 		#Get the datetime timestamps
	mosaictimelist = [dt.datetime.strptime(date, "%Y-%m-%d %H:%M:%S") for date in mosaiclatlondata.datetime.values]
	mosaiclatlondata['datetime'] = mosaictimelist
	mosaiclats = mosaiclatlondata['lat']                             			#Get the latutides
	mosaiclons = mosaiclatlondata['lon']				                        #Get the longitudes
	
	
	
	locs = []
	
	for i in timeswrfdatetime:
		if i in mosaictimelist:
			mosaiclat = mosaiclatlondata.loc[mosaiclatlondata.datetime == i, 'lat'].values[0]
			mosaiclon = mosaiclatlondata.loc[mosaiclatlondata.datetime == i, 'lon'].values[0]
		
			print('Timestamp = ',i)
			print('MOSAiC latitude = ',mosaiclat)
			print('MOSAiC longitude = ',mosaiclon)
	
			latindex,lonindex = getclosest(mosaiclat,mosaiclon,wrflats,wrflons) 	      #Get the latindex and lonindex using the getclosest function
					
			
			print('WRF latitude = ',wrflats[latindex][0])
			print('WRF longitude = ',wrflons[lonindex][0])
			
			#locs.extend(repeat((lonindex[0][0],lonindex[1][0]),6))
			locs.extend((lonindex[0][0],lonindex[1][0]))
			
			#First one: south_north (j in ncview)
			#Second one: west_east (i in ncview)
			
			
	#print(locs)
	#print(locs.shape)
	#np.save('/home/WUR/barte035/WRFChem/Python/WarmAirIntrusion/locs.npy',locs)
	np.save('/home/WUR/barte035/WRFChem/Python/WarmAirIntrusion/locs'+appendixname+'.npy',locs)

'''
        float T2(Time, south_north, west_east) ;
                T2:FieldType = 104 ;
                T2:MemoryOrder = "XY " ;
                T2:description = "TEMP at 2 M" ;
                T2:units = "K" ;
                T2:stagger = "" ;
                T2:coordinates = "XLONG XLAT XTIME" ;
'''

def make_netcdf_file(filename,appendixname):		
	wrffile = filename
	#wrffile = '/archive/ESG/barte035/MOSAiC/WarmAirIntrusion/input_61levels/wrfout/wrfout_warmairintrusion_chem_d03_2020-04-05_00:00:00'
	ds = nc.Dataset(wrffile,'r')
	timeswrf = ds.variables['XTIME'][:]
	
	locs2 = np.load('/home/WUR/barte035/WRFChem/Python/WarmAirIntrusion/locs'+appendixname+'.npy')
	#locs2 = np.load('/home/WUR/barte035/WRFChem/Python/WarmAirIntrusion/locs.npy')
	
	print(locs2)
	print(locs2.shape)
	print(timeswrf.shape)
	
	netcdf_names = list(ds.variables.keys())
	netcdf_names.remove('Times')
	netcdf_names.remove('XTIME')
	netcdf_names.remove('P_TOP')
	
#	netcdf_names = ['o3','T2','SST','XLAT','XLONG']

	fn = '/home/WUR/barte035/WRFChem/Python/WarmAirIntrusion/track'+appendixname+'.nc'
	#os.remove(fn)
	dsfn = nc.Dataset(fn, 'w', format='NETCDF4')
	time = dsfn.createDimension('time',None)
	bio_emissions_dimension_stag = dsfn.createDimension('bio_emissions_dimension_stag',159)
	bottom_top = dsfn.createDimension('bottom_top',60)
	bottom_top_stag = dsfn.createDimension('bottom_top_stag',61)
	soil_layers_stag = dsfn.createDimension('soil_layers_stag',4)
	
	#Here we extract from the WRFoutput MOSAiC location and save it as a different netcdf file.
	for varname in netcdf_names:	
		print(ds.variables[varname].shape)		
		if ds.variables[varname].ndim == 3:
			value = dsfn.createVariable(varname, 'f4', ('time',))
			for i in range(len(timeswrf)):
				value[i] = ds.variables[varname][i,locs2[i*2],locs2[i*2+1]]				
		if ds.variables[varname].ndim == 4 and ds.variables[varname].shape[1] == 159:
			value = dsfn.createVariable(varname, 'f4', ('time','bio_emissions_dimension_stag',))
			for i in range(len(timeswrf)):
				value[i,:] = ds.variables[varname][i,:,locs2[i*2],locs2[i*2+1]]
		if ds.variables[varname].ndim == 4 and ds.variables[varname].shape[1] == 60:
			value = dsfn.createVariable(varname, 'f4', ('time','bottom_top',))
			for i in range(len(timeswrf)):
				value[i,:] = ds.variables[varname][i,:,locs2[i*2],locs2[i*2+1]]
		if ds.variables[varname].ndim == 4 and ds.variables[varname].shape[1] == 61:
			value = dsfn.createVariable(varname, 'f4', ('time','bottom_top_stag',))
			for i in range(len(timeswrf)):
				value[i,:] = ds.variables[varname][i,:,locs2[i*2],locs2[i*2+1]]
		if ds.variables[varname].ndim == 4 and ds.variables[varname].shape[1] == 4:
			value = dsfn.createVariable(varname, 'f4', ('time','soil_layers_stag',))
			for i in range(len(timeswrf)):
				value[i,:] = ds.variables[varname][i,:,locs2[i*2],locs2[i*2+1]]
				

filename = '/archive/ESG/barte035/MOSAiC/WarmAirIntrusion/input_61levels/wrfout/wrfout_warmairintrusion_chem_v9_d03_2020-04-05_00:00:00'
appendixname = '_v9_d03'
make_locs_array(filename,appendixname)	
make_netcdf_file(filename,appendixname)
