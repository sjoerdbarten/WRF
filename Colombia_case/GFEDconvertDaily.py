"""
emissions unit from GFED --> kgDM m-2 month-1
emissions factor unit --> gr<species> kg-1DM
to convert: DM emissions * emissions factor
	unit: gr<species> m-2 month-1
unit for WRFChem: microgram <species> m-2 s-1 (for aerosols) mole <species> km-2 h-1 (for gas)
	
	1 gr(CO) m-2 month-1 = CO (gr)* mrCO (mol/gr) /(1 m2 / 1000000 km2/m2) /(1month* 31 days/month)
	
	aerosol --> PM25 gr m-2 month-1 = 1000000 (ug/gr) /(1 (m2/m2) * 31*24*60*60 (second/month))
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import netCDF4 as nc
from netCDF4 import Dataset
from mpl_toolkits import basemap
import os

GFEDdatamap="/home/WUR/barte035/WRFChem/GFED/fire_emissions_v4_R1_1293/data"
FINNdatamap="/home/WUR/barte035/WRFChem/FINN/src"
WRFChemmap="/home/WUR/barte035/WRFChem/WPS"

sources = 'SAVA','BORF','TEMF','DEFO','PEAT','AGRI'
emiss_factor = os.path.join(GFEDdatamap,'GFED4_Emission_Factors.txt')
#geodata = os.path.join(WRFChemmap,'geo_em.d01.nc')
#geodata = os.path.join(WRFChemmap,'geo_em.d02.nc')
GFED_data = os.path.join(GFEDdatamap,'GFED4.1s_2014.hdf5')

"""
#Read in emission factors
"""
EFs = np.zeros((41, 6)) # 41 species, 6 sources
k = 0
f = open(emiss_factor)
while 1:
	line = f.readline()
	if line == "":
		break
	
	if line[0] != '#':
		contents = line.split()
		#species.append(contents[0])
		EFs[k,:] = contents[1:]
		k += 1
f.close()



EF_CO = EFs[2,:]
EF_NOx = EFs[7,:]
EF_PM25 = EFs[9,:]
EF_OC = EFs[12,:]
EF_BC = EFs[13,:]
EF_SO2 = EFs[14,:]
EF_C2H6 = EFs[15,:]
EF_CH3OH = EFs[16,:]
EF_C2H5OH = EFs[17,:]
EF_C3H8 = EFs[18,:]
EF_C2H4 = EFs[20,:]
EF_C3H6 = EFs[21,:]
EF_C10H16 = EFs[23,:]
EF_toluene = EFs[27,:]
EF_bigene = EFs[28,:]
EF_bigalk = EFs[29,:]
EF_CH2O = EFs[30,:]
EF_NH3 = EFs[33,:]
EF_CH3COOH = EFs[37,:]
EF_MEK = EFs[38,:]

"""
#read in GFED data
"""

g = h5py.File(GFED_data, 'r')
xlat = g['lat'][:]
xlon = g['lon'][:]

##CO emission global
CO_emissions = np.zeros((720, 1440))
NOx_emissions = np.zeros((720, 1440))
PM25_emissions = np.zeros((720, 1440))
OC_emissions = np.zeros((720, 1440))
BC_emissions = np.zeros((720, 1440))
SO2_emissions = np.zeros((720, 1440))
C2H6_emissions = np.zeros((720, 1440))
CH3OH_emissions = np.zeros((720, 1440))
C2H5OH_emissions = np.zeros((720, 1440))
C3H8_emissions = np.zeros((720, 1440))
C2H4_emissions = np.zeros((720, 1440))
C3H6_emissions = np.zeros((720, 1440))
C10H16_emissions = np.zeros((720, 1440))
toluene_emissions = np.zeros((720, 1440))
bigene_emissions = np.zeros((720, 1440))
bigalk_emissions = np.zeros((720, 1440))
CH2O_emissions = np.zeros((720, 1440))
NH3_emissions = np.zeros((720, 1440))
CH3COOH_emissions = np.zeros((720, 1440))
MEK_emissions = np.zeros((720, 1440))

string = '/emissions/01/DM' # read in DM emissions for January only
DM_emissions = g[string][:]

"""
regrid
"""
#global array
xlat_globe = np.flipud(xlat[:,0]) #because the lat is not exceeding, it needs to be flipped
xlon_globe = xlon[0,:]
##define colombia array

"""
write a netCDF file
"""

loop_day = 1
loop_day1 = 0
loop_day2 = 0

loop_3h_num = 1

domain = 1

for loop_day in range(1,32):
		if domain == 1:
			geodata = os.path.join(WRFChemmap,'geo_em.d01.nc')
			mf = nc.Dataset(geodata)
			lat_col = mf.variables['XLAT_M'][0,:,:]
			lon_col = mf.variables['XLONG_M'][0,:]
		if domain == 2:
			geodata = os.path.join(WRFChemmap,'geo_em.d02.nc')
			mf = nc.Dataset(geodata)
			lat_col = mf.variables['XLAT_M'][0,:,:]
			lon_col = mf.variables['XLONG_M'][0,:]
		if loop_day == (0,9,1):
			loop_day1 = 0
		elif loop_day == 10:
			loop_day1 = 1
		elif loop_day == 20:
			loop_day1 = 2
		elif loop_day == 30:
			loop_day1 = 3
		if loop_day == 0 or loop_day == 10 or loop_day == 20 or loop_day == 30:
			loop_day2 = 0	
		if loop_day == 1 or loop_day == 11 or loop_day == 21 or loop_day == 31:
			loop_day2 = 1
		if loop_day == 2 or loop_day == 12 or loop_day == 22:
			loop_day2 = 2
		if loop_day == 3 or loop_day == 13 or loop_day == 23:
			loop_day2 = 3
		if loop_day == 4 or loop_day == 14 or loop_day == 24:
			loop_day2 = 4
		if loop_day == 5 or loop_day == 15 or loop_day == 25:
			loop_day2 = 5		
		if loop_day == 6 or loop_day == 16 or loop_day == 26:
			loop_day2 = 6		
		if loop_day == 7 or loop_day == 17 or loop_day == 27:
			loop_day2 = 7
		if loop_day == 8 or loop_day == 18 or loop_day == 28:
			loop_day2 = 8		
		if loop_day == 9 or loop_day == 19 or loop_day == 29:
			loop_day2 = 9		
	
		# read in the day fraction
		d = '/emissions/01/daily_fraction/day_'+str(loop_day)
		if loop_day == 10:
			d = '/emissions/01/daily_fraction/day_'+str(loop_day1)+str(loop_day2)
		day_frac = g[d][:]
	
		print 'domain='+str(domain),'loopday='+str(loop_day)
		
		for source in range(6):
			# read in the fractional contribution of each source
			string = '/emissions/01/partitioning/DM_'+sources[source]
			contribution = g[string][:]
			specs = ['CO','NOx','PM25','OC','BC','SO2','C2H6','CH3OH','C2H5OH','C3H8','C2H4','C3H6','C10H16','toluene','bigene','bigalk','CH2O','NH3','CH3COOH','MEK']
			for k in range(len(specs)):
				if specs[k]=='CO':
					CO_emissions += DM_emissions* day_frac*contribution*EF_CO[source]/28.01*10e6/(31*24) #mole CO per km2 per hours
				elif specs[k]=='NOx':
					NOx_emissions += DM_emissions* day_frac*contribution*EF_NOx[source]/30.01*10e6/(31*24) #NOx is represented as NO (mole NO per km2 per hours)
				elif specs[k]=='PM25':
					PM25_emissions += DM_emissions* day_frac*contribution*EF_PM25[source]*1000000/(31*24*60*60) #microgr PM25 per m2 per second
				elif specs[k]=='OC':
					OC_emissions += DM_emissions* day_frac*contribution*EF_OC[source]*1000000/(31*24*60*60) #microgr OC per m2 per second
				elif specs[k]=='BC':
					BC_emissions += DM_emissions* day_frac*contribution*EF_BC[source]*1000000/(31*24*60*60) #microgr BC per m2 per second
				elif specs[k]=='SO2':
					SO2_emissions += DM_emissions* day_frac*contribution*EF_SO2[source]/64.02*10e6/(31*24) #mole SO2 per km2 per hours
				elif specs[k]=='C2H6':
					C2H6_emissions += DM_emissions* day_frac*contribution*EF_C2H6[source]/30.07*10e6/(31*24) #mole C2H6 per km2 per hours
				elif specs[k]=='CH3OH':
					CH3OH_emissions += DM_emissions* day_frac*contribution*EF_CH3OH[source]/32.04*10e6/(31*24) #mole CH3OH per km2 per hours
				elif specs[k]=='C2H5OH':
					C2H5OH_emissions += DM_emissions* day_frac*contribution*EF_C2H5OH[source]*46.07/10e6/(31*24) #mole C2H5OH per km2 per hours
				elif specs[k]=='C3H8':
					C3H8_emissions += DM_emissions* day_frac*contribution*EF_C3H8[source]/44.1*10e6/(31*24) #mole C3H8 per km2 per hours
				elif specs[k]=='C2H4':
					C2H4_emissions += DM_emissions* day_frac*contribution*EF_C2H4[source]/28.05*10e6/(31*24) #mole C2H4 per km2 per hours
				elif specs[k]=='C3H6':
					C3H6_emissions += DM_emissions* day_frac*contribution*EF_C3H6[source]/42.08*10e6/(31*24) #mole C3H6 per km2 per hours
				elif specs[k]=='C10H16':
					C10H16_emissions += DM_emissions* day_frac*contribution*EF_C10H16[source]/136.24*10e6/(31*24) #mole C10H16 per km2 per hours
				elif specs[k]=='toluene':
					toluene_emissions += DM_emissions* day_frac*contribution*EF_toluene[source]/92.14*10e6/(31*24) #toluene is C7H8(mole per km2 per hours)
				elif specs[k]=='bigene':
					bigene_emissions += DM_emissions* day_frac*contribution*EF_bigene[source]/56*10e6/(31*24) #mole butene per km2 per hours
				elif specs[k]=='bigalk':
					bigalk_emissions += DM_emissions* day_frac*contribution*EF_bigalk[source]/58*10e6/(31*24) #mole butane per km2 per hours
				elif specs[k]=='CH2O':
					CH2O_emissions += DM_emissions* day_frac*contribution*EF_CH2O[source]/30.03*10e6/(31*24) #CH2O is formaldehyde
				elif specs[k]=='NH3':
					NH3_emissions += DM_emissions* day_frac*contribution*EF_NH3[source]/17.03*10e6/(31*24) #mole NH3 per km2 per hours
				elif specs[k]=='CH3COOH':
					CH3COOH_emissions += DM_emissions* day_frac*contribution*EF_CH3COOH[source]/60.05*10e6/(31*24) #mole CH3COOH per km2 per hours
				else:
					MEK_emissions += DM_emissions* day_frac*contribution*EF_MEK[source]/72.11*10e6/(31*24) #mole MEK per km2 per hours
		
		#print NOx_emissions.shape
		#print np.max(NOx_emissions)
		
		if domain == 1:
			CO_emissions_colombia = basemap.interp(np.flipud(CO_emissions), xlon_globe, xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			NOx_emissions_colombia = basemap.interp(np.flipud(NOx_emissions), xlon_globe, xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			PM25_emissions_colombia = basemap.interp(np.flipud(PM25_emissions), xlon_globe, xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			OC_emissions_colombia = basemap.interp(np.flipud(OC_emissions), xlon_globe, xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			BC_emissions_colombia = basemap.interp(np.flipud(BC_emissions), xlon_globe, xlat_globe,lon_col, lat_col, order=1).reshape((1,1,99,99))
			SO2_emissions_colombia = basemap.interp(np.flipud(SO2_emissions), xlon_globe, xlat_globe,lon_col, lat_col, order=1).reshape((1,1,99,99))
			C2H6_emissions_colombia = basemap.interp(np.flipud(C2H6_emissions), xlon_globe,xlat_globe, lon_col, lat_col,order=1).reshape((1,1,99,99))
			CH3OH_emissions_colombia = basemap.interp(np.flipud(CH3OH_emissions), xlon_globe,xlat_globe, lon_col, lat_col,order=1).reshape((1,1,99,99))
			C2H5OH_emissions_colombia = basemap.interp(np.flipud(C2H5OH_emissions),xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			C3H8_emissions_colombia = basemap.interp(np.flipud(C3H8_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			C2H4_emissions_colombia = basemap.interp(np.flipud(C2H4_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			C3H6_emissions_colombia = basemap.interp(np.flipud(C3H6_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			C10H16_emissions_colombia = basemap.interp(np.flipud(C10H16_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			toluene_emissions_colombia = basemap.interp(np.flipud(toluene_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			bigene_emissions_colombia = basemap.interp(np.flipud(bigene_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			bigalk_emissions_colombia = basemap.interp(np.flipud(bigalk_emissions), xlon_globe,xlat_globe,lon_col,lat_col,order=1).reshape((1,1,99,99))
			CH2O_emissions_colombia = basemap.interp(np.flipud(CH2O_emissions), xlon_globe,xlat_globe, lon_col, lat_col,order=1).reshape((1,1,99,99))
			NH3_emissions_colombia = basemap.interp(np.flipud(NH3_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			CH3COOH_emissions_colombia = basemap.interp(np.flipud(CH3COOH_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,99,99))
			MEK_emissions_colombia = basemap.interp(np.flipud(MEK_emissions), xlon_globe, xlat_globe,lon_col, lat_col, order=1).reshape((1,1,99,99))
	
		if domain == 2:
			CO_emissions_colombia = basemap.interp(np.flipud(CO_emissions), xlon_globe, xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			NOx_emissions_colombia = basemap.interp(np.flipud(NOx_emissions), xlon_globe, xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			PM25_emissions_colombia = basemap.interp(np.flipud(PM25_emissions), xlon_globe, xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			OC_emissions_colombia = basemap.interp(np.flipud(OC_emissions), xlon_globe, xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			BC_emissions_colombia = basemap.interp(np.flipud(BC_emissions), xlon_globe, xlat_globe,lon_col, lat_col, order=1).reshape((1,1,183,183))
			SO2_emissions_colombia = basemap.interp(np.flipud(SO2_emissions), xlon_globe, xlat_globe,lon_col, lat_col, order=1).reshape((1,1,183,183))
			C2H6_emissions_colombia = basemap.interp(np.flipud(C2H6_emissions), xlon_globe,xlat_globe, lon_col, lat_col,order=1).reshape((1,1,183,183))
			CH3OH_emissions_colombia = basemap.interp(np.flipud(CH3OH_emissions), xlon_globe,xlat_globe, lon_col, lat_col,order=1).reshape((1,1,183,183))
			C2H5OH_emissions_colombia = basemap.interp(np.flipud(C2H5OH_emissions),xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			C3H8_emissions_colombia = basemap.interp(np.flipud(C3H8_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			C2H4_emissions_colombia = basemap.interp(np.flipud(C2H4_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			C3H6_emissions_colombia = basemap.interp(np.flipud(C3H6_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			C10H16_emissions_colombia = basemap.interp(np.flipud(C10H16_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			toluene_emissions_colombia = basemap.interp(np.flipud(toluene_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			bigene_emissions_colombia = basemap.interp(np.flipud(bigene_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			bigalk_emissions_colombia = basemap.interp(np.flipud(bigalk_emissions), xlon_globe,xlat_globe,lon_col,lat_col,order=1).reshape((1,1,183,183))
			CH2O_emissions_colombia = basemap.interp(np.flipud(CH2O_emissions), xlon_globe,xlat_globe, lon_col, lat_col,order=1).reshape((1,1,183,183))
			NH3_emissions_colombia = basemap.interp(np.flipud(NH3_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			CH3COOH_emissions_colombia = basemap.interp(np.flipud(CH3COOH_emissions), xlon_globe,xlat_globe, lon_col, lat_col, order=1).reshape((1,1,183,183))
			MEK_emissions_colombia = basemap.interp(np.flipud(MEK_emissions), xlon_globe, xlat_globe,lon_col, lat_col, order=1).reshape((1,1,183,183))

		CO_emissions = 0
		NOx_emissions = 0
		PM5_emissions = 0
		OC_emissions = 0
		BC_emissions = 0
		SO2_emissions = 0
		C2H6_emissions = 0
		CH3OH_emissions = 0
		C2H5OH_emissions = 0
		C3H8_emissions = 0
		C2H4_emissions = 0
		C2H4_emissions = 0
		C3H6_emissions = 0
		C10H16_emissions = 0
		toluene_emissions = 0
		bigene_emissions = 0
		bigalk_emissions = 0
		CH2O_emissions = 0
		NH3_emissions = 0
		CH3COOH_emissions = 0
		MEK_emissions = 0
			
		FINN = Dataset(os.path.join(FINNdatamap,'wrffirechemi_d0'+str(domain)+'_2014-01-'+str(loop_day1)+str(loop_day2)+'_00:00:00'),'r')
		#we need FINN data to have fire size and faction of land type. These data will be used on plume rise module
		GFED = Dataset('wrffirechemi_d0'+str(domain)+'_2014-01-'+str(loop_day1)+str(loop_day2)+'_00:00:00','w',format='NETCDF4_CLASSIC')
		#WRF-Chem version 3.2 only allows netcdf version 3 file format with classic type.
	
		print 'wrffirechemi_d0'+str(domain)+'_2014-01-'+str(loop_day1)+str(loop_day2)+'_00:00:00 has been made'
		
		##global attributes
		GFED.Title = "GFED emissions for WRF"
		GFED.History = "Created on feb 2018"
		GFED.Author = "Sjoerd Barten"
		GFED.START_DATE = "2014-01-01_00:00:00"
		GFED.SIMULATION_START_DATE = "2014-01_00:00:00"
		if domain == 1:
			GFED.WEST_EAST_GRID_DIMENSION = int(99)
			GFED.SOUTH_NORTH_GRID_DIMENSION = int(99)
			GFED.DX = float(27000.0)
			GFED.DY = float(27000.0)
			GFED.WEST_EAST_PATCH_END_UNSTAG = int(99)
			GFED.WEST_EAST_PATCH_END_STAG = int(100)
			GFED.SOUTH_NORTH_PATCH_END_UNSTAG = int(99)
			GFED.SOUTH_NORTH_PATCH_END_STAG = int(100)
			GFED.GRID_ID = int(1)
			GFED.I_PARENT_START = int(1)
			GFED.J_PARENT_START = int(1)
			GFED.PARENT_GRID_RATIO = int(1)
			GFED.createDimension('west_east',99)
			GFED.createDimension('south_north',99)
		if domain == 2:
			GFED.WEST_EAST_GRID_DIMENSION = int(183)
			GFED.SOUTH_NORTH_GRID_DIMENSION = int(183)
			GFED.DX = float(9000.0)
			GFED.DY = float(9000.0)
			GFED.WEST_EAST_PATCH_END_UNSTAG = int(183)
			GFED.WEST_EAST_PATCH_END_STAG = int(184)
			GFED.SOUTH_NORTH_PATCH_END_UNSTAG = int(183)
			GFED.SOUTH_NORTH_PATCH_END_STAG = int(184)
			GFED.GRID_ID = int(2)
			GFED.I_PARENT_START = int(22)
			GFED.J_PARENT_START = int(22)
			GFED.PARENT_GRID_RATIO = int(3)
			GFED.createDimension('west_east',183)
			GFED.createDimension('south_north',183)
		GFED.BOTTOM_TOP_GRID_DIMENSION = int(61)
		GFED.GRIDTYPE = "C"
		GFED.DIFF_OPT = int(0)
		GFED.KM_OPT = int(4)
		GFED.DAMP_OPT = int(0)
		GFED.DAMPCOEF = float(0.0)
		GFED.KHDIF = float(0.0)
		GFED.KVDIF = float(0.0)
		GFED.MP_PHYSICS = int(10)
		GFED.RA_LW_PHYSICS = int(1)
		GFED.RA_SW_PHYSICS = int(1)
		GFED.SF_SFCLAY_PHYSICS = int(1)
		GFED.SF_SURFACE_PHYSICS = int(2)
		GFED.BL_PBL_PHYSICS = int(1)
		GFED.CU_PHYSICS = int(5)
		GFED.SURFACE_INPUT_SOURCE = int(3)
		GFED.SST_UPDATE = int(0)
		GFED.GRID_FDDA = int(0)
		GFED.GFDDA_INTERVAL_M = int(0)
		GFED.GFDDA_END_H = int(0)
		GFED.GRID_SFDDA = int(0)
		GFED.SGFDDA_INTERVAL_M = int(0)
		GFED.SGFDDA_END_H = int(0)
		GFED.WEST_EAST_PATCH_START_UNSTAG = int(1)
		GFED.WEST_EAST_PATCH_START_STAG = int(1)
		GFED.SOUTH_NORTH_PATCH_START_UNSTAG = int(1)
		GFED.SOUTH_NORTH_PATCH_START_STAG = int(1)
		GFED.BOTTOM_TOP_PATCH_START_UNSTAG = int(1)
		GFED.BOTTOM_TOP_PATCH_END_UNSTAG = int(60)
		GFED.BOTTOM_TOP_PATCH_START_STAG = int(1)
		GFED.BOTTOM_TOP_PATCH_END_STAG = int(61)
		GFED.PARENT_ID = int(1)
		GFED.DT = float(60.0)
		GFED.CEN_LAT = FINN.getncattr('CEN_LAT')
		GFED.CEN_LON = FINN.getncattr('CEN_LON')
		GFED.TRUELAT1 = FINN.getncattr('TRUELAT1')
		GFED.TRUELAT2 = FINN.getncattr('TRUELAT2')
		GFED.MOAD_CEN_LAT = FINN.getncattr('MOAD_CEN_LAT')
		GFED.STAND_LON = FINN.getncattr('STAND_LON')
		GFED.POLE_LAT = FINN.getncattr('POLE_LAT')
		GFED.POLE_LON = FINN.getncattr('POLE_LAT')
		GFED.GMT = FINN.getncattr('GMT')
		GFED.JULYR = int(2014)
		GFED.JULDAY = int(loop_day)
		GFED.MAP_PROJ = int(3)
		GFED.MMINLU = "USGS"
		GFED.NUM_LAND_CAT = int(24)
		GFED.ISWATER = int(16)
		GFED.ISLAKE = int(-1)
		GFED.ISICE = int(24)
		GFED.ISURBAN = int(1)
		GFED.ISOILWATER = int(14)
	
		GFED.createDimension('DateStrLen',19)
		GFED.createDimension('emissions_zdim_stag',1)
		GFED.createDimension('Time', 1)
		Times = GFED.createVariable('Times','S1',('Time','DateStrLen'))			

		Times.units = "secs since 1970-01-01 00:00:00"
		Times.long_name = "synthesized time coordinate from Times(time)"
		Times._CoordinateAxisType = "Times"
	
		Times [:] = ['2','0','1','4','-','0','1','-',''+str(loop_day1)+'',''+str(loop_day2)+'','_','0','0',':','0','0',':','0','0']
	
		co_emiss = GFED.createVariable('ebu_co', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		co_emiss.MemoryOrder = 'XYZ'
		co_emiss.description = 'Carbon Monoxide Emissions'
		co_emiss.units = 'mole km-2 hr-1'
		co_emiss.stagger = 'Z'
		co_emiss.FieldType = 104

		#NOx is defined as NO since NO and NO2 are rapidly interconverted in the atmosphere (Akagi et al, 2011)
		no_emiss = GFED.createVariable('ebu_no', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		no_emiss.MemoryOrder = 'XYZ'
		no_emiss.description = 'Nitrogen Monoxide Emissions'
		no_emiss.units = 'mole km-2 hr-1'
		no_emiss.stagger = 'Z'
		no_emiss.FieldType = 104

		so2_emiss = GFED.createVariable('ebu_so2', 'f4',('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		so2_emiss.MemoryOrder = 'XYZ'
		so2_emiss.description = 'Sulphur Dioxide Emissions'
		so2_emiss.units = 'mole km-2 hr-1'
		so2_emiss.stagger = 'Z'
		so2_emiss.FieldType = 104

		bigalk_emiss = GFED.createVariable('ebu_bigalk','f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		bigalk_emiss.MemoryOrder = 'XYZ'
		bigalk_emiss.description = 'Alkanes with >4 carbon Emissions'
		bigalk_emiss.units = 'mole km-2 hr-1'
		bigalk_emiss.stagger = 'Z'
		bigalk_emiss.FieldType = 104

		bigene_emiss = GFED.createVariable('ebu_oli', 'f4',('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		bigene_emiss.MemoryOrder = 'XYZ'
		bigene_emiss.description = 'Alkenes with >4 carbon Emissions'
		bigene_emiss.units = 'mole km-2 hr-1'
		bigene_emiss.stagger = 'Z'
		bigene_emiss.FieldType = 104	

		c2h4_emiss = GFED.createVariable('ebu_ol2', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		c2h4_emiss.MemoryOrder = 'XYZ'
		c2h4_emiss.description = 'Ethylen Emissions'
		c2h4_emiss.units = 'mole km-2 hr-1'
		c2h4_emiss.stagger = 'Z'
		c2h4_emiss.FieldType = 104

		c2h5oh_emiss = GFED.createVariable('ebu_eth', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		c2h5oh_emiss.MemoryOrder = 'XYZ'
		c2h5oh_emiss.description = 'Ethanol Emissions'
		c2h5oh_emiss.units = 'mole km-2 hr-1'
		c2h5oh_emiss.stagger = 'Z'
		c2h5oh_emiss.FieldType = 104

		#c2h6_emiss = GFED.createVariable('ebu_c2h6', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		#c2h6_emiss.MemoryOrder = 'XYZ'
		#c2h6_emiss.description = 'Ethane Emissions'
		#c2h6_emiss.units = 'mole km-2 hr-1'
		#c2h6_emiss.stagger = 'Z'
		#c2h6_emiss.FieldType = 104

		c3h8_emiss = GFED.createVariable('ebu_hc3', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		c3h8_emiss.MemoryOrder = 'XYZ'
		c3h8_emiss.description = 'Propane Emissions'
		c3h8_emiss.units = 'mole km-2 hr-1'
		c3h8_emiss.stagger = 'Z'
		c3h8_emiss.FieldType = 104

		c3h6_emiss = GFED.createVariable('ebu_olt', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		c3h6_emiss.MemoryOrder = 'XYZ'
		c3h6_emiss.description = 'Propene Emissions'
		c3h6_emiss.units = 'mole km-2 hr-1'
		c3h6_emiss.stagger = 'Z'
		c3h6_emiss.FieldType = 104

		ch2o_emiss = GFED.createVariable('ebu_ch2o', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		ch2o_emiss.MemoryOrder = 'XYZ'
		ch2o_emiss.description = 'Formaldehyde Emissions'
		ch2o_emiss.units = 'mole km-2 hr-1'
		ch2o_emiss.stagger = 'Z'
		ch2o_emiss.FieldType = 104

		#ch3oh_emiss = GFED.createVariable('ebu_in_ch3oh', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		#ch3oh_emiss.MemoryOrder = 'XYZ'
		#ch3oh_emiss.description = 'Methanol Emissions'
		#ch3oh_emiss.units = 'mole km-2 hr-1'
		#ch3oh_emiss.stagger = 'Z'
		#ch3oh_emiss.FieldType = 104
	
		mek_emiss = GFED.createVariable('ebu_mek', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		mek_emiss.MemoryOrder = 'XYZ'
		mek_emiss.description = 'Methyl Ethyl Keton Emissions'
		mek_emiss.units = 'mole km-2 hr-1'
		mek_emiss.stagger = 'Z'
		mek_emiss.FieldType = 104
	
		toluene_emiss = GFED.createVariable('ebu_tol', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		toluene_emiss.MemoryOrder = 'XYZ'
		toluene_emiss.description = 'Lumped Aromatics Emissions'
		toluene_emiss.units = 'mole km-2 hr-1'
		toluene_emiss.stagger = 'Z'
		toluene_emiss.FieldType = 104

		nh3_emiss = GFED.createVariable('ebu_nh3', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		nh3_emiss.MemoryOrder = 'XYZ'
		nh3_emiss.description = 'Ammonia Emissions'
		nh3_emiss.units = 'mole km-2 hr-1'
		nh3_emiss.stagger = 'Z'
		nh3_emiss.FieldType = 104
	
		#c10h16_emiss = GFED.createVariable('ebu_in_c10h16', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		#c10h16_emiss.MemoryOrder = 'XYZ'
		#c10h16_emiss.description = 'Terpenes Emissions'
		#c10h16_emiss.units = 'mole km-2 hr-1'
		#c10h16_emiss.stagger = 'Z'
		#c10h16_emiss.FieldType = 104
	
		ch3cooh_emiss = GFED.createVariable('ebu_ora2', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		ch3cooh_emiss.MemoryOrder = 'XYZ'
		ch3cooh_emiss.description = 'Acetic Acid Emissions'
		ch3cooh_emiss.units = 'mole km-2 hr-1'
		ch3cooh_emiss.stagger = 'Z'
		ch3cooh_emiss.FieldType = 104
	
		oc_emiss = GFED.createVariable('ebu_oc', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		oc_emiss.MemoryOrder = 'XYZ'
		oc_emiss.description = 'Organic Carbon Emissions'
		oc_emiss.units = 'ug m-2 s-1'
		oc_emiss.stagger = 'Z'
		oc_emiss.FieldType = 104
	
		pm25_emiss = GFED.createVariable('ebu_pm25', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		pm25_emiss.MemoryOrder = 'XYZ'
		pm25_emiss.description = 'PM25 Emissions'
		pm25_emiss.units = 'ug m-2 s-1'
		pm25_emiss.stagger = 'Z'
		pm25_emiss.FieldType = 104	

		bc_emiss = GFED.createVariable('ebu_bc', 'f4', ('Time','emissions_zdim_stag', 'south_north', 'west_east'))
		bc_emiss.MemoryOrder = 'XYZ'
		bc_emiss.description = 'Black Carbon Emissions'
		bc_emiss.units = 'ug m-2 s-1'
		bc_emiss.stagger = 'Z'
		bc_emiss.FieldType = 104
	
		FIRESIZE_AGEF = GFED.createVariable('FIRESIZE_AGEF', 'f4', ('Time','south_north', 'west_east'))
		FIRESIZE_AGEF.MemoryOrder = 'XY'
		FIRESIZE_AGEF.description = 'mean firesize of extra tropical forest from FINN'
		FIRESIZE_AGEF.units = 'm2'
		FIRESIZE_AGEF.stagger = ' '
		FIRESIZE_AGEF.FieldType = 104
	
		FIRESIZE_AGGR = GFED.createVariable('FIRESIZE_AGGR', 'f4', ('Time','south_north', 'west_east'))
		FIRESIZE_AGGR.MemoryOrder = 'XY'
		FIRESIZE_AGGR.description = 'mean firesize of grassland from FINN'
		FIRESIZE_AGGR.units = 'm2'
		FIRESIZE_AGGR.stagger = ' '
		FIRESIZE_AGGR.FieldType = 104
	
		FIRESIZE_AGSV = GFED.createVariable('FIRESIZE_AGSV', 'f4', ('Time','south_north', 'west_east'))
		FIRESIZE_AGSV.MemoryOrder = 'XY'
		FIRESIZE_AGSV.description = 'mean firesize of savana from FINN'
		FIRESIZE_AGSV.units = 'm2'
		FIRESIZE_AGSV.stagger = ' '
		FIRESIZE_AGSV.FieldType = 104

		FIRESIZE_AGTF = GFED.createVariable('FIRESIZE_AGTF', 'f4', ('Time','south_north', 'west_east'))
		FIRESIZE_AGTF.MemoryOrder = 'XY'
		FIRESIZE_AGTF.description = 'mean firesize of tropical forest from FINN'
		FIRESIZE_AGTF.units = 'm2'
		FIRESIZE_AGTF.stagger = ' '
		FIRESIZE_AGTF.FieldType = 104

		MEAN_FCT_AGEF = GFED.createVariable('MEAN_FCT_AGEF', 'f4', ('Time','south_north', 'west_east'))
		MEAN_FCT_AGEF.MemoryOrder = 'XY'
		MEAN_FCT_AGEF.description = 'mean fraction of extra tropical forest from FINN'
		MEAN_FCT_AGEF.units = 'm2'
		MEAN_FCT_AGEF.stagger = ' '
		MEAN_FCT_AGEF.FieldType = 104

		MEAN_FCT_AGGR = GFED.createVariable('MEAN_FCT_AGGR', 'f4', ('Time','south_north', 'west_east'))
		MEAN_FCT_AGGR.MemoryOrder = 'XY'
		MEAN_FCT_AGGR.description = 'mean fraction of grassland from FINN'
		MEAN_FCT_AGGR.units = 'm2'
		MEAN_FCT_AGGR.stagger = ' '
		MEAN_FCT_AGGR.FieldType = 104

		MEAN_FCT_AGSV = GFED.createVariable('MEAN_FCT_AGSV', 'f4', ('Time','south_north', 'west_east'))
		MEAN_FCT_AGSV.MemoryOrder = 'XY'
		MEAN_FCT_AGSV.description = 'mean fraction of savana from FINN'
		MEAN_FCT_AGSV.units = 'm2'
		MEAN_FCT_AGSV.stagger = ' '
		MEAN_FCT_AGSV.FieldType = 104

		MEAN_FCT_AGTF = GFED.createVariable('MEAN_FCT_AGTF', 'f4', ('Time','south_north', 'west_east'))
		MEAN_FCT_AGTF.MemoryOrder = 'XY'
		MEAN_FCT_AGTF.description = 'mean fraction of tropical forest from FINN'
		MEAN_FCT_AGTF.units = 'm2'
		MEAN_FCT_AGTF.stagger = ' '
		MEAN_FCT_AGTF.FieldType = 104

		co_emiss[:,:,:,:] = CO_emissions_colombia
		no_emiss[:,:,:,:] = NOx_emissions_colombia
		pm25_emiss[:,:,:,:] = PM25_emissions_colombia
		oc_emiss[:,:,:,:] = OC_emissions_colombia
		bc_emiss[:,:,:,:] = BC_emissions_colombia
		so2_emiss[:,:,:,:] = SO2_emissions_colombia
		#c2h6_emiss[:,:,:,:] = C2H6_emissions_colombia
		#ch3oh_emiss[:,:,:,:] = CH3OH_emissions_colombia
		c2h5oh_emiss[:,:,:,:] = C2H5OH_emissions_colombia+C2H6_emissions_colombia
		c3h8_emiss[:,:,:,:] = C3H8_emissions_colombia+CH3OH_emissions_colombia
		c2h4_emiss[:,:,:,:] = C2H4_emissions_colombia
		c3h6_emiss[:,:,:,:] = C3H6_emissions_colombia
		#c10h16_emiss[:,:,:,:] = C10H16_emissions_colombia
		toluene_emiss[:,:,:,:] = toluene_emissions_colombia
		bigene_emiss[:,:,:,:] = bigene_emissions_colombia+C10H16_emissions_colombia
		bigalk_emiss[:,:,:,:] = bigalk_emissions_colombia
		ch2o_emiss[:,:,:,:] = CH2O_emissions_colombia
		nh3_emiss[:,:,:,:] = NH3_emissions_colombia
		ch3cooh_emiss[:,:,:,:] = CH3COOH_emissions_colombia
		mek_emiss[:,:,:,:] = MEK_emissions_colombia
		FIRESIZE_AGEF[:,:,:] = FINN.variables['FIRESIZE_AGEF']
		FIRESIZE_AGGR[:,:,:] = FINN.variables['FIRESIZE_AGGR']
		FIRESIZE_AGSV[:,:,:] = FINN.variables['FIRESIZE_AGSV']
		FIRESIZE_AGTF[:,:,:] = FINN.variables['FIRESIZE_AGTF']
		MEAN_FCT_AGEF[:,:,:] = FINN.variables['MEAN_FCT_AGEF']
		MEAN_FCT_AGGR[:,:,:] = FINN.variables['MEAN_FCT_AGGR']
		MEAN_FCT_AGSV[:,:,:] = FINN.variables['MEAN_FCT_AGSV']
		MEAN_FCT_AGTF[:,:,:] = FINN.variables['MEAN_FCT_AGTF']

		GFED.close()
		
loop_day += 1
		
