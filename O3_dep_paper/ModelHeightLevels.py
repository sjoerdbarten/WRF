print 'This script gets the height of the WRF levels in meters. Use heightlevels(wrffile,time,x,y,botlevel,toplevel) where wrffile is the path to the wrffile, time is the time from start of simulation in hours, x is the x gridbox number, y is the y gridbox number, botlevel is the bottom level for which you want to calculate the height for and toplevel is the top level for which you want to calculate the height for.'

import netCDF4 as nc
import math

def heightlevels(wrffile,time,x,y,botlevel,toplevel):
	ncfile = nc.Dataset(wrffile,'r')
	p = ncfile['P'][time,botlevel:toplevel,y,x]
	pb = ncfile['PB'][time,botlevel:toplevel,y,x]
	T = ncfile['T'][time,botlevel:toplevel,y,x]+273.15
	g = 9.81
	Rstar = 8314.
	M = 28.9647
	
	wrf_p = p+pb
	
	print(wrf_p)	
	
	for i in range(botlevel,toplevel,1):
		h = (math.log(wrf_p[0]/wrf_p[i])*Rstar*T[i])/(g*M)
		print(h)
		
	return(h)

heightlevels('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_nudgedBL_fixeddep_d01_2008-08-10_00:00:00',0,150,150,0,33)
