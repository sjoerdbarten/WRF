"This script is to create budgets of anthropogenic, biogenic en biomass burning emissions"

import netCDF4 as nc
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import os
np.set_printoptions(threshold=2)

ncfile = nc.Dataset('/home/WUR/barte035/WRFChem/WRFV3/run/wrfout_d01_2014-01-01_00:00:00','r')

daysim = 30
timesim = daysim*24
areacell = (27000.*27000.)*1e-6 #km2
N_IC = 250. #moles emitted per lightning strike IC
N_CG = 250. #moles emitted per lightning strike CG

#ANTHROPOGENIC CALCULATIONS
NOAnthro = np.zeros(shape=(timesim+1,1,99,99))
for i in range(0,timesim):
	NOAnthro[i] = ncfile.variables['E_NO'][i,0,:,:]*14.0067*areacell #g N hr-1
NOAnthro_Total = np.sum(NOAnthro) #g N simulation-1
NOAnthro_Hourly = NOAnthro_Total/timesim #g N hr-1 avg
NOAnthro_Daily = NOAnthro_Hourly*24 #g N day-1
NOAnthro_Yearly = NOAnthro_Daily*1e-9*365	#Gg N yr-1

NO2Anthro = np.zeros(shape=(timesim+1,1,99,99))
for i in range(0,timesim):
	NO2Anthro[i] = ncfile.variables['E_NO2'][i,0,:,:]*14.0067*areacell #g N hr-1
NO2Anthro_Total = np.sum(NO2Anthro) #g N simulation-1
NO2Anthro_Hourly = NO2Anthro_Total/timesim #g N hr-1 avg
NO2Anthro_Daily = NO2Anthro_Hourly*24 #g N day-1
NO2Anthro_Yearly = NO2Anthro_Daily*1e-9*365	#Gg N yr-1

#SOIL NOX CALCULATIONS
NOBio = np.zeros(shape=(timesim+1,1,99,99))
NO2Bio = np.zeros(shape=(timesim+1,1,99,99))
for i in range(0,timesim):
	NOBio[i] = ncfile.variables['EBIO_NO'][i,:,:]*14.0067*areacell #g N hr-1
	NO2Bio[i] = ncfile.variables['EBIO_NO2'][i,:,:]*14.0067*areacell #g N hr-1
NOBio_Total = np.sum(NOBio) #g N simulation-1
NOBio_Hourly = NOBio_Total/timesim #g N hr-1 avg
NOBio_Daily = NOBio_Hourly*24 #g N day-1
NOBio_Yearly = NOBio_Daily*1e-9*365	#Gg N yr-1
NO2Bio_Total = np.sum(NO2Bio) #g N simulation-1
NO2Bio_Hourly = NO2Bio_Total/timesim #g N hr-1 avg
NO2Bio_Daily = NO2Bio_Hourly*24 #g N day-1
NO2Bio_Yearly = NO2Bio_Daily*1e-9*365	#Gg N yr-1

#LIGHTNING NOX CALCULATIONS	
IC_LNOx = ncfile.variables['IC_FLASHCOUNT'][timesim,:,:] #flashes
CG_LNOx = ncfile.variables['CG_FLASHCOUNT'][timesim,:,:] #flashes
IC_LNOx_Total = sum(IC_LNOx)*N_IC*14.0067 #g N simulation-1
CG_LNOx_Total = sum(CG_LNOx)*N_IC*14.0067 #g N simulation-1
LNOx_Total = IC_LNOx_Total+CG_LNOx_Total #g N simulation-1
LNOx_Hourly = LNOx_Total/timesim #g N hr-1 avg
LNOx_Daily = LNOx_Hourly*24 #g N day-1
LNOx_Yearly = LNOx_Daily*365.*1e-9 #Gg N yr-1

'''
#Biomass burning
NOFire = np.zeros(shape=(timesim+1,1,99,99))
NO2Fire = np.zeros(shape=(timesim+1,1,99,99))
for i in range(0,timesim):
	NOFire[i] = ncfile.variables['ebu_no'][timesim,1,:,:]*14.0067*areacell #g N hr-1 (OR ebu_in_no)
	NO2Fire[i] = ncfile.variables['ebu_no2'][timesim,1,:,:]*14.0067*areacell #g N hr-1
NOFire_Total = np.sum(NOFire) #g N simulation-1
NOFire_Daily = NOFire_Total/daysim #g N day-1
NOFire_Yearly = NOFire_Daily*1e-9*365	#Gg N yr-1
NO2Fire_Total = np.sum(NO2Fire) #g N simulation-1
NO2Fire_Daily = NO2Fire_Total/daysim #g N day-1
NO2Fire_Yearly = NO2Fire_Daily*1e-9*365	#Gg N yr-1
'''

#Biomass burning from input files
GFEDpath = '/home/WUR/barte035/WRFChem/GFED/CBMZdaily/'
GFEDfiles = [os.path.join(GFEDpath,filename) for filename in os.listdir(GFEDpath) if filename.startswith('wrffire')]
NOFire = np.zeros(shape=(daysim,1,99,99))
for i in range (0,daysim):
	GFEDfile = nc.Dataset(GFEDfiles[i],'r')
	NOFire[i] = GFEDfile.variables['ebu_in_no'][0,0,:,:]*14.0067*areacell*24 #g N day-1
NOFire_Total = np.sum(NOFire) #g N simulation-1
NOFire_Hourly = NOFire_Total/timesim #g N hr-1 avg
NOFire_Daily = NOFire_Hourly*24 #g N day-1
NOFire_Yearly = NOFire_Daily*1e-9*365	#Gg N yr-1



#PRINT STATEMENTS
print 'Total anthropogenic N flux (NO) =',NOAnthro_Total,'g N over',daysim,'days.'
print 'Average daily anthropogenic N flux (NO) =',NOAnthro_Daily,'g N day-1'
print 'Estimated anthropogenic N budget (NO) =',NOAnthro_Yearly,'Gg N yr-1'
print ''

print 'Total anthropogenic N flux (NO2) =',NO2Anthro_Total,'g N over',daysim,'days.'
print 'Average daily anthropogenic N flux (NO2) =',NO2Anthro_Daily,'g N day-1'
print 'Estimated anthropogenic N budget (NO2) =',NO2Anthro_Yearly,'Gg N yr-1'
print ''
print 'Total anthropogenic N flux (NOx) =',NO2Anthro_Total+NOAnthro_Total,'g N over',daysim,'days.'
print 'Average daily anthropogenic N flux (NOx) =',NO2Anthro_Daily+NOAnthro_Daily,'g N day-1'
print 'Estimated anthropogenic N budget (NOx) =',NO2Anthro_Yearly+NOAnthro_Yearly,'Gg N yr-1'
print ''

print 'Total lightning N simulation duration =',LNOx_Total,'g N over',daysim,'days.'
print 'Average daily lightning N flux =',LNOx_Daily,'g N day-1'
print 'Estimated lightning N budget =',LNOx_Yearly,'Gg N yr-1'
print ''

print 'Total biogenic N flux (NO) =',NOBio_Total,'g N over',daysim,'days.'
print 'Average daily biogenic N flux (NO) =',NOBio_Daily,'g N day-1'
print 'Estimated biogenic N budget (NO) =',NOBio_Yearly,'Gg N yr-1'
print ''
print 'Total biogenic N flux (NO2) =',NO2Bio_Total,'g N over',daysim,'days.'
print 'Average daily biogenic N flux (NO2) =',NO2Bio_Daily,'g N day-1'
print 'Estimated biogenic N budget (NO2) =',NO2Bio_Yearly,'Gg N yr-1'
print ''
print 'Total biogenic N flux (NOx) =',NO2Bio_Total+NOBio_Total,'g N over',daysim,'days.'
print 'Average daily biogenic N flux (NOx) =',NO2Bio_Daily+NOBio_Daily,'g N day-1'
print 'Estimated biogenic N budget (NOx) =',NO2Bio_Yearly+NOBio_Yearly,'Gg N yr-1'
print ''

print 'Total biomass burning N flux (NO) =',NOFire_Total,'g N over',daysim,'days.'
print 'Average daily biomass burning N flux (NO) =',NOFire_Daily,'g N day-1'
print 'Estimated biomass burning N budget (NO) =',NOFire_Yearly,'Gg N yr-1'
print ''
