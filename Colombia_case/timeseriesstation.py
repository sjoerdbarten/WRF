import netCDF4 as nc
from pylab import *
import csv
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab
from sklearn.metrics import mean_squared_error,r2_score
from math import sqrt
from scipy.stats import pearsonr


'''
=IF(AND(U2>0.26,U2<0.76),1,2)

STATION:
BAR-Parque Las Aguas 	= 1	O3
BEL-U.S.Buenaventura 	= 2	NO, NO2, O3
C. Alto Rendimiento 	= 14	CO, O3
Cabecera 		= 13	CO
CAL-Corp. Lasallista 	= 3	O3
Carvajal - Sevillana 	= 15	CO, O3
Centro 			= 16	CO, NO, NO2
Ciudadela 		= 17	CO
Compartir 		= 18	O3
Fontibon 		= 19	CO, O3
Guaymaral 		= 20	NO, NO2, O3
ITA-Casa Justicia 	= 4	NO, NO2
ITA-Col. Concejo 	= 5	O3
Kennedy 		= 21	CO, NO, NO2
Las Ferias 		= 22	CO, O3
MED-Museo Antioquia 	= 6	CO
MED-Politecnico JIC 	= 7	NO, NO2
MED-Politecnico JIC (S) = 8	NO
MED-UN Fac. Minas 	= 9	NO, NO2
MED-UN Nucleo Volador 	= 10	NO, NO2, O3
MED-Univ. Medellin (S) 	= 11	O3
MED-Villahermosa 	= 12	O3
MinAmbiente 		= 23	O3
Pance 			= 24	O3
Puente Aranda 		= 25	CO, NO, NO2, O3
San Christobal 		= 26	O3
Suba 			= 27	NO, NO2, O3
Tunal 			= 28	CO, NO, NO2, O3
Universidad del Valle 	= 29	NO2, O3
Usaquen 		= 30	O3

POLLUTANT:
CO 	= 1 (mg/m3)
NO 	= 2 (ug/m3)
NO2 	= 3 (ug/m3)
O3 	= 4 (ug/m3)

UNITS:
ug/m3 	= 1
mg/m3 	= 2

PIXEL GROUPS:
yloc = 42, xloc = 21: 24				Cali
yloc = 43, xloc = 21: 29				Cali
yloc = 43, xloc = 22: 18				Cali
yloc = 48, xloc = 34: 14,15,19,21,22,23,25,26,28	Bogota
yloc = 49, xloc = 34: 20,27,30				Bogota
yloc = 56, xloc = 27: 3,4,5,7,8				Medellin
yloc = 57, xloc = 27: 2,6,9,10,11,12			Medellin
yloc = 57, xloc = 28: 1					Medellin
yloc = 60, xloc = 39: 13,16,17				Buacaramanga
'''

print 'This script reads StationData from 30 ground-stations in Colombia and plots it with WRF-Chem data for the appropriate pixels, use timeseriesstation(variable,city,figsave) to plot any variable at any city. Variables: CO, NO, NO2, O3, NOX. Cities: Cali, Bogota, Bucaramanga, Medellin, BELUS (most rural station). Figsave = 1 saves figures.'

def timeseriesstation(pollutant,city,figsave):
	txtfile = '/archive/ESG/barte035/Colombia/StationData/SISAIRE_201401_CO_NO_NO2_O3_1h.txt'
	csvfile = '/archive/ESG/barte035/Colombia/StationData/SISAIRE_201401_CO_NO_NO2_O3_1h.csv'
	wrffile = '/home/WUR/hoek071/wrfout_d01_2014-01-01_00:00:00'

	StationData = pd.read_csv(csvfile)
	ncfile = nc.Dataset(wrffile,'r')

	daysim = 31
	#timesim = daysim*24-5
	timesim = len(ncfile.variables['XTIME']) #Times or XTIME
	spinup = 24
	
	avg2 = np.zeros(timesim)
	avgday2 = np.zeros(timesim)
	avgnight2 = np.zeros(timesim)
	std2 = np.zeros(timesim)
	stdday2 = np.zeros(timesim)
	stdnight2 = np.zeros(timesim)
	
	#CONVERT STATIONDATA UNIT TO WRF UNIT
	if pollutant == 'CO':
		plotvar = 'co'
		StationData['Concentration'] = (24.45*(StationData['Concentration']*1000.))/28.0101	#mg/m3 to ppb
		pollegend = 'CO'
	if pollutant == 'NO':
		plotvar = 'no'
		StationData['Concentration'] = (24.45*(StationData['Concentration']))/30.00610		#ug/m3 to ppb
		pollegend = 'NO'
	if pollutant == 'NO2':
		plotvar = 'no2'
		StationData['Concentration'] = (24.45*(StationData['Concentration']))/46.00550		#ug/m3 to ppb
		pollegend = 'NO$_2$'
	elif pollutant == 'O3':
		plotvar = 'o3'
		StationData['Concentration'] = (24.45*(StationData['Concentration']))/47.99820		#ug/m3 to ppb
		pollegend = 'O$_3$'

	a = StationData.loc[StationData['Pollutant'] == pollutant]
	
	if pollutant != 'NOX':
		if city == 'Cali':
			for i in range(spinup,timesim,1):
				avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([24,29,18])]
				avg2[i] =  mean(avg1.values[:,3])
				avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([24,29,18]) & a['DagNacht2'].str.match('Dag')]
				avgday2[i] = mean(avgday1.values[:,3])
				avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([24,29,18]) & a['DagNacht2'].str.match('Nacht')]
				avgnight2[i] = mean(avgnight1.values[:,3])
				std2[i] = std(avg1.values[:,3])
				stdday2[i] = std(avgday1.values[:,3])
				stdnight2[i] = std(avgnight1.values[:,3])
			if plotvar == 'o3':
				plotvar = (ncfile.variables[plotvar][:,0,42,21]*1000.+ncfile.variables[plotvar][:,0,43,21]*1000.+ncfile.variables[plotvar][:,0,43,22]*1000.)/3
				st1 = a.loc[a['StationCode'] == 24]
				st2 = a.loc[a['StationCode'] == 29]
				st3 = a.loc[a['StationCode'] == 18]
				leg1 = 'Pance'
				leg2 = 'Universidad del Valle'
				leg3 = 'Compartir'
				legarray = ['WRF-Chem',leg1,leg2,leg3]
				plotoption = 3
			if plotvar == 'no2':
				plotvar = ncfile.variables[plotvar][:,0,43,21]*1000.
				st1 = a.loc[a['StationCode'] == 29]
				leg1 = 'Universidad del Valle'
				legarray = ['WRF-Chem',leg1]
				plotoption = 1

		if city == 'Bucaramanga':
			for i in range(spinup,timesim,1):
				avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([13,16,17])]
				avg2[i] =  mean(avg1.values[:,3])
				avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([13,16,17]) & a['DagNacht2'].str.match('Dag')]
				avgday2[i] = mean(avgday1.values[:,3])
				avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([13,16,17]) & a['DagNacht2'].str.match('Nacht')]
				avgnight2[i] = mean(avgnight1.values[:,3])
				std2[i] = std(avg1.values[:,3])
				stdday2[i] = std(avgday1.values[:,3])
				stdnight2[i] = std(avgnight1.values[:,3])
			if plotvar == 'co':
				plotvar = ncfile.variables[plotvar][:,0,60,39]*1000.
				st1 = a.loc[a['StationCode'] == 13]
				st2 = a.loc[a['StationCode'] == 16]
				st3 = a.loc[a['StationCode'] == 17]
				leg1 = 'Cabecera'
				leg2 = 'Centro'
				leg3 = 'Ciudadela'
				legarray = ['WRF-Chem',leg1,leg2,leg3]
				plotoption = 3
			if plotvar == 'no':
				plotvar = ncfile.variables[plotvar][:,0,60,39]*1000.
				st1 = a.loc[a['StationCode'] == 16]
				leg1 = 'Centro'
				legarray = ['WRF-Chem',leg1]
				plotoption = 1
			if plotvar == 'no2':
				plotvar = ncfile.variables[plotvar][:,0,60,39]*1000.
				st1 = a.loc[a['StationCode'] == 16]
				leg1 = 'Centro'
				legarray = ['WRF-Chem',leg1]
				plotoption = 1					

		if city == 'Bogota':
			for i in range(spinup,timesim,1):
				avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([14,15,19,20,21,22,23,25,26,27,28,30])]
				avg2[i] = mean(avg1.values[:,3])
				avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([14,15,19,20,21,22,23,25,26,27,28,30]) & a['DagNacht2'].str.match('Dag')]
				avgday2[i] = mean(avgday1.values[:,3])
				avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([14,15,19,20,21,22,23,25,26,27,28,30]) & a['DagNacht2'].str.match('Nacht')]			
				avgnight2[i] = mean(avgnight1.values[:,3])
				std2[i] = std(avg1.values[:,3])
				stdday2[i] = std(avgday1.values[:,3])
				stdnight2[i] = std(avgnight1.values[:,3])
			if plotvar == 'co':		
				plotvar = ncfile.variables[plotvar][:,0,48,34]*1000.
				st1 = a.loc[a['StationCode'] == 14]
				st2 = a.loc[a['StationCode'] == 15]
				st3 = a.loc[a['StationCode'] == 19]
				st4 = a.loc[a['StationCode'] == 21]
				st5 = a.loc[a['StationCode'] == 22]
				st6 = a.loc[a['StationCode'] == 25]
				st7 = a.loc[a['StationCode'] == 28]
				leg1 = 'C. Alto Rendimiento'
				leg2 = 'Carvajal - Sevillana'
				leg3 = 'Fontibon'
				leg4 = 'Kennedy'
				leg5 = 'Las Ferias'
				leg6 = 'Puente Aranda'
				leg7 = 'Tunal'
				legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5,leg6,leg7]	
				plotoption = 7	
			if plotvar == 'o3':
			 	plotvar = (ncfile.variables[plotvar][:,0,48,34]*1000.*8+ncfile.variables[plotvar][:,0,49,34]*1000.*3)/11
				#r = np.array([ncfile.variables[plotvar][:,0,48,34]*1000.,ncfile.variables[plotvar][:,0,48,34]*1000.,ncfile.variables[plotvar][:,0,48,34]*1000.,ncfile.variables[plotvar][:,0,48,34]*1000.,ncfile.variables[plotvar][:,0,48,34]*1000.,ncfile.variables[plotvar][:,0,48,34]*1000.,ncfile.variables[plotvar][:,0,48,34]*1000.,ncfile.variables[plotvar][:,0,48,34]*1000.,ncfile.variables[plotvar][:,0,49,34]*1000.,ncfile.variables[plotvar][:,0,49,34]*1000.,ncfile.variables[plotvar][:,0,49,34]*1000.])
				#plotvar = np.mean(r,axis=0)
				#plotvarstd = np.std(r,axis=0)
				st1 = a.loc[a['StationCode'] == 14]
				st2 = a.loc[a['StationCode'] == 15]
				st3 = a.loc[a['StationCode'] == 19]
				st4 = a.loc[a['StationCode'] == 20]
				st5 = a.loc[a['StationCode'] == 22]
				st6 = a.loc[a['StationCode'] == 23]
				st7 = a.loc[a['StationCode'] == 25]
				st8 = a.loc[a['StationCode'] == 26]
				st9 = a.loc[a['StationCode'] == 27]
				st10 = a.loc[a['StationCode'] == 28]
				st11 = a.loc[a['StationCode'] == 30]
				leg1 = 'C. Alto Rendimiento'
				leg2 = 'Carvajal - Sevillana'
				leg3 = 'Fontibon'
				leg4 = 'Guaymaral'
				leg5 = 'Las Ferias'
				leg6 = 'MinAmbiente'
				leg7 = 'Puente Aranda'
				leg8 = 'San Christobal'
				leg9 = 'Suba'
				leg10 = 'Tunal'
				leg11 = 'Usaquen'
				legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5,leg6,leg7,leg8,leg9,leg10,leg11]
				plotoption = 11
			if plotvar == 'no':
			 	plotvar = (ncfile.variables[plotvar][:,0,48,34]*1000.*3+ncfile.variables[plotvar][:,0,49,34]*1000.*2)/5
				st1 = a.loc[a['StationCode'] == 20]
				st2 = a.loc[a['StationCode'] == 21]
				st3 = a.loc[a['StationCode'] == 25]
				st4 = a.loc[a['StationCode'] == 27]
				st5 = a.loc[a['StationCode'] == 28]
				leg1 = 'Guaymaral'
				leg2 = 'Kennedy'
				leg3 = 'Puente Aranda'
				leg4 = 'Suba'
				leg5 = 'Tunal'
				legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5]
				plotoption = 5		
			if plotvar == 'no2':
			 	plotvar = (ncfile.variables[plotvar][:,0,48,34]*1000.*3+ncfile.variables[plotvar][:,0,49,34]*1000.*2)/5
				st1 = a.loc[a['StationCode'] == 20]
				st2 = a.loc[a['StationCode'] == 21]
				st3 = a.loc[a['StationCode'] == 25]
				st4 = a.loc[a['StationCode'] == 27]
				st5 = a.loc[a['StationCode'] == 28]
				leg1 = 'Guaymaral'
				leg2 = 'Kennedy'
				leg3 = 'Puente Aranda'
				leg4 = 'Suba'
				leg5 = 'Tunal'
				legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5]
				plotoption = 5

		if city == 'Medellin':
			for i in range(spinup,timesim,1):
				avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([3,4,5,7,8,2,6,9,10,11,12,1])]
				avg2[i] = mean(avg1.values[:,3])
				avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([3,4,5,7,8,2,6,9,10,11,12,1]) & a['DagNacht2'].str.match('Dag')]
				avgday2[i] = mean(avgday1.values[:,3])
				avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([3,4,5,7,8,2,6,9,10,11,12,1]) & a['DagNacht2'].str.match('Nacht')]			
				avgnight2[i] = mean(avgnight1.values[:,3])
				std2[i] = std(avg1.values[:,3])
				stdday2[i] = std(avgday1.values[:,3])
				stdnight2[i] = std(avgnight1.values[:,3])
			if plotvar == 'co':
				plotvar = ncfile.variables[plotvar][:,0,57,27]*1000.
				st1 = a.loc[a['StationCode'] == 6]
				leg1 = 'MED-Museo Antioquia'
				legarray = ['WRF-Chem',leg1]
				plotoption = 1
			if plotvar == 'o3':
				plotvar = (ncfile.variables[plotvar][:,0,57,27]*1000.*4+ncfile.variables[plotvar][:,0,57,28]*1000.+ncfile.variables[plotvar][:,0,56,27]*1000.*2)/7
				st1 = a.loc[a['StationCode'] == 1]
				st2 = a.loc[a['StationCode'] == 2]
				st3 = a.loc[a['StationCode'] == 3]
				st4 = a.loc[a['StationCode'] == 5]
				st5 = a.loc[a['StationCode'] == 10]
				st6 = a.loc[a['StationCode'] == 11]
				st7 = a.loc[a['StationCode'] == 12]
				leg1 = 'BAR-Parque Las Aguas'
				leg2 = 'BEL-U.S.Buenaventura'
				leg3 = 'CAL-Corp. Lasallista'
				leg4 = 'ITA-Col. Concejo'
				leg5 = 'MED-UN Nucleo Volador'
				leg6 = 'MED-Univ. Medellin (S)'
				leg7 = 'MED-Villahermosa'
				legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5,leg6,leg7]
				plotoption = 7
			if plotvar == 'no':
				plotvar = (ncfile.variables[plotvar][:,0,57,27]*1000.*3+ncfile.variables[plotvar][:,0,56,27]*1000.*3)/6
				st1 = a.loc[a['StationCode'] == 2]
				st2 = a.loc[a['StationCode'] == 4]
				st3 = a.loc[a['StationCode'] == 7]
				st4 = a.loc[a['StationCode'] == 8]
				st5 = a.loc[a['StationCode'] == 9]
				st6 = a.loc[a['StationCode'] == 10]
				leg1 = 'BEL-U.S.Buenaventura'
				leg2 = 'ITA-Casa Justicia'
				leg3 = 'MED-Politecnico JIC'
				leg4 = 'MED-Politecnico JIC (S)'
				leg5 = 'MED-UN Fac. Minas'
				leg6 = 'MED-UN Nucleo Volador'
				legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5,leg6]
				plotoption = 6		
			if plotvar == 'no2':
				plotvar = (ncfile.variables[plotvar][:,0,57,27]*1000.*3+ncfile.variables[plotvar][:,0,56,27]*1000.*2)/5
				st1 = a.loc[a['StationCode'] == 2]
				st2 = a.loc[a['StationCode'] == 4]
				st3 = a.loc[a['StationCode'] == 7]
				st4 = a.loc[a['StationCode'] == 9]
				st5 = a.loc[a['StationCode'] == 10]
				leg1 = 'BEL-U.S.Buenaventura'
				leg2 = 'ITA-Casa Justicia'
				leg3 = 'MED-Politecnico JIC'
				leg4 = 'MED-UN Fac. Minas'
				leg5 = 'MED-UN Nucleo Volador'
				legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5]
				plotoption = 5
				
		if city == 'BELUS':
			for i in range(spinup,timesim,1):
				avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([2])]
				avg2[i] = mean(avg1.values[:,3])
				avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([2]) & a['DagNacht2'].str.match('Dag')]
				avgday2[i] = mean(avgday1.values[:,3])
				avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([2]) & a['DagNacht2'].str.match('Nacht')]			
				avgnight2[i] = mean(avgnight1.values[:,3])
				std2[i] = std(avg1.values[:,3])
				stdday2[i] = std(avgday1.values[:,3])
				stdnight2[i] = std(avgnight1.values[:,3])
			if plotvar == 'o3':
				plotvar = ncfile.variables[plotvar][:,0,57,27]*1000.
				st1 = a.loc[a['StationCode'] == 2]
				leg1 = 'BEL-U.S.Buenaventura'
				legarray = ['WRF-Chem',leg1]
				plotoption = 1
			if plotvar == 'no':
				plotvar = ncfile.variables[plotvar][:,0,57,27]*1000.
				st1 = a.loc[a['StationCode'] == 2]
				leg1 = 'BEL-U.S.Buenaventura'
				legarray = ['WRF-Chem',leg1]
				plotoption = 1		
			if plotvar == 'no2':
				plotvar = ncfile.variables[plotvar][:,0,57,27]*1000.
				st1 = a.loc[a['StationCode'] == 2]
				leg1 = 'BEL-U.S.Buenaventura'
				legarray = ['WRF-Chem',leg1]
				plotoption = 1
							
	if pollutant == 'NOX':
		plotvar = 'nox'
		pollegend = 'NO$_x$'
		n = StationData.loc[StationData['Pollutant'] == 'NO']
		n['Concentration'] = (24.45*(StationData['Concentration']))/30.00610		#ug/m3 to ppb
		m = StationData.loc[StationData['Pollutant'] == 'NO2']
		m['Concentration'] = (24.45*(StationData['Concentration']))/46.00550		#ug/m3 to ppb
		for i in range(spinup,timesim,1):
			if city == 'Bucaramanga':
				avgn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([16])]
				avgm1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([16])]
				avgdayn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([16]) & n['DagNacht2'].str.match('Dag')]
				avgdaym1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([16]) & m['DagNacht2'].str.match('Dag')]
				avgnightn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([16]) & n['DagNacht2'].str.match('Nacht')]			
				avgnightm1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([16]) & m['DagNacht2'].str.match('Nacht')]			
				if avgn1.values[:,3] > 0 and avgm1.values[:,3] > 0:
					avgo1 = avgn1.values[:,3]+avgm1.values[:,3]
					avg2[i] = mean(avgo1)
					std2[i] = std(avgo1)
				if avgdayn1.values[:,3] > 0 and avgdaym1.values[:,3] > 0:
					avgdayo1 = avgdayn1.values[:,3]+avgdaym1.values[:,3]
					avgday2[i] = mean(avgdayo1)
					stdday2[i] = std(avgdayo1)
				if avgnightn1.values[:,3] > 0 and avgnightm1.values[:,3] > 0:
					avgnighto1 = avgnightn1.values[:,3]+avgnightm1.values[:,3]
					avgnight2[i] = mean(avgnighto1)
					stdnight2[i] = std(avgnighto1)	
			if city == 'Bogota':
				avgn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([20,21,25,27,28])]
				avgm1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([20,21,25,27,28])]
				avgdayn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([20,21,25,27,28]) & n['DagNacht2'].str.match('Dag')]
				avgdaym1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([20,21,25,27,28]) & m['DagNacht2'].str.match('Dag')]
				avgnightn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([20,21,25,27,28]) & n['DagNacht2'].str.match('Nacht')]			
				avgnightm1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([20,21,25,27,28]) & m['DagNacht2'].str.match('Nacht')]				
				avgj1 = np.zeros(min([avgn1.shape[0],avgm1.shape[0]]))
				avgdayj1 = np.zeros(min([avgdayn1.shape[0],avgdaym1.shape[0]]))
				avgnightj1 = np.zeros(min([avgnightn1.shape[0],avgnightm1.shape[0]]))
				for j in range(0,min([avgn1.shape[0],avgm1.shape[0]]),1):
					if avgn1.values[j,3] > 0 and avgm1.values[j,3] > 0:
						avgj1[j] = avgn1.values[j,3]+avgm1.values[j,3]
					if avgn1.values[0,21] == 'Dag':
						if avgdayn1.values[j,3] > 0 and avgdaym1.values[j,3] > 0:
							avgdayj1[j] = avgdayn1.values[j,3]+avgdaym1.values[j,3]
					if avgn1.values[0,21] == 'Nacht':
						if avgnightn1.values[j,3] > 0 and avgnightm1.values[j,3] > 0:
							avgnightj1[j] = avgnightn1.values[j,3]+avgnightm1.values[j,3]
				avg2[i] = mean(avgj1)
				std2[i] = std(avgj1)
				if avgn1.values[0,21] == 'Dag':
					avgday2[i] = mean(avgdayj1)
					stdday2[i] = std(avgdayj1)
				if avgn1.values[0,21] == 'Nacht':
					avgnight2[i] = mean(avgnightj1)
					stdnight2[i] = std(avgnightj1)
			if city == 'Medellin':
				avgn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([2,4,7,9,10])]
				avgm1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([2,4,7,9,10])]
				avgdayn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([2,4,7,9,10]) & n['DagNacht2'].str.match('Dag')]
				avgdaym1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([2,4,7,9,10]) & m['DagNacht2'].str.match('Dag')]
				avgnightn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([2,4,7,9,10]) & n['DagNacht2'].str.match('Nacht')]			
				avgnightm1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([2,4,7,9,10]) & m['DagNacht2'].str.match('Nacht')]			
				avgj1 = np.zeros(min([avgn1.shape[0],avgm1.shape[0]]))
				avgdayj1 = np.zeros(min([avgdayn1.shape[0],avgdaym1.shape[0]]))
				avgnightj1 = np.zeros(min([avgnightn1.shape[0],avgnightm1.shape[0]]))
				for j in range(0,min([avgn1.shape[0],avgm1.shape[0]]),1):
					if avgn1.values[j,3] > 0 and avgm1.values[j,3] > 0:
						avgj1[j] = avgn1.values[j,3]+avgm1.values[j,3]
					if avgn1.values[0,21] == 'Dag':
						if avgdayn1.values[j,3] > 0 and avgdaym1.values[j,3] > 0:
							avgdayj1[j] = avgdayn1.values[j,3]+avgdaym1.values[j,3]
					if avgn1.values[0,21] == 'Nacht':
						if avgnightn1.values[j,3] > 0 and avgnightm1.values[j,3] > 0:
							avgnightj1[j] = avgnightn1.values[j,3]+avgnightm1.values[j,3]
				avg2[i] = mean(avgj1)
				std2[i] = std(avgj1)
				if avgn1.values[0,21] == 'Dag':
					avgday2[i] = mean(avgdayj1)
					stdday2[i] = std(avgdayj1)
				if avgn1.values[0,21] == 'Nacht':
					avgnight2[i] = mean(avgnightj1)
					stdnight2[i] = std(avgnightj1)
			if city == 'BELUS':
				avgn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([2])]
				avgm1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([2])]
				avgdayn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([2]) & n['DagNacht2'].str.match('Dag')]
				avgdaym1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([2]) & m['DagNacht2'].str.match('Dag')]
				avgnightn1 = n.loc[(n['DateConvert2'] == i) & n['StationCode'].isin([2]) & n['DagNacht2'].str.match('Nacht')]			
				avgnightm1 = m.loc[(m['DateConvert2'] == i) & m['StationCode'].isin([2]) & m['DagNacht2'].str.match('Nacht')]			
				if avgn1.values[:,3] > 0 and avgm1.values[:,3] > 0:
					avgo1 = avgn1.values[:,3]+avgm1.values[:,3]
					avg2[i] = mean(avgo1)
					std2[i] = std(avgo1)
				if avgdayn1.values[:,3] > 0 and avgdaym1.values[:,3] > 0:
					avgdayo1 = avgdayn1.values[:,3]+avgdaym1.values[:,3]
					avgday2[i] = mean(avgdayo1)
					stdday2[i] = std(avgdayo1)
				if avgnightn1.values[:,3] > 0 and avgnightm1.values[:,3] > 0:
					avgnighto1 = avgnightn1.values[:,3]+avgnightm1.values[:,3]
					avgnight2[i] = mean(avgnighto1)
					stdnight2[i] = std(avgnighto1)
			
					
		if city == 'Bucaramanga':
			plotvar = ncfile.variables['no'][:,0,60,39]*1000.+ncfile.variables['no2'][:,0,60,39]*1000.
			stn1 = n.loc[n['StationCode'] == 16]
			stm1 = m.loc[m['StationCode'] == 16]
			leg1 = 'Centro'
			legarray = ['WRF-Chem',leg1]
			plotoption = 1
			stotime1 = np.zeros(1000)
			stoconc1 = np.zeros(1000)
			for i in range(spinup,timesim,1):
				stn1row = stn1.loc[stn1['DateConvert2'] == i]
				stm1row = stm1.loc[stm1['DateConvert2'] == i]
				if (stn1row.values[:,3] > 0 and stm1row.values[:,3] > 0):
					stoconc1[i] = stn1row.values[:,3]+stm1row.values[:,3]
					stotime1[i] = stm1row.values[:,2]
		if city == 'Bogota':
			plotvar = ((ncfile.variables['no'][:,0,48,34]*1000.*3+ncfile.variables['no'][:,0,49,34]*1000.*2)/5)+((ncfile.variables['no2'][:,0,48,34]*1000.*3+ncfile.variables['no2'][:,0,49,34]*1000.*2)/5)
			stn1 = n.loc[n['StationCode'] == 20]
			stn2 = n.loc[n['StationCode'] == 21]
			stn3 = n.loc[n['StationCode'] == 25]
			stn4 = n.loc[n['StationCode'] == 27]
			stn5 = n.loc[n['StationCode'] == 28]
			stm1 = m.loc[m['StationCode'] == 20]
			stm2 = m.loc[m['StationCode'] == 21]
			stm3 = m.loc[m['StationCode'] == 25]
			stm4 = m.loc[m['StationCode'] == 27]
			stm5 = m.loc[m['StationCode'] == 28]
			leg1 = 'Guaymaral'
			leg2 = 'Kennedy'
			leg3 = 'Puente Aranda'
			leg4 = 'Suba'
			leg5 = 'Tunal'
			legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5]
			plotoption = 5
			stotime1 = np.zeros(1000)
			stoconc1 = np.zeros(1000)
			stotime2 = np.zeros(1000)
			stoconc2 = np.zeros(1000)
			stotime3 = np.zeros(1000)
			stoconc3 = np.zeros(1000)
			stotime4 = np.zeros(1000)
			stoconc4 = np.zeros(1000)
			stotime5 = np.zeros(1000)
			stoconc5 = np.zeros(1000)
			for i in range(spinup,timesim,1):
				stn1row = stn1.loc[stn1['DateConvert2'] == i]
				stm1row = stm1.loc[stm1['DateConvert2'] == i]
				if (stn1row.values[:,3] > 0 and stm1row.values[:,3] > 0):
					stoconc1[i] = stn1row.values[:,3]+stm1row.values[:,3]
					stotime1[i] = stm1row.values[:,2]
				stn2row = stn2.loc[stn2['DateConvert2'] == i]
				stm2row = stm2.loc[stm2['DateConvert2'] == i]
				if (stn2row.values[:,3] > 0 and stm2row.values[:,3] > 0):
					stoconc2[i] = stn2row.values[:,3]+stm2row.values[:,3]
					stotime2[i] = stm2row.values[:,2]
				stn3row = stn3.loc[stn3['DateConvert2'] == i]
				stm3row = stm3.loc[stm3['DateConvert2'] == i]
				if (stn3row.values[:,3] > 0 and stm3row.values[:,3] > 0):
					stoconc3[i] = stn3row.values[:,3]+stm3row.values[:,3]
					stotime3[i] = stm3row.values[:,2]
				stn4row = stn4.loc[stn4['DateConvert2'] == i]
				stm4row = stm4.loc[stm4['DateConvert2'] == i]
				if (stn4row.values[:,3] > 0 and stm4row.values[:,3] > 0):
					stoconc4[i] = stn4row.values[:,3]+stm4row.values[:,3]
					stotime4[i] = stm4row.values[:,2]
				stn5row = stn5.loc[stn5['DateConvert2'] == i]
				stm5row = stm5.loc[stm5['DateConvert2'] == i]
				if (stn5row.values[:,3] > 0 and stm5row.values[:,3] > 0):
					stoconc5[i] = stn5row.values[:,3]+stm5row.values[:,3]
					stotime5[i] = stm5row.values[:,2]
		if city == 'Medellin':
			plotvar = ((ncfile.variables['no'][:,0,57,27]*1000.*3+ncfile.variables['no'][:,0,56,27]*1000.*2)/5)+((ncfile.variables['no2'][:,0,57,27]*1000.*3+ncfile.variables['no2'][:,0,56,27]*1000.*2)/5)		
			stn1 = n.loc[n['StationCode'] == 2]
			stn2 = n.loc[n['StationCode'] == 4]
			stn3 = n.loc[n['StationCode'] == 7]
			stn4 = n.loc[n['StationCode'] == 9]
			stn5 = n.loc[n['StationCode'] == 10]
			stm1 = m.loc[m['StationCode'] == 2]
			stm2 = m.loc[m['StationCode'] == 4]
			stm3 = m.loc[m['StationCode'] == 7]
			stm4 = m.loc[m['StationCode'] == 9]
			stm5 = m.loc[m['StationCode'] == 10]
			leg1 = 'BEL-U.S.Buenaventura'
			leg2 = 'ITA-Casa Justicia'
			leg3 = 'MED-Politecnico JIC'
			leg4 = 'MED-UN Fac. Minas'
			leg5 = 'MED-UN Nucleo Volador'
			legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5]
			plotoption = 5
			stotime1 = np.zeros(1000)
			stoconc1 = np.zeros(1000)
			stotime2 = np.zeros(1000)
			stoconc2 = np.zeros(1000)
			stotime3 = np.zeros(1000)
			stoconc3 = np.zeros(1000)
			stotime4 = np.zeros(1000)
			stoconc4 = np.zeros(1000)
			stotime5 = np.zeros(1000)
			stoconc5 = np.zeros(1000)
			for i in range(spinup,timesim,1):
				stn1row = stn1.loc[stn1['DateConvert2'] == i]
				stm1row = stm1.loc[stm1['DateConvert2'] == i]
				if (stn1row.values[:,3] > 0 and stm1row.values[:,3] > 0):
					stoconc1[i] = stn1row.values[:,3]+stm1row.values[:,3]
					stotime1[i] = stm1row.values[:,2]
				stn2row = stn2.loc[stn2['DateConvert2'] == i]
				stm2row = stm2.loc[stm2['DateConvert2'] == i]
				if (stn2row.values[:,3] > 0 and stm2row.values[:,3] > 0):
					stoconc2[i] = stn2row.values[:,3]+stm2row.values[:,3]
					stotime2[i] = stm2row.values[:,2]
				stn3row = stn3.loc[stn3['DateConvert2'] == i]
				stm3row = stm3.loc[stm3['DateConvert2'] == i]
				if (stn3row.values[:,3] > 0 and stm3row.values[:,3] > 0):
					stoconc3[i] = stn3row.values[:,3]+stm3row.values[:,3]
					stotime3[i] = stm3row.values[:,2]
				stn4row = stn4.loc[stn4['DateConvert2'] == i]
				stm4row = stm4.loc[stm4['DateConvert2'] == i]
				if (stn4row.values[:,3] > 0 and stm4row.values[:,3] > 0):
					stoconc4[i] = stn4row.values[:,3]+stm4row.values[:,3]
					stotime4[i] = stm4row.values[:,2]
				stn5row = stn5.loc[stn5['DateConvert2'] == i]
				stm5row = stm5.loc[stm5['DateConvert2'] == i]
				if (stn5row.values[:,3] > 0 and stm5row.values[:,3] > 0):
					stoconc5[i] = stn5row.values[:,3]+stm5row.values[:,3]
					stotime5[i] = stm5row.values[:,2]
		
		if city == 'BELUS':
			plotvar = ncfile.variables['no'][:,0,57,27]*1000.+ncfile.variables['no2'][:,0,57,27]*1000.
			stn1 = n.loc[n['StationCode'] == 2]
			stm1 = m.loc[m['StationCode'] == 2]
			leg1 = 'BEL-U.S.Buenaventura'
			legarray = ['WRF-Chem',leg1]
			plotoption = 1
			stotime1 = np.zeros(1000)
			stoconc1 = np.zeros(1000)
			for i in range(spinup,timesim,1):
				stn1row = stn1.loc[stn1['DateConvert2'] == i]
				stm1row = stm1.loc[stm1['DateConvert2'] == i]
				if (stn1row.values[:,3] > 0 and stm1row.values[:,3] > 0):
					stoconc1[i] = stn1row.values[:,3]+stm1row.values[:,3]
					stotime1[i] = stm1row.values[:,2]

	'''	
	#SELECT YLOC,XLOC
	if station in [24]:
		yloc = 42
		xloc = 21
		st1 = a.loc[a['StationCode'] == 24]
		leg1 = 'Pance'
		legarray = ['WRF-Chem',leg1]
		for i in range(spinup,timesim,1):
			avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([24])]
			avg2[i] =  mean(avg1.values[:,3])
			avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([24]) & a['DagNacht2'].str.match('Dag')]
			avgday2[i] = mean(avgday1.values[:,3])
			avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([24]) & a['DagNacht2'].str.match('Nacht')]
			avgnight2[i] = mean(avgnight1.values[:,3])
	if station in [29]:
		yloc = 43
		xloc = 21
		st1 = a.loc[a['StationCode'] == 29]
		leg1 = 'Universidad del Valle'
		legarray = ['WRF-Chem',leg1]
		for i in range(spinup,timesim,1):
			avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([29])]
			avg2[i] =  mean(avg1.values[:,3])
			avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([29]) & a['DagNacht2'].str.match('Dag')]
			avgday2[i] = mean(avgday1.values[:,3])
			avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([29]) & a['DagNacht2'].str.match('Nacht')]
			avgnight2[i] = mean(avgnight1.values[:,3])
	if station in [18]:
		yloc = 42
		xloc = 21
		st1 = a.loc[a['StationCode'] == 18]
		leg1 =  'Compartir'
		legarray = ['WRF-Chem',leg1]
		for i in range(spinup,timesim,1):
			avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([18])]
			avg2[i] =  mean(avg1.values[:,3])
			avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([18]) & a['DagNacht2'].str.match('Dag')]
			avgday2[i] = mean(avgday1.values[:,3])
			avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([18]) & a['DagNacht2'].str.match('Nacht')]
			avgnight2[i] = mean(avgnight1.values[:,3])
	if station in [14,15,19,21,22,23,25,26,28]:
		yloc = 48
		xloc = 34
		st1 = a.loc[a['StationCode'] == 14]
		st2 = a.loc[a['StationCode'] == 15]
		st3 = a.loc[a['StationCode'] == 19]
		st4 = a.loc[a['StationCode'] == 21]
		st5 = a.loc[a['StationCode'] == 22]
		st6 = a.loc[a['StationCode'] == 23]
		st7 = a.loc[a['StationCode'] == 25]
		st8 = a.loc[a['StationCode'] == 26]
		st9 = a.loc[a['StationCode'] == 28]
		leg1 = 'C. Alto Rendimiento'
		leg2 = 'Carvajal - Sevillana'
		leg3 = 'Fontibon'
		leg4 = 'Kennedy'
		leg5 = 'Las Ferias'
		leg6 = 'MinAmbiente'
		leg7 = 'Puente Aranda'
		leg8 = 'San Christobal'
		leg9 = 'Tunal'
		legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5,leg6,leg7,leg8,leg9]
		for i in range(spinup,timesim,1):
			avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([14,15,19,21,22,23,25,26,28])]
			avg2[i] = mean(avg1.values[:,3])
			avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([14,15,19,21,22,23,25,26,28]) & a['DagNacht2'].str.match('Dag')]
			avgday2[i] = mean(avgday1.values[:,3])
			avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([14,15,19,21,22,23,25,26,28]) & a['DagNacht2'].str.match('Nacht')]			
			avgnight2[i] = mean(avgnight1.values[:,3])
	if station in [20,27,30]:
		yloc = 49
		xloc = 34
		st1 = a.loc[a['StationCode'] == 20]
		st2 = a.loc[a['StationCode'] == 27]
		st3 = a.loc[a['StationCode'] == 30]
		leg1 = 'Guaymaral'
		leg2 = 'Suba'
		leg3 = 'Usaquen'
		legarray = ['WRF-Chem',leg1,leg2,leg3]
		for i in range(spinup,timesim,1):
			avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([20,27,30])]    839         plt.plot([0,100000],[0,100000],'r-')
			avg2[i] =  mean(avg1.values[:,3])
			avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([20,27,30]) & a['DagNacht2'].str.match('Dag')]
			avgday2[i] = mean(avgday1.values[:,3])
			avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([20,27,30]) & a['DagNacht2'].str.match('Nacht')]
			avgnight2[i] = mean(avgnight1.values[:,3])
	if station in [3,4,5,7,8]:
		yloc = 56
		xloc = 27
		st1 = a.loc[a['StationCode'] == 3]
		st2 = a.loc[a['StationCode'] == 4]
		st3 = a.loc[a['StationCode'] == 5]
		st4 = a.loc[a['StationCode'] == 7]
		st5 = a.loc[a['StationCode'] == 8]
		leg1 = 'CAL-Corp. Lasallista'
		leg2 = 'ITA-Casa Justicia'
		leg3 = 'ITA-Col. Concejo'
		leg4 = 'MED-Politecnico JIC'
		leg5 = 'MED-Politecnico JIC (S)'
		legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5]
		for i in range(spinup,timesim,1):
			avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([3,4,5,7,8])]
			avg2[i] =  mean(avg1.values[:,3])
			avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([3,4,5,7,8]) & a['DagNacht2'].str.match('Dag')]
			avgday2[i] = mean(avgday1.values[:,3])
			avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([3,4,5,7,8]) & a['DagNacht2'].str.match('Nacht')]
			avgnight2[i] = mean(avgnight1.values[:,3])
	if station in [2,6,9,10,11,12]:
		yloc = 57
		xloc = 27
		st1 = a.loc[a['StationCode'] == 2]
		st2 = a.loc[a['StationCode'] == 6]
		st3 = a.loc[a['StationCode'] == 9]
		st4 = a.loc[a['StationCode'] == 10]
		st5 = a.loc[a['StationCode'] == 11]
		st6 = a.loc[a['StationCode'] == 12]
		leg1 = 'BEL-U.S.Buenaventura'
		leg2 = 'MED-Museo Antioquia'
		leg3 = 'MED-UN Fac. Minas'
		leg4 = 'MED-UN Nucleo Volador'
		leg5 = 'MED-Univ. Medellin (S)'
		leg6 = 'MED-Villahermosa'
		legarray = ['WRF-Chem',leg1,leg2,leg3,leg4,leg5,leg6]
		for i in range(spinup,timesim,1):
			avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([2,6,9,10,11,12])]
			avg2[i] =  mean(avg1.values[:,3])
			avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([2,6,9,10,11,12]) & a['DagNacht2'].str.match('Dag')]
			avgday2[i] = mean(avgday1.values[:,3])
			avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([2,6,9,10,11,12]) & a['DagNacht2'].str.match('Nacht')]
			avgnight2[i] = mean(avgnight1.values[:,3])
	if station in [1]:
		yloc = 57
		xloc = 28
		st1 = a.loc[a['StationCode'] == 1]
		leg1 =  'BAR-Parque Las Aguas'
		legarray = ['WRF-Chem',leg1]
		for i in range(spinup,timesim,1):
			avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([1])]
			avg2[i] =  mean(avg1.values[:,3])
			avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([1]) & a['DagNacht2'].str.match('Dag')]
			avgday2[i] = mean(avgday1.values[:,3])
			avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([1]) & a['DagNacht2'].str.match('Nacht')]
			avgnight2[i] = mean(avgnight1.values[:,3])
	if station in [13,16,17]:
		yloc = 60
		xloc = 39
		st1 = a.loc[a['StationCode'] == 13]

		st2 = a.loc[a['StationCode'] == 16]
		st3 = a.loc[a['StationCode'] == 17]
		leg1 = 'Cabacera'
		leg2 = 'Centro'
		leg3 = 'Ciudadela'
		legarray = ['WRF-Chem',leg1,leg2,leg3]
		for i in range(spinup,timesim,1):
			avg1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([13,16,17])]
			avg2[i] =  mean(avg1.values[:,3])
			avgday1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([13,16,17]) & a['DagNacht2'].str.match('Dag')]
			avgday2[i] = mean(avgday1.values[:,3])
			avgnight1 = a.loc[(a['DateConvert2'] == i) & a['StationCode'].isin([13,16,17]) & a['DagNacht2'].str.match('Nacht')]
			avgnight2[i] = mean(avgnight1.values[:,3])
	'''
	
	max1 = max(plotvar)
	max2 = max(avg2)
	max3 = max([max1,max2])
	
	timeplotvar = np.arange(0,timesim,1)/24.

	#for i in range(0,plotvarstd.shape[0],1):
	#	if i not in [20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740]:
	#		plotvarstd[i] = 0

	
	plt.plot(timeplotvar,plotvar,linewidth=3,color='black',zorder=3)
	if plotoption == 1:
		if pollutant != 'NOX':
			plt.plot(st1['DateConvert2']/24., st1['Concentration'], 'ro',markersize=3,markeredgewidth=0.0)
		if pollutant == 'NOX':
			plt.plot(stotime1/24., stoconc1, 'ro',markersize=3,markeredgewidth=0.0)
	if plotoption == 3:
		plt.plot(st1['DateConvert2']/24., st1['Concentration'], 'ro',markersize=3,markeredgewidth=0.0)
		plt.plot(st2['DateConvert2']/24., st2['Concentration'], 'bo',markersize=3,markeredgewidth=0.0)
		plt.plot(st3['DateConvert2']/24., st3['Concentration'], 'go',markersize=3,markeredgewidth=0.0)
	if plotoption == 5:
		if pollutant != 'NOX':
			plt.plot(st1['DateConvert2']/24., st1['Concentration'], 'ro',markersize=3,markeredgewidth=0.0)
			plt.plot(st2['DateConvert2']/24., st2['Concentration'], 'bo',markersize=3,markeredgewidth=0.0)
			plt.plot(st3['DateConvert2']/24., st3['Concentration'], 'go',markersize=3,markeredgewidth=0.0)
			plt.plot(st4['DateConvert2']/24., st4['Concentration'], 'co',markersize=3,markeredgewidth=0.0)
			plt.plot(st5['DateConvert2']/24., st5['Concentration'], 'mo',markersize=3,markeredgewidth=0.0)
		if pollutant == 'NOX':
			plt.plot(stotime1/24., stoconc1, 'ro',markersize=3,markeredgewidth=0.0)
			plt.plot(stotime2/24., stoconc2, 'bo',markersize=3,markeredgewidth=0.0)
			plt.plot(stotime3/24., stoconc3, 'go',markersize=3,markeredgewidth=0.0)
			plt.plot(stotime4/24., stoconc4, 'co',markersize=3,markeredgewidth=0.0)
			plt.plot(stotime5/24., stoconc5, 'mo',markersize=3,markeredgewidth=0.0)
	if plotoption == 6:
		plt.plot(st1['DateConvert2']/24., st1['Concentration'], 'ro',markersize=3,markeredgewidth=0.0)
		plt.plot(st2['DateConvert2']/24., st2['Concentration'], 'bo',markersize=3,markeredgewidth=0.0)
		plt.plot(st3['DateConvert2']/24., st3['Concentration'], 'go',markersize=3,markeredgewidth=0.0)
		plt.plot(st4['DateConvert2']/24., st4['Concentration'], 'co',markersize=3,markeredgewidth=0.0)
		plt.plot(st5['DateConvert2']/24., st5['Concentration'], 'mo',markersize=3,markeredgewidth=0.0)
		plt.plot(st6['DateConvert2']/24., st6['Concentration'], 'yo',markersize=3,markeredgewidth=0.0)
	if plotoption == 7:
		plt.plot(st1['DateConvert2']/24., st1['Concentration'], 'ro',markersize=3,markeredgewidth=0.0)
		plt.plot(st2['DateConvert2']/24., st2['Concentration'], 'bo',markersize=3,markeredgewidth=0.0)
		plt.plot(st3['DateConvert2']/24., st3['Concentration'], 'go',markersize=3,markeredgewidth=0.0)
		plt.plot(st4['DateConvert2']/24., st4['Concentration'], 'co',markersize=3,markeredgewidth=0.0)
		plt.plot(st5['DateConvert2']/24., st5['Concentration'], 'mo',markersize=3,markeredgewidth=0.0)
		plt.plot(st6['DateConvert2']/24., st6['Concentration'], 'yo',markersize=3,markeredgewidth=0.0)
		plt.plot(st7['DateConvert2']/24., st7['Concentration'], 'ko',markersize=3,markeredgewidth=0.0)
	if plotoption == 11:
		plt.plot(st1['DateConvert2']/24., st1['Concentration'], 'ro',markersize=3,markeredgewidth=0.0)
		plt.plot(st2['DateConvert2']/24., st2['Concentration'], 'bo',markersize=3,markeredgewidth=0.0)
		plt.plot(st3['DateConvert2']/24., st3['Concentration'], 'go',markersize=3,markeredgewidth=0.0)
		plt.plot(st4['DateConvert2']/24., st4['Concentration'], 'co',markersize=3,markeredgewidth=0.0)
		plt.plot(st5['DateConvert2']/24., st5['Concentration'], 'mo',markersize=3,markeredgewidth=0.0)
		plt.plot(st6['DateConvert2']/24., st6['Concentration'], 'yo',markersize=3,markeredgewidth=0.0)
		plt.plot(st7['DateConvert2']/24., st7['Concentration'], 'ko',markersize=3,markeredgewidth=0.0)
		plt.plot(st8['DateConvert2']/24., st8['Concentration'], marker='o',c='0.80',linestyle='None',markersize=3,markeredgewidth=0.0)
		plt.plot(st9['DateConvert2']/24., st9['Concentration'], marker='o',c='0.60',linestyle='None',markersize=3,markeredgewidth=0.0)
		plt.plot(st10['DateConvert2']/24., st10['Concentration'], marker='o',c='0.40',linestyle='None',markersize=3,markeredgewidth=0.0)
		plt.plot(st11['DateConvert2']/24., st11['Concentration'], marker='o',c='0.20',linestyle='None',markersize=3,markeredgewidth=0.0)			
	#plt.plot(ncfile.variables['ho'][:,0,yloc,xloc]*100000000.,linewidth=3)
	#plt.plot(ncfile.variables['iso'][:,0,yloc,xloc]*100-00.,linewidth=3)
	#plt.plot((np.convolve(abs(plotvar-avg2), np.ones((24,))/24, mode='valid')),linewidth=3) # Here we make a moving average
	#plt.plot(timeplotvar,avg2,linewidth=1,color='black',zorder=2)
	plt.legend(legarray,fontsize=7)
	plt.xlim([spinup/24, daysim])
	plt.xticks(np.arange(spinup/24, daysim-(6/24), 3))
	plt.xlabel('Time since start of simulation [days]')
	plt.ylabel(''+str(pollegend)+' mixing ratio [ppb]')
	if figsave == 1:
        	plt.savefig('/home/WUR/barte035/WRFChem/Python/Figures/'+str(city)+''+str(pollutant)+'.png',bbox_inches='tight')
	plt.show()


	'''
	plt.plot(plotvar,linewidth=3)
	plt.plot(st1['DateConvert2'], st1['Concentration'], 'ro',markersize=3,markeredgewidth=0.0)
	plt.plot(st2['DateConvert2'], st2['Concentration'], 'bo',markersize=3,markeredgewidth=0.0)
	plt.plot(st3['DateConvert2'], st3['Concentration'], 'go',markersize=3,markeredgewidth=0.0)
	plt.plot(st4['DateConvert2'], st4['Concentration'], 'co',markersize=3,markeredgewidth=0.0)
	plt.plot(st5['DateConvert2'], st5['Concentration'], 'mo',markersize=3,markeredgewidth=0.0)
	plt.plot(st6['DateConvert2'], st6['Concentration'], 'yo',markersize=3,markeredgewidth=0.0)
	plt.plot(st7['DateConvert2'], st7['Concentration'], 'ko',markersize=3,markeredgewidth=0.0)
	plt.plot(st8['DateConvert2'], st8['Concentration'], marker='o',c='0.66',linestyle='None',markersize=3,markeredgewidth=0.0)
	plt.plot(st9['DateConvert2'], st9['Concentration'], marker='o',c='0.33',linestyle='None',markersize=3,markeredgewidth=0.0)
	plt.plot(st10['DateConvert2'], st10['Concentration'], marker='o',c='0.33',linestyle='None',markersize=3,markeredgewidth=0.0)
	plt.plot(st11['DateConvert2'], st1['Concentration'], marker='o',c='0.33',linestyle='None',markersize=3,markeredgewidth=0.0)	
	#if station in [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,23,25,26,27,28,30]: # If more then one station is presetent, we plot all the stations (even if data is not available).
	#	plt.plot(st2['DateConvert2'], st2['Concentration'], 'bo',markersize=3,markeredgewidth=0.0)
	#	plt.plot(st3['DateConvert2'], st3['Concentration'], 'go',markersize=3,markeredgewidth=0.0)
	#if station in [2,3,4,5,6,7,8,9,10,11,12,14,15,19,21,22,23,25,26,28]:
	#	plt.plot(st4['DateConvert2'], st4['Concentration'], 'co',markersize=3,markeredgewidth=0.0)
	#	plt.plot(st5['DateConvert2'], st5['Concentration'], 'mo',markersize=3,markeredgewidth=0.0)
	#if station in [2,6,9,10,11,12,14,15,19,21,22,23,25,26,28]:
	#	plt.plot(st6['DateConvert2'], st6['Concentration'], 'yo',markersize=3,markeredgewidth=0.0)
	#if station in [14,15,19,21,22,23,25,26,28]:
	#	plt.plot(st7['DateConvert2'], st7['Concentration'], 'ko',markersize=3,markeredgewidth=0.0)
	#	plt.plot(st8['DateConvert2'], st8['Concentration'], marker='o',c='0.66',linestyle='None',markersize=3,markeredgewidth=0.0)
	#	plt.plot(st9['DateConvert2'], st9['Concentration'], marker='o',c='0.33',linestyle='None',markersize=3,markeredgewidth=0.0)
	#plt.plot(ncfile.variables['ho'][:,0,yloc,xloc]*100000000.,linewidth=3)
	#plt.plot(ncfile.variables['iso'][:,0,yloc,xloc]*100-00.,linewidth=3)
	#plt.plot((np.convolve(abs(plotvar-avg2), np.ones((24,))/24, mode='valid')),linewidth=3) # Here we make a moving average
	plt.legend(legarray,fontsize=6)
	plt.xlim([spinup, timesim])
	plt.xticks(np.arange(spinup, timesim, 72))
	plt.xlabel('Time since start of simulation [hours]')
	plt.ylabel(''+str(pollegend)+' concentration [ppb]')
	plt.show()
	'''
	
	#REMOVE TIMES WHERE NO OBSERVATIONAL DATA IS AVAILABLE (GET THE VALUE OF 1)
	plotvarday = np.zeros(timesim)
	plotvarnight = np.zeros(timesim)
	#plotvarstdday = np.zeros(timesim)
	#plotvarstdnight = np.zeros(timesim)
	for i in range(spinup,timesim,1):
		plotvarday[i] = plotvar[i]
		plotvarnight[i] = plotvar[i]
		#plotvarstdday[i] = plotvarstd[i]
		#plotvarstdnight[i] = plotvarstd[i]
	for i in range(spinup,timesim,1):
		if pollutant != 'NOX':
			if avg2[i] == 1:
				avg2[i] = nan
				std2[i] = nan
				plotvar[i] = nan
				#plotvarstd[i] = nan
			if avgday2[i] == 1:
				avgday2[i] = nan
				stdday2[i] = nan
				plotvarday[i] = nan
				#plotvarstdday[i] = nan
			if avgnight2[i] == 1:
				avgnight2[i] = nan
				stdnight2[i] = nan
				plotvarnight[i] = nan
				#plotvarstdnight[i] = nan
		if pollutant == 'NOX':
			if avg2[i] == 0:
				avg2[i] = nan
				std2[i] = nan
				plotvar[i] = nan
				#plotvarstd[i] = nan
			if avgday2[i] == 0:
				avgday2[i] = nan
				stdday2[i] = nan
				plotvarday[i] = nan
				#plotvarstdday[i] = nan
			if avgnight2[i] == 0:
				avgnight2[i] = nan
				stdnight2[i] = nan
				plotvarnight[i] = nan
				#plotvarstdnight[i] = nan
			
	avg2 = avg2[~np.isnan(avg2)]
	std2 = std2[~np.isnan(std2)]
	avgday2 = avgday2[~np.isnan(avgday2)]
	stdday2 = stdday2[~np.isnan(stdday2)]
	avgnight2 = avgnight2[~np.isnan(avgnight2)]
	stdnight2 = stdnight2[~np.isnan(stdnight2)]
	plotvar = plotvar[~np.isnan(plotvar)]
	plotvarday = plotvarday[~np.isnan(plotvarday)]
	plotvarnight = plotvarnight[~np.isnan(plotvarnight)]
	#plotvarstd = plotvarstd[~np.isnan(plotvarstd)]
	#plotvarstdday = plotvarstdday[~np.isnan(plotvarstdday)]
	#plotvarstdnight = plotvarstdnight[~np.isnan(plotvarstdnight)]
	
	#stdwrfmean = np.nanmean(plotvarstd)
	#stdwrfmeanday = np.nanmean(plotvarstdday)
	#stdwrfmeannight = np.nanmean(plotvarstdnight)
	stdobsmean = np.nanmean(std2)
	stdobsmeanday = np.nanmean(stdday2)
	stdobsmeannight = np.nanmean(stdnight2)	
	
	for i in range(0,stdday2.shape[0],1):
	#for i in range(0,plotvarstdday.shape[0],1):
		if i not in [20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740]:
			stdday2[i] = 0
			#plotvarstd[i] = 0
	for i in range(0,stdnight2.shape[0],1):
	#for i in range(0,plotvarstdnight.shape[0],1):
		if i not in [20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740]:
			stdnight2[i] = 0
			#plotvarstd[i] = 0
	
	'''
	plt.scatter(avg2, plotvar, c='r')
	plt.plot([0,100000],[0,100000],'b-')
	plt.xlim([0, max3])
	plt.ylim([0, max3])
	plt.ylabel(''+str(pollegend)+' concentration WRF-Chem [ppb]')
	plt.xlabel(''+str(pollegend)+' concentration observations [ppb]')	
	plt.show()
	'''
	
	#plt.scatter(avgday2, plotvarday, c='y')
	#plt.scatter(avgnight2, plotvarnight, c='b')
	plt.errorbar(plotvarday, avgday2, stdday2, linestyle='None',markersize=5,marker='o',c='y',capsize=0,alpha=0.7)
	plt.errorbar(plotvarnight, avgnight2, stdnight2, linestyle='None',markersize=5,marker='o',c='b',capsize=0,alpha=0.7)
	plt.plot([0,100000],[0,100000],'r-')
	plt.xlim([0, max3])
	plt.ylim([0, max3])
	#plt.xlim([0, 40])
	#plt.ylim([0, 40])
	plt.xlabel(''+str(pollegend)+' mixing ratio WRF-Chem [ppb]')
	plt.ylabel(''+str(pollegend)+' mixing ratio observations [ppb]')
	plt.legend(['1:1 line','Day','Night'],fontsize=7)
	if figsave == 1:
        	plt.savefig('/home/WUR/barte035/WRFChem/Python/Figures/'+str(city)+''+str(pollutant)+'scatter.png',bbox_inches='tight')	
	plt.show()	
	  
	#STATISTICS
	meanobs = np.mean(avg2)
	meanwrf = np.mean(plotvar)
	rmse = sqrt(mean_squared_error(avg2, plotvar))
	nrmse = rmse/np.mean(avg2)
	pearsonrr = pearsonr(avg2, plotvar)
	r2 = (r2_score(avg2,plotvar))
	MBi = np.zeros(avg2.shape[0])
	MBEi = np.zeros(avg2.shape[0])
	for i in range(0,avg2.shape[0],1):
		MBi[i] = (plotvar[i]-avg2[i])
		MBEi[i] = ((plotvar[i]-avg2[i])/avg2[i])
	MB = np.sum(MBi)/avg2.shape[0]
	MBE = np.sum(MBEi)/avg2.shape[0]

	meanobsday = np.mean(avgday2)
	meanwrfday = np.mean(plotvarday)	
	rmseday = sqrt(mean_squared_error(avgday2, plotvarday))
	nrmseday = rmseday/np.mean(avgday2)
	pearsonrrday = pearsonr(avgday2, plotvarday)
	r2day = (r2_score(avgday2,plotvarday))
	MBiday = np.zeros(avgday2.shape[0])
	MBEiday = np.zeros(avgday2.shape[0])
	for i in range(0,avgday2.shape[0],1):
		MBiday[i] = (plotvarday[i]-avgday2[i])
		MBEiday[i] = ((plotvarday[i]-avgday2[i])/avgday2[i])
	MBday = np.sum(MBiday)/avgday2.shape[0]
	MBEday = np.sum(MBEiday)/avgday2.shape[0]

	meanobsnight = np.mean(avgnight2)
	meanwrfnight = np.mean(plotvarnight)	
	rmsenight = sqrt(mean_squared_error(avgnight2, plotvarnight))
	nrmsenight = rmsenight/np.mean(avgnight2)
	pearsonrrnight = pearsonr(avgnight2, plotvarnight)
	r2night = (r2_score(avgnight2,plotvarnight))
	MBinight = np.zeros(avgnight2.shape[0])
	MBEinight = np.zeros(avgnight2.shape[0])
	for i in range(0,avgnight2.shape[0],1):
		MBinight[i] = (plotvarnight[i]-avgnight2[i])
		MBEinight[i] = ((plotvarnight[i]-avgnight2[i])/avgnight2[i])
	MBnight = np.sum(MBinight)/avgnight2.shape[0]
	MBEnight = np.sum(MBEinight)/avgnight2.shape[0]

	print ''
	print 'Mean observations',meanobs
	print 'STD observations',stdobsmean
	print 'Mean WRF-Chem',meanwrf
	#print 'STD WRF-Chem',stdwrfmean
	#print 'MB =',MB
	print 'MBE =',MB
	print 'RMSE =',rmse
	print 'NRMSE =',nrmse
	print 'Pearson R =',pearsonrr[0]
	print '2-tailed P =',pearsonrr[1]
	print 'R-square =',r2

	print ''
	print 'Mean observations (day)',meanobsday
	print 'STD observations (day)',stdobsmeanday
	print 'Mean WRF-Chem (day)',meanwrfday
	#print 'STD WRF-Chem (day)',stdwrfmeanday
	#print 'MB (day) =',MBday
	print 'MBE (day) =',MBday
	print 'RMSE (day) =',rmseday
	print 'NRMSE (day) =',nrmseday
	print 'Pearson R (day) =',pearsonrrday[0]
	print '2-tailed P (day) =',pearsonrrday[1]
	print 'R-square (day) =',r2day

	print ''
	print 'Mean observations (night)',meanobsnight
	print 'STD observations (night)',stdobsmeannight
	print 'Mean WRF-Chem (night)',meanwrfnight
	#print 'STD WRF-Chem (night)',stdwrfmeannight
	#print 'MB (night) =',MBnight
	print 'MBE (night) =',MBnight
	print 'RMSE (night) =',rmsenight
	print 'NRMSE (night) =',nrmsenight
	print 'Pearson R (night) =',pearsonrrnight[0]
	print '2-tailed P (night) =',pearsonrrnight[1]
	print 'R-square (night) =',r2night

	print ''
	print '(1)',meanwrfday
	print '(2)',meanwrfnight
	print '(3)',meanobsday
	print '(4)',meanobsnight
	print '(5)',stdobsmeanday
	print '(6)',stdobsmeannight
	print '(7)',MBday
	print '(8)',MBnight
	print '(9)',nrmseday,'2dec'
	print '(10)',nrmsenight,'2dec'
	print '(11)',pearsonrrday[0],'2dec'
	print '(12)',pearsonrrnight[0],'2dec'
	
	plotvardiurnal = plotvar[spinup:plotvar.shape[0]]
	avg2diurnal = avg2[spinup:avg2.shape[0]]

	daysim2 = timesim/24

	plotvardiurnalmean = np.zeros(24)
	avg2diurnalmean = np.zeros(24)
	plotvardiurnalstd = np.zeros(24)
	avg2diurnalstd = np.zeros(24)
	for i in range(0,24,1):
		globals()["plotvardiurnal" + str(i)] = np.zeros(28) #was 29
		globals()["avg2diurnal" + str(i)] = np.zeros(28) #was 29
		for j in range(0,28,1):	#was 29
			globals()["plotvardiurnal" + str(i)][j] = plotvardiurnal[(j*24)+i]
			globals()["avg2diurnal" + str(i)][j] = avg2diurnal[(j*24)+i]
		plotvardiurnalmean[i] = np.mean(globals()["plotvardiurnal" + str(i)])
		plotvardiurnalstd[i] = np.std(globals()["plotvardiurnal" + str(i)])
		avg2diurnalmean[i] = np.mean(globals()["avg2diurnal" + str(i)])
		avg2diurnalstd[i] = np.std(globals()["avg2diurnal" + str(i)])
		
	temp = [plotvardiurnalmean[0],plotvardiurnalmean[1],plotvardiurnalmean[2],plotvardiurnalmean[3],plotvardiurnalmean[4]]
	plotvardiurnalmean[0:19] = plotvardiurnalmean[5:24]
	plotvardiurnalmean[19:24] = temp
	temp = [plotvardiurnalstd[0],plotvardiurnalstd[1],plotvardiurnalstd[2],plotvardiurnalstd[3],plotvardiurnalstd[4]]
	plotvardiurnalstd[0:19] = plotvardiurnalstd[5:24]
	plotvardiurnalstd[19:24] = temp
	temp = [avg2diurnalmean[0],avg2diurnalmean[1],avg2diurnalmean[2],avg2diurnalmean[3],avg2diurnalmean[4]]
	avg2diurnalmean[0:19] = avg2diurnalmean[5:24]
	avg2diurnalmean[19:24] = temp
	temp = [avg2diurnalstd[0],avg2diurnalstd[1],avg2diurnalstd[2],avg2diurnalstd[3],avg2diurnalstd[4]]
	avg2diurnalstd[0:19] = avg2diurnalstd[5:24]
	avg2diurnalstd[19:24] = temp
				
	plt.plot(np.arange(24),plotvardiurnalmean,color='black',linewidth=2)
	plt.plot(np.arange(24),avg2diurnalmean,color='red',linewidth=2)
	plt.fill_between(np.arange(24),plotvardiurnalmean-plotvardiurnalstd,plotvardiurnalmean+plotvardiurnalstd,alpha=0.2,color='black')
	plt.fill_between(np.arange(24),avg2diurnalmean-avg2diurnalstd,avg2diurnalmean+avg2diurnalstd,alpha=0.2,color='red')	
	plt.axvline(x=6,linestyle='--',color='grey')
	plt.axvline(x=18,linestyle='--',color='grey')
	plt.xlim(0,23)
	#plt.ylim(0,100) #for NO, comment for other vars!
	plt.xticks(np.arange(0,23.1,6))
	plt.axvspan(xmin=0.,xmax=6.,ymin=0.97,ymax=1,facecolor='blue',alpha=0.3)
	plt.axvspan(xmin=6.,xmax=18.,ymin=0.97,ymax=1,facecolor='yellow',alpha=0.3)
	plt.axvspan(xmin=18.,xmax=24.,ymin=0.97,ymax=1,facecolor='blue',alpha=0.3)
	plt.ylabel(''+str(pollegend)+' mixing ratio [ppb]')
	plt.xlabel('Local time [hours]')
	plt.legend(['WRF-Chem','Observations'],fontsize=7)
	if figsave == 1:
        	plt.savefig('/home/WUR/barte035/WRFChem/Python/Figures/'+str(city)+''+str(pollutant)+'diurnal.png',bbox_inches='tight')
	plt.show()

	return(plotvardiurnalmean,avg2diurnalmean)
