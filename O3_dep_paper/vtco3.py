import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def vtcfluxes(wrffile,savefig):
	wrf = nc.Dataset(wrffile,'r')
	print(wrf)
	arrvar = ['VTC','VTCO','VTCDMS','VTCODMS','VTC2DMS','VTCCO2','VTCOCO2','VTC2CO2','VTCCH4']
	#arrvar = ['VTC','VTCO','VTCDMS','VTCODMS','VTC2DMS','VTCCO2','VTCOCO2','VTC2CO2','VTCCH4','DEP_VEL']
	legendlist = ['b-','b:','g-','g:','g--','y-','y:','y--','c-','r-']
	plt.figure(figsize=(15,6))	
	#all to cm/s (from m/s)
	for i in range(0,len(arrvar)):
		plotvar = wrf[arrvar[i]][1:,35,121]*100
		plt.plot(plotvar,legendlist[i],label=arrvar[i])
	plt.legend(loc='best')
	plt.xlabel('Time [days]')
	plt.ylabel('Piston velocity [cm s$^{-1}$]')
	plt.xlim(0,len(plotvar))
	plt.xticks(np.arange(0.,len(plotvar),24.),np.arange(0.,len(plotvar)/24.,1.))
	if savefig == 1:
		plt.savefig('Figures/pistonvelocitytimeseries',dpi=300)
	plt.show()
	
	arrvar2 = ['VTC','VTCDMS','VTCCO2','VTCCH4']
	#arrvar2 = ['VTC','VTCDMS','VTCCO2','VTCCH4','DEP_VEL']
	legendlist2 = ['b','g','y','c','r']
	polyfit = 2
	plt.figure(figsize=(15,6))
	#wind = wrf['UST'][:,35,121]
	wind = (wrf['U_PHY'][1:,0,35,121]**2+wrf['V_PHY'][1:,0,35,121]**2)**0.5
	for i in range(0,len(arrvar2)):
		plotvar = wrf[arrvar2[i]][1:,35,121]*100
		plt.scatter(wind,plotvar,c=legendlist2[i],label=arrvar2[i])
		plt.plot(np.unique(wind),np.poly1d(np.polyfit(wind,plotvar,polyfit))(np.unique(wind)),c=legendlist2[i])
		print(arrvar2[i]+' polyfitvar: '+str(np.polyfit(wind,plotvar,polyfit)))
	plt.legend(loc='best')
	plt.xlabel('Wind speed [m s$^{-1}$]')
	plt.ylabel('Piston velocity [cm s$^{-1}$]')
	plt.xlim(0,max(wind)*1.01)
	plt.ylim(0,max(wrf['VTCCO2'][1:,35,121])*100*1.01)
	#plt.ylim(0,max(wrf['DEP_VEL'][1:,35,121])*100*1.01)
	if savefig == 1:
		plt.savefig('Figures/correlationwindpistonvelocity',dpi=300)
	plt.show()
	
	plt.figure(figsize=(15,6))
	temp = (wrf['T_PHY'][1:,0,35,121])
	for i in range(0,len(arrvar2)):
		plotvar = wrf[arrvar2[i]][1:,35,121]*100
		plt.scatter(temp,plotvar,c=legendlist2[i],label=arrvar2[i])
		plt.plot(np.unique(temp),np.poly1d(np.polyfit(temp,plotvar,polyfit))(np.unique(temp)),c=legendlist2[i])
		print(arrvar2[i]+' polyfitvar: '+str(np.polyfit(temp,plotvar,polyfit)))
	plt.legend(loc='best')
	plt.xlabel('Sea surface temperature [K]')
	plt.ylabel('Piston velocity [cm s$^{-1}$]')
	print(min(temp))
	print(temp)
	plt.xlim(min(temp[np.nonzero(temp)])-0.5,max(temp)+0.5)
	plt.ylim(0,max(wrf['VTCCO2'][1:,35,121])*100*1.01)
	#plt.ylim(0,max(wrf['DEP_VEL'][1:,35,121])*100*1.01)
	if savefig == 1:
		plt.savefig('Figures/correlationwindpistonvelocitytemp',dpi=300)
	plt.show()
	
	a04co2 = 2073.1
	a05co2 = 125.62
	a06co2 = 3.6276
	a07co2 = 0.043219
	a01dms = 2674.0
	a02dms = -147.1
	a03dms = 3.726
	a04dms = -0.038
	
	scwcco2 = np.zeros(len((wrf['T_PHY'][1:,0,35,121])))
	scwcdms = np.zeros(len((wrf['T_PHY'][1:,0,35,121])))
	sc_correctionco2 = np.zeros(len((wrf['T_PHY'][1:,0,35,121])))
	sc_correctiondms = np.zeros(len((wrf['T_PHY'][1:,0,35,121])))
	print(scwcco2)
	
	scwcco2 = a04co2-a05co2*(wrf['T_PHY'][1:,0,35,121]-273.15)+a06co2*(wrf['T_PHY'][1:,0,35,121]-273.15)**2-a07co2*(wrf['T_PHY'][1:,0,35,121]-273.15)**3 #!schmidt number
	scwcdms = a01dms+a02dms*(wrf['T_PHY'][1:,0,35,121]-273.15)+a03dms*(wrf['T_PHY'][1:,0,35,121]-273.15)**2+a04dms*(wrf['T_PHY'][1:,0,35,121]-273.15)**3 #!schmidt #
	print(scwcco2)

	sc_correctionco2 = (scwcco2/660.)**0.5
	sc_correctiondms = (scwcdms/660.)**0.5
	
	plt.figure(figsize=(15,6))
	for i in range(0,len(arrvar2)):
		plotvar = wrf[arrvar2[i]][1:,35,121]*100
		if arrvar2[i] == 'VTCCO2':
			plotvar = plotvar*sc_correctionco2
		if arrvar2[i] == 'VTCDMS':
			plotvar = plotvar*sc_correctiondms
		plt.scatter(wind,plotvar,c=legendlist2[i],label=arrvar2[i])
		plt.plot(np.unique(wind),np.poly1d(np.polyfit(wind,plotvar,polyfit))(np.unique(wind)),c=legendlist2[i])
		print(arrvar2[i]+' polyfitvar: '+str(np.polyfit(wind,plotvar,polyfit)))
	plt.legend(loc='best')
	plt.xlabel('Wind speed [m s$^{-1}$]')
	plt.ylabel('Piston velocity normalized for Schmidt Number 660 [cm s$^{-1}$]')
	plt.xlim(0,max(wind)*1.01)
	#plt.xlim(0,max((wrf['U'][:,0,35,121]**2*wrf['V'][:,0,35,121]**2)**0.5)*1.01)
	plt.ylim(0,max(wrf['VTCCO2'][1:,35,121])*100*1.01*max(sc_correctionco2))
	if savefig == 1:
		plt.savefig('Figures/correlationwindpistonvelocitysccorrection',dpi=300)
	plt.show()


	
