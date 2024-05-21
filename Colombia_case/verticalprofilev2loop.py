import numpy as np
from matplotlib import pyplot as pl


print 'This script loops over verticalprofilev2 and plots average profiles of any variable over 3 locations. verticalprofilev2loop(plotvar,area,figsave2).'

plotvar = 'nox'

execfile("verticalprofilev2.py")
pressurebogota,concentrationbogota,tempbogota = verticalprofilev2(plotvar,'Bogota',1,0)
pressureorinoco,concentrationorinoco,temporinoco = verticalprofilev2(plotvar,'Orinoco',1,0)
pressureamazon,concentrationamazon,tempamazon = verticalprofilev2(plotvar,'Amazone',1,0)
	
pressurebogotamean = np.mean(pressurebogota,axis=0)
concentrationbogotamean = np.mean(concentrationbogota,axis=0)
errorbogota = np.std(concentrationbogota,axis=0)
tempbogotamean = np.mean(tempbogota,axis=0)
pressureorinocomean = np.mean(pressureorinoco,axis=0)
concentrationorinocomean = np.mean(concentrationorinoco,axis=0)
errororinoco = np.std(concentrationorinoco,axis=0)
temporinocomean = np.mean(temporinoco,axis=0)
pressureamazonmean = np.mean(pressureamazon,axis=0)
concentrationamazonmean = np.mean(concentrationamazon,axis=0)
erroramazon = np.std(concentrationamazon,axis=0)
tempamazonmean = np.mean(tempamazon,axis=0)

heightbogotamean = np.zeros(pressurebogotamean.shape[0])
heightorinocomean = np.zeros(pressurebogotamean.shape[0])
heightamazonmean = np.zeros(pressurebogotamean.shape[0])

#for i in range(0,pressurebogotamean.shape[0]):
#	heightbogotamean[i] = -np.log(pressurebogotamean[i]/pressurebogotamean[0])*(8.3143*tempbogotamean[i])/(0.02896*9.81)
#	heightorinocomean[i] = -np.log(pressureorinocomean[i]/pressureorinocomean[0])*(8.3143*temporinocomean[i])/(0.02896*9.81)
#	heightamazonmean[i] = -np.log(pressureamazonmean[i]/pressureamazonmean[0])*(8.3143*tempamazonmean[i])/(0.02896*9.81)

#print max(concentrationbogotamean),max(pressurebogotamean),errorbogota[0]
print pressurebogotamean[49],pressureorinocomean[49],pressureamazonmean[49]

plt.plot(concentrationbogotamean,pressurebogotamean,color='red',label='Bogota',zorder=10,linewidth=2)
plt.plot(concentrationorinocomean,pressureorinocomean,color='green',label='Orinoco',zorder=10,linewidth=2)
plt.plot(concentrationamazonmean,pressureamazonmean,color='blue',label='Amazon',zorder=10,linewidth=2)
#plt.plot(concentrationbogotamean,heightbogotamean,color='red',label='Bogota',zorder=10,linewidth=2)
#plt.plot(concentrationorinocomean,heightorinocomean,color='green',label='Orinoco',zorder=10,linewidth=2)
#plt.plot(concentrationamazonmean,heightamazonmean,color='blue',label='Amazon',zorder=10,linewidth=2)
plt.fill_betweenx(pressurebogotamean,concentrationbogotamean-errorbogota,concentrationbogotamean+errorbogota,alpha=0.2,color='red')
plt.fill_betweenx(pressureorinocomean,concentrationorinocomean-errororinoco,concentrationorinocomean+errororinoco,alpha=0.2,color='green')
plt.fill_betweenx(pressureamazonmean,concentrationamazonmean-erroramazon,concentrationamazonmean+erroramazon,alpha=0.2,color='blue')
#plt.fill_betweenx(heightbogotamean,concentrationbogotamean-errorbogota,concentrationbogotamean+errorbogota,alpha=0.2,color='red')
#plt.fill_betweenx(heightorinocomean,concentrationorinocomean-errororinoco,concentrationorinocomean+errororinoco,alpha=0.2,color='green')
#plt.fill_betweenx(heightamazonmean,concentrationamazonmean-erroramazon,concentrationamazonmean+erroramazon,alpha=0.2,color='blue')
ax = plt.axes()
plt.gca().invert_yaxis()
plt.xlabel('NO$_x$ mixing ratio [ppb]')
plt.ylabel('Pressure [hPa]')
#plt.ylabel('Height above surface [m]')
plt.xlim([0, 11])
plt.yscale('log')
plt.ylim([1000,50])
#plt.ylim([0,15000])
plt.yticks([1000,900,800,700,600,500,400,300,200,100,50])
ax.set_yticks([1000,900,800,700,600,500,400,300,200,100,50])
ax.set_yticklabels(["1000","900","800","700","600","500","400","300","200","100","50"])
plt.legend(['Bogota','Orinoco','Amazon'],fontsize=9,loc=1)
plt.text(1,60,'(d)',fontsize=9)
plt.show()
