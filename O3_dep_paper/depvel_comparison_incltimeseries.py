from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib import pyplot
from sklearn.metrics import mean_squared_error
from math import sqrt
from scipy.stats import gaussian_kde
from sklearn.metrics import r2_score
import matplotlib.colors as colors
import pandas as pd
import datetime as dt

wrflat= np.load('/home/WUR/barte035/WRFChem/Python/Arrays/longitudes.npy')
wrflon=	np.load('/home/WUR/barte035/WRFChem/Python/Arrays/latitudes.npy')
wrfseaice= np.load('/home/WUR/barte035/WRFChem/Python/Arrays/seaice.npy')
depvel_run1= np.load('/home/WUR/barte035/WRFChem/Python/Arrays/depvel_run2_new.npy')		#CHANGED THIS FROM DEPVEL_RUN1 (DEFAULT) TO DEPVEL_RUN2 (NUDGED)
depvel_run3= np.load('/home/WUR/barte035/WRFChem/Python/Arrays/depvel_run3_new.npy')
depvel_fixeddep= np.load('/home/WUR/barte035/WRFChem/Python/Arrays/depvel_fixeddep_new.npy')
timearr = np.load('/home/WUR/barte035/WRFChem/Python/Arrays/Surfaceo3/timearr.npy')

wrfdata_fullfile = nc.Dataset('/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_chem_nudgedBL_nofixeddep_d01_2008-08-10_00:00:00','r')
#LU_INDEX 16 = ocean, 24 = ice
lu_index = wrfdata_fullfile.variables['LU_INDEX'][:,:,:]

#depvel_run1 = np.where(lu_index == 16, np.where(depvel_run1 > 0.0006, depvel_run1-0.0008, depvel_run1), depvel_run1)
depvel_run1 = np.where(lu_index == 16, depvel_fixeddep, depvel_run1)
#depvel_run1 = np.where(lu_index == 24, np.where(depvel_run1 > 0.001, depvel_run1-0.001, depvel_run1), depvel_run1)
depvel_run1 = np.where(lu_index == 24, np.where(depvel_run1 > 0.0003, depvel_fixeddep, depvel_run1), depvel_run1)
meandepvel_run1=np.mean(depvel_run1,axis=0)*100.
meandepvel_run3=np.mean(depvel_run3,axis=0)*100.
meandepvel_fixeddep=np.mean(depvel_fixeddep,axis=0)*100.

#To mask deposition to land
meandepvel_run1[meandepvel_run1 > 0.05] = np.nan
meandepvel_run3[meandepvel_run3 > 0.05] = np.nan
meandepvel_fixeddep[meandepvel_fixeddep > 0.05] = np.nan
meandepvel_run1[meandepvel_run1 > 0.048] = 0.048

depvel_run1_copy = depvel_run1
depvel_run1_copy[depvel_run1_copy < 0.00042] = 0.00042

markerlat = 57.7831
markerlon = -32.0054

markerlat2 = 57.6524
markerlon2 = 20.4429

markerlat3 = 69.3729
markerlon3 = -169.563

print(wrfseaice)
print(wrfseaice.shape)

#Plot mean deposition velocities and total flux
fig,axes = plt.subplots(1,3,figsize=(20,5.6), gridspec_kw = {'wspace':0.34, 'hspace':0.005})

m1 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[0])
m1.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
lons,lats = m1(wrflat,wrflon,inverse=False)
m1.contour(lons,lats,wrfseaice,[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
m1.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m1.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
m1.scatter(markerlon,markerlat,s=95,marker='^',zorder=1003,color='red',latlon=True,linewidth=0)
#m1.scatter(markerlon2,markerlat2,s=95,marker='*',zorder=1003,color='red',latlon=True,linewidth=0)
#m1.scatter(markerlon3,markerlat3,s=95,marker='h',zorder=1003,color='red',latlon=True,linewidth=0)
im = m1.imshow(meandepvel_run1,cmap='winter',interpolation='none',zorder=1000, vmin=0.030,vmax=0.050)
#im = m1.imshow(meandepvel_fixeddep,cmap='viridis',interpolation='none',zorder=1000, vmin=0.045,vmax=0.055)

m2 = Basemap(projection='npstere',boundinglat=57.3,lon_0=0.,resolution='l',ax=axes[1])
m2.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m2.contour(lons,lats,wrfseaice,[0.5,1.5],colors='w',zorder=1002,linewidths=[2,2])
m2.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m2.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,False],zorder=1001)
m2.scatter(markerlon,markerlat,s=95,marker='^',zorder=1003,color='green',latlon=True,linewidth=0)
#m2.scatter(markerlon2,markerlat2,s=95,marker='*',zorder=1003,color='green',latlon=True,linewidth=0)
#m2.scatter(markerlon3,markerlat3,s=95,marker='h',zorder=1003,color='green',latlon=True,linewidth=0)
im2 = m2.imshow(meandepvel_run3,cmap='Wistia',interpolation='none',zorder=1000, vmin=0.008,vmax=0.018)

cax = fig.add_axes([axes[0].get_position().x1+0.008,axes[0].get_position().y0+0.038,0.012,axes[0].get_position().y1-axes[0].get_position().y0-0.08])
cbar = fig.colorbar(im,cax=cax,extend='both',ticks=[0.030,0.035,0.040,0.045,0.050])
cbar.ax.tick_params(labelsize=14,width=2)
cbar.set_label('Mean deposition velocity [cm s$^{-1}$]',fontsize=14)
cax2 = fig.add_axes([axes[1].get_position().x1+0.008,axes[1].get_position().y0+0.038,0.012,axes[1].get_position().y1-axes[1].get_position().y0-0.08])
cbar2 = fig.colorbar(im2,cax=cax2,extend='both',ticks=[0.008,0.013,0.018])
cbar2.ax.tick_params(labelsize=14,width=2)
cbar2.set_label('Mean deposition velocity [cm s$^{-1}$]',fontsize=14)

ax=axes[2]
print(min(depvel_run3[24:,20,59]*100.),max(depvel_run3[24:,20,59]*100.),np.mean(depvel_run3[24:,20,59]*100.))
ax.plot(timearr[24:],depvel_run1[24:,20,59]*100.,color='red',zorder=3,label='NUDGED',linewidth=2)
ax.plot(timearr[24:],depvel_run3[24:,20,59]*100.,color='green',zorder=3,label='COAREG',linewidth=2)
#ax.plot(timearr[24:],depvel_run1[24:,11,173]*100.,color='red',zorder=3,linewidth=2,linestyle='dashed')
#ax.plot(timearr[24:],depvel_run3[24:,11,173]*100.,color='green',zorder=3,linewidth=2,linestyle='dashed')
#ax.plot(timearr[24:],depvel_run1[24:,200,110]*100.,color='red',zorder=3,linewidth=2,linestyle='dotted')
#ax.plot(timearr[24:],depvel_run3[24:,200,110]*100.,color='green',zorder=3,linewidth=2,linestyle='dotted')
ax.grid(zorder=1)
ax.set_xlim([timearr[24],timearr[-1]])
ax.set_ylim([0,0.06])
ax.tick_params(axis="both", width=1, length=10)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.set_yticks([0,0.01,0.02,0.03,0.04,0.05,0.06])
ax.set_yticklabels([0,0.01,0.02,0.03,0.04,0.05,0.06], fontsize=12)
ax.set_xticklabels(['11-Aug','14-Aug','17-Aug','20-Aug','23-Aug','26-Aug','29-Aug','01-Sep','04-Sep','07-Sep'],fontsize=12,rotation=45,ha='right')
ax.legend(loc=1,fontsize=13)
ax.set_ylabel('Deposition velocity [cm s$^{-1}$]',fontsize=14)

#CHANGE COLORBARS AND ADD TWO MORE LOCATIONS TO (C)

xtext,ytext = m1(-133.,48.)
axes[0].set_title("NUDGED",fontsize=18)
axes[1].set_title("COAREG",fontsize=18)
axes[0].text(xtext,ytext,s='(a)',fontsize=14,zorder=1010)
axes[1].text(xtext,ytext,s='(b)',fontsize=14,zorder=1010)
axes[2].text(timearr[34],0.056,s='(c)',fontsize=14,zorder=1010)
plt.savefig('Figures/DEPVELplots/DEPVELMEANANDFLUX_new',dpi=600)
plt.show()
