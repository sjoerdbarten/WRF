import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import csv
import netCDF4 as nc
from os.path import dirname, join as pjoin
import scipy.io as sio
import pandas as pd
from matplotlib.lines import Line2D

def makemapofstations():
    #DEFINE LON/LAT in array
    barrow = [71.3230,-156.6114]
    storhofdi = [63.400,-20.288]
    summit = [72.5800018311,-38.4799995422]
    ahtari = [62.583333,24.183333]
    bredkalen = [63.85,15.333333]
    esrange = [67.883333,21.066667]
    karasjok = [69.466667,25.216667]
    karvatn = [62.783333,8.883333]
    lerwick = [60.13922,-1.185319]
    oulanka = [66.320278,29.401667]
    pallas = [67.973333333,24.116111111]
    tustervatn = [65.833333,13.916667]
    villum = [81.6,-16.67]
    vindeln = [64.25,19.766667]
    virolahti = [60.526667,27.686111]
    zeppelin = [78.90715,11.88668]    #Ny-Alesund
    ascos = [87.5178,-7.882] #rough estimate [87.4,-6.0]
    hurdal = [60.372386,11.078142]
    whitehorse = [60.718609,-135.049193]
    yellowknife = [62.45207,-114.364]
    normanwells = [65.27926,-126.813]
    fortliard = [60.23583,-123.467]
    inuvik = [68.36005,-133.727]
    denalinp = [63.72,-148.97]
    alert = [82.4991455078,-62.3415260315]
    eureka = [79.99,-85.94]
    sodankyla = [67.3638000488,26.630399704]
    scoresbysund = [70.4848,-21.9512] #ittoqqortoormiit
    #thule = [76.5166702271,-68.7666702271] #pittufik
    nyalesund = [78.923576355,11.9236602783] #this is downhill, zeppelin is uphill at ny-alesund
    resolute = [74.71,-94.97] #resolute airport   
    

    lons = [barrow[1],storhofdi[1],summit[1],ahtari[1],bredkalen[1],esrange[1],karasjok[1],karvatn[1],lerwick[1],oulanka[1],pallas[1],tustervatn[1],villum[1],vindeln[1],virolahti[1],zeppelin[1],ascos[1],hurdal[1],whitehorse[1],yellowknife[1],normanwells[1],fortliard[1],inuvik[1],alert[1],denalinp[1]]
    lats = [barrow[0],storhofdi[0],summit[0],ahtari[0],bredkalen[0],esrange[0],karasjok[0],karvatn[0],lerwick[0],oulanka[0],pallas[0],tustervatn[0],villum[0],vindeln[0],virolahti[0],zeppelin[0],ascos[0],hurdal[0],whitehorse[0],yellowknife[0],normanwells[0],fortliard[0],inuvik[0],alert[0],denalinp[0]]    
    text = ['BRW','ICE','SUM','AHT','BRE','ESR','KAS','KRV','SIS','OUX','PAL','TUV','VIL','VDL','VIR','ZEP','ASC','HUR','WHI','YEL','NOR','FOR','INU','ALT','DEN']
    lons_sondes = [ascos[1],alert[1],eureka[1],sodankyla[1],scoresbysund[1],nyalesund[1],summit[1],resolute[1],lerwick[1]]
    lats_sondes = [ascos[0],alert[0],eureka[0],sodankyla[0],scoresbysund[0],nyalesund[0],summit[0],resolute[0],lerwick[0]]
    text_sondes = ['ASC','ALT','EUK','SDK','SCO','NYA','SUM','RSL','SIS']
    
    #Define basemap
    geopath = '/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/geo_em.d01.nc'
    geo = nc.Dataset(geopath,'r')
    wrflat = geo.variables['XLAT_M'][0,:,:]
    wrflon = geo.variables['XLONG_M'][0,:,:]
    we  = geo.variables['XLAT_M'].shape[2]
    sn  = geo.variables['XLAT_M'].shape[1]
    lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
    lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
    dx_dom = 30000
    wrffile = '/archive/ESG/barte035/MOSAiC/o3_analysis/CorrectFiles/wrfout/wrfout_nudgedBL_fixeddep_d01_2008-08-10_00:00:00'
    ds = nc.Dataset(wrffile,'r')
    
    #shippos_ascos = sio.loadmat("/home/WUR/barte035/WRFChem/o3_analysis_DATA/shippos.mat")
    #shippos = shippos_ascos["shippos"]

    shippos_ascos = pd.read_csv("/home/WUR/barte035/WRFChem/o3_analysis_DATA/OdenWeatherStation_ascii.txt", sep=' ', skiprows=3, header=None, skipinitialspace=True).drop([0])
    shippos_ascos = shippos_ascos[[0,7,8]]
    shippos_ascos = shippos_ascos.loc[(shippos_ascos[0] >= 228.41666666667) & (shippos_ascos[0] <= 247.)]
    lon_path = shippos_ascos[8].values
    lat_path = shippos_ascos[7].values
   
    fig, ax = plt.subplots(figsize=(12,12))
    m = Basemap(projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0,resolution='l')
    lonswrf,latswrf = m(ds.variables['XLONG'][0,:,:],ds.variables['XLAT'][0,:,:],inverse=False)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.1,1.5],colors='w',zorder=2,alpha=0.1)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.2,1.5],colors='w',zorder=2,alpha=0.2)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.3,1.5],colors='w',zorder=2,alpha=0.3)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.4,1.5],colors='w',zorder=2,alpha=0.4)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.5,1.5],colors='w',zorder=2,alpha=0.5)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.6,1.5],colors='w',zorder=2,alpha=0.6)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.7,1.5],colors='w',zorder=2,alpha=0.7)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.8,1.5],colors='w',zorder=2,alpha=0.8)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[0.9,1.5],colors='w',zorder=2,alpha=0.9)
    m.contourf(lonswrf,latswrf,ds.variables['SEAICE'][0,:,:],[1.0,1.5],colors='w',zorder=2,alpha=1.0)   
    m.contourf(lonswrf,latswrf,ds.variables['SNOWC'][0,:,:],[0.5,1.5],colors='w',zorder=2,alpha=0.9)
    x_path, y_path = m(lon_path, lat_path)
    m.plot(x_path,y_path, linewidth=3,color='black',zorder=4)
    m.drawcoastlines(linewidth=0.3)
    m.drawcountries(linewidth=0.3)
    m.fillcontinents(color='peru',lake_color='cornflowerblue')
    m.drawparallels(np.arange(-80.,81.,10.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='cornflowerblue')
    ax = plt.gca()
#    for iii in range(0,len(lons_sondes)):
#        x2,y2 = m(lons_sondes[iii],lats_sondes[iii])
#	if text_sondes[iii] != 'ASC':
#     		plt.plot(x2,y2,marker='s',markersize=12,color='red')
#	if text_sondes[iii] == 'ASC':
#     		plt.plot(x2,y2,marker='s',markersize=12,color='red',zorder=3,label='Radiosondes')	
#	if text_sondes[iii] != 'NYA':
#		plt.text(x2+0.1e6,y2,text_sondes[iii],color='black',weight="bold")
#	if text_sondes[iii] == 'NYA':
#		plt.text(x2+0.1e6,y2,'NYA/ZEP',color='black',weight="bold")
    for ii in range(0,len(lons)):
        x,y = m(lons[ii],lats[ii])
        #if text[ii] != 'ASC':
        #	plt.plot(x,y,marker='o',markersize=10,color='green')
        if text[ii] in ['ASC','ALT','VIL','ZEP','BRW','SUM']:
            plt.plot(x,y,marker='o',markersize=10,color='green')
        if text[ii] in ['DEN','ESR','KAS','INU','SIS','PAL','ICE','YEL']:
            plt.plot(x,y,marker='o',markersize=10,color='magenta')		
        if text[ii] in ['AHT','BRE','FOR','HUR','KRV','NOR','OUX','TUV','VDI', 'VIR','WHI']:
            plt.plot(x,y,marker='o',markersize=10,color='cyan')		
        #if text[ii] == 'ASC':
        #	plt.plot(x,y,marker='o',markersize=10,color='green',zorder=3,label='Surface Ozone')
        if text[ii] not in ['ESR','BRE','OUX','PAL']:
            plt.text(x+0.1e6,y,text[ii],color='black',weight="bold")
        if text[ii] == 'ESR':
            plt.text(x-0.32e6,y-0.05e6,text[ii],color='black',weight="bold")
        if text[ii] == 'PAL':
            plt.text(x-0.32e6,y,text[ii],color='black',weight="bold")
        if text[ii] in ['OUX','BRE']:
            plt.text(x+0.1e6,y-0.05e6,text[ii],color='black',weight="bold")
#    legend_elements = [Line2D([0],[0],marker='o',markersize=10,color='green',label='Surface Ozone',lw=0),
#    		       Line2D([0],[0],marker='s',markersize=12,color='red',label='Radiosondes',lw=0)]
    legend_elements = [Line2D([0],[0],marker='o',markersize=10,color='green',label='High Arctic',lw=0),
    			Line2D([0],[0],marker='o',markersize=10,color='magenta',label='Remote',lw=0),
			Line2D([0],[0],marker='o',markersize=10,color='cyan',label='Terrestrial',lw=0),]
    ax.legend(handles=legend_elements,loc=2,numpoints=1)
    fig.tight_layout()
    plt.savefig('Figures/2008_stationmap_nosondes_hires.png', format='png', dpi=500)
    plt.show()
    
makemapofstations()
