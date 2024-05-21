#""" plot_orbit_EUR.py
#
#takes:
#- an 1D [norbit x nscanl x ngr_pix] array with data to be plotted
#- two 2D [norbit x nscanl x ngr_pix,4] array with corner lat/lon informations for the pixels
#- a minimum and maximum value for the color scale of the plot
#- a figure title and label

#and plots:
#- the orbits on the map

#"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from matplotlib import cm
import os

#select dataset
choose_dataset = 11
path = "/home/WUR/barte035/WRFChem/OMI/DATA/"

if choose_dataset == 1: #THIS ONE HAS VALUES (HORIZONTAL???)
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T002200_o50339_fitB_v1.nc")
if choose_dataset == 2:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T020100_o50340_fitB_v1.nc")
if choose_dataset == 3:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T034000_o50341_fitB_v1.nc")
if choose_dataset == 4:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T051900_o50342_fitB_v1.nc")		
if choose_dataset == 5:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T065800_o50343_fitB_v1.nc")
if choose_dataset == 6:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T083600_o50344_fitB_v1.nc")
if choose_dataset == 7:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T101500_o50345_fitB_v1.nc")
if choose_dataset == 8:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T115400_o50346_fitB_v1.nc")			
if choose_dataset == 9:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T133300_o50347_fitB_v1.nc")
if choose_dataset == 10:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T151200_o50348_fitB_v1.nc")
if choose_dataset == 11: #THIS ONE HAS VALUES
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T165100_o50349_fitB_v1.nc")
if choose_dataset == 12: #THIS ONE HAS VALUES
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T183000_o50350_fitB_v1.nc")
if choose_dataset == 13:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T200900_o50351_fitB_v1.nc")
if choose_dataset == 14:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T214700_o50352_fitB_v1.nc")
if choose_dataset == 15:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140101T232600_o50353_fitB_v1.nc")
if choose_dataset == 16: #THIS ONE HAS VALUES (HORIZONTAL???)
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T010500_o50354_fitB_v1.nc")
if choose_dataset == 17:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T024400_o50355_fitB_v1.nc")
if choose_dataset == 18:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T042300_o50356_fitB_v1.nc")
if choose_dataset == 19:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T060200_o50357_fitB_v1.nc")
if choose_dataset == 20:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T074100_o50358_fitB_v1.nc")
if choose_dataset == 21:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T092000_o50359_fitB_v1.nc")
if choose_dataset == 22:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T105900_o50360_fitB_v1.nc")
if choose_dataset == 23:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T123700_o50361_fitB_v1.nc")
if choose_dataset == 24:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T141600_o50362_fitB_v1.nc")
if choose_dataset == 25: #THIS ONE HAS VALUES
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T155500_o50363_fitB_v1.nc")
if choose_dataset == 26: #THIS ONE HAS VALUES
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T173400_o50364_fitB_v1.nc")
if choose_dataset == 27:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T191300_o50365_fitB_v1.nc")
if choose_dataset == 28:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T205200_o50366_fitB_v1.nc")
if choose_dataset == 29:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140102T223100_o50367_fitB_v1.nc")
if choose_dataset == 30: #THIS ONE HAS VALUES (HORIZONTAL???)
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T001000_o50368_fitB_v1.nc")
if choose_dataset == 31:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T014900_o50369_fitB_v1.nc")
if choose_dataset == 32:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T032700_o50370_fitB_v1.nc")
if choose_dataset == 33:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T050600_o50371_fitB_v1.nc")
if choose_dataset == 34:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T064500_o50372_fitB_v1.nc")
if choose_dataset == 35:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T082400_o50373_fitB_v1.nc")
if choose_dataset == 36:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T100300_o50374_fitB_v1.nc")
if choose_dataset == 37:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T114200_o50375_fitB_v1.nc")
if choose_dataset == 38:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T132100_o50376_fitB_v1.nc")
if choose_dataset == 39:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T150000_o50377_fitB_v1.nc")
if choose_dataset == 40: #THIS ONE HAS VALUES
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T163900_o50378_fitB_v1.nc")
if choose_dataset == 41: #THIS ONE HAS VALUES
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T181700_o50379_fitB_v1.nc")
if choose_dataset == 42:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T195600_o50380_fitB_v1.nc")
if choose_dataset == 43:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T213500_o50381_fitB_v1.nc")
if choose_dataset == 44:
	OMI_file = os.path.join(path,"QA4ECV_L2_NO2_OMI_20140103T231400_o50382_fitB_v1.nc")							
	
#Access file and variables (2D arrays: [nscanl; ngrpix])
ds       = nc.Dataset(OMI_file,'r')
prod     = ds.groups['PRODUCT']
geo      = ds.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['GEOLOCATIONS']
indata = ds.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['INPUT_DATA']
detres = ds.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS']
errflag = prod.variables['processing_error_flag'][0,:,:]
amf_trop = prod.variables['amf_trop'][0,:,:]
amf_geo = detres.variables['amf_geo'][0,:,:]
vcd_trop = prod.variables['tropospheric_no2_vertical_column'][0,:,:]
clats    = geo.variables['latitude_bounds'][0,:,:,:]
clons    = geo.variables['longitude_bounds'][0,:,:,:]
zenithflag = geo.variables['solar_zenith_angle'][0,:,:]
snowiceflag = indata.variables['snow_ice_flag'][0,:,:]
amfflag = amf_trop/amf_geo
cloudflag = detres.variables['cloud_radiance_fraction_no2'][0,:,:]

#Select valid pixels (converted to 1D [nscanl x ngrpix] array)
#NO FILTERS
#amf_trop_sel = amf_trop[np.where((errflag == 0))]
#vcd_trop_sel = vcd_trop[np.where((errflag == 0))] #parr_rs
#clats_sel    = clats[np.where((errflag == 0))] #clats_rs
#clons_sel    = clons[np.where((errflag == 0))] #clons_rs
#NO AMF FLAG
#amf_trop_sel = amf_trop[np.where((errflag == 0) & (zenithflag < 80) & (cloudflag <= 0.5) & ((snowiceflag < 10) | (snowiceflag == 255)))]
#vcd_trop_sel = vcd_trop[np.where((errflag == 0) & (zenithflag < 80) & (cloudflag <= 0.5) & ((snowiceflag < 10) | (snowiceflag == 255)))]
#clats_sel = clats[np.where((errflag == 0) & (zenithflag < 80) & (cloudflag <= 0.5) & ((snowiceflag < 10) | (snowiceflag == 255)))]
#clons_sel = clons[np.where((errflag == 0) & (zenithflag < 80) & (cloudflag <= 0.5) & ((snowiceflag < 10) | (snowiceflag == 255)))]
#ALL FILTERS
cloudflag_sel = cloudflag[np.where((errflag == 0) & (zenithflag < 80) & (amfflag > 0.2) & (cloudflag <= 1) & ((snowiceflag < 10) | (snowiceflag == 255)))]
amf_trop_sel = amf_trop[np.where((errflag == 0) & (zenithflag < 80) & (amfflag > 0.2) & (cloudflag <= 0.2) & ((snowiceflag < 10) | (snowiceflag == 255)))]
vcd_trop_sel = vcd_trop[np.where((errflag == 0) & (zenithflag < 80) & (amfflag > 0.2) & (cloudflag <= 0.2) & ((snowiceflag < 10) | (snowiceflag == 255)))]
clats_sel = clats[np.where((errflag == 0) & (zenithflag < 80) & (amfflag > 0.2) & (cloudflag <= 1) & ((snowiceflag < 10) | (snowiceflag == 255)))]
clons_sel = clons[np.where((errflag == 0) & (zenithflag < 80) & (amfflag > 0.2) & (cloudflag <= 1) & ((snowiceflag < 10) | (snowiceflag == 255)))]

def plot_update_range(bmap):
    #lot = [(-10,40),(-10,60),(25,60),(25,40)]
    lot = [(-83,-5),(-83,16),(-60,16),(-60,-5)]
    #lot = [(-82.36575317382812,-4.975730895996094),(-82.36575317382812,15.24386596679),(-60.02581787109,15.24386596679),(-60.02581787109,-4.975730895996094)]
    xs = []
    ys = []
    for i in range(0,len(lot)-1):
        xs.append(lot[i][0])
        xs.append(lot[i+1][0])
        ys.append(lot[i][1])
        ys.append(lot[i+1][1])

    xs.append(lot[0][0])
    xs.append(lot[-1][0])
    ys.append(lot[0][1])
    ys.append(lot[-1][1])
    bmap.plot(xs,ys,latlon=True,color='k',zorder=5)

#def plot_orbit_EUR(orbarray,corner_lats,corner_lons,r_min,r_max,figtitle,figlabel):
def plot_orbit_Colombia(parr_rs,clats_rs,clons_rs,r_min,r_max,figtitle,figlabel):
    #parr = orbarray[:,5:26]
    #clats = corner_lats[:,5:26]
    #clons = corner_lons[:,5:26]
    
    #parr_rs = parr.reshape(parr.shape[0]*parr.shape[1]*parr.shape[2])
    #clats_rs = clats.reshape(clats.shape[0]*clats.shape[1]*clats.shape[2],3)
    #clons_rs = clons.reshape(clons.shape[0]*clons.shape[1]*clons.shape[2],3)
    
    fig,ax = plt.subplots(figsize=(16,6))
    m = Basemap(projection='cyl',llcrnrlat=-4.975730895996094,urcrnrlat=15.24386596679,llcrnrlon=-82.36575317382812,urcrnrlon=-60.02581787109,resolution='l')
    m.drawcoastlines(zorder=3)
    m.drawcountries(zorder=3)
    m.drawmapboundary(fill_color='lightgray')
    parallels = np.arange(-5.,15.,5.)
    meridians = np.arange(-90.,-60,5.)
    m.drawparallels(parallels,labels=[False,True,True,False])
    m.drawmeridians(meridians,labels=[True,False,False,True])
    plot_update_range(m)
    m.fillcontinents(color='lightgray',lake_color='lightgray') 
    
    x,y = m(clons_rs,clats_rs)
    pols = zip(x,y)
    pols = np.swapaxes(pols,1,2)
    
    ###Select plotting variable###
    #NO2 columns
    #coll = PolyCollection(pols,array=parr_rs,cmap=cm.jet,norm=LogNorm(),edgecolors='none',zorder=2,alpha=0.8)
    coll = PolyCollection(pols,array=parr_rs,cmap=cm.jet,edgecolors='none',zorder=2,alpha=0.8)
    #difference plots
    #coll = PolyCollection(pols,array=parr_rs,cmap=cm.bwr,edgecolors='none',zorder=2,alpha=0.8)
    #Air Mass Factors
    #coll = PolyCollection(pols,array=parr_rs,cmap=cm.jet,edgecolors='none',zorder=2,alpha=0.8)
    coll.set_clim([r_min,r_max])
    
    ax.add_collection(coll)
        
    c = fig.colorbar(coll,ax=ax)
    c.set_label(figlabel,size=16)
    #NO2 columns
    #c.set_ticks([1,2,3,4,6,8,11,15,20])
    #c.set_ticklabels([1,2,3,4,6,8,11,15,20])
    #Air Mass Factors
    c.set_ticks([0.,0.3,0.6,0.9,1.2,1.5,1.8])
    c.set_ticklabels([0.,0.3,0.6,0.9,1.2,1.5,1.8])
    #AMF differences & VCD_rec/VCD_QA differences
    #c.set_ticks([-1.,-0.5,0.,0.5,1.])
    #c.set_ticklabels([-1.,-0.5,0.,0.5,1.])
    #c.set_ticks([-0.5,0.,0.5])
    #c.set_ticklabels([-0.5,0.,0.5])
    #AMF ratios
    #c.set_ticks([0.5,0.75,1.,1.25,1.5])
    #c.set_ticklabels([0.5,0.75,1.,1.25,1.5])
    #NO2 VCD differences
    #c.set_ticks([-7,-5,-3,-1,1,3,5,7])
    #c.set_ticklabels([-7,-5,-3,-1,1,3,5,7])
    #c.set_ticks([-5,-3,-1,1,3,5])
    #c.set_ticklabels([-5,-3,-1,1,3,5])
    
    plt.title(figtitle,size=18,weight='bold')
    plt.show()


#plot_update_range(Basemap)
#plot_orbit_EUR(orbarray=vcd_trop_sel,corner_lats=clats_sel,corner_lons=clons_sel,r_min=0,r_max=1,figtitle='test',figlabel='NO2')

    
    
    
