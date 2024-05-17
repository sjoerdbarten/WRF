# This script loads the LAND, SEA and NOx variables and makes a timeseries figure.

# 1. Input
# 1.1 Load packages
import os
import csv
import matplotlib.dates as mdates
import wrf
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.transforms import Bbox as mtBbox
np.bool = np.bool_

# 1.2 Set paths
home       = os.path.expanduser('~')     # windows proof
datapath   = '%s/atmmod/WRF/run/'%home

# 1.3 User Settings -- This is where you need to modify things
domain = 'd01'
(yy,mm,dd,HH) = (2018,6,5,0)
datapath_obs = '%s/atmmod/ConcentrationData/'%home


# 1.4 Load observations from file
# --- Observations are from 2016-05-09 00:00h to 2006-05-12 00:00h
stations   = ('Vredepeel','Zegveld','Wekerom')
NOx_obs    = dict()                                                            # initialise an empty dictionary to store the observation data in
for station in stations:                                                       # loop over 3 stations
    obsfile = '%s/%s_NOX.csv'%(datapath_obs,station)       # assign filename for that station
    NOx_obs[station] = list()                                                  # create an empty list in the dictionary to store that station's data in
    with open(obsfile,newline='') as csvfile:                                  # open the file
        filereader = csv.reader(csvfile,delimiter=' ')                         # read the file
        for row in filereader:                                                 # read line by line
           NOx_obs[station].append(float(row[0]))                              # convert this line to a float and append it in the (initially empty) list
    NOx_obs[station] = np.asarray(NOx_obs[station])                                    # convert from list to numpy array

# 1.5 Load variables from wrfoutfile
infilename = '%s/wrfout_%s_%4d-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)
infile     = nc.Dataset(infilename,'r')
name_project = infile.getncattr('TITLE')
xlat       = wrf.getvar(infile,'XLAT')
xlon       = wrf.getvar(infile,'XLONG')

#cart_proj = get_cartopy(rh2

# --- Define station coordinates and find indices into xlat/xlon ---
(lat_1, lon_1)  = (51.5405  , 5.8531)   # Vredepeel 
(lat_2, lon_2)  = (52.134596, 4.808560) # Zegveld-De Meije
(lat_3, lon_3)  = (52.104891, 5.723901) # Wekerom
(ix_1,iy_1)     = wrf.to_np(wrf.ll_to_xy(infile,lat_1,lon_1))
(ix_2,iy_2)     = wrf.to_np(wrf.ll_to_xy(infile,lat_2,lon_2))
(ix_3,iy_3)     = wrf.to_np(wrf.ll_to_xy(infile,lat_3,lon_3))
print('Vredepeel:        (iy,ix) = (%2d,%2d)'%(iy_1,ix_1))
print('Zegveld de Meije: (iy,ix) = (%2d,%2d)'%(iy_2,ix_2))
print('Wekerom:          (iy,ix) = (%2d,%2d)'%(iy_3,ix_3))

times           = wrf.to_np(wrf.getvar(infile,'Times'      , wrf.ALL_TIMES))
NOx             = wrf.to_np(wrf.getvar(infile,'NOx'        , wrf.ALL_TIMES))                    # ppbv NOx = NO + NO2
LND_tracer      = wrf.to_np(wrf.getvar(infile,'LAND_tracer', wrf.ALL_TIMES))*1000.              # ppbv LND_tracer is emitted at a constant rate from all land cells
SEA_tracer      = wrf.to_np(wrf.getvar(infile,'SEA_tracer' , wrf.ALL_TIMES))*1000.              # ppbv SEA_tracer is emitted at a constant rate from all sea  cells
ntimes,nz,ny,nx = NOx.shape

# 2. --- Processing ----
#NOx                                      = (NOx-0.0115)*1000*1.35                   # ppb --> ug/m3 with background correction
NOx                                       = (NOx-0.0000)*1000*1.35                   # ppb --> ug/m3 without background correction
NOx_mod = dict()
LND_mod = dict()
SEA_mod = dict()
for station in stations:
    if station == 'Vredepeel':
       (iz,iy,ix) = (0,iy_1,ix_1)
    elif station == 'Zegveld':
       (iz,iy,ix) = (0,iy_2,ix_2)
    elif station == 'Wekerom':
       (iz,iy,ix) = (0,iy_3,ix_3)
    else:
        print('Warning: station not found!')
        (iz,iy,ix) = (-1,-1,-1)
    NOx_mod[station] = NOx       [:,iz,iy,ix]  # ug/m3
    LND_mod[station] = LND_tracer[:,iz,iy,ix]  # ppbv
    SEA_mod[station] = SEA_tracer[:,iz,iy,ix]  # ppbv
    

# 3. Output
COLORS = {'Vredepeel':'r', 'Zegveld':'g', 'Wekerom':'b'}
pp      = PdfPages('plt_timeSeries_Land_Sea_NOx.pdf') 
FS      = 10 # FontSize

if True:
    f1 = plt.figure(1,figsize=(15,9))
    f1.clf()
    ax = f1.add_subplot(111)
    for station in stations:
        ax.plot(times,SEA_mod[station],COLORS[station])
    ax.set_xlabel('Time',fontsize=FS)
    ax.set_ylabel('SEA tracer @Stations (ppbv)',fontsize=FS)
    ax.tick_params(axis='both', which='major', labelsize=FS)
    ax.tick_params(axis='both', which='minor', labelsize=FS)
    ax.set_title(name_project[-21:])
    ax.legend(stations)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%b-%Y_%H'))
    pp.savefig()

if True:
    f2 = plt.figure(2,figsize=(15,9))
    f2.clf()
    ax = f2.add_subplot(111)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%b-%Y_%H'))
    for station in stations:
        ax.plot(times,NOx_obs[station][:67],'-' ,marker='*',color=COLORS[station])
        ax.plot(times,NOx_mod[station]     ,'--',color=COLORS[station])
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.xaxis.set_tick_params(rotation = 0, labelsize = 8)
    ax.set_xbound(times[0], times[-1])
    ax.set_xlabel('Time',fontsize=FS)
    ax.set_ylabel('NOx (ug/m3)',fontsize=FS)
    ax.tick_params(axis='both', which='major', labelsize=FS)
    ax.tick_params(axis='both', which='minor', labelsize=FS)
    ax.set_title('WRF: '+name_project[-21:])
    ax.legend(('OBS Vredepeel','MOD Vredepeel','OBS Zegveld','MOD Zegveld','OBS Wekerom','MOD Wekerom'))
    pp.savefig()

pp.close()    
 
