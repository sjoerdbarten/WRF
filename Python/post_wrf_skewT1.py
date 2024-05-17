#   Example script to produce plots for a WRF real-data run, with the ARW coordinate dynamics option.
#   Plot SkewT's at a number of locations

# 1. Input

import os
import wrf
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units
np.bool = np.bool_

#User settings
##############

home       = os.path.expanduser('~')    # set home directory 
datapath   = '%s/atmmod/WRF/run'%home # set full datapath

# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,0)

# which domain is selected
domain = "d02"

#Select the time instances and model level to visualize

# step time in timesteps of wrfout 
# (i.e . if value is 6 every sixth time stamp in wrfout is plotted)
 
step_times =  1  

(station, station_lat, station_lon) = ('Cabauw (NL)', 51.93, 4.93)
Stations            = ('Cabauw (NL)',)
coords              = {'Cabauw (NL)':(51.93,4.93)}


# 1.4 Load variables from wrfoutfile

infilename = '%s/wrfout_%s_%4d-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)
infile     = nc.Dataset(infilename,'r')

name_project = infile.getncattr('TITLE')
xlat       = wrf.getvar(infile,'XLAT')
xlon       = wrf.getvar(infile,'XLONG')
times      = wrf.getvar(infile,'Times'      , wrf.ALL_TIMES)

plot_times = range(0,len(times),step_times)

inittime = times[0]
init_str = np.datetime_as_string(inittime.values,unit='s',timezone='UTC')

pp                  = PdfPages('plt_skewT1.pdf') 
for it in plot_times:

    curtime = times[it]
    cur_str = np.datetime_as_string(curtime.values,unit='s',timezone='UTC')
    
    print('Processing time: ',cur_str)

    Tc_wrf     = wrf.getvar(infile,'tc'      ,timeidx=it).values *units.degC                            # T in C
    Td_wrf     = wrf.getvar(infile,'td'      ,timeidx=it).values *units.degC                           # dew point temperature
    p_wrf      = wrf.getvar(infile,'pressure',timeidx=it).values *units.hPa                           # grid point pressure
    uvm        = wrf.getvar(infile,'uvmet'   ,timeidx=it)                            # umet and vmet averaged to mass points
    u_wrf      = uvm[0,:,:,:].values*1.94386 *units.knots                                                # This is a 4D array where, uvm(0,:,:,:) is umet, and uvm(1,:,:,:) is vmet, and this function rotate winds to earth coord.
    v_wrf      = uvm[1,:,:,:].values*1.94386  *units.knots                                               # extract u and v from uvm array, and turn wind into kts

    for station in Stations:
        lat_station,lon_station = coords[station]     # LOOP through above 20 station locations and
                                                      # plot a skewT if location is inside model domain

        # Get ij point in model domain for station
        # ; loc(0) is south-north (y) and loc(1) is west-east (x)
        (ix,iy)                 = wrf.to_np(wrf.ll_to_xy(infile,lat_station,lon_station))
      
      
      
        # Now we have the pressure, temperature and dew point temperature in the whole domain
        # Select one vertical column , t =0 , x=30, y=30
        p           = p_wrf[:,iy,ix]
        T           = Tc_wrf[:,iy,ix]
        Td          = Td_wrf[:,iy,ix]
        u           = u_wrf[:,iy,ix]     
        v           = v_wrf[:,iy,ix]
	
        lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
        lfc_pressure, lfc_temperature = mpcalc.lfc(p, T, Td)
        el_pressure, el_temperature = mpcalc.el(p, T, Td)
        prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
        cape, cin = mpcalc.cape_cin(p, T, Td,prof)
 
        fig = plt.figure(figsize=(10, 10))
        skew = SkewT(fig)

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
        skew.plot(p, T, 'r', linewidth=3,label='T')
        skew.plot(p, Td, 'b', linewidth=3,label='Td')
        skew.plot_barbs(p, u, v)
        skew.ax.set_ylabel('hPa')
        skew.ax.set_xlabel('Degree C')
#
# set the axes accorng to the temperature and 
# dew point temperature profile
#

        skew.ax.set_xlim([min(Td.magnitude)-10.,max(T.magnitude)+10.])	

#
# plot the lifting condensation level
#
        skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black',label='LCL')


#
# plot the profile of a parcel
#
        prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
        skew.plot(p, prof, 'k', linewidth=2,label='parcel profile')
       
        skew.ax.legend(loc='lower left')
# Add the relevant special lines

        skew.plot_dry_adiabats(color='purple',label='dry adiobats')
        skew.plot_moist_adiabats(color='saddlebrown',label='moist adiabats')
        skew.plot_mixing_lines(label='mixing lines (values in g/kg)')

        skew.ax.legend(loc='lower left')
        mixing_ratio = np.array([0.0004, 0.001, 0.002, 0.004, 0.007, 0.01,
                                     0.016, 0.024, 0.032]).reshape(-1, 1)
        td = mpcalc.dewpoint(mpcalc.vapor_pressure(600.*units.mbar, mixing_ratio))
       
        for mx,dewp in zip(mixing_ratio[:,0],td[:,0]):
             skew.ax.text(dewp,600.,str(1000.*mx),va='bottom',ha='center',color='g',clip_on=True)
        skew.ax.set_title(name_project[-21:]+': Skew-T log-P diagram \n')

        skew.ax.text(0,1,'Init: '+init_str,va='bottom',transform=skew.ax.transAxes)
        skew.ax.text(1,1,'valid: '+cur_str,ha='right',va='bottom',transform=skew.ax.transAxes)
        textstr= '\n'.join((
    'P_par=%.1f hPa' % (p[0].magnitude ),
    'T_LCL=%.1f hPa' % (lcl_pressure.magnitude),
    'P_LCL=%.1f degree C' % (np.nan_to_num(lcl_temperature.magnitude,nan=-999.)),
    'P_LFC=%.1f hPa' % (np.nan_to_num(lfc_pressure.magnitude,nan=-999.)),
    'T_LFC=%.1f degree C' % (np.nan_to_num(lfc_temperature.magnitude,nan=-999.)),
    'P_EL=%.2f hPa' % (np.nan_to_num(el_pressure.magnitude,nan=-999.)),
    'T_EL=%.2f degree C' % (np.nan_to_num(el_temperature.magnitude,nan=-999.)),
    'CAPE=%.2f J/kg' % (np.nan_to_num(cape.magnitude,nan=-999.)),
    'CIN=%.2f J/kg' % (np.nan_to_num(cin.magnitude,nan=-999.))
    )) 
    
        props = dict(boxstyle='round', facecolor='wheat')

# place a text box in upper left in axes coords
        skew.ax.text(0.05, 0.95, textstr, transform=skew.ax.transAxes,
        verticalalignment='top', bbox=props)

        pp.savefig()
        plt.close(fig)
pp.close()
infile.close()
