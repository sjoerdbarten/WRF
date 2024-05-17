import os
import numpy as np
import netCDF4 as nc
from scipy.interpolate import griddata
from itertools import product
import datetime as dt
from dateutil.rrule import rrule, HOURLY

def msg(string,type=2):
    import sys
    if type == 1: print(string),
    if type == 2: print(string)
    sys.stdout.flush()
    return    

### User-specific settings ----------------------------------------------------------------------------------------------------
domains          = ['d01','d02','d03']                        # domains for which to prepare wrfchemi-files
start_date       = '2018-06-05 00:00:00'
end_date         = '2018-06-07 18:00:00'
### End User-specific settings ------------------------------------------------------------------------------------------------

### Other settings
 
sdate            = dt.datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S")
edate            = dt.datetime.strptime(end_date, "%Y-%m-%d %H:%M:%S")

homedir         = os.getenv("HOME")
MACCpath         = homedir+'/atmmod/EmissionData' #%homedir                   # location where the MACC Emission data files are located
geo_em_path      = homedir+'/atmmod/WPS'       #%homedir                      # location where the geo_em.dNN.nc files are located
wrfchemipath     = homedir+'/atmmod/WRF/run'                                    # location where the wrfchemifiles will be written to

tracers          = ['NOx','PM10','SO2','PM25','BC','CO','NO2','NO']                  # tracers to put in the wrfchemi-files, they will appear in the wrfchemi-file (in the same order)
tracers_extra    = ['LAND_tracer','SEA_tracer','tracer_stack_k01','tracer_stack_k05','tracer_stack_k15','tracer_stack_k20','VOC']
Nsnaps           = 13 

outunits         = {'PM25': 'kg/km2/hr', 'CO':'mole/km2/hr', 'PM10':'kg/km2/hr','CH4':'mole/km2/hr','SO2':'mole/km2/hr','NMVOC':'kg/km2/hr','NH3':'mole/km2/hr','NOx':'mole/km2/hr','BC':'kg/km2/hr','NO2':'mole/km2/hr','NO':'mole/km2/hr'}
outunits_extra   = {'LAND_tracer':'mole/km2/hr','SEA_tracer':'mole/km2/hr','tracer_stack_k01':'mole/km2/hr','tracer_stack_k05':'mole/km2/hr','tracer_stack_k15':'mole/km2/hr','tracer_stack_k20':'mole/km2/hr','VOC':'mole/km2/hr'}
    
emis_cat_names   = [['01: Energy industries                            '],   #  0        Stationary: 0,1,2,3,4,         10,11,12
                    ['02: Non-industrial combustion                    '],   #  1        Traffic:              5,6,7,8,9
                    ['34: Industry                                     '],   #  2
                    ['05: Fossil fuel production and distribution      '],   #  3 
                    ['06: Solvent and other product use                '],   #  4
                    ['71: Road transport, exhaust, gasoline            '],   #  5
                    ['72: Road transport, exhaust, diesel              '],   #  6
                    ['73: Road transport, exhaust, LPG and natural gas '],   #  7
                    ['74: Road transport, gasoline evaporation         '],   #  8
                    ['75: Road transport, tyre, brake and road wear    '],   #  9
                    ['08: Non-road transport                           '],   # 10  # Laurie: this includes shipping emissions
                    ['09: Waste                                        '],   # 11
                    ['10: Agriculture                                  '] ]  # 12     

#           f_EC    f_EC    f_OC    f_OC
#           avg     std     avg     std     # over countries
PM25TOBC = np.asarray([[0.071, 0.017, 0.088, 0.099],  #  1
                       [0.167, 0.036, 0.507, 0.257],  #  2
#                      [0.031, 0.021, 0.046, 0.012],  #  3
#                      [0.028, 0.029, 0.088, 0.082],  #  4
                       [0.030, 0.025, 0.067, 0.047],  # 34
                       [0.630, 0.000, 0.130, 0.000],  #  5
                       [0.050, 0.000, 0.900, 0.000],  #  6
                       [0.300, 0.000, 0.530, 0.000],  # 71
                       [0.720, 0.000, 0.210, 0.000],  # 72
                       [0.250, 0.113, 0.442, 0.200],  # 73
                       [0.000, 0.000, 0.000, 0.000],  # 74
                       [0.010, 0.000, 0.230, 0.000],  # 75
                       [0.419, 0.081, 0.406, 0.041],  #  8
                       [0.249, 0.078, 0.593, 0.062],  #  9
                       [0.095, 0.078, 0.780, 0.097]]) # 10
           
            
#--- input

if not 'progress' in globals(): progress = list()

if not 'TimeProfilesLoaded' in progress:

    TPfilename    = '%s/time_profiles_2011.nc'%MACCpath
    TPfile        = nc.Dataset(TPfilename,'r')     
    country_code  = TPfile.variables['country_code' ][:]
    country_name  = TPfile.variables['country_name' ][:]  
    emis_cat      = TPfile.variables['emis_cat'     ][:]     
    emis_cat_code = TPfile.variables['emis_cat_code'][:]
    emis_cat_name = TPfile.variables['emis_cat_name'][:]
    time          = TPfile.variables['time'         ][:]         
    date          = TPfile.variables['date'         ][:]         
    time_factors  = TPfile.variables['time_factors' ][:]
    TPfile.close()
    
    dates = [dt.datetime(sdate.year,x[1],x[2,],x[3],x[4],x[5]) for x in date]
    
    country_tmp = [x[0].decode('utf-8')+x[1].decode('utf-8')+x[2].decode('utf-8') for x in country_code]
    
    country_ids        = {k:v for v,k in enumerate(country_tmp)}
    ncountries         = len(country_ids.keys())
           
    progress.append('TimeProfilesLoaded')


if not 'MACCloaded' in progress:
    #--- init
    msg('Loading gridded MACC files: ',type=1)
    E_macc                     = dict()
    for tracer in (selectedtracers for selectedtracers in tracers if selectedtracers not in ['BC','NO2','NO']): 
        print(tracer)
        MACCncfilename         = '%s/MACC_%s.nc'%(MACCpath,tracer)
        msg('%s...'%MACCncfilename,type=1)
        MACCncfile             = nc.Dataset(MACCncfilename,'r')
        E_macc[tracer]         = MACCncfile.variables['E_%s'%tracer][:] # kg per area per hr or mole per area per hour
        if not 'lonmacc' in globals():
            lonmacc            = MACCncfile.variables['Longitude'][:]
            latmacc            = MACCncfile.variables['Latitude' ][:]
            country            = MACCncfile.variables['Country'  ][:]
            LONmacc,LATmacc    = np.meshgrid(lonmacc,latmacc)
            ny,nx              = LONmacc.shape
        MACCncfile.close()
    E_macc['NO2']              = E_macc['NOx'][:]*.05
    E_macc['NO']               = E_macc['NOx'][:]*.95    
    E_macc['BC']               = E_macc['PM25'][:] * 0
    for isnap in range(Nsnaps):
        E_macc['BC'][isnap] = E_macc['PM25'][isnap] * (PM25TOBC[isnap,0] + PM25TOBC[isnap,2])
    progress.append('MACCloaded') 
    
    
    
if not 'geo_em_loaded' in progress:
   msg('Loading lat/lon from geo_em...',type=1)
   LATwrf                       = dict()
   LONwrf                       = dict()
   GridResolution               = dict()
   for domain in domains:
       geofilename              = '%s/geo_em.%s.nc'%(geo_em_path,domain)
       geofile                  = nc.Dataset(geofilename,'r')       
       LATwrf[domain]           = geofile.variables['XLAT_M' ][0,:] # centerpoints
       LONwrf[domain]           = geofile.variables['XLONG_M'][0,:] # centerpoints
       GridResolution[domain]   = geofile.getncattr('DX')
       geofile.close()
   progress.append('geo_em_loaded')
   msg('done.')
                    
if not 'MACCinterpolatedToWRF' in progress:
    Xwrf                        = dict()
    Ywrf                        = dict()
    E_wrf                       = dict()
    Country_wrf                 = dict()
    for domain in domains:
        nlatwrf,nlonwrf         = LONwrf[domain].shape
        Country_wrf[domain]     = griddata((LONmacc.flatten(),LATmacc.flatten()),country.flatten() , (LONwrf[domain],LATwrf[domain]),method='nearest')
        if GridResolution[domain] > 7000.:
            print('XLONG_M')
	    # Note: d04 has a higher resolution than NEI --? interpolate or nearest neighbor
            #       d03 has the same resolution as NEI, but may be shifted in x or y direction --> nearest neighbour regridding
            #       d02 has a similar resolution as MACC --> nearest neighbor interpolation with griddata
            #       d01 has a coarser resolution than MACC  --> regrid cell by cell
            #E_wrf[domain,tracer][isnap,:,:] = griddata((LONmacc.flatten(),LATmacc.flatten()),E_macc[tracer][isnap].flatten() , (LONwrf[domain],LATwrf[domain]),method='nearest')
            msg('Interpolating %s..'%(domain),type=1)
            for tracer in tracers:
                E_wrf[domain,tracer] = np.zeros((Nsnaps,nlatwrf,nlonwrf))
            
            Xbl                 = LONwrf[domain]*np.nan # boundary left
            Xbr                 = LONwrf[domain]*np.nan # boundary right
            Ix                  = np.arange(1,nlonwrf)
            Xbl[:,Ix]           = (LONwrf[domain][:,Ix]+LONwrf[domain][:,Ix-1])/2.
            Ix                  = np.arange(nlonwrf-1)
            Xbr[:,Ix]           = (LONwrf[domain][:,Ix]+LONwrf[domain][:,Ix+1])/2.
                                
            Xbl[:, 0]           = LONwrf[domain][:, 0] - (Xbr[:, 0]-LONwrf[domain][:, 0])
            Xbr[:,-1]           = LONwrf[domain][:,-1] + (LONwrf[domain][:,-1]-Xbl[:,-1])
                                
            Ybb                 = LATwrf[domain]*np.nan # boundary bottom
            Ybt                 = LATwrf[domain]*np.nan # boundary top
            Iy                  = np.arange(1,nlatwrf)
            Ybb[Iy,:]           = (LATwrf[domain][Iy,:]+LATwrf[domain][Iy-1,:])/2.
            Iy                  = np.arange(nlatwrf-1)
            Ybt[Iy,:]           = (LATwrf[domain][Iy,:]+LATwrf[domain][Iy+1,:])/2.
                                
            Ybb[ 0,:]           = LATwrf[domain][ 0,:] - (Ybt[ 0,:]-LATwrf[domain][ 0,:])
            Ybt[-1,:]           = LATwrf[domain][-1,:] + (LATwrf[domain][-1,:]-Ybb[-1,:])
            
            where               = np.where
            t0                  = dt.datetime.now()
            for iy,ix in product(range(nlatwrf),range(nlonwrf)):
                Iy,Ix           = where( (LONmacc > Xbl[iy,ix]) & (LONmacc <= Xbr[iy,ix]) & (LATmacc > Ybb[iy,ix]) & (LATmacc <= Ybt[iy,ix]) )
                if len(Ix)>0:
                    for tracer in tracers:
                        E_wrf[domain,tracer][:,iy,ix] = E_macc[tracer][:,Iy.min():Iy.max()+1,Ix.min():Ix.max()+1].mean(axis=2).mean(axis=1)
            
                if ((np.remainder(iy,10) == 0) & (ix == nlonwrf-1)):
                    t1           = dt.datetime.now()
                    delt           = t1-t0
                    msg('iy=%3d of %3d; %7.4f seconds'%(iy,nlatwrf,delt.microseconds/1.e6)) # v0: 0.39 sec; v1 (Ilat outside loop): 0.007 sec
                    t0 = dt.datetime.now()
        elif GridResolution[domain] <= 7000.:
            msg('Interpolating %s..'%(domain),type=1)            
            for tracer in tracers:
                E_wrf[domain,tracer] = np.zeros((Nsnaps,nlatwrf,nlonwrf))
                for isnap in range(Nsnaps):
                    E_wrf[domain,tracer][isnap,:,:] = griddata((LONmacc.flatten(),LATmacc.flatten()),E_macc[tracer][isnap].flatten() , (LONwrf[domain],LATwrf[domain]),method='nearest')
        else: 
            msg('Resolution not found: %s'%domain)
    progress.append('MACCinterpolatedToWRF')
    msg('done.')
                                      
                  
if not 'wrfchemi_written' in progress:
    
    for domain,hour in product(domains,rrule(freq=HOURLY,dtstart=sdate,until=edate)):
 
        print('Writing wrfchemifile for %s %s '%(domain,hour))
     
        nlat,nlon                 = LATwrf[domain].shape 
        nz                        =  15

        TP                        = np.ones((1,Nsnaps,nlat,nlon))
	
        It                        = np.where(np.array(dates) == hour)
        It                        = It[0]
	
        for ic in range(ncountries):
            I                     = np.where(Country_wrf[domain] == ic)
            Ilat                  = I[0]
            Ilon                  = I[1]
            
            TP_country            = np.tile(time_factors[It,:,ic,np.newaxis,np.newaxis],(1,1,nlat,nlon)) # This doesn't seem to work on Enthought Python for Windows
            TP[:,:,Ilat,Ilon]     = TP_country[:,:,Ilat,Ilon]
        
        geofilename               = '%s/geo_em.%s.nc'%(geo_em_path,domain)
        geofile                   = nc.Dataset(geofilename,'r')
        wrfchemifilename          = '%s/wrfchemi_%s_%s'%(wrfchemipath,domain,hour.strftime("%Y-%m-%d_%H:%M:%S"))
        
        wrfchemifile              = nc.Dataset(wrfchemifilename, 'w',format = 'NETCDF3_CLASSIC' )

        wrfchemifile.createDimension('Time',        None) # will be 1
        wrfchemifile.createDimension('bottom_top' , nz)
        wrfchemifile.createDimension('south_north', nlat)
        wrfchemifile.createDimension('west_east'  , nlon)
        wrfchemifile.createDimension('DateStrLen' , 19)

     
        ncTimes                   = wrfchemifile.createVariable('Times'               ,'c',dimensions=('Time','DateStrLen'))
        ncXLAT                    = wrfchemifile.createVariable('XLAT'                ,'d',dimensions=('south_north','west_east'))
        ncXLONG                   = wrfchemifile.createVariable('XLONG'               ,'d',dimensions=('south_north','west_east'))
        ncTimes[:]                = [hour.strftime("%Y-%m-%d_%H:%M:%S")]
                                     
        ncXLAT.FieldType          = 104.
        ncXLAT.MemoryOrder        = 'XY'
        ncXLAT.description        = 'LATITUDE, SOUTH IS NEGATIVE'
        ncXLAT.units              = 'degree_north'
        ncXLAT.stagger            = ' ' 
        ncXLAT.coordinates        = 'XLONG XLAT'        
        ncXLAT[:]                 = LATwrf[domain]
        
        ncXLONG.FieldType         = 104.
        ncXLONG.MemoryOrder       = 'XY'
        ncXLONG.description       = 'LONGITUDE, WEST IS NEGATIVE'
        ncXLONG.units             = 'degree_north'
        ncXLONG.stagger           = ' ' 
        ncXLONG.coordinates       = 'XLONG XLAT' 
        ncXLONG[:]                = LONwrf[domain]
        
        for tracer in tracers:
            ncE_tracer            = wrfchemifile.createVariable('E_%s'%tracer          ,'d',dimensions=('Time','bottom_top','south_north','west_east'),fill_value = 0.)
            ncE_tracer.FieldType    = 104.
            ncE_tracer.MemoryOrder  = 'XYZ'
            ncE_tracer.description  = '%s emissions from stationary sources'%tracer
            ncE_tracer.units        = outunits[tracer]
            ncE_tracer.stagger      = ' ' 
            ncE_tracer.coordinates  = 'XLONG XLAT'   
            Isnap                   = np.asarray((0,1,2,3,4,5,6,7,8,9,10,11,12)) #         Stationary: 0,1,2,3,4,         10,11,12
            if (tracer == 'SO2'):
              print('We are doing SO2 now!')
              E_wrf[domain,tracer][10,:,:] = E_wrf[domain,tracer][10,:,:]*0.1 
            ncE_tracer[0,0,:,:]     = (TP[0,Isnap,:,:] * E_wrf[domain,tracer][Isnap,:,:]).sum(axis=0)

        xland = geofile.variables['LANDMASK'][:]
       
        for tracer in tracers_extra:
            ncE_tracer            = wrfchemifile.createVariable('E_%s'%tracer          ,'d',dimensions=('Time','bottom_top','south_north','west_east'),fill_value = 0.)
            ncE_tracer.FieldType    = 104.
            ncE_tracer.MemoryOrder  = 'XYZ'
            ncE_tracer.description  = '%s emissions from stationary sources'%tracer
            ncE_tracer.units        = outunits_extra[tracer]
            ncE_tracer.stagger      = ' ' 
            ncE_tracer.coordinates  = 'XLONG XLAT'   
            if tracer ==  'LAND_tracer':
                ncE_tracer[0,0,:,:] = xland[0,:,:]*1.0
            elif tracer ==  'SEA_tracer':
                ncE_tracer[0,0,:,:] = (1.-xland[0,:,:])*1.0
            elif tracer == 'tracer_stack_k01':
                ncE_tracer[0, 0,15,15] = 10.
            elif  tracer == 'tracer_stack_k05':
                ncE_tracer[0, 4,15,15] = 10.
            elif  tracer == 'tracer_stack_k15':
                ncE_tracer[0, 9,15,15] = 10.
            elif  tracer == 'tracer_stack_k20':   # Note: nz = 15, so emissions to the 20th level are not included in the wrfchemi file
                ncE_tracer[0,14,15,15] = 10.
            elif  tracer == 'VOC':
                ncE_tracer[0,0,:,:] = 0.
            else: 
                print('wrong tracer')
           
        
        
        wrfchemifile.TITLE        = "OUTPUT FROM WRFCHEMI V4.1"    
        wrfchemifile.CEN_LAT      = geofile.getncattr('CEN_LAT')
        wrfchemifile.CEN_LON      = geofile.getncattr('CEN_LON')
        wrfchemifile.TRUELAT1     = geofile.getncattr('TRUELAT1')
        wrfchemifile.TRUELAT2     = geofile.getncattr('TRUELAT2') 
        wrfchemifile.MOAD_CEN_LAT = geofile.getncattr('MOAD_CEN_LAT')
        wrfchemifile.STAND_LON    = geofile.getncattr('STAND_LON')
        wrfchemifile.POLE_LAT     = geofile.getncattr('POLE_LAT') 
        wrfchemifile.POLE_LON     = geofile.getncattr('POLE_LON')
#       wrfchemifile.GMT          = geofile.getncattr('GMT')  
#       wrfchemifile.JULYR        = geofile.getncattr('JULYR')
#       wrfchemifile.JULDAY       = geofile.getncattr('JULDAY')
        wrfchemifile.MAP_PROJ     = geofile.getncattr('MAP_PROJ')
        wrfchemifile.MMINLU       = geofile.getncattr('MMINLU')
        wrfchemifile.NUM_LAND_CAT = geofile.getncattr('NUM_LAND_CAT')
        wrfchemifile.ISWATER      = geofile.getncattr('ISWATER')
        wrfchemifile.ISLAKE       = geofile.getncattr('ISLAKE')   
        wrfchemifile.ISICE        = geofile.getncattr('ISICE')    
        wrfchemifile.ISURBAN      = geofile.getncattr('ISURBAN')     
        wrfchemifile.ISOILWATER   = geofile.getncattr('ISOILWATER')
        
        wrfchemifile.close()
    progress.append('wrfchemi_written')
    print ('The wrfchemi files are written to %s.'%wrfchemipath)
