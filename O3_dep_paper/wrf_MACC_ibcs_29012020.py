#wrf_macc_ibcs.pyip
""" Routines that read MACC/CAMS 3D fields output and process it to format
usable to add as initial or boundary conditions in WRF. 

Authors : Ingrid Super (based on wrf_ct_ibcs.py)
Last changes: 19-12-2014 
Further modified by Laurens Ganzeveld, February 2019

#############################################################################
## ESSENTIAL INFORMATION
##
## This script is now set-up for reading in MACC files at a time-resolution
## of 6 hours, 4 timesteps per day! See specific edits compared to the code
## used to process the more recent CAMS files written out with dt=3hr.
## 
## RUN FIRST THE make_1x1_molefractions(dates....
## FUNCTION FOR 1 DAY EXTRA COMPARED TO THE wrfbdy_d* & wrfinput_d* FILES
## TO BE ABLE TO SUCCESSFULLY THEN 1) CONVERT (DON FORGET TO CHANGE DIRECTORY
## NAMES (unconverted/../), 2) EMPTY THE FILES AND THEN WRITE THE wrfbdy_d*
## AND wrfinput_d* FILES. THE PROCESSING SHOULD ALSO BE PROPERLY ENDED TO 
## SECURE WRITING THE UPDATED FILES PROPERLY!!!
#############################################################################
"""
import datetime as dtm
import os
import netCDF4 as nc
import numpy as np
import sys

# DONT FORGET TO CREATE ALL DIRECTORIES (also output/inverted/...)

inpath_macc = '/Storage/ganzeveld/models/WRF/input/EU/CAMS-MACC/2003/'
## MAQ_lg_20190411+ don forget to modify the path settings: unconverted vs ../directory after copying the unconverted files to ../
outpath_macc = '/Storage/ganzeveld/models/WRF/output/inverted/EU/2003/'
inpath_wrfinput = '/Storage/ganzeveld/models/WRF/input/EU/2003/'

#define start and end date as in macc input
#find the correct dates by typing the following to get normal dates from proleptic gregorian:
#>ncdump -v time -t [filename.nc]
#start_date=dtm.datetime(2016,12,20,0,0)
#end_date=dtm.datetime(2016,12,30,0,0)

# EU 2003 summer heatwave
#start_date=dtm.datetime(2003,7,31,0,0)
#end_date=dtm.datetime(2003,8,15,0,0)

# ASCOS campaign
#start_date=dtm.datetime(2009,2,01,0,0)
#end_date=dtm.datetime(2009,2,10,0,0)

# O3 study
#start_date=dtm.datetime(2008,8,10,0,0)
#end_date=dtm.datetime(2008,9,7,0,0)

# Warm Air Intrusion
start_date=dtm.datetime(2020,4,13,0,0)
end_date=dtm.datetime(2020,4,15,0,0)

#note that as of August 2014 the MACC system has changed (includes chemistry)
#the resolution has increased (~0.7 deg.); the only difference of importance
#here is that the hybrid A coefficient is already multiplied by 1000 (P(0))
#so this needs to be switched off here;
#units have changed from mole/mole to kg/kg, so use molar mass ratio correction

def calc_pres_fields(input,tstep,start_date):
    mf=nc.Dataset(input)
    zlen=len(mf.variables['level'])

    # determine the number of levels in the MACC files to properly use the actual a and b coefficients
    # in determining the height of these levels for the further interpolation. Check always carefully the 
    # pressure fields in the 3d_molefractions file which should be around 1000 hPa.

    nlevstart=int(60-zlen+1)

    # reading in file that contains the a and b parameters for the 60 lev version of the ECMWF/CAMS/MACC system
    file='/Storage/ganzeveld/models/WRF/input/CAMS-MACC/ab_ECMWF_L60.txt' #path
    hyab=np.genfromtxt(file,skip_header=2)
    acoef=hyab[:,0]
    bcoef=hyab[:,1]

    #print 'acoef: array size and values'
    acoef.shape #to evaluate the size of the array acoef...
    #print acoef
    #print bcoef

    xlen=len(mf.dimensions['longitude'])
    ylen=len(mf.dimensions['latitude'])
    
    P0=(mf.variables['sp'][:])
    Psurf=(mf.variables['sp'][:])

    print 'calc_pres_fields: zlen & nlevstart: ',zlen,nlevstart

    infield=zlen*[Psurf[tstep,:]]
    infield=np.asarray(infield)

    #the sigma pressure levs have a pressure of:
    #P=A*P(0)+B*P(s) where P(0) is a constant (1000 hPa) and P(s) is the surface pressure
    #A and B are coefficients for each layer

    ps=[]

    #print infield

    #after August 2014
    for k in range(zlen):
        kk=k+nlevstart # MAQ_lg_20190420+ corrected for the actual number of layers found in the CAMS/MACC file
        ps.append(acoef[kk]+bcoef[kk]*infield[k,:,:])
    #before August 2014
    #else:
    #for k in range(zlen):
    #    ps.append(acoef[k]*P0+bcoef[k]*infield[k,:,:])

    #these are pressures at half levs, similar to the eta levs in wrfinput

    return ps

print 'End calculation pressure levs'

def merge_netcdf_fields(input,dsname,tstep,massr,start_date,addtracers=False,okdebug=True):
    """ Merge TM5 resolutions, adapted from original hdf using code from tm5tools """
    from maptools import select_map

    m,nx,ny=select_map('Global Cylinder') # for global image plus merge_fields
    readfile=nc.Dataset(input)
    
    #print(nx,ny)
    #ny=10

    print readfile

    # get dimensions of dataset in file, also get coordinates of grid
    lonsmacc=readfile.variables['longitude'][:]
    im=len(readfile.dimensions['longitude'])
    for j in range(im):
       if lonsmacc[j]>180:
          lonsmacc[j]=lonsmacc[j]-360.
    jm=len(readfile.dimensions['latitude'])
    dx=abs(lonsmacc[1]-lonsmacc[0])
    dy=abs(readfile.variables['latitude'][1]-readfile.variables['latitude'][0])
    xbeg=lonsmacc[0]

    # MAQ_lg_20180119+ using the first line for the assignment of latitude puts the domain at the right location in global file
    #   this must be further checked
    #if start_date>dtm.datetime(2014,8,1,0,0):
    ybeg=readfile.variables['latitude'][-1]
    #else:
    #   ybeg=readfile.variables['latitude'][0]

    lm=len(readfile.dimensions['level'])

    if dsname=='sp':
       infield=calc_pres_fields(input=input,tstep=tstep,start_date=start_date)
       infield=np.asarray(infield)
       infield=infield[:,::-1,:] #reversion necessary because latitudinal direction is inverted
    if dsname=='co':
       infield=readfile.variables[dsname][:]*massr
       #infield=infield[tstep,:,:,:]
       infield=infield[tstep,:,::-1,:] #reversion necessary because latitudinal direction is inverted
    if dsname=='go3':
       infield=readfile.variables[dsname][:]*massr
       #infield=infield[tstep,:,:,:]
       infield=infield[tstep,:,::-1,:] #reversion necessary because latitudinal direction is inverted
    if dsname=='so2':
       infield=readfile.variables[dsname][:]*massr
       #infield=infield[tstep,:,:,:]
       infield=infield[tstep,:,::-1,:] #reversion necessary because latitudinal direction is inverted
    if dsname=='nox':
       infield=readfile.variables[dsname][:]*massr
       #infield=infield[tstep,:,:,:]
       infield=infield[tstep,:,::-1,:] #reversion necessary because latitudinal direction is inverted
    if dsname=='hcho':
       infield=readfile.variables[dsname][:]*massr
       #infield=infield[tstep,:,:,:]
       infield=infield[tstep,:,::-1,:] #reversion necessary because latitudinal direction is inverted

    readfile.close()

    target=np.zeros((lm,ny,360),dtype=np.float32)    

    # construct lats and lons arrays
    lons=xbeg+np.arange(im)*dx+dx/2.
    lats=ybeg+np.arange(jm)*dy+dy/2.

    # determine which cells to take from the global image
    istart=180+xbeg
    ie=istart+im*dx
    js=90+ybeg
    je=js+jm*dy
    #print(ybeg,js,je,jm,lats)


    # put field to global grid, all elements outside lon and lat are masked, set order=0 to subsample instead of interpolate
    outfield=[]
    #in contrast to CT, MACC data is from low to high pressure; to get correct outfields, turn it around!

    print 'routine merge_netcdf_fields, nx & ny: ',nx,ny
    print 'routine merge_netcdf_fields, infield.shape: ',infield.shape

    for l in (range(lm))[::-1]: 
        outfield.append(m.transform_scalar(infield[l,:,:],lons,lats,nx,ny,order=0))
    outfield=np.array(outfield)
    # grab non-masked and put into target array
    target[:,js:je,istart:ie]=outfield[:,js:je,istart:ie]

    return target

print 'End merge_netcdf_fields'

def make_1x1_molefractions(dates=(start_date,end_date),inpath=inpath_macc,outpath=outpath_macc):
    #    sys.path.append('/home/WUR/ganze004/soft/python/lib')
    #    sys.path.append('/Storage/ganzeveld/models/WRF/python')

    sys.path.append('/Users/visse218/Code/carbontracker/python')
    sys.path.append('/Users/visse218/Code/pythonlib/trunk')

    import tm5tools
    from ct_netcdf import CT_CDF, std_savedict
    import timedate  # this is a module from the pythonlib code from Wouter
    
    dectime0=tm5tools.date2num(dtm.datetime(2000,1,1))
    (sd,ed)=dates
    dtw=dtm.timedelta(days=1)
    dt = dtm.timedelta(hours=24)
    days = timedate.timegen(sd,ed,dtw)

    print 'make_1x1_molefractions: parameters days and dt:'
    print days, dt

    for day in days:
        saveas = os.path.join(outpath,'3d_molefractions_1x1_%s.nc'%day.strftime('%Y%m%d'))
        if os.path.exists(saveas):
            print 'Skipping existing file: %s' % saveas
            continue
        else:
            print 'Creating new file: %s' % saveas

        print 'Start reading file for day: '
	print day
        for i,hour in enumerate(range(0,24,6)):

            print 'Start reading file. loop i,hour: '
	    print i

            dd = day + dtm.timedelta(hours=hour)
            ed = dd + dt
### here comes the use of tm5tools function called mergefields, which purpose is to merge the data from the different resolutions on to 1x1 grid - we rewrote the function into new one merge_netcdf_fields included as well in this file
            #ncfile = input_macc
            #print ncfile 

## MAQ_lg_20190411+ modified to deal with the 6hr resolution MACC data with time for point #1 = 00.00 
### MACC files starting at 3:00, with a temporal resolution of 3hrs, require coding as outcomment now below to secure 
### a proper construction of the 3d_molefraction file. This is to secure start at 0:00 (which is in the file of the previous day)
            #if hour==0:
            #   new_day=day-dtm.timedelta(days=1)
            #   tstep=3  # last timestep of previous day: 7 for a temporal resolution of 3hrs (8 steps) and 3 for 6hr time resolution 
	                 # file of the previous day
            #else:
            new_day=day
            tstep=i-1

            print 'Start reading file: # of timesteps '
	    print tstep

            ncfile = os.path.join(inpath,'MACC_%s.nc'%(new_day.strftime('%Y%m%d')))
            #ncfile = os.path.join(inpath,'CAMS_20140101_merged.nc')
		
	    print 'reading MACC input file: '
            print ncfile
                
            # MACC using the varxxx codes
            # retrieval of limited selection of tracers (CO, SO2, NOx, O3, HCHO); 
            # param; CH4=62,S O2=122, CO=123, HCHO=124, NOx=129, O3=203 

            co = merge_netcdf_fields(input=ncfile,dsname='co',tstep=tstep,massr=1.,start_date=start_date)
            o3 = merge_netcdf_fields(input=ncfile,dsname='go3',tstep=tstep,massr=1.,start_date=start_date)
            so2 = merge_netcdf_fields(input=ncfile,dsname='so2',tstep=tstep,massr=1.,start_date=start_date)
	    no2 = merge_netcdf_fields(input=ncfile,dsname='nox',tstep=tstep,massr=1.,start_date=start_date)
	    hcho = merge_netcdf_fields(input=ncfile,dsname='hcho',tstep=tstep,massr=1.,start_date=start_date)

            # Create 3d pressure field, press=var134
            press = merge_netcdf_fields(input=ncfile,dsname='sp',tstep=tstep,massr=1.,start_date=start_date)
            nlev=press.shape[0]

            print 'finished creating 3d pressure field'

            if hour == 0:
            
                saveas = os.path.join(outpath,'3d_molefractions_1x1_%s.nc'%day.strftime('%Y%m%d'))

                # Create NetCDF output file
                #
                ncf=CT_CDF(saveas,'create')

                dimgrid=ncf.AddLatLonDim()
                dimlevs=(ncf.def_dim('levs',nlev),)
                dimdate=ncf.AddDateDim()
                dimidateformat=ncf.AddDateDimFormat()

                #
                # save NetCDF data
                #

                savedict=std_savedict.copy()
                savedict['name']='levs'
                savedict['long_name']='levs'
                savedict['values']=range(1,nlev+1)
                savedict['actual_range']=(1,nlev)
                savedict['dims']=dimlevs
                savedict['units']='lev'
                ncf.AddData(savedict)
            
            else:
                pass  
            # extra lines added for relevant species (AV)
            
            savedict=std_savedict.copy()
            savedict['name']='o3'
            savedict['values']=o3.tolist()
            savedict['dims']=dimdate+dimlevs+dimgrid
            savedict['units']='mol mol-1'
            savedict['long_name']='mole_fraction_of_ozone_in_air'
            savedict['standard_name']='mole_fraction_of_ozone_in_air'
            savedict['count']=i
            ncf.AddData(savedict)

            savedict=std_savedict.copy()
            savedict['name']='so2'
            savedict['values']=so2.tolist()
            savedict['dims']=dimdate+dimlevs+dimgrid
            savedict['units']='mol mol-1'
            savedict['long_name']='mole_fraction_of_sulphur_dioxide_in_air'
            savedict['standard_name']='mole_fraction_of_sulphur_dioxide_in_air'
            savedict['count']=i
            ncf.AddData(savedict)

            savedict=std_savedict.copy()
            savedict['name']='no2'
            savedict['values']=no2.tolist()
            savedict['dims']=dimdate+dimlevs+dimgrid
            savedict['units']='mol mol-1'
            savedict['long_name']='mole_fraction_of_nitric_dioxide_in_air'
            savedict['standard_name']='mole_fraction_of_nitric_dioxide_in_air'
            savedict['count']=i
            ncf.AddData(savedict)

            savedict=std_savedict.copy()
            savedict['name']='co'
            savedict['values']=co.tolist()
            savedict['dims']=dimdate+dimlevs+dimgrid
            savedict['units']='mol mol-1'
            savedict['long_name']='mole_fraction_of_carbon_monoxide_in_air'
            savedict['standard_name']='mole_fraction_of_carbon_monoxide_in_air'
            savedict['count']=i
            ncf.AddData(savedict)
            
            savedict=std_savedict.copy()
            savedict['name']='hcho'
            savedict['values']=hcho.tolist()
            savedict['dims']=dimdate+dimlevs+dimgrid
            savedict['units']='mol mol-1'
            savedict['long_name']='mole_fraction_of_formaldehyde_in_air'
            savedict['standard_name']='mole_fraction_of_formaldehyde_in_air'
            savedict['count']=i
            ncf.AddData(savedict)
            
            #savedict=std_savedict.copy()
            #savedict['name']='hno3'
            #savedict['values']=hno3.tolist()
            #savedict['dims']=dimdate+dimlevs+dimgrid
            #savedict['units']='mol mol-1'
            #savedict['long_name']='mole_fraction_of_nitric_acid_in_air'
            #savedict['standard_name']='mole_fraction_of_nitric_acid_in_air'
            #savedict['count']=i
            #ncf.AddData(savedict)
            
            #savedict=std_savedict.copy()
            #savedict['name']='iso'
            #savedict['values']=iso.tolist()
            #savedict['dims']=dimdate+dimlevs+dimgrid
            #savedict['units']='mol mol-1'
            #savedict['long_name']='mole_fraction_of_isoprene_in_air'
            #savedict['standard_name']='mole_fraction_of_isoprene_in_air'
            #savedict['count']=i
            #ncf.AddData(savedict)
            
            #savedict=std_savedict.copy()
            #savedict['name']='pan'
            #savedict['values']=pan.tolist()
            #savedict['dims']=dimdate+dimlevs+dimgrid
            #savedict['units']='mol mol-1'
            #savedict['long_name']='mole_fraction_of_peroxyacyl_nitrates_in_air'
            #savedict['standard_name']='mole_fraction_of_peroxyacyl_nitrates_in_air'
            #savedict['count']=i
            #ncf.AddData(savedict)
            
            savedict=std_savedict.copy()
            savedict['name']='pressure'
            savedict['values']=press.tolist()
            savedict['dims']=dimdate+dimlevs+dimgrid
            savedict['units']='Pa'
            savedict['long_name']='pressure_at_center_levs'
            savedict['standard_name']='air pressure'
            savedict['count']=i
            ncf.AddData(savedict)

            savedict=ncf.StandardVar(varname='idate')

## MAQ_lg_20190411+ modified to deal with the 6hr resolution MACC data
            # MAQ_lg_20190412+ modified definition of time-axis where in the orginal code dt/2 =12 was added
            #savedict['values']=(dd+dt/2).timetuple()[0:6]
            savedict['values']=(dd).timetuple()[0:6]

            savedict['dims']=dimdate+dimidateformat
            savedict['count']=i
            ncf.AddData(savedict)

            savedict=ncf.StandardVar(varname='date')

## MAQ_lg_20190411+ modified to deal with the 6hr resolution MACC data
            # MAQ_lg_20190412+ modified definition of time-axis where in the orginal code dt/2 =12 was added
            #savedict['values']=tm5tools.date2num(dd+dt/2)-dectime0
            savedict['values']=tm5tools.date2num(dd)-dectime0

            savedict['dims']=dimdate
            savedict['count']=i
            ncf.AddData(savedict)
 
        ncf.close()
    return None

print 'End merge_netcdf fields'

# MAQ_av_20180326+
def convert_1x1_molefractions(dates=(start_date,end_date),outpath=outpath_macc):
    """The 1x1_molefractions files should be converted from kg/kg to mol/mol using 
    the molar mass of air and the molar mass of the species"""
    Ma = 28.97e-3 #kg/mol

## MAQ_lg_20190411+ modified to deal with the MACC data: reduced number of tracers and no NO and NO2
    #specs = ['o3','so2','no2','co','hcho','dms','sulf']

    specs = ['o3','so2','no2','co','hcho']
    #specs = ['o3','no','no2','co','hcho','hno3','pan']

    #Mspec = [48.e-3,64.06e-3,46.01e-3,28.01e-3,30.03e-3,62.13e-3,96.06e-3] #kg/mol

    Mspec = [48.e-3,30.01e-3,46.01e-3,28.01e-3,30.03e-3] #kg/mol
    #Mspec = [48.e-3,30.01e-3,46.01e-3,28.01e-3,30.03e-3,63.01e-3,121.05e-3] #kg/mol
    
    mf_files = [os.path.join(outpath,filename) for filename in os.listdir(outpath) if filename.startswith('3d_molefractions_1x1_')]
    for each_path in mf_files:
        mf = nc.Dataset(each_path,mode='r+',format="NETCDF4")
        for i in range(0,len(specs)):
            var = mf.variables[specs[i]][:]
            var[:,:,:,:] *= Ma / Mspec[i]
            mf.variables[specs[i]][:] = var
        mf.close()
        print "Done converting the tracers in %s"%each_path
# MAQ_av_20180326-

def zero_tracers_in_wrf_ibcs(inpath,tracer_list=None,okdebug=True):
    """ Routine that zeroes the fields for supplied list of tracers in the WRF initial and boundary conditions - wrfinput & wrfbdy files"""
    if tracer_list == None: 
        print "No tracers supplied, please revise tracer_list in the input parameters"
        return None

    input_paths = [os.path.join(inpath,filename) for filename in os.listdir(inpath) if filename.startswith('wrfinput_')]
    for each_path in input_paths:
        if os.path.exists(each_path):
            ncfile = nc.Dataset(each_path, mode='r+')
            for each_tracer in tracer_list:
                if each_tracer in ncfile.variables.keys():
                    ncobj = ncfile.variables[each_tracer]
                    zero_tracer = np.zeros(ncobj.shape,float)
                    ncobj[:] = zero_tracer
                else:
                    if okdebug: print 'Tracer %s skipped - not in the input file.'%each_tracer
                    continue
            ncfile.close()
            print "Done zeroing the tracers in %s"%each_path
        else:
            if okdebug: print "Skipping %s - file path doesn't exist."%each_path
            continue
    bdy_path = os.path.join(inpath,'wrfbdy_d01')
    end_list = ['BXS', 'BXE', 'BYS', 'BYE', 'BTXS', 'BTXE', 'BTYS', 'BTYE']
    if os.path.exists(bdy_path):
        ncfile = nc.Dataset(bdy_path, mode='r+')
        for each_tracer in tracer_list:
            for each_ending in end_list:
                new_tracer = '_'.join([each_tracer,each_ending])                    
                if new_tracer in ncfile.variables.keys():
                    ncobj = ncfile.variables[new_tracer]
                    zero_tracer = np.zeros(ncobj.shape,float)
                    ncobj[:] = zero_tracer
                else:
                    if okdebug: print 'Tracer %s skipped - not in the bdy file.'%new_tracer
                    continue
        ncfile.close()
        print "Done zeroing the tracers in %s"%bdy_path
    else:
        if okdebug: print "Skipping %s - file path doesn't exist."%bdy_path
        pass
    return None

def add_bg_to_wrfinput(inpath_wrf,inpath_macc=outpath_macc,tracer_names=False,okdebug=True):
    """Routine that interpolates bgCO2 to the wrfinput and wrfbdy files as initial and boundary conditions"""
    from pylab import logical_and,logical_or,where
    from scipy.interpolate import griddata, interp1d
### loading the wrf input files one by one
    input_paths = [os.path.join(inpath_wrf,filename) for filename in os.listdir(inpath_wrf) if filename.startswith('wrfinput_')]

### calculating the initial conditions first!
    for each_path in input_paths:
        ncfile = nc.Dataset(each_path,mode='r+')
        init_times = ncfile.variables['Times'][:]
        itime = dtm.datetime(int(''.join(init_times[0][:4])),int(''.join(init_times[0][5:7])),int(''.join(init_times[0][8:10])),int(''.join(init_times[0][11:13])),int(''.join(init_times[0][14:16])),int(''.join(init_times[0][17:])))
### finding the correct maccfile and index for the time at wrfinput
        macc_filename = os.path.join(inpath_macc,'3d_molefractions_1x1_%s.nc'%itime.date().strftime('%Y%m%d'))
	print "add_bg_to_wrf: macc_filename: " 
	print macc_filename
        if not os.path.exists(macc_filename):
            if okdebug: print "Skipping %s - file path to MACC data doesn't exist."%macc_filename
            raise IOError,'MACC file with input data not available.'
        else:
            maccfile = nc.Dataset(macc_filename)
            macc_time = maccfile.variables['idate'][:]
	    
	    print 'add_bg_to_wrfinput, macc_time'
	    print macc_time
	    print macc_time[0][3]

            macc_dt = dtm.datetime(macc_time[1][0],macc_time[1][1],macc_time[1][2],macc_time[1][3],macc_time[1][4],macc_time[1][5])-dtm.datetime(macc_time[0][0],macc_time[0][1],macc_time[0][2],macc_time[0][3],macc_time[0][4],macc_time[0][5])
            imacc = None
            for i_macc in range(macc_time.shape[0]):
	    
	        print 'itime'
		print itime, i_macc
     
                #if ~logical_and(itime>=dtm.datetime(macc_time[i_macc][0],macc_time[i_macc][1],macc_time[i_macc][2],macc_time[i_macc][3],macc_time[i_macc][4],macc_time[i_macc][5])-macc_dt/2,itime< dtm.datetime(macc_time[i_macc][0],macc_time[i_macc][1],macc_time[i_macc][2],macc_time[i_macc][3],macc_time[i_macc][4],macc_time[i_macc][5])+macc_dt/2):
                #    print 'check logical'
                #    continue
                #else:
                if okdebug: "Wrf itime is found at macc_time with index %s"%i_macc
                imacc = i_macc
                break
            if imacc == None:
                if okdebug: "Wrf itime was not found in the macc file."
                raise ValueError,'Unknown time index'
### at this point we have the correct macc index for the required wrf input time
        print "Loading wrf grid data..."
### we need to load wrf grid and pressure levs, interpolate macc data to the wrf grid
        wrflat = ncfile.variables['XLAT'][:]
        wrflon = ncfile.variables['XLONG'][:]
        nx = wrflat.shape[2]
        ny = wrflat.shape[1]
        nt = wrflat.shape[0]

        lonmin = wrflon.min()-3
        lonmax = wrflon.max()+3
        latmin = wrflat.min()-3
        latmax = wrflat.max()+3
# this differs a bit from the matlab script mostly because python starts to count from 0 and matlab from 1
        Xi = np.arange(0,nx)
        Yi = np.arange(0,ny)
        [Xi,Yi] = np.meshgrid(Xi,Yi)
        Ixy_ic     = where((Xi>=0).flatten())[0].reshape((ny,nx))
      

#calculating WRF pressure fields at initial time
        eta = ncfile.variables['ZNU'][:]
        nz = eta.shape[1]
        ptop = ncfile.variables['P_TOP'][:]

        MUB = ncfile.variables['MUB'][:]
        MUP_ic = ncfile.variables['MU'][:]
        MU_ic = np.tile((MUB.flatten()[Ixy_ic].reshape(MUP_ic.shape)+MUP_ic )[:,np.newaxis,:,:] ,(1,nz,1,1))        
        outp_ic = np.tile(eta.reshape(1,nz,1,1),(nt,1,ny,nx)) * MU_ic + ptop ### gives 10x higher pressure?!?! check it out
# converting 2D lat lon to 3D (nz, ny, nz)
        wrflon_3d = np.tile(wrflon.flatten()[Ixy_ic].reshape(Ixy_ic.shape)[np.newaxis,:,:],(nz,1,1))
        wrflat_3d = np.tile(wrflat.flatten()[Ixy_ic].reshape(Ixy_ic.shape)[np.newaxis,:,:],(nz,1,1))

#interpolation from MACC to WRF for initial conditions
        VAR_ic = np.ones((nt,nz,ny,nx))*np.nan
        print "Making initial conditions for time %s"%itime.isoformat().replace('T','_')
#### load the MACC grid and find the start index and number of cells needed
        macclat = maccfile.variables['latitude'][:]
        macclon = maccfile.variables['longitude'][:]
        Ilon = where(logical_and(macclon.flatten()>lonmin,macclon.flatten()<lonmax))[0]
        Ilat = where(logical_and(macclat.flatten()>latmin,macclat.flatten()<latmax))[0]
        Lon1 = Ilon[0]
        Lat1 = Ilat[0]
        LonN = Ilon.size
        LatN = Ilat.size        
        p1 = 0
        pN = 36

        inlat = maccfile.variables['latitude'][Lat1:Lat1+LatN]
        inlon = maccfile.variables['longitude'][Lon1:Lon1+LonN]
        inp = maccfile.variables['pressure'][imacc,p1:p1+pN,Lat1:Lat1+LatN,Lon1:Lon1+LonN]

# put all variables in common dimensions (Z, LAT, LON)
        inlon = np.tile(inlon[np.newaxis,np.newaxis,:],(pN,LatN,1))
        inlat = np.tile(inlat[np.newaxis,:,np.newaxis],(pN,1,LonN))
# interpolating MACC data to WRF grid
        in_nx = inp.shape[2]
        in_ny = inp.shape[1]
        in_nz = inp.shape[0]
        incoord = np.array([inlon[0].flatten(),inlat[0].flatten()]).T
        outcoord = np.array([wrflon_3d[0].flatten(),wrflat_3d[0].flatten()]).T
        Ptmp = np.ones((in_nz,ny,nx))*np.nan

        for i,tracer_name in enumerate(tracer_names):
            inC = maccfile.variables[tracer_name][imacc,p1:p1+pN,Lat1:Lat1+LatN,Lon1:Lon1+LonN]
            outCtmp = np.ones((in_nz,ny,nx))*np.nan
# first a horizontal interpolation at each lev of macc
            for iz in range(in_nz):
                Cin = inC[iz].flatten()
                Pin = inp[iz].flatten()
#### ORIGINAL SCRIPT HAS LINEAR INTERPOLATION, doesn't work on maunaloa??
                outCtmp[iz] = griddata(incoord,Cin,outcoord,method='nearest').reshape(ny,nx)
                Ptmp[iz] = griddata(incoord,Pin,outcoord,method='nearest').reshape(ny,nx)

            print "Done with horizontal interpolation %s"%tracer_name
# secondly, vertical interpolation 
            outC = np.ones((nz,ny,nx))*np.nan
            for iy in range(ny):
                for ix in range(nx): #interp1d requires the pressure values to be increasing, so we reverse the arrays
                    Pin = Ptmp[:,iy,ix][::-1]
                    Pout = outp_ic[0,:,iy,ix][::-1]
                    Cin = outCtmp[:,iy,ix][::-1]
                    valid_Pout = where(logical_and(Pout<Pin[-1],Pout>=Pin[0]))[0]
                    invalid_Pout = where(Pout>=Pin[-1])[0]
                    invalid_Pout2 = where(Pout<Pin[0])[0]
                    f_vert = interp1d(Pin,Cin,kind='nearest') ## could be linear, works 
                    reverse_Cout = f_vert(Pout[valid_Pout])    
# we fix the values near surface of WRF outside the pressure of MACC to be equal to MACC surface values and reverse back the rest of the interpolated values
                    valid_Pout=np.asarray(valid_Pout)
                    outC[nz-1-invalid_Pout,iy,ix] = outCtmp[0,iy,ix]
                    outC[nz-1-invalid_Pout2,iy,ix] = 0
                    outC[nz-1-valid_Pout,iy,ix] = reverse_Cout
            print "Done with vertical interpolation %s"%tracer_name
            if where(np.isnan(outC))[0].size>0 :
                print "Nan found in output %s fields, please revise script!"%tracer_name
                raise ValueError("NaN found in output fields")
            VAR_ic[0] = outC

	    print 'outC'
            print outC

            tracer=ncfile.variables[tracer_name]
            tracer[:]=VAR_ic*1e6 #ppm
        ncfile.close()
        maccfile.close()
        print "Initial conditions written in %s"%each_path
    return None

def add_bg_to_wrfbdy(inpath_wrf,inpath_macc=outpath_macc,tracer_names=False,okdebug=True):
    """Routine that interpolates bgCO2 to the wrfinput and wrfbdy files as initial and boundary conditions"""
    from pylab import logical_and,logical_or,where,date2num
    from scipy.interpolate import griddata, interp1d

    ncfile = nc.Dataset(os.path.join(inpath_wrf,'wrfinput_d01'))
    bdyfile = nc.Dataset(os.path.join(inpath_wrf,'wrfbdy_d01'),mode='r+')
    bdy_times = bdyfile.variables['Times'][:]
    btimes = [dtm.datetime(int(''.join(bdy_times[b][:4])),int(''.join(bdy_times[b][5:7])),int(''.join(bdy_times[b][8:10])),int(''.join(bdy_times[b][11:13])),int(''.join(bdy_times[b][14:16])),int(''.join(bdy_times[b][17:]))) for b in range(bdy_times.shape[0])]
    dt = (date2num(btimes[-1])-date2num(btimes[-2]))*86400 # interval between calls in seconds - needed for the tendencies
    btimes.append(btimes[-1]+(btimes[-1]-btimes[-2])) # add one extra time step needed to calculate the tendency of the last time step
    print "Loading wrf grid data..."
### load the wrf data
    wrflat = ncfile.variables['XLAT'][:]
    wrflon = ncfile.variables['XLONG'][:]
    eta = ncfile.variables['ZNU'][:]
    ptop = ncfile.variables['P_TOP'][:]

    lonmin = wrflon.min()-3
    lonmax = wrflon.max()+3
    latmin = wrflat.min()-3
    latmax = wrflat.max()+3

    nx = wrflat.shape[2]
    
    ny = wrflat.shape[1]
    nt = len(btimes)
    nz = eta.shape[1]
    nb = len(bdyfile.dimensions['bdy_width'])

# this differs a bit from the matlab script mostly because python starts to count from 0 and matlab from 1
    Xi = np.arange(0,nx)
    Yi = np.arange(0,ny)
    [Xi,Yi] = np.meshgrid(Xi,Yi)
    Ixy = {}
#    Ixy['ic']     = where((Xi>=0).flatten())[0].reshape((ny,nx))
    Ixy['BXS']    = where((Xi<nb).flatten())[0].reshape((ny,nb))
    Ixy['BXE']    = where((Xi>=nx-nb).flatten())[0].reshape((ny,nb))
    Ixy['BYS']    = where((Yi<nb).flatten())[0].reshape((nb,nx))
    Ixy['BYE']    = where((Yi>=ny-nb).flatten())[0].reshape((nb,nx))    

    bdys = ['BXS','BXE','BYS','BYE']
    tdys = ['BTXS','BTXE','BTYS','BTYE']
#calculating WRF pressure fields for the entire period
    print "Calculating WRF pressure fields.."
    MUB = ncfile.variables['MUB'][0]
    MUP = {}
    MU = {}
    outp={}
    for i in range(len(bdys)):
        bdy = bdys[i]
        MUP[bdy] = np.zeros((bdyfile.variables['MU_'+bdy].shape[0]+1,bdyfile.variables['MU_'+bdy].shape[1],bdyfile.variables['MU_'+bdy].shape[2]))
        MUP[bdy][:-1] = bdyfile.variables['MU_'+bdy][:]
        MUP[tdys[i]] = bdyfile.variables['MU_'+tdys[i]][:]
        MUP[bdy][-1] = MUP[bdy][nt-1]+MUP[tdys[i]][nt-2]*dt 
#pressure state is not available at the last time step, therefore we calculate that one from the tendencies

    for i in range(len(bdys)):
        bdy = bdys[i]
        if bdy in ['BXS','BXE']:
            aMU = MUB.flatten()[Ixy[bdy]]
            MU[bdy] = np.tile((np.tile(np.transpose(aMU,(1,0))[np.newaxis,:,:],(nt,1,1))+MUP[bdy])[:,:,np.newaxis,:],(1,1,nz,1))
            nh = ny
            nttmp = nt
        elif bdy in ['BYS','BYE']:
            aMU = MUB.flatten()[Ixy[bdy]]
            MU[bdy] = np.tile((np.tile(aMU[np.newaxis,:,:],(nt,1,1))+MUP[bdy])[:,:,np.newaxis,:],(1,1,nz,1))
            nh = nx
            nttmp = nt      
        outp[bdy]= np.tile(eta.reshape(1,1,nz,1),(nttmp,nb,1,nh)) * MU[bdy] + ptop

# converting 2D lat lon to 3D (nz, ny, nz)
    wrflon_3d = {}
    wrflat_3d = {}
    for i in range(len(bdys)):
        bdy = bdys[i]
        wrflon_3d[bdy] = np.tile(wrflon.flatten()[Ixy[bdy]][np.newaxis,:,:],(nz,1,1))
        wrflat_3d[bdy] = np.tile(wrflat.flatten()[Ixy[bdy]][np.newaxis,:,:],(nz,1,1))
### initialize variables to calculate and write
    VAR = {}
    for i in range(len(bdys)):
        bdy = bdys[i]
        tdy = tdys[i]
        if bdy in ['BXS','BXE']:
            VAR[bdy] = np.ones((nt,nb,nz,ny))*np.nan
            VAR[tdy] = np.ones((nt,nb,nz,ny))*np.nan
        elif bdy in ['BYS','BYE']:
            VAR[bdy] = np.ones((nt,nb,nz,nx))*np.nan
            VAR[tdy] = np.ones((nt,nb,nz,nx))*np.nan
### loop over times
    for j,tracer in enumerate(tracer_names):
        for it in range(nt):
            itime = btimes[it]
            print "Starting calculations for time step %s at %s"%(itime.isoformat(),str(dtm.datetime.utcnow()))
### finding the correct maccfile and index for the time at wrfinput
            macc_filename = os.path.join(inpath_macc,'3d_molefractions_1x1_%s.nc'%itime.date().strftime('%Y%m%d'))

            print 'add_bg_to_wrfbdy, MACC filename'
	    print macc_filename

            if not os.path.exists(macc_filename):
                if okdebug: print "Skipping %s - file path to MACC/CAMS data doesn't exist."%macc_filename
                raise IOError,'MACC file with input data not available.'
            else:
                maccfile = nc.Dataset(macc_filename)
                macc_time = maccfile.variables['idate'][:]
                macc_dt = dtm.datetime(macc_time[1][0],macc_time[1][1],macc_time[1][2],macc_time[1][3],macc_time[1][4],macc_time[1][5])-dtm.datetime(macc_time[0][0],macc_time[0][1],macc_time[0][2],macc_time[0][3],macc_time[0][4],macc_time[0][5])
                imacc = None
                print 'add_bg_to_wrfbdy, itime, macc_time and macc_dt'
	        print itime,macc_time,macc_dt
                for i_macc in range(macc_time.shape[0]):

                    print 'add_bg_to_wrfbdy: i_macc, itime, condition met?'
	            print i_macc,itime,~logical_and(itime>=dtm.datetime(macc_time[i_macc][0],macc_time[i_macc][1],macc_time[i_macc][2],macc_time[i_macc][3],macc_time[i_macc][4],macc_time[i_macc][5])-macc_dt*0,itime<dtm.datetime(macc_time[i_macc][0],macc_time[i_macc][1],macc_time[i_macc][2],macc_time[i_macc][3],macc_time[i_macc][4],macc_time[i_macc][5])-macc_dt*0+macc_dt/2)
		    print 'add_bg_to_wrfbdy, dtm.datetime'
		    print dtm.datetime(macc_time[i_macc][0],macc_time[i_macc][1],macc_time[i_macc][2],macc_time[i_macc][3],macc_time[i_macc][4],macc_time[i_macc][5])-macc_dt*0,dtm.datetime(macc_time[i_macc][0],macc_time[i_macc][1],macc_time[i_macc][2],macc_time[i_macc][3],macc_time[i_macc][4],macc_time[i_macc][5])-macc_dt*0+macc_dt/2
		    print imacc
		    
                    if ~logical_and(itime>=dtm.datetime(macc_time[i_macc][0],macc_time[i_macc][1],macc_time[i_macc][2],macc_time[i_macc][3],macc_time[i_macc][4],macc_time[i_macc][5])-macc_dt*0,itime<dtm.datetime(macc_time[i_macc][0],macc_time[i_macc][1],macc_time[i_macc][2],macc_time[i_macc][3],macc_time[i_macc][4],macc_time[i_macc][5])-macc_dt*0+macc_dt/2):
                        print 'Continue i_macc loop'
                        continue
                    else:
                        if okdebug: print "Wrf itime is found at macc_time with index %s"%i_macc
                        imacc = i_macc
			print 'add_bg_to_wrfbdy, okdebug'
			print i_macc
                        break
                if imacc == None:
                    if okdebug: print "Wrf itime was not found in the maccfile."
                    raise ValueError,'Unknown time index'
### at this point we have the correct macc index for the required wrf time step
            print "Loading MACC grid data..."

#### load the macc grid and find the start index and number of cells needed
            macclat = maccfile.variables['latitude'][:]
            macclon = maccfile.variables['longitude'][:]
            Ilon = where(logical_and(macclon.flatten()>lonmin,macclon.flatten()<lonmax))[0]
            Ilat = where(logical_and(macclat.flatten()>latmin,macclat.flatten()<latmax))[0]
            Lon1 = Ilon[0]
            Lat1 = Ilat[0]
            LonN = Ilon.size
            LatN = Ilat.size        
            p1 = 0
            pN = 36

            inlat = maccfile.variables['latitude'][Lat1:Lat1+LatN]
            inlon = maccfile.variables['longitude'][Lon1:Lon1+LonN]
            inp = maccfile.variables['pressure'][imacc,p1:p1+pN,Lat1:Lat1+LatN,Lon1:Lon1+LonN]
            inC = maccfile.variables[tracer][imacc,p1:p1+pN,Lat1:Lat1+LatN,Lon1:Lon1+LonN]
		
            inlon = np.tile(inlon[np.newaxis,np.newaxis,:],(pN,LatN,1))
            inlat = np.tile(inlat[np.newaxis,:,np.newaxis],(pN,1,LonN))
            print "Interpolating macc to WRF grid.."

# interpolating macc data to WRF grid
            for i in range(len(bdys)):
                bdy = bdys[i]
                tdy = tdys[i]

                in_nx = inp.shape[2]
                in_ny = inp.shape[1]
                in_nz = inp.shape[0]				
                incoord = np.array([inlon[0].flatten(),inlat[0].flatten()]).T
                outcoord = np.array([wrflon_3d[bdy][0].flatten(),wrflat_3d[bdy][0].flatten()]).T
                Ptmp = np.ones((in_nz,wrflon_3d[bdy].shape[1],wrflon_3d[bdy].shape[2]))*np.nan
                outCtmp = np.ones((in_nz,wrflon_3d[bdy].shape[1],wrflon_3d[bdy].shape[2]))*np.nan
# first a horizontal interpolation at each level of macc
                for iz in range(in_nz):
                    Zin = inC[iz].flatten()
                    Pin = inp[iz].flatten()
#### ORIGINAL SCRIPT HAS LINEAR INTERPOLATION, doesn't work on maunaloa??
                    outCtmp[iz] = griddata(incoord,Zin,outcoord,method='nearest').reshape(wrflon_3d[bdy].shape[1],wrflon_3d[bdy].shape[2])
                    Ptmp[iz] = griddata(incoord,Pin,outcoord,method='nearest').reshape(wrflon_3d[bdy].shape[1],wrflon_3d[bdy].shape[2])

                print "Done with horizontal interpolation"
# secondly, vertical interpolation 
            
                outC = np.ones((wrflon_3d[bdy].shape[1],nz,wrflon_3d[bdy].shape[2]))*np.nan
                #print wrflon_3d[bdy].shape
                for iy in range(wrflon_3d[bdy].shape[1]):
                    for ix in range(wrflon_3d[bdy].shape[2]): #interp1d requires the pressure values to be increasing, so we reverse the arrays
                        Pin = Ptmp[:,iy,ix][::-1]
                        #print outp[bdy].shape
                        #print it,iy,ix
                        if bdy in ['BXS','BXE']:
                           Pout = outp[bdy][it,ix,:,iy][::-1]
                        elif bdy in ['BYS','BYE']:
                           Pout = outp[bdy][it,iy,:,ix][::-1]
                        Cin = outCtmp[:,iy,ix][::-1]
                        valid_Pout = where(logical_and(Pout<Pin[-1],Pout>=Pin[0]))[0]
                        invalid_Pout = where(Pout>=Pin[-1])[0]
                        invalid_Pout2 = where(Pout<Pin[0])[0]
                        f_vert = interp1d(Pin,Cin,kind='nearest') ## could be linear, works 
                        reverse_Cout = f_vert(Pout[valid_Pout])    
# we fix the values near surface of WRF outside the pressure of macc to be equal to macc surface values and reverse back the rest of the interpolated values
                        valid_Pout=np.asarray(valid_Pout)
                        outC[iy,nz-1-invalid_Pout,ix] = outCtmp[0,iy,ix]
                        outC[iy,nz-1-invalid_Pout2,ix] = 0
                        outC[iy,nz-1-valid_Pout,ix] = reverse_Cout
                print "Done with vertical interpolation"
                if where(np.isnan(outC))[0].size>0 :
                    print "Nan found in output %s fields, please revise script!"%tracer_name
                    raise ValueError("NaN found in output fields")
                if bdy in ['BXS','BXE']: 
                    VAR[bdy][it] = np.transpose(outC,(2,1,0))
                elif bdy in ['BYS','BYE']:
                    VAR[bdy][it] = np.transpose(outC,(0,1,2))
                 
### calculate tendencies
        print "Calculating tendencies.."
        for i in range(len(bdys)):
            bdy = bdys[i]
            tdy = tdys[i]
            VAR[tdy] = (VAR[bdy][1:nt]-VAR[bdy][:nt-1])/dt

### for 'BXE', 'BYE' we need to reverse the direction of the boundary grid
        VAR['BXE'] = VAR['BXE'][:,::-1,:,:]
        VAR['BTXE'] = VAR['BTXE'][:,::-1,:,:]
        VAR['BYE'] = VAR['BYE'][:,::-1,:,:]
        VAR['BTYE'] = VAR['BTYE'][:,::-1,:,:]
    
#writing the bc and tendencies in the bdy

        for key in VAR.keys():
            varname = tracer+'_'+key
            tracern = bdyfile.variables[varname] 
            tracern[:] = VAR[key][:nt-1]*1e6 #ppm
    bdyfile.close()
    ncfile.close()
    maccfile.close()
    print "Boundary conditions written in %s"%os.path.join(inpath_wrf,'wrfbdy_d01')

    return None

# First make the 3d_molefractions_1x1_ file (outcommenting the calls of the other functions)
# writing the 3d_molecfractions_1x1_ file to a subdirectory unconverted
# then run the conversion and indeed create the wrfbdy file using only ONCE the convert subtoutine!!
# DONT FORGET TO COPY THE 3d_molefractions FILES FROM THE DIR. UNCONVERTED TO THE DIRECTORY ../
# AND TO CHANGE THE DIRECTORY NAME ABOVE OF USING THE 3d_molefractions FILE

#make_1x1_molefractions(dates=(start_date,end_date),inpath=inpath_macc,outpath=outpath_macc)
#convert_1x1_molefractions(dates=(start_date,end_date),outpath=outpath_macc)

## MAQ_lg_20190411+ modified to deal with the limited no of tracers in MACC, no NO and NO2, instead NOx represented by NO2
#tracers=['o3','so2','no2','co','hcho','dms','sulf']

#tracers=['O3','NO','NO2','CO']
tracers=['o3','so2','no2','co','hcho']

#tracers=['o3','no2','no','co','hcho','hno3','pan']

zero_tracers_in_wrf_ibcs(inpath_wrfinput,tracer_list=tracers,okdebug=True)
add_bg_to_wrfinput(inpath_wrf=inpath_wrfinput,inpath_macc=outpath_macc,tracer_names=tracers,okdebug=True)
#add_bg_to_wrfbdy(inpath_wrf=inpath_wrfinput,inpath_macc=outpath_macc,tracer_names=tracers,okdebug=True)
