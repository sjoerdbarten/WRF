import numpy as np
import netCDF4 as nc


# User settings -----------
domains              = ('d01','d02','d03')#,'d04')
# End user settings

WRFpath              = '../WRF/run'     #%homedir
tracers              = ('NOx','PM10','SO2','PM25','NO','NO2','O3','LAND_tracer','SEA_tracer','BC','Rn','CO','tracer_stack_k01','tracer_stack_k05','tracer_stack_k15','tracer_stack_k20','VOC','tracer_18','tracer_19','tracer_20','tracer_ens')
bdys                 = ('BXS','BXE','BYS','BYE')
btdys                = ('BTXS','BTXE','BTYS','BTYE')


ibc = dict()
ibc['NOx']           = 0.0030   # NOx   ppm
ibc['PM10']          = 20.0     # PM10  ug/m3     
ibc['SO2']           = 0.0014   # SO2   ppm    
ibc['PM25']          =  5.0     # PM25  ug/m3     
ibc['NO']            = 0.0003   # NO    ppm     ~ 4 ug/m3
ibc['NO2']           = 0.0027   # NO2   ppm     ~ 20 ug/m3
ibc['O3_1']          = 0.0300   # O3    ppm     
ibc['O3_2']          = 0.1000   # O3    ppm     
ibc['O3_3']          = 0.3000   # O3    ppm     
ibc['LAND_tracer']   = 0.0000   # LAND_tracer
ibc['SEA_tracer']    = 0.0000   # SEA_tracer 
ibc['BC']            = 0.2      # BC    ug/m3     ~0.5 ug/m3
ibc['Rn']            = 0.5      # Rn    Bq/m3     ~0.5 Bq/m3
ibc['CO']            = 0.1500   # CO    ppm
ibc['tracer_stack_k01']     = 0.0000   # tracer_13  
ibc['tracer_stack_k05']     = 0.0000   # tracer_14  
ibc['tracer_stack_k15']     = 0.0000   # tracer_15  
ibc['tracer_stack_k20']     = 0.0000   # tracer_16     # Note: in chem kemit = nz = 15, so emissions to the 20th level are not included in the wrfchemi file
ibc['VOC']     = 0.0001   # tracer_17  
ibc['tracer_18']     = 0.0001   # tracer_18  
ibc['tracer_19']     = 0.0001   # tracer_19  
ibc['tracer_20']     = 0.0001   # tracer_20  
ibc['tracer_ens']    = 0.0001   # tracer_ens 

# wrfinput
for domain in domains:
    print('ic:'),
    icf              = nc.Dataset('%s/wrfinput_%s'%(WRFpath,domain),'r+')
    for tracer in tracers:
        print('%s '%tracer),
        var          = icf.variables[tracer]
        shape        = var.shape
        if tracer != 'O3':
            var[:]       = np.ones(shape)*ibc[tracer]
        else:
            print(var.dimensions)
            var[:]           = np.zeros(shape)
            var[:,  :10,:,:] = ibc['O3_1']
            var[:,10:20,:,:] = ibc['O3_2']
            var[:,20:  ,:,:] = ibc['O3_3']

    icf.close()   
print('done.')

# wrfbdy
bcfilename = '%s/wrfbdy_d01'%WRFpath
print('bc: '+bcfilename)
bcf                  = nc.Dataset('%s/wrfbdy_d01'%WRFpath,  'r+')
for tracer in tracers:
    print('%s '%tracer)
    for bdy in bdys:
        var          = bcf.variables['%s_%s'%(tracer,bdy)]
        shape        = var.shape
#      print(shape)
        if tracer != 'O3':
            var[:]   = np.ones(shape) * ibc[tracer]
        else:
            var[:,:,  :10,:] = ibc['O3_1']
            var[:,:,11:20,:] = ibc['O3_2']
            var[:,:,20:  ,:] = ibc['O3_3']
    for btdy in btdys:
        var          = bcf.variables['%s_%s'%(tracer,btdy)]
        shape        = var.shape
        var[:]       = np.zeros(shape) 
            
bcf.close()
print('done.')
print('The tracer concentrations have been updated in the wrfinput_dnn.nc and wrfbdy_d01.nc files on %s'%WRFpath)


  
