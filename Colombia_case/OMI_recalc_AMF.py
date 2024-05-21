"""OMI_recalc_AMF.py

This script takes:
- ods: an OMI netCDF dataset
- wds: a WRF-Chem netCDF dataset
- wrf_lat/wrf_lon: pixel lat/lon data from the WRF-Chem cell of interest

checks if:
- the search distance is within the limit
- the scanline yields valid data

calculates:
- an interpolated NO2 mixing ratio profile at a detailed vertical grid 
  based on TM5 pressure
- the sub-columns at the detailed pressure grid and subsequently for the
  TM5 layers
- the adjusted AMF based on the WRF-Chem sub-column profile evaluated at
  TM5 pressure
- the modified tropospheric VCD based on the adjusted AMF

and outputs:
- the adjusted tropospheric AMF
- the adjusted tropospheric VCD
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy import interpolate

def recalc_AMF(ods,wds,scanl,gr_pix,wrf_lat,wrf_lon,time,debug_lev):
    
    ### Extract OMI variables at the pixel ###
    prod = ods.groups['PRODUCT']
    detres = ods.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS']
    
    P_TM5_surf = prod.variables['tm5_surface_pressure'][0,scanl,gr_pix]
    trop_lev   = prod.variables['tm5_tropopause_layer_index'][0,scanl,gr_pix]
    SCD        = detres.variables['scd_no2'][0,scanl,gr_pix]
    VCD_strat  = detres.variables['stratospheric_no2_vertical_column'][0,scanl,gr_pix]
    AMF_strat  = detres.variables['amf_strat'][0,scanl,gr_pix]
    VCD_trop   = prod.variables['tropospheric_no2_vertical_column'][0,scanl,gr_pix]
    AMF_trop   = prod.variables['amf_trop'][0,scanl,gr_pix]
    AMF_total  = prod.variables['amf_total'][0,scanl,gr_pix]
    AK         = prod.variables['averaging_kernel'][0,scanl,gr_pix,0:trop_lev+1]
    P_TM5      = prod.variables['tm5_pressure_level_a'][:,1] + prod.variables['tm5_pressure_level_b'][:,1] * P_TM5_surf * 100
    
    #Abort program if a processing flag is raised
    #err_flag = prod.variables['processing_error_flag'][0,:,:]
    #if err_flag[scanl,gr_pix] > 0.:
    #    print 'Processing error raised'
    #    AMF_trop_adj = np.nan
    #    VCD_trop_adj = np.nan
    #    dist = np.nan
    #    return AMF_trop_adj, VCD_trop_adj
    
    ### Extract WRF-Chem variables ###
    WRF_latitude  = wds.variables['XLAT'][0,:,:]
    WRF_longitude = wds.variables['XLONG'][0,:,:]
    
    #Find indices of WRF pixel
    wrf_latlon_tup = np.where((WRF_latitude == wrf_lat) & (WRF_longitude == wrf_lon))
    lat_ind        = wrf_latlon_tup[0][0]
    lon_ind        = wrf_latlon_tup[1][0]
    
    wrf_no2  = wds.variables['no2'][time,:,lat_ind,lon_ind]
    wrf_Pres = wds.variables['P'][time,:,lat_ind,lon_ind] + wds.variables['PB'][time,:,lat_ind,lon_ind]
    wrf_Psfc = wds.variables['PSFC'][time,lat_ind,lon_ind]
    
    ### dP-based calculation of sub-columns ###
    
    #Interpolate WRF-Chem NO2 column to the TM5 grid
    #max_P_surf = max(wrf_Psfc,P_TM5_surf * 100.)
    interp     = interpolate.interp1d(np.append(wrf_Psfc,wrf_Pres),np.append(wrf_no2[0],wrf_no2),bounds_error=False,fill_value = 0.) #MAQ_AV_13112017+ changed from 0 / wrf_no2[0] 
    #interp    = interpolate.interp1d(wrf_Pres,wrf_no2) #MAQ_AV_13112017+ changed from 0 / wrf_no2[0] 
    #print np.append(max_P_surf,wrf_Pres)
    #print np.append(max_P_surf,P_TM5[0:trop_lev+1])
    
    num_points = 50
    #P_det = np.linspace(wrf_Psfc,P_TM5[0],num_points,endpoint=False)
    P_det = np.linspace(P_TM5_surf * 100.,P_TM5[0],num_points,endpoint=False)
    for i in range(0,P_TM5[0:trop_lev+1].shape[0]-1):
        x = np.linspace(P_TM5[i],P_TM5[i+1],num_points,endpoint=False)
        P_det = np.append(P_det,x)
    
    no2_det = interp(P_det)
    
    #Plot NO2 profiles
    if debug_lev >= 100.:
        plt.figure()
        plt.plot(wrf_no2*1000.,wrf_Pres / 100., 'kx-',label='$[NO_2] \ | \ P_{WRF-Chem}$')
        plt.plot(no2_det*1000.,P_det / 100.,'r-',label='$[NO_2] \ | \ P_{det}$')
        plt.axhline(y=P_TM5[trop_lev]/100.,c='b',label=r'$P_{tropopause,TM5}$')
        plt.ylim(50,1050)
        plt.gca().invert_yaxis()
        plt.xlabel(u'[NO$_2$] [ppb]')
        plt.ylabel('Pressure [hPa]')
        plt.title(r'[NO$_2$] profile')
        plt.legend(loc='upper right')
        plt.grid(True)
        plt.show()
    
    #Initialize required constants
    Ma = 0.02896
    N_A = 6.022e23
    A = 27000.*27000.
    g = 9.81
    t = 18
    
    #Calculate sub-columns - WRF
    if debug_lev >= 50.:
        Ps = np.zeros((wrf_Pres.shape[0]+1))
        Ps[0] = wrf_Psfc
        for i in range(1,Ps.shape[0]):
            Ps[i] = wrf_Pres[i-1] - Ps[i-1] + wrf_Pres[i-1]

        dP = Ps[0:-1] - Ps[1:]

        m_c = dP * A / g
        n_a = m_c / Ma

        psc_no2_wrf = n_a * 1.e-6 * wrf_no2 * N_A / (A * 1.e4)

    #Calculate sub-columns at interpolated pressure
    #dP_det = P_det - np.append(P_det[1:],P_TM5[trop_lev])
    dP_det = P_det - np.append(P_det[1:],wrf_Pres[-1])
    
    m_c_det = dP_det * A / g
    n_a_det = m_c_det / Ma
    
    psc_no2_det = n_a_det * 1.e-6 * no2_det * N_A / (A * 1.e4)
    
    #Sum sub-columns 
    psc_no2_TM5 = np.zeros((P_TM5[0:trop_lev+1].shape[0]))
    for i in range(0,psc_no2_TM5.shape[0]):
        psc_no2_TM5[i] = np.sum(psc_no2_det[num_points*i:num_points*i+num_points])
    
    #Check mass conservation
    if debug_lev >= 50.:
        WRF_sum_total = np.sum(psc_no2_wrf[i] for i in range(psc_no2_wrf.shape[0]) if wrf_Pres[i] >= P_TM5[trop_lev])
        TM5_sum_total = np.sum(psc_no2_wrf[0:trop_lev+1])
        print 'Mass conversion check (1 means equal NO2 column): ' + str(WRF_sum_total / TM5_sum_total)
        WRF_sum_lower = np.sum(psc_no2_wrf[i] for i in range(psc_no2_wrf.shape[0]) if wrf_Pres[i] > 88500.)
        if np.nanmax(P_TM5) > 88500.:
            TM5_sum_lower = np.sum(psc_no2_TM5[i] for i in range(psc_no2_TM5.shape[0]) if P_TM5[i] > 88500.)
            #det_sum_lower = np.sum(psc_no2_det[i] for i in range(psc_no2_det.shape[0]) if P_det[i] > 90000.)
            print 'Mass conversion check lowermost 115hPa: ' + str(WRF_sum_lower / TM5_sum_lower)
    
    #Plot VCD profiles
    if debug_lev >= 100.:
        plt.figure()
        plt.plot(psc_no2_wrf,wrf_Pres / 100.,'kx-',label=r'$\vec{x}_{WRF-Chem}$')
        plt.plot(psc_no2_det,P_det / 100.,'r-',label=r'$\vec{x}_{WRF-Chem \ | \ P_{det} }$')
        plt.plot(psc_no2_TM5,P_TM5[0:trop_lev+1] / 100.,'r.--',label=r'$\vec{x}_{WRF-Chem \ | \ P_{TM5} }$')
        plt.axhline(y=P_TM5[trop_lev]/100.,c='b',label=r'$P_{tropopause,TM5}$')
        plt.ylim(50,1050)
        plt.gca().invert_yaxis()
        plt.xlabel(u'NO$_2$ VCD [molec. cm$^{-2}$]')
        plt.ylabel('Pressure [hPa]')
        plt.title(r'NO$_2$ VCD profile')
        plt.legend(loc='upper right')
        plt.grid(True)
        plt.show()
        
    #Plot cumulative VCD profiles
    if debug_lev >= 100.:
        cum_psc_no2_wrf = np.cumsum(psc_no2_wrf)
        cum_psc_no2_TM5 = np.cumsum(psc_no2_TM5)
        cum_psc_no2_det = np.cumsum(psc_no2_det)
        
        print cum_psc_no2_wrf[-1], cum_psc_no2_TM5[trop_lev]
        
        plt.figure()
        plt.plot(cum_psc_no2_wrf,wrf_Pres / 100.,'kx-',label=r'$\vec{x}_{WRF-Chem}$')
        plt.plot(cum_psc_no2_TM5,P_TM5[0:trop_lev+1] / 100.,'r.--',label=r'$\vec{x}_{WRF-Chem \ | \ P_{TM5} }$')
        plt.plot(cum_psc_no2_det,P_det / 100.,'b--',label=r'$\vec{x}_{WRF-Chem \ | \ P_{det} }$')
        plt.axhline(y=P_TM5[trop_lev]/100.,c='b',label=r'$P_{tropopause,TM5}$')
        plt.ylim(50,1050)
        plt.gca().invert_yaxis()
        plt.xlabel(u'NO$_2$ VCD [molec. cm$^{-2}$]')
        plt.ylabel('Pressure [hPa]')
        plt.title(r'Cumulative NO$_2$ VCD profile')
        plt.legend(loc='upper left')
        plt.grid(True)
        plt.show()
    
    if debug_lev >= 100.:
        fig,axes=plt.subplots(1,2,figsize=(16,6))
        axes[0].plot(psc_no2_wrf/1.e15,wrf_Pres / 100.,'k.-',label=r'$\vec{x}_{WRF-Chem}$')
        #axes[0].plot(psc_no2_det,P_det / 100.,'r-',label=r'$\vec{x}_{WRF-Chem \ | \ P_{det} }$')
        axes[0].plot(psc_no2_TM5/1.e15,P_TM5[0:trop_lev+1] / 100.,'r.-',label=r'$\vec{x}_{WRF-Chem \ | \ P_{TM5} }$')
        #axes[0].axhline(y=P_TM5[trop_lev]/100.,c='b',label=r'$P_{tropopause,TM5}$')
        axes[0].set_ylim(50,1050)
        axes[0].invert_yaxis()
        axes[0].set_xlabel(u'$\hat{x}_{l,NO_2}$ [10$^{15}$ molec. cm$^{-2}$]',size=16)
        axes[0].set_ylabel('Pressure [hPa]',size=16)
        axes[0].set_title(r'NO$_2$ VCD profile',size=16)
        axes[0].legend(loc='upper right')
        axes[0].annotate("a)", xy=(0.92, 0.05), xycoords="axes fraction",weight='bold')
        axes[0].grid(True)
        
        axes[1].plot(cum_psc_no2_wrf/1.e15,wrf_Pres / 100.,'k.-',label=r'$\vec{x}_{WRF-Chem}$')
        axes[1].plot(cum_psc_no2_TM5/1.e15,P_TM5[0:trop_lev+1] / 100.,'r.-',label=r'$\vec{x}_{WRF-Chem \ | \ P_{TM5} }$')
        #axes[1].plot(cum_psc_no2_det,P_det / 100.,'b--',label=r'$\vec{x}_{WRF-Chem \ | \ P_{det} }$')
        #axes[1].axhline(y=P_TM5[trop_lev]/100.,c='b',label=r'$P_{tropopause,TM5}$')
        axes[1].set_ylim(50,1050)
        axes[1].invert_yaxis()
        axes[1].set_xlabel(u'$\hat{x}_{l,NO_2}$ [10$^{15}$ molec. cm$^{-2}$]',size=16)
        axes[1].set_ylabel('Pressure [hPa]',size=16)
        axes[1].set_title(r'Cumulative NO$_2$ VCD profile',size=16)
        axes[1].legend(loc='upper left')
        axes[1].annotate("b)", xy=(0.92, 0.05), xycoords="axes fraction",weight='bold')
        axes[1].grid(True)
        plt.show()
    
    #Re-calculate AMF according to Appendix D in Boersma et al. (2016)
    #MAQ_AV_281117+ added tropospheric AK calculation according to PSD
    AK_trop = AK * (AMF_total / AMF_trop)
    
    AMF_trop_adj = AMF_trop * np.sum(AK_trop * psc_no2_TM5) / np.sum(psc_no2_TM5)
    if debug_lev >= 50.:
        print 'Initial AMF = ' + str(AMF_trop) + ', Adjusted AMF = ' + str(AMF_trop_adj)
    
    #psc_no2_TM5[0:4] *= 1.2
    #AMF_trop_sens = AMF_trop * np.sum(AK * psc_no2_TM5) / np.sum(psc_no2_TM5)
    #if debug_lev >= 100.:
    #    print 'Sensitivity of PBL NO2, re-calculated AMF = ' + str(AMF_trop_sens)
    
    #Re-calculate tropospheric VCD
    SCD_strat = VCD_strat * AMF_strat
    VCD_trop_adj = (SCD - SCD_strat) / AMF_trop_adj
    VCD_trop_rec = (SCD - SCD_strat) / AMF_trop
    if debug_lev >= 50.:
        print 'Initial VCD_trop = ' + str(VCD_trop) + ', Adjusted VCD_trop = ' + str(VCD_trop_adj)
        print 'Initial VCD_trop = ' + str(VCD_trop) + ', Re-calculated VCD_trop = ' + str(VCD_trop_rec)
    
    #return AMF_trop, AMF_trop_adj, VCD_trop, VCD_trop_adj, VCD_trop_rec
    return AMF_trop_adj,VCD_trop_adj
