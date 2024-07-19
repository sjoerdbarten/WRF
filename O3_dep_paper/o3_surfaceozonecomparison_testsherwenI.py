import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import datetime as dt
from sklearn.metrics import mean_squared_error
from math import sqrt
from scipy.stats import gaussian_kde
from matplotlib.lines import Line2D

def mean_error(set_1,set_2):
	me = np.nanmean(set_1-set_2)
	return me

def mean_abs_error(set_1,set_2):
	mae = np.nanmean(np.abs(set_1-set_2))
        return mae

def root_mean_squared_error(set_1,set_2):
	rms = sqrt(mean_squared_error(set_1, set_2))
	return rms

def loaddata(o3datafile,wrfdatafile,dailyaverage,histplot,statistics):
	stationdata = pd.read_csv(o3datafile)
	stationdata.index = stationdata['Unnamed: 0']
	stationdata.index = pd.to_datetime(stationdata.index)
	stationdata.drop(['Unnamed: 0'],axis=1)
	bf1 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-10_00:00:00','r')
	bf2 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-13_01:00:00','r')
	bf3 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-21_01:00:00','r')
	bf4 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-31_01:00:00','r')
	cf1 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_coaregv2_d01_2008-08-10_00:00:00','r')
	cf2 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_coaregv2_d01_2008-08-15_01:00:00','r')
	df1 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chem_nudgedBL_coareg_d01_2008-08-10_00:00:00','r')
	df2 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chem_nudgedBL_coareg_d01_2008-08-11_01:00:00','r')
	df3 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chem_nudgedBL_coareg_d01_2008-08-14_01:00:00','r')
	wrfdata_fixeddep = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_fixeddep_d01_2008-08-10_00:00:00','r')
	wrfdata_nudgedBL_fixeddep = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_nudgedBL_fixeddep_d01_2008-08-10_00:00:00','r')
	wrfdata_nudgedBL_nofixeddep = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chem_nudgedBL_nofixeddep_d01_2008-08-10_00:00:00','r')
	sh1 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chem_nudgedBL_coareg_sherweniodide_d01_2008-08-10_00:00:00','r')
	sh2 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chem_nudgedBL_coareg_sherweniodide_d01_2008-08-11_01:00:00','r')
		
	#baserun = standard wrf, run1 = adjusted deposition to snow and ice, run2 = coareg routine I=f(T)
	#wrfo3data_base = wrfdata['o3'][:,0,:,:]*1000. #surface ozone to ppb
	wrfo3data_base = np.concatenate((bf1['o3'][0:73,0,:,:]*1000.,bf2['o3'][0:192,0,:,:]*1000.,bf3['o3'][0:240,0,:,:]*1000.,bf4['o3'][:,0,:,:]*1000.),axis=0)
	wrfo3data_run1 = wrfdata_fixeddep['o3'][:,0,:,:]*1000. #surface ozone to ppb
	wrfo3data_run2 = np.concatenate((cf1['o3'][0:121,0,:,:]*1000.,cf2['o3'][:,0,:,:]*1000.),axis=0) #surface ozone to ppb
	wrfo3data_run3 = wrfdata_nudgedBL_nofixeddep['o3'][:,0,:,:]*1000. #surface ozone to ppb
	#wrfo3data_run4 = wrfdata_nudgedBL_fixeddep['o3'][:,0,:,:]*1000. #surface ozone to ppb
	wrfo3data_run5 = np.concatenate((df1['o3'][0:24,0,:,:]*1000.,df2['o3'][0:73,0,:,:]*1000.,df3['o3'][:,0,:,:]*1000.),axis=0)
	wrfo3data_run4 = np.concatenate((sh1['o3'][0:25,0,:,:]*1000.,sh2['o3'][:,0,:,:]*1000.),axis=0)
		
	#subsetting data based on wrf output
	wrfstarttime = dt.datetime(2008,8,10,00,00,00)
	maxtimestep = wrfo3data_base.shape[0]
	wrfo3data_base = wrfo3data_base[0:maxtimestep]
	wrfo3data_run1 = wrfo3data_run1[0:maxtimestep]
	wrfo3data_run2 = wrfo3data_run2[0:maxtimestep]
	wrfo3data_run3 = wrfo3data_run3[0:maxtimestep]	
	wrfo3data_run4 = wrfo3data_run4[0:maxtimestep]	
	wrfo3data_run5 = wrfo3data_run5[0:maxtimestep]	
	timearr = np.array([wrfstarttime + dt.timedelta(hours=i) for i in xrange(maxtimestep)])
	stationdata = stationdata[timearr[0]:timearr[-1]]
	stationdatacams = stationdata[::6]
	
	#import CAMS data
	camsdata = np.load('/home/WUR/barte035/WRFChem/o3_analysis_DATA/CAMS-MACC/MACCnpo3array.npy')
	camstime = np.load('/home/WUR/barte035/WRFChem/o3_analysis_DATA/CAMS-MACC/MACCnptimearray.npy')
	
	#get wrf locations
	loc_barrow = [188,96] #[71.3230,-156.6114]
    	loc_storhofdi = [30,89] #[63.400,-20.288]
    	loc_summit = [74,84] #[72.5800018311,-38.4799995422]
    	loc_ahtari = [30,166] #[62.583333,24.183333]
    	loc_bredkalen = [29,150] #[63.85,15.333333]
    	loc_esrange = [46,154] #[67.883333,21.066667]
    	loc_karasjok = [54,157] #[69.466667,25.216667]
    	loc_karvatn = [22,140] #[62.783333,8.883333]
    	loc_lerwick = [11,122] #[60.13922,-1.185319]
    	loc_oulanka = [46,168] #[66.320278,29.401667]
    	loc_pallas = [48,158] #[67.973333333,24.116111111]
    	loc_tustervatn = [36,146] #[65.833333,13.916667]
    	loc_villum = [93,115] #[81.6,-16.67]
    	loc_vindeln = [32,157] #[64.25,19.766667]
    	loc_virolahti = [25,176] #[60.526667,27.686111]
    	loc_zeppelin = [85,132] #[78.90715,11.88668]    
    	loc_ascos = [114,123] #[87.4,-6.0] #rough, NEEDS TIME ADAPTIVE LOCATION
	loc_whitehorse = [202,46] #[60.718609,-135.049193]
    	loc_yellowknife = [167,29] #[62.45207,-114.364]
    	loc_normanwells = [180,49] #[65.27926,-126.813]
    	loc_fortliard = [186,30] #[60.23583,-123.467]
    	loc_inuvik = [180,65] #[68.36005,-133.727]
    	loc_denalinp = [209,73] #[63.72,-148.97]
    	loc_alert = [111,99] #[82.4991455078,-62.3415260315]
	loc_hurdal = [14,146] #[60.372386,11.078142]
	
	#baserun vars
	wrf_barrow_base = wrfo3data_base[:,loc_barrow[0],loc_barrow[1]]
	wrf_storhofdi_base = wrfo3data_base[:,loc_storhofdi[0],loc_storhofdi[1]]
	wrf_summit_base = wrfo3data_base[:,loc_summit[0],loc_summit[1]]
	wrf_ahtari_base = wrfo3data_base[:,loc_ahtari[0],loc_ahtari[1]]
	wrf_bredkalen_base = wrfo3data_base[:,loc_bredkalen[0],loc_bredkalen[1]]
	wrf_esrange_base = wrfo3data_base[:,loc_esrange[0],loc_esrange[1]]
	wrf_karasjok_base = wrfo3data_base[:,loc_karasjok[0],loc_karasjok[1]]
	wrf_karvatn_base = wrfo3data_base[:,loc_karvatn[0],loc_karvatn[1]]
	wrf_lerwick_base = wrfo3data_base[:,loc_lerwick[0],loc_lerwick[1]]
	wrf_oulanka_base = wrfo3data_base[:,loc_oulanka[0],loc_oulanka[1]]
	wrf_pallas_base = wrfo3data_base[:,loc_pallas[0],loc_pallas[1]]
	wrf_tustervatn_base = wrfo3data_base[:,loc_tustervatn[0],loc_tustervatn[1]]
	wrf_villum_base = wrfo3data_base[:,loc_villum[0],loc_villum[1]]
	wrf_vindeln_base = wrfo3data_base[:,loc_vindeln[0],loc_vindeln[1]]
	wrf_virolahti_base = wrfo3data_base[:,loc_virolahti[0],loc_virolahti[1]]
	wrf_zeppelin_base = wrfo3data_base[:,loc_zeppelin[0],loc_zeppelin[1]]
	wrf_ascos_base = wrfo3data_base[:,loc_ascos[0],loc_ascos[1]]
	wrf_whitehorse_base = wrfo3data_base[:,loc_whitehorse[0],loc_whitehorse[1]]
	wrf_yellowknife_base = wrfo3data_base[:,loc_yellowknife[0],loc_yellowknife[1]]
	wrf_normanwells_base = wrfo3data_base[:,loc_normanwells[0],loc_normanwells[1]]
	wrf_fortliard_base = wrfo3data_base[:,loc_fortliard[0],loc_fortliard[1]]
	wrf_inuvik_base = wrfo3data_base[:,loc_inuvik[0],loc_inuvik[1]]
	wrf_denalinp_base = wrfo3data_base[:,loc_denalinp[0],loc_denalinp[1]]
	wrf_alert_base = wrfo3data_base[:,loc_alert[0],loc_alert[1]]
	wrf_hurdal_base = wrfo3data_base[:,loc_hurdal[0],loc_hurdal[1]]
	
	#run1 vars
	wrf_barrow_run1 = wrfo3data_run1[:,loc_barrow[0],loc_barrow[1]]
	wrf_storhofdi_run1 = wrfo3data_run1[:,loc_storhofdi[0],loc_storhofdi[1]]
	wrf_summit_run1 = wrfo3data_run1[:,loc_summit[0],loc_summit[1]]
	wrf_ahtari_run1 = wrfo3data_run1[:,loc_ahtari[0],loc_ahtari[1]]
	wrf_bredkalen_run1 = wrfo3data_run1[:,loc_bredkalen[0],loc_bredkalen[1]]
	wrf_esrange_run1 = wrfo3data_run1[:,loc_esrange[0],loc_esrange[1]]
	wrf_karasjok_run1 = wrfo3data_run1[:,loc_karasjok[0],loc_karasjok[1]]
	wrf_karvatn_run1 = wrfo3data_run1[:,loc_karvatn[0],loc_karvatn[1]]
	wrf_lerwick_run1 = wrfo3data_run1[:,loc_lerwick[0],loc_lerwick[1]]
	wrf_oulanka_run1 = wrfo3data_run1[:,loc_oulanka[0],loc_oulanka[1]]
	wrf_pallas_run1 = wrfo3data_run1[:,loc_pallas[0],loc_pallas[1]]
	wrf_tustervatn_run1 = wrfo3data_run1[:,loc_tustervatn[0],loc_tustervatn[1]]
	wrf_villum_run1 = wrfo3data_run1[:,loc_villum[0],loc_villum[1]]
	wrf_vindeln_run1 = wrfo3data_run1[:,loc_vindeln[0],loc_vindeln[1]]
	wrf_virolahti_run1 = wrfo3data_run1[:,loc_virolahti[0],loc_virolahti[1]]
	wrf_zeppelin_run1 = wrfo3data_run1[:,loc_zeppelin[0],loc_zeppelin[1]]
	wrf_ascos_run1 = wrfo3data_run1[:,loc_ascos[0],loc_ascos[1]]
	wrf_whitehorse_run1 = wrfo3data_run1[:,loc_whitehorse[0],loc_whitehorse[1]]
	wrf_yellowknife_run1 = wrfo3data_run1[:,loc_yellowknife[0],loc_yellowknife[1]]
	wrf_normanwells_run1 = wrfo3data_run1[:,loc_normanwells[0],loc_normanwells[1]]
	wrf_fortliard_run1 = wrfo3data_run1[:,loc_fortliard[0],loc_fortliard[1]]
	wrf_inuvik_run1 = wrfo3data_run1[:,loc_inuvik[0],loc_inuvik[1]]
	wrf_denalinp_run1 = wrfo3data_run1[:,loc_denalinp[0],loc_denalinp[1]]
	wrf_alert_run1 = wrfo3data_run1[:,loc_alert[0],loc_alert[1]]
	wrf_hurdal_run1 = wrfo3data_run1[:,loc_hurdal[0],loc_hurdal[1]]

	#run2 vars
	wrf_barrow_run2 = wrfo3data_run2[:,loc_barrow[0],loc_barrow[1]]
	wrf_storhofdi_run2 = wrfo3data_run2[:,loc_storhofdi[0],loc_storhofdi[1]]
	wrf_summit_run2 = wrfo3data_run2[:,loc_summit[0],loc_summit[1]]
	wrf_ahtari_run2 = wrfo3data_run2[:,loc_ahtari[0],loc_ahtari[1]]
	wrf_bredkalen_run2 = wrfo3data_run2[:,loc_bredkalen[0],loc_bredkalen[1]]
	wrf_esrange_run2 = wrfo3data_run2[:,loc_esrange[0],loc_esrange[1]]
	wrf_karasjok_run2 = wrfo3data_run2[:,loc_karasjok[0],loc_karasjok[1]]
	wrf_karvatn_run2 = wrfo3data_run2[:,loc_karvatn[0],loc_karvatn[1]]
	wrf_lerwick_run2 = wrfo3data_run2[:,loc_lerwick[0],loc_lerwick[1]]
	wrf_oulanka_run2 = wrfo3data_run2[:,loc_oulanka[0],loc_oulanka[1]]
	wrf_pallas_run2 = wrfo3data_run2[:,loc_pallas[0],loc_pallas[1]]
	wrf_tustervatn_run2 = wrfo3data_run2[:,loc_tustervatn[0],loc_tustervatn[1]]
	wrf_villum_run2 = wrfo3data_run2[:,loc_villum[0],loc_villum[1]]
	wrf_vindeln_run2 = wrfo3data_run2[:,loc_vindeln[0],loc_vindeln[1]]
	wrf_virolahti_run2 = wrfo3data_run2[:,loc_virolahti[0],loc_virolahti[1]]
	wrf_zeppelin_run2 = wrfo3data_run2[:,loc_zeppelin[0],loc_zeppelin[1]]
	wrf_ascos_run2 = wrfo3data_run2[:,loc_ascos[0],loc_ascos[1]]
	wrf_whitehorse_run2 = wrfo3data_run2[:,loc_whitehorse[0],loc_whitehorse[1]]
	wrf_yellowknife_run2 = wrfo3data_run2[:,loc_yellowknife[0],loc_yellowknife[1]]
	wrf_normanwells_run2 = wrfo3data_run2[:,loc_normanwells[0],loc_normanwells[1]]
	wrf_fortliard_run2 = wrfo3data_run2[:,loc_fortliard[0],loc_fortliard[1]]
	wrf_inuvik_run2 = wrfo3data_run2[:,loc_inuvik[0],loc_inuvik[1]]
	wrf_denalinp_run2 = wrfo3data_run2[:,loc_denalinp[0],loc_denalinp[1]]
	wrf_alert_run2 = wrfo3data_run2[:,loc_alert[0],loc_alert[1]]
	wrf_hurdal_run2 = wrfo3data_run2[:,loc_hurdal[0],loc_hurdal[1]]

	#run3 vars
	wrf_barrow_run3 = wrfo3data_run3[:,loc_barrow[0],loc_barrow[1]]
	wrf_storhofdi_run3 = wrfo3data_run3[:,loc_storhofdi[0],loc_storhofdi[1]]
	wrf_summit_run3 = wrfo3data_run3[:,loc_summit[0],loc_summit[1]]
	wrf_ahtari_run3 = wrfo3data_run3[:,loc_ahtari[0],loc_ahtari[1]]
	wrf_bredkalen_run3 = wrfo3data_run3[:,loc_bredkalen[0],loc_bredkalen[1]]
	wrf_esrange_run3 = wrfo3data_run3[:,loc_esrange[0],loc_esrange[1]]
	wrf_karasjok_run3 = wrfo3data_run3[:,loc_karasjok[0],loc_karasjok[1]]
	wrf_karvatn_run3 = wrfo3data_run3[:,loc_karvatn[0],loc_karvatn[1]]
	wrf_lerwick_run3 = wrfo3data_run3[:,loc_lerwick[0],loc_lerwick[1]]
	wrf_oulanka_run3 = wrfo3data_run3[:,loc_oulanka[0],loc_oulanka[1]]
	wrf_pallas_run3 = wrfo3data_run3[:,loc_pallas[0],loc_pallas[1]]
	wrf_tustervatn_run3 = wrfo3data_run3[:,loc_tustervatn[0],loc_tustervatn[1]]
	wrf_villum_run3 = wrfo3data_run3[:,loc_villum[0],loc_villum[1]]
	wrf_vindeln_run3 = wrfo3data_run3[:,loc_vindeln[0],loc_vindeln[1]]
	wrf_virolahti_run3 = wrfo3data_run3[:,loc_virolahti[0],loc_virolahti[1]]
	wrf_zeppelin_run3 = wrfo3data_run3[:,loc_zeppelin[0],loc_zeppelin[1]]
	wrf_ascos_run3 = wrfo3data_run3[:,loc_ascos[0],loc_ascos[1]]
	wrf_whitehorse_run3 = wrfo3data_run3[:,loc_whitehorse[0],loc_whitehorse[1]]
	wrf_yellowknife_run3 = wrfo3data_run3[:,loc_yellowknife[0],loc_yellowknife[1]]
	wrf_normanwells_run3 = wrfo3data_run3[:,loc_normanwells[0],loc_normanwells[1]]
	wrf_fortliard_run3 = wrfo3data_run3[:,loc_fortliard[0],loc_fortliard[1]]
	wrf_inuvik_run3 = wrfo3data_run3[:,loc_inuvik[0],loc_inuvik[1]]
	wrf_denalinp_run3 = wrfo3data_run3[:,loc_denalinp[0],loc_denalinp[1]]
	wrf_alert_run3 = wrfo3data_run3[:,loc_alert[0],loc_alert[1]]
	wrf_hurdal_run3 = wrfo3data_run3[:,loc_hurdal[0],loc_hurdal[1]]

	#run4 vars
	wrf_barrow_run4 = wrfo3data_run4[:,loc_barrow[0],loc_barrow[1]]
	wrf_storhofdi_run4 = wrfo3data_run4[:,loc_storhofdi[0],loc_storhofdi[1]]
	wrf_summit_run4 = wrfo3data_run4[:,loc_summit[0],loc_summit[1]]
	wrf_ahtari_run4 = wrfo3data_run4[:,loc_ahtari[0],loc_ahtari[1]]
	wrf_bredkalen_run4 = wrfo3data_run4[:,loc_bredkalen[0],loc_bredkalen[1]]
	wrf_esrange_run4 = wrfo3data_run4[:,loc_esrange[0],loc_esrange[1]]
	wrf_karasjok_run4 = wrfo3data_run4[:,loc_karasjok[0],loc_karasjok[1]]
	wrf_karvatn_run4 = wrfo3data_run4[:,loc_karvatn[0],loc_karvatn[1]]
	wrf_lerwick_run4 = wrfo3data_run4[:,loc_lerwick[0],loc_lerwick[1]]
	wrf_oulanka_run4 = wrfo3data_run4[:,loc_oulanka[0],loc_oulanka[1]]
	wrf_pallas_run4 = wrfo3data_run4[:,loc_pallas[0],loc_pallas[1]]
	wrf_tustervatn_run4 = wrfo3data_run4[:,loc_tustervatn[0],loc_tustervatn[1]]
	wrf_villum_run4 = wrfo3data_run4[:,loc_villum[0],loc_villum[1]]
	wrf_vindeln_run4 = wrfo3data_run4[:,loc_vindeln[0],loc_vindeln[1]]
	wrf_virolahti_run4 = wrfo3data_run4[:,loc_virolahti[0],loc_virolahti[1]]
	wrf_zeppelin_run4 = wrfo3data_run4[:,loc_zeppelin[0],loc_zeppelin[1]]
	wrf_ascos_run4 = wrfo3data_run4[:,loc_ascos[0],loc_ascos[1]]
	wrf_whitehorse_run4 = wrfo3data_run4[:,loc_whitehorse[0],loc_whitehorse[1]]
	wrf_yellowknife_run4 = wrfo3data_run4[:,loc_yellowknife[0],loc_yellowknife[1]]
	wrf_normanwells_run4 = wrfo3data_run4[:,loc_normanwells[0],loc_normanwells[1]]
	wrf_fortliard_run4 = wrfo3data_run4[:,loc_fortliard[0],loc_fortliard[1]]
	wrf_inuvik_run4 = wrfo3data_run4[:,loc_inuvik[0],loc_inuvik[1]]
	wrf_denalinp_run4 = wrfo3data_run4[:,loc_denalinp[0],loc_denalinp[1]]
	wrf_alert_run4 = wrfo3data_run4[:,loc_alert[0],loc_alert[1]]
	wrf_hurdal_run4 = wrfo3data_run4[:,loc_hurdal[0],loc_hurdal[1]]

	#run5 vars
	wrf_barrow_run5 = wrfo3data_run5[:,loc_barrow[0],loc_barrow[1]]
	wrf_storhofdi_run5 = wrfo3data_run5[:,loc_storhofdi[0],loc_storhofdi[1]]
	wrf_summit_run5 = wrfo3data_run5[:,loc_summit[0],loc_summit[1]]
	wrf_ahtari_run5 = wrfo3data_run5[:,loc_ahtari[0],loc_ahtari[1]]
	wrf_bredkalen_run5 = wrfo3data_run5[:,loc_bredkalen[0],loc_bredkalen[1]]
	wrf_esrange_run5 = wrfo3data_run5[:,loc_esrange[0],loc_esrange[1]]
	wrf_karasjok_run5 = wrfo3data_run5[:,loc_karasjok[0],loc_karasjok[1]]
	wrf_karvatn_run5 = wrfo3data_run5[:,loc_karvatn[0],loc_karvatn[1]]
	wrf_lerwick_run5 = wrfo3data_run5[:,loc_lerwick[0],loc_lerwick[1]]
	wrf_oulanka_run5 = wrfo3data_run5[:,loc_oulanka[0],loc_oulanka[1]]
	wrf_pallas_run5 = wrfo3data_run5[:,loc_pallas[0],loc_pallas[1]]
	wrf_tustervatn_run5 = wrfo3data_run5[:,loc_tustervatn[0],loc_tustervatn[1]]
	wrf_villum_run5 = wrfo3data_run5[:,loc_villum[0],loc_villum[1]]
	wrf_vindeln_run5 = wrfo3data_run5[:,loc_vindeln[0],loc_vindeln[1]]
	wrf_virolahti_run5 = wrfo3data_run5[:,loc_virolahti[0],loc_virolahti[1]]
	wrf_zeppelin_run5 = wrfo3data_run5[:,loc_zeppelin[0],loc_zeppelin[1]]
	wrf_ascos_run5 = wrfo3data_run5[:,loc_ascos[0],loc_ascos[1]]
	wrf_whitehorse_run5 = wrfo3data_run5[:,loc_whitehorse[0],loc_whitehorse[1]]
	wrf_yellowknife_run5 = wrfo3data_run5[:,loc_yellowknife[0],loc_yellowknife[1]]
	wrf_normanwells_run5 = wrfo3data_run5[:,loc_normanwells[0],loc_normanwells[1]]
	wrf_fortliard_run5 = wrfo3data_run5[:,loc_fortliard[0],loc_fortliard[1]]
	wrf_inuvik_run5 = wrfo3data_run5[:,loc_inuvik[0],loc_inuvik[1]]
	wrf_denalinp_run5 = wrfo3data_run5[:,loc_denalinp[0],loc_denalinp[1]]
	wrf_alert_run5 = wrfo3data_run5[:,loc_alert[0],loc_alert[1]]
	wrf_hurdal_run5 = wrfo3data_run5[:,loc_hurdal[0],loc_hurdal[1]]
	
	#CAMS data
	cams_barrow = camsdata[:,loc_barrow[0],loc_barrow[1]]
	cams_storhofdi = camsdata[:,loc_storhofdi[0],loc_storhofdi[1]]
	cams_summit = camsdata[:,loc_summit[0],loc_summit[1]]
	cams_ahtari = camsdata[:,loc_ahtari[0],loc_ahtari[1]]
	cams_bredkalen = camsdata[:,loc_bredkalen[0],loc_bredkalen[1]]
	cams_esrange = camsdata[:,loc_esrange[0],loc_esrange[1]]
	cams_karasjok = camsdata[:,loc_karasjok[0],loc_karasjok[1]]
	cams_karvatn = camsdata[:,loc_karvatn[0],loc_karvatn[1]]
	cams_lerwick = camsdata[:,loc_lerwick[0],loc_lerwick[1]]
	cams_oulanka = camsdata[:,loc_oulanka[0],loc_oulanka[1]]
	cams_pallas = camsdata[:,loc_pallas[0],loc_pallas[1]]
	cams_tustervatn = camsdata[:,loc_tustervatn[0],loc_tustervatn[1]]
	cams_villum = camsdata[:,loc_villum[0],loc_villum[1]]
	cams_vindeln = camsdata[:,loc_vindeln[0],loc_vindeln[1]]
	cams_virolahti = camsdata[:,loc_virolahti[0],loc_virolahti[1]]
	cams_zeppelin = camsdata[:,loc_zeppelin[0],loc_zeppelin[1]]
	cams_ascos = camsdata[:,loc_ascos[0],loc_ascos[1]]
	cams_whitehorse = camsdata[:,loc_whitehorse[0],loc_whitehorse[1]]
	cams_yellowknife = camsdata[:,loc_yellowknife[0],loc_yellowknife[1]]
	cams_normanwells = camsdata[:,loc_normanwells[0],loc_normanwells[1]]
	cams_fortliard = camsdata[:,loc_fortliard[0],loc_fortliard[1]]
	cams_inuvik = camsdata[:,loc_inuvik[0],loc_inuvik[1]]
	cams_denalinp = camsdata[:,loc_denalinp[0],loc_denalinp[1]]
	cams_alert = camsdata[:,loc_alert[0],loc_alert[1]]
	cams_hurdal = camsdata[:,loc_hurdal[0],loc_hurdal[1]]

	loclist = ['barrow','storhofdi','summit','ahtari','bredkalen','esrange','karasjok','karvatn','lerwick','oulanka','pallas','tustervatn','villum','vindeln','virolahti','zeppelin','ascos','whitehorse','yellowknife','normanwells','fortliard','inuvik','denalinp','alert','hurdal']
	plotlist_base = [wrf_barrow_base,wrf_storhofdi_base,wrf_summit_base,wrf_ahtari_base,wrf_bredkalen_base,wrf_esrange_base,wrf_karasjok_base,wrf_karvatn_base,wrf_lerwick_base,wrf_oulanka_base,wrf_pallas_base,wrf_tustervatn_base,wrf_villum_base,wrf_vindeln_base,wrf_virolahti_base,wrf_zeppelin_base,wrf_ascos_base,wrf_whitehorse_base,wrf_yellowknife_base,wrf_normanwells_base,wrf_fortliard_base,wrf_inuvik_base,wrf_denalinp_base,wrf_alert_base,wrf_hurdal_base]
	plotlist_run1 = [wrf_barrow_run1,wrf_storhofdi_run1,wrf_summit_run1,wrf_ahtari_run1,wrf_bredkalen_run1,wrf_esrange_run1,wrf_karasjok_run1,wrf_karvatn_run1,wrf_lerwick_run1,wrf_oulanka_run1,wrf_pallas_run1,wrf_tustervatn_run1,wrf_villum_run1,wrf_vindeln_run1,wrf_virolahti_run1,wrf_zeppelin_run1,wrf_ascos_run1,wrf_whitehorse_run1,wrf_yellowknife_run1,wrf_normanwells_run1,wrf_fortliard_run1,wrf_inuvik_run1,wrf_denalinp_run1,wrf_alert_run1,wrf_hurdal_run1]
	plotlist_run2 = [wrf_barrow_run2,wrf_storhofdi_run2,wrf_summit_run2,wrf_ahtari_run2,wrf_bredkalen_run2,wrf_esrange_run2,wrf_karasjok_run2,wrf_karvatn_run2,wrf_lerwick_run2,wrf_oulanka_run2,wrf_pallas_run2,wrf_tustervatn_run2,wrf_villum_run2,wrf_vindeln_run2,wrf_virolahti_run2,wrf_zeppelin_run2,wrf_ascos_run2,wrf_whitehorse_run2,wrf_yellowknife_run2,wrf_normanwells_run2,wrf_fortliard_run2,wrf_inuvik_run2,wrf_denalinp_run2,wrf_alert_run2,wrf_hurdal_run2]
	plotlist_run3 = [wrf_barrow_run3,wrf_storhofdi_run3,wrf_summit_run3,wrf_ahtari_run3,wrf_bredkalen_run3,wrf_esrange_run3,wrf_karasjok_run3,wrf_karvatn_run3,wrf_lerwick_run3,wrf_oulanka_run3,wrf_pallas_run3,wrf_tustervatn_run3,wrf_villum_run3,wrf_vindeln_run3,wrf_virolahti_run3,wrf_zeppelin_run3,wrf_ascos_run3,wrf_whitehorse_run3,wrf_yellowknife_run3,wrf_normanwells_run3,wrf_fortliard_run3,wrf_inuvik_run3,wrf_denalinp_run3,wrf_alert_run3,wrf_hurdal_run3]
	plotlist_run4 = [wrf_barrow_run4,wrf_storhofdi_run4,wrf_summit_run4,wrf_ahtari_run4,wrf_bredkalen_run4,wrf_esrange_run4,wrf_karasjok_run4,wrf_karvatn_run4,wrf_lerwick_run4,wrf_oulanka_run4,wrf_pallas_run4,wrf_tustervatn_run4,wrf_villum_run4,wrf_vindeln_run4,wrf_virolahti_run4,wrf_zeppelin_run4,wrf_ascos_run4,wrf_whitehorse_run4,wrf_yellowknife_run4,wrf_normanwells_run4,wrf_fortliard_run4,wrf_inuvik_run4,wrf_denalinp_run4,wrf_alert_run4,wrf_hurdal_run4]
	plotlist_run5 = [wrf_barrow_run5,wrf_storhofdi_run5,wrf_summit_run5,wrf_ahtari_run5,wrf_bredkalen_run5,wrf_esrange_run5,wrf_karasjok_run5,wrf_karvatn_run5,wrf_lerwick_run5,wrf_oulanka_run5,wrf_pallas_run5,wrf_tustervatn_run5,wrf_villum_run5,wrf_vindeln_run5,wrf_virolahti_run5,wrf_zeppelin_run5,wrf_ascos_run5,wrf_whitehorse_run5,wrf_yellowknife_run5,wrf_normanwells_run5,wrf_fortliard_run5,wrf_inuvik_run5,wrf_denalinp_run5,wrf_alert_run5,wrf_hurdal_run5]
	camslist = [cams_barrow,cams_storhofdi,cams_summit,cams_ahtari,cams_bredkalen,cams_esrange,cams_karasjok,cams_karvatn,cams_lerwick,cams_oulanka,cams_pallas,cams_tustervatn,cams_villum,cams_vindeln,cams_virolahti,cams_zeppelin,cams_ascos,cams_whitehorse,cams_yellowknife,cams_normanwells,cams_fortliard,cams_inuvik,cams_denalinp,cams_alert,cams_hurdal]
	
	if statistics == True:
		#HA = high arctic
		#RE = remote
		#DI = diurnal
		filters = ["HA","zRE","HA","zzDI","zzDI","zzDI","zzDI","zzDI","zRE","zzDI","zzDI","zzDI","HA","zzDI","zzDI","HA","HA","zzDI","zzDI","zzDI","zzDI","zRE","zzDI","HA","zzDI"]
		loclist_abb = ['BRW','ICE','SUM','AHT','BRE','ESR','KAS','KRV','SIS','OUX','PAL','TUV','VIL','VDI','VIR','ZEP','ASC','WHI','YEL','NOR','FOR','INU','DEN','ALT','HUR']
	
	bias_base = np.zeros(len(loclist))
	mae_base = np.zeros(len(loclist))
	#rmse_base = np.zeros(len(loclist))
	bias_fixdep = np.zeros(len(loclist))
	mae_fixdep = np.zeros(len(loclist))
	#rmse_fixdep = np.zeros(len(loclist))
	bias_coareg = np.zeros(len(loclist))
	mae_coareg = np.zeros(len(loclist))
	#rmse_coareg = np.zeros(len(loclist))
	bias_nudged = np.zeros(len(loclist))
	mae_nudged = np.zeros(len(loclist))
	#rmse_nudged = np.zeros(len(loclist))
	bias_nofixdep_nudged = np.zeros(len(loclist))
	mae_nofixdep_nudged = np.zeros(len(loclist))
	#rmse_nofixdep_nudged = np.zeros(len(loclist))
	bias_fixdep_nudged = np.zeros(len(loclist))
	mae_fixdep_nudged = np.zeros(len(loclist))
	#rmse_fixdep_nudged = np.zeros(len(loclist))
	bias_coareg_nudged = np.zeros(len(loclist))
	mae_coareg_nudged = np.zeros(len(loclist))
	#rmse_coareg_nudged = np.zeros(len(loclist))
	bias_cams = np.zeros(len(loclist))
	mae_cams = np.zeros(len(loclist))
	#rmse_cams = np.zeros(len(loclist))
	
	for i in range(len(loclist)):
		filterr = ~np.isnan(plotlist_base[i]-stationdata[loclist[i]+'_o3'])
		
		me = mean_error(plotlist_base[i],stationdata[loclist[i]+'_o3'])
		mae = mean_abs_error(plotlist_base[i],stationdata[loclist[i]+'_o3'])
		#rmse = root_mean_squared_error(plotlist_base[i][filterr][~np.isnan(plotlist_base[i][filterr])],stationdata[loclist[i]+'_o3'].values[filterr][~np.isnan(plotlist_base[i][filterr])])
		mer1 = mean_error(plotlist_run1[i],stationdata[loclist[i]+'_o3'])
		maer1 = mean_abs_error(plotlist_run1[i],stationdata[loclist[i]+'_o3'])
		#rmser1 = root_mean_squared_error(plotlist_run1[i][filterr][~np.isnan(plotlist_run1[i][filterr])],stationdata[loclist[i]+'_o3'].values[filterr][~np.isnan(plotlist_base[i][filterr])])
		mer2 = mean_error(plotlist_run2[i],stationdata[loclist[i]+'_o3'])
		maer2 = mean_abs_error(plotlist_run2[i],stationdata[loclist[i]+'_o3'])
		#rmser2 = root_mean_squared_error(plotlist_run2[i][filterr][~np.isnan(plotlist_run2[i][filterr])],stationdata[loclist[i]+'_o3'].values[filterr][~np.isnan(plotlist_base[i][filterr])])
		mer3 = mean_error(plotlist_run3[i],stationdata[loclist[i]+'_o3'])
		maer3 = mean_abs_error(plotlist_run3[i],stationdata[loclist[i]+'_o3'])
		#rmser3 = root_mean_squared_error(plotlist_run3[i][filterr][~np.isnan(plotlist_run3[i][filterr])],stationdata[loclist[i]+'_o3'].values[filterr][~np.isnan(plotlist_base[i][filterr])])
		mer4 = mean_error(plotlist_run4[i],stationdata[loclist[i]+'_o3'])
		maer4 = mean_abs_error(plotlist_run4[i],stationdata[loclist[i]+'_o3'])
		#rmser4 = root_mean_squared_error(plotlist_run4[i][filterr][~np.isnan(plotlist_run4[i][filterr])],stationdata[loclist[i]+'_o3'].values[filterr][~np.isnan(plotlist_base[i][filterr])])
		mer5 = mean_error(plotlist_run5[i],stationdata[loclist[i]+'_o3'])
		maer5 = mean_abs_error(plotlist_run5[i],stationdata[loclist[i]+'_o3'])
		#rmser5 = root_mean_squared_error(plotlist_run5[i][filterr][~np.isnan(plotlist_run5[i][filterr])],stationdata[loclist[i]+'_o3'].values[filterr][~np.isnan(plotlist_base[i][filterr])])
		mecams = mean_error(camslist[i],stationdatacams[loclist[i]+'_o3'])
		maecams = mean_abs_error(camslist[i],stationdatacams[loclist[i]+'_o3'])
		#rmsecams = root_mean_squared_error(camslist[i][filterr][~np.isnan(camslist[i][filterr])],stationdatacams[loclist[i]+'_o3'].values[filterr][~np.isnan(camslist[i][filterr])])
		
		bias_base[i] = me
		mae_base[i] = mae
		#rmse_base[i] = rmse
		bias_fixdep[i] = mer1
		mae_fixdep[i] = maer1
		#rmse_fixdep[i] = rmser1
		bias_coareg[i] = mer2
		mae_coareg[i] = maer2
		#rmse_coareg[i] = rmser2
		bias_nofixdep_nudged[i] = mer3
		mae_nofixdep_nudged[i] = maer3
		#rmse_nofixdep_nudged[i] = rmser3
		bias_fixdep_nudged[i] = mer4
		mae_fixdep_nudged[i] = maer4
		#rmse_fixdep_nudged[i] = rmser4
		bias_coareg_nudged[i] = mer5
		mae_coareg_nudged[i] = maer5
		#rmse_coareg_nudged[i] = rmser5
		bias_cams[i] = mecams
		mae_cams[i] = maecams
		#rmse_cams[i] = rmsecams

		plt.figure(figsize=(15,8))
		plt.plot(timearr,stationdata[loclist[i]+'_o3'],'ko',label='Observartions')
		plt.plot(timearr,plotlist_base[i],'r-',label='Base')
		#plt.plot(timearr,plotlist_run1[i],'g-',label='Run 1')
		#plt.plot(timearr,plotlist_run2[i],'b-',label='Run 2')
		plt.plot(timearr,plotlist_run3[i],color='orange',linestyle='-',label='Nudged')
		plt.plot(timearr,plotlist_run4[i],'m-',label='COAREGSherwen')
		plt.plot(timearr,plotlist_run5[i],'c-',label='NudgedCoareg')
		plt.plot(camstime,camslist[i],marker='x',linewidth=0,color='blue',markersize=8,label='MACC')
		plt.ylim([0,60])
		plt.ylabel('Surface ozone mixing ratio [ppb]')
		plt.title(loclist[i])
		if loclist[i] == 'Summit':
			plt.legend(loc=4)
		else:
			plt.legend(loc=1)
		plt.text(timearr[6],57,"Bias/MAE Base         = "+str(round(me,2))+"/"+str(round(mae,2))+" ppb")
		plt.text(timearr[6],55,"Bias/MAE Nudged       = "+str(round(mer3,2))+"/"+str(round(maer3,2))+" ppb")
		plt.text(timearr[6],53,"Bias/MAE NudgedCoareg = "+str(round(mer5,2))+"/"+str(round(maer5,2))+" ppb")
		plt.text(timearr[6],51,"Bias/MAE MACC         = "+str(round(mecams,2))+"/"+str(round(maecams,2))+" ppb")
		plt.text(timearr[6],51,"Bias/MAE COAREGSherwen  = "+str(round(mer4,2))+"/"+str(round(maer4,2))+" ppb")
		
		'''
		plt.text(timearr[6],55,'MAE base = '+str(round(mae,2))+' ppb')
		plt.text(timearr[6],57,'Bias base = '+str(round(me,2))+' ppb')		
		plt.text(timearr[6],51,'MAE run1 = '+str(round(maer1,2))+' ppb')
		plt.text(timearr[6],53,'Bias run1 = '+str(round(mer1,2))+' ppb')		
		plt.text(timearr[6],47,'MAE run2 = '+str(round(maer2,2))+' ppb')
		plt.text(timearr[6],49,'Bias run2 = '+str(round(mer2,2))+' ppb')		
		plt.text(timearr[5*24],55,'MAE run3 = '+str(round(maer3,2))+' ppb')
		plt.text(timearr[5*24],57,'Bias run3 = '+str(round(mer3,2))+' ppb')		
		plt.text(timearr[5*24],51,'MAE run4 = '+str(round(maer4,2))+' ppb')
		plt.text(timearr[5*24],53,'Bias run4 = '+str(round(mer4,2))+' ppb')		
		plt.text(timearr[5*24],47,'MAE run5 = '+str(round(maer5,2))+' ppb')
		plt.text(timearr[5*24],49,'Bias run5 = '+str(round(mer5,2))+' ppb')		
		plt.text(timearr[5*24],43,'MAE MACC = '+str(round(maecams,2))+' ppb')
		plt.text(timearr[5*24],45,'Bias MACC = '+str(round(mecams,2))+' ppb')
		'''
		
		plt.savefig('Figures/Surfaceplots/Timeseries/SherwenIvsMcDonaldI/'+loclist[i],dpi=300)
		#plt.show()
		plt.close()
		
		if statistics == True:
			#filter out nan for boxplot
			filterr = ~np.isnan(plotlist_base[i]-stationdata[loclist[i]+'_o3'])
			#data = [(plotlist_base[i]-stationdata[loclist[i]+'_o3'])[filterr],(plotlist_run1[i]-stationdata[loclist[i]+'_o3'])[filterr],(plotlist_run2[i]-stationdata[loclist[i]+'_o3'])[filterr],(plotlist_run3[i]-stationdata[loclist[i]+'_o3'])[filterr],(plotlist_run4[i]-stationdata[loclist[i]+'_o3'])[filterr],(plotlist_run5[i]-stationdata[loclist[i]+'_o3'])[filterr],(camslist[i]-stationdatacams[loclist[i]+'_o3'])[filterr]]
			data = [(plotlist_base[i]-stationdata[loclist[i]+'_o3'])[filterr],(plotlist_run3[i]-stationdata[loclist[i]+'_o3'])[filterr],(plotlist_run5[i]-stationdata[loclist[i]+'_o3'])[filterr],(camslist[i]-stationdatacams[loclist[i]+'_o3'])[filterr],(plotlist_run4[i]-stationdata[loclist[i]+'_o3'])[filterr]]
		
			fig, ax = plt.subplots()
			ax.boxplot(data)
			plt.hlines(y=0,xmin=0,xmax=len(data)+0.5,linestyles='--')
			#plt.vlines(x=3.5,ymin=-100,ymax=100,colors='k',linestyles='solid',linewidth=3)
			ax.set_ylim([-40,40])
			plt.title(loclist[i])
			ax.set_ylabel('WRF - Obs [ppb]',fontsize=13)
			ax.set_xticklabels(['Base','Nudged','NudgedCOAREG','MACC','SherwenCOAREG'],fontsize=13)
			plt.savefig('Figures/Surfaceplots/Boxplot/SherwenIvsMcDonaldI/'+loclist[i]+'_boxplot',dpi=300)
			#plt.show()
			plt.close()
						
			if i == len(loclist)-1:
				#df = pd.DataFrame(data=[bias_base,bias_fixdep,bias_coareg,bias_nofixdep_nudged,bias_fixdep_nudged,bias_coareg_nudged,bias_cams,
				#			     mae_base,mae_fixdep,mae_coareg,mae_nofixdep_nudged,mae_fixdep_nudged,mae_coareg_nudged,mae_cams],
				#	          index=['bias_base','bias_fixdep','bias_coareg','bias_nofixdep_nudged','bias_fixdep_nudged','bias_coareg_nudged','bias_MACC',
				#			      'mae_base','mae_fixdep','mae_coareg','mae_nofixdep_nudged','mae_fixdep_nudged','mae_coareg_nudged','mae_MACC'],
				#		  columns=loclist_abb).T
				df = pd.DataFrame(data=[bias_base,bias_nofixdep_nudged,bias_coareg_nudged,bias_cams,
							     mae_base,mae_nofixdep_nudged,mae_coareg_nudged,mae_cams],
					          index=['bias_Base','bias_Nudged','bias_NudgedCOAREG','bias_MACC','bias_SherwenCOAREG'
							      'mae_Base','mae_Nudged','mae_NudgedCOAREG','mae_MACC','bias_SherwenCOAREG'],
						  columns=loclist_abb).T
				df['loctype'] = filters
				loctypes = ['HA','zRE','zzDI']
				markertypes = ['o','^','s']
				#colortypes = ['red','green','blue','orange','magenta','cyan','black']
				colortypes = ['red','orange','cyan','blue','magenta']
				plooptext = ['bias','mae']
				
				for p in [0,1]:
					df_sorted = df.sort(['loctype'], ascending=[True])
					
					print(df_sorted)
					
					#runstext = [plooptext[p]+'_base',plooptext[p]+'_fixdep',plooptext[p]+'_coareg',plooptext[p]+'_nofixdep_nudged',plooptext[p]+'_fixdep_nudged',plooptext[p]+'_coareg_nudged',plooptext[p]+'_MACC']
					runstext = [plooptext[p]+'_Base',plooptext[p]+'_Nudged',plooptext[p]+'_NudgedCOAREG',plooptext[p]+'_MACC',plooptext[p]+'_SherwenCOAREG']
					counter=1.
				
					fig, ax = plt.subplots()
					for idz, z in enumerate(loctypes):
						for idy, y in enumerate(runstext):
							plotvalues = df_sorted[y][df_sorted['loctype']==z].values
							xvalues = np.arange(counter,counter+len(plotvalues),1)
							plt.plot(xvalues,plotvalues,linewidth=0,marker=markertypes[idz],markersize=6,color=colortypes[idy],markeredgecolor=colortypes[idy],markeredgewidth=0.0)
						counter=counter+len(plotvalues)
						if idz == 0 or idz == 1:
							plt.vlines(x=counter-0.5,ymin=-100,ymax=100,linewidth=2)
					plt.hlines(y=0,xmin=0,xmax=df_sorted.shape[0]+2,linewidth=1,linestyle=':')
					ax.set_xticks(np.arange(1,df_sorted.shape[0]+1))
					ax.set_xticklabels(df_sorted.index,rotation=45,ha='right')
					ax.set_xlim([0.5,counter-0.5])
					if p == 0:
						ax.set_ylim([-18,12])
					if p == 1:
						ax.set_ylim([0,18])
					ax.set_ylabel(plooptext[p]+' [ppb]',fontsize=13)
					if plooptext[p] == 'mae':
						custom_lines = [Line2D([0], [0], color='red', lw=0, marker='o', markersize=6),
               							Line2D([0], [0], color='orange', lw=0, marker='o', markersize=6),
                						Line2D([0], [0], color='cyan', lw=0, marker='o', markersize=6),
								Line2D([0], [0], color='blue', lw=0, marker='o', markersize=6),
								Line2D([0], [0], color='magenta', lw=0, marker='o', markersize=6)]
						ax.legend(custom_lines, ['Base', 'Nudged', 'NudgedCOAREG', 'MACC', 'SherwenCOAREG'],loc=1)
					plt.savefig('Figures/Surfaceplots/Boxplot/SherwenIvsMcDonaldI/'+plooptext[p]+'_AllstationsCAMS',dpi=300)
					plt.show()
				

			
		if histplot == True:
			fig,ax=plt.subplots(1,1,figsize=(10,10))
			ax.hist(stationdata[loclist[i]+'_o3'],bins=60,range=(0.,60.),normed=True,cumulative=False,histtype='step',alpha=1,color='black',label='Observations')		     
			ax.hist(plotlist_base[i],bins=60,range=(0.,60.),normed=True,cumulative=False,histtype='step',alpha=1,color='red',label='Base')
			#ax.hist(plotlist_run1[i],bins=60,range=(0.,60.),normed=True,cumulative=False,histtype='step',alpha=1,color='green',label='Run 1')
			#ax.hist(plotlist_run2[i].ravel(),bins=60,range=(0.,60.),normed=True,cumulative=False,histtype='step',alpha=1,color='blue',label='Run 2')
			ax.hist(plotlist_run3[i].ravel(),bins=60,range=(0.,60.),normed=True,cumulative=False,histtype='step',alpha=1,color='orange',label='Nudged')
			ax.hist(plotlist_run4[i].ravel(),bins=60,range=(0.,60.),normed=True,cumulative=False,histtype='step',alpha=1,color='magenta',label='SherwenCOAREG')
			ax.hist(plotlist_run5[i].ravel(),bins=60,range=(0.,60.),normed=True,cumulative=False,histtype='step',alpha=1,color='cyan',label='NudgedCOAREG')
			ax.hist(camslist[i].ravel(),bins=60,range=(0.,60.),normed=True,cumulative=False,histtype='step',alpha=1,color='blue',label='MACC')
			ax.yaxis.set_label_position("right")
			ax.yaxis.tick_right()
			ax.set_xlim([0,60])
			ax.set_ylim([0,0.2])
			ax.set_xlabel('Surface ozone mixing ratio [ppb]',fontsize=13)
			ax.set_ylabel('Frequency [-]',fontsize=13)
			ax.legend(loc='upper right')
			plt.title(loclist[i])
			plt.savefig('Figures/Surfaceplots/Histogram/SherwenIvsMcDonaldI/HIST_inclnudged_cams'+loclist[i],dpi=150)
			#plt.show()
			plt.close()				
		
		if dailyaverage == True:
			station_daily = np.zeros(int(np.trunc(len(plotlist_base[i])/24.)))
			base_daily = np.zeros(int(np.trunc(len(plotlist_base[i])/24.)))
			base_run1 = np.zeros(int(np.trunc(len(plotlist_base[i])/24.)))
			base_run2 = np.zeros(int(np.trunc(len(plotlist_base[i])/24.)))
			base_run3 = np.zeros(int(np.trunc(len(plotlist_base[i])/24.)))
			base_run4 = np.zeros(int(np.trunc(len(plotlist_base[i])/24.)))
			base_run5 = np.zeros(int(np.trunc(len(plotlist_base[i])/24.)))
			camsdailyavg = np.zeros(int(np.trunc(len(plotlist_base[i])/24.)))
			timearr_day = np.zeros(int(np.trunc(len(plotlist_base[i])/24.)))
						
			for day in range(0,int(np.trunc(len(plotlist_base[i])/24.))):
				timearr_day[day] = day
				station_daily[day] = np.nanmean(stationdata[loclist[i]+'_o3'][(day*24):(day*24)+23])
				base_daily[day] = np.nanmean(plotlist_base[i][(day*24):(day*24)+23])
				base_run1[day] = np.nanmean(plotlist_run1[i][(day*24):(day*24)+23])
				base_run2[day] = np.nanmean(plotlist_run2[i][(day*24):(day*24)+23])
				base_run3[day] = np.nanmean(plotlist_run3[i][(day*24):(day*24)+23])
				base_run4[day] = np.nanmean(plotlist_run4[i][(day*24):(day*24)+23])
				base_run5[day] = np.nanmean(plotlist_run5[i][(day*24):(day*24)+23])
				camsdailyavg[day] = np.nanmean(camslist[i][(day*24):(day*24)+23])
				
			plt.figure(figsize=(15,8))
			plt.plot(timearr_day,station_daily,'ko',label='Observartions')
			plt.plot(timearr_day,base_daily,'r-',label='Base')
			#plt.plot(timearr_day,base_run1,'g-',label='Run 1')
			#plt.plot(timearr_day,base_run2,'b-',label='Run 2')
			plt.plot(timearr_day,base_run3,color='orange',linestyle='-',label='Nudged')
			plt.plot(timearr_day,base_run4,'m-',label='SherwenCOAREG')
			plt.plot(timearr_day,base_run5,'c-',label='NudgedCOAREG')
			plt.plot(timearr_day,camsdailyavg,'b-',label='MACC')
			plt.ylim([0,60])
			plt.title(loclist[i])
			plt.legend(loc=1)
			plt.savefig('Figures/Surfaceplots/DailyAvg/SherwenIvsMcDonaldI/'+loclist[i]+'_DAILY_inclnudged_cams',dpi=300)
			plt.close()

		'''
		fig, axes = plt.subplots(2, 2)
		ax1 = axes[0,0]
		ax2 = axes[0,1]
		ax3 = axes[1,0]
		ax4 = axes[1,1]
		ax1.plot(timearr,stationdata[loclist[i]+'_o3'],'bo',label='Observartions')
		ax1.plot(timearr,plotlist_base[i],'k-',label='Baserun')
		ax1.plot(timearr,plotlist_run1[i],'g-',label='Run 1')
		#ax2.plot(timearr,plotlist_run2[i],'y-',label='Run 2')
		ax1.set_ylim([0,60])
		ax1.set_title(loclist[i])
		ax1.legend(loc=2)
		x,y=stationdata[loclist[i]+'_o3'][~np.isnan(stationdata[loclist[i]+'_o3']-plotlist_base[i])],plotlist_base[i][~np.isnan(stationdata[loclist[i]+'_o3']-plotlist_base[i])]
		xy = np.vstack([x,y])
		z = gaussian_kde(xy)(xy)
		idx = z.argsort()
		x, y, z = x[idx], y[idx], z[idx]
		ax2.scatter(x, y, c=z, cmap='viridis' ,s=5, edgecolor='')
		ax2.scatter(stationdata[loclist[i]+'_o3'],plotlist_base[i])
		ax2.plot([0,100],[0,100],'r-',label='1:1 line')
		ax2.set_xlim([0,60])
		ax2.set_ylim([0,60])
		ax2.set_xlabel('Observed O3 mixing ratio [ppb]')
		ax2.set_ylabel('Modelled O3 mixing ratio [ppb]')
		ax2.text(2,55,'MAE = '+str(round(mae,2))+' ppb')
		ax2.text(2,57,'ME = '+str(round(me,2))+' ppb')
		x,y=stationdata[loclist[i]+'_o3'][~np.isnan(stationdata[loclist[i]+'_o3']-plotlist_run1[i])],plotlist_run1[i][~np.isnan(stationdata[loclist[i]+'_o3']-plotlist_run1[i])]
		xy = np.vstack([x,y])
		z = gaussian_kde(xy)(xy)
		idx = z.argsort()
		x, y, z = x[idx], y[idx], z[idx]
		ax3.scatter(x, y, c=z, cmap='viridis' ,s=5, edgecolor='')
		ax3.plot([0,100],[0,100],'r-',label='1:1 line')
		ax3.set_xlim([0,60])
		ax3.set_ylim([0,60])
		ax3.set_xlabel('Observed O3 mixing ratio [ppb]')
		ax3.set_ylabel('Modelled O3 mixing ratio [ppb]')
		ax3.text(2,55,'MAE = '+str(round(mae,2))+' ppb')
		ax3.text(2,57,'ME = '+str(round(me,2))+' ppb')
		#plt.savefig('Figures/Surfaceplots/'+loclist[i])
		plt.show()
		'''
		
#loaddata('/home/WUR/barte035/WRFChem/o3_analysis_DATA/o3_station_data_new.csv','/archive/ESG/barte035/MOSAiC/o3_analysis/wrfout_polar_chem_5days')
loaddata('/home/WUR/barte035/WRFChem/o3_analysis_DATA/o3_station_data_new.csv','/lustre/backup/WUR/ESG/barte035/wrfout_chemdt10_d01_2008-08-10_00:00:00',True,True,True)
