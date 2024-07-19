from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import datetime as dt
from sklearn.metrics import mean_squared_error
from math import sqrt


#34 launches in total:
#Chronological order of launches (time in UTC):
# 2008-08-12 17:16 Summit
# 2008-08-13 11:00 Lerwick
# 2008-08-13 11:00 Sodankyla
# 2008-08-13 11:10 NyAlesund
# 2008-08-13 23:21 Resolute
# 2008-08-15 11:00 Scoresbysund
# 2008-08-20 10:54 NyAlesund
# 2008-08-20 11:00 Lerwick
# 2008-08-20 11:00 Sodankyla
# 2008-08-20 16:57 Summit
# 2008-08-20 23:15 Alert
# 2008-08-20 23:15 Eureka
# 2008-08-22 17:36 ASCOS
# 2008-08-22 23:00 Scoresbysund
# 2008-08-25 11:39 ASCOS
# 2008-08-26 07:36 Sodankyla
# 2008-08-27 10:54 NyAlesund
# 2008-08-27 11:00 Lerwick
# 2008-08-27 11:21 ASCOS
# 2008-08-27 16:20 Summit
# 2008-08-27 23:15 Alert
# 2008-08-27 23:15 Eureka
# 2008-08-28 10:00 Scoresbysund
# 2008-08-29 11:41 ACSOS
# 2008-08-31 11:40 ASCOS
# 2008-09-03 11:00 Lerwick
# 2008-09-03 11:00 Sodankyla
# 2008-09-03 11:08 NyAlesund
# 2008-09-03 19:21 Summit
# 2008-09-03 23:18 Alert
# 2008-09-04 11:00 Scoresbysund
# 2008-09-04 17:15 ASCOS
# 2008-09-04 23:28 Resolute
# 2008-09-06 11:26 ASCOS

def sondecomparison(sondedatamap,wrfdatafile):
	bf1 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-10_00:00:00','r')
	bf2 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-13_01:00:00','r')
	bf3 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-21_01:00:00','r')
	bf4 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chemdt15_d01_2008-08-31_01:00:00','r')
	cf1 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_coaregv2_d01_2008-08-10_00:00:00','r')
	cf2 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_coaregv2_d01_2008-08-15_01:00:00','r')	
	df1 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chem_nudgedBL_coareg_d01_2008-08-10_00:00:00','r')
	df2 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chem_nudgedBL_coareg_d01_2008-08-11_01:00:00','r')
	df3 = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chem_nudgedBL_coareg_d01_2008-08-14_01:00:00','r')
	#wrfdata = nc.Dataset(wrfdatafile,'r')
	wrfdata_o3 = np.concatenate((bf1['o3'][0:73,:,:,:],bf2['o3'][0:192,:,:,:],bf3['o3'][0:240,:,:,:],bf4['o3'][:,:,:,:]),axis=0)
	wrfdata_p = np.concatenate((bf1['P'][0:73,:,:,:],bf2['P'][0:192,:,:,:],bf3['P'][0:240,:,:,:],bf4['P'][:,:,:,:]),axis=0)
	wrfdata_pb = np.concatenate((bf1['PB'][0:73,:,:,:],bf2['PB'][0:192,:,:,:],bf3['PB'][0:240,:,:,:],bf4['PB'][:,:,:,:]),axis=0)
	wrfdata_fixeddep = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_fixeddep_d01_2008-08-10_00:00:00','r')
	wrfdata_nudgedBL_nofixeddep = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_chem_nudgedBL_nofixeddep_d01_2008-08-10_00:00:00','r')
	wrfdata_nudgedBL_fixeddep = nc.Dataset('/lustre/backup/WUR/ESG/barte035/wrfout_nudgedBL_fixeddep_d01_2008-08-10_00:00:00','r')
	wrfdata_coareg_o3 = np.concatenate((cf1['o3'][0:121,:,:,:],cf2['o3'][:,:,:,:]),axis=0)
	wrfdata_coareg_p = np.concatenate((cf1['P'][0:121,:,:,:],cf2['P'][:,:,:,:]),axis=0)
	wrfdata_coareg_pb = np.concatenate((cf1['PB'][0:121,:,:,:],cf2['PB'][:,:,:,:]),axis=0)
	wrfdata_nudgedBL_coareg_o3 = np.concatenate((df1['o3'][0:24,:,:,:],df2['o3'][0:73,:,:,:],df3['o3'][:,:,:,:]),axis=0)
	wrfdata_nudgedBL_coareg_p = np.concatenate((df1['P'][0:24,:,:,:],df2['P'][0:73,:,:,:],df3['P'][:,:,:,:]),axis=0)
	wrfdata_nudgedBL_coareg_pb = np.concatenate((df1['PB'][0:24,:,:,:],df2['PB'][0:73,:,:,:],df3['PB'][:,:,:,:]),axis=0)

	#subsetting data based on wrf output
	wrfstarttime = dt.datetime(2008,8,10,00,00,00)
	timearr = np.array([wrfstarttime + dt.timedelta(hours=i) for i in xrange(wrfdata_coareg_o3.shape[0])])	

	#get wrf locations
	loc_sodankyla = [48,162] 	#[67.37,26.63]
    	loc_ascos = [114,123] 		#[87.4,-6.0] #rough, NEEDS TIME ADAPTIVE LOCATION
    	loc_summit = [74,84] 		#[72.5800018311,-38.4799995422]
    	loc_lerwick = [11,122] 		#[60.13922,-1.185319]
	loc_alert = [111,99]		#[82.50,-62.33]
	loc_eureka = [121,86]		#[79.99,-85.94]
	loc_nyalesund =	[85,132]	#[78.92,11.92] #same as zeppelin (makes sense)
	loc_resolute = [129,67]		#[74.71,-94.97]
	loc_scoresbysund = [57,97]	#[70.4848,-21.9512]

	wrf_sodankyla = wrfdata_o3[:,:,loc_sodankyla[0],loc_sodankyla[1]]*1000. #to ppb
	wrf_ascos = wrfdata_o3[:,:,loc_ascos[0],loc_ascos[1]]*1000. #to ppb
	wrf_summit = wrfdata_o3[:,:,loc_summit[0],loc_summit[1]]*1000. #to ppb
	wrf_lerwick = wrfdata_o3[:,:,loc_lerwick[0],loc_lerwick[1]]*1000. #to ppb
	wrf_alert = wrfdata_o3[:,:,loc_alert[0],loc_alert[1]]*1000. #to ppb
	wrf_eureka = wrfdata_o3[:,:,loc_eureka[0],loc_eureka[1]]*1000. #to ppb
	wrf_nyalesund = wrfdata_o3[:,:,loc_nyalesund[0],loc_nyalesund[1]]*1000. #to ppb
	wrf_resolute = wrfdata_o3[:,:,loc_resolute[0],loc_resolute[1]]*1000. #to ppb
	wrf_scoresbysund = wrfdata_o3[:,:,loc_scoresbysund[0],loc_scoresbysund[1]]*1000. #to ppb	
	wrfp_sodankyla = wrfdata_p[:,:,loc_sodankyla[0],loc_sodankyla[1]]+wrfdata_pb[:,:,loc_sodankyla[0],loc_sodankyla[1]]
	wrfp_ascos = wrfdata_p[:,:,loc_ascos[0],loc_ascos[1]]+wrfdata_pb[:,:,loc_ascos[0],loc_ascos[1]]
	wrfp_summit = wrfdata_p[:,:,loc_summit[0],loc_summit[1]]+wrfdata_pb[:,:,loc_summit[0],loc_summit[1]]
	wrfp_lerwick = wrfdata_p[:,:,loc_lerwick[0],loc_lerwick[1]]+wrfdata_pb[:,:,loc_lerwick[0],loc_lerwick[1]]
	wrfp_alert = wrfdata_p[:,:,loc_alert[0],loc_alert[1]]+wrfdata_pb[:,:,loc_alert[0],loc_alert[1]]
	wrfp_eureka = wrfdata_p[:,:,loc_eureka[0],loc_eureka[1]]+wrfdata_pb[:,:,loc_eureka[0],loc_eureka[1]]
	wrfp_nyalesund = wrfdata_p[:,:,loc_nyalesund[0],loc_nyalesund[1]]+wrfdata_pb[:,:,loc_nyalesund[0],loc_nyalesund[1]]
	wrfp_resolute = wrfdata_p[:,:,loc_resolute[0],loc_resolute[1]]+wrfdata_pb[:,:,loc_resolute[0],loc_resolute[1]]
	wrfp_scoresbysund = wrfdata_p[:,:,loc_scoresbysund[0],loc_scoresbysund[1]]+wrfdata_pb[:,:,loc_scoresbysund[0],loc_scoresbysund[1]]
	
	#run1
	wrf_sodankyla_run1 = wrfdata_fixeddep['o3'][:,:,loc_sodankyla[0],loc_sodankyla[1]]*1000. #to ppb
	wrf_ascos_run1 = wrfdata_fixeddep['o3'][:,:,loc_ascos[0],loc_ascos[1]]*1000. #to ppb
	wrf_summit_run1 = wrfdata_fixeddep['o3'][:,:,loc_summit[0],loc_summit[1]]*1000. #to ppb
	wrf_lerwick_run1 = wrfdata_fixeddep['o3'][:,:,loc_lerwick[0],loc_lerwick[1]]*1000. #to ppb
	wrf_alert_run1 = wrfdata_fixeddep['o3'][:,:,loc_alert[0],loc_alert[1]]*1000. #to ppb
	wrf_eureka_run1 = wrfdata_fixeddep['o3'][:,:,loc_eureka[0],loc_eureka[1]]*1000. #to ppb
	wrf_nyalesund_run1 = wrfdata_fixeddep['o3'][:,:,loc_nyalesund[0],loc_nyalesund[1]]*1000. #to ppb
	wrf_resolute_run1 = wrfdata_fixeddep['o3'][:,:,loc_resolute[0],loc_resolute[1]]*1000. #to ppb
	wrf_scoresbysund_run1 = wrfdata_fixeddep['o3'][:,:,loc_scoresbysund[0],loc_scoresbysund[1]]*1000. #to ppb	
	wrfp_sodankyla_run1 = wrfdata_fixeddep['P'][:,:,loc_sodankyla[0],loc_sodankyla[1]]+wrfdata_fixeddep['PB'][:,:,loc_sodankyla[0],loc_sodankyla[1]]
	wrfp_ascos_run1 = wrfdata_fixeddep['P'][:,:,loc_ascos[0],loc_ascos[1]]+wrfdata_fixeddep['PB'][:,:,loc_ascos[0],loc_ascos[1]]
	wrfp_summit_run1 = wrfdata_fixeddep['P'][:,:,loc_summit[0],loc_summit[1]]+wrfdata_fixeddep['PB'][:,:,loc_summit[0],loc_summit[1]]
	wrfp_lerwick_run1 = wrfdata_fixeddep['P'][:,:,loc_lerwick[0],loc_lerwick[1]]+wrfdata_fixeddep['PB'][:,:,loc_lerwick[0],loc_lerwick[1]]
	wrfp_alert_run1 = wrfdata_fixeddep['P'][:,:,loc_alert[0],loc_alert[1]]+wrfdata_fixeddep['PB'][:,:,loc_alert[0],loc_alert[1]]
	wrfp_eureka_run1 = wrfdata_fixeddep['P'][:,:,loc_eureka[0],loc_eureka[1]]+wrfdata_fixeddep['PB'][:,:,loc_eureka[0],loc_eureka[1]]
	wrfp_nyalesund_run1 = wrfdata_fixeddep['P'][:,:,loc_nyalesund[0],loc_nyalesund[1]]+wrfdata_fixeddep['PB'][:,:,loc_nyalesund[0],loc_nyalesund[1]]
	wrfp_resolute_run1 = wrfdata_fixeddep['P'][:,:,loc_resolute[0],loc_resolute[1]]+wrfdata_fixeddep['PB'][:,:,loc_resolute[0],loc_resolute[1]]
	wrfp_scoresbysund_run1 = wrfdata_fixeddep['P'][:,:,loc_scoresbysund[0],loc_scoresbysund[1]]+wrfdata_fixeddep['PB'][:,:,loc_scoresbysund[0],loc_scoresbysund[1]]

	#run2
	wrf_sodankyla_run2 = wrfdata_coareg_o3[:,:,loc_sodankyla[0],loc_sodankyla[1]]*1000. #to ppb
	wrf_ascos_run2 = wrfdata_coareg_o3[:,:,loc_ascos[0],loc_ascos[1]]*1000. #to ppb
	wrf_summit_run2 = wrfdata_coareg_o3[:,:,loc_summit[0],loc_summit[1]]*1000. #to ppb
	wrf_lerwick_run2 = wrfdata_coareg_o3[:,:,loc_lerwick[0],loc_lerwick[1]]*1000. #to ppb
	wrf_alert_run2 = wrfdata_coareg_o3[:,:,loc_alert[0],loc_alert[1]]*1000. #to ppb
	wrf_eureka_run2 = wrfdata_coareg_o3[:,:,loc_eureka[0],loc_eureka[1]]*1000. #to ppb
	wrf_nyalesund_run2 = wrfdata_coareg_o3[:,:,loc_nyalesund[0],loc_nyalesund[1]]*1000. #to ppb
	wrf_resolute_run2 = wrfdata_coareg_o3[:,:,loc_resolute[0],loc_resolute[1]]*1000. #to ppb
	wrf_scoresbysund_run2 = wrfdata_coareg_o3[:,:,loc_scoresbysund[0],loc_scoresbysund[1]]*1000. #to ppb	
	wrfp_sodankyla_run2 = wrfdata_coareg_p[:,:,loc_sodankyla[0],loc_sodankyla[1]]+wrfdata_coareg_pb[:,:,loc_sodankyla[0],loc_sodankyla[1]]
	wrfp_ascos_run2 = wrfdata_coareg_p[:,:,loc_ascos[0],loc_ascos[1]]+wrfdata_coareg_pb[:,:,loc_ascos[0],loc_ascos[1]]
	wrfp_summit_run2 = wrfdata_coareg_p[:,:,loc_summit[0],loc_summit[1]]+wrfdata_coareg_pb[:,:,loc_summit[0],loc_summit[1]]
	wrfp_lerwick_run2 = wrfdata_coareg_p[:,:,loc_lerwick[0],loc_lerwick[1]]+wrfdata_coareg_pb[:,:,loc_lerwick[0],loc_lerwick[1]]
	wrfp_alert_run2 = wrfdata_coareg_p[:,:,loc_alert[0],loc_alert[1]]+wrfdata_coareg_pb[:,:,loc_alert[0],loc_alert[1]]
	wrfp_eureka_run2 = wrfdata_coareg_p[:,:,loc_eureka[0],loc_eureka[1]]+wrfdata_coareg_pb[:,:,loc_eureka[0],loc_eureka[1]]
	wrfp_nyalesund_run2 = wrfdata_coareg_p[:,:,loc_nyalesund[0],loc_nyalesund[1]]+wrfdata_coareg_pb[:,:,loc_nyalesund[0],loc_nyalesund[1]]
	wrfp_resolute_run2 = wrfdata_coareg_p[:,:,loc_resolute[0],loc_resolute[1]]+wrfdata_coareg_pb[:,:,loc_resolute[0],loc_resolute[1]]
	wrfp_scoresbysund_run2 = wrfdata_coareg_p[:,:,loc_scoresbysund[0],loc_scoresbysund[1]]+wrfdata_coareg_pb[:,:,loc_scoresbysund[0],loc_scoresbysund[1]]
	
	#run3
	wrf_sodankyla_run3 = wrfdata_nudgedBL_nofixeddep['o3'][:,:,loc_sodankyla[0],loc_sodankyla[1]]*1000. #to ppb
	wrf_ascos_run3 = wrfdata_nudgedBL_nofixeddep['o3'][:,:,loc_ascos[0],loc_ascos[1]]*1000. #to ppb
	wrf_summit_run3 = wrfdata_nudgedBL_nofixeddep['o3'][:,:,loc_summit[0],loc_summit[1]]*1000. #to ppb
	wrf_lerwick_run3 = wrfdata_nudgedBL_nofixeddep['o3'][:,:,loc_lerwick[0],loc_lerwick[1]]*1000. #to ppb
	wrf_alert_run3 = wrfdata_nudgedBL_nofixeddep['o3'][:,:,loc_alert[0],loc_alert[1]]*1000. #to ppb
	wrf_eureka_run3 = wrfdata_nudgedBL_nofixeddep['o3'][:,:,loc_eureka[0],loc_eureka[1]]*1000. #to ppb
	wrf_nyalesund_run3 = wrfdata_nudgedBL_nofixeddep['o3'][:,:,loc_nyalesund[0],loc_nyalesund[1]]*1000. #to ppb
	wrf_resolute_run3 = wrfdata_nudgedBL_nofixeddep['o3'][:,:,loc_resolute[0],loc_resolute[1]]*1000. #to ppb
	wrf_scoresbysund_run3 = wrfdata_nudgedBL_nofixeddep['o3'][:,:,loc_scoresbysund[0],loc_scoresbysund[1]]*1000. #to ppb	
	wrfp_sodankyla_run3 = wrfdata_nudgedBL_nofixeddep['P'][:,:,loc_sodankyla[0],loc_sodankyla[1]]+wrfdata_nudgedBL_nofixeddep['PB'][:,:,loc_sodankyla[0],loc_sodankyla[1]]
	wrfp_ascos_run3 = wrfdata_nudgedBL_nofixeddep['P'][:,:,loc_ascos[0],loc_ascos[1]]+wrfdata_nudgedBL_nofixeddep['PB'][:,:,loc_ascos[0],loc_ascos[1]]
	wrfp_summit_run3 = wrfdata_nudgedBL_nofixeddep['P'][:,:,loc_summit[0],loc_summit[1]]+wrfdata_nudgedBL_nofixeddep['PB'][:,:,loc_summit[0],loc_summit[1]]
	wrfp_lerwick_run3 = wrfdata_nudgedBL_nofixeddep['P'][:,:,loc_lerwick[0],loc_lerwick[1]]+wrfdata_nudgedBL_nofixeddep['PB'][:,:,loc_lerwick[0],loc_lerwick[1]]
	wrfp_alert_run3 = wrfdata_nudgedBL_nofixeddep['P'][:,:,loc_alert[0],loc_alert[1]]+wrfdata_nudgedBL_nofixeddep['PB'][:,:,loc_alert[0],loc_alert[1]]
	wrfp_eureka_run3 = wrfdata_nudgedBL_nofixeddep['P'][:,:,loc_eureka[0],loc_eureka[1]]+wrfdata_nudgedBL_nofixeddep['PB'][:,:,loc_eureka[0],loc_eureka[1]]
	wrfp_nyalesund_run3 = wrfdata_nudgedBL_nofixeddep['P'][:,:,loc_nyalesund[0],loc_nyalesund[1]]+wrfdata_nudgedBL_nofixeddep['PB'][:,:,loc_nyalesund[0],loc_nyalesund[1]]
	wrfp_resolute_run3 = wrfdata_nudgedBL_nofixeddep['P'][:,:,loc_resolute[0],loc_resolute[1]]+wrfdata_nudgedBL_nofixeddep['PB'][:,:,loc_resolute[0],loc_resolute[1]]
	wrfp_scoresbysund_run3 = wrfdata_nudgedBL_nofixeddep['P'][:,:,loc_scoresbysund[0],loc_scoresbysund[1]]+wrfdata_nudgedBL_nofixeddep['PB'][:,:,loc_scoresbysund[0],loc_scoresbysund[1]]

	#run4
	wrf_sodankyla_run4 = wrfdata_nudgedBL_fixeddep['o3'][:,:,loc_sodankyla[0],loc_sodankyla[1]]*1000. #to ppb
	wrf_ascos_run4 = wrfdata_nudgedBL_fixeddep['o3'][:,:,loc_ascos[0],loc_ascos[1]]*1000. #to ppb
	wrf_summit_run4 = wrfdata_nudgedBL_fixeddep['o3'][:,:,loc_summit[0],loc_summit[1]]*1000. #to ppb
	wrf_lerwick_run4 = wrfdata_nudgedBL_fixeddep['o3'][:,:,loc_lerwick[0],loc_lerwick[1]]*1000. #to ppb
	wrf_alert_run4 = wrfdata_nudgedBL_fixeddep['o3'][:,:,loc_alert[0],loc_alert[1]]*1000. #to ppb
	wrf_eureka_run4 = wrfdata_nudgedBL_fixeddep['o3'][:,:,loc_eureka[0],loc_eureka[1]]*1000. #to ppb
	wrf_nyalesund_run4 = wrfdata_nudgedBL_fixeddep['o3'][:,:,loc_nyalesund[0],loc_nyalesund[1]]*1000. #to ppb
	wrf_resolute_run4 = wrfdata_nudgedBL_fixeddep['o3'][:,:,loc_resolute[0],loc_resolute[1]]*1000. #to ppb
	wrf_scoresbysund_run4 = wrfdata_nudgedBL_fixeddep['o3'][:,:,loc_scoresbysund[0],loc_scoresbysund[1]]*1000. #to ppb	
	wrfp_sodankyla_run4 = wrfdata_nudgedBL_fixeddep['P'][:,:,loc_sodankyla[0],loc_sodankyla[1]]+wrfdata_nudgedBL_fixeddep['PB'][:,:,loc_sodankyla[0],loc_sodankyla[1]]
	wrfp_ascos_run4 = wrfdata_nudgedBL_fixeddep['P'][:,:,loc_ascos[0],loc_ascos[1]]+wrfdata_nudgedBL_fixeddep['PB'][:,:,loc_ascos[0],loc_ascos[1]]
	wrfp_summit_run4 = wrfdata_nudgedBL_fixeddep['P'][:,:,loc_summit[0],loc_summit[1]]+wrfdata_nudgedBL_fixeddep['PB'][:,:,loc_summit[0],loc_summit[1]]
	wrfp_lerwick_run4 = wrfdata_nudgedBL_fixeddep['P'][:,:,loc_lerwick[0],loc_lerwick[1]]+wrfdata_nudgedBL_fixeddep['PB'][:,:,loc_lerwick[0],loc_lerwick[1]]
	wrfp_alert_run4 = wrfdata_nudgedBL_fixeddep['P'][:,:,loc_alert[0],loc_alert[1]]+wrfdata_nudgedBL_fixeddep['PB'][:,:,loc_alert[0],loc_alert[1]]
	wrfp_eureka_run4 = wrfdata_nudgedBL_fixeddep['P'][:,:,loc_eureka[0],loc_eureka[1]]+wrfdata_nudgedBL_fixeddep['PB'][:,:,loc_eureka[0],loc_eureka[1]]
	wrfp_nyalesund_run4 = wrfdata_nudgedBL_fixeddep['P'][:,:,loc_nyalesund[0],loc_nyalesund[1]]+wrfdata_nudgedBL_fixeddep['PB'][:,:,loc_nyalesund[0],loc_nyalesund[1]]
	wrfp_resolute_run4 = wrfdata_nudgedBL_fixeddep['P'][:,:,loc_resolute[0],loc_resolute[1]]+wrfdata_nudgedBL_fixeddep['PB'][:,:,loc_resolute[0],loc_resolute[1]]
	wrfp_scoresbysund_run4 = wrfdata_nudgedBL_fixeddep['P'][:,:,loc_scoresbysund[0],loc_scoresbysund[1]]+wrfdata_nudgedBL_fixeddep['PB'][:,:,loc_scoresbysund[0],loc_scoresbysund[1]]

	#run5
	wrf_sodankyla_run5 = wrfdata_nudgedBL_coareg_o3[:,:,loc_sodankyla[0],loc_sodankyla[1]]*1000. #to ppb
	wrf_ascos_run5 = wrfdata_nudgedBL_coareg_o3[:,:,loc_ascos[0],loc_ascos[1]]*1000. #to ppb
	wrf_summit_run5 = wrfdata_nudgedBL_coareg_o3[:,:,loc_summit[0],loc_summit[1]]*1000. #to ppb
	wrf_lerwick_run5 = wrfdata_nudgedBL_coareg_o3[:,:,loc_lerwick[0],loc_lerwick[1]]*1000. #to ppb
	wrf_alert_run5 = wrfdata_nudgedBL_coareg_o3[:,:,loc_alert[0],loc_alert[1]]*1000. #to ppb
	wrf_eureka_run5 = wrfdata_nudgedBL_coareg_o3[:,:,loc_eureka[0],loc_eureka[1]]*1000. #to ppb
	wrf_nyalesund_run5 = wrfdata_nudgedBL_coareg_o3[:,:,loc_nyalesund[0],loc_nyalesund[1]]*1000. #to ppb
	wrf_resolute_run5 = wrfdata_nudgedBL_coareg_o3[:,:,loc_resolute[0],loc_resolute[1]]*1000. #to ppb
	wrf_scoresbysund_run5 = wrfdata_nudgedBL_coareg_o3[:,:,loc_scoresbysund[0],loc_scoresbysund[1]]*1000. #to ppb	
	wrfp_sodankyla_run5 = wrfdata_nudgedBL_coareg_p[:,:,loc_sodankyla[0],loc_sodankyla[1]]+wrfdata_nudgedBL_coareg_pb[:,:,loc_sodankyla[0],loc_sodankyla[1]]
	wrfp_ascos_run5 = wrfdata_nudgedBL_coareg_p[:,:,loc_ascos[0],loc_ascos[1]]+wrfdata_nudgedBL_coareg_pb[:,:,loc_ascos[0],loc_ascos[1]]
	wrfp_summit_run5 = wrfdata_nudgedBL_coareg_p[:,:,loc_summit[0],loc_summit[1]]+wrfdata_nudgedBL_coareg_pb[:,:,loc_summit[0],loc_summit[1]]
	wrfp_lerwick_run5 = wrfdata_nudgedBL_coareg_p[:,:,loc_lerwick[0],loc_lerwick[1]]+wrfdata_nudgedBL_coareg_pb[:,:,loc_lerwick[0],loc_lerwick[1]]
	wrfp_alert_run5 = wrfdata_nudgedBL_coareg_p[:,:,loc_alert[0],loc_alert[1]]+wrfdata_nudgedBL_coareg_pb[:,:,loc_alert[0],loc_alert[1]]
	wrfp_eureka_run5 = wrfdata_nudgedBL_coareg_p[:,:,loc_eureka[0],loc_eureka[1]]+wrfdata_nudgedBL_coareg_pb[:,:,loc_eureka[0],loc_eureka[1]]
	wrfp_nyalesund_run5 = wrfdata_nudgedBL_coareg_p[:,:,loc_nyalesund[0],loc_nyalesund[1]]+wrfdata_nudgedBL_coareg_pb[:,:,loc_nyalesund[0],loc_nyalesund[1]]
	wrfp_resolute_run5 = wrfdata_nudgedBL_coareg_p[:,:,loc_resolute[0],loc_resolute[1]]+wrfdata_nudgedBL_coareg_pb[:,:,loc_resolute[0],loc_resolute[1]]
	wrfp_scoresbysund_run5 = wrfdata_nudgedBL_coareg_p[:,:,loc_scoresbysund[0],loc_scoresbysund[1]]+wrfdata_nudgedBL_coareg_pb[:,:,loc_scoresbysund[0],loc_scoresbysund[1]]
	
	#24 unique timesteps to analyze
	timelist = [dt.datetime(2008,8,12,17,00,00),	# 2008-08-12 17:16 Summit
		    dt.datetime(2008,8,13,11,00,00),	# 2008-08-13 11:00 Lerwick # 2008-08-13 11:00 Sodankyla # 2008-08-13 11:10 NyAlesund
		    dt.datetime(2008,8,13,23,00,00),	# 2008-08-13 23:21 Resolute
		    dt.datetime(2008,8,15,11,00,00),	# 2008-08-15 11:00 Scoresbysund
		    dt.datetime(2008,8,20,11,00,00),	# 2008-08-20 10:54 NyAlesund # 2008-08-20 11:00 Lerwick # 2008-08-20 11:00 Sodankyla
		    dt.datetime(2008,8,20,17,00,00),	# 2008-08-20 16:57 Summit
		    dt.datetime(2008,8,20,23,00,00),	# 2008-08-20 23:15 Alert # 2008-08-20 23:15 Eureka
		    dt.datetime(2008,8,22,18,00,00),	# 2008-08-22 17:36 ASCOS
		    dt.datetime(2008,8,22,23,00,00),	# 2008-08-22 23:00 Scoresbysund
		    dt.datetime(2008,8,25,12,00,00),	# 2008-08-25 11:39 ASCOS
		    dt.datetime(2008,8,26,8,00,00),	# 2008-08-26 07:36 Sodankyla
		    dt.datetime(2008,8,27,11,00,00),	# 2008-08-27 10:54 NyAlesund # 2008-08-27 11:00 Lerwick # 2008-08-27 11:21 ASCOS
		    dt.datetime(2008,8,27,16,00,00),	# 2008-08-27 16:20 Summit
		    dt.datetime(2008,8,27,23,00,00),	# 2008-08-27 23:15 Alert # 2008-08-27 23:15 Eureka
		    dt.datetime(2008,8,28,10,00,00),	# 2008-08-28 10:00 Scoresbysund
		    dt.datetime(2008,8,29,11,00,00),	# 2008-08-29 11:41 ACSOS
		    dt.datetime(2008,8,31,11,00,00),	# 2008-08-31 11:40 ASCOS
		    dt.datetime(2008,9,3,11,00,00),	# 2008-09-03 11:00 Lerwick # 2008-09-03 11:00 Sodankyla # 2008-09-03 11:08 NyAlesund
		    dt.datetime(2008,9,3,19,00,00),	# 2008-09-03 19:21 Summit
		    dt.datetime(2008,9,3,23,00,00),	# 2008-09-03 23:18 Alert
		    dt.datetime(2008,9,4,11,00,00),	# 2008-09-04 11:00 Scoresbysund
		    dt.datetime(2008,9,4,17,00,00),	# 2008-09-04 17:15 ASCOS
		    dt.datetime(2008,9,4,23,00,00),	# 2008-09-04 23:28 Resolute
		    dt.datetime(2008,9,6,11,00,00)]	# 2008-09-06 11:26 ASCOS
		    
	timelist_locs = [['Summit'],
			 ['Lerwick','Sodankyla','NyAlesund'],
			 ['Resolute'],
			 ['Scoresbysund'],
			 ['Lerwick','Sodankyla','NyAlesund'],
			 ['Summit'],
			 ['Alert','Eureka'],
			 ['ASCOS'],
			 ['Scoresbysund'],
			 ['ASCOS'],
			 ['Sodankyla'],
			 ['Lerwick','NyAlesund','ASCOS'],
			 ['Summit'],
			 ['Alert','Eureka'],
			 ['Scoresbysund'],
			 ['ACSOS'],
			 ['ACSOS'],
			 ['Lerwick','Sodankyla','NyAlesund'],
			 ['Summit'],
			 ['Alert'],
			 ['Scoresbysund'],
			 ['ASCOS'],
			 ['Resolute'],
			 ['ASCOS']]
			 
	
	##Get array of lists of pressure, o3 and int(numsondes)
	sum_sondep,sum_sondeo3,sum_numsondes = readdata('Summit','/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes/Summit/',115)
	eur_sondep,eur_sondeo3,eur_numsondes = readdata('Eureka','/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes/Eureka/',67)
	alt_sondep,alt_sondeo3,alt_numsondes = readdata('Alert','/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes/Alert/',64)
	nya_sondep,nya_sondeo3,nya_numsondes = readdata('NyAlesund','/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes/NyAlesund/',152)
	sco_sondep,sco_sondeo3,sco_numsondes = readdata('Scoresbysund','/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes/Scoresbysund/',144)
	sod_sondep,sod_sondeo3,sod_numsondes = readdata('Sodankyla','/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes/Sodankyla/',144)
	ler_sondep,ler_sondeo3,ler_numsondes = readdata('Lerwick','/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes/Lerwick/',45)
	res_sondep,res_sondeo3,res_numsondes = readdata('Resolute','/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes/Resolute/',52)
	asc_sondep,asc_sondeo3,asc_numsondes = readdata('ASCOS','/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes/ASCOS/',1)
		
	#Plot WRF and radiosondes per location and timestep:
	plotdata('Summit',wrf_summit,wrfp_summit,wrf_summit_run1,wrfp_summit_run1,wrf_summit_run2,wrfp_summit_run2,wrf_summit_run3,wrfp_summit_run3,wrf_summit_run4,wrfp_summit_run4,wrf_summit_run5,wrfp_summit_run5,timearr,timelist,timelist_locs,sum_sondep,sum_sondeo3,sum_numsondes)
	plotdata('Eureka',wrf_eureka,wrfp_eureka,wrf_eureka_run1,wrfp_eureka_run1,wrf_eureka_run2,wrfp_eureka_run2,wrf_eureka_run3,wrfp_eureka_run3,wrf_eureka_run4,wrfp_eureka_run4,wrf_eureka_run5,wrfp_eureka_run5,timearr,timelist,timelist_locs,eur_sondep,eur_sondeo3,eur_numsondes)
	plotdata('Alert',wrf_alert,wrfp_alert,wrf_alert_run1,wrfp_alert_run1,wrf_alert_run2,wrfp_alert_run2,wrf_alert_run3,wrfp_alert_run3,wrf_alert_run4,wrfp_alert_run4,wrf_alert_run5,wrfp_alert_run5,timearr,timelist,timelist_locs,alt_sondep,alt_sondeo3,alt_numsondes)
	plotdata('NyAlesund',wrf_nyalesund,wrfp_nyalesund,wrf_nyalesund_run1,wrfp_nyalesund_run1,wrf_nyalesund_run2,wrfp_nyalesund_run2,wrf_nyalesund_run3,wrfp_nyalesund_run3,wrf_nyalesund_run4,wrfp_nyalesund_run4,wrf_nyalesund_run5,wrfp_nyalesund_run5,timearr,timelist,timelist_locs,nya_sondep,nya_sondeo3,nya_numsondes)
	plotdata('Scoresbysund',wrf_scoresbysund,wrfp_scoresbysund,wrf_scoresbysund_run1,wrfp_scoresbysund_run1,wrf_scoresbysund_run2,wrfp_scoresbysund_run2,wrf_scoresbysund_run3,wrfp_scoresbysund_run3,wrf_scoresbysund_run4,wrfp_scoresbysund_run4,wrf_scoresbysund_run5,wrfp_scoresbysund_run5,timearr,timelist,timelist_locs,sco_sondep,sco_sondeo3,sco_numsondes)
	plotdata('Sodankyla',wrf_sodankyla,wrfp_sodankyla,wrf_sodankyla_run1,wrfp_sodankyla_run1,wrf_sodankyla_run2,wrfp_sodankyla_run2,wrf_sodankyla_run3,wrfp_sodankyla_run3,wrf_sodankyla_run4,wrfp_sodankyla_run4,wrf_sodankyla_run5,wrfp_sodankyla_run5,timearr,timelist,timelist_locs,sod_sondep,sod_sondeo3,sod_numsondes)
	plotdata('Lerwick',wrf_lerwick,wrfp_lerwick,wrf_lerwick_run1,wrfp_lerwick_run1,wrf_lerwick_run2,wrfp_lerwick_run2,wrf_lerwick_run3,wrfp_lerwick_run3,wrf_lerwick_run4,wrfp_lerwick_run4,wrf_lerwick_run5,wrfp_lerwick_run5,timearr,timelist,timelist_locs,ler_sondep,ler_sondeo3,ler_numsondes)
	plotdata('Resolute',wrf_resolute,wrfp_resolute,wrf_resolute_run1,wrfp_resolute_run1,wrf_resolute_run2,wrfp_resolute_run2,wrf_resolute_run3,wrfp_resolute_run3,wrf_resolute_run4,wrfp_resolute_run4,wrf_resolute_run5,wrfp_resolute_run5,timearr,timelist,timelist_locs,res_sondep,res_sondeo3,res_numsondes)
	plotdata('ASCOS',wrf_ascos,wrfp_ascos,wrf_ascos_run1,wrfp_ascos_run1,wrf_ascos_run2,wrfp_ascos_run2,wrf_ascos_run3,wrfp_ascos_run3,wrf_ascos_run4,wrfp_ascos_run4,wrf_ascos_run5,wrfp_ascos_run5,timearr,timelist,timelist_locs,asc_sondep,asc_sondeo3,asc_numsondes)

def readdata(location,mypath,skiplines):
	filelist = [f for f in listdir(mypath) if isfile(join(mypath, f))]
	fileloc = np.zeros(len(filelist))
	sondep = np.empty((len(filelist),),dtype=object)
	sondeo3 = np.empty((len(filelist),),dtype=object)

	if location == 'Summit':
		for i in range(len(filelist)):
			fileloc = pd.read_csv(mypath+filelist[i], skiprows=skiplines, sep='\s+')
			fileloc = fileloc.drop(fileloc.index[0]).astype(float)
			sondep[i] = (fileloc['Press']*100).tolist()		#Pressure in Pa
			sondeo3[i] = (fileloc['O3Mix']*1000).tolist()		#Mixing ratio in ppb	
	if location == 'Alert':
		for i in range(len(filelist)):
			fileloc = pd.read_csv(mypath+filelist[i], skiprows=skiplines, sep='\s+', header=None)
			sondep[i] = (fileloc[1]*100).tolist()
			sondeo3[i] = (fileloc[15]*1000).tolist()
	if location in ['Eureka','NyAlesund','Scoresbysund','Sodankyla']:
		for i in range(len(filelist)):
			fileloc = pd.read_csv(mypath+filelist[i], skiprows=skiplines, sep='\s+', header=None)
			if location == 'Sodankyla' and filelist[i] == 'so080826.b07':
				fileloc = pd.read_csv(mypath+filelist[i], skiprows=58, sep='\s+', header=None)
			sondep[i] = (fileloc[0]*100).tolist()						#Pressure in Pa
			sondeo3[i] = (((fileloc[6]*0.001)/(fileloc[0]*100))*1e9).tolist()		#Mixing ratio in ppb from partial pressure						
	if location in ['Lerwick','Resolute']:
		for i in range(len(filelist)):
			fileloc = pd.read_csv(mypath+filelist[i], skiprows=skiplines)
			if location == 'Lerwick':
				fileloc = fileloc.drop(['LevelCode'],axis=1)
				fileloc = fileloc.dropna(how='any').astype(float)
			sondep[i] = (fileloc['Pressure']*100).tolist()								#Pressure in Pa
			sondeo3[i] = (((fileloc['O3PartialPressure']*0.001)/(fileloc['Pressure']*100))*1e9).tolist()		#Mixing ratio in ppb from partial pressure	
	if location == 'ASCOS':
		#overwrite fileloc,sondep,sondeo
		fileloc = np.zeros(7)
		sondep = np.empty((7,),dtype=object)
		sondeo3 = np.empty((7,),dtype=object)
		fileloc = pd.read_csv(mypath+filelist[0],skiprows=skiplines)
		fileloc = fileloc.drop(fileloc.index[0]).astype(float)
		for i in range(7):
			if i == 0:
				sondep[i] = (fileloc['Pressure']*100).tolist()
				fileloc['Ozone'] = fileloc['Ozone'].replace(99.9,np.nan)
				sondeo3[i] = (((fileloc['Ozone']*0.001)/(fileloc['Pressure']*100))*1e9).tolist()
			if i > 0:
				sondep[i] = (fileloc['Pressure.'+str(i)]*100).tolist()
				fileloc['Ozone.'+str(i)] = fileloc['Ozone.'+str(i)].replace(99.9,np.nan)
				sondeo3[i] = (((fileloc['Ozone.'+str(i)]*0.001)/(fileloc['Pressure.'+str(i)]*100))*1e9).tolist()
			#sondep[i] = [x for x in sondep[i] if ~np.isnan(x)]
			#sondep[i] = [x for x in sondeo3[i] if ~np.isnan(x)]
			#sondeo3[i] = [x for x in sondep[i] if ~np.isnan(x)]
			#sondeo3[i] = [x for x in sondeo3[i] if ~np.isnan(x)]
			
		filelist = [1,2,3,4,5,6,7] #to reset numsondes to 7			
	return sondep,sondeo3,len(filelist)
					
def plotdata(loc,wrf_o3,wrf_p,wrf_o3_run1,wrf_p_run1,wrf_o3_run2,wrf_p_run2,wrf_o3_run3,wrf_p_run3,wrf_o3_run4,wrf_p_run4,wrf_o3_run5,wrf_p_run5,timearr_wrf,timelist_allsondes,timelist_locs_allsondes,sonde_p,sonde_o3,numsondes):
	a = 0
	for i in range(len(timelist_locs_allsondes)):
		if loc in timelist_locs_allsondes[i] and timelist_allsondes[i] in timearr_wrf:
			print('Processing sonde ',loc,timelist_allsondes[i])
			timestep = np.where(timearr_wrf == timelist_allsondes[i])[0][0]
			wrfo3 = wrf_o3[timestep,:]
			wrfp = wrf_p[timestep,:]
			wrfo3_run1 = wrf_o3_run1[timestep,:]
			wrfp_run1 = wrf_p_run1[timestep,:]
			wrfo3_run2 = wrf_o3_run2[timestep,:]
			wrfp_run2 = wrf_p_run2[timestep,:]
			wrfo3_run3 = wrf_o3_run3[timestep,:]
			wrfp_run3 = wrf_p_run3[timestep,:]
			wrfo3_run4 = wrf_o3_run4[timestep,:]
			wrfp_run4 = wrf_p_run4[timestep,:]
			wrfo3_run5 = wrf_o3_run5[timestep,:]
			wrfp_run5 = wrf_p_run5[timestep,:]
			sondeo3 = sonde_o3[a]
			sondep = sonde_p[a]
			a = a+1
			
			plt.plot(wrfo3,wrfp,'ro',label='Base')
			#plt.plot(wrfo3_run1,wrfp_run1,'go',label='Run 1')
			#plt.plot(wrfo3_run2,wrfp_run2,'bo',label='Run 2')
			plt.plot(wrfo3_run3,wrfp_run3,color='orange',marker='o',linewidth=0,label='Nudged')
			#plt.plot(wrfo3_run4,wrfp_run4,'mo',label='Run 4')
			plt.plot(wrfo3_run5,wrfp_run5,'co',label='NudgedCOAREG')
			plt.plot(sondeo3,sondep,color='black',label='Observations')
			plt.gca().invert_yaxis()
			plt.yscale('log')
			plt.ylim([105000,30000])
			plt.yticks([100000,90000,80000,70000,60000,50000,40000,30000],[1000,900,800,700,600,500,400,300])
			plt.xlim([0,100])
			plt.xlabel('Ozone mixing ratio [ppb]')
			plt.ylabel('Pressure [hPa]')
			plt.legend(loc=2)
			plt.title(loc+'   '+timelist_allsondes[i].strftime("%m/%d/%Y, %H:%M:%S"))
			plt.savefig('Figures/Sondeplots/'+loc+'_final_'+timelist_allsondes[i].strftime("%Y%m%d%H%M%S"))
			plt.show()
		
#sondecomparison('/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes','/archive/ESG/barte035/MOSAiC/o3_analysis/wrfout_polar_chem_5days')
sondecomparison('/home/WUR/barte035/WRFChem/o3_analysis_DATA/sondes','/lustre/backup/WUR/ESG/barte035/wrfout_chemdt10_d01_2008-08-10_00:00:00')
