;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   Plot data in xy-plot 
 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

begin
 ; What kind of plot do you prefer?
   type                                      = "newpdf"
   wks                                       = gsn_open_wks(type,"plt_TimeSeries_NOx") 
 ; The WRF ARW input file. This needs to have a ".nc" appended, so just do it.
   fileID                                    = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")
   
 ; --- Define station coordinates and find indices into xlat/xlon ---
   iz                                        = 0
   lat_1                                     = 51.5405 ; Vredepeel 
   lon_1                                     =  5.8531 ; Vredepeel 
   loc_1                                     = wrf_user_latlon_to_ij(fileID, lat_1, lon_1)
   ix_1                                      = loc_1(1)
   iy_1                                      = loc_1(0)
   print("Vredepeel: (iy,ix) = (" +iy_1 + ", " + ix_1 + ")") 

   lat_2                                     = 52.134596 ; Zegveld-De Meije
   lon_2                                     =  4.808560 ; Zegveld-De Meije
   loc_2                                     = wrf_user_latlon_to_ij(fileID, lat_2, lon_2)
   ix_2                                      = loc_2(1)
   iy_2                                      = loc_2(0)
   print("Zegveld: (iy,ix) = (" +iy_2 + ", " + ix_2 + ")") 

   lat_3                                     = 52.104891 ; Wekerom
   lon_3                                     =  5.723901 ; Wekerom
   loc_3                                     = wrf_user_latlon_to_ij(fileID, lat_3, lon_3)
   ix_3                                      = loc_3(1)
   iy_3                                      = loc_3(0)
   print("Wekerom: (iy,ix) = (" +iy_3 + ", " + ix_3 + ")") 

 ; --- Input ---
 ; Observations are from 2016-05-09 00:00h to 2006-05-12 00:00h
    NOx_Vredepeel_obs = asciiread("./ConcentrationData/Vredepeel_NOX.csv",-1,"float")
    NOx_Zegveld_obs = asciiread("./ConcentrationData/Zegveld_NOX.csv",-1,"float")
    NOx_Wekerom_obs = asciiread("./ConcentrationData/Wekerom_NOX.csv",-1,"float") 
   times                                     = wrf_times_c(fileID->Times,0)
   NOx                                       = wrf_user_getvar(fileID,"NOx", -1)             ; NOx        (nz,ny,nx) ppbv (95% * NO@30 kg/mol and 5% NO2 @46 kg/mol)
   LND_tracer                                = wrf_user_getvar(fileID,"LAND_tracer", -1)*1000; LND_tracer (nz,ny,nx) ppmv 
   LND_tracer@units                          = "ppbv"
   SEA_tracer                                = wrf_user_getvar(fileID,"SEA_tracer" , -1)*1000; SEA_tracer (nz,ny,nx) ppmv 
   SEA_tracer@units                          = "ppbv"
   ntimes                                    = dimsizes(times)                               ; number of times in the file
   mdims                                     = getfilevardimsizes(fileID,"P")                ; get some dimension sizes for the file
   nd                                        = dimsizes(mdims)                               ; 0:nt=73; 1:nz=30; 2:ny=30; 3:nx=30;
   time                                      = ispan(0,ntimes-1,1)
   
   ; print (times(0))

 ; --- Processing ----
   ntimes                                    = dimsizes(times)
 ; NOx                                       = (NOx-0.0115)*1000*1.35                        ; ppb --> ug/m3 with background correction
   NOx                                       = (NOx-0.0000)*1000*1.35                        ; ppb --> ug/m3 without background correction
   NOx@units                                 = "ug/m3"
   NOx_ts_1                                  = NOx(:,iz,iy_1,ix_1)
   NOx_ts_2                                  = NOx(:,iz,iy_2,ix_2)
   NOx_ts_3                                  = NOx(:,iz,iy_3,ix_3)
   LND_ts_1                                  = LND_tracer(:,iz,iy_1,ix_1)
   LND_ts_2                                  = LND_tracer(:,iz,iy_2,ix_2)
   LND_ts_3                                  = LND_tracer(:,iz,iy_3,ix_3)
   SEA_ts_1                                  = SEA_tracer(:,iz,iy_1,ix_1)
   SEA_ts_2                                  = SEA_tracer(:,iz,iy_2,ix_2)
   SEA_ts_3                                  = SEA_tracer(:,iz,iy_3,ix_3)

 ; --- Output ----
   opts                                      = True
 ; opts@xyLineColor                          = "blue"
 ; opts@xyLineThicknessF                     = 5   
   opts@pmLegendDisplayMode                  = "Always"                                      ; turn on legend
   opts@pmLegendSide                         = "Top"                                         ; Change location of 
   opts@pmLegendParallelPosF                 = 0.25                                          ; positive moves units right
   opts@pmLegendOrthogonalPosF               = -0.35                                         ; negative moves legend down
   opts@pmLegendWidthF                       = 0.15                                          ; Change width and
   opts@pmLegendHeightF                      = 0.18                                          ; height of legend.
   opts@lgPerimOn                            = False                                         ; turn off box around
   opts@lgLabelFontHeightF                   = .02                                           ; label font height
   opts@tmYLFormat                           = "f"                                           ; remove trailing ".0"                                          
   opts@tmXTOn                               = True
   opts@tmXUseBottom                         = True
   opts@tmXBLabelFontHeightF                 = 0.01
   opts@tmXBLabelAngleF                      = 45
 
 ; --------------------------------------------------
 ; The time_axis_label function adds additional
 ; resources to "res" to produce nicely-formatted 
 ; time labels on X axis. This function only works 
 ; if you have a time "units" recognized by the 
 ; cd_calendar function.
 ; --------------------------------------------------
   
 
   tickopts                                  = True
   tickopts@ttmAxis                          = "XB" ; x-axis bottom
   tickopts@ttmFormat                        = "%D-%c %H:%M"
 ; plotopts@ttmNumticks                      = 7
   tickopts@ttmMajorStride                   = 6
   tickopts@ttmMinorStride                   = 2
   
   time_axis_labels(times,opts,tickopts)
 ; time_axis_labels(time,  resplot, restick)
   
   opts_NOx                                  = opts
   opts_NOx@tiXAxisString                    = "time into simulation (h)"
   opts_NOx@tiYAxisString                    = "NOx (ug/m3)"
   opts_NOx@xyExplicitLegendLabels           = (/" Wekerom (obs)", " Wekerom (WRF)", " Zegveld (obs)", " Zegveld (WRF)", " Vredepeel (obs)", "Vredepeel (WRF)" /)        ; create explicit labels; inverse order
   opts_NOx@xyLineColors                     = (/"blue", "blue", "green", "green", "red", "red" /)                                       ; change line color; inverse order
   opts_NOx@xyMonoDashPattern                = False
   opts_NOx@xyDashPatterns                   = (/ 16, 0, 16, 0, 16, 0 /)                      ; inverse order
   opts_NOx@xyLineThicknesses                = (/ 4.0, 6.0, 4.0, 6.0, 4.0, 6.0/)              ; line thickeness; inverse order
   
   opts_LNDSEA                               = opts
   opts_LNDSEA@tiXAxisString                 = "time into simulation (h)"
   opts_LNDSEA@tiYAxisString                 = "SEA tracer @Land stations (ppbv)"
   opts_LNDSEA@xyExplicitLegendLabels        = (/"SEA Wekerom","SEA Zegveld","SEA Vredepeel"/)  ; create explicit labels
   opts_LNDSEA@xyLineColors                  = (/"green", "red",  "blue"/)                    ; change line color
   opts_LNDSEA@xyMonoDashPattern             = False
   opts_LNDSEA@xyDashPatterns                = (/ 0, 0, 0 /)                                  ; Make 1st 3 curves solid, last 3 dashed
   opts_LNDSEA@xyLineThicknesses             = (/  5.0, 5.0, 5.0, 5.0, 5.0, 5.0/)             ; line thickeness

   print("Discarding the first 24 hours of the run as spin-off...")
   its                                       = 25
   ite                                       = ntimes-1
   print(dimsizes(NOx_Vredepeel_obs))
   NOx_obs_1                                 = runave(NOx_Vredepeel_obs( its:ite),5,0)
   NOx_obs_2                                 = runave(NOx_Zegveld_obs(its:ite),5,0)
   NOx_obs_3                                 = runave(NOx_Wekerom_obs(its:ite),5,0)
   print (dimsizes(NOx_ts_1))
   NOx_mod_1                                 = NOx_ts_1(its:ite)
   NOx_mod_2                                 = NOx_ts_2(its:ite)
   NOx_mod_3                                 = NOx_ts_3(its:ite)
   
   print("Vredepeel: (NOx,obs,NOx,mod) = " + avg(NOx_obs_1) + "," + avg(NOx_mod_1) + ")" ); dNOX = 13.4 ug/m3 --> 0.0099 ppm |
   print("Zegveld: (NOx,obs,NOx,mod) = " + avg(NOx_obs_2) + "," + avg(NOx_mod_2) + ")" ); dNOx = 17.7 ug/m3 --> 0.0131 ppm |
   print("Wekerom: (NOx,obs,NOx,mod) = " + avg(NOx_obs_3) + "," + avg(NOx_mod_3) + ")" ); dNOx =  3.6 ug/m3 --> 0.0026 ppm |--> avg 0.0085 ppm too high
   print(SEA_ts_1)
   plot                                      = gsn_csm_xy(wks, times(its:ite), (/SEA_ts_3(its:ite), SEA_ts_2(its:ite), SEA_ts_1(its:ite)/), opts_LNDSEA)
   plot                                      = gsn_csm_xy(wks, times(its:ite), (/NOx_obs_3, NOx_mod_3, NOx_obs_2, NOx_mod_2, NOx_obs_1, NOx_mod_1 /), opts_NOx)  ; inverse order
   print("Figures written to plt_TimeSeries_NOx.pdf.")
   print("Observed concentrations are in $HOME/WRF/ConcentrationData/LML2006_05.xls, sheet GRAPHS")
 end
