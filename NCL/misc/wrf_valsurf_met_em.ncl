;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
;
; The WRF ARW input file.  
; This needs to have a ".nc" appended, so just do it.
  a = addfile("/home/WUR/dubbi001/atmos2020/WRF/run/wrfout_d01_2020-02-07_00:00:00.nc","r")
  b = addfile("/home/WUR/dubbi001/atmos2020/WPS/met_em.d01.2020-02-07_00:00:00.nc","r")

; We generate plots, but what kind do we prefer?
;  type = "x11"
  type = "newpdf"
; type = "ps"
; type = "ncgm"

  wks = gsn_open_wks(type,"plt_Val1")
  gsn_define_colormap(wks,"hotres")

; Set some basic resources
  res = True
  res@MainTitle                   = "REAL-TIME WRF"

  pltres = True
  mpres = True
  mpres@mpDataBaseVersion = "HighRes"
  mpres@mpNationalLineThicknessF = 3.0 ; interior boundaries 
  mpres@mpGeophysicalLineThicknessF = 3.0 ; lines separating land/ocean 
  mpres@mpGeophysicalLineColor = "Black"


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; What times and how many time steps are in the data set?
  times  = wrf_user_list_times(a)  ; get times in the file
  ntimes = dimsizes(times)         ; number of times in the file

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


 do it = 0,2 ; ntimes-1,1            ; TIME LOOP

    print("Working on time: " + times(it) )
    res@TimeLabel = times(it)   ; Set Valid time to use on plots


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need        

  ;  slp_met_em = b->msl(it,::-1,:)
  ;  lon_era5 = b->longitude
  ;  lat_era5 = b->latitude(::-1)
  ;  lon_wrf = a->XLONG(it,:,:)
  ;  lat_wrf = a->XLAT(it,:,:)
   
 ;   reg_wrf = rgrid2rcm_Wrap(lat_era5,lon_era5,slp_met_em,lat_wrf,lon_wrf,0)
    reg_wrf = b->PMSL(0,:,:)
    print(reg_wrf)
    slp = wrf_user_getvar(a,"slp",it)    ; slp
      wrf_smooth_2d( slp, 3 )            ; smooth slp
    if ( it .eq. 0 ) then                
      tc = wrf_user_getvar(a,"tc",it)    ; 3D tc 
      td = wrf_user_getvar(a,"td",it)    ; 3D td 
      u  = wrf_user_getvar(a,"ua",it)    ; 3D U at mass points
      v  = wrf_user_getvar(a,"va",it)    ; 3D V at mass points
      tc2 = tc(0,:,:)                    ; Use lowest T at time zero
      td2 = td(0,:,:)                    ; Use lowest Td at time zero
      u10 = u(0,:,:)                     ; Use lowest level at time 0
      v10 = v(0,:,:)
      rh2x = wrf_user_getvar(a,"rh2",it)  ; RH at 2m, mass point
    else
      tc2 = wrf_user_getvar(a,"T2",it)   ; T2 in Kelvin
         tc2 = tc2-273.16                  ; T2 in C
	 tc2@units = "C"
      td2 =  wrf_user_getvar(a,"td2",it) ; Td2 in C
      u10 = wrf_user_getvar(a,"U10",it)  ; u at 10 m, mass point
      v10 = wrf_user_getvar(a,"V10",it)  ; v at 10 m, mass point
      rh2x = wrf_user_getvar(a,"rh2",it)  ; RH at 2m, mass point
    end if

    tf2 = 1.8*tc2+32.                    ; Turn temperature into Fahrenheit
      tf2@description = "Surface Temperature"
      tf2@units = "F"
    td_f = 1.8*td2+32.                   ; Turn temperature into Fahrenheit
      td_f@description = "Surface Dew Point Temp" 
      td_f@units = "F"
    u10 = u10*1.94386                ; Turn wind into knots
    v10 = v10*1.94386
      u10@units = "kts"
      v10@units = "kts"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; Plotting options for T                
      opts = res                         
      opts@cnFillOn = True  
     ; opts@ContourParameters = (/ 8., 35., 1./)
      opts@gsnSpreadColorEnd = -3  ; End third from the last color in color map
      contour_tc = wrf_contour(a,wks,tc2,opts)
      delete(opts)


    ; Plotting options for Td
      opts = res         
      opts@cnFillOn = True 
      opts@cnLinesOn = True
      opts@cnLineLabelsOn = True
      opts@ContourParameters = (/ 5., 30., 1./) 
      opts@cnLineLabelBackgroundColor = -1
      opts@gsnSpreadColorEnd = -3  ; End third from the last color in color map
      contour_td = wrf_contour(a,wks,td2,opts)
      delete(opts)

    ; Plotting options for RH2
      opts = res         
      opts@cnFillOn = True 
      opts@cnLinesOn = True
      opts@cnLineLabelsOn = True
      opts@ContourParameters = (/ 30., 100., 2.5/) 
      opts@cnLineLabelBackgroundColor = -1
      opts@gsnSpreadColorEnd = -3  ; End third from the last color in color map
      contour_rh = wrf_contour(a,wks,rh2x,opts)
      delete(opts)

    ; Plotting options for SLP                     
      opts = res         
      opts@cnLineColor = "Blue"
      opts@cnHighLabelsOn = True
      opts@cnLowLabelsOn = True
      opts@ContourParameters = (/ 980., 1030., 2.5 /)
      opts@cnLineLabelBackgroundColor = -1
      opts@gsnContourLineThicknessesScale = 2.0
      contour_psl = wrf_contour(a,wks,reg_wrf,opts)
      delete(opts)

    ; Plotting options for Wind Vectors                 
      opts = res         
      opts@FieldTitle = "Wind"       ; overwrite Field Title
      opts@NumVectors = 47           ; density of wind barbs
      vector = wrf_vector(a,wks,u10,v10,opts)
      delete(opts)
  

    ; MAKE PLOTS                                       
      plot = wrf_map_overlays(a,wks,(/contour_tc,contour_psl,vector/),pltres,mpres)
     ; plot = wrf_map_overlays(a,wks,(/contour_rh,vector/),pltres,mpres)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  end do        ; END OF TIME LOOP

end
