
;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   Interpolating to specified pressure levels

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
;
; The WRF ARW input file.  
; This needs to have a ".nc" appended, so just do it.
  a = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")

; We generate plots, but what kind do we prefer?
;  type = "x11"
 type = "newpdf"
; type = "ps"
; type = "ncgm"
  wks = gsn_open_wks(type,"plt_PressureLevel1")

; Set some Basic Plot options
  ARWres = True
  ARWres@MainTitle                   = "REAL-TIME WRF"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; What times and how many time steps are in the data set?
  times  = wrf_user_list_times(a)  ; get times in the file
  ntimes = dimsizes(times)         ; number of times in the file

; The specific pressure levels that we want the data interpolated to.
  pressure_levels = (/ 850., 700., 500., 300./)   ; pressure levels to plot
  nlevels         = dimsizes(pressure_levels)     ; number of pressure levels

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  do it = 0,ntimes-1,2             ; TIME LOOP

    print("Working on time: " + times(it) )
    ARWres@TimeLabel = times(it)   ; Set Valid time to use on plots

    mpres = True                   ; Create map background
    map = wrf_map(wks,a,mpres)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need        

    tc = wrf_user_getvar(a,"tc",it)        ; T in C
    u  = wrf_user_getvar(a,"ua",it)        ; u averaged to mass points
    v  = wrf_user_getvar(a,"va",it)        ; v averaged to mass points
    p  = wrf_user_getvar(a, "pressure",it) ; pressure is our vertical coordinate
    z  = wrf_user_getvar(a, "z",it)        ; grid point height
    rh = wrf_user_getvar(a,"rh",it)        ; relative humidity

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    do level = 0,nlevels-1                 ; LOOP OVER LEVELS

      pressure = pressure_levels(level)

      tc_plane = wrf_user_intrp3d(tc,p,"h",pressure,0.,False)
      z_plane  = wrf_user_intrp3d( z,p,"h",pressure,0.,False)
      rh_plane = wrf_user_intrp3d(rh,p,"h",pressure,0.,False)
      u_plane  = wrf_user_intrp3d( u,p,"h",pressure,0.,False)
      v_plane  = wrf_user_intrp3d( v,p,"h",pressure,0.,False)

      spd     = (u_plane*u_plane + v_plane*v_plane)^(0.5) ; m/sec
      spd@description = "Wind Speed"
      spd@units = "m/s"
      u_plane = u_plane*1.94386     ; kts
      v_plane = v_plane*1.94386     ; kts
      u_plane@units = "kts"
      v_plane@units = "kts"


      ; Plotting options for T                
        opts = ARWres                          
        opts@cnLineColor = "Red"
        opts@ContourParameters = (/ 5.0 /)
        opts@cnInfoLabelOrthogonalPosF = 0.07  ; offset second lable information
        opts@gsnContourLineThicknessesScale = 2.0
        contour_tc = wrf_contour(a,wks,tc_plane,opts)
        delete(opts)


      ; Plotting options for RH                
        opts = ARWres                          
        opts@cnFillOn = True  
        opts@ContourParameters = (/ 10., 90., 10./)
        opts@cnFillColors = (/"White","White","White", \
                              "White","Chartreuse","Green",\
                              "Green3","Green4", \
                              "ForestGreen","PaleGreen4"/)
        contour_rh = wrf_contour(a,wks,rh_plane,opts)
        delete(opts)


      ; Plotting options for Wind Speed                
        opts = ARWres                          
        opts@cnLineColor = "MediumSeaGreen"
        opts@ContourParameters = (/ 10. /)
        opts@cnInfoLabelOrthogonalPosF = 0.07  ; offset second lable information
        opts@gsnContourLineThicknessesScale = 3.0
        contour_spd = wrf_contour(a,wks,spd,opts)
        delete(opts)


      ; Plotting options for Wind Vectors                 
        opts = ARWres          
        opts@FieldTitle = "Wind"   ; overwrite Field Title
        opts@NumVectors = 47       ; wind barb density
        vector = wrf_vector(a,wks,u_plane,v_plane,opts)
        delete(opts)


      ; Plotting options for Geopotential Heigh
        opts_z = ARWres                          
        opts_z@cnLineColor = "NavyBlue"
        opts_z@gsnContourLineThicknessesScale = 3.0


      ; MAKE PLOTS                                       

        if ( pressure .eq. 850 ) then   ; plot temp, rh, height, wind barbs
          opts_z@ContourParameters = (/ 20.0 /)
          contour_height = wrf_contour(a,wks,z_plane,opts_z)
          wrf_map_overlay(wks,map,(/contour_rh,contour_tc,contour_height, \
                                    vector/),True)
        end if

        if ( pressure .eq. 700 ) then   ; plot temp, height, wind barbs
          opts_z@ContourParameters = (/ 30.0 /)
          contour_height = wrf_contour(a,wks, z_plane,opts_z)
          wrf_map_overlay(wks,map,(/contour_tc,contour_height, \
                                    vector/),True)
        end if

        if ( pressure .eq. 500 ) then   ; plot temp, height, wind barbs
          opts_z@ContourParameters = (/ 60.0 /)
          contour_height = wrf_contour(a,wks, z_plane,opts_z)
          wrf_map_overlay(wks,map,(/contour_tc,contour_height, \
                                    vector/),True)
        end if

        if ( pressure .eq. 300 ) then   ; plot windspeed, height, wind barbs
          opts_z@ContourParameters = (/ 60.0 /)
          contour_height = wrf_contour(a,wks, z_plane,opts_z)
          wrf_map_overlay(wks,map,(/contour_spd,contour_height, \
                                    vector/),True)
        end if
        delete(opts_z)

    end do      ; END OF LEVEL LOOP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  end do        ; END OF TIME LOOP

end
