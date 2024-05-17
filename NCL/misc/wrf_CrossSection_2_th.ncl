
;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   Plot data on a cross section
;   This script will plot data at a set angle through a specified point
;   Add some label info to the Y-axis

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
;
; The WRF ARW input file.  
; This needs to have a ".nc" appended, so just do it.
  a = addfile("$HOME/WRF/WRFV3/run/wrfout_d02_2006-05-09_00:00:00.nc","r")


; We generate plots, but what kind do we prefer?
;  type = "x11"
  type = "newpdf"
; type = "ps"
; type = "ncgm"
  wks = gsn_open_wks(type,"plt_CrossSection_2_th")


; Set some basic resources
  res = True
  res@MainTitle = "REAL-TIME WRF"
  res@Footer = False
  
  pltres = True


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  times  = wrf_user_getvar(a,"times",-1) ; get times in the file
  ntimes = dimsizes(times)          ; number of times in the file
  FirstTime = True

  mdims = getfilevardimsizes(a,"P") ; get some dimension sizes for the file
  nd = dimsizes(mdims)

;---------------------------------------------------------------

  do it = 0,ntimes-1,3                  ; TIME LOOP

    print("Working on time: " + times(it) )
    res@TimeLabel = times(it)           ; Set Valid time to use on plots

    th  = wrf_user_getvar(a,"th",it)     ; Th in K
;    rh = wrf_user_getvar(a,"rh",it)      ; relative humidity
      z   = wrf_user_getvar(a, "z",it)     ; grid point height
      u  = wrf_user_getvar(a,"ua",it)    ; 3D U at mass points
      v  = wrf_user_getvar(a,"va",it)    ; 3D V at mass points
      w  = wrf_user_getvar(a,"wa",it)    ; 3D V at mass points

      Utot = sqrt(u^2+v^2)


    if ( FirstTime ) then                ; get height info for labels
      zmin = 0.
      zmax = max(z)/1000.
      nz   = floattoint(zmax/2 + 1)
      FirstTime = False
    end if

;---------------------------------------------------------------

; Plot a cross session that run south-north through the middle of the plot
; For this we need a pivot point and a angle

;                   |
;       angle=0 is  |
;                   |
; 

        angle = 0
        plane = new(2,float)
        plane = (/ mdims(nd-1)/2, mdims(nd-2)/2 /)    ; pivot point is center of domain (x,y)
        opts = False                                  ; start and end points not specified

        rh_plane = wrf_user_intrp3d(rh,z,"v",plane,angle,opts)
        tc_plane = wrf_user_intrp3d(tc,z,"v",plane,angle,opts)

        dim = dimsizes(rh_plane)                      ; Find the data span - for use in labels
        zspan = dim(0)

      ; Options for XY Plots
        opts_xy                         = res
        opts_xy@tiYAxisString           = "Height (km)"
        opts_xy@cnMissingValPerimOn     = True
        opts_xy@cnMissingValFillColor   = 0
        opts_xy@cnMissingValFillPattern = 11
        opts_xy@tmYLMode                = "Explicit"
        opts_xy@tmYLValues              = fspan(0,zspan,nz)                    ; Create tick marks
        opts_xy@tmYLLabels              = sprintf("%.1f",fspan(zmin,zmax,nz))  ; Create labels
        opts_xy@tiXAxisFontHeightF      = 0.020
        opts_xy@tiYAxisFontHeightF      = 0.020
        opts_xy@tmXBMajorLengthF        = 0.02
        opts_xy@tmYLMajorLengthF        = 0.02
        opts_xy@tmYLLabelFontHeightF    = 0.015
        opts_xy@PlotOrientation         = tc_plane@Orientation
        

      ; Plotting options for RH
        opts_rh = opts_xy
        opts_rh@ContourParameters       = (/ 10., 90., 10. /)
        opts_rh@pmLabelBarOrthogonalPosF = -0.07
        opts_rh@cnFillOn                = True
        opts_rh@cnFillColors            = (/"White","White","White", \
                                            "White","Chartreuse","Green", \
                                            "Green3","Green4", \
                                            "ForestGreen","PaleGreen4"/)

      ; Plotting options for Temperature
        opts_tc = opts_xy
        opts_tc@cnInfoLabelOrthogonalPosF = 0.00
        opts_tc@ContourParameters  = (/ 5. /)




    ; Plotting options for Wind Vectors                 
      opts_wind = opts_xy         
      opts_wind@FieldTitle = "Wind"       ; overwrite Field Title
      opts_wind@NumVectors = 47           ; density of wind barbs
      vector = wrf_vector(a,wks,u10,v10,opts)
      delete(opts)



      ; Get the contour info for the rh and temp
        contour_tc = wrf_contour(a,wks,tc_plane,opts_tc)
        contour_rh = wrf_contour(a,wks,rh_plane,opts_rh)


      ; MAKE PLOTS         
        plot = wrf_overlays(a,wks,(/contour_rh,contour_tc/),pltres)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  end do        ; END OF TIME LOOP

end
