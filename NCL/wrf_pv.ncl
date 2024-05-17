
;   Example script to produce Vorticity plots from WRF ARW model data
;   Novemner 2008

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
;
; The WRF ARW input file.  
; This needs to have a ".nc" appended, so just do it.
  a = addfile("./wrfout_d01_2018-06-05_00:00:00.nc","r")

; We generate plots, but what kind do we prefer?
  type = "x11"
; type = "newpdf"
; type = "ps"
; type = "ncgm"
  wks = gsn_open_wks(type,"plt_pv")

; Set some basic resources
  res = True
  res@MainTitle                   = "REAL-TIME WRF"

  pltres = True
  mpres = True


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; What times and how many time steps are in the data set?
  times  = wrf_user_getvar(a,"times",-1) ; get times in the file
  ntimes = dimsizes(times)         ; number of times in the file

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  do it = 0,ntimes-1,2             ; TIME LOOP

    print("Working on time: " + times(it) )
    res@TimeLabel = times(it)   ; Set Valid time to use on plots


  ; Get the data
    pvo   = wrf_user_getvar(a,"pvo",it)
    avo   = wrf_user_getvar(a,"avo",it)
    p     = wrf_user_getvar(a,"pressure",it)


  ; Interpolate to pressure
    pv_plane = wrf_user_intrp3d(pvo,p,"h",300.,0,False)
    av_plane = wrf_user_intrp3d(avo,p,"h",500.,0,False)


    ; Plotting options 
      opts = res                         
      opts@cnFillOn = True  
      opts@gsnSpreadColorEnd = -3  ; End third from the last color in color map
      opts@ContourParameters = (/ 0., 100., 10./)
      contour_a = wrf_contour(a,wks,av_plane,opts)
      opts@ContourParameters = (/ 0., 5., 1./)
      contour = wrf_contour(a,wks,pv_plane,opts)
      delete(opts)

    ; MAKE PLOTS                                       
      plot = wrf_map_overlays(a,wks,(/contour_a/),pltres,mpres)
      plot = wrf_map_overlays(a,wks,(/contour/),pltres,mpres)


  end do        ; END OF TIME LOOP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end
