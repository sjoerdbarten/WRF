
;   Example script to produce dbz plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   November 2008

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
  wks = gsn_open_wks(type,"plt_dbz1")
  gsn_define_colormap(wks,"WhViBlGrYeOrReWh")       ; Overwrite the standard color map

; Set some basic resources
  res = True
  res@MainTitle                   = "REAL-TIME WRF"

  pltres = True
  mpres = True
  mpres@mpDataBaseVersion = "HighRes"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Which times and how many time steps are in the data set?
  times = wrf_user_getvar(a,"times",-1)  ; get all times in the file
  ntimes = dimsizes(times)         ; number of times in the file

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  do it = 0,ntimes-1,3               ; TIME LOOP
 
     print("Working on time: " + times(it) )
     res@TimeLabel = times(it)   ; Set Valid time to use on plots


; First get the variables we will need        
; Both dbz and mdbz will be calculated using constant intercept parameters
; (as in early Reisner-2), with values of 8x10^6, 2x10^7, and 4x10^6 m^-4,
; for rain, snow, and graupel, respectively.

     mdbz = wrf_user_getvar(a,"mdbz",it) 
     dbz = wrf_user_getvar(a,"dbz",it) 

     opts = res                         
     opts@cnFillOn = True  
     opts@ContourParameters = (/ 5., 75., 5./)
     contour = wrf_contour(a,wks,dbz(1,:,:),opts)     ; plot only lowest level
     plot = wrf_map_overlays(a,wks,(/contour/),pltres,mpres)
 
     contour = wrf_contour(a,wks,mdbz,opts)
     plot = wrf_map_overlays(a,wks,(/contour/),pltres,mpres)
 

   end do        ; END OF TIME LOOP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
