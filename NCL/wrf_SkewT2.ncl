
;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   Plot SkewT's at a number of x/y locations

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;load "./WRFUserARW.ncl"

begin
;
; The WRF ARW input file.  
; This needs to have a ".nc" appended, so just do it.
  a = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")

;  We generate plots, but what kind do we prefer?
;  type = "x11"
   type = "newpdf"
; type = "ps"
; type = "ncgm"
  wks = gsn_open_wks(type,"plt_SkewT2")

  gsn_define_colormap(wks,"WhViBlGrYeOrReWh")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  What times and how many time steps are in the data set?

  times  = wrf_user_list_times(a)  ; get times in the file
  ntimes = dimsizes(times)         ; number of times in the file

; This is the big loop over all of the time periods to process.

  FirstTime = True

  do it = 0,ntimes-1
   print("Working on time " + it )
   time = it

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need        

    tc  = wrf_user_getvar(a,"tc",time)            ; T in C
    td  = wrf_user_getvar(a,"td",time)            ; dew point temperature
    p   = wrf_user_getvar(a, "pressure",time)     ; grid point pressure
    z   = wrf_user_getvar(a, "z",time)            ; grid point height
    uvm = wrf_user_getvar(a,"uvmet",time)         ; umet and vmet averaged to mass points
						  ; This is a 4D array where
						  ; uvm(0,:,:,:) is umet, and
						  ; uvm(1,:,:,:) is vmet, and
						  ; This function rotate winds to earth coord.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;       ************************************************************

    u  =  uvm(0,:,:,:)*1.94386
    v  =  uvm(1,:,:,:)*1.94386
    dims_v = dimsizes(v)
    x_str = dims_v(2) / 6
    y_str = dims_v(1) / 6

    do ip = 0, 4   ; correspond to 5 x/y locations below

      loc_x = x_str * (ip + 1)
      loc_y = y_str * (ip + 1)

      ip_locs = "Sounding at x=" +loc_x+ " ; y=" +loc_y+ " ; Time = " +times(it)

  ;      Define a few skew-T plotting options
         skewtOpts          = True
         skewtOpts@DrawHeightScale   = True    ; plot height scale on side
         skewtOpts@DrawHeightScaleFt = False   ; plot height scale in km
         skewtOpts@DrawStandardAtm   = True    ; draw standard atm on plot
         skewtOpts@vpXF              = 0.12    ; controls off-set from left
         skewtOpts@vpYF              = 0.87    ; controls off-set from top
         skewtOpts@vpWidthF          = 0.75    ; controls size of plot
         skewtOpts@vpHeightF         = 0.75    ; controls size of plot
         skewtOpts@DrawFahrenheit    = False   ; use deg C scale
         skewtOpts@tiMainFontHeightF = 0.015   ; change height of main title
        ;skewtOpts@DrawColLine       = False   ; draw lines in black
         skewtOpts@DrawColAreaFill   = True    ; color on background plot
        ;skewtOpts@DrawColAreaColor  = "Green" ; final color may depend on the color table used
         skewtOpts@DrawColAreaColor  = 53      ; Light Green for WhViBlGrYeOrReWh color table
         skewtOpts@PrintOpts         = False   ; do not print options out


  ;      Get the skew-T background
         skewtOpts@tiMainString = ip_locs
         skewt_bkgd             = skewT_BackGround (wks, skewtOpts)
         draw (skewt_bkgd)


  ;      Draw the skew-T plot
         dataOpts           = True
         dataOpts@Parcel    = 1
         dataOpts@WspdWdir  = False  ; wind speed and dir [else: u,v]
         dataOpts@HspdHdir  = True   ; wind speed and dir [else: u,v]
         dataOpts@PlotWindH = False  ; plot wind barbs at h lvls [pibal; special]
         skewT_data = skewT_PlotData(wks, skewt_bkgd, p(:,loc_y, loc_x), \
                                                     tc(:,loc_y, loc_x), \
                                                     td(:,loc_y, loc_x), \
                                                     z(:,loc_y, loc_x), \
                                                     u(:,loc_y, loc_x), \
                                                     v(:,loc_y, loc_x), \
                                                     dataOpts)

  ;      Close the frame
         frame(wks)

         delete(skewtOpts)
         delete(dataOpts)
         delete(skewT_data)
         delete(skewt_bkgd)

    end do


  ;       ************************************************************

   delete(dims_v)
   FirstTime = False

  end do  ;       end of the time loop

end
