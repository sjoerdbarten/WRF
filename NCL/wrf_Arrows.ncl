
;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   Plot data on a cross section
;   This script will plot data from a a given point A to point B
;   Vertical coordinate is pressure

;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;load "WRFUserARW.ncl"



load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"


begin
;
; The WRF ARW input file.  
; This needs to have a ".nc" appended, so just do it.
  a = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")

; We generate plots, but what kind do we prefer?
  type = "newpdf"
  wks = gsn_open_wks(type,"Arrows")



;MAKE SURE THE STARTING POINT IS ALWAYS LESS THAN THE END POINT!!


;FILL IN THINGS HERE: ==========================================================

xloc = 12
yloc = 18

angle = 90.  ;keep angle between 0 and 180.

 ;                   |
 ;                   |
 ;       angle=0 is  |, angle=90 is ------
 ;                   |
 ;                   |

;END FILL IN THINGS ============================================================


plane = new(2,float)
plane = (/ xloc-1, yloc-1 /)    ; pivot point is center of domain

print("center point: xloc = " + xloc + " , yloc = " + yloc )





gsn_define_colormap(wks,"ncl_default")
gsn_draw_colormap(wks)


; Set some Basic Plot options
  ARWres = True
  ARWres@MainTitle          = "REAL-TIME WRF"
  ARWres@vpWidthF           = .9               ; overwrite basic plot size
  ARWres@vpHeightF          = 1.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; What times and how many time steps are in the data set?
  FirstTime = True
  FirstTimeMap = True
  times  = wrf_user_getvar(a,"times",-1) ;wrf_user_list_times(a)  ; get times in the file
  ntimes = dimsizes(times)         ; number of times in the file

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


xlat = wrf_user_getvar(a,"XLAT",0)
xlon = wrf_user_getvar(a,"XLONG",0)
ter = wrf_user_getvar(a,"HGT",0)

 

  do it = 0,ntimes-1,3               ; TIME LOOP

    print("Working on time: " + times(it) )
    ARWres@TimeLabel = times(it)   ; Set Valid time to use on plots

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need        

    th = wrf_user_getvar(a,"th",it)        ; theta
    rh = wrf_user_getvar(a,"rh",it)        ; relative humidity
    z  = wrf_user_getvar(a, "z",it)        ; grid point height
    p  = wrf_user_getvar(a, "pressure",it) ; grid point height
    u  = wrf_user_getvar(a,"ua",it)    ; 3D U at mass points
    v  = wrf_user_getvar(a,"va",it)    ; 3D V at mass points
    w  = wrf_user_getvar(a,"wa",it)    ; 3D V at mass points

    if (FirstTime) then
       zmin = 0.
       zmax = 5.
       nz = floattoint(zmax+1)
    end if


dimsth = dimsizes(th)
nlev = dimsth(0)
Utot = new((/dimsth(0),dimsth(1),dimsth(2)/),float)
angle_wind = new((/dimsth(0),dimsth(1),dimsth(2)/),float)
new_angle = new((/dimsth(0),dimsth(1),dimsth(2)/),float)

do lev = 0,(nlev-1)
    Utot(lev,:,:) = sqrt(u(lev,:,:)^2+v(lev,:,:)^2)
    angle_wind(lev,:,:) = (360/(2*3.141592654))*atan2(u(lev,:,:),v(lev,:,:))+180
    new_angle(lev,:,:) = 180 + angle_wind(lev,:,:) - angle
end do

    U_plane_axis = Utot * cos(new_angle*3.141592654/180)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




    do ip = 1, 1      	; we are just going to do one crossection for now

  ;   Specify the start and end points of the cross section              
  ;   For this case, set opts in wrf_user_intrp3d to True
        opts = False


if ((angle.ge.0).AND.(angle.lt.45)) then
   X_plane = wrf_user_intrp2d(xlat,plane,angle,opts)
   X_desc = "latitude"
else if ((angle.lt.135)) then
   X_plane = wrf_user_intrp2d(xlon,plane,angle,opts)
   X_desc = "longitude"
else if ((angle.ge.135).AND.(angle.le.180)) then
   X_plane = wrf_user_intrp2d(xlat,plane,angle,opts)
   X_desc = "latitude"
end if
end if
end if



        ; Interpolate data vertically (in z)

        rh_plane = wrf_user_intrp3d(rh,z,"v",plane,angle,opts)
        th_plane = wrf_user_intrp3d(th,z,"v",plane,angle,opts)
        U_plane = wrf_user_intrp3d(U_plane_axis,z,"v",plane,angle,opts)
        w_plane = wrf_user_intrp3d(w,z,"v",plane,angle,opts)	


    if (FirstTime) then
      ; THIS IS NEEDED FOR Y LABLES - ALWAYS DO
      zz = wrf_user_intrp3d(z,z,"v",plane,angle,False)
      b = ind(zz(:,0) .gt. zmax*1000. )
      zmax_pos = b(0) -1
      if (abs(zz(zmax_pos,0)-zmax*1000.) .lt. abs(zz(zmax_pos+1,0)-zmax*1000.) ) then
         zspan = b(0) - 1
      else
         zspan = b(0)
      end if
      delete(zz)
      delete(b)
   

      print("Y labels set to: 0 - " + zmax + " km")
      FirstTime = False
      ; END OF ALWAYS DO
    end if

      ; X-axis lables
      dimsX = dimsizes(X_plane)
      xmin  = X_plane(0)
      xmax  = X_plane(dimsX(0)-1)
      print(xmax)
      print(xmin)
      xspan = dimsX(0)-1
      nx    = floattoint( (xmax-xmin)/2 + 1)



  ; Options for XY Plots
        opts_xy                         = ARWres
        opts_xy@tiXAxisString           = X_desc
        opts_xy@tiYAxisString           = "Height (km)"
        opts_xy@cnMissingValPerimOn     = True
        opts_xy@cnMissingValFillColor   = 0
        opts_xy@cnMissingValFillPattern = 11
	
	opts_xy@tmXTOn                  = False
        opts_xy@tmYROn                  = False
        opts_xy@tmXBMode                = "Explicit"
	print(xspan)
	print(nx)
        opts_xy@tmXBValues              = fspan(0,xspan,10)                    ; Create tick marks
        opts_xy@tmXBLabels              = sprintf("%.1f",fspan(xmin,xmax,10))  ; Create labels
        opts_xy@tmXBLabelFontHeightF    = 0.015


        opts_xy@tmYLMode                = "Explicit"
        opts_xy@tmYLValues              = fspan(0,zspan,nz)
	opts_xy@tmYLLabels              = sprintf("%.1f",fspan(zmin,zmax,nz))  ; Create labels
        opts_xy@tiXAxisFontHeightF      = 0.020
        opts_xy@tiYAxisFontHeightF      = 0.020
        opts_xy@tmXBMajorLengthF        = 0.02
        opts_xy@tmYLMajorLengthF        = 0.02
        opts_xy@tmYLLabelFontHeightF    = 0.015
        opts_xy@PlotOrientation         = th_plane@Orientation





  ; Plotting options for RH
        opts_rh = opts_xy
        opts_rh@ContourParameters       = (/ 10., 90., 10. /)
        opts_rh@cnFillOn                = True
        opts_rh@cnFillColors            = (/"White", \
                                            "White","DarkOliveGreen1","Chartreuse","Green", \
                                            "Green1","Green3","Green4", \
                                            "ForestGreen","PaleGreen4"/)

  ; Plotting options for Temperature
        opts_th = opts_xy
        opts_th@ContourParameters       = (/ 286.,308.,1. /)
        opts_th@gsnSpreadColors         = True
        opts_th@cnFillOn                = True


  ; Plotting options for Wind Vectors                 
       opts_wind = opts_xy         
       opts_wind@NumVectors = 42           ; density of wind barbs
       ;opts_wind@vcRefMagnitudeF = 1.
       ;opts_wind@vcRefLengthF = 0.045
       opts_wind@vcGlyphStyle = "LineArrow"
       opts_wind@vcLineArrowHeadMinSizeF = 0.005
       opts_wind@vcLineArrowHeadMaxSizeF = 0.005


  ; Get the contour info for the rh and temp
        contour_th = wrf_contour(a,wks,th_plane(0:zmax_pos,:),opts_th)
        contour_rh = wrf_contour(a,wks,rh_plane(0:zmax_pos,:),opts_rh)
	wrf_smooth_2d(w_plane(0:zmax_pos,:),3)
        windvector = wrf_vector(a,wks,U_plane(0:zmax_pos,:),100.*w_plane(0:zmax_pos,:),opts_wind)


  ; MAKE PLOTS         
        if (FirstTimeMap) then
	   lat_plane = wrf_user_intrp2d(xlat,plane,angle,opts)
	   lon_plane = wrf_user_intrp2d(xlon,plane,angle,opts)
	   dimsX = dimsizes(lat_plane)
	   mpres = True
	    mpres@mpDataBaseVersion = "HighRes"
            mpres@mpNationalLineThicknessF = 3.0 ; interior boundaries 
            mpres@mpGeophysicalLineThicknessF = 3.0 ; lines separating land/ocean 
            mpres@mpGeophysicalLineColor = "Black"

	   pltres = True
	   pltres@FramePlot = False
	   optsM = ARWres
	   optsM@NoHeaderFooter = True
	   optsM@cnFillOn = True
	   optsM@lbTitleOn = False
	   contour = wrf_contour(a,wks,ter,optsM)
	   plot = wrf_map_overlays(a,wks,(/contour/),pltres,mpres)
	   lnres = True
	   lnres@gsLineThicknessF = 3.0
	   lnres@gsLineColor = "Red"
	   do ii = 0,dimsX(0)-2
	      gsn_polyline(wks,plot,(/lon_plane(ii),lon_plane(ii+1)/),(/lat_plane(ii),lat_plane(ii+1)/),lnres)
	   end do
	   frame(wks)
	   delete(lon_plane)
	   delete(lat_plane)
	   pltres@FramePlot = True
	end if   

        plot = wrf_overlays(a,wks,(/contour_th,windvector/),pltres)


  ; Delete options and fields, so we don't have carry over
        delete(opts_th)
        delete(opts_rh)
        delete(opts_wind)
        delete(th_plane)
        delete(rh_plane)
	delete(U_plane)

    end do  ; make next cross section

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    FirstTimeMap = False
  end do        ; END OF TIME LOOP

end
