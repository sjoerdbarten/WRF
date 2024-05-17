;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   Plot data on a cross section
;   This script will plot data at a set angle through a specified point
 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
 
begin
 ; The WRF ARW input file. 
 ; This needs to have a ".nc" appended, so just do it.
   fileID                                    = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")
  
 ; What kind of plots do you prefer?
 ; type                                      = "x11"
   type                                      = "newpdf"
 ; type                                      = "ps"
 ; type                                      = "ncgm"
   wks                                       = gsn_open_wks(type,"plt_CrossSection_NOx")
  
 ; Set some basic resources
   opts                                      = True
   opts@Footer                               = False
   pltopts                                   = True
   FirstTime                                 = True
 
 ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   times                                     = wrf_user_getvar(fileID,"times",-1)                 ; get times in the file
   ntimes                                    = dimsizes(times)                                    ; number of times in the file
   mdims                                     = getfilevardimsizes(fileID,"P")                     ; get some dimension sizes for the file
   nd                                        = dimsizes(mdims)                                    ; 0:nt=73; 1:nz=30; 2:ny=30; 3:nx=30;
 ; ---------------------------------------------------------------
   do it = 2,ntimes-1,2                                                                           ; TIME LOOP
 
      print("Working on time: " + times(it) )
      opts@TimeLabel                         = times(it)                                          ; Set Valid time to use on plots
 
    ; First get the variables yoy will need       
      NOx                                    = wrf_user_getvar(fileID,"NOx", it); *1000           ; NOx        (nz,ny,nx) ppbv (95% * NO@30 kg/mol and 5% NO2 @46 kg/mol)
      NOx = (NOx-0.0115)*1000
      NOx@units                              = "ppbv"
      z                                      = wrf_user_getvar(fileID,"z",it)                     ; grid point height (nz,ny,nx)
      izmax                                  = 18
      NOx_adjusted                           = NOx(0:izmax+1,:,:)
      z_adjusted                             = z(0:izmax+1,:,:)
      if ( FirstTime ) then   
         zmin                                = 0.
      	 zmax                                = max(z_adjusted(izmax,:,:))/1000.	                  ; set maximum height at 2 km
      	 nz                                  = floattoint(zmax*2 + 1)                             ; number of ticklabels
	 FirstTime 		             = False
      end if
    ; ---------------------------------------------------------------
    ; Plot a cross session that runs West-East through the middle of the plot
    ; For this we need a pivot point and a angle
    ;                   |
    ;       angle=0 is  |     angle = 90 = -----
    ;                   |
      angle                                  = 90
      plane                                  = new(2,float)
      plane                                  = (/ mdims(nd-1)/2, mdims(nd-2)/2 /)                 ; pivot point is center of domain (x,y); mdims(nd-1)=x; mdims(nd-2)=y
      opts_intrp                             = False                                              ; start and end points not specified
      NOx_plane                              = wrf_user_intrp3d(NOx_adjusted,z_adjusted,"v",plane,angle,opts_intrp)
      dim                                    = dimsizes(NOx_plane)                                ; Find the data span - for use in labels
      zspan                                  = dim(0)
      if ( FirstTime ) then
         zz                                  = wrf_user_intrp3d(z,z,"v",plane,angle,opts_intrp)
         b                                   = ind(zz(:,0) .gt. zmax*1000. )
         zmax_pos                            = b(0) - 1
         if ( abs(zz(zmax_pos,0)-zmax*1000.) .lt. abs(zz(zmax_pos+1,0)-zmax*1000.) ) then
            zspan                            = b(0) - 1
         else
            zspan                            = b(0)
         end if
         delete(zz)
         delete(b)
         FirstTime = False
      end if
 
    ; Options for XY Plots        
      opts_xy                                = opts
      opts_xy@tiYAxisString                  = "Height (km)"
      opts_xy@tiXAxisString                  = "West-East (gridpoint)"
      opts_xy@cnMissingValPerimOn            = True
      opts_xy@cnMissingValFillColor          = 0
      opts_xy@cnMissingValFillPattern        = 11
      opts_xy@tmYLMode                       = "Explicit"
      opts_xy@tmYLValues                     = fspan(0,zspan,nz)                                   ; Create tick marks
      opts_xy@tmYLLabels                     = sprintf("%.1f",fspan(zmin,zmax,nz))                 ; Create labels
      opts_xy@tiXAxisFontHeightF             = 0.020
      opts_xy@tiYAxisFontHeightF             = 0.020
      opts_xy@tmXBMajorLengthF               = 0.02
      opts_xy@tmYLMajorLengthF               = 0.02
      opts_xy@tmYLLabelFontHeightF           = 0.015
      opts_xy@PlotOrientation                = NOx_plane@Orientation
  
    ; Plotting options for NOx_tracer
      opts_NOx                               = opts_xy
      opts_NOx@ContourParameters             = (/ 0., 50, 2. /)
      opts_NOx@cnInfoLabelOrthogonalPosF     = -0.10
      opts_NOx@cnFillOn                      = True
      opts_NOx@cnLinesOn                     = False
      opts_NOx@MainTitle                     = "NOx concentration"
  
    ; Get the contour info for the rh and temp
      contour_NOx                            = wrf_contour(fileID,wks,NOx_plane,opts_NOx)
    ; MAKE PLOTS        
      plot                                   = wrf_overlays(fileID,wks,(/contour_NOx/),pltopts)
 
   end do        ; END OF TIME LOOP
   print("Figures are written to plt_CrossSection_NOx.pdf")

end  
