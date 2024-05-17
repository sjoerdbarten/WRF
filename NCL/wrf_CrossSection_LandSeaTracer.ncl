;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   Plot data on a cross section
;   This script will plot data at a set angle through a specified point
;   Michiel van der Molen, Wageningen U, March 2017

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
 ;
 ; The WRF ARW input file.  
 ; This needs to have a ".nc" appended, so just do it.
   fileID                                    = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")

 ; What kind of plots do you prefer?
 ; type                                      = "x11"
 ; type                                      = "ps"
 ; type                                      = "ncgm"
   type                                      = "newpdf"
   wks                                       = gsn_open_wks(type,"plt_CrossSection_LANDSEA")      ; open WorkStation on which to draw graphics

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
      opts@TimeLabel                          = times(it)                                         ; Set Valid time to use on plots

    ; First get the variables you will need        
      SEA_tracer                              = wrf_user_getvar(fileID,"SEA_tracer", it) * 1000   ; SEA_tracer        (nz,ny,nx) ppm --> ppb
      SEA_tracer@units                        = "ppbv"
      LND_tracer                              = wrf_user_getvar(fileID,"LAND_tracer",it) * 1000   ; LND_tracer        (nz,ny,nx) ppm --> ppb
      LND_tracer@units                        = "ppbv"

      z                                       = wrf_user_getvar(fileID,"z",it)                    ; grid point height (nz,ny,nx)
      izmax                                   = 18                                               ; ~ 2 km
      LND_tracer_adjusted                     = LND_tracer(0:izmax+1,:,:)
      SEA_tracer_adjusted                     = SEA_tracer(0:izmax+1,:,:)
      z_adjusted                              = z(0:izmax+1,:,:)
      if ( FirstTime ) then                                                                       ; get height info for labels
         zmin                                 = 0.
         zmax                                 = max(z_adjusted(izmax,:,:))/1000.                  ; set maximum height at 5km
         nz                                   = floattoint(zmax*2 + 1)                            ; number of ticklabels
         FirstTime                            = False
      end if 
  
  ; ---------------------------------------------------------------
  ; Plot a cross session that runs West-East through the middle of the plot 
  ; For this we need a pivot point and a angle

  ;                   |
  ;       angle=0 is  |     angle = 90 = -----
  ;                   |
  ; 
      angle                               = 90
      plane                               = new(2,float)
      plane                               = (/ mdims(nd-1)/2, mdims(nd-2)/2 /)                 ; pivot point is center of domain (x,y); mdims(nd-1)=x; mdims(nd-2)=y
      opts_intrp                          = False                                              ; start and end points not specified
      SEA_plane                           = wrf_user_intrp3d(SEA_tracer_adjusted,z_adjusted,"v",plane,angle,opts_intrp)
      LND_plane                           = wrf_user_intrp3d(LND_tracer_adjusted,z_adjusted,"v",plane,angle,opts_intrp)
      dim                                 = dimsizes(SEA_plane)                                ; Find the data span - for use in labels
      zspan                               = dim(0)
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
      opts_xy                             = opts
      opts_xy@tiYAxisString               = "Height (km)" 
      opts_xy@tiXAxisString               = "West-East (gridpoint)" 
      opts_xy@cnMissingValPerimOn         = True
      opts_xy@cnMissingValFillColor       = 0
      opts_xy@cnMissingValFillPattern     = 11
      opts_xy@tmYLMode                    = "Explicit"
      opts_xy@tmYLValues                  = fspan(0,zspan,nz)                                  ; Create tick marks
      opts_xy@tmYLLabels                  = sprintf("%.1f",fspan(zmin,zmax,nz))                ; Create labels
      opts_xy@tiXAxisFontHeightF          = 0.020
      opts_xy@tiYAxisFontHeightF          = 0.020
      opts_xy@tmXBMajorLengthF            = 0.02
      opts_xy@tmYLMajorLengthF            = 0.02
      opts_xy@tmYLLabelFontHeightF        = 0.015
      opts_xy@PlotOrientation             = SEA_plane@Orientation 

    ; Plotting options for SEA_tracer
      opts_SEA                            = opts_xy
      opts_SEA@ContourParameters          = (/ 0., 3.0, 0.1 /)
      opts_SEA@cnInfoLabelOrthogonalPosF  = -0.10
      opts_SEA@cnFillOn                   = True
      opts_SEA@cnLinesOn		    = False
      opts_SEA@MainTitle                  = "SEA_tracer"
              
      opts_LND                            = opts_xy
      opts_LND@ContourParameters          = (/ 0., 3.0, 0.1 /)
      opts_LND@cnInfoLabelOrthogonalPosF  = -0.10
      opts_LND@cnFillOn                   = True
      opts_SEA@cnLinesOn		    = False
      opts_LND@MainTitle                  = "LAND_tracer"
 
    ; Get the contour info for the rh and temp
      contour_SEA                         = wrf_contour(fileID,wks,SEA_plane,opts_SEA)
      contour_LND                         = wrf_contour(fileID,wks,LND_plane,opts_LND)

    ; MAKE PLOTS         
      plot                                = wrf_overlays(fileID,wks,(/contour_SEA/),pltopts)
      plot                                = wrf_overlays(fileID,wks,(/contour_LND/),pltopts)

    ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   end do        ; END OF TIME LOOP
   print("Figures are written to plt_CrossSection_LANDSEA.pdf")
end
