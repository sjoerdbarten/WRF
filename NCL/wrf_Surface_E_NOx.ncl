
;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   Plot data on a cross section
;   This script will plot data at a set angle through a specified point

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
 ;
 ; The WRF ARW input file.  
 ; This needs to have a ".nc" appended, so just do it.
   fileID             	              = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")


 ; What kind of plots do you prefer?
 ; type             	              = "x11"
 ; type             	              = "ps"
 ; type             	              = "ncgm"
   type             	              = "newpdf"
   wks                                = gsn_open_wks(type,"plt_Surface_E_NOx")


 ; Set some basic options
   opts              	              = True
   opts@Footer      	              = False
   opts@gsnSpreadColors               = True 
   opts@PlotLevelID                   = "Surface" 
   pltopts         	              = True
   mapopts                            = True
   mapopts@mpDataBaseVersion          ="Highres"
   mapopts@mpNationalLineThicknessF   = 3.0                                                ; interior boundaries 
   mapopts@mpGeophysicalLineThicknessF= 3.0                                                ; lines separating land/ocean 
   mapopts@mpGeophysicalLineColor     = "Black" 

 ; --- get times and dimensions
   times                              = wrf_user_getvar(fileID,"times",-1)                 ; get times in the file
   ntimes                             = dimsizes(times)                                    ; number of times in the file
   mdims                              = getfilevardimsizes(fileID,"P")                     ; get some dimension sizes for the file
   nd                                 = dimsizes(mdims)                                    ; 0:nt=73; 1:nz=30; 2:ny=30; 3:nx=30;

 ; --- Start time loop
   do it = 2,ntimes-1,2                                                                    ; TIME LOOP   
   
    ; --- Set plot options  
      print("Working on time: " + times(it) )
      opts@TimeLabel                  = times(it)                                          ; Set Valid time to use on plots
    
      opts_E_NOx                      = opts
      opts_E_NOx@cnFillOn             = True
      opts_E_NOx@MainTitle            = "E_NOx (emissions in mole/km2/hr)"
      opts_E_NOx@cnLevelSelectionMode = "ExplicitLevels"                                   ; set explicit contour levels 
      opts_E_NOx@cnLevels             = (/ 1.00e-01, 1.78e-01, 3.16e-01, 5.62e-01, 1.00e0, 1.78e0, 3.16e0, 5.62e0, 1.00e1, 1.78e1, 3.16e1, 5.62e1, 1.00e2, 1.78e2, 3.16e2, 5.62e2, 1.00e3, 1.78e3, 3.16e3, 5.62e3 /)

      opts_NOx                        = opts
      opts_NOx@cnFillOn               = True
      opts_NOx@MainTitle              = "NOx (surface concentration in ug/m3)"

      opts_wv                         = opts
      opts_wv@FieldTitle              = "Wind"                                             ; overwrite Field Title
      opts_wv@NumVectors              = 50                                                 ; density of wind barbs
      opts_wv@vcGlyphStyle            ="WindBarb"

    ; --- Get the variables we need        
      u                               = wrf_user_getvar(fileID,"ua",    it)                ; 3D U at mass points
      v                               = wrf_user_getvar(fileID,"va",    it)                ; 3D V at mass points
      E_NOx                           = wrf_user_getvar(fileID,"E_NOx", it)                ; E_NOx      (mole/km2/hr)
      NOx                             = wrf_user_getvar(fileID,"NOx",   it)                ; NOx        (nz,ny,nx) ppmv (95% * NO@30 kg/mol and 5% NO2 @46 kg/mol)
      NOx                             = (NOx-0.0115)*1000                                  ; NOx convert ppmv to ppbv. (2018-03: I do not remember why the offset of 0.0115 ppm is removed.
      NOx@units                       = "ppbv"

    ; --- MAKE PLOTS         
      contour_E_NOx                   = wrf_contour(     fileID,wks,E_NOx(0,:,:),            opts_E_NOx)
      contour_NOx                     = wrf_contour(     fileID,wks,  NOx(0,:,:),            opts_NOx)
      vector                          = wrf_vector(      fileID,wks,    u(0,:,:),v(0,:,:),   opts_wv)
      plot                            = wrf_map_overlays(fileID,wks,  contour_E_NOx,         pltopts,mapopts)
      plot                            = wrf_map_overlays(fileID,wks,(/contour_NOx,vector/),  pltopts,mapopts) 
 
   end do                                                                                  ; END OF TIME LOOP
   print ("Figures written to plt_Surface_E_NOx.pdf.")
   print ("Open PDF in two-page view to compare the emission field with the concentration field.")
end                                                                                        ; end script
