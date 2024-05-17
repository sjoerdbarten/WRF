
;   Sample script to create and plot PW for WRF-ARW model output
;   November 2009

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
  wks = gsn_open_wks(type,"plt_PW")
  gsn_define_colormap(wks,"gsdtol")

; Set some basic resources
  res = True
  res@MainTitle                   = "REAL-TIME WRF"
  res@Footer = False

  pltres = True
  mpres = True

  gas_const = 287. ; J/K/kg
  Cp        = 1004. ; J/K/kg

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; What times and how many time steps are in this data set?
  times = wrf_user_getvar(a,"times",-1)  ; get times in the file
  ntimes = dimsizes(times)         ; number of times in the file

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  do it = 0,ntimes-1,3                       ; TIME LOOP

    print("Working on time: " + times(it) )
    res@TimeLabel = times(it)      ; Set Valid time to use on plots

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need        

    pres = wrf_user_getvar(a,"pres",it)                ; pressure in pa
    ph = wrf_user_getvar(a,"PH",it)                    ; pert. geopt
    phb = wrf_user_getvar(a,"PHB",it)                  ; base geopt
    th = wrf_user_getvar(a,"theta",it)                 ; theta (T+300)
    qvapor = wrf_user_getvar(a,"QVAPOR",it)    
  
; Calc height, pressure and virtual temperature
    height = (ph + phb)/9.81                           ; height at full levels
    temp = th * (pres/100000) ^ (gas_const/Cp)
    vtemp = (1 + 0.61*qvapor) * temp                   ; virtual temp
  
    zdiff=height(0,:,:)
    zdiff=0.
    pw_sfc_ptop=height(0,:,:)
    pw_sfc_ptop=0.

    dim = dimsizes(pres)
    do k=0,dim(0)-1
      zdiff(:,:)= (height(k+1,:,:) - height(k,:,:))
      pw_sfc_ptop(:,:) =  pw_sfc_ptop(:,:) + ((pres(k,:,:)/(gas_const * vtemp(k,:,:))) * qvapor(k,:,:) * zdiff(:,:))
    end do
    pw_sfc_ptop@description = "PW"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    opts = res                          
    opts@cnFillOn = True  
    ;opts@lbTitleOn = False
    opts@lbLabelBarOn = True
    ;opts@lbLabelsOn = False
    ;opts@pmLabelBarDisplayMode = "NoCreate"
    contour  = wrf_contour(a,wks,pw_sfc_ptop,opts)
    delete(opts)
   
    ; MAKE PLOTS                                       
    plot = wrf_map_overlays(a,wks,(/contour/),pltres,mpres)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  end do  ;       end of the time loop

end
