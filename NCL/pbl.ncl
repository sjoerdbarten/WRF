load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "WRFOptions.ncl"     ; set basic plot options here
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
  
  nzlevels=45
  ntimes=66
  Ricrit=0.25
  ; calculate RI profiles for middle domain
  yloc=15
  xloc=15

  Ri=new((/nzlevels/),"float")
  pblh= new((/ntimes/),"float")
  time=new((/ntimes/),"float")
  
  
; load WRF model forecasts
 a = addfile("./wrfout_d01_2018-06-05_00:00:00.nc","r")
  
  colorplot=new(1,logical)
  
  colorplot= True

  wks_type="ps"
  wks_type@wkOrientation= "landscape"
  wks = gsn_open_wks(wks_type ,"WRFpblheight")
  
  do kk=0,ntimes-1
    time(kk)=kk
  end do
  
  do it = 0, ntimes-1
  p1   = wrf_user_getvar(a, "pres",it)     ; pressure is our vertical coordinate
  z   = wrf_user_getvar(a, "height",it)     ; grid point height
  qv1  = wrf_user_getvar(a,"QVAPOR",it)*1000. ; Qv      
  tc1  = wrf_user_getvar(a,"tc",it)     ; T in C
  u1  = wrf_user_getvar(a,"U",it)     ; 
  v1  = wrf_user_getvar(a,"V",it)     ; 
  

  theta1=(tc1+273.15)*(p1/p1(0,yloc,xloc))^-0.285-273.15

 
  do iz=0,nzlevels-1
  Ri(iz)=9.81/300.*z(iz,xloc,yloc)*(theta1(iz,xloc,yloc)-theta1(0,xloc,yloc))/(v1(iz,xloc,yloc)^2+u1(iz,xloc,yloc)^2)
  end do

  iz=0
; Find PBL depth
  do while (Ri(iz).lt.Ricrit)
    iz=iz+1
  end do
  pblh(it)=z(iz,xloc,yloc)-(z(iz,xloc,yloc)-z(iz-1,xloc,yloc))*(Ri(iz)-Ricrit)/(Ri(iz)-Ri(iz-1))

  end do
  pblh2=(/time,pblh/)
  
;Edit MS 19-3-2014  
;New: 
  
  matrixtoascii=(/time,pblh/)
  
 
  ftmx="21f12.4"
  opt=True
  opt@tspace=40
 
  opt@title="Time,PBLheight(m)"
  opt@fout="pblheight.txt"
  write_matrix(transpose(matrixtoascii),ftmx,opt)
  print("***Successfully finished pbl.ncl script***")


;Orig:  
;  asciiwrite("pblheight.txt",(/pblh2/))
;End Edit MS 19-3-2014
  resources                        = True
  resources@xyLineThicknesses   = (/4.,4.,4./) ; Define line thicknesses
;                                               ; (1.0 is the default).
  resources@tiMainString    = "PBL height"  ; Title for the XY plot
  resources@tiXAxisString   = "time (h)"    ; Label for the X axis
  resources@tiYAxisString   = "pbl height (m)"    ; Label for the Y axis
  resources@tiMainFont      = "Helvetica" ; Font for title
  resources@tiXAxisFont     = "Helvetica" ; Font for X axis label
  resources@tiYAxisFont     = "Helvetica" ; Font for Y axis label

  resources@xyMarkLineModes = (/"Lines","Lines","Lines","Markers"/)
  resources@xyMarkers       = (/0,0,0,2/)   ; (none, dot, asterisk)
  resources@xyMarkerColor   = 255           ; Marker color
  resources@xyMarkerSizeF   = 0.03        ; Marker size (default
                                          ; is 0.01)

  plot = gsn_xy(wks,time,pblh,resources)

end
