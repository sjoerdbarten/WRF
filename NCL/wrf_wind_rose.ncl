;############################################################
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/wind_rose.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

  f = addfile("./wrfout_d01_2018-06-05_00:00:00.nc","r")

  uvm10 = wrf_user_getvar(f,"uvmet10",-1) 
  u10 = uvm10(0,:,:,:)
  v10 = uvm10(1,:,:,:)
                                  
  
; Pick one grid point
  
  loc = wrf_user_ll_to_ij(f,5., 52., True)    ; fill with correct lon ,lat  change  -87.37, 25.65
  mx = loc(0) - 1
  ny = loc(1) - 1

  print(loc)
  lat1  = f->XLAT(0,ny,mx)
  lon1  = f->XLONG(0,ny,mx)
  res               = True
  res@tiMainString  = "WRF: All Times: grid point ["+sprintf("%5.2f",lat1)+" , "+sprintf("%5.2f",lon1) +"]"
  wspd1= ndtooned( sqrt(u10(:,ny,mx)^2 + v10(:,ny,mx)^2) )
  wdir1= ndtooned( atan2(u10(:,ny,mx),v10(:,ny,mx))/0.01745329 +180. )

  numPetals =  16   
  circFr    = 10. 
  spdBounds         = (/ 10., 20., 30., 40. /)
  colorBounds       = (/ "blue", "green", "yellow", "red" /)

  wks = gsn_open_wks("pdf","rose")
  plot = WindRoseColor (wks,wspd1,wdir1,numPetals,circFr,spdBounds,colorBounds,res)

end
