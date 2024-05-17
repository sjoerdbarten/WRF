
;   Calculate lat, lon naar yloc and xloc

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

 a = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")



; Note, coordinates should be in decimales, not in minutes seconds etc!!
  

  loc = wrf_user_latlon_to_ij(a, 51.971, 4.93)
  print("Yloc/Xloc location is: " +loc)
  

 
end
