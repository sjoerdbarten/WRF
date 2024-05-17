load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "WRFOptions.ncl"     ; set basic plot options here
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"


begin


; load WRF model forecasts
a = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")

total_ypoints = 30
total_xpoints = 30

xlat = wrf_user_getvar(a,"XLAT",0)
xlon = wrf_user_getvar(a,"XLONG",0)
ter = wrf_user_getvar(a,"HGT",0)

matrixtoascii=new((/5,total_ypoints*total_xpoints/),float)

teller = 0
do pointsy = 0,total_ypoints-1
   do pointsx = 0,total_xpoints-1
      matrixtoascii(:,teller) = (/pointsy,pointsx,xlat(pointsy,pointsx),xlon(pointsy,pointsx),ter(pointsy,pointsx)/)
      teller = teller + 1
   end do
end do


  ftmx="5f12.4"
  opt=True
  opt@tspace=40
  opt@title="yloc,xloc,latitude,longitude,height"
  opt@fout="petersheight.txt"

  write_matrix(transpose(matrixtoascii),ftmx,opt)

  print("***Successfully finished peter_hgt.ncl script***")


end
