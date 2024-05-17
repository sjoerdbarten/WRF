load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "WRFOptions.ncl"     ; set basic plot options here
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

; load WRF model forecasts
 a = addfile("./wrfout_d01_2018-06-05_00:00:00.nc","r")


;===========END OF USER MODIFICATION=========================================================================
  
  dummy=a->XLAT(0,:,:)
  dummydims=dimsizes(dummy)
  xix=new(dummydims(0)*dummydims(1),float)
  xiy=new(dummydims(0)*dummydims(1),float)
  xlatitude=new(dummydims(0)*dummydims(1),float)
  xlongitude=new(dummydims(0)*dummydims(1),float)
  
  counter=0
  do ix =0,dummydims(0)-1
    do iy=0,dummydims(1)-1
  
    latitude=a->XLAT(0,iy,ix)
    longitude=a->XLONG(0,iy,ix)
    xix(counter)=ix*1.
    xiy(counter)=iy*1.
    xlatitude(counter)=latitude
    xlongitude(counter)=longitude
    counter=counter+1
    end do
  end do
 
  matrixtoascii=(/xix(0:counter-1),xiy(0:counter-1),xlatitude(0:counter-1),xlongitude(0:counter-1)/)
 
  ftmx="4f12.4"
  opt=True
  opt@tspace=40
 
  opt@title="ix, iy, longitude, latitude"
  opt@fout="coorc.txt"
  write_matrix(transpose(matrixtoascii),ftmx,opt)
  print("Successfully finished script make_coorcs.ncl")
   
 end
