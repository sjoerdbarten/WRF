load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "WRFOptions.ncl"     ; set basic plot options here
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

; load WRF model forecasts
 a = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")

  nmodeltimesteps=66
  nmodellevels=45
 
  yloc=15
  xloc=10
  
  lat = a->XLAT(0,yloc,xloc)
  lon = a->XLONG(0,yloc,xloc)
  print(lat)
  print(lon)

;===========END OF USER MODIFICATION=========================================================================
  
  times=new(nmodeltimesteps+1,float)
  klev=new(nmodellevels,float)
  
  do kk=0,nmodeltimesteps
    times(kk)=kk
  end do

  do kk=1,nmodellevels
    klev(kk-1)=kk
  end do
  
  matrixtoascii=new((/9,(nmodellevels+1)*nmodeltimesteps+nmodeltimesteps/),float)
    
  ; Calculate PBL depth
  do k=0,nmodeltimesteps-1
  Ux = a->U(k,:,yloc,xloc)
  Vx = a->V(k,:,yloc,xloc)
  Qx = a->QVAPOR(k,:,yloc,xloc)
  Qcx = a->QCLOUD(k,:,yloc,xloc)
  zx = wrf_user_getvar(a, "z",k)     ; grid point height
  pressx = wrf_user_getvar(a, "pressure",k)
  Tx = wrf_user_getvar(a, "tc",k)
  Wx = wrf_user_getvar(a,"wa",k)
  
  latitude=a->XLAT(:,yloc,xloc)
  longitude=a->XLONG(:,yloc,xloc)
  ;Make header ervy hour
  matrixtoascii(0,k*(nmodellevels+1)+k)=k
  matrixtoascii(1:8,k*(nmodellevels+1)+k)=-999
  
  matrixtoascii(0,k*(nmodellevels+1)+k+1:k*(nmodellevels+1)+nmodellevels+k)=klev*1.
  matrixtoascii(1,k*(nmodellevels+1)+k+1:k*(nmodellevels+1)+nmodellevels+k)=zx(:,yloc,xloc)
  matrixtoascii(2,k*(nmodellevels+1)+k+1:k*(nmodellevels+1)+nmodellevels+k)=pressx(:,yloc,xloc)
  matrixtoascii(3,k*(nmodellevels+1)+k+1:k*(nmodellevels+1)+nmodellevels+k)=Ux
  matrixtoascii(4,k*(nmodellevels+1)+k+1:k*(nmodellevels+1)+nmodellevels+k)=Vx
  matrixtoascii(5,k*(nmodellevels+1)+k+1:k*(nmodellevels+1)+nmodellevels+k)=Wx(:,yloc,xloc)
  matrixtoascii(6,k*(nmodellevels+1)+k+1:k*(nmodellevels+1)+nmodellevels+k)=Tx(:,yloc,xloc)
  matrixtoascii(7,k*(nmodellevels+1)+k+1:k*(nmodellevels+1)+nmodellevels+k)=Qx
  matrixtoascii(8,k*(nmodellevels+1)+k+1:k*(nmodellevels+1)+nmodellevels+k)=Qcx
  end do
  
  ftmx="9f14.6"
  opt=True
  opt@tspace=40
 
  opt@title="IZ, Z(m),P(hPa),U(m/s),V(m/s),W(m/s),T(degC),Q(kg/kg),Qcloud(kg/kg)"
  opt@fout="3D.txt"
  write_matrix(transpose(matrixtoascii),ftmx,opt)
  print("***Successfully finished make3D.ncl script***")
   
 end
