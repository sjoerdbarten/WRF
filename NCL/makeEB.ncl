load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "WRFOptions.ncl"     ; set basic plot options here
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

; Calculate PBLdepth

  function pbldepth(u,v,z,T)
  begin
  ricrit=new(1,float)
  pblheight=new(1,float)
  ribulk=new(41,float)
  ricrit=0.2
   do i=0,28
    ribulk(i)=9.81/300.*z(i)*(T(i)-T(0))/(u(i)^2+v(i)^2)
   end do
   i=0
   do while (ribulk(i).lt.ricrit)
    i=i+1
   end do
   pblheight=z(i-1)+(z(i)-z(i-1))/(ribulk(i)-ribulk(i-1))*(ricrit-ribulk(i-1))
  return(pblheight)
  end


begin

; load WRF model forecasts
  a = addfile("./wrfout_d02_2018-06-05_00:00:00.nc","r")

  nmodeltimesteps=66
  
  yloc=15
  xloc=15

;===========END OF USER MODIFICATION=========================================================================
  
  times=new(nmodeltimesteps+1,float)
 
  do kk=0,nmodeltimesteps
    times(kk)=kk
  end do
  
;  plot = new(16,graphic)
;  colorplot=new(1,logical)
;  colorplot = True
;  resources=True
;  resources@gsnDraw         = False
;  resources@gsnFrame        = False
;  if (colorplot.eq.True) 
;    gsn_define_colormap(wks,"cosam12") 
;  end if

   
  ; Calculate PBL depth

  pbl1 = a->PBLH(:,yloc,xloc)
  pbl2=new(nmodeltimesteps+1,float)
  do k=0,nmodeltimesteps
  U1 = a->U(k,:,yloc,xloc)
  V1 = a->V(k,:,yloc,xloc)
  geop1 = a->PH(k,:,yloc,xloc)
  geo1 = a->PHB(k,:,yloc,xloc)
  pp1 = a->P(k,:,yloc,xloc)
  ppb1 = a->PB(k,:,yloc,xloc)
  p1=pp1+ppb1
  z1=(geo1+geop1)/9.81
  TEMP1 = wrf_user_getvar(a,"tc",k)    ; 3D tc 
  T1=TEMP1(:,yloc,xloc)
  T1P=(T1+273.15)*(100000/p1)^0.286
  pbl2(k)=pbldepth(U1,V1,z1,T1P)
  end do
  
  hfx = a->HFX(:,yloc,xloc)
  tsk = a->TSK(:,yloc,xloc)-273.15
  GRDFLX = a->GRDFLX(:,yloc,xloc)
  lfx = a->LH(:,yloc,xloc)
  ust = a->UST(:,yloc,xloc)
  T2= a->T2(:,yloc,xloc) -273.15 
  U10 = a->U10(:,yloc,xloc)
  V10 = a->V10(:,yloc,xloc)
  swdown = a->SWDOWN(:,yloc,xloc)
  alb = a->ALBEDO(:,yloc,xloc)
  glw = a->GLW(:,yloc,xloc)
  gsw = swdown*(1.-alb)
  luindex= a->LU_INDEX(:,yloc,xloc)
  Q2m=a->Q2(:,yloc,xloc)
  surfpress=a->PSFC(:,yloc,xloc)
  soilmoist=a->SMOIS(:,0,yloc,xloc)
  soiltemp=a->TSLB(:,0,yloc,xloc)
  albedo=a->ALBEDO(:,yloc,xloc)
  emissiv=a->EMISS(:,yloc,xloc)
  latitude=a->XLAT(:,yloc,xloc)
  longitude=a->XLONG(:,yloc,xloc)
     
  matrixtoascii=(/times(0:nmodeltimesteps),surfpress(0:nmodeltimesteps),T2(0:nmodeltimesteps),Q2m(0:nmodeltimesteps),U10(0:nmodeltimesteps),V10(0:nmodeltimesteps),ust(0:nmodeltimesteps),swdown(0:nmodeltimesteps),gsw(0:nmodeltimesteps),glw(0:nmodeltimesteps),tsk(0:nmodeltimesteps),hfx(0:nmodeltimesteps),lfx(0:nmodeltimesteps),GRDFLX(0:nmodeltimesteps),pbl1(0:nmodeltimesteps),pbl2(0:nmodeltimesteps),luindex(0:nmodeltimesteps),soiltemp(0:nmodeltimesteps),soilmoist(0:nmodeltimesteps),albedo(0:nmodeltimesteps),emissiv(0:nmodeltimesteps),longitude(0:nmodeltimesteps),latitude(0:nmodeltimesteps)/)
 
  ftmx="23f12.4"
  opt=True
  opt@tspace=40
  opt@title="Time,Surfpress(Pa),T2m(C),Q2m(kg/kg),U10(m/s),V10(m/s),Ustar(m/s),SWdown(W/m2),Net_SW(W/m2),GLW(W/m2),H(W/m2),LvE(W/m2),G(W/m2),luindex(-),soiltemp(K),soil_moisture(-),albedo(-),emissivity(-),longitude,latitude"
  opt@fout="EB.txt"
  write_matrix(transpose(matrixtoascii),ftmx,opt)
  print("***Successfully finished script makeEB.ncl***")
   
 end
