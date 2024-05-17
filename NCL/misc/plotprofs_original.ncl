load "WRFOptions.ncl"     ; set basic plot options here
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
; load radiosondes De Bilt
;  sonde1=asciiread("/home/gert-jan/WRF/WRFV2/run/Results/BSIKDOY2972007/DATA/200710230000.TXT",(/100,7/),"float")
;  sonde2=asciiread("/home/gert-jan/WRF/WRFV2/run/Results/BSIKDOY2972007/DATA/200710231200.TXT",(/100,7/),"float")
;  sonde3=asciiread("/home/gert-jan/WRF/WRFV2/run/Results/BSIKDOY2972007/DATA/200710240000.TXT",(/100,7/),"float")
;  sonde4=asciiread("/home/gert-jan/WRF/WRFV2/run/Results/BSIKDOY2972007/DATA/200710241200.TXT",(/100,7/),"float")
;  sonde5=asciiread("/home/gert-jan/WRF/WRFV2/run/Results/BSIKDOY2972007/DATA/200710250000.TXT",(/100,7/),"float")
;;  sonde6=asciiread("/home/gert-jan/WRF/WRFV2/run/Results/BSIKDOY2972007/DATA/200710251200.TXT",(/100,7/),"float")
;  sonde7=asciiread("/home/gert-jan/WRF/WRFV2/run/Results/BSIKDOY2972007/DATA/200710260000.TXT",(/100,7/),"float")
;  sonde8=asciiread("/home/gert-jan/WRF/WRFV2/run/Results/BSIKDOY2972007/DATA/200710261200.TXT",(/100,7/),"float")
;  sondes=(/sonde1,sonde2,sonde3,sonde4,sonde5,sonde6,sonde7/)
;  print(sondes)
  sondetimes=(/0,12,24,36,48,60,72,84/)
  sondecounter=0
  
  
; load Cabauw profiles
;  data=asciiread("/home/gert-jan/WRF/WRFV2/run/Results/BSIKDOY2972007/DATA/cesar_tower_meteo_200710_lb1_v10_t10.asc",(/24000,53/),"float")
;  zdata=(/200.,140.,80.,40.,20.,10.,2./)
;  zdatawind=(/200.,140.,80.,40.,20.,10./)
;  Tobs=new((/48*6,7/),"float")
;  Qobs=new((/48*6,7/),"float")
;  Uobsvect=new((/48*6,6/),"float")
;  uobs=new((/48*6,6/),"float")
;  vobs=new((/48*6,6/),"float")
;  Udirobs=new((/48*6,6/),"float")

;  Tobs=data(22*24*6:24*24*6-1,3:9)
;  do iz=0,6
;    Tobs(:,iz)=Tobs(:,iz)+9.81/1005.*zdata(iz)
;  end do
;  TdObs=data(22*24*6:24*24*6-1,10:16)
;  Qobs=610.7*10^(7.5*TdObs/(237.3+TdObs))/1013.E2*0.622*1000.
;  Uobsvect=data(22*24*6:24*24*6-1,17:22)
;  Udirobs=data(22*24*6:24*24*6-1,41:46)

;  do itt=0,2*6*24-1
;  do iz=0,5
;      if (Udirobs(itt,iz)>0.and.Udirobs(itt,iz).le.90)
;      uobs(itt,iz)=Uobsvect(itt,iz)*sin(Udirobs(itt,iz)*3.1415927/180.)*-1.
;      vobs(itt,iz)=Uobsvect(itt,iz)*cos(Udirobs(itt,iz)*3.1415927/180.)*-1.
;     end if
;      if (Udirobs(itt,iz).gt.90.and.Udirobs(itt,iz).le.180)
;      uobs(itt,iz)=Uobsvect(itt,iz)*cos((Udirobs(itt,iz)-90.)*3.1415927/180.)*-1.
;      vobs(itt,iz)=Uobsvect(itt,iz)*sin((Udirobs(itt,iz)-90.)*3.1415927/180.)
;      end if
;      if (Udirobs(itt,iz).gt.180.and.Udirobs(itt,iz).le.270)
 ;     uobs(itt,iz)=Uobsvect(itt,iz)*cos((Udirobs(itt,iz)-180.)*3.1415927/180.)
 ;     vobs(itt,iz)=Uobsvect(itt,iz)*sin((Udirobs(itt,iz)-180.)*3.1415927/180.)
;      end if
;      if (Udirobs(itt,iz).gt.270.and.Udirobs(itt,iz).le.360)
;      uobs(itt,iz)=Uobsvect(itt,iz)*cos((Udirobs(itt,iz)-270.)*3.1415927/180.)
;      vobs(itt,iz)=Uobsvect(itt,iz)*sin((Udirobs(itt,iz)-270.)*3.1415927/180.)*-1.
;      end if
;  end do
;  end do

; load WRF model forecasts
  a = addfile("$HOME/WRF/PRACT/WRF/PAR/WRFV3/run/wrfout_d01_2016-12-27_12:00:00_YSU_WSM6.nc","r")
  b = addfile("$HOME/WRF/PRACT/WRF/PAR/WRFV3/run/wrfout_d01_2016-12-27_12:00:00_YSU_WSM6.nc","r")
  c = addfile("$HOME/WRF/PRACT/WRF/PAR/WRFV3/run/wrfout_d01_2016-12-27_12:00:00_YSU_WSM6.nc","r")
  
  times=new(72,float)
  colorplot=new(1,logical)
  
  colorplot= True

  wks_type="ps"
  wks_type@wkOrientation= "landscape"
  wks = gsn_open_wks(wks_type ,"WRFprofs_fogdec2016")
  
  do kk=0,72-1
    times(kk)=kk+12.
  end do
  
  do it = 0, 60-1
  
  p1   = wrf_user_getvar(a, "p",it)     ; pressure is our vertical coordinate
  p2   = wrf_user_getvar(b, "p",it)     ; pressure is our vertical coordinate
  p3   = wrf_user_getvar(c, "p",it)     ; pressure is our vertical coordinate
  z   = wrf_user_getvar(a, "height",it)     ; grid point height
  qv1  = wrf_user_getvar(a,"QVAPOR",it)*1000. ; Qv      
  qv2  = wrf_user_getvar(b,"QVAPOR",it)*1000. ; Qv
  qv3  = wrf_user_getvar(c,"QVAPOR",it)*1000. ; Qv
  
  ql1  = wrf_user_getvar(a,"QCLOUD",it)*1000. ; Qv      
  ql2  = wrf_user_getvar(b,"QCLOUD",it)*1000. ; Qv
  ql3  = wrf_user_getvar(c,"QCLOUD",it)*1000. ; Qv
  
  tc1  = wrf_user_getvar(a,"tc",it)     ; T in C
  tc2  = wrf_user_getvar(b,"tc",it)     ; T in C
  tc3  = wrf_user_getvar(c,"tc",it)     ; T in C
  
  td1  = wrf_user_getvar(a,"td",it)     ; Tdew in C
  td2  = wrf_user_getvar(b,"td",it)     ; Tdew in C
  td3  = wrf_user_getvar(c,"td",it)     ; Tdew in C
  
  tdewdep1=tc1-td1
  tdewdep2=tc2-td2
  tdewdep3=tc3-td3

  u1  = wrf_user_getvar(a,"U",it)     ; 
  u2  = wrf_user_getvar(b,"U",it)     ; 
  u3  = wrf_user_getvar(c,"U",it)     ; 
  u=(/u1,u2,u3/)

  v1  = wrf_user_getvar(a,"V",it)     ; 
  v2  = wrf_user_getvar(b,"V",it)     ; 
  v3  = wrf_user_getvar(c,"V",it)     ; 
  v=(/v1,v2,v3/)
  
  Uvect=(/sqrt(u1(:,150,150)^2+v1(:,150,150)^2),sqrt(u2(:,150,150)^2+v2(:,150,150)^2),sqrt(u3(:,150,150)^2+v3(:,150,150)^2)/)
  Udir=new((/3,33/),"float")
  rd2=45./atan(1.)

  do nmod=0,2
       if (nmod.eq.0)
        vcomp=v1
	ucomp=u1
       end if
       if (nmod.eq.1)
	vcomp=v2
	ucomp=u2
       end if 
       if (nmod.eq.2)
        vcomp=v3
	ucomp=u3
       end if
       do iz=0,32
       Udir(nmod,iz)=atan2(ucomp(iz,150,150),vcomp(iz,150,150))*rd2 +180.
;       if (Udir(nmod,iz).lt.0)
;	  Udir(nmod,iz)=Udir(nmod,iz)+360.
;       end if
      end do
  end do

  theta1=(tc1+273.15)*(p1/p1(0,150,150))^-0.285-273.15
  theta2=(tc2+273.15)*(p2/p2(0,150,150))^-0.285-273.15
  theta3=(tc3+273.15)*(p3/p3(0,150,150))^-0.285-273.15

  theta=(/theta1,theta2,theta3/)
  ql=(/ql1,ql2,ql3/)
  tdewdep=(/tdewdep1,tdewdep2,tdewdep3/)
;  print(theta)
  
  resources                    = True          ; Indicate you want to
                                               ; set some resources.
  plot = new(6,graphic)
  plottemp=new(2,graphic)

  resources@gsnDraw         = False
  resources@gsnFrame        = False
  if (colorplot.eq.True) 
    gsn_define_colormap(wks,"cosam12") 
  end if


;  resources@xyLineColors        = (/20,100,200/)    ; Define line colors.
  resources@xyDashPatterns       = (/0,0,0/)    ; Define line colors.
  resources@xyLineThicknesses   = (/6.,6.,6./) ; Define line thicknesses
;                                               ; (1.0 is the default).
					       
  resources@tiMainString    = "Temperature  "+sprinti("%2.0i",it)+" h. after 200710230000Z"  ; Title for the XY plot
  resources@tiXAxisString   = "Pot temp."    ; Label for the X axis
  resources@tiYAxisString   = "z (m)"    ; Label for the Y axis
  resources@tiMainFont      = "Helvetica" ; Font for title
  resources@tiXAxisFont     = "Helvetica" ; Font for X axis label
  resources@tiYAxisFont     = "Helvetica" ; Font for Y axis label
  resources@tiMainFontHeightF= 0.03
  resources@xyMarkLineModes = (/"Lines"/)
  resources@trYMaxF         = 500
  resources@trXMaxF         = 25
  resources@trXMinF         = 0
  resources@xyLineColors    = (/13,8,5/)
  resources@tmEqualizeXYSizes = True
  resources@tmYLLabelFontHeightF= 0.04
  resources@tiXAxisFontHeightF= 0.03
  resources@tiYAxisFontHeightF= 0.03
  resources@pmLegendDisplayMode="Always"
  resources@pmLegendSide="Top"
  resources@xyExplicitLegendLabels= (/"YSU","MYJ","MRF"/)
   
  plot(0) = gsn_xy(wks,theta(:,:,150,150),z(0:15,150,150),resources)
  resources@xyMarkLineModes = (/"Markers"/)
  resources@xyMarkerSizeF   = 0.03        ; Marker size (default
                                          ; is 0.01)
  resources@xyMarkerThicknessF = 3
  resources@pmLegendDisplayMode="Never"
  
    ; new graph
  resources@tiMainString    = "Spec Hum"  ; Title for the XY plot
  resources@tiXAxisString   = "q (g/kg)"    ; Label for the X axis
  resources@tiYAxisString   = "z (m)"    ; Label for the Y axis
  resources@tiMainFont      = "Helvetica" ; Font for title
  resources@tiXAxisFont     = "Helvetica" ; Font for X axis label
  resources@tiYAxisFont     = "Helvetica" ; Font for Y axis label
  resources@xyMarkLineModes = (/"Lines"/)
  resources@trXMaxF         = 12
  resources@trXMinF         = 0
  qv=(/qv1,qv2,qv3/)
  plot(1) = gsn_xy(wks,qv(:,:,150,150),z(0:15,150,150),resources)
  
  resources@xyMarkLineModes = (/"Markers"/)
  resources@xyMarkerSizeF   = 0.03        ; Marker size (default
                                          ; is 0.01)
  
;    if (it.eq.sondetimes(sondecounter))
;      resources@xyMarkLineModes = (/"Markers"/)
;      resources@xyMarker        = 4
;      resources@xyMarkerSizeF   = 0.03        ; Marker size (default
 ;                                         ; is 0.01)
 ;     spechumsond=610.7*10^(7.5*sondes(sondecounter,:,3)/10/(237.3+sondes(sondecounter,:,1)/10))*0.622/(sondes(sondecounter,:,0)*100)*1000.
 ;     plotsound=gsn_xy(wks,spechumsond,sondes(sondecounter,:,1),resources)
 ;    
 ;    resources@xyMarker        = 0
 ; end if

  resources@tiMainString    = "liquid water content"  ; Title for the XY plot
  resources@tiXAxisString   = "q L (g/kg)"    ; Label for the X axis
  resources@tiYAxisString   = "z (m)"    ; Label for the Y axis
  resources@tiMainFont      = "Helvetica" ; Font for title
  resources@tiXAxisFont     = "Helvetica" ; Font for X axis label
  resources@tiYAxisFont     = "Helvetica" ; Font for Y axis label
  resources@trXMaxF         = 1
  resources@trXMinF         = 0
  resources@xyMarkLineModes = (/"Lines"/)
;  resources@pmLegendDisplayMode="Always"
;  resources@xyExplicitLegendLabels= (/"+ profiler","o sounding","x tower"/)
  plot(2) = gsn_xy(wks,ql(:,:,150,150),z(0:15,150,150),resources)
  resources@xyMarkLineModes = (/"Markers"/)
  resources@xyMarkerSizeF   = 0.03        ; Marker size (default
                                          ; is 0.01)
;  resources@pmLegendDisplayMode="Never"
 ; plot3= gsn_xy(wks,uobs(it*6,:),zdatawind,resources)
  
;    if (it.eq.sondetimes(sondecounter))
;      resources@xyMarkLineModes = (/"Markers"/)
;      resources@xyMarker        = 4
;      resources@xyMarkerSizeF   = 0.03        ; Marker size (default
;                                          ; is 0.01)
      
;      udirsond=new(201,float)
;      vdirsond=new(201,float)
;      do iz=0,40
;      if (sondes(sondecounter,iz,6)>0.and.sondes(sondecounter,iz,6).le.90)
;      udirsond(iz)=sondes(sondecounter,iz,5)*0.1*sin(sondes(sondecounter,iz,6)*3.1415927/180.)*-1.
;      end if
;      if (sondes(sondecounter,iz,6).gt.90.and.sondes(sondecounter,iz,6).le.180)
;      udirsond(iz)=sondes(sondecounter,iz,5)*0.1*cos((sondes(sondecounter,iz,6)-90)*3.1415927/180.)*-1.
;      end if
;      if (sondes(sondecounter,iz,6).gt.180.and.sondes(sondecounter,iz,6).le.270)
;      udirsond(iz)=sondes(sondecounter,iz,5)*0.1*cos((sondes(sondecounter,iz,6)-180)*3.1415927/180.)
;      end if
; ;     if (sondes(sondecounter,iz,6).gt.270.and.sondes(sondecounter,iz,6).le.360)
;      udirsond(iz)=sondes(sondecounter,iz,5)*0.1*cos((sondes(sondecounter,iz,6)-270)*3.1415927/180.)
;      end if
;      end do
;      print(udirsond)
;      plotsound=gsn_xy(wks,udirsond,sondes(sondecounter,0:40,1),resources)

     
 ;    resources@xyMarker        = 0
 ; end if

  resources@tiMainString    = "Dew point depression"  ; Title for the XY plot
  resources@tiXAxisString   = "T-Td (K)"    ; Label for the X axis
  resources@tiYAxisString   = "z (m)"    ; Label for the Y axis
  resources@tiMainFont      = "Helvetica" ; Font for title
  resources@tiXAxisFont     = "Helvetica" ; Font for X axis label
  resources@tiYAxisFont     = "Helvetica" ; Font for Y axis label
  resources@trXMaxF         = 5
  resources@trXMinF         = 0
  resources@xyMarkLineModes = (/"Lines"/)
  plot(3) = gsn_xy(wks,tdewdep(:,:,150,150),z(0:15,150,150),resources)
  resources@xyMarkLineModes = (/"Markers"/)
  resources@xyMarkerSizeF   = 0.03        ; Marker size (default
                                          ; is 0.01)
  ;plot4= gsn_xy(wks,vobs(it*6,:),zdatawind,resources)
  ;;
;
;     if (it.eq.sondetimes(sondecounter))
;      resources@xyMarkLineModes = (/"Markers"/)
;      resources@xyMarker        = 4
;      resources@xyMarkerSizeF   = 0.03        ; Marker size (default
;                                          ; is 0.01)
;					  
;      do iz=0,40
;      if (sondes(sondecounter,iz,6)>0.and.sondes(sondecounter,iz,6).le.90)
;      vdirsond(iz)=sondes(sondecounter,iz,5)*0.1*cos(sondes(sondecounter,iz,6)*3.1415927/180.)*-1.
;      end if
;      if (sondes(sondecounter,iz,6).gt.90.and.sondes(sondecounter,iz,6).le.180)
;      vdirsond(iz)=sondes(sondecounter,iz,5)*0.1*sin((sondes(sondecounter,iz,6)-90)*3.1415927/180.)
;      end if
;      if (sondes(sondecounter,iz,6).gt.180.and.sondes(sondecounter,iz,6).le.270)
;      vdirsond(iz)=sondes(sondecounter,iz,5)*0.1*sin((sondes(sondecounter,iz,6)-180)*3.1415927/180.)
;      end if
; ;     if (sondes(sondecounter,iz,6).gt.270.and.sondes(sondecounter,iz,6).le.360)
;      vdirsond(iz)=sondes(sondecounter,iz,5)*0.1*sin((sondes(sondecounter,iz,6)-270)*3.1415927/180.)*-1.
;      end if
;      end do
 ;     plotsound=gsn_xy(wks,vdirsond,sondes(sondecounter,0:40,1),resources)
 ;    
  ;   resources@xyMarker        = 0
;  end if


  
  resources@tiMainString    = "Wind magnitude"  ; Title for the XY plot
  resources@tiXAxisString   = "U_vect (m/s)"    ; Label for the X axis
  resources@tiYAxisString   = "z (m)"    ; Label for the Y axis
  resources@tiMainFont      = "Helvetica" ; Font for title
  resources@tiXAxisFont     = "Helvetica" ; Font for X axis label
  resources@tiYAxisFont     = "Helvetica" ; Font for Y axis label
  resources@trXMaxF         = 18
  resources@trXMinF         = 0
  resources@xyMarkLineModes = (/"Lines"/)
  plot(4) = gsn_xy(wks,Uvect,z(0:15,150,150),resources)
  
  resources@xyMarkLineModes = (/"Markers"/)
  resources@xyMarkerSizeF   = 0.03        ; Marker size (default
                                          ; is 0.01)
 ; plot3= gsn_xy(wks,Uobsvect(it*6,:),zdatawind,resources)
  
  

  resources@tiMainString    = "Direction"  ; Title for the XY plot
  resources@tiXAxisString   = "dir (o)"    ; Label for the X axis
  resources@tiYAxisString   = "z (m)"    ; Label for the Y axis
  resources@tiMainFont      = "Helvetica" ; Font for title
  resources@tiXAxisFont     = "Helvetica" ; Font for X axis label
  resources@tiYAxisFont     = "Helvetica" ; Font for Y axis label
  resources@trXMaxF         = 360
  resources@trXMinF         = 0
  resources@xyMarker        = 0
  resources@xyMarkLineModes = (/"Lines"/)
  plot(5) = gsn_xy(wks,Udir,z(0:15,150,150),resources)

  resources@xyMarkLineModes = (/"Markers"/)
  resources@xyMarkerSizeF   = 0.03        ; Marker size (default
                                          ; is 0.01)
  ;plot4= gsn_xy(wks,Udirobs(it*6,:),zdata,resources)
  

;  if (it.lt.12)
;    resources@xyMarker        = 2
;    plot5= gsn_xy(wks,profiler(it*33:it*33+32,3),profiler(it*33:it*33+32,1),resources)
;    overlay(plot(5),plot5)
;    resources@xyMarker        = 0
;  end if
;  if (it.ge.12)
;    resources@xyMarker        = 2
;    plot5= gsn_xy(wks,profilerdir(2*(it-12),3:26),profileralti(3,3:26),resources)
;    overlay(plot(5),plot5)
;    resources@xyMarker        = 0
;  end if 

;  if (it.eq.sondetimes(sondecounter))
;      resources@xyMarkLineModes = (/"Markers"/)
;      resources@xyMarker        = 4
;      resources@xyMarkerSizeF   = 0.03        ; Marker size (default
;                                          ; is 0.01)
;     plotsound=gsn_xy(wks,sondes(sondecounter,:,6),sondes(sondecounter,:,1),resources)
;     
;     sondecounter=sondecounter+1
;     resources@xyMarker        = 0
;  end if
  
  resout = True
  gsn_panel(wks,plot,(/2,3/),resout)
  txres=True
  txres@txFontHeightF=0.02
  gsn_text_ndc(wks,"+ profiler; o sounding; x tower",0.3,0.1,txres)

;  legend=create "Legend" legendClass wks
;  "lgLineLabelStrings": (/"YSU","MYJ","MRF"/)
;  end create
;  draw(legend)

  delete(resources)

  end do
  
end
