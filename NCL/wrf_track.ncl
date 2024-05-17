;********************************************************
; Plot storm stracks from wrfout files.
;********************************************************
;
; JUN-18-2005 
; So-Young Ha (MMM/NCAR)
; SEP-01-2006
; Slightly modified by Mary Haley to add some extra comments.
; ===========================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"


begin

; DATES
  ;date = (/1512,1600,1612,1700,1712,1800,1812,1900/)
  ;ndate = dimsizes(date)

 ; sdate = sprinti("%4.0i",date)

; Experiment name (for legend)
  EXP = (/"EXP_I"/)                ; (/"EXP_I","EXP_II","EXP_III"/)
  nexp = dimsizes(EXP)

; To get lat/lon info.

 a = addfile("/lustre/scratch/WUR/ESG/dries002/EMS_WRF_workshop/Emily/answers/wrfout_d01_2005-07-13_12:00:00","r")
  lat2d = a->XLAT(0,:,:)
  lon2d = a->XLONG(0,:,:)
  dimll = dimsizes(lat2d)
  nlat  = dimll(0)
  mlon  = dimll(1)

; Sea Level Pressure
  slp = wrf_user_getvar(a,"slp",0)
  dims = dimsizes(slp)

; Array for track

  timex  = wrf_user_list_times(a)
  tdims = dimsizes(timex)
  ;tdims = tdimsx(0)
  time = new(tdims,string)
  imin = new(tdims,integer)
  jmin = new(tdims,integer)
  smin = new(tdims,integer)

slpmin=new(tdims,float)
; =======
;  ndate
; =======

  do ifs=0,tdims-1
    time(ifs) = timex(ifs)
    print(timex(ifs))
    slp2d = wrf_user_getvar(a,"slp",ifs)

; We need to convert 2-D array to 1-D array to find the minima.
    slp1d     = ndtooned(slp2d)
    smin(ifs) = minind(slp1d)
    slpmin(ifs) = min(slp2d)

; Convert the index for 1-D array back to the indeces for 2-D array.
    minijx = min(ind(slp1d.eq.min(slp2d)))
    minij     = ind_resolve(minijx(0),dims)
    imin(ifs) = minij(0,0)
    jmin(ifs) = minij(0,1)

    print(time(ifs)+" : "+min(slp2d)+" ("+imin(ifs)+","+jmin(ifs)+")")

  end do
;
lines = time + sprintf("%8.2f",slpmin) + sprintf("%8.2f",imin) + sprintf("%8.2f",jmin)
asciiwrite("minslp.txt", lines)

; Graphics section

  wks=gsn_open_wks("pdf","track")              ; Open PS file.
  gsn_define_colormap(wks,"BlGrYeOrReVi200")  ; Change color map.

  res                     = True
  res@gsnDraw             = False             ; Turn off draw.
  res@gsnFrame            = False             ; Turn off frame advance.
  res@gsnMaximize         = True              ; Maximize plot in frame.

  res@tiMainString = ""       ; Main title

  WRF_map_c(a,res,0)                          ; Set up map resources
                                              ;    (plot options)
  plot = gsn_csm_map(wks,res)                 ; Create a map.

; Set up resources for polymarkers.
  gsres                = True
  gsres@gsMarkerIndex  = 16                  ; filled dot
  ;gsres@gsMarkerSizeF = 0.005               ; default - 0.007
  cols                  = (/5,160,40/)

; Set up resources for polylines.
  res_lines                      = True
  res_lines@gsLineThicknessF     = 3.           ; 3x as thick

  dot  = new(tdims,graphic)    ; Make sure each gsn_add_polyxxx call
  line = new(tdims,graphic)    ; is assigned to a unique variable.

; Loop through each date and add polylines to the plot.
  do i = 0,tdims-7,6
     res_lines@gsLineColor           = cols(0)
     xx=(/lon2d(imin(i),jmin(i)),lon2d(imin(i+6),jmin(i+6))/)
     yy=(/lat2d(imin(i),jmin(i)),lat2d(imin(i+6),jmin(i+6))/)
     line(i) = gsn_add_polyline(wks,plot,xx,yy,res_lines)
  end do

  lon1d = ndtooned(lon2d)
  lat1d = ndtooned(lat2d)

; Loop through each date and add polymarkers to the plot.
  do i = 0,tdims-1,6
     print("dot:"+lon1d(smin(i))+","+lat1d(smin(i)))
     gsres@gsMarkerColor  = cols(0)
     dot(i)=gsn_add_polymarker(wks,plot,lon1d(smin(i)),lat1d(smin(i)),gsres)
  end do

; Date (Legend)
  txres               = True
  txres@txFontHeightF = 0.015
  txres@txFontColor   = cols(0)

  txid1 = new(tdims,graphic)
  txid2 = new(tdims,graphic)
; Loop through each date and draw a text string on the plot.
  do i = 0, tdims-1,12
     txres@txJust = "CenterRight"
     ix = smin(i) - 4
     print("Eye:"+ix)
     if(i.eq.1) then
        txres@txJust = "CenterLeft"
        ix = ix + 8
     end if
     txid1(i) = gsn_add_text(wks,plot,tostring(i),lon1d(ix),lat1d(ix),txres)
     txid2(i) = gsn_add_text(wks,plot,sprintf("%6.2f",slpmin(i)),lon1d(ix)+2,lat1d(ix),txres)
  end do
  

; Add marker and text for legend. (Or you can just use "pmLegend" instead.)
  txres@txJust = "CenterLeft"

  txid2 = new(nexp,graphic)
  pmid2 = new(nexp,graphic)
  do i = 0,nexp-1
     gsres@gsMarkerColor  = cols(i)
     txres@txFontColor    = cols(i)
     ii = ((/129,119,109/))  ; ilat
     jj = ((/110,110,110/))  ; jlon
     ji = ii*mlon+jj         ; col x row
    ; pmid2(i) = gsn_add_polymarker(wks,plot,lon1d(ji(i)),lat1d(ji(i)),gsres)
   ;  txid2(i) = gsn_add_text(wks,plot,EXP(i),lon1d(ji(i)+5),lat1d(ji(i)),txres)
  end do

  draw(plot)
  frame(wks)
  
end
