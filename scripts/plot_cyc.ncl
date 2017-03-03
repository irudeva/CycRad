begin
;************************************************
; read in netCDF file
;************************************************
yr = 2005
t = 1282
  slpyr = yr-1
  a = addfile("/Users/Irina/work/DATA/ERAint/erain.mslp."+slpyr+".nc","r")
  slp = short2flt(a->msl(t,:,:))
  slp1 =0.01*slp
  copy_VarMeta(slp, slp1)



  b = addfile("../output/cyc."+yr+".nc","r")
  var = b->cyclone(t,:,:)

  wks = gsn_open_wks("png","../output/cyc.t"+t)                  ; send graphics to PNG file

res = True
res@gsnDraw              = False             ; do not draw the plot
res@gsnFrame             = False             ; do not advance the frame

res@cnFillOn            = True        ; Turn on contour fill
res@cnLinesOn           = False       ; Turn off contour lines
res@cnFillPalette       = "BlueDarkRed18"  ; Set color map

res@tiMainString        = "t = " + t

plot1 = gsn_csm_contour_map(wks,slp1,res)

res@cnFillOn            = False        ; Turn on contour fill
res@cnLinesOn           = True       ; Turn off contour lines
res@cnMonoFillColor     = True         ; Use single pattern color
res@cnLineColor         = "green"      ; Turn off contour lines

plot2 = gsn_csm_contour(wks,var,res)         ; create a default plot


overlay(plot1,plot2)



draw(plot1)
frame(wks)


end
