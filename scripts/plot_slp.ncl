begin
;************************************************
; read in netCDF file
;************************************************
  a = addfile("/Users/Irina/work/DATA/ERAint/erain.mslp.2005.nc","r")
;************************************************
; read in zonal winds
;************************************************
  var = short2flt(a->msl(1287,:,:))
  printVarSummary(var)
  var1 =0.01*var
  copy_VarMeta(var, var1)
  printVarSummary(var1)
;************************************************
; create default plot
;************************************************
  wks = gsn_open_wks("png","../output/slp.t1287")                  ; send graphics to PNG file
  plot = gsn_csm_contour_map(wks,var1,False)         ; create a default plot
end
