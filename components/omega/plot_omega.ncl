;======================================================================
; mpas_1.ncl
;
; Concepts illustrated:
;   - Plotting MPAS data
;   - Plotting unstructured data
;   - Plotting data with missing lat/lon coordinates
;   - Using cnFillPalette to assign a color palette to contours
;======================================================================
;
; These files are loaded by default in NCL V6.2.0 and newer
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;---Read data from MPAS Grid
    f  = addfile("OmegaSphereMesh.nc","r")
    lonCell = f->lonCell
    latCell = f->latCell
    lonEdge = f->lonEdge
    latEdge = f->latEdge

    f1  = addfile("SWoutput.nc","r")

;---Convert to degrees from radians
    RAD2DEG = 180.d0/(4.0*atan(1.0))
    lonCell = lonCell*RAD2DEG
    latCell = latCell*RAD2DEG
    lonEdge = lonEdge*RAD2DEG
    latEdge = latEdge*RAD2DEG

;---Get time dimension
    fdims = getfiledimsizes(f)
    nt = fdims(1) 
    print(nt)

;****************************
;****************************

;   do n = 0,nt-1
    do n = nt-1,nt-1

    print(n)
;---levels for contours
  nlevels   = 254
; levels    = fspan(min(sp),max(sp),nlevels)
; levels    = fspan(0,1.5,nlevels)
; levels    = fspan(-4,0.2,nlevels)
; levels    = 10^levels
    
;---Start the graphics
    if ( n .lt. 10 ) then
      ifile = "mpas_0" + n 
    else
      ifile = "mpas_" + n 
    end if
    wks = gsn_open_wks("png",ifile)

    ;sp = f->kineticEnergyCell(n,:,0)
    ;sp = f->layerThickness(0,:,0)
    ;sp = f1->LayerThickness(:,0)
    sp = f1->NormalVelocity(:,0)
    levels    = fspan(min(sp),max(sp),nlevels)

    res                      = True              ; Plot mods desired.

    res@gsnMaximize          = True              ; Maximize plot

    res@mpDataBaseVersion     = "MediumRes"
;   res@mpLandFillColor       = "wheat2"
;   res@mpOceanFillColor      = "transparent"   ; no fill
;   res@mpGridAndLimbOn       = False
;   res@mpOutlineOn           = True
;   res@mpFillDrawOrder       = "PostDraw"

    res@cnFillMode           = "RasterFill"      ; turn raster on      
;   res@cnFillMode           = "AreaFill"      ; turn raster on      
    res@cnFillOn             = True              ; color plot desired
    res@cnLinesOn            = True              ; turn off contour lines
    res@cnLineLabelsOn       = False             ; turn off contour labels

;   res@cnFillPalette        = cmap(60:208,:)       ; Don't use white
    res@cnLevelSelectionMode  = "ExplicitLevels"
    res@cnLevels              = levels
;   res@cnFillPalette         = "BlGrYeOrReVi200"
    res@cnFillPalette         = "WhiteBlueGreenYellowRed"
;   res@cnFillPalette         = "hotres"
;   res@cnFillPalette         = "GMT_ocean"
    res@cnLinesOn            = False             ; turn off contour lines

;   res@cnLevelSelectionMode = "ManualLevels"
;   res@cnMinLevelValF       = 0
;   res@cnMaxLevelValF       = 0.5
;   res@cnLevelSpacingF      = 0.005
    res@lbBoxLinesOn         = False

  res@lbLabelStrings        = sprintf("%.2f",levels)  ; format decimal places of labels
  res@lbLabelStride         = 25                ; Over 250 levels, so only plot every 25th level

  res@pmLabelBarWidthF      = 0.75              ; make labelbar slightly longer
  res@pmLabelBarHeightF     = 0.08              ;   and narrower
  res@pmLabelBarOrthogonalPosF = 0.1            ; move labelbar down

    res@gsnAddCyclic         = False

;   res@sfXArray             = lonCell
;   res@sfYArray             = latCell
    res@sfXArray             = lonEdge
    res@sfYArray             = latEdge
    res@gsnAddCyclic         = False

    res@vpXF = 0.0 
    res@vpYF = 0.6 
    res@vpWidthF = 1.0 
    res@vpHeightF = 0.6 

    tt = (n+1) * 5
    if ( tt .lt. 10 ) then
      res@tiMainString         = "Kinetic energy on cell at day 00" + (n+1)*5
    else if ( (tt .ge. 10) .and. (tt .lt. 100)) then
      res@tiMainString         = "Kinetic energy on cell at day 0" + (n+1)*5
    else
      res@tiMainString         = "Kinetic energy on cell at day " + (n+1)*5
    end if
    end if

    plot = gsn_csm_contour_map(wks,sp,res)  

    sfile = ifile+".png"
    afile = "trim_"+ifile+".png"

    system("convert "+" "+sfile+" -trim   "+" "+afile)

    end do

    exit

end
