;+
; NAME
;   ml_griddataWF
;
; Quick hack for WFIRST to allow turning off SDSS flags.
;
; PURPOSE:
;   Grid irregularly-sampled data onto a regular grid.  Assumes
;   input is in the form of 3 1-dimensional vectors containing
;   the x coordinates, y coordinates, and fluxes respectively.
;   Can also optionally take in a vector of inverse variance values
;   with which to construct an inverse variance arry for the final
;   output grid.
;
;   Note that these individual-fiber inverse variances are NOT
;   used in determining the weights for the actual combination step.
;   This is because they overwhelm the radial distance weighting
;   (since ivar ~ 1/flux), effectively downranking central pixels 
;   compared to wings of profile.  This causes >10% errors in total 
;   flux conservation and biases the effective PSF substantially.
;   The effect is worse the smaller the intrinsic size of the object.
;
;   Inverse variances are thus only used for their zeroes- anything
;   with zero inverse variance is ignored, while anything else is
;   combined assuming constant inverse variance.
;
;   Note that the output grid does NOT account for covariance.
;
;   Algorithm used is a flux-conserving version of Shepards
;   method with exponential inverse distance weighting within a 
;   limiting radius.  The algorithm is extremely similar to that
;   used by the CALIFA survey described by Sanchez+, with the addition
;   of inverse variances in the weight function.
;
;   Note that this routine is intended to be run on a single
;   wavelength at a time; it does not treat wavelength dependance.
;
; CALLING SEQUENCE:
;   image=ml_griddata(x,y,f,dim_out,rlim,sigma,[scale=,ivar=,invarimg=,maskvec=,maskimg=])
;
; INPUTS:
;   x: 1-d vector containing X coordinates of fiber centers
;      in output pixel coordinates
;   y: 1-d vector containing Y coordinates of fiber centers
;      in output pixel coordinates
;   f: 1-d vector containing flux in each fiber at the given location
;   dim_out: Array giving the [x,y] dimensions of output image
;   rlim: Scalar fixing boundary limit radius for influence of a given
;         fiber, in output pixel units.
;   sigma: Scalar defining gaussian sigma for the spatial weight
;         function, in output pixel units.
;
; OPTIONAL INPUTS:
;   scale: Scaling factor to be applied to the output.  Required
;          as input flux units tend to be per fiber area, while
;          output flux units should be per spaxel.
;          1.0 by default
;   ivar: 1-d vector containing inverse variances corresponding to f
;   maskvec: 1-d vector containing MANGA_DRPPIXFLAG mask values for f
;   \noflag: Do not do any flagging
;
; OUTPUTS:
;   image: 2d image array containing data interpolated onto regular grid.
;
; OPTIONAL OUTPUTS:
;   invarimg: 2d image array containing inverse variance of output image.
;   maskimg: 2d image array containing output MANGA_DRPPIXFLAG values
;
; PROCEDURES CALLED:
;
; EXAMPLES:
;   image=ml_griddata(x,y,f,[40,40],1.6/pixscale,0.7/pixscale,scale=0.43,ivar=ivar, $
;     invarimg=invarimg,maskvec=maskvec,maskimg=maskimg)
;
; BUGS: 
;
; REVISION HISTORY:
;   06-Dec-2013  David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;       First written based on older version without ivar or scaling,
;       incorporating speed tweaks by Maryna Tsybulska.
;   07-Apr-2014  Substantially retooled inverse variance calculations.
;       Now ONLY uses this information to do rejections, does not combine
;       spectra inverse variance weights with spatial kernel weights, because
;       this overwhelms the spatial weighting. (D. Law)
;   27-Apr-2014  Overhauled to include mask calculations (D. Law)
;   23-Apr-2015  Revise to reduce calls to bitmask routines, providing
;       major speedup (R. D'Souza, D. Law)
;-

function ml_griddatawf, x, y, f, dim_out, rlim, sigma, scale=scale, ivar=ivar, invarimg=invarimg, maskvec=maskvec, maskimg=maskimg,noflag=noflag

; Dimensions
ntot = n_elements(f) ; Number of total samples
all_dim = [dim_out, ntot]; Output X size, Y size, total samples

; Default scale multiplier is 1.0
if (keyword_set(scale)) then scale=scale $
else scale=1.0

; Default inverse variance vector has everything
; set to a constant of 1
if (keyword_set(ivar)) then ivar=ivar $
else ivar=replicate(1.,n_elements(f))

; Default mask vector has everything equal to 0
if (keyword_set(maskvec)) then maskvec=maskvec $
else maskvec=replicate(0L,n_elements(f))

; Output images:
fimg = dblarr(dim_out)
invarimg=dblarr(dim_out)
maskimg=lonarr(dim_out)

; Safety check that x,y,f have the same length
nx=n_elements(x)
ny=n_elements(y)
nf=n_elements(f)
ni=n_elements(ivar)
nm=n_elements(maskvec)
if ((nx ne ny) or (nx ne nf) or (nx ne ni) or (nx ne nm)) then begin
  splog,'WARNING!  x,y,f,ivar,maskvec do not have the same length.'
  return,fimg
endif

; X and Y output pixel coordinate arrays
arr_xcoord = dindgen(dim_out[0])
arr_ycoord = dindgen(dim_out[1])

; Calculate a 3d array of weights for all locations in the image
; for all input elements.  Weight array combines a term describing
; distance from the input element and a term describing the inverse
; variance of the input element.
arr_weights = dblarr(all_dim)

for i=0,ntot-1 do begin
   ; Array defining radii in output grid away from the fiber location
   ; Initialize to rlim+1 (i.e. something greater than rlim
   ; radius criterion within which we actually care about the radii)
   arr_radius = replicate(double(rlim+1), dim_out)
   ; Figure out the region of influence in which we actually need to calculate radii
   xmin = x[i] - rlim > 0
   xmax = x[i] + rlim < (dim_out[0]-1)
   ymin = y[i] - rlim > 0
   ymax = y[i] + rlim < (dim_out[1]-1)
  
   ; Calculate actual radii in this region of influence
   dim_c = [n_elements(arr_xcoord[xmin:xmax]),n_elements(arr_ycoord[ymin:ymax])]
   arr_radius[xmin:xmax,ymin:ymax] = sqrt(rebin((arr_xcoord[xmin:xmax]-x[i])^2,dim_c) $
                              + rebin(transpose((arr_ycoord[ymin:ymax]-y[i])^2),dim_c))
   ; Calculate weights in the region where radius < rlim
   tocalc = where(arr_radius le rlim)
   ; Weights are the exponential falloff of influence with
   ; increasing distance from the input fiber location.  Things with ivar=0
   ; will be given zero weight later on, keep them non-zero here
   ; so that we can just track which input elements affect which
   ; output pixels
   arr_weights[tocalc+i*n_elements(arr_radius)] = exp(-0.5/sigma^2*arr_radius[tocalc]^2)
  
endfor
; Figure out the normalization matrix- sum of arr_weights
; Safety case for where there is only 1 exposure, so no 3rd dimension
if (ntot eq 1) then matr_norm=arr_weights $
; Sum over the 3rd dimension of arr_weights.  First make sure that
; any input element that has zero inverse variance contributions nothing
; to the normalization sum.  Do this by taking a logical AND between ivar
; and 1 (which will give 0 where ivar=0, and 1 elsewhere) and
; recasting this as a 3d array of the correct dimensions so that it
; can simply be multiplied by the arr_weights.
else matr_norm = total(rebin(reform((ivar and 1),1,1,ntot),dim_out[0],dim_out[1],ntot)*arr_weights,3)

; Flag where the normalization matrix is zero; there is no good data here
nodata=where(matr_norm eq 0,nnodata)
; We don't want to divide by zero where there is no data; set the normalization
; matrix to 1 in these cases
if (nnodata gt 0) then matr_norm[nodata]=1.

; Set up pixel flags if so desired
if (~keyword_set(noflag)) then begin
  ; read and store all the MANGA flags
  flagdeadfiber=sdss_flagval('MANGA_DRP3PIXMASK','DEADFIBER')
  flaglocov=sdss_flagval('MANGA_DRP3PIXMASK','LOWCOV')
  flagnocov=sdss_flagval('MANGA_DRP3PIXMASK','NOCOV')
  flagnouse=sdss_flagval('MANGA_DRP3PIXMASK','DONOTUSE')
  flag3dreject = sdss_flagval('MANGA_DRP2PIXMASK','3DREJECT')

  ; What spatial elements have entirely zero arr_weights?  These are the
  ; output pixels with no coverage.  (Note that this isn't the same as where
  ; matr_norm is zero, because matr_norm can be zero even in the middle
  ; of the ifu bundle if fibers had zero inverse variance.)

  temp=fix(not(total(arr_weights,3)))
  mask_nocov=temp*flagnocov
  ; Everything with no coverage gets added to the 'DONOTUSE' flag too
  mask_dnu=temp*flagnouse
endif

ngood=lonarr(dim_out)
ndead=lonarr(dim_out)
nbad=lonarr(dim_out)

; Apply the weights to calculate the output flux and inverse variance
; arrays by summing over all input elements.  Make sure to multiply
; by (ivar[i] and 1) to zero out contributions from input elements
; with zero inverse variance (i.e., bad values).
for i=0,ntot-1 do begin
   alpha=arr_weights[*,*,i]*(ivar[i] and 1) / matr_norm
   fimg += f[i] * alpha
   ; Combine inverse variances by using formula
   ; sigma^2=Sum(alpha^2*sigma_i^2)
   ; where sigma = 1./sqrt(ivar)
   ; Note that ivarimg here defined is the *variance* array
   ; not the *inverse variance* array- we'll flip it later.
   if (ivar[i] ne 0) then invarimg += alpha*alpha/ivar[i]

   if (~keyword_set(noflag)) then begin
     ; Work out pixel counts for the masks
     ; Count pixels with good values
     if (ivar[i] ne 0) then ngood+=fix(arr_weights[*,*,i] and 1)
     ; Count pixels with '3DREJECT' set
     if ((maskvec[i] and flag3dreject) ne 0) then nbad+=fix(arr_weights[*,*,i] and 1)
     ; Count pixels with 'DEADFIBER' set
     if ((maskvec[i] and flagdeadfiber) ne 0) then ndead+=fix(arr_weights[*,*,i] and 1)
   endif
endfor

if (~keyword_set(noflag)) then begin
  ; Set stuff with ngood < 30% *max(ngood) to LOWCOV
  mask_lowcov=lonarr(dim_out)
  index=where(ngood le 0.3*max(ngood),nindex)
  if (nindex gt 0) then mask_lowcov[index]=flaglocov
  ; Everything with LOWCOV set also gets DONOTUSE set
  if (nindex gt 0) then mask_dnu[index]=flagnouse
  ; Set DEADFIBER masks
  mask_dead=(ndead and 1)*flagdeadfiber

  ; Combine all bitmasks
  maskimg = mask_nocov
  maskimg = maskimg or mask_lowcov
  maskimg = maskimg or mask_dead
  maskimg = maskimg or mask_dnu
endif

; Flip from variance to real inverse variance
invarimg=1./invarimg
; Inverse variance is zero where no data
if (nnodata gt 0) then invarimg[nodata]=0.

; Account for an areal correction factor (if applicable).
; This would be used, for instance, to transform from input
; fluxes per fiber (i.e., in an area 3.14*r_fib^2) to output
; fluxes per spatial element (spaxel)
fimg *= scale
invarimg /= scale*scale

return, fimg
end
