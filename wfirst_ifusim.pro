; v5 (July 23 2015)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Resample psf cube to match to our wavelength sampling
function resample_psf_cube,psfcube,psfwave,wave

nx=(size(psfcube))[1]
ny=(size(psfcube))[2]
nz=(size(psfcube))[3]
nwave=n_elements(wave)

xgrid=indgen(nx)
ygrid=indgen(ny)
zgrid=indgen(nz)

; Compute the locations of 'wave' in channel units
zvals=interpol(zgrid,psfwave,wave)

; Linear interpolation from the input cube
newcube=interpolate(psfcube,xgrid,ygrid,zvals,/grid)

return,newcube
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Given a footprint mask, make a 'mesh' version
; for overlay purposes
function makemesh,mask,nx,ny,filename,color=rcol

if (~keyword_set(rcol)) then rcol='green'

openw,lun,filename,/get_lun

preamble0='# Region file'
preamble1='global color='+rcol+' dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
preamble2='image'

printf,lun,preamble0
printf,lun,preamble1
printf,lun,preamble2

Xsize=(size(mask))[1]
Ysize=(size(mask))[2]

mesh=mask
mesh[*,*]=0.

imin=0
imax=max(mask)
n=imax-imin+1

for i=0,n-1 do begin
  index=where(mask eq i,nindex)
  ; Split index into row, column indices
  s = SIZE(mask)
  ncol = s(1)
  col = index MOD ncol
  row = index / ncol
  ; Loop to see if they are edges
  for j=0,nindex-1 do begin
    thesum=total(mask[col[j]-1:col[j]+1,row[j]-1:row[j]+1])
    if (thesum ne 9*i) then mesh[col[j],row[j]]=1
  endfor
endfor

; Loop over rows
for i=0,ny-1 do begin
  ; Starting mask id
  idstart=i*nx
  ; Ending mask id
  idstop=(i+1)*nx-1
  ; Starting box
  index=where(mask eq idstart)
  s = SIZE(mask)
  ncol = s(1)
  col = index MOD ncol
  row = index / ncol
  ; Identify lower left corner by distance from actual corner
  dist=sqrt((0.-col)^2+(0.-row)^2)
  junk=min(dist,minat)
  llx=col[minat]+0.5
  lly=row[minat]+0.5
  ; Ending box
  index=where(mask eq idstop)
  s = SIZE(mask)
  ncol = s(1)
  col = index MOD ncol
  row = index / ncol
  ; Identify lower right corner by distance from actual corner
  dist=sqrt((Xsize-col)^2+(0.-row)^2)
  junk=min(dist,minat)
  lrx=col[minat]+1.5
  lry=row[minat]+0.5
  ; Write the line
  printf,lun,strcompress('line('+string(llx)+','+string(lly)+','+string(lrx)+','+string(lry)+')',/remove_all)
endfor

; Do the top row
i=ny-1
  ; Starting mask id
  idstart=i*nx
  ; Ending mask id
  idstop=(i+1)*nx-1
  ; Starting box
  index=where(mask eq idstart)
  s = SIZE(mask)
  ncol = s(1)
  col = index MOD ncol
  row = index / ncol
  ; Identify top left corner by distance from actual corner
  dist=sqrt((0.-col)^2+(Ysize-row)^2)
  junk=min(dist,minat)
  tlx=col[minat]+0.5
  tly=row[minat]+1.5
  ; Ending box
  index=where(mask eq idstop)
  s = SIZE(mask)
  ncol = s(1)
  col = index MOD ncol
  row = index / ncol
  ; Identify lower right corner by distance from actual corner
  dist=sqrt((Xsize-col)^2+(Ysize-row)^2)
  junk=min(dist,minat)
  trx=col[minat]+1.5
  try=row[minat]+1.5
  ; Write the line
  printf,lun,strcompress('line('+string(tlx)+','+string(tly)+','+string(trx)+','+string(try)+')',/remove_all)

; Loop over columns
for i=0,nx-1 do begin
  ; Starting mask id
  idstart=i
  ; Ending mask id
  idstop=(ny-1)*nx+i
  ; Starting box
  index=where(mask eq idstart)
  s = SIZE(mask)
  ncol = s(1)
  col = index MOD ncol
  row = index / ncol
  ; Identify lower left corner by distance from actual corner
  dist=sqrt((0.-col)^2+(0.-row)^2)
  junk=min(dist,minat)
  llx=col[minat]+0.5
  lly=row[minat]+0.5
  ; Ending box
  index=where(mask eq idstop)
  s = SIZE(mask)
  ncol = s(1)
  col = index MOD ncol
  row = index / ncol
  ; Identify top left corner by distance from actual corner
  dist=sqrt((0.-col)^2+(Ysize-row)^2)
  junk=min(dist,minat)
  tlx=col[minat]+0.5
  tly=row[minat]+1.5
  ; Write the line
  printf,lun,strcompress('line('+string(llx)+','+string(lly)+','+string(tlx)+','+string(tly)+')',/remove_all)
endfor

; Do the final column
i=nx-1
  ; Starting mask id
  idstart=i
  ; Ending mask id
  idstop=(ny-1)*nx+i
  ; Starting box
  index=where(mask eq idstart)
  s = SIZE(mask)
  ncol = s(1)
  col = index MOD ncol
  row = index / ncol
  ; Identify lower right corner by distance from actual corner
  dist=sqrt((Xsize-col)^2+(0.-row)^2)
  junk=min(dist,minat)
  lrx=col[minat]+1.5
  lry=row[minat]+0.5
  ; Ending box
  index=where(mask eq idstop)
  s = SIZE(mask)
  ncol = s(1)
  col = index MOD ncol
  row = index / ncol
  ; Identify top right corner by distance from actual corner
  dist=sqrt((Xsize-col)^2+(Ysize-row)^2)
  junk=min(dist,minat)
  trx=col[minat]+1.5
  try=row[minat]+1.5
  ; Write the line
  printf,lun,strcompress('line('+string(lrx)+','+string(lry)+','+string(trx)+','+string(try)+')',/remove_all)

close,lun
free_lun,lun

; DRL- note some problems with box size???
; Last box in each row is too big?  Prob b/c not integer number
; of pixels per slice, being icky at the edge...

return,mesh
end




pro wfirst_ifusim,run=run

if (keyword_set(run)) then outdir=strcompress('out/'+run+'/',/remove_all) $
else outdir='out/'

spawn,strcompress('mkdir -p '+outdir)

; Read in parameters
exposures=yanny_readone(concat_dir(getenv('WFIRST_DIR'),'ifupars.par'),'IFUPARS',hdr=hdr)
nexp=n_elements(exposures)

spawn,strcompress('cp '+concat_dir(getenv('WFIRST_DIR'),'ifupars.par') + ' ' + outdir)
spawn,strcompress('cp '+concat_dir(getenv('WFIRST_DIR'),'wfirst_ifusim.pro') + ' ' + outdir)

; Constants
radpas=4.8481e-6; Radians per arcsecond
c=2.9979e10; Speed of light (cm/s)
hc=1.986e-16; Planck constant times speed of light (erg cm)
k=1.381e-16; Boltzmann constant (erg/K)
teldiam=float(yanny_par(hdr,'teldiam'))
A=!PI*teldiam*teldiam/4.; Area of primary (cm^2)

; Type of IFU (currently only 'slicer')
ifutype=strtrim(yanny_par(hdr,'ifutype'),2)

; For convenience define some colors for the regions files for
; each exposure
rcol=['green','red','cyan','magenta','yellow','blue']
; All blue after exposure 6
if (nexp gt 6) then rcol=[[rcol],replicate('blue',nexp-6)]

; Detector and other parameters
det_dkcurr=float(yanny_par(hdr,'det_dkcurr'))
det_rn=float(yanny_par(hdr,'det_rn'))
ifu_trans=float(yanny_par(hdr,'ifu_trans'))

; Slice setup
det_pixscale=float(yanny_par(hdr,'det_pixscale'))
slicewidth=float(yanny_par(hdr,'slicewidth'))
nslice=fix(yanny_par(hdr,'nslices'))
ninslice=fix(float(yanny_par(hdr,'slicelength'))/det_pixscale)
; FOV size in arcsec
XfovAS=float(yanny_par(hdr,'slicelength'))
YfovAS=nslice*slicewidth
; Input grid size is FOV + 1 arcsec for padding
XgridAS=XfovAS+1.
YgridAS=YfovAS+1.

; Read in wavelength solution in microns
wavefile=yanny_par(hdr,'wavesol')
data=read_ascii(concat_dir(getenv('WFIRST_DIR'),wavefile))
data1=data.field1
wave=data1[1,*]/1e4
nwave=n_elements(wave)
; Define dlam which is width of each channel in Angstroms
dwave=fltarr(nwave)
dwave[0]=(wave[1]-wave[0])*1e4
for i=1,nwave-1 do dwave[i]=(wave[i]-wave[i-1])*1e4


; Input image is a star (blackbody) on a 0.015 arcsec grid
inppixsize=0.015; arcsec (use this because it means the pixel scales
; of 0.075'' or 0.15'' are exactly 5 and 10 sim units respectively)
Xgrid=XgridAS/inppixsize
Ygrid=YgridAS/inppixsize
inpcube=fltarr(fix(Xgrid),fix(Ygrid),nwave)
inpcube[Xgrid/2,Ygrid/2,*]=1.

; Read in spectrum from a file
spectype=yanny_par(hdr,'spectype')
specfile=yanny_par(hdr,'spectrum')
if (spectype eq 'fits') then begin
  inp=readfits(concat_dir(getenv('WFIRST_DIR'),specfile))
  inpwave=inp[0,*]
  inpflux=inp[1,*]
  fluxnew=interpol(inpflux,inpwave,wave)
endif
if (spectype eq 'txt') then begin
  inp=read_ascii(concat_dir(getenv('WFIRST_DIR'),specfile))
  data1=inp.field1
  inpwave=data1[0,*]/1e4
  inpflux=data1[1,*]
  fluxnew=interpol(inpflux,inpwave,wave)
endif
; Stick it into the cube
for i=0,nwave-1 do $
  inpcube[*,*,i]=inpcube[*,*,i]*fluxnew[i]
writefits,concat_dir(outdir,'inpcube.fits'),inpcube

; Energy to photon conversion for the cube
energy=hc/(wave/1e4)
; We're using a rough guess at A, but it doesn't matter
; since it will divide out later
factor=A/energy*dwave
; Convert cube to units of photons/s/channel
for i=0,nwave-1 do inpcube[*,*,i]=inpcube[*,*,i]*factor[i]


; Blur input cube by PSF.
; Assume diffraction limited Gaussian
;obscube=inpcube
;fwhm=1.22*wave*1e-6/2.36 * 180./!PI*3600. ; fwhm in arcsec, 2.36m primary
;for i=0,nwave-1 do begin
;  obscube[*,*,i]=filter_image(inpcube[*,*,i],fwhm_gaussian=fwhm[i]/inppixsize)
;endfor
;writefits,'out/blurcube.fits',obscube

; Blur input cube by PSF.
; Assume diffraction limited airy function
;obscube=inpcube
;airmin=1.22*wave*1e-6/2.36 * 180./!PI*3600. ; airy minimum in arcsec, 2.36m primary
;for i=0,nwave-1 do begin
;  psf=airy(Xgrid,Ygrid,Xgrid/2,Ygrid/2,airmin[i]/inppixsize,1.)
;  obscube[*,*,i]=convolve(inpcube[*,*,i],psf)
;endfor
;writefits,concat_dir(outdir,'blurcube.fits'),obscube

; Blur input cube by PSF
; Assume WebbPSF input
; Read in the PSF reference cube
psffile1=yanny_par(hdr,'psfcube')
psffile2=yanny_par(hdr,'psfwave')
psfcube=readfits(concat_dir(getenv('WFIRST_DIR'),psffile1),psfhdr)
; Read in the PSF reference wavelengths in microns
psfwave=readfits(concat_dir(getenv('WFIRST_DIR'),psffile2))
; Interpolate PSF cube to desired spatial sampling,
; conserve total flux
psf_pscl=float(fxpar(psfhdr,'PIXELSCL'))/inppixsize
newxsize=fix((size(psfcube))[1]*psf_pscl)
newysize=fix((size(psfcube))[2]*psf_pscl)
psfcube=congrid(psfcube,newxsize,newysize,(size(psfcube))[3],/interp,/center)/psf_pscl/psf_pscl
; Interpolate PSF cube to desired wavelength sampling
psfcube=resample_psf_cube(psfcube,psfwave,wave)
; Trim PSF cube to have dimension <= input object cube
if (newxsize gt fix(Xgrid)) then psfcube=psfcube[(newxsize-fix(Xgrid))/2:(newxsize-fix(Xgrid))/2+fix(Xgrid)-1,*,*]
if (newysize gt fix(Ygrid)) then psfcube=psfcube[*,(newysize-fix(Ygrid))/2:(newysize-fix(Ygrid))/2+fix(Ygrid)-1,*]
; New cube is convolved at each wavelength
obscube=inpcube
for i=0,nwave-1 do $
  obscube[*,*,i]=convolve(inpcube[*,*,i],psfcube[*,*,i])
writefits,concat_dir(outdir,'blurcube.fits'),obscube

; Slice it up
nspec=nslice*ninslice
xarr=fltarr(nspec*nexp)
yarr=fltarr(nspec*nexp)
; Define a footprint mask for each spectrum element for each exposure
specmask=intarr(Xgrid,Ygrid,nexp)
specmask[*,*,*]=-1
print,'Making footprint masks'
for p=0,nexp-1 do begin
  for j=0,nslice-1 do begin
    for k=0,ninslice-1 do begin
      index=j*ninslice+k
      ; X position relative to array center in arcsec
      thisx=k*det_pixscale-XfovAS/2.+exposures[p].xdither
      thisy=j*slicewidth-YfovAS/2.+exposures[p].ydither
      xarr[index+nspec*p]=thisx
      yarr[index+nspec*p]=thisy

      xstart=(thisx-det_pixscale/2.)/inppixsize+Xgrid/2.
      xstop=(thisx+det_pixscale/2.)/inppixsize+Xgrid/2.
      ystart=(thisy-slicewidth/2.)/inppixsize+Ygrid/2.
      ystop=(thisy+slicewidth/2.)/inppixsize+Ygrid/2.

      if ((xstart ge 0)and(xstop lt Xgrid)and(ystart ge 0)and(ystop lt Ygrid)) then $
        specmask[xstart:xstop,ystart:ystop,p]=index
    endfor
  endfor

  ; Make a mesh overlay for this footprint
  thismask=specmask[*,*,p]
  mtemp=makemesh(thismask,ninslice,nslice,strcompress(concat_dir(outdir,'mesh_exp'+string(p)+'.reg'),/remove_all),color=rcol[p])
endfor

; Define zodiacal light
zodi=mrdfits(concat_dir(getenv('WFIRST_DIR'),'RefFiles/zodi_medium.fits'),1)
; Wavelength in microns
zodi_wave=zodi.wavelength
; Surface brightness in MJy/sr
zodi_sb=zodi.sb
; Convert to erg/s/cm2/Hz/arcsec^2
zodi_sb=zodi_sb*(1e-17)*2.3504*1e-11
; Convert to on-sky spatial element
; erg/s/Hz/spaxel seen by telescope
zodi_sb=zodi_sb*slicewidth*det_pixscale*A; Inconsistency on telescope area, but
; it doesn't matter cuz we'll divide it out again immediately later.
; Convert to erg/s/Angstrom/spaxel
zodi_sb=zodi_sb*c/(zodi_wave*1e-4)/(zodi_wave*1e-4)*1e-8
; Resample onto wavelength grid
zodi_flux=interpol(zodi_sb,zodi_wave,wave)
; Energy to photon conversion and integrate over
; width of each spectral element to photon/s/channel
zodi_flux=zodi_flux/energy*dwave

; Thermal background
thermal=mrdfits(concat_dir(getenv('WFIRST_DIR'),'RefFiles/thermal.fits'),0)
thermal_w=thermal[0,*]; Wave in Angstroms
thermal_f=thermal[1,*]; erg/s/cm2/Ang/arcsec^2
; Resample onto wavelength grid
thermal_flux=interpol(thermal_f,thermal_w/1e4,wave)
; Convert to erg/s/Ang/spaxel
thermal_flux=thermal_flux*slicewidth*det_pixscale*A
; Energy to photon conversion and integrate over
; width of each spectral element to photon/s/channel
thermal_flux=thermal_flux/energy*dwave

; Define effective area in cm2: Area multiplied by throughput
eafile=yanny_par(hdr,'effarea')
data=read_ascii(concat_dir(getenv('WFIRST_DIR'),eafile),data_start=3)
data1=data.field1
ea_wave=data1[0,*]/1e4
ea_val=data1[1,*]
; Interpolate to effective area in cm2
effarea=interpol(ea_val,ea_wave,wave)*1e4
tput_all=effarea/A

; Get sliced spectra and add noise sources
print,'Making ',nspec*nexp,' spectra'
; Make ncovar simulated spectra with noise
; (ncovar >1 helps constrain covariance effects)
ncovar=fix(yanny_par(hdr,'ncovar'))
flux=fltarr(nspec*nexp,nwave,ncovar)
; Also make ivar spectrum, and a truth spectrum (no noise)
ivar=fltarr(nspec*nexp,nwave); Ivar of the simulated spectra
flux1=fltarr(nspec*nexp,nwave); The truth spectrum (no noise)

for p=0,nexp-1 do begin
  ; Background spectrum in counts/channel
  ; from zodiacal and dark current
  bgspec=(zodi_flux*tput_all + thermal_flux*tput_all + det_dkcurr)*exposures[p].exptime
  ;bgspec=(zodi_flux*tput_all + det_dkcurr)*exposures[p].exptime
  for q=0,nspec-1 do begin
    ; Where does this spatial element cover?
    index=where(specmask eq q,nindex)
    xyindex=array_indices(specmask,index)
    temp=fltarr(nindex,nwave)
    ; Sum fluxes in this footprint
    for i=0,nindex-1 do temp[i,*]=obscube[xyindex[0,i],xyindex[1,i],*]
    ; Science spectrum in counts/channel
    scispec=total(temp,1)*exposures[p].exptime*tput_all

    ; Noise vector based on total signal plus readnoise
    noisevec=sqrt(bgspec+scispec+det_rn*det_rn)
    ; Random realization of resulting noise spectrum
    thisspec=randomn(systime_seed,nwave,ncovar)
    ; Add in science spectrum
    for i=0,ncovar-1 do thisspec[*,i]=thisspec[*,i]*noisevec+scispec

    ; Flux calibrate the spectra and the ivar to units of 1e-17 erg/s/cm2/Ang
    for i=0,ncovar-1 do $
      flux[q+p*nspec,*,i]=thisspec[*,i]/exposures[p].exptime/tput_all/(A/energy*dwave)*1e17
    flux1[q+p*nspec,*]=scispec/exposures[p].exptime/tput_all/(A/energy*dwave)*1e17
    temp=noisevec/exposures[p].exptime/tput_all/(A/energy*dwave)*1e17
    ivar[q+p*nspec,*]=1.D/(temp*temp)
  endfor
endfor

; Reconstruct cube1 from the first exposure only
print,'Making simple cube'
simpcube=fltarr(ninslice,nslice,nwave,ncovar)
simpcube1=fltarr(ninslice,nslice,nwave)
simpivar=fltarr(ninslice,nslice,nwave)
for j=0,nslice-1 do begin
  for k=0,ninslice-1 do begin
    index=j*ninslice+k
    for i=0,ncovar-1 do simpcube[k,j,*,i]=flux[index,*,i]
    simpcube1[k,j,*]=flux1[index,*]
    simpivar[k,j,*]=ivar[index,*]
  endfor
endfor
writefits,concat_dir(outdir,'simpcube.fits'),simpcube

;simp_sum=fltarr(nwave)
;simp_sum1=fltarr(nwave)
;simp_sigma=fltarr(nwave)
;aper=0.25
;aperx=round(aper/det_pixscale)
;apery=round(aper/slicewidth)
; Note that this isn't matched to the circular extraction from later cube!
;for i=0,nwave-1 do simp_sum[i]=total(simpcube[17:21,10:11,i])
;for i=0,nwave-1 do simp_sum1[i]=total(simpcube1[17:21,10:11,i])
;for i=0,nwave-1 do simp_sigma[i]=sqrt(total((1./simpivar)[17:21,10:11,i]))
;for i=0,nwave-1 do simp_sum[i]=total(simpcube[19-aperx:19+aperx,10-apery:10+apery,i])
;for i=0,nwave-1 do simp_sum1[i]=total(simpcube1[19-aperx:19+aperx,10-apery:10+apery,i])
;for i=0,nwave-1 do simp_sigma[i]=sqrt(total((1./simpivar)[19-aperx:19+aperx,10-apery:10+apery,i]))

; sigma in a circular aperture should be about 0.9 of sigma in a
; square aperture

; Reconstruct cube2
outppixsize=0.05; arcsec
XOsize=fix(Xgrid*inppixsize/outppixsize)
YOsize=fix(Ygrid*inppixsize/outppixsize)
; ncovar realizations of recovered noised data cube
reccube=fltarr(XOsize,YOsize,nwave,ncovar)
; 1 realization of perfect SNR data cube
reccube1=fltarr(XOsize,YOsize,nwave)
; 1 realization of the noise cube (perfect)
recivar=fltarr(XOsize,YOsize,nwave)

; Scaling factor is the ratio between a spatial
; element on sky and our output spaxel size
scalefac=(outppixsize*outppixsize)/(slicewidth*det_pixscale)
rdiffrac=fltarr(nwave)
print,'Making full cube'
for i=0,nwave-1 do begin
  print,i
  ; Perfect fluxes of elements at this wave
  theseflux1=flux1[*,i]
  ; Ivar of all elements at this wave
  theseivar=ivar[*,i]
  ; Define grid
  thesex=XOsize/2.+xarr/outppixsize
  thesey=YOsize/2.+yarr/outppixsize
  ; Choose cutoff distance of a slice width, scale distance of a det pixel
  rdiffrac[i]=wave[i]*1e-4*1.22/teldiam*180./!PI*3600.; rough diffraction limit in arcsec
  rlim=slicewidth/outppixsize
  rscale=det_pixscale/outppixsize/2.

  ; Noise+science cubes
  for j=0,ncovar-1 do begin
    theseflux=flux[*,i,j]
    thisimg=ml_griddata(thesex,thesey,theseflux,[XOsize,YOsize],rlim,rscale,scale=scalefac,ivar=theseivar,invarimg=invarimg,/noflag)
    reccube[*,*,i,j]=thisimg
    if (j eq 0) then recivar[*,*,i]=invarimg
  endfor

  ; Infinite SNR truth cube
  thisimg1=ml_griddata(thesex,thesey,theseflux1,[XOsize,YOsize],rlim,rscale,scale=scalefac,/noflag)
  reccube1[*,*,i]=thisimg1
endfor
writefits,concat_dir(outdir,'reccube.fits'),reccube

; Determine S/N by aperture extraction
finalspec=fltarr(nwave,ncovar)
finalspec1=fltarr(nwave)
idealrms=fltarr(nwave)
aperrad=fltarr(nwave)
measfwhm=fltarr(nwave)
print,'Extracting spectra'
; Define aperture for spectral extraction
nreff=float(yanny_par(hdr,'nreff'))
for i=0,nwave-1 do begin
  ; Fit 2d gaussian model to the noiseless image to get good fits
  ; and do aperture phot on this noiseless image too
  slice1=reccube1[*,*,i]
  model=gauss2dfit(slice1,psf,/tilt)
  xcen=psf[4]
  ycen=psf[5]
  measfwhm[i]=(psf[2]+psf[3])/2.*2.35

  if (nreff ge 0) then begin
    aperrad[i]=measfwhm[i]/2.35*nreff ; nreff effective radii
  endif else begin
    aperrad[i]= -nreff/outppixsize
  endelse

  ; Aperture photometry on the infinite SNR cube
  aper,slice1,xcen,ycen,theflux,theeflux,thesky,theesky,1.,aperrad[i],/flux,/nan,setskyval=0.,/silent
  finalspec1[i]=theflux

  ; Aperture photometry on the noise+science cubes
  for j=0,ncovar-1 do begin
    slice=reccube[*,*,i,j]
    aper,slice,xcen,ycen,theflux,theeflux,thesky,theesky,1.,aperrad[i],/flux,/nan,setskyval=0.,/silent
    finalspec[i,j]=theflux
  endfor

  ; Aperture photometry on the sigma^2 image
  slice_sig2=1./recivar[*,*,i]
  aper,slice_sig2,xcen,ycen,thess,theess,theskyss,theeskyss,1.,aperrad[i],/flux,/nan,setskyval=0.,/silent
  idealrms[i]=sqrt(thess);*(det_pixscale*slicewidth/outppixsize/outppixsize)
endfor
; SNR is ratio of truth signal to ideal noise
snr=finalspec1/idealrms

; Need to account for covariance
; Do this by comparing ideal noise to actual noise realizations after resampling
; using a moving boxcar
realrms=fltarr(nwave,ncovar)
rmsratio=fltarr(nwave,ncovar)
boxsize=25
for j=0,ncovar-1 do begin
  temp=finalspec1-finalspec[*,j]
  for i=boxsize,nwave-boxsize-1 do realrms[i,j]=sqrt((moment(temp[i-boxsize:i+boxsize]))[1]) 
  ; Average ratio of realrms/finalrms is the covariance correction factor
  rmsratio[*,j]=realrms[*,j]/idealrms
endfor
allratio=total(rmsratio,2)/float(ncovar)
covarfac=(moment(allratio[boxsize:nwave-boxsize-1]))[0]
; Apply to derived SNR
snr=snr/covarfac


; Write some intermediary files
counts_spec=finalspec1/(1./exposures[0].exptime/tput_all/(A/energy*dwave)*1e17)
counts_noise=idealrms/(1./exposures[0].exptime/tput_all/(A/energy*dwave)*1e17)*covarfac
writefits,concat_dir(outdir,'counts_spec.fits'),counts_spec
writefits,concat_dir(outdir,'counts_noise.fits'),counts_noise

; Write out to FITS file
snrout=fltarr(2,nwave)
snrout[0,*]=wave
snrout[1,*]=snr
writefits,concat_dir(outdir,'snr.fits'),snrout

    set_plot,'ps'
    origp=!p
    origx=!x
    origy=!y
    !p.multi=[0,1,2]
    !p.thick=0 & !x.thick=0 & !y.thick=0 & !p.charthick=0
    plotname=concat_dir(outdir,'snr.ps')
    dfpsplot, plotname, /color
    loadct,39

plot,wave,snr,xtitle='Wavelength (microns)',ytitle='SNR',title='Supernova z=0.5, 500 sec',thick=2,xthick=2,ythick=2,charthick=2,charsize=1.5,xrange=[0.5,2.1],xstyle=1

    dfpsclose
    ; Convert to PDF and remove old PS file
    spawn, strcompress('ps2pdf '+plotname+' '+ml_strreplace(plotname,'.ps','.pdf'))
    spawn, strcompress('rm -f '+plotname)
    !p.multi=0
    !p=origp
    !x=origx
    !y=origy
set_plot,'x'

openw,lun, concat_dir(outdir,'results.txt'), /get_lun, width=250
printf,lun,'# wave aperture counts noise snr'
for i=0,nwave-1 do $
  printf,lun,wave[i]*1e4,aperrad[i]*outppixsize,counts_spec[i],counts_noise[i],snr[i]

close,lun
free_lun,lun

stop
return
end
