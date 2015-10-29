pro makepsfcube

; Assemble WebbPSF frames into a reference psf cube
; so that all in one file
indir='~/Downloads/webbpsf/'
files=file_search(indir,'wfirstifu*')

nfiles=n_elements(files)
frame0=readfits(files[0],hdr)
nx=(size(frame0))[1]
ny=(size(frame0))[2]

pixelscl=fxpar(hdr,'PIXELSCL')

psfcube=fltarr(nx,ny,nfiles)
psfwave=fltarr(nfiles)

for i=0,nfiles-1 do begin
  new=readfits(files[i],hdr)
  thiswave=float(fxpar(hdr,'WAVE0'))*1e6; microns
  psfcube[*,*,i]=new
  psfwave[i]=thiswave
endfor

mkhdr,hdrcube,psfcube
fxaddpar,hdrcube,'PIXELSCL',pixelscl

writefits,'RefFiles/psfcube.fits',psfcube,hdrcube
writefits,'RefFiles/psfwave.fits',psfwave

stop

return
end
