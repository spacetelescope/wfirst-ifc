; Return a 2d airy function array
; Give it an array size, function center
; FWHM (distance to first minimum), and total flux.
; Be warned that total flux scaling may not work if substantial
; flux falls off of the field.
; obsairy is obscured airy disk for a given eta obscuration

function airy,xsize,ysize,xcen,ycen,min1rad,totalflux,eta=eta,obsairy=obsairy

if (~keyword_set(eta)) then eta=0.3

; min1rad is the radius in pixel at which to have the first
; zero of the airy function.  Note that beselj(radius,1) has
; first zero where radius is 3.831705

img=fltarr(xsize,ysize)
img2=fltarr(xsize,ysize)

for i=0,xsize-1 do begin
  for j=0,ysize-1 do begin
    radius=sqrt((i-xcen)^2+(j-ycen)^2)
    ; Scale radius to correct situation
    radius=radius/min1rad*3.831705
    ; Safety case, problems with radius=0
    if (radius lt 1e-4) then radius=1e-4
    value=(beselj(radius,1)/radius)^2
    img[i,j]=value

    ; Obscured airy function
    value2=(beselj(radius,1)/radius - eta*beselj(radius*eta,1)/radius)^2
    img2[i,j]=value2
  endfor
endfor

; Scale to desired total
sum=total(img)
img=img*totalflux/sum

sum2=total(img2)
obsairy=img2*totalflux/sum2

return,img
end
