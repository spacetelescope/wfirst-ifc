function makeblackbody,mag=mag,T=T,wave=inpwave

; Default source has 19th magnitude Vega at 5556 Angstroms
if (~keyword_set(mag)) then mag=19.00

; Default temperature
if (~keyword_set(T)) then T=5000.;Kelvin

norm=3.39e-9 ;erg/s/cm2/Angstrom for Vega at 5556 Angstroms
norm=norm*10.^(-mag/2.5)
inpwave=(findgen(16000)+4000.)/1e8; Every Angstrom from 4000 to 2 microns, in cm
  PI=3.1415926536
  radpas=4.8481e-6; Radians per arcsecond
  c=2.9979e10; Speed of light (cm/s)
  hc=1.986e-16; Planck constant times speed of light (erg cm)
  k=1.381e-16; Boltzmann constant (erg/K)
  A=PI*236.*236./4.; Area of primary (cm^2)
flux=2*hc*c/(inpwave^5)/(exp(hc/(inpwave*k*T))-1.D)

temp=abs(inpwave-5556.e-8)
junk=min(temp,minat)

; Normalize flux
flux=flux/flux[minat]*norm
; convert wavelength to microns
inpwave=inpwave*1e4

return,flux
end
