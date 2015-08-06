; Source has 21th magnitude Vega at 5556 Angstroms
norm=3.39e-9 ;erg/s/cm2/Angstrom for Vega at 5556 Angstroms
mag=20.78
norm=norm*10.^(-mag/2.5)
inpwave=(findgen(15000)+5556.)/1e8; Every Angstrom from 5556 to 2 microns, in cm
  PI=3.1415926536
  radpas=4.8481e-6; Radians per arcsecond
  c=2.9979e10; Speed of light (cm/s)
  hc=1.986e-16; Planck constant times speed of light (erg cm)
  k=1.381e-16; Boltzmann constant (erg/K)
  A=PI*236.*236./4.; Area of primary (cm^2)
T=5000.; kelvin
flux=2*hc*c/(inpwave^5)/(exp(hc/(inpwave*k*T))-1.D)
flux=flux/flux[0]*norm
