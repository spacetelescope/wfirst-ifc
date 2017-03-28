pro TRplots

rootdir='./'

; Test 1
; Law result
temp=read_ascii(concat_dir(rootdir,'Run23/results_Run23.txt'),data_start=1)
t1_law=temp.field1
; Rubin result
temp=read_ascii(concat_dir(rootdir,'Run23/rubinfeb12_blackbody_R07.txt'),data_start=1)
t1_rubin=temp.field1
; Klaus result
t1_klaus=mrdfits(concat_dir(rootdir,'Run23/klaus_signaltonoise_test1_feb19h.fits'),1)

; Test 2
; Law result
temp=read_ascii(concat_dir(rootdir,'Run24/results_Run24.txt'),data_start=1)
t2_law=temp.field1
; Rubin result
temp=read_ascii(concat_dir(rootdir,'Run24/rubinfeb12_mid-z_single_R07.txt'),data_start=1)
t2_rubin=temp.field1
; Klaus result
t2_klaus=mrdfits(concat_dir(rootdir,'Run24/klaus_signaltonoise_test2_feb19h.fits'),1)

; Test 3
; Law result
temp=read_ascii(concat_dir(rootdir,'Run25/results_Run25.txt'),data_start=1)
t3_law=temp.field1
; Rubin result
temp=read_ascii(concat_dir(rootdir,'Run25/rubinfeb12_mid-z_merged_R07.txt'),data_start=1)
t3_rubin=temp.field1
; Klaus result
t3_klaus=mrdfits(concat_dir(rootdir,'Run25/klaus_signaltonoise_test3_feb19h.fits'),1)

; Test 4
; Law result
temp=read_ascii(concat_dir(rootdir,'Run26/results_Run26.txt'),data_start=1)
t4_law=temp.field1
; Rubin result
temp=read_ascii(concat_dir(rootdir,'Run26/rubinfeb12_high-z_R07.txt'),data_start=1)
t4_rubin=temp.field1
; Klaus result
t4_klaus=mrdfits(concat_dir(rootdir,'Run26/klaus_signaltonoise_test4_feb19h.fits'),1)

; Resampling
t1_rubinres=interpol(t1_rubin[5,*],t1_rubin[0,*],t1_law[0,*])
r1_rubin=t1_rubinres/t1_law[4,*]
t1_klausres=interpol(t1_klaus.sn,t1_klaus.wavelength*1e4,t1_law[0,*])
r1_klaus=t1_klausres/t1_law[4,*]
t2_rubinres=interpol(t2_rubin[5,*],t2_rubin[0,*],t2_law[0,*])
r2_rubin=t2_rubinres/t2_law[4,*]
t2_klausres=interpol(t2_klaus.sn,t2_klaus.wavelength*1e4,t2_law[0,*])
r2_klaus=t2_klausres/t2_law[4,*]
t3_rubinres=interpol(t3_rubin[5,*],t3_rubin[0,*],t3_law[0,*])
r3_rubin=t3_rubinres/t3_law[4,*]
t3_klausres=interpol(t3_klaus.sn,t3_klaus.wavelength*1e4,t3_law[0,*])
r3_klaus=t3_klausres/t3_law[4,*]
t4_rubinres=interpol(t4_rubin[5,*],t4_rubin[0,*],t4_law[0,*])
r4_rubin=t4_rubinres/t4_law[4,*]
t4_klausres=interpol(t4_klaus.sn,t4_klaus.wavelength*1e4,t4_law[0,*])
r4_klaus=t4_klausres/t4_law[4,*]
ravg_klaus=(r1_klaus+r2_klaus+r3_klaus+r4_klaus)/4.
ravg_rubin=(r1_rubin+r2_rubin+r3_rubin+r4_rubin)/4.

set_plot,'ps'
  origp=!p
  origx=!x
  origy=!y
  !p.multi=[0,1,2]
  !p.thick=0 & !x.thick=0 & !y.thick=0 & !p.charthick=0
dfpsplot, 'wfirst_trplots.ps', /color,ysize=7., xsize=5.
loadct,39

xtemp=[0.7,1.0]

plot,t1_law[0,*]/1e4,t1_law[4,*],xrange=[0.4,2.01],/xstyle,yrange=[0,140],/ystyle,charsize=1.5,xtitle='Wavelength (microns)',ytitle='SNR',xthick=5,ythick=5,charthick=5,thick=5,/nodata
oplot,t1_rubin[0,*]/1e4,t1_rubin[5,*],thick=7,color=80
oplot,t1_klaus.wavelength,t1_klaus.sn,thick=7,color=250
oplot,t1_law[0,*]/1e4,t1_law[4,*],thick=7
oplot,xtemp,[60,60],thick=7
oplot,xtemp,[40,40],thick=7,color=250
oplot,xtemp,[20,20],thick=7,color=80
xyouts,1.1,60,'LAW',charthick=5
xyouts,1.1,40,'PANDEIA',color=250,charthick=5
xyouts,1.1,20,'RUBIN',color=80,charthick=5


plot,t2_law[0,*]/1e4,t2_law[4,*],xrange=[0.4,2.01],/xstyle,yrange=[0,6.99],/ystyle,charsize=1.5,xtitle='Wavelength (microns)',ytitle='SNR',xthick=5,ythick=5,charthick=5,thick=5,/nodata
oplot,t2_rubin[0,*]/1e4,t2_rubin[5,*],thick=7,color=80
oplot,t2_klaus.wavelength,t2_klaus.sn,thick=7,color=250
oplot,t2_law[0,*]/1e4,t2_law[4,*],thick=7

plot,t3_law[0,*]/1e4,t3_law[4,*],xrange=[0.4,2.01],/xstyle,yrange=[0,16],/ystyle,charsize=1.5,xtitle='Wavelength (microns)',ytitle='SNR',xthick=5,ythick=5,charthick=5,thick=5,/nodata
oplot,t3_rubin[0,*]/1e4,t3_rubin[5,*],thick=7,color=80
oplot,t3_klaus.wavelength,t3_klaus.sn,thick=7,color=250
oplot,t3_law[0,*]/1e4,t3_law[4,*],thick=7

plot,t4_law[0,*]/1e4,t4_law[4,*],xrange=[0.4,2.01],/xstyle,yrange=[0,1.5],/ystyle,charsize=1.5,xtitle='Wavelength (microns)',ytitle='SNR',xthick=5,ythick=5,charthick=5,thick=5,/nodata
oplot,t4_rubin[0,*]/1e4,t4_rubin[5,*],thick=7,color=80
oplot,t4_klaus.wavelength,t4_klaus.sn,thick=7,color=250
oplot,t4_law[0,*]/1e4,t4_law[4,*],thick=7

dfpsclose
  ; Convert to PDF and remove old PS file
  ;spawn, strcompress('ps2pdf '+pname + ' ' + pnamepdf)
  !p.multi=0
  !p=origp
  !x=origx
  !y=origy
cleanplot



set_plot,'ps'
  origp=!p
  origx=!x
  origy=!y
  !p.multi=[0,1,2]
  !p.thick=0 & !x.thick=0 & !y.thick=0 & !p.charthick=0
dfpsplot, 'wfirst_example.ps', /color,ysize=7., xsize=5.
loadct,39

wave=readfits(concat_dir(rootdir,'Run25/snr.fits'))
wave=wave[0,*]
ideal=readfits(concat_dir(rootdir,'Run25/finalspec1.fits'))
real=readfits(concat_dir(rootdir,'Run25/finalspec.fits'))
real=real[*,0]

;plot,wave,real
plot,wave,real,xrange=[0.4,2.01],/xstyle,yrange=[0,0.005],/ystyle,charsize=1.5,xtitle='Wavelength (microns)',ytitle='Flux (mJy)',xthick=5,ythick=5,charthick=5,thick=5,/nodata
oplot,wave,ideal-0.0005,thick=5,color=250
oplot,wave,real,thick=7

dfpsclose
  ; Convert to PDF and remove old PS file
  ;spawn, strcompress('ps2pdf '+pname + ' ' + pnamepdf)
  !p.multi=0
  !p=origp
  !x=origx
  !y=origy
cleanplot


set_plot,'ps'
  origp=!p
  origx=!x
  origy=!y
  !p.multi=[0,1,2]
  !p.thick=0 & !x.thick=0 & !y.thick=0 & !p.charthick=0
dfpsplot, 'wfirst_example2.ps', /color,ysize=7., xsize=5.
loadct,39

wave=readfits(concat_dir(rootdir,'Run24/snr.fits'))
wave=wave[0,*]
ideal=readfits(concat_dir(rootdir,'Run24/finalspec1.fits'))
real=readfits(concat_dir(rootdir,'Run24/finalspec.fits'))
real=real[*,0]

;plot,wave,real
plot,wave,real,xrange=[0.4,2.01],/xstyle,yrange=[0,0.005],/ystyle,charsize=1.5,xtitle='Wavelength (microns)',ytitle='Flux (mJy)',xthick=5,ythick=5,charthick=5,thick=5,/nodata
oplot,wave,ideal-0.0005,thick=5,color=250
oplot,wave,real,thick=7

dfpsclose
  ; Convert to PDF and remove old PS file
  ;spawn, strcompress('ps2pdf '+pname + ' ' + pnamepdf)
  !p.multi=0
  !p=origp
  !x=origx
  !y=origy
cleanplot




stop
return
end

