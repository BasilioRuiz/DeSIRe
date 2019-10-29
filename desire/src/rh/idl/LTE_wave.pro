pro LTE_wave,lam0,deltalamini,Nlambda,step,filename

 ;Spectrum example Fe I 6173
 ;LTE_wave,6173.334,-0.5,100L,10e-3,'Fe6173_spectrum_10mA.wave'
 ;lam0=6173.334
 ;deltalamini=-0.5
 ;Nlambda = 100L
 ;step=10e-3
 ;filename='Fe6173_spectrum_10mA.wave'
 
 lambda = dindgen(nlambda)*step+lam0+deltalamini
 lambda=lambda/1e1  ;[nm]
 openw,1,/XDR,filename
 writeu,1,nlambda,airtovacuum(lambda)
 close,1
 
 return
 end