function LOWPASS, spectrum, CutOff, filter=filter
;+
; NAME:	
;       LOWPASS
; PURPOSE:
;       Apply a low-pass filter to spectrum
; CATEGORY:
; CALLING SEQUENCE:
;       newspect = LOWPASS( spectrum, Cutoff [,filter=filter] )
; INPUTS:
;       SPECTRUM  --  Fourier spectrum
;       CUTOFF    --  Cut Off frequency (float between 0.0 ... 1.0)
; KEYWORD PARAMETERS:
;       FILTER    --  The Low-Pass filter itself
; OUTPUTS:
;       NEWSPECT  --  Filtered power spectrum
; COMMON BLOCKS:
; NOTES:
;       
; MODIFICATION HISTORY:
;       Adopted from Pete Nisenson, Han Uitenbroek, March 1993
;-
  c = [0.074, 0.302, 0.233, 0.390]

  sz = size(spectrum)  &  nx = sz(1)
  if (sz(0) eq 1) then begin
    filter = fltarr(nx)
    r = findgen(nx/2)
    if (nx mod 2) then $
      r = [r, (nx-1)/2, reverse(r)] / (CutOff*nx) $
    else $
      r = [r, reverse(r)] / (CutOff*nx)
  endif else begin
    ny = sz(2)
    filter = fltarr(nx, ny)
    r0 = sqrt(nx^2 + ny^2)
    r  = dist(nx, ny) / (CutOff*r0)
  endelse

  index = where(r le 1.5)
  r1 = (1.0 - r(index)^2)
  r2 = r1^2
  r3 = r2 * r1

  filter(index) = (c[0] + r1*c[1] + r2*c[2] + r3 * c[3]) > 0.0

  return, spectrum*filter
end
