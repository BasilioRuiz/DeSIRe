PRO writeBF, atom, fileName

  openw, unit, fileName, /GET_LUN
  printf, unit, FORMAT='(A, //A)', "# Photoionization rates from TopBase", $
   "# j   i  alpha [m^2]    Nlambda     Wavel. Dep.   lamb_min [nm]"

  Nrad = atom.Nline + atom.Ncont
  FOR kr=atom.Nline, Nrad-1 DO BEGIN
    rt = atom.transition(kr)

    CASE (rt.shape) OF
      3: shapeStr = "HYDROGENIC"
      4: shapeStr = " EXPLICIT "
    ENDCASE

    printf, unit, FORMAT='(/"#   ", A)', atom.labels(rt.i)
    printf, unit, FORMAT='(I3, X, I3, 3X, E10.4, 5X, I3, 7X, A, 5X, F6.1)', $
     rt.j, rt.i, rt.strength, rt.Nlambda, shapeStr, rt.lambdaMin

    IF (rt.shape EQ 4) THEN BEGIN
      lambda=*rt.lambda_ptr
      alpha=*rt.alpha_ptr 
      FOR la=0, rt.Nlambda-1 DO $
       printf, unit, FORMAT='(X, F9.3, 5X, E10.4)', lambda[la], alpha[la]
    ENDIF
  ENDFOR

  free_lun, unit
END

PRO adjustAlpha, atom

  cLight   = 2.99792458E+08
  hPlanck  = 6.626176E-34
  ERydb    = 2.17991E-18
  one_eV   = 1.60219E-19

  MEGABARN_TO_M2 = 1.0E-22
  HYDROGENIC = 3
  EXPLICIT   = 4

  openw, scriptFile, 'script.txt', /GET_LUN

  Nrad = atom.Nline + atom.Ncont
  FOR kr=atom.Nline, Nrad-1 DO BEGIN
    rt = atom.transition(kr)
    alpha0 = rt.strength
    shape  = rt.shape
    Nlamb  = rt.Nlambda  &  lambdaMin = rt.lambdaMin
    lambda=*rt.lambda_ptr
    alpha=*rt.alpha_ptr

    choice = 1  &  plotNew = 0
    WHILE (choice NE 0) DO BEGIN
      IF (rt.shape EQ EXPLICIT) THEN BEGIN
        plot, lambda, alpha, XTITLE='Wavelength [nm]', $
         YTITLE='Cross-section [m!U2!N]', /YLOG
        oplot, lambda, alpha0 * (lambda/rt.lambda0)^3, COLOR=175B
        oplot, lambda, alpha, PSYM=5, COLOR=200B
      ENDIF ELSE BEGIN
        lambda = rt.lambdaMin + (rt.lambda0 - rt.lambdaMin) * $
         findgen(rt.Nlambda) / (rt.Nlambda - 1)
        alpha_hydr = alpha0 * (lambda/rt.lambda0)^3
        plot, lambda, alpha_hydr, XTITLE='Wavelength [nm]', $
         YTITLE='Cross-section [m!U2!N]]', /YLOG
        oplot, lambda, alpha_hydr, PSYM=5, COLOR=200B
      ENDELSE
      xyouts, 0.15, 0.9, atom.labels(rt.i), /NORM, CHARSIZE=1.2

      IF (shape EQ EXPLICIT) THEN BEGIN
        IF (plotNew) THEN $
         oplot, lambdaNew, alphaNew, COLOR=150B

        read, 'Choose modification, (-1=quit, 0=nothing, 1=resample, ' + $
         '2=fft, 3=hydrogenic) > ', choice
        choice = -1 > choice < 3
      ENDIF ELSE BEGIN
        IF (plotNew) THEN BEGIN
          alpha_hydr = alpha0 * (lambdaNew/rt.lambda0)^3
          oplot, lambdaNew, alpha_hydr, COLOR=150B
          oplot, lambdaNew, alpha_hydr, PSYM=5
        ENDIF
        read, 'Choose modification, (-1=quit, 0=nothing, 1=resample) > ', $
         choice
        choice = -1 > choice < 1
      ENDELSE
      CASE (choice) OF
        -1: GOTO, break
        0:

        1:  BEGIN
          read, 'Give Nlamb > ', Nlamb
          Nlamb = Nlamb > 2
          IF (shape EQ HYDROGENIC) THEN BEGIN
            read, 'Give lambdaMin [nm] > ', lambdaMin
            lambdaMin = lambdaMin < (rt.lambda0 - 1.0)
            lambdaNew = lambdaMin + (rt.lambda0 - lambdaMin) * $
             findgen(Nlamb) / (Nlamb - 1)
          ENDIF ELSE BEGIN
            sample = rt.Nlambda / Nlamb
            lambdaNew = lambda(sample*findgen(Nlamb))
            alphaNew  = alpha(sample*findgen(Nlamb))
            lambdaMin = lambdaNew(Nlamb-1)
          ENDELSE
          plotNew = 1

          printf, scriptFile, $
           FORMAT='(A1, A20, A1, 2X, "RESAMPLE", 2X, I4)', $
           "'", atom.labels(rt.i), "'", Nlamb
        END

	2: BEGIN
          read, 'Give Nlamb and Rlowpass [0.0, 1.0] > ', Nlamb, Rlowpass
          Nlamb = Nlamb > 2
          Rlowpass = 0.001 > Rlowpass < 0.999

          ;; --- Apply a low pass filter in energy space -- -------------;;

          Nfft = 1
          WHILE (Nfft LT rt.Nlambda) DO Nfft = Nfft * 2
          E = (hPlanck * cLight) / lambda
          E_ft = E(0) + (E(rt.Nlambda-1) - E(0)) * findgen(Nfft) / $
           float(Nfft - 1)
          tabinv, E, E_ft, Eeff
          aa = linear(alpha, Eeff)
 
          ;; --- Symmetrize to avoid aliasing --           ------------- ;;
 
          ft_alpha = fft([reverse(aa), aa], -1)
          alpha_lp = fft(lowpass(ft_alpha, rLowPass), 1)

          lambdaNew = lambda[0] + (lambda[rt.Nlambda-1] - lambda[0]) * $
           findgen(Nlamb)/(Nlamb-1)
          tabinv, (cLight * hPlanck)/E_ft, lambdaNew, leff
          alphaNew = linear(float(alpha_lp(Nfft:*)), leff)
          plotNew = 1

          printf, scriptFile, $
           FORMAT='(A1, A20, A1, 2X, "FFT", 2X, I4, 2X, F7.5)', $
           "'", atom.labels(rt.i), "'", Nlamb, Rlowpass
        END
        3: BEGIN
          shape = HYDROGENIC
          read, 'Give Nlamb and lambdaMin [nm] > ', Nlamb, lambdaMin
          Nlamb = Nlamb > 2
          lambdaMin = lambdaMin < (rt.lambda0 - 1.0)

          lambdaNew = lambdaMin + (rt.lambda0 - lambdaMin) * $
           findgen(Nlamb) / (Nlamb - 1)
          plotNew = 1

          printf, scriptFile, $
           FORMAT='(A1, A20, A1, 2X, "HYDROGENIC", 2X, I4, 2X, F7.1)', $
           "'", atom.labels(rt.i), "'", Nlamb, lambdaMin
        END
      ENDCASE
    ENDWHILE
    rt.shape = shape
    rt.Nlambda = Nlamb
    rt.lambdaMin = lambdaMin
    IF ((rt.shape EQ EXPLICIT)  AND  (plotNew)) THEN BEGIN
      lambdaNew=*rt.lambda_ptr ;, /SET
      alphaNew=*rt.alpha_ptr   ;, /SET
    ENDIF
    atom.transition(kr) = rt
  ENDFOR

break:
  free_lun, scriptFile

  writebf, atom, 'bf.dat'
END
