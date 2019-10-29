PRO write1datmosmod, atmos, FILENAME=filename

  CM_TO_M = 1.0E-2
  KM_TO_M = 1.0E+3
  G_TO_KG = 1.0E-3

  openw, unit, filename, /GET_LUN

  printf, unit, atmos.atmosID

  CASE (atmos.scale) OF
    'MASS_SCALE': printf, unit, FORMAT='(2X,"Mass scale")'
    'TAU500_SCALE': printf, unit, FORMAT='(2X,"Tau scale")'
    'GEOMETRIC_SCALE': printf, unit, FORMAT='(2X,"Height scale")'
  ENDCASE

  printf, unit, $
   FORMAT='("*",/"* lg g",/F6.2,/"*",/"* Ndep",/I4,/"*")', $
   alog10(atmos.gravitation), atmos.Ndep


  CASE (atmos.scale) OF
    'MASS_SCALE': BEGIN
      printf, unit, $
       FORMAT='("*lg column mass", 5X,"Temperature",8X,"Ne",9X,"V",14X,"Vturb")'
      
      FOR k=0, atmos.Ndep-1 DO $
       printf, unit, FORMAT='(5E15.6)', $
               alog10(atmos.cmass[k] * CM_TO_M^2 / G_TO_KG), $
               atmos.T[k], atmos.n_elec[k] * CM_TO_M^3 , $
               atmos.v[k]/KM_TO_M, atmos.vturb[k]/KM_TO_M

    END
    'TAU500_SCALE': BEGIN
      printf, unit, $
       FORMAT='("*lg tau 500    ", 5X,"Temperature",8X,"Ne",9X,"V",14X,"Vturb")'
      
      FOR k=0, atmos.Ndep-1 DO $
       printf, unit, FORMAT='(5E15.6)', $
               alog10(atmos.tau500[k]), atmos.T[k], $
               atmos.n_elec[k] * CM_TO_M^3 , $
               atmos.v[k]/KM_TO_M, atmos.vturb[k]/KM_TO_M

    END
    'GEOMETRIC_SCALE': BEGIN
      printf, unit, $
       FORMAT='("*height [km]   ", 5X,"Temperature",8X,"Ne",9X,"V",14X,"Vturb")'

      FOR k=0, atmos.Ndep-1 DO $
       printf, unit, FORMAT='(E17.8, 4E15.6)', $
               atmos.height[k], atmos.T[k], atmos.n_elec[k] * CM_TO_M^3 , $
               atmos.v[k]/KM_TO_M, atmos.vturb[k]/KM_TO_M
    END
  ENDCASE

  free_lun, unit
END
