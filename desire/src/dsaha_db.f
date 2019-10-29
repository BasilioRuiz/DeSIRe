      FUNCTION dSAHA_db(THETA,CHI,dU1,dU2)
c	calcula la detivada respecto a t, del logaritmo de la
C     SAHA-EGGERT EQUATION
c	la derivada del logaritmo de saha respecto a pe es -1/pe
      real*4 THETA,CHI,dU1,dU2
      real*8 dSAHA_DB


      dSAHA_db=1.d0*(du2-du1)+(theta/5040.d0)*(2.5d0+chi*THETA*dlog(10.d0))
      
      RETURN
      END   
