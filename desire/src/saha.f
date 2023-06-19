      FUNCTION SAHA(THETA,CHI,U1,U2,PE) 
C     SAHA-EGGERT EQUATION  
      implicit real*4 (a-h,o-z)
      SAHA1=U2*10.**(9.0805126-THETA*CHI)/(U1*PE*THETA**2.5) 
      if(SAHA1 .lt. 1.e38)then
        SAHA=SAHA1
      else
        SAHA=1.e38
      endif  
      RETURN
      END   

      FUNCTION SAHAdouble(THETA,CHI,U1,U2,PE) 
      	implicit real*4 (a-h,o-z)
      real*8 THETA,PE,SAHAdouble
      real*4 CHI,U1,U2
      real*8 CHId,U1d,U2d
      
      CHId=dble(CHI)
      U1d=dble(U1)
      U2d=dble(U2)

      SAHAdouble=U2d*10.d0**(9.0805126d0-THETA*CHId)/(U1d*PE*THETA**2.5d0) 

      RETURN
      END   

      FUNCTION SAHAdoubleinv(THETA,CHI,U1,U2,PE) 
      	implicit real*4 (a-h,o-z)
      real*8 THETA,PE,SAHAdoubleinv
      real*4 CHI,U1,U2
      real*8 CHId,U1d,U2d
      
      CHId=dble(CHI)
      U1d=dble(U1)
      U2d=dble(U2)
      
      SAHAdoubleinv=8.3078262d-10*10.d0**(THETA*CHId)*THETA**2.5d0*PE*U1d/U2d

      RETURN
      END   
