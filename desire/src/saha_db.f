      FUNCTION SAHA_DB(THETA,CHI,U1,U2,PE) 
C     SAHA-EGGERT EQUATION  
      real*4 THETA,CHI,U1,U2,PE
      real*8 SAHA_DB,SAHA1,theta25
      
      theta25=(theta*1.d0)**2.5d0
      SAHA1=U2*10.d0**(9.0805126d0-THETA*CHI*1.d0)/(U1*PE*theta25) 

      if(SAHA1 .lt. 1.d38)then
        SAHA_DB=SAHA1
      else
        SAHA_DB=1.d38
      endif  

      RETURN
      END   


