      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)

      implicit real*4 (a-h,o-z)
      include 'PARAMETER'  !mfitmax

      real*4 U(mfitmax,mfitmax),W(mfitmax),V(mfitmax,mfitmax),
     &       B(mfitmax),X(mfitmax),TMP(mfitmax)


      DO 12 J=1,N

        S=0.

c        IF(W(J).NE.0.)THEN
        IF(abs(W(J)).GT.1.e-18)THEN

          DO 11 I=1,M

            S=S+U(I,J)*B(I)
c            print*,'svbksb 17',i,j,U(I,J),B(i),w(j)

11        CONTINUE

          S=S/W(J)

        ENDIF

        TMP(J)=S

12    CONTINUE

      DO 14 J=1,N

        S=0.

        DO 13 JJ=1,N

          S=S+V(J,JJ)*TMP(JJ)
c          print*,'svbksb 36',j,jj,V(J,JJ),TMP(JJ)

13      CONTINUE

        X(J)=S

14    CONTINUE

      RETURN

      END

