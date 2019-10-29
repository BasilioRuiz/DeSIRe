        subroutine marqcoeferr(y,ndata,a,m,mfit,beta)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER' !por kn,kl y kld
        parameter (mmax=16*kn+2,kld4=4*kld,kldn=mmax*kld4)        
        dimension y(ndata),beta(*)
        dimension dyda(mmax),a(*),ymod(kld4),ymodobs(kld4),dvcal(kldn)
        dimension ymodobs_RH(kld4)
        real*4 sigreal(kld4)
	integer m(*),npos(kld4)
	
	common/uvesalida/ymodobs,ymodobs_RH
        common/derivchi/dvcal  ! deriv. chi^2 para los errores
        common/primera2/ntotal,ntotal4,ists
        common/posiciones/npos 
	common/sigrealchi/sigreal,chireal,sumsq

      do 12 j=1,mfit
         beta(j)=0.
12    continue

c      print*,'marqcoeferr 23 a(1)=',a(1)
      call fperfil2err(m,a,ymod,dvcal)
      do i=1,ndata
c         ymodobs(i)=ymod(npos(i))
         ymodobs(i)=ymodobs_RH(npos(i))
      end do

      sig0=0.	!Diferencia cuadratica entre perfiles obs. y cal.

      do i=1,ndata
         dy=(y(i)-ymodobs(i))/sigreal(i)
         sig0=sig0+dy*dy
      enddo

      do 15 i=1,ndata

	 do j=1,mfit
	    jj=(j-1)*ntotal4+npos(i)
	    dyda(j)=dvcal(jj)
	 end do

         do 14 j=1,mfit
            ww=dyda(j)/sigreal(i)
            beta(j)=beta(j)+ww*ww
14       continue

15    continue

      do 16 j=1,mfit
         beta(j)=beta(j)/sig0
16    continue


	
      return
      end
