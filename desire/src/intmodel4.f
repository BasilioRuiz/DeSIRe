c esta rutina interpola en tau el modelo ttau,tt,ppe...y da la salida en
c ttau,tt,ppe...por medio de un polinomio de grado ngrado

	subroutine intmodel4(ngrado,n,tau,ttau,tt,ppe,mmic,hh,vvz,gg,ffi,zz,ppg,rrho)
	implicit real*4 (a-h,o-z)

	include 'PARAMETER'  !por kt
	real ttau(*),tt(*),ppe(*),hh(*),MMIC(*),VVZ(*),gg(*),ffi(*)
	real zz(*),ppg(*),rrho(*)
	real tau(*),t(kt),pe(kt),h(kt),vz(kt),MIC(kt),G(kt),FI(kt)
	real z(kt),pg(kt),rho(kt)
	real xa(11),ya(11)

	num=n
c interpolamos
	n2=int(ngrado/2)
	
	do i=1,n
	       CALL LOCATE(TTAU,NUM,TAU(i),J)
	       n3=j-n2-1
               if(n3.lt.0)n3=0
               if(n3+ngrado+1.gt.num)n3=num-ngrado-1
	       do k=1,ngrado+1
		     xa(k)=ttau(n3+k)
	       end do
	       do k=1,ngrado+1
		     ya(k)=tt(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(i),T(i),ERROR)
	       do k=1,ngrado+1
		     ya(k)=alog(ppe(n3+k))
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(i),pe(i),ERROR)
	       
	       do k=1,ngrado+1
		     ya(k)=hh(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(i),h(i),ERROR)
	       do k=1,ngrado+1
	             ya(k)=vvz(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(i),vz(i),ERROR)
	       do k=1,ngrado+1
		     ya(k)=mmic(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(i),mic(i),ERROR)
	       do k=1,ngrado+1
		     ya(k)=gg(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(i),g(i),ERROR)
	       do k=1,ngrado+1
		     ya(k)=ffi(n3+k)
	       end do
               do k=2,ngrado+1
                   dya=ya(k)-ya(k-1)
                   if(abs(dya).ge.350)ya(k)=ya(k)-360*int(dya/350.)
               end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(i),fi(i),ERROR)
	       do k=1,ngrado+1
		  ya(k)=zz(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(i),z(i),ERROR)
	       do k=1,ngrado+1
		  ya(k)=alog(ppg(n3+k))
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(i),pg(i),ERROR)
	       do k=1,ngrado+1
		  ya(k)=alog(rrho(n3+k))
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(i),rho(i),ERROR)
	       
	    end do
            do i=2,n
               dfi=fi(i)-fi(i-1)
               if(abs(dfi).ge.350)fi(i)=fi(i)-360.*int(dfi/360.)
            end do 
        
	    do i=1,n
	        ttau(i)=tau(i)
	        tt(i)=t(i)
		ppe(i)=exp(pe(i))
	        mmic(i)=mic(i)
	        hh(i)=h(i)
	        vvz(i)=vz(i)
	        gg(i)=g(i)
	        ffi(i)=fi(i)
	        zz(i)=z(i)
	        ppg(i)=exp(pg(i))
	        rrho(i)=exp(rho(i))
	    end do

	return
	end
