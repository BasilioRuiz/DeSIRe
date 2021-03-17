c taulinea4sub
c il=0 :pasa de una atmosfera modelo en direccion z a una
c       atmosfera en la linea de vision
c il=1 :transforma una atmosfera en la linea de vision a una atmosfera
c       en direccion z
c ntau : numero de puntos en tau
c tau  : log10(tau)
c vx   : Omega(equator)*R= tangential velocity at equator (in Km/s)

        subroutine taulinea4sub(il,cth,deglon,deglat,vx,vy,atmos,z,pg,ro,ntau)
        implicit real*4 (a-h,o-z)
        include 'PARAMETER'  !por kt
        real*4 atmos(*),sqrtb,deglon,deglat,logcth,vx,vy,cth
        real*4 Bx,By,Bz,Bx2,By2,Bz2,z(*),pg(*),ro(*)
        real*4 tau(kt),t(kt),p(kt),vmic(kt),h(kt),vz(kt),gamma(kt),phi(kt)
        real*4 ttau(kt)
        integer il,ntau

        pi=3.14159265
        pi180=0.01745329250000000
        epsilon=1.e-4
        logcth=alog10(abs(cth))
c       print*,'logcth longitude latitude=',logcth,deglon,deglat
        
c changing gamma a phi to radian or degrees
        if(il.eq.0)then
            do i=1,ntau
	        atmos(i+6*ntau)=atmos(i+6*ntau)*pi180 !input gedrees radians (Local reference system)
	        atmos(i+7*ntau)=atmos(i+7*ntau)*pi180 
            enddo
        endif
        if(abs(cth-1.).le.epsilon)then
           if(il.eq.1)then
              do i=1,ntau
	             atmos(i+6*ntau)=atmos(i+6*ntau)/pi180 !input gedrees radians (Local reference system)
	             atmos(i+7*ntau)=atmos(i+7*ntau)/pi180 
              enddo
           endif
           return
        end if   

       logcth=alog10(abs(cth))     
	do i=1,ntau
	   ttau(i)=atmos(i)
	   t(i)=atmos(i+ntau)
	   p(i)=atmos(i+2*ntau)
	   vmic(i)=atmos(i+3*ntau)
	   h(i)=atmos(i+4*ntau)
	   vz(i)=atmos(i+5*ntau)
	   gamma(i)=atmos(i+6*ntau)
	   phi(i)=atmos(i+7*ntau)
	end do
		
	clon=cos(deglon*pi180)
	slon=sin(deglon*pi180)
	clat=cos(deglat*pi180)
	slat=sin(deglat*pi180)	
	      
	if(il.eq.0)then

	   do i=1,ntau
	      tau(i)=ttau(i)-logcth
	      cgaz=cos(gamma(i))
	      sgaz=sin(gamma(i))
	      cfiz=cos(phi(i))
	      sfiz=sin(phi(i))
	      
	      Bx=sgaz*cfiz
	      By=sgaz*sfiz
	      Bz=cgaz
	      
c turn around X axis latitude angle	      
	      Bx1= Bx
	      By1= By*clat+Bz*slat
	      Bz1=-By*slat+Bz*clat
	      
c	      vx1= vx
c	      vy1= vy*clat+vzi*slat
c	      vz1=-vy*slat+vzi*clat
	      
c turn around Y axis longitude angle	      	      
	      Bx2= Bx1*clon+Bz1*slon
	      By2= By1
	      Bz2=-Bx1*slon+Bz1*clon
	      
c	      vx2= vx1*clon+vz1*slon
c	      vy2= vy1
c	      vz2=-vx1*slon+vz1*clon
	      
c	      vz(i)=vx*clat*slon+vz(i)*cth !changing sign vz>0 downward vx>0 positive X, right of the disk  
c              vz(i)=-vz2   !changing sign vz>0 downward vx>0 positive X, right of the disk 

c evaluating new gamma and phi
              cga=Bz2      
c	      cga=-sth*cfiz*sgaz+cth*cgaz
	      if(cga.gt.1.)then
		     cga=1.-epsilon/100.
	      endif
	      if(cga.lt.-1.)then
		     cga=-1.+epsilon/100.
	      endif
	      	      
	      sga=sqrtb(1.e0-cga*cga) !always positive

c	      if(sgaz.lt.0)sga=-sga
	      if(abs(sga).le.epsilon/10000000.)then
	         phi(i)=0.             !atencion, cambio para HORST
	      else
		 gamma(i)=acos(cga)
	         cfi=Bx2/sga
	         sfi=By2/sga
	         phi(i)=atan2(sfi,cfi)
              end if

              ggg=gamma(i)
	      if(ggg.lt.-2.*pi)ggg=ggg-2.*pi*int(ggg/(2.*pi))
	      if(ggg.gt.-2.*pi.and.ggg.lt.-pi)ggg=2.*pi+ggg
	      if(ggg.lt.0.and.ggg.gt.-pi)ggg=-ggg
	      if(ggg.ge.2.*pi)ggg=ggg-2.*pi*int(ggg/(2.*pi))
	      if(ggg.gt.pi)ggg=2.*pi-ggg
              gamma(i)=ggg
	    end do

c	    print*,'taulinea4sub line 117: going out after transform to LoS', gamma(1)
c	    call quitasaltos2(ntau,pi,phi)
	    call quitasaltos(ntau,pi,phi)

	   else
            do i=1,ntau
	       tau(i)=ttau(i)+logcth
c	       vz(i)=(vz(i)-vx*sth)/cth
c               vzi=vz(i)               !keep the sign because I'm resolving vz
c               if(i.eq.1)print*,'taulinea4sub 127',gamma(i),cos(gamma(i)),sin(gamma(i))

	       cga=cos(gamma(i))
	       sga=sin(gamma(i))
	       cfi=cos(phi(i))
	       sfi=sin(phi(i))
c	       cgaz=sth*cfi*sga+cth*cga
	     	      
	       Bx=sga*cfi
	       By=sga*sfi
	       Bz=cga
	      
c turn around Y axis -longitude angle	      	      
	       Bx1= Bx*clon-Bz*slon
	       By1= By
	       Bz1= Bx*slon+Bz*clon	      
	      
c turn around X axis -latitude angle	      
	       Bx2= Bx1
	       By2= By1*clat-Bz1*slat
	       Bz2= By1*slat+Bz1*clat
c	       vz(i)=(vz(i)-vx*clat*slon)/cth
c	       vz(i)=(vz(i)-vx*slon-vy*slat*clon)/cth	       
c evaluating new gamma and phi
           cgaz=Bz2      

	       if(cgaz.gt.1.)cgaz=1.-epsilon/100.
	       if(cgaz.lt.-1.)cgaz=-1.+epsilon/100.
	       sgaz=sqrt(1.e0-cgaz*cgaz)	       
           gamma(i)=acos(cgaz)/pi180    !in degrees (LOS reference system)

	       if(abs(sgaz).le.epsilon/1000000.)then
	          phi(i)=atmos(i+7*ntau)   !fi !atencion, cambio para HORST
	       else
              cfiz=Bx2/sgaz
	          sfiz=By2/sgaz                
	          phi(i)=atan2(sfiz,cfiz)/pi180 !in degrees (LOS reference system)
           end if

           ggg=gamma(i)
           if(ggg.lt.-360.)ggg=ggg-360.*int(ggg/360.)
           if(ggg.gt.-360..and.ggg.lt.-180.)ggg=360.+ggg
           if(ggg.lt.0.and.ggg.gt.-180.)ggg=-ggg
           if(ggg.ge.360.)ggg=ggg-360.*int(ggg/360.)
           if(ggg.gt.180.)ggg=360.-ggg
           gamma(i)=ggg
        end do
c	    call quitasaltos2(ntau,180.,phi)
        call quitasaltos(ntau,180.,phi)
      end if

c         interpolamos
c	  call intmodel4(3,ntau,ttau,tau,t,p,vmic,h,vz,gamma,phi,z,pg,ro)

	    do i=1,ntau
	       atmos(i)=tau(i)
c	       atmos(i+ntau)=t(i)
c	       atmos(i+2*ntau)=p(i)
c	       atmos(i+3*ntau)=vmic(i)
	       atmos(i+4*ntau)=h(i)
	       atmos(i+5*ntau)=vz(i)
 	       atmos(i+6*ntau)=gamma(i)
	       atmos(i+7*ntau)=phi(i)
	    end do
	    
	return
	end
c ----------------------------------------------------------------------------------
c taulinea5sub
c il=0 :pasa gamma y fi en radianes de atmosfera en Local reference frame to
c       atmosfera en la linea de vision
c il=1 :transforma una atmosfera en la linea de vision a una atmosfera Local RF
c ntau : numero de puntos en tau
c
c 10/10/20 epm: Esta subrutina es llamada por write_atmos_RH() la cual pasa
c los arrays 'gamma_rad' y 'phi_rad' de arriba a abajo, es decir, desde
c ntau hasta 1 ya que es como RH lo espera.

	subroutine taulinea5sub(il,deglat,deglon,gamma_rad,phi_rad,ntau)
        implicit real*4 (a-h,o-z)
	include 'PARAMETER'  !por kt
	real*4 gamma_rad(*),phi_rad(*),sqrtb,deglon,deglat
	real*4 Bx,By,Bz,Bx2,By2,Bz2
	integer il,ntau

	pi=3.14159265
	pi180=pi/180.
	epsilon=1.e-5
        		
	clon=cos(deglon*pi180)
	slon=sin(deglon*pi180)
	clat=cos(deglat*pi180)
	slat=sin(deglat*pi180)	
	      
	if(il.eq.0)then
	   do i=1,ntau
	      cgaz=cos(gamma_rad(i))
	      sgaz=sin(gamma_rad(i))
	      cfiz=cos(phi_rad(i))
	      sfiz=sin(phi_rad(i))
	      
	      Bx=sgaz*cfiz
	      By=sgaz*sfiz
	      Bz=cgaz
	      
c turn around X axis latitude angle	      
	      Bx1= Bx
	      By1= By*clat+Bz*slat
	      Bz1=-By*slat+Bz*clat
	      	      
c turn around Y axis longitude angle	      	      
	      Bx2= Bx1*clon+Bz1*slon
	      By2= By1
	      Bz2=-Bx1*slon+Bz1*clon
	      
c evaluating new gamma and phi
              cga=Bz2      
	      if(cga.gt.1.)then
		  cga=1.-epsilon/100.
	      endif
	      if(cga.lt.-1.)then
		  cga=-1.+epsilon/100.
	      endif
	      	      
	      sga=sqrtb(1.e0-cga*cga) !always positive

	      if(abs(sga).le.epsilon/10000000.)then
	         phi_rad(i)=0.             !atencion, cambio para HORST
	      else
		 gamma_rad(i)=acos(cga)
	         cfi=Bx2/sga
	         sfi=By2/sga
	         phi_rad(i)=atan2(sfi,cfi)
              end if
	    end do

	  else
	  
            do i=1,ntau
	       cga=cos(gamma_rad(i))
	       sga=sin(gamma_rad(i))
	       cfi=cos(phi_rad(i))
	       sfi=sin(phi_rad(i))
	     	      
	       Bx=sga*cfi
	       By=sga*sfi
	       Bz=cga
	      
c turn around Y axis -longitude angle	      	      
	       Bx1= Bx*clon-Bz*slon
	       By1= By
	       Bz1= Bx*slon+Bz*clon	      
	      
c turn around X axis -latitude angle	      
	       Bx2= Bx1
	       By2= By1*clat-Bz1*slat
	       Bz2= By1*slat+Bz1*clat
	      
c evaluating new gamma and phi
               cgaz=Bz2                     
	       if(cgaz.gt.1.)cgaz=1.-epsilon/100.
	       if(cgaz.lt.-1.)cgaz=-1.+epsilon/100.
	       sgaz=sqrt(1.e0-cgaz*cgaz)
	       
               gamma_rad(i)=acos(cgaz)    

	       if(abs(sgaz).gt.epsilon/1000000.)then
                  cfiz=Bx2/sgaz
	          sfiz=By2/sgaz                
	          phi_rad(i)=atan2(sfiz,cfiz)
               end if
	    end do

	  end if


	    
	return
	end



           

