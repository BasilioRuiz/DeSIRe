c taulinea3
c il=0 :pasa de una atmosfera modelo en direccion z a una
c       atmosfera en la linea de vision
c il=1 :transforma una atmosfera en la linea de vision a una atmosfera
c       en direccion z
c ntau : numero de puntos en tau
c tau  : log10(tau)
	subroutine taulinea4(il,imodel2,deglon,deglat,vx,vy,atmos,z1,pg1,ro1,z2,pg2,ro2,ntau)
	include 'PARAMETER'   !para kt
	parameter (kt8=8*kt+2)
	implicit real*4 (a-h,o-z)
	integer ihe,il,ntau,imodel2
	real*4 atmos(*),atmos1(kt8),atmos2(kt8),cth,vx,vy,deglon,deglat
        real*4 z1(*),pg1(*),ro1(*),z2(*),pg2(*),ro2(*)
        
        common/truecth/cth !para pasarselo a RH
	pi=3.14159265
	cth=cos(deglon*pi/180.)*cos(deglat*pi/180.)
	
c        print*,'taulinea4 19 V gamma fi=', atmos(5*ntau+10),atmos(6*ntau+10),atmos(7*ntau+10)	    
	do i=1,8*ntau+2
           atmos1(i)=atmos(i)
	end do
	if(imodel2 .ne. 0)then
	   do i=1,8*ntau+2
	      atmos2(i)=atmos(i+8*ntau+2)
c	      print*,'taulinea3 i atmos1 atmos2',ntau,i,atmos1(i),atmos2(i) 
	   end do
	endif   
        
	if(cth.ne.1. .and. il.eq. 0)then
           print*,'We suppose input/output models are radial: transforming to LOS'
	end if
	if(cth.ne.1. .and. il.eq. 1)then
           print*,'Transforming back to LOCAL reference frame'
	end if

c	print*,'estoy en taulinea4  37',il,cth,deglon,deglat,vx,atmos1(1),z1(1),pg1(1),ro1(1),ntau
c        print*,'taulinea4 38 V gamma fi=', atmos1(5*ntau+10),atmos1(6*ntau+10),atmos1(7*ntau+10)
c        print*,'taulinea4 39 cth=',cth,'il= ',il
c        if(cth.gt.0.9999998)cth=0.99
c        print*,'taulinea4 41 cth=',cth,'il= ',il
        call taulinea4sub(il,cth,deglon,deglat,vx,vy,atmos1,z1,pg1,ro1,ntau)
c        print*,'taulinea4 39 V gamma fi=', atmos1(5*ntau+10),atmos1(6*ntau+10),atmos1(7*ntau+10)	

        if(imodel2 .ne. 0)call taulinea4sub(il,cth,deglon,deglat,vx,vy,atmos2,z2,pg2,ro2,ntau)

	do i=1,8*ntau+2
           atmos(i)=atmos1(i)
	end do
c	print*,'taulinea4 47 V gamma fi=',il, atmos(5*ntau+10),atmos(6*ntau+10),atmos(7*ntau+10)	
	
	if(imodel2 .ne. 0)then
	   do i=1,8*ntau+2
              atmos(i+8*ntau+2)=atmos2(i)
	   end do
	endif   

	return
	end

