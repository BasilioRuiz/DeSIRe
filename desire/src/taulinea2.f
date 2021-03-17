c taulinea2
c il=0 :pasa de una atmosfera modelo en direccion z a una
c       atmosfera en la linea de vision
c il=1 :transforma una atmosfera en la linea de vision a una atmosfera
c       en direccion z
c ntau : numero de puntos en tau
c tau  : log10(tau)
        subroutine taulinea2(il,cth,ihe,vx,atmos,ntau)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        parameter (kt8=8*kt+2)
        integer ihe,il,ntau
        real*4 atmos(*),atmos1(kt8),atmos2(kt8),cth,vx
        character*40 msg
c       real*4 gammal1(kt),gammal2(kt),phil1(kt),phil2(kt)
c       common/LOS/ gammal1,phil1,gammal2,phil2  !para ASP: subsir

        do i=1,8*ntau+2
           atmos1(i)=atmos(i)
           atmos2(i)=atmos(i+8*ntau+2)
        end do

c	if(cth.ne.1.)then
	if(abs(cth-1) .lt. 1.e-4)then 
           write(msg,*)cth
           call error(KWARN,'taulinea2','The atmospheres are assumed to be'
     &     //         ' line of sight.\n However, hydrostatic equilibrium is'
     &     //         ' computed along the vertical\n using the heliocentric'
     &     //         ' angle provided (cth = '//trim(msg)//')')
	end if

c	call taulinea(il,cth,ihe,vx,atmos1,ntau,gammal1,phil1)  ! gammal1???
c	call taulinea(il,cth,ihe,vx,atmos2,ntau,gammal2,phil2)

        call taulinea(il,cth,ihe,vx,atmos1,ntau)
        call taulinea(il,cth,ihe,vx,atmos2,ntau)

	do i=1,8*ntau+2
           atmos(i)=atmos1(i)
           atmos(i+8*ntau+2)=atmos2(i)
	end do

	return
	end

