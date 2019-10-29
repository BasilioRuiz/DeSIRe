c leepesos lee un fichero (vobs) de datos de pesos para los parametros de stokes, pesos
c el formato es el mismo que el de los perfiles observados
	
	subroutine leepesos(vobs,ist,pesos)
        implicit real*4 (a-h,o-z)

	include 'PARAMETER'  !por kld
	real*4 pesos(*),dl,si(kld),sq(kld),su(kld),sv(kld)
	integer ist(4)
	character vobs*(*)
	character*100 control
	common/canal/icanal
	common/nombrecontrol/control

	ican=56
	open(ican,file=vobs,status='old',err=991)
	print*,'Reading weights from file: ',trim(vobs)
	k=0
	n0=-1
        i=0
	do while(k.lt.kld)
	   k=k+1
	   read(ican,*,err=992,end=990)a,dl,si(k),sq(k),su(k),sv(k)
	end do
990     ntl=i
        nli=k-1
	call sfromiquv(ist,nli,si,sq,su,sv,pesos)  !cuelga de leeuve2
	close(ican)
	return

991	print*,' '
        print*,'The file containing the weights does NOT exist. We do not introduce weights!!!'
        print*,trim(vobs)
        print*,'______________________________________________________________________________'
c	open(icanal,file=control,fileopt='eof')
c        write(icanal,*) 'STOP: The file containing the observed or the stray light profiles does NOT exist:'
c        write(icanal,*) trim(vobs)
c	close(icanal)
c	stop
        return

992	print*,' '
	print*,'STOP: Incorrect format in the file containing the weights:'
        print*,trim(vobs)
        print*,'______________________________________________________________________________'
c	open(icanal,file=control,fileopt='eof')
c	write(icanal,*) 'STOP: Incorrect format in the file containing the observed or stray light profiles'
c        write(icanal,*) trim(vobs)
c	close(icanal)

	stop


	end
