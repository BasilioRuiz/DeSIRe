c write_atmos_RH write the model atmosphere in RH format

c Calling a C Function from Fortran
c To call a C function from a Fortran program, ensure that the C function's 
c name is spelled the way the Fortran compiler expects it to be. When you control 
c the name of the C function, the simplest solution is to give it a name that 
c consists of lowercase letters with a terminal underscore. For example, the following C function:
c int fromfort_() {...}
c could be declared in a Fortran program as follows:
c external FROMFORT

        subroutine write_atmos_RH(ID_model,RH_model,RH_mag,atmosLG,ntau)
        implicit real*4 (a-h,o-z)
	include 'PARAMETER'   !por kt,kn,kl,kld
	
	character RH_model*(*),RH_mag*(*),ID_model*(*)
	real*4 atmosLG(*),T,elec_dens,vz,vmic !,logtau
        integer ntau,HEon
c	real*8 degrad,Bs,gs,fis
c	real*8 vector(3*ntau),vector2(3*ntau)
c	real*8 degrad,x,epsil
	real*8 log_mass
	real*8 mass_column(ntau)      !,mass_columnold(ntau)
	real*8 gravity,dlog10g
	real*4 pg1(kt),z1(kt),ro1(kt),pg2(kt),z2(kt),ro2(kt),z1i,z0,zn
	real*4 gamma_rad(kt),phi_rad(kt)
	real*4 deglat,deglon
	common/HEonoff/HEon  !HEon=1 HE, HEon=0 Non HE
	common/zetas/pg1,z1,ro1,pg2,z2,ro2    !para el calculo de z,pg y ro en amp22 y amp22err
	common/latlong/deglat,deglon

c        print*,'write_atmos_RH campo=',ntau,atmosLG(4*ntau+1),atmosLG(6*ntau+1),atmosLG(7*ntau+1)
c        stop
c	epsil=1.d-6
	bolr=1.3806488e-16
c	degrad=3.1415926535897932385d0/180.d0
c       degrad=0.017453293005625408d0
	gravity=27413.847972d0      !cgs
	dlog10g=dlog10(gravity)    !4.43797
	ican=53
	epsilon=1.e-4

	open(ican,file=RH_model)
	write(ican,*)ID_model(1:20)

        write(ican,*)'Mass scale'
        write(ican,'("* log g [cm s^-2]")')
	write(ican,*)dlog10g            !4.43797  !log g
	write(ican,'("* Ndep")')
	write(ican,*)ntau
	write(ican,'("*")')	
        write(ican,'("*lgcolumnMass Temperature Ne V Vturb")')
        
c        print*,'write_atmos_RH ',dlog10g,ntau,HEon,10**atmosLG(ntau),10**atmosLG(1)
c         print*,'write_atmos_RH ',HEon
c         stop
        
        if(HEon .eq.0)then
           print*,'We are running out of Hydrostatic equilibrium'
           taun=10**atmosLG(ntau)
           rhon=atmosLG(10*ntau+2+ntau)
           zn=atmosLG(8*ntau+2+ntau)
	   delta_tau=(10**atmosLG(ntau-1)-taun)       !chnaging sign (-d(tau))
           delta_z=(zn-atmosLG(8*ntau+2+ntau-1))*1.e5
	   romean=(rhon+atmosLG(10*ntau+2+ntau-1))/2.
	   absorp=delta_tau/delta_z/romean
c	   print*,'write_atmos_RH  delta_z=',delta_z,'romean=',romean,'absorp=',absorp,'T=',atmosLG(2*ntau),'pe=',atmosLG(3*ntau)
           mass_column(ntau)=taun/absorp
c           mass_columnold(ntau)=mass_column(ntau)
           rho0=rhon
           z0=zn
           log_mass= dlog10(mass_column(ntau))
           T=atmosLG(2*ntau)
           elec_dens=atmosLG(3*ntau)/bolr/T
           vz=-atmosLG(6*ntau)*1.e-5
           vmic=atmosLG(4*ntau)*1.e-5
           write(ican,100)log_mass,T,elec_dens,vz,vmic      !WARNING Vz=0 --> vz*0.
	   do k=ntau-1,1,-1
	      T=atmosLG(k+ntau)
              elec_dens=atmosLG(k+2*ntau)/bolr/T
              vz=-atmosLG(k+5*ntau)*1.e-5
              vmic=atmosLG(k+3*ntau)*1.e-5
	      rho1=atmosLG(10*ntau+2+k)
	      z1i=atmosLG(8*ntau+2+k)
	      delta_z=(z1i-z0)*1.e5
	      x=rho1/rho0-1.0
	      if(x .gt. epsil)then
	         mass_column(k)=mass_column(k+1)-rho0*delta_z*x/dlog(1.d0+x)
	      else
	         mass_column(k)=mass_column(k+1)-rho0*delta_z/(1.d0-0.5d0*x)
	      endif
              rho0=rho1 
              z0=z1i
c              mass_columnold(i)=mass_column(i)
              log_mass= dlog10(mass_column(k))
              write(ican,100)log_mass,T,elec_dens,vz,vmic      !WARNING Vz=0 --> vz*0.
c              if(k .eq. 301)print*,'stop 1 en write_atmos_RH',log_mass,T,elec_dens,vz,vmic
c              if(k .eq. 301)stop
	  end do
       else
c          print*,'We are in Hydrostatic equilibrium'
	  do k=ntau,1,-1
             T=atmosLG(k+ntau)
             elec_dens=atmosLG(k+2*ntau)/bolr/T
             vz=-atmosLG(k+5*ntau)*1.e-5
             vmic=atmosLG(k+3*ntau)*1.e-5
c             print*,'write_atmos_RH pe',k,atmosLG(k+2*ntau)
             log_mass= dlog10(atmosLG(9*ntau+2+k)*1.d0)-dlog10g
c         print*,'write_atmos_RH log_mass=',atmosLG(9*ntau+2+k)*1.d0,'T=',atmosLG(k+ntau),'pe=',atmosLG(k+2*ntau)
c             print*,'write_atmos_RH 2',log_mass,T,elec_dens,vz,vmic
             write(ican,100)log_mass,T,elec_dens,vz,vmic      !WARNING Vz=0 --> vz*0.
c             if(k .eq. 300)print*,' 300 en write_atmos_RH',log_mass,T,elec_dens,vz,vmic
c             if(k .eq. 301)print*,' 301 en write_atmos_RH',log_mass,T,elec_dens,vz,vmic
c             if(k .eq. 301)stop
          end do
       end if   
c       stop
c       print*,'write_atmos_RH 113 tau=36',atmosLG(36+ntau)
       
        close(ican)

        do i=1,ntau
          gamma_rad(i)=atmosLG(6*ntau+i)
          phi_rad(i)=atmosLG(7*ntau+i)
        end do
        if(abs(deglat) .gt. epsilon .or. abs(deglon) .gt. epsilon)then
           call taulinea5sub(1,deglat,deglon,gamma_rad,phi_rad,ntau) !cambiamos de LoS a local (Z)
        endif
        
        open(ican,file=RH_mag)
c        print*,'Voy a escribir en ',RH_mag,' en write_atmos_RH el campo '
        do i=ntau,1,-1
c            print*,'write_atmos_RH campo',i,atmosLG(4*ntau+i),atmosLG(6*ntau+i), atmosLG(7*ntau+i)
c            write(ican,*) atmosLG(4*ntau+i)*1.d-4,atmosLG(6*ntau+i), atmosLG(7*ntau+i)
             write(ican,*)atmosLG(4*ntau+i)*1.d-4,gamma_rad(i),phi_rad(i)  
        end do
        close(ican)
        
        
100     format(1x,f13.9,1x,f10.2,1x,1pe12.5,1x,e11.4,1x,e11.4,1x)
	return	

	end
	

