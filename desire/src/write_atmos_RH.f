c_____________________________________________________________________________
c
c write_atmos_RH
c write the model atmosphere in RH format
c
c 11/11/20 epm: Save the model in memory instead of file.
c_____________________________________________________________________________

        subroutine write_atmos_RH(ID_model,RH_model,RH_mag,atmosLG,ntau)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        character RH_model*(*),RH_mag*(*),ID_model*(*)
        character geoscale
        integer*4 ntau,nB,HEon
        integer*4 codes(100)
        real*4    deglat,deglon
        real*4    atmosLG(*)
        real*4    T(kt),elec_dens(kt),vz(kt),vmic(kt)
        real*4    B(kt),gamma(kt),phi(kt)
        real*8    gravity,dlog10g
        real*8    mass_column(ntau)
        real*8    log_mass(kt)

        common/HEonoff/HEon  !HEon=1 HE, HEon=0 Non HE
        common/latlong/deglat,deglon
c       Usamos este common compartido con write_atmos_RH_tau() para salvar
c       las variables en el (en vez de usar save) y ahorrar algo de memoria.
        common/saveatmos/log_mass,T,elec_dens,vz,vmic,B,gamma,phi

        epsil=1.e-6
        bolr=1.3806488e-16
        gravity=27413.847972d0   !cgs
        dlog10g=dlog10(gravity)  !4.43797

c       ican=53
c       RH_model='Temporary.atmos'
c       open(ican,file=RH_model)
c       write(ican,*)ID_model(1:20)
c       write(ican,*)'Mass scale'
c       write(ican,'("* log g [cm s^-2]")')
c       write(ican,*)dlog10g
c       write(ican,'("* Ndep")')
c       write(ican,*)ntau
c       write(ican,'("*")')
c       write(ican,'("*lgcolumnMass Temperature Ne V Vturb")')

        call toascii(trim(ID_model),100,codes)
        ncodes=len_trim(ID_model)
        geoscale='M'

c       Le damos la vuelta a los arrays porque SIR va desde lo mas profundo
c       a la superficie y RH va al reves, de la superficie a lo profundo.

        i=1
        if(HEon.eq.0)then
           call error(KWARN,'write_atmos_RH',
     &                'We are running out of hydrostatic equilibrium')
           taun=10**atmosLG(ntau)
           rhon=atmosLG(10*ntau+2+ntau)
           zn=atmosLG(8*ntau+2+ntau)
           delta_tau=(10**atmosLG(ntau-1)-taun)  !changing sign (-d(tau))
           delta_z=(zn-atmosLG(8*ntau+2+ntau-1))*1.e5
           romean=(rhon+atmosLG(10*ntau+2+ntau-1))/2.
           absorp=delta_tau/delta_z/romean
           mass_column(ntau)=taun/absorp
           rho0=rhon
           z0=zn
           T(i)=atmosLG(2*ntau)
           elec_dens(i)=atmosLG(3*ntau)/bolr/T(i)
           vz(i)=-atmosLG(6*ntau)*1.e-5
           vmic(i)=atmosLG(4*ntau)*1.e-5
           log_mass(i)= dlog10(mass_column(ntau))
c          write(ican,100)log_mass(i),T(i),elec_dens(i),vz(i),vmic(i)
           i=i+1
           do k=ntau-1,1,-1
              T(i)=atmosLG(k+ntau)
              elec_dens(i)=atmosLG(k+2*ntau)/bolr/T(i)
              vz(i)=-atmosLG(k+5*ntau)*1.e-5
              vmic(i)=atmosLG(k+3*ntau)*1.e-5
              rho1=atmosLG(10*ntau+2+k)
              z1i=atmosLG(8*ntau+2+k)
              delta_z=(z1i-z0)*1.e5
              x=rho1/rho0-1.0
              if(x.gt.epsil)then
                 mass_column(k)=mass_column(k+1)-rho0*delta_z*x/dlog(1.d0+x)
              else
                 mass_column(k)=mass_column(k+1)-rho0*delta_z/(1.d0-0.5d0*x)
              endif
              rho0=rho1
              z0=z1i
              log_mass(i)= dlog10(mass_column(k))
c             write(ican,100)log_mass(i),T(i),elec_dens(i),vz(i),vmic(i)
              i=i+1
           end do

        else
           do k=ntau,1,-1
              T(i)=atmosLG(k+ntau)
              elec_dens(i)=atmosLG(k+2*ntau)/bolr/T(i)
              vz(i)=-atmosLG(k+5*ntau)*1.e-5
              vmic(i)=atmosLG(k+3*ntau)*1.e-5
              log_mass(i)= dlog10(dble(atmosLG(9*ntau+2+k)))-dlog10g
c             write(ican,100)log_mass(i),T(i),elec_dens(i),vz(i),vmic(i)
              i=i+1
           end do
        end if

c       close(ican)

c       Escritura del campo magnetico.

        nB=0
        i=1
        do k=ntau,1,-1
           B(i)=atmosLG(4*ntau+k)*1.e-4
           gamma(i)=atmosLG(6*ntau+k)
           phi(i)=atmosLG(7*ntau+k)
           if(B(i).gt.0) nB=ntau
           i=i+1
        end do

c       RH_mag='Temporary_field.mod'
c       if(RH_mag.ne.'none')then
c          open(ican,file=RH_mag)
c          do i=1,ntau
c             write(ican,*)B(i),gamma(i),phi(i)
c          end do
c          close(ican)
c       end if

c       Pasamos todo a RH.

        call siratmos(codes,ncodes,geoscale,dlog10g,ntau,
     &                log_mass,T,elec_dens,vz,vmic,nB,B,gamma,phi)

100     format(1x,f13.9,1x,f10.2,1x,1pe12.5,1x,e11.4,1x,e11.4,1x)

        return
        end

c_____________________________________________________________________________
c
c write_atmos_RH_tau
c write the model atmosphere in RH format in logtau scale (needs H populations)
c
c 11/11/20 epm: Save the model in memory instead of file.
c_____________________________________________________________________________

        subroutine write_atmos_RH_tau(ID_model,RH_model,RH_mag,atmosLG,
     &                                ntau,natmos)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        character RH_model*(*),RH_mag*(*),ID_model*(*)
        character geoscale
        integer*4 ntau,nB,natmos
        integer*4 codes(100)
        real*4    deglat,deglon
        real*4    atmosLG(*)
        real*4    T(kt),elec_dens(kt),vz(kt),vmic(kt)
        real*4    B(kt),gamma(kt),phi(kt)
        real*4    hydro1(6,kt),hydro2(6,kt)
        real*4    nh1(kt),nh2(kt),nh3(kt),nh4(kt),nh5(kt),nh6(kt)
        real*8    gravity,dlog10g
        real*8    log_tau(kt)

        character*40 msg

        save nh1,nh2,nh3,nh4,nh5,nh6

        common/hydroge_populations/hydro1,hydro2  !float (6,kt), H populations
        common/latlong/deglat,deglon
c       Usamos este common compartido con write_atmos_RH() para salvar
c       las variables en el (en vez de usar save) y ahorrar algo de memoria.
        common/saveatmos/log_tau,T,elec_dens,vz,vmic,B,gamma,phi

        epsil=1.e-6
        bolr=1.3806488e-16
        gravity=27413.847972d0   !cgs
        dlog10g=dlog10(gravity)  !4.43797

c       ican=53
c       RH_model='Temporary.atmos'
c       open(ican,file=RH_model)
c       write(ican,*)ID_model(1:20)
c       write(ican,*)'Tau scale'
c       write(ican,'("* log g [cm s^-2]")')
c       write(ican,*)dlog10g
c       write(ican,'("* Ndep")')
c       write(ican,*)ntau
c       write(ican,'("*")')
c       write(ican,'("*lgtau Temperature Ne V Vturb")')

        call toascii(trim(ID_model),100,codes)
        ncodes=len_trim(ID_model)
        geoscale='T'

c       Le damos la vuelta a los arrays porque SIR va desde lo mas profundo
c       a la superficie y RH va al reves, de la superficie a lo profundo.

        i=1
        do k=ntau,1,-1
           T(i)=atmosLG(k+ntau)
           elec_dens(i)=atmosLG(k+2*ntau)/bolr/T(i)
           vz(i)=-atmosLG(k+5*ntau)*1.e-5
           vmic(i)=atmosLG(k+3*ntau)*1.e-5
           log_tau(i)= dble(atmosLG(k))
c          write(ican,100)log_tau(i),T(i),elec_dens(i),vz(i),vmic(i)
           i=i+1
        end do

c       write(ican,'("*")')
c       write(ican,'("* Hydrogen populations ")')
c       write(ican,'("*    nh(1)        nh(2)        n(h3)        nh(4)'
c    &  //         '        nh(5)        np ")')
        if(natmos.eq.1)then
           i=1
           do k=ntau,1,-1
              nh1(i)=hydro1(1,k)
              nh2(i)=hydro1(2,k)
              nh3(i)=hydro1(3,k)
              nh4(i)=hydro1(4,k)
              nh5(i)=hydro1(5,k)
              nh6(i)=hydro1(6,k)
              i=i+1
c             write(ican,101)(hydro1(j,k),j=1,6)
           end do
        else if(natmos.eq.2)then
           i=1
           do k=ntau,1,-1
              nh1(i)=hydro2(1,k)
              nh2(i)=hydro2(2,k)
              nh3(i)=hydro2(3,k)
              nh4(i)=hydro2(4,k)
              nh5(i)=hydro2(5,k)
              nh6(i)=hydro2(6,k)
              i=i+1
c             write(ican,101)(hydro2(j,k),j=1,6)
           end do
        else
           write(msg,'(a,i2)') 'Wrong number of atmospheres = ',natmos
           call error(KSTOP,'write_atmos_RH_tau',msg)
        end if

c       close(ican)

c       Escritura del campo magnetico.

        nB=0
        i=1
        do k=ntau,1,-1
           B(i)=atmosLG(4*ntau+k)*1.e-4
           gamma(i)=atmosLG(6*ntau+k)
           phi(i)=atmosLG(7*ntau+k)
           if(B(i).gt.0) nB=ntau
           i=i+1
        end do

c       RH_mag  ='Temporary_field.mod'
c       if(RH_mag.ne.'none')then
c          open(ican,file=RH_mag)
c          do i=1,ntau
c             write(ican,*)B(i),gamma(i),phi(i)
c          end do
c          close(ican)
c       end if

c       Pasamos todo a RH.
        call siratmos(codes,ncodes,geoscale,dlog10g,ntau,
     &                log_tau,T,elec_dens,vz,vmic,nB,B,gamma,phi)
        call siratmosnh(nh1,nh2,nh3,nh4,nh5,nh6)

100     format(1x,f13.9,1x,f10.2,1x,1pe12.5,1x,e11.4,1x,e11.4,1x)
101     format(6(1x,e12.5))

        return
        end
