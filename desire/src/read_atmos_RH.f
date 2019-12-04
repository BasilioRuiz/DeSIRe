c
c read_atmos_RH reads the model atmosphere in RH (Multi) format
c and runs the program "convert" to reevaluate the Hydrostatic Equilibrium
c using RH subroutines
c natmos=1: First component, natmos=2: second component
c
c 19/06/19 epm: New routines to run RH executables as functions.
c ________________________________________________________________________

      subroutine read_atmos_RH(model_multi,natmos,atmos,B,gamma,phi,z,pg,rho,
     +ntau,Vmac,fill,strayfac)

      implicit real*4 (a-h,o-z)
      include 'PARAMETER'   !por kt,kn,kl,kld

      character model_multi*(*)
      real*4 atmos(*)
      real*4 Bi
      real*4 pg(*),z(*),rho(*)
      real*4 B(*),gamma(*),phi(*)
      real*4 Vmac,fill,strayfac
      real*4 rhomedio(kt),rho1(kt)      !,pg1(kt)
c     19/06/19 epm: we don't need 'tau(kt)' anymore.
c     real*4 tau(kt)
c     real*4 wgt,abu,ei1,ei2,pmusum,asum
c     real*4 cgases,avog
      real*4 bol,gravity
      integer ntau,natmos,nn,HEon
      integer ndepths,istatus
      character*100 ruta,ruta2,ruta3,filewavename
      character*80 men1,men2,men3
      real*4 alog10mu
      real*4 tsi,psi,pgsi,pesomediosi,rhosi,kappasi

c     19/06/19 epm: Variables to be read from a file written by 'conversion'.
c     real*4 tau_RH(kt), mass(kt), z1(kt), T(kt), nele(kt), vz(kt), Vturb(kt)
c     19/06/19 epm: Double precision for the variables of conversion().
      real*8 tau_RH(kt), mass(kt), z1(kt), T(kt), nele(kt), vz(kt), Vturb(kt)

c     19/06/19 epm: New functions and variables to run RH applications.
      integer*4 conversion      !, system
      integer*4 atmos_nspace
      integer*8 cpu1, cpu2, wall1, wall2

      common/nohaycampo/nohaycampo1,nohaycampo2
      common/ruta/ruta,ruta2,ruta3,filewavename
      common/contornopg/ncontpg,pg01,pg02
      common/contornoro/ro01,ro02
      common/alog10mu/alog10mu             !alog10mu=alog10(xmu) from desire                            
      common/HEonoff/HEon  !HEon=1 HE, HEon=0 Non HE

c     19/06/19 epm: New functions to run RH applications.
      write(0, *) ''
      write(0, *) '$ CALLING conversion from "librh.a"'
      write(0, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      call chrono(cpu1, wall1)
c     istatus=system(ruta3)
      istatus = conversion(atmos_nspace, kt, tau_RH, mass, z1, T, nele, vz,
     &                     Vturb)
      call chrono(cpu2, wall2)
      write(0, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(0, '(1x,a,8x,i8)') 'conversion: status = ', istatus
      write(0, '(1x,a,i8)') 'conversion: CPU  time (ms) = ', cpu2 - cpu1
      write(0, '(1x,a,i8)') 'conversion: Wall time (ms) = ', wall2 - wall1
      write(0, *) ''
      if (istatus .ne. 0) then
         print*,'*** ERROR #1 in read_atmos_RH.f: '
         print*,'Arrays for the converted model with insufficient memory.'
         print'(1x,a,i6,a,$)','Allocated ', kt, ' bytes but '
         print'(i6,a)',atmos_nspace, ' are needed.'
         stop ' '
      end if
      ndepths = atmos_nspace

c     19/06/19 epm: We don't need this anymore.
c     We have writen in "atmos.txt.out":
c     ndepths
c     tau(i)  mass(i)  z(i)  T(i)  ne(i)  vz(i)  Vturb(i)
c     open(50,file='atmos.txt.out')
c     read(50,*)ndepths
c     do i=1,ndepths
c        read(50,*)tau_RH(i),mass(i),z1(i),T(i),nele(i),vz(i),Vturb(i) !MKSA
c        print*,'leyendo conversion ',i,tau_RH(i),mass(i),z1(i),T(i),nele(i),vz(i),Vturb(i)
c        tau(i)=alog10(tau_RH(i))
c      end do 
c      print*,'read_atmos ntau=',ndepths,'   tau(1)=',tau(1)
c      close(50)

      bol=1.3806488e-16         !erg/K
      gravity=2.7414e+4         !cm/s^2
c     cgases=83145100.
c     avog=6.023e23
c     print*,bol,gravity,ndepths
      ntau=ndepths
      if (ntau.gt.kt) then
         men1='STOP: The model has more depth points than the current limit kt.'
         men2='      Decrease the number of grid points or change the PARAMETER file.'
         men3=' '
         call mensaje(2,men1,men2,men3)
         stop
      end if

      nn=0
      if(natmos .eq. 2)nn=8*ntau+2

      do i=2,ntau
         deltaz=real((z1(i-1)-z1(i)))*1.e2  !cm
         if(deltaz .lt. 1)deltaz=1.
         rhomedio(i-1)=1.e-1*real(mass(i)-mass(i-1))/deltaz
      end do
      rho1(1)=1.5*rhomedio(1)-0.5*rhomedio(2)
      do i=2,ntau
         rho1(i)=0.5*rhomedio(i-1)+0.5*rhomedio(i)
      end do  
      Bi=0
      do i=1,ntau
         k=ntau-i+1
c        19/06/19 epm: Now we are not using 'tau(k)' but 'tau_RH(k)'.
c        atmos(i+nn)=tau(k)-alog10mu              !log tau
         atmos(i+nn)=real(dlog10(tau_RH(k)))-alog10mu   !log tau
         atmos(i+nn+ntau)=real(T(k))    !T
         atmos(i+nn+2*ntau)=real(nele(k)*T(k))*bol*1.e-6  !pe (cgs conversion e-/m^3 a e/cm^3)
         atmos(i+nn+3*ntau)=real(Vturb(k))*1.e2 !micro
         atmos(i+nn+4*ntau)=B(k)          !B
         Bi=Bi+B(k)
         atmos(i+nn+5*ntau)=-real(vz(k))*1.e2   !Los Velocity cm/s redshifted positive
         atmos(i+nn+6*ntau)=gamma(k)      !inclination phi(k)
         atmos(i+nn+7*ntau)=phi(k)        !azimuth
         z(i)=real(z1(k))*1.e-3                 !Km
         rho(i)=rho1(k)
         if(HEon .eq. 1)then
            pg(i)=real(mass(k))*gravity*1.e-1
         else
            tsi=real(T(k))
            psi=real(nele(k)*T(k))*bol*1.e-6      !pe (cgs conversion e-/m^3 a e/cm^3)
            call PgmufromPeT_isub(tsi,psi,pgsi,pesomediosi,rhosi,kappasi)
            pg(i)=pgsi
         end if
c        print*,'read_atmos_RH 105 ',i,atmos(i+nn),atmos(i+nn+4*ntau),atmos(i+nn+6*ntau),atmos(i+nn+7*ntau)
c        print*,'read_atmos_RH 105 ',i,mass(i),pg(i),rho(i),gravity
      end do

c        if(ncontpg.eq.1 .and. natmos .eq. 1)pgcontorno=pg01
c        if(ncontpg.eq.1 .and. natmos .eq. 2)pgcontorno=pg02
c        if(ncontpg.ne.1)pgcontorno=gravity*mass(1)*1.e-1
c        if(ncontpg.ne.1) then
c          print*,'check boundary condition in Pg'
c          print*,'at read_atmos_RH'
c          stop
c        endif
c        pg1(1)=pgcontorno
c        do i=2,ntau
c          pg1(i)=pg1(i-1)+gravity*(mass(i)-mass(i-1))*1.e-1
c        end do
c        do i=1,ntau
c          pg(i)=pg1(ntau-i+1)
c        end do
        if(natmos .eq. 1 .and. Bi .gt. 0)nohaycampo1=1
        if(natmos .eq. 2 .and. Bi .gt. 0)nohaycampo2=1

c     pmusum=0.0
c     asum=0.0
c     do j=1,92
c        jj=j
c        call neldatb(jj,0.,wgt,abu,ei1,ei2)
c        pmusum=pmusum+wgt*abu
c        asum=asum+abu
c     end do

c     print*,'Vmac= ',Vmac,'Filling Factor=',fill,'Stray light=',strayfac 
c     print*,'read_atmos atmos(1)',atmos(1)
      atmos(8*ntau+1+nn)=Vmac
      atmos(8*ntau+2+nn)=fill
      atmos(16*ntau+5)=strayfac

      return
      end
