c departures
c calls RH to evaluate the departure coefficients
c RH was modified to write departure coefficients in ASCII 
c The modified subroutine was  pops_xdr.c that now writes
c the file:    test_pops_xdr_2.txt  (lines from 109 to 125)
c and          exit RH              (line 129)
c
c it is called by fperfil2.f
c Basilio 12/01/2018
c
c 15/05/19 epm: New routines to run RH executables as functions.
c 31/07/19 epm: New routine to pass Barklem data to RH.
c 04/04/20 epm: Cosine of heliocentric angle as input to solveray().
c 05/05/20 epm: Output on screen centralized.
c 07/07/20 epm: New routine to pass Kurucz data to RH.
c 08/08/20 epm: Use common/niveles instead of passing arguments.
c               New routine to pass abundance values to RH.
c 10/10/20 epm: New routine to pass the wavetable to RH.
c_____________________________________________________________________________

      subroutine departures(label_ID_model,RH_model,
     &                      RH_magneticfield,natmos,atmosLG,ntau,
     &                      ntotal_lines,beta1dbl,beta2dbl,stok_RH)

      implicit real*4 (a-h,o-z)
      include 'PARAMETER'
      parameter (kt8=8*kt+2,kt16=16*kt+5,kt11=11*kt+2,kld4=4*kld,kt12=11*kt+3)
c     06/06/19 epm: Array size for 'departure'.
      parameter (kt500=500*kt)
c     10/06/19 epm: Arrays for the spectrum.
      parameter (n_wave_NLTE=50000)
c     08/08/20 epm: Number of abundances for RH.
      parameter (narh=99)

c     Departure coefficients.
      integer nlevels,natmos
      integer nlow_i(kl),nup_i(kl),linea_nlte(kl)
      integer ntau,ntotal_lines
      real*8 bol
      real*4 stok_RH(kld4)  !,tau5RH(kt),tauRH_step(kt)
c     real*4 beta1(kl,kt),beta2(kl,kt),departure(kt500)
      real*8 beta1dbl(kl,kt),beta2dbl(kl,kt)
c     06/06/19 epm: Double precision for departure coefficients.
      real*8 departure(kt500)
      character*100 label_ID_model,RH_model,RH_magneticfield
      character*2 atom_nlte(kl),atomi
      integer istatus,istatus1,istatus2
      integer nspect

      real*4 tau(kt),T(kt),Pe(kt),Pg(kt),z(kt),ro(kt)
      real*4 Vmac,fill,strayfac
c     real*4 atmosnew(kt16)
c     real*4 pe1_change(kt),pe2_change(kt),pg1_change(kt),pg2_change(kt)
c     real*4 tauoriginal(kt)

c     Para la atmosfera.
      real*4 atmosLG(kt12)
c     real*4 pg1(kt),z1(kt),ro1(kt),pg2(kt),z2(kt),ro2(kt),z1i,z0,zn

c     15/05/19 epm: New functions and variables to run RH applications.
      integer*4 rhf1d, solveray
      integer*4 atmos_nspace, atom_nlevel, spectrum_nspect
      integer*8 cpu1, cpu2, wall1, wall2
      real*8    wave_NLTE_vacuum(n_wave_NLTE)
      real*8    SI_NLTE(n_wave_NLTE), SQ_NLTE(n_wave_NLTE)
      real*8    SU_NLTE(n_wave_NLTE), SV_NLTE(n_wave_NLTE)
      real*8    xmu
      real*4    alog10mu

c     31/07/19 brc: Barklem data (read from common/brklm).
      integer*4 ntotblends, atom_arr(kl), istage_arr(kl)
      real*4    alfa_arr(kl), sigma_arr(kl), wave_arr(kl)

c     08/08/20 epm: Abundance values (read from common/abundances).
      real*4 abu(narh)

c     05/05/20 epm: Command line flags and error message.
      integer*4 flagv, flagquiet
      character*200 msg

c     real*4 atmos1LGold(kt12)
c     common/atmos1LGold/atmos1LGold
c     common/ichange/ichange

      data ivez/0/

c     data (pe1_change(i), i=1,kt)/kt*1.0/
c     data (pe2_change(i), i=1,kt)/kt*1.0/
c     data (pg1_change(i), i=1,kt)/kt*1.0/
c     data (pg2_change(i), i=1,kt)/kt*1.0/

      common/istatus12/istatus1,istatus2

c     31/07/19 epm: Use this common to pass Barklem data to RH.
      common/brklm/ntotblends,atom_arr,istage_arr,alfa_arr,sigma_arr,wave_arr
c     04/04/20 epm: Logaritm of cosine of heliocentric angle.
      common/alog10mu/alog10mu   !alog10mu=alog10(xmu) from desire.f
c     05/05/20 epm: Common with the command line flags.
      common/commandline/flagv,flagquiet
c     08/08/20 epm: Line data from spectra().
      common/niveles/nlow_i,nup_i,atom_nlte,linea_nlte
c     08/08/20 epm: Use this common to pass abundance values to RH.
      common/abundances/abu
c     11/11/20 brc: Common from desire.f to pass H populations.
      common/scalemodelRH/imassortau  !0=logtau, 1=mass column

c     common/zetas/pg1,z1,ro1,pg2,z2,ro2  !para el calculo de z,pg y ro en amp22 y amp22err

c     common/NLTEdepartures/beta1,beta2
c     common/NLTEprofiles/stok_RH
c     common/atmosSIRfromRH/atmosnew,pe1_change,pe2_change,pg1_change,pg2_change
c     common/tauoriginal/tauoriginal

      ivez=ivez+1
      bol=1.3806488d-16  !erg/s
      istatus=0

      if(imassortau .eq. 0)then
          call write_atmos_RH_tau(label_ID_model,RH_model,RH_magneticfield,
     &                            atmosLG,ntau,natmos)
      else
          call write_atmos_RH(label_ID_model,RH_model,RH_magneticfield,
     &                        atmosLG,ntau)
      end if

      Vmac=atmosLG(8*ntau+1)
      fill=atmosLG(8*ntau+2)
      strayfac=atmosLG(11*ntau+3)

c     Passing data to RH.
      if (ivez .eq. 1) then
c        31/07/19 epm: New routine to pass Barklem data to RH.
         call sirbarklem(ntotblends, atom_arr, istage_arr, nlow_i, nup_i,
     &                   alfa_arr, sigma_arr, wave_arr)
c        07/07/20 epm: New routine to pass Kurucz data to RH.
         call lines_kurucz()
c        08/08/20 epm: New routine to pass abundance values to RH.
         call sirabundances(narh, abu)
c        10/10/20 epm: New routine to pass the wavetable to RH.
         call lines_wave()
      end if

c     15/05/19 epm: New functions to run RH applications.
      call error(KLINE,'','@CALLING rhf1d')
      if (flagv .eq. 1) then
         call error(KLINE,'',
     &        '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
         call chrono(cpu1, wall1)
      end if
      istatus1 = rhf1d(atmos_nspace, atom_nlevel, kt500, departure)
      if (flagv.eq.1) then
         call chrono(cpu2, wall2)
         write(msg,'(a,i8,a,i8,a)') 'rhf1d: CPU  time (ms) = ', cpu2 - cpu1,
     &                           '\n rhf1d: Wall time (ms) = ', wall2 - wall1
         call error(KLINE,'',msg)
         call error(KLINE,'',
     &        '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
      end if
      if (istatus1 .eq. -1) then
         write(msg,'(a,i6,a,i6,a)') 'Departure coefficients array with'
     &   //    ' insufficient memory.\n Allocated ', kt500, ' bytes'
     &   //    ' but ', atmos_nspace*atom_nlevel,' are needed'
         call error(KSTOP,'departures',msg)
      end if
      if (atmos_nspace .ne. ntau) then
         call error(KSTOP,'departures','The value of ntau (departures.f)'
     &   //         ' and atmos_nspace (rhf1d.c) is different')
      end if
      ntau = atmos_nspace
      nlevels = atom_nlevel

c     15/05/19 epm: Al convertir 'rhf1d' en funcion rhf1d() el error recogido
c     no es el de la llamada al sistema (que implicaba un aborto de RH por un
c     problema de convergencia) sino el retorno de la funcion rhf1d() que
c     ahora es 0 si correcto o -1 si array con memoria insuficiente. Por lo
c     tanto comento todo el trozo siguiente ya que es inutil.
c
c     ichange=0
c     18/06/19 brc: Si RH no converge, quizas es un problema de no equilibrio
c     hidrostatico. Rellamamos a RH poniendo antes la atmosfera en equilibrio.
c     if (istatus1 .ne. 0) then
c        call error(KWARN,'departures','Error running rhf1d:'
c    &   //         ' check "keyword.input" or change input model.\n'
c    &   //         ' Running HE in RH to avoid RH-nonconvergence')
c        do k=1,ntau
c           tau(k)=atmosLG(k)
c           T(k)=atmosLG(ntau+k)
c           Pe(k)=atmosLG(2*ntau+k)
c           Z(k)=atmosLG(8*ntau+2+k)
c           Pg(k)=atmosLG(9*ntau+2+k)
c           ro(k)=atmosLG(10*ntau+2+k)
c        end do
c        call equisubmu(ntau,tau,T,Pe,Pg,Z,ro)
c        do k=1,ntau
c           atmosLG(2*ntau+k)=Pe(k)
c           atmosLG(8*ntau+2+k)=Z(k)
c           atmosLG(9*ntau+2+k)=Pg(k)
c           atmosLG(10*ntau+2+k)=ro(k)
c        end do
c        if(imassortau .eq. 0)call write_atmos_RH_tau(label_ID_model,RH_model,
c    &                             RH_magneticfield,atmosLG,ntau,natmos)
c        if(imassortau .eq. 1)call write_atmos_RH(label_ID_model,RH_model,
c    &                             RH_magneticfield,atmosLG,ntau)
c        call error(KLINE,'','@RECALLING rhf1d')
c        if (flagv .eq. 1) then
c           call error(KLINE,'',
c    &           '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
c           call chrono(cpu1, wall1)
c        end if
c        istatus1 = rhf1d(atmos_nspace, atom_nlevel, kt500, departure)
c        if (flagv.eq.1) then
c           call chrono(cpu2, wall2)
c           write(msg,'(a,i8,a,i8,a)') 'rhf1d: CPU  time (ms) = ', cpu2 - cpu1,
c    &                              '\n rhf1d: Wall time (ms) = ', wall2 - wall1
c           call error(KLINE,'',msg)
c           call error(KLINE,'',
c    &           '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
c        end if
c        if (istatus1 .ne. 0) then
c           call error(KWARN,'departures','Persistent error running rhf1d.\n'
c    &      //         ' Trying to continue.....')
c        end if
c     end if

c     27/05/19 epm: New functions to run RH applications.
      if (istatus1 .eq. 0) then
         if (flagv .eq. 1) then
            call error(KLINE,'','@CALLING solveray')  !only print it with '-v'
            call error(KLINE,'',
     &           '------------------------------------------------------------')
            call chrono(cpu1, wall1)
         end if
         xmu = 10**alog10mu   !input to solveray
         istatus2 = solveray(spectrum_nspect, n_wave_NLTE, wave_NLTE_vacuum,
     &                       SI_NLTE, SQ_NLTE, SU_NLTE, SV_NLTE, xmu)
         if (flagv .eq. 1) then
            call chrono(cpu2, wall2)
            write(msg,'(a,i8,a,i8,a)') 'solveray: CPU  time (ms) = ', cpu2 - cpu1,
     &                              '\n solveray: Wall time (ms) = ', wall2 - wall1
            call error(KLINE,'',msg)
            call error(KLINE,'',
     &           '------------------------------------------------------------')
         end if
         if (istatus2 .eq. -1) then
            write(msg,'(a,i6,a,i6,a)') 'Arrays for the spectrum with'
     &      //    ' insufficient memory.\n Allocated ', n_wave_NLTE, ' bytes'
     &      //    ' but ', spectrum_nspect,' are needed'
            call error(KSTOP,'departures',msg)
         end if
         nspect = spectrum_nspect
      endif

      if (istatus1 .eq. 0) then
         atomi=' '
         do iline=1,ntotal_lines
            if (linea_nlte(iline) .ne. 0) then   !only if NLTE line
c              06/06/19 epm: We don't need read_populations() anymore.
c              if (atom_nlte(iline) .ne. atomi) then   !only if unread atom
c                 atomi=atom_nlte(iline)
c                 call read_populations(ntau,atomi,nlevels,nlow_i,nup_i,departure)
c              end if
               do j=1,ntau
                  j1=(j-1)*nlevels
                  k1=ntau-j+1          !RH's atmosphere model runs from upper layers to lower ones
                  beta1dbl(iline,k1)=departure(j1+nlow_i(iline)+1)
                  beta2dbl(iline,k1)=departure(j1+nup_i(iline)+1)
                  if(beta1dbl(iline,k1) .lt. 1.d-10) beta1dbl(iline,k1)=1.d-10
                  if(beta2dbl(iline,k1) .lt. 1.d-10) beta2dbl(iline,k1)=1.d-10
                  if(beta1dbl(iline,k1) .gt. 1.d10)  beta1dbl(iline,k1)=1.d10
                  if(beta2dbl(iline,k1) .gt. 1.d10)  beta2dbl(iline,k1)=1.d10
               end do
            end if
         end do

c        Reading the NLTE spectrum.
         call read_spectrum(nspect, n_wave_NLTE, wave_NLTE_vacuum,
     &                      SI_NLTE, SQ_NLTE, SU_NLTE, SV_NLTE, stok_RH)
      end if

      return
      end

c-----------------------------------------------------------------------------

        subroutine read_populations(ntau,atom_nlte_i,nlevels,nlow_i,nup_i,departure)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        character*2 atom_nlte_i
        character*100 file_pops
        integer nlow_i(*),nup_i(*),ntau,itotalmax
        real*4 departure(*),dep_ii
        data nlevelsmax/500/

        ican=57

        itotalmax=nlevelsmax*ntau  
        if(atom_nlte_i(2:2) .eq. " ")then 
           file_pops='pops.'//atom_nlte_i(1:1)//'.out.txt'
        else
           file_pops='pops.'//atom_nlte_i//'.out.txt'
        endif

        open(ican,file=file_pops,status='old',err=999)

        read(ican,*,err=999)i1,i2,dd
        itotal=1
        do while(itotal.lt.itotalmax)
           read(ican,*,end=221,err=999)i1,i2,dd
           itotal=itotal+1
        end do
221     nlevels=itotal/ntau
        close(ican)

        open(ican,file=file_pops,status='old',err=999)
        do ii=1,itotal
           read(ican,*,err=999)i1,i2,dep_ii
           if(dep_ii .lt. 1.e-10)dep_ii=1.e-10
           if(dep_ii .gt. 1.e10)dep_ii=1.e10
           departure(ii)=dep_ii
        end do
        close(ican)

        return

999     call error(KSTOP,'read_populations','Error reading '//file_pops)

        end

c-----------------------------------------------------------------------------

        subroutine read_spectrum(nspect, n_wave_NLTE, wave_NLTE_vacuum,
     +                           SI_NLTE, SQ_NLTE, SU_NLTE, SV_NLTE, stok_RH)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
c       parameter (n_wave_NLTE=50000)
        parameter (kt8=8*kt+2,kt16=16*kt+5,kt11=11*kt+2,kl4=4*kl,kld4=4*kld)

        integer nspect, n_wave_NLTE
        integer nlam_NLTE,iprimera,iLTE,istage,nlow,nup,ikk !,nxx
        real*8  wave_NLTE(n_wave_NLTE),wave_NLTE_vacuum(n_wave_NLTE)
        real*8  SI_NLTE(n_wave_NLTE),SQ_NLTE(n_wave_NLTE),SU_NLTE(n_wave_NLTE),SV_NLTE(n_wave_NLTE) 
        integer nlam_LTE,igrado
        real*8  wave_LTE(kld)
        real*4  SI_LTE(kld),SQ_LTE(kld),SU_LTE(kld),SV_LTE(kld),continuo
        real*4  stok_RH(*)
        real*8  x0,x1,x2,x,y0,y1,y2,y,delta2          !,conhsra_RH,cc_RH
        real*4  continuohNLTEarr(kld)

        real*4  stok(kld4),sigd(kld4)                  !,conhsra,ccc
        real*4  dlamda0(kl),dlongd(kld)       !,con_i(kl),dlamda00(kl)
        integer nlin(kl),npas(kl)
        integer ist(4),nble(kl)

        integer mult(2)
        real*4  tam(2),loggf,wlengt1,zeff,alfa,sigma,energy

        character design(2)*1,atom*2

        real*8  wlengtdbl, wavedbl_arr(kl)
        
        common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/observaciones/stok,sigd
        common/ldeo/dlongd,dlamda0
        common/wave_LTEnew/wave_LTE
        common/numero_LTE/nlam_LTE
        common/continuosNLTEarr/continuohNLTEarr
        common/wavearrdble/wavedbl_arr

        data iprimera/0/

        iprimera=iprimera+1

        nlam_NLTE = nspect

        call vacuum_to_air(nlam_NLTE,wave_NLTE_vacuum,wave_NLTE) !wave_NLTE is air wavelength

        if(iprimera .eq. 1)then
c          Datos de la linea.
           ifrec=0
           ikk=1
c          Leemos la primera componente de cada blend.
c          iln=numero de la linea,ntl num.total (sin blends).
           do iln=1,ntl
              wlengtdbl=wavedbl_arr(ikk)
              ikk=ikk+nble(iln)

              do ii=1,npas(iln)
                 ifrec=ifrec+1
                 wave_LTE(ifrec)=wlengtdbl+dlongd(ifrec)*1.d-3 
              end do
            end do
            nlam_LTE=ifrec
        end if

c       Interpolating.
        iLTE=1
        iii=0
        do iln=1,ntl
           do i=1,npas(iln)
              iii=iii+1
              x=wave_LTE(iii)
              continuo=continuohNLTEarr(iii)           !*x*x*1.e-28/30.
              if( continuo .lt. 1.e-30)continuo=1.e-30

              call locating(wave_NLTE,nlam_NLTE,x,iLTE)
              if(iLTE.eq.1)iLTE=2
              ii0=iLTE-1
              x0=wave_NLTE(ii0)
              ishift=0
              ii1=iLTE+ishift
              do while(dabs(wave_NLTE(ii1)-x0) .lt. 1.e-8 )
                 ishift=ishift+1
                 ii1=iLTE+ishift
                 if(ii1 .eq. n_wave_NLTE)then
                    call error(KSTOP,'read_spectrum',
     &                         'Problems interpolating wave NLTE?')
                 end if
              end do
              x1=wave_NLTE(ii1)
              ishift=0
              ii2=ii1+1+ishift
              do while(dabs(wave_NLTE(ii2)-x1) .lt. 1.e-8 )
                 ishift=ishift+1
                 ii2=ii1+ishift
                 if(ii2 .eq. n_wave_NLTE)then
                    call error(KSTOP,'read_spectrum',
     &                         'Problems interpolating wave NLTE?')
                 end if
              end do
              x2=wave_NLTE(ii2)

              delta2=(x2-x0)/(x1-x0)
              igrado=2
              if( delta2 .gt. 2.5 .or. delta2 .lt. 1.2)igrado=1

              y0=SI_NLTE(ii0)
              y1=SI_NLTE(ii1)
              y2=SI_NLTE(ii2)

              if(igrado.eq.2)call INTERPOLATING2(x0,x1,x2,y0,y1,y2,x,y)
              if(igrado.eq.1)call INTERPOLATING(x1,x2,y1,y2,x,y)

              SI_LTE(iii)=real(y)/continuo

              y0=SQ_NLTE(ii0)
              y1=SQ_NLTE(ii1)
              y2=SQ_NLTE(ii2)
              if(igrado .eq.2)call INTERPOLATING2(x0,x1,x2,y0,y1,y2,x,y)
              if(igrado .eq.1)call INTERPOLATING(x1,x2,y1,y2,x,y)
              SQ_LTE(iii)=real(y)/continuo

              y0=SU_NLTE(ii0)
              y1=SU_NLTE(ii1)
              y2=SU_NLTE(ii2)
              if(igrado .eq.2)call INTERPOLATING2(x0,x1,x2,y0,y1,y2,x,y)
              if(igrado .eq.1)call INTERPOLATING(x1,x2,y1,y2,x,y)
              SU_LTE(iii)=real(y)/continuo

              y0=SV_NLTE(ii0)
              y1=SV_NLTE(ii1)
              y2=SV_NLTE(ii2)
              if(igrado .eq.2)call INTERPOLATING2(x0,x1,x2,y0,y1,y2,x,y)
              if(igrado .eq.1)call INTERPOLATING(x1,x2,y1,y2,x,y)
              SV_LTE(iii)=real(y)/continuo
           end do
        end do

        call sfromiquv(ist,nlam_LTE,SI_LTE,SQ_LTE,SU_LTE,SV_LTE,stok_RH)

        return
        end

c-----------------------------------------------------------------------------

        subroutine INTERPOLATING(x1,x2,y1,y2,x,y)

        implicit real*4 (a-h,o-z)
        real*8 x,x1,x2,y,y1,y2
        real*8 deltax

        deltax=(x2-x1)
        if(dabs(deltax) .lt. 1.d-5)then
           y=y1
        else
           y=y1+(x-x1)*(y2-y1)/deltax
        endif

        return
        end

c-----------------------------------------------------------------------------

        subroutine INTERPOLATING2(x1,x2,x3,y1,y2,y3,x,y)

        implicit real*4 (a-h,o-z)
        real*8 x,x1,x2,x3,y,y1,y2,y3
        real*8 x21,x13,x32,a,b,c

        x21=x2-x1
        x13=x1-x3
        x32=x3-x2

        if(dabs(x13) .lt. 1.d-5)then
           y=y2
        else
           c=((y2-y1)/x21-(y3-y2)/x32)/x13
           b=(y3-y2)/x32-c*(x3+x2)
           a=y2-x2*(b+c*x2)
           y=a+x*(b+c*x)
        endif

        return
        end

c-----------------------------------------------------------------------------

        subroutine locating(wave_NLTE,nlam_NLTE,x,iLTE)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        real*8 wave_NLTE(*),x,ww
        integer iLTE,nlam_NLTE,ii
        character*300 msg

        ii=iLTE
        ww=wave_NLTE(ii)

        if(dabs(ww-x) .lt. 1.d-5)then
           iLTE=ii
        else if(ww .ge. x)then
           do while (wave_NLTE(ii) .gt. x  .and. ii.ge.2)
              ii=ii-1
           end do
           iLTE=ii
        else if(ww .lt. x)then
           do while (wave_NLTE(ii) .lt. x  .and. ii.lt.nlam_NLTE)
              ii=ii+1
           end do
           iLTE=ii-1
        endif

        if(iLTE .lt. 1 )iLTE=1
        if(iLTE .gt. (nlam_NLTE-1) )iLTE=nlam_NLTE-1
        if(iLTE .gt. 1 .and. iLTE .lt.  nlam_NLTE)then
           if( wave_NLTE(iLTE) .gt. x .or. wave_NLTE(iLTE+1) .lt. x)then
              write(msg,'(3(a,f10.3),a)') 'Wavelength is not'
     &        //    ' correctly located for interpolation.\n'
     &        //    ' Looking for (SIR)', x, ' between'
     &        //    ' (RH)', wave_NLTE(iLTE),' and ', wave_NLTE(iLTE+1), '\n'
     &        //    ' Possibly the SIR grid file is larger than'
     &        //    ' (or displaced respect to)\n the additional'
     &        //    ' LTE grid file for RH'
              call error(KSTOP,'locating',msg)
           end if
        end if

        return
        end

c-----------------------------------------------------------------------------

        subroutine vacuum_to_air(nwave,w_vac,w_air)

        implicit real*4 (a-h,o-z)
        integer nwave,i
        real*8 w_vac(*),w_air(*),refrax,ss,ww

c       12/06/19 epm: 'w_vac' entra en nm, 'w_air' sale en Angstrom.

        do i=1,nwave
           ww=w_vac(i)*10.d0   !nm --> Angstrom
           ss=1.d8/ww/ww
           refrax = 1.d0 + 0.0000834254d0 + 0.02406147d0 / (130.d0 - ss) + 0.00015998d0 / (38.9d0 - ss)
           w_air(i)=ww/refrax
        end do

        return
        end

c-----------------------------------------------------------------------------

        subroutine air_to_vacuum(nwave,w_air,w_vac)

        implicit real*4 (a-h,o-z)
        integer nwave,i
        real*8 w_vac(*),w_air(*),refrax,ss,ww

c       12/06/19 epm: 'w_air' entra en Angstrom, 'w_vac' sale en nm.

        do i=1,nwave
           ww=w_air(i)
           ss=1.d8/ww/ww
           refrax=1.d0+8.336624212083d-5+2.408926869968d-2/(130.1065924522d0 - ss) +
     &            1.599740894897d-4 / (38.92568793293d0 - ss)
           w_vac(i)=ww*refrax/10.d0  !Angstrom --> nm
        end do

        return
        end

c ----------------------------------------------------------------------------

        FUNCTION w_vac_fn(w_air)

        real*8 w_vac_fn,w_air,refrax,ss

c       12/06/19 epm: 'w_air' entra en Angstrom, 'w_vac_fn' sale en Angstrom.

        ss=1.d8/w_air/w_air
        refrax=1.d0+8.336624212083d-5+2.408926869968d-2/(130.1065924522d0 - ss) +
     &         1.599740894897d-4 / (38.92568793293d0 - ss)
        w_vac_fn=w_air*refrax

        return
        end
