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
c ________________________________________________________________________

      subroutine departures(label_ID_model,RH_model,
     +RH_magneticfield,natmos,atmosLG,ntau,
     +ntotal_lines,linea_nlte,atom_nlte,nlow_i,nup_i,beta1,beta2,stok_RH)  !,pe_dep

      implicit real*4 (a-h,o-z)
      include 'PARAMETER'   !por kt,kn,kl,kld
      parameter (kt8=8*kt+2,kt16=16*kt+5,kt11=11*kt+2,kld4=4*kld,kt12=11*kt+3) 
c     03/06/19 epm: Array size for 'departure'.
      parameter (kt500=500*kt)
c     10/06/19 epm: Arrays for the spectrum.
      parameter (n_wave_NLTE=50000)

c     Departure coefficients.
      integer nlevels,natmos
      integer linea_nlte(kl)
      integer nlow_i(kl),nup_i(kl),ntau,ntotal_lines
      real*8 bol
      real*4 stok_RH(kld4)   !,tau5RH(kt),tauRH_step(kt)
c     real*4 beta1(kl,kt),beta2(kl,kt),departure(kt500)  !,pe_dep(kt)
      real*4 beta1(kl,kt),beta2(kl,kt)                   !,pe_dep(kt)
c     06/06/19 epm: Double precision for departure coefficients.
      real*8 departure(kt500)
      character*100 label_ID_model,RH_model,RH_magneticfield !BRC-RH Jun 20 2017
      character*2 atom_nlte(kl),atomi
      integer istatus,istatus1,istatus2
      integer nspect

      real*4 tau(kt),T(kt),Pe(kt),Pg(kt),z(kt),ro(kt)
      real*4 Vmac,fill,strayfac
c     real*4 atmosnew(kt16)
c     real*4 pe1_change(kt),pe2_change(kt),pg1_change(kt),pg2_change(kt)
c     real*4 tauoriginal(kt)
      character*100 vprint,ruta,ruta2,ruta3,filewavename

c     Para la atmosfera.
      real*4 atmosLG(kt12)
c     real*4 pg1(kt),z1(kt),ro1(kt),pg2(kt),z2(kt),ro2(kt),z1i,z0,zn

c     31/07/19 brc: New variables in a common to pass Barklem data to RH.
      integer*4 ntotblends, atom_arr(kl), istage_arr(kl)
      real*4    alfa_arr(kl), sigma_arr(kl), wave_arr(kl)

c     15/05/19 epm: New functions and variables to run RH applications.
      integer*4 rhf1d, solveray, system
      integer*4 atmos_nspace, atom_nlevel, spectrum_nspect
      integer*8 cpu1, cpu2, wall1, wall2
      real*8    wave_NLTE_vacuum(n_wave_NLTE)
      real*8    SI_NLTE(n_wave_NLTE), SQ_NLTE(n_wave_NLTE)
      real*8    SU_NLTE(n_wave_NLTE), SV_NLTE(n_wave_NLTE)

c     real*4 atmos1LGold(kt12)
c     common/atmos1LGold/atmos1LGold
c     common/ichange/ichange

      data ivez/0/

c     data (pe1_change(i), i=1,kt)/kt*1.0/
c     data (pe2_change(i), i=1,kt)/kt*1.0/
c     data (pg1_change(i), i=1,kt)/kt*1.0/
c     data (pg2_change(i), i=1,kt)/kt*1.0/

      common/vprint/vprint
      common/ruta/ruta,ruta2,ruta3,filewavename
      common/istatus12/istatus1,istatus2
c     31/07/19 brc: New common to pass Barklem data to RH.
      common/brklm/ntotblends,atom_arr,istage_arr,alfa_arr,sigma_arr,wave_arr

c     common/zetas/pg1,z1,ro1,pg2,z2,ro2    !para el calculo de z,pg y ro en amp22 y amp22err

c     common/NLTEdepartures/beta1,beta2
c     common/NLTEprofiles/stok_RH
c     common/atmosSIRfromRH/atmosnew,pe1_change,pe2_change,pg1_change,pg2_change  
c     common/tauoriginal/tauoriginal

      ivez=ivez+1
      bol=1.3806488d-16  !erg/s
      istatus=0 
c     istatus =system("rm spectrum*")

      call write_atmos_RH(label_ID_model,RH_model,RH_magneticfield,atmosLG,ntau)  !BRC-RH Jun 19 2017

c     11/06/19 epm: QUEDA PENDIENTE SUPRIMIR ESTA LLAMADA HASTA COMPROBAR QUE
c     ES NECESARIA. PARA ESO HAY QUE EJECUTAR CON REDISTRIBUCION PARCIAL.
      istatus = system("rm PRD*.dat")

      Vmac=atmosLG(8*ntau+1)
      fill=atmosLG(8*ntau+2)
      strayfac=atmosLG(11*ntau+3)

c     31/07/19 epm: New routine to pass Barklem data to RH.
      call desirelines(ntotblends, atom_arr, istage_arr, nlow_i, nup_i,
     &                 alfa_arr, sigma_arr, wave_arr)

c     15/05/19 epm: New functions to run RH applications.
      write(0, *) ''
      write(0, *) '$ CALLING rhf1d from "librh.a"'
      write(0, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      call chrono(cpu1, wall1)
c     istatus1 = system(ruta)
      istatus1 = rhf1d(atmos_nspace, atom_nlevel, kt500, departure)
      call chrono(cpu2, wall2)
      write(0, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(0, '(1x,a,8x,i8)') 'rhf1d: status = ', istatus1
      write(0, '(1x,a,i8)') 'rhf1d: CPU  time (ms) = ', cpu2 - cpu1
      write(0, '(1x,a,i8)') 'rhf1d: Wall time (ms) = ', wall2 - wall1
      write(0, *) ''
      if (istatus1 .ne. 0) then
         print*,'*** ERROR #1 in departures.f: '
         print*,'Departure coefficients array with insufficient memory.'
         print'(1x,a,i6,a,$)','Allocated ', kt500, ' bytes but '
         print'(i6,a)',atmos_nspace*atom_nlevel, ' are needed.'
         stop ' '
      end if
      if (atmos_nspace .ne. ntau) then
         print*,'*** ERROR #2 in departures.f: '
         print'(1x,a,$)','The value of "ntau" (departures.f)'
         print*,'and "atmos_nspace" (rhf1d.c) is different.'
         stop ' '
      end if
      ntau = atmos_nspace
      nlevels = atom_nlevel

c     ichange=0
c     18/06/19 brc: Si RH no converge, quizas es un problema de no equilibrio
c     hidrostatico. Rellamamos a RH poniendo antes la atmosfera en equilibrio.
      if (istatus1 .ne. 0) then
         print*,'Error running rhf1d: check keyword.input or change input model'
         print*,'Running HE in DeSIRe (to avoid RH-nonconvergence)!'
         do k=1,ntau 
            tau(k)=atmosLG(k)
            T(k)=atmosLG(ntau+k)
            Pe(k)=atmosLG(2*ntau+k)
c           Z(k)=atmosLG(8*ntau+2+k)
            Pg(k)=atmosLG(9*ntau+2+k)
c           ro(k)=atmosLG(10*ntau+2+k)
         end do
         call equisubmu(ntau,tau,T,Pe,Pg,Z,ro)
         do k=1,ntau
            atmosLG(2*ntau+k)=Pe(k)
            atmosLG(8*ntau+2+k)=Z(k)
            atmosLG(9*ntau+2+k)=Pg(k)
            atmosLG(10*ntau+2+k)=ro(k)
         end do
         call write_atmos_RH(label_ID_model,RH_model,RH_magneticfield,atmosLG,ntau)
         write(0, *) ''
         write(0, *) '$ RECALLING rhf1d from "librh.a"'
         write(0, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         call chrono(cpu1, wall1)
c        istatus1 = system(ruta)
         istatus1 = rhf1d(atmos_nspace, atom_nlevel, kt500, departure)
         call chrono(cpu2, wall2)
         write(0, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write(0, '(1x,a,8x,i8)') 'rhf1d: status = ', istatus1
         write(0, '(1x,a,i8)') 'rhf1d: CPU  time (ms) = ', cpu2 - cpu1
         write(0, '(1x,a,i8)') 'rhf1d: Wall time (ms) = ', wall2 - wall1
         write(0, *) ''
         if (istatus1 .ne. 0) print*, 'Error running rhf1d: check keyword.input or change input model'
      end if

c     if (ivez .eq. 1. .and. istatus1 .ne. 0) stop 'STOP Error running rhf1d: check keyword.input or change input model'
c     if (istatus1 .ne. 0) stop 'STOP Error running rhf1d: check keyword.input or change input model'

c     27/05/19 epm: New functions to run RH applications.
      if (istatus1 .eq. 0) then
         write(0, *) ''
         write(0, *) '$ CALLING solveray from "librh.a"'
         write(0, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         call chrono(cpu1, wall1)
c        istatus2 = system(ruta2)
         istatus2 = solveray(spectrum_nspect, n_wave_NLTE, wave_NLTE_vacuum,
     &                       SI_NLTE, SQ_NLTE, SU_NLTE, SV_NLTE)
         call chrono(cpu2, wall2)
         write(0, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write(0, '(1x,a,8x,i8)') 'solveray: status = ', istatus2
         write(0, '(1x,a,i8)') 'solveray: CPU  time (ms) = ', cpu2 - cpu1
         write(0, '(1x,a,i8)') 'solveray: Wall time (ms) = ', wall2 - wall1
         write(0, *) ''
         if (istatus2 .ne. 0) then
            print*,'*** ERROR #3 in departures.f: '
            print*,'Arrays for the spectrum with insufficient memory.'
            print'(1x,a,i6,a,$)','Allocated ', n_wave_NLTE, ' bytes but '
            print'(i6,a)',spectrum_nspect, ' are needed.'
            stop ' '
         end if
         nspect = spectrum_nspect
      endif

c     if (istatus2 .ne. 0) print*,'Error running solveray: check keyword.input'
c     if (ivez .eq. 1. .and. istatus2 .ne. 0) stop 'STOP Error running solveray: check keyword.inputor change input model '

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
                  beta1(iline,k1)=real(departure(j1+nlow_i(iline)+1))
                  beta2(iline,k1)=real(departure(j1+nup_i(iline)+1))
                  if(beta1(iline,k1) .lt. 1.e-10) beta1(iline,k1)=1.e-10
                  if(beta2(iline,k1) .lt. 1.e-10) beta2(iline,k1)=1.e-10
                  if(beta1(iline,k1) .gt. 1.e10)  beta1(iline,k1)=1.e10
                  if(beta2(iline,k1) .gt. 1.e10)  beta2(iline,k1)=1.e10
               end do
            end if
         end do

c        Reading the NLTE spectrum.
c        call read_spectrum(stok_RH)
         call read_spectrum(nspect, n_wave_NLTE, wave_NLTE_vacuum,
     &                      SI_NLTE, SQ_NLTE, SU_NLTE, SV_NLTE, stok_RH)
      end if

      return
      end

c --------------------------------------------------------------------------------

        subroutine read_populations(ntau,atom_nlte_i,nlevels,nlow_i,nup_i,departure)
        implicit real*4 (a-h,o-z)

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
        print*,' reading pops from ',file_pops

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

999     print*,'STOP: Incorrect format in the file containing departures.'
        stop

        end

c --------------------------------------------------------------------------------

c       subroutine read_spectrum(stok_RH)
        subroutine read_spectrum(nspect, n_wave_NLTE, wave_NLTE_vacuum,
     +                           SI_NLTE, SQ_NLTE, SU_NLTE, SV_NLTE, stok_RH)
        implicit real*4 (a-h,o-z)
        include 'PARAMETER'   !por kt,kn,kl,kld

c       parameter (n_wave_NLTE=50000)
        parameter (kt8=8*kt+2,kt16=16*kt+5,kt11=11*kt+2,kl4=4*kl,kld4=4*kld)

        integer nspect, n_wave_NLTE
        integer nlam_NLTE,iprimera,iLTE,nxx,istage,nlow,nup,ikk
        real*8 wave_NLTE(n_wave_NLTE),wave_NLTE_vacuum(n_wave_NLTE)
        real*8 SI_NLTE(n_wave_NLTE),SQ_NLTE(n_wave_NLTE),SU_NLTE(n_wave_NLTE),SV_NLTE(n_wave_NLTE) 
        integer nlam_LTE,igrado
        real*8 wave_LTE(kld)
        real*4 SI_LTE(kld),SQ_LTE(kld),SU_LTE(kld),SV_LTE(kld),continuo
        real*4 stok_RH(*)
        character*100 ruta,ruta2,ruta3,filewavename
        real*8 x0,x1,x2,x,y0,y1,y2,y,delta2          !,conhsra_RH,cc_RH
        real*4 continuohNLTEarr(kld)

        real*4 stok(kld4),sigd(kld4)                  !,conhsra,ccc
        real*4 dlamda0(kl),dlongd(kld),dlamda00(kl)   !,con_i(kl)
        integer nlin(kl),npas(kl)
        integer ist(4),nble(kl)

        integer mult(2)
        real tam(2),loggf,wlengt1,zeff,alfa,sigma,energy
        character design(2)*1,atom*2

        common/ruta/ruta,ruta2,ruta3,filewavename
        common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/observaciones/stok,sigd
        common/ldeo/dlongd,dlamda0
        common/ldonew/dlamda00
        common/wave_LTEnew/wave_LTE
        common/numero_LTE/nlam_LTE
        common/continuosNLTEarr/continuohNLTEarr

        data iprimera/0/

        iprimera=iprimera+1

c       10/06/19 epm: We don't need this anymore.
c       print*,' reading RH spectrum from ',filewavename
c       open(1,file=filewavename)
c       read(1,*)nlam_NLTE
c       if (nlam_NLTE .gt.n_wave_NLTE) then
c          print*,'Please change n_wave_NLTE in read_spectrum (departures.f) because is smaller than ',nlam_NLTE
c          stop
c       end if
c       do iNLTE=1,nlam_NLTE
c          read(1,*)wave_NLTE_vacuum(iNLTE),SI_NLTE(iNLTE),SQ_NLTE(iNLTE),SU_NLTE(iNLTE),SV_NLTE(iNLTE)
c       end do
c       close(1)

        nlam_NLTE = nspect

        call vacuum_to_air(nlam_NLTE,wave_NLTE_vacuum,wave_NLTE) !wave_NLTE is air wavelength

        if(iprimera .eq. 1)then
c          datos de la linea
           ifrec=0
           ikk=1
           do iln=1,ntl        !iln=numero de la linea,ntl num.total
              nxx=nlin(ikk)
              call leelineasiii(nxx,atom,istage,wlengt1,zeff,energy,        !BRC 11/01/2018   
     &                         loggf,mult,design,tam,alfa,sigma,nlow,nup)     !BRC 11/01/2018
              dlamda00(iln)=wlengt1
              ikk=ikk+nble(iln)

c              ccc=conhsra(wlengt1)*conhsra_RH(wlengt1)
c              cc_RH=conhsra_RH(wlengt1)
c              con_i(iln)=conhsra(wlengt1)*real(conhsra_RH(wlengt1))           !*real(cc_RH)
c              print*,'departures conhsra ',conhsra(wlengt1)*real(conhsra_RH(wlengt1))

              do ii=1,npas(iln)
                 ifrec=ifrec+1
                 wave_LTE(ifrec)=wlengt1*1.d0+dlongd(ifrec)*1.d-3  
              end do          
            end do  
            nlam_LTE=ifrec
        end if   
     
c interpolating
        iLTE=1
        iii=0
      do iln=1,ntl
c          d00=dlamda00(iln)*1.e-14            
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
                  print*,'departures: Problems interpolating wave NLTE?'
                  stop
                end if  
             enddo                  
                   x1=wave_NLTE(ii1)
                   ishift=0
             ii2=ii1+1+ishift
             do while(dabs(wave_NLTE(ii2)-x1) .lt. 1.e-8 )
                ishift=ishift+1
                ii2=ii1+ishift
c               print*,'departures 250 ishift=',ishift
                if(ii2 .eq. n_wave_NLTE)then
                  print*,'departures: Problems interpolating wave NLTE?'
                  stop
                end if  
             enddo  
                   x2=wave_NLTE(ii2)
                   
                   delta2=(x2-x0)/(x1-x0)
                   igrado=2
                   if( delta2 .gt. 2.5 .or. delta2 .lt. 1.2)igrado=1
               
                   y0=SI_NLTE(ii0)
                   y1=SI_NLTE(ii1)
                   y2=SI_NLTE(ii2)  
                   
                   if(igrado .eq.2)call INTERPOLATING2(x0,x1,x2,y0,y1,y2,x,y) 
                   if(igrado .eq.1)call INTERPOLATING(x1,x2,y1,y2,x,y)

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
           
c ---------------------------------------------------------------------------------------------------

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

c ---------------------------------------------------------------------------------------------------

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
c      	   print*,'departures 312',c,x21,x13,x32,x1,x2,x3
      	   
           return
           end

c ---------------------------------------------------------------------------------------------------

           subroutine locating(wave_NLTE,nlam_NLTE,x,iLTE)
           implicit real*4 (a-h,o-z)

           real*8 wave_NLTE(*),x,ww
           integer iLTE,nlam_NLTE,ii
           
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
	       print*,'Wavelength is not correctly located for interpolation'
	       print*,'WARNING at locating (departure)' 
	       print*,'Looking for (SIR)',x,' between (RH)',wave_NLTE(iLTE),' and ',wave_NLTE(iLTE+1)
	       print*,'Possibly the SIR grid file is larger than (or displaced respect to) the additional LTE grid file for RH'
               stop
	     end if
           end if         	   
           return
           end

c ---------------------------------------------------------------------------------------------------

        subroutine vacuum_to_air(nwave,w_vac,w_air)
        implicit real*4 (a-h,o-z)

        integer nwave,i
        real*8 w_vac(*),w_air(*),refrax,ss,ww

        do i=1,nwave
           ww=w_vac(i)*10.d0   ! w_vac in nm, ww in Angstrom
           ss=1.d8/ww/ww
           refrax = 1.d0 + 0.0000834254d0 + 0.02406147d0 / (130.d0 - ss) + 0.00015998d0 / (38.9d0 - ss)
           w_air(i)=ww/refrax
        end do

        return
        end

c ---------------------------------------------------------------------------------------------------

        subroutine air_to_vacuum(nwave,w_air,w_vac)
        implicit real*4 (a-h,o-z)

        integer nwave,i
        real*8 w_vac(*),w_air(*),refrax,ss,ww

        do i=1,nwave
           ww=w_air(i)       !wvac in A
           ss=1.d8/ww/ww
           refrax=1.d0+8.336624212083d-5+2.408926869968d-2/(130.1065924522d0 - ss) + 
     &           1.599740894897d-4 / (38.92568793293d0 - ss)
                 w_vac(i)=ww*refrax    
              end do 
        return
        end

c ------------------------------------------------------------------------------------------------------          

          FUNCTION w_vac_fn(w_air)

        real*8 w_vac_fn,w_air,refrax,ss
           ss=1.d8/w_air/w_air  !wvac in A
           refrax=1.d0+8.336624212083d-5+2.408926869968d-2/(130.1065924522d0 - ss) + 
     &           1.599740894897d-4 / (38.92568793293d0 - ss)
           w_vac_fn=w_air*refrax   

        return
        end

c ---------------------------------------------------------------------------------------------------
c
c tau5RHfromtau5SIR evaluates optical depth scale of a model running RH (conversion)
c        subroutine tau5RHfromtau5SIR(tau5RH)
c        implicit real*4 (a-h,o-z)
c	include 'PARAMETER'   !por kt,kn,kl,kld
c	
c        real*4 tau_RHi,massi,z1i,Ti,nei,vzi,Vturbi
c	real*4 tau5RH(kt)
c	character*100 ruta,ruta2,ruta3,filewavename
c	common/ruta/ruta,ruta2,ruta3,filewavename
c	      
c        print*,'$ CALLING conversion from  ',ruta3   
c        istatus=system(ruta3)
c        if ( istatus .ne. 0 ) stop 'Error running conversion: check keyword.input' 
c
c        open(50,file='atmos.txt.out')
c        read(50,*)ndepths
c        do i=1,ndepths
c         read(50,*)tau_RHi,massi,z1i,Ti,nei,vzi,Vturbi !MKSA
c          tau5RH(ndepths-i+1)=alog10(tau_RHi)
c        end do 
c        close(50)
c
c	return	
c
c	end
