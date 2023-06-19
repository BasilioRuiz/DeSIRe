c
c File          _______ : desire.f
c Description   _______ : Stokes Inversion based on Response functions
c                         considering NLTE atomic populations evaluated
c                         with RH.
c                         Takes into account the position over the disk:
c                         it supposes that input/out models are in the Local
c                         Reference Frame (LRF). To sinthesize/invert
c                         input/output models are transformed to the Line
c                         of Sight Reference Frame (LoS).
c Project       _______ : DeSIRe
c Creation date _______ : 22/05/18
c Author        _______ : Basilio Ruiz Cobo (brc@iac.es)
c                         Han Uytenbroek, Carlos Quintero, David Orozco
c
c Versions <n.month> where n increases with the year starting in 1 for 2019:
c
c 1.00 (20/05/19): Starting version with system calls (2019 = year 1).
c 1.09 (09/09/19): First version with executables converted to functions.
c 1.10 (10/10/19): Barklem coefficients through memory from SIR to RH.
c 1.11 (11/11/19): RH new code, SIR depurated.
c 1.12 (12/12/19): RH updated.
c 2.03 (03/03/20): Change in amp2.f (azimuth threshold modified).
c 2.04 (04/04/20): Avoid disk writing, RH updated.
c 2.05 (05/05/20): Output on screen centralized.
c 2.06 (06/06/20): New integration methods.
c 2.07 (07/07/20): Kurucz data through memory from SIR to RH.
c 2.08 (08/08/20): Abundance values from SIR to RH.
c 2.09 (09/09/20): Source function re-evaluated.
c 2.10 (10/10/20): Wavetable and some keywords through memory from SIR to RH.
c 2.11 (11/11/20): Atmospheric model and magnetic field through memory.
c 3.03 (03/03/21): Some corrections about weights, chi**2 not getting better,
c                  and allowing not even-spaced models in optical depth.
c 3.05 (05/05/21): Some bugs corrected. Input and output are LoS models.
c                  New continuum evaluation. New option without normalization.
c 3.06 (06/06/21): Han's release memory. DCs in stimulated emission.
c 3.11 (11/11/21): JJ coupling. PARAMETER adjusted, align common data.
c 4.03 (03/03/22): Some includes added and gfortran flags changed.
c 4.04 (04/04/22): Makefiles for M1. Some Bezier routines changed.
c 4.12 (12/12/22): Option (in dtrol file) to write errors and H populations.
c 5.03 (03/03/23): Option (in dtrol file) to write models in RH format.
c                  Option (in dtrol file) to use numerical functions.
c 5.06 (06/06/23): Numerical RF implementation. Change in the chi**2 values
c                  written in file. Offset in velocity (see PARAMETER). New
c                  uncertainties determination.
c
c_____________________________________________________________________________
c
c     DIMENSIONES
      implicit real*4 (a-h,o-z)
      include 'PARAMETER'

      parameter (kt16=16*kt+5)
      parameter (kl4=4*kl,kld4=4*kld)
      parameter (kldt=kld*kn,kldt4=4*kldt)
      parameter (mmax=16*kn+2,kldn=mmax*kld4)
      parameter (kn16=16*kn+4,kt11=11*kt+2,kt12=11*kt+3)

c     PARA LA MALLA
      integer nlin(kl),nlins(kl4),npas(kl),npass(kl4),npos(kld4),nposi(kld),indice(kl)
      integer nlincheck(kl),imaya,ici,iauto
      integer ist(4),nble(kl),nlinsn(kl),il
      real*4  dlamda0(kl),dlamda(kld),pist(4)
      real*4  dlamdas(kld4),cth,vx,vy
      integer m(18),HEon,numerical(18),numericalor(18)
      integer icanal,ican,ncifno,ncifnoOUT,ixt1,ixt2

c     PARA LA ATMOSFERA
      integer icalerr,ncontpg,ipgmag
      integer ntau,ierror,nrays
      integer startj,NEW_J,OLD_J
      real*4  atmos(kt16),atmosr(kt16),atmosrlin(kt16),atmosr_old(kt16)
      real*4  atmosout(kt16),atmoslin(kt16),tau(kt),atmosb(kt16),errores(kt16)
      real*4  atmos1err(8*kt),atmos2err(8*kt)
      real*4  pg1(kt),z1(kt),ro1(kt),T1(kt),Pe1(kt)
      real*4  pg2(kt),z2(kt),ro2(kt),T2(kt),Pe2(kt)
      real*4  pg1b(kt),z1b(kt),ro1b(kt),deglat,deglon
      real*4  pg2b(kt),z2b(kt),ro2b(kt)
      real*4  pg01,pg02,ro01,ro02,alog10mu
      real*4  hydro1(6,kt),hydro2(6,kt),hydro1_conv(6,kt),hydro2_conv(6,kt)

c     PARA EL PERFIL
      real*4  stok(kld4),sig(kld4),stray(kld4),sigd(kld4) !sigd es copia de sig
      real*4  sigescala(kld4),deltastokes_i(kld4)
      real*4  pesos(kld4),ymodobs_RH(kld4)
      real*4  scal(kld4),smax(kl,4),vic(kl),perfil(kld4)
      integer npasobs(kl),numberfrec(4),ntl,icontinuo_normz
      real*4  dlamdaobs(kld),dlamdacheck(kld)

c     PARA LA INVERSION
      integer mreadi2,mreadi3,mreadi4,meves,iratio(4)
      real*4  covar(mfitmax,mfitmax),alpha(mfitmax,mfitmax),beta(mfitmax)
      real*4  mreadr2,mreadr3,mreadr4

c     PARA LOS COEFICIENTES DE ALEJAMIENTO
      integer iRH1,iRH2
      integer iRH1_store(101),iRH2_store(101)
      real*4  rnlte_th

c     PARA LOS ERRORES
      real*4  maximo
      real*4  tup(kt),peup(kt),pgup(kt),zup(kt),roup(kt)
      real*4  tdown(kt),pedown(kt),pgdown(kt),zdown(kt),rodown(kt)
      real*4  errz1(kt),errpg1(kt),errro1(kt)
      real*4  errz2(kt),errpg2(kt),errro2(kt)    

c     PARA EL FICHERO LOG
      integer nguvec(500),posicionciclo(20),mvec(500,18)   !,mvecmax(18)
      real*4  addvec(500),snchivec(500),chprintvec(500),amac1vec(500),amac2vec(500)
      real*4  amic1vec(500),amic2vec(500),fill2vec(500),contrastevec(500)
      real*4  porcienvector(500)
      real*4  snpr(20),chisnpr(20)

c     CADENAS
      character cha1*2,cha2*2,chag*1,chamod*4,chaerr*4,chalog*4,ca(18)*3
      character chaper*4
      character chaweight*4,chahyd*4,chaRH*4,chaRHmag*4  
      character chasnychi*4
      character*100 uveobs,uveout
      character*100 modelout00,modelout1,modelout2
      character*100 modelin1,modelin2,mod1inic,mod2inic
      character*100 modelin2ini,coorfile,extensioncoor,modelout1b
      character*100 mod_hyd_in1,mod_hyd_in2,mod_hyd_out1,mod_hyd_out2
      character*100 malla,control,controlb,fcontrol,modelerr1,modelerr2,difusa
      character*100 pesos_file
      character*100 snychi
      character*100 filtro,psfper,extensionfiltro,extensionperfil
      character*100 mreadc2,mreadc4
      character*100 nomlineas,fichabun
      character*100 mas_or_tau
      character*100 linlarga
      character*70  cartel

c     RH
      character*100 RH_model,RH_magneticfield,label_ID_model

c     10/05/19 epm: To measure time.
      integer*8 cpu1, cpu2, wall1, wall2

c     21/06/19 epm: Command line.
      character*100 arg
      integer*4 flagv, flagquiet, flagtime, flagnolog

c     05/05/20 epm: Error messages.
      character*20  msg1,msg2
      character*400 msg

c     COMUNES
      common/responde2/ist,ntau,ntl,nlin,npas,nble
      common/ldeo/dlamda,dlamda0
      common/smalla/ntls,nlins,npass
      common/smalla1/dlamdas
      common/atmosfera/atmos
      common/atmosferaout/atmosout
      common/iamplioold/iamplioold
      common/tol/tol
      common/filfac/fill
      common/alamda0/alamda0
      common/uvesalida/scal,ymodobs_RH
      common/observaciones/stok,sigd
      common/mu/cth  !para equisub en amp22
      common/sigrealchi/sig,chireal,sumsq
      common/cambio/ncambiono,ncambioprec  !para splines
      common/ifiltro/ifiltro
      common/filtro/filtro
      common/contraste/contr
      common/repeticion/factorrep  !factor de diag para evitar repeticiones
      common/nciclos/nciclos       !para marquardt2, leeuve3, leemodi222
      common/posiciones/npos
      common/calerr/icalerr  !si calculo errores=1 else =0 (para amp2)
      common/difusa/stray
c     El common siguiente es para pasarle a marquardt2 los indices iniciales y 
c     finales de las perturbaciones a gamma y fi aditivas.
c     common/ifies/iga1,ifi11,iga2,ifi22 
      common/ifies/ipa1,ipa11,iga1,ifi11,ipa2,ipa22,iga2,ifi22
      common/ieliminofrec/ielimino  !para el print de la S/N en marquardt2
      common/ivez/ivez
      common/ficlineas/nomlineas  !para leelineas
      common/fichabun/fichabun    !para leeabun
      common/numberfrec/numberfrec,iratio  !para marqcoef2
      common/nommod/modelin1,modelin2      !para escribeFR
      common/atmosr/atmosr                 !para escribeFR
      common/tau/tau
      common/iautomatico/iauto  !se pasa a fperfil2
      common/iprimeravez/iprimeravez  !OJO: se cambia en fperfil2
      common/nspl_lin/nspl_lin !para seleccionar interp.lineal=1 o splines=0
                               !o penalty (sol. regularizada =3 o 2, 5 o 4)
      common/contornopg/ncontpg,pg01,pg02
      common/contornoro/ro01,ro02
      common/primerchi/snn,chisn
      common/zetas/pg1,z1,ro1,pg2,z2,ro2  !para el calculo de z,pg y ro en amp22 y amp22err
      common/pgmag/ipgmag
c     17/06/19 brc: Common repetido con valores absolutos.
c     common/mu/cthabs
      common/anguloheliocent/xmu  !for equisubmu and related routines
                                  !to evaluate equ. hidr. in the z direction,
                                  !not along the LOS
      common/alog10mu/alog10mu    !alog10mu=alog10(xmu) from desire
      common/thresholdNLTE/rnlte_th
      common/RHnames/label_ID_model,RH_model,RH_magneticfield
      common/iRH/iRH1,iRH2  !eq 1 we call RH, 0 don't call RH
      common/HEonoff/HEon  !HEon=1 HE, HEon=0 Non HE
      common/istatus12/istatusdep1,istatusdep2
      common/variablePSF/psfper
      common/numerical/numerical
      common/latlong/deglat,deglon
      common/integrationmethod/intemethod  !0=herm_int,1=herm,2=bzr3,3=bzr3log 
      common/scalemodelRH/imassortau  !integer    0=logtau, 1=mass column para departures
      common/hydroge_populations/hydro1,hydro2  !float (6,kt),H populations para departures
      common/new_evaluation/new_evaluation  !to nodos (at blends2) to evaluate uncertainties if new_evaluation=2
      common/atmoserr/atmos1err,atmos2err  !from fperfil2 for error evaluation
      common/signaltonoise/sn  !to evaluate_errors
      common/iwrtRHformat/iwriteRHformat  !=1 --> writes models on disk

c     05/05/20 epm: Common to save command line flags.
      common/commandline/flagv,flagquiet
      
c     05/05/21 brc: Common to set on/off the continuum normalization.
      common/icontnormalz/icontinuo_normz  !icontinuo_normz=0 input/output profiles are in SI units
c     23/05/23 brc: Common to define chi**2 scale (calculated  in marqcoef2).
      common/escalachi/escalachi,sigescala,chinew
c     06/06/23 brc: Common to know if there is a difference between observed and synthesized lines wave grid
c     common/iFP/iFP  !=1 in case there is a diference between observed and synthesized lines wave grid
      common/deltastokes_i/deltastokes_i
      
      data NEW_J/3/, OLD_J/4/  !see enum solution in "rh.h"

c_____________________________________________________________________________

c     10/05/19 epm: Initial time.
      call chrono(cpu1, wall1)

c     21/06/19 epm: Read the command line.
      fcontrol  = ''
      flagv     = 0
      flagquiet = 0
      flagtime  = 0
      flagnolog = 0
      i = 1
      do
         call get_command_argument(i, arg) !0 = executable, 1 = first argument
         if (len_trim(arg) .eq. 0) exit    ! no more arguments, leave the loop
         if (trim(arg).eq.'-h' .or. trim(arg).eq.'-help') then
            print*,''
            print*,'desire [flags] <control-file>'
            print*,''
            print*,'-h|help    : print this help message'
            print*,'-n|nolog   : do not create a log file'
            print*,'-q|quiet   : do not show SIR messages on console'
            print*,'-t|time    : show CPU time when quiet option'
            print*,'-v|verbose : show RH messages on console'
c           print*,'-trace=n   : write the program trace (n=depth [1-3])'
            print*,''
            stop
         else if (trim(arg).eq.'-n' .or. trim(arg).eq.'-nolog') then
            flagnolog = 1
         else if (trim(arg).eq.'-q' .or. trim(arg).eq.'-quiet') then
            flagquiet = 1
         else if (trim(arg).eq.'-t' .or. trim(arg).eq.'-time') then
            flagtime = 1
         else if (trim(arg).eq.'-v' .or. trim(arg).eq.'-verbose') then
            flagv = 1
c        else if (arg(1:7).eq.'-trace=') then
c           read(arg(8:len_trim(arg)),'(i3)',err=665)flagtrace
         else if (arg(1:1).eq.'-') then
            print*,'desire: invalid option ',trim(arg)
            print*,"Try 'desire -help' for more information"
            stop
         else
            fcontrol = arg
         end if
         i = i + 1
      end do
c     04/04/20 epm: Pass some command line flags to RH.
      call sirflags(flagv)

      call error(KLINE,'','')
      call error(KLINE,'',' __________________________________________________________ ')
      call error(KLINE,'','|                                                          |')
      CARTEL=             '|            DeSIRe version 5.06  (06/Jun/2023)            |'
      call error(KLINE,'',CARTEL)
      call error(KLINE,'','|__________________________________________________________|')
      call error(KLINE,'','')

      if (len_trim(fcontrol) .eq. 0) then
         write(*,*) ''
         write(*,'(a,$)') ' Control file: '
         read(*,'(a)') fcontrol
      end if

c     Inicializa algunas variables
c     ntau : numero de puntos en tau
c     tau  : log10(tau)
c     ntl  : numero total de lineas
c     nlin : indice de cada linea
c     npas : numero de puntos en cada linea
c     dlamda:cada delta de l.d.o. en ma

      cotaminima=-0.99  !valor minimo de cualquier p. Stokes
                        !cualquier Stokes menor de cotaminima tendra una sigma=1.e15

      vx=0.0            !vx=1.95e5 velocidad rotacion fotosfera solar ecuador

      icanal=24         !canal del fichero log
      ican=23           !canal de lectura del fichero de control

c     Se define el fichero log
      chalog='.log'
      chasnychi='.chi'
      controlb=fcontrol

      call quitaex(controlb)
      call concatena(controlb,chalog,control)
      call concatena(controlb,chasnychi,snychi)

      linlarga='------------------------------------------------------------------------------'

c     Se leen las condiciones de ejecucion

c     do i=1,18
c        mvecmax(i)=0  !para los errores
c     end do
      iwriteRHformat=0
      nciclos=1
      ici=0
      ll=0
      ilast=0  !contador #iteracion que se salta por < umbral
      rnlte_thOLD=11.
      new_evaluation=1

      do while (ici.lt.nciclos)  !main loop............................

         ivez=0
         ici=ici+1
         iprimeravez=0
         ncambiono=0    !cambio el numero de nodos
         ncambioprec=0  !cambio la precision
         if(ici .gt. 1)rnlte_thOLD=rnlte_th
         open(ican,file=fcontrol,status='old',err=666)

         rciclos   =mreadr3(ican,ici,1.)     !numero de ciclos (repeticiones)
         uveobs    =mreadc2(ican,ici)        !fichero de entrada perfiles
         difusa    =mreadc2(ican,ici)        !fichero entrada  luz difusa
         filtro    =mreadc2(ican,ici)        !fichero con transformada de PSF
         malla     =mreadc2(ican,ici)        !malla,(npas=n. puntos por linea)
         nomlineas =mreadc2(ican,ici)        !fichero con parametros atomicos
         fichabun  =mreadc2(ican,ici)        !fichero con abundancias
         modelin1  =mreadc2(ican,ici)        !nombre del modelo 1 anterior 
         modelin2  =mreadc2(ican,ici)        !nombre del modelo 2 anterior
         coorfile  =mreadc2(ican,ici)        !nombre del fichero con coordenadas X/Rsun, Y/Rsun, Vrot(cm/s)
         pist(1)   =mreadr3(ican,ici,1.)     !i (peso de i en ajuste)
         pist(2)   =mreadr3(ican,ici,1.)     !q
         pist(3)   =mreadr3(ican,ici,1.)     !u
         pist(4)   =mreadr3(ican,ici,1.)     !v
         iauto     =mreadi3(ican,ici,0)      !seleccion automatica de nodos
         m(1)      =mreadi3(ican,ici,0)      !nodos en t1
         m(2)      =mreadi3(ican,ici,0)      !nodos en pe1
         m(3)      =mreadi3(ican,ici,0)      !nodos en micro1
         m(4)      =mreadi3(ican,ici,0)      !   "   "    h1
         m(5)      =mreadi3(ican,ici,0)      !   "   "    vz1
         m(6)      =mreadi3(ican,ici,0)      !   "   "    gamma1
         m(7)      =mreadi3(ican,ici,0)      !   "   "    fi1
         m(8)      =mreadi3(ican,ici,0)      !   "   "    macro1
         m(9)      =mreadi3(ican,ici,0)      !   "   "    t2
         m(10)     =mreadi3(ican,ici,0)      !   "   "    pe2
         m(11)     =mreadi3(ican,ici,0)      !   "   "    micro2
         m(12)     =mreadi3(ican,ici,0)      !   "   "    h2
         m(13)     =mreadi3(ican,ici,0)      !   "   "    vz2
         m(14)     =mreadi3(ican,ici,0)      !   "   "    gamma2
         m(15)     =mreadi3(ican,ici,0)      !   "   "    fi2
         m(16)     =mreadi3(ican,ici,0)      !   "   "    macro2
         m(17)     =mreadi3(ican,ici,0)      !   "   "    f.f.2
         m(18)     =mreadi3(ican,ici,0)      !luz difusa
         sn        =mreadr3(ican,ici,1000.)  !segnal ruido estimada para i
         contr     =mreadr3(ican,ici,-1.)    !contraste Ic1/Ic2
         tol       =mreadr3(ican,ici,1.e-4)  !tolerancia inversion
         alamda0   =mreadr3(ican,ici,1.e-3)  !factor diagonal inicial
         nspl_lin  =mreadi3(ican,ici,0)      !0=splines, 1=lineal (regularizada splines(2 o 4) o lineal(3 o 5)
         pg01      =mreadr3(ican,ici,0.)     !condi. cont. pg atm. 1/negative value= density
         pg02      =mreadr3(ican,ici,0.)     !condi. cont. pg atm. 2/negative value= density
         rnlte_th  =mreadr3(ican,ici,10.)    !NLTE thresh. (0=NLTE,>10 LTE)10 (LTEdefault
         nrays     =mreadi3(ican,ici,3)      !number of rays
         intemethod=mreadi3(ican,ici,0)      !integration method 0=herm_int,1=herm,2=bzr3,3=bzr3log
         new_evaluation=mreadi4(ican,ici,1)  !optional new evaluation and output writting
                                             !0= no output, no new evaluation (in last cycle 0=1)
                                             !1= writting output model + profile (default)
                                             !2= writting output model + profile + errors
                                             !3= writting output model + profile + errors + H population
                                             !4= 3 option + model in RH format 
         mas_or_tau=mreadc4(ican,ici,'m')    !optional RH models in mass column or logtau 

         imassortau=0  !RH reads tau models (when imassortau=1 RH reads mass column models)
         if(mas_or_tau(1:1) .eq. 'm' .or. mas_or_tau(1:1) .eq. 'M')imassortau=1

         if(nrays.le.0)nrays=3  !default value
         nciclos=int(rciclos)
         if(new_evaluation.eq.0 .and. ici.eq.nciclos)new_evaluation=1
         umbralchi=(rciclos-nciclos)*10.
c        inew=1  !in RH: No using old values of background op & J. In SIR evaluating .per after every cycle 
c        if(rnlte_th.lt.0)then
c           inew=0  !in RH: Using old values of background op & J. In SIR NO evaluating .per after every cycle 
c           rnlte_th=-1*rnlte_th
c        end if
         maxnumnodos=0
         numcalculation=0
         do inum=1,18
            numerical(inum)=0  !RF for parameter i numerically evaluated in NLTE
            if(m(inum).ge.1000)then
               m(inum)=m(inum)-1000
               numerical(inum)=1
               numcalculation=1
            endif
            if(m(inum).le.-1000)then
               m(inum)=m(inum)+1000
c              numerical(inum)=1  !in this case we are NOT evaluating RF because the perturbation i
c                                 !taken equal to the other atmosphere
            endif
            if(m(inum).gt.maxnumnodos)maxnumnodos=m(inum)
         end do
         if(maxnumnodos.gt.kn)then
            write(msg1,*) kn
            call error(KSTOP,'desire','The number of nodes in some variable'
     &      //         ' is larger than kn = '//trim(adjustl(msg1))//'\n'
     &      //         ' Decrease the number of nodes or change the PARAMETER'
     &      //         ' file')
         end if

         call error(KLINE,'','')
         call extension3(coorfile,extensioncoor,n)
         if(extensioncoor(1:n).eq.'coor')then
            call coor(coorfile,deglat,deglon,vx,vy)
            write(msg,*)'Vx(cm/s) =',Vx,' Vy(cm/s) =',Vy
            call error(KLITE,'',msg)
         else if(extensioncoor(1:n).eq.'coorP')then
            call coorP(coorfile,deglat,deglon,vx,vy)
            write(msg,*)'Vx(cm/s) =',Vx,' Vy(cm/s) =',Vy
            call error(KLITE,'',msg)
         else
            call error(KSTOP,'desire','Coordinates file must have a'
     &      //         ' .coor or a .coorP extension\n'
     &      //         ' .coor  contains Xc/Rsun Yc/Rsun Vx[cm/s] Vy[cm/s]\n'
     &      //         ' .coorP contains Xc/Rsun Yc/Rsun P B0')
         endif
         call error(KPARA,'',linlarga)

         HEon=1  ! HEon=1 means HE, HEon=0 meaning OUT of HE
         if(abs(m(2))+abs(m(10)).ne.0 .and. nciclos.gt.0)HEon=0

         imodel2=0
         if(modelin2.ne.'' .or. modelin2.ne.' ')imodel2=1

c        pi=3.14159265
c        cth=cos(deglon*pi/180.)*cos(deglat*pi/180.)
c        xmu=abs(cth)
c        if(xmu .gt. 1.00)xmu=1.00
c        alog10mu=alog10(xmu)
         cth=1.0     !to work in LRF   
         cthabs=1.0  !to work in LRF               
         xmu=1.0     !to work in LRF
         alog10mu=0.
c        04/04/20 epm: El coseno del angulo heliocentrico pasa por un common
c        a departures() y alli pasa por argumento de entrada a solveray().
c        open(1,file='ray.input')
c        write(1,*)xmu
c        write(1,*)0
c        close(1)

         ipgmag=0  !we do not consider magnetic pressure terms

c        ---------------------------------------------------------------

         if(rnlte_th .lt. 0.)rnlte_th=0.

         ncontpg=0
         ro01=0.
         ro02=0.
         if(pg01.gt.0. .or. pg02.gt.0. .and. abs(m(2))+abs(m(10)).eq.0)then
            if(pg02.lt.0.)pg02=pg01
            if(pg01.lt.0.)pg01=pg02
            ncontpg=1
         end if
         if(pg01.lt.0. .or. pg02.lt.0. .and. abs(m(2))+abs(m(10)).eq.0)then
            if(pg02.lt.0.)pg02=pg01
            if(pg01.lt.0.)pg01=pg02
            pg01=abs(pg01)
            pg02=abs(pg02)
            ro01=pg01
            ro02=pg02
            ncontpg=-1
         end if

c        Cadenas para el output
         ca(1) =' t1'
         ca(2) =' p1'
         ca(3) =' m1'
         ca(4) =' h1'
         ca(5) =' v1'
         ca(6) =' g1'
         ca(7) =' f1'
         ca(8) =' M1'
         ca(9) =' t2'
         ca(10)=' p2'
         ca(11)=' m2'
         ca(12)=' h2'
         ca(13)=' v2'
         ca(14)=' g2'
         ca(15)=' f2'
         ca(16)=' M2'
         ca(17)=' ff'
         ca(18)=' st'

         if(nciclos.eq.0)then  !para sintesis uso los 4 parametros (LRB)
            do i=1,4
               iratio(i)=0  !iratio(Q)=1 -> se invierte Q/I, etc
               ist(i)=1
               if(pist(i).lt.0)iratio(i)=1
               pist(i)=1
            enddo
            iratio(1)=0  !I no se puede sintetizar como I/I
         else
            do i=1,4
               ist(i)=0
               iratio(i)=0  !iratio(Q)=1 -> se invierte Q/I, etc
               if(pist(i).lt.0)iratio(i)=1
               pist(i)=abs(pist(i))
               if(pist(i).gt.0)ist(i)=1
            end do
         endif
         if(iratio(2).eq.1.or.iratio(3).eq.1.or.iratio(4).eq.1)iratio(1)=0

c        Definimos el fichero con los parametros atomicos
         ilineas=0
         ilineas=meves(nomlineas,100)
         if(ilineas.eq.0)then
            call error(KSTOP,'desire','Atomic parameters file is not defined'
     &      //         ' in the control file')
         endif

c        Definimos el fichero con las abundancias
         iabun=0
         iabun=meves(fichabun,100)
         if(iabun.eq.0)then
            call error(KSTOP,'desire','Abundance file is not defined in the'
     &      //         ' control file')
         endif

c        Definimos los nombres de los ficheros de salida
         chamod='.mod'
         chaerr='.err'
         chaper='.per'
         chaweight='.wht'
         chahyd='.hyd'
         chaRH='.atm'
         chaRHmag='.mag'

         pesos_file=uveobs
         call quitaex3(pesos_file,ixt3)
         call concatena(pesos_file,chaweight,pesos_file)
         
c        Determining if the observed file ends on .per or .spc
         call extension3(uveobs,extensionperfil,n)
         icontinuo_normz=1  !in this case input/output profiles are normalized to HSRA QS continuum
         if(extensionperfil(1:n).ne. 'per')then
            icontinuo_normz=0 !in this case input/output profiles are in SI units
            chaper='.'//extensionperfil(1:3)
         endif         
         
         ncifno=1
         ncifnoOUT=1

         chag='_'
         cha1(1:2)='  '
         cha2(1:2)='  '
         modelout1=modelin1

         if(imassortau .eq. 0)then 
            mod_hyd_in1=modelout1
            call quitaex3(mod_hyd_in1,ix0h)
         end if
         
         call quitaex2(modelout1,ixt1) 
         if(ici+ixt1-1.ge.9)ncifnoOUT=2
         if(ici+ixt1-1.ge.10)ncifno=2

         if(imodel2 .eq. 1)then
            modelout2=modelin2 
            modelin2ini=modelin2
            if(imassortau .eq. 0)then 
               mod_hyd_in2=modelout2
               call quitaex3(mod_hyd_in2,ix2h)
            end if   
            call quitaex2(modelout2,ixt2)
         end if  
         
         if(ici.eq.1.or.ici.eq.0)then
            mod1inic=modelout1
            call chachi(1+ixt1,cha1,ncifno)
            if(imassortau .eq. 0)call concatena(mod_hyd_in1,chahyd,mod_hyd_in1) 
            if(imodel2 .eq. 1)then
               mod2inic=modelout2
               call chachi(1+ixt2,cha2,ncifno)
               if(imassortau .eq. 0)call concatena(mod_hyd_in2,chahyd,mod_hyd_in2) 
            end if   
         else
            if(modelout1.eq.mod1inic)then
               call chachi(ici-1+ixt1,cha1,ncifno)
               call concatena(modelout1,chag,modelin1)
               call concatena(modelin1,cha1,modelin1)
               if(imassortau .eq. 0)call concatena(modelin1,chahyd,mod_hyd_in1) 
               call concatena(modelin1,chamod,modelin1)
               call chachi(ici+ixt1,cha1,ncifnoOUT)
            else
               call chachi(1+ixt1,cha1,ncifno)
            endif
            if(imodel2.eq.1)then
               if(modelout2.eq.mod2inic )then
                  call chachi(ici-1+ixt2,cha2,ncifno)
                  call concatena(modelout2,chag,modelin2)
                  call concatena(modelin2,cha2,modelin2)
                  if(imassortau .eq. 0)call concatena(modelin2,chahyd,mod_hyd_in2)                   
                  call concatena(modelin2,chamod,modelin2)
                  call chachi(ici+ixt2,cha2,ncifnoOUT)
               else
                  call chachi(1+ixt2,cha2,ncifno)               
               endif
            endif
         endif

         modelout00=modelout1  !requiered to build the output name (RH format) when new_evaluation=4 and nciclos=0
         call concatena(modelout1,chag,modelout1) !"guess_"
         call concatena(modelout1,cha1,modelout1) !"guess_1"
         
         mod_hyd_out1=modelout1  
         call concatena(mod_hyd_out1,chahyd,mod_hyd_out1) 
         if(imodel2.eq.1 .and. imassortau.eq.0)then 
            mod_hyd_out2=modelout2  ! "guess" without ".mod'
            call concatena(mod_hyd_out2,chahyd,mod_hyd_out2)         
         end if           
         call concatena(modelout1,chaper,uveout)
         call concatena(modelout1,chaerr,modelerr1)
         modelout1b=modelout1
         call concatena(modelout1,chamod,modelout1) 
         
         if(imodel2.eq.1)then
            call concatena(modelout2,chag,modelout2)
            call concatena(modelout2,cha2,modelout2)
            call concatena(modelout2,chaerr,modelerr2)
            call concatena(modelout2,chamod,modelout2)
         endif
         if(contr.lt.0)then
            contr=-1.e6
         else
            contr=(contr-1.)/(contr+1.)
         end if

         ifiltro=0
         ifiltro=meves(filtro,100)

         if(ifiltro.eq.1)then
            call extension3(filtro,extensionfiltro,n)
            if(extensionfiltro(1:n) .eq. extensionperfil(1:n))then
               psfper=filtro
               ifiltro=2
            endif
         endif

         istray=meves(difusa,100)

c        -----------------------------------------------------
c        Leemos la malla, el perfil, el modelo y la luz difusa
c        -----------------------------------------------------
         imodel1=meves(modelin1,100)
         imodel2=meves(modelin2,100)
         if(imodel1.eq.1.and.imodel2.eq.1)ierror=0  !existen ambos
         if(imodel1.eq.0.and.imodel2.eq.1)ierror=1  !no existe mod1 pero si mod2
         if(imodel1.eq.1.and.imodel2.eq.0)ierror=2  !si existe mod1 pero no mod2
         if(imodel1.eq.0.and.imodel2.eq.0)ierror=3  !no existe ninguno de los 2 modelos

         if(ierror.eq.1 .or. ierror.eq.3)then
            call error(KSTOP,'desire','Model 1 does not exist')
         endif

         if(rnlte_th .lt. 10)then
            if(rnlte_thOLD .ge. 10)then
               startj=NEW_J
            else
               startj=OLD_J
            end if
            label_ID_model=modelin1
c           11/11/20 epm: El modelo de atmosferas y el campo pasan por memoria.
c           call read_keyword_input_RH('ATMOS_FILE',RH_model)
c           call read_keyword_input_RH('STOKES_INPUT',RH_magneticfield)
c           10/10/20 epm: Pass some keywords to RH.
            call sirkeywords(nrays,startj)
         end if 

         if(ici.eq.1)call leeabun(0)
         if(ici.eq.1)call leemallab(malla,ntl,nli,nlin,npas,dlamda,nble)

c        #1) Transformation of model format from multi-RH to SIR format
c        if(multiformat1.ne.0 .and. ici.lt.2)then
c        ... 11/11/20 epm: Suprimo todo lo relativo a multiformat ...
c        end if  !#1) end if(multiformat1.ne.0 .and. ici.lt.2)

         if(ici.eq.1)call leemodi222(0,modelin1,modelin2,atmos,ntau,ierror,
     &                               z1,pg1,ro1,z2,pg2,ro2)

         do i=1,ntau
            tau(i)=atmos(i)
            T1(i)=atmos(ntau+i)
            Pe1(i)=atmos(2*ntau+i)
            T2(i)=atmos(9*ntau+2+i)
            Pe2(i)=atmos(10*ntau+2+i)               
         end do
         if(m(2).eq.0 .and. ici.eq.1)then
            if(ncontpg.eq.0) call equisubmu(ntau,tau,T1,Pe1,pg1,z1,ro1)
            if(ncontpg.ne.0) call equisubmu_cont(ntau,tau,T1,Pe1,pg01,pg1,z1,ro1)
            do i=1,ntau
               atmos(2*ntau+i)=Pe1(i)
            end do
         end if
         if(imodel2.eq.1 .and. m(10).eq.0 .and. ici.eq.1)then
            if(ncontpg.eq.0)call equisubmu(ntau,tau,T2,Pe2,pg2,z2,ro2)
            if(ncontpg.ne.0)call equisubmu_cont(ntau,tau,T2,Pe2,pg02,pg2,z2,ro2)
            do i=1,ntau
               atmos(10*ntau+2+i)=Pe2(i)
            end do
         end if

         if(imassortau.eq. 0)then
            call lee_hydro(mod_hyd_in1,ihydro1,ntau,tau,T1,Pe1,hydro1) !reading H populations model1 if exist (ihydro1=1)
            if(imodel2 .eq. 1)call lee_hydro(mod_hyd_in2,ihydro2,ntau,tau,T2,Pe2,hydro2) !reading H populations model2 if exist (ihydro2=1)
         end if

         if(ierror.eq.2)then  !no podemos invertir el modelo 2
            do jj=9,17
               m(jj)=0
            end do
         endif
         if(ici.eq. 1)then
            if(atmos(16*ntau+5).gt.1e-5 .or. m(18).ne.0)then
               if(istray.ne.0)then
                  call leeuveobsindic(difusa,ist,ntl,nliobs,nlin,npasobs,dlamdaobs,nble,stray)

                  ntlcheck=ntl  !para comprobar que difusa tiene las mismas lambdas que perf. obs.
                  do l=1,nliobs
                     dlamdacheck(l)=dlamdaobs(l)
                  end do
                  do l=1,ntl
                     nlincheck(l)=nlin(l)
                  end do

c              Por si no tenemos el fichero con la malla en sintesis:
c              (si lo tenemos, estas variables se machacan luego y no pasa nada)
                  nli=nliobs
                  do i=1,ntl
                     npas(i)=npasobs(i)
                  end do
                  do i=1,nli
                     nposi(i)=i
                     dlamda(i)=dlamdaobs(i)
                  end do

               else
                  if(atmos(16*ntau+5).gt.0.)then
                     call error(KWARN,'desire','As the stray light file does'
     &            //         ' not exist,\n the stray light factor is being'
     &            //         ' changed to zero')
                  end if
                  atmos(16*ntau+5)=0.
                  m(18)=0  !LB (si no casca, porque trato de invertir una variable = cero)
               end if
            end if
            do i=1,16*ntau+5
              atmosb(i)=atmos(i)
            end do
         end if
c        -------------------------
c        Para el calculo de las FR
c        -------------------------
         if(nciclos.eq.0)then
            m(2)=1
         end if
         if(nciclos.lt.0)then  !para calculo de las FR mnodos=ntau
c           if(kn.lt.ntau)then
c              call error(KSTOP,'desire','The number of depth points in the'
c    &         //         ' model is larger than the current kn value\n'
c    &         //         ' Decrease the number of depth points or change'
c    &         //         ' the PARAMETER file')
c           end if
            nRF=ntau
            if(kn.lt.ntau)nRF=kn

            if(nciclos.eq.-1)then
               do i=1,15
                  if(m(i).ne.0)m(i)=nRF  !previously here said m(i)=ntau
               enddo
               if(m(1).ne.0 .and. m(2) .eq. 0)m(2)=nRF
               if(m(8).ne.0)m(8)=1  !la macro solo puede tener 1 nodo
            endif
            iauto=0
         end if

c        ---------------------------------------------------------------
         if(ici.eq.1)then
            imaya=meves(malla,100)       !(LRB)

            if(nciclos.lt.1)then
               uveout=uveobs
               if(imaya.eq.1)then  !si tenemos malla
                  call leemallab(malla,ntl,nli,nlin,npas,dlamda,nble)
                  nliobs=nli
                  do i=1,ntl
                     npasobs(i)=npas(i)
                  end do
                  do i=1,nli
                     nposi(i)=i
                     dlamdaobs(i)=dlamda(i)
                  end do
               else
                  if(istray.eq.1)then
                     call error(KWARN,'desire','The wavelength grid is being'
     &               //         ' generated automatically\n from the stray'
     &               //         ' light profile. No blends are considered')
                  else
                     call error(KSTOP,'desire','In synthesis mode, you must'
     &               //         ' specify the wavelength grid')
                  endif
               endif
            else
               if(imaya.eq.0)then
                  call error(KWARN,'desire','The wavelength grid is being'
     &            //         ' generated automatically\n from the observed'
     &            //         ' profiles. No blends are considered')

                  call leeuveobsindic(uveobs,ist,ntl,nliobs,nlin,npasobs,
     &                             dlamdaobs,nble,stok)
                  nli=nliobs
                  do i=1,ntl
                     npas(i)=npasobs(i)
                  end do
                  do i=1,nli
                     nposi(i)=i
                     dlamda(i)=dlamdaobs(i)
                  end do
               else
                  call leeuveobsindic(uveobs,ist,ntl,nliobs,indice,npasobs,
     &                             dlamdaobs,nble,stok)
                  call leemalla2(malla,ntl,nli,nliobs,nlin,npasobs,npas,
     &                        dlamdaobs,dlamda,nble,nposi,indice)  
c                 nli    = number of frequencies for each line in grid
c                 nliobs = number of frequencies for each line in observed profiles
c                 iFP=0 !=1 if Fabry-Perot= non equispaced observed wavelength grid
c                 if(nli .ne. nliobs)then
c                    if(nli .lt. nliobs)then
c                       print*,'desire line 855: the number of frequencies points in the grid is lower than in the observed profile'
c                       stop
c                    end if
c                    iFP=1
c                 endif 
               endif
            endif
         end if
         
         if(atmos(16*ntau+5).gt.0.)then  !comprobaciones varias, esto es si tenemos difusa
            if(ntlcheck.ne.ntl)then
               call error(KSTOP,'desire','The number of lines'
     &         //         ' in the files containing the observed\n'
     &         //         ' and stray light profiles are not equal')
            endif

            do l=1,ntl
               if (l.eq.1)then
                  numm=1
               else
                  numm=numm+nble(l-1)
               endif
               if(nlin(numm).ne.nlincheck(l))then
                  call error(KSTOP,'desire','The order of the lines in the'
     &            //         ' stray light file is not equal\n to that in'
     &            //         ' the file containing the observed profiles')
               endif
            enddo

            do l=1,nli
               if(abs(dlamda(l)-dlamdacheck(l)).gt.1.)then
                  write(msg,*)'Wavelength #',l,':',dlamda(l),dlamdacheck(l)
                  call error(KSTOP,'desire','The wavelengths in the file'
     &            //         ' containing the stray light profile differ\n'
     &            //         ' by more than 1 mA from those generated by the'
     &            //         ' wavelength grid\n'
     &            //         msg)
               endif
            enddo

         endif

c        Print these values on screen.
         write(msg,*)'Number of wavelengths in the wavelength grid  : ',nli
         call error(KLITE,'',msg)
         write(msg,*)'Number of wavelengths in the observed profiles: ',nliobs
         call error(KLITE,'',msg)
         call error(KLINE,'','')

         if(nciclos.ge.1 .or. nciclos.eq.0)then
            if(new_evaluation .ge. 2)then
               if(ierror.eq.0)then
                  call error(KPARA,'',
     &                 'Output models  : '//modelout1(1:20)//' and '//modelout2(1:20)//'\n'
     &            //  ' Uncertainties  : '//modelerr1(1:20)//' and '//modelerr2(1:20)//'\n'
     &            //  ' Output profiles: '//uveout(1:20)//' (2 components)')
               else if(ierror.eq.2)then
                  call error(KPARA,'',
     &                 'Output model   : '//modelout1(1:50)//'\n'
     &            //  ' Uncertainties  : '//modelerr1(1:50)//'\n'
     &            //  ' Output profiles: '//uveout(1:43)//' (1 component)')
               else
                  call error(KSTOP,'desire','Models not specified or unexistent')
               end if
               if(new_evaluation .eq. 4)then
                  if(nciclos .eq. 0)modelout1b=modelout00
                  call concatena(modelout1b,chaRH,RH_model)
                  call concatena(modelout1b,chaRHmag,RH_magneticfield) 
               end if
            else if(new_evaluation .eq. 1 .or. ici .eq. nciclos)then
                if(ierror.eq.0)then
                  call error(KPARA,'',
     &                 'Output models  : '//modelout1(1:20)//' and '//modelout2(1:20)//'\n'
     &            //  ' Output profiles: '//uveout(1:20)//' (2 components)')
               else if(ierror.eq.2)then
                  call error(KPARA,'',
     &                 'Output model   : '//modelout1(1:50)//'\n'
     &            //  ' Output profiles: '//uveout(1:43)//' (1 component)')
               else
                  call error(KSTOP,'desire','Models not specified or unexistent')
               end if           
            end if
         end if

c        ---------------------
c        Sentencias de control
c        ---------------------
         if(ici .eq. 1)then  

            if(ierror.ne.0.and.contr.gt.-1.e-6)then
               contr=-1.e-6
               call error(KWARN,'desire','No contrast calculation can be done'
     &      //         ' without two models')
            end if

            nfrecmax=0
            ntlblends=0
            do i=1,ntl
               nfrecmax=max0(nfrecmax,npas(i))  !numero max. de frec.
               nlinsn(i)=nlin(ntlblends+1)
               do j=1,nble(i)
                  ntlblends=ntlblends+1
               end do
            end do

            if (ntlblends.gt.kl) then
               write(msg,'(a,i3,a,i3,a)') 'The number of lines (',ntlblends,
     &            ') is larger than the current limit kl = ',kl,
     &            '\n Decrease this number or change the PARAMETER file'
               call error(KSTOP,'desire',msg)
            end if
            if (nli.gt.kld) then
               write(msg,'(a,i5,a,i5,a)') 'The number of wavelengths (',nli,
     &            ') is larger than the current limit kld = ',kld,
     &            '\n Decrease the number of wavelengths or change the PARAMETER file'
               call error(KSTOP,'desire',msg)
            end if
            if (nfrecmax.gt.kld) then
               write(msg,'(a,i5,a)') 'There is a line with more wavelengths'
     &         //    ' than the current limit kld = ',kld,'\n Decrease the'
     &         //    ' number of wavelengths or change the PARAMETER file'
               call error(KSTOP,'desire',msg)
            end if
            if (ntau.gt.kt) then
               write(msg,'(a,i4,a)') 'The model has more depth points than the'
     &         //    ' current limit kt = ',kt,'\n Decrease the number of grid'
     &         //    ' points or change the PARAMETER file'
               call error(KSTOP,'desire',msg)
            end if
            
c           05/05/21 brc: New routine for avoiding multiple line readings.
c           Reads all the spectral information and make some calculations.
c           Sends the information via common to blends2, blendscon2, departures
c           and fperfil2. Input via common: responde2
            call spectra()
            if(nciclos .ge. 1. and. icontinuo_normz .eq. 0) then !input profiles in SI
                call renormaliza(ist,ntl,npas,stok)
            end if

c           Definimos nlins,ntls,npass,dlamda0s,dlamdas
            ntls=0
            nb=0
            k4=0
            k4c=0
            do i=1,ntl
               vic(i)=1.  !inicializo i del continuo
               do ii=1,4
                  smax(i,ii)=1./sn
               end do
            end do

            ii=0
            vmax=-2.
            vmin=2.
            nfrecI=0 !number of frequencies in I
            do i=1,4
               do j=1,ist(i)
                  iinli=ii*nli
                  ii=ii+1
                  k3=0
                  k3c=0
                  jj=0
                  do k=1,ntl
                     do jble=1,nble(k)
                        nb=nb+1
                        jj=jj+1
                        nlins(nb)=nlin(jj)
                     end do
                     ntls=ntls+1
                     npass(ntls)=npas(k)
                     do l=1,npasobs(k)
                        k3=k3+1
                        k4=k4+1
                        npos(k4)=nposi(k3)+iinli
                        sk4=stok(k4)
                        pesos(k4)=1.00

                        if(sk4.gt.cotaminima)then
                           if(i.eq.1 .and. nciclos.gt.0)then
                              xmax=-1.
                              nfrecI=nfrecI+1
                              do lll=1,npasobs(k)
                                 if(xmax.le.stok(k4-l+lll))xmax=stok(k4-l+lll)
                              end do
                              vic(k)=xmax  !i continuo
                              sss=abs( vic(k)-sk4 )
                           else
                              sss=abs(sk4)
                           end if
                           if(sss.gt.smax(k,i))smax(k,i)=sss
                        end if

                     end do  !fin del do en frecuencias obs
                     do l=1,npas(k)
                        k3c=k3c+1
                        k4c=k4c+1
                        dlamdas(k4c)=dlamda(k3c)
                     end do  !fin del do en frecuencias malla

                  end do  !fin del do en lineas
               end do  !fin del do en j
            end do  !fin del do en i

            if(nciclos.gt.0)call leepesos(pesos_file,ist,pesos)

            maximo=smax(1,1)
            do k=1,ntl
               if(smax(k,1).gt.maximo)maximo=smax(k,1)
            enddo

            do k=1,ntl
               if(smax(k,1).lt..01 .and. nciclos.gt.0)then
                  write(msg,'(a,i2)') 'Weight is being changed for the line ',k
                  call error(KWARN,'desire',msg)
                  smax(k,1)=maximo
               endif
            enddo
            
         end if !end if(ici .eq. 1)then

         nfrecs=k4

c        Para controlar el maximo numero de nodos en modo automatico ("*" en m(i))
         if(m(8).gt.1)m(8)=1
         if(m(16).gt.1)m(16)=1
         if(m(17).gt.1)m(17)=1

         if(iauto.eq.1)then
            n_1000=0
            n_det_min=0
            do i=1,18
               if(m(i).eq.1000)then
                  n_1000=n_1000+1
               else
                  if(m(i).gt.0)n_det_min=n_det_min+m(i)
               end if
            end do
            ndd1=min(mfitmax,nfrecs)
            nmax1=nint((ndd1-n_det_min)/float((n_1000+2)))

            if(m(1).eq.1000)m(1)=2*nmax1  !temperatura (doble num de nodos)
            do i=2,7
               if(m(i).eq.1000)m(i)=nmax1
            end do
            if(m(9).eq.1000)m(9)=2*nmax1  !temperatura (doble num de nodos)
            do i=10,15
               if(m(i).eq.1000)m(i)=nmax1
            end do
            kv=0
         else
            do i=1,18
               if(m(i).eq.1000)then
                  call error(KSTOP,'desire',"'*' is not allowed in non"
     &            //         " automatic mode")
               end if
            end do
         end if

         if(nciclos.gt.0)call nodos2aut(ntau,0,ierror,ici,m)  !en la version LBR estaba sin comentar

         mfit=0
         do i=1,18
            if(m(i).gt.0)mfit=mfit+m(i)
         end do

         if(mfit.gt.mfitmax .and. nciclos .gt. 0)then
            write(msg1,*) mfit
            write(msg2,*) mfitmax
            call error(KSTOP,'desire','The number of free parameters'
     &      //         ' ('//trim(adjustl(msg1))//') is larger than'
     &      //         ' mfitmax = '//trim(adjustl(msg2))//'\n'
     &      //         ' Decrease the number of free parameters or change'
     &      //         ' the PARAMETER file')
         end if

c        Escribe en atmosr los valores de atmos en los nodos
         if(nciclos.ne.0)call comprime2(ntau,m,atmos,atmosr)
     
         nfree=nfrecs-mfit

         write(msg,'(a,i8,1x,i6,1x,i8)')
     &         'Number-of-observ-freqs/fitable-params/free-params:  ',
     &         nfrecs,mfit,nfree
         call error(KPARA,'',msg)

         if(mfit.gt.nfrecs .and. nciclos.gt.0)then
            write(msg,'(a,i2,a,i3,a)')
     &            'The number of free parameters (',mfit,
     &            ')\n is larger than the number of observables (',nfrecs,')'
            call error(KSTOP,'desire',msg)
         end if

         if(ici.eq.1)then
            k40=0
            ntotal4=0
            ielimino=0
            sigmamax=1.e15

            do i=1,4
               numberfrec(i)=0
               do j=1,ist(i)
                  indicei=0
                  sigc=1./(sn*sqrt(pist(i)))
                  do k=1,ntl
                     numberfrec(i)=numberfrec(i)+npasobs(k) !umero de ldo para cada stokes invertido (0 si no se invierte)
                     sigc2=sigc*sqrt(smax(k,i)/vic(k))
                     do l=1,npasobs(k)
                        k40=k40+1
                        indicei=indicei+1
                        if(pesos(k40).gt.1.e-10)pesos_inv=1./pesos(k40)
                        if(pesos(k40).le.1.e-10)pesos_inv=1.e10
                        sig(k40)=sigc2
                        sigescala(k40)=vic(k)/sn
                        deltastokes_i(k40)=sigescala(k40) !inicializado: en marqcoef2 se calcula el maximo entre abs(dy) y este valor a cada frecuencia
                        if(stok(k40).lt.cotaminima)then
                           sigx=sigmamax
                           sig(k40)=sigx
                           sigescala(k40)=sigx
                           ielimino=ielimino+1
                        end if
                        sig(k40)=sig(k40)*pesos_inv
                        sigd(k40)=sig(k40)
                     end do
                     do l=1,npas(k)
                        ntotal4=ntotal4+1
                     end do
                  end do
                  if(contr.gt.-1.e5)then
                     if(i.eq.1)then
                        sig(k40)=4.*sig(k40-1)/sqrt(float(k3)*pist(i))
                     else
                        sig(k40)=sigmamax/sqrt(pist(i))
                     end if
                     sigd(k40)=sig(k40)
                  end if
               end do
            end do

            if(ielimino.gt.0)then
               write(msg,'(a,f5.1,a,i3,a)') 'Being smaller than ',cotaminima,
     &            ', ',ielimino,' data points are not considered for inversion'
               call error(KWARN,'desire',msg)
            endif
         end if
         if(nciclos.gt.0)then
            if(iauto.eq.0)then
               call error(KLINE,'',linlarga)
               write(msg,*) "ite  DE     s/n      chi**2    Mac1   "
     &         //           "Mac2   mic1   mic2   fill2 Ic1/Ic2 stray"
               call error(KLITE,'',msg)
               call error(KLINE,'',linlarga)
            end if
            if(iauto.eq.1)then
               call error(KLINE,'',linlarga)
               write(msg,2000) "ite","DE","s/neqv","chi**2",(ca(jj),jj=1,18)
               call error(KLITE,'',msg)
               call error(KLINE,'',linlarga)
            end if
         end if

         chireal_store=chireal
         sumsq_store=sumsq
         icorrijo=0
               
c        ===============================================================
c                        Empieza el ciclo iterativo
c        ===============================================================

         if(nciclos.eq.0 .and. new_evaluation.eq.4)iwriteRHformat=1
         factorrep=1.
         diag=-1.0
         if(ici.eq.1)chi0=1.e20
         varchi=1.0
         isigo=1
         it=0
         iamplio=0
         ngu=0
         icalerr=0  !=0 pues de momento no vamos a calcular errores
         iRH1=0
         iRH2=0
         if(rnlte_th.lt.10.)IRH1=1
         if(rnlte_th.lt.10. .and. imodel2.eq.1)IRH2=1
         deltaT1_rel=0.
         deltaT2_rel=0.
         iamplioold=0
         
         if(it .eq. 0 .and. ici .eq. 1)then
            chirealgood  = 1.e8
            escalachigood= 1.
            chinewgood   = 1.e8
            sumsqgood    = 1.e8
            snchigood    = 1.e8
            chprintgood  = 1.e8
         end if   
                  
         if(ici.gt.1 .and. chprint.le.umbralchi)then
            isigo=0
            ilast=ilast+1  !contador #iteracion que se salta por < umbral
            write(msg,'(a,e10.3,a,e10.3)')
     &            'Skipping this cycle because chi = ',chprint,
     &            ' < threshold = ',  umbralchi
            call error(KLINE,'',msg)
         else
            if(ici .gt. 1)then
               write(msg,'(a,e10.3,a,e10.3)')
     &               'Doing this new cycle because chi = ',chprint,
     &               ' > threshold = ',  umbralchi
               call error(KLINE,'',msg)
            end if  
         end if

         do while(isigo.eq.1)
            it=it+1
            if(it.gt.1 .and. new_evaluation.ge.3)then
               call copia_hydro(ntau,hydro1_conv,hydro1,imodel2,hydro2_conv,hydro2)
            end if
            if(it.gt.1 .and. icorrijo.eq.0)then           
               iRH1=0
               iRH2=0
            end if
            icorrijo=0
            if(it.gt.2)then
               if(iamplio.eq.0 .and. iRH1_store(it-1).eq.1
     &            .and. iRH1_store(it-2).eq.0 .and. rnlte_th.lt.10.)IRH1=1
               if(iamplio.eq.0 .and. imodel2.eq.1 .and. iRH2_store(it-1).eq.1
     &            .and. iRH2_store(it-2).eq.0 .and. rnlte_th.lt.10.)IRH2=1
            end if
            if(it.gt.1)then
               if(iamplio.eq.1 .and. deltaT1_rel.gt.rnlte_th
     &            .and. rnlte_th.lt.10.)IRH1=1
               if(iamplio.eq.1 .and. imodel2.eq.1 .and. deltaT2_rel.gt.rnlte_th
     &            .and. rnlte_th.lt.10.)IRH2=1
            end if
            
            iRH1_store(it)=IRH1
            iRH2_store(it)=IRH2

            diag0=diag

            if(nciclos.gt.0)then
               do inod=1,m(1)
                  atmosr_old(inod)=atmosr(inod)
               end do
               mm=m(1)+m(2)+m(3)+m(4)+m(5)+m(6)+m(7)+m(8)+m(9)+m(10)
               do inod=mm+1,mm+m(11)
                  atmosr_old(inod)=atmosr(inod)
               end do
            end if

            nn=8*ntau+2
            call marquardt2(stok,sig,nfrecs,atmosr,m,mfit,
     &                      covar,alpha,chisq,diag,iamplio,beta)

            if(nciclos.gt.0)then
               deltaT1_rel=0.
               do inod=1,m(1)
                  deltaT1_reli=abs((atmosr(inod)-atmosr_old(inod))/atmosr_old(inod))
                  if(deltaT1_reli.gt.deltaT1_rel)deltaT1_rel=deltaT1_reli
               end do
               if(imodel2.eq.1)then
                  deltaT2_reli=0.
                  do inod=mm+1,mm+m(11)
                     deltaT21_reli=abs((atmosr(inod)-atmosr_old(inod))/atmosr_old(inod))
                     if(deltaT2_reli.gt.deltaT2_rel)deltaT2_rel=deltaT2_reli
                  end do
               end if
            end if

            if(nciclos.le.0)then
               if(nciclos.eq.0 .and. rnlte_th .ge.10)
     &            call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,scal)
               if(nciclos.eq.0 .and. rnlte_th .lt.10)
     &            call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,
     &                         ymodobs_RH)
c              10/05/19 epm: Measure time before existing the program.
               call chrono(cpu2, wall2)
               call error(KLINE,'','')
               write(msg,'(a,i8,a,i8)') 'CPU  time (ms) = ', cpu2 - cpu1,
     &                               '\n Wall time (ms) = ', wall2 - wall1
               call error(KPARA,'',msg)

               stop  !end of synthesis..................................
            end if

            nfree=nfrecs-mfit

            if(it.eq.1)diag0=alamda0
            if(iRH1 .eq.1 .or. iRH2.eq.1)chi00=chi0
            varchi=(chi0-chisq)/(chi0+chisq)
            varchiLTE=varchi
            if(iRH1 .eq.1 .or. iRH2.eq.1)varchi=(chi00-chisq)/(chi00+chisq)
            call diagonal(diag0,diag,isigo,iamplio,varchi,it)

            if(isigo. eq. 0 .and. iRH1 .eq. 0 .and. varchi .gt. 1.e-3 .and. numcalculation .eq. 1)then            
               isigo=1
               icorrijo=1
               iRH1=1  
               iRH2=1
            end if

            if(istatusdep1*istatusdep2.ne.0)iamplio=0
            if(istatusdep1*istatusdep2.ne.0)then
               call error(KWARN,'desire','Setting old atmosphere due to'
     &         //         ' non RH convergence')
            end if
            istatusdep1=0
            istatusdep2=0

c           Copiamos la primera sintesis del primer ciclo en perfil.
            if(ici .eq. 1 .and. it.eq.1)then
               if(rnlte_th.lt.10.)then
                  do i=1,k40  !numero de longitudes de onda
                     perfil(i)=ymodobs_RH(i)
                  end do              
               else
                  do i=1,k40  !numero de longitudes de onda
                     perfil(i)=scal(i)
                  end do
               end if 
               chi0=chisq  
               chiw=sumsq/float(nfrecs-ielimino)
               snchi=1./sqrt(chiw)

               add=-20.
               if(diag0.gt.0.d0)add=alog10(diag0)
               chprint=chireal/float(nfree-ielimino)
            end if
            
            if(iamplio.eq.1)then  !converge
            
               chirealgood=chireal
               escalachigood=escalachi
               chinewgood=chinew
               sumsqgood=sumsq
               snchigood=snchi
               chprintgood=chprint
              
               ngu=ngu+1
               chi0=chisq

               chireal_store=chireal
               sumsq_store=sumsq

c              Copiamos la atmosfera de salida (via common) atmosout en atmos
c              para no machacar el common y la duplicamos en atmoslin.
c              La inversion se realiza en la linea de vision, en radianes, pero
c              los modelos de entrada y salida se escriben en disco en grados.
               do i=1,16*ntau+5
                  atmos(i)=atmosout(i)
                  atmosb(i)=atmosout(i)
               end do

               do i=1,ntau
                  pg1b(i)=pg1(i)
                  pg2b(i)=pg2(i)
                  z1b(i)=z1(i)
                  z2b(i)=z2(i)
                  ro1b(i)=ro1(i)
                  ro2b(i)=ro2(i)
               end do
               
              if(it.gt.1 .and. new_evaluation.ge.3)then
                 call copia_hydro(ntau,hydro1,hydro1_conv,imodel2,hydro2,hydro2_conv)
              end if

c             Si amplio copio en la variable perfil el perfil sintetico.
              if(rnlte_th.lt.10.)then
                  do i=1,k40  !numero de longitudes de onda
                     perfil(i)=ymodobs_RH(i)
                  end do              
              else
                  do i=1,k40  !numero de longitudes de onda
                     perfil(i)=scal(i)
                  end do
              end if    

               fill2=atmosb(16*ntau+4)
               amac1=atmosb(8*ntau+1)
               amac2=atmosb(16*ntau+3)
               amic1=atmosb(3*ntau+1)*1.e-5
               amic2=atmosb(11*ntau+3)*1.e-5
               porcien=atmosb(16*ntau+5)

               contraste=contr
               if(nlin(ntlblends).eq.0)contraste=scal(nliobs)
               contraste=(1+contraste)/(1.-contraste)

               call comprime2(ntau,m,atmos,atmosr)
               do j=1,mfit
                  atmosrlin(j)=atmosr(j)
               end do 
               
               if(ngu.lt.25)then
                  cotavar=exp((-42+float(ngu))/3.)-1.e-4
                  if(cotavar.gt.0.0033)cotavar=0.0033
               else
                  if(ngu.lt.50)cotavar=(ngu-24)*.003
                  if(ngu.ge.50)cotavar=(ngu-49)*.08
                  if(ngu.eq.50)then
                     call error(KWARN,'desire','More iterations will not be'
     &               //         ' permitted unless significant variations\n'
     &               //         ' of chi**2 occur. If you still want to use'
     &               //         ' the same nodes, restart\n'
     &               //         ' the inversion using the output models')
                  end if
               end if

               if(abs(varchiLTE).le.cotavar*300.)then
                  isigo=0
                  write(msg,'(a,e10.3,a,e10.3)')
     &                  'Stop iterating because chi variation = ',varchiLTE,
     &                  ' <= ',  cotavar*300.
                  call error(KLINE,'',msg)
               end if
               if(ngu.eq.100)then
                  isigo=0
                  call error(KWARN,'desire','100 iterations. If you want to'
     &            //         ' continue, restart the inversion\n and use the'
     &            //         ' output models of this cycle')
               end if

               chiw=sumsq/float(nfrecs-ielimino)
               snchi=1./sqrt(chiw)

               add=-20.
               if(diag0.gt.0.d0)add=alog10(diag0)
               chprint=chireal/float(nfree-ielimino)

               ll=ll+1
               nguvec(ll)=ngu
               addvec(ll)=add
               snchivec(ll)=snchi
               chprintvec(ll)=chprint
               amac1vec(ll)=amac1
               amac2vec(ll)=amac2
               amic1vec(ll)=amic1
               amic2vec(ll)=amic2
               fill2vec(ll)=fill2
               contrastevec(ll)=contraste
               porcienvector(ll)=porcien
               posicionciclo(ici)=ll

               kv=0
               do i=1,18
                  mvec(ll,i)=m(i)
               end do

               if(iauto.eq.0) then
                  write(msg,1002)ngu,add,snchi,chprint,amac1,amac2,
     &                           amic1,amic2,fill2,contraste,porcien
                  call error(KLITE,'',msg)
               end if
               if(iauto.eq.1) then
                  write(msg,2002)ngu,add,snchi,chprint,(m(jj),jj=1,18)
                  call error(KLITE,'',msg)
               end if

            end if  !end if(iamplio.eq.1)

         end do  !end do while(isigo.eq.1)

c        ===============================================================
c                        Acaba el ciclo iterativo
c        ===============================================================

         chireal=chireal_store
         sumsq=sumsq_store

         do i=1,16*ntau+5
            atmos(i)=atmosb(i)
            atmosout(i)=atmosb(i)
         end do
         do j=1,mfit
            atmosr(j)=atmosrlin(j)
         end do 

         if(ngu.eq.0)then
            if(chprint .gt. umbralchi)then
               call error(KWARN,'desire','Finishing without changing'
     &         //         ' initial model')
            end if
         end if  !03/03/21 brc: termino el if(ngu.eq.0)

         if(rnlte_th.ge.10. .and. new_evaluation.ne.0)then
            call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,perfil)
            call leemodi222(1,modelout1,modelout2,atmosb,ntau,ierror,
     &                      z1,pg1,ro1,z2,pg2,ro2)
            if(new_evaluation .ge. 2)then
               if(ierror .eq. 2)then !solo atmosfera 1
                  do i=1,8*ntau+2
                     errores(i)=abs(atmosb(i)) !initializing with the atmosphere value
                  end do
                  do j=1,8
                     ini1=j*ntau
                     mj=m(j)
                     if(j .eq. 2)mj=m(1)
                     call evaluate_errors(ntau,mj,ini1,ini1,atmos1err,errores)
                  end do               
               else
                  do i=1,16*ntau+5  !2 componentes
                     errores(i)=abs(atmosb(i)) !initializing with the atmosphere value
                  end do
                  do j=1,8
                     ini1=j*ntau
                     ini2=(8+j)*ntau+2
                     mj=m(j)
                     mj8=m(j+8)
                     if(j .eq. 2)mj=m(1)
                     if(j .eq. 2)mj8=m(9)
                     call evaluate_errors(ntau,mj,ini1,ini1,atmos1err,errores)
                     call evaluate_errors(ntau,mj8,ini1,ini2,atmos2err,errores)
                  end do 
               end if            
               call leemodi222(1,modelerr1,modelerr2,errores,
     &                         ntau,ierror,z1*0,z1*0,z1*0,z2*0,z2*0,z2*0)
            end if
         end if

c        ---------------------------------------------------------
c        Hacemos una sintesis NLTE
         iwriteRHformat=0
         if(rnlte_th.lt.10. .and. new_evaluation.ne.0)then
               if(new_evaluation .eq. 4)iwriteRHformat=1
               ievaluo=1          
               IRH1=1
               if(imodel2.eq.1)IRH2=1
               
               iamplio=1
               nciclosor=nciclos
c              Desactivamos las FR numericas en la sintesis final pero
c              guardamos la variable para reactivarla mas adelante.               
               do inum=1,18
                  numericalor(inum)=numerical(inum)
                  numerical(inum)=0
               end do
               diagor=diag
  
               startj=NEW_J
               call sirkeywords(nrays,startj)             
               
               call comprime2(ntau,m,atmos,atmosr)
               call marquardt2(stok,sig,nfrecs,atmosr,m,mfit,
     &                      covar,alpha,chisq,12345.5,iamplio,beta)


               call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,ymodobs_RH)
               call leemodi222(1,modelout1,modelout2,atmos,ntau,ierror,
     &                      z1,pg1,ro1,z2,pg2,ro2)
     
               if(new_evaluation .ge. 2)then
                  if(new_evaluation .ge. 3)then
                    if(imassortau .eq. 0)then
                       call write_hydro(mod_hyd_out1,ntau,tau,hydro1)                  !writing H populations model1 if exist (ihydro1=1)
                       if(imodel2 .eq. 1)call write_hydro(mod_hyd_out2,ntau,tau,hydro2)!writing H populations model2 if exist (ihydro2=1)
                    end if
                  end if
                  if(ierror .eq. 2)then !solo atmosfera 1
                     do i=1,8*ntau+2
                        errores(i)=abs(atmos(i)) !initializing with the atmosphere value, used as a threshold (uncertainties < 2*abs(atmos(i)))
                     end do
                     j=1 !Temperature
                     ini1=j*ntau
                     call evaluate_errors(ntau,m(1),ini1,ini1,atmos1err,errores)
                     j=2 !Pe
                     do i=1,ntau
                        tup(i)=atmos(ntau+i)+errores(ntau+i)
                        if(tup(i).gt. 2*atmos(ntau+i))tup(i)=2*atmos(ntau+i)
                        tdown(i)=atmos(ntau+i)-errores(ntau+i)
                        if(tdown(i) .lt. 0.5*atmos(ntau+i))tdown(i)=0.5*atmos(ntau+i)
                        if(tdown(i) .lt. 2000.)tdown(i)=2000.
                        peup(i)=1.05*atmos(2*ntau+i)    !initialization
                        pedown(i)=0.95*atmos(2*ntau+i)  !initialization
                     end do   
                     call equisubmu(ntau,tau,tup,peup,pgup,zup,roup)
                     call equisubmu(ntau,tau,tdown,pedown,pgdown,zdown,rodown)
                     do i=1,ntau
                        errorpe1=abs(peup(i)-atmos(2*ntau+i))
                        errorpe2=abs(pedown(i)-atmos(2*ntau+i))
                        if(errorpe2 .gt. errorpe1)errorpe1=errorpe2
                        errores(2*ntau+i)=errorpe1  !uncertainties of Pe
                        errorz1=abs(zup(i)-z1(i))
                        errorz2=abs(zdown(i)-z1(i))
                        if(errorz2 .gt. errorz1)errorz1=errorz2
                        errz1(i)=errorz1  !uncertainties of Pg
                        errorpg1=abs(pgup(i)-pg1(i))
                        errorpg2=abs(pgdown(i)-pg1(i))
                        if(errorpg2 .gt. errorpg1)errorpg1=errorpg2
                        errpg1(i)=errorpg1  !uncertainties of Pg
                        errorro1=abs(roup(i)-ro1(i))
                        errorro2=abs(rodown(i)-ro1(i))
                        if(errorro2 .gt. errorro1)errorro1=errorro2
                        errro1(i)=errorro1  !uncertainties of density                                                            
                     end do
                   
                     do j=3,8
                        ini1=j*ntau
c                       Errors as input is the value of the parameter.
c                       As output contains the uncertainties.
                        call evaluate_errors(ntau,m(j),ini1,ini1,atmos1err,errores)
                     end do               
                  else
                     do i=1,16*ntau+5  !2 componentes
c                       Initializing with the atmosphere value, used as
c                       threshols (uncertainties < 2*abs(atmos(i))).
                        errores(i)=abs(atmos(i))
                     end do
                     do j=1,8
                        ini1=j*ntau
                        ini2=(8+j)*ntau+2
c                       Errors as input is the value of the parameter.
c                       As output contains the uncertainties.
                        mj=m(j)
                        mj8=m(j+8)
                        if(j .eq. 2)mj=m(1)
                        if(j .eq. 2)mj8=m(9)                       
                        call evaluate_errors(ntau,mj,ini1,ini1,atmos1err,errores)
                        call evaluate_errors(ntau,mj8,ini1,ini2,atmos2err,errores)
                     end do 
                  end if            
                  call leemodi222(1,modelerr1,modelerr2,errores,ntau,
     &                            ierror,errz1,errpg1,errro1,z2*0,z2*0,z2*0)
               end if

               diag=diagor
               nciclos=nciclosor
               do inum=1,18
                  numerical(inum)=numericalor(inum)
               end do
 
               if(iauto.eq.0) then
                  write(msg,1002)ngu,add,snchi,chprint,amac1,amac2,
     &                           amic1,amic2,fill2,contraste,porcien
                  call error(KLITE,'',msg)
               end if
               if(iauto.eq.1) then
                  write(msg,2002)ngu,add,snchi,chprint,(m(jj),jj=1,18)
                  call error(KLITE,'',msg)
               end if

         end if  !end if (rnlte_th.lt.10. .and. new_evaluation.ne.0)

         close(ican)

         if(ilast.lt.1)then
            snpr(ici)=snn         !viene del common de la iterarion anterior
            chisnpr(ici)=chisn    !idem
         else
            snpr(ici)=snchi       !ultima iteracion
            chisnpr(ici)=chprint  !idem
         end if

         call error(KLINE,'','')
         write(msg,'(a,i2,a)') '============================ End of cycle ',
     &         ici,' ================================='
         call error(KPARA,'',msg)

      end do   !end main loop...........................................

c     Write the ".chi" file.
c     05/05/20 epm: To run in parallel, we need this information.
      open(78,file=snychi)
      write(78,*) snchigood,escalachigood/chirealgood,chinewgood
      close(78)

c     Write the ".log" file.
      if(flagnolog.eq.0)then
         open(icanal,file=control)
         if(iauto.eq.0)write(icanal,1000) "ite","DE","s/neqv","chi**2","Mac1",
     &                       "Mac2","mic1","mic2","fill2","Ic1/Ic2","stray"
         if(iauto.eq.1)write(icanal,2000) "ite","DE","s/neqv","chi**2",
     &                       (ca(kk),kk=1,18)

         write(icanal,786)0,snpr(1),chisnpr(1)

         do i=1,ll
            if(iauto.eq.0)write(icanal,1002)nguvec(i),addvec(i),snchivec(i),
     &                          chprintvec(i),amac1vec(i),amac2vec(i),
     &                          amic1vec(i),amic2vec(i),fill2vec(i),
     &                          contrastevec(i),porcienvector(i)
            if(iauto.eq.1)write(icanal,2002)nguvec(i),addvec(i),snchivec(i),
     &                          chprintvec(i),(mvec(i,jj),jj=1,18)

            do j=1,nciclos-1
               if(i.eq.posicionciclo(j))then
                  write(icanal,*)''
                  if(iauto.eq.0)write(icanal,1000)"ite","DE","s/neqv","chi**2",
     &                                "Mac1","Mac2","mic1","mic2","fill2",
     &                                "Ic1/Ic2","stray"
                  if(iauto.eq.1)write(icanal,2000)"ite","DE","s/neqv","chi**2",
     &                                (ca(kk),kk=1,18)

                  write(icanal,786)0,snpr(j+1),chisnpr(j+1)
               end if
            end do
         end do
         close(icanal)
      end if

      call error(KLINE,'',' __________________________________________________________ ')
      call error(KLINE,'','|                                                          |')
      call error(KLINE,'',CARTEL)
      call error(KLINE,'','|__________________________________________________________|')
      call error(KLINE,'','')

c     10/05/19 epm: Measure time before existing the program.
      call chrono(cpu2, wall2)
      write(msg,'(a,i8,a,i8)') 'CPU  time (ms) = ', cpu2 - cpu1,
     &                      '\n Wall time (ms) = ', wall2 - wall1
      call error(KPARA,'',msg)
      if (flagquiet.eq.1 .and. flagtime.eq.1) then
         write(*,2003)'s/n=',snchivec(ll),'chi**2=',chprintvec(ll),
     &                'CPU=',cpu2 - cpu1,'Wall=',wall2 - wall1
      end if

c     Formats.
786   format(1x,i3,6x,1pe9.2,1x,e10.3)
1000  format(2x,a2,2x, a2,5x,a3,6x,a6,2x,  a6,1x,a6,3x,  a4,3x,a4,3x,
     &       a5,1x,a7,1x,a5)
1002  format(1x,i3,1x,f4.0,1x,1pe9.2,1x, e10.3,1x, 0pf6.3,5(1x, f6.3),1x,f6.3)
2000  format(2x,a2,2x, a2,5x,a3,6x,a6,3x,18(a3))
2002  format(1x,i3,1x,f4.0,1x,1pe9.2,1x, e10.3, 18(i3))
2003  format(1x,a,1pe9.2,3x,a,e10.3,3x,a,i8,3x,a,i8)

c     Happy ending.
      stop

c     07/07/20 epm: Error reading the command line.
665   call error(KSTOP,'desire','Invalid argument: '//arg)

c     21/06/19 epm: Error opening the control file.
666   call error(KSTOP,'desire','File not found: '//fcontrol)

      end

c ______________________________________________________________________

      subroutine diagonal(diag0,diag,isigo,iamplio,varchi,it)

      implicit real*4(a-h,o-z)
      include 'PARAMETER'
      real*4 diagon(100)
      character diver*11,espacio*35
      character*100 msg
      common/iteradiagonal/itt
      common/repeticion/factorrep !factor de diag para evitar repeticiones
      data ir/1/

      diver=' increases '
      espacio='---------------------------------'

      if(it.eq.1)then
         ir=1
         itt=0
      endif

      if(ir.eq.1)diagon(ir)=diag0
      ir=ir+1
      if(ir.eq.100)ir=6
      diagon(ir)=diag

      if(abs(diag).lt.1.e-15)then
         isigo=0
         iamplio=0
         return
      else if(diag.gt.diag0)then
         isigo=1
         iamplio=0
         itt=itt+1
      else if(diag.le.diag0)then
         isigo=1
         iamplio=1
         itt=0
      end if

      if(ir.gt.5)then
         dir02=abs(diagon(ir)-diagon(ir-2))
         dir13=abs(diagon(ir-1)-diagon(ir-3))
         if(dir02 .lt. 1.e-8 .or. dir13 .lt. 1.e-8)then
            factorrep=0.
         else
            factorrep=1.
         end if
      end if

      if(ir.gt.7)then
         dir01=abs(diagon(ir)-diagon(ir-1))
         dir03=abs(diagon(ir)-diagon(ir-3))
         dir04=abs(diagon(ir)-diagon(ir-4))
         dir05=abs(diagon(ir)-diagon(ir-5))

         if(dir01 .lt. 1.e-12 .and.dir02 .lt. 1.e-12 .and.
     &      dir03 .lt. 1.e-12 .and.dir04 .lt. 1.e-12 .and.
     &      dir05 .lt. 1.e-12 )factorrep=1.
      end if

      if(diag0.gt.0)then
         if(diag.gt.diag0)then
            varchi=1.0
            add=-20.
            if(diag0.gt.0.0)add=alog10(diag0)
            write(msg,1100)add,diver,espacio
            call error(KLITE,'',msg)
            if(itt.ge.4)then
               isigo=0
               iamplio=0
            end if
         else
           if(varchi .lt. 1.e-3)then
               isigo=0
               iamplio=0
           end if   
         end if
      end if

1100  format(5x,f4.0,1x,a11,a)
      return
      end

c_____________________________________________________________________________
c
c rutina convierte
c Reescribe el perfil de luz difusa en las longitudes de onda de la 
c malla
c entradas: s     - perfil de luz difusa en las ldo observadas
c           nposi - array de enteros que contiene las posiciones
c                  de la malla correspondientes a las de las ldo
c                   observadas
c           nliobs- numero total de longitudes de onda observadas
c           nli   - numero total de longitudes de onda de la malla
c salidas:  stray - (via common) perfil de luz difusa en las ldo de la malla
c Basilio 24 Junio 1994
c_____________________________________________________________________________

       subroutine convierte(s,nposi,nliobs,nli)

       include 'PARAMETER'
       parameter (kld4=4*kld)
       real*4 stray(kld4),s(*)
       integer nposi(*)
       common/difusa/stray

       do i=1,nliobs
          stray(nposi(i))=s(i)
       end do

       return
       end

c_____________________________________________________________________________

       subroutine coor(coorfile,deglat,deglon,Vx,Vy)

       implicit real*4 (a-h,o-z)
       include 'PARAMETER'
       character*(*) coorfile
       real*4 Xc,Yc,Vx,Vy,deglat,deglon
       integer icoor,ican

       ican=1
       icoor=meves(coorfile,100)
       if(icoor.eq.1)then
          open(ican,file=coorfile,status='old',err=760)
c         Xc, Yc=position/Rsun, Vx, VY horizontal vel cm/s(1.95e5)
          read(ican,*,err=770,END=770)Xc,Yc,Vx,Vy
          close(ican)
c         03/03/21 epm: De momento solo funciona con valores cero.
          if ((abs(Xc)+abs(Yc)+abs(Vx)+abs(Vy)) .gt. 0.001) then
             call error(KWARN,'coor','The synthesis/inversion is done in line'
     &       //               ' of sight reference system (i.e. mu=1)')
          end if
          Vx = 0
          Vy = 0
          deglat = 0
          deglon = 0
c         if(abs(Xc) .lt. 0.0003)then
c            if(Xc .lt. 0.)Xc=-0.0003
c            if(Xc .ge. 0.)Xc=0.0003
c         endif
c         if(abs(Yc) .lt. 0.0001)then
c            if(Yc .lt. 0.)Yc=-0.0001
c            if(Yc .ge. 0.)Yc=0.0001
c         endif
c         call LatLonfromXcYc(Xc,Yc,deglat,deglon)
          return
       end if

760    call error(KSTOP,'coor','The coordinates file does not exist\n'
     & //               ' File: '//coorfile)

770    call error(KSTOP,'coor','Incorrect format in the coordinates file\n'
     & //               ' File: '//coorfile)

       return
       end

c_____________________________________________________________________________

        subroutine coorP(coorfile,deglat,deglon,Vx,Vy)

c Reading a file with extension .coorP that contains Xc,Yc, P & B0
c evaluates the longitud (deglon), latitud (deglat) and the components
c Vx & Vy of the rotational velocity in a Local reference frame which X axis is
c parallel to the horizontal.
c INPUT:
c     Xc & Yc coordinates of the observation point in units of solar radius
c     P & B0 angles defining the solar rotation axis in degrees
c     You will need to look for the angles P,B0 of the
c     Sun for the observation date
c     you can use the IDL procedure /home/brc/idlsoft.dir/PandB0.pro
c     or http://bass2000.obspm.fr/ephem.php
c OUTPUT:
c     deglat & deglon: latitude and longitude of the observation point (in radians)
c     Vx, Vy velocity in cm/s

       implicit real*4 (a-h,o-z)
       include 'PARAMETER'
       character*(*) coorfile
       real*4 Xc,Yc,Vx,Vy,deglat,deglon,P,B0
       integer icoor,ican

       ican=1
       icoor=meves(coorfile,100)
       if(icoor.eq.1)then
          open(ican,file=coorfile,status='old',err=860)
c         Xc, Yc=position/Rsun, P, B0 solar parameters in degrees
          read(ican,*,err=870,END=870)Xc,Yc,P,B0
          close(ican)
c         03/03/21 epm: De momento solo funciona con valores cero.
          if ((abs(Xc)+abs(Yc)+abs(Vx)+abs(Vy)) .gt. 0.001) then
             call error(KWARN,'coor','The synthesis/inversion is done in line'
     &       //               ' of sight reference system (i.e. mu=1)')
          end if
          Vx = 0
          Vy = 0
          deglat = 0
          deglon = 0
c         call LatLonfromXcYcPB(Xc,Yc,P,B0,deglat,deglon,Vx,Vy)
          return
       end if

860    call error(KSTOP,'coorP','The coordinates file does not exist\n'
     & //               ' File: '//coorfile)

870    call error(KSTOP,'coorP','Incorrect format in the coordinates file\n'
     & //               ' File: '//coorfile)

       return
       end

c_____________________________________________________________________________

       subroutine LatLonfromXcYc(Xc,Yc,deglat,deglon)

       implicit real*4 (a-h,o-z)
       real*4 Xc,Yc,deglat,deglon,xxc,epsi,pi180,rc 

       epsi=1.e-5
       pi180=57.2957795131

       rc=sqrt(Xc*Xc+Yc*Yc)
       if(1.-rc .lt. epsi)then
         Xc=Xc*(1.-epsi)/rc
         Yc=Yc*(1.-epsi)/rc
       endif
       deglat=asin(Yc)
       xxc=Xc/cos(deglat)
       if(xxc .lt. (-1.))xxc=-1.+epsi
       if(xxc .gt. (1.))xxc=1.-epsi
       deglon=asin(xxc)*pi180
       deglat=deglat*pi180

       return
       end

c_____________________________________________________________________________

c Evaluates the horizontal component (Vx,Vy) of the
c tangential differential rotation of the SUN for
c a point (Xp,Yp) in units of the solar radius  over the projected disk
c (Xp>0 to the rightside= West, Yp>0 to the upside= North)

       subroutine LatLonfromXcYcPB(Xc,Yc,P,B0,deglat,deglon,Vx,Vy)

       implicit real*4 (a-h,o-z)
       real*4 Xc,Yc,P, B0,deglat,deglon,Vx,Vy
       real*4 Rc,Zc,sinlat,epsi,pi180 
       real*4 cos2theta,cos4theta,Omega,Velocity

       epsi=1.e-5
       pi180=57.2957795131

       Rc=sqrt(Xc*Xc+Yc*Yc)
       if(Rc .gt. 1.0)then    !the point should be inside the disk
          Xc=Xc/Rc
          Yc=Yc/Rc
       endif

       Zc=sqrt(1.-Xc*Xc-Yc*Yc)  

       call LatLonfromXcYc(Xc,Yc,deglat,deglon)

c We need to rotate a angle P counter clock side (to the east) around Z axis
c and later rotate around X axis to the South
       sinlat=-cos(B0/pi180)*sin(P/pi180)*Xc+ 
     & cos(B0/pi180)*cos(P/pi180)*Yc +sin(B0/pi180)*Zc

c Following Kuker & Rudiger (2008)
c http://iopscience.iop.org/article/10.1088/1742-6596/118/1/012029/meta
c the rotation velocity is
c Omega=(14.05 -1.492*cos^2(theta) -2.606*cos^4(theta)) degrees/day with theta=colatitude
c so cos(theta)=sinlat

       cos2theta=sinlat*sinlat
       cos4theta=cos2theta*cos2theta

       Omega=14.05 -1.492*cos2theta -2.606*cos4theta   ! degrees/day
       Omega=Omega/pi180/86400.               ! radians/second
       Velocity=695700.*Omega                          ! Km/s

c this velocity is parallel to the equator so, we need to rotate the reference
c System, an angle P clock side (to the west) around Z axis a angle P
       Vx=Velocity*cos(P/pi180)*1.e5
       Vy=Velocity*sin(P/pi180)*1.e5

       return
       end

c_____________________________________________________________________________

c renormaliza divide los perfiles observados por el continuo RH del HSRA

        subroutine renormaliza(ist,ntl,npas,stok)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        real*4 stok(*),si(kld),sq(kld),su(kld),sv(kld)
        real*4 continuohNLTEarr(kld)
        integer ntl,k,n,npas(*),ist(*)
        common/continuosNLTEarr/continuohNLTEarr 

        n=0
        do i=1,ntl
           do j=1,npas(i)
              n=n+1
           end do
        end do      
        call iquvfroms(ist,n,si,sq,su,sv,stok)
        k=0
        do i=1,ntl
           do j=1,npas(i)
              k=k+1
              con=continuohNLTEarr(k)
              si(k)=si(k)/con
              sq(k)=sq(k)/con
              su(k)=su(k)/con
              sv(k)=sv(k)/con                    
           end do  
        end do           
        call sfromiquv(ist,n,si,sq,su,sv,stok)   

        return
        end

c_____________________________________________________________________________   

        subroutine evaluate_errors(ntau,mm,ini1,ini2,atmoserr,errores)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        integer ntau,mm
        real*4 atmoserr(*),errores(*)
        common/signaltonoise/sn !signal to noise input in main program
        
        if(mm .ge. 0)then
            do i=1,ntau
               ileo=ini1+i
               iescribo=ini2+i
               errfac=atmoserr(ileo)
               if(errfac .gt. 2*errores(iescribo))errfac=2*errores(iescribo)
               errores(iescribo)=errfac
            end do
        end if         

        return
        end

c_____________________________________________________________________________   
