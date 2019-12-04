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
c Versions:
c
c 1.00 (20/05/19): Version base con llamadas al sistema (2019 = año 1).
c 1.09 (09/09/19): Primera version con ejecutables convertidos en funciones.
c 1.10 (10/10/19): Coeficientes de Barklem por argumento desde SIR a RH.
c 1.11 (11/11/19): Codigo RH nuevo y codigo SIR depurado.
c
c_____________________________________________________________________________
c
c     DIMENSIONES
      implicit real*4 (a-h,o-z)

      include 'PARAMETER'

      parameter (kt16=16*kt+5,na=92)         
      parameter (kl4=4*kl,kld4=4*kld)
      parameter (kldt=kld*kn,kldt4=4*kldt)
      parameter (mmax=16*kn+2,kldn=mmax*kld4)
c     parameter (mfitmax=200)  !incluido en PARAMETER
      parameter (kn16=16*kn+4,kt11=11*kt+2,kt12=11*kt+3)

c     PARA LA MALLA
      integer nlin(kl),nlins(kl4),npas(kl),npass(kl4),npos(kld4),nposi(kld),indice(kl)
      integer nlincheck(kl),imaya,ici,iauto
      integer ist(4),nble(kl),nlinsn(kl),il
      real*4 dlamda0(kl),dlamda(kld),pist(4)
      real*4 dlamdas(kld4),cth,vx,vy
      integer m(18), HEon, numerical(18),numericalor(18)
      integer icanal,ican,ncifno,ncifnoOUT,ixt1,ixt2
      real abu(na)

c     PARA LA ATMOSFERA
      integer icalerr,ncontpg,ipgmag,multiformat1,multiformat2
      integer ntau,ierror,nrays,ndepths
      real*4 atmos(kt16),atmosr(kt16),atmosrlin(kt16),atmosr_old(kt16)
      real*4 atmosout(kt16),atmoslin(kt16),tau(kt)
c     real*4 atmosoutold(kt16)
      real*4 pg1(kt),z1(kt),ro1(kt),T1(kt),Pe1(kt)
      real*4 pg2(kt),z2(kt),ro2(kt),T2(kt),Pe2(kt)
      real*4 B1(kt),gamma1(kt),phi1(kt)
      real*4 B2(kt),gamma2(kt),phi2(kt)
      real*4 B1_RH(kt),gamma1_RH(kt),phi1_RH(kt)
      real*4 B2_RH(kt),gamma2_RH(kt),phi2_RH(kt)
      real*4 pg1b(kt),z1b(kt),ro1b(kt),deglat,deglon
      real*4 pg2b(kt),z2b(kt),ro2b(kt)
      real*4 pg01,pg02,ro01,ro02,alog10mu
      real*4 Vmac1,fill1,strayfac,Vmac2,fill2    
      real*4 depth(kt),temp(kt),Nelec(kt),Vz(kt),Vtur(kt)
      real*4 tauoriginal(kt)

c     PARA EL PERFIL
      real*4 stok(kld4),sig(kld4),stray(kld4),sigd(kld4) !sigd es copia de sig
      real*4 pesos(kld4),ymodobs_RH(kld4)
      real*4 scal(kld4),smax(kl,4),vic(kl),perfil(kld4)
      integer npasobs(kl),numberfrec(4),ntl
      real*4 dlamdaobs(kld),dlamdacheck(kld)

c     PARA LA INVERSION
      integer mreadi2,mreadi3,meves,iratio(4)
      real*4 covar(mfitmax,mfitmax),alpha(mfitmax,mfitmax),beta(mfitmax)
      real*4 mreadr2,mreadr3                 !,ileoNLTE

c     PARA LOS COEFICIENTES DE ALEJAMIENTO
      integer iRH1,iRH2
      integer iRH1_store(101),iRH2_store(101)
      real*4 rnlte_th

c     PARA LOS ERRORES
c     real*4 dvcal(kldn),x(mfitmax)   !,taur(kt)
c     real*4 errores(mfitmax)
c     integer ifiable(kn16),indi(kn16)
      real*4 maximo

c     PARA EL FICHERO LOG
      integer nguvec(500),posicionciclo(20),mvec(500,18),mvecmax(18)
      real*4  addvec(500),snchivec(500),chprintvec(500),amac1vec(500),amac2vec(500)
      real*4  amic1vec(500),amic2vec(500),fill2vec(500),contrastevec(500)
      real*4  porcienvector(500)
      real*4  chisnpr(20),snpr(20)

c     CADENAS
      character*100 uveobs,uveout,modelout1,modelout2,modelin1,modelin2,mod1inic,mod2inic
      character*100 modelin2ini,coorfile,extensioncoor,extensionmodel
      character*100 model_multi1,model_multi2
      character*101 Bmodel_multi1,Bmodel_multi2
      character*100 malla,control,controlb,fcontrol,modelerr1,modelerr2,difusa
      character*100 pesos_file
      character*100 snychi
      character*100 filtro,psfper,extensionfiltro
      character cha1*2,cha2*2,chag*1,chamod*4,chaper*4,chaerr*4,chalog*4,ca(18)*3
      character chaweight*4
      character chasnychi*4
      character*100 mreadc2
      character*100 nomlineas,fichabun
      character*80 men1,men2,men3
      character*120 linlarga
      character*70 cartel
      character*100 vprint,ruta,ruta2,ruta3,filewavename
      character*4 cth_char

c     RH
      character*100 RH_model,RH_abundance,RH_magneticfield,label_ID_model !BRC-RH Jun 20 2017
      character*100 starting_j      !,RH_rays,RH_starting_j
c     character*100 key1,key1back,key2back
c     real*4 pe1_change(kt),pe2_change(kt),pg1_change(kt),pg2_change(kt)
c     real*4 atmosnew(kt16)

c     10/05/19 epm: To measure time.
      integer*8 cpu1, cpu2, wall1, wall2

c     21/06/19 epm: Command line.
      character*100 arg
      integer*4 flagv

      data iprimeravez_ruta/0/
      data nraysini/3/
c     data (pe1_change(i), i=1,kt)/kt*1.0/
c     data (pe2_change(i), i=1,kt)/kt*1.0/
c     data (pg1_change(i), i=1,kt)/kt*1.0/
c     data (pg2_change(i), i=1,kt)/kt*1.0/

c     COMUNES
      common/responde2/ist,ntau,ntl,nlin,npas,nble
      common/ldeo/dlamda,dlamda0
      common/smalla/ntls,nlins,npass
      common/smalla1/dlamdas
      common/atmosfera/atmos
      common/atmosferaout/atmosout 
c     common/atmosferaoutold/atmosoutold 
      common/iamplioold/iamplioold
      common/tol/tol
      common/filfac/fill
      common/alamda0/alamda0
      common/uvesalida/scal,ymodobs_RH
      common/observaciones/stok,sigd
      common/mu/cth                   !para equisub en amp22
      common/sigrealchi/sig,chireal,sumsq
      common/canal/icanal
      common/nombrecontrol/control
      common/cambio/ncambiono,ncambioprec          !para splines
      common/ifiltro/ifiltro
      common/filtro/filtro
c     common/nlte/nlte
      common/contraste/contr
      common/repeticion/factorrep !factor de diag para evitar repeticiones
      common/nciclos/nciclos   !para marquardt2, leeuve3, leemodi22
c     common/derivchi/dvcal  ! deriv. chi^2 para los errores
      common/posiciones/npos 
c     common/errores/errores
      common/calerr/icalerr     !si calculo errores=1 else =0 (para amp2)
      common/difusa/stray
c     common/fiable/ifiable,indi
c     El common siguiente es para pasarle a marquardt2 los indices iniciales y 
c     finales de las perturbaciones a gamma y fi aditivas
c     common/ifies/iga1,ifi11,iga2,ifi22 
      common/ifies/ipa1,ipa11,iga1,ifi11,ipa2,ipa22,iga2,ifi22
      common/ieliminofrec/ielimino !para el print de la S/N en marquardt2
      common/ivez/ivez
      common/ficlineas/nomlineas   !para leelineasii
      common/fichabun/fichabun     !para atmdatb
      common/numberfrec/numberfrec,iratio     !para marqcoef2
      common/nommod/modelin1,modelin2       !para escribeFR
      common/atmosr/atmosr                  !para escribeFR
      common/tau/tau
      common/iautomatico/iauto ! se pasa a fperfil2fperfil2
      common/iprimeravez/iprimeravez !se cambia en fperfil2
      common/nspl_lin/nspl_lin !para seleccionar interp.lineal=1 o splines=0
                               !o penalty (sol. regularizada =3 o 2, 5 o 4)
      common/contornopg/ncontpg,pg01,pg02
      common/contornoro/ro01,ro02
      common/primerchi/chisn,snn
c     common/departcoef/beta1,beta2
      common/zetas/pg1,z1,ro1,pg2,z2,ro2    !para el calculo de z,pg y ro en amp22 y amp22err
      common/pgmag/ipgmag
c     17/06/19 brc: common repetido con valores absolutos.
c     common/mu/cthabs
      common/anguloheliocent/xmu    !for equisubmu and related routines
                                    !to evaluate equ. hidr. in the z direction,
                                    !not along the LOS
      common/alog10mu/alog10mu             !alog10mu=alog10(xmu) from desire
      common/thresholdNLTE/rnlte_th
      common/RHnames/label_ID_model,RH_model,RH_magneticfield
      common/vprint/vprint
      common/iRH/iRH1,iRH2  ! eq 1 we call RH, 0 don't call RH
c     common/ruta/ruta,ruta2,filewavename
      common/ruta/ruta,ruta2,ruta3,filewavename
      common/HEonoff/HEon  !HEon=1 HE, HEon=0 Non HE
      common/istatus12/istatusdep1,istatusdep2
      common/variablePSF/psfper
      common/numerical/numerical
      common/latlong/deglat,deglon
c     common/atmosSIRfromRH/atmosnew,pe1_change,pe2_change,pg1_change,pg2_change 
c     common/tauoriginal/tauoriginal
c     common/ileoNLTE/ileoNLTE

c     22/08/19 epm: Common to save command line flags.
      common/commandline/flagv

c     EXTERNAS
      external mreadi2,mreadr2,mreadc2,meves,mreadi3,mreadr3

c     10/05/19 epm: Initial time.
      call chrono(cpu1, wall1)

c ntau : numero de puntos en tau
c tau  : log10(tau)
c ntl  : numero total de lineas
c nlin : indice de cada linea
c npas : numero de puntos en cada linea
c dlamda:cada delta de l.d.o. en ma

c_______________________________________________________________________

c     10/10/19 epm: Print removing any leading and trailing blanks.

c     21/06/19 epm: Read the command line.
      fcontrol = ''
      flagv = 0
      i = 1
      do
         call get_command_argument(i, arg) !0 = executable, 1 = first argument
         if (len_trim(arg) .eq. 0) exit    ! no more arguments, leave the loop
         if (trim(arg) .eq. '-h') then
            print*,''
            print*,'desire [flags] <control-file>'
            print*,''
            print*,'-h : print this help message'
            print*,'-v : show all the messages on console'
            print*,''
            stop
         else if (trim(arg) .eq. '-v') then
            flagv = 1
         else
            fcontrol = arg
         end if
         i = i + 1
      end do

      print*, ''
      print*, ' __________________________________________________________ '
      print*, '|                                                          |'
      CARTEL=' |            DeSIRe version 1.11  (11/Nov/2019)            |'
      print'(a)',CARTEL
      print*, '|__________________________________________________________|'
      print*, ''
      print*, ''

      if (len_trim(fcontrol) .eq. 0) then
         write(*,'(a,$)') ' Control file: '
         read(*,'(a)') fcontrol
      end if

c     Inicializa algunas variables

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

c	open(icanal,file=control,fileopt='eof')    
      linlarga='------------------------------------------------------------------------------'

c 	Se leen las condiciones de ejecucion

	do i=1,18
           mvecmax(i)=0  !para los errores
        end do
c        iprimeravez=0
        nciclos=1
	ici=0
	ll=0
	do while (ici.lt.nciclos)
        ivez=0
        ici=ici+1 
        iprimeravez=0

c        iprimeravez=ici-1
	ncambiono=0	!cambio el numero de nodos
	ncambioprec=0	!cambio la precision
      open(ican,file=fcontrol,status='old',err=666)
	
	nciclos =mreadi2(ican,ici)        !numero de ciclos (repeticiones)
	uveobs  =mreadc2(ican,ici)        !fichero de entrada perfiles
	difusa  =mreadc2(ican,ici)        !fichero entrada  luz difusa
	filtro  =mreadc2(ican,ici)        !fichero con transformada de PSF
        malla   =mreadc2(ican,ici)        !malla,(npas=n. puntos por linea)
	nomlineas=mreadc2(ican,ici)       !fichero con parametros atomicos
	fichabun =mreadc2(ican,ici)       !fichero con abundancias
        modelin1 =mreadc2(ican,ici)       !nombre del modelo 1 anterior 
        modelin2 =mreadc2(ican,ici)       !nombre del modelo 2 anterior
        coorfile =mreadc2(ican,ici)       !nombre del fichero con coordenadas X/Rsun, Y/Rsun, Vrot(cm/s)
        pist(1)  =mreadr3(ican,ici,1.)        !i (peso de i en ajuste)
        pist(2)  =mreadr3(ican,ici,1.)        !q
        pist(3)  =mreadr3(ican,ici,1.)        !u
        pist(4)  =mreadr3(ican,ici,1.)        !v
        iauto   =mreadi3(ican,ici,0)        !seleccion automatica de nodos
	m(1)    =mreadi3(ican,ici,0)        !nodos en t1
	m(2)    =mreadi3(ican,ici,0)        !nodos en pe1
	m(3)    =mreadi3(ican,ici,0)        !nodos en micro1
	m(4)    =mreadi3(ican,ici,0)        !   "   "    h1
	m(5)    =mreadi3(ican,ici,0)        !   "   "    vz1
	m(6)    =mreadi3(ican,ici,0)        !   "   "    gamma1
	m(7)    =mreadi3(ican,ici,0)        !   "   "    fi1
	m(8)    =mreadi3(ican,ici,0)        !   "   "    macro1
	m(9)    =mreadi3(ican,ici,0)        !   "   "    t2
	m(10)    =mreadi3(ican,ici,0)       !   "   "    pe2
	m(11)    =mreadi3(ican,ici,0)       !   "   "    micro2
	m(12)    =mreadi3(ican,ici,0)       !   "   "    h2
	m(13)    =mreadi3(ican,ici,0)       !   "   "    vz2
	m(14)    =mreadi3(ican,ici,0)       !   "   "    gamma2
	m(15)    =mreadi3(ican,ici,0)       !   "   "    fi2
	m(16)    =mreadi3(ican,ici,0)       !   "   "    macro2
	m(17)    =mreadi3(ican,ici,0)       !   "   "    f.f.2
        m(18)    =mreadi3(ican,ici,0)       !luz difusa
	sn      =mreadr3(ican,ici,1000.)    !segnal ruido estimada para i
	contr    =mreadr3(ican,ici,-1.)     !contraste Ic1/Ic2
	tol     =mreadr3(ican,ici,1.e-4)    !tolerancia inversion 
	alamda0 =mreadr3(ican,ici,1.e-3)    !factor diagonal inicial 
        nspl_lin =mreadi3(ican,ici,0)	    !0=splines, 1=lineal (regularizada splines(2 o 4) o lineal(3 o 5)
        pg01    =mreadr3(ican,ici,0.)       !condi. cont. pg atm. 1/negative value= density
        pg02    =mreadr3(ican,ici,0.)       !condi. cont. pg atm. 2/negative value= density
        rnlte_th =mreadr3(ican,ici,10.)     !NLTE thresh. (0=NLTE,>10 LTE)10 (LTE)= default
        nrays=  mreadi3(ican,ici,3)         !number of rays
        ruta    = mreadc2(ican,ici)
        if(nrays .eq. 0)nrays=nraysini
        if(ici.eq.1)nraysini=nrays
c        inew=1 !in RH: No using old values of background op & J. In DeSIRe evaluating .per after every cycle 
c        if(rnlte_th .lt. 0)then
c          inew=0 !in RH: Using old values of background op & J. In DeSIRe NO evaluating .per after every cycle 
c          rnlte_th=-1*rnlte_th
c        end if
        do inum=1,18
          numerical(inum)=0         !RF for parameter i numerically evaluated in NLTE
          if(m(inum) .ge. 1000)then
             m(inum)=m(inum)-1000
             numerical(inum)=1
          endif
          if(m(inum) .le. -1000)then
             m(inum)=m(inum)+1000
c             numerical(inum)=1 !in this case we are NOT evaluating RF because the perturbation i
c                               !taken equal to the other atmosphere
          endif
        end do  
        
        call extension3(coorfile,extensioncoor,n)
        print*,''
        if(extensioncoor(1:n) .eq. 'coor')then
c           print*,'DeSIRe: extension coorfile=',extensioncoor(1:n)
           call coor(coorfile,deglat,deglon,vx,vy)
           print*,'Vx=',Vx, 'cm/s Vy=',Vy,' cm/s'
        else if(extensioncoor(1:n) .eq. 'coorP')then
c           print*,'DeSIRe: extension coorfile=',extensioncoor(1:n)
           call coorP(coorfile,deglat,deglon,vx,vy)
           print*,'Vx=',Vx, 'cm/s Vy=',Vy,' cm/s'
        else
           print*,'DeSIRe: Coordinates file must have a .coor or a .coorP extension'
           print*,'.coor  contains Xc/Rsun Yc/Rsun Vx[cm/s] Vy[cm/s]'
           print*,'.coorP contains Xc/Rsun Yc/Rsun P B0'
           stop
        endif
        
        HEon=1          ! HEon=1 means HE, HEon=0 meaning OUT of HE
        if(abs(m(2))+abs(m(10)) .ne. 0 .and. nciclos .gt. 0)HEon=0
        
        multiformat1=0
        multiformat2=0
        call extension3(modelin1,extensionmodel,n)
        if(extensionmodel(1:n) .eq. 'multi' .or. extensionmodel(1:n) .eq. 'atmos')then
           multiformat1=1
           model_multi1=modelin1
           modelout1=modelin1  
	   call quitaex2(modelout1,ixt1)
           call concatena(modelout1,'.mod',modelin1) 
        end if  
       imodel2=0
        if(modelin2 .ne. ''.or.modelin2 .ne. ' ')then
           imodel2=1
           call extension3(modelin2,extensionmodel,n)
           if(extensionmodel(1:n) .eq. 'multi'.or. extensionmodel(1:n) .eq. 'atmos')then
              multiformat2=1
              model_multi2=modelin2
              modelout2=modelin2  
	      call quitaex2(modelout2,ixt1)
              call concatena(modelout2,'.mod',modelin2) 
           end if
        endif   
            
c defining the path from   solveray from the path for  rhf1d   
        iprimeravez_ruta=iprimeravez_ruta+1
        if(iprimeravez_ruta.eq.1)then
           islash=0
           isalgo=0
           call busco(' ',ruta,1,islashb)
           do while(islash .lt. islashb-1 .and. isalgo .ne.1)
              call busco('/',ruta(islash+1:100),1,islash1) 
              if(islash+islash1 .ge. islashb-1 )then
                 isalgo=1
              else
                 islash=islash+islash1
              end if
           end do
           ruta2=ruta(1:islash)//'solveray'
           ruta3=ruta(1:islash)//'conversion'           
        end if
        
        pi=3.14159265
	cth=cos(deglon*pi/180.)*cos(deglat*pi/180.)
        
        open(1,file='ray.input')
          xmu=abs(cth)
          if(xmu .gt.1.00)xmu=1.00
	  write(1,*)xmu
	  write(1,*)0
	close(1)
	
	write(cth_char, '(f4.2)' )xmu
	filewavename='spectrum_'//cth_char(1:4)//'.asc'
	
	alog10mu=alog10(xmu)
	cthabs=1.0     !to work in LRF
	xmu=1.0        !to work in LRF
        
        vprint='v'
        ipgmag=0  ! we do not consider magnetic pressure terms
c ____________________________________________________________________

        if(rnlte_th .lt. 0.)rnlte_th=0.
      	if(vprint .ne. 'q'.and.vprint.ne.'Q')print*,'  '
        if(vprint .ne. 'q'.and.vprint.ne.'Q')print*,trim(linlarga)
	
c	xmu=abs(cth)
c	if(cth .ne. 1.)then
c	   print*,'Resulting models are stratified versus line of sight'
c	   print*,'We will use heliocentric angle value only for hydrostatic equilibrium equation'
c           cth=1.  !por defecto, no pasamos de la linea de vision a z
c	endif
c        cthabs=abs(cth)

        ncontpg=0
        ro01=0.
        ro02=0.
        if(pg01.gt.0. .or. pg02.gt.0. .and. abs(m(2))+abs(m(10)).eq.0)then
           print*,'Boundary condition in gas pressure'
	   if(pg02.lt.0.)pg02=pg01
	   if(pg01.lt.0.)pg01=pg02
           ncontpg=1
        end if
        if(pg01.lt.0. .or. pg02.lt.0. .and. abs(m(2))+abs(m(10)).eq.0)then
           print*,'Boundary condition in density'
	   if(pg02.lt.0.)pg02=pg01
	   if(pg01.lt.0.)pg01=pg02
           pg01=abs(pg01)
           pg02=abs(pg02)
           ro01=pg01
           ro02=pg02
           ncontpg=-1
        end if
	if(ncontpg.eq.0 .and. abs(m(2))+abs(m(10)) .eq.0)then
           print*,'Boundary condition in electronic pressure'
	   print*,'Better results are expected setting boundary condition in gas pressure/density'
	end if

c cadenas para el output
        ca(1)=' t1'      
        ca(2)=' p1'      
        ca(3)=' m1'      
        ca(4)=' h1'      
        ca(5)=' v1'      
        ca(6)=' g1'      
        ca(7)=' f1'      
        ca(8)=' M1'      
        ca(9)=' t2'      
        ca(10)=' p2'      
        ca(11)=' m2'      
        ca(12)=' h2'      
        ca(13)=' v2'      
        ca(14)=' g2'      
        ca(15)=' f2'      
        ca(16)=' M2'      
        ca(17)=' ff'      
        ca(18)=' st'     


	if(nciclos.eq.0)then   !para sintesis uso los 4 parametros (LRB)
          do i=1,4
            iratio(i)=0     !iratio(Q)=1 -> se invierte Q/I, etc
            ist(i)=1
            if(pist(i).lt.0)iratio(i)=1
            pist(i)=1
          enddo
          iratio(1)=0  !I no se puede sintetizar como I/I
	else
          do i=1,4
             ist(i)=0
             iratio(i)=0     !iratio(Q)=1 -> se invierte Q/I, etc
             if(pist(i).lt.0)iratio(i)=1
             pist(i)=abs(pist(i))
             if(pist(i).gt.0)ist(i)=1
          end do
	  if(ist(1).eq.0)print*,'Stokes I is not considered in this inversion'
	  if(ist(2).eq.0)print*,'Stokes Q is not considered in this inversion'
	  if(ist(3).eq.0)print*,'Stokes U is not considered in this inversion'
	  if(ist(4).eq.0)print*,'Stokes V is not considered in this inversion'
	endif
        if(iratio(2).eq.1.or.iratio(3).eq.1.or.iratio(4).eq.1)iratio(1)=0

c definimos el fichero con los parametros atomicos
        ilineas=0
	ilineas=meves(nomlineas,100)
        if(ilineas.eq.0)then 
            print*,''
            print*,'STOP: Atomic parameters file is not defined in the control file'
            stop ' '
c		men1='The default LINES file is being used for supplying atomic parameters.'
c		print*,trim(men1)
c                nomlineas='~/sir/default/LINES'
        endif

c definimos el fichero con las abundancias
        iabun=0
	iabun=meves(fichabun,100)
        if(iabun.eq.0)then 
            print*,''
            print*,'STOP: Abundance file is not defined in the control file'
            stop ' '
c		men1='The default ABUNDANCES file is being used for supplying abundances.'
c		print*,trim(men1)
c                fichabun='~/sir/default/ABUNDANCES'
        endif

c definimos los nombres de los ficheros de salida
	chamod='.mod'
	chaerr='.err'
	chaper='.per'
	chaweight='.wht'

        ncifno=1
c        if(ici.ge.10)ncifno=2 
        ncifnoOUT=1
        if(ici.gt.9)ncifnoOUT=2
        if(ici.gt.10)ncifno=2 
        chag='_'
	cha1(1:2)='  '
	cha2(1:2)='  '
	modelout1=modelin1
	modelout2=modelin2
	modelin2ini=modelin2
	pesos_file=uveobs
  
	call quitaex2(modelout1,ixt1)
        call quitaex2(modelout2,ixt2)
        call quitaex3(pesos_file,ixt3)
        call concatena(pesos_file,chaweight,pesos_file)
     
        if(ici.eq.1.or.ici.eq.0)then
             mod1inic=modelout1
             mod2inic=modelout2
	     call chachi(1+ixt1,cha1,ncifno)
	     call chachi(1+ixt2,cha2,ncifno)	     
        else
             if(modelout1.eq.mod1inic)then 
	        call chachi(ici-1+ixt1,cha1,ncifno)
                call concatena(modelout1,chag,modelin1)
                call concatena(modelin1,cha1,modelin1)
                call concatena(modelin1,chamod,modelin1)
c	        call chachi(ici+ixt1,cha1,ncifno)
 	        call chachi(ici+ixt1,cha1,ncifnoOUT)
             else
	        call chachi(1+ixt1,cha1,ncifno)
             endif
             if(imodel2.eq.1)then
	        if(modelout2.eq.mod2inic )then
		  call chachi(ici-1+ixt2,cha2,ncifno)
		  call concatena(modelout2,chag,modelin2)
		  call concatena(modelin2,cha2,modelin2)
		  call concatena(modelin2,chamod,modelin2)
		  call chachi(ici+ixt2,cha2,ncifnoOUT)        	        
		else
		  call chachi(1+ixt2,cha2,ncifno)
		endif
	     endif	
        endif

        call concatena(modelout1,chag,modelout1)
        call concatena(modelout1,cha1,modelout1)
        call concatena(modelout1,chaper,uveout)
        call concatena(modelout1,chaerr,modelerr1)
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
        if(ifiltro.eq.0)then 
            men1='PSF NOT taken into account'
            print*,trim(men1)
        endif
	
	if(ifiltro.eq.1)then
	     men1='The PSF will be read from: '//filtro(1:40)
	     print*,trim(men1)
             call extension3(filtro,extensionfiltro,n)
             if(extensionfiltro(1:n) .eq. 'per')then
		psfper=filtro
		ifiltro=2
		print*,'Reading PSF variable with wave from ',trim(psfper)
	     endif	
	endif		

	istray=meves(difusa,100)
        if(istray.eq.0)then
            men1='NO stray light file has been specified.'
            print*,trim(men1)
        endif
 
        if(nspl_lin.eq. 0)print*,'Splines interpolation between nodes '
        if(nspl_lin.eq. 1)print*,'Linear interpolation between nodes '
        if(nspl_lin.eq. 2)print*,'Light regularized solution (Splines interpolation between nodes)'
	if(nspl_lin.eq. 3)print*,'Light regularized solution (Linear interpolation between nodes)'
	if(nspl_lin.eq. 4)print*,'Regularized solution (Splines interpolation between nodes)'
	if(nspl_lin.eq. 5)print*,'Regularized solution (Linear interpolation between nodes)'

c       ------------------------------------------------------
c       leemos la malla, el perfil ,el modelo y la luz difusa:
c       ------------------------------------------------------
        imodel1=meves(modelin1,100)
  	imodel2=meves(modelin2,100)
        if(imodel1.eq.1.and.imodel2.eq.1)ierror=0 !existen ambos
	if(imodel1.eq.0.and.imodel2.eq.1)ierror=1 !no existe mod1 pero si mod2 
	if(imodel1.eq.1.and.imodel2.eq.0)ierror=2 !si existe mod1 pero no mod2
  	if(imodel1.eq.0.and.imodel2.eq.0)ierror=3 !no existe ninguno de los 2 modelos
  	
        if(ierror.eq.2)then
          men1='     WARNING: Model 2 does NOT exist. Only model 1 is considered'
          print*,trim(men1)
	endif
	if(ierror.eq.1 )then
          print*,''
          print*,'STOP: Model 1 does NOT exist'
          stop ' '
	endif
	
  	label_ID_model=modelin1 
  	call read_keyword_input_RH('ATMOS_FILE',RH_model)           !BRC-RH Jun 19 2017
c 	call read_keyword_input_RH('ABUND_FILE',RH_abundance)       !BRC-RH Jun 19 2017
c 	call read_keyword_input_RH('NRAYS',RH_rays)
 	if(iprimeravez_ruta .eq. 1 .or. nrays.ne.nraysini)call write_keyword_input_RH('NRAYS',nrays)                  !BRC-RH Jul 3 2018	 
c        call read_keyword_input_RH('STARTING_J',RH_starting_j)
c        STARTING_J      ='NEW_J'  
        if(iprimeravez_ruta .eq. 1)starting_j='NEW_J'
        if(iprimeravez_ruta .gt. 1)starting_j='OLD_J'
        if(iprimeravez_ruta .le. 2)call write_keyword_input_RHstr('STARTING_J',starting_j)                  !BRC-RH May 2 2019	 
c 	call write_keyword_input_RH('NRAYS',nrays)                  !BRC-RH Jul 3 2018	 
  	call read_keyword_input_RH('STOKES_INPUT',RH_magneticfield) !BRC-RH Jun 20 2017  

c       if(nciclos.eq. -10)then  
            call read_keyword_input_RH('ABUND_FILE',RH_abundance)
            call leeabun(0,abu) 
            call write_abun_RH(abu,RH_abundance)
            call leemallab(malla,ntl,nli,nlin,npas,dlamda,nble)
c           call lines_kurucz(nomlineas)

           if(multiformat1 .ne. 0 .and. ici.lt.2)then          !multiformat1=1 --> input model en multi format 
           !if the input is in multi format  we tranform a mass column
           !we have to write the input file name in a file with the name in keyword.input:  RH_model
           !for instance we are copying FALCtau.atmos to Temporary.atmos
           call read_model_atmos(model_multi1,ndepths,depth,temp,Nelec,Vz,Vtur)
           ntau=ndepths          
           print*,'copying ',trim(model_multi1),' in multi format in file: ',trim(RH_model)
           call copying(model_multi1,RH_model)  
           
           Bmodel_multi1='B'//model_multi1
           print*,'Reading magnetic field from file: ',trim(Bmodel_multi1),' from surface to bottom'
           do i=1,ndepths
              B1_RH(i)=0.
              gamma1_RH(i)=0.
              phi1_RH(i)=0.
              b1(i)=0.
              gamma1(i)=0.
              phi1(i)=0.
           end do
           Vmac1=0.
           fill1=1.
           strayfac=0.
           open(50,file=Bmodel_multi1,status='old',err=991)
           read(50,*,err=771,end=991)Vmac1,fill1,strayfac
           read(50,*,err=771)B1i,gamma1i,phi1i !si hay una sola linea valores constantes
           do i=1,ndepths
              B1_RH(i)=B1i*1.e-4
              gamma1_RH(i)=gamma1i*3.14159265/180.
              phi1_RH(i)=phi1i*3.14159265/180.
              b1(i)=B1i
              gamma1(i)=gamma1i
              phi1(i)=phi1i
           end do

	   do i=2,ndepths
	      read(50,*,err=991,end=991)B1i,gamma1i,phi1i
	      B1_RH(i)=B1i*1.e-4
              gamma1_RH(i)=gamma1i*3.14159265/180.
              phi1_RH(i)=phi1i*3.14159265/180.	 
              b1(i)=B1i
              gamma1(i)=gamma1i
              phi1(i)=phi1i
           end do 
           goto 991
771           print*,'The model ',trim(Bmodel_multi1),' does not exist or has a bad format'
              print*,'It should contain two or more rows with 3 values each:'
              print*,'Vmac[km/s], filling factor[0-1], stray [%] in the first row' 
              print*,'B1i[G], gamma1i[degrees], phi1i[degrees] in the second [and later rows]'
              stop
991       close(50)  
             open(51,file=RH_magneticfield)
             do i=1,ntau
                write(51,*)B1_RH(i),gamma1_RH(i),phi1_RH(i)
             end do
             close(51)
c working in Local reference frame (no LoS)          
             call read_atmos_RH(RH_model,1,atmos,b1,gamma1,phi1,z1,pg1,ro1,ntau,Vmac1,fill1,strayfac)
             tau1=atmos(1)
             tau1=nint(tau1*100.)/100.
             taun=atmos(ntau)
             step=(tau1-taun)/float(ntau-1)
             step=nint(step*1000.)/1000.
	     taun=tau1-float(ntau-1)*step 
	  
             do i=1,ntau
                tauoriginal(i)=tau1-(i-1)*step
             end do 
             call interpolate1(tauoriginal,atmos,ntau,1,z1,pg1,ro1)
             call leemodi222(1,modelin1,modelin2,atmos,ntau,2,
     &                       z1,pg1,ro1,z1,pg1,ro1)

c          call read_atmos_RH(RH_model,1,atmos,b1,gamma1,phi1,z1,pg1,ro1,ntau,Vmac1,fill1,strayfac)
c          tau1=atmos(1)
c          tau1=nint(tau1*100.)/100.
c          taun=atmos(ntau)
c          step=(tau1-taun)/float(ntau-1)
c          step=nint(step*1000.)/1000.
c	  taun=tau1-float(ntau-1)*step 
c          do i=1,ntau
c             tauoriginal(i)=tau1-(i-1)*step
c          end do 
c          call interpolate1(tauoriginal,atmos,ntau,1,z1,pg1,ro1)
c          call leemodi222(1,modelin1,modelin2,atmos,ntau,2,
c     &                       z1,pg1,ro1,z1,pg1,ro1)

c And now with atmosphere 2:       
          if(multiformat2 .eq. 1)then
             call read_model_atmos(model_multi2,ndepths,depth,temp,Nelec,Vz,Vtur)
             ntau=ndepths          
             print*,'copying ',trim(model_multi2),' in multi format in file: ',trim(RH_model)
             call copying(model_multi2,RH_model)
             
             Bmodel_multi2='B'//model_multi2
             print*,'Reading magnetic field (second component) from file: ',trim(Bmodel_multi2)
             do i=1,ndepths
                B2_RH(i)=0.
                gamma2_RH(i)=0.
                phi2_RH(i)=0.
                b2(i)=0.
                gamma2(i)=0.
                phi2(i)=0.
             end do 
             Vmac2=0.
             fill2=0.
             open(50,file=Bmodel_multi2,status='old',err=992)
             read(50,*,err=772,end=992)Vmac2,fill2,strayfac2
             read(50,*,err=772)B2i,gamma2i,phi2i !si ahy una sola linea valores constantes
             do i=1,ndepths
                B2_RH(i)=B2i*1.e-4
                gamma2_RH(i)=gamma2i*3.14159265/180.
                phi2_RH(i)=phi2i*3.14159265/180.
                b2(i)=B2i
                gamma2(i)=gamma2i
                phi2(i)=phi12
            end do

	    do i=2,ndepths
	       read(50,*,err=992,end=992)B2i,gamma2i,phi2i
	       B2_RH(i)=B2i*1.e-4
               gamma2_RH(i)=gamma2i*3.14159265/180.
               phi2_RH(i)=phi2i*3.14159265/180.	 
               b2(i)=B2i
               gamma2(i)=gamma2i
               phi2(i)=phi2i
           end do 
772           print*,'The model ',trim(Bmodel_multi2),' does not exist or has a bad format'
              print*,'It should contain two or more rows with 3 values each:'
              print*,'Vmac[km/s], filling factor[0-1], stray [%] in the first row' 
              print*,'B1i[G], gamma1i[degrees], phi1i[degrees] in the second [and later rows]'
              stop                      
992       close(50)  
             open(51,file=RH_magneticfield)
             do i=1,ntau
                write(51,*)B2_RH(i),gamma2_RH(i),phi2_RH(i)
             end do
             close(51)
          
             call read_atmos_RH(RH_model,1,atmos,b2,gamma2,phi2,z2,pg2,ro2,ntau,Vmac2,fill2,strayfac)
             tau1=atmos(1)
             tau1=nint(tau1*100.)/100.
             taun=atmos(ntau)
             step=(tau1-taun)/float(ntau-1)
             step=nint(step*1000.)/1000.
	     taun=tau1-float(ntau-1)*step 
             do i=1,ntau
                tauoriginal(i)=tau1-(i-1)*step
             end do 
             call interpolate1(tauoriginal,atmos,ntau,2,z2,pg2,ro2)
             call leemodi222(1,modelin1,modelin2,atmos,ntau,0,
     &                       z1,pg1,ro1,z2,pg2,ro2)  
     
             end if
           end if 
c           stop
c        endif   
c          call read_atmos_RH(RH_model,1,atmos,b2,gamma2,phi2,z2,pg2,ro2,ntau,Vmac2,fill2,strayfac)
c          tau1=atmos(1)
c          tau1=nint(tau1*100.)/100.
c          taun=atmos(ntau)
c          step=(tau1-taun)/float(ntau-1)
c          step=nint(step*1000.)/1000.
c	  taun=tau1-float(ntau-1)*step 
c          do i=1,ntau
c             tauoriginal(i)=tau1-(i-1)*step
c          end do 
c          call interpolate1(tauoriginal,atmos,ntau,2,z2,pg2,ro2)
c          call leemodi222(1,modelin1,modelin2,atmos,ntau,0,
c     &                       z1,pg1,ro1,z2,pg2,ro2)          
c          end if
c        end if

        if(imodel2 .eq. 0)then ierror=1
        call leemodi222(0,modelin1,modelin2,atmos,ntau,ierror,z1,pg1,ro1,z2,pg2,ro2)

  	if(m(2) .eq. 0 )then
           do i=1,ntau
	       tau(i)=atmos(i)
	       T1(i)=atmos(ntau+i)
	       Pe1(i)=atmos(2*ntau+i)
	   end do  
	   print*,'Running HE in DeSIRe (model 1)!'
           if(ncontpg.eq.0) call equisubmu(ntau,tau,T1,Pe1,pg1,z1,ro1)
	   if(ncontpg.ne.0) call equisubmu_cont(ntau,tau,T1,Pe1,pg01,pg1,z1,ro1)
	   do i=1,ntau
	      atmos(2*ntau+i)=Pe1(i)
	   end do 
  	end if
        if(imodel2.eq.1 .and. m(10) .eq. 0 .and. nciclos .ge. 1)then
            print*,'Running HE in DeSIRe (model 2)!'
            do i=1,ntau
               tau(i)=atmos(8*ntau+2+i)
	       T2(i)=atmos(9*ntau+2+i)
	       Pe2(i)=atmos(10*ntau+2+i)
	    end do     
            if(ncontpg.eq.0)call equisubmu(ntau,tau,T2,Pe2,pg2,z2,ro2)
            if(ncontpg.ne.0)call equisubmu_cont(ntau,tau,T2,Pe2,pg02,pg2,z2,ro2)
            do i=1,ntau
	       atmos(10*ntau+2+i)=Pe2(i)	          
	    end do   		
  	end if

c	----------------------------------------------
c       Calculo de la atmosfera en la linea de vision:
c	----------------------------------------------
        il=0
        call taulinea4(il,imodel2,deglon,deglat,vx,vy,atmos,z1,pg1,ro1,z2,pg2,ro2,ntau)  

c        do i=1,16*ntau+5
c	   atmosoutold(i)=atmos(i)
c	end do        
        
c  	call read_keyword_input_RH('ATMOS_FILE',RH_model)           !BRC-RH Jun 19 2017
c 	call read_keyword_input_RH('ABUND_FILE',RH_abundance)       !BRC-RH Jun 19 2017
c 	call read_keyword_input_RH('NRAYS',RH_rays)
c	call write_keyword_input_RH('NRAYS',nrays)                  !BRC-RH Jul 3 2018	 	
c  	call read_keyword_input_RH('STOKES_INPUT',RH_magneticfield) !BRC-RH Jun 20 2017  

c        call leeabun(0,abu) 
c        call write_abun_RH(abu,RH_abundance)
c  	label_ID_model=modelin1 

	if(ierror.eq.2)then    !no podemos invertir el modelo 2!!
            do jj=9,17
               m(jj)=0
            end do
        endif

        
        if(atmos(16*ntau+5).gt.1e-5.or.m(18).ne.0)then
          if(istray.ne.0)then
            call leeuveobsindic(difusa,ist,ntl,nliobs,nlin,npasobs,dlamdaobs,nble,stray)  !(LRB)
	      
	    ntlcheck=ntl  !para comprobar que difusa tiene las mismas lambdas que perf. obs.
	    do l=1,nliobs
	      dlamdacheck(l)=dlamdaobs(l)
	    enddo
	    do l=1,ntl
	      nlincheck(l)=nlin(l)
	    enddo

c	    -----------------------------------------------------------------
c           Por si no tenemos el fichero con la malla en sintesis:
c	    (si lo tenemos, estas variables se machacan luego y no pasa nada)
            nli=nliobs
            do i=1,ntl
               npas(i)=npasobs(i)
            end do
            do i=1,nli
               nposi(i)=i
               dlamda(i)=dlamdaobs(i)
            end do
c	    -----------------------------------------------------------------

           else
	     men1='   WARNING: The STRAY LIGHT FACTOR is being changed to ZERO'
             if(atmos(16*ntau+5).gt.0.)print*,trim(men1)
             atmos(16*ntau+5)=0.
             m(18)=0  !LB (si no casca, porque trato de invertir una variable = cero)
           end if
        end if

c -----------------------------------------------------------------
c Para el calculo de las FR
c -----------------------------------------------------------------
          if(nciclos.eq.0)then
            m(2)=1
          end if  
	  if(nciclos.lt.0)then !para calculo de las FR mnodos=ntau
            men1=' STOP: The number of depth points in the model' 
            men2=' is larger than the current kn value' 
            men3=' Decrease the number of depth points or change the PARAMETER file.' 

c            if(kn.lt.ntau)call mensaje(4,men1,men2,men3)
            nRF=ntau
            if(kn.lt.ntau)nRF=kn
 
            if(nciclos.eq.-1)then
              do i=1,15
                if(m(i).ne.0)m(i)=nRF !; previously here said m(i)=ntau
              enddo
              if(m(1).ne.0 .and. m(2) .eq. 0)m(2)=nRF
              if(m(8).ne.0)m(8)=1 !la macro solo puede tener 1 nodo
            endif  
            iauto=0
          end if
c -----------------------------------------------------------------
	
        imaya=meves(malla,100)       !(LRB)

     	if(nciclos.lt.1)then
           uveout=uveobs
	   if(imaya.eq.1)then   !si tenemos malla
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
		 men1='The wavelength grid is being generated automatically from the'
		 men2='stray light profile. NO blends are considered.'
                 print*,trim(men1)
                 print*,trim(men2)
              else
	         men1='STOP: In synthesis mode, you MUST specify the wavelength grid'
		 call mensaje(1,men1,men2,men3)
	      endif
	   endif
	     
           print*,'Output profiles: ',trim(uveout)
                    
           if(ici.eq.1)then
c	     open(icanal,file=control,access='append')
	     open(icanal,file=control)
	        write(icanal,*)trim(linlarga )
	     close(icanal)
	   endif
       
	else 	
           if(imaya.eq.0)then
	     men1='The wavelength grid is being generated automatically from the'
	     men2='observed profiles. NO blends are considered.'
             print*,trim(men1)
             print*,trim(men2)
	           
             call leeuveobsindic(uveobs,ist,ntl,nliobs,nlin,npasobs,dlamdaobs,nble,stok)

             nli=nliobs
             do i=1,ntl
                npas(i)=npasobs(i)
             end do
             do i=1,nli
                nposi(i)=i
                dlamda(i)=dlamdaobs(i)
             end do
           else
             call leeuveobsindic(uveobs,ist,ntl,nliobs,indice,npasobs,dlamdaobs,nble,stok)
	     call leemalla2(malla,ntl,nli,nliobs,nlin,npasobs,npas,
     &                             dlamdaobs,dlamda,nble,nposi,indice)
           endif

	endif
	
	if(atmos(16*ntau+5).gt.0.)then   !comprobaciones varias !esto es si tenemos difusa

	  if(ntlcheck.ne.ntl)then
	    men1='STOP: The number of lines in the files containing the observed and'
	    men2='      stray light profiles are not equal.'
	    call mensaje(2,men1,men2,men3)
          endif

	  do l=1,ntl
	    if (l.eq.1)then
		numm=1
	    else
		numm=numm+nble(l-1)  
	    endif
	    if(nlin(numm).ne.nlincheck(l))then
	      men1='STOP: The order of the lines in the stray light file is not'
	      men2='      equal to that in the file containing the observed profiles.'
	      call mensaje(2,men1,men2,men3)
	    endif
	  enddo
	    
          do l=1,nli
	    if(abs(dlamda(l)-dlamdacheck(l)).gt.1.)then
	      men1='STOP: The wavelengths in the file containing the stray light profile differ by'
	      men2='      more than 1. mA from those generated by the wavelength grid:'
	      print*,'      Wavelength #: ',l,':',dlamda(l),dlamdacheck(l)
	      call mensaje(2,men1,men2,men3)
	    endif
	  enddo
	
	endif

	if(nciclos.ge.1)then
           if(ierror.eq.0)then
	      men1='Output models  : '//modelout1(1:20)//' and '//modelout2(1:20)
	      men2='Uncertainties  : '//modelerr1(1:20)//' and '//modelerr2(1:20)
	      men3='Output profiles: '//uveout(1:20)//' (2 components)'
	      print*,trim(men1)
	      print*,trim(men2)
	      print*,trim(men3)
		
           else if (ierror.eq.2)then
	      men1='Output model   : '//modelout1(1:50)
	      men2='Uncertainties  : '//modelerr1(1:50)
	      men3='Output profiles: '//uveout(1:43)//' (1 component)'

	      print*,trim(men1)
	      print*,trim(men2)
	      print*,trim(men3)

	   else 
              men1=' STOP: Models not specified or unexistent. '
	      call mensaje(1,men1,men2,men3)
	      stop
           end if
        end if

        print*,trim(linlarga)
	if(ici.eq.1)then
c	   open(icanal,file=control,access='append')
	   open(icanal,file=control)	   
	   write(icanal,*) trim(linlarga)
	endif 

c	---------------------
c     Sentencias de control
c	---------------------
 
        if(ierror.ne.0.and.contr.gt.-1.e-6)then
           contr=-1.e-6
	   men1=' No contrast calculation can be done without two models!'
           print*,trim(men1)
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
           print*,' '
           write(*,'(a34,i2,a38,i2)') 'STOP: The number of lines (',ntlblends,') is larger than the current limit kl=',kl
	     print*,'      Decrease this number or change the PARAMETER file.'
           print*,' '
           print*,'______________________________________________________________________________'
           stop 
        end if
	if (nli.gt.kld)then
           print*,' ' 
           write(*,'(a33,i4,a40,i4)') 'STOP: The number of wavelengths (',nli,') is larger than the current limit kld= ',kld
	     print*,'      Decrease the number of wavelengths or change the PARAMETER file.'
           print*,' '
           print*,'______________________________________________________________________________'
           stop 
        end if
	if (nfrecmax.gt.kld) then
           men1='STOP: There is a line with more wavelengths than the current limit kld'
	   men2='      Decrease the number of wavelengths or change the PARAMETER file.'
	   call mensaje(2,men1,men2,men3)
        end if
        if (ntau.gt.kt) then
	   men1='STOP: The model has more depth points than the current limit kt.'
	   men2='      Decrease the number of grid points or change the PARAMETER file.'
	   call mensaje(2,men1,men2,men3)
        end if

c       Definimos nlins,ntls,npass,dlamda0s,dlamdas:
	ntls=0
        nb=0
	k4=0
	k4c=0
	do i=1,ntl
	   vic(i)=1.     !inicializo i del continuo
           do ii=1,4
	     smax(i,ii)=1./sn
           end do
	end do

	ii=0
	vmax=-2.
	vmin=2.
        areaVazul=0.
	areaVroja=0.
	ipolaridad=1
	do i=1,4
	   do j=1,ist(i)
              iinli=ii*nli
              ii=ii+1
	      k3=0
	      k3c=0
	      suma=0
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
    	               if(i.eq.1.and.nciclos.gt.0)then
                           xmax=-1.
                           do lll=1,npasobs(k)
                             if(xmax.le.stok(k4-l+lll))xmax=stok(k4-l+lll)
                           end do
	                   vic(k)=xmax	!i continuo
                           sss=abs( vic(k)-sk4 )
	               else
	                   sss=abs(sk4)
	               end if
	               if(sss.gt.smax(k,i))smax(k,i)=sss
                     end if

	         end do !fin del do en frecuencias obs
                 do l=1,npas(k)
	            k3c=k3c+1
	            k4c=k4c+1
	            dlamdas(k4c)=dlamda(k3c)
	         end do !fin del do en frecuencias malla

	      end do 	!fin del do en lineas
	   end do       !fin del do en j
	end do          !fin del do en i
	
	if(nciclos.gt.0)call leepesos(pesos_file,ist,pesos)     
	
	maximo=smax(1,1)
	do k=1,ntl
            if(smax(k,1).gt.maximo)maximo=smax(k,1)
        enddo

	do k=1,ntl
           if(smax(k,1).lt..01.and.nciclos.gt.0)then
              print*,'Weight is being changed for the line',k
              smax(k,1)=maximo
           endif
        enddo     

	nfrecs=k4

c Para controlar el maximo numero de nodos en modo automatico ("*" en m(i))
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
           
           if(m(1).eq.1000)m(1)=2*nmax1 !temperatura (doble num de nodos)
           do i=2,7
              if(m(i).eq.1000)m(i)=nmax1
           end do
           if(m(9).eq.1000)m(9)=2*nmax1 !temperatura (doble num de nodos)
           do i=10,15
              if(m(i).eq.1000)m(i)=nmax1
           end do
           kv=0
           if(n_1000.gt.0.or.iauto.eq.1)then
             do i=1,18
                if(m(i).ne.0)then
                   print*,'maximun number of nodes for variable ',i,'= ',m(i)
                end if
             end do
           end if
	else
          do i=1,18
             if(m(i).eq.1000)then
               print*,'* is not allowed in non automatic mode'
               stop
             end if
          end do
   	end if

 	if(nciclos.gt.0)call nodos2aut(ntau,0,ierror,ici,m)      !en la version LBR estaba sin comentar
 
	mfit=0
	do i=1,18
           if(m(i).gt.0)mfit=mfit+m(i)
	end do

	if(mfit.gt.mfitmax .and. nciclos .gt. 0)then 
           print*,' ' 
           print*, 'STOP: The number of free parameters ',mfit
           print*, '      is larger than the current limit mfitmax',mfitmax
	   print*,'       Decrease the number of free parameters.'
           print*,' ' 
           stop
	endif
 
c       Escribe en atmosr los valores de atmos en los nodos:
	if(nciclos.ne.0)call comprime2(ntau,m,atmos,atmosr)

	nfree=nfrecs-mfit

      print*,' '
c     print*,'Number-of-observ-freqs/fitable-params/free-params: ',nfrecs,mfit,nfree	
      write(*,'(a,i8,1x,i6,1x,i8)') ' Number-of-observ-freqs/fitable-params/free-params: ',nfrecs,mfit,nfree
      print*,'______________________________________________________________________________'
      print*,' '
	
        if(mfit.gt.nfrecs.and.nciclos.gt.0)then
           print*,' ' 
           write(*,'(a37,i2,a27)') 'STOP: The number of free parameters (',mfit,') is larger than the number'
           write(*,'(a22,i3,a2)') '      of observables (',nfrecs,').'
           print*,' '
           print*,'______________________________________________________________________________'
           stop
        end if

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
                 numberfrec(i)=numberfrec(i)+npasobs(k) !numero de ldo para
                                             !cada stokes invertido 
                                             !es cero si no se invierte 

 	         sigc2=sigc*sqrt(smax(k,i)/vic(k))
	         do l=1,npasobs(k)
	            k40=k40+1
                    indicei=indicei+1

		    if(pesos(k40).gt.1.e-10)pesos_inv=1./pesos(k40)
                    if(pesos(k40).le.1.e-10)pesos_inv=1.e10
                    sig(k40)=sigc2    
                    
                    if(stok(k40).lt.cotaminima)then
                       sigx=sigmamax
                       sig(k40)=sigx
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
            write(*,'(i3,a32,f5.1,a34)') ielimino,' data points, being smaller than',cotaminima,', are NOT considered for inversion'
	endif

	if(nciclos.gt.0)then
	     if(iauto.eq.0)write( *,1000)"it","DE","s/n","chi**2","Mac1","Mac2","mic1","mic2"
     &                ,"fill2","Ic1/Ic2","stray"
             if(iauto.eq.1)write(*,2000)"ite","DE","s/neqv","chi**2",
     &                                   (ca(jj),jj=1,18)

        end if


c _____________________________________________________________________________
c                        empieza el ciclo iterativo
c _____________________________________________________________________________

        factorrep=1.
	diag=-1.0
	chi0=1.e20
	varchi=1.0
	isigo=1
	it=0
	iamplio=0
	ngu=0
        icalerr=0   !=0 pues de momento no vamos a calcular errores
        iRH1=0
        iRH2=0
c        if(rnlte_th .lt. 10. .and. ici.eq. 1)IRH1=1
c        if(rnlte_th .lt. 10. .and. imodel2.eq.1 .and. ici.eq. 1)IRH2=1
        if(rnlte_th .lt. 10. )IRH1=1
        if(rnlte_th .lt. 10. .and. imodel2.eq.1 )IRH2=1
        deltaT1_rel=0.
        deltaT2_rel=0.
        iamplioold=0
   	    
	do while(isigo.eq.1)
	    it=it+1
c            if(ici .gt. 1 )then
c              iRH1=0
c	      iRH2=0	
c	    endif  
	    if(it .gt. 1)then
	      iRH1=0
	      iRH2=0	
	    endif  
	    if(it .gt. 2)then
	      if(iamplio .eq. 0 .and. iRH1_store(it-1).eq.1 .and. iRH1_store(it-2).eq.0 .and. rnlte_th .lt. 10.)IRH1=1   
	      if(iamplio .eq. 0 .and. imodel2.eq.1 .and. iRH2_store(it-1).eq.1 .and. iRH2_store(it-2).eq.0 .and. rnlte_th .lt. 10.)IRH2=1
	    end if
	    if(it .gt. 1)then      
	      if(iamplio .eq. 1 .and. deltaT1_rel .gt. rnlte_th .and. rnlte_th .lt. 10.)IRH1=1
	      if(iamplio .eq. 1 .and. imodel2.eq.1 .and. deltaT2_rel .gt. rnlte_th .and. rnlte_th .lt. 10.)IRH2=1
	    end if
c	    print*,'DeSIRe 1276 it',it,iRH1

     	    iRH1_store(it)=IRH1
	    iRH2_store(it)=IRH2
	    
	    diag0=diag
	    
	    if(nciclos.gt.0)then
	       do inod=1,m(1)
	          atmosr_old(inod)=atmosr(inod)
	       enddo 
	       mm=m(1)+m(2)+m(3)+m(4)+m(5)+m(6)+m(7)+m(8)+m(9)+m(10)
	       do inod=mm+1,mm+m(11)
	          atmosr_old(inod)=atmosr(inod)
	       enddo 
	    endif   
	    	    
	    nn=8*ntau+2

	    call marquardt2(stok,sig,nfrecs,atmosr,m,mfit,
     &		         covar,alpha,chisq,diag,iamplio,beta)

c     	    iRH1_store(it)=IRH1
c	    iRH2_store(it)=IRH2
            if(nciclos.gt.0)then
               deltaT1_rel=0.
     	       do inod=1,m(1)
     	         deltaT1_reli=abs((atmosr(inod)-atmosr_old(inod))/atmosr_old(inod))
     	         if( deltaT1_reli .gt. deltaT1_rel)deltaT1_rel=deltaT1_reli
	       enddo 
	       if(imodel2.eq.1)then	    
	          deltaT2_reli=0.
	          do inod=mm+1,mm+m(11)
     	             deltaT21_reli=abs((atmosr(inod)-atmosr_old(inod))/atmosr_old(inod))
     	             if( deltaT2_reli .gt. deltaT2_rel)deltaT2_rel=deltaT2_reli
	          enddo 
	       endif 
	    endif   

	    if(nciclos.le.0)then
		if(nciclos.eq.0 .and. rnlte_th .ge.10)call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,scal)
		if(nciclos.eq.0 .and. rnlte_th .lt.10)call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,ymodobs_RH)
c           10/05/19 epm: Measure time before existing the program.
            call chrono(cpu2, wall2)
            print*,''
            print'(1x,a,i8)','CPU  time (ms) = ', cpu2 - cpu1
            print'(1x,a,i8)','Wall time (ms) = ', wall2 - wall1
            print*,''
            stop
	    end if

	    nfree=nfrecs-mfit
  
	    if(it.eq.1)diag0=alamda0
	
	    varchi=(chi0-chisq)/(chi0+chisq)

	    call diagonal(diag0,diag,isigo,iamplio,varchi,it)
            if(istatusdep1*istatusdep2 .ne. 0)iamplio=0
            if(istatusdep1*istatusdep2 .ne. 0)print*,'setting OLD atmosphere due to non RH convergence !!!!!!'
            istatusdep1=0
            istatusdep2=0
	    
	    if(iamplio.eq.1)then  !converge
		ngu=ngu+1
	        chi0=chisq

c copiamos la atmosfera de salida (via common) atmosout en atmos para no 
c machacar el common
c y la duplicamos en atmoslin para pasarla al centro del disco con taulinea4
c la inversion se realiza en la linea de vision pero los modelos de entrada 
c y salida estan en el centro del disco
c                 print*,'desire 1126 V gamma fi=', atmosout(5*ntau+10),atmosout(6*ntau+10),atmosout(7*ntau+10)
c	        do i=1,16*ntau+5
c	           atmosoutold(i)=atmos(i)
c	        end do
c	        iamplioold=1
	        do i=1,16*ntau+5
	           atmos(i)=atmosout(i)
	        end do


	        do i=1,ntau
	           pg1b(i)=pg1(i)
	           pg2b(i)=pg2(i)
	           z1b(i)=z1(i)
	           z2b(i)=z2(i)
	           ro1b(i)=ro1(i)
	           ro2b(i)=ro2(i)
	        end do

                do i=1,k40   !numero de longitudes de onda
                   perfil(i)=scal(i)
                end do


	        fill2=atmos(16*ntau+4)
	        amac1=atmos(8*ntau+1)
	        amac2=atmos(16*ntau+3)
	        amic1=atmos(3*ntau+1)*1.e-5
	        amic2=atmos(11*ntau+3)*1.e-5
                porcien=atmos(16*ntau+5)

                contraste=contr
	        if(nlin(ntlblends).eq.0)contraste=scal(nliobs)
                contraste=(1+contraste)/(1.-contraste)

               	call comprime2(ntau,m,atmos,atmosr)

	        if(ngu.lt.25)then
	           cotavar=exp((-42+float(ngu))/3.)-1.e-4
	           if(cotavar.gt.0.0033)cotavar=0.0033
                else 
                   if(ngu.lt.50)cotavar=(ngu-24)*.003
                   if(ngu.ge.50)cotavar=(ngu-49)*.08
                   if(ngu.eq.50)then
                      print*,'More iterations will not be permitted unless significant variations of chi**2 occur.'
                      print*,'If you still want to use the same nodes, restart the inversion using the output models.'
                   end if            
                end if
                

	        if(varchi.le.cotavar)isigo=0
                if(ngu.eq.100)then 
                   isigo=0
                   print*,'100 iterations. If you want to continue, restart the inversion and'
                   print*,'use the output models of this cycle'                 
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
				
c		if(isigo.eq.1)posicionciclo(ici)=ll
		posicionciclo(ici)=ll

				
                kv=0
                do i=1,18
                   mvec(ll,i)=m(i)
                   if(mvecmax(i).lt.m(i))mvecmax(i)=m(i)
                end do
                
	        if(iauto.eq.0 .and. vprint .ne. 'q'.and.vprint.ne.'Q')
     &	        write(*,1002)ngu,add,snchi,chprint,amac1,amac2,amic1,amic2,
     &                        fill2,contraste,porcien

	        if(iauto.eq.1.and. vprint .ne. 'q'.and.vprint.ne.'Q')
     &	        write(*,2002)ngu,add,snchi,chprint,(m(jj),jj=1,18)

	    end if

c            stop

	end do 
c _____________________________________________________________________________
c                        acaba el ciclo iterativo
c _____________________________________________________________________________



	if(vprint .ne. 'q'.and.vprint.ne.'Q')print*,' '
	if(vprint .ne. 'q'.and.vprint.ne.'Q')
     &              write(*,423)' ============================ End of cycle ',ici,' ================================='
	if(vprint .ne. 'q'.and.vprint.ne.'Q')print*,' '

	do i=1,16*ntau+5
	    atmoslin(i)=atmos(i)
	end do

	if(ngu.eq.0)then
             ierror=0 

	     call leemodi222(0,modelin1,modelin2,atmoslin,ntau,ierror,
     &                       z1,pg1,ro1,z2,pg2,ro2)
	     call leemodi222(1,modelout1,modelout2,atmoslin,ntau,ierror,
     &                       z1,pg1,ro1,z2,pg2,ro2)
     	    
     	    call comprime2(ntau,mvecmax,atmoslin,atmosr)
            do j=1,mfit
                atmosrlin(j)=atmosr(j)
            end do
    	    do i=1,16*ntau+5
	       atmos(i)=atmoslin(i)
	    end do
	    if(rnlte_th .lt. 10. )print*,'Finishing without changing initial model (NLTE) writing ',ymodobs_RH(1),ymodobs_RH(2)
	    if(rnlte_th .lt. 10. )call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,ymodobs_RH)
	    if(rnlte_th .ge. 10. )print*,'Finishing without changing initial model (LTE) writing ',scal(1),scal(2)
	    if(rnlte_th .ge. 10. )call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,scal)
	    call taulinea4(0,imodel2,deglon,deglat,vx,vy,atmos,z1,pg1,ro1,z2,pg2,ro2,ntau) 
	    do i=1,16*ntau+5
	       atmoslin(i)=atmos(i)
	    end do

 	 else
 	 
c	 call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,perfil)
         if(rnlte_th .ge. 10. )then
            call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,perfil)
         endif   
c         call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,ymodobs_RH)
      
         mfit=0

         if(mvecmax(1).gt.0)mfit=mfit+mvecmax(1)
         ipa1=mfit       !ipa1 es el indice anterior a la gamma comp. 1
         if(mvecmax(2).gt.0)mfit=mfit+mvecmax(2)
         ipa11=mfit
	 do i=3,5
           if(mvecmax(i).gt.0)mfit=mfit+mvecmax(i)
	 end do
         iga1=mfit    !iga1 es el indice anterior a la gamma comp. 1
         if(mvecmax(6).gt.0)mfit=mfit+mvecmax(6)
         if(mvecmax(7).gt.0)mfit=mfit+mvecmax(7)
         ifi11=mfit   !ifi11 es el indice ultimo de la fi comp. 1

	 do i=8,9
           if(mvecmax(i).gt.0)mfit=mfit+mvecmax(i)
	 end do
	 ipa2=mfit         !ipa2 es el indice anterior a la gamma comp. 2
	 if(mvecmax(10).gt.0)mfit=mfit+mvecmax(10)
	 ipa22=mfit
	 do i=11,13
           if(mvecmax(i).gt.0)mfit=mfit+mvecmax(i)
	 end do	 
         iga2=mfit    !iga2 es el indice anterior a la gamma comp. 2
         if(mvecmax(14).gt.0)mfit=mfit+mvecmax(14)
         if(mvecmax(15).gt.0)mfit=mfit+mvecmax(15)
         ifi22=mfit   !ifi22 es el indice ultimo de la fi comp. 1
	 do i=16,18
           if(mvecmax(i).gt.0)mfit=mfit+mvecmax(i)
	 end do

c         call comprime2(ntau,mvecmax,atmoslin,atmosr)
c         do j=1,mfit
c            atmosrlin(j)=atmosr(j)
c	 end do

c            IRH1=1
c	    if(imodel2.eq.1)IRH2=1
c	    iamplio=1
c	    call marquardt2(stok,sig,nfrecs,atmosrlin,m,mfit,
c     &		         covar,alpha,chisq,diag,iamplio,beta)
                 
	 do j=1,mfit
            atmosrlin(j)=atmosr(j)
	 end do

c	 if(iauto .eq. 1)then 
	    do i=1,ntau
	       tau(i)=atmos(i)
	       T1(i)=atmos(ntau+i)
	       Pe1(i)=atmos(2*ntau+i)
	    end do   

	    if(ncontpg.eq.0 .and. m(2) .eq. 0) call equisubmu(ntau,tau,T1,Pe1,pg1,z1,ro1)
	    if(ncontpg.ne.0 .and. m(2) .eq. 0) call equisubmu_cont(ntau,tau,T1,Pe1,pg01,pg1,z1,ro1)
	    do i=1,ntau
	       atmos(i)=tau(i)
	       atmos(ntau+i)=T1(i)
	       atmos(2*ntau+i)=Pe1(i)
	    end do 
	    
	    if(imodel2.eq.1)then
               do i=1,ntau
               	  tau(i)=atmos(8*ntau+2+i)
	          T2(i)=atmos(9*ntau+2+i)
	          Pe2(i)=atmos(10*ntau+2+i)
	       end do     
               if(ncontpg.eq.0 .and. m(10) .eq. 0)call equisubmu(ntau,tau,T2,Pe2,pg2,z2,ro2)
               if(ncontpg.ne.0 .and. m(10) .eq. 0)call equisubmu_cont(ntau,tau,T2,Pe2,pg02,pg2,z2,ro2)
               do i=1,ntau
               	  atmos(8*ntau+2+i)=tau(i)
	          atmos(9*ntau+2+i)=T2(i)
	          atmos(10*ntau+2+i)=Pe2(i)	          
	       end do   
            end if   
c         end if

c-----------------------------------------------------------------     
c hacemos una sintesis NLTE
c         if(ici .eq. nciclos .and. rnlte_th .lt. 10. )then
          ievaluo=0
c          if(rnlte_th .lt. 10.  .and. inew .eq. 1)ievaluo=1
c          if(rnlte_th .lt. 10.  .and. ici .eq. nciclos)ievaluo=1
c    	    do i=1,16*ntau+5
c	       atmos(i)=atmoslin(i)
c	    end do


	  call taulinea4(1,imodel2,deglon,deglat,vx,vy,atmos,z1,pg1,ro1,z2,pg2,ro2,ntau) 

          call leemodi222(1,modelout1,modelout2,atmos,ntau,ierror,
     &                       z1,pg1,ro1,z2,pg2,ro2)
          
          if(rnlte_th .lt. 10.)ievaluo=1

          if(ievaluo.eq. 1)then
            IRH1=1
	    if(imodel2.eq.1)IRH2=1
	    iamplio=1
	    nciclosor=nciclos
	    diagor=diag
c	    if(ici.eq.nciclosor)then
	       do inum=1,18
                  numericalor(inum)=numerical(inum)
                  numerical(inum)=0
               end do
c            endif  
c	    iprimeravezor=iprimeravez
c	    iprimeravez=0
	    diag=-1
	    nciclos=0
	    
            if(m(2) .eq. 0 )then
               do i=1,ntau
	          tau(i)=atmos(i)
	          T1(i)=atmos(ntau+i)
	          Pe1(i)=atmos(2*ntau+i)
	       end do  
               if(ncontpg.eq.0) call equisubmu(ntau,tau,T1,Pe1,pg1,z1,ro1)
	       if(ncontpg.ne.0) call equisubmu_cont(ntau,tau,T1,Pe1,pg01,pg1,z1,ro1)
	       do i=1,ntau
	          atmos(2*ntau+i)=Pe1(i)
	       end do 
c	       call taulinea4(0,imodel2,deglon,deglat,vx,vy,atmos,z1,pg1,ro1,z2,pg2,ro2,ntau) 
  	    end if
            if(imodel2.eq.1 .and. m(10) .eq. 0 .and. nciclos .ge. 1)then
               print*,'Running HE in DeSIRe (model 2)!'
               do i=1,ntau
                  tau(i)=atmos(8*ntau+2+i)
	          T2(i)=atmos(9*ntau+2+i)
	          Pe2(i)=atmos(10*ntau+2+i)
	       end do     
               if(ncontpg.eq.0)call equisubmu(ntau,tau,T2,Pe2,pg2,z2,ro2)
               if(ncontpg.ne.0)call equisubmu_cont(ntau,tau,T2,Pe2,pg02,pg2,z2,ro2)
               do i=1,ntau
	          atmos(10*ntau+2+i)=Pe2(i)	          
	       end do   		
  	    end if
  	    
            call taulinea4(0,imodel2,deglon,deglat,vx,vy,atmos,z1,pg1,ro1,z2,pg2,ro2,ntau) 
               
	    call marquardt2(stok,sig,nfrecs,atmosrlin,m,mfit,
     &		         covar,alpha,chisq,diag,iamplio,beta)
	    
            call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,ymodobs_RH)
c            iprimeravez=iprimeravezor
            diag=diagor
            nciclos=nciclosor
c            if(ici.eq.nciclosor)then
	       do inum=1,18
                  numerical(inum)=numericalor(inum)
               end do
c            endif  
          end if         
          
c          ievaluo=0
c fin de intrusion-----------------------------                   
c          if(ievaluo.eq. 1)then
c            print*,'Evaluating the Final NLTE profiles (after every cycle) using RH'
c	    do i=1,ntau
c	       pg1(i)=pg1b(i)
c	       pg2(i)=pg2b(i)
c	       z1(i)=z1b(i)
c	       z2(i)=z2b(i)
c	       ro1(i)=ro1b(i)
c	       ro2(i)=ro2b(i)
c	    end do
c              call comprime2(ntau,mvecmax,atmos,atmosrlin)
c              IRH1=1
c	      if(imodel2.eq.1)IRH2=1
c	      iamplio=1
c	      nciclosor=nciclos
c	      diagor=diag
c	      iprimeravezor=iprimeravez
c	      iprimeravez=0
c	      diag=-1
c	      nciclos=0
c	      
cc    	      do i=1,16*ntau+5
cc	         atmos(i)=atmoslin(i)
cc	      end do
c	      
c	      call marquardt2(stok,sig,nfrecs,atmosrlin,m,mfit,
c     &		         covar,alpha,chisq,diag,iamplio,beta)
c     
cc              do i=1,16*ntau+5
cc	         atmoslin(i)=atmos(i)
cc	      end do
c     
c              nciclos=nciclosor
c              diag=diagor
c              iprimeravez=iprimeravezor
c              
c              print*,'writing output profiles at file:',   uveout
c              call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,ymodobs_RH)
cc             call taulinea2(1,cth,ihemi,vx,atmos,ntau)
cc             call taulinea3(1,dglon,deglat,vx,atmos,ntau)
c
cc              call taulinea4(1,imodel2,deglon,deglat,vx,vy,atmos,z1,pg1,ro1,z2,pg2,ro2,ntau) 
cc	      call leemodi222(1,modelout1,modelout2,atmos,ntau,ierror,
cc     &                       z1,pg1,ro1,z2,pg2,ro2)
c	    else
c	      IRH1=0
c	      if(imodel2.eq.1)IRH2=0
c	    end if
cc	    call taulinea4(1,imodel2,deglon,deglat,vx,vy,atmos,z1,pg1,ro1,z2,pg2,ro2,ntau) 
c
cc	    call leemodi222(1,modelout1,modelout2,atmos,ntau,ierror,
cc    &                       z1,pg1,ro1,z2,pg2,ro2)
c
cc-----------------------------------------------------------------     
         end if
c               
c         icalerr=1   !=1 pues vamos a calcular errores
c         iauto=0
c         
c         do i=1,16*ntau+5
c	    atmos(i)=atmoslin(i)  !go into fperfil2err through common common/atmosfera/atmos
c	 end do
c	 	    
c     	 call marquarderr(stok,nfrecs,atmosr,mvecmax,mfit)
c                                                       
c         do j=1,iga1
c            x(j)=atmosrlin(j)*(1.+sqrt(abs(errores(j))))
c	 end do
c         do j=iga1+1,ifi11
c            x(j)=atmosrlin(j)+sqrt(abs(errores(j)))
c	 end do
c         do j=ifi11+1,iga2
c            x(j)=atmosrlin(j)*(1.+sqrt(abs(errores(j))))
c	 end do
c         do j=iga2+1,ifi22
c            x(j)=atmosrlin(j)+sqrt(abs(errores(j)))
c	 end do
c         do j=ifi22+1,mfit
c            x(j)=atmosrlin(j)*(1.+sqrt(abs(errores(j))))
c	 end do
c	  	 
c         do i=1,ntau                      !nuevo para el calculo errores
c	    pg1b(i)=pg1(i)
c	    pg2b(i)=pg2(i)
c	    z1b(i)=z1(i)
c	    z2b(i)=z2(i)
c	    ro1b(i)=ro1(i)
c	    ro2b(i)=ro2(i)
c         end do
c
c	 call amp2err(ntau,mvecmax,atmoslin,x)
c	 call taulinea4(1,imodel2,deglon,deglat,vx,vy,atmoslin,z1b,pg1b,ro1b,z2b,pg2b,ro2b,ntau) 
c
c	 do i=ntau+1,16*ntau+5
c	    atmos(i)=abs(atmos(i)-atmoslin(i))     !errores
c	 end do
c	 do i=1,ntau                      !nuevo para el calculo errores
c	    z1(i)=abs(z1(i)-z1b(i))
c	    z2(i)=abs(z2(i)-z2b(i))
c            pg1(i)=abs(pg1(i)-pg1b(i))
c	    pg2(i)=abs(pg2(i)-pg2b(i))
c            ro1(i)=abs(ro1(i)-ro1b(i))
c	    ro2(i)=abs(ro2(i)-ro2b(i))
c         end do
c
cc escribimos los errores
c	 call leemodi222(1,modelerr1,modelerr2,atmos,ntau,ierror,
c     &                       z1,pg1,ro1,z2,pg2,ro2)
c
c	     
c         close(13)

	 close(ican)

         chisnpr(ici)=chisn
         snpr(ici)=snn
	
	 open(78,file=snychi,access='append')
           write(78,*) snchi,chprint,sumsq*sn*sn/float(nfree-ielimino)
         close(78)
	 end do	!fin del do en el numero de ciclos (ici)
	 close(icanal)

c	open(icanal,file=control,access='append')
	open(icanal,file=control)
        if(iauto.eq.0)write(icanal,1000)"ite","DE","s/neqv","chi**2","Mac1","Mac2","mic1","mic2"
     &                ,"fill2","Ic1/Ic2","stray"



        if(iauto.eq.1)write(icanal,2000)"ite","DE","s/neqv","chi**2",(ca(kk),kk=1,18)

             write(icanal,786)0,chisnpr(1),snpr(1)
!esta es la correcta

        do i=1,ll
c          print*,'desire 1912 escribo en log ',i,nguvec(i),addvec(i)
 	  if(iauto.eq.0)write(icanal,1002)nguvec(i),addvec(i),snchivec(i),chprintvec(i),amac1vec(i),
     &                      amac2vec(i),amic1vec(i),amic2vec(i),
     &                      fill2vec(i),contrastevec(i),porcienvector(i)
 	  if(iauto.eq.1) write(icanal,2002)nguvec(i),addvec(i),snchivec(i),chprintvec(i),
     & (mvec(i,jj),jj=1,18)

	  do j=1,nciclos-1   
             if(i.eq.posicionciclo(j))then     
c               print*,'desire 1931 escribo en log ',i,j,posicionciclo(j)
	       if(iauto.eq.0)write(icanal,1000)"ite","DE","s/neqv","chi**2","Mac1","Mac2",
     &                "mic1","mic2","fill2","Ic1/Ic2","stray"
               if(iauto.eq.1)write(icanal,2000)"ite","DE","s/neqv","chi**2",
     &                          (ca(kk),kk=1,18)

               write(icanal,786)0,chisnpr(j+1),snpr(j+1)

	     endif
	  enddo
	enddo
	close(icanal)

      print*, ' __________________________________________________________ '
      print*, '|                                                          |'
      print'(a)',CARTEL
      print*, '|__________________________________________________________|'
      print*, ''

c     10/05/19 epm: Measure time before existing the program.
      call chrono(cpu2, wall2)
      print'(1x,a,i8)','CPU  time (ms) = ', cpu2 - cpu1
      print'(1x,a,i8)','Wall time (ms) = ', wall2 - wall1
      print*,''

c     Formats.
786   format(1x,i3,6x,1pe9.2,1x,e10.3)
1000  format(2x,a2,2x, a2,5x,a3,6x,a6,2x,  a6,1x,a6,3x,  a4,3x,a4,3x,
     &       a5,1x,a7,1x,a5)
1002  format(1x,i3,1x,f4.0,1x,1pe9.2,1x, e10.3,1x, 0pf6.3,5(1x, f6.3),1x,f6.3)
2000  format(2x,a2,2x, a2,5x,a3,6x,a6,3x,18(a3))
2002  format(1x,i3,1x,f4.0,1x,1pe9.2,1x, e10.3, 18(i3))
423   format(a,i2,a)

c     Happy ending.
      stop

c     21/06/19 epm: Error opening control file.
666   print*,'STOP: File "',trim(fcontrol),'" not found.'
      stop ' '

      end

c _______________________________________________________________________

      subroutine diagonal_FUNCIONA_PEOR(diag0,diag,isigo,iamplio,varchi,it)

      implicit real*4(a-h,o-z)
      real*4 diagon(7),diagonold(7) 
      character diver*11,espacio*35

      common/iteradiagonal/itt
      common/repeticion/factorrep !factor de diag para evitar repeticiones
      data ir/1/

      diver=' increases '
      espacio=' _________________________________ '

      if(it.eq.1)then
         ir=1
         itt=0

         diagonold(1)=diag0
         diagon(1)=diag0
         do j=2,7 
            diagonold(j)=diag
         end do
      end if
      ir=it
      if(ir .gt. 6)ir=6
      do j=1,ir
         diagonold(j)=diagon(j)
      end do 
      do j=1,ir-1
         diagon(j)=diagonold(j+1)
      end do   
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

      if(ir.gt.5.)then
         dir02=abs(diagon(ir)-diagon(ir-2))
         dir13=abs(diagon(ir-1)-diagon(ir-3))
         if(dir02 .lt. 1.e-8 .or. dir13 .lt. 1.e-8)then
            factorrep=0.
         else
            factorrep=1.
         end if 
      end if

c     if(ir.eq.7)then
c        dir01=abs(diagon(ir)-diagon(ir-1))
c        dir03=abs(diagon(ir)-diagon(ir-3))
c        dir04=abs(diagon(ir)-diagon(ir-4))
c        dir05=abs(diagon(ir)-diagon(ir-5))
c
c        if(dir01 .lt. 1.e-12 .and.dir02 .lt. 1.e-12 .and.
c     &     dir03 .lt. 1.e-12 .and.dir04 .lt. 1.e-12 .and.
c     &     dir05 .lt. 1.e-12 )factorrep=1.
c     end if

      if(diag0.gt.0)then
         if(diag.gt.diag0)then
            varchi=1.0
            add=-20.
            if(diag0.gt.0.0)add=alog10(diag0)
            write(*,1100)add,diver,espacio
            if(itt.ge.7)then
               isigo=0
               iamplio=0
            end if
         end if
      end if

      return
1100  format(5x,f4.0,1x,a11,a52)
      end

c _______________________________________________________________________

      subroutine diagonal(diag0,diag,isigo,iamplio,varchi,it)

      implicit real*4(a-h,o-z)
      real*4 diagon(100) 
      character diver*11,espacio*35

      common/iteradiagonal/itt
      common/repeticion/factorrep !factor de diag para evitar repeticiones
      data ir/1/

      diver=' increases '
      espacio=' _________________________________ '

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

      if(ir.gt.5.)then
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
            write(*,1100)add,diver,espacio
            if(itt.ge.7)then
               isigo=0
               iamplio=0
            end if
         end if
      end if

      return
1100  format(5x,f4.0,1x,a11,a52)
      end

c _______________________________________________________________
c rutina convierte
c reescribe el perfil de luz difusa en las longitudes de onda de la 
c malla
c entradas: s     - perfil de luz difusa en las ldo observadas
c           nposi - array de enteros que contiene las posiciones
c	            de la malla correspondientes a las de las ldo
c                   observadas
c           nliobs- numero total de longitudes de onda observadas
c           nli   - numero total de longitudes de onda de la malla
c salidas:  stray - (via common) perfil de luz difusa en las ldo de la malla
c Basilio 24 Junio 1994
c _______________________________________________________________

        subroutine convierte(s,nposi,nliobs,nli)

	include 'PARAMETER' !para kld
	parameter (kld4=4*kld)           
	real*4 stray(kld4),s(*)
        integer nposi(*)
        common/difusa/stray

        do i=1,nliobs
           stray(nposi(i))=s(i)
  	end do

	return
        end 

c _______________________________________________________________
c "amp2err" construye la atmosfera completa a partir de la antigua y 
c de las perturbaciones multiplicativas nuevas (atmosfera reducida) 
c 'atmos' es la atmosfera antigua completa a la salida es la atm. nueva
c 'pert'  es la atmosfera antigua reducida
c 'atmosr' es la atmosfera perturbada reducida
c a la salida atmos contendra la nueva atmosfera perturbada en todos los
c puntos, como una interpolacion por splines cubicos de la atmosfera
c perturbada en los nodos (reducida)
c
c    i variable   i  variable
c    -   tau1     -    tau2   profundidad optica a 5000 /AA
c    1   t1       9    t2     temperatura en ambos modelos (k)
c    2   p1      10    p2     presion electronica (dinas/cm**2)
c    3   mic1    11    mic2   microturbulencia (cm/s)
c    4   h1      12    h2     campo magnetico (G)
c    5   v1      13    v2     velocidad eje z (cm/s)
c    6   g1      14    g2     gamma (radianes)
c    7   f1      15    f2     fi (radianes)
c    8   mac1    16    mac2   macroturbulencia (en Km/s)
c    -   ff1     17    ff2    factor de llenado (ff1=1-ff2)
c                18    %      peso de la luz difusa
c trabajo siempre con el ff del segundo modelo
c m(i) es el numero de nodos de la varible i. 
c Asi si m(i)=0 la variable i no se modifica
c     si m(i)=1 la variable i se modifica mediante un factor mult. cte.
c     si m(i)=2 la variable i se modifica mediante un factor mult. lineal
c     .........
c     si m(i)=-1 la variable i se modifica igual que la misma variable 
c               de la otra atmosfera (es decir como i+/-8)
c En mdata(1-18) guardo los indices anteriores a la variable i (atm. ampliada)
c
c Basilio 22-3-93 
c Basilio 23-3-93 (modificacion para contemplar el caso m(i)=-1)
c Basilio y Jose Carlos 9-1-95 (modificacion perturbaciones aditivas en fi)
c Basilio y Jose Carlos 6-2-95 (modificacion perturbaciones aditivas en gamma)
c Basilio y Jose Carlos 9-1-96 (modificacion eliminacion de cotas para errores)
c _______________________________________________________________________

	subroutine amp2err(ntau,m,atmos,atmosr)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER'  !para kt
	integer m(*),mdata(18)
	real*4 atmos(*),atmosr(*)
	real*4 x(kt),y(kt),yy(kt),pert(14*kt+4),f(kt,kt)
	real*4 tau(kt),t1(kt),p1(kt),tnew1(kt),pnew1(kt)
	real*4 t2(kt),p2(kt),tnew2(kt),pnew2(kt)
	real*4 pg1(kt),z1(kt),ro1(kt),b1(kt),gam1(kt)
        real*4 pg2(kt),z2(kt),ro2(kt),b2(kt),gam2(kt)
        character*26 var(18) 
        integer icalerr,ncontpg,ipgmag
        real*4 prec,pg01,pg02,ro01,ro02
        
        common/preciso/prec      
        common/calerr/icalerr !si calculo errores=1 else =0
	common/zetas/pg1,z1,ro1,pg2,z2,ro2
        common/contornopg/ncontpg,pg01,pg02
        common/contornoro/ro01,ro02
        common/pgmag/ipgmag

	epsilon=1.e-2

c precision equilibrio hidrostatico en tanto por uno (necesaria para equisubmu)
        prec=1.e-3  ! 
     
	call comprime2(ntau,m,atmos,pert) !atmosfera antigua reducida
	
        do i=1,16
	   mdata(i)=i*ntau+2*int(i/9)  ! indi. anteri. a la var. i (ampliada)
	end do
        mdata(17)=mdata(16)+1
        mdata(18)=mdata(17)+1
   
	do i=1,ntau
	   tau(i)=atmos(i)
           t1(i)=atmos(ntau+i)
	   p1(i)=atmos(2*ntau+i)	!inicializamos la presion
           t2(i)=atmos(9*ntau+2+i)
	   p2(i)=atmos(10*ntau+2+i)	
	end do

	kred=0		!indice reducido
        kamp=ntau	!indice ampliado (los ntau puntos de tau1)

	cota=1.e10
	cotapres=2.
        cotafi=1.e10 
        cota1=-1.e10
        cota2=1.e10

	do i=1,18	!do en grupos de varibles (1=t,2=p,...etc)

           ntau2=ntau
           if(i.eq.8.or.i.eq.16.or.i.eq.17.or.i.eq.18)ntau2=1!mac1,mac2,ff2,%

           if(m(i).eq.1)then	            !si pert. constante sumo
              kred=kred+1
              if(i.eq.7.or.i.eq.15.or.i.eq.6.or.i.eq.14)then
                 y1=atmosr(kred)-pert(kred)
  	         if(y1.lt.-cotafi)y1=-cotafi   !acoto inferiormente
	         if(y1.gt.cotafi)y1=cotafi     !acoto superiormente
      	      else
                 if(abs(pert(kred)).lt.1.e-10)goto 999   !por aqui nunca va a pasar!!!!
	         y1=(atmosr(kred)/pert(kred))-1.    !perturbacion multiplicativa
  	         if(y1.lt.-cota)y1=-cota   !acoto inferiormente
	         if(y1.gt.cota)y1=cota     !acoto superiormente
	         y1=y1*pert(kred) 	            !perturbaciona aditiva
              end if

              if(i.eq.17)then        !si es el f.f
                 varfill=(1./atmos(kamp+1))-1.
                 varfill=varfill+y1  
	         atmos(kamp+1)=1./(1.+ varfill)
              else if(i.eq.18)then        !si es el %
                 varpercen=atmos(kamp+1)/(100.-atmos(kamp+1))
                 varpercen=varpercen+y1
	         atmos(kamp+1)=100.*varpercen/(1.+ varpercen)	     
	      else
                 do j=1,ntau2
                    atmos(kamp+j)=atmos(kamp+j)+y1
                 end do
              end if 
  
           else if(m(i).gt.1)then

              mm=(ntau-1)/(m(i)-1)   !espaciado entre nodos 
              kred1=kred+m(i)           !indice del ultimo nodo de cada variable
              do j=1,m(i)
                 kred=kred+1
                 jj=(j-1)*mm+1	               !indice de tau en los nodos
              	 x(j)=atmos(jj)	               !tau en los nodos 
                 if(i.eq.7.or.i.eq.15.or.i.eq.6.or.i.eq.14)then
                    y(j)=atmosr(kred)-pert(kred)
  	            if(y(j).lt.-cotafi)y(j)=-cotafi   !acoto inferiormente
	            if(y(j).gt.cotafi)y(j)=cotafi     !acoto superiormente
	         else if (i.eq.2.or.i.eq.10) then
                    y1=atmosr(kred)-pert(kred)        !perturbaciona aditiva en la presion e
  	            if(y1.lt.-cotapres)y1=-cotapres   !acoto inferiormente
	            if(y1.gt.cotapres)y1=cotapres     !acoto superiormente
                    y(j)=y1*pert(kred1)               !escala con la presion en el ultimo nodo
      	         else
	            y(j)=(atmosr(kred)/pert(kred))-1.  !pert. multiplicativa
	            if(y(j).lt.cota1)y(j)=cota1 !acoto inferiormente
	            if(y(j).gt.cota2)y(j)=cota2  !acoto superiormente
	            y(j)=y(j)*pert(kred)     !perturb. aditiva en los nodos
                 end if
	      end do
    
	      if(ntau2.ne.ntau)then
		  print*,'a la macro o al ff se les asigna mas de 1 variable? '
	          stop
              end if

	      call splines22(x,y,m(i)-2,ntau,tau,yy,f)
              
                 do j=1,ntau2
                    atmos(kamp+j)=atmos(kamp+j)+yy(j)
                 end do    
c              end if                 
	   end if	     
	   kamp=kamp+ntau2
	   if(i.eq.8)kamp=kamp+ntau+1	!los ntau puntos de tau2 y el de ff1	
	end do

c en caso de que no se corrija la presion
c ponemos las presiones en equilibrio hidrostatico con las temperaturas
c	vart1=0.0
c	vart2=0.0
        do i=1,ntau
           tnew1(i)=atmos(ntau+i)
           tnew2(i)=atmos(9*ntau+2+i)
	end do

        do i=1,ntau
           pnew1(i)=atmos(2*ntau+i)
           pnew2(i)=atmos(10*ntau+2+i)
	   b1(i)=atmos(4*ntau+i)
           b2(i)=atmos(12*ntau+2+i)
	   gam1(i)=atmos(6*ntau+i)
           gam2(i)=atmos(14*ntau+2+i)
	end do


c ----------------------------------------------------------------------------
c coloco la p1 y p2 en equilibrio hidrostatico

	   if(m(2).eq.0) then	   
	        if(ncontpg.eq.0)call equisubmu(ntau,tau,tnew1,p1,pg1,z1,ro1)
                if(ncontpg.eq.-1)then
                   pg01=ro01*83145100.*tnew1(ntau)/1.302
                   pg1(ntau)=pg01
                   ro1(ntau)=ro01
                   call pgpefromrho(tnew1(ntau),ro1(ntau),p1(ntau),pg1(ntau))
                endif    
	        if(ncontpg.ne.0)then
                   if(ipgmag.ne.1)call
     & equisubmu_cont(ntau,tau,tnew1,p1,pg01,pg1,z1,ro1)
                   if(ipgmag.eq.1)call
     & equisubmu_contmag(ntau,tau,tnew1,p1,pg01,pg1,z1,ro1,b1,gam1)
                end if
	        do i=1,ntau
                   atmos(ntau+i)=tnew1(i)
                   atmos(2*ntau+i)=p1(i)
	        end do	        
	   else
	        do i=1,ntau
                   atmos(ntau+i)=tnew1(i)
                   p1(i)= atmos(2*ntau+i)                  
                end do
                call pgzrofrompetau(ntau,tau,tnew1,p1,pg1,z1,ro1)
	   end if

	     if(m(10).eq.0) then
	        if(ncontpg.eq.0)call equisubmu(ntau,tau,tnew2,p2,pg2,z2,ro2)

                if(ncontpg.eq.-1)then
                   pg02=ro02*83145100.*tnew2(ntau)/1.302
                   pg2(ntau)=pg02
                   ro2(ntau)=ro02
                   call pgpefromrho(tnew2(ntau),ro2(ntau),p2(ntau),pg2(ntau))
                endif   

	        if(ncontpg.ne.0)then
                   if(ipgmag.ne.1)call
     & equisubmu_cont(ntau,tau,tnew2,p2,pg02,pg2,z2,ro2)
                   if(ipgmag.eq.1)call
     & equisubmu_contmag(ntau,tau,tnew2,p2,pg02,pg2,z2,ro2,b2,gam2)
                end if
                do i=1,ntau
                   atmos(9*ntau+2+i)=tnew2(i)
                   atmos(10*ntau+2+i)=p2(i)
	        end do
	     else
                do i=1,ntau
                   atmos(9*ntau+2+i)=tnew2(i)
                   p2(i)= atmos(10*ntau+2+i)
	        end do
	        call pgzrofrompetau(ntau,tau,tnew2,p2,pg2,z2,ro2)
	     end if
c ----------------------------------------------------------------------------


c Si m(i)=-1 se toma la variable de la otra atmosfera
	do i=1,16   !el ff2 no puede tener -1 (en ese caso se toma como 0) 
           ntau2=ntau
           if(i.eq.8.or.i.eq.16)ntau2=1 !si mac1,mac2 o ff2
           if(m(i).eq.-1)then
              ii=i-8
              if(ii.le.0)ii=i+8
              do j=1,ntau2 
		 atmos(mdata(i)+j)=atmos(mdata(ii)+j)
              end do
           end if
        end do 

c en cualquier caso los ff tienen que ser complementarios
        atmos(8*ntau+2)=1.-atmos(16*ntau+4)
 
	return

999	var(1)='Temperatura (modelo 1)    '
        var(2)='Presion elec.  (modelo 1) '
        var(3)='V. microturb.   (modelo 1)'
        var(4)='Campo magnetico (modelo 1)'
        var(5)='Velocidad  (modelo 1)     '
        var(6)='Inclinacion mag.(modelo 1)'
        var(7)='Azimuth campo (modelo 1)  '
        var(8)='V. macroturb.   (modelo 1)'
        var(9)='Temperatura (modelo 2)    '
        var(10)='Presion elec.  (modelo 2) '
        var(11)='V. microturb.   (modelo 2)'
        var(12)='Campo magnetico (modelo 2)'
        var(13)='Velocidad  (modelo 2)     '
        var(14)='Inclinacion mag.(modelo 2)'
        var(15)='Azimuth campo (modelo 2)  '
        var(16)='V. macroturb.   (modelo 2)'
        var(17)='Factor llenado (modelo 2) '
        var(18)='Peso % de la luz difusa '
       
        print*,' '
        print*,'********************************************************************'
        print*,'se esta tratando de invertir la variable '
        print*,var(i)
        print*,'y su valor inicial es 0 en algun punto'
        print*,'o inicializas de otra manera o no inviertas dicha variable'
        print*,'********************************************************************'
        print*,' '
        stop 'en amp2err'

	end

c _______________________________________________________________________

       subroutine coor(coorfile,deglat,deglon,Vx,Vy)
       character*(*) coorfile
       real*4 Xc,Yc,Vx,Vy,deglat,deglon
       integer icoor,ican
 
       ican=1
       icoor=meves(coorfile,100)
       if(icoor.eq. 1)then
          open(ican,file=coorfile)
          read(ican,*,err=770,END=770)Xc,Yc,Vx,Vy     !Xc, Yc=position/Rsun, Vx,VY horizontal vel cm/s(1.95e5)
          if(abs(xc) .lt. 0.0003)then
             if(xc .lt. 0.)xc=-0.0003
             if(xc .ge. 0.)xc=0.0003
          endif
c          if(abs(yc) .lt. 0.0001)then
c             if(yc .lt. 0.)yc=-0.0001
c             if(yc .ge. 0.)yc=0.0001
c          endif
          call LatLonfromXcYc(Xc,Yc,deglat,deglon)
          return
770       print*,'Incorrect format in the coordinates file',trim(coorfile)
          print*,'STOP'  
          stop
       else  
          print*,'Incorrect format in the coordinates file',trim(coorfile)
          print*,'STOP'
          stop
       end if  
       end

c _______________________________________________________________________

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

       character*(*) coorfile
       real*4 Xc,Yc,Vx,Vy,deglat,deglon,P,B0
       integer icoor,ican
 
       ican=1
       icoor=meves(coorfile,100)
       if(icoor.eq. 1)then
          open(ican,file=coorfile)
          read(ican,*,err=870,END=870)Xc,Yc,P,B0     !Xc, Yc=position/Rsun, P, B0 solar parameters in degrees
          call LatLonfromXcYcPB(Xc,Yc,P,B0,deglat,deglon,Vx,Vy)
          return
870       print*,'Incorrect format in the coordinates file',trim(coorfile)
          print*,'STOP'  
          stop
       else  
          print*,'Incorrect format in the coordinates file',trim(coorfile)
          print*,'STOP'
          stop
       end if  
       end

c ________________________________________________________________________ 

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

c _______________________________________________________________________

c evaluates the horizontal component (Vx,Vy) of the
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
     
c       print*,'deglat =',deglat, asin(sinlat)*pi180

c Following Kuker & Rudiger (2008)
c  http://iopscience.iop.org/article/10.1088/1742-6596/118/1/012029/meta
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

c       print*,' Vx=',Vx,' cm/s Vy=',Vy,' cm/s'       
       
       return
       end 
c _______________________________________________________________________
