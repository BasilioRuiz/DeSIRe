c rutina leelineas
c datos de la linea

c Se supone que el fichero nomlineas tiene una cabecera con DOS
c lineas de menos de 80 caracteres

c       iln es el numero de la linea dentro del fichero lineas
c       atom es el simbolo del atomo,istage el estado de ionizacion,
c       wlengt la l.d.o. en angstroms
c       zeff es la correccion multiplicativa al damping (gamma6) 
c       energy potencial de excitacion del nivel inferior en ev.
c       loggf=log10(peso estadistico (g)* fuerza oscilador (f))
c       mult(i) da el spin total (spin=.5*float(mult(i)-1)),
c       design(i) da el momento angular orbital (ej:S=0,P=1,D=2,F=3...),
c       tam(i) da el momento angular total.
c       code1 representa mediante letras minusculas los momentos angulares
c       orbitales semienteros, asi p=1/2, f=3/2, h=5/2,k=7/2,m=9/2,o=11/2
c                              r=13/2,t=15/2, u=17/2,v=19/2,w=21/2
c
c  Basilio 22 Junio 1994
c  Basilio 23 Junio 1994 modificacion oam semienteros
c  Basilio 11 Enero 2018 including lower (nlow) and upper level (nup), default=0
c ...................................................................

        subroutine leelineas(ixx,iln,atom,istage,wlengtdbl,zeff,energy,loggf,
     &                       mult,design,tam,alfa,sigma,nlow,nup,gLSd1,gLSd2)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        character   atom*2
        character   msg*20,mensaje*33
        character   nomlineas*100,linea*200
        character*1 design(2)
        character*1 code1(21)
        integer*4   ixx,iln,istage,mult(2),nlow,nup
        integer*4   ich1,ich2,cpl1,cpl2
        real*4      gLSd1,gLSd2
        real*4      tam(2),loggf,wlengt,zeff,alfa,sigma,rlow,rup
        real*8      wlengtdbl,real8leido

c       10/10/21 epm: Los numeros cuanticos salen por argumento de read_*()
c       y cargan el common aqui declarado. Cuidado con los nombres porque
c       Fortran no distingue entre mayusculas y minusculas.
        real*4      dL(2),dS(2)                         !LS
        real*4      dL1(2),dS1(2),dJ1(2),dele(2),dK(2)  !JK
        real*4      j1(2),l1(2),j2(2),l2(2)             !JJ
c       04/04/21 epm: Common para albergar informacion del acoplamiento.
        integer*4   nlines_lte,cpl_low(kl),cpl_up(kl)
        real*4      qn1_low(kl),qn1_up(kl),qn2_low(kl),qn2_up(kl),
     &              qn3_low(kl),qn3_up(kl),qn4_low(kl),qn4_up(kl),
     &              qn5_low(kl),qn5_up(kl),qn6_low(kl),qn6_up(kl)
        common/cpldata/nlines_lte,cpl_low,cpl_up,
     &                 qn1_low,qn1_up,qn2_low,qn2_up,qn3_low,qn3_up,
     &                 qn4_low,qn4_up,qn5_low,qn5_up,qn6_low,qn6_up

        common/ficlineas/nomlineas
        common/numerito/num

        data nlines_lte/0/   !ademas lo hace estatico
        data code1/'p','1','f','2','h','3','k','4','m','5','o',
     &  '6','r','7','t','8','u','9','v','0','w'/ 

        ican=1
        mensaje=' containing the atomic parameters'

        call cabecera(ican,nomlineas,mensaje,ifail)
        if(ifail.eq.1)then
           call error(KSTOP,'leelineas','The file containing the atomic'
     &     //         ' parameters does not exist')
        end if

c inicializamos algunos argumentos de salida
        do i=1,2  !2 niveles, low & up
           mult(i)=0
           design(i)=char(0)
        end do

c contamos las lineas
        num=0
        nxx=0
        do while(nxx.ne.iln)
           read(ican,'(a)',err=997)linea
           call busco('=',linea,1,igu)
           read (linea(1:igu-1),*,err=998) nxx
           num=num+1
        end do
        close(ican)
        call cabecera(ican,nomlineas,mensaje,ifail)

c saltamos las num-1 primeras lineas
        do i=1,num-1
           read(ican,'(a)',err=997)linea
        end do

c leemos la linea como una cadena
        read(ican,'(a)',err=998)linea
        close(ican)

c buscamos el '=' y leemos nxx y atom
        call busco('=',linea,1,igu)
        read (linea(1:igu-1),*,err=998) nxx
        atom=linea(igu+1:igu+2)

c buscamos el proximo blanco
        call busco(' ',linea,igu,ibl1)

        istage=nteroleido(ibl1,ibl2,linea)

        if(istage.gt.2)then
           write(msg,*) istage
           call error(KSTOP,'leelineas','No ionization stages higher'
     &     //         ' than 2 can be handled.\n The ionization'
     &     //         ' stage '//trim(adjustl(msg))//' has been found'
     &     //         ' in the file\n containing the atomic parameters')
        end if
        wlengtdbl=real8leido(ibl2,ibl3,linea)
        zeff=realleido(ibl3,ibl30,linea)
        energy=realleido(ibl30,ibl4,linea)
        loggf=realleido(ibl4,ibl5,linea)

c buscamos el proximo '-'
        call busco('-',linea,ibl5,iguion)

c buscamos si existe un parentesis antes del guion. El entero esta a la
c izquierda del corchete
        call busco('(',linea,ibl5,iopenp1)
        call busco('{',linea,ibl5,iopenll1)
        call busco('(',linea,iguion+1,iopenp2)
        call busco('{',linea,iguion+1,iopenll2)
        if(iopenp1 .gt. iguion .and. iopenp1 .lt. 200)then
           iopenp2=iopenp1
           iopenp1=200
        end if
        if(iopenll1 .gt. iguion .and. iopenll1 .lt. 200)then
           iopenll2=iopenll1
           iopenll1=200
        end if

c first term
        if(iopenp1 .eq. 200 .and. iopenll1 .eq. 200)then
           call read_normal(linea,1,ibl5,200,iguion,mult(1),design(1),tam(1),
     &                      ich1,dL(1),dS(1))
           gLSd1=-1.
           cpl1=0 !term 1 LS
        else
           if(iopenp1 .gt. ibl5 .and. iopenp1.lt. iguion)then  !existe un ( antes del -
              call read_parent(linea,ibl5,iopenp1,mult(1),design(1),tam(1),
     &                         ich1,gLSd1,dL1(1),dS1(1),dJ1(1),dele(1),dK(1))
              cpl1=1 !term 1 JK
           else if(iopenll1 .gt. ibl5 .and. iopenll1.lt. iguion)then
              call read_parentJJ(linea,ibl5,iopenll1,tam(1),
     &                           ich1,gLSd1,j1(1),l1(1),j2(1),l2(1))
              cpl1=2 !term 1 JJ
           else
              call error(KSTOP,'leelineas','The file containing the atomic'
     &        //         ' parameters is not written correctly')
           end if
        end if

c second term
        if(iopenp2 .eq. 200 .and. iopenll2 .eq. 200)then
           call read_normal(linea,2,ibl5,200,iguion,mult(2),design(2),tam(2),
     &                      ich2,dL(2),dS(2))
           gLSd2=-1.
           cpl2=0 !term 2 LS
        else
           if(iopenp2 .gt. ibl5 .and. iopenp2 .gt. iguion .and. iopenp2 .lt. 200)then  !existe un ( despues del -
              call read_parent(linea,ibl5,iopenp2,mult(2),design(2),tam(2),
     &                         ich2,gLSd2,dL1(2),dS1(2),dJ1(2),dele(2),dK(2))
              cpl2=1 !term 2 JK
           else if(iopenll2 .gt. ibl5 .and. iopenll2 .gt. iguion .and. iopenll2 .lt. 200)then
              call read_parentJJ(linea,ibl5,iopenll2,tam(2),
     &                           ich2,gLSd2,j1(2),l1(2),j2(2),l2(2))
              cpl2=2 !term 2 JJ
           else
              call error(KSTOP,'leelineas','The file containing the atomic'
     &        //         ' parameters is not written correctly')
           end if
        end if

        alfa=realleido(ich2+5,ial1,linea)
        sigma=realleido(ial1,ial2,linea)
        nlow=0
        nup=0
        rlow=realleido(ial2,ial3,linea)
        nlow=nint(rlow)
        rup=realleido(ial3,ial4,linea)
        nup=nint(rup)

c       04/04/21 epm: Guardo en un common la informacion a pasar a RH.
        if (ixx .eq. 1) nlines_lte=0  !garantiza la correcta re-inicializacion

        if (nlow.eq.0 .and. nup.eq.0) then  !solo lineas NLTE
           nlines_lte=nlines_lte+1

           cpl_low(nlines_lte) = cpl1
           cpl_up(nlines_lte)  = cpl2

           if (cpl1.eq.0) then      !LS
              qn1_low(nlines_lte) = tam(1)
              qn2_low(nlines_lte) = dL(1)
              qn3_low(nlines_lte) = dS(1)
              qn4_low(nlines_lte) = 0
              qn5_low(nlines_lte) = 0
              qn6_low(nlines_lte) = 0
           else if (cpl1.eq.1) then !JK
              qn1_low(nlines_lte) = tam(1)
              qn2_low(nlines_lte) = dL1(1)
              qn3_low(nlines_lte) = dS1(1)
              qn4_low(nlines_lte) = dJ1(1)
              qn5_low(nlines_lte) = dele(1)
              qn6_low(nlines_lte) = dK(1)
           else if (cpl1.eq.2) then !JJ
              qn1_low(nlines_lte) = tam(1)
              qn2_low(nlines_lte) = j1(1)
              qn3_low(nlines_lte) = l1(1)
              qn4_low(nlines_lte) = j2(1)
              qn5_low(nlines_lte) = l2(1)
              qn6_low(nlines_lte) = 0
           end if
           if (cpl2.eq.0) then      !LS
              qn1_up(nlines_lte)  = tam(2)
              qn2_up(nlines_lte)  = dL(2)
              qn3_up(nlines_lte)  = dS(2)
              qn4_up(nlines_lte)  = 0
              qn5_up(nlines_lte)  = 0
              qn6_up(nlines_lte)  = 0
           else if (cpl2.eq.1) then !JK
              qn1_up(nlines_lte)  = tam(2)
              qn2_up(nlines_lte)  = dL1(2)
              qn3_up(nlines_lte)  = dS1(2)
              qn4_up(nlines_lte)  = dJ1(2)
              qn5_up(nlines_lte)  = dele(2)
              qn6_up(nlines_lte)  = dK(2)
           else if (cpl2.eq.2) then !JJ
              qn1_up(nlines_lte)  = tam(2)
              qn2_up(nlines_lte)  = j1(2)
              qn3_up(nlines_lte)  = l1(2)
              qn4_up(nlines_lte)  = j2(2)
              qn5_up(nlines_lte)  = l2(2)
              qn6_up(nlines_lte)  = 0
           end if
        end if

        return

c Mensajes de error

998     call error(KSTOP,'leelineas','The file containing the atomic'
     &  //         ' parameters is not written correctly')

997     write(msg,*) iln
        call error(KSTOP,'leelineas','The line '//trim(adjustl(msg))
     &  //         ' appearing in the profiles and/or wavelength grid\n does'
     &  //         ' not exist in the file containing the atomic parameters')

100     format(1x,a20,i2,1x,a59)
110     format(1x,a15,i2,1x,a43)

        return
        end

c ___________________________________________________________________

        subroutine read_normal(linea,iparte,ibl5,iopen1,iguion,
     &                         mult1,design1,tam1,ich1,dL,dS)

        include 'PARAMETER'
        character linea*(*),design1*1
        integer ibl5,iopen1,iguion,mult1,ibar1,ntero1,iint1,ich1,iparte,ii
        real*4 tam1,dL,dS
        character*1 code1(21)

        data code1/'p','1','f','2','h','3','k','4','m','5','o',
     &  '6','r','7','t','8','u','9','v','0','w'/ 

        ii=iguion
        if(iparte.eq.2)ii=200
        if(iopen1.lt.ii)then
           read (linea(iopen1-2:iopen1-1),'(i2)',err=698)mult1
           call busco('/',linea,iopen1,ibar1) 
           read (linea(iopen1+1:ibar1-1),'(i2)',err=698)ntero1
           design1=code1(ntero1)
        else             !no tenemos corchete a la izquierda del -
c el entero empieza 7 characteres antes del guion
          call busco('-',linea,ibl5,iguion) 
          iint1=iguion-7
          if(iparte.eq.2)iint1=iguion+1
          read (linea(iint1:iint1+1),'(i2)',err=698)mult1
          dS=(mult1-1)/2.
          ich1=iint1+2 
          design1=linea(ich1:ich1)
          dL=-1.
          if(design1.eq.'S')dL=0.
          if(design1.eq.'P')dL=1.
          if(design1.eq.'D')dL=2.
          if(design1.eq.'F')dL=3.
          if(design1.eq.'G')dL=4.
          if(design1.eq.'H')dL=5.
          if(design1.eq.'I')dL=6.
          if(design1.eq.'K')dL=7.
          if(design1.eq.'L')dL=8.
          if(design1.eq.'M')dL=9.
          if(design1.eq.'N')dL=10.
          if(design1.eq.'O')dL=11.
          if(design1.eq.'R')dL=12.
          if(design1.eq.'T')dL=13.
          if(design1.eq.'U')dL=14.
          if(design1.eq.'V')dL=15.
          if(design1.eq.'W')dL=16.
          if (dL.lt.0) then
             call error(KSTOP,'read_normal','The file containing the atomic'
     &       //         ' parameters is not written correctly')
          end if
        end if
        read (linea(ich1+1:ich1+4),'(f4.1)',err=698)tam1

        return

698     call error(KSTOP,'read_normal','The file containing the atomic'
     &  //         ' parameters is not written correctly')

        return
        end

c ___________________________________________________________________

        subroutine read_parent(linea,ibl5,iopenp1,mult1,design1,tam1,
     &                         ich2,gLSd,dL1,dS1,dJ1,dele,dK)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        character linea*(*),design1*1,dL1ch*1,delech*1
        integer ibl5,iopenp1,mult1,iclosep1,d2S1,ich2
        real*4 tam1,dL1,dS1,dJ1,dele,dK
        real*4 gLSd,JKcoupling

        call busco(')',linea,ibl5,iclosep1)
        read (linea(iopenp1+1:iopenp1+1),'(i2)',err=798)d2S1 !2*S1+1
c       read (linea(iopenp1+1:iopenp1+1),'(i2)')d2S1
        dS1=(d2S1-1)/2.
        dL1=-1.
        dL1ch=linea(iopenp1+2:iopenp1+2)  !L1
        if(dL1ch.eq.'S')dL1=0.
        if(dL1ch.eq.'P')dL1=1.
        if(dL1ch.eq.'D')dL1=2.
        if(dL1ch.eq.'F')dL1=3.
        if(dL1ch.eq.'G')dL1=4.
        if(dL1ch.eq.'H')dL1=5.
        if(dL1ch.eq.'I')dL1=6.
        if(dL1ch.eq.'K')dL1=7.
        if(dL1ch.eq.'L')dL1=8.
        if(dL1ch.eq.'M')dL1=9.
        if(dL1ch.eq.'N')dL1=10.
        if(dL1ch.eq.'O')dL1=11.
        if(dL1ch.eq.'R')dL1=12.
        if(dL1ch.eq.'T')dL1=13.
        if(dL1ch.eq.'U')dL1=14.
        if(dL1ch.eq.'V')dL1=15.
        if(dL1ch.eq.'W')dL1=16.
        if (dL1.lt.0) then
           call error(KSTOP,'read_parent','The file containing the atomic'
     &     //         ' parameters is not written correctly')
        end if
        read(linea(iopenp1+3:iclosep1-1),'(f4.1)',err=798)dJ1
        delech=linea(iclosep1+1:iclosep1+1)
        dele=-1.
        if(delech.eq.'s')dele=0.
        if(delech.eq.'p')dele=1.
        if(delech.eq.'d')dele=2.
        if(delech.eq.'f')dele=3.
        if(delech.eq.'g')dele=4.
        if(delech.eq.'h')dele=5.
        if(delech.eq.'i')dele=6.
        if(delech.eq.'k')dele=7.
        if(delech.eq.'l')dele=8.
        if(delech.eq.'m')dele=9.
        if(delech.eq.'n')dele=10.
        if(delech.eq.'o')dele=11.
        if(delech.eq.'q')dele=12.
        if(delech.eq.'r')dele=13.
        if(delech.eq.'t')dele=15.
        if (dele.lt.0) then
           call error(KSTOP,'read_parent','The file containing the atomic'
     &     //         ' parameters is not written correctly')
        end if

        read (linea(iclosep1+2:iclosep1+2),'(i1)',err=798)mult1
        design1=linea(iclosep1+3:iclosep1+3)
        dK=-1.
        if(design1.eq.'S')dK=0.
        if(design1.eq.'P')dK=1.
        if(design1.eq.'D')dK=2.
        if(design1.eq.'F')dK=3.
        if(design1.eq.'G')dK=4.
        if(design1.eq.'H')dK=5.
        if(design1.eq.'I')dK=6.
        if(design1.eq.'K')dK=7.
        if(design1.eq.'L')dK=8.
        if(design1.eq.'M')dK=9.
        if(design1.eq.'N')dK=10.
        if(design1.eq.'O')dK=11.
        if(design1.eq.'Q')dK=12.
        if(design1.eq.'p')dK=0.5
        if(design1.eq.'f')dK=1.5
        if(design1.eq.'h')dK=2.5
        if(design1.eq.'k')dK=3.5
        if(design1.eq.'m')dK=4.5
        if(design1.eq.'o')dK=5.5
        if(design1.eq.'r')dK=6.5
        if(design1.eq.'t')dK=7.5
        if(design1.eq.'u')dK=8.5
        if(design1.eq.'v')dK=9.5
        if(design1.eq.'w')dK=10.5
        if (dK.lt.0) then
           call error(KSTOP,'read_parent','The file containing the atomic'
     &     //         ' parameters is not written correctly')
        end if

        read(linea(iclosep1+4:iclosep1+6),'(f4.1)',err=798)tam1
        gLSd=JKcoupling(dL1,dS1,dJ1,dele,dK,tam1)

        ich2=iclosep1+3

        return

798     call error(KSTOP,'read_parent','The file containing the atomic'
     &  //         ' parameters is not written correctly')

        return
        end

c ___________________________________________________________________

        subroutine read_parentJJ(linea,ibl5,iopenp1,tam1,ich2,gLSd,
     &                           j1,l1,j2,l2)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        character linea*(*)  
        character*1 l1ch
        integer ibl5,iopenp1,iclosep1,ich2
        real*4 tam1,j1,l1,j2,l2
        real*4 gLSd,JJcoupling

        call busco('}',linea,ibl5,iclosep1)

c reading the first real after {, =j1
        call nobusco(' ',linea,iopenp1+1,inbl1)  !first no-blank char after {
        j1=read1real(linea,inbl1,ich_ifin)

c reading the next character, = l1 
        call nobusco(' ',linea,ich_ifin+1,i2)
        l1ch=linea(i2:i2)
        l1=-1.
        if(l1ch.eq.'s' .or. l1ch.eq.'S')l1=0.
        if(l1ch.eq.'p' .or. l1ch.eq.'P')l1=1.
        if(l1ch.eq.'d' .or. l1ch.eq.'D')l1=2.
        if(l1ch.eq.'f' .or. l1ch.eq.'F')l1=3.
        if(l1ch.eq.'g' .or. l1ch.eq.'G')l1=4.
        if(l1ch.eq.'h' .or. l1ch.eq.'H')l1=5.
        if(l1ch.eq.'i' .or. l1ch.eq.'I')l1=6.
        if(l1ch.eq.'k' .or. l1ch.eq.'K')l1=7.
        if(l1ch.eq.'l' .or. l1ch.eq.'L')l1=8.
        if(l1ch.eq.'m' .or. l1ch.eq.'M')l1=9.
        if(l1ch.eq.'n' .or. l1ch.eq.'N')l1=10.
        if(l1ch.eq.'o' .or. l1ch.eq.'O')l1=11.
        if(l1ch.eq.'q' .or. l1ch.eq.'Q')l1=12.
        if(l1ch.eq.'r' .or. l1ch.eq.'R')l1=13.
        if(l1ch.eq.'t' .or. l1ch.eq.'T')l1=15.
        if (l1.lt.0) then
           call error(KSTOP,'read_parentJJ','The file containing the atomic'
     &     //         ' parameters is not written correctly')
        end if

c reading the next real = j2
        call nobusco(' ',linea,i2+1,i3)
        j2=read1real(linea,i3,ich_ifin)

c reading the next character, = l2
        call nobusco(' ',linea,ich_ifin+1,i2)
        l1ch=linea(i2:i2)
        l2=-1.
        if(l1ch.eq.'s' .or. l1ch.eq.'S')l2=0.
        if(l1ch.eq.'p' .or. l1ch.eq.'P')l2=1.
        if(l1ch.eq.'d' .or. l1ch.eq.'D')l2=2.
        if(l1ch.eq.'f' .or. l1ch.eq.'F')l2=3.
        if(l1ch.eq.'g' .or. l1ch.eq.'G')l2=4.
        if(l1ch.eq.'h' .or. l1ch.eq.'H')l2=5.
        if(l1ch.eq.'i' .or. l1ch.eq.'I')l2=6.
        if(l1ch.eq.'k' .or. l1ch.eq.'K')l2=7.
        if(l1ch.eq.'l' .or. l1ch.eq.'L')l2=8.
        if(l1ch.eq.'m' .or. l1ch.eq.'M')l2=9.
        if(l1ch.eq.'n' .or. l1ch.eq.'N')l2=10.
        if(l1ch.eq.'o' .or. l1ch.eq.'O')l2=11.
        if(l1ch.eq.'q' .or. l1ch.eq.'Q')l2=12.
        if(l1ch.eq.'r' .or. l1ch.eq.'R')l2=13.
        if(l1ch.eq.'t' .or. l1ch.eq.'T')l2=15.  
        if (l2.lt.0) then
           call error(KSTOP,'read_parentJJ','The file containing the atomic'
     &     //         ' parameters is not written correctly')
        end if

c reading the next real (after }= tam1 i.e. J
        call nobusco(' ',linea,iclosep1+1,i2)
        tam1=read1real(linea,i2,ich2)
        gLSd=JJcoupling(j1,l1,j2,l2,tam1)

        return

798     call error(KSTOP,'read_parentJJ','The file containing the atomic'
     &  //         ' parameters is not written correctly')

        return
        end

c ___________________________________________________________________
c
c real8leido     - lee un real*8 si esta entre dos blancos.
c
c entradas:
c     ini        - es la posicion de un blanco anterior al real 
c     linea      - es la cadena a leer
c salidas:
c     ifi        - es la posicion del primer blanco posterior al real
c     real8leido - es el real*8 buscado 
c
c Basilio -22 de Junio 1994
c ...................................................................

        real*8 function real8leido(ini,ifi,linea)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        character linea*(*) 
        character str*13

        str='             '
c buscamos el primer 'no blanco' siguiente a ini 
        call nobusco(' ',linea,ini,inbl1) 
c buscamos el proximo 'blanco'  
        call busco(' ',linea,inbl1,ifi)
c leemos el real contenido en linea(inbl1:ifi-1)        
        str(1:ifi-inbl1)=linea(inbl1:ifi-1)
        read(str,'(d15.7)',err=998)real8leido

        return

c Mensajes de error

998     call error(KSTOP,'real8leido','The file containing the atomic'
     &  //         ' parameters is not written correctly')

        return
        end

c ___________________________________________________________________
c
c realleido      - lee un real si esta entre dos blancos.
c
c entradas:
c     ini        - es la posicion de un blanco anterior al real 
c     linea      - es la cadena a leer
c salidas:
c     ifi        - es la posicion del primer blanco posterior al real
c     realleido  - es el real buscado 
c
c Basilio -22 de Junio 1994
c ...................................................................

        real*4 function realleido(ini,ifi,linea)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        character linea*(*) 
        character str*13

        str='             '
c buscamos el primer 'no blanco' siguiente a ini 
        call nobusco(' ',linea,ini,inbl1) 
c buscamos el proximo 'blanco'  
        call busco(' ',linea,inbl1,ifi)
c leemos el real contenido en linea(inbl1:ifi-1)        
        str(1:ifi-inbl1)=linea(inbl1:ifi-1)
        read (str,'(f13.5)',err=998) realleido

        return

998     call error(KSTOP,'realleido','The file containing the atomic'
     &  //         ' parameters is not written correctly')

        return
        end

c ___________________________________________________________________
c
c nteroleido     - lee un entero si esta entre dos blancos.
c
c entradas:
c     ini        - es la posicion de un blanco anterior al entero 
c     linea      - es la cadena a leer
c salidas:
c     ifi        - es la posicion del primer blanco posterior al entero
c     nteroleido - es el entero buscado 
c
c Basilio -22 de Junio 1994
c ...................................................................

        integer*4 function nteroleido(ini,ifi,linea)

        implicit real*4 (a-h,o-z)
        character linea*(*)

c buscamos el primer 'no blanco' siguiente a ini
        call nobusco(' ',linea,ini,inbl1)
c buscamos el proximo 'blanco'
        call busco(' ',linea,inbl1,ifi)
c leemos el real contenido en linea(inbl1:ifi-1)

        read (linea(inbl1:ifi-1),'(i1)') nteroleido

        return
        end

c _______________________________________________________________

       real*4 function read1real(linea,ich_ini,ich_ifin) 

       implicit real*4 (a-h,o-z)
       integer ich_ini,ich_ifin,i,iascii
       character linea*(*)

       i=ich_ini
       iascii=iachar(linea(i:i))
       do while (iascii .ge. 46 .and. iascii  .le. 57 .and. iascii .ne. 47)
          i=i+1
          iascii=iachar(linea(i:i))
       end do
       ich_ifin=i-1
       n= ich_ifin-ich_ini+1
       if(n .eq. 1)read (linea(ich_ini:ich_ifin),'(f1.1)')read1real
       if(n .eq. 2)read (linea(ich_ini:ich_ifin),'(f2.1)')read1real      
       if(n .eq. 3)read (linea(ich_ini:ich_ifin),'(f3.1)')read1real      
       if(n .eq. 4)read (linea(ich_ini:ich_ifin),'(f4.1)')read1real
       if(n .eq. 5)read (linea(ich_ini:ich_ifin),'(f5.1)')read1real 

       return
       end 

c _______________________________________________________________
c
c busco          - da la posicion de un character que se encuentre
c                  a partir de una posicion dada en un string  
c entradas:
c     ch         - es el character buscado  
c     ini        - es la posicion a partir de la cual comienza la
c                  busqueda
c salidas:
c     ifi        - es la posicion buscada 
c
c Basilio -22 de Junio 1994
c ................................................................

        subroutine busco(ch,linea,ini,ifi)
        implicit real*4 (a-h,o-z)

        character ch*1,linea*(*)
        lon=len(linea)
        if(ini .lt. 1)ini=1
        ifi=ini
        if(lon .gt. 0)then
           i=ini
           do while(linea(i:i).ne.ch.and.i.lt.lon)
              i=i+1
           end do
           ifi=i
        end if

        return
        end

c _______________________________________________________________
c
c nobusco        - da la primera posicion en que no aparece
c                  un character dado, comenzando la busqueda
c                  a partir duna posicion dada de un string  
c entradas:
c     ch         - es el character (NO) buscado  
c     ini        - es la posicion a partir de la cual comienza la
c                  busqueda
c salidas:
c     ifi        - es la posicion buscada 
c
c Basilio -22 de Junio 1994
c ................................................................

        subroutine nobusco(ch,linea,ini,ifi)
        implicit real*4 (a-h,o-z)

        character ch*1,linea*(*)
        lon=len(linea)

        i=ini
        if(i .lt. 1)i=1
        if(i .gt. lon)i=lon

        do while(linea(i:i).eq.ch.and.i.lt.lon)
           i=i+1
        end do
        ifi=i

        return
        end

c _______________________________________________________________
c
c buscon         - da la posicion de un character que se encuentre
c                  a partir de una posicion dada en un string  
c entradas:
c     ch(1:n)    - es el substring buscado 
c
c     ini        - es la posicion a partir de la cual comienza la
c                  busqueda
c salidas:
c     ifi        - es la posicion buscada 
c
c Basilio -22 de Junio 1994
c ................................................................

        subroutine buscon(ch,n,linea,ini,ifi,ifound)
        implicit real*4 (a-h,o-z)

        character ch*(*),linea*(*)
        integer n,ini,ifi,l,ifound
        i=ini
        l=ifi
        ifound=0
        if(i .lt. 1)i=1

        do while(linea(i:i+n-1).ne.ch(1:n).and.i.le.l-n)
           i=i+1
        end do
        if(i .le. l-n)then 
          ifound=1
          ini=i
        end if

        return
        end

c _______________________________________________________________
c
c nobuscorev     - da la ultima posicion en que no aparece
c                  un character dado, comenzando la busqueda
c                  a partir duna posicion dada de un string  
c entradas:
c     ch         - es el character (NO) buscado  
c     ini        - es la posicion a partir de la cual comienza la
c                  busqueda
c salidas:
c     ifi        - es la posicion buscada 
c
c Basilio -22 de Junio 1994
c ................................................................

        subroutine nobuscorev(ch,linea,ini,ifi)
        character ch*1,linea*200 

        i=ini
        do while(linea(i:i).eq.ch.and.i.gt.1)
           i=i-1
        end do
        ifi=i

        return
        end

c _______________________________________________________________
c
c calculates Lande factor g of a level using jk-coupling, (j1l-coupling) following
c Egidio Landi & Landolfi (2005) Atomic Polarization, page 77
c example, for line Fe I 15652.873:  gJ1l=JKcoupling(2,2.5,4.5,3,3.5,4)

        real*4 function JKcoupling(L1,S1,J1,l,K,J)
        real*4 L1,S1,J1,l,K,J,gamma_Eg

        JKcoupling=1.+gamma_Eg(J,0.5,K)+gamma_Eg(J,K,0.5)*gamma_Eg(K,J1,l)*gamma_Eg(J1,S1,L1)

        return
        end

c _______________________________________________________________

c calculates Lande factor g of a level using jj-coupling, following
c Egidio Landi & Landolfi (2005) Atomic Polarization, page 77
c example, for line Si I 8550.3532:  gjj=JJcoupling(1.5,1,0.5,0.,1.)=1.167

        real*4 function JJcoupling(j1,l1,j2,l2,J)
        real*4 j1,l1,j2,l2,J,gamma_Eg

        JJcoupling=1.+gamma_Eg(J,j1,j2)*gamma_Eg(j1,0.5,l1)+gamma_Eg(J,j2,j1)*gamma_Eg(j2,0.5,l2)

        return
        end

c _______________________________________________________________

        real*4 function gamma_Eg(A,B,C)

        include 'PARAMETER'
        real*4 A,B,C

        if(abs(A) .lt. 0.1)then
           call error(KSTOP,'gamma_Eg','Error in evaluation for JK coupling'
     &     //         ' because one angular momentun is 0\n Revise the atomic'
     &     //         ' configuration of a JK case in the atomic data file')
        else  
           gamma_Eg=(A*(A+1)+B*(B+1)-C*(C+1))/(2*A*(A+1))
        end if

        return
        end

c _______________________________________________________________
