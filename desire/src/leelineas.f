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

        subroutine leelineas(iln,atom,istage,wlengtdbl,zeff,energy,
     &                       loggf,mult,design,tam,alfa,sigma,nlow,nup,gLSd1,gLSd2)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        real*4 dS1,dL1,dele1,dK,JKcoupling,gLSd1,gLSd2  !ESTO gLSd1,gLSd2 TIENE QUE SALIR
        character*1 dL1ch,dele1ch
        integer mult(2),nlow,nup
        real*4 tam(2),loggf,wlengt,zeff,alfa,sigma,rlow,rup
        real*8 wlengtdbl,real8leido
        character design(2)*1,atom*2,linea*200  
        character*1 code1(21)
        character*100 nomlineas
        character*33 mensaje
        character*20 msg
        common/ficlineas/nomlineas
        common/numerito/num
        data code1/'p','1','f','2','h','3','k','4','m','5','o'
     & ,'6','r','7','t','8','u','9','v','0','w'/        

        ican=1
        mensaje=' containing the atomic parameters' 

        call cabecera(ican,nomlineas,mensaje,ifail)
        if(ifail.eq.1)goto 999

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
        atom=linea(igu+1:igu+3)

c buscamos el proximo blanco 
        call busco(' ',linea,igu,ibl1) 

        istage=nteroleido(ibl1,ibl2,linea)

        if(istage.gt.2)then
           write(msg,'(i1)') istage
           call error(KSTOP,'leelineas','No ionization stages higher'
     &     //         ' than 2 can be handled.\n'
     &     //         ' The ionization stage '//trim(msg)//' has been found'
     &     //         ' in the file\n containing the atomic parameters')
        endif
        wlengtdbl=real8leido(ibl2,ibl3,linea)
        zeff=realleido(ibl3,ibl30,linea)
        energy=realleido(ibl30,ibl4,linea)
        loggf=realleido(ibl4,ibl5,linea)

c buscamos el proximo '-'. 
        call busco('-',linea,ibl5,iguion) 

c buscamos si existe un parentesis antes del guion.El entero esta a la izquierda del corchete
        call busco('(',linea,ibl5,iopenp1) 
        call busco('(',linea,iguion+1,iopenp2) 

        if(iopenp1 .eq. 200)then 
           call busco('[',linea,ibl5,iopen1) !'No existen (: Version clasica'
           iparte=1
           call read_normal(linea,iparte,ibl5,iopen1,iguion,mult(1),design(1),tam(1),ich1)  
           call busco('[',linea,iguion,iopen2)          
           iparte=2
           call read_normal(linea,iparte,ibl5,iopen2,iguion,mult(2),design(2),tam(2),ich2)         
           gLSd1=-1.
           gLSd2=-1.
        else
           if(iopenp1 .gt. ibl5 .and. iopenp1.lt. iguion .and. iopenp2 .eq. 200 )then   !existe un ( solo antes del -
              call read_parent(linea,ibl5,iopenp1,mult(1),design(1),tam(1),ich1,gLSd1)
              call busco('[',linea,iguion,iopen2)          
              iparte=2
              call read_normal(linea,iparte,ibl5,iopen2,iguion,mult(2),design(2),tam(2),ich2)
              gLSd2=-1.
           end if
           if(iopenp1 .gt. ibl5 .and. iopenp1.gt. iguion .and. iopenp1.lt.200)then !existe un ( solo despues del -
              call busco('[',linea,ibl5,iopen1) 
              iparte=1
              call read_normal(linea,iparte,ibl5,iopen1,iguion,mult(1),design(1),tam(1),ich1)  
              call read_parent(linea,ibl5,iopenp1,mult(2),design(2),tam(2),ich2,gLSd2) 
              gLSd1=-1.
           end if
           if(iopenp2 .gt. iopenp1 .and. iopenp2 .lt. 200)then !'existen dos ( uno antes y otro despues del -'
              call read_parent(linea,ibl5,iopenp1,mult(1),design(1),tam(1),ich1,gLSd1)
              call read_parent(linea,iguion,iopenp2,mult(2),design(2),tam(2),ich2,gLSd2)            
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

        return

c Mensajes de error

999     call error(KSTOP,'leelineas','The file containing the atomic'
     &  //         ' parameters does not exist')

998     call error(KSTOP,'leelineas','The file containing the atomic'
     &  //         ' parameters is not written correctly')

997     write(msg,'(i2)') iln
        call error(KSTOP,'leelineas','The line '//trim(msg)//' appearing'
     &  //         ' in the profiles and/or wavelength grid\n does not exist'
     &  //         ' in the file containing the atomic parameters')

100     format(1x,a20,i2,1x,a59)
110     format(1x,a15,i2,1x,a43)

        return
        end

c----------------------------------------------------------------------------------        

        subroutine read_normal(linea,iparte,ibl5,iopen1,iguion,mult1,design1,tam1,ich1)          

        include 'PARAMETER'
        character linea*(*),design1*1
        integer ibl5,iopen1,iguion,mult1,ibar1,ntero1,iint1,ich1,iparte,ii
        real*4 tam1
        character*1 code1(21)
        data code1/'p','1','f','2','h','3','k','4','m','5','o'
     & ,'6','r','7','t','8','u','9','v','0','w'/ 
     
        ii=iguion
        if(iparte.eq.2)ii=200
        if(iopen1.lt.ii)then
           read (linea(iopen1-2:iopen1-1),'(i2)',err=698)mult1
           call busco('/',linea,iopen1,ibar1) 
           read (linea(iopen1+1:ibar1-1),'(i2)',err=698)ntero1
           design1=code1(ntero1)
        else             !no tenemos corchete a la izquierda del -
c El entero empieza 7 characteres antes del guion
          call busco('-',linea,ibl5,iguion) 
          iint1=iguion-7
          if(iparte.eq.2)iint1=iguion+1
          read (linea(iint1:iint1+1),'(i2)',err=698)mult1
          ich1=iint1+2 
          design1=linea(ich1:ich1)
        end if
        read (linea(ich1+1:ich1+4),'(f4.1)',err=698)tam1
        return
        
698     call error(KSTOP,'read_normal','The file containing the atomic'
     &  //         ' parameters is not written correctly (normal term)')

        return
        end 

c----------------------------------------------------------------------------------        

        subroutine read_parent(linea,ibl5,iopenp1,mult1,design1,tam1,ich2,gLSd)
        
        include 'PARAMETER'
        character linea*(*),design1*1,dL1ch*1,dele1ch*1
        integer ibl5,iopenp1,mult1,iclosep1,d2S1,ich2
        real*4 tam1,gLSd,dS1,dL1,dJ1,dele1,dK,dJ,JKcoupling
        
        call busco(')',linea,ibl5,iclosep1) 
        read (linea(iopenp1+1:iopenp1+1),'(i2)',err=798)d2S1 !2*S1+1
c       read (linea(iopenp1+1:iopenp1+1),'(i2)')d2S1
        dS1=(d2S1-1)/2.
        dL1=0.
        dL1ch=linea(iopenp1+2:iopenp1+2)    !L1
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
        read(linea(iopenp1+3:iclosep1-1),'(f4.1)',err=798)dJ1
        dele1ch=linea(iclosep1+1:iclosep1+1)
        if(dele1ch.eq.'s')dele1=0.
        if(dele1ch.eq.'p')dele1=1.
        if(dele1ch.eq.'d')dele1=2.
        if(dele1ch.eq.'f')dele1=3.
        if(dele1ch.eq.'g')dele1=4.
        if(dele1ch.eq.'h')dele1=5.
        if(dele1ch.eq.'i')dele1=6.
        if(dele1ch.eq.'k')dele1=7.
        if(dele1ch.eq.'l')dele1=8.
        if(dele1ch.eq.'m')dele1=9.
        if(dele1ch.eq.'n')dele1=10.
        if(dele1ch.eq.'o')dele1=11.
        if(dele1ch.eq.'q')dele1=12.
        if(dele1ch.eq.'r')dele1=13.
        if(dele1ch.eq.'t')dele1=15.
              
        read (linea(iclosep1+2:iclosep1+2),'(i1)',err=798)mult1
        design1=linea(iclosep1+3:iclosep1+3)  
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
                                 
        read(linea(iclosep1+4:iclosep1+6),'(f4.1)',err=798)tam1          
        dJ=tam1             
        gLSd=JKcoupling(dL1,dS1,dJ1,dele1,dK,dJ)
        
        ich2=iclosep1+3
        return
        
798     call error(KSTOP,'read_parent','The file containing the atomic'
     &  //         ' parameters is not written correctly (normal term)')

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

c Mensajes de error

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

c ................................................................
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

c ................................................................
c
c calculates Lande factor g of a level using J-K coupling, (J-l1 coupling) following
c Egidio Landi & Landolfi (2005) Atomic Polarization, page 77
c example, for line Fe I 15652.873:  gJ1l=JKcoupling(2,2.5,4.5,3,3.5,4)

        real*4 function JKcoupling(L1,S1,J1,l,K,J)
        real*4 L1,S1,J1,l,K,J,gamma_Eg
        
        JKcoupling=1.+gamma_Eg(J,0.5,K)+gamma_Eg(J,K,0.5)*gamma_Eg(K,J1,l)*gamma_Eg(J1,S1,L1)

        return
        end

c ................................................................

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

c-----------------------------------------------------------------

        
