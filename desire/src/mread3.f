c Lecturas con valores defecto por argumento
c ------------------------------------------
c lee una linea cuyos primeros characteres son un string hasta el simbolo :
c luego vienen una
c serie de n enteros (con n menor o igual a 15).
c La separacion entre los enteros se supone realizada por los simolos , o ;
c Se ignoran los blancos 
c Esta function retorna el entero
c que ocupa el orden iciclo (el icilcloesimo vamos), en caso de que iciclo sea
c mayor que n retorna el ultimo entero.
c Ignora todo aquello que se encuentre a la derecha del simbolo !
c Basilio -Enero 1993

        integer function mreadi3(ican,iciclo,idefault)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        character espacio*80,dato*300,linea*300,dat*13,d*1
        integer n(50),idefault,ican,iciclo
        character*100 msg

        read(ican,'(a)') linea
        long=len(linea)

        i=1
        do while(linea(i:i).ne.':' .and. i.lt.long)
           i=i+1
        end do
        lcom=i

        i=lcom+1
        do while(linea(i:i).ne.'!' .and. i.lt.long)
           i=i+1
        end do
        long=i

        espacio(1:lcom)=linea(1:lcom)
        leng=long-lcom
        dato(1:leng-1)=linea(lcom+1:long-1)
        dato(leng:leng)=';'
        ivar=1
        ii=0
        dat(1:13)='             '

        do 11 i=1,leng
           d=dato(i:i)
           if(d.eq.'.')then
              call error(KSTOP,'mreadi3','An integer is expected in one of'
     &        //         ' the fields of the control file.\n'
     &        //         ' A real has been found instead')
           end if
           if(d.eq.' ')goto 11
           if(d.ne.',' .and. d.ne.';')then
              ii=ii+1
              dat(ii:ii)=d
           else
              if(ii.eq.0)goto 11
              ii=0
              if(dat(1:1).eq.'*')then
                 n(ivar)=kn
              else
                read(dat,'(i13)',err=666) n(ivar)
              end if
              ivar=ivar+1
              dat(1:13)='             '
           end if
11      continue
        if(ivar.eq.1)then
           mreadi3=idefault  !si no hay nada escrito pongo el default
           goto 25
        endif
        ivar=ivar-1
        if(iciclo.ge.ivar)then
           mreadi3=n(ivar)
        else
           mreadi3=n(iciclo)
        end if

25      write(msg,*) espacio(1:lcom),mreadi3
        call error(KLITE,'',msg)

        return

666     call error(KSTOP,'mreadi3','Wrong integer format in one of the fields'
     &  //         ' of the control file')

        return
        end

c ____________________________________________________________________________
c
c lee una linea cuyos primeros characteres son un string (hasta el simbolo : )
c y luego vienen una
c serie de n real*4 (con n menor o igual a 15).
c La separacion entre los real*4 se supone realizada por los simolos , o ;
c Se ignoran los blancos 
c Esta function retorna el real*4
c que ocupa el orden iciclo (el icilcloesimo vamos), en caso de que iciclo sea
c mayor que n retorna el ultimo real*4.
c Ignora todo aquello que se encuentre a la derecha del simbolo !
c Basilio -Enero 1993

        real*4 function mreadr3(ican,iciclo,rdefault)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        integer ican,iciclo
        real*4 x(50),rdefault
        character espacio*80,dato*300,linea*300,dat*13,d*1
        character*100 msg

        read(ican,'(a)') linea
        long=len(linea)

        i=1
        do while(linea(i:i).ne.':' .and. i.lt.long)
           i=i+1
        end do
        lcom=i 

        i=lcom+1
        do while(linea(i:i).ne.'!' .and. i.lt.long)
           i=i+1
        end do
        long=i

        espacio(1:lcom)=linea(1:lcom)
        leng=long-lcom
        dato(1:leng-1)=linea(lcom+1:long-1)
        dato(leng:leng)=';'
        ivar=1
        ii=0 
        dat(1:13)='             '
        do 11 i=1,leng
           d=dato(i:i)
           if(d.eq.' ')goto 11 
           if(d.ne.',' .and. d.ne.';')then
              ii=ii+1
              dat(ii:ii)=d
           else
              if(ii.eq.0)goto 11
c             Si es un entero lo pasamos a real
              j=1
              do while(dat(j:j).ne.'.' .and. dat(j:j).ne.'e' .and.
     &                 dat(j:j).ne.'d' .and. j.le.ii)
                 j=j+1
              end do
              if(j.gt.ii)dat(ii+1:ii+1)='.'

              ii=0
c             08/08/20 epm: Primero comprobamos que leemos un real.
c             Luego volvemos a leer porque el formato "fw.d" acepta 1.e-2
c             pero no 1e-2 (sin el punto detras del 1).
              read(dat,'(f13.5)',err=676) x(ivar)
              read(dat,*) x(ivar)
              ivar=ivar+1
              dat(1:13)='             '
           end if
11      continue

        if(ivar.eq.1)then
           mreadr3=rdefault  !si no hay nada escrito pongo el default
           goto 26
        endif

        ivar=ivar-1
        if(iciclo.ge.ivar)then
           mreadr3=x(ivar)
        else
           mreadr3=x(iciclo)
        end if

        if(lcom.gt.30)then
           lcom=30
           espacio(lcom:lcom)=':'
        end if

c       09/09/20 epm: Pasamos por argumento el numero de decimales a mostrar.
c26     iw=id+7  !iw >= id+7, pues 1pe<iw>.<id> se representa +1.ddddE+ee
c       write(fmt,'(a,i2.2,a1,i1,a1)') "(1x,a30,1x,1pe",iw,".",id,")"
c       write(msg,fmt)espacio(1:lcom),mreadr3
c       call error(KLITE,'',msg)
c
c       09/09/20 epm: Finalmente no implemento lo anterior, mejor defino el
c       formato en funcion del valor del numero real.
26      if(int(mreadr3).gt.1 .and. int(mreadr3).lt.10)then
           write(msg,'(1x,a30,1x,f9.6)')espacio(1:lcom),mreadr3
        else
           write(msg,'(1x,a30,1x,1pe9.2)')espacio(1:lcom),mreadr3
        end if
        call error(KLITE,'',msg)

        return

676     call error(KSTOP,'mreadr3','Wrong real format in one of the fields'
     &  //         ' of the control file')

        return
        end
