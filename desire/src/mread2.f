c Lecturas normales
c -----------------
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

        integer function mreadi2(ican,iciclo)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        character espacio*80,dato*300,linea*300,dat*13,d*1
        integer n(50),ican,iciclo
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
              call error(KSTOP,'mreadi2','An integer is expected in one of'
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
                 n(ivar)=1000
              else
                read(dat,'(i13)',err=666) n(ivar)
              end if
              ivar=ivar+1
              dat(1:13)='             '
           end if
11      continue
        if(ivar.eq.1)then
           mreadi2=0  !si no hay nada escrito lo hago cero. LRB
           goto 25
        endif
        ivar=ivar-1
        if(iciclo.ge.ivar)then
           mreadi2=n(ivar)
        else
           mreadi2=n(iciclo)
        end if

25      write(msg,*) espacio(1:lcom),mreadi2
        call error(KLITE,'',msg)

        return

666     call error(KSTOP,'mreadi2','Wrong integer format in one of the fields'
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

        real*4 function mreadr2(ican,iciclo)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        integer ican,iciclo
        real*4 x(50)
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
           mreadr2=0.  !si no hay nada escrito lo hago cero, LRB
           goto 26
        endif

        ivar=ivar-1
        if(iciclo.ge.ivar)then
           mreadr2=x(ivar)
        else
           mreadr2=x(iciclo)
        end if

        if(lcom.gt.30)then
           lcom=30
           espacio(lcom:lcom)=':'
        end if

26      write(msg,100)espacio(1:lcom),mreadr2
        call error(KLITE,'',msg)
100     format(1x,a30,1x,1pe9.2)

        return

676     call error(KSTOP,'mreadr2','Wrong real format in one of the fields'
     &  //         ' of the control file')

        return
        end

c_____________________________________________________________________________
c
c lee una linea cuyos primeros characteres son un string (hasta el simbolo :)
c y luego vienen una
c serie de n cadenas o palabras (con n menor o igual a 15).
c La separacion entre las cadenas se supone realizada por los simolos , o ;
c Se ignoran los blancos 
c Esta function retorna la cadena
c que ocupa el orden iciclo (el icilcloesimo vamos), en caso de que iciclo sea
c mayor que n retorna la ultima cadena.
c Ignora todo aquello que se encuentre a la derecha del simbolo !
c Basilio -Enero 1993

        character*100 function mreadc2(ican,iciclo)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        integer ican,iciclo
        character*100 x(100)
        character titulo*100,dato*400,linea*400,dat*100,d*1,titprint*150

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

        titulo(1:lcom)=linea(1:lcom)
        leng=long-lcom
        dato(1:leng-1)=linea(lcom+1:long-1)
        dato(leng:leng)=';'
        ivar=1
        ii=0
        dat(1:100)='                                                                                                    '
        do 11 i=1,leng
           d=dato(i:i)
           if(d.eq.' ')goto 11
           if(d.ne.',' .and. d.ne.';')then
              ii=ii+1
              dat(ii:ii)=d
           else
              if(ii.eq.0)goto 11
              x(ivar)(1:100)=dat(1:100)
              ii=0
              ivar=ivar+1
              dat(1:100)='                                                                                                    '
           end if
11      continue
        ivar=ivar-1
        if(ivar.ne.0)then
           if(iciclo.ge.ivar)then
              mreadc2(1:100)=x(ivar)(1:100)
           else
              mreadc2(1:100)=x(iciclo)(1:100)
           end if
        else
           mreadc2(1:100)='                                                                                                    '

        endif

        titprint=titulo(1:lcom)//' '//mreadc2(1:48)
        call error(KLINE,'',titprint)

        return
        end

c_____________________________________________________________________________

        subroutine mreadmalla(ican,numblends,n,ini,paso,fin,ierror,blanco)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        real*4 x(50)
        character dato*300,linea*300,dat*13,d*1
        integer n(*),blanco,ican,numblends,ierror
        real*4 ini,paso,fin

        ierror=0
        blanco=0
        read(ican,'(a)',end=999) linea
        longg=len(linea)

        if(longg.eq.0 .or. longg.eq.1)then  !la linea esta vacia, o no esta bien escrita
           blanco=1
           return
        endif

        lcom=0

        i=lcom+1
        do while(linea(i:i).ne.':' .and. i.lt.longg)
           i=i+1
        end do
        if(i.eq.longg)then  !la linea esta vacia, o no esta bien escrita
           blanco=1
           return
        endif
        long=i

        leng=long-lcom
        dato(1:leng-1)=linea(lcom+1:long-1)
        dato(leng:leng)=';'
        ivar=1
        ii=0
        dat(1:13)='             '

        do 11 i=1,leng
           d=dato(i:i)
           if(d.eq.'.')then
              call error(KSTOP,'mreadmalla','Incorrect format in the file'
     &        //         ' containing the wavelength grid')
           end if
           if(d.eq.' ')goto 11
           if(d.ne.',' .and. d.ne.';')then
              ii=ii+1
              dat(ii:ii)=d
           else
              if(ii.eq.0)goto 11
              ii=0
              read(dat,'(i13)') n(ivar)
              ivar=ivar+1
              dat(1:13)='             '
           end if
11      continue
        if(ivar.eq.1)then
           call error(KSTOP,'mreadmalla','Incorrect format in the file'
     &     //         ' containing the wavelength grid')
        endif
        ivar=ivar-1
        numblends=ivar

c       Ahora leemos las lambdas:

        leng=longg-long
        dato(1:leng-1)=linea(long+1:longg-1)
        dato(leng:leng)=';'
        ivar=1
        ii=0
        dat(1:13)='             '

        do 12 i=1,leng
           d=dato(i:i)
           if(d.eq.' ')goto 12
           if(d.ne.',' .and. d.ne.';')then
              ii=ii+1
              dat(ii:ii)=d
           else
              if(ii.eq.0)goto 12
c             si es un entero lo pasamos a real:
              j=1
              do while(dat(j:j).ne.'.' .and. dat(j:j).ne.'e' .and.
     &                 dat(j:j).ne.'d' .and. j.le.ii)
                 j=j+1
              end do
              if(j.gt.ii)dat(ii+1:ii+1)='.'

              ii=0
              read(dat,'(f13.5)') x(ivar)
              ivar=ivar+1
              dat(1:13)='             '
           end if
12      continue

        if(ivar.gt.4)then
           call error(KWARN,'mreadmalla','There are more than the'
     &     //         ' three items required to determine the wavelength\n'
     &     //         ' range in the file containing the wavelength grid.'
     &     //         ' The first three are adopted')
        endif
        if(ivar.lt.4)then  !falta algun numero
           call error(KSTOP,'mreadmalla','Some parameter required to'
     &     //         ' determine the wavelength range is missing in the\n'
     &     //         ' file containing the wavelength grid. Check the values'
     &     //         ' after the : symbol')
        endif

        ini=x(1)
        paso=x(2)
        fin=x(3)
        
        return

999     ierror=1

        return
        end
