        subroutine cabecera(ican,nombre,mensaje,ifail)

        implicit real*4 (a-h,o-z)
        character nombre*(*),mensaje*(*)
        character*80 cab

        ifail=0
        ilin=0
        icha=0

        open(ican,file=nombre,status='old',err=999)

        do while(icha.le.19.and.ilin.eq.0)
           read(ican,'(a)',end=991)cab
           icha=icha+1 
           if(cab(3:6).eq.'____'.or.cab(3:6).eq.'----'.or.
     &        cab(3:6).eq.'++++'.or.cab(3:6).eq.'===='.or.
     &        cab(3:6).eq.'....'.or.cab(3:6).eq.',,,,'.or.
     &        cab(3:6).eq.'::::'.or.cab(3:6).eq.';;;;'.or.
     &        cab(3:6).eq.'""""'.or.cab(3:6).eq.'****'.or.
     &        cab(3:6).eq.'////'.or.cab(3:6).eq.'$$$$'.or.
     &        cab(3:6).eq.'####'.or.cab(3:6).eq.'!!!!'.or.
     &        cab(3:6).eq.'>>>>'.or.cab(3:6).eq.'<<<<')ilin=1
        end do

        if(icha.ge.20)then  !there is no header!!
             goto 991
        endif

        return

991     close(ican)         !no hay cabecera!!!!
        open(ican,file=nombre,status='old',err=999)

        return

999     ifail=1
        return

        end
