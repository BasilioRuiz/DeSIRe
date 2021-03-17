        subroutine chachi(num,cha,nc)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

c       nc: Dimension de "cha" en la rutina que llama a "chachi".
c       "cha" se rellena por la izquierda con blancos.
c       La rutina chachi convierte el entero num en un character cha    

        character*(*) cha
        integer num,k,j,nc

        if (nc.gt.10) then
           call error(KSTOP,'chachi','The number of characters should be'
     &     //         ' lower than 10')
        end if        
        
        kk=10**nc 
        k=num-(num/kk)*kk      !k es el entero formado por las ultimas
                               !nc cifras de num
        do i=nc-1,0,-1
           j=10**i
           n=k/j
           k=k-j*n
           cha(nc-i:nc-i)=char(n+48)
        end do

        return
        end

c *****************************************************************************

c Quita lo que haya detras del punto '.'
c Ejemplo: "model_1.mod" --> cha="model_1"

        subroutine quitaex0(cha)

        implicit real*4 (a-h,o-z)
        character*(*) cha
        character*100 chacho

        long=len(cha)

        i=1

        do while(i.lt.long.and.cha(i:i).ne.'.')
           chacho(i:i)=cha(i:i)
           i=i+1
        end do
        
        do j=i,long
           chacho(j:j)=' '
        end do

        cha(1:long)=chacho(1:long)

        return
        end

c *****************************************************************************

c Quita lo que haya detras del punto '.' o del '_', lo que encuentre primero
c Ejemplo: "model_1.mod" --> cha="model"

        subroutine quitaex(cha)

        implicit real*4 (a-h,o-z)
        character*(*) cha
        character*100 chacho
        integer long,i,j

        long=len(cha)

        i=1

        do while(i.lt.long.and.cha(i:i).ne.'.'.and.cha(i:i).ne.'_')
           chacho(i:i)=cha(i:i)
           i=i+1
        end do
        
        do j=i,long
           chacho(j:j)=' '
        end do

        cha(1:long)=chacho(1:long)

        return
        end

c *****************************************************************************

c Quita lo que haya detras del punto '.' o del '_', lo que encuentre primero
c Ademas devuelve el numero entero que hay detras del '_'
c Ejemplo: "model_1.mod" --> cha="model", ixt=1

        subroutine quitaex2(cha,ixt)

        implicit real*4 (a-h,o-z)
        integer ixt,long,i,j
        character*(*) cha
        character*100 chacho
        
        long=len(cha)

        i=1

        do while(i.lt.long.and.cha(i:i).ne.'.'.and.cha(i:i).ne.'_')
           chacho(i:i)=cha(i:i)
           i=i+1
        end do

        if(cha(i:i).eq.'_')then
             ixt=ichar(cha(i+1:i+1))-48
        else
             ixt=0
        endif
        
        do j=i,long
           chacho(j:j)=' '
        end do

        cha(1:long)=chacho(1:long)

        return
        end

c *****************************************************************************

c Quita lo que haya detras del punto '.'
c Ademas devuelve el numero entero que hay detras del '_'
c Ejemplo: "model_1.mod" --> cha="model_1", ixt=1

        subroutine quitaex3(cha,ixt)

        implicit real*4 (a-h,o-z)
        integer ixt,long,i,j
        character*(*) cha
        character*100 chacho    

        long=len(cha)

        i=1

        do while(i.lt.long.and.cha(i:i).ne.'.')
           chacho(i:i)=cha(i:i)
           i=i+1
        end do

        if(cha(i:i).eq.'_')then
             ixt=ichar(cha(i+1:i+1))-48
        else
             ixt=0
        endif
        
        do j=i,long
           chacho(j:j)=' '
        end do

        cha(1:long)=chacho(1:long)

        return
        end

c *****************************************************************************

        subroutine extension3(cha,ext,n)  ! returns the extension of a file after . ot -

        implicit real*4 (a-h,o-z)
        character*(*) cha,ext
        integer long,i,n,j
        
        long=len(cha)

        i=1
        do while(i.lt.long.and.cha(i:i).ne.'.')
           i=i+1
        end do
        j=i+1
        do while(j.lt.long.and.cha(j:j).ne.' ')
           j=j+1
        end do
        
        ext(1:j-i)=cha(i+1:j)
        n=j-i-1

        return
        end

c *****************************************************************************

        subroutine concatena(cha1,cha2,cha)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        character*(*) cha1
        character*(*) cha2
        character*(*) cha
        integer long1, long2, long3,i,n1,n2

        long1 = len(cha1)
        long2 = len(cha2)
        long3 = len(cha)

        i=1
        do while(i.lt.long1.and.cha1(i:i).ne.' ')
           i=i+1
        end do
        if(i .eq. long1 .and. cha1(i:i).ne.' ')i=i+1
        n1=i-1

        if (long3.le.n1) then
           call error(KSTOP,'concatena','Two strings cannot be merged in one'
     &     //         ' shorter')
        end if
        
        i=1
        do while(i.lt.long2.and.cha2(i:i).ne.' ')
           i=i+1
        end do
        if(i .eq. long2 .and. cha2(i:i).ne.' ')i=i+1
        n2=i-1
       
        cha(1:n1)=cha1(1:n1)
        if (n1+n2.le.long3) then
           cha(n1+1:n1+n2) = cha2(1:n2)
           cha(n1+n2+1:long3) = ' ' 
        else
           cha(n1+1:long3) = cha2(1:long3-n1)
           call error(KWARN,'concatena','The second string '//trim(cha)
     &     //         ' has been cut back')
        end if

        return
        end

c *****************************************************************************

        subroutine quitaend0(cha,n) !quita la extension de blancos al final

        implicit real*4 (a-h,o-z)
        integer i,n 
        character*(*) cha

        long=len(cha)

        i=long

        do while(i.gt.1.and.cha(i:i).eq.' ')
           i=i-1
        end do
        n=i

        return
        end

c *****************************************************************************

        subroutine copying(model_multi1,RH_model)
        
        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        character*(*) model_multi1,RH_model
        character*100 modelcp1,modelcp2
        integer n,n1,istatus,system
        
        modelcp1='cp '//model_multi1
        call quitaend0(modelcp1,n)          
        modelcp1=modelcp1(1:n)//' '
        call quitaend0(RH_model,n1)
        modelcp2=modelcp1(1:n+1)//RH_model(1:n1)
        istatus=system(modelcp2(1:n+1+n1))
        if(istatus .ne. 0)then
          call error(KSTOP,'copying',
     &               'Error copying '//trim(model_multi1)//' in '//RH_model)
        end if  

        return
        end

c *****************************************************************************

c 10/10/20 epm: Subrutina para convertir un string en un array de enteros
c con los codigos ASCII correspondientes. Es usada para pasar strings a C,
c que como no es recomendado, lo hacemos pasando sus codigos ASCII.
c
c string = (input) cadena que queremos convertir.
c nmax   = (input) maximo numero de caracteres a convertir
c codes  = (output) array con los enteros ASCII de la cadena.

      subroutine toascii(string, nmax, codes)

      implicit none

      character ascii_table(126)
      character string*(*), car*1
      integer*4 nmax, nlen, i, j
      integer*4 codes(nmax)

c     Caracteres ASCII del 1 al 126 auqnue los 31 primeros no son imprimibles.
      data ascii_table /' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
     &                  ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
     &                  ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
     &                  ' ',' ','!','"','#','$','%','&',"'",'(',
     &                  ')','*','+',',','-','.','/','0','1','2',
     &                  '3','4','5','6','7','8','9',':',';','<',
     &                  '=','>','?','@','A','B','C','D','E','F',
     &                  'G','H','I','J','K','L','M','N','O','P',
     &                  'Q','R','S','T','U','V','W','X','Y','Z',
     &                  '[','\',']','^','_','`','a','b','c','d',
     &                  'e','f','g','h','i','j','k','l','m','n',
     &                  'o','p','q','r','s','t','u','v','w','x',
     &                  'y','z','{','|','}','~'/

      nlen = len(string)
      do i = 1, nlen
         if (i .gt. nmax) exit
         car = string(i:i)

c        Siguiendo estrictamente Fortran 77.
c        do j = 1, 126
c           if (car .eq. ascii_table(j)) then
c              codes(i) = j
c              exit
c           end if
c        end do

c        A partir de Fortran 95.
         codes(i) = iachar(car)

      end do

      return
      end

c *****************************************************************************
