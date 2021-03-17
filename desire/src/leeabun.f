c El fichero de abundancias debe contener una cabecera
c el numero de lineas de la cabecera no debe exceder de 20 y
c cada linea no debe tener mas de 70 caracters ( si tiene
c mas no se escribiran al modificar e fichero)
c La cabecera debe acabar necesariamente con una linea que utilice alguno
c de los simbolos (ver el 'if' de lectura en el 'do while')

        subroutine leeabun(iesc)           !iesc=1 escribo

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        parameter (na=92,narh=99)
        character mensaje*26,fichabun*100
        integer*4 in
        real*4    aa,abu(narh)   

        common/fichabun/fichabun
        common/abundances/abu

        ican=1
        mensaje=' containing the abundances' 
        call cabecera(ican,fichabun,mensaje,ifail)
        if(ifail.eq.1)goto 999

c       Los 92 elementos.
        do i=1,na
           read(ican,*,err=998)in,aa
           abu(in)=aa  !10.**(aa-12.)
        end do

c       Los 7 elementos finales que espera RH.
        do i=na+1,narh
           abu(i)=-7.96
        end do

        close(1)
        return

998     call error(KSTOP,'leeabun','The file containing the abundances'
     &  //         ' is not written correctly\n File: '//fichabun)

999     call error(KSTOP,'leeabun','The file containing the abundances'
     &  //         ' does not exist\n File: '//fichabun)

        return
        end
