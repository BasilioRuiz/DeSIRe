
      subroutine mensaje(numero,mensaje1,mensaje2,mensaje3)
      implicit real*4 (a-h,o-z)

      character(*) mensaje1,mensaje2,mensaje3
      integer numero
c     common/canal/icanal

c     10/10/19 epm: Prints removing any leading and trailing blanks.

      print*,''
      print*,trim(mensaje1)
      if (numero.ge.2) print*,trim(mensaje2)
      if (numero.ge.3) print*,trim(mensaje3)
      print*,''
      print*,'_______________________________________________________________________________'

c     open(icanal,file=control,fileopt='eof')
c     write(icanal,*) ''
c     write(icanal,*) trim(mensaje1)
c     if (numero.ge.2) write(icanal,*) trim(mensaje2)
c     if (numero.ge.3) write(icanal,*) trim(mensaje3)
c     write(icanal,*) ''
c     write(icanal,*) '_______________________________________________________________________________'
c     close(icanal)

      stop
      end
