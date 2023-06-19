c
c File          _______ : interface_fns.f
c Description   _______ : Interface functions.
c Project       _______ : DeSIRe
c Creation date _______ : 12/12/22
c Author        _______ : epm@iac.es
c


c_____________________________________________________________________________
c
c  Function: template()
c
c  Set a value as a function(x, y).
c  @param  x value in the x axis.
c  @param  y value in the y axis.
c  @param  f established value.
c  @return return code (0 if successful).
c_____________________________________________________________________________

      integer*4 function template( x, y, f )

      implicit none
      real*4 x, y, f

      f = x * y

      template = 0
      return
      end


c_____________________________________________________________________________
c
c  Function: Pg_T_Pe()
c
c  Calculate the gaseous pressure (Pg) from the temperature (T) and
c  the electronic  pressure (Pe).
c  @param  x temperature.
c  @param  y logarithm of the electronic pressure.
c  @param  f logarithm of the gaseous pressure.
c  @return return code (0 if successful).
c_____________________________________________________________________________

      integer*4 function Pg_T_Pe( x, y, f )

      implicit none
      character s1*20, s2*20
      real*4    x, y, f
      real*4    T, pg, pe, logpe, pp(10)

      T = x      !axis X = T
      logpe = y  !axis Y = log(Pe)

      pe = 10**(logpe)
      call gasccota(T, pe, pg, pp)
      f = alog10(pg)  !table value = log(Pg)

c     By definition, NaN is not equal to anything, even itself.
      if (f .ne. f) then
         write(s1, *) T
         write(s2, *) logpe
         write(*,*) ' ** NaN with  T = '//trim(adjustl(s1))
     &   //         '  and  log(Pe) = '//trim(adjustl(s2))
         Pg_T_Pe = -1
      else
         Pg_T_Pe = 0
      end if

      return
      end


c_____________________________________________________________________________
