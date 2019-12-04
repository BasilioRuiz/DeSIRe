c
c rnorma rutina que reescala las funciones con el paso de la integral
c ________________________________________________________________________

      subroutine rnorma(ntau,cont,r)

      implicit real*4 (a-h,o-z)
      include 'PARAMETER'   !para kt
      parameter (aln10=2.3025851)
      real*4 tau(kt),taue(kt),r(*)
      real*4 deltae(kt),deltai(kt),delt2i(kt)      !,tauRH_step(kt)

      common/segunda/tau,taue,deltae,deltai,delt2i
c     common/tauRH_step/tauRH_step

      paso= (tau(2)-tau(1))*aln10/cont
      do i=1,ntau
         r(i)=r(i)*taue(i)*paso   !*tauRH_step(i)
      end do

      return
      end
