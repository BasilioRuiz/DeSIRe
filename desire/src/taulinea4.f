c
c taulinea4
c il=0 :pasa de una atmosfera modelo en direccion z a una
c       atmosfera en la linea de vision
c il=1 :transforma una atmosfera en la linea de vision a una atmosfera
c       en direccion z
c ntau :numero de puntos en tau
c tau  :log10(tau)
c ________________________________________________________________________

      subroutine taulinea4(il,imodel2,deglon,deglat,vx,vy,atmos,z1,pg1,ro1,z2,pg2,ro2,ntau)

      implicit real*4 (a-h,o-z)
      include 'PARAMETER'
      parameter (kt8=8*kt+2)
      integer il,ntau,imodel2
      real*4 atmos(*),atmos1(kt8),atmos2(kt8),cth,vx,vy,deglon,deglat
      real*4 z1(*),pg1(*),ro1(*),z2(*),pg2(*),ro2(*)

      common/truecth/cth !para pasarselo a RH

      pi=3.14159265
      cth=cos(deglon*pi/180.)*cos(deglat*pi/180.)

      do i=1,8*ntau+2
         atmos1(i)=atmos(i)
      end do
      if(imodel2 .ne. 0)then
         do i=1,8*ntau+2
            atmos2(i)=atmos(i+8*ntau+2)
         end do
      endif

c     if(abs(cth-1) .gt. 1.e-4  .and. il.eq. 0)then
c        call error(KWARN,'taulinea4','We suppose input/output models'
c    &   //         ' are radial. LOS is considered (cth = 1)')
c     end if
c     if(abs(cth-1) .gt. 1.e-4  .and. il.eq. 1)then
c        call error(KWARN,'taulinea4','We suppose input/output models'
c    &   //         ' are radial. LOS is considered (cth = 1)')
c     end if

      call taulinea4sub(il,cth,deglon,deglat,vx,vy,atmos1,z1,pg1,ro1,ntau)

      if(imodel2 .ne. 0)call taulinea4sub(il,cth,deglon,deglat,vx,vy,atmos2,z2,pg2,ro2,ntau)

      do i=1,8*ntau+2
         atmos(i)=atmos1(i)
      end do

      if(imodel2 .ne. 0)then
         do i=1,8*ntau+2
            atmos(i+8*ntau+2)=atmos2(i)
         end do
      endif

      return
      end
