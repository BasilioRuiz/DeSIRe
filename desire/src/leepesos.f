c leepesos lee un fichero (vobs) de datos de pesos para los parametros de stokes, pesos
c el formato es el mismo que el de los perfiles observados

      subroutine leepesos(vobs,ist,pesos)

      implicit real*4 (a-h,o-z)
      include 'PARAMETER'
      real*4 pesos(*),dl,si(kld),sq(kld),su(kld),sv(kld)
      integer ist(4)
      character vobs*(*)

      ican=56
      open(ican,file=vobs,status='old',err=991)
      k=0
      n0=-1
      i=0
      do while(k.lt.kld)
         k=k+1
         read(ican,*,err=992,end=990)a,dl,si(k),sq(k),su(k),sv(k)
      end do
990   ntl=i
      nli=k-1
      call sfromiquv(ist,nli,si,sq,su,sv,pesos)  !cuelga de leeuve2
      close(ican)
      return

991   call error(KWARN,'leepesos','The file containing weights does not'
     & //        ' exist: '//trim(vobs)//'\n We do not introduce weights')
      return

992   call error(KSTOP,'leepesos','Incorrect format in the file containing'
     & //        ' the weights\n File: '//vobs)

      return
      end
