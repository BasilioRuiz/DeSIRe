c leeuveobsindic lee un fichero de datos de parametros de stokes :vobs

        subroutine leeuveobsindic(vobs,ist,ntl,nli,nlin,npas,dl,nble,stok)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        real*4 stok(*),dl(*),si(kld),sq(kld),su(kld),sv(kld)
        integer npas(*),ist(*),nble(*),nlin(*),nli
        character vobs*(*)
        common/contraste/contr

        ican=56

        open(ican,file=vobs,status='old',err=991)
        k=0
        n0=-1
        i=0
        do while(k.lt.kld)
           k=k+1
           read(ican,*,err=992,end=990)a,dl(k),si(k),sq(k),su(k),sv(k)
           n1=nint(a)
           if(n1.ne.n0)then
             j=1
             i=i+1

             if(i.gt.kl) then  !si el numero de lineas es mayor que kl
                 call error(KSTOP,'leeuveobsindic','The number of lines in'
     &           //         ' the observed/stray light profiles is larger'
     &           //         ' than\n the current limit. Decrease this number'
     &           //         ' or change the PARAMETER file')
             end if

             nlin(i)=n1  !para eliminar la malla si es necesario
           end if
           n0=n1
           npas(i)=j
           j=j+1
        end do
990     ntl=i
        nli=k-1  !este es el numero de longitudes de onda

        if(nli.gt.kld)then
           call error(KSTOP,'leeuveobsindic','The number of wavelengths in'
     &     //         ' the observed/stray light profiles is larger than\n'
     &     //         ' the current limit kld. Decrease this number or'
     &     //         ' change the PARAMETER file')
        end if

        do j=1,ntl
           nble(j)=1  !si no tenemos malla, no tenemos blends
           do jj=j+1,ntl-1
             if(nlin(j).eq.nlin(jj+1))goto 993
           end do
        end do

        if(contr.gt.-1.e5)then
           nli=nli+1
           ntl=ntl+1
           npas(ntl)=1
           si(nli)=contr
           sq(nli)=0.
           su(nli)=0.
           sv(nli)=0.
           dl(nli)=0.
        end if

        call sfromiquv(ist,nli,si,sq,su,sv,stok) !cuelga de leeuve2

        close(ican)
        return

c       Mensajes de error:

991     call error(KSTOP,'leeuveobsindic','The file containing the'
     &  //         ' observed/stray light profiles does not exist\n'
     &  //         ' File: '//vobs)

992     call error(KSTOP,'leeuveobsindic','Incorrect format in the file'
     &  //         ' containing the observed/stray light profiles\n'
     &  //         ' File: '//vobs)

993     call error(KSTOP,'leeuveobsindic','The wavelength samples of a given'
     &  //         ' line have to be written contiguously\n in the file'
     &  //         ' containing the observed/stray light profiles')

        return
        end

c-----------------------------------------------------------------------------

c leeI lee un fichero de datos de parametros de stokes :devuleve solo Stokes I
        
        subroutine leeI(vobs,ntl,nli,npas,dl,si)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        real*4 dl(*),si(*)
        integer npas(*)
        character vobs*(*)
        common/contraste/contr

        ican=56

        open(ican,file=vobs,status='old',err=791)
        k=0
        n0=-1
        i=0
        do while(k.lt.kld)
           k=k+1
           read(ican,*,err=792,end=790)a,dl(k),si(k)  !,sq,su,sv
           n1=nint(a)
           if(n1.ne.n0)then
              j=1
              i=i+1
           endif
           n0=n1
           npas(i)=j
           j=j+1
        end do
790     ntl=i
        nli=k-1

        if(contr.gt.-1.e5)then
           nli=nli+1
           ntl=ntl+1
           npas(ntl)=1
           si(nli)=contr
           dl(nli)=0.
        end if

        close(ican)
        return

791     call error(KSTOP,'leeI','The file containing the PSF in .per format'
     &  //         ' does not exist\n File: '//vobs)

792     call error(KSTOP,'leeI','Incorrect format in the file containing'
     &  //         ' the PSF in .per format\n File: '//vobs)

        return
        end
