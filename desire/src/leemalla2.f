c rutina leemalla2
c lee el fichero de control de lineas y longitudes de onda:'mallaobs'
c ntau : numero de puntos en tau
c tau  : log10(tau)
c ntl  : numero total de lineas
c nlin : indice de cada linea y de cada blend
c npas : numero de puntos en cada linea
c dlamda:cada delta de l.d.o. en ma
c nble :numero de blends de cada linea

        subroutine leemalla2(mallaobs,ntl,nli,nliobs,nlin,npasobs,
     &                       npas,dlamdaobs,dlamda,nble,nposi,indice)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        integer npasobs(*),npas(kl),nble(kl),nd(kl),nposi(*),nlin(kl),ndd(kl,kl),nblee(kl)
        real*4 dlamdaobs(*),dlamda(kld),dini,dipa,difi,dinii(kl),dipaa(kl),difii(kt)
        integer indice(*),nddd(50),blanco      !,epsiloni
        character mallaobs*(*),mensajito*31
        character msg1*20,msg2*20
        common/nciclos/nciclos
        common/contraste/contr

        if(contr.gt.-1.e5.and.nciclos.gt.0)then
           ntl=ntl-1
           nliobs=nliobs-1
        end if

        ican=57
        mensajito=' containing the wavelength grid'
        call cabecera(ican,mallaobs,mensajito,ifail)

        if(ifail.eq.1)then
           call error(KSTOP,'leemalla2','The file containing the wavelength'
     &     //         ' grid does not exist\n File: '//mallaobs)
        end if

        jj=0
        k=0
        nli=0
        numlin=0
        do i=1,1000
           call mreadmalla(ican,ntonto,nddd,ddinii,ddipaa,ddifii,ierror,blanco)
           if(ierror.eq.1)goto 8   !si he llegado al final del fichero
           if(blanco.ne.1)then
              numlin=numlin+1
              if(numlin.gt.kl.or.ntonto.gt.kl)then
                 write(msg1,*) min(numlin,ntonto)
                 write(msg2,*) kl
                 call error(KSTOP,'leemalla2','The number of lines in the'
     &           //         ' wavelength grid ('//trim(adjustl(msg1))//') is'
     &           //         ' larger than kl = '//trim(adjustl(msg2))//'\n'
     &           //         ' Decrease the number of lines or change'
     &           //         ' the PARAMETER file')
              end if
              nblee(numlin)=ntonto
              dinii(numlin)=ddinii
              dipaa(numlin)=ddipaa
              difii(numlin)=ddifii

              do j=1,nblee(numlin)
                 ndd(numlin,j)=nddd(j)
              end do
           end if   !el de blanco.ne.1
        end do

8       close(ican)
        ntl=numlin

        if(ntl.eq.0)then
           call error(KSTOP,'leemalla2','Incorrect format in the file'
     &     //         ' containing the wavelength grid\n File: '//mallaobs)
        end if

c       Ahora las ordenamos como en los perfiles:

        do i=1,ntl     !indice en los perfiles 
           icheck=0
           do l=1,numlin   !indice en la malla
               if(indice(i).eq.ndd(l,1).and.icheck.lt.1)then   !he encontrado una
                   icheck=icheck+1
                   nble(i)=nblee(l)
                   do ll=1,nblee(l)
                       nd(ll)=ndd(l,ll)
                   end do
                   dini=dinii(l)
                   dipa=dipaa(l)
                   difi=difii(l)

                   do j=1,nble(i)
                      jj=jj+1
                      nlin(jj)=nd(j)
                   end do
                   ntlblends=jj
                   pasmin=1.e30
                   primero=dlamdaobs(k+1) 
                   ultimo=dlamdaobs(k+npasobs(i)) 
                   do j=1,npasobs(i)-1
                        k=k+1
                        pasoobs=dlamdaobs(k+1)-dlamdaobs(k) 
                        pasmin=min(pasmin,pasoobs) !cambio amin0 por min
                   end do

                   k=k+1
c                  if(dipa.gt.pasmin)dipa=pasmin      !provisional!!!!
                   if(dini.gt.primero)dini=primero
                   if(difi.lt.ultimo)difi=ultimo
                   ntimes=int(pasmin/dipa)        !numero de veces que la red fina divide a la obs 
c                  dipa=pasmin/ntimes

                   n1=nint((primero-dini)/dipa)  !numero de puntos anteriores
                   n2=nint((difi-ultimo)/dipa)   !nuero de puntos posteriores
                   dini=primero-n1*dipa
                   difi=ultimo+n2*dipa

                   npas(i)=int((difi-dini)/dipa)+1
                   if(10*( (difi-dini)/dipa+1 -int( (difi-dini)/dipa+1 ) ).gt..5) npas(i)=npas(i)+1
                   do j=1,npas(i)
                      nli=nli+1
                      dlamda(nli)=dini+(j-1)*dipa
                   end do
              end if
           end do  !fin del do en l

           if(icheck.eq.0)then
              call error(KSTOP,'leemalla2','There are lines in the'
     &        //         ' observed/stray light profiles which\n do not'
     &        //         ' appear in the wavelength grid.\n Check also the'
     &        //         ' wavelength grid for missing symbols')
           end if
        end do  !fin del do en i

c       k=1
c       epsiloni=dipa/2.
c       do i=1,nliobs
c           dd=dlamdaobs(i)
c           do while(k.lt.nli.and.dd.gt.dlamda(k)+epsiloni)
c              k=k+1
c           end do
c           do while(k.lt.nli.and.dd.lt.dlamda(k)-epsiloni)
c              k=k+1
c           end do
c           do while(k.lt.nli.and.dd.gt.dlamda(k)+epsiloni)
c              k=k+1
c           end do
c           nposi(i)=k
c       end do

        do i=1,nli
           nposi(i)=i
        end do

        if(contr.gt.-1.e5.and.nciclos.gt.0)then
           ntl=ntl+1
           ntlblends=ntlblends+1
           nlin(ntlblends)=0
           npas(ntl)=1
           nble(ntl)=1
           nli=nli+1
           nliobs=nliobs+1
           dlamda(nli)=0.
           nposi(nliobs)=nli
        end if

        return
        end
