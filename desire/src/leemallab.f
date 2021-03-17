c rutina leemallab
c lee el fichero de control de lineas y longitudes de onda:'mallaobs'
c ntau : numero de puntos en tau
c tau  : log10(tau)
c ntl  : numero total de lineas
c nlin : indice de cada linea
c npas : numero de puntos en cada linea
c dlamda:cada delta de l.d.o. en ma
c nble :numero de blends de cada linea

        subroutine leemallab(mallaobs,ntl,nli,nlin,npas,dlamda,nble)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        integer nlin(*),npas(*),nble(*),nd(kl),blanco,ntl,nli,ierror
        real*4 dlamda(*),dini,dipa,difi
        character mallaobs*(*)
        character*31 mensajito

        ican=57
        mensajito=' containing the wavelength grid'
        call cabecera(ican,mallaobs,mensajito,ifail)

        if(ifail.eq.1)then
           call error(KSTOP,'leemallab','The file containing the wavelength'
     &     //         ' grid does not exist\n File: '//mallaobs)
        end if

        jj=0
        k=0 
        nli=0
        numlin=0
        do i=1,1000  
           call mreadmalla(ican,ntonto,nd,dini,dipa,difi,ierror,blanco)
           if(ierror.eq.1)goto 8   !si he llegado al final del fichero
           if(blanco.ne.1)then
              numlin=numlin+1   

              if(numlin.gt.kl.or.jj+ntonto.gt.kl)then !comprobamos numero de lineas
                 call error(KSTOP,'leemallab','The number of lines in the'
     &           //         ' wavelength grid is larger than the\n current'
     &           //         ' limit. Decrease this number or change the'
     &           //         ' PARAMETER file')
              end if

              nble(numlin)=ntonto
              do j=1,nble(numlin)
                 jj=jj+1
                 nlin(jj)=nd(j)
              end do

              npas(numlin)=int((difi-dini)/dipa)+1
              if(10*( (difi-dini)/dipa+1 -int( (difi-dini)/dipa+1 ) ).gt..5) npas(numlin)=npas(numlin)+1

              if(nli+npas(numlin).gt.kld)then  !comprobamos numero longitudes onda
                 call error(KSTOP,'leemallab','The wavelength grid has more'
     &           //         ' points than the current limit kld.\n'
     &           //         ' Decrease the number of wavelengths or'
     &           //         ' change the PARAMETER file')
              end if

              do j=1,npas(numlin)
                 nli=nli+1
                 dlamda(nli)=dini+(j-1)*dipa
              end do
           end if   !el de blanco.ne.1
        end do

8       close(ican)
        ntl=numlin

        if(ntl.eq.0)then
           call error(KSTOP,'leemallab','Incorrect format in the file'
     &     //         ' containing the wavelength grid\n File: '//mallaobs)
        end if

        return
        end
