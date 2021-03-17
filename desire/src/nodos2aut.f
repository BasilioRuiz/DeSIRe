c rutina nodos2aut

        subroutine nodos2aut(ntau,iauto,ierror,ici,mnodos)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        integer ndiv(kt),mnodos(*),ntau,ierror,ici
        character*4 vari(18)
        character*100 msg

        data vari/'t1  ','p1  ','mic1','h1  ','vz1 ','gam1','fi1 ','mac1',
     &            't2  ','p2  ','mic2','h2  ','vz2 ','gam2','fi2 ','mac2',
     &            'f.f.','  % '/

c       Calculamos todos los divisores de ntau
        ndiv(1)=-1      !el primer divisor (perturbacion = otra variable)
        ndiv(2)=0       !segundo divisor   (no hay perturbacion)
        ndiv(3)=1       !el tercero        (perturbacion constante)
        ndiv(4)=2       !los 2 extremos    (perturbacion lineal)
        k=4             !ya llevamos 4
        do i=3,ntau     !es i-1 divisor de ntau-1?
           m=( (ntau-1)/(i-1) )*(i-1)
           if(m.eq.(ntau-1))then
              k=k+1
              ndiv(k)=i
           end if
        end do
        num=k

c       ----------------------
c       Comprobaciones varias:
c       ----------------------
        if(ierror.eq.2)then    !no podemos tomar el parametro del otro modelo
           do j=1,8
              if(mnodos(j).eq.-1)then
                 mnodos(j)=0
                 call error(KWARN,'nodos2aut','Nodes for '//vari(j)
     &           //         ' changed to 0 (there is only one component)')
              endif
           enddo
        endif

        if(mnodos(8).gt.1)then
           mnodos(8)=1
           call error(KWARN,'nodos2aut','The number of nodes for the'
     &     //         ' macroturbulence (model 1) should be\n'
     &     //         ' smaller or equal than 1. One node is being adopted')
        end if
        if(mnodos(16).gt.1)then
           mnodos(16)=1
           call error(KWARN,'nodos2aut','The number of nodes for the'
     &     //         ' macroturbulence (model 2) should be\n'
     &     //         ' smaller or equal than 1. One node is being adopted')
        end if
        if(mnodos(17).gt.1.or.mnodos(17).lt.0)then
           mnodos(17)=1
           call error(KWARN,'nodos2aut','The number of nodes for'
     &     //         ' the filling factor has to be either\n'
     &     //         ' 0 or 1. One node is being adopted')
        end if
        if(mnodos(18).gt.1.or.mnodos(18).lt.0)then
           mnodos(18)=1
           call error(KWARN,'nodos2aut','The number of nodes for the'
     &     //         ' stray light factor has to be either\n'
     &     //         ' 0 or 1. One node is being adopted')
        end if
c       ----------------------

        if(iauto.eq.0)then
c          Buscamos el divisor menor o igual mas cercano a cada mnodos
           do j=1,18
              mold=mnodos(j)
              if(mnodos(j).lt.-1)mnodos(j)=-1
              k=1
              do while( ndiv(k).le.mnodos(j) .and. k.le.num)
                 k=k+1
              end do
              mnodos(j)=ndiv(k-1)
              if(mnodos(j).ne.mold)then
                 write(msg,'(a,a,a,i2)') 'Nodes for ',vari(j),' changed to ',
     &                 mnodos(j)
                 call error(KWARN,'nodos2aut',msg)
              end if
           end do

        else  !seleccion automatica de nodos

           if(mnodos(1).ge.1)mnodos(1)=ndiv(ici+3)
           if(mnodos(2).ge.1)mnodos(2)=1  !solo permito una cte
           if(mnodos(3).ge.1)mnodos(3)=1

           do i=4,7
              if(mnodos(i).ge.1)mnodos(i)=ndiv(ici+1)
           enddo

           do j=1,8
              write(msg,'(a,a,a,i2)') 'Nodes for ',vari(j),': ',mnodos(j)
              call error(KLINE,'',msg)
           enddo

           if(ierror.eq.0)then
              if(mnodos(9).ge.1)mnodos(9)=ndiv(ici+2)
              if(mnodos(10).ge.1)mnodos(10)=1  !solo permito una cte
              if(mnodos(11).ge.1)mnodos(11)=1
              do i=12,15
                 if(mnodos(i).ge.1)mnodos(i)=ndiv(ici+1)
              enddo
              do j=9,17
                 write(msg,'(a,a,a,i2)') 'Nodes for ',vari(j),': ',mnodos(j)
                 call error(KLINE,'',msg)
              enddo
           endif
           write(msg,'(a,a,a,i2)') 'Nodes for ',vari(18),': ',mnodos(18)
           call error(KPARA,'',msg)
        endif

        return
        end
