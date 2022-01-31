c leemodi222
c iesc=0 :lee dos atmosferas modelo contando puntos en tau de la 1
c iesc=1 :escribe dos atmosferas modelo
c ntau : numero de puntos en tau
c tau  : log10(tau)
c ierror=0 O.K. ierror=1 el modelo 1 no existe, ierror=2 el modelo 2 no existe
c ierror=3 no existen los modelos 1 ni 2
C WARNING! leemodi222 reads angular variables written in degrees but the output is in radians
C WARNING! leemodi222 writes angular variables in degrees but the input is in radians

        subroutine leemodi222(iesc,model1,model2,atmos,ntau,ierror,
     &                        z1,pg1,ro1,z2,pg2,ro2)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        real*4 atmos(*),vmac1,vmac2
        real*4 tau(kt),t(kt),pe(kt),mic(kt),h(kt),vz(kt),gam(kt),phi(kt)
        real*4 pg1(*),z1(*),ro1(*)
        real*4 pg2(*),z2(*),ro2(*)
        character model1*(*),model2*(*),a*1
        character*20 msg1,msg2
        common/nohaycampo/nohaycampo1,nohaycampo2
        common/nciclos/nciclos   !para marquardt, leeuve3, leemodi2
        data ivez/0/

        epsilon=2.e-5
        ican=52
        dtor=1.7453293005625408e-2 ! degrees to radians
        if(iesc.eq.0)then   !reading

           ivez=ivez+1
c          leemos para contar las lineas
           if(ivez .eq. 1)then
              open(ican,file=model1,status='old',err=800)
              i=0
              read(ican,*,err=802)vmac1,fill1,peso1
              do while(i.lt.kt+2)
                 i=i+1
c                err=200 por si hay lineas en blanco
                 read(ican,*,end=200, err=200)a
              end do
200           ntau=i-1
c             por si al final del fichero hay un caracter blanco
              if(a .eq. ' ')then
                 ntau=i-2
              end if  
              close(ican)
           end if
           if(ntau. gt. kt)then
              write(msg1,*) ntau
              write(msg2,*) kt
              call error(KSTOP,'leemodi222','The number of depth points in'
     &        //         ' the model ('//trim(adjustl(msg1))//') is larger'
     &        //         ' than kt = '//trim(adjustl(msg2))//'\n'
     &        //         ' Decrease the number of depth points or change'
     &        //         ' the PARAMETER file')
           end if

c          leemos el modelo 1
           if(ierror.ne.1)then  
              open(ican,file=model1,status='old',err=800)
              read(ican,*,err=802)vmac1,fill1,peso1
              do i=1,ntau
                 read(ican,*)tau(i),t(i),pe(i),mic(i),
     &           h(i),vz(i),gam(i),phi(i),z1(i),pg1(i),ro1(i)
              end do
              close(ican)

              if(abs(fill1-1.0).gt.1.e-6 .and. ierror.eq.2)then 
                 call error(KWARN,'leemodi222','The filling factor is being'
     &           //         ' changed to 1.0')
                 fill1=1.
              endif
              if(fill1.lt.(-1.*epsilon).or.fill1.gt.1.+epsilon)then
                 call error(KSTOP,'leemodi222','For model 1, the filling'
     &           //         ' factor is outside the interval [0,1]')
              end if
              atmos(8*ntau+1)=vmac1
              atmos(8*ntau+2)=fill1
              if(peso1.lt.(-1.*epsilon).or.peso1.gt.100.+epsilon) then
                 call error(KSTOP,'leemodi222','The stray light of model 1'
     &           //         ' is outside the interval [0,100]')
              endif           
              atmos(16*ntau+5)=peso1  !% de luz difusa
              nohaycampo1=0
              do i=1,ntau
                 atmos(i)=tau(i)
                 atmos(i+ntau)=t(i)
                 atmos(i+2*ntau)=pe(i)
                 atmos(i+3*ntau)=mic(i)
                 atmos(i+4*ntau)=h(i)
                 if(abs(h(i)).gt.1.)nohaycampo1=1
                 atmos(i+5*ntau)=vz(i)
                 atmos(i+6*ntau)=gam(i)*dtor  !changing to radians
                 atmos(i+7*ntau)=phi(i)*dtor  !changing to radians
                 if(i .gt. 1)then
                    if((z1(i) - z1(i-1)) .lt. 1.)z1(i)=z1(i-1)+1.
                 end if 
                 if(pg1(i) .lt. 0)pg1(i)=abs(pg1(i))
                 if(ro1(i) .lt. 0)ro1(i)=abs(ro1(i))
              end do
              paso=atmos(2)-atmos(1)
              if(paso.ge.0)then
                 call error(KSTOP,'leemodi222','The models supplied are not'
     &           //         ' ordered with decreasing log(tau)')
              end if
           end if

c          leemos el modelo 2
           if(ierror.ne.2)then  
              if(ivez .eq. 1)then
                 open(ican,file=model2,status='old',err=801)
                 i=0
                 read(ican,*,err=802)vmac2,fill2,peso2
                 do while(i.lt.kt+2)
                    i=i+1
                    read(ican,*,end=221, err=221)a
                 end do
221              ntau2=i-1
                 if(a .eq. ' ')ntau=i-2
                 close(ican)
              end if
              if(ntau2.ne.ntau)then
                 call error(KSTOP,'leemodi222','The initial models are'
     &           //         ' discretized in different spatial grids')
              end if

              open(ican,file=model2,status='old',err=900)
              read(ican,*,err=802)vmac2,fill2,peso2
              do i=1,ntau
                 read(ican,*,end=222,err=802)tau(i),t(i),pe(i),mic(i),
     &           h(i),vz(i),gam(i),phi(i),z2(i),pg2(i),ro2(i)
222           end do  
              close(ican)

              if(abs(1-(fill1+fill2)).gt.epsilon)then
                 fill2=1.-fill1
                 write(msg1,*)fill2
                 call error(KWARN,'leemodi222','For model 2, the filling'
     &           //         ' factor is taken to be '//adjustl(msg1))
              end if
              atmos(16*ntau+3)=vmac2
              atmos(16*ntau+4)=fill2

              nohaycampo2=0
              do i=1,ntau
                 atmos(i+8*ntau+2)=tau(i)
                 atmos(i+9*ntau+2)=t(i)
                 atmos(i+10*ntau+2)=pe(i)
                 atmos(i+11*ntau+2)=mic(i)
                 atmos(i+12*ntau+2)=h(i)
                 if(abs(h(i)).gt.1.)nohaycampo2=1
                 atmos(i+13*ntau+2)=vz(i)
                 atmos(i+14*ntau+2)=gam(i)*dtor  !changing to radians
                 atmos(i+15*ntau+2)=phi(i)*dtor  !changing to radians
                 if(i .gt. 1)then
                    if((z2(i) - z2(i-1)) .lt. 1.)z2(i)=z2(i-1)+1.
                 end if 
                 if(pg2(i) .lt. 0)pg2(i)=abs(pg2(i))
                 if(ro2(i) .lt. 0)ro2(i)=abs(ro2(i))
              end do
              if(abs(1-(fill1+fill2)).gt.epsilon)then
                 fill2=1.-fill1
                 write(msg1,*)fill2
                 call error(KWARN,'leemodi222','For model 2, the filling'
     &           //         ' factor is taken to be '//adjustl(msg1))
              end if
              paso2=tau(2)-tau(1)
              if(paso2.ge.0)then
                 call error(KSTOP,'leemodi222','The models supplied are not'
     &           //         ' ordered with decreasing log(tau)')
              end if

              do i=2,ntau-1
                 paso1=atmos(i+1)-atmos(i)
                 paso2=atmos(8*ntau+i+3)-atmos(8*ntau+i+2)
                 if(abs(paso2-paso1).gt.epsilon)then
                    call error(KSTOP,'leemodi222','The same spatial grid has'
     &              //         ' to be used to discretize models 1 and 2')
                 end if
              end do

           end if !if(ierror.ne.2)

        else

c          escribimos los modelos
           peso=atmos(16*ntau+5)

           open(ican,file=model1,err=803)
           write(ican,*)atmos(8*ntau+1),atmos(8*ntau+2),peso
           do i=1,ntau
              if(atmos(i+ntau) .lt. 500)atmos(i+ntau)=500.            !T
              if(atmos(i+ntau) .gt. 9.999e4)atmos(i+ntau)=9.999e4     !T
              gamma_deg=atmos(i+6*ntau)/dtor  !inclination in degrees
              phi_deg=atmos(i+7*ntau)/dtor    !azimuth in degrees
              write(ican,100)(atmos(i+j*ntau),j=0,5),gamma_deg,phi_deg,z1(i),pg1(i),ro1(i)
           end do
           close(ican)
                
           if(ierror.eq.0)then
              open(ican,file=model2,err=803)
              write(ican,*)atmos(16*ntau+3),atmos(16*ntau+4),peso
              do i=1,ntau
                 if(atmos(i+9*ntau+2) .lt. 500)atmos(i+9*ntau+2)=500.        !T
                 if(atmos(i+9*ntau+2) .gt. 9.999e4)atmos(i+9*ntau+2)=9.999e4 !T
                   gamma_deg=atmos(i+14*ntau+2)/dtor  !inclination in degrees
                   phi_deg=atmos(i+15*ntau+2)/dtor    !azimuth in degrees                 
                 write(ican,100)(atmos(i+j*ntau+2),j=8,13),gamma_deg,phi_deg,z2(i),pg2(i),ro2(i)
              end do
              close(ican)
           end if

        end if

100     format(1x,f8.4,1x,f11.2,1x,1pe12.5,1x,8(e11.4,1x))
        return


800     call error(KSTOP,'leemodi222','The file containing model 1'
     &  //         ' does not exist')

801     call error(KSTOP,'leemodi222','The file containing model 2'
     &  //         ' does not exist')

802     call error(KSTOP,'leemodi222','Incorrect format in the file(s)'
     &  //         ' containing the model(s)')

803     call error(KSTOP,'leemodi222','Unable to open the file model 1 or 2')


900     call error(KWARN,'leemodi222','Is the file containing model 2'
     &  //         ' unexistent?\n File: '//model2)

        if(ierror.eq.0)then
           ierror=2
        else
           ierror=3
        end if 

        return  
        end
