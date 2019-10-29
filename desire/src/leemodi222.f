c leemodi222
c iesc=0 :lee dos atmosferas modelo contando puntos en tau de la 1
c iesc=1 :escribe dos atmosferas modelo
c ntau : numero de puntos en tau
c tau  : log10(tau)
c ierror=0 O.K. ierror=1 el modelo 1 no existe, ierror=2 el modelo 2 no existe
c ierror=3 no existen los modelos 1 ni 2

	subroutine leemodi222(iesc,model1,model2,atmos,ntau,ierror,
     &                        z1,pg1,ro1,z2,pg2,ro2)
        implicit real*4 (a-h,o-z)

	include 'PARAMETER'  !por kt
	real*4 atmos(*),a(8),vmac1,vmac2
	real*4 tau(kt),t(kt),pe(kt),mic(kt),h(kt),vz(kt),gam(kt),phi(kt)
        real*4 pg1(*),z1(*),ro1(*)
        real*4 pg2(*),z2(*),ro2(*)
        character model1*(*),model2*(*)
	character*100 control
	character*80 men1
	common/canal/icanal
	common/nombrecontrol/control
	common/nohaycampo/nohaycampo1,nohaycampo2
        common/nciclos/nciclos   !para marquardt, leeuve3, leemodi2

	epsilon=2.e-5	
	ican=52
	if(iesc.eq.0)then

c leemos el modelo 1
           if(ierror.ne.1)then  
	     open(ican,file=model1,status='old',err=991)
	     ntau=0
	     read(ican,*,err=999)vmac1,fill1,peso1
	     do while(ntau.lt.kt)
	        ntau=ntau+1
	        read(ican,*,end=221,err=999)tau(ntau),t(ntau),pe(ntau),mic(ntau),
     &	        h(ntau),vz(ntau),gam(ntau),phi(ntau),z1(ntau),pg1(ntau),ro1(ntau)
	     end do
221	     ntau=ntau-1
	     close(ican)
	     if(abs(fill1-1.0).ne.0 .and. ierror.eq.2)then 
	        men1='     WARNING: The FILLING FACTOR is being changed to 1.0'
	        print*,men1
	        fill1=1.
	     endif
	     if(fill1.lt.(-1.*epsilon).or.fill1.gt.1.+epsilon)then
	        men1='STOP: For model 1, the filling factor is outside the interval [0,1] '
	        call mensaje(1,men1,men1,men1)
             end if
	     atmos(8*ntau+1)=vmac1
	     atmos(8*ntau+2)=fill1
	     if(peso1.lt.(-1.*epsilon).or.peso1.gt.100.+epsilon) then
	      men1='STOP: In model 1, the stray light factor is outside the interval [0,100]'
	      call mensaje(1,men1,men1,men1)
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
                atmos(i+6*ntau)=gam(i)
                atmos(i+7*ntau)=phi(i)
	        if(i .gt. 1)then
	          if((z1(i) - z1(i-1)) .lt. 1.)z1(i)=z1(i-1)+1.
	        end if 
	        if(pg1(i) .lt. 0)pg1(i)=abs(pg1(i))
	        if(ro1(i) .lt. 0)ro1(i)=abs(ro1(i))
	     end do
	     paso=atmos(2)-atmos(1)
	     if(paso.ge.0)then
	        men1='STOP: The models supplied are NOT ordered with DECREASING log(tau).'    
	        call mensaje(1,men1,men1,men1)
	     end if

             do i=2,ntau-1
	        paso1=atmos(i+1)-atmos(i)
                if(abs(paso1-paso).gt.epsilon)then
	           if(nciclos.gt.0)then
		      men1='STOP: Model 1 is NOT equally spaced. For inversion, an evenly spaced grid is required'
	              call mensaje(1,men1,men1,men1)
	           endif
                end if
	     end do
	       
           end if

c leemos el modelo 2
           if(ierror.ne.2)then  
	     open(ican,file=model2,status='old',err=988)
	     ntau2=0
	     read(ican,*,err=999)vmac2,fill2,peso2
	     do while(ntau2.lt.kt)
	       ntau2=ntau2+1
	       read(ican,*,end=222,err=999)tau(ntau2),t(ntau2),pe(ntau2),mic(ntau2),
     &	       h(ntau2),vz(ntau2),gam(ntau2),phi(ntau2),z2(ntau2),pg2(ntau2),ro2(ntau2)
             end do
222	     ntau2=ntau2-1   
             close(ican)
             if(ntau2.ne.ntau)go to 993
	     if(abs(1-(fill1+fill2)).gt.epsilon)then
                fill2=1.-fill1
                print*,'     WARNING: For model 2, the filling factor is taken to be ',fill2
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
                atmos(i+14*ntau+2)=gam(i)
                atmos(i+15*ntau+2)=phi(i)
	        if(i .gt. 1)then
	          if((z2(i) - z2(i-1)) .lt. 1.)z2(i)=z2(i-1)+1.
	        end if 
	        if(pg2(i) .lt. 0)pg2(i)=abs(pg2(i))
	        if(ro2(i) .lt. 0)ro2(i)=abs(ro2(i))
	     end do
	     if(abs(1-(fill1+fill2)).gt.epsilon)then
                fill2=1.-fill1
                print*,'     WARNING: For model 2, the filling factor is taken to be ',fill2
	     end if
	     paso2=tau(2)-tau(1)
	     if(paso2.ge.0)then
	        men1='STOP: The models supplied are NOT ordered with DECREASING log(tau).'    
	        call mensaje(1,men1,men1,men1)
	     end if

             do i=2,ntau-1
	        paso1=atmos(i+1)-atmos(i)
                paso2=atmos(8*ntau+i+3)-atmos(8*ntau+i+2)
                if(abs(paso2-paso1).gt.epsilon)then
	           men1='STOP: The same spatial grid has to be used to discretize models 1 and 2.'       
	           call mensaje(1,men1,men1,men1)
                end if
	     end do
	     
             return   

993	     if(ntau2.ne.ntau)then
  	       men1='STOP: The initial models are discretized in DIFFERENT SPATIAL GRIDS.'             
	       call mensaje(1,men1,men1,men1)
             end if
           end if
         
	else
	
c escribimos los modelos
	   peso=atmos(16*ntau+5)

	   open(ican,file=model1,err=990)
           goto 800 
990        print*,' '
	   print*,'WARNING: The file containing model 1 does NOT exist'
           ierror=1
	   goto 989

800	   write(ican,*)atmos(8*ntau+1),atmos(8*ntau+2),peso
	   do i=1,ntau
	      if(atmos(i+ntau) .lt. 500)atmos(i+ntau)=500.            !T
              if(atmos(i+ntau) .gt. 9.999e4)atmos(i+ntau)=9.999e4     !T
	      write(ican,100)(atmos(i+j*ntau),j=0,7),z1(i),pg1(i),ro1(i)
	   end do
	   close(ican)
          	
989	   if(ierror.eq.0)then
	      open(ican,file=model2,err=988)
	      write(ican,*)atmos(16*ntau+3),atmos(16*ntau+4),peso
	      do i=1,ntau
	      	 if(atmos(i+9*ntau+2) .lt. 500)atmos(i+9*ntau+2)=500.        !T
                 if(atmos(i+9*ntau+2) .gt. 9.999e4)atmos(i+9*ntau+2)=9.999e4 !T
	         write(ican,100)(atmos(i+j*ntau+2),j=8,15),z2(i),pg2(i),ro2(i)
	      end do
	      close(ican)
	   end if
        end if 

100     format(1x,f8.4,1x,f11.2,1x,1pe12.5,1x,8(e11.4,1x))
	return

991     men1='STOP: The file containing model 1 does NOT exist'
	call mensaje(1,men1,men1,men1)


988 	print*,'WARNING: Is the file containing model 2, ',model2,', unexistent?? '

        if(ierror.eq.0)then
           ierror=2
        else
           ierror=3
        end if 
	return	


999     men1='STOP: (leemodi222) Incorrect format in the file(s) containing the model(s).'
	call mensaje(1,men1,men1,men1)

	end
