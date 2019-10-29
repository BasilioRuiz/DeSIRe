        subroutine interpolate1(tauoriginal,atmos,ntau,imodel,z,pg,ro)
        
        implicit real*4 (a-h,o-z)
	include 'PARAMETER'   !por kt,kn,kl,kld    
	
	real*4 atmos(*),z(*),pg(*),ro(*),tauoriginal(*)
	real*4 taunew(ntau),ti,tau1,xa(3),ya(3)
	real*4 tau(ntau),atmosnew(16*kt+2)
	real*4 znew(kt),pgnew(kt),ronew(kt)
	integer ntau,ngrado,n2,num,n,imodel,irow,nn
	
	num=ntau !number of original depth points
	ngrado=2
	n2=int(ngrado/2)
	n=ntau   !number of final depth points
	nn=(8*ntau+2)*(imodel-1) !=0 for imodel=1, =8*ntau+2 for imodel=2

	do i=1,num
	   tau(i)=atmos(nn+i)	   
	end do
	
c	tau1=tau(1)
c	step=(tau1-tau(n))/float(n-1)
c	step=nint(step*1000.)/1000.
c	tau1=tau(n)+float(n-1)*step
	
	do i=1,n
c	   taunewi=tau1-step*(i-1)	
	   taunewi=tauoriginal(i)
           atmosnew(i)=taunewi	
           
	   call locate(tau,num,taunewi,j)
	   n3=j-n2-1
           if(n3.lt.0)n3=0
c           print*,'j n3 n2 tau  ',j,n3,n2,tau(n3+1),tau(n3+2),tau(n3+3)
           if(j .le. 0)ngrado=1
           if(n3 .gt. 0)ngrado=2
           if(j .ge. num-1)ngrado=1
c           if(j .ge. num)ngrado=0

           if(n3+ngrado+1.gt.num)n3=num-ngrado-1

	   do k=1,ngrado+1
	      xa(k)=tau(n3+k)
           end do
           irow=1                          !T
           do k=1,ngrado+1
	      ya(k)=atmos(nn+irow*ntau+n3+k)   
           end do	
c           if(ngrado .eq. 1)ya(2)=ya(1)+(ya(2)-ya(1))*(ntau-i+1)/float(ntau)
	   call polint(xa,ya,ngrado+1,taunewi,ti,error) 
	   atmosnew(nn+irow*ntau+i)=ti
            
           irow=2                           !log(pe)
	   do k=1,ngrado+1
	      ya(k)=alog(abs(atmos(nn+irow*ntau+n3+k)))  
c	      print*,'interpoate1 ',i,k,xa(k),ya(k)
	   end do
c	   if(ngrado .eq. 1)ya(2)=ya(1)+(ya(2)-ya(1))*(ntau-i+1)/float(ntau)

	   call polint(xa,ya,ngrado+1,taunewi,ti,error)
	   atmosnew(nn+irow*ntau+i)=exp(ti)
c	   print*,'interpoate1 ',i,taunewi,ti,exp(ti)	   
	   
	   do irow=3,7                     !mic,h,vz,g,phi
	      do k=1,ngrado+1
	         ya(k)=atmos(nn+irow*ntau+n3+k)   
              end do
c              if(ngrado .eq. 1)ya(2)=ya(1)+(ya(2)-ya(1))*(ntau-i+1)/float(ntau)

	      call polint(xa,ya,ngrado+1,taunewi,ti,error) 
	      atmosnew(nn+irow*ntau+i)=ti
           end do

	   do k=1,ngrado+1
	      ya(k)=z(n3+k)                         !z
           end do	
c           if(ngrado .eq. 1)ya(2)=ya(1)+(ya(2)-ya(1))*(ntau-i+1)/float(ntau)

	   call polint(xa,ya,ngrado+1,taunewi,ti,error) 
	   znew(i)=ti
	   
	   do k=1,ngrado+1
	      ya(k)=alog(abs(pg(n3+k)))             !pg
           end do	
c           if(ngrado .eq. 1)ya(2)=ya(1)+(ya(2)-ya(1))*(ntau-i+1)/float(ntau)
	   call polint(xa,ya,ngrado+1,taunewi,ti,error) 
	   pgnew(i)=exp(ti)
	   
           do k=1,ngrado+1
	      ya(k)=alog(abs(ro(n3+k)))             !rho
           end do	
c           if(ngrado .eq. 1)ya(2)=ya(1)+(ya(2)-ya(1))*(ntau-i+1)/float(ntau)
	   call polint(xa,ya,ngrado+1,taunewi,ti,error) 
	   ronew(i)=exp(ti)
	end do   
	do i=1,n
	   do kk=1,8
	      atmos(nn+i+(kk-1)*ntau)=atmosnew(nn+i+(kk-1)*ntau)
	   end do
	   z(i)=znew(i)
	   pg(i)=pgnew(i)
	   ro(i)=ronew(i)
	end do
	do i=2,n
	   deltaz=z(i)-z(i-1)
	   if(deltaz .lt. 1.)deltaz=1
	   z(i)=z(i-1)+deltaz
	enddo
	

	return
	end
	