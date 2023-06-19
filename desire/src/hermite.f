	subroutine hermite(s,ds,etall,n,svec,nn,bt,tk,pk,hk,
     &           vk,gk,fk,mk,rt,rp,rh,rv,rg,rf,rm,mnodos,icontorno)

        implicit real*4 (a-h,o-z)
	include 'PARAMETER'
        parameter (kt16=16*kt, kfac=5, ktnew=kfac*(kt-1)+1)
        parameter (dl10=2.3025851)

	integer n,nn
	real*4 svec(4)
        integer*4 mnodos(*)

c internas
        real*4 sol(4,kt),tau(kt),taue(kt)  
        real*4 deltae(kt),deltai(kt),delt2i(kt)
        real*4 ds(kt), etall(4,4,kt)
        real*4 s(kt)
        real*4 OO(4,4,kt)
        real*4 tk(kt16),pk(kt16),hk(kt16),vk(kt16),gk(kt16),fk(kt16),mk(kt16)
        real*4 bt(kt)
        real*4 rt(4,kt),rp(4,kt),rh(4,kt),rv(4,kt),rg(4,kt),rf(4,kt),rm(4,kt)

c interpolacion
        real*4 xa(11),ya(11)
        real*4 taunew(ktnew),tauenewxx(ktnew)
        real*4 etallnew(4,4,ktnew),deltainew(ktnew),delt2inew(ktnew)
        real*4 etalnew(4,4,ktnew),OOnew(4,4,ktnew)
        real*4 solnew(4,ktnew),snew(ktnew),dsnew(ktnew)
        integer iold(kt)
		
        common/segunda/tau,taue,deltae,deltai,delt2i
                        
        if(abs(etall(1,1,icontorno)).le.50.)then
           interpola=0  
           nnew=n
           do j=1,nnew
              taunew(j)=tau(j)
              xx=taue(j)*dl10
              deltainew(j)=deltai(j)
              delt2inew(j)=delt2i(j)
              do jj=1,4
                 do kk=1,4             
                    etallnew(jj,kk,j)=etall(jj,kk,j)
                    etalnew(jj,kk,j)=etall(jj,kk,j)*xx
                 enddo
              end do 
              snew(j)=s(j)
              dsnew(j)=ds(j)
           enddo   
	       icontornonew=icontorno
        else                           !intepolating if(abs(etall(1,1,icontorno)).gt.10.)then
           icontornonew=kfac*(icontorno-1)+1
           interpola=1

           do j=1,n
              iold(j)=(j-1)*kfac+1
           enddo

            nnew=kfac*(n-1)+1
            steptaunew=(tau(1)-tau(2))/(kfac*1.)
            do j=1,nnew
               taunew(j)=tau(1)-(j-1)*steptaunew
               tauenewxx(j)=10.**taunew(j)*dl10
            end do  
            do j=2,nnew
               deltainew(j)=(taunew(j)-taunew(j-1))/2.
               delt2inew(j)=deltainew(j)*deltainew(j)/3.
            end do
            deltainew(1)=deltainew(2)
            delt2inew(1)=delt2inew(2)
            ngrado=2
            n2=int(ngrado/2)
            do i=1,nnew
               taunewi=taunew(i)
               tauenewxxi=tauenewxx(i)
               jx=(i-1)/kfac+1
               n3=jx-n2-1
               if(n3.lt.0)n3=0
               if(n3+ngrado+1.gt.nnew)n3=nnew-ngrado-1
               do k=1,ngrado+1
                  xa(k)=tau(n3+k)
                  ya(k)=s(n3+k)	               
               enddo
	           call polint(xa,ya,ngrado+1,taunewi,snew(i),dy)
	       
	           do k=1,ngrado+1
	              ya(k)=ds(n3+k)	               
	           enddo
               call polint(xa,ya,ngrado+1,taunewi,dsnew(i),dy)
	       
 	           do jj=1,4
	              do kk=1,4             
	                 do k=1,ngrado+1
	                    ya(k)=etall(jj,kk,n3+k)	               
	                 enddo
	                 call polint(xa,ya,ngrado+1,taunewi,yyy,dy)
	                 etallnew(jj,kk,i)=yyy
	                 etalnew(jj,kk,i)=yyy*tauenewxxi
	                 do k=1,ngrado+1
	                    ya(k)=oo(jj,kk,n3+k)	               
	                 enddo
	                 call polint(xa,ya,ngrado+1,taunewi,yyy,dy)
	                 oonew(jj,kk,i)=yyy
	              enddo
	           end do  
             enddo   
         endif

         call inicializanew(ktnew,icontornonew,etalnew,snew,dsnew,oonew,solnew,ierror)

         call numden(icontornonew,nnew,ktnew,taunew,deltainew,delt2inew,etalnew,snew,dsnew,OOnew,solnew,svec)
      
         if( interpola .eq. 0)then
            do i=1,n
               do jj=1,4
                  sol(jj,i)=solnew(jj,i)
	              do kk=1,4 
                     OO(jj,kk,i)=oonew(jj,kk,i)	         
	              enddo
	          enddo   
           enddo
        else
           do i=1,n
              ii=iold(i)
              do jj=1,4           
                 sol(jj,i)=solnew(jj,ii)
	             do kk=1,4 
                    OO(jj,kk,i)=oonew(jj,kk,ii)	         
                 enddo
              enddo   
           enddo
        end if
      
        if(ierror.eq.1)then
           call error(KSTOP,'hermite','Error in the integration of the RTE')
        endif
	
        do i=1,n
           sol(1,i) = sol(1,i)-s(i)
        enddo

        if(mnodos(1).ne.0)call delta1(n,tk,bt,sol,OO,etall,rt)
        if(mnodos(1).ne.0.or.mnodos(2).ne.0.or.mnodos(4).ne.0.or.
     &     mnodos(6).ne.0)call delta2(n,pk,sol,OO,rp)
        if(mnodos(3).ne.0)call delta2(n,mk,sol,OO,rm)
        if(mnodos(4).ne.0)call delta2(n,hk,sol,OO,rh)
        if(mnodos(5).ne.0)call delta2(n,vk,sol,OO,rv)
        if(mnodos(6).ne.0)call delta2(n,gk,sol,OO,rg)
        if(mnodos(7).ne.0)call delta2(n,fk,sol,OO,rf)

        return
        end

c ___________________________________________________________________
       subroutine inicializanew(ktt,icontorno,etal,s,ds,oo,sol,ierror)

       implicit real*4 (a-h,o-z)

       include 'PARAMETER' !por kt
       parameter (kt16=16*kt)
       real*4 etal(4,4,ktt),inveta(4,4),s(*),ds(*),oo(4,4,ktt),sol(4,ktt)

       do k=1,icontorno
	  do j=1,4
	     do i=1,4
	        inveta(i,j)=etal(i,j,k)
	     end do
	  end do

	  call matinx2(inveta,ierror)	!calcula la matriz inversa 4*4

	  do i=1,4
	     sol(i,k)=inveta(i,1)*ds(k)
	  end do
	  sol(1,k)=sol(1,k)+s(k)
       end do

	do l=2,icontorno+1   !aprox. de difusion para el op. evolucion
	   do j=1,4
              do i=1,4
                OO(i,j,l)=0.
	      enddo
            enddo       

	    OO(1,1,l-1)=0.   !(sol(1,l)-smed)/(sol(1,i-1)-smed)
	enddo

       return
       end

c-------------------------------------------------------------------
       subroutine numden(icontorno,n,kt,tau,deltai,delt2i,etal,s,ds,OO,sol,svec)
       implicit real*4 (a-h,o-z)
              
       integer icontorno,n,i,ii,jj,k
       real*4 etal(4,4,kt),etal2(4,4,kt),detal(4,4,kt)
       real*4 s(kt), ds(kt),den(4), num(4,4), den2(4,4)
       real*4 sol(4,kt),OO(4,4,kt),tau(kt)
       real*4 deltai(kt),delt2i(kt),svec(4)
       
       call deriva4cuad(etal,detal,n,tau)

       call etal2sub(icontorno,n,kt,etal,etal2) 
       
       do i=icontorno+1,n
           delta=deltai(i)
           delt2=delt2i(i)

	   do jj=1,4
              e0=etal(jj,1,i-1)
	      e1=etal(jj,1,i)
              s00=s(i-1)
              s01=s(i)
	      den(jj) =-( e0*s00 + e1*s01 )*
     &	                 delta + sol(jj,i-1)
	      den(jj) =-(detal(jj,1,i-1)*s00 - detal(jj,1,i)*s01  +
     &                  e0*ds(i-1) -  e1*ds(i) +
     &                  etal2(jj,1,i-1)*s00 - etal2(jj,1,i)*s01 )*
     &                  delt2 + den(jj)

	      do ii=1,4
	         den2(ii,jj)=etal(ii,jj,i-1)*delta+
     & delt2*(detal(ii,jj,i-1)+etal2(ii,jj,i-1))
	         den(jj)=den(jj) + (etal(jj,ii,i-1)*delta+
     & delt2*(detal(jj,ii,i-1)+etal2(jj,ii,i-1)))*sol(ii,i-1)
	         num(ii,jj) = - etal(ii,jj,i)*delta + 
     & delt2*(detal(ii,jj,i)+etal2(ii,jj,i))
	      enddo
	      num(jj,jj) = 1.0 + num(jj,jj)
	      den2(jj,jj) = 1.0 + den2(jj,jj)
	   enddo
	   call matinx(num)

	   do jj=1,4
              suma=0.
	      do ii=1,4
	         suma = suma + num(jj,ii)*den(ii)
	         sumaoo= 0.e0
	         do k=1,4
	            sumaoo = sumaoo + num(ii,k)*den2(k,jj)
	         enddo
	         OO(ii,jj,i-1) = sumaoo
	      enddo
	      sol(jj,i) = suma
	   enddo
	enddo
	
	do ii=1,4
	   svec(ii) = sol(ii,n)
	enddo
	
	do jj=1,4
	   do ii=1,4
	      OO(ii,jj,n) = 0.e0
	   enddo
	   OO(jj,jj,n) = 1.e0
	enddo
	
	do i=n-2,1,-1
	   do jj=1,4
	      do kk=1,4
                 suma=0.e0
	         do k=1,4
	            suma=suma+OO(kk,k,i+1)*OO(k,jj,i)
	         enddo
	         den2(kk,jj) = suma
              enddo
	   enddo
	   do jj=1,4
	      do kk=1,4
	         OO(kk,jj,i)=den2(kk,jj)
	      enddo
	   enddo
	enddo
       
       return
       end
c----------------------------------------------------------------       
      subroutine etal2sub(icontorno,n,kt,etal,etal2) 
      integer icontorno,n,i,ii,jj,k
      real*4 etal(4,4,kt),etal2(4,4,kt),ssuma
      
	do i=icontorno,n
	   do jj=1,4
	      do ii=1,4
                ssuma=0.
	         do k=1,4
	            ssuma = ssuma +  etal(ii,k,i)*etal(k,jj,i)
	         enddo
	         etal2(ii,jj,i) = ssuma
	      enddo
	   enddo
	enddo
	
       return
       end
c----------------------------------------------------------------       
       
       
	subroutine hermite_no_interpol(s,ds,etall,n,svec,nn,bt,tk,pk,hk,
     &           vk,gk,fk,mk,rt,rp,rh,rv,rg,rf,rm,mnodos,icontorno)

        implicit real*4 (a-h,o-z)
	include 'PARAMETER' !por kt
        parameter (kt16=16*kt)
        parameter (dl10=2.3025851)

	integer n,nn
	real*4 svec(4)
        integer*4 mnodos(*)

c internas
        real*4 sol(4,kt),tau(kt),taue(kt),xx
        real*4 deltae(kt),deltai(kt),delt2i(kt)
	real*4 ds(kt), detal(4,4,kt), etal(4,4,kt), etall(4,4,kt)
	real*4 etal2(4,4,kt)
	real*4 s(kt), den(4), num(4,4), den2(4,4)
	real*4 OO(4,4,kt)
        real*4 tk(kt16),pk(kt16),hk(kt16),vk(kt16),gk(kt16),fk(kt16),mk(kt16)
	real*4 bt(kt)
	real*4 rt(4,kt),rp(4,kt),rh(4,kt),rv(4,kt),rg(4,kt),rf(4,kt),rm(4,kt)

        common/segunda/tau,taue,deltae,deltai,delt2i


c se calcula la derivada de la matriz de absorcion

	do i=1,n
	   xx = taue(i)*dl10
	   do jj=1,4
	      do kk=1,4
	         etal(jj,kk,i) = etall(jj,kk,i)*xx
	      enddo
	   enddo
  	end do 	
	call deriva4cuad(etal,detal,n,tau)


c se integra
        call inicializa(icontorno,etal,s,ds,oo,sol,ierror)
	if(ierror.eq.1)return
	  
	do i=icontorno,n
	   do jj=1,4
	      do ii=1,4
                ssum=0.
	         do k=1,4
	            ssum = ssum +  etal(ii,k,i)*etal(k,jj,i)
	         enddo
	         etal2(ii,jj,i) = ssum
	      enddo
	   enddo
	enddo

 	do i=icontorno+1,n
           delta=deltai(i)
           delt2=delt2i(i)
	   do jj=1,4
              e0=etal(jj,1,i-1)
	      e1=etal(jj,1,i)
              s00=s(i-1)
              s01=s(i)
	      den(jj) =-( e0*s00 + e1*s01 )*
     &	                 delta + sol(jj,i-1)
	      den(jj) =-(detal(jj,1,i-1)*s00 - detal(jj,1,i)*s01  +
     &                  e0*ds(i-1) -  e1*ds(i) +
     &                  etal2(jj,1,i-1)*s00 - etal2(jj,1,i)*s01 )*
     &                  delt2 + den(jj)

	      do ii=1,4
	         den2(ii,jj)=etal(ii,jj,i-1)*delta+
     & delt2*(detal(ii,jj,i-1)+etal2(ii,jj,i-1))
	         den(jj)=den(jj) + (etal(jj,ii,i-1)*delta+
     & delt2*(detal(jj,ii,i-1)+etal2(jj,ii,i-1)))*sol(ii,i-1)
	         num(ii,jj) = - etal(ii,jj,i)*delta + 
     & delt2*(detal(ii,jj,i)+etal2(ii,jj,i))
	      enddo
	      num(jj,jj) = 1.0 + num(jj,jj)
	      den2(jj,jj) = 1.0 + den2(jj,jj)
	   enddo
     
	   call matinx(num)

	   do jj=1,4
              sum=0.
	      do ii=1,4
	         sum = sum + num(jj,ii)*den(ii)
	         sumoo= 0.e0
	         do k=1,4
	            sumoo = sumoo + num(ii,k)*den2(k,jj)
	         enddo
	         OO(ii,jj,i-1) = sumoo
	      enddo
	      sol(jj,i) = sum
	   enddo
	enddo

	do ii=1,4
	   svec(ii) = sol(ii,n)
	enddo

	do jj=1,4
	   do ii=1,4
	      OO(ii,jj,n) = 0.e0
	   enddo
	   OO(jj,jj,n) = 1.e0
	enddo
	
	do i=n-2,icontorno,-1
	   do jj=1,4
	      do kk=1,4
                 sum=0.e0
	         do k=1,4
	            sum=sum+OO(kk,k,i+1)*OO(k,jj,i)
	         enddo
	         den2(kk,jj) = sum
              enddo
	   enddo
	   do jj=1,4
	      do kk=1,4
	         OO(kk,jj,i)=den2(kk,jj)
	      enddo
	   enddo
	enddo
	do i=1,icontorno-1
           do jj=1,4
	      do kk=1,4
	         OO(kk,jj,i)=0.
	      enddo
	   enddo
	end do

	do i=1,n
	   sol(1,i) = sol(1,i)-s(i)
	enddo

	if(mnodos(1).ne.0)call DELTA1(n,tk,bt,sol,OO,etall,rt)
	if(mnodos(1).ne.0.or.mnodos(2).ne.0.or.mnodos(4).ne.0.or.
     &     mnodos(6).ne.0)call DELTA2(n,pk,sol,OO,rp)
	if(mnodos(3).ne.0)call DELTA2(n,mk,sol,OO,rm)
	if(mnodos(4).ne.0)call DELTA2(n,hk,sol,OO,rh)
	if(mnodos(5).ne.0)call DELTA2(n,vk,sol,OO,rv)
	if(mnodos(6).ne.0)call DELTA2(n,gk,sol,OO,rg)
	if(mnodos(7).ne.0)call DELTA2(n,fk,sol,OO,rf)

	return
	end
