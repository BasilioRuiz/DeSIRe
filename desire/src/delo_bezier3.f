	subroutine delo_bezier3(s,ds,etall,n,svec,nn,bt,tk,pk,hk,
     &           vk,gk,fk,mk,rt,rp,rh,rv,rg,rf,rm,mnodos)

        implicit real*4 (a-h,o-z)
	include 'PARAMETER' !por kt
        parameter (kt16=16*kt)
        parameter (dl10=2.3025851)

	integer n,nn,i,j,k,kr,ierror,jj
	real*4 svec(4)
        integer*4 mnodos(*)

c internas
        real*4 sol(4,kt),tau(kt),taue(kt),s0(kt)
        real*4 deltae(kt),deltai(kt),delt2i(kt)
	real*4 ds(kt), etall(4,4,kt)  
	real*4 s(kt)
	real*4 OO(4,4,kt),suma1ij
        real*4 tk(kt16),pk(kt16),hk(kt16),vk(kt16),gk(kt16),fk(kt16),mk(kt16)
	real*4 bt(kt)
	real*4 rt(4,kt),rp(4,kt),rh(4,kt),rv(4,kt),rg(4,kt),rf(4,kt),rm(4,kt)
	
c para la integracion
        real*8 taunui(kt)
        real*8 xx(kt),etadb(kt),Ck(kt),etaii
        real*8 Sp(4,kt),solp(4,kt),deltabar(kt),solpk1(4)
        real*8 K4(4,4,kt),K4der(4,4,kt)
        real*8 K4_k(4,4),K4_k1(4,4),K42_k(4,4),K42_k1(4,4)
        real*8 e_k(4),e_k1(4),ebar_k(4,4),ebar_k1(4,4),fkvec(4),suma1,suma2  
        real*8 Spder(4,kt),Sgotvec(4)
        real*8 dk,dk3,dk_3,Exk,alphak,betak,gammak,epsik
        real*8 Ek(kt),Fk1(kt)
        real*8 Ak(4,4),Bk(4),Xk(4),Bmat(4,4),OOk(4,4),OOkr(4,4,kt)!,C(4),OOkr8(4,4),svecr8(4)
        real*8 OO1k_1(4,4),OO1k(4,4),OOfromkto1(4,4,kt)
        real*8 Bmatinv(4,4),Akinv(4,4)

        common/segunda/tau,taue,deltae,deltai,delt2i

c We need to reverse the S, taue & etall arrays and integrate from n to 1 
	do k=1,n
	   kr=n-k+1
           xx(k)=dble(taue(kr))                !x in paper
           etaii=dble(etall(1,1,kr))
           etadb(k)=etaii
           do j=1,4
              do i=1,4
                 K4_k(i,j)=dble(etall(i,j,kr))/etaii
                 K4(i,j,k)=K4_k(i,j)
              end do
              K4(j,j,k)=K4(j,j,k)-1.d0
              Sp(j,k)=dble(S(kr))*K4_k(j,1) !changing S_SIR to S_Delo: dI/dtau=K(I-S_SIR)  dI/dtau=KI-S_Delo
           end do                           !S_Delo=K*S_SIR
        end do
               
        call bezier_quadratic_coeff(xx,etadb,n,Ck) 
        
        do k=1,n-1
           deltabar(k)=(xx(k+1)-xx(k))/3.*(etadb(k)+etadb(k+1)+Ck(k)) !tau_(k+1)-tau_k
        end do

        taunui(1)=dble(taue(n)) !superficie
        do k=1,n-1
           taunui(k+1)=taunui(k)+deltabar(k)
        end do
        
c evaluation of the coefficients of eq 29-30-31 of  
c de la Cruz & Piskunov (2013) ApJ 764, 33        
c https://iopscience.iop.org/article/10.1088/0004-637X/764/1/33/pdf       

         solp(1,n)=dble(S(1))
         do i=2,4
           solp(i,n)=0.
         end do 
         
c        call derivaarray(taunui,Sp,n,Spder)                  !using classical derivative of the Source function
         call bezier_cubic_deriv_array4(taunui,Sp,n,Spder)    !using Bezier3 derivative of the Source function
c         call derivaarray44(taunui,K4,n,K4der)               !using classical derivative of the Asorption matrix
c         call bezier_cubic_deriv_array44(taunui,K4,n,K4der)  !using Bezier3 derivative of the Asorption matrix
         call bezier_cubic_deriv_absormodf(taunui,K4,n,K4der) !using Bezier3 derivative of the Asorption matrix (ONLY 6 values)
                                                                                                          
        do i=1,4
           do j=1,4
             K4_k(i,j)=K4(i,j,n)
           end do
        end do  
        call matrizmult(K4_k,K4_k,K42_k)

        do k=n-1,1,-1  !para obtener Ik a partir de Ik+1
           dk=deltabar(k)
           dk3=dk*dk*dk
           dk_3=dk/3.d0

           do j=1,4
              solpk1(j)=solp(j,k+1)           !I(k+1),Q(k+1)...V(k+1)
              do i=1,4
                 K42_k1(i,j)=K42_k(i,j)       ! K(k+1)*K(k+1)
                 K4_k1(i,j)=K4_k(i,j)         ! K(k+1)     
                 K4_k(i,j)=K4(i,j,k)          ! K(k)    
              end do
           end do    
c           call matrizmult(K4_k,K4_k,K42_k)   ! K(k)*K(k)
           call absormodfsquare(K4_k,K42_k)
           
c e_k vector of equation 29. First part of E vector definition.  
           do i=1,4
             suma1=0.d0
             suma2=0.d0
             do j=1,4
                 suma1= suma1 +K4_k(i,j) *Sp(j,k)      !K_k*S_k
                 suma2= suma2 +K4_k1(i,j)*Sp(j,k+1)    !K_k+1*S_k+1
             end do
             e_k(i)=   dk_3*(suma1+Spder(i,k))  +Sp(i,k)
             e_k1(i)= -dk_3*(suma2+Spder(i,k+1))+Sp(i,k+1)
           end do

c e^bar_k matrix of equation 29. Appears in the second part of E vector definition.  
           do j=1,4
              do i=1,4
                 ebar_k(i,j) =-(dk_3*(K42_k(i,j) +K4der(i,j,k)  +K4_k(i,j)) +K4_k(i,j) )
                 ebar_k1(i,j)= (dk_3*(K42_k1(i,j)+K4der(i,j,k+1)+K4_k1(i,j))-K4_k1(i,j))
              end do
           end do 		         
            
c F vector of equation 30
           call matrizmultvec(ebar_k1,solpk1,fkvec)  !right term of equation 30 [ ]I_k+1
           do i=1,4
              fkvec(i)=e_k1(i)+fkvec(i)
           end do   

c S gotica (S-K*I)(k+1) vector of equation 31
           call matrizmultvec(K4_k1,-solpk1,Sgotvec)
           do i=1,4
              Sgotvec(i)=Sp(i,k+1)+Sgotvec(i)
           end do   
                     
           if(dk .gt. 6.d-2)then    
              Exk=dexp(-dk)
              alphak=(-6.d0+dk*(6.d0+dk*(-3.d0+dk))+6.d0*Exk)/dk3
              betak=(6.d0+(-6.d0-dk*(6.d0+dk*(3.d0+dk)))*Exk)/dk3
              gammak=3.d0*(6.d0+(dk-4.d0)*dk-2.d0*(dk+3.d0)*Exk)/dk3
              epsik=3.d0*(2.d0*dk-6.d0+(6.d0+dk*(dk+4.d0))*Exk)/dk3
           else
              Exk=1.-dk*(1.d0-dk*(0.5d0-0.166667d0*dk))
              alphak=dk*(0.25d0-dk*(0.05d0-dk*(1.d0/120.d0-dk/840.d0)))
              betak= dk*(0.25d0-dk*(0.2d0-dk*(1.d0/12.d0-dk/42.d0)))
              gammak=dk*(0.25d0-dk*(0.1d0-dk*(1.d0/40.d0-dk/210.d0)))
              epsik= dk*(0.25d0-dk*(3.d0/20.d0-dk*(1.d0/20.d0-dk/84.d0)))
           end if 
           
c A matrix of equation 31 and B vector of equation 31
c Bmat= matrix multiplying I_k+1  in eq 31 (and 30, 29)
          do j=1,4
             do i=1,4
                Ak(i,j)=alphak*K4_k(i,j)-gammak*ebar_k(i,j)
                Bmat(i,j)=epsik*ebar_k1(i,j)-betak*K4_k1(i,j)
             end do
             Ak(j,j)=1.+Ak(j,j)
             Bmat(j,j)=Bmat(j,j)+Exk
             Bk(j)=Exk*solpk1(j)+alphak*Sp(j,k)+betak*Sgotvec(j)+gammak*e_k(j)+epsik*fkvec(j)    
          end do   
                   
          call matinx3(Ak,Akinv,ierror)
c           call linearsystem(Ak,Bk,Xk)
          call matrizmultvec(Akinv,Bk,Xk)
          do i=1,4
             solp(i,k)=Xk(i)
             sol(i,n-k+1)=real(Xk(i))
          end do  
                  
          call matrizmult(Akinv,Bmat,OOk)  

          do i=1,4
             do j=1,4
                OOkr(i,j,k+1)=OOk(i,j)       !O(k,k+1): I(tau_k)=O(tau_k,tau_k+1)*I(tau_k+1) from k=n-1 to k=1
             end do
          end do  
                 
        end do  !do en alturas (k)
        
        do j=1,4
           do i=1,4
              OOfromkto1(i,j,1)=0.d0            !O(1,1)  1-->1
           end do
           OOfromkto1(j,j,1)=1.d0 
        end do     
        do i=1,4
           do j=1,4
              OOfromkto1(i,j,2)=OOkr(i,j,2)     !O(1,2)  2-->1
           end do
        end do        
        do k=3,n
           do i=1,4
             do j=1,4
                OOk(i,j)=OOkr(i,j,k)               ! O(k-1,k) 
                OO1k_1(i,j)=OOfromkto1(i,j,k-1)    ! O(1,k-1)
             end do
           end do  
           call matrizmult(OO1k_1,OOk,OO1k)   ! O(1,k-1)*O(k-1,k)=O(1,k)
           do i=1,4
             do j=1,4
                OOfromkto1(i,j,k)=OO1k(i,j)  !O(1,k)
             end do
           end do  
        end do
        
        do ii=1,4
	   svec(ii) = sol(ii,n)
	enddo
        do k=1,n
           do j=1,4
              do i=1,4
                 OO(i,j,k)=real(OOfromkto1(i,j,n-k+1))  !now O(tau_k,tau_k-1): I(tau_k)=O(tau_k,tau_k-1)*I(tau_k-1) from k=n-1 to k=1
              end do   
           end do   
        end do
		
	if(ierror.eq.1)then
           call error(KSTOP,'delo_bezier3',
     &                'Error in the integration of the RTE')
	endif
	
	do i=1,n
	   s0(i) = 0.e0
	   sol(1,i) = sol(1,i)-s(i)
	enddo
	
	
	if(mnodos(1).ne.0)call delta1(n,tk,bt,sol,OO,etall,rt)
	if(mnodos(1).ne.0.or.mnodos(2).ne.0.or.mnodos(4).ne.0.or.
     &     mnodos(6).ne.0)call delta1(n,pk,s0,sol,OO,etall,rp)
	if(mnodos(3).ne.0)call delta1(n,mk,s0,sol,OO,etall,rm)
	if(mnodos(4).ne.0)call delta1(n,hk,s0,sol,OO,etall,rh)
	if(mnodos(5).ne.0)call delta1(n,vk,s0,sol,OO,etall,rv)
	if(mnodos(6).ne.0)call delta1(n,gk,s0,sol,OO,etall,rg)
	if(mnodos(7).ne.0)call delta1(n,fk,s0,sol,OO,etall,rf)

	return
	end

c ___________________________________________________________________
       subroutine inicializa(icontorno,etal,s,ds,oo,sol,ierror)

       implicit real*4 (a-h,o-z)

       include 'PARAMETER' !por kt
       parameter (kt16=16*kt)
       real*4 etal(4,4,kt),inveta(4,4),s(*),ds(*),oo(4,4,kt),sol(4,kt)

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
