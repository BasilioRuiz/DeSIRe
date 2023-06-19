c Similar to delo_bezier3 but integrating in a log(tau_nu) scale
c instead to in tau_nu scale (as delo_bezier3 does)
c delo_bezier3  supposes that Sgotic=S-(K/eta-1)*I is a cubic function of tau_nu between two tau_nu points
c delo_bezier3log supposes tau*exp(-tau)*Sgotic is a cubic function of tau_nu between two tau_nu points

	subroutine delo_bezier3log(s,ds,etall,n,svec,nn,bt,tk,pk,hk,
     &           vk,gk,fk,mk,rt,rp,rh,rv,rg,rf,rm,mnodos)

        implicit real*4 (a-h,o-z)
	include 'PARAMETER'
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
        real*8 taunui(kt),logtaunui(kt),taui0,taui1
        real*8 xx(kt),etadb(kt),Ck(kt),etaii
        real*8 Sp(4,kt),solp(4,kt),deltabar(kt),solpk1(4)
        real*8 K4(4,4,kt),K4der(4,4,kt)
        real*8 K4_k(4,4),K4_k1(4,4),K42_k(4,4),K42_k1(4,4)
        real*8 e_k(4),e_k1(4),ebar_k(4,4),ebar_k1(4,4),fkvec(4),suma1,suma2  
        real*8 Spder(4,kt),Sgotvec(4)
        real*8 dk,h,Exk,alphak,betak,gammak,epsik
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

        taunui(1)=dble(taue(n)) !surface
        logtaunui(1)=dlog(taunui(1))
        do k=1,n-1
           taunui(k+1)=taunui(k)+deltabar(k)
           logtaunui(k+1)=dlog(taunui(k+1))
        end do
        
c evaluation of the coefficients of eq 29-30-31 of  
c de la Cruz & Piskunov (2013) ApJ 764, 33        
c https://iopscience.iop.org/article/10.1088/0004-637X/764/1/33/pdf       
         solp(1,n)=dble(S(1))
         do i=2,4
           solp(i,n)=0.
         end do 
        
c Now derivatives are to respect logtaunui    
c         call derivaarray(logtaunui,Sp,n,Spder)                  !using classical derivative of the Source function
          call bezier_cubic_deriv_array4(logtaunui,Sp,n,Spder)    !using Bezier3 derivative of the Source function
c         call derivaarray44(logtaunui,K4,n,K4der)               !using classical derivative of the Asorption matrix
c         call bezier_cubic_deriv_array44(logtaunui,K4,n,K4der)  !using Bezier3 derivative of the Asorption matrix
          call bezier_cubic_deriv_absormodf(logtaunui,K4,n,K4der) !using Bezier3 derivative of the Asorption matrix (ONLY 6 values)
                                                          
        do i=1,4
           suma1=0.d0
           do j=1,4
              K4_k(i,j)=K4(i,j,n)
              suma1= suma1 +K4_k(i,j)*Sp(j,n)
           end do
           e_k(i)=suma1-Spder(i,n)  ! antes suma1+Spder(i,n) 
        end do  
c        call matrizmult(K4_k,K4_k,K42_k)
        call absormodfsquare(K4_k,K42_k)
        
        do i=1,4
           do j=1,4
              ebar_k(i,j)=-(-K4der(i,j,n)+taunui(n)*(K42_k(i,j)+K4_k(i,j)))  !antes -(K4der(i,j,n)
           end do
        end do         

        do k=n-1,1,-1  !para obtener Ik a partir de Ik+1
           dk=deltabar(k)
           taui0=taunui(k)
           taui1=taunui(k+1)
           h=logtaunui(k+1)-logtaunui(k)

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
c e^bar_k matrix of equation 29. Appears in the second part of E vector definition. 
           do i=1,4
             suma1=0.d0
             do j=1,4
                 suma1= suma1 +K4_k(i,j) *Sp(j,k)      !K_k*S_k
                 ebar_k1(i,j)=ebar_k(i,j)
                 ebar_k(i,j)=-(-K4der(i,j,k)+taui0*(K42_k(i,j)+K4_k(i,j))) !antes -(K4der(i,j,k)+
             end do
             e_k1(i)=e_k(i)
             e_k(i)= taui0*suma1-Spder(i,k)      !antes   taui0*suma1+Spder(i,k)       
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
            else
               Exk=1.-dk*(1.d0-dk*(0.5d0-0.166667d0*dk))
            end if 
            alphak=taui0*h*(0.5d0+h/12.d0*(1.d0-taui0))
            betak =taui1*h*(0.5d0-h/12.d0*(1.d0-taui1))*Exk
            gammak=taui0*h*h/12.d0
            epsik=-taui1*h*h/12.d0*Exk     
           
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
                           
           call matinx3 (Ak,Akinv,ierror)
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
           call error(KSTOP,'delo_bezier3log',
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
		
c______________________________________________________________________________________________________
c Similar to delo_bezier3 but integrating in a log(tau_nu) scale
c instead to in tau_nu scale (as delo_bezier3 does)
c delo_bezier3  supposes that Sgotic=S-(K/eta-1)*I is a cubic function of tau_nu between two tau_nu points
c delo_bezier3log supposes tau*exp(-tau)*Sgotic is a cubic function of tau_nu between two tau_nu points

	subroutine delo_bezier3logcont(s,ds,etall,n,svec,nn,bt,tk,pk,hk,
     &           vk,gk,fk,mk,rt,rp,rh,rv,rg,rf,rm,mnodos)

        implicit real*4 (a-h,o-z)
	include 'PARAMETER' !por kt
        parameter (kt16=16*kt)
        parameter (dl10=2.3025851)

	integer*4 n,nn,i,j,k,kr,ierror,jj
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
        integer*4 kcontorno
        real*8 taunui(kt),logtaunui(kt),taui0,taui1
        real*8 xx(kt),etadb(kt),Ck(kt),etaii
        real*8 Sp(4,kt),solp(4,kt),deltabar(kt),solpk1(4)
        real*8 K4(4,4,kt),K4der(4,4,kt)
        real*8 K4_k(4,4),K4_k1(4,4),K42_k(4,4),K42_k1(4,4)
        real*8 e_k(4),e_k1(4),ebar_k(4,4),ebar_k1(4,4),fkvec(4),suma1,suma2  
        real*8 Spder(4,kt),Sgotvec(4)
        real*8 dk,h,Exk,alphak,betak,gammak,epsik
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

        taunui(1)=dble(taue(n)) !surface
        logtaunui(1)=dlog(taunui(1))
        kcontorno=n
        do k=1,n-1
           taunui(k+1)=taunui(k)+deltabar(k)
           logtaunui(k+1)=dlog(taunui(k+1))
           if(logtaunui(k+1) .gt. 10. .and. kcontorno .eq. n)kcontorno=k+1
        end do
        if (kcontorno .le. 5)kcontorno=5
c        kcontorno=n  !BORRAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
c evaluation of the coefficients of eq 29-30-31 of  
c de la Cruz & Piskunov (2013) ApJ 764, 33        
c https://iopscience.iop.org/article/10.1088/0004-637X/764/1/33/pdf     
        do k=n,kcontorno,-1
           solp(1,k)=dble(S(n-k+1))
           sol(1,n-k+1)=S(n-k+1)
           do i=2,4
             solp(i,k)=0.d0
             sol(i,n-k+1)=0.
           end do 
        end do   
        
c Now derivatives are to respect logtaunui    
c         call derivaarray(logtaunui,Sp,n,Spder)                  !using classical derivative of the Source function
c          call bezier_cubic_deriv_array4(logtaunui,Sp,n,Spder)    !using Bezier3 derivative of the Source function
          call bezier_cubic_deriv_array4ini(logtaunui,Sp,1,kcontorno,Spder)    !using Bezier3 derivative of the Source function
c         call derivaarray44(logtaunui,K4,n,K4der)               !using classical derivative of the Asorption matrix
c         call bezier_cubic_deriv_array44(logtaunui,K4,n,K4der)  !using Bezier3 derivative of the Asorption matrix
c          call bezier_cubic_deriv_absormodf(logtaunui,K4,n,K4der) !using Bezier3 derivative of the Asorption matrix (ONLY 6 values)
          call bezier_cubic_deriv_absormodfini(logtaunui,K4,1,kcontorno,K4der) !using Bezier3 derivative of the Asorption matrix (ONLY 6 values)
                                                          
        do i=1,4
           suma1=0.d0
           do j=1,4
              K4_k(i,j)=K4(i,j,n)
              suma1= suma1 +K4_k(i,j)*Sp(j,n)
           end do
           e_k(i)=suma1+Spder(i,n)
        end do  
c        call matrizmult(K4_k,K4_k,K42_k)
        call absormodfsquare(K4_k,K42_k)
        
        do i=1,4
           do j=1,4
              ebar_k(i,j)=-(K4der(i,j,n)+taunui(n)*(K42_k(i,j)+K4_k(i,j)))
           end do
        end do         

        do k=kcontorno-1,1,-1  !para obtener Ik a partir de Ik+1
           dk=deltabar(k)
           taui0=taunui(k)
           taui1=taunui(k+1)
           h=logtaunui(k+1)-logtaunui(k)

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
c e^bar_k matrix of equation 29. Appears in the second part of E vector definition. 
           do i=1,4
             suma1=0.d0
             do j=1,4
                 suma1= suma1 +K4_k(i,j) *Sp(j,k)      !K_k*S_k
                 ebar_k1(i,j)=ebar_k(i,j)
                 ebar_k(i,j)=-(K4der(i,j,k)+taui0*(K42_k(i,j)+K4_k(i,j)))
             end do
             e_k1(i)=e_k(i)
             e_k(i)= taui0*suma1+Spder(i,k)              
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
            else
               Exk=1.-dk*(1.d0-dk*(0.5d0-0.166667d0*dk))
            end if 
            alphak=taui0*h*(0.5d0+h/12.d0*(1.d0-taui0))
            betak =taui1*h*(0.5d0-h/12.d0*(1.d0-taui1))*Exk
            gammak=taui0*h*h/12.d0
            epsik=-taui1*h*h/12.d0*Exk     
           
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
                           
           call matinx3 (Ak,Akinv,ierror)
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
        do k=3,kcontorno  !previously 3, n
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
	do k=1,kcontorno !previously 1, n
           do j=1,4
              do i=1,4
                 OO(i,j,k)=real(OOfromkto1(i,j,n-k+1))  !now O(tau_k,tau_k-1): I(tau_k)=O(tau_k,tau_k-1)*I(tau_k-1) from k=n-1 to k=1
              end do   
           end do   
        end do
        do k=kcontorno+1,n 
           do j=1,4
              do i=1,4
                 OO(i,j,k)=0.
              end do   
           end do   
        end do
  		
	if(ierror.eq.1)then
           call error(KSTOP,'delo_bezier3logcont',
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
		
c______________________________________________________________________________________________________
