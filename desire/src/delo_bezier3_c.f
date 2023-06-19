c---------------------------------------------------------------------------    
    	 subroutine delo_bezier3_c(s,ds,etall,n,svec,nn,bt,tk,pk,
     &           vk,mk,rt,rp,rv,rm)
    
c Assumes that source function behaves as cubic splines in tau_nu
c instead as linear in tau_nu as Delo does.     

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'   !por kt
        parameter (kt16=16*kt)
        parameter (dl10=2.3025851)

        integer n,nn,nni
        real*4 svec

c internas
        real*4 sol(kt),tau(kt),taue(kt)
        real*4 deltae(kt),deltai(kt),delt2i(kt)
        real*4 ds(*), detal(kt), etal(kt), etall(*)
        real*4 etal2(kt)
        real*4 s(*), den, num, den2
        real*4 OO(kt),OOr(kt)
        real*4 tk(kt),pk(kt),vk(kt),mk(kt)
        real*4 bt(kt)
        real*4 rt(*),rp(*),rv(*),rm(*)
c para la integracion
        real*8 taunui(kt),bezier_cubic_deriv
        real*8 xx(kt),etadb(kt),Ck(kt)
        real*8 Sp(kt),solp(kt),deltabar(kt)
        real*8 dk,dk3,Exk,alphak,betak,gammak,epsik
        real*8 Ek(kt),Fk1(kt)
c para la interpolacion de prueba
c        real*4 xi(2000), yi(2000), Eb(2000), Fb(2000)
                
c       deltae(i)=(taue(i)-taue(i+1))
c       luego    -deltae(i-1)=taue(i)-taue(i-1)
c       deltai(i) = (tau(i)-tau(i-1))/2.
c       delt2i(i) = delta(i)*delta(i)/3. 
        common/segunda/tau,taue,deltae,deltai,delt2i

c se calcula donde poner el contorno
     	icontorno=1
c se calcula la derivada de la matriz de absorcion
	call derivacuad(etal,detal,n,tau)

c se integra
        call inicializa_c(icontorno,etal,s,ds,oo,sol,ierror)

        do i=1,n
           xx(i)=dble(taue(n-i+1))                !x in paper
           etadb(i)=dble(etall(n-i+1))            ! S in paper
        end do

        call bezier_quadratic_coeff(xx,etadb,n,Ck) 
        
        do i=1,n-1
           deltabar(i)=(xx(i+1)-xx(i))/3.*(etadb(i)+etadb(i+1)+Ck(i)) !tau_(k+1)-tau_k
        end do

        taunui(1)=dble(taue(n)) !superficie
        do i=1,n-1
           taunui(i+1)=taunui(i)+deltabar(i)
        end do
        
c evaluation of the coefficients of eq 20 of  
c de la Cruz & Piskunov (2013) ApJ 764, 33        
c https://iopscience.iop.org/article/10.1088/0004-637X/764/1/33/pdf
c We need to reverse the S array and integrate from n to 1        
        do i=1,n
           Sp(i)=dble(S(n-i+1))            ! S in paper
        end do
        ipcontorno=n-icontorno+1
        do i=n,ipcontorno,-1
           solp(i)=sol(n-i+1)
        end do 
        call bezier_cubic_coeff(taunui,Sp,n,Ek,Fk1)
        
        do k=ipcontorno+1,n
           OOr(k)=1.
        end do   
        do k=ipcontorno-1,1,-1  !para obtener Ik a partir de Ik+1
           dk=deltabar(k)
           dk3=dk*dk*dk
           
           if(dk .gt. 6.e-2)then           
              Exk=dexp(-dk)
              alphak=(-6.+dk*( 6.+dk*(-3.+dk))+6.*Exk)/dk3
              betak=(6.+(-6.-dk*(6.+dk*(3.+dk)))*Exk)/dk3
              gammak=3*(6+(dk-4.)*dk-2.*(dk+3.)*Exk)/dk3
              epsik=3*(2*dk-6.+(6+dk*(dk+4.))*Exk)/dk3
           else
              Exk=1.-dk*(1.-dk*(0.5-0.166667*dk))
              alphak=dk*(0.25-dk*(0.05-dk*(1./120-dk/840.)))
              betak= dk*(0.25-dk*(0.2-dk*(1./12.-dk/42.)))
              gammak=dk*(0.25-dk*(0.1-dk*(1./40.-dk/210.)))
              epsik= dk*(0.25-dk*(3./20.-dk*(1./20.-dk/84.)))
           end if   
           solp(k)=solp(k+1)*Exk+alphak*Sp(k)+betak*Sp(k+1)+
     &             gammak*Ek(k)+epsik*Fk1(k)
           OOr(k+1)=real(Exk)    !hay que darle la vuelta
           sol(n-k+1)=real(solp(k))
        end do
        do i=1,n
          OO(i)=OOr(n-i+1)
        end do  

        svec = real(solp(1))
	do i=n-2,1,-1
	   OO(i)=OO(i+1)*OO(i)
	enddo

c calculamos las f.resp
        do i=1,n-1
           ss=sol(i) - s(i)
           ooo=OO(i)
           rt(i)=ooo*(tk(i)*ss-etall(i)*bt(i))
           rp(i)=ooo*pk(i)*ss
           rv(i)=ooo*vk(i)*ss
           rm(i)=ooo*mk(i)*ss
        enddo
	    rt(n)=0.
	    rp(n)=0.
	    rv(n)=0.
	    rm(n)=0.
        
	     return
	     end    
    
    

       