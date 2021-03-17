	subroutine hermite_c(s,ds,etall,n,svec,nn,bt,tk,pk,
     &           vk,mk,rt,rp,rv,rm)

        implicit real*4 (a-h,o-z)
	include 'PARAMETER'   !por kt
        parameter (kt16=16*kt, kfac=10, ktnew=kfac*(kt-1)+1)
        parameter (dl10=2.3025851)

	integer n,nn
	real*4 svec

c internas
        real*4 sol(kt),tau(kt),taue(kt)
        real*4 deltae(kt),deltai(kt),delt2i(kt)
	real*4 ds(*), etall(*)
	real*4 s(*)
	real*4 OO(kt)
        real*4 tk(kt),pk(kt),vk(kt),mk(kt)
	real*4 bt(kt)
	real*4 rt(*),rp(*),rv(*),rm(*)
	
c interpolacion
        real*4 xa(11),ya(11)
	real*4 taunew(ktnew),tauenewxx(ktnew)
	real*4 etallnew(ktnew),deltainew(ktnew),delt2inew(ktnew)
	real*4 etalnew(ktnew),OOnew(ktnew)
	real*4 solnew(ktnew),snew(ktnew),dsnew(ktnew)
	integer iold(kt)	

        common/segunda/tau,taue,deltae,deltai,delt2i

c se calcula donde poner el contorno
        icontorno=1
	do i=1,n-1
           taunu=abs(etall(i))*deltae(i)
           if(taunu.gt.5.)icontorno=i
        end do

        if(abs(etall(icontorno)).le.10.)then
           interpola=0  
           nnew=n
           do j=1,nnew
              taunew(j)=tau(j)
              xx=taue(j)*dl10
              deltainew(j)=deltai(j)
              delt2inew(j)=delt2i(j)         
	      etallnew(j)=etall(j)
	      etalnew(j)=etall(j)*xx
	      snew(j)=s(j)
	      dsnew(j)=ds(j)
	    enddo   
	    icontornonew=icontorno
         else              
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
               do k=1,ngrado+1
	          ya(k)=etall(n3+k)	               
	       enddo
	       call polint(xa,ya,ngrado+1,taunewi,yyy,dy)
	       etallnew(i)=yyy
	       etalnew(i)=yyy*tauenewxxi
	       do k=1,ngrado+1
	          ya(k)=oo(n3+k)	               
	       enddo
	       call polint(xa,ya,ngrado+1,taunewi,yyy,dy)
	       oonew(i)=yyy
	    enddo   
         endif 
 
         call inicializa_c(icontornonew,etalnew,snew,dsnew,oonew,solnew,ierror)

         call numden_c(icontornonew,nnew,ktnew,taunew,deltainew,delt2inew,etalnew,snew,dsnew,OOnew,solnew,svec)

         if( interpola .eq. 0)then
            do i=1,n
               sol(i)=solnew(i)
               OO(i)=oonew(i)	         
            enddo
         else
            do i=1,n
               ii=iold(i)
               sol(i)=solnew(ii)
               OO(i)=oonew(ii)	         
            enddo
         end if
 
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
c ________________________________________________________________________
       subroutine inicializa_c(icontorno,etal,s,ds,oo,sol,ierror)

       implicit real*4 (a-h,o-z)

       include 'PARAMETER' !por kt
       parameter (kt16=16*kt)
       real*4 etal(*),inveta,s(*),ds(*),oo(*),sol(*)

       do k=1,icontorno
	 inveta=0.
	 if(abs(etal(k)) .gt. 1.e-30)then inveta=1./etal(k)
	 sol(k)=inveta*ds(k)+s(k)
	 OO(k)=0.
       end do

       return
       end

c-------------------------------------------------------------------
       subroutine numden_c(icontorno,n,kt,tau,deltai,delt2i,etal,s,ds,OO,sol,svec)
       implicit real*4 (a-h,o-z)
              
       integer icontorno,n,i
       real*4 etal(kt),etal2(kt),detal(kt)
       real*4 s(kt), ds(kt),den, num, den2
       real*4 sol(kt),OO(kt),tau(kt)
       real*4 deltai(kt),delt2i(kt),svec
           
       	do i=icontorno,n
  	   etal2(i) = etal(i)*etal(i)
	enddo
	call derivacuad(etal,detal,n,tau)

 	do i=icontorno+1,n
           delta=deltai(i)
           delt2=delt2i(i)
           e0=etal(i-1)
	   e1=etal(i)
           s00=s(i-1)
           s01=s(i)
	   den =-( e0*s00 + e1*s01 )*
     &	                 delta + sol(i-1)
	   den =-(detal(i-1)*s00 - detal(i)*s01  +
     &                  e0*ds(i-1) -  e1*ds(i) +
     &                  etal2(i-1)*s00 - etal2(i)*s01 )*
     &                  delt2 + den

	   den2=e0*delta+
     & delt2*(detal(i-1)+etal2(i-1))
	   den=den + (e0*delta+
     & delt2*(detal(i-1)+etal2(i-1)))*sol(i-1)
	   num= - e1*delta +
     & delt2*(detal(i)+etal2(i))
	   num = 1.0 + num
	   den2 = 1.0 + den2
        
	   OO(i-1) = den2/num
	   sol(i) = den/num
	enddo

	svec = sol(n)
	OO(n) = 1.0

	do i=n-2,1,-1
	   OO(i)=OO(i+1)*OO(i)
	enddo
            
       return
       end
c----------------------------------------------------------------       

