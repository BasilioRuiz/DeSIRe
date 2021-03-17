c	gasb calcula las preiones parciales y sus derivadas con
c	respectoa la temperatura y la presion.es una modificacion
c	de la rutina gas.

c	dp es la derivada de p respecto a t, ddp respecto a pe.

      subroutine gasb_db(theta_r4,pe_r4,p_r4,dp_r4,ddp_r4)
c	...............................................................
      implicit integer*4 (i-n)
      implicit real*8 (a-h, o-z)
      include 'PARAMETER'
      parameter (ncontr=28)
        real*4 theta_r4,pe_r4,t_r4
        real*4 p_r4(99),dp_r4(99),ddp_r4(99)
        real*4 weight_r4,weighti_r4,alfaii_r4,chi1i_r4,chi2i_r4
        real*4 cmol(91),dcmol(91)
        real*4 alfai_r4(ncontr),chi1_r4(ncontr),chi2_r4(ncontr)
        real*4 u0_r4(ncontr),u1_r4(ncontr),u2_r4(ncontr)
        real*4 du0_r4(ncontr),du1_r4(ncontr),du2_r4(ncontr)

        real*8 alfai(ncontr)
	real*8 p(99),dp(99),ddp(99)
	real*8 abundance(ncontr),totalbundance,XA(ncontr),weightarr(ncontr)
	real*8 musuminv,mu,elec_gram,uma,Kboltz,R,ne,rho
	
	t_r4=5040./theta_r4
        t=5040.d0/theta_r4
        pe=1.d0*pe_r4

      if(t .gt. 1.5d4)then
         uma=1.660538921d-24
         Kboltz=1.38064852d-16
         R=8.3144598d7                   !cte de los gases 
         musuminv=0.
         elec_gram=0.
         do i=1,ncontr
            p(i)=0.                      ! p(h)/p(h')
	    dp(i)=0.	                 !dlf1/f1
	    ddp(i)=0.	                 !ddlf1/f1
      	    iii=i
    	    call neldatb(iii,0.,weighti_r4,alfaii_r4,chi1i_r4,chi2i_r4)
    	    chi1_r4(i)=chi1i_r4
    	    chi2_r4(i)=chi2i_r4
    	    weightarr(i)=weighti_r4*1.d0
    	    abundance(i)=weighti_r4*alfaii_r4*1.d0
    	    totalbundance=totalbundance+abundance(i)
    	end do
    	do i=1,ncontr
    	   abundance(i)=abundance(i)/totalbundance  !grams element i/grams material= Xi
    	   XA(i)=abundance(i)/weightarr(i)             !mol element i/grams material 
    	   musuminv=musuminv+(i+1)*XA(i)
    	   elec_gram=elec_gram+i*XA(i)
       	end do
       	mu=1./musuminv                              !mean molecular weight
       	elec_gram=elec_gram/uma
       	ne=pe/Kboltz/t 
       	rho=ne/elec_gram
       	pg=rho*R*t/mu
       	
        p(84)=pg                          ! gas pressure 
          dp(84)=0.
	  ddp(84)=1./pe                   !derivative of logarithm
        p(85)= XA(1)/uma*rho*Kboltz*T     ! p(h')
	  dp(85)=0.		
	  ddp(85)=1./pe			
        p(86)=1.                          ! p(h+)/p(h')	
          dp(86)=0.
	  ddp(86)=0.
        p(87)=0.                          ! p(h-)/p(h')
	  dp(87)=0.
	  ddp(87)=0.
        p(88)=0.                          ! p(h2+)/p(h')
	  dp(88)=0.		
	  ddp(88)=0.		
        p(89)=0.                          ! p(h2)/p(h')
	  dp(89)=0.		
	  ddp(89)=0.	                  !ddlf5/f5
        p(90)=pe/p(85)                    ! pe/p(h') 
	  dp(90)=0.		
	  ddp(90)=0.		

        p(91)=pe/(Kboltz*T)               ! n(e)=pe/kt
	  dp(91)=-1./t
	  ddp(91)=1./pe
      else
      
      if(pe.le.0)then
         call error(KWARN,'gasb_db','Negative values of the electron pressure'
     &   //         ' have been found.\n These are being changed to 1.e-10')
         pe=1.d-10
         g4=0.
         g5=0.
         dg4=0.
         ddg4=0.
         dg5=0.
         ddg5=0.
      else
         call molecb(theta_r4,cmol,dcmol)
         do i=1,2
            call acota(cmol(i),-30.,30.)
         end do
         g4=pe*10.d0**(cmol(1))
         dg4=dcmol(1)*dlog(10.d0)
         ddg4=1./pe
         g5=pe*10.d0**(cmol(2))
         dg5=dcmol(2)*dlog(10.d0)
	 ddg5=1./pe
      end if


c ahora leo las funciones de particion u0,u1,u2 y sus derivadas
      do 5 i=1,ncontr
      	  iii=i
          call neldatb(iii,0.,weight_r4,alfai_r4(i),
     *                 chi1_r4(i),chi2_r4(i))
          alfai(i)=alfai_r4(i)*1.d0
5     end do

      do 4 i=1,ncontr
      	  iii=i
     	  call nelfctb(iii,t_r4,u0_r4(iii),u1_r4(iii),u2_r4(iii),
     *du0_r4(iii),du1_r4(iii),du2_r4(iii))
4     end do      
      g2=saha_db(theta_r4,chi1_r4(1),u0_r4(1),u1_r4(1),pe_r4)   ! p(h+)/p(h)
	dg2=dsaha_db(theta_r4,chi1_r4(1),du0_r4(1),du1_r4(1))   !derivo log(g2) con t
	ddg2=-1./pe 
                           !derivada log(g2) con pe
      p(92)=g2
	dp(92)=dg2
	ddp(92)=ddg2
      g3=saha_db(theta_r4,0.754,1.,u0_r4(1),pe_r4)   ! p(h)/p(h-)

c      call acota_db(g3,1.d-30,1.d30)
      call acota_db(g3,1.d-38,1.d38)

      g3=1.d0/g3                              ! p(h-)/p(h) 
	dg3=-1.*dsaha_db(theta_r4,0.754,0.,du0_r4(1))
	ddg3=1./pe
      p(93)=g3
	dp(93)=dg3
	ddp(93)=ddg3
      g1=0.
	dlg1=0.       !las dl son derivadas no de log(g) sino de g
	ddlg1=0.
	
c      print*,'gasb_db 80 p(93)=',p(93),	ncontr

      do 1 i=2,ncontr
      a=saha_db(theta_r4,chi1_r4(i),u0_r4(i),u1_r4(i),pe_r4)
	da=dsaha_db(theta_r4,chi1_r4(i),du0_r4(i),du1_r4(i))
	dda=-1./pe
      b=saha_db(theta_r4,chi2_r4(i),u1_r4(i),u2_r4(i),pe_r4)
	dlb=b*dsaha_db(theta_r4,chi2_r4(i),du1_r4(i),du2_r4(i))
	ddlb=-b/pe
      c=1.+a*(1.+b)
c      call acotasig_db(c,1.d-20,1.d20)
      call acotasig_db(c,1.d-38,1.d38)

      p(i)=alfai(i)/c   ! p/ph' for neutral he,li, ...
	dp(i)=-(a*da*(1.+b)+a*dlb)/c

	ddp(i)=-(a*dda*(1.+b)+a*ddlb)/c
        ss1=(1.+2.*b)
c        call acotasig_db(ss1,1.d-20,1.d20)
        call acotasig_db(ss1,1.d-38,1.d38)
   
	ss=p(i)*a*ss1
        dss=dp(i)+da+2.*dlb/ss1
c        print*,'gasb_db 163 ',dp(i),da,ss1,dlb,dlb/ss1

	ddss=ddp(i)+dda+2.*ddlb/ss1

c	g1=g1+p(i)*a*ss1
	g1=g1+ss
	dlg1=dlg1+ss*dss             !ojo estas no son derivadas del log
c	print*,'gasb_db 168 ',ss,dss,dlg1
1	ddlg1=ddlg1+ss*ddss
        a=1.+g2+g3
	dla=g2*dg2+g3*dg3
	ddla=g2*ddg2+g3*ddg3


	if(g5.lt.1.d-35)then
           call error(KSTOP,'gasb_db','The electronic pressure is too small.'
     &     //         ' Check the atmospheric models')
	end if

c      call acotasig_db(g5,1.d-20,1.d20)
      call acotasig_db(g5,1.d-38,1.d38)

      b=2.*(1.+g2/g5*g4)
	dlb=(b-2.)*(dg2-dg5+dg4)
	ddlb=(b-2.)*(ddg2-ddg5+ddg4)
      c=g5
	dlc=dg5*g5
	ddlc=ddg5*g5
      d=g2-g3
	dld=g2*dg2-g3*dg3
	ddld=g2*ddg2-g3*ddg3
      e=g2/g5*g4
	de=dg2-dg5+dg4
	dde=ddg2-ddg5+ddg4
	dle=e*de
	ddle=e*dde

c	call acotasig_db(a,1.d-15,1.d15)
c	call acotasig_db(b,1.d-15,1.d15)
c	call acotasig_db(c,1.d-15,1.d15)
c	call acotasig_db(d,1.d-15,1.d15)
c	call acotasig_db(e,1.d-15,1.d15)
	call acotasig_db(a,1.d-25,1.d25)
	call acotasig_db(b,1.d-25,1.d25)
	call acotasig_db(c,1.d-25,1.d25)
	call acotasig_db(d,1.d-25,1.d25)
	call acotasig_db(e,1.d-25,1.d25)

      c1=c*b**2+a*d*b-e*a**2
	dlc1=dlc*b*b+(c*2.*b+a*d)*dlb+dla*(d*b-2.*e*a)+dld*a*b-dle*a*a
	ddlc1=ddlc*b*b+(c*2.*b+a*d)*ddlb+ddla*(d*b-2.*e*a)+ddld*a*b-ddle*a*a
		
      c2=2.*a*e-d*b+a*b*g1
	dlc2=dla*(2.*e+b*g1)+dlb*(a*g1-d)-dld*b+dle*2.*a+a*b*dlg1
	ddlc2=ddla*(2.*e+b*g1)+ddlb*(a*g1-d)-ddld*b+ddle*2.*a+a*b*ddlg1
      c3=-(e+b*g1)
	dlc3=-dle-dlb*g1-b*dlg1
	ddlc3=-ddle-ddlb*g1-b*ddlg1
           

c      call acotasig_db(c1,1.d-15,1.d15)
      call acotasig_db(c1,1.d-25,1.d25)

      f1=0.5*c2/c1
c      print*,'gasb_db 220',f1,c1,c2

c	dlf1=.5*(dlc2/c1-(dlc1*c2)/(c1*c1))
      dc1=dlc1/c1  !!!!!!!!
      dc2=dlc2/c2
      dlf1=f1*(dc2-dc1)  

c	ddlf1=.5*(ddlc2/c1-(ddlc1*c2)/(c1*c1))
      ddc1=ddlc1/c1  !!!!!!!!
      ddc2=ddlc2/c2
      ddlf1=f1*(ddc2-ddc1)             


	dlf1=-dlf1+sign(1.d0,c1)*(2.*f1*dlf1-dlc3/c1+dlc1*c3/(c1*c1))
     *	/(2.*dsqrt(f1**2-c3/c1))

	ddlf1=-ddlf1+sign(1.d0,c1)*(2.*f1*ddlf1-ddlc3/c1+ddlc1*c3/(c1*c1))
     *	/(2.*dsqrt(f1**2-c3/c1))

      f1=-f1+sign(1.d0,c1)*dsqrt(f1**2-c3/c1)
c        print*,'gasb_db 238',f1,c3,c1,f1**2-c3/c1

      f5=(1.-a*f1)/b
	if(dabs(f5).lt.1.d-30)then
	   dlf5=0.
	   ddlf5=0.
	   df5=0.
	   ddf5=0.
	else
	   dlf5=(-dla*f1-a*dlf1)/b-((1.-a*f1)*dlb)/(b*b)
	   ddlf5=(-ddla*f1-a*ddlf1)/b-((1.-a*f1)*ddlb)/(b*b)
	   df5=dlf5/f5
	   ddf5=ddlf5/f5
	end if
        f4=e*f5
	if(dabs(f4).lt.1.d-30)then
	   dlf4=0.
	   ddlf4=0.
	   df4=0.
	   ddf4=0.
	else
	   dlf4=f5*dle+e*dlf5
	   ddlf4=f5*ddle+e*ddlf5
	   df4=dlf4/f4
	   ddf4=ddlf4/f4
	end if
      

      f3=g3*f1
	dlf3=f3*dg3+g3*dlf1
	ddlf3=f3*ddg3+g3*ddlf1
      f2=g2*f1
	dlf2=f2*dg2+g2*dlf1
	ddlf2=f2*ddg2+g2*ddlf1

	if(abs(f1).lt.1.d-30)then
	   divf1=0.
	   ddivf1=0.
	else
	   divf1=dlf1/f1
	   ddivf1=ddlf1/f1
	end if
	dlf200=dg2+divf1
	ddlf200=ddg2+ddivf1


        fe=f2-f3+f4+g1
	if(abs(fe).lt.1.e-30)then
	   dlfe=0.
	   ddlfe=0.
	   dfe=0.
	   ddfe=0.
	else
	   dlfe=dlf2-dlf3+dlf4+dlg1
	   ddlfe=ddlf2-ddlf3+ddlf4+ddlg1
	   dfe=dlfe/fe
	   ddfe=ddlfe/fe
	end if

c      call acotasig(f1,1.e-20,1.)
       call acotasig_db(f2,1.d-25,1.d0)
c      call acotasig(f3,1.e-20,1.)
c      call acotasig(f4,1.e-20,1.)
c      call acotasig(f5,1.e-20,1.)
      call acotasig_db(fe,1.d-25,2.d0)
      
      phtot=pe/fe
	dlphtot=-phtot*dfe
	ddlphtot=(1.-ddlfe*phtot)/fe
	dphtot=-dfe
	ddphtot=1./pe-ddfe

        if(f5.gt.1.d-4) goto 2
          	const6=g5/pe*f1**2
		dconst6=dg5+2.*divf1		!dlf1/f1
		ddconst6=ddg5-1./pe+2.*ddivf1	!ddlf1/f1
      		const7=f2-f3+g1
c                call acotasig(const7,1.e-15,2.)
		dlconst7=dlf2-dlf3+dlg1
		ddlconst7=ddlf2-ddlf3+ddlg1
		dconst7=dlconst7/const7
		ddconst7=ddlconst7/const7

     		 do 3 i=1,5
      			f5=phtot*const6
			df5=dphtot+dconst6
			ddf5=ddphtot+ddconst6
		      	f4=e*f5
			df4=de+df5
			ddf4=dde+ddf5
		      	fe=const7+f4
                        call acotasig_db(fe,1.d-25,2.d0)
c                        call acotasig(f4,1.e-20,1.)
c                        call acotasig(f5,1.e-20,1.)
			dfe=(dlconst7+df4*f4)/fe
			ddfe=(ddlconst7+ddf4*f4)/fe
			phtot=pe/fe
		        dphtot=-dfe
3       		ddphtot=1./pe-ddfe

		dlf5=df5*f5
		ddlf5=ddf5*f5
		dlf4=df4*f4
		ddlf4=ddf4*f4
		dlfe=dfe*fe
		ddlfe=ddfe*fe
		dlphtot=dphtot*phtot
		ddlphtot=ddphtot*phtot
2     pg=pe*(1.+(f1+f2+f3+f4+f5+0.1014)/fe)
        
        dlpg=pe*(dlf1+dlf2+dlf3+dlf4+dlf5)/fe-(pg-pe)*(dlfe/fe)
	ddlpg=pe*(ddlf1+ddlf2+ddlf3+ddlf4+ddlf5)/fe-(pg-pe)*(ddlfe/fe)
	ddlpg=ddlpg+pg/pe
      p(1)=f1   ! p(h)/p(h')
	dp(1)=divf1		!dlf1/f1
	ddp(1)=ddivf1		!ddlf1/f1

      call acotasig_db(pg,1.d-25,1.d25)
      p(84)=pg   ! gas pressure
        dp(84)=dlpg/pg
	ddp(84)=ddlpg/pg
	if(t .gt. 2.d4)then
	   f2=1.0
	   dlf200=0.
	   ddlf200=0.
c	   phtot=4.e-2
c	   dphtot=0.
c	   ddphtot=1./pe
	end if   
      p(85)=phtot   ! p(h')
	dp(85)=dphtot			!dlphtot/phtot
	ddp(85)=ddphtot			!ddlphtot/phtot
      p(86)=f2   ! p(h+)/p(h')	
        dp(86)=dlf200
	ddp(86)=ddlf200
      p(87)=f3   ! p(h-)/p(h')
	if(abs(f3).lt.1.d-30)then
	   divf3=0.
	   ddivf3=0.
	else
	   divf3=dlf3/f3
	   ddivf3=ddlf3/f3
	end if

	dp(87)=divf3
	ddp(87)=ddivf3
      p(88)=f4   ! p(h2+)/p(h')
	dp(88)=df4		!dlf4/f4
	ddp(88)=ddf4		!ddlf4/f4
      p(89)=f5   ! p(h2)/p(h')
	dp(89)=df5		!dlf5/f5
	ddp(89)=ddf5		!ddlf5/f5
      p(90)=fe   ! pe/p(h')
	dp(90)=dfe		!dlfe/fe
	ddp(90)=ddfe		!ddlfe/fe

      p(91)=pe/(1.38054e-16*t) ! n(e)=pe/kt
	dp(91)=-1./t
	ddp(91)=1./pe
c      print*,'gasb fin'
      do i=1,ncontr
         p_r4(i)=real(p(i))
         dp_r4(i)=real(dp(i))
         ddp_r4(i)=real(ddp(i))
      end do
      do i=84,93
         p_r4(i)=real(p(i))
         dp_r4(i)=real(dp(i))
         ddp_r4(i)=real(ddp(i))
      end do

      end if

c      print*,'gasb_db 432 p(93)=',p(93),p_r4(93)
      
      return
      end
c	...............................................................
