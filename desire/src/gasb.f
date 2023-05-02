c=============================================================================
c     gasb calcula las preiones parciales y sus derivadas con
c     respectoa la temperatura y la presion.es una modificacion
c     de la rutina gas.
c
c     dp es la derivada de p respecto a t, ddp respecto a pe.
c=============================================================================

      subroutine gasb(theta,pe,p,dp,ddp)
      implicit real*4 (a-h, o-z)
      include 'PARAMETER'
      parameter (ncontr=28)
      dimension cmol(91),alfai(ncontr),chi1(ncontr),chi2(ncontr),
     *u0(ncontr),u1(ncontr),u2(ncontr)
      dimension du0(ncontr),du1(ncontr),du2(ncontr),dcmol(91)
      real*4 p(99),dp(99),ddp(99)
      real*4 abundance(ncontr),totalbundance,XA(ncontr),weightarr(ncontr)
      real*4 musuminv,mu,elec_gram,uma,Kboltz,R,ne,rho

      t=5040./theta

      if(t .gt. 1.5e4)then
         uma=1.660538921e-24
         Kboltz=1.38064852e-16
         R=8.3144598e7                   !cte de los gases 
         musuminv=0.
         elec_gram=0.
         do i=1,ncontr
            p(i)=0.                      ! p(h)/p(h')
            dp(i)=0.                         !dlf1/f1
            ddp(i)=0.                         !ddlf1/f1
            iii=i
            call neldatb(iii,0.,weighti,alfaii,chi1i,chi2i)
            weightarr(i)=weighti
            abundance(i)=weighti*alfaii
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
         ddp(89)=0.                          !ddlf5/f5
         p(90)=pe/p(85)                    ! pe/p(h') 
         dp(90)=0.                
         ddp(90)=0.                
         p(91)=pe/(Kboltz*T)               ! n(e)=pe/kt
         dp(91)=-1./t
         ddp(91)=1./pe
      else
      
      if(pe.le.0)then
         call error(KWARN,'gasb','Negative values of the electron pressure'
     &   //         ' have been found\n These are being changed to 1.e-10')
         pe=1.e-10
         g4=0.
         g5=0.
         dg4=0.
         ddg4=0.
         dg5=0.
         ddg5=0.
      else
         call molecb(theta,cmol,dcmol)
         do i=1,2
            call acota(cmol(i),-30.,30.)
         end do
         g4=pe*10.**(cmol(1))
         dg4=dcmol(1)*alog(10.)
         ddg4=1./pe
         g5=pe*10.**(cmol(2))
         dg5=dcmol(2)*alog(10.)
         ddg5=1./pe
      end if
      
c evaluating partition functions u0,u1,u2 and its derivatives
       do i=1,ncontr
         iii=i
         call neldatb(iii,0.,weight,alfai(i),chi1(i),chi2(i))
         call nelfctb(iii,t,u0(iii),u1(iii),u2(iii),du0(iii),du1(iii),
     *du2(iii))
       end do
      
      g2=saha(theta,chi1(1),u0(1),u1(1),pe)   ! p(h+)/p(h)
      dg2=dsaha(theta,chi1(1),du0(1),du1(1)) !derivo log(g2) con t
      ddg2=-1./pe 
                           !derivada log(g2) con pe
      p(92)=g2
      dp(92)=dg2
      ddp(92)=ddg2
      g3=saha(theta,0.754,1.,u0(1),pe)   ! p(h)/p(h-)

      call acota(g3,1.e-30,1.e30)
      g3=1.0/g3                              ! p(h-)/p(h) 
      dg3=-1.*dsaha(theta,0.754,0.,du0(1))
      ddg3=1./pe
      p(93)=g3
      dp(93)=dg3
      ddp(93)=ddg3
      g1=0.
      dlg1=0.       !las dl son derivadas no de log(g) sino de g
      ddlg1=0.

      do i=2,ncontr
         a=saha(theta,chi1(i),u0(i),u1(i),pe)
         da=dsaha(theta,chi1(i),du0(i),du1(i))
         dda=-1./pe
         b=saha(theta,chi2(i),u1(i),u2(i),pe)
         dlb=b*dsaha(theta,chi2(i),du1(i),du2(i))
         ddlb=-b/pe
         c=1.+a*(1.+b)
         call acotasig(c,1.e-20,1.e20)
         p(i)=alfai(i)/c   ! p/ph' for neutral he,li, ...
         dp(i)=-(a*da*(1.+b)+a*dlb)/c
         ddp(i)=-(a*dda*(1.+b)+a*ddlb)/c
         ss1=(1.+2.*b)
         call acotasig(ss1,1.e-20,1.e20)
         ss=p(i)*a*ss1
         dss=dp(i)+da+2.*dlb/ss1
         ddss=ddp(i)+dda+2.*ddlb/ss1
         g1=g1+ss
         dlg1=dlg1+ss*dss             !ojo estas no son derivadas del log
         ddlg1=ddlg1+ss*ddss
      end do
      a=1.+g2+g3
      dla=g2*dg2+g3*dg3
      ddla=g2*ddg2+g3*ddg3

      if(g5.lt.1.e-35)then
          call error(KSTOP,'gasb','The electronic pressure is too small.'
     &    //         ' Check the atmospheric models')
      end if

      call acotasig(g5,1.e-20,1.e20)
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

        call acotasig(a,1.e-15,1.e15)
        call acotasig(b,1.e-15,1.e15)
        call acotasig(c,1.e-15,1.e15)
        call acotasig(d,1.e-15,1.e15)
        call acotasig(e,1.e-15,1.e15)

      c1=c*b**2+a*(d*b-e*a)
        dlc1=dlc*b*b+(c*2.*b+a*d)*dlb+dla*(d*b-2.*e*a)+dld*a*b-dle*a*a
        ddlc1=ddlc*b*b+(c*2.*b+a*d)*ddlb+ddla*(d*b-2.*e*a)+ddld*a*b-ddle*a*a
      c2=2.*a*e-d*b+a*b*g1
        dlc2=dla*(2.*e+b*g1)+dlb*(a*g1-d)-dld*b+dle*2.*a+a*b*dlg1
        ddlc2=ddla*(2.*e+b*g1)+ddlb*(a*g1-d)-ddld*b+ddle*2.*a+a*b*ddlg1
      c3=-(e+b*g1)
        dlc3=-dle-dlb*g1-b*dlg1
        ddlc3=-ddle-ddlb*g1-b*ddlg1
           
      call acotasig(c1,1.e-15,1.e15)

      f1=0.5*c2/c1
      dc1=dlc1/c1  !!!!!!!!
      dc2=dlc2/c2
      dlf1=f1*(dc2-dc1)  

      ddc1=ddlc1/c1  !!!!!!!!
      ddc2=ddlc2/c2
      ddlf1=f1*(ddc2-ddc1)             

      dlf1=-dlf1+sign(1.,c1)*(2.*f1*dlf1-dlc3/c1+dlc1*c3/(c1*c1))
     *        /(2.*sqrt(f1**2-c3/c1))
      ddlf1=-ddlf1+sign(1.,c1)*(2.*f1*ddlf1-ddlc3/c1+ddlc1*c3/(c1*c1))
     *        /(2.*sqrt(f1**2-c3/c1))

      f1=-f1+sign(1.,c1)*sqrt(f1**2-c3/c1)

      f5=(1.-a*f1)/b
        if(abs(f5).lt.1.e-30)then
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
        if(abs(f4).lt.1.e-30)then
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

        if(abs(f1).lt.1.e-30)then
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

        call acotasig(f2,1.e-20,1.)
        call acotasig(fe,1.e-15,2.)
      
        phtot=pe/fe
        dlphtot=-phtot*dfe
        ddlphtot=(1.-ddlfe*phtot)/fe
        dphtot=-dfe
        ddphtot=1./pe-ddfe

        if(f5.le.1.e-4)then
            const6=g5/pe*f1**2
            dconst6=dg5+2.*divf1                !dlf1/f1
            ddconst6=ddg5-1./pe+2.*ddivf1        !ddlf1/f1
            const7=f2-f3+g1
            dlconst7=dlf2-dlf3+dlg1
            ddlconst7=ddlf2-ddlf3+ddlg1
            dconst7=dlconst7/const7
            ddconst7=ddlconst7/const7

            do i=1,5
                f5=phtot*const6
                df5=dphtot+dconst6
                ddf5=ddphtot+ddconst6
                f4=e*f5
                df4=de+df5
                ddf4=dde+ddf5
                fe=const7+f4
                call acotasig(fe,1.e-15,2.)
                dfe=(dlconst7+df4*f4)/fe
                ddfe=(ddlconst7+ddf4*f4)/fe
                phtot=pe/fe
                dphtot=-dfe
                ddphtot=1./pe-ddfe
            end do

            dlf5=df5*f5
            ddlf5=ddf5*f5
            dlf4=df4*f4
            ddlf4=ddf4*f4
            dlfe=dfe*fe
            ddlfe=ddfe*fe
            dlphtot=dphtot*phtot
            ddlphtot=ddphtot*phtot
         end if       
         pg=pe*(1.+(f1+f2+f3+f4+f5+0.1014)/fe)
        
        dlpg=pe*(dlf1+dlf2+dlf3+dlf4+dlf5)/fe-(pg-pe)*(dlfe/fe)
        ddlpg=pe*(ddlf1+ddlf2+ddlf3+ddlf4+ddlf5)/fe-(pg-pe)*(ddlfe/fe)
        ddlpg=ddlpg+pg/pe

        p(1)=f1   ! p(h)/p(h')
        dp(1)=divf1                !dlf1/f1
        ddp(1)=ddivf1                !ddlf1/f1

        call acotasig(pg,1.e-20,1.e20)
        p(84)=pg   ! gas pressure
        dp(84)=dlpg/pg
        ddp(84)=ddlpg/pg
        if(t .gt. 2.e4)then
           f2=1.0
           dlf200=0.
           ddlf200=0.
        end if   
        p(85)=phtot   ! p(h')
        dp(85)=dphtot                        !dlphtot/phtot
        ddp(85)=ddphtot                        !ddlphtot/phtot
        p(86)=f2   ! p(h+)/p(h')        
        dp(86)=dlf200
        ddp(86)=ddlf200
        p(87)=f3   ! p(h-)/p(h')
        if(abs(f3).lt.1.e-30)then
           divf3=0.
           ddivf3=0.
        else
           divf3=dlf3/f3
           ddivf3=ddlf3/f3
        end if

        dp(87)=divf3
        ddp(87)=ddivf3
        p(88)=f4   ! p(h2+)/p(h')
        dp(88)=df4                !dlf4/f4
        ddp(88)=ddf4                !ddlf4/f4
        p(89)=f5   ! p(h2)/p(h')
        dp(89)=df5                !dlf5/f5
        ddp(89)=ddf5                !ddlf5/f5
        p(90)=fe   ! pe/p(h')
        dp(90)=dfe                !dlfe/fe
        ddp(90)=ddfe                !ddlfe/fe
        p(91)=pe/(1.38054e-16*t) ! n(e)=pe/kt
        dp(91)=-1./t
        ddp(91)=1./pe

      end if

      return
      end


c=============================================================================
c     gasc calcula las presiones parciales
c=============================================================================

      subroutine gasc(t,pe,pg,pp)
      parameter (ncontr=28)
      dimension cmol(91),alfai(ncontr),chi1(ncontr),chi2(ncontr),
     *u0(ncontr),u1(ncontr),u2(ncontr)
        real du0,du1,du2,dcmol(91)
        real pp(*),p(28)

      theta=5040./t
      call molecb(theta,cmol,dcmol)
      g4=pe*10.**(cmol(1))
      g5=pe*10.**(cmol(2))
c evaluating partition functions u0,u1,u2 and its derivatives
      do i=1,ncontr
         iii=i
         call neldatb(iii,0.,weight,alfai(i),chi1(i),chi2(i))
         call nelfctb(iii,t,u0(iii),u1(iii),u2(iii),du0,du1,du2)
      end do 
      g2=saha(theta,chi1(1),u0(1),u1(1),pe)   ! p(h+)/p(h)
      g3=1./saha(theta,0.754,1.,u0(1),pe)     ! p(h-)/p(h)
      pp(10)=g3
      g1=0.

       do i=2,ncontr
          a=saha(theta,chi1(i),u0(i),u1(i),pe)
          b=saha(theta,chi2(i),u1(i),u2(i),pe)
          if(a .gt. 1.e30)then
             p(i)=0.
             if(b .lt. 1.e30)g1=g1+alfai(i)*((1.+2.*b)/(1.+b))
             if(b .ge. 1.e30)g1=g1+2.*alfai(i)
          else
             if(b .ge. 1.e30)g1=g1+2.*alfai(i)
             if(b .lt. 1.e30)then
                if(a*b .gt. 1.e35)then
                   p(i)=0.
                   g1=g1+alfai(i)*((1./2.+2.*b)/(1.+b))
                else
                   c=1.+a*(1.+b)
                   p(i)=alfai(i)/c   ! p/ph' for neutral he,li, ...
                   g1=g1+p(i)*a*(1.+2.*b)
                endif   
             endif   
          end if  
       end do

        pp(2)=p(2)
        pp(3)=p(6)
        pp(4)=p(11)
        pp(5)=p(12)
 
        a=1.+g2+g3
        b=2.*(1.+g2/g5*g4)
        c=g5
        d=g2-g3
        e=g2/g5*g4
        c1=c*b**2+a*(d*b-e*a)
        if(c1 .gt. 1.e33)c1=1.e33 
        c2=2.*a*e-d*b+a*b*g1        
        c3=-(e+b*g1)
        f1=0.5*c2/c1
        f1=-f1+sign(1.,c1)*sqrt(f1**2-c3/c1)
        f5=(1.-a*f1)/b
        f4=e*f5
        f3=g3*f1
        f2=g2*f1
        fe=f2-f3+f4+g1
        phtot=pe/fe

        if(f5.le.1.e-4)then
           const6=g5/pe*f1**2
           const7=f2-f3+g1
           do i=1,5
             f5=phtot*const6
             f4=e*f5
             fe=const7+f4
             phtot=pe/fe
           end do
        end if

        pg=pe*(1.+(f1+f2+f3+f4+f5+0.1014)/fe)
        pp(1)=f1   ! p(h)/p(h')
        pp(6)=f2   ! p(h+)/p(h')        
        pp(7)=f5   ! p(h2)/p(h')
        pp(8)=fe   ! pe/p(h')
        pp(9)=pe/(1.38054e-16*t) ! n(e)=pe/kt

      return
      end


c=============================================================================
c     gasccota es igual que gasc con control de fe < 0 (if de linea 524)
c=============================================================================

      subroutine gasccota(t,pe,pg,pp)
      parameter (ncontr=28)
      dimension cmol(91),alfai(ncontr),chi1(ncontr),chi2(ncontr),
     *u0(ncontr),u1(ncontr),u2(ncontr)
        real du0,du1,du2,dcmol(91)
        real pp(*),p(28)

      theta=5040./t
      call molecb(theta,cmol,dcmol)
      g4=pe*10.**(cmol(1))
      g5=pe*10.**(cmol(2))
c evaluating partition functions u0,u1,u2 and its derivatives
      do i=1,ncontr
         iii=i
         call neldatb(iii,0.,weight,alfai(i),chi1(i),chi2(i))
         call nelfctb(iii,t,u0(iii),u1(iii),u2(iii),du0,du1,du2)
      end do 
      g2=saha(theta,chi1(1),u0(1),u1(1),pe)   ! p(h+)/p(h)
      g3=1./saha(theta,0.754,1.,u0(1),pe)     ! p(h-)/p(h)
      pp(10)=g3
      g1=0.

       do i=2,ncontr
          a=saha(theta,chi1(i),u0(i),u1(i),pe)
          b=saha(theta,chi2(i),u1(i),u2(i),pe)
          if(a .gt. 1.e30)then
             p(i)=0.
             if(b .lt. 1.e30)g1=g1+alfai(i)*((1.+2.*b)/(1.+b))
             if(b .ge. 1.e30)g1=g1+2.*alfai(i)
          else
             if(b .ge. 1.e30)g1=g1+2.*alfai(i)
             if(b .lt. 1.e30)then
                if(a*b .gt. 1.e35)then
                   p(i)=0.
                   g1=g1+alfai(i)*((1./2.+2.*b)/(1.+b))
                else
                   c=1.+a*(1.+b)
                   p(i)=alfai(i)/c   ! p/ph' for neutral he,li, ...
                   g1=g1+p(i)*a*(1.+2.*b)
                endif   
             endif   
          end if  
       end do

        pp(2)=p(2)
        pp(3)=p(6)
        pp(4)=p(11)
        pp(5)=p(12)
 
        a=1.+g2+g3
        b=2.*(1.+g2/g5*g4)
        c=g5
        d=g2-g3
        e=g2/g5*g4
        c1=c*b**2+a*(d*b-e*a)
        if(c1 .gt. 1.e33)c1=1.e33 
        c2=2.*a*e-d*b+a*b*g1        
        c3=-(e+b*g1)
        f1=0.5*c2/c1
        f1=-f1+sign(1.,c1)*sqrt(f1**2-c3/c1)
        f5=(1.-a*f1)/b
        f4=e*f5
        f3=g3*f1
        f2=g2*f1
        fe=f2-f3+f4+g1
        if(fe .le. 0)then !this is because f3 (H- population is too high --> 0)
           f3=1.e-20
           fe=f2+f4+g1
        end if   
        phtot=pe/fe

        if(f5.le.1.e-4)then
           const6=g5/pe*f1**2
           const7=f2-f3+g1
           do i=1,5
             f5=phtot*const6
             f4=e*f5
             fe=const7+f4
             phtot=pe/fe
           end do
        end if

        pg=pe*(1.+(f1+f2+f3+f4+f5+0.1014)/fe)
        pp(1)=f1   ! p(h)/p(h')
        pp(6)=f2   ! p(h+)/p(h')        
        pp(7)=f5   ! p(h2)/p(h')
        pp(8)=fe   ! pe/p(h')
        pp(9)=pe/(1.38054e-16*t) ! n(e)=pe/kt

      return
      end
