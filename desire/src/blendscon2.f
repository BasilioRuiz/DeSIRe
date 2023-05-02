c_______________________________________________________________
c blendscon2 (rutina del programa inversion)
c igual que la rutina fstokes0 pero teniendo en cuenta blends
c calcula las funciones respuesta de I
c para un angulo heliocentrico de coseno mu. se supone lte,
c el maximo numero de puntos en tau es 64
c Basilio Ruiz (23 julio 93)
c ______________________________________________________________
c
c ist=1 (i); =2 (q); =3 (u); =4 (v)

         subroutine blendscon2(atmos,inten,rt,rp,rv,rm,mnodos,beta1,beta2,atmoserr)
         implicit real*4 (a-h,o-z)
         include 'PARAMETER'

c para las funciones respuesta
         real*4 rt(*),rp(*),rv(*),rm(*)
         real*4 grt(kt),grp(kt),grv(kt),grm(kt)
         real*4 x(kt)

c para los perfiles 
         real*4 inten(*),icont,continuoharr(kld)  !con_i(kl)
         real*4 wcon_i1(kl),wcon_i2(kl)

c para la atmosfera
         real*4 atmos(*),atmoserr(*)
         real*4 grtmax(kt),grpmax(kt),grmmax(kt),grvmax(kt)
         real*4 tau(kt),t(kt),pe(kt),vtur(kt),vz(kt)
         real*4 taue(kt),vof(kt),logpe(kt)
         real*4 continuoh    
         integer mnodos(*)

c para la matriz de absorcion y sus derivadas
         real*4 dab(kt),tk(kt),pk(kt),vk(kt),mk(kt)
         real*4 dabtot(kt,kld)
         real*4 tktot(kt,kld),pktot(kt,kld),mktot(kt,kld),vktot(kt,kld)
         real*4 www,dyt(kt),dyp(kt),alpha(kt)
 
c para la malla
         real*4 dlongd(kld),dlamda0(kl)
         real*4 wlengt,wlengt1,lambda,dlamda,wvac,wc,c
         integer nlin(kl),npas(kl),ist(4),nble(kl)

c para los parametros atomicos y coeficientes de absorcion
         integer atomic_number,atom_arr(kl),istage_arr(kl)
         real*4 alfa_arr(kl),sigma_arr(kl),wave_arr(kl)
         real*8 wavedbl_arr(kl),wlengtdbl,wlengt1dbl
         real*4 wlengt1_arr(kl),gf_arr(kl),energy_arr(kl)
         real*4 chi10_arr(kl),chi20_arr(kl),abu_arr(kl)
         integer*4 nel_arr(kl)         
         real*4 lambda_arr(kl),croot_arr(kl),wvac_arr(kl),weight_arr(kl)
         real loggf,nair,mvdop,meta0,ma
         real*4 y(kt),table(kt,16)
         real*4 ck5(kt),dk5(kt),ddk5(kt),zeff
         real*4 ckappa,ckappa5,dkappa,dkappa5,ddkappa,ddkappa5
         real*4 alfa,sigma 
         integer nlow,nup

c para el damping        
         real*4 chydro_arr(kl),xmu1_arr(kl),xmu2_arr(kl),xmu3_arr(kl)
         real*4 vv_arr(kl),beta_arr(kl)         
         
c para la emisividad y la funcion fuente
         real*8 beta1(kl,kt),beta2(kl,kt),blow,bup !,bratio
         real*8 wwwdbl,a_plck,c_plck,exalfa_plck,bpdob
         real*4 emisividad(kt,kld),emisividad_dt(kt,kld)
         real*4 bp(kt),bt(kt),bp0(kt),bt0(kt),dbp(kt)
         real*4 source(kt),source_dt(kt)
         
c para las presiones parciales
         integer ivar(10)
         real pg(99),dpg(99),ddpg(99),pi(10),dpi(10),ddpi(10)
         real pt(kt,10),dpt(kt,10),ddpt(kt,10),abu,piis
         real*4 pgas(kt),dpgas(kt),ddpgas(kt)      

c para la inclusion de RP en RT
         real*4 ax(kt),bx(kt),cx(kt),dx(kt),d1x(kt),d2x(kt),fx(kt)
         real*4 px(kt),qx(kt),rx(kt),sx(kt),tx(kt),wx(kt,kt)
         real*4 kac,kat,kap

c para hermite_c
         real*4 deltae(kt),deltai(kt),delt2i(kt)

c para barklem
         real*4 bol,pir,v0,melectron,mhidrogeno,xmasaproton,avo,uma
         real*4 gas,coc2,coc3
         real*8 saha_db,dsaha_db
         real*8 u12_db,u23_db,u33_db
         real*8 du12_db,du23_db,du33_db
         real*8 ddu12_db,ddu23_db,ddu33_db

c lugares comunes de memoria
         common/responde2/ist,ntau,ntl,nlin,npas,nble
         common/ldeo/dlongd,dlamda0
         common/brklm/ntotal_lines,atom_arr,istage_arr,alfa_arr,sigma_arr,wave_arr
         common/wavearrdble/wavedbl_arr
         common/datosatom/wvac_arr,weight_arr,croot_arr,wlengt1_arr,lambda_arr
         common/datosdamping/chydro_arr,xmu1_arr,xmu2_arr,xmu3_arr,vv_arr,beta_arr
         common/loggfarr/gf_arr,energy_arr
         common/pot_ion/chi10_arr,chi20_arr,nel_arr,abu_arr
         common/piis/piis
         common/yder/y,dyt,dyp,alpha  !coef. de abs. del continuo y su der. t,p
         common/segunda/tau,taue,deltae,deltai,delt2i
         common/offset/voffset  !para respuestas
         common/iautomatico/iautomatico
c        common/mu/cth                              !esto esta puesto a 1 en sir
         common/anguloheliocent/xmu 
c        common/continuos/con_i
         common/continuosarr/continuoharr
         
         data iprimera/0/

c nble es el numero de componentes de cada linea
         c=2.99792458e+10 	!vel. de la luz en cm/seg
         piis=1./sqrt(3.1415926)
         g=xmu*2.7414e+4	 !gravedad cm/s^2 en fotosfera solar   
         avog=6.023e23
         bol=1.3807e-16         !erg/s
         pir=3.1415926
         v0=1e6                 !cm/s
         melectron=9.1094e-24
         mhidrogeno=1.67442e-24
         xmasaproton=1.6526e-24
         avo=6.023e23
         bohr=0.0529177249e-7   !cm
         uma=1.660540e-24
         gas=8.31451e7         !constante de los gases en cgs
         coc3=1.212121          ! cociente polarizabilidad del H2 con HI
         coc2=.3181818          ! idem para el He

         ntotal=0
         do i=1,ntl
           do j=1,npas(i)
		      ntotal=ntotal+1
           end do
	     end do

	     if(iprimera.eq.0)then
	        do i=1,ntau
	           tau(i)=atmos(i)
	           taue(i)=10.**(tau(i))
	        end do
	        do i=1,ntau-1
              deltae(i)=taue(i)-taue(i+1)
            end do 
	        do i=2,ntau
              deltai(i)=(tau(i)-tau(i-1))/2.0
	          delt2i(i)=deltai(i)*deltai(i)/3.0
            end do
	        paso=tau(1)-tau(2)
        end if

        iprimera=iprimera+1
        
        x(1)=g*(tau(1)-tau(2))*2.3025851         
        do i=1,ntau-1                             
           x(i+1)=g*(tau(i)-tau(i+1))*2.3025851   
        end do                                    
        
c leemos la atmosfera
        do i=1,ntau
           grtmax(i)=0.  !inicializamos para calcular el maximo en lamda para cada tau
           grpmax(i)=0.    
           grmmax(i)=0.
           grvmax(i)=0.
           t(i)=atmos(i+ntau)
           pe(i)=atmos(i+2*ntau)
           logpe(i)=alog(pe(i))
           vtur(i)=atmos(i+3*ntau)
           vz(i) =atmos(i+5*ntau)       !velocidad lines de vision
           vof(i)=vz(i)-voffset
        end do

c calculamos el peso molecular medio pmu
	    pmusum=0.0
        asum=0.0
        do i=1,92
  	       ii=i
	       call neldatb(ii,0.,wgt,abu,ei1,ei2)
	       pmusum=pmusum+wgt*abu
           asum=asum+abu
        end do

c calculo las presiones parciales (pg) y las derivadas de sus logaritmos
c con respecto a la t (dpg) y, con pe (ddpg)
        ivar(1)=1
        ivar(2)=2
        ivar(3)=6
        ivar(4)=11
        ivar(5)=12
        ivar(6)=86
        ivar(7)=89
        ivar(8)=90
        ivar(9)=91
        ivar(10)=93
	
        do i=1,ntau
           ps=pe(i)
           ts=t(i)
           theta=5040./ts
           call gasb_db(theta,ps,pg,dpg,ddpg)
	    
           do j=1,10
              k=ivar(j)
              pt(i,j)=pg(k)
              dpt(i,j)=dpg(k)
              ddpt(i,j)=ddpg(k)
              pi(j)=pg(k)
              dpi(j)=dpg(k)
              ddpi(j)=ddpg(k)
           end do

           pgas(i)=pg(84)      !presion gaseosa
           dpgas(i)=dpg(84)    !der. log(presion gaseosa) respectp a T
           ddpgas(i)=ddpg(84)  !der. log(presion gaseosa) respectp a pe
           call kappach(5.e-5,ts,ps,pi,dpi,ddpi,ck5(i),dk5(i),ddk5(i))
           cc=avog/pmusum
           kac=ck5(i)*cc
           kat=dk5(i)*cc
           kap=ddk5(i)*cc
           tauk=taue(i)/2./kac/kac

           fx(i)=dpg(84)*pg(84)
           cx(i)=1./(ddpg(84)*pg(84))
           bx(i)=kap/(ddpg(84)*pg(84))
           ax(i)=kat-kap*dpg(84)/ddpg(84)
           d1x(i)=x(i)*tauk
           d2x(i)=x(i)*tauk
           dx(i)=1.+d1x(i)*bx(i)
        end do

        rx(ntau)=0. 
        sx(ntau)=0. 
        tx(ntau)=0. 
        do i=1,ntau-1
           rx(i)=(1.- d2x(i+1)*bx(i+1))/dx(i)
           sx(i)=(d2x(i+1)/dx(i))*ax(i+1)
           tx(i)=(d1x(i)/dx(i))*ax(i)
           px(i)=-cx(i)*(tx(i)+fx(i))
        end do

        do i=2,ntau
           qx(i)=-cx(i-1)*(sx(i-1)+tx(i)*rx(i-1))
        end do
        qx(1)=0.

        do i=1,ntau-1
           do j=1,i+1
              wx(i,j)=0.
           end do
           do j=ntau-1,ntau
              wx(i,j)=0.
           end do
           do j=i+2,ntau-2
              r1=1.
              do k=i,j-2
                 r1=r1*rx(k)
              end do
              wxi=-(sx(j-1)+tx(j)*rx(j-1))*r1 
              wx(i,j)=cx(i)*wxi
              
           end do
        end do

c datos de la linea
        ikk0=0
        ikk1=0
        ixx=0
        do 999 iln=1,ntl	!iln=numero de la linea,ntl num.total
           do i=1,npas(iln)
               do k=1,ntau
                  dabtot(k,i)=0.
                  tktot(k,i)=0.
                  pktot(k,i)=0.
                  vktot(k,i)=0.
                  mktot(k,i)=0.
               end do
           end do
           do ible=1,nble(iln)
              ixx=ixx+1
              nxx=nlin(ixx) 

              energy=energy_arr(ixx)
              gf=gf_arr(ixx)
              istage=istage_arr(ixx)
              alfa=alfa_arr(ixx)
              sigma=sigma_arr(ixx)
              wlengt=wave_arr(ixx)
              wlengtdbl=wavedbl_arr(ixx)
              
              if(ible.eq.1)wlengt1=wlengt
              if(ible.eq.1)wlengt1dbl=wlengtdbl
              dlamda0(iln)=wlengt1
c             continuoh=con_i(iln) 

c parametros atomicos
             chi10=chi10_arr(ixx)
             chi20=chi20_arr(ixx)
             nel=nel_arr(ixx)
             abu=abu_arr(ixx)

             wvac=wvac_arr(ixx)     
             weight=weight_arr(ixx)   
             croot= croot_arr(ixx) 
             wc=wvac/c               !inverso de la frecuencia (segundos)

c calculo de terminos constantes.
             weinv=1./weight
             dlo=4.6686e-5*wvac*wvac !separac.de l.d.o=dlo*h*(m1*g1-m2*g2)
             crad=.22233/wvac        !(l.d.o.*coef.clasico ensanch. natural)
             eta00=1.49736e-2*gf*abu*wvac !para el coef. de absor.de la linea

c chydro (gamma6)coef.van der waals ensanch.colisional
             chydro=chydro_arr(ixx)
             xmu1=xmu1_arr(ixx)
             xmu2=xmu2_arr(ixx)
             xmu3=xmu3_arr(ixx)  
             vv=vv_arr(ixx)
             beta=beta_arr(ixx)        

c estratificacion en tau 5000
c genero una tabla lineal en logaritmo de tau a l.d.o.=5000.angs.
c desde tauini a taufin, con ntau valores.
             lambda=lambda_arr(ixx)           !igual que  wwwdbl pero en real*4 (para kappach)
             wwwdbl=wlengtdbl*1.d-8           !Atencion cojo la ldo central de cada linea
             a_plck=1.43880/wwwdbl            !for Planck function
             c_plck=1.1910627d-5/(wwwdbl**5)  !for Planck function

             do 71 i=1,ntau
                ps=pe(i)
                ts=t(i)
                theta=5040./ts
                do j=1,10
                   k=ivar(j)
                   pi(j)=pt(i,j)
                   dpi(j)=dpt(i,j)
                   ddpi(j)=ddpt(i,j)
                   pg(k)=pi(j)
                   dpg(k)=dpi(j)
                   ddpg(k)=ddpi(j)
                end do

c calculo el coeficiente de absorcion del continuo por cm**3
c (ckappa) y sus derivadas con respecto a t y pe:dkappa,
c ddkappa y lo divido por dicho coef. evaluado a 5000.a
c calculo el ckappa para lambda=(5000 a=5.e-5 cm),(ckappa5)
               call kappach(lambda,ts,ps,pi,dpi,ddpi,ckappa,
     &                   dkappa,ddkappa)
     
               ckappa5=ck5(i)
               dkappa5=dk5(i)
               ddkappa5=ddk5(i)

               ckappa=ckappa/ckappa5
               dkappa=(dkappa-ckappa*dkappa5)/ckappa5
               ddkappa=(ddkappa-ckappa*ddkappa5)/ckappa5
               dkappa5=dkappa5/ckappa5
               ddkappa5=ddkappa5/ckappa5

c calculo vdop=(delta doppler*c/l.d.o.) y su derivada
               vdop=sqrt(croot*ts+vtur(i)**2) !cm/seg.
               vdop2=sqrt((2.*gas*ts)/weight+vtur(i)**2)
               dvdop=croot/vdop/2.
               dvdop2=gas/(vdop2*weight)
               mvdop=vtur(i)/vdop

c calculo el coeficiente absorcion linea en cada tau y sus derivadas	      
c calculo eta0=(coef.absorcion linea/coef. absorcion continuo)
c necesito calcular la fraccion de atomos del elemento en el
c el nivel de la transicion respecto al numero total de atomos
c del elemento en cualquier estado de ionizacion.asi, tengo que
c utilizar las ecuacione de saha y boltzmann.
c llamare u12 al cociente de poblaciones entre el estado de ioniza
c cion 1 (neutro) y 2.igualmente u23 entre los iones 2 y 3
c necesito llamar a nelfctb (definida en atmdatb con un 'entry')
c para calcular las funciones de particion y sus derivadas a cada
c temperatura
               call nelfctb(nel,ts,u1,u2,u3,du1,du2,du3)

c corrijo los potenciales de ionizacion chi10 y chi20 con
c un termino proporcional a la raiz cubica de la densidad de e-
               rcu=(pg(91))**(1./3.)
               chi1=chi10-6.96e-7*rcu
               chi2=chi20-1.1048e-6*rcu

               u12_db=saha_db(theta,chi1,u1,u2,ps)        
               du12_db=u12_db*dsaha_db(theta,chi1,du1,du2) !derivada de u12 con t
               ddu12_db=-1.*u12_db/ps        !    "    "   "  con pe
               u23_db=saha_db(theta,chi2,u2,u3,ps)      !n3/n2
               du23_db=u23_db*dsaha_db(theta,chi2,du2,du3) !derivada de u23 con t
               ddu23_db=-1.*u23_db/ps                !    "    "   "  con pe
               u33_db=1.+u12_db*(1.+u23_db)                !(n1+n2+n3)/n1
               du33_db=du12_db*(1.+u23_db)+u12_db*du23_db
               ddu33_db=ddu12_db*(1.+u23_db)+u12_db*ddu23_db
               u12=real(u12_db)
               du12=real(du12_db)
               ddu12=real(ddu12_db)
               u23=real(u23_db)
               du23=real(du23_db)
               ddu23=real(ddu23_db)
               u33=real(u33_db)
               du33=real(du33_db)
               ddu33=real(ddu33_db)
           
               eta0=eta00*1.e1**(-theta*energy)/(u1*vdop*u33*ckappa5)
               deta0=eta0*(alog(10.)*theta/ts*energy-real(du1_db)-dvdop/vdop-real(du33_db/u33_db))
               deta0=deta0-eta0*dkappa5
               ddeta0=eta0*real(-ddu33_db/u33_db)
               ddeta0=ddeta0-eta0*ddkappa5	  
 
               if(istage.eq.2)then
                  eta0=eta0*u1/u2*u12
                  deta0=deta0*u1/u2*u12+eta0*(du1-du2+du12/u12)
                  ddeta0=ddeta0*u1/u2*u12+eta0*(ddu12/u12)
               end if

c introduzco la correccion por emision estimulada
              blow=beta1(ixx,i)
              blows=real(blow)
              bup=beta2(ixx,i)
              c3e=1.4388/(ts*wvac)
              corre=1.-bup*exp(-c3e)/blow
              dcorre=(corre-1.)*c3e/ts
 
              deta0=deta0*corre+eta0*dcorre
              ddeta0=ddeta0*corre
              eta0=eta0*corre
              meta0=-eta0*mvdop/vdop

c calculo el damping 'a' mediante Unsold
c despreciaremos la dependencia del damping con t o pe en el calcu
c lo de las funciones respuesta.
c pg(1)=p(h)/p(h'),pg(90)=pe/p(h'),pg(91)=densidad de e-
c pg(2)=p(he)/p(h'),pg(89)=p(h2)/p(h')
c   a=(chydro*(pg(1)/pg(90)*pg(91))*t(i)**0.3*((.992093+weinv)**.3
c  &  +.6325*pg(2)/pg(1)*(.2498376+weinv)**.3+.48485*pg(89)/pg(1)*
c  &  (.4960465+weinv)**.3)+crad)/(12.5663706*vdop)
               if(abs(alfa).lt.1.e-25.or.abs(sigma).lt.1.e-25)then
                  if(pg(1).gt.1.e-20)then !EVITAMOS DAMPING NULO POR PG(1)=0
	                 aj=chydro*(pg(1)/pg(90)*pg(91))*ts**0.3
     	             ai=(.992093+weinv)**.3+.6325*pg(2)/pg(1)*(.2498376+weinv)**.3+
     &	                 .48485*pg(89)/pg(1)*(.4960465+weinv)**.3
	                 a=(aj*ai+crad)/(12.5663706*vdop)
                     ma=-a*mvdop/vdop
	                 daj=aj*(dpg(1)-dpg(90)+dpg(91)+.3/ts) !derivada de aj cont
	                 ddaj=aj*(ddpg(1)-ddpg(90)+ddpg(91))  !derivada de aj con p
	                 dai=.6325*pg(2)/pg(1)*(.2498376+weinv)**.3*(dpg(2)-dpg(1))+
     &	                 .48485*pg(89)/pg(1)*(.4960465+weinv)**.3*(dpg(89)-dpg(1))
	                 ddai=.6325*pg(2)/pg(1)*(.2498376+weinv)**.3*(ddpg(2)-ddpg(1))+
     &	                 .48485*pg(89)/pg(1)*(.4960465+weinv)**.3*(ddpg(89)-ddpg(1))

	                 da=(ai*daj+aj*dai)/(12.5663706*vdop)-a*dvdop/vdop
	                 dda=(ai*ddaj+aj*ddai)/(12.5663706*vdop)
                  else  !CORRECCION PARA EVITAR DAMPING NULO POR PG(1)=0
                     aj=chydro*(1./pg(90)*pg(91))*ts**0.3
                     ai=.6325*pg(2)*(.2498376+weinv)**.3+
     &	                .48485*pg(89)*(.4960465+weinv)**.3
	                 a=(aj*ai+crad)/(12.5663706*vdop)
	                 ma=-a*mvdop/vdop
	
	                 daj=aj*(-dpg(90)+dpg(91)+.3/ts) !derivada de aj cont
	                 ddaj=aj*(-ddpg(90)+ddpg(91))  !derivada de aj con p
                     dai=.6325*pg(2)*(.2498376+weinv)**.3*dpg(2)+
     &	                 .48485*pg(89)*(.4960465+weinv)**.3*dpg(89)
	                 ddai=.6325*pg(2)*(.2498376+weinv)**.3*ddpg(2)+
     &	                .48485*pg(89)*(.4960465+weinv)**.3*ddpg(89)

	                 da=(ai*daj+aj*dai)/(12.5663706*vdop)-a*dvdop/vdop
	                 dda=(ai*ddaj+aj*ddai)/(12.5663706*vdop)
                  end if !CORRECCION PARA EVITAR DAMPING NULO POR PG(1)=0
               else

c       calculo del damping 'a' mediante BARKLEM.       
c       d son las derivadas totales respecto a la temperatura.
c       dd son las derivadas totales respecto a la presion.
                  dam=beta*(pg(91)/pg(90))*(pg(1)*xmu1**(-vv)+coc2*pg(2)*(xmu2**(-vv))
     &                +coc3*pg(89)*xmu3**(-vv))*ts**vv
        
                  a=(1./(4.*pir))*(crad/vdop+dam/vdop2) 
                  ma=-a*mvdop/vdop         

                  ddam=dam*(dpg(91)-dpg(90))+beta*(pg(91)/pg(90))*(xmu1**(-vv)*pg(1)*
     &                 dpg(1)+coc2*xmu2**(-vv)*pg(2)*dpg(2)+coc3*pg(89)*dpg(89)
     &                 *xmu3**(-vv))*ts**(vv)+vv*(dam/ts)

                  dddam=dam*(ddpg(91)-ddpg(90))+beta*(pg(91)/pg(90))*(xmu1**(-vv)*
     &                  pg(1)*ddpg(1)+coc2*xmu2**(-vv)*pg(2)*ddpg(2)+coc3*pg(89)*ddpg(89)*
     &                  xmu3**(-vv))

                  da=(crad/(4.*pir))*(-dvdop/(vdop**2.))+(1./(4.*pir))*(ddam*vdop2-
     &               dam*dvdop2)/(vdop2**2.)
                  dda=(1./(4.*pir))*(dddam/vdop2)
               endif   

c calculo la funcion de planck en lambda y su derivada con t
c dtplanck" es la derivada de la f. de planck con temperatura.
               exalfa_plck=dexp(a_plck/ts)
               if(ible .eq. 1)then 
                  bpdob=c_plck/(exalfa_plck-1.d0)
                  bp0(i)=real(bpdob)
                  bt0(i)=real(bpdob*bpdob*a_plck*exalfa_plck/(dble(ts)*dble(ts)*c_plck))
               end if
            
               bpdob=c_plck*bup/(blow*exalfa_plck-bup)
               bp(i)=real(bpdob)
               bt(i)=real(bpdob*bpdob*a_plck*exalfa_plck*blow/(dble(ts)*dble(ts)*c_plck*bup))

               eta0=eta0*blows
               deta0=deta0*blows
               ddeta0=ddeta0*blows
               meta0=meta0*blows

               y(i)=ckappa              !ckappa/ckappa5
               dyt(i)=dkappa
               dyp(i)=ddkappa
               table(i,1)=dkappa        !derivada con t
               table(i,2)=ddkappa       !    "     "  pe
               table(i,3)=eta0    !kap. linea/kcont.5000
               table(i,4)=deta0   !derivada con t
               table(i,5)=ddeta0  !derivada con pe
               table(i,6)=wc*vdop !despl. doppler en l.d.o(cm)
               table(i,7)=vz(i)/vdop !velocidad eje z (u.doppler)
               table(i,8)=dkappa5       !derivada con t de kappacon5
               table(i,9)=ddkappa5      !  "       "  pe      "
               table(i,10)=a      !damping
               table(i,11)=da   !derivada del damping con t
               table(i,12)=dda    !derivada del damping con p
               table(i,13)=dvdop/vdop   !derivada de log(vdop) con t
               table(i,14)=mvdop/vdop   !   "      "    "      con mic
               table(i,15)=meta0        !derivada de eta0 con la micro
               table(i,16)=ma           !derivada de a    con la micro
71	         continue            !do en log(tau)

c        call derivacuad(bp,dbp,ntau,tau)
             amaxim=0.
             do i=1,ntau
                if(table(i,10).ge.amaxim)amaxim=table(i,10)
             end do

            if(amaxim.gt.3.)then
               do i=1,ntau
                  table(i,10)=(3.*table(i,10))/amaxim
                  table(i,11)=(3.*table(i,11))/amaxim
                  table(i,12)=(3.*table(i,12))/amaxim
                  table(i,16)=(3.*table(i,16))/amaxim
               end do
            end if
    
c inicializamos la emisividad a todas las frecuencias con la emisividad del continuo
        if(ible .eq. 1)then !evaluo la emisividad del continuo
           do j=1,ntau
              emj=bp0(j)*y(j)
              demj=bt0(j)*y(j)+bp0(j)*dyt(j)
              do i=1,npas(iln)
                 emisividad(j,i)=emj
                 emisividad_dt(j,i)=demj
              end do
           end do
        end if 

c muestreo en lambda , calculo las l.d.o. de cada punto
        do 10 i=1,npas(iln)
	       ikk=ikk0+i
	       if(ible.eq.nble(iln).and.i.eq.npas(iln))ikk0=ikk0+npas(iln)
	       dlamda=(dlongd(ikk)+(wlengt1-wlengt)*1.e3)*1.e-11  !en cm.

c calculo etar,etal,etap y sus derivadas
	   do 101 j=1,ntau     !do en tau
	      a=table(j,10)	!damping
	      dldop=table(j,6) !desplazamiento doppler en cm
	      dvcam=table(j,7) !campo de velocidades unidades doppler
	      v=dlamda/dldop-dvcam !(l.d.o.+c.vel.) en unidades doppler
	      if(ible.eq.1)then
	      	t0=y(j)		!coef. absor. continuo/c.a.c 5000
	      	t1=table(j,1)     !derivada de t0 con t
	      	t2=table(j,2)     !  "         "      pe
	      else
	        t0=0.
	        t1=0.
	        t2=0.
	      end if
	      t3=table(j,3)     !coef. absor. linea/c.a.c 5000
	      t4=table(j,4)     !derivada de t3 con t
	      t5=table(j,5)     !  "         "   "  pe
	      t15=table(j,15)   !  "         "   "  micro
	      t11=table(j,11)   !derivada del damping con t
	      t12=table(j,12)   !  "       "     "     "  pe
	      t16=table(j,16)   !  "       "     "     "  micro
	      t13=table(j,13)   !derivada de log(vdop) con t
	      t14=table(j,14)   !  "      "      "      "  micro
	
	      call mvoigtc(a,v,t13,t14,wc,dldop,etar,vetar,
     &	                   ettar,ettvr,ettmr)

              dabtot(j,i)= dabtot(j,i)+t0 +t3*etar
              tktot(j,i) = tktot(j,i) +t1 +t4*etar + t3*(ettar*t11+ettvr)
              mktot(j,i) = mktot(j,i)     +t15*etar+ t3*(ettar*t16+ettmr)
              pktot(j,i) = pktot(j,i) +t2 +t5*etar + t3* ettar*t12 
              vktot(j,i) = vktot(j,i)              + t3* vetar
                           
c calculamos la emisividad
              emisividad(j,i)=emisividad(j,i)+bp(j)*t3*etar
              emisividad_dt(j,i)=emisividad_dt(j,i)+bt(j)*t3*etar+bp(j)*(t4*etar+t3*(ettar*t11+ettvr))
          
101	     continue !fin do en tau(estamos aun dentro del do en lamda)
10	  continue	!fin del do en lambda(pasos)
	end do   !fin del do en blends
c	................................................................

         do 9 i=1,npas(iln)
	        do k=1,ntau
	           dabk=dabtot(k,i)
	           tkk=tktot(k,i)
	           pk(k)=pktot(k,i)
               vk(k)=vktot(k,i)
	           mk(k)=mktot(k,i)
               sourcek=emisividad(k,i)/dabk 
               source_dt(k)=(emisividad_dt(k,i)-sourcek*tkk)/dabk
               dab(k)=dabk
               tk(k)=tkk
               source(k)=sourcek
	       end do

           call derivacuad(source,dbp,ntau,tau)
          
	       call hermite_c(source,dbp,dab,ntau,icont,kt,source_dt,tk,pk,vk,mk,
     &             grt,grp,grv,grm)	
     
c introducimos el paso en tau para pasar la integral sobre las f. resp.
c a sumatorio y normalizmos por el continuo
c las perturbaciones son relativas por tanto las funciones respuesta equi
c valentes en los nodos son multiplicadas por el valor del parametro
c en el nodo (en 'nodos')
c !ojo las perturbaciones son relativas a los parametros en z no a linea vision

        ikk1=ikk1+1
        continuoh=continuoharr(ikk1)

        if(mnodos(1).ne.0)then
           if(mnodos(2).eq.0)then
	          do kk=1,ntau  !introduccion de grp en grt
                 suma=0.
                 do kj=1,kk-2
                    suma=suma+wx(kj,kk)*grp(kj)
                 end do
                 correc=grp(kk)*px(kk)+suma
	             if(kk.gt.1)correc=correc+grp(kk-1)*qx(kk)
                 grt(kk)=grt(kk)+correc
              end do
           end if
           call rnorma(ntau,continuoh,grt)
           if(iautomatico.ne.1.or.mnodos(1).eq.1)then
              call nodos(grt,ntau,tau,t,mnodos(1),grtmax) !al final de blends2
              do ie=1,ntau
                 atmoserr(ntau+ie)=grtmax(ie)
              end do	
           end if    
        end if
        if(mnodos(2).ne.0)then
           call rnorma(ntau,continuoh,grp,grpmax)
           if(iautomatico.ne.1.or.mnodos(2).eq.1)then
              call nodos(grp,ntau,tau,pe,mnodos(2),grpmax)    !OJO nodosp escala con la pe en el ultimo nodo 
              do ie=1,ntau
                 atmoserr(2*ntau+ie)=grpmax(ie)
              end do  
           end if   
        end if                                    
        if(mnodos(3).ne.0)then
           call rnorma(ntau,continuoh,grm)
           if(iautomatico.ne.1.or.mnodos(3).eq.1)then
              call nodos(grm,ntau,tau,vtur,mnodos(3),grmmax) 
              do ie=1,ntau
                 atmoserr(3*ntau+ie)=grmmax(ie)
              end do  
           end if   
        end if
        if(mnodos(5).ne.0)then
          call rnorma(ntau,continuoh,grv)
          if(iautomatico.ne.1.or.mnodos(5).eq.1)then
             call nodos(grv,ntau,tau,vof,mnodos(5),grvmax)
             do ie=1,ntau
                atmoserr(5*ntau+ie)=grvmax(ie)
             end do  
          end if             
        end if

c las f. respuesta salen ordenadas en longitud de onda,perfil y tau
c rt(tau1:i(l1,l2,...);tau2:i......)
	      do j=1,mnodos(1)
	         iktt=(j-1)*ntotal+ikk1
		     rt(iktt)=grt(j)
	      end do   
	      do j=1,mnodos(2)
	         iktt=(j-1)*ntotal+ikk1
		     rp(iktt)=grp(j)
	      end do   
	      do j=1,mnodos(3)
	         iktt=(j-1)*ntotal+ikk1
		     rm(iktt)=grm(j)
	      end do   
	      do j=1,mnodos(5)
	         iktt=(j-1)*ntotal+ikk1
		     rv(iktt)=grv(j)
	      end do   

	     inten(ikk1)=icont/continuoh
9 	continue	!fin del do en lambda
999	continue	!fin del do en lineas

	return
	end
