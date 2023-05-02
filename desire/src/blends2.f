c_______________________________________________________________
c blends2 (rutina del programa de inversion 2 componentes DOS)
c se basa en la rutina "blends0" asi pues esta disenyado
c igual que la rutina "fstokes0" pero teniendo en cuenta blends
c calcula las funciones respuesta de los perfiles de stokes
c para un angulo heliocentrico de coseno mu. se supone lte,
c se basa en el programa "zeemanlines" de Wittmann y usa sus rutinas
c modificadas y en el "zeeman" de Rees usando el metodo DELO para integrar
c el maximo numero de puntos en tau es 64
c Basilio Ruiz 23-3-93
c Basilio y Jose Carlos 9-1-95 para pert. aditivas en fi
c Basilio y Jose Carlos 6-2-95 para pert. aditivas en gamma
c _________________________________________________________________
c
c ist=1 (i); =2 (q); =3 (u); =4 (v)

        subroutine blends2(atmos,stok,rt,rp,rh,rv,rg,rf,rm,mnodos,beta1,beta2,atmoserr)
        
        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        parameter (kt4=4*kt,kt16=16*kt,kt7=7*kt)
        parameter (kld4=4*kld) 
        real*4 melectron,mhidrogeno

c para las funciones respuesta
        real*4 rt(*),rp(*),rh(*),rv(*),rg(*),rf(*),rm(*)
        real*4 rt4(4,kt),rp4(4,kt),rh4(4,kt),rv4(4,kt),rg4(4,kt),rf4(4,kt)
        real*4 rm4(4,kt)
        real*4 grt(kt),grp(kt),grh(kt),grv(kt),grg(kt),grf(kt),grm(kt)
        real*4 x(kt)
        
c para los perfiles 
        real*4 stok(*),svec(4),continuoharr(kld)  !con_i(kl)

c para la atmosfera
        real*4 atmos(*),atmoserr(*)
        real*4 grtmax(kt),grpmax(kt),grmmax(kt),grhmax(kt)
        real*4 grvmax(kt),grgmax(kt),grfmax(kt)
        real*4 tau(kt),t(kt),pe(kt),vtur(kt),h(kt),vz(kt)
        real*4 vof(kt)
        real*4 taue(kt),gamma(kt),phi(kt),agamma(kt),aphi(kt)
        real*4 continuoh    
        integer mnodos(*)
        
c para la matriz de absorcion y sus derivadas
        real*4 dab(kt16),tk(kt16),pk(kt16),hk(kt16),vk(kt16),gk(kt16)
        real*4 fk(kt16),mk(kt16)
        real*4 dabtot(kt7,kld),gktot(kt7,kld),fktot(kt7,kld)
        real*4 tktot(kt7,kld),pktot(kt7,kld),hktot(kt7,kld)
        real*4 vktot(kt7,kld),mktot(kt7,kld)
        real*4 dabb(16),gkb(16),fkb(16),tkb(16),pkb(16),hkb(16),vkb(16)
        real*4 mkb(16)
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
        real*4 loggf,nair,mvdop,meta0,ma
        real*4 y(kt),table(kt,16)
        real*4 ck5(kt),dk5(kt),ddk5(kt),zeff
        real*4 ckappa,ckappa5,dkappa,dkappa5,ddkappa,ddkappa5
        
c para la emisividad y la funcion fuente
        real*8 beta1(kl,kt),beta2(kl,kt),blow,bup!,bratio
        real*8 wwwdbl,a_plck,c_plck,exalfa_plck,bpdob
        real*4 emisividad(kt,kld),emisividad_dt(kt,kld)
        real*4 bp(kt),bt(kt),bp0(kt),bt0(kt),dbp(kt)
        real*4 source(kt),source_dt(kt)
        integer linea_all_nlte(kl),lineaallnlte   !1= si algún blend de la linea iln es nlte       
c para el damping        
         real*4 chydro_arr(kl),xmu1_arr(kl),xmu2_arr(kl),xmu3_arr(kl)
         real*4 vv_arr(kl),beta_arr(kl)
        
c para el patron zeeman
        parameter (mc=20)       !numero maximo de componentes zeeman
        integer mult(2),ji(2),jf(2)
        real*4 tam(2),abu
        real*4 dlp(mc),dll(mc),dlr(mc),sp(mc),sl(mc),sr(mc)
        real*4 dlp_arr(mc,kl),dll_arr(mc,kl),dlr_arr(mc,kl)
        real*4 sp_arr(mc,kl),sl_arr(mc,kl),sr_arr(mc,kl)
        integer np_arr(kl),nl_arr(kl),nr_arr(kl)

c para las presiones parciales
        integer ivar(10)
        real*4 pg(99),dpg(99),ddpg(99),pi(10),dpi(10),ddpi(10)
        real*4 pt(kt,10),dpt(kt,10),ddpt(kt,10)
        real*4 pgas(kt),dpgas(kt),ddpgas(kt)      !,ro(kt),ck5_ro(kt)

c para la inclusion de RP en RT
        real*4 ax(kt),bx(kt),cx(kt),dx(kt),d1x(kt),d2x(kt),fx(kt)
        real*4 px(kt),qx(kt),rx(kt),sx(kt),tx(kt),wx(kt,kt)
        real*4 kac,kat,kap
        
c para hermite
        integer*4 intemethod  !integration method 0=herm_int,1=herm,2=bzr3,3=bzr3log
        real*4 deltae(kt),deltai(kt),delt2i(kt)
        real*8 saha_db,dsaha_db
        real*8 u12_db,u23_db,u33_db
        real*8 du12_db,du23_db,du33_db
        real*8 ddu12_db,ddu23_db,ddu33_db

c 05/05/20 epm: Error message.
        character*100 msg
             
c lugares comunes de memoria
        common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/brklm/ntotal_lines,atom_arr,istage_arr,alfa_arr,sigma_arr,wave_arr
        common/wavearrdble/wavedbl_arr
        common/datosatom/wvac_arr,weight_arr,croot_arr,wlengt1_arr,lambda_arr
        common/datosdamping/chydro_arr,xmu1_arr,xmu2_arr,xmu3_arr,vv_arr,beta_arr
        common/ldeo/dlongd,dlamda0  ! dlongd es cada delta de l.d.o. en ma
        common/loggfarr/gf_arr,energy_arr
        common/piis/piis
        common/yder/y,dyt,dyp,alpha !coef.abs.cont.y su der.t,p 
        common/segunda/tau,taue,deltae,deltai,delt2i
        common/offset/voffset  !para respuestas
        common/iautomatico/iautomatico
c       common/mu/cth   !este esta puesto a 1 en sir
        common/anguloheliocent/xmu 
        common/pgmag/ipgmag
c       common/continuos/con_i
        common/integrationmethod/intemethod  !integration method 0=herm_int,1=herm,2=bzr3,3=bzr3log 
        common/pot_ion/chi10_arr,chi20_arr,nel_arr,abu_arr
        common/patronzeeman/dlp_arr,dll_arr,dlr_arr,sp_arr,sl_arr,sr_arr,np_arr,nl_arr,nr_arr  
        common/componente_nlte/linea_all_nlte
        common/continuosarr/continuoharr
                    
        data iprimera/0/
                 
c nble es el numero de componentes de cada linea
        c=2.99792458e+10        !vel. de la luz en cm/seg
        piis=1./sqrt(3.1415926)
        g=xmu*2.7414e+4         !gravedad cm/s^2 en fotosfera solar   
        avog=6.023e23
        epsilon=.95  !cota para corregir RB por presiones Pmag<epsilon*Ptot

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
        epsilon2=epsilon/(1.-epsilon)  ! Pmag< epsilon2 * Pg
        coc3=1.212121          ! cociente polarizabilidad del H2 con HI
        coc2=.3181818          ! idem para el He

        ntotal=0
        do i=1,ntl
           do j=1,npas(i)
                ntotal=ntotal+1
           end do
        end do
        ists=ist(1)+ist(2)+ist(3)+ist(4)
        ntotal4=ntotal*ists
        
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
           grhmax(i)=0. 
           grvmax(i)=0.
           grgmax(i)=0.    
           grfmax(i)=0.
                            
           t(i)=atmos(i+ntau)
           pe(i)=atmos(i+2*ntau)
           
           if (pe(i).lt.0.) then
              write(msg,'(a,i4,a,1pe10.3)') ' pe(',i,')= ',pe(i)
              call error(KSTOP,'blends2','The electronic pressure turns out'
     &        //         ' to be negative at some optical depth\n'
     &        //         msg)
           end if
           vtur(i)=atmos(i+3*ntau)
           h(i)=atmos(i+4*ntau)
           vz(i) =atmos(i+5*ntau)       !velocidad linea de vision
           vof(i)=vz(i)-voffset
           gamma(i)=atmos(i+6*ntau)           
           phi(i)=atmos(i+7*ntau)
c          agamma(i)=tan(gamma(i)/2.0)
c          aphi(i)=tan(phi(i)/4.0)
           aphi(i)=1. !no qeremos pertb. relativas en fi sino aditivas
           agamma(i)=1. !no qeremos pertb. relativas en gamma sino aditivas
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
            call gasb(theta,ps,pg,dpg,ddpg)

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

c       call densidad(ntau,tau,t,pe,ck5_ro,pgas,ro) !ck5_ro=kap/ro
        ikk0=0
        ikk1=0
        ixx=0
        do 999 iln=1,ntl        !iln=numero de la linea,ntl num.total
           lineaallnlte=linea_all_nlte(iln) !=1 si alguna de las componentes es nlte
           do i=1,npas(iln)
              do k=1,ntau*7
                 dabtot(k,i)=0.
                 tktot(k,i)=0.
                 pktot(k,i)=0.
                 hktot(k,i)=0.
                 gktot(k,i)=0.
                 fktot(k,i)=0.
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

c subniveles zeeman
             np=np_arr(ixx)
             nl=nl_arr(ixx)
             nr=nr_arr(ixx)
        
        do i=1,np   
           dlp(i)=dlp_arr(i,ixx)
           sp(i)=sp_arr(i,ixx)
        end do   
        do i=1,nl
           dll(i)=dll_arr(i,ixx)
           sl(i)=sl_arr(i,ixx)
        end do 
        do i=1,nr
           dlr(i)=dlr_arr(i,ixx)
           sr(i)=sr_arr(i,ixx)
        end do 
       
c estratificacion en tau 5000
c genero una tabla lineal en logaritmo de tau a l.d.o.=5000.angs.
c desde tauini a taufin, con ntau valores.
        lambda=lambda_arr(ixx)
       
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

c calculo el damping 'a' mediante WITTMANN.
        if(abs(alfa).lt.1.e-25.or.abs(sigma).lt.1.e-25)then

           if(pg(1).gt.1.e-20)then !!!!!!!! CORRECCION PARA EVITAR DAMPING NULO POR PG(1)=0
        
              aj=chydro*(pg(1)/pg(90)*pg(91))*ts**0.3
              ai=(.992093+weinv)**.3+.6325*pg(2)/pg(1)*(.2498376+weinv)
     &           **.3+.48485*pg(89)/pg(1)*(.4960465+weinv)**.3
              a=(aj*ai+crad)/(12.5663706*vdop)
              ma=-a*mvdop/vdop
        
              daj=aj*(dpg(1)-dpg(90)+dpg(91)+.3/ts) !derivada de aj cont
              ddaj=aj*(ddpg(1)-ddpg(90)+ddpg(91)) !derivada de aj con p

              dai=.6325*pg(2)/pg(1)*(.2498376+weinv)**.3*(dpg(2)-dpg(1))
     &            +.48485*pg(89)/pg(1)*(.4960465+weinv)**.3*(dpg(89)
     &             -dpg(1))
              ddai=.6325*pg(2)/pg(1)*(.2498376+weinv)**.3*(ddpg(2)
     &             -ddpg(1))+.48485*pg(89)/pg(1)*(.4960465+weinv)**
     &             .3*(ddpg(89)-ddpg(1))

              da=(ai*daj+aj*dai)/(12.5663706*vdop)-a*dvdop/vdop
              dda=(ai*ddaj+aj*ddai)/(12.5663706*vdop)

           else !!!!!!!! CORRECCION PARA EVITAR DAMPING NULO POR PG(1)=0

              aj=chydro*(1./pg(90)*pg(91))*ts**0.3
              ai=.6325*pg(2)*(.2498376+weinv)**.3+
     &      .48485*pg(89)*(.4960465+weinv)**.3
              a=(aj*ai+crad)/(12.5663706*vdop)
              ma=-a*mvdop/vdop
        
              daj=aj*(-dpg(90)+dpg(91)+.3/ts) !derivada de aj cont
              ddaj=aj*(-ddpg(90)+ddpg(91)) !derivada de aj con p
              dai=.6325*pg(2)*(.2498376+weinv)**.3*dpg(2)+
     &      .48485*pg(89)*(.4960465+weinv)**.3*dpg(89)
              ddai=.6325*pg(2)*(.2498376+weinv)**.3*ddpg(2)+
     &       .48485*pg(89)*(.4960465+weinv)**.3*ddpg(89)

	    da=(ai*daj+aj*dai)/(12.5663706*vdop)-a*dvdop/vdop
	    dda=(ai*ddaj+aj*ddai)/(12.5663706*vdop)

           end if!!!!!!!! CORRECCION PARA EVITAR DAMPING NULO POR PG(1)=0

       else

c      calculo del damping 'a' mediante BARKLEM.       
c      d son las derivadas totales respecto a la temperatura.
c      dd son las derivadas totales respecto a la presion.
           dam=beta*(pg(91)/pg(90))*(pg(1)*xmu1**(-vv)+coc2*pg(2)*
     &        (xmu2**(-vv))+coc3*pg(89)*xmu3**(-vv))*ts**vv
        
           a=(1./(4.*pir))*(crad/vdop+dam/vdop2) 
           ma=-a*mvdop/vdop         

           ddam=dam*(dpg(91)-dpg(90))+beta*(pg(91)/pg(90))*(xmu1**(-vv)*
     &       pg(1)*dpg(1)+coc2*xmu2**(-vv)*pg(2)*dpg(2)+coc3*pg(89)*
     &       dpg(89)*xmu3**(-vv))*ts**(vv)+vv*(dam/ts)

           dddam=dam*(ddpg(91)-ddpg(90))+beta*(pg(91)/pg(90))*(xmu1**
     &           (-vv)*pg(1)*ddpg(1)+coc2*xmu2**(-vv)*pg(2)*ddpg(2)+
     &           coc3*pg(89)*ddpg(89)*xmu3**(-vv))

           da=(crad/(4.*pir))*(-dvdop/(vdop**2.))+(1./(4.*pir))*(ddam*
     &        vdop2-dam*dvdop2)/(vdop2**2.)

           dda=(1./(4.*pir))*(dddam/vdop2)

       endif

c calculo la funcion de planck en lambda y su derivada con t
c dtplanck" es la derivada de la f. de planck con temperatura.
                exalfa_plck=dexp(a_plck/ts)
                if(ible .eq. 1 .and. lineaallnlte .eq. 1)then 
                   bpdob=c_plck/(exalfa_plck-1.d0)
                   bp0(i)=real(bpdob)
                   bt0(i)=real(bpdob*bpdob*a_plck*exalfa_plck/(dble(ts)*dble(ts)*c_plck))
                end if
                
c                blow=beta1(ixx,i)
c                blows=real(blow)
c                bup=beta2(ixx,i)
c                bratio=blow/bup
                
                if(lineaallnlte .eq. 1)then
                   bpdob=c_plck*bup/(blow*exalfa_plck-bup)
                   bp(i)=real(bpdob)
                   bt(i)=real(bpdob*bpdob*a_plck*exalfa_plck*blow/(dble(ts)*dble(ts)*c_plck*bup))
                else
                   bpdob=c_plck/(exalfa_plck-1.d0)
                   source(i)=real(bpdob)
                   source_dt(i)=real(bpdob*bpdob*a_plck*exalfa_plck/(dble(ts)*dble(ts)*c_plck)) 
                end if

                eta0=eta0*blows
                deta0=deta0*blows
                ddeta0=ddeta0*blows
                meta0=meta0*blows

                y(i)=ckappa              !ckappa/ckappa5
                dyt(i)=dkappa
                dyp(i)=ddkappa
c                table(i,1)=dkappa        !derivada con t
c                table(i,2)=ddkappa       !    "     "  pe
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
71      continue            !do en log(tau)
 
c        call derivacuad(bp,dbp,ntau,tau)  

        amaxim=3.
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
        if(ible .eq. 1 .and. lineaallnlte .eq. 1)then !evaluo la emisividad del continuo
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
           dlamda=real((dble(dlongd(ikk))+(wlengt1dbl-wlengtdbl)*1.d3)*1.d-11)  !en cm.

c calculo etar,etal,etap y sus derivadas
           do 101 j=1,ntau     !do en tau
              k=(j-1)*7+1
              a=table(j,10)     !damping
              dldop=table(j,6) !desplazamiento doppler en cm
              dvcam=table(j,7) !campo de velocidades unidades doppler
              v=dlamda/dldop-dvcam !(l.d.o.+c.vel.) en unidades doppler
              hh=h(j)/dldop
              if(ible.eq.1)then
                t0=y(j)         !coef. absor. continuo/c.a.c 5000
                t1=dyt(j)       !derivada de t0 con t
                t2=dyp(j)       !  "         "      pe
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
              sg=sin(gamma(j))
              cg=cos(gamma(j))
              s2g=sg*sg
              c2g=1.+cg*cg
              scg=sg*cg
              sf=sin(2.*phi(j))
              cf=cos(2.*phi(j))
        
              call mvoigt(nr,dlr,sr,a,v,hh,t13,t14,wc,dldop,etar,vetar,
     &              getar,ettar,ettvr,esar,vesar,gesar,essar,essvr,
     &              ettmr,essmr)
              call mvoigt(nl,dll,sl,a,v,hh,t13,t14,wc,dldop,etal,vetal,
     &              getal,ettal,ettvl,esal,vesal,gesal,essal,essvl,
     &              ettml,essml)
              call mvoigt(np,dlp,sp,a,v,hh,t13,t14,wc,dldop,etap,vetap,
     &              getap,ettap,ettvp,esap,vesap,gesap,essap,essvp,
     &              ettmp,essmp)

c calculamos la matriz de absorcion
              tm=.5*(etar+etal)
              tn=.5*(etar-etal)
              sm=.5*(esar+esal)
              sn=.5*(esar-esal)
              tpm=.5*(etap-tm)
              spm=.5*(esap-sm)

              fi=t0 + t3 * .5 * ( etap*s2g + tm*c2g )
              fq=t3 * tpm * s2g * cf
              fu=t3 * tpm * s2g * sf
              fv=t3 * tn * cg

              fq1=t3 * spm * s2g * cf 
              fu1=t3 * spm * s2g * sf 
              fv1=t3 * sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,dabb)
              do iii=1,7
                 jjj=k+iii-1
                 dabtot(jjj,i)=dabtot(jjj,i)+dabb(iii)
              end do
              
c calculamos la emisividad
              if(lineaallnlte .eq. 1)then
                 emisividad(j,i)=emisividad(j,i)+bp(j)*t3*etap
                 emisividad_dt(j,i)=emisividad_dt(j,i)+bt(j)*t3*etap+bp(j)*(t4*etap+t3*(ettap*t11+ettvp))
              end if

c calculamos la derivada de la matriz de absorcion con gamma
              if(mnodos(6).ne.0)then

              fi=  2.* t3 * tpm * scg
              fq= fi * cf
              fu= fi * sf
              fv=-t3 * tn * sg

              fi1= 2.* t3 * spm * scg
              fq1=fi1 * cf 
              fu1=fi1 * sf 
              fv1=-t3 * sn * sg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,gkb)
              do iii=1,7
                 jjj=k+iii-1
                 gktot(jjj,i)=gktot(jjj,i)+gkb(iii)
              end do

              end if

c calculamos la derivada de la matriz de absorcion con phi
              if(mnodos(7).ne.0)then

              fi=  0.
              fq= -2.* t3 * tpm * s2g * sf
              fu=  2.* t3 * tpm * s2g * cf
              fv=  0.

              fq1=-2.* t3 * spm * s2g * sf 
              fu1= 2.* t3 * spm * s2g * cf 
              fv1= 0.

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,fkb)
              do iii=1,7
                 jjj=k+iii-1
                 fktot(jjj,i)=fktot(jjj,i)+fkb(iii)
              end do

              end if

c calculamos la derivada de la matriz de absorcion con t

              detar=t4*etar+t3*(ettar*t11+ettvr) !cambio el signo +ettvr
              detal=t4*etal+t3*(ettal*t11+ettvl)
              detap=t4*etap+t3*(ettap*t11+ettvp)
              desar=t4*esar+t3*(essar*t11+essvr)
              desal=t4*esal+t3*(essal*t11+essvl)
              desap=t4*esap+t3*(essap*t11+essvp)
              tm=.5*(detar+detal)
              tn=.5*(detar-detal)
              sm=.5*(desar+desal)
              sn=.5*(desar-desal)
              tpm=.5*(detap-tm)
              spm=.5*(desap-sm)

              fi=t1 + .5 * ( detap*s2g + tm*c2g )
              fq=tpm * s2g * cf
              fu=tpm * s2g * sf
              fv=tn * cg

              fq1=spm * s2g * cf 
              fu1=spm * s2g * sf 
              fv1=sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,tkb)
              do iii=1,7
                 jjj=k+iii-1
                 tktot(jjj,i)=tktot(jjj,i)+tkb(iii)
              end do

c calculamos la derivada de la matriz de absorcion con mic
              if(mnodos(3).ne.0)then

              detar=t15*etar+t3*(ettar*t16+ettmr)!cambio el signo +ettmr
              detal=t15*etal+t3*(ettal*t16+ettml)
              detap=t15*etap+t3*(ettap*t16+ettmp)
              desar=t15*esar+t3*(essar*t16+essmr)
              desal=t15*esal+t3*(essal*t16+essml)
              desap=t15*esap+t3*(essap*t16+essmp)
              tm=.5*(detar+detal)
              tn=.5*(detar-detal)
              sm=.5*(desar+desal)
              sn=.5*(desar-desal)
              tpm=.5*(detap-tm)
              spm=.5*(desap-sm)

              fi=.5 * ( detap*s2g + tm*c2g )
              fq=tpm * s2g * cf
              fu=tpm * s2g * sf
              fv=tn * cg

              fq1=spm * s2g * cf 
              fu1=spm * s2g * sf 
              fv1=sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,mkb)
              do iii=1,7
                 jjj=k+iii-1
                 mktot(jjj,i)=mktot(jjj,i)+mkb(iii)
              end do
              end if

c calculamos la derivada de la matriz de absorcion con pe
              if(mnodos(1).ne.0.or.mnodos(2).ne.0.or.mnodos(4).ne.0.or.
     &        mnodos(6).ne.0)then

              detar=t5*etar+t3*(ettar*t12) 
              detal=t5*etal+t3*(ettal*t12)
              detap=t5*etap+t3*(ettap*t12)
              desar=t5*esar+t3*(essar*t12)
              desal=t5*esal+t3*(essal*t12)
              desap=t5*esap+t3*(essap*t12)
              tm=.5*(detar+detal)
              tn=.5*(detar-detal)
              sm=.5*(desar+desal)
              sn=.5*(desar-desal)
              tpm=.5*(detap-tm)
              spm=.5*(desap-sm)

              fi=t2 + .5 * ( detap*s2g + tm*c2g )
              fq=tpm * s2g * cf
              fu=tpm * s2g * sf
              fv=tn * cg

              fq1=spm * s2g * cf 
              fu1=spm * s2g * sf 
              fv1=sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,pkb)
              do iii=1,7
                 jjj=k+iii-1
                 pktot(jjj,i)=pktot(jjj,i)+pkb(iii)
              end do
              end if

c calculamos la derivada de la matriz de absorcion con el campo
              if(mnodos(4).ne.0)then

              tm=.5*(getar+getal)
              tn=.5*(getar-getal)
              sm=.5*(gesar+gesal)
              sn=.5*(gesar-gesal)
              tpm=.5*(getap-tm)
              spm=.5*(gesap-sm)

              fi=t3 * .5 * ( getap*s2g + tm*c2g )
              fq=t3 * tpm * s2g * cf
              fu=t3 * tpm * s2g * sf
              fv=t3 * tn * cg

              fq1=t3 * spm * s2g * cf 
              fu1=t3 * spm * s2g * sf 
              fv1=t3 * sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,hkb)
              do iii=1,7
                 jjj=k+iii-1
                 hktot(jjj,i)=hktot(jjj,i)+hkb(iii)
              end do
              end if

c calculamos la derivada de la matriz de absorcion con la velocidad
              if(mnodos(5).ne.0)then

              tm=.5*(vetar+vetal)
              tn=.5*(vetar-vetal)
              sm=.5*(vesar+vesal)
              sn=.5*(vesar-vesal)
              tpm=.5*(vetap-tm)
              spm=.5*(vesap-sm)

              fi=t3 * .5 * ( vetap*s2g + tm*c2g )
              fq=t3 * tpm * s2g * cf
              fu=t3 * tpm * s2g * sf
              fv=t3 * tn * cg

              fq1=t3 * spm * s2g * cf 
              fu1=t3 * spm * s2g * sf 
              fv1=t3 * sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,vkb)
              do iii=1,7
                 jjj=k+iii-1
                 vktot(jjj,i)=vktot(jjj,i)+vkb(iii)
              end do
              
              end if
        
101          continue !fin do en tau(estamos aun dentro del do en lamda)
10        continue      !fin del do en lambda(pasos)
        end do   !fin del do en blends
   
        if(lineaallnlte .ne. 1)call derivacuad(source,dbp,ntau,tau)
        
        do 9 i=1,npas(iln)
          call matabs2(dabtot,i,ntau,dab,1)       
c se calcula donde poner el contorno 
           icontorno=1
           do j=1,ntau-1
              if(abs(dab(16*j-15))*deltae(j).gt.5.)icontorno=j
           end do
           if(mnodos(1).ne.0)call matabs2(tktot,i,ntau,tk,icontorno)
           if(mnodos(1).ne.0.or.mnodos(2).ne.0.or.mnodos(4).ne.0.or.
     &        mnodos(6).ne.0)call matabs2(pktot,i,ntau,pk,icontorno)
           if(mnodos(3).ne.0)call matabs2(mktot,i,ntau,mk,icontorno)
           if(mnodos(4).ne.0)call matabs2(hktot,i,ntau,hk,icontorno)
           if(mnodos(5).ne.0)call matabs2(vktot,i,ntau,vk,icontorno)
           if(mnodos(6).ne.0)call matabs2(gktot,i,ntau,gk,icontorno)
           if(mnodos(7).ne.0)call matabs2(fktot,i,ntau,fk,icontorno)
           
           if(lineaallnlte .eq. 1)then
              do j=1,icontorno-1
                 source(j)=bp0(j)
                 source_dt(j)=bt0(j)
              end do
              do j=icontorno,ntau
                 k=(j-1)*7+1
                 source(j)=emisividad(j,i)/dabtot(k,i) 
                 source_dt(j)=(emisividad_dt(j,i)-source(j)*tktot(k,i))/dabtot(k,i)
              end do
              call derivacuad(source,dbp,ntau,tau)
           end if
           if(intemethod .eq. 0) then
             call hermite(source,dbp,dab,ntau,svec,kt,source_dt,tk,pk,hk,
     &              vk,gk,fk,mk,rt4,rp4,rh4,rv4,rg4,rf4,rm4,mnodos,icontorno)
           else if (intemethod .eq. 1) then
              call hermite_no_interpol(source,dbp,dab,ntau,svec,kt,source_dt,tk,pk,hk,
     &              vk,gk,fk,mk,rt4,rp4,rh4,rv4,rg4,rf4,rm4,mnodos,icontorno)
           else if (intemethod .eq. 2) then
              call delo_bezier3(source,dbp,dab,ntau,svec,kt,source_dt,tk,pk,hk,
     &               vk,gk,fk,mk,rt4,rp4,rh4,rv4,rg4,rf4,rm4,mnodos)
           else if (intemethod .eq. 3) then
              call delo_bezier3log(source,dbp,dab,ntau,svec,kt,source_dt,tk,pk,hk,
     &                vk,gk,fk,mk,rt4,rp4,rh4,rv4,rg4,rf4,rm4,mnodos)
           else if (intemethod .eq. 4) then
              call delo_bezier3logcont(source,dbp,dab,ntau,svec,kt,source_dt,tk,pk,hk,
     &                vk,gk,fk,mk,rt4,rp4,rh4,rv4,rg4,rf4,rm4,mnodos)
           end if

c introducimos el paso en tau para pasar la integral sobre las f. resp.
c a sumatorio y normalizmos por el continuo

        ikk1=ikk1+1
        continuoh=continuoharr(ikk1)
        ikk4=ikk1-ntotal
        ikk5=ikk1-ntotal
        do isv=1,4
           ikk5=ikk5+ntotal
           if(ist(isv).eq.1)then
              ikk4=ikk4+ntotal

              if(mnodos(2).ne.0.or.mnodos(1).ne.0 .or. ipgmag .eq. 1)then
                 do kk=1,icontorno-1
                    grp(kk)=0.
                 end do                     
                 do kk=icontorno,ntau
                    grp(kk)=rp4(isv,kk)
                 end do 
              end if

              if(mnodos(1).ne.0)then
                 do kk=1,icontorno-1
                    grt(kk)=0.
                 end do 
                 do kk=icontorno,ntau   !introduccion de grp en grt
                    if(mnodos(2).eq.0)then
                       suma=0.
                       do kj=1,kk-2
                          suma=suma+wx(kj,kk)*grp(kj)
                       end do
                       correc=grp(kk)*px(kk)+suma
                       if(kk.gt.1)correc=correc+grp(kk-1)*qx(kk)
                       grt(kk)=rt4(isv,kk)+correc
                    else
                       grt(kk)=rt4(isv,kk)
                    end if
                 end do 
                 call rnorma(ntau,continuoh,grt)
                 if(iautomatico.ne.1.or.mnodos(1).eq.1)then
                    call nodos(grt,ntau,tau,t,mnodos(1),grtmax)
                    do ie=1,ntau
                       atmoserr(ntau+ie)=grtmax(ie)
                    end do
                 end if   
              end if
     
              if(mnodos(2).ne.0)then
                 call rnorma(ntau,continuoh,grp)
                 if(iautomatico.ne.1.or.mnodos(2).eq.1)then
                    call nodos(grp,ntau,tau,pe,mnodos(2),grpmax)
                    do ie=1,ntau
                       atmoserr(2*ntau+ie)=grpmax(ie)
                    end do  
                 end if
              end if

              if(mnodos(3).ne.0)then
                 do kk=1,icontorno-1
                    grm(kk)=0.
                 end do
                 do kk=icontorno,ntau
                    grm(kk)=rm4(isv,kk)
                 end do 
                 call rnorma(ntau,continuoh,grm)
                 if(iautomatico.ne.1.or.mnodos(3).eq.1)then
                   call nodos(grm,ntau,tau,vtur,mnodos(3),grmmax) 
                   do ie=1,ntau
                       atmoserr(3*ntau+ie)=grmmax(ie)
                   end do  
                 end if  
              end if

              if(mnodos(4).ne.0)then
                 do kk=1,icontorno-1
                   grh(kk)=0.
                 end do
                 if(mnodos(2).eq.0.and.ipgmag.eq.1)then
                    do kk=icontorno,ntau
                       correc=sin(gamma(kk))
                       correc=correc*correc*h(kk)/4./3.1415926  !B*sin^2(g)/4./pi
                       corrijorb=0.
                       if(correc*h(kk)/2..lt.epsilon2*pgas(kk))corrijorb=1.
                       correc=correc*corrijorb*cx(kk)
                       grh(kk)=rh4(isv,kk)-correc*rp4(isv,kk)
                    end do 
                 else
                    do kk=icontorno,ntau
                       grh(kk)=rh4(isv,kk)
                    end do 
                 end if
                 call rnorma(ntau,continuoh,grh)
                 if(iautomatico.ne.1.or.mnodos(4).eq.1)then
                    call nodos(grh,ntau,tau,h,mnodos(4),grhmax)
                    do ie=1,ntau
                       atmoserr(4*ntau+ie)=grhmax(ie)
                    end do  
                 end if 
             end if
                 
             if(mnodos(5).ne.0)then
                 do kk=1,icontorno-1
                   grv(kk)=0.
                 end do
                 do kk=icontorno,ntau
                    grv(kk)=rv4(isv,kk)
                 end do 
                 call rnorma(ntau,continuoh,grv)
                 if(iautomatico.ne.1.or.mnodos(5).eq.1)then
                    call nodos(grv,ntau,tau,vof,mnodos(5),grvmax)
                    do ie=1,ntau
                       atmoserr(5*ntau+ie)=grvmax(ie)
                    end do  
                 end if                     
              end if

              if(mnodos(6).ne.0)then
                 do kk=1,icontorno-1
                   grg(kk)=0.
                 end do
                 if(mnodos(2).eq.0.and.ipgmag.eq.1)then
                   do kk=icontorno,ntau
                      correc=sin(gamma(kk))*h(kk)*h(kk)/4./3.1415926  !B^2*sin(g)/4./pi
                      corrijorg=0.
                      if(correc*sin(gamma(kk))/2..lt.epsilon2*pgas(kk))corrijorg=1.
                      correc=correc*corrijorg*cos(gamma(k))*cx(kk)
                      grg(kk)=rg4(isv,kk)-correc*rp4(isv,kk)
                   end do 
                 else
                   do kk=icontorno,ntau
c                     grg(kk)=rg4(isv,kk)*2./(1.0+agamma(kk)*agamma(kk))
                      grg(kk)=rg4(isv,kk)
                   end do 
                 end if
                 call rnorma(ntau,continuoh,grg)
                 if(iautomatico.ne.1.or.mnodos(6).eq.1)then
                    call nodos(grg,ntau,tau,agamma,mnodos(6),grgmax)
                    do ie=1,ntau
                       atmoserr(6*ntau+ie)=grgmax(ie)
                    end do                 
                 end if   
              end if
              
              if(mnodos(7).ne.0)then
                 do kk=1,icontorno-1
                   grf(kk)=0.
                 end do
                 do kk=icontorno,ntau
c                   grf(kk)=rf4(isv,kk)*4./(1.+aphi(kk)*aphi(kk))
                    grf(kk)=rf4(isv,kk)
                 end do
                 call rnorma(ntau,continuoh,grf)
                 if(iautomatico.ne.1.or.mnodos(7).eq.1)then
                    call nodos(grf,ntau,tau,aphi,mnodos(7),grfmax)
                    do ie=1,ntau
                       atmoserr(7*ntau+ie)=grfmax(ie)
                    end do                 
                 end if
              end if

c las f. respuesta salen ordenadas en longitud de onda,perfil y tau
c rt(tau1:i(l1,l2,...),q(l1,..),....v(l1....);tau2:i......)
              do j=1,mnodos(1)
                 iktt=(j-1)*ntotal4+ikk4
                 rt(iktt)=grt(j)
              end do   
              do j=1,mnodos(2)
                 iktt=(j-1)*ntotal4+ikk4
                 rp(iktt)=grp(j)
              end do   
              do j=1,mnodos(3)
                 iktt=(j-1)*ntotal4+ikk4
                 rm(iktt)=grm(j)
              end do   
              do j=1,mnodos(4)
                 iktt=(j-1)*ntotal4+ikk4 !ojo las perturbaciones son relativas a
                 rh(iktt)=grh(j)         !los parametros en z no a linea vision
              end do   
              do j=1,mnodos(5)
                 iktt=(j-1)*ntotal4+ikk4
                 rv(iktt)=grv(j)
              end do   
              do j=1,mnodos(6)
                 iktt=(j-1)*ntotal4+ikk4
                 rg(iktt)=grg(j)
              end do   
              do j=1,mnodos(7)
                 iktt=(j-1)*ntotal4+ikk4
                 rf(iktt)=grf(j)
              end do  
 
              stok(ikk4)=svec(isv)/continuoh
              end if
            end do

9         continue      !fin del do en lambda
999     continue        !fin del do en lineas


        return
        end

c __________________________________________________________________________

c rutina nodos calcula las funciones respuesta equivalentes en los nodos

        subroutine nodos(gr,n,tau,yy,no,grmax)
        
        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        real*4 gr(*),g(kt),tau(*),x(kt),y(kt),yy(kt),zz(kt),f(kt,kt)
        real*4 grmax(*)
        common/new_evaluation/new_evaluation  !from desire.f to evaluate uncertainties if new_evaluation=2
        
        if(no.le.0)return
        
c        if(new_evaluation .ge. 2)then
           do i=1,n
              if(abs(gr(i)) .gt. grmax(i))grmax(i)=abs(gr(i))
           end do
c        end if 

        if(no.eq.1)then   !si perturbacion constante

           sumagr=0.
           yymedio=0.
 
           do i=1,n
c              sumagr=sumagr+gr(i)*yy(i)  !pert. multiplicativa constante
               sumagr=sumagr+gr(i)        !pert. aditiva constante
               yymedio=yymedio+yy(i)      !necesario si pert. aditiva constante
           end do
c          gr(1)=sumagr                  !pert. multiplicativa constante
           gr(1)=sumagr*yymedio/float(n)  !necesario si pert. aditiva constante: escalamos al velor medio del parametro
        else

           do i=1,n
              zz(i)=yy(i)       !guardamos la variable para no perderla
           end do

           m=(n-1)/(no-1)   !n es el numero de taus
                            !m es el numero de pasos en cada subintervalo
           do i=1,no        !no es el numero de nodos total i el indice del nodo
              j=(i-1)*m+1   !j el indice en tau
              x(i)=tau(j)   !x(i) es el valor de tau en el nodo
              y(i)=zz(j)    !y(i) es el valor del parametro en el nodo
           end do

           call splines22(x,y,no-2,n,tau,zz,f)

           do i=1,no
              ggg=0.
              do j=1,n
                 ggg=ggg+gr(j)*f(j,i)
              end do
              g(i)=ggg
           end do
           do i=1,no
              gr(i)=g(i)*y(i)
           end do
        end if 

        return
        end

c __________________________________________________________________________


c matabs rutina que llena la matriz de absorcion

        subroutine matabs(fi,fq,fu,fv,fq1,fu1,fv1,dab)

        implicit real*4 (a-h,o-z)

        real*4 dab(*)

        dab(1)=fi
        dab(2)=fq
        dab(3)=fu
        dab(4)=fv
        dab(5)=fu1
        dab(6)=fq1
        dab(7)=fv1

        return
        end

c__________________________________________________________________________

c matabs2 rutina que construye la matriz de absorcion a partir de los 7

        subroutine matabs2(dab7,ifrec,ntau,dab,icontorno)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        parameter (kt7=7*kt)

        real*4 dab7(kt7,kld)
        real*4 dab(*)

        do jj=icontorno,ntau
           j=16*(jj-1)
           j7=7*(jj-1)

           fi=dab7(j7+1,ifrec)
           fq=dab7(j7+2,ifrec)

           fu=dab7(j7+3,ifrec)
           fv=dab7(j7+4,ifrec)
           fu1=dab7(j7+5,ifrec)
           fq1=dab7(j7+6,ifrec)
           fv1=dab7(j7+7,ifrec)

           dab(j+1)=fi
           dab(j+2)=fq
           dab(j+3)=fu 
           dab(j+4)=fv

           dab(j+5)=fq
           dab(j+6)=fi
           dab(j+7)=-fv1
           dab(j+8)=fu1 

           dab(j+9)=fu
           dab(j+10)=fv1
           dab(j+11)=fi
           dab(j+12)=-fq1

           dab(j+13)=fv
           dab(j+14)=-fu1
           dab(j+15)=fq1
           dab(j+16)=fi
        end do

        return
        end

c _________________________________________________________________

        subroutine densidad(ntau,tau,t,pe,kap,pg,ro)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        real*4 tau(*),t(*),pe(*),kap(kt),pg(*),ro(*),pp(10),kac,d1(10),d2,d3
        real wgt,abu,ei1,ei2
        common/constantes/g,avog        !gravedad,n. avogadro/pmu
        common/anguloheliocent/xmu 
c       common/mu/cth                   !esto esta puesto a 1 en sir
      
        g=2.7414e+4*xmu         !gravedad cm/s^2 en fotosfera solar
        avog=6.023e23
        bol=1.380662e-16      !cte boltzmann cgs

        do i=1,10
           d1(i)=0
        end do

c calculamos el peso molecular medio pmu
        pmusum=0.0
        asum=0.0
        do i=1,92
           ii=i
           call neldatb(ii,0.,wgt,abu,ei1,ei2)
           asum=asum+abu                
           pmusum=pmusum+wgt*abu
        end do

        do i=1,ntau
           tsi=t(i)
           psi=pe(i)
           call gasc(tsi,psi,psg,pp)
           pg(i)=psg

           call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d3)

c  pp(8) es la abundancia electronica =pe/p(h') 
           pesomedio=pmusum/(asum+pp(8)) !peso molec. medio

           kap(i)=kac*avog/pmusum
           ro(i)=pg(i)*pesomedio/bol/t(i)/avog 
        end do
 
        return
        end
