c fperfil2
c dada una atmosfera (2 componentes) "atmos" que entra por el common
c y dada una perturbacion "atmosr" esta rutina construye la nueva 
c atmosfera (mediante amp2), calcula los nuevos perfiles "scal" y las nuevas
c fr (blends0), convoluciona los perfiles y las fr con la macro y 
c la psf (deconv y deconv2) y calcula las derivadas de la chi^2 que 
c salen por "dscal"
c
c Basilio 23-3-93
c Modificacion luz difusa Basilio 20-6-94
c ________________________________________________________________________
        subroutine fperfil2(m,atmosr,scal,dscal,scal_RH)

        implicit real*4 (a-h,o-z)
 
        include 'PARAMETER' !incluidng vthresh  the threshold for velocity
        parameter (kt8=8*kt+2,kt16=16*kt+5,kt11=11*kt+2,kt12=11*kt+3)
        parameter (kl4=4*kl)    !numero maximo de lineas
        parameter (kld4=4*kld)
        parameter (kldt=kld*kn,kldt4=4*kldt)

c para la malla
        integer nlin(kl),npas(kl),ist(4),m(*),mdata(18) !,mnod1(8),mnod2(8)
        integer numerical(18)
        integer mnod1tot(8),mnod2tot(8),npos(kld4),ntl,ntls
        integer mprimera(18)
        integer nlins(kl4),npass(kl4),nble(kl),ifiltro,nlam_LTE,nRFpoints
        real*4 dlamda0(kl),dlamda0s(kl4),dlamda(kld),dlamdas(kld4)

c para la atmosfera
        integer ntau
        real*4 atmosr(*),atmos(kt16),atmostry(kt16),vmac1,vmac2
        real*4 atmos1(kt8),atmos2(kt8),vof(kt),gam1(kt),fi1(kt)
        real*4 atmos1err(kt8),atmos2err(kt8)
        real*4 atmos1LG(kt12),atmos2LG(kt12) 

c para los perfiles y f. respuesta
        real*4 scal(kld4),scal1(kld4),scal2(kld4),dscal(*),sin01(kld4),sin02(kld4),stray(kld4)
        real*4 scal_corregida(kld4)
        real*4 stok(kld4),difer(kld4),sig(kld4),diff_NLTE_LTE(kld4)
        real*4 rt1(kldt4),rh1(kldt4),rv1(kldt4),rg1(kldt4),rf1(kldt4)
        real*4 rp1(kldt4),rm1(kldt4)
        real*4 rt2(kldt4),rh2(kldt4),rv2(kldt4),rg2(kldt4),rf2(kldt4)
        real*4 rp2(kldt4),rm2(kldt4)
        integer icalerr
        real*4 pg1(kt),z1(kt),ro1(kt)
        real*4 pg2(kt),z2(kt),ro2(kt)
        real*4 stok_RH_1(kld4),stok_RH_2(kld4)
        real*4 scal_RH(kld4)
        
        character*100 RH_model,RH_magneticfield,label_ID_model !BRC-RH Jun 20 2017

c departure coefficients
        integer iRH1,iRH2
        integer ntotal_lines,store_ntot
        integer atomic_number,atom_arr(kl),istage_arr(kl)
        real*4 alfa_arr(kl),sigma_arr(kl),wave_arr(kl)
        integer icallingRH1,icallingRH2
        real*8 beta1_1(kl,kt), beta2_1(kl,kt)  
        real*8 beta1_2(kl,kt), beta2_2(kl,kt)  
        real*4 hydro1(6,kt),hydro2(6,kt),popH(6),hydroRH(6*kt)
        real*4 betaH1(6,kt),betaH2(6,kt)
        save betaH1,betaH2  !conserva los valores para la siguiente entrada
        
c comunes
        common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/ldeo/dlamda,dlamda0
        common/smalla/ntls,nlins,npass
        common/smalla1/dlamdas
        common/smalla0/dlamda0s        
        common/atmosfera/atmos
        common/atmosferaout/atmostry
        common/iamplioold/iamplioold
        common/nohaycampo/nohaycampo1,nohaycampo2
        common/primera2/ntotal,ntotal4,ists
        common/nciclos/nciclos   !del principal
        common/difusa/stray
        common/numtotblends/ntlblends
        common/iautomatico/iautomatico
        common/observaciones/stok,sig
        common/posiciones/npos 
        common/mprimera/mprimera ! para guardar m inicial para iteaciones posteriores 
        common/ndata/ndata !para fperfil2 y automatico
        common/offset/voffset  !para respuestas
        common/iprimeravez/iprimeravez
        common/calerr/icalerr !si calculo errores=1 else =0
        common/ifiltro/ifiltro
        common/iRH/iRH1,iRH2  ! eq 1 we call RH, 0 don't call RH
        common/RHnames/label_ID_model,RH_model,RH_magneticfield
        common/zetas/pg1,z1,ro1,pg2,z2,ro2    !para el calculo de z,pg y ro en amp22 y amp22err
        common/departcoef_atm1/beta1_1,beta2_1 
        common/departcoef_atm2/beta1_2,beta2_2 
        common/alog10mu/alog10mu             !alog10mu=alog10(xmu) from desire

        common/numerical/numerical
        common/numero_LTE/nlam_LTE
        common/brklm/ntotal_lines,atom_arr,istage_arr,alfa_arr,sigma_arr,wave_arr
        common/hydroge_populations/hydro1,hydro2    !from deSIRe float (6,kt),H populations 
        common/scalemodelRH/imassortau  !integer    0=logtau, 1=mass column from deSIRe
        common/atmoserr/atmos1err,atmos2err  !to desire for error evaluation

        data store_ntot/1/
        data idiff/0/

        epsilon=2.e-5 !precision para comparar reales
        
        if(iprimeravez .eq. 0)then
c datos de la linea
            ifrec=0
            do iln=1,ntl        !iln=numero de la linea,ntl num.total
               do ii=1,npas(iln)
                 ifrec=ifrec+1
               end do          
            end do  
            nlam_LTE=ifrec
        end if  
        
        nohaycampo=nohaycampo1+nohaycampo2
        if(nohaycampo .eq. 0)nRFpoints=nlam_LTE
        if(nohaycampo .ne. 0)nRFpoints=(ist(1)+ist(2)+ist(3)+ist(4))*nlam_LTE

        icallingRH1=iRH1
        icallingRH2=iRH2

        if(iprimeravez.eq.0)then  !The first iteration for every cycle iprimeravez=0
             ixx=0
             ikk0=0
             iii=0

             do iln=1,ntl
                do ible=1,nble(iln)
                   ixx=ixx+1
                   do j=1,ntau
                      beta1_1(ixx,j)=1.0d0
                      beta2_1(ixx,j)=1.0d0
                      beta1_2(ixx,j)=1.0d0
                      beta2_2(ixx,j)=1.0d0            
                   end do
                end do   
             end do
             ntotal_lines=ixx
             store_ntot=ntotal_lines
             do iH=1,6
                do j=1,ntau
                    betaH1(iH,j)=1.00
                    betaH2(iH,j)=1.00
                end do
             end do
             do i=1,18
                mprimera(i)=m(i) ! guardamos los nodos de entrada
             end do 
           
              if(nohaycampo.eq.0)then
             ntls=ntl
             do i=1,ntl
                npass(i)=npas(i)
             end do
             do i=2,4
                ist(i)=0
             end do
             ist(1)=1
          end if
          ntotal=0
          iii=0
          do i=1,ntl
             ntotal=ntotal+npas(i)
             do ii=1,nble(i)
                iii=iii+1 
             end do
             ntlblends=iii
          end do
          ists=ist(1)+ist(2)+ist(3)+ist(4)
          ntotal4=ntotal*ists
          voffset=-vthresh 
        end if
        iprimeravez=iprimeravez+1

        nli=0
        do i=1,ntl
           nli=nli+npas(i)
        end do
        
c copiamos la atmosfera antigua en atmostry; "amp2" escribira en 
c atmostry la nueva atmosfera
        do i=1,16*ntau+5
           atmostry(i)=atmos(i)
        end do
        if(nciclos.ne.0.and.iprimeravez.ne.1)call amp2(ntau,m,atmostry,atmosr)
    
        do i=1,16
           mdata(i)=i*ntau+2*int(i/9)  ! indi. anterior a la var. i (ampliada)
        end do
        mdata(17)=mdata(16)+1
        mdata(18)=mdata(17)+1

c la macro y los f.f son
        vmac1=atmostry(mdata(8)+1)
        vmac2=atmostry(mdata(16)+1)
        fill2=atmostry(mdata(17)+1)
        fill1=1.0-fill2
        peso2=atmostry(mdata(18)+1)/100.
        peso1=1.-peso2

        do i=1,ntotal4
           scal1(i)=0.
           scal2(i)=0.
           sin01(i)=0.
           sin02(i)=0.
        end do
        
        if(icalerr.eq.1)then
           do i=1,8
              mnod1tot(i)=m(i)
              mnod2tot(i)=m(i+8)
           end do
        else
           do i=1,8
              if(mprimera(i).gt.1.and.iautomatico.eq.1)then
                 mnod1tot(i)=ntau
              else
                 mnod1tot(i)=mprimera(i)
              end if
              if(mprimera(i+8).gt.1.and.iautomatico.eq.1)then
                 mnod2tot(i)=ntau
              else
                 mnod2tot(i)=mprimera(i+8)
              end if
           end do
        end if

        ntotal_lines=store_ntot
c dividimos la atmosfera en las 2 componentes
        do i=1,ntau
           atmos1(i)=atmostry(i)                      !in SIR case, tau is in the LoS Ref. Frame
           atmos1LG(i)=atmostry(i)+alog10mu           !tau correction: RH input in Local Ref Frame
           
           atmos2(i)=atmostry(i+8*ntau+2)
           atmos2LG(i)=atmostry(i+8*ntau+2)+alog10mu
        end do
        do i=ntau+1,8*ntau+2
           atmos1(i)=atmostry(i)
           atmos1LG(i)=atmostry(i)            !atmos1LG(8*ntau+1)=Vmac1 , atmos1LG(8*ntau+2)=fill1
           atmos2(i)=atmostry(i+8*ntau+2)
           atmos2LG(i)=atmostry(i+8*ntau+2)   !atmos2LG(8*ntau+1)=Vmac2 , atmos2LG(8*ntau+2)=fill2
        end do
        do i=1,ntau
           atmos1LG(8*ntau+2+i)=z1(i)
           atmos1LG(9*ntau+2+i)=pg1(i)
           atmos1LG(10*ntau+2+i)=ro1(i)
           atmos2LG(8*ntau+2+i)=z2(i)
           atmos2LG(9*ntau+2+i)=pg2(i)
           atmos2LG(10*ntau+2+i)=ro2(i)
        end do  
        atmos1LG(11*ntau+3)=atmostry(mdata(18)+1)  !stray
        atmos2LG(11*ntau+3)=atmostry(mdata(18)+1)  !stray (copied)

c calculamos las funciones respuesta y el perfil observado para cada atmosfera
c solo en el caso de que fill1.ne.0 se hace lo siguiente
c ************************* atmosfera 1 ***************************************
        k=0        !inicializo el indice de las derivadas
        if(abs(fill1).gt.epsilon)then
           if(nohaycampo1.eq.0.and.ist(1).ne.0)then
              numer=0
              numer=numerical(1)+numerical(3)+numerical(5)!T1, mic1, vz1
              if(icallingRH1.eq.1)then 
                 if(iprimeravez .ge. 2 .and. imassortau .eq. 0)then
                    do i=1,ntau
                       call hpopulations(atmos1(i+ntau),atmos1(i+2*ntau),popH)
                       do j=1,6
                          hydro1(j,i)=popH(j)*betaH1(j,i)  !going out through common to write_atmos_RH_tau call by departures and main
                       end do 
                    end do
                 end if
                 if(abs(atmos1LG(5*ntau+1)) .lt. 1)atmos1LG(5*ntau+1)=1
                 call departures(label_ID_model,RH_model,RH_magneticfield,1,atmos1LG,
     &                           ntau,ntotal_lines,beta1_1,beta2_1,stok_RH_1)

                 if(imassortau .eq. 0)then         
                    call rhhpop(ntau,hydroRH)
                    do i=1,ntau
                       do j=1,6
                          betaH1(j,i)=hydroRH((i-1)*6+j)/hydro1(j,i)
                       end do
                    end do
                 end if   
              end if   
              call blendscon2(atmos1,scal1,rt1,rp1,rv1,rm1,
     &                        mnod1tot,beta1_1,beta2_1,atmos1err)
              if(icallingRH1.eq.1)then !debe estar aquí despues deblndescon2 porque dlamda0 se define en blendscon2
                 if(vmac1.gt.0. .or. ifiltro.ge.1)then
                    do j=ntotal+1,ndata
                       stok_RH_1(j)=0.
                    end do
                    call deconv(stok_RH_1,1,ntl,npas,dlamda0,dlamda,vmac1)  !aqui convolucionamos el perfil NLTE (sin campo)
                 end if 
              end if

c              if(numer .lt. 1)then
                   if(vmac1.gt.0. .or. ifiltro.ge.1)then
                      do j=1,ntotal
                         sin01(j)=scal1(j)
                      end do
                      do j=ntotal+1,ndata
                         scal1(j)=0.
                         sin01(j)=0.
                      end do
                      call deconv(scal1,1,ntl,npas,dlamda0,dlamda,vmac1)!aqui convolucionamos el perfil LTE (sin campo)
                      call deconv2(sin01,1,ntl,npas,dlamda0,dlamda,vmac1)!y aqui con la derivada de la macro
                   end if 
c              else  !in this case does not exist scal1 nor sin01
              if(numer .ge. 1  .and. icallingRH1.eq.1)then
                 call numericalsub_con(1,atmos1LG,mnod1tot,             !llamamos a la numerica: entra el perfil NLTE ya convolucinado
     &                             stok_RH_1,vmac1,ifiltro,ntotal,ndata,
     &                             rt1,rv1,rm1) 
c                 do j=1,ntotal
c                    scal1(j)=stok_RH_1(j)                               !sobre-escribimos el perfil NLTE sobre el LTE en caso numerico
c                 end do   
              end if  
           else  !o sea si hay campo1 
              numer=0
              do ii=1,7
                 numer=numer+numerical(ii) !T1, pe1, mic1, h1, vz1, g1, f1
              end do 
              if(icallingRH1.eq.1)then   
                 if(iprimeravez .ge. 2 .and. imassortau .eq. 0)then
                    do i=1,ntau
                       call hpopulations(atmos1(i+ntau),atmos1(i+2*ntau),popH)
                       do j=1,6
                          hydro1(j,i)=popH(j)*betaH1(j,i)
                       end do 
                    end do
                 end if

                 call departures(label_ID_model,RH_model,
     &                   RH_magneticfield,1,atmos1LG,ntau,ntotal_lines,
     &                   beta1_1,beta2_1,stok_RH_1)

                 if(imassortau .eq. 0)then       
                    call rhhpop(ntau,hydroRH)
                    do i=1,ntau
                       do j=1,6
                          betaH1(j,i)=hydroRH((i-1)*6+j)/hydro1(j,i)
                       end do
                    end do
                 end if
              end if   

              call blends2(atmos1,scal1,rt1,rp1,rh1,rv1,rg1,rf1,rm1,
     &                     mnod1tot,beta1_1,beta2_1,atmos1err)
c Dado que dlamda0 entra en el common via blends2 la sentencia siguiente 
c no puede colocarse antes de la llamada a blends2  
            k1=0
            do i=1,4
               do j=1,ist(i)
                  do klin=1,ntl
                     k1=k1+1
                     dlamda0s(k1)=dlamda0(klin)             
                  end do
               end do
            end do
            if(icallingRH1.eq.1)then                             
               if(vmac1.gt.0 .or. ifiltro.ge.1)then  
                  call deconv(stok_RH_1,1,ntls,npass,dlamda0s,dlamdas,vmac1)!convolucionamos los perfiles NLTE
               end if
            end if    
c calculamos la convolucion con el perfil gaussiano y la psf
c calculamos la convolucion con la derivada del perfil gaussiano y la psf
            if(vmac1.gt.0 .or. ifiltro.ge.1)then
               do j=1,ntotal4
                  sin01(j)=scal1(j)
               end do
               call deconv(scal1,1,ntls,npass,dlamda0s,dlamdas,vmac1)
               call deconv2(sin01,1,ntls,npass,dlamda0s,dlamdas,vmac1)
            end if 
            if(numer .ge. 1 .and. icallingRH1.eq.1 )then
               call numericalsub(1,atmos1LG,mnod1tot,
     &                             stok_RH_1,vmac1,ifiltro,ntotal,ndata,
     &                             rt1,rp1,rh1,rv1,rg1,rf1,rm1)   
            end if 
          end if !fin del if sobre el campo (linea 272 else at 322)
        end if !fin del if sobre el filling factor: atmosfera 1 (linea 271)

c ************************* atmosfera 2 ***************************************
        if(abs(fill2).gt.epsilon)then
           if(nohaycampo2.eq.0.and.ist(1).ne.0)then
              numer=0
              if(icallingRH2.eq.1)then 
                 numer=numerical(9)+numerical(11)+numerical(13) !T2, mic2, vz2
                 if(iprimeravez .ge. 2 .and. imassortau .eq. 0)then
                    do i=1,ntau
                       call hpopulations(atmos2(i+ntau),atmos2(i+2*ntau),popH)
                       do j=1,6
                          hydro2(j,i)=popH(j)*betaH2(j,i)  !going out through common to write_atmos_RH_tau call by departures and main
                       end do 
                    end do
                 end if
                 if(abs(atmos2LG(5*ntau+1)) .lt. 1)atmos2LG(5*ntau+1)=1
                 call departures(label_ID_model,RH_model,RH_magneticfield,2,atmos2LG,
     &                           ntau,ntotal_lines,beta1_2,beta2_2,stok_RH_2) 
                 if(imassortau .eq. 0)then  
                    call rhhpop(ntau,hydroRH)
                    do i=1,ntau
                       do j=1,6
                          betaH2(j,i)=hydroRH((i-1)*6+j)/hydro2(j,i)
                       end do
                    end do
                 end if
                 if(vmac2.gt.0. .or. ifiltro.ge.1)then
                    do j=ntotal+1,ndata
                       stok_RH_2(j)=0.
                    end do
                    call deconv(stok_RH_2,1,ntl,npas,dlamda0,dlamda,vmac2)
                 end if  
              end if             
              call blendscon2(atmos2,scal2,rt2,rp2,rv2,rm2,mnod2tot,beta1_2,beta2_2,atmos2err)
              if(numer .ge. 1)then
                 call error(KSTOP,'fperfil2','No numerical RF are implemented'
     &           //         ' for second component')
              end if   
c              if(numer .ge. 1)call numericalsub_con(2,atmos2LG,mnod2tot,stok_RH_2,vmac2,ifiltro,ntotal,ndata,rt2,rv2,rm2)              
              if(vmac2.gt.0 .or. ifiltro.ge.1)then
                 do j=1,ntotal
                    sin02(j)=scal2(j)
                 end do
                 do j=ntotal+1,ndata
                   scal2(j)=0.
                   sin02(j)=0.
                 end do

                 call deconv(scal2,1,ntl,npas,dlamda0,dlamda,vmac2)
                 call deconv2(sin02,1,ntl,npas,dlamda0,dlamda,vmac2)
              end if
           else  !o sea si hay campo2
              if(icallingRH2.eq.1)then 
                 if(iprimeravez .ge. 2 .and. imassortau .eq. 0)then
                    do i=1,ntau
                       call hpopulations(atmos2(i+ntau),atmos2(i+2*ntau),popH)
                       do j=1,6
                          hydro2(j,i)=popH(j)*betaH2(j,i)  !going out through common to write_atmos_RH_tau call by departures and main
                       end do 
                    end do
                 end if
                 call departures(label_ID_model,RH_model,RH_magneticfield,2,atmos2LG,
     &                           ntau,ntotal_lines,beta1_2,beta2_2,stok_RH_2) 
                 if(imassortau .eq. 0)then
                    call rhhpop(ntau,hydroRH)
                    do i=1,ntau
                       do j=1,6
                          betaH2(j,i)=hydroRH((i-1)*6+j)/hydro2(j,i)
                       end do
                    end do
                 end if 
                 if(vmac2.gt.0 .or. ifiltro.ge.1)call deconv(stok_RH_2,1,ntls,npass,dlamda0s,dlamdas,vmac2)
              end if
              call blends2(atmos2,scal2,rt2,rp2,rh2,rv2,rg2,rf2,rm2,mnod2tot,beta1_2,beta2_2,atmos2err)
c Dado que dlamda0 entra en el common via blends2 la sentencia siguiente 
c no puede colocarse antes de la llamada a blends2  
              k1=0
              do i=1,4
                 do j=1,ist(i)
                    do klin=1,ntl
                       k1=k1+1
                       dlamda0s(k1)=dlamda0(klin)
                    end do
                 end do
              end do
c calculamos la convolucion con la derivada del perfil gaussiano y la psf
              if(vmac2.gt.0 .or. ifiltro.ge.1)then
                 do i=1,ntotal4
                    sin02(i)=scal2(i)
                 end do
                 call deconv(scal2,1,ntls,npass,dlamda0s,dlamdas,vmac2)
                 call deconv2(sin02,1,ntls,npass,dlamda0s,dlamdas,vmac2)
              end if
           end if
        end if 

c *********************mezclamos las dos atmosferas****************************
c y los perfiles de salida seran (teniendo en cuenta la luz difusa)
        if(abs(fill2).gt.epsilon .or. peso2 .gt. 0)then ! i.e. if we have 2 atmospheres
           do j=1,ntotal4
              scal(j)=peso1*(scal1(j)*fill1+scal2(j)*fill2)+peso2*stray(j)
           end do        
           if(nlin(ntlblends).eq.0)then
              scal(nli)=(scal1(nli)-scal2(nli))/2.
              scal_RH(nli)=(stok_RH_1(nli)-stok_RH_2(nli))/2.
           end if    
           if( icallingRH1 .eq. 1 .or.  icallingRH2 .eq. 1 )then
              if(nlin(ntlblends).eq.0)then
                 scal_RH(nli)=(stok_RH_1(nli)-stok_RH_2(nli))/2.
              end if   
              idiff=1
              do i=1,ntotal4
                 scal_RH(i)=peso1*(stok_RH_1(i)*fill1+stok_RH_2(i)*fill2)+peso2*stray(i)
                 diff_NLTE_LTE(i)=scal_RH(i)-scal(i)
              end do
           end if
        else ! in case of only one atmosphere
           do j=1,ntotal4
              scal(j)=scal1(j)
           end do  
           if( icallingRH1 .eq. 1)then
              idiff=1
              do i=1,ntotal4
                 scal_RH(i)=stok_RH_1(i)
                 diff_NLTE_LTE(i)=scal_RH(i)-scal(i)
              end do
           end if
        end if
        do i=1,ntotal4
           scal_corregida(i)=scal(i)
        end do
        if( idiff .eq. 1) then !in anytime NLTE has been evaluated
           do i=1,ntotal4
             scal_corregida(i)=scal(i)+diff_NLTE_LTE(i) !synthetic profile with SIR +new/old difference NLTE
           end do
        end if

C version que estaba funcionando: para no modificar scal defino scal_corregida
c        do j=1,ndata
c           difer(j)=(stok(j)-scal(npos(j)))/sig(j)/sig(j)
c        end do
c stok son las observaciones
         do j=1,ndata
            difer(j)=(stok(j)-scal_corregida(npos(j)))/sig(j)/sig(j)
         end do
         
c *******************reevaluamos el numero de nodos y escribimos el vector derivada
        vmac10=vmac1
        if(numer .ge. 1)vmac10=-10. !to avoid convolution of (already convolved) numerical RF with Mac and/or filter 
        if(abs(fill1).gt.epsilon)then
           if(nohaycampo1.eq.0.and.ist(1).ne.0)then

              if(nlin(ntlblends).eq.0)then
                 do i=nli,mnod1tot(1)*ntotal,ntotal               
                    rt1(i)=rt1(i)/2./fill1
                 end do
                 do i=nli,mnod1tot(2)*ntotal,ntotal               
                    rp1(i)=rp1(i)/2./fill1
                 end do
              end if
            
              fill=fill1*peso1
              call sub1(mnod1tot(1),difer,npos,vmac10,fill,rt1,k,
     &                  dscal,atmos1(ntau+1),mprimera(1)) !tempe.
              call sub1(mnod1tot(2),difer,npos,vmac10,fill,rp1,k,            
     &                  dscal,atmos1(2*ntau+1),mprimera(2)) !presi. 
              call sub1(mnod1tot(3),difer,npos,vmac10,fill,rm1,k,
     &                  dscal,atmos1(3*ntau+1),mprimera(3)) !micro.
              call cero1(mnod1tot(4),ntotal4,k,dscal)           !campo 

              if(mnod1tot(5).gt.0)then
                 do i=1,ntau
                    vof(i)=atmos1(5*ntau+i)-voffset
                 end do
              end if

              call sub1(mnod1tot(5),difer,npos,vmac10,fill,rv1,k,
     &                  dscal,vof,mprimera(5)) !veloc.
              call cero1(mnod1tot(6),ntotal4,k,dscal)           !inclinacion 
              call cero1(mnod1tot(7),ntotal4,k,dscal)           !azimuth 
              call mult1(mnod1tot(8),ntotal4,k,fill,vmac10,sin01,dscal) !macro
              do i=1,8
                 m(i)=mnod1tot(i)   !nuevo numero de nodos
              end do

          else  !o sea si hay campo1 

              if(nlin(ntlblends).eq.0)then
                 do i=nli,mnod1tot(1)*ntotal4,ntotal4               
                    rt1(i)=rt1(i)/2./fill1
                 end do
                 do i=nli,mnod1tot(2)*ntotal4,ntotal4               
                    rp1(i)=rp1(i)/2./fill1
                 end do
              end if
        
c calculamos la convolucion de las funciones respuesta
              fill=fill1*peso1

              call sub2(mnod1tot(1),difer,npos,vmac10,fill,dlamda0s,ist,rt1,k,
     &                  dscal,atmos1(ntau+1),mprimera(1)) !tempe.
     
              call sub2(mnod1tot(2),difer,npos,vmac10,fill,dlamda0s,ist,rp1,k,  
     &                  dscal,atmos1(2*ntau+1),mprimera(2)) !presi.
              call sub2(mnod1tot(3),difer,npos,vmac10,fill,dlamda0s,ist,rm1,k,
     &                  dscal,atmos1(3*ntau+1),mprimera(3)) !micro.
              call sub2(mnod1tot(4),difer,npos,vmac10,fill,dlamda0s,ist,rh1,k,
     &                  dscal,atmos1(4*ntau+1),mprimera(4)) !campo

              if(mnod1tot(5).gt.0)then
                 do i=1,ntau
                    vof(i)=atmos1(5*ntau+i)-voffset
                 end do
              end if
              
              call sub2(mnod1tot(5),difer,npos,vmac10,fill,dlamda0s,ist,rv1,k,
     &                  dscal,vof,mprimera(5)) !veloc.

              if(mnod1tot(6).gt.0)then
                 do i=1,ntau
                    gam1(i)=1.
                 end do
              end if
              
              call sub2(mnod1tot(6),difer,npos,vmac10,fill,dlamda0s,ist,rg1,k,
     &                  dscal,gam1,mprimera(6)) !incli.
     
              if(mnod1tot(7).gt.0)then
                 do i=1,ntau
                    fi1(i)=1.
                 end do
              end if

              call sub2(mnod1tot(7),difer,npos,vmac10,fill,dlamda0s,ist,rf1,k,
     &                  dscal,fi1,mprimera(7)) !azimuth.
              
              call mult1(mnod1tot(8),ntotal4,k,fill,vmac10,sin01,dscal)   !macro

              do i=1,8
                 m(i)=mnod1tot(i)   !nuevo numero de nodos
              end do
           end if
        end if

c solo en el caso de que fill2.ne.0 se hace lo siguiente
c ************************* atmosfera 2 ***************************************

        if(abs(fill2).gt.epsilon)then
           if(nohaycampo2.eq.0.and.ist(1).ne.0)then

              if(nlin(ntlblends).eq.0)then
                 do i=nli,mnod2tot(1)*ntotal,ntotal               
                    rt2(i)=-rt2(i)/2./fill2
                 end do
                 do i=nli,mnod2tot(2)*ntotal,ntotal               
                    rp2(i)=-rp2(i)/2./fill2
                 end do
              end if
              fill=fill2*peso1

              call sub1(mnod2tot(1),difer,npos,vmac2,fill,rt2,k,
     &                  dscal,atmos2(ntau+1),mprimera(9)) !tempe.
              call sub1(mnod2tot(2),difer,npos,vmac2,fill,rp2,k,              
     &                  dscal,atmos2(2*ntau+1),mprimera(10)) !presi.
              call sub1(mnod2tot(3),difer,npos,vmac2,fill,rm2,k,
     &                  dscal,atmos2(3*ntau+1),mprimera(11)) !micro.
              call cero1(mnod2tot(4),ntotal4,k,dscal)           !campo 

              if(mnod2tot(5).gt.0)then
                 do i=1,ntau
                    vof(i)=atmos2(5*ntau+i)-voffset
                 end do
              end if

              call sub1(mnod2tot(5),difer,npos,vmac2,fill,rv2,k,
     &                  dscal,vof,mprimera(13)) !veloc.

              call cero1(mnod2tot(6),ntotal4,k,dscal)            !inclinacion  
              call cero1(mnod2tot(7),ntotal4,k,dscal)            !azimuth 
              call mult1(mnod2tot(8),ntotal4,k,fill,vmac2,sin02,dscal) !macro
                  do i=1,8
                     m(i+8)=mnod2tot(i)   !nuevo numero de nodos
                  end do

           else  !o sea si hay campo2 

c Dado que dlamda0 entra en el common via blends2 la sentencia siguiente 
c no puede colocarse antes de la llamada a blends2  
              k1=0
              do i=1,4
                 do j=1,ist(i)
                    do klin=1,ntl
                       k1=k1+1
                       dlamda0s(k1)=dlamda0(klin)
                    end do
                 end do
              end do

              if(nlin(ntlblends).eq.0)then
                 do i=nli,mnod2tot(1)*ntotal4,ntotal4               
                    rt2(i)=-rt2(i)/2./fill2
                 end do
                 do i=nli,mnod2tot(2)*ntotal4,ntotal4               
                    rp2(i)=-rp2(i)/2./fill2
                 end do
              end if

c calculamos la convolucion de las funciones respuesta
              fill=fill2*peso1

              call sub2(mnod2tot(1),difer,npos,vmac2,fill,dlamda0s,ist,rt2,k,
     &                  dscal,atmos2(ntau+1),mprimera(9)) !tempe.
              call sub2(mnod2tot(2),difer,npos,vmac2,fill,dlamda0s,ist,rp2,k, !ojo llamo a sub2
     &                  dscal,atmos2(2*ntau+1),mprimera(10)) !presi.
              call sub2(mnod2tot(3),difer,npos,vmac2,fill,dlamda0s,ist,rm2,k,
     &                  dscal,atmos2(3*ntau+1),mprimera(11)) !micro.
              call sub2(mnod2tot(4),difer,npos,vmac2,fill,dlamda0s,ist,rh2,k,
     &                  dscal,atmos2(4*ntau+1),mprimera(12)) !campo

              if(mnod2tot(5).gt.0)then
                 do i=1,ntau
                    vof(i)=atmos2(5*ntau+i)-voffset
                 end do
              end if

              call sub2(mnod2tot(5),difer,npos,vmac2,fill,dlamda0s,ist,rv2,k,
     &                  dscal,vof,mprimera(13)) !veloc.

              if(mnod2tot(6).gt.0)then
                 do i=1,ntau
                    gam1(i)=1.
                 end do
              end if
              call sub2(mnod2tot(6),difer,npos,vmac2,fill,dlamda0s,ist,rg2,k,
     &                  dscal,gam1,mprimera(14)) !incli.

              if(mnod2tot(7).gt.0)then
                 do i=1,ntau
                    fi1(i)=1.
                 end do
              end if
              call sub2(mnod2tot(7),difer,npos,vmac2,fill,dlamda0s,ist,rf2,k,
     &                  dscal,fi1,mprimera(15)) !azimuth.
              call mult1(mnod2tot(8),ntotal4,k,fill,vmac2,sin02,dscal)   !macro

              do i=1,8
                 m(i+8)=mnod2tot(i)   !nuevo numero de nodos
              end do

           end if

c la variable es (1-fill2)/fill2 . Su derivada respecto a fill2 es -1/fill2^2
c como tenemos modf. multiplicativas el factor es fill2*(fill2-1)
   
           if(m(17).eq.1)then
              factor=fill2*(fill2-1.)*peso1
              do j=1,ntotal4
                 k=k+1 
                 dscal(k)=factor*(scal2(j)-scal1(j)) 
              end do
           end if

        end if

c calculamos la funcion respuesta al peso de la luz difusa
c la variable es peso2/(1.-peso2). La derivada de peso2 respecto 
c a la variable es (1.-peso2)^2 Y teniendo en cuenta pert. mult.
c el factor es peso2*(1.-peso2)
        if(m(18).eq.1)then
           factor=peso1*peso2
           do j=1,ntotal4
              k=k+1 
              dscal(k)=factor*(stray(j)-scal1(j)*fill1-scal2(j)*fill2)
           end do
        end if

        if(nciclos.ne.0)call comprime2(ntau,m,atmostry,atmosr) 

        return
        end

c _____________________________________________________________________________

        subroutine sub1(mi,difer,npos,vmac,fill,rt,k,dscal,t,mp)

        implicit real*4 (a-h,o-z) 
        include 'PARAMETER' !por kt,kn,kl y kld
        parameter (kl4=4*kl)    !numero maximo de lineas
        parameter (kld4=4*kld)

        integer ist(4),nlin(kl),npas(kl),nble(kl),npos(*),ifiltro
        real*4 dlamda(kld),dlamda0(kl)
        real*4 rt(*),dscal(*),difer(*),t(*)
        common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/ldeo/dlamda,dlamda0
        common/primera2/ntotal,ntotal4,ists
        common/iautomatico/iautomatico
        common/ifiltro/ifiltro

        if(iautomatico.eq.1.and.mp.gt.1.and.mi.gt.1)call automatico(mi,mp,ntotal,difer,npos,rt,t) 

        do i=1,mi
           npun=ntotal*(i-1)+1
           if(vmac .gt. -5)then
              if(vmac.gt.0 .or. ifiltro .ge. 1)then
c El argumento rt(npun) pasa el array a partir de ese indice.
                 call deconv(rt(npun),1,ntl,npas,dlamda0,dlamda,vmac)
              end if
           end if   
        end do 

        do i=1,mi
           npun=ntotal*(i-1)+1
           do j=npun,npun+ntotal-1
              k=k+1 
              dscal(k)=rt(j)*fill      !t1(l1,l2...lntot4),t2(l1,l2...
           end do
           do j=ntotal+1,ntotal4
              k=k+1 
              dscal(k)=0.
           end do
        end do                   

        return
        end

c _____________________________________________________________________________

        subroutine sub2(mi,difer,npos,vmac,fill,dlamda0s,ist,rt,k,dscal,t,mp)

        implicit real*4 (a-h,o-z) 
        include 'PARAMETER' !por kt,kn,kl y kld

        parameter (kl4=4*kl)    !numero maximo de lineas
        parameter (kld4=4*kld)

        integer ist(4),nlins(kl4),npass(kl4),npos(*),ifiltro
        real*4 dlamda0s(*),dlamdas(kld4)
        real*4 rt(*),dscal(*),difer(*),t(*)
        common/smalla/ntls,nlins,npass
        common/smalla1/dlamdas
        common/primera2/ntotal,ntotal4,ists
        common/iautomatico/iautomatico
        common/ifiltro/ifiltro

        if(iautomatico.eq.1.and.mp.gt.1.and.mi.gt.1)call automatico(mi,mp,ntotal4,difer,npos,rt,t) 

        do i=1,mi
           npun=ntotal4*(i-1)+1
           if(vmac .gt. -5)then
              if(vmac.gt.0 .or. ifiltro .ge. 1)then
c             El argumento rt(npun) pasa el array a partir de ese indice.
                 call deconv(rt(npun),1,ntls,npass,dlamda0s,dlamdas,vmac)
              end if 
          end if     
        end do 

        do i=1,mi
           npun=ntotal4*(i-1)+1
           do j=npun,npun+ntotal4-1
              k=k+1 
              dscal(k)=rt(j)*fill      !t1(l1,l2...lntot4),t2(l1,l2...
           end do
        end do        
        return
        end
        
c _____________________________________________________________________________

        subroutine cero1(mi,ntotal4,k,dscal)
        implicit real*4 (a-h,o-z)

        real*4 dscal(*)

        do i=1,mi
           do j=1,ntotal4
              k=k+1 
              dscal(k)=0.
           end do
        end do                   

        return
        end

c _____________________________________________________________________________

        subroutine mult1(mi,ntotal4,k,fill,vmac,sinp,dscal)
        implicit real*4 (a-h,o-z)

        real*4 fill,vmac,sinp(*),dscal(*) 

        if(mi.eq.1)then
           prod=fill*vmac
           do j=1,ntotal4
              k=k+1 
              dscal(k)=sinp(j)*prod           
           end do
        end if

        return
        end 

c _____________________________________________________________________________

        FUNCTION atomic_number(symbol)
c        atomic Number returns the atomic number corresponding to a given symbol
        integer atomic_number,iel
        character*2 symbol  

        CHARACTER*2 ATM(92)/'H','HE','LI','BE','B','C','N','O','F','NE',
     *'NA','MG','AL','SI','P','S','CL','AR','K','CA','SC','TI','V','CR',
     *'MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR',
     *'RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN',
     *'SN','SB','TE','I','XE','CS','BA','LA','CE','PR','ND','PM',
     *'SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA','W',
     *'RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN',
     *'FR','RA','AC','TH','PA','U'/
     
        CHARACTER*2 ATM2(92)/'H','He','Li','Be','B','C','N','O','F','Ne',
     *'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     *'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     *'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     *'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm',
     *'Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W',
     *'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
     *'Fr','Ra','Ac','Th','Pa','U'/

       CHARACTER*2 ATM3(92)/'H','He','Li','Be','B','C','N','o','f','ne',
     *'na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr',
     *'mn','fe','co','ni','cu','zn','ga','ge','as','se','br','kr',
     *'rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in',
     *'sn','sb','te','i','xe','cs','ba','la','ce','pr','nd','pm',
     *'sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w',
     *'re','os','ir','pt','au','hg','tl','pb','bi','po','at','rn',
     *'fr','ra','ac','th','pa','u'/
     
       iel=1
       do while(iel .lt. 92 .and. symbol.ne.ATM(iel)
     &   .and. symbol.ne.ATM2(iel) .and. symbol.ne.ATM3(iel) )
         iel=iel+1
       end do
       atomic_number=iel
       
       return
       
      end
