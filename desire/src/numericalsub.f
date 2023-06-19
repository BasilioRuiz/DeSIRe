c Numerical Response Functions
c it is call from fperfil2 if nodos >1000 --> then nodos=nodos-1000
c calls RH to evaluate RF numerically
c numericalsub evaluates RF for the Stokes case
c numericalsub_con evaluates RF to t, vz & mic for the nonStokes case
c Basilio 15/12/2022
c ________________________________________________________________________

        subroutine numericalsub(natmos,atmosLG,mnod,stok_RHinput,vmac,ifiltro,ntotal,ndata,rt,rp,rh,rvz,rg,rf,rmic)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        parameter (kt8=8*kt+2,kt16=16*kt+5,kt11=11*kt+2,kld4=4*kld,kt12=11*kt+3) 
        parameter (kl4=4*kl)    !numero maximo de lineas
        parameter (aln10=2.3025851, bol=1.3806488d-16)

        integer natmos,mnod(*)
        integer nlam_LTE,numerical(18)
        integer ntotal_lines
        real*4  stok_RHinput(kld4),stok_RH(kld4)  
        real*4  rt(*),rp(*),rh(*),rvz(*),rg(*),rf(*),rmic(*)        
        real*4  dlamda0(kl),dlamda0s(kl4),dlamda(kld),dlamdas(kld4)
        real*4  delta,deltaT,deltavz,deltamic
        character*100 label_ID_model,RH_model,RH_magneticfield
        character*40  msg1,msg2
        
c para la atmosfera
        real*4  atmosLG(*),atmosLGpert(kt11+1)
        real*4  atmosLGr(kt11+1),atmosLGrpert(kt11+1)
        real*8  beta(kl,kt)
        integer nlin(kl),npas(kl),ist(4),nble(kl),nfrecuencies
        integer ntls,nlins(kl4),npass(kl4)
        integer atomic_number,atom_arr(kl),istage_arr(kl)
        real*4  alfa_arr(kl),sigma_arr(kl),wave_arr(kl)

        common/istatus12/istatus1,istatus2
        common/numero_LTE/nlam_LTE
        common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/smalla/ntls,nlins,npass
        common/ldeo/dlamda,dlamda0
        common/smalla1/dlamdas
        common/smalla0/dlamda0s  
        common/numerical/numerical
        common/RHnames/label_ID_model,RH_model,RH_magneticfield
        common/brklm/ntotal_lines,atom_arr,istage_arr,alfa_arr,sigma_arr,wave_arr
        common/ndata/ndatosobs !viene de nmarqcoef2: numero de frecuencias incluyendo Stokes activos
        
        call comprime3(ntau,mnod,atmosLG,atmosLGr) !atmosLGr es la atmosfera en los nodos
        isum=0
        if(natmos .eq. 1)then
           do inodo=1,8
              isum=isum+mnod(inodo)
           end do 
        end if   
        mtotalnodos=isum 
c        do i=1, mtotalnodos              ! initialzed NO perturbation
c            atmosLGrpert(i)=atmosLGr(i)  ! atmosphere at nodes + perturbation (multiplicative)
c        end do
            
        deltaT=4.    ! T perturbation (K) at nodes  
        deltaPe=1.e-3 ! Factor of Pe perturbation at nodes
        deltamic=1.e4 ! micro perturbation (cm/s) at nodes 
        deltaB=5.    ! B perturbation (G)
        deltavz=5.e3  ! Vz perturbation (cm/s) at nodes 
        deltagamma=2.5 ! inclination perturbation (degrees) at nodes 
        deltaphi=5.  ! azimuth perturbation (degrees) at nodes
        nfrecuencies=(ist(1)+ist(2)+ist(3)+ist(4))*nlam_LTE
        inodosum=0
        numer=0
        do ii=1,7
           numer=numer+numerical(ii)
        end do      
        do i=1,8
           if(numerical(i).eq.1)then
              write(msg1, '(a)') 'Evaluating RF to'
              write(msg2, '(a,i2,a)') 'numerically at ',mnod(i),' nodes'
              if(i.eq.1)call error(KLINE,'',trim(msg1)//' T '//trim(msg2))
              if(i.eq.2)call error(KLINE,'',trim(msg1)//' Pe '//trim(msg2))
              if(i.eq.3)call error(KLINE,'',trim(msg1)//' micro '//trim(msg2))
              if(i.eq.4)call error(KLINE,'',trim(msg1)//' magnetic field strength '//trim(msg2))
              if(i.eq.5)call error(KLINE,'',trim(msg1)//' Vz '//trim(msg2))
              if(i.eq.6)call error(KLINE,'',trim(msg1)//' inclination '//trim(msg2))
              if(i.eq.7)call error(KLINE,'',trim(msg1)//' azimuth '//trim(msg2))

              do inodo=1,mnod(i) !do in nodes to introduce aditive perturbation in T 
                 do k=1, mtotalnodos                 !reset 
                     atmosLGrpert(k)=atmosLGr(k)     !atmosphere at nodes + perturbation (multiplicative)
                 end do
                 do itau=1,8*ntau+3
                    atmosLGpert(itau)=atmosLG(itau)  !reset 
                 end do
                 if(i .eq. 1 )delta=deltaT
                 if(i .eq. 2 )delta=deltaPe
                 if(i .eq. 3 )delta=deltamic
                 if(i .eq. 4 )delta=deltaB
                 if(i .eq. 5 )delta=deltavz
                 if(i .eq. 6 )delta=deltagamma
                 if(i .eq. 7 )delta=deltaphi
                 if(numer .eq. 0)then
                    call error(KSTOP,'numericalsub','No free parameters with numerical RF')
                 end if   
                 atmosnodo= atmosLGr(inodosum+inodo)
                 atmosnodo_d=atmosnodo/delta 
                 if(i .ne. 2)then
                     atmosLGrpert(inodosum+inodo)=atmosnodo+delta  !perturbation at nodes
                 else
                     atmosLGrpert(inodosum+inodo)=atmosnodo*(1+delta)  !perturbation Pe at nodes
                 end if
                 call amp3(ntau,mnod,atmosLGpert,atmosLGrpert)                 !perturbation at every tau-grid
c                 atmosLGrpert(inodosum+inodo)=atmosnodo         !reset of perturbation at nodes   
c calling RH to synthetize the profile, after that we substract form non-perturbated one to calculate RF 
                 call departures(label_ID_model,RH_model,RH_magneticfield
     &             ,natmos,atmosLGpert,ntau,ntotal_lines,beta,beta,stok_RH) 
  
                 if(vmac.gt.0. .or. ifiltro.ge.1)then
                    call deconv(stok_RH,1,ntls,npass,dlamda0s,dlamdas,vmac)
                 end if 

                 if (i .eq. 1 )then
                    do ifrec=1,nfrecuencies
                       iktt=(inodo-1)*nfrecuencies+ifrec
                       rt(iktt)=(stok_RH(ifrec)-stok_RHinput(ifrec))*atmosnodo_d 
                    end do
                 end if 
                 
                 if (i .eq. 2 )then
                    do ifrec=1,nfrecuencies
                       iktt=(inodo-1)*nfrecuencies+ifrec
                       rp(iktt)=(stok_RH(ifrec)-stok_RHinput(ifrec))/delta  
                    end do
                 end if
                 if (i .eq. 3 )then
                    do ifrec=1,nfrecuencies
                       iktt=(inodo-1)*nfrecuencies+ifrec
                       rmic(iktt)=(stok_RH(ifrec)-stok_RHinput(ifrec))*atmosnodo_d 
                    end do
                 end if
                 if (i .eq. 4)then
                    do ifrec=1,nfrecuencies
                       iktt=(inodo-1)*nfrecuencies+ifrec
                       rh(iktt)=(stok_RH(ifrec)-stok_RHinput(ifrec))*atmosnodo_d
                    end do
                 end if   
                 if (i .eq. 5)then
                    do ifrec=1,nfrecuencies
                       iktt=(inodo-1)*nfrecuencies+ifrec
                       rvz(iktt)=(stok_RH(ifrec)-stok_RHinput(ifrec))*atmosnodo_d
                    end do
                 end if 
                 if (i .eq. 6)then
                    do ifrec=1,nfrecuencies
                       iktt=(inodo-1)*nfrecuencies+ifrec
                       rg(iktt)=(stok_RH(ifrec)-stok_RHinput(ifrec))/delta !aditiva
                    end do
                 end if   
                 if (i .eq. 7)then
                    do ifrec=1,nfrecuencies
                       iktt=(inodo-1)*nfrecuencies+ifrec
                       rf(iktt)=(stok_RH(ifrec)-stok_RHinput(ifrec))/delta  !aditiva
                    end do
                 end if  
              end do
           end if
           inodosum=inodosum+mnod(i)
        end do  

        return 
       
        end

c ________________________________________________________________________

        subroutine numericalsub_con(natmos,atmosLG,mnod,stok_RHinput,vmac,ifiltro,ntotal,ndata,rt,rvz,rmic)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        parameter (kt8=8*kt+2,kt16=16*kt+5,kt11=11*kt+2,kld4=4*kld,kt12=11*kt+3) 
        parameter (aln10=2.3025851, bol=1.3806488d-16)

        integer natmos,mnod(*)
        integer nlam_LTE,numerical(18)
        integer ntotal_lines
        real*4  stok_RHinput(kld4),stok_RH(kld4),rt(*),rvz(*),rmic(*)   
        real*4  dlamda0(kl),dlamda(kld)
        real*4  delta,deltaT,deltavz,deltamic
        character*100 label_ID_model,RH_model,RH_magneticfield
        character*40  msg1,msg2
        
c para la atmosfera
        real*4  atmosLG(*),atmosLGpert(kt16)
        real*4  atmosLGr(kt11+3),atmosLGrpert(kt11+3)
        real*8  beta(kl,kt)
        integer nlin(kl),npas(kl),ist(4),nble(kl),nfrecuencies
        integer atomic_number,atom_arr(kl),istage_arr(kl)
        real*4  alfa_arr(kl),sigma_arr(kl),wave_arr(kl)

        common/istatus12/istatus1,istatus2
        common/numero_LTE/nlam_LTE
        common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/ldeo/dlamda,dlamda0
        common/numerical/numerical
        common/RHnames/label_ID_model,RH_model,RH_magneticfield
        common/brklm/ntotal_lines,atom_arr,istage_arr,alfa_arr,sigma_arr,wave_arr

        
        call comprime3(ntau,mnod,atmosLG,atmosLGr)
        isum=0
        if(natmos .eq. 1)then
           do inodo=1,8
              isum=isum+mnod(inodo)
           end do 
        end if   
        mtotalnodos=isum 
        do i=1, mtotalnodos
            atmosLGrpert(i)=atmosLGr(i)  
        end do
            
        deltaT=10.    ! T perturbation (K) at nodes  
        deltamic=1.e4 ! micro perturbation (cm/s) at nodes  
        deltavz=1.e4  ! Vz perturbation (cm/s) at nodes 
        nfrecuencies=(ist(1)+ist(2)+ist(3)+ist(4))*nlam_LTE
        inodosum=0
        numer=numerical(1)+numerical(3)+numerical(5)       
        do i=1,8
           if(numerical(i).eq.1)then
              write(msg1, '(a)') 'Evaluating RF to'
              write(msg2, '(a,i2,a)') 'numerically at ',mnod(i),' nodes'
              if(i.eq.1)call error(KLINE,'',trim(msg1)//' T '//trim(msg2))
              if(i.eq.3)call error(KLINE,'',trim(msg1)//' micro '//trim(msg2))
              if(i.eq.5)call error(KLINE,'',trim(msg1)//' Vz '//trim(msg2))

              do inodo=1,mnod(i) !do in nodes to introduce aditive perturbation in T 
                 do itau=1,8*ntau+3
                    atmosLGpert(itau)=atmosLG(itau)
                 end do
                 if(i .eq. 1 )delta=deltaT
                 if(i .eq. 3 )delta=deltamic
                 if(i .eq. 5 )delta=deltavz
                 if(numer .eq. 0)then
                    call error(KSTOP,'numericalsub_con',
     &                 'In non-magnetic numerical RF available only for T, micro or Vz')
                 end if  
                 atmosnodo= atmosLGr(inodosum+inodo)
                 atmosnodo_d=atmosnodo/delta 
                 atmosLGrpert(inodosum+inodo)=atmosnodo+delta  !perturbation at nodes
                 call amp3(ntau,mnod,atmosLGpert,atmosLGrpert)                 !perturbation at every tau-grid
                 atmosLGrpert(inodosum+inodo)=atmosnodo         !reset of perturbation at nodes   
c calling RH to synthetize the profile, after that we substract form non-perturbated one to calculate RF 
                 call departures(label_ID_model,RH_model,RH_magneticfield
     &             ,natmos,atmosLGpert,ntau,ntotal_lines,beta,beta,stok_RH) 
                 if(vmac.gt.0. .or. ifiltro.ge.1)then
                     do j=ntotal+1,ndata
                       stok_RH(j)=0.
                    end do
                    call deconv(stok_RH,1,ntl,npas,dlamda0,dlamda,vmac)
                 end if 
                 if (i .eq. 1 )then
                    do ifrec=1,nfrecuencies
                       iktt=(inodo-1)*nfrecuencies+ifrec
                       rt(iktt)=(stok_RH(ifrec)-stok_RHinput(ifrec))*atmosnodo_d 
                    end do
                 end if 
                 if (i .eq. 3 )then
                    do ifrec=1,nfrecuencies
                       iktt=(inodo-1)*nfrecuencies+ifrec
                       rmic(iktt)=(stok_RH(ifrec)-stok_RHinput(ifrec))*atmosnodo_d 
                    end do
                 end if
                 if (i .eq. 5)then
                    do ifrec=1,nfrecuencies
                       iktt=(inodo-1)*nfrecuencies+ifrec
                       rvz(iktt)=(stok_RH(ifrec)-stok_RHinput(ifrec))*atmosnodo_d
                    end do
                 end if    
              end do
           end if
           inodosum=inodosum+mnod(i)
        end do  
 
        return 
       
        end

