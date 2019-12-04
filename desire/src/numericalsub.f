c numerical
c calls RH to evaluate RF numerically
c Basilio 19/06/2019
c ________________________________________________________________________
	subroutine numericalsub(label_ID_model,RH_model,RH_magneticfield,natmos,atmosLG,
     +mnod,stok_RHinput,vmac1,ifiltro,ntotal,ndata,rt)
     	
     	implicit real*4 (a-h,o-z)
     	include 'PARAMETER'   !por kt,kn,kl,kld
	parameter (kt8=8*kt+2,kt16=16*kt+5,kt11=11*kt+2,kld4=4*kld,kt12=11*kt+3) 
	
        integer natmos,mnod(*)
        integer nlam_LTE
        real*8 bol
        real*4 stok_RHinput(kld4),stok_RH(kld4),rt(*)   
        real*4 dlamda0(kl),dlamda(kld)
        character*100 label_ID_model,RH_model,RH_magneticfield !BRC-RH Jun 20 2017
        integer istatus,system

        real*4 tau(kt),T(kt),Pe(kt),Pg(kt),z(kt),ro(kt)
        character*100 vprint,ruta,ruta2,ruta3,filewavename
        
c para la atmosfera
        real*4 atmosLG(kt12),atmosLGpert(kt12)
        integer nlin(kl),npas(kl),ist(4),nble(kl),nfrecuencies
        data ivez/0/
        
        common/vprint/vprint
        common/ruta/ruta,ruta2,ruta3,filewavename
        common/istatus12/istatus1,istatus2
        common/numero_LTE/nlam_LTE
       	common/responde2/ist,ntau,ntl,nlin,npas,nble
       	common/ldeo/dlamda,dlamda0


        ivez=ivez+1
        bol=1.3806488d-16  !erg/s
        istatus=0 
        
        print*,'entering numerical to evaluate RF to T numerically'
 
        deltaT=1. ! a perturbation of 100 K at nodes     
        nfrecuencies=(ist(1)+ist(2)+ist(3)+ist(4))*nlam_LTE

c We perturb the atmospheres at the nodes and run RH for each perturbed model      
        if(mnod(1) .eq. 1)then  !constant perturbation
           Tmean=0
           do k=1,ntau
              tau1=atmosLG(k)      !logtau
              atmosLGpert(k)=tau1  !logtau
              tau(k)=tau1
              T1=atmosLG(ntau+k)+deltaT
              atmosLGpert(ntau+k)=T1
              T(k)=T1
              Pe(k)=atmosLG(2*ntau+k)
              Pg(k)=atmosLG(9*ntau+2+k)
              Tmean=Tmean+T1
           end do
           Tmean=Tmean/deltaT/float(ntau)
           call equisubmu(ntau,tau,T,Pe,Pg,Z,ro)
           
	   do k=1,ntau
	      atmosLGpert(2*ntau+k)=Pe(k)
	      do kk=3,7 
	        atmosLGpert(kk*ntau+k)=atmosLG(kk*ntau+k)
	      end do  
	      atmosLGpert(8*ntau+2+k)=Z(k)
	      atmosLGpert(9*ntau+2+k)=Pg(k)
	      atmosLGpert(10*ntau+2+k)=ro(k)
	   end do 
	   atmosLGpert(11*ntau+2+3)=atmosLG(11*ntau+2+3)
	   
           call write_atmos_RH(label_ID_model,RH_model,RH_magneticfield,atmosLGpert,ntau)
c           istatus = system("rm PRD*.dat")
           print*,'$ CALLING RH from  ',ruta    !(1:40)
           istatus1 =system(ruta)
           print*,'$ CALLING solveray from  ',ruta2 
           istatus2 =system(ruta2)
      	   call read_spectrum(stok_RH)
      	   
      	   if(vmac1.gt.0. .or. ifiltro.ge.1)then
	      do j=ntotal+1,ndata
                 stok_RH(j)=0.
	      end do
	      call deconv(stok_RH,1,ntl,npas,dlamda0,dlamda,vmac1)
            end if 
      	    
      	   do ifrec=1,nfrecuencies
      	     rt(ifrec)=(stok_RH(ifrec)-stok_RHinput(ifrec))*Tmean
c      	     print*,'numericalsub ',ifrec,stok_RH(ifrec),stok_RHinput(ifrec),rt(ifrec)
      	   end do
        endif
        if(mnod(1) .gt. 1)then
          stop
        endif
        
        return 
       
        end
