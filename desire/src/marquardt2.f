        subroutine marquardt2(y,sig,ndata,a,mnodos,mfit,
     &                        covar,alpha,chisq,alamda,iamplio,beta)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'
        parameter (kld4=4*kld,nodos=18)

c       dimension y(ndata),sig(ndata),a(*),mnodos(nodos),v(mmx,mmx),w(mmx),
c    &  covar(mfit,mfit),alpha(mfit,mfit),atry(mmx),beta(mfit),da(mmx)
        dimension y(ndata),sig(ndata),a(*),mnodos(*),v(mfitmax,mfitmax),
     &  w(mfitmax),covar(mfitmax,mfitmax),alpha(mfitmax,mfitmax),
     &  atry(mfitmax),beta(mfitmax),da(mfitmax)

        real*4 sigreal(kld4)
        integer iRH1,iRH2
        real*4 rnlte_th
        integer mnodosold(nodos)
        character*100 msg

        common/tol/tol
        common/alamda0/alamda0
        common/ochisq/ochisq
        common/sigrealchi/sigreal,chireal,sumsq
        common/repeticion/factorrep !factor de lambda para evitar repeticiones
        common/nciclos/nciclos  !del principal
c       El common siguiente viene del principal y es para pasarle a marquardt2
c       los indices iniciales y finales de las perturbaciones a gamma y fi
c       aditivas.
        common/ifies/ipa1,ipa11,iga1,ifi11,ipa2,ipa22,iga2,ifi22
        common/ieliminofrec/ielimino  !(del ppal) para el print de la S/N
        common/primerchi/snn,chisn
        common/iRH/iRH1,iRH2  !eq 1 we call RH, 0 don't call RH
        common/thresholdNLTE/rnlte_th

        data it/0/
        xpi=3.14159265

        if(alamda.lt.0.)then
           alamda=alamda0
           call marqcoef2(y,sig,ndata,a,mnodos,mfit,alpha,beta,chisq)
           if(nciclos.lt.1)return
           iRH1=0
           iRH2=0
           ochisq=chisq
           chip=sumsq/float(ndata-ielimino)
           snn=1./sqrt(chip)
           chisn=chireal/float(ndata-ielimino-mfit)
           write(msg,786)0,snn,chisn
           call error(KLITE,'',msg)
        endif

786     format(1x,i3,6x,1pe9.2,1x,e10.3)

c       ::::::::::::::::::::::::hasta aqui solo para la primera iteracion

        if(alamda .gt. 12345.49 .and. alamda .lt. 12345.51)then
           call marqcoef2(y,sig,ndata,a,mnodos,mfit,covar,da,chisq)
           return
        endif

        it=it+1

        do j=1,mfit
           atry(j)=a(j)
           do k=1,mfit
              covar(j,k)=alpha(j,k)
           end do
           covar(j,j)=alpha(j,j)*(1.+alamda)
           da(j)=beta(j)
        end do

        call svdmatriz2(covar,mfit,mnodos,da,tol,v,w)

        datTmax=0

        do j=1,ipa1  !nodos en T atm 1
           atry(j)=a(j)*(1.+da(j))
           if(it .eq.1 .and. abs(da(j)) .gt. datTmax)datTmax=abs(da(j))
        end do
        if(it .eq.1 .and. datTmax .gt. rnlte_th .and. rnlte_th .lt. 10. )iRH1=1
        do j=ipa1+1,ipa11                 !nodos en Pe atm 1
           if(da(j) .gt. .25) da(j)=.25
           if(da(j) .lt. -0.25) da(j)=-0.25
           atry(j)=a(j)*(1.+da(j))        !escala con la Pe en el ultimo nodo
        end do
        do j=ipa11+1,iga1                 !nodos en mic,H,Vz atm 1
           atry(j)=a(j)*(1.+da(j))
        end do
        do j=iga1+1,ifi11                 !nodos en gamma y fi atm 1
           atry(j)=a(j)+da(j)
        end do
        do j=ifi11+1,ipa2                 !resto nodos atm1 y nodos T atm2
           atry(j)=a(j)*(1.+da(j))
        end do
        do j=ipa2+1,ipa22                 !nodos en Pe atm 2
           atry(j)=a(j)*(1.+da(j))        !escala con la Pe en el ultimo nodo
        end do
        do j=ipa22+1,iga2                 !nodos en mic,H,Vz atm 2
           atry(j)=a(j)*(1.+da(j))
        end do
        do j=iga2+1,ifi22                 !nodos en gamma y fi atm 2
           atry(j)=a(j)+da(j)
        end do
        do j=ifi22+1,mfit                 !resto nodos atm2
           atry(j)=a(j)*(1.+da(j))
        end do

        mfitold=mfit
        do j=1,18
           mnodosold(j)=mnodos(j)
        end do
        ipa1old=ipa1
        ipa2old=ipa2
        ipa11old=ipa11
        ipa22old=ipa22
        iga1old=iga1
        iga2old=iga2
        if11old=ifi11
        if22old=ifi22

        call marqcoef2(y,sig,ndata,atry,mnodos,mfit,covar,da,chisq)    

        if(chisq.lt.ochisq)then

           if(abs(factorrep).gt.1.e-3)then
              if(alamda.gt.1.e-4)then
                 alamda=0.1*alamda
              else
                 alamda=0.5*alamda
              end if
           end if
           ochisq=chisq

           do j=1,mfit
              do k=1,mfit
                 alpha(j,k)=covar(j,k)
              end do
              beta(j)=da(j)
              a(j)=atry(j)

              iamplio=1
           end do

        else

           if(alamda.le.1.e-3)then
              alamda=100.*alamda
           else
              if(alamda.lt.1.e3)then
                 alamda=10.*alamda
              else
                 alamda=2.*alamda
              end if
           end if

           chisq=ochisq
           iamplio=0

           mfit=mfitold
           do j=1,18
              mnodos(j)=mnodosold(j)
           end do
           ipa1=ipa1old
           ipa2=ipa2old
           ipa11=ipa11old
           ipa22=ipa22old
           iga1=iga1old
           iga2=iga2old
           if11=ifi11old
           if22=ifi22old

        endif

        return
        end
