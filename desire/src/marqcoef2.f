         subroutine marqcoef2(y,sig,ndata,a,m,mfit,alpha,beta,chisq)

         implicit real*4 (a-h,o-z)
         include 'PARAMETER'  !por kn,kl,kld,mfitmax
         parameter (kn16=16*kn+4)         !
         parameter (mmax=16*kn+2,kld4=4*kld,kldn=mmax*kld4)        
         dimension y(ndata),sig(ndata),alpha(mfitmax,mfitmax),beta(*)
         dimension dyda(mfitmax),a(*),ymod(kld4),ymodobs(kld4),dvcal(kldn)
         dimension ymodobs_RH(kld4),diff_NLTE_LTE(kld4)      !,ynew(kld4)
         real*4 scal_RH(kld4)
         integer m(*),npos(kld4),numberfrec(4),iratio(4)
         real*4 sigreal(kld4),stokesratio(kld4),sigescala(kld4)
         real*4 deltastokes_i(kld4)
         real*4 dcastigo(mfitmax),ddcastigo(mfitmax)
         integer iRH1,iRH2            !,ileoNLTE
        
         common/uvesalida/stokesratio,ymodobs_RH
         common/sigrealchi/sigreal,chireal,sumsq
         common/nciclos/nciclos   !del principal
         common/derivchi/dvcal  ! deriv. chi^2 para los errores
         common/primera2/ntotal,ntotal4,ists
         common/posiciones/npos 
         common/numberfrec/numberfrec,iratio     !para marqcoef2
         common/ndata/ndatosobs !para fperfil2 y automatico
         common/ifies/ipa1,ipa11,iga1,ifi11,ipa2,ipa22,iga2,ifi22
         common/nspl_lin/nspl_lin !para seleccionar penalty2 (4,5); penalty1(2 ,3) o no(0,1)
         common/iRH/iRH1,iRH2 
         common/storing_diff/diff_NLTE_LTE
         common/storing_NLTE/scal_RH
         common/escalachi/escalachi,sigescala,chinew
         common/deltastokes_i/deltastokes_i
         data idiff/0/
         data ivez/0/
        
         ivez=ivez+1
         if(ivez .eq. 1)escalachi=0
         chisq=0.
         chireal=0.
         sumsq=0.
         chinew=0.

         ndatosobs=ndata
         icastigo=0
         if(nspl_lin .eq. 2 .or. nspl_lin .eq. 3)icastigo=1
         if(nspl_lin .eq. 4 .or. nspl_lin .eq. 5)icastigo=2

         call fperfil2(m,a,ymod,dvcal,scal_RH)

         mfit=0
         if(m(1).gt.0)mfit=mfit+m(1)
         ipa1=mfit       !ipa1 es el indice anterior a la gamma comp. 2 
         if(m(2).gt.0)mfit=mfit+m(2)
         ipa11=mfit
         do i=3,5
           if(m(i).gt.0)mfit=mfit+m(i)
         end do
         iga1=mfit    !iga1 es el indice anterior a la gamma comp. 1
         if(m(6).gt.0)mfit=mfit+m(6)
         if(m(7).gt.0)mfit=mfit+m(7)
         ifi11=mfit   !ifi11 es el indice ultimo de la fi comp. 1
         do i=8,9
           if(m(i).gt.0)mfit=mfit+m(i)
         end do
         ipa2=mfit         !ipa2 es el indice anterior a la gamma comp. 2
         if(m(10).gt.0)mfit=mfit+m(10)
         ipa22=mfit
         do i=11,13
           if(m(i).gt.0)mfit=mfit+m(i)
         end do	
         iga2=mfit    !iga2 es el indice anterior a la gamma comp. 2
         if(m(14).gt.0)mfit=mfit+m(14)
         if(m(15).gt.0)mfit=mfit+m(15)
         ifi22=mfit   !ifi22 es el indice ultimo de la fi comp. 1
         do i=16,18
           if(m(i).gt.0)mfit=mfit+m(i)
         end do
         do j=1,mfit
            do  k=1,j
               alpha(j,k)=0.
            end do
            beta(j)=0.
         end do 

c =======================================================================
c escribimos la FR
      if(nciclos .lt. 0)then
        i=0
        iratio(1)=0 !no permitimos escribir la FR de I/I 
        do ii=1,4
           do iii=1,numberfrec(ii)
              i=i+1
              do j=1,mfit
                 jj=(j-1)*ntotal4+i
                 jji=(j-1)*ntotal4+iii
                 if(iratio(ii).eq.1)dvcal(jj)=
     &             ( dvcal(jj) - dvcal(jji)*ymod(i)/ymod(iii) ) /ymod(iii)
              end do
           end do
        end do
        call escribeFR(m,ndata,ntotal4,npos,dvcal)
      end if
c =======================================================================
c calculamos el castigo (y sus derivadas) por picos en las estratificaciones
      if(icastigo.eq.1) call penalty(mfit,m,a,castigo,dcastigo,ddcastigo)
      if(icastigo.eq.2) call penalty2(mfit,m,a,castigo,dcastigo)
      wcast=1.  !pesocastigo

      do i=1,ndata
         ymodobs(i)=ymod(npos(i))        !synthetic profile with SIR
      end do
      
       if( iRH1 .eq. 1 .or.  iRH2 .eq. 1 )then
         idiff=1
         do i=1,ndata
             ymodobs_RH(i)=scal_RH(npos(i))  !synthetic profile with RH
             diff_NLTE_LTE(i)=(ymodobs_RH(i)-ymodobs(i)) 
         end do
      end if
      
      if( idiff .eq. 1) then !in anytime NLTE has been evaluated
         do i=1,ndata
            ymodobs(i)=ymodobs(i)+diff_NLTE_LTE(i) !synthetic profile with SIR +new/old difference NLTE
         end do
      end if
      
      i=0
      do ii=1,4
         do iii=1,numberfrec(ii)
            i=i+1
            sss=sig(i)
            do j=1,mfit
              jj=(j-1)*ntotal4+npos(i)
              jji=(j-1)*ntotal4+npos(iii)
              if(iratio(ii).eq.1)then
                 dyda(j)=(dvcal(jj)-dvcal(jji)*ymodobs(i)/ymodobs(iii))/ymodobs(iii)
              else
                 dyda(j)=dvcal(jj)
              end if
            end do
            sig2i=1./sss/sss
            sig2reali=1./(sigreal(i)*sigreal(i))  
            if(iratio(ii).ne.1)stokesratio(i)=ymodobs(i)
            if(iratio(ii).eq.1)stokesratio(i)=ymodobs(i)/ymodobs(iii)
            dy=y(i)-stokesratio(i)
            do j=1,mfit
               wt=dyda(j)*sig2i
               do k=1,j
                 alpha(j,k)=alpha(j,k)+wt*dyda(k)
               end do
	           beta(j)=beta(j)+dy*wt
            end do
            if(sss.lt.1.e14)then 
               chisq=chisq+dy*dy*sig2i
               chireal=chireal+dy*dy*sig2reali 
               if(ivez .eq. 1)escalachi=escalachi+sigescala(i)*sigescala(i)*sig2reali 
               sumsq=sumsq+dy*dy
               if(ii.eq.1)chinew=chinew+dy*dy/y(i) 
               if(ii.ne.1)then
                  cota=abs(y(i))
                  if(cota .lt. sigescala(i))cota=sigescala(i)
                  chinew=chinew+dy*dy/cota
               end if
               if(abs(dy) .gt. sigescala(i))then
                  deltastokes_i(i)=abs(dy)
               else
                  deltastokes_i(i)=sigescala(i)
               end if   
            end if   
         end do
      end do
      
c------------------------------castigo 1------------------------------------------------
        if(icastigo .eq.1)then
           castigo=1.+castigo*wcast   !multiplicativamente

           do j=1,mfit
	          alpha(j,1)=alpha(j,1)*castigo+beta(j)*dcastigo(1)*wcast+
     &             beta(1)*dcastigo(j)*wcast+chisq*ddcastigo(j)*wcast
	          do k=2,j
	             alpha(j,k)=alpha(j,k)*castigo+beta(j)*dcastigo(k)*wcast
     &                                     +beta(k)*dcastigo(j)*wcast
              end do
           end do
           do j=1,mfit
              beta(j)=beta(j)*castigo+chisq*dcastigo(j)*wcast
           end do
           chisq=chisq*castigo
           chireal=chireal*castigo
        end if
c------------------------------castigo 2------------------------------------------------
        if(icastigo .eq.2)then
           castigo=1.+castigo*wcast   !multiplicativamente

           do j=1,mfit
	         do k=1,j
	            alpha(j,k)=alpha(j,k)*castigo+beta(j)*dcastigo(k)*wcast+
     &    beta(k)*dcastigo(j)*wcast+chisq*dcastigo(j)*dcastigo(k)*wcast
             end do
          end do
          do j=1,mfit
              beta(j)=beta(j)*castigo+chisq*dcastigo(j)*wcast
          end do
          chisq=chisq*castigo
          chireal=chireal*castigo
        end if
c---------------------------------------------------------------------------------------

        do  j=2,mfit
           do k=1,j-1
              alpha(k,j)=alpha(j,k)
           end do
        end do

      return
      end
