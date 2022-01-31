c leeuve2 lee un fichero de datos de parametros de stokes :vobs

        subroutine leeuve2(idis,vobs,ist,ntl,nlin,npas,dl,stok)
        implicit real*4 (a-h,o-z)

        include 'PARAMETER'
        real*4 stok(*),dl(*),si(kld),sq(kld),su(kld),sv(kld)
        real*4 continuohNLTEarr(kld)
        integer npas(*),nlin(*),ist(*),idis,icontinuo_normz
        character vobs*(*)
c       05/05/21 brc: Common to set on/off the continuum normalization in leeuve2
        common/icontnormalz/icontinuo_normz !icontinuo_normz=0 input/output profiles are in SI units
        common/continuosNLTEarr/continuohNLTEarr 
        
        ican=56

        if(idis.eq.0)then
            open(ican,file=vobs,status='old',err=991)
            k=0
            if(icontinuo_normz .eq. 1)then
               do i=1,ntl
                  do j=1,npas(i)
                     k=k+1
                     read(ican,*,err=992)a,dl(k),si(k),sq(k),su(k),sv(k)
                     nlin(i)=nint(a)
                  end do
               end do
            else
                do i=1,ntl
                  do j=1,npas(i)
                     k=k+1
                     con=continuohNLTEarr(k)
                     read(ican,*,err=992)a,dl(k),sii,sqi,suu,svv
                     si(k)=sii/con
                     sq(k)=sqq/con
                     su(k)=suu/con
                     sv(k)=svv/con                     
                     nlin(i)=nint(a)
                  end do
               end do
            end if   
            call sfromiquv(ist,k,si,sq,su,sv,stok)
        else
           open(ican,file=vobs)
           k=0
           do i=1,ntl
              do j=1,npas(i)
                 k=k+1
              end do
           end do
           call iquvfroms(ist,k,si,sq,su,sv,stok)
           k=0
           if(icontinuo_normz .eq. 1)then           
              do i=1,ntl
                 do j=1,npas(i)
                    k=k+1
                    write(ican,993)nlin(i),dl(k),si(k),sq(k),su(k),sv(k)
                 end do
              end do
           else
              do i=1,ntl
                 do j=1,npas(i)
                    k=k+1
                    con=continuohNLTEarr(k)                     
                    write(ican,993)nlin(i),dl(k),si(k)*con,sq(k)*con,su(k)*con,sv(k)*con
                 end do
              end do           
           end if
       end if
       close(ican)
       return

c      Mensajes de error.
991    call error(KSTOP,'leeuve2','The file containing the observed/stray'
     & //         ' light profiles does not exist\n File: '//vobs)

992    call error(KSTOP,'leeuve2','Incorrect format in the file containing'
     & //         ' the observed/stray light profiles\n File: '//vobs)

c      Formato de escritura.
993    format(1x,i5,1x,f11.4,1x,4(e14.7,1x))

       end

c______________________________________________________________________
c sfromiquv guarda en stok los vectores si,sq,su,sv segun ist(4)
c si por ejemplo ist=1,0,1,0 stok contiene si y a continuacion su

	subroutine sfromiquv(ist,n,si,sq,su,sv,stok)
        implicit real*4 (a-h,o-z)

	integer ist(*),n
	real*4 si(*),sq(*),su(*),sv(*),stok(*)

	ks=0
	if(ist(1).eq.1)then
	   do i=1,n
	      ks=ks+1
	      stok(ks)=si(i)
	   end do
	end if
	if(ist(2).eq.1)then
	   do i=1,n
	      ks=ks+1
	      stok(ks)=sq(i)
	   end do
	end if
	if(ist(3).eq.1)then
	   do i=1,n
	      ks=ks+1
	      stok(ks)=su(i)
	   end do
	end if
	if(ist(4).eq.1)then
	   do i=1,n
	      ks=ks+1
	      stok(ks)=sv(i)
	   end do
	end if

	return
	end
c______________________________________________________________________
c iquvfroms divide stok en los vectores si,sq,su,sv segun ist(4)
c si por ejemplo ist=1,0,1,0 sq=0 sv=0

	subroutine iquvfroms(ist,n,si,sq,su,sv,stok)
        implicit real*4 (a-h,o-z)

	integer ist(*),n
	real*4 si(*),sq(*),su(*),sv(*),stok(*)

	        ks=0
	        if(ist(1).eq.1)then
	           do i=1,n
	              ks=ks+1
	              si(i)=stok(ks)
	           end do
	        else
	           do i=1,n
	              si(i)=0.
	           end do
	        end if
	        if(ist(2).eq.1)then
	           do i=1,n
	              ks=ks+1
	              sq(i)=stok(ks)
	           end do
	        else
	           do i=1,n
	              sq(i)=0.
	           end do
	        end if
	        if(ist(3).eq.1)then
	           do i=1,n
	              ks=ks+1
	              su(i)=stok(ks)
	           end do
	        else
	           do i=1,n
	              su(i)=0.
	           end do
	        end if
	        if(ist(4).eq.1)then
	           do i=1,n
	              ks=ks+1
	              sv(i)=stok(ks)
	           end do
	        else
	           do i=1,n
	              sv(i)=0.
	           end do
	        end if
	
	return
	end
