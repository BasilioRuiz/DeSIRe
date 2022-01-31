c DECONV convoluciona vin con una gaussiana de "velocidad" s
c y con el filtro leido del fichero filtro (col. 1 ldo, col2 trans.)
c da el resultado en vin
c isigno (=1 convolucion),(=-1 deconvolucion),(0= return)
c n es el numero de puntos
c el vector a deconvolucionar entra por medio del common vobservado

      subroutine deconv(vin,isigno,nlins,npass,dlamda0s,dlamdas,s)
      
      implicit real*4 (a-h,o-z)
      include 'PARAMETER'

      real*4 vin(*),dlamda0s(*),dlamdas(*)
      real*4 frec,ex,cota,paso,sigma,pi,c,s
      real*4 v1(m1),v2(m1),expo(m1)
      integer npass(*),ifiltro,nlins,isigno,linea
      character*20 msg1,msg2,msg3
      common/ifiltro/ifiltro
      data ivez/0/

      ntot=m1  !antes ponia ntot=kld

	c=3.e5		!c en km/seg luego "s vendra en km/seg"
	pi=3.141590    !dlamda0 vendra dada en a
	cota=26.
        
	if(ivez.eq.0.and.ifiltro.eq.1)then
          ivez=1
          call psf(ntot,nlins,npass,dlamda0s,dlamdas) !da la transf. en ftrans
        end if 
        if(ivez.eq.0.and.ifiltro.eq.2)then
          ivez=1
          call psf_variablePSF(ntot,nlins,npass,dlamda0s,dlamdas)
	end if

	k=0
	if(isigno.eq.0.and.ifiltro.eq.0)return
	if(abs(s).lt.1.e-15 .and.ifiltro.eq.0)return

c descomponemos vin en datos para cada linea: v1
	k1=0
	kk=0
   
c rellenamos hasta 'ntot' con ceros
	do j=1,nlins
	  n2=npass(j)
	  i1=(ntot-n2)/2
	  i2=i1+n2

          if(n2.gt.ntot)then
             write(msg1,*) j
             write(msg2,*) n2
             write(msg3,*) ntot
             call error(KSTOP,'deconv','Line '//trim(adjustl(msg1))//' has'
     &       //         ' more frequencies ('//trim(adjustl(msg2))//') than'
     &       //         ' m1 = '//trim(adjustl(msg3))//'\n'
     &       //         ' Decrease the number of frequencies or change'
     &       //         ' the PARAMETER file')
          end if

	  if(n2.gt.1)then
	    sigma=1.e3*dlamda0s(j)*s/c
	    paso=abs((dlamdas(k1+2)-dlamdas(k1+1)))
	    if(paso.lt.1.e-20)paso=1.
	    do i=1,n2
		k1=k1+1
		v1(i)=vin(k1)		!*1.e-13
	    end do

	    do i=1,i1
	       v2(i)=v1(1)
	    end do

	    do i=i2+1,ntot
	       v2(i)=v1(n2)
	    end do

	    do i3=i1+1,i2
		v2(i3)=v1(i3-i1)
	    end do

	    do i=1,ntot/2
		frec=(i-1)/(paso*ntot)
	        ex=2.*pi*pi*sigma*sigma*frec*frec
		if(ex.gt.cota)ex=cota
	        expo(i)=exp(-1.*isigno*ex)
	    end do
	    do i=ntot/2+1,ntot
		frec=(ntot-(i-1))/(paso*ntot)
	        ex=2.*pi*pi*sigma*sigma*frec*frec
		if(ex.gt.cota)ex=cota
	        expo(i)=exp(-1.*isigno*ex)
	    end do

	    linea=j     
c	    if(n2.gt.1)call deconvsub2(v2,ntot,expo,linea)
	    if(n2.gt.1 .and.ifiltro.ne.1)call deconvsub2(v2,ntot,expo,linea)
	    if(n2.gt.1 .and.ifiltro.eq.1)call deconvsub2(v2,ntot,expo,1)
	    
	    do i=i1+1,i2
		kk=kk+1
		vin(kk)=v2(i)		
	    end do
	  else
	    do i=1,n2
	       kk=kk+1
	       k1=k1+1
	    end do
	  end if
	end do

	return
	end
c_______________________________________________________________________
c psf lee la PSF del fichero filtro se supone que contiene dos columnas
c la primera con la ldo en mA centrada en el maximo
c la segunda la psf (no necesita estar normalizada en area)
c escribe su transformada en ftrans 
c da los valores solo para el muestreo definido por dlamda para  
c cada linea 
c num=numero de puntos del fichero filtro
c ntot=numero de puntos que quiero

	subroutine psf(ntot,nlin,npas,dlamda0,dlamda)

	implicit real*4 (a-h,o-z)
	include 'PARAMETER'

        parameter (m2=2*m1,nmx=m1)
	real*4 dlamda0(*),dlamda(*)
	real x(nmx),y(nmx)
	
	real*4 entrada(m2)
	real*4 ftrans(4*kl,m2)
	real xa(11),ya(11)
	integer npas(*)

        character*100 filtro

	common/filtro/filtro
	common/ftransformada/ftrans

        ngrado=2        !interpolo con parabolas
	n2=int(ngrado/2)

c leemos el fichero filtro
	open(55,file=filtro,status='old',err=100)
	do ii=1,nmx
           read(55,*,end=127,err=101)x(ii),y(ii)
	end do
127	num=ii-1

        if(num.eq.nmx-1)then
           call error(KWARN,'psf','The PSF file is being truncated (it has'
     &     //         ' more than 401 wavelengths)')
        end if


c interpolamos
	k1=0 
        kcount=0  
        ntot2=ntot/2
	do j=1,nlin
	   n=npas(j)
	   paso=abs((dlamda(k1+2)-dlamda(k1+1)))
	   k1=k1+n
	  	   
	   do i=1,ntot
              if(i.le.ntot2)then
                 x0=(i-1)*paso
              else
                 x0=(i-1-ntot)*paso
              end if
        
              entrada(2*i)=0.

	      if(x0.lt.x(1).or.x0.gt.x(num))then
		 entrada(2*i-1)=0.
              else
	         call locate(x,num,x0,jj)  !jj es el indice anterior
	         n3=jj-n2-1
		 if(n3.lt.0)n3=0
                 if(n3+ngrado+1.gt.num)n3=num-ngrado-1
	         do k=1,ngrado+1
		    xa(k)=x(n3+k)
	         end do
	         do k=1,ngrado+1
		    ya(k)=y(n3+k)
	         end do
	         call polint(xa,ya,ngrado+1,x0,salida,dy)
                 entrada(2*i-1)=salida
	      end if
           end do
	   
 	   call four1(entrada,ntot,1)

           area=entrada(1)
	   do i=1,2*ntot
	      ftrans(j,i)=entrada(i)/area
	   end do 
 
	end do !do en lineas
        close(55)      

        return

100     call error(KSTOP,'psf','The file containing the PSF does not exist\n'
     &  //         ' File: '//filtro)

101     call error(KSTOP,'psf','The file containing the PSF is not written'
     &  //         ' properly\n File: '//filtro)

        end
c_______________________________________________________________________
c DECONVSUB2 (incluida en DECONV) multiplica las transformadas de
c Fourier de la gaussiana, de la PSF y del perfil y anti-transforma

c rutina deconvsub2(data,n,datb,isigno)
c isigno=1 convolucion, -1 deconvolucion
c data : vector a convolucionar
c datb : transformada del vector con/deconvolucionador

	subroutine deconvsub2(data,n,datb,linea)
	
        implicit real*4 (a-h,o-z)
	include 'PARAMETER'

	parameter(m2=2*m1)
	real*4 data(*),datb(*),entrada(m2)
	real*4 ftrans(4*kl,m2)
	integer ifiltro,n,linea
	common/ftransformada/ftrans
	common/ifiltro/ifiltro

	if(n.le.1)return
	if(n.ne.m1)then
	   call error(KSTOP,'deconvsub2','Error: n not equal to m1')
	end if
        
        do i=1,n
           entrada(2*i-1)=data(i)
           entrada(2*i)=0.
        end do
 	call four1(entrada,n,1)
	
	k=0
	if (ifiltro.ge.1)then
	   do i=1,n
              preal1=entrada(k+1)*datb(i)
              pimag1=entrada(k+2)*datb(i)
              preal2=ftrans(linea,k+1)
              pimag2=ftrans(linea,k+2)
              entrada(k+1)=preal1*preal2-pimag1*pimag2  !parte real
              entrada(k+2)=preal1*pimag2+preal2*pimag1  !parte imaginaria
              k=k+2
           end do 
	else
           do i=1,n
              entrada(k+1)=entrada(k+1)*datb(i) !parte real
              entrada(k+2)=entrada(k+2)*datb(i) !parte imaginaria
              k=k+2
           end do 
	end if
	   
	call four1(entrada,n,-1)

c normalizamos dividiendo por el numero de puntos
        xn=float(n) 

        do i=1,n
           data(i)=entrada(2*i-1)/xn !parte real
        end do 

	return
	end
	
c********************************************************
c	Subrutina FOUR1 del NUMERICAL RECIPES
	
	subroutine four1(data,nn,isigno)
	implicit real*4 (a-h,o-z)
	real*8 wr,wi,wpr,wpi,wtemp,theta
	dimension data(*)

	n=2*nn
	j=1
	do 11 i=1,n,2

	  if (j.gt.i) then
	    tempr=data(j)
	    tempi=data(j+1)
	    data(j)=data(i)
	    data(j+1)=data(i+1)
	    data(i)=tempr
	    data(i+1)=tempi
	  endif
	  m=n/2
 1	  if ((m.ge.2).and.(j.gt.m)) then
	    j=j-m
	    m=m/2
	    goto 1
	  end if
	  j=j+m
 11	continue

	mmax=2
 2	if (n.gt.mmax) then

	istep=2*mmax
	theta=6.28318530717959d0/(isigno*mmax)
	wpr=-2.d0*dsin(0.5d0*theta)**2
	wpi=dsin(theta)
	wr=1.d0
	wi=0.d0
	do 13 m=1,mmax,2
	  do 12 i=m,n,istep
	    j=i+mmax
	    tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
	    tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
	    data(j)=data(i)-tempr
	    data(j+1)=data(i+1)-tempi
	    data(i)=data(i)+tempr
	    data(i+1)=data(i+1)+tempi
 12	  continue
	  wtemp=wr
	  wr=wr*wpr-wi*wpi+wr
	  wi=wi*wpr+wtemp*wpi+wi
 13	continue
	mmax=istep
	goto 2
	endif

	return
	end
c------------------------------------------------------------------------------------------

c psf_variablePSF lee la PSF del fichero filtro.per 
c con el mismo formato que un perfil
c es decir con 6 columnas de las cuales 
c (col. 1 indice de la linea
c col2=ldo en mA centrada en el maximo
c col3 PSF centrada en el maximo
c se ignoran la columnas 4 a 6)
c escribe su transformada en ftrans 
c da los valores solo para el muestreo definido por dlamda para  
c cada linea 
c num=numero de puntos del fichero filtro
c ntot=numero de puntos que quiero

	subroutine psf_variablePSF(ntot,nlin,npas,dlamda0,dlamda)
	
        implicit real*4 (a-h,o-z)
	include 'PARAMETER'

        parameter (m2=2*m1,nmx=4*kld)
	real*4 dlamda0(*),dlamda(*)
	real x(nmx),y(nmx),xx(nmx),yy(nmx)
	
	real*4 entrada(m2)
	real*4 ftrans(4*kl,m2)
	real xa(11),ya(11)
	integer npas(*)
        integer ntlpsf,nlipsf,npaspsf(kl)
        
        character*100 filtro,psfper
        character*20  msg1,msg2

	common/filtro/filtro
	common/ftransformada/ftrans
	common/variablePSF/psfper

        ngrado=2        !interpolo con parabolas
	n2=int(ngrado/2)

c leemos el fichero filtro
        call leeI(psfper,ntlpsf,nlipsf,npaspsf,xx,yy)!subroutine in leeuveobsindic
                                                     !the PSF files should contain indx,lam, psf
                                                     !ntlpsf total number of lines (no Stokes)
                                                     !nlipsf total number of frecuencies
                                                     !npaspsf array of frecuencies of each line (no Stokes)
        nstokes= nlin/ntlpsf 

        if( nstokes .ne.1 .and. nstokes .ne.2 .and. 
     &      nstokes .ne.3 .and. nstokes .ne.4 )then
           write(msg1,*) ntlpsf
           write(msg2,*) nlin
           call error(KSTOP,'psf_variablePSF','PSF file should have the'
     &     //         ' same number ('//trim(adjustl(msg1))//') of'
     &     //         ' spectral lines\n than the .per file'
     &     //         ' ('//trim(adjustl(msg2))//')')
        end if
        
c interpolating
	k1=0 
        ntot2=ntot/2
        j=0
        do istokes=1,nstokes
        kk=0
	do ji=1,ntlpsf
	   j=j+1
	   num=npaspsf(ji)
	   
	   do ii=1,num
	      kk=kk+1
	      x(ii)=xx(kk)
	      y(ii)=yy(kk)
           end do
	   n=npas(j)
	   paso=abs((dlamda(k1+2)-dlamda(k1+1)))
	   k1=k1+n
	   do i=1,ntot
              if(i.le.ntot2)then
                 x0=(i-1)*paso
              else
                 x0=(i-1-ntot)*paso
              end if
        
              entrada(2*i)=0.

	      if(x0.lt.x(1).or.x0.gt.x(num))then
		 entrada(2*i-1)=0.
              else
	         call locate(x,num,x0,jj)  !jj es el indice anterior
	         n3=jj-n2-1
		 if(n3.lt.0)n3=0
                 if(n3+ngrado+1.gt.num)n3=num-ngrado-1
	         do k=1,ngrado+1
		    xa(k)=x(n3+k)
	         end do
	         do k=1,ngrado+1
		    ya(k)=y(n3+k)
	         end do
	         call polint(xa,ya,ngrado+1,x0,salida,dy)
                 entrada(2*i-1)=salida
	      end if
           end do
	   
 	   call four1(entrada,ntot,1)

           area=entrada(1)
	   do i=1,2*ntot
	      ftrans(j,i)=entrada(i)/area
	   end do 
        end do !do in lines
	end do !do in stokes

	return

	end
c_______________________________________________________________________

