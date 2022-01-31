c DECONV2 convoluciona vin con la derivada de una gaussiana de
c "velocidad" s con respecto a s 
c  y con el filtro leido del fichero "filter.dat"(col. 1 ldo, col2 trans.)
c da el resultado en vin
c isigno (=1 convolucion),(=-1 deconvolucion),(0= return)
c n es el numero de puntos

      subroutine deconv2(vin,isigno,nlins,npass,dlamda0s,dlamdas,s)

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

	c=3.e5	       !c en km/seg luego "s vendra en km/seg"
	pi=3.141590    !dlamda0 vendra dada en a
	cota=26.
        
	if(isigno.eq.0)return

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
             call error(KSTOP,'deconv2','Line '//trim(adjustl(msg1))//' has'
     &       //         ' more wavelengths ('//trim(adjustl(msg2))//') than'
     &       //         ' m1 = '//trim(adjustl(msg3))//'\n'
     &       //         ' Decrease the number of wavelengths or change'
     &       //         ' the PARAMETER file')
          end if

	  if(n2.gt.1)then
	    sigma=1.e3*dlamda0s(j)*s/c
	    paso=abs((dlamdas(k1+2)-dlamdas(k1+1)))

	    if(paso.lt.1.e-20)paso=1.
	    do i=1,npass(j)
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
		expo(i)=-expo(i)*2.*ex/s   !esto es por la derivada
	    end do
	    do i=ntot/2+1,ntot
		frec=(ntot-(i-1))/(paso*ntot)
	        ex=2.*pi*pi*sigma*sigma*frec*frec
		if(ex.gt.cota)ex=cota
	        expo(i)=exp(-1.*isigno*ex)
		expo(i)=-expo(i)*2.*ex/s   !esto es por la derivada
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

