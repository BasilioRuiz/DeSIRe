c PgmufromPeT_isub calcula Pg,pesomedio,kappa y rho desde Pe y T evaluados a una altura dada
c              siendo kappa coeficiente absor. del continuo a 5000 por gramo
c Basilio Ruiz 30/08/2016
c _______________________________________________________________
	subroutine PgmufromPeT_isub(T,Pe,Pg,pesomedio,rho,kappa)
    
	implicit real*4 (a-h,o-z)

	parameter (cgases=83145100.)
	real*4 T,Pg,Pe,pesomedio,rho,kappa
    
	real*4 wgt,abu,ei1,ei2,pp(10),d1(10)
	real*4 kac,asum,pmsum
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
        common/precisoitera/precitera 
        common/nmaxitera/nmaxitera
        
        precitera=1.e-5
        nmaxitera=250
        
        g=2.7414e+4
        avog=6.023e23
                      
        do i=1,10
           d1(i)=0
        end do
        
c calculamos el peso molecular medio pesomedio
	pmusum=0.0
        asum=0.0
	do i=1,92
  	   ii=i
	   call neldatb(ii,0.,wgt,abu,ei1,ei2)
	   pmusum=pmusum+wgt*abu
           asum=asum+abu
	end do
c	print*,'pmusum asum=',pmusum,asum

c	Pe=Pg*1.e-4 !inicializamos 
c	call pefrompg11(T,Pg,Pe)
c	print*,'en PemufromPgT_isub.f T Pg Pe=',T,Pg,Pe

	call gasc(T,Pe,Pg,pp)
	
c	print*,'pp=',pp
	
	call kappach(5.e-5,T,Pe,pp,d1,d1,kac,d2,d2)
	pesomedio=pmusum/(asum+pp(8)) !peso molec. medio
	
c	print*,'pesomedio=',pesomedio
        kappa=kac*avog/pmusum         !coeficiente absorcion continuo a 5000 por gramo
        rho=pesomedio*Pg/T/cgases
        
        
           
       return
       end 

