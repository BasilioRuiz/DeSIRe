c************************************************
		program PemufromPgT_i
c************************************************

c PemufromPgT_i calcula Pe,pesomedio,kappa y rho desde Pg y T evaluados a una altura dada
c              siendo kappa coeficiente absor. del continuo a 5000 por gramo
c             escribe como salida el fichero output_PemufromPgT_i
c Basilio Ruiz 30/08/2016
c _______________________________________________________________


        real*4 T,Pg,Pe,pesomedio,rho,kappa
        real*4 thermo(14)
        character*100 fichabun
        common/fichabun/fichabun
        
c        fichabun='THEVENIN'
              
        print*,'write abundances filename abundfile, T and Pg'
        read*,fichabun,T,Pg
        
        call PemufromPgT_isub(T,Pg,Pe,pesomedio,rho,kappa)
        call thermosub(T,Pe,thermo)
        
c        print*,'Pe pesomedio rho kappa',Pe,pesomedio,rho,kappa
c        print*,'ne nH nH+  ',thermo(5:7)
c        print*,'grd_ion alpha delta  ',thermo(8:10)
c        print*,'c_v c_p c_s nabl_ad ',thermo(11:14)
        
        
	open(2,file='output_PemufromPgT_i')
	write(2,1002) 'Pe','pesomed','rho','kappa','  ne   ','  nH   ','   nH+  ',
     &	'   grd_ion','   alpha ','   delta ','   c_v   ','  c_p   ',
     &  '   c_s  ','  nabl_ad'
        
        write(2,1001)Pe,pesomedio,rho,kappa,thermo(5:14)
        
        close(2)
        
1001     format(1x,14(1pe10.3,1x))
1002     format(1x,14(a10,1x))	

        end 

