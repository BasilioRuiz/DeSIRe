c spectra
c reads all the spectral information and make some calculations
c sends the information via common to blends2, blendscon2, departures, fperfil2
c input via common: responde2
c ________________________________________________________________________

         subroutine spectra( )

         implicit real*4 (a-h,o-z)
         include 'PARAMETER'

         parameter (kt8=8*kt+2,kt16=16*kt+5,kt11=11*kt+2,kt12=11*kt+3)
         parameter (kld4=4*kld)

c        para la malla
         integer nlin(kl),npas(kl),ist(4)
         integer npos(kld4),ntl,ntls,nble(kl)
         real*4  dlamda0(kl),dlamda(kld)

c        departure coefficients
         character atom_nlte(kl)*2
         integer   nlow_i(kl),nup_i(kl),linea_nlte(kl),linea_all_nlte(kl)
         integer   ntotal_lines
         integer   nlow,nup
         integer   atomic_number,atom_arr(kl),istage_arr(kl)
         real*4    alfa_arr(kl),sigma_arr(kl),wave_arr(kl)

c        para el continuo
         real*4 continuohNLTEarr(kld), continuoharr(kld)
         real*8 conhsra_RH, conhsra     

c        para los parametros atomicos
         real*4 lambda,nair,wlengt1_arr(kl),loggf,gf_arr(kl),energy_arr(kl)
         real*4 lambda_arr(kl),croot_arr(kl),wvac_arr(kl),weight_arr(kl)
         real*8 wlengtdbl, wavedbl_arr(kl), elemcode_arr(kl),wlengt1dbl,wcondbl

c        para el damping        
         real*4 chydro_arr(kl),xmu1_arr(kl),xmu2_arr(kl),xmu3_arr(kl)
         real*4 vv_arr(kl),beta_arr(kl)

c        para el patron zeeman
         parameter (mc=20)  !numero maximo de componentes zeeman
         character atom*2
         character design(2)*1
         integer   mult(2),ji(2),jf(2)
         integer   nel_arr(kl)
         real*4    tam(2),abu,melectron,mhidrogeno
         real*4    dlp(mc),dll(mc),dlr(mc),sp(mc),sl(mc),sr(mc)
         real*4    dlp_arr(mc,kl),dll_arr(mc,kl),dlr_arr(mc,kl)
         real*4    sp_arr(mc,kl),sl_arr(mc,kl),sr_arr(mc,kl)
         integer   np_arr(kl),nl_arr(kl),nr_arr(kl),icontinuo_normz
         real*4    chi10_arr(kl),chi20_arr(kl),abu_arr(kl)        
         real*4    gLS1,gLS2

c        otras variables de los commons
         real*8    loggf8_arr(kl)

         common/responde2/ist,ntau,ntl,nlin,npas,nble
         common/brklm/ntotal_lines,atom_arr,istage_arr,alfa_arr,sigma_arr,wave_arr
         common/wavearrdble/wavedbl_arr
         common/loggfarr/gf_arr,energy_arr
         common/ldeo/dlamda,dlamda0
         common/patronzeeman/dlp_arr,dll_arr,dlr_arr,sp_arr,sl_arr,sr_arr,np_arr,nl_arr,nr_arr
         common/datosatom/wvac_arr,weight_arr,croot_arr,wlengt1_arr,lambda_arr
         common/datosdamping/chydro_arr,xmu1_arr,xmu2_arr,xmu3_arr,vv_arr,beta_arr
         common/continuosarr/continuoharr
         common/continuosNLTEarr/continuohNLTEarr         
         common/niveles/nlow_i,nup_i,linea_nlte,atom_nlte
         common/pot_ion/chi10_arr,chi20_arr,nel_arr,abu_arr
         common/elemcode/elemcode_arr
         common/loggfarr8/loggf8_arr
         common/componente_nlte/linea_all_nlte !para blends2 y blendscon2

         c=2.99792458e+10        !vel. de la luz en cm/seg
         piis=1./sqrt(3.1415926)
         avog=6.023e23
         bol=1.3807e-16          !erg/s
         pir=3.1415926
         melectron=9.1094e-24
         mhidrogeno=1.67442e-24
         xmasaproton=1.6526e-24
         avo=6.023e23
         bohr=0.0529177249e-7    !cm
         uma=1.660540e-24
         gas=8.31451e7           !constante de los gases en cgs
         v0=1e6                  !cm/s
         coc3=1.212121           !cociente polarizabilidad del H2 con HI
         coc2=.3181818           !idem para el He

         ixx=0
         ikk0=0
         iii=0

         do iln=1,ntl
            linea_all_nlte(iln)=0   !1= si algun blend de la linea iln es nlte
            do ible=1,nble(iln)
               ixx=ixx+1
               nxx=nlin(ixx)
               call leelineas(ixx,nxx,atom,istage,wlengtdbl,zeff,energy,
     &                        loggf,mult,design,tam,alfa,sigma,nlow,nup,gLS1,gLS2)

               atom_arr(ixx)=atomic_number(atom)
               elemcode_arr(ixx)=dble(atom_arr(ixx))+(istage-1)/100.0d0
               energy_arr(ixx)=energy
               gf_arr(ixx)=1.e1**(loggf)
               loggf8_arr(ixx)=dble(loggf)
               nlow_i(ixx)=nlow
               nup_i(ixx)=nup
               atom_nlte(ixx)="  "
               linea_nlte(ixx)=0
               if(nlow .ne. 0 .or. nup .ne. 0)then
                  linea_nlte(ixx)=ixx   !indice de la componente nlte (0 si es lte)
                  atom_nlte(ixx)=atom
                  linea_all_nlte(iln)=1 !1= si algun blend de la linea iln es nlte
               end if

               istage_arr(ixx)=istage
               alfa_arr(ixx)=alfa
               sigma_arr(ixx)=sigma
               wavedbl_arr(ixx)=wlengtdbl
               wlengt=real(wlengtdbl)
               wave_arr(ixx)=wlengt  
               if(ible.eq.1)wlengt1=wlengt
               if(ible.eq.1)wlengt1dbl=wlengtdbl
               dlamda0(iln)=wlengt1

c parametros atomicos (estamos dentro de los do en lineas y en blends)
c llamo a atmdat que devuelve weight (peso molecular),abu (abunda
c cia), chi1,chi2 (pot.de ionizacion del atomo neutro e ion),
c u1,u2,u3 (funciones de particion atomo neutro,ion,ion2).
               call atmdatb(atom,0.,nel,weight,abu,chi10,chi20,u1,u2,u3,
     &                      du1,du2,du3)
               abu_arr(ixx)=abu
               chi10_arr(ixx)=chi10
               chi20_arr(ixx)=chi20
               nel_arr(ixx)=nel
c calculo la l.d.o. en el vacio (wvac), con la function refrax
               nair=1.0004
               if(wlengt.ge.1800.)nair=refrax(wlengt*1.e-4,15.,760.,0.)

c ATENCION en blends2 existe lambda (escalar) y 
c y dlamda=(dlongd(ikk)+(wlengt1-wlengt)*1.e3)*1.e-11 (escalar)
c mientras aqui dlamda(ikk)=dlongd(ikk).
c cambiar en la linea 725 de blends2 dlamda por otra variable y
c comprobar que todos los dlongd se han convertido a dlamda.
               lambda=wlengt*1.e-8  !l.d.o. en el aire en cm.
               lambda_arr(ixx)=lambda  !enviar a blends2 cargar antes de linea 465 para kappach
               wvac=nair*lambda        !l.d.o. en el vacio en cm.
               weinv=1./weight
               croot=1.66286e+8*weinv  !2r/m  para anchura doppler(cm**2/s**2)
               wvac_arr(ixx)=wvac      !enviar a blends2 caragar antes de linea 559  construi ahi wc=wvac/c
               weight_arr(ixx)=weight  !constuir  weinv=1./weight en blends2 antes de la linea 573 (donde se invocan alfa y sigma
               croot_arr(ixx)=croot    !enviar a blends2 cargar antes de linea 359  
               wlengt1_arr(ixx)= wlengt1  !enviar a blends2 cargar antes de linea  648
               dlo=4.6686e-5*wvac*wvac !separac.de l.d.o=dlo*h*(m1*g1-m2*g2)

c calculo chydro (gamma6)coef.van der waals ensanch.colisional
c corregido con zeff (termino 'semiempirico' caca de la vaca)
c (se supone despreciable el debido a stark)
c si el damping asi calculado es mayor que 3 lo reescalaremos

               if(abs(alfa).lt.1.e-25.or.abs(sigma).lt.1.e-25)then
                  chi1=chi10
                  if(istage.eq.2)chi1=chi20
                  eupper=energy+1.2398539e-4/wvac
                  ediff1=amax1(chi1-eupper-chi10*float(istage-1),1.0)

                  ediff2=amax1(chi1-energy-chi10*float(istage-1),3.0)
                  chydro=lambda*1.e1**(+.4*alog10(1./ediff1**2-1./ediff2**2)-
     &               12.213)*5.34784e+3
                  chydro=chydro*zeff
                  if(istage.eq.2)chydro=chydro*1.741
               else
c xmu1 es la masa reducida del hidrogeno y el atomo. Los indices 2
c y 3 hacen referencia al helio neutro y al hidrogeno molecular.
                  xmu1=uma*(1.008*weight)/(1.008+weight)
                  xmu2=uma*(4.0026*weight)/(4.0026+weight)
                  xmu3=uma*(2.016*weight)/(2.016+weight)

                  arr=2.-alfa*.5-1.
                  gammaf=1.+(-.5748646+(.9512363+(-.6998588+(.4245549-
     &               .1010678*arr)*arr)*arr)*arr)*arr
                  vv=(1.-alfa)/2.
 
                  beta=lambda*2*(4./pir)**(alfa/2.)*gammaf*(v0**alfa)*sigma*
     &               ((8.*bol/pir)**vv)
               end if
               chydro_arr(ixx)=chydro
               xmu1_arr(ixx)=xmu1
               xmu2_arr(ixx)=xmu2
               xmu3_arr(ixx)=xmu3   
               vv_arr(ixx)=vv
               beta_arr(ixx)=beta

c subniveles zeeman
c calculo la parte entera "ji" y la parte fraccionaria "jf" del
c momento angular total (tam), necesarios para zeeman
               do i=1,2
                  ji(i)=int(tam(i))
                  jf(i)=int(10*(tam(i)-ji(i)))
               end do

               if(gLS1.gt.0 .or. gLS2.gt.0)then
                  call zeeman_jk(mc,mult,design,tam,ji,jf,dlo,np,nl,nr,dlp,dll,dlr,
     &                           sp,sl,sr,gLS1,gLS2)
               else
                  call zeeman(mc,mult,design,tam,ji,jf,dlo,np,nl,nr,dlp,dll,dlr,
     &                        sp,sl,sr)
               endif

               do i=1,np  !cargar de la misma manera en blends2
                  dlp_arr(i,ixx)=dlp(i)
                  sp_arr(i,ixx)=sp(i)
               end do   
               do i=1,nl
                  dll_arr(i,ixx)=dll(i)
                  sl_arr(i,ixx)=sl(i)
               end do 
               do i=1,nr
                  dlr_arr(i,ixx)=dlr(i)
                  sr_arr(i,ixx)=sr(i)
               end do 
               np_arr(ixx)=np
               nl_arr(ixx)=nl
               nr_arr(ixx)=nr
            end do  !end do en blends
            do ipas=1,npas(iln)
               iii=iii+1
               wcondbl=wlengt1dbl+dlamda(iii)*1.d-3
               continuohNLTEarr(iii)=real(conhsra_RH(wcondbl)) 
               continuoharr(iii)=real(conhsra(wcondbl))
            end do 
         end do  !end do en lineas     
         ntotal_lines=ixx

         return
         end
