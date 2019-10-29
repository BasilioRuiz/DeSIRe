c lines_kurucz transfrom an atomic parameters file (SIR FORMAT)
c to a Kurucz format file
c ...................................................................

        subroutine lines_kurucz(nomlineas)
        
 	include 'PARAMETER'       
        character*100 nomlineas,kurpath,kurfile,a
        integer mult(2),nlow,nup,istage,ist(4),ntau
        integer nlin(kl),npas(kl),ntl,nble(kl)
        real tam(2),loggf,wlengt,zeff,alfa,sigma,energy
        character design(2)*1,atom*2,linea*100  
        common/responde2/ist,ntau,ntl,nlin,npas,nble
        
c        print*,'lines_kurucz ntl nlin=',ntl,npas(1),nlin(1),npas(1)
c        stop
        
        call read_keyword_input_RH('KURUCZ_DATA',kurpath)
        open(33,file=kurpath)
        read(33,*)kurfile
c        print*,'lines_kurucz : voy a escribir',nomlineas,'en ',kurfile
        close(33)
        
c counting lines in nomlineas       
        open(33,file=nomlineas,err=901) 
         ilines=1
         do while (ilines .lt. kl)
            read(33,*,end=900)a
            ilines=ilines+1
         end do   
900     close(33)   

        kurfile='TEST.'//kurfile !Borrar esta línea

        print*,'lines_kurucz: hemos leido:',ilines,' lineas'
        print*,'vamos a escribir en (prueba)',kurfile
        open(33,file=kurfile)
        
        ixx=0
        do iln=1,ntl
	   do ible=1,nble(iln)
	      ixx=ixx+1
              nxx=nlin(ixx) 
	      call leelineasiii(nxx,atom,istage,wlengt,zeff,energy,
     &                       loggf,mult,design,tam,alfa,sigma,nlow,nup)
c              print*,ixx,atom,istage,wlengt,nlow,nup
              if(nlow .eq. 0 .and. nup .eq. 0)then !LTE case
                 call write_kurucz(atom,istage,wlengt,energy,
     &                       loggf,mult,design,tam)
              end if
           end do
        end do   
        close(33)
        
        return
901     print*,'The file ', nomlineas,' does not exist'
        stop
        end
        
        subroutine write_kurucz(atom,istage,wlengt,energy,
     &                       loggf,mult,design,tam) 
     
        real*4 tam(2),loggf,wlengt,energy, wavenm, elemcode
        real*8 wlengt_air,w_vac,wlengt_vac
        real*8 aloggamma_rad,aloggamma_4,aloggamma_6
        
        real*4 eVtocm_1,energy_low,energy_up,hyper
        integer mult(2),istage,atomic_number,iso,nonlteindex
        character design(2)*1,atom*2,linea*100  
        character*10 label1,label2
        character*1 code1(21),aString,zero  
        character*4 nada
        parameter (eVtocm_1=8065.5443)   

        wavenm=wlengt/10.
        elemcode=atomic_number(atom)*1.0+(istage-1)/100.  !element number + charge/100.
        energy_low=energy*eVtocm_1
c        energy_up=energy_low+1e8/wlengt
        
        write(aString, 787)mult(1)
        label1=aString//design(1)
        write(aString, 787)mult(2)
        label2=aString//design(2)
        wlengt_air=dble(wlengt)
        wlengt_vac=w_vac(wlengt_air)*1.d-8 ! vacumm ldo in cm 
        
        energy_up=energy_low+1.e0/wlengt_vac
                
        aloggamma_rad=dlog(0.22233d0/wlengt_vac)
        
        aloggamma_4=-6.00d0
        aloggamma_6=-7.60d0
        nonlteindex=0
        iso=0
        hyper=0.0
        zero='0'
        nada='    '

        write(*,786)wavenm,loggf,elemcode,energy_low,tam(1),label1,
     &  energy_up,tam(2),label2,aloggamma_rad,aloggamma_4,aloggamma_6,
     &  nada,nonlteindex,nonlteindex,iso,hyper,iso,hyper,
     &  iso,iso,zero,zero,zero,zero,
     &  iso,nada,iso,iso,iso
     
        write(33,786)wavenm,loggf,elemcode,energy_low,tam(1),label1,
     &  energy_up,tam(2),label2,aloggamma_rad,aloggamma_4,aloggamma_6,
     &  nada,nonlteindex,nonlteindex,iso,hyper,iso,hyper,
     &  iso,iso,zero,zero,zero,zero,
     &  iso,nada,iso,iso,iso
        
786    FORMAT(F10.4,F7.3,F6.2,F12.3,F5.2,1X,A10,
     & F12.3,F5.2,1X,A10,3F6.2,
     & A4,2I2,I3,F6.3,I3,F6.3,
     & 2I5,1X,A1,A1,1X,A1,A1,
     & i1,A3,2I5,I6)  
787    FORMAT(I1)  
c+++++++++++^^^^^^^++++++^^^^^^^^^^^^+++++^++++++++++^^^^^^^^^^^^+++++^++++++++++c      
c  630.1500 -0.718 26.00   45333.872  2.0 s6D)5s e5D   29469.022  2.0 5Dsp3P z5P  7.93 -5.46 -7.58K94  0 0  0 0.000  0 0.000    0    0           1503 1835      
c123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
c        10        20        30        40         50        60        70        80        90       100       120       120       130       140        150     160      
      
c Line files with names of the form GF* or HY* have the following 160 column format:
c1                                                                             80
c+++++++++++^^^^^^^++++++^^^^^^^^^^^^+++++^++++++++++^^^^^^^^^^^^+++++^++++++++++
c   800.7110  0.116 27.00   45924.947  3.5 (3F)5s e2F   33439.661  4.5 (3F)4p y2G
c    wl(nm)  log gf elem      E(cm-1)   J   label        E'(cm-1)   J'   label'  
c                   code                                                         
c                        [char*28 level descriptor  ][char*28 level descriptor  ]
c                                                                                
ccontinuing                                                                      
c81                                                                           160
c^^^^^^++++++^^^^^^++++^^++^^^++++++^^^++++++^^^^^+++++^+^+^+^+++^^^^^+++++^^^^^^
c  8.19 -5.38 -7.59K88  0 0 59-2.584 59 0.000  104  -77F6 -5 0    1140 1165     0
c  log   log   log ref NLTE iso log iso  log     hyper  F F'    eveglande     iso
c Gamma Gamma Gamma   level hyper f iso frac   shift(mK)     ^    oddglande shift
c  rad  stark  vdW    numbers                    E    E'     ^abc  (x1000)   (mA)
c                                                         I*1^char*3
c                                                           codes
c
cFORMAT(F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,A10,
c3F6.2,A4,2I2,I3,F6.3,I3,F6.3,2I5,1X,A1,A1,1X,A1,A1,i1,A3.2I5,I6) 
c 1 wavelength (nm)  air above 200 nm   F11.4
c 2 log gf  F7.3
c 3 element code = element number + charge/100.  F6.2
c 4 first energy level in cm-1   F12.3
c        (if allowed, with same parity as ground state) 
c        (negative energies are predicted or extrapolated)
c 5 J for first level   F5.1
c   blank for legibility   1X
c 6 label field for first level   A10
c 7 second energy level in cm-1   F12.3
c        (if allowed, with parity opposite first level) 
c        (negative energies are predicted or extrapolated)
c 8 J for second level   F5.1
c   blank for legibility   1X
c 9 label field for second level   A10
c10 log of radiative damping constant, Gamma Rad  F6.2 or F6.3
c11 log of stark damping constant/electron number. Gamma Stark  F6.2 or F6.3
c12 log of van der Waals damping constant/neutral hydrogen number, 
c       Gamma van der Waals   F6.2 or F6.3
c13 reference that can be expanded in subdirectory LINES   A4  
c14 non-LTE level index for first level   I2
c15 non-LTE level index for second level   I2
c16 isotope number   I3
c17 hyperfine component log fractional strength  F6.3
c18 isotope number  (for diatomics there are two and no hyperfine)   I3
c19 log isotopic abundance fraction   F6.3
c20 hyperfine shift for first level in mK to be added to E  I5
c21 hyperfine shift for second level in mK to be added to E'  I5
c   the symbol "F" for legibilty   1X
c22 hyperfine F for the first level    I1
c23 note on character of hyperfine data for first level: z none, ? guessed  A1
c   the symbol "-" for legibility    1X
c24 hyperfine F' for the second level  I1
c25 note on character of hyperfine data for second level: z none, ? guessed  A1
c26 1-digit code, sometimes for line strength classes   I1
c27 3-character code such as AUT for autoionizing    A3  
c28 lande g for the even level times 1000   I5
c29 lande g for the odd level times 1000   I5
c30 isotope shift of wavelength in mA 

     
     
        return
        end
