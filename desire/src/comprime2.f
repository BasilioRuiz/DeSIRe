c "comprime2" extrae los valores de las atmosferas en los nodos
c y los devuelve en atmosr
c
c    i variable   i  variable
c    -   tau1     -    tau2   profundidad optica a 5000 /AA
c    1   t1       9    t2     temperatura en ambos modelos (k)
c    2   p1      10    p2     presion electronica (dinas/cm**2)
c    3   mic1    11    mic2   microturbulencia (cm/s)
c    4   h1      12    h2     campo magnetico (G)
c    5   v1      13    v2     velocidad eje z (cm/s)
c    6   g1      14    g2     gamma (radianes)
c    7   f1      15    f2     fi (radianes)
c    8   mac1    16    mac2   macroturbulencia (en Km/s)
c    -   ff1     17    ff2    factor de llenado (ff1=1-ff2)
c                18    %stray
c trabajo siempre con el ff del segundo modelo
c
c Basilio 22-3-93 
c Basilio  1-4-93 modifico la variable ff2 a (1-ff2)/ff2
c Basilio y Jose Carlos 9-1-95 pertu. aditiva para fi (no mult. a tan(fi))
c Basilio y Jose Carlos 6-2-95 pertu. aditiva para gamma (no mult. a tan(gamma))
c
c _____________________________________________________________
        subroutine comprime2(ntau,m,atmos,atmosr)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER' !incluidng vthresh: the threshold for velocity

        real*4 atmos(*),atmosr(*)
        integer m(*)
        common/offset/voffset  !para respuestas
        common/ivez/ivez

        if(ivez.eq.0)then
           ivez=1
c           vmin=20.e5
c           if(m(13).ne.0.or.m(5).ne.0)vmin=1.e20
c           do i=1,ntau              
c              if(atmos(2+13*ntau+i).lt.vmin.and.m(13).ne.0)
c     &           vmin=atmos(2+13*ntau+i)
c              if(atmos(5*ntau+i).lt.vmin.and.m(5).ne.0)vmin=atmos(5*ntau+i)
c           end do
c           voffset=vmin-vthresh    !cm/s
           voffset=-vthresh    !cm/s
        end if

        kred=0     !indice reducido
        kamp=ntau  !indice ampliado (los ntau puntos de tau1)

        do i=1,18  !do en grupos de variables (1=t,2=p,...etc)
           ntau2=ntau
           if(i.eq.8.or.i.eq.16.or.i.eq.17.or.i.eq.18)ntau2=1  !si mac1,mac2 o ff22,%
           if(m(i).eq.1)then  !si pert. constante promedio la atm.
              kred=kred+1
              sum=0.
              do ii=1,ntau2
                 sum=sum+atmos(kamp+ii) 
              end do
              atmosr(kred)=sum/float(ntau2)
              if(i.eq.17)atmosr(kred)=(1.e0-atmosr(kred))/atmosr(kred)  !ff
              if(i.eq.18)atmosr(kred)=atmosr(kred)/(100.-atmosr(kred))  !%
              if(i.eq.5.or.i.eq.13)atmosr(kred)=atmosr(kred)-voffset
           else if(m(i).gt.1)then
              mm=(ntau-1)/(m(i)-1)   !espaciado entre nodos 
              do j=1,m(i)
                 kred=kred+1
                 jj=kamp+(j-1)*mm+1
                 atmosr(kred)=atmos(jj)
                 if(i.eq.5.or.i.eq.13)atmosr(kred)=atmosr(kred)-voffset                                      
              end do
           end if     
           kamp=kamp+ntau2
              if(i.eq.8)kamp=kamp+ntau+1  !los ntau puntos de tau2 y el de ff1
       end do
       
       return
       end


c _____________________________________________________________
c
c "comprime3" = comprime2 for 1 atmosphere
c _____________________________________________________________

        subroutine comprime3(ntau,m,atmos,atmosr)

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        real*4 atmos(*),atmosr(*)
        integer m(*)
        common/offset/voffset  !para respuestas
        common/ivez/ivez

        if(ivez.eq.0)then
           ivez=1
c           vmin=20.e5
c           if(m(13).ne.0.or.m(5).ne.0)vmin=1.e20
c           do i=1,ntau              
c              if(atmos(2+13*ntau+i).lt.vmin.and.m(13).ne.0)
c     &           vmin=atmos(2+13*ntau+i)
c              if(atmos(5*ntau+i).lt.vmin.and.m(5).ne.0)vmin=atmos(5*ntau+i)
c           end do
c          voffset=vmin-vthresh     !cm/s
           voffset=-vthresh
        end if

        kred=0     !indice reducido
        kamp=ntau  !indice ampliado (los ntau puntos de tau1)

        do i=1,8   !do en grupos de variables (1=t,2=p,...etc)
           ntau2=ntau
           if(i.eq.8)ntau2=1  !si mac1,mac2 o ff22,%
           if(m(i).eq.1)then  !si pert. constante promedio la atm.
              kred=kred+1
              sum=0.
              do ii=1,ntau2
                 sum=sum+atmos(kamp+ii) 
              end do
              atmosr(kred)=sum/float(ntau2)
              if(i.eq.5)atmosr(kred)=atmosr(kred)-voffset
           else if(m(i).gt.1)then
              mm=(ntau-1)/(m(i)-1)   !espaciado entre nodos 
              do j=1,m(i)
                 kred=kred+1
                 jj=kamp+(j-1)*mm+1
                 atmosr(kred)=atmos(jj)
                 if(i.eq.5)atmosr(kred)=atmosr(kred)-voffset                                    
              end do
           end if     
           kamp=kamp+ntau2
           if(i.eq.8)kamp=kamp+ntau+1	!los ntau puntos de tau2 y el de ff1
        end do
        return
        end


c _____________________________________________________________
