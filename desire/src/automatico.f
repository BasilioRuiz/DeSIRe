        subroutine automatico(mi,mp,ntotal4,difer,npos,rt,t)

        implicit real*4 (a-h,o-z) 
        include 'PARAMETER'
        parameter (kld4=4*kld)

        integer npos(*),ntotal4,nodosposibles(kt),icalerr
        real*4 difer(*),rt(*),derivada(kt),t(*)
        real*4 x(kt),z(kt),f(kt,kt),tau(kt)

        common/nodosposibles/nummax,nodosposibles
        common/ndata/ndata !para fperfil2 y automatico
        common/tau/tau
        common/calerr/icalerr !si calculo errores=1 else =0

        data iprimeraveza/0/

c calculamos los posibles nodos 
       ntau=mi
        if(iprimeraveza.eq.0)then 
           iprimeraveza=1
           k=2
           nodosposibles(1)=1
           nodosposibles(2)=2
           do j=3,ntau
              resto=(ntau-1)-(j-1)*((ntau-1)/(j-1))
              if(abs(resto).lt.1.e-4)then
                 k=k+1
                 nodosposibles(k)=j
              end if
           end do 
           nummax=k
        end if

c calculamos la derivada de la chi^2
        do i=1,mi
           s=0.
           do j=1,ndata
              k=(i-1)*ntotal4+npos(j)
              s=s+difer(j)*rt(k)
           end do
           derivada(i)=s
        end do

c buscamos el numero de nodos optimo
        if(icalerr.eq.0)then
           if(mp.gt.1)then
             call criterio(ntau,derivada,mi,mp)
           else
             mi=1
           end if
        end if

        if(mi.gt.1)then
            m=(ntau-1)/(mi-1) !n es el numero de taus
                              !m es el numero de pasos en cada subintervalo
            do i=1,mi         !mi es el numero de nodos total i el indice del nodo
               j=(i-1)*m+1    !j el indice en tau
               x(i)=tau(j)    !x(i) es el valor de tau en el nodo
            end do

            call splines22(x,x,mi-2,ntau,tau,z,f)
        end if
           
        call nodos_sub(rt,ntau,mi,f,ntotal4,t)
  
        return
        end

c_____________________________________________________________________________

        subroutine nodos_sub(rt,n,no,f,ntotal4,t)

        implicit real*4 (a-h,o-z) 
        include 'PARAMETER'

        integer no
        real*4 rt(*),f(kt,kt),gr(kt),t(*),y(kt)

        if(no.le.0)return

        if(no.eq.1)then   !si perturbacion constante
           ymedio=0
           do i=1,n
              ymedio=ymedio+t(i)
           end do
           ymedio=ymedio/float(n)

           do ilam=1,ntotal4
              sumart=0.
              do i=1,n
                 kk=(i-1)*ntotal4+ilam
                 sumart=sumart+rt(kk)
              end do
              rt(ilam)=sumart*ymedio 
           end do
        else

           m=(n-1)/(no-1)  !n es el numero de taus
                           !m es el numero de pasos en cada subintervalo
           do i=1,no       !no es el numero de nodos total i el indice del nodo
              j=(i-1)*m+1  !j el indice en tau
              y(i)=t(j)    !y(i) es el valor del parametro en el nodo
           end do

           do ilam=1,ntotal4
              do i=1,no
                 ggg=0.
                 do j=1,n
                    kk=(j-1)*ntotal4+ilam
                    ggg=ggg+rt(kk)*f(j,i)
                 end do
                 gr(i)=ggg
              end do
              do i=1,no
                 kk=(i-1)*ntotal4+ilam
                 rt(kk)=gr(i)*y(i)
              end do
           end do
        end if

        return
        end
        
c_____________________________________________________________________________

        subroutine criterio(ntau,derivada,mi,mp)

        implicit real*4 (a-h,o-z) 
        include 'PARAMETER'
        parameter (kld4=4*kld)

        integer nodosposibles(kt),paso
        real*4 derivada(*)

        common/nodosposibles/nummax,nodosposibles

c       data sesgo1/10./

        a0=abs(derivada(1))/2.
        a1=derivada(1)/2.
        do i=2,ntau-1
           a0=a0+abs(derivada(i))
           a1=a1+derivada(i)
        end do
        a0=a0+abs(derivada(ntau))/2.
        a1=a1+derivada(ntau)/2.

        a1=a1/a0

        if(abs(a1-1.) .lt. 1.e-4)then
          mi=1
          return
        end if

        amax=abs(a1)
        imax=1
c       sesgo1=10. !factor para penalizar un numero elevado de nodos
c       sesgo1=sesgo1*.9 !factor para penalizar un numero elevado de nodos
        sesgo1=0. 
        sesgo2=1.0 !factor para penalizar un numero elevado de nodos

        nummax2=1
        i=2  
        do while(i.lt.nummax.and.nodosposibles(i).le.mp)
           nummax2=i
           i=i+1 
        end do   

c nuevo criterio (refuerzo el caso 2) integrando desde la parte significativa
c distinta de cero y poniendo el nodo en la mitad
c calculo el maximo y el minimo de derivada
        dmax=-a0
        dmin=a0
        idmax=0
        idmin=0
        do i=1,ntau
           if(derivada(i).gt.dmax)then
             dmax=derivada(i)
             idmax=i
           end if   
           if(derivada(i).lt.dmin)then
             dmin=derivada(i)
             idmin=i
           end if  
        end do
        if(dmax*dmin.ge.0)then
          mi=1
          return
        end if

        i00n=min(idmax,idmin) !extremo (max o min) inferior
        i00x=max(idmax,idmin) !extremo superior


c buscamos el primer corte despues del maximo         
         i11=i00n
         do while(derivada(i11)*derivada(i11+1).gt.0.and.i11.lt.ntau-1)
           i11=i11+1
         end do

         aijmax=-10.*a0
         aijmin=10*a0

c al caso de los 2 nodos le echamos de comer a parte
        aij=derivada(1)/2.
        do j=2,i11-1
          aij=aij+derivada(j)
        end do
        aij=(aij+derivada(i11)/2.)/a0
        if(aij.ge.aijmax)aijmax=aij                     
        if(aij.le.aijmin)aijmin=aij
        aij=derivada(i11+1)/2.
        do j=i11+2,ntau-1
           aij=aij+derivada(j)
        end do
        aij=(aij+derivada(ntau)/2.)/a0
        if(aij.ge.aijmax)aijmax=aij                     
        if(aij.le.aijmin)aijmin=aij
        a=abs(aijmax-aijmin)*sesgo1
        if(abs(a-amax) .lt. abs(amax)*1.e-5)then
           mi=2
           return
        end if
        if(a.gt.amax)then
           amax=a
           imax=2
       end if

        do i=2,nummax2 
           nod=nodosposibles(i)
           paso=nint((ntau-1)/float(nod))
           aijmax=-10.*a0
           aijmin=10*a0
           do j=1,nod
              aij=derivada((j-1)*paso)/2.
              do k=(j-1)*paso+1,j*paso-1
                 aij=aij+derivada(k)
              end do
              aij=aij+derivada(j*paso)/2.
              aij=aij/a0
              if(aij.ge.aijmax)aijmax=aij               
              if(aij.le.aijmin)aijmin=aij
           end do
           a=abs(aijmax-aijmin)
           if(a.gt.sesgo2*amax)then !!!!!!!!!!!!!!!!!!!!!!
              amax=a
              imax=i
           end if
        end do  

        mi=nodosposibles(imax)

        return
        end
