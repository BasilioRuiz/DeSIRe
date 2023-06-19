c Solve by svd the problem A*x=B being x & B vector of 4 componentes and A a 4x4 matrix
c matrix A is destroyed in the output

        subroutine linearsystem(a,b,x)

	integer n, j, i
        real*8 a(4,4),b(4),w(4),v(4,4),x(4)
	real*8 wt,ww,wi(2),tol,wmax,wmin

        call svd4(a,w,v)

        n=4 
        tol=1.d-8
c calculamos el maximo
	      wmax=0.d0
	      do j=1,4
	         if(w(j).gt.wmax)wmax=w(j)
	      end do

c anulamos los terminos menores que el max. por tol.
	      wmin=wmax*tol
	      do j=1,4
	         if(w(j).lt.wmin)w(j)=0
	      end do	

	call SVBKSBr8(a,w,v,b,x)

	return
	end

c-----------------------------------------------
      subroutine svd4(da,w,v)

      real*8 rv1(4),w(4),v(4,4)
      real*8 g,anorm,scale,s,f,h,ab,c,z,y,x,da(4,4)
      integer i,j,k,l,m,n
      
      g=0.0d0

      scale=0.0d0

      anorm=0.0d0

      m=4
      n=4

      do 25 i=1,n

        l=i+1

        rv1(i)=scale*g

        g=0.0d0

        s=0.0d0

        scale=0.0

        if (i.le.m) then

          do 11 k=i,m

            scale=scale+dabs(da(k,i))

11        continue
          if (dabs(scale).gt.1.d-20) then

            do 12 k=i,m

              da(k,i)=da(k,i)/scale

              s=s+da(k,i)*da(k,i)

12          continue

            f=da(i,i)

            g=-sign(dsqrt(s),f)

            h=f*g-s

            da(i,i)=f-g

            if (i.ne.n) then

              do 15 j=l,n

                s=0.0d0

                do 13 k=i,m

                  s=s+da(k,i)*da(k,j)

13              continue

                f=s/h

                do 14 k=i,m

                  da(k,j)=da(k,j)+f*da(k,i)

14              continue

15            continue

            endif

            do 16 k= i,m

              da(k,i)=scale*da(k,i)

16          continue

          endif

        endif

        w(i)=scale *g

        g=0.0d0

        s=0.0d0

        scale=0.0d0

        if ((i.le.m).and.(i.ne.n)) then

          do 17 k=l,n

            scale=scale+dabs(da(i,k))

17        continue

          if (dabs(scale).gt.1.d-20) then

            do 18 k=l,n

              da(i,k)=da(i,k)/scale

              s=s+da(i,k)*da(i,k)

18          continue

            f=da(i,l)

            g=-sign(dsqrt(s),f)

            h=f*g-s

            da(i,l)=f-g

            do 19 k=l,n

              rv1(k)=da(i,k)/h

19          continue

            if (i.ne.m) then

              do 23 j=l,m

                s=0.0d0

                do 21 k=l,n

                  s=s+da(j,k)*da(i,k)

21              continue

                do 22 k=l,n

                  da(j,k)=da(j,k)+s*rv1(k)

22              continue

23            continue

            endif

            do 24 k=l,n

              da(i,k)=scale*da(i,k)

24          continue

          endif

        endif

	ab=dabs(w(i))+dabs(rv1(i))
        anorm=max(anorm,ab)

25    continue

      do 32 i=n,1,-1

        if (i.lt.n) then
          if (dabs(g).gt.1.d-20) then

            do 26 j=l,n
              v(j,i)=(da(i,j)/da(i,l))/g 
26          continue

            do 29 j=l,n

              s=0.0

              do 27 k=l,n

                s=s+da(i,k)*v(k,j)

27            continue

              do 28 k=l,n

                v(k,j)=v(k,j)+s*v(k,i)
28            continue

29          continue

          endif

          do 31 j=l,n

            v(i,j)=0.0

            v(j,i)=0.0

31        continue

        endif

        v(i,i)=1.0

        g=rv1(i)

        l=i

32    continue

      do 39 i=n,1,-1

        l=i+1

        g=w(i)*1.d0

        if (i.lt.n) then

          do 33 j=l,n

            da(i,j)=0.0d0

33        continue

        endif
        if (dabs(g).gt.1.d-20) then

          g=1.0d0/g

          if (i.ne.n) then

            do 36 j=l,n

              s=0.0

              do 34 k=l,m

                s=s+da(k,i)*da(k,j)

34            continue

              f=(s/da(i,i))*g

              do 35 k=i,m

                da(k,j)=da(k,j)+f*da(k,i)

35            continue

36          continue

          endif

          do 37 j=i,m

            da(j,i)=da(j,i)*g

37        continue

        else

          do 38 j= i,m

            da(j,i)=0.0

38        continue

        endif

        da(i,i)=da(i,i)+1.0

39    continue

      do 49 k=n,1,-1

        do 48 its=1,30

          do 41 l=k,1,-1

            nm=l-1
            if ( dabs(rv1(l)) .lt. 1.d-20 )  go to 2
            if ( abs(w(nm)) .lt. 1.e-15)  go to 1

41        continue

1         c=0.0d0

          s=1.0d0

          do 43 i=l,k

            f=s*rv1(i)
            if ( dabs(f) .gt. 1.d-20) then

              g=w(i)*1.d0

              h=dsqrt(f*f+g*g)

              w(i)=real(h)

              h=1.0d0/h

              c= (g*h)

              s=-(f*h)

              do 42 j=1,m

                y=da(j,nm)

                z=da(j,i)

                da(j,nm)=(y*c)+(z*s)

                da(j,i)=-(y*s)+(z*c)

42            continue

            endif

43        continue

2         z=w(k)*1.d0

          if (l.eq.k) then

            if (z.lt.0.0) then

              w(k)=-real(z)

              do 44 j=1,n

                v(j,k)=-v(j,k)

44            continue

            endif

            go to 3

          endif

          if (its.eq.50) then
	      w(1)=1.
	      do i=2,n
	         w(i)=0.
	      end do
	      goto 999
c	      pause 'no convergence in 50 iterations SVDCMP'
          end if

          x=w(l)*1.d0

          nm=k-1

          y=w(nm)*1.d0

          g=rv1(nm)

          h=rv1(k)
          
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)

          g=dsqrt(f*f+1.0)

          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x

          c=1.0d0

          s=1.0d0

          do 47 j=l,nm

            i=j+1

            g=rv1(i)

            y=w(i)*1.d0

            h=s*g

            g=c*g

            z=dsqrt(f*f+h*h)

            rv1(j)=z

            c=f/z

            s=h/z

            f= (x*c)+(g*s)

            g=-(x*s)+(g*c)

            h=y*s

            y=y*c

            do 45 nm=1,n

              x=v(nm,j)

              z=v(nm,i)

              v(nm,j)= (x*c)+(z*s)

              v(nm,i)=-(x*s)+(z*c)
                           
45          continue

            z=dsqrt(f*f+h*h)

            w(j)=z
            if (dabs(z).gt.1.d-20) then

              z=1.0d0/z

              c=f*z

              s=h*z

            endif

            f= (c*g)+(s*y)

            x=-(s*g)+(c*y)

            do 46 nm=1,m

              y=da(nm,j)

              z=da(nm,i)

              da(nm,j)= (y*c)+(z*s)

              da(nm,i)=-(y*s)+(z*c)

46          continue

47        continue

          rv1(l)=0.0

          rv1(k)=f

          w(k)=x

48      continue

3       continue

49    continue

       
999      return

      end

c ------------------------------------------------------------
c evaluates the square of modified (K/eta-1) absroption matrix
      SUBROUTINE absormodfsquare(K,K2)
      real*8 K(4,4),K2(4,4)
      real*8 a,b,c,d,e,f,a2,b2,c2,d2,e2,f2
      
      a=K(1,2)
      b=K(1,3)
      c=K(1,4)
      d=K(2,3)
      e=K(2,4)
      f=K(3,4)
      
      a2=a*a
      b2=b*b
      c2=c*c
      e2=e*e
      f2=f*f
      
      K2(1,1)=a2+b2+c2
      K2(2,2)=a2-d2-e2
      K2(3,3)=b2-d2-f2
      K2(4,4)=c2-e2-f2

      K2(2,1)=d*b+e*c
      K2(1,2)=-K2(2,1)
      
      K2(3,1)=c*f-d*a
      K2(1,3)=-K2(3,1)
      
      K2(1,4)=a*e+b*f
      K2(4,1)=-K2(1,4)
      
      K2(2,3)=a*b-e*f
      K2(3,2)=K2(2,3)
      
      K2(2,4)=a*c+d*f
      K2(4,2)=K2(2,4)
      
      K2(3,4)=b*c-d*e
      K2(4,3)=K2(3,4)

      return
      end
c ------------------------------------------------------------
      SUBROUTINE matrizmult(A,B,C)
      real*8 A(4,4),B(4,4),C(4,4),suma
      integer i,j,k
      
      do i=1,4
         do j=1,4
            suma=0.
            do k=1,4
               suma=suma+A(i,k)*B(k,j)
            end do
            C(i,j)=suma
         end do
      end do
      return
      end
c ------------------------------------------------------------
      SUBROUTINE matrizmultvec(A,B,C)     
      real*8 A(4,4),B(4),C(4),suma
      integer i,j
      
      do i=1,4
         suma=0.
         do j=1,4
            suma=suma+A(i,j)*B(j)
         end do
         C(i)=suma
      end do
      return
      end      
      
c-------------------------------------------------------------
      SUBROUTINE SVBKSBr8(U,W,V,B,X)

      real*8 U(4,4),W(4),V(4,4),B(4),X(4),TMP(4)
      real*8 S
      integer J,JJ
      
      DO 12 J=1,4
         S=0.
         IF(dabs(W(J)).GT.1.d-18)THEN
            DO 11 I=1,4
               S=S+U(I,J)*B(I)
11          CONTINUE
            S=S/W(J)
         ENDIF
         TMP(J)=S
12    CONTINUE
      DO 14 J=1,4
         S=0.
         DO 13 JJ=1,4
            S=S+V(J,JJ)*TMP(JJ)
13       CONTINUE
         X(J)=S
14    CONTINUE

      RETURN
      END

c-----------------------------------------------
c expansion of exp(-x) valid X > 0.04 
      REAL*8 FUNCTION DEXPD(X)  
      REAL*8 X
      
      if(x .gt. 0.04)then
        DEXPD=dexp(-x)
      else
        DEXPD=1.d0-x*(1.d0-x*(0.5d0-0.16666666666666667d0*x))
      end if
      return
      end

c------------------------------------------------------      
        real*8 function bezier_cubic_deriv(x1,x2,x3,y1,y2,y3)

        include 'PARAMETER'
        real*8 x1,x2,x3,y1,y2,y3,yder
        real*8 hk,hk_1,dkplus,dkminus,dkprod,alphak
        real*8 abhk,abhk_1,eps

        eps=dabs(x3-x1)*1.d-8
        hk=x3-x2
        hk_1=x2-x1   
        abhk=dabs(hk)
        abhk_1=dabs(hk_1)
        if(abhk .lt. eps .and. abhk_1 .lt. eps)then 
           call error(KSTOP,'bezier_cubic_deriv',
     &                'Stop at bezier_cubic_deriv because x1=x2=x3')
        endif
        if(abhk .ge. eps .and. abhk_1 .ge. eps) then 
            alphak=(1.+hk/(hk+hk_1))/3.       
            dkplus=(y3-y2)/hk
            dkminus=(y2-y1)/hk_1
        else 
            if(abhk .ge. eps) then 
               alphak=0. 
               dkplus=(y3-y2)/hk
               dkminus=1.
               if(dkplus .lt. -eps)dkminus=-1.
            endif
            if(abhk_1 .ge. eps) then 
               alphak=1. 
               dkminus=(y2-y1)/hk_1  
               dkplus=1.
               if(dkminus .lt. -eps)dkplus=-1.
            endif
        endif
                   
        dkprod=dkplus*dkminus
        yder=0.
        if(dkprod.gt.0.)yder=dkprod/(alphak*dkplus+(1.-alphak)*dkminus)
        bezier_cubic_deriv=yder
        return
        end    

c-----------------------------------------------
        real*8 function bezier_cubic_deriv2(hk_1,hk,y21,y32)

        include 'PARAMETER'
        real*8 hk_1,hk,y21,y32,yder
        real*8 dkplus,dkminus,dkprod,alphak
        real*8 abhk,abhk_1,eps
 
c hk_1=x21=x2-x1 hk=x32=x3-x2 
c y21=y2-y1 y32=y3-y2

        abhk=dabs(hk)
        abhk_1=dabs(hk_1)
        eps=(abhk+abhk_1)*1.d-8
        yder=0.
        
        if(abhk .lt. eps .and. abhk_1 .lt. eps)then 
           call error(KSTOP,'bezier_cubic_deriv2',
     &                'Stop at bezier_cubic_deriv2 because x1=x2=x3')
        endif
        if(abhk .ge. eps .and. abhk_1 .ge. eps) then 
            alphak=(1.+(hk/(hk+hk_1)))/3.       
            dkplus=y32/hk
            dkminus=y21/hk_1
            dkprod=dkplus*dkminus
            if(dkprod.gt.0.)yder=dkprod/(alphak*dkplus+(1.-alphak)*dkminus)           
        else 
            if(abhk .ge. eps)yder=y32/hk 
            if(abhk_1 .ge. eps)yder=y21/hk_1
        endif
                   
        bezier_cubic_deriv2=yder
        return
        end          
        
c----------------------------------------------------------------------------
c bezier_cubic_deriv_array
c It evaluates the derivatives of one array Sk(nk) 
c It uses function bezier_cubic_deriv

        subroutine bezier_cubic_deriv_array(xk,Sk,nk,SKder)  
        real*8 xk(*),Sk(*),Skder(*)
        integer j,nk
        real*8 xki_1,xki,xki1
        real*8 yki_1,yki,yki1
        real*8 bezier_cubic_deriv2,    bezier_cubic_deriv
              
        xki=xk(1) 
        xki1=xk(2) 
        yki=Sk(1) 
        yki1=Sk(2) 
        xki_1=xki  
        yki_1=yki   
c        Skder(1)=bezier_cubic_deriv(xki_1,xki,xki1,yki_1,yki,yki1)    
        Skder(1)=bezier_cubic_deriv2(xki-xki_1,xki1-xki,yki-yki_1,yki1-yki) 
             
        do j=2,nk-1
           xki_1=xki 
           xki=xki1
           xki1=xk(j+1) 
           yki_1=yki 
           yki=yki1
           yki1=Sk(j+1) 
c           Skder(j)=bezier_cubic_deriv(xki_1,xki,xki1,yki_1,yki,yki1)   
           Skder(j)=bezier_cubic_deriv2(xki-xki_1,xki1-xki,yki-yki_1,yki1-yki) 
        end do
        Skder(nk)=Skder(nk-1)
        return
        end  
c-----------------------------------------------
c bezier_cubic_deriv_array4
c It evaluates the derivatives of one array Kk(4,4,nk) 
c It uses function bezier_cubic_deriv

        subroutine bezier_cubic_deriv_array4(xk,Sk,nk,SKder)  
        real*8 xk(nk),Sk(4,nk),Skder(4,nk)
        integer i,j,k,nk
        real*8 x21,x32,y21,y32
        real*8 bezier_cubic_deriv2
                          
        x21=0.
        x32=xk(2)-xk(1) 
        y21=0.
        do j=1,4
           y32=Sk(j,2)-Sk(j,1)              
           Skder(j,1)=bezier_cubic_deriv2(x21,x32,y21,y32)  
        end do           
                                             
        do k=2,nk-1
          x21=x32
          x32=xk(k+1)-xk(k)
          do j=1,4
               y32=Sk(j,k+1)-Sk(j,k)
               y21=Sk(j,k)-Sk(j,k-1)
               Skder(j,k)=bezier_cubic_deriv2(x21,x32,y21,y32)    
          end do  
        end do
        do j=1,4
           Skder(j,nk)=Skder(j,nk-1)
        end do
        return
        end         
c-----------------------------------------------
c bezier_cubic_deriv_array4ini
c It evaluates the derivatives of one array Kk(4,4,nk) 
c It uses function bezier_cubic_deriv

        subroutine bezier_cubic_deriv_array4ini(xk,Sk,ni,nk,SKder)  
        real*8 xk(nk),Sk(4,nk),Skder(4,nk)
        integer i,j,k,nk,ni
        real*8 x21,x32,y21,y32
        real*8 bezier_cubic_deriv2
        
        do k=1,ni-1
           do j=1,4
              Skder(j,k)=0.
           end do  
        end do
        if(ni .le. 1)then                        
          x21=0.
          x32=xk(2)-xk(1) 
          y21=0.
          do j=1,4
             y32=Sk(j,2)-Sk(j,1)              
             Skder(j,1)=bezier_cubic_deriv2(x21,x32,y21,y32)  
          end do                                                       
          do k=2,nk-1
            x21=x32
            x32=xk(k+1)-xk(k)
            do j=1,4
               y32=Sk(j,k+1)-Sk(j,k)
               y21=Sk(j,k)-Sk(j,k-1)
               Skder(j,k)=bezier_cubic_deriv2(x21,x32,y21,y32)    
            end do  
          end do
        else
           x32=xk(ni)-xk(ni-1) 
           do k=ni,nk-1
            x21=x32
            x32=xk(k+1)-xk(k)
            do j=1,4
               y32=Sk(j,k+1)-Sk(j,k)
               y21=Sk(j,k)-Sk(j,k-1)
               Skder(j,k)=bezier_cubic_deriv2(x21,x32,y21,y32)    
            end do  
          end do  
        end if
                   
        do j=1,4
           Skder(j,nk)=Skder(j,nk-1)
        end do
        return
        end         

c-----------------------------------------------
c bezier_cubic_deriv_array44 
c It evaluates the derivatives of one array Kk(4,4,nk) 
c It uses function bezier_cubic_deriv

        subroutine bezier_cubic_deriv_array44(xk,Kk,nk,KKder)  
        real*8 xk(nk),Kk(4,4,nk),Kkder(4,4,nk)
        integer i,j,k,nk
        real*8 x21,x32,y21,y32
        real*8 bezier_cubic_deriv2
                          
        x21=0.
        x32=xk(2)-xk(1) 
        y21=0.
        do j=1,4
           do i=1,4 
              y32=Kk(i,j,2)-Kk(i,j,1)              
c              Kkder(i,j,k)=bezier_cubic_deriv(xki_1,xki,xki1,yki_1,yki,yki1)  
              Kkder(i,j,k)=bezier_cubic_deriv2(x21,x32,y21,y32)  
           end do  
        end do           
                                             
        do k=2,nk-1
          x21=x32
          x32=xk(k+1)-xk(k)
          do j=1,4
            do i=1,4
               y32=Kk(i,j,k+1)-Kk(i,j,k)
               y21=Kk(i,j,k)-Kk(i,j,k-1)
               Kkder(i,j,k)=bezier_cubic_deriv2(x21,x32,y21,y32)    
            end do
          end do  
        end do
        do j=1,4
           do i=1,4
              Kkder(i,j,nk)=Kkder(i,j,nk-1)
           end do  
        end do
        return
        end    

c-----------------------------------------------
c bezier_cubic_deriv_absormodf  
        subroutine bezier_cubic_deriv_absormodf(xk,Kk,nk,KKder)  
        real*8 xk(nk),Kk(4,4,nk),Kkder(4,4,nk)
        real*8 y(6,nk),yder(6,nk)
        integer i,j,k,nk,l
        real*8 x21,x32,y21,y32
        real*8 bezier_cubic_deriv2
                          
        x21=0.
        x32=xk(2)-xk(1) 
        y21=0.
        
        do k=1,nk
           y(1,k)=Kk(1,2,k) !etaQ 
           y(2,k)=Kk(1,3,k) !etaU
           y(3,k)=Kk(1,4,k) !etaV
           y(4,k)=Kk(2,3,k) !rhoV
           y(5,k)=Kk(2,4,k) !-rhoU
           y(6,k)=Kk(3,4,k) !rhoQ
        end do
        
        do l=1,6
           y32=y(l,2)-y(l,1)
           yder(l,1)=bezier_cubic_deriv2(x21,x32,y21,y32)  
        end do
        
        do k=2,nk-1
          x21=x32
          x32=xk(k+1)-xk(k)       
          do l=1,6
             y32=y(l,k+1)-y(l,k)         
             y21=y(l,k)-y(l,k-1)                     
             yder(l,k)=bezier_cubic_deriv2(x21,x32,y21,y32)  
          end do           
        end do
        do l=1,6
           yder(l,nk)=yder(l,nk-1)
        end do   
        
        do k=1,nk
           Kkder(1,1,k)=0.
           Kkder(1,2,k)=yder(1,k)
           Kkder(1,3,k)=yder(2,k)
           Kkder(1,4,k)=yder(3,k)
           
           Kkder(2,1,k)=yder(1,k)
           Kkder(2,2,k)=0.
           Kkder(2,3,k)=yder(4,k)
           Kkder(2,4,k)=yder(5,k)
                
           Kkder(3,1,k)=yder(2,k)
           Kkder(3,2,k)=-yder(4,k)
           Kkder(3,3,k)=0.
           Kkder(3,4,k)=yder(6,k)
                  
           Kkder(4,1,k)=yder(3,k)
           Kkder(4,2,k)=-yder(5,k)
           Kkder(4,3,k)=-yder(6,k)
           Kkder(4,4,k)=0.        
        end do   
        
        return
        end    

c----------------------------------------------------------------------
c bezier_cubic_deriv_absormodfini  
          subroutine bezier_cubic_deriv_absormodfini(xk,Kk,ni,nk,KKder)  
        real*8 xk(nk),Kk(4,4,nk),Kkder(4,4,nk)
        real*8 y(6,nk),yder(6,nk)
        integer i,j,k,nk,l,ni
        real*8 x21,x32,y21,y32
        real*8 bezier_cubic_deriv2
               
        do k=1,ni-1       
          do j=1,4
             do i=1,4
                Kkder(i,j,k)=0.
             end do 
          end do             
        end do 
        
        if(ni .le. 1)then
           x21=0.
           y21=0.
        else   
           x21=xk(ni)-xk(ni-1) 
        end if
        x32=xk(ni+1)-xk(ni) 
        
        if(ni .gt. 1)then
           do k=ni-1,nk
             y(1,k)=Kk(1,2,k) !etaQ 
             y(2,k)=Kk(1,3,k) !etaU
             y(3,k)=Kk(1,4,k) !etaV
             y(4,k)=Kk(2,3,k) !rhoV
             y(5,k)=Kk(2,4,k) !-rhoU
             y(6,k)=Kk(3,4,k) !rhoQ
          end do
        else
          do k=ni,nk
             y(1,k)=Kk(1,2,k) !etaQ 
             y(2,k)=Kk(1,3,k) !etaU
             y(3,k)=Kk(1,4,k) !etaV
             y(4,k)=Kk(2,3,k) !rhoV
             y(5,k)=Kk(2,4,k) !-rhoU
             y(6,k)=Kk(3,4,k) !rhoQ
          end do
        end if
        
        do l=1,6
           if(ni .gt. 1)then
              y21=y(l,ni)-y(l,ni-1)
           end if
           y32=y(l,ni+1)-y(l,ni)
           yder(l,ni)=bezier_cubic_deriv2(x21,x32,y21,y32)  
        end do
        
        do k=ni+1,nk-1
          x21=x32
          x32=xk(k+1)-xk(k)       
          do l=1,6
             y32=y(l,k+1)-y(l,k)         
             y21=y(l,k)-y(l,k-1)                     
             yder(l,k)=bezier_cubic_deriv2(x21,x32,y21,y32)  
          end do           
        end do
        do l=1,6
           yder(l,nk)=yder(l,nk-1)
        end do   
        
        do k=ni,nk
           Kkder(1,1,k)=0.
           Kkder(1,2,k)=yder(1,k)
           Kkder(1,3,k)=yder(2,k)
           Kkder(1,4,k)=yder(3,k)
           
           Kkder(2,1,k)=yder(1,k)
           Kkder(2,2,k)=0.
           Kkder(2,3,k)=yder(4,k)
           Kkder(2,4,k)=yder(5,k)
                
           Kkder(3,1,k)=yder(2,k)
           Kkder(3,2,k)=-yder(4,k)
           Kkder(3,3,k)=0.
           Kkder(3,4,k)=yder(6,k)
                  
           Kkder(4,1,k)=yder(3,k)
           Kkder(4,2,k)=-yder(5,k)
           Kkder(4,3,k)=-yder(6,k)
           Kkder(4,4,k)=0.        
        end do   
        
        return
        end    
        
c-----------------------------------------------
c It evaluates Ek & Fk1 cubic-bezier control points
c It uses function bezier_cubic_deriv

        subroutine bezier_cubic_coeff(xk,yk,nk,Ek,Fk1)  
        real*8 xk(*),yk(*),Ek(*),Fk1(*)
        integer i,j,nk,n,jold
        real*8 xki_1,xki,xki1,xki2,xkiold,eps
        real*8 yki_1,yki,yki1,yki2
        real*8 hk3
        real*8 yder,yder1,bezier_cubic_deriv
        
        eps=1.d-8
                   
        do j=1,nk
             if(j .eq. 1)then 
                xki=xk(j) 
                xki1=xk(j+1) 
                xki2=xk(j+2)
                yki=yk(j) 
                yki1=yk(j+1) 
                yki2=yk(j+2)
                xki_1=xki  
                yki_1=yki
             endif 
             if(j .ge. 2 .and. j .le. nk-2)then 
                xki_1=xk(j-1) 
                xki=xk(j)
                xki1=xk(j+1) 
                xki2=xk(j+2)
                yki_1=yk(j-1) 
                yki=yk(j) 
                yki1=yk(j+1) 
                yki2=yk(j+2)
             endif   
             if(j .ge. 2 .and. j .gt. nk-2)then 
                xki_1=xk(j-1) 
                xki=xk(j) 
                xki1=xk(j+1) 
                yki_1=yk(j-1) 
                yki=yk(j) 
                yki1=yk(j+1) 
                xki2=xki1
                yki2=yki1
             endif   
                  
             if (j .eq. 1)then   
                yder=bezier_cubic_deriv(xki_1,xki,xki1,yki_1,yki,yki1)
                yder1=bezier_cubic_deriv(xki,xki1,xki2,yki,yki1,yki2)
             else
                yder=yder1
                yder1=bezier_cubic_deriv(xki,xki1,xki2,yki,yki1,yki2)
             endif   
             hk3=(xki1-xki)/3.
             Ek(j)=yki+hk3*yder
             Fk1(j)=yki1-hk3*yder1

        end do
        
        return
        end  

c-----------------------------------------------
c It evaluates C quadratic-bezier control point
c It uses function bezier_cubic_deriv

        subroutine bezier_quadratic_coeff(xk,yk,nk,Ck)  
        real*8 xk(*),yk(*),Ck(*)
        integer i,j,nk,n,jold
        real*8 xki_1,xki,xki1,xki2,xkiold,eps
        real*8 yki_1,yki,yki1,yki2
        real*8 hk,C0,C1
        real*8 yder,yder1,bezier_cubic_deriv
        
        eps=1.d-8
                   
        do j=1,nk
             if(j .eq. 1)then 
                xki=xk(j) 
                xki1=xk(j+1) 
                xki2=xk(j+2)
                yki=yk(j) 
                yki1=yk(j+1) 
                yki2=yk(j+2)
                xki_1=xki  
                yki_1=yki
             endif 
             if(j .ge. 2 .and. j .le. nk-2)then 
                xki_1=xk(j-1) 
                xki=xk(j)
                xki1=xk(j+1) 
                xki2=xk(j+2)
                yki_1=yk(j-1) 
                yki=yk(j) 
                yki1=yk(j+1) 
                yki2=yk(j+2)
             endif   
             if(j .ge. 2 .and. j .gt. nk-2)then 
                xki_1=xk(j-1) 
                xki=xk(j) 
                xki1=xk(j+1) 
                yki_1=yk(j-1) 
                yki=yk(j) 
                yki1=yk(j+1) 
                xki2=xki1
                yki2=yki1
             endif   
                  
             if (j .eq. 1)then   
                yder=bezier_cubic_deriv(xki_1,xki,xki1,yki_1,yki,yki1)
                yder1=bezier_cubic_deriv(xki,xki1,xki2,yki,yki1,yki2)
             else
                yder=yder1
                yder1=bezier_cubic_deriv(xki,xki1,xki2,yki,yki1,yki2)
             endif   
             hk=(xki1-xki)
             C0=yki-hk*yder   !antes yki+hk*yder     !!!!Atencion: cambio de signo porque en el input la red en tau va de arriba a abajo...
             C1=yki1+hk*yder1 !antes yki1-hk*yder1   !!!!Atencion: cambio de signo porque en el input la red en tau va de arriba a abajo...
             Ck(j)=(C0+C1)/2.

        end do
        
        return
        end  

c-----------------------------------------------
c interpolates using cubic bezier an array yk defined at xk and returns y (i.e. values at x)
c nk=n_elements(xk), n=n_elements(x)
c uses function bezier_cubic_deriv

        subroutine bezier_cubic_interp(xk,yk,nk,x,y,n,Eb,Fb)  
        real*8 xk(*),yk(*),x(*),y(*),Eb(*),Fb(*)
        integer i,j,nk,n,jold
        real*8 xki_1,xki,xki1,xki2,xkiold,eps
        real*8 yki_1,yki,yki1,yki2
        real*8 yder,yder1,hk3,Ebi,Fbi,bezier_cubic_deriv
        real*8 u,u1

        jold=0
        xkiold=-1000.
        eps=1.e-8
                   
        do i=1,n
          call locate(xk,nk,x(i),j)
          if(j .gt. nk-1) j=nk-1
          if(j .lt. 1)j=1
          
          if(j .ne. jold)then
             if(j .eq. 1)then 
                xki=xk(j) 
                xki1=xk(j+1) 
                xki2=xk(j+2)
                yki=yk(j) 
                yki1=yk(j+1) 
                yki2=yk(j+2)
                xki_1=xki  
                yki_1=yki
             endif 
             if(j .ge. 2 .and. j .le. nk-2)then 
                xki_1=xk(j-1) 
                xki=xk(j)
                xki1=xk(j+1) 
                xki2=xk(j+2)
                yki_1=yk(j-1) 
                yki=yk(j) 
                yki1=yk(j+1) 
                yki2=yk(j+2)
             endif   
             if(j .ge. 2 .and. j .gt. nk-2)then 
                xki_1=xk(j-1) 
                xki=xk(j) 
                xki1=xk(j+1) 
                yki_1=yk(j-1) 
                yki=yk(j) 
                yki1=yk(j+1) 
                xki2=xki1
                yki2=yki1
             endif   
                  
             if (dabs(xki_1-xkiold) .gt. eps)then   
                yder=bezier_cubic_deriv(xki_1,xki,xki1,yki_1,yki,yki1)
                yder1=bezier_cubic_deriv(xki,xki1,xki2,yki,yki1,yki2)
             else
                yder=yder1
                yder1=bezier_cubic_deriv(xki,xki1,xki2,yki,yki1,yki2)
             endif   
             hk3=(xki1-xki)/3.
             Ebi=yki+hk3*yder
             Fbi=yki1-hk3*yder1
          endif     
          u=(x(i)-xki)/(xki1-xki)
          u1=1.-u 
          y(i)=u1*u1*u1*yki+u*u*u*yki1+3.*u*u1*u1*Ebi+3.*u*u*u1*Fbi
          Eb(i)=Ebi
          Fb(i)=Fbi
          jold=j
          xkiold=xki
        end do
        
        return
        end  

c-----------------------------------------------

c derivapruebaarray
       subroutine derivaarray(x,y,n,yder)
       real*8 x(*),y(*),yder(*)
       integer n,i
       
       do i=2,n-1
         yder(i)=((y(i+1)-y(i-1)))/(x(i+1)-x(i-1))
       end do
       yder(1)=(y(2)-y(1))/(x(2)-x(1))
       yder(n)=(y(n)-y(n-11))/(x(n)-x(n-1))
       return
       end

c-----------------------------------------------
c derivapruebaarray4
       subroutine derivaarray44(x,y44,n,yder44)
       real*8 x(*),y44(4,4,n),yder44(4,4,n),yi(n),yider(n)
       integer n,i,j,k
       
       do i=1,4
          do j=1,4
             do k=1,n
                yi(k)=y44(i,j,k)
             end do
             call derivaarray(x,yi,n,yider)
             do k=1,n
                yder44(i,j,k)=yider(k)
             end do
         end do
       end do  
       return
       end   