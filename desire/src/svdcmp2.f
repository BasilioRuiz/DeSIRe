      subroutine svdcmp2(a,m,n,mp,np,w,v)

      implicit real*4 (a-h,o-z)
      include 'PARAMETER'   !mfitmax

      real*4 a(mfitmax,mfitmax),w(mfitmax),v(mfitmax,mfitmax)
      real*8 rv1(mfitmax)
      real*8 g,anorm,scale,s,f,h,ab,c,z,y,x,da(mfitmax,mfitmax)
      integer i,j,k,l,m
      
      g=0.0d0

      scale=0.0d0

      anorm=0.0d0

      do i=1,m
          do j=1,n
             da(i,j)=a(i,j)*1.d0
          enddo
      enddo

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
c          if (scale.ne.0.0) then
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

        w(i)=real(scale *g)

        g=0.0d0

        s=0.0d0

        scale=0.0d0

        if ((i.le.m).and.(i.ne.n)) then

          do 17 k=l,n

            scale=scale+dabs(da(i,k))

17        continue
c          if (scale.ne.0.0) then
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

	ab=abs(w(i))+abs(rv1(i))
        anorm=max(anorm,ab)

25    continue

      do 32 i=n,1,-1

        if (i.lt.n) then
c          if (g.ne.0.0) then
          if (dabs(g).gt.1.d-20) then

            do 26 j=l,n

              v(j,i)=real( (da(i,j)/da(i,l))/g )
26          continue

            do 29 j=l,n

              s=0.0

              do 27 k=l,n

                s=s+da(i,k)*v(k,j)

27            continue

              do 28 k=l,n

                v(k,j)=v(k,j)+real(s)*v(k,i)

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
c        if (g.ne.0.0) then
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
c            if ((dabs(rv1(l))+anorm).eq.anorm)  go to 2
            if ( dabs(rv1(l)) .lt. 1.d-20 )  go to 2
c            if ((abs(w(nm))*1.d0+anorm).eq.anorm)  go to 1
            if ( abs(w(nm)) .lt. 1.e-15)  go to 1

41        continue

1         c=0.0d0

          s=1.0d0

          do 43 i=l,k

            f=s*rv1(i)
c            if ((dabs(f)+anorm).ne.anorm) then
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

              v(nm,j)= real((x*c)+(z*s))

              v(nm,i)=real(-(x*s)+(z*c))
              
c              vnmj=(x*c)+(z*s)
c              vnmi=-(x*s)+(z*c)
c             if(vnmj+1 .eq. vnmj)print*,'svdcmp2 vnmj nm es nan!!!!',vnmj
c              if(vnmi+1 .eq. vnmi)print*,'svdcmp2 vnmi nm es nan!!!!',vnmi
c              if(vnmj+1 .eq. vnmj)stop
c             if(vnmi+1 .eq. vnmi)stop
c              if(v(nm,j)+1. .eq. v(nm,j))print*,'svdcmp2 v(nm,j) nm es nan!!!!',v(nm,j)
c              if(v(nm,i)+1. .eq. v(nm,i))print*,'svdcmp2 v(nm,i) nm es nan!!!!',v(nm,i)
c              if(v(nm,j)+1. .eq. v(nm,j))stop
c              if(v(nm,i)+1. .eq. v(nm,i))stop
              
45          continue

            z=dsqrt(f*f+h*h)

            w(j)=real(z)
c            if (z.ne.0.0) then
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

          w(k)=real(x)

48      continue

3       continue

49    continue

999     do i=1,m
          do j=1,n
             a(i,j)=real(da(i,j))
          enddo
       enddo  
       
      return

      end

