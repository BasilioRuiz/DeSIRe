	subroutine deriva4cuad(eta,deta,n,x)
        implicit real*4 (a-h,o-z)

c	include 'PARAMETER'  !solo por kt	
	integer n
	real*4 eta(4,4,n), deta(4,4,n), y(n), dy(n),x(*)
	integer indi(7),indj(7)
	data indi/1,1,1,1,2,2,3/
	data indj/1,2,3,4,3,4,4/
	

	do i=1,7
	   do ii=1,n 
	      y(ii) = eta(indi(i),indj(i),ii)	      	   
c	      print*,'deriva4cuad 15',ii,y(ii),indi(i),indj(i)
	   enddo

	   call derivacuad(y,dy,n,x)

	   do ii=1,n
c	      print*,'deriva4cuad 21',ii,dy(ii),indi(i),indj(i),ii
	      deta(indi(i),indj(i),ii) = dy(ii)
	   enddo
	enddo

	do i=1,n
	   deta(2,2,i) = deta(1,1,i)
	   deta(3,3,i) = deta(1,1,i)
	   deta(4,4,i) = deta(1,1,i)
	   deta(2,1,i) = deta(1,2,i)
	   deta(3,1,i) = deta(1,3,i)
	   deta(4,1,i) = deta(1,4,i)
	   deta(3,2,i) = -deta(2,3,i)
	   deta(4,2,i) = -deta(2,4,i)
	   deta(4,3,i) = -deta(3,4,i)
	enddo

	return
	end
