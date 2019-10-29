
	subroutine chachi(num,cha,nc)
        implicit real*4 (a-h,o-z)

c	nc: Dimension de "cha" en la rutina que llama a "chachi".
c	"cha" se rellena por la izquierda con blancos.
c	La rutina chachi convierte el entero num en un character cha	

	character*(*) cha
	integer num,k,j,nc

	if (nc.gt.10) then   !Esto no sirve para nada
	   print*,'El numero de caracteres debe ser menor que 10'
	   stop
	end if
        
	
        kk=10**nc 
	k=num-(num/kk)*kk      !k es el entero formado por las ultimas
			       !nc cifras de num
	do i=nc-1,0,-1
	   j=10**i
	   n=k/j
	   k=k-j*n
	   cha(nc-i:nc-i)=char(n+48)
	end do

	return
	end
c *****************************************************************************

	subroutine quitaex0(cha) !quita la extension detras del punto
        implicit real*4 (a-h,o-z)

	character*(*) cha
	character*100 chacho

	long=len(cha)

	i=1

	do while(i.lt.long.and.cha(i:i).ne.'.')
	   chacho(i:i)=cha(i:i)
	   i=i+1
	end do
	
	do j=i,long
	   chacho(j:j)=' '
	end do

	cha(1:long)=chacho(1:long)

	return
	end
c *****************************************************************************

	subroutine quitaex(cha)
        implicit real*4 (a-h,o-z)

	character*(*) cha
	character*100 chacho
	integer long,i,j

	long=len(cha)

	i=1

	do while(i.lt.long.and.cha(i:i).ne.'.'.and.cha(i:i).ne.'_')
	   chacho(i:i)=cha(i:i)
	   i=i+1
	end do
	
	do j=i,long
	   chacho(j:j)=' '
	end do

	cha(1:long)=chacho(1:long)

	return
	end



*******************************************************************************
	subroutine quitaex2(cha,ixt)
        implicit real*4 (a-h,o-z)

	integer ixt,long,i,j
	character*(*) cha
	character*100 chacho
        

	long=len(cha)

	i=1

	do while(i.lt.long.and.cha(i:i).ne.'.'.and.cha(i:i).ne.'_')
	   chacho(i:i)=cha(i:i)
	   i=i+1
	end do

	if(cha(i:i).eq.'_')then
             ixt=ichar(cha(i+1:i+1))-48
	else
   	     ixt=0
	endif
	
	do j=i,long
	   chacho(j:j)=' '
	end do

	cha(1:long)=chacho(1:long)

	return
	end

************************************************************
	subroutine quitaex3(cha,ixt)
        implicit real*4 (a-h,o-z)

	integer ixt,long,i,j
	character*(*) cha
	character*100 chacho
        

	long=len(cha)

	i=1

	do while(i.lt.long.and.cha(i:i).ne.'.')
	   chacho(i:i)=cha(i:i)
	   i=i+1
	end do

	if(cha(i:i).eq.'_')then
             ixt=ichar(cha(i+1:i+1))-48
	else
   	     ixt=0
	endif
	
	do j=i,long
	   chacho(j:j)=' '
	end do

	cha(1:long)=chacho(1:long)

	return
	end
*******************************************************
	subroutine extension3(cha,ext,n)  ! returns the extension of a file after . ot -
        implicit real*4 (a-h,o-z)

	character*(*) cha,ext
        integer long,i,n,j
        
	long=len(cha)

	i=1
	do while(i.lt.long.and.cha(i:i).ne.'.')
	   i=i+1
	end do
	j=i+1
	do while(j.lt.long.and.cha(j:j).ne.' ')
	   j=j+1
	end do
	
	ext(1:j-i)=cha(i+1:j)
	n=j-i-1
c        print*,'cha input=',cha
c        print*,'cha output=',ext
c        print*,'n=',n

	return
	end
*******************************************************		   	   
	subroutine concatena(cha1,cha2,cha)
        implicit real*4 (a-h,o-z)

	character*(*) cha1
	character*(*) cha2
	character*(*) cha
	integer long1, long2, long3,i,n1,n2

	long1 = len(cha1)
	long2 = len(cha2)
	long3 = len(cha)

	i=1
	do while(i.lt.long1.and.cha1(i:i).ne.' ')
	   i=i+1
	end do
	if(i .eq. long1 .and. cha1(i:i).ne.' ')i=i+1
	n1=i-1

	if (long3.le.n1) then !tampoco hace nada (dimens. se definen en el principal)
	   print*,'No puedo concatenar dos nombres en uno mas corto '
	   print*,'que el primero'
	   return
	end if

	
	i=1
	do while(i.lt.long2.and.cha2(i:i).ne.' ')
	   i=i+1
	end do
	if(i .eq. long2 .and. cha2(i:i).ne.' ')i=i+1
	n2=i-1
       
	cha(1:n1)=cha1(1:n1)
	if (n1+n2.le.long3) then
	   cha(n1+1:n1+n2) = cha2(1:n2)
	   cha(n1+n2+1:long3) = ' ' 
	else
	   cha(n1+1:long3) = cha2(1:long3-n1)
	   print*,'Hemos cortado la segunda cadena de caracteres: ',cha
	end if

	return
	end
	subroutine quitaend0(cha,n) !quita la extension de blancos al final
        implicit real*4 (a-h,o-z)

        integer i,n 
	character*(*) cha

	long=len(cha)

	i=long

	do while(i.gt.1.and.cha(i:i).eq.' ')
	   i=i-1
	end do
c	print*,'quitaend0:',cha(1:i),'--- i=',i
	n=i

	return
	end
	
	subroutine copying(model_multi1,RH_model)
	
	character*(*) model_multi1,RH_model
	character*100 modelcp1,modelcp2
	integer n,n1,istatus,system
	
        modelcp1='cp '//model_multi1
        call quitaend0(modelcp1,n)          
        modelcp1=modelcp1(1:n)//' '
        call quitaend0(RH_model,n1)
        modelcp2=modelcp1(1:n+1)//RH_model(1:n1)
        istatus=system(modelcp2(1:n+1+n1))
        if(istatus .ne. 0)then
          print*,'Error copying ',model_multi1,' in ',RH_model
          stop
        end if  
	return
	end