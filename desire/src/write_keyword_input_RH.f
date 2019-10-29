c write_keyword_input_RH.f  writes the RH input file keyword.input ; BRC Jun 19 2017
c looks for prompt (must be follwed by "=" and then by filename_prompt in keyword.input
c then writes keyword2.input changing the content of prompt by the integer icontent

        subroutine write_keyword_input_RH(prompt,icontent)
        implicit real*4 (a-h,o-z)
        character*100 nombre,linea,nombre2
        character*(*) prompt
        character*100 filename_prompt
        integer len_pr,ican,num,nxx,iline,ifi,i,ieq,ieq1,ieq2
        integer icontent
        
        nombre="keyword.input"
        nombre2="keyword2.input"
        ican=50
        open(ican,file=nombre,status='old',err=999)
        open(ican+1,status='new',file=nombre2)
c contamos las lineas
        num=0
        nxx=0
        iline=1
        ifi=1
        i=100
        len_pr=len(prompt)

        do while(iline.ne.2000 .and. i.eq.100)
           read(ican,'(a)',err=997,END=997)linea 
           call busco ('#',linea,ifi,ieq) 
           if(ieq .eq. 100)then               !if there isn't any #
             call busco (' ',linea,ifi,ieq2) 
             if(ieq2 .ne. 100)then
                 i=1
                 do while(linea(i:i+len_pr).ne. prompt .and.i.lt.100-len_pr)
                    i=i+1
                 end do
                 do while(linea(i:i).ne.'='.and.i.lt.100)
                   i=i+1
                end do
             endif   
           endif  
           if(ieq .ne. 100. or. i.eq.100)write(ican+1,101)linea
           iline=iline+1
        end do
        ini=i
        ifi=ini
        call busco('=',linea,ifi,ieq)
        call nobusco(' ',linea,ieq+1,ieq1)
        call busco(' ',linea,ieq1+2,ieq2) 

c        write(ican+1,102)' NRAYS =',icontent
        write(ican+1,102)prompt,' =',icontent

        do while(iline.ne.2000)
           read(ican,'(a)',END=996)linea 
           write(ican+1,101)linea
           iline=iline+1
        end do
        
996     close (ican+1)
        close (ican)
c        print*,'renombro keyword2.input'
        call system("mv keyword2.input keyword.input")
        return

999     print*,' '
        print*,'STOP: The file keyword.input does NOT exist.'
        print*,' '
        print*,'_______________________________________________________________________________'
        stop
        
        return
997     print*,'I do not find NRAYS'
        call system("rm keyword2.input")
        return
101     FORMAT(A78)
102     FORMAT(A7,A2,i3)
        return
        end
c write_keyword_input_RH.f  writes the RH input file keyword.input ; BRC Jun 19 2017
c looks for prompt (must be follwed by "=" and then by filename_prompt in keyword.input
c then writes keyword2.input changing the content of prompt by the integer icontent

        subroutine write_keyword_input_RHstr(prompt,stringcontent)
        implicit real*4 (a-h,o-z)
        character*100 nombre,linea,nombre2
        character*(*) prompt,stringcontent
        character*100 filename_prompt
        integer len_pr,ican,num,nxx,iline,ifi,i,ieq,ieq1,ieq2
c        integer icontent
        
        
        nombre="keyword.input"
        nombre2="keyword2.input"
        ican=50
        open(ican,file=nombre,status='old',err=999)
        open(ican+1,status='new',file=nombre2)
c contamos las lineas
        num=0
        nxx=0
        iline=1
        ifi=1
        i=100
        len_pr=len(prompt)

        do while(iline.ne.2000 .and. i.eq.100)
           read(ican,'(a)',err=997,END=997)linea 
           call busco ('#',linea,ifi,ieq) 
           if(ieq .eq. 100)then               !if there isn't any #
             call busco (' ',linea,ifi,ieq2) 
             if(ieq2 .ne. 100)then
                 i=1
                 do while(linea(i:i+len_pr).ne. prompt .and.i.lt.100-len_pr)
                    i=i+1
                 end do
                 do while(linea(i:i).ne.'='.and.i.lt.100)
                   i=i+1
                end do
             endif   
           endif  
           if(ieq .ne. 100. or. i.eq.100)write(ican+1,101)linea
           iline=iline+1
        end do
        ini=i
        ifi=ini
        call busco('=',linea,ifi,ieq)
        call nobusco(' ',linea,ieq+1,ieq1)
        call busco(' ',linea,ieq1+2,ieq2) 

c        write(ican+1,102)' NRAYS =',icontent
        write(ican+1,102)prompt,' = ',stringcontent

        do while(iline.ne.2000)
           read(ican,'(a)',END=996)linea 
           write(ican+1,101)linea
           iline=iline+1
        end do
        
996     close (ican+1)
        close (ican)
c        print*,'renombro keyword2.input'
        call system("mv keyword2.input keyword.input")
        return

999     print*,' '
        print*,'STOP: The file keyword.input does NOT exist.'
        print*,' '
        print*,'_______________________________________________________________________________'
        stop
        
        return
997     print*,'I do not find NRAYS'
        call system("rm keyword2.input")
        return
101     FORMAT(A78)
102     FORMAT(A12,A3,A12)
        return
        end
