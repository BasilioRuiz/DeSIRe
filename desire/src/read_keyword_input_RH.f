c read_keyword_input_RH.f  reads the RH input file keyword.input ; BRC Jun 19 2017
c looks for prompt (must be follwed by "=" and then by filename_prompt in keyword.input

        subroutine read_keyword_input_RH(prompt,filename_prompt)
        implicit real*4 (a-h,o-z)
        character*100 nombre,linea
        character*(*) prompt,filename_prompt
        integer len_pr,ican,num,nxx,iline,ifi,i,ieq,ieq1,ieq2
        
        nombre="keyword.input"
        ican=50
        open(ican,file=nombre,status='old',err=999)
        
c contamos las lineas
        num=0
        nxx=0
        iline=1
        ifi=0
        i=100
        len_pr=len(prompt)

        do while(iline.ne.2000 .and. i.eq.100)
           read(ican,'(a)',err=997,END=997)linea 
           call busco ('#',linea,ifi,ieq) 
           if(ieq .eq. 100)then
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
           iline=iline+1
        end do
        ini=i
        ifi=ini
        call busco('=',linea,ifi,ieq)
        call nobusco(' ',linea,ieq+1,ieq1)
        call busco(' ',linea,ieq1+2,ieq2) 

        filename_prompt=linea(ieq1:ieq2-1)
        print*,''
        print*,prompt,' = ',trim(filename_prompt)

        close (ican)
        return

999     print*,''
        print*,'STOP: The file "keyword.input" does NOT exist.'
        stop ' '

997     print*,''
        print*,'STOP: "keyword.input" does not contain ',trim(prompt)
        stop ' '

        end
