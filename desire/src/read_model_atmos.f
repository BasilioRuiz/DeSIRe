c read_keyword_input_RH.f  reads the RH input file keyword.input ; BRC Jun 19 2017
c looks for prompt (must be follwed by "=" and then by filename_prompt in keyword.input
c ________________________________________________________________________

        subroutine read_model_atmos(name,ndepth,depth,temp,Nelec,Vz,Vtur)
        include 'PARAMETER'
        implicit real*4 (a-h,o-z)
        character*(*) name
        character*100 linea
        integer ndepth,ini,ifi,ifound
        real*4 depth(kt),temp(kt),Nelec(kt),Vz(kt),Vtur(kt)

        ican=50
        open(ican,file=name,status='old',err=997)

c       Contamos las lineas.
        iline=1
        ifound=0
        ifi=20

        do while(iline.lt.100 .and. ifound .eq. 0)
           read(ican,'(a)',err=997)linea 
           ini=1
           call buscon('Ndep',4,linea,ini,ifi,ifound)
           iline=iline+1
        end do
        if(iline .eq. 100)then
           close(ican)
           open(ican,file=name,status='old',err=997)
           iline=1
           ifound=0
           ifi=20
           do while(iline.lt.100 .and. ifound .eq. 0)
              read(ican,'(a)',err=997)linea 
              ini=1
              call buscon('NDEP',4,linea,ini,ifi,ifound)
              iline=iline+1
           end do
        end if
        if(iline .eq. 100)then 
          print*,'The string Ndep (nor NDEP) does not appear inside file ',name
          print*,'STOP at read_model_atmos'
          stop
        end if

        read(ican,'(i4)',err=997)ndepth

        iline=1
        ifound=0
        ifi=100
        do while(iline.lt.100 .and. ifound .eq. 0)
           read(ican,'(a)',err=997)linea 
           ini=1
           call buscon('Temp',4,linea,ini,ifi,ifound)
           iline=iline+1
        end do
        if(iline .eq. 100)then
           close(ican)
           open(ican,file=name,status='old',err=997)
           iline=1
           ifound=0
           ifi=100
           do while(iline.lt.100 .and. ifound .eq. 0)
              read(ican,'(a)',err=997)linea 
              ini=1
              call buscon('TEMP',4,linea,ini,ifi,ifound)
              iline=iline+1
           end do
        end if

        if(iline .eq. 100)then 
          print*,'The string Temp (nor TEMP) does not appear inside file ',name
          print*,'STOP at read_model_atmos'
          stop
        end if
        do idep=1,ndepth
           read(ican,*,err=997)depth(idep),temp(idep),Nelec(idep),Vz(idep),Vtur(idep)
        end do
        close (ican)
        return

997     print*,' '
        print*,'STOP: error reading ',name
        stop
        return
        end
