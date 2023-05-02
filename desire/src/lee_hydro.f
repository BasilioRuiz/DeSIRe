c read H populations model1 if exist (ihidro1=1)

         subroutine lee_hydro(mod_hyd_in1,ihydro1,ntau,tau,t,pe,hydro1)

         implicit real*4 (a-h,o-z)
         include 'PARAMETER'
         logical there
         real*4  hydro1(6,kt),tau(*),t(*),pe(*)
         character*(*) mod_hyd_in1
         character*1 cabecera  !ignora el resto de la linea

c        11/11/20 brc: Si se suministran las poblaciones de H desde fuera
c        (tipicamente "observed.hyd" para sintesis y "guess.hyd" para
c        inversion) entonces se toman dichas poblaciones.
         inquire(file=mod_hyd_in1,exist=there)
         if(there)then
            ihydro1=1
            call error(KPARA,'','Using logtau scale in RH with'
     &      //                  ' H populations file: '//mod_hyd_in1)
            open(8,file=mod_hyd_in1,status='old',err=900)
            read(8,'(a)')cabecera
            do i=1,ntau
               read(8,*,err=901)tauj,(hydro1(j,i),j=1,6)
            end do
            close(8)
         else
            ihydro1=0
            call error(KPARA,'','Evaluating LTE H populations to use'
     &      //         ' logtau scale in RH')
            call calculate_hydro(mod_hyd_in1,ntau,tau,t,pe,hydro1)
         end if

         return

900      call error(KSTOP,'lee_hydro','Error opening H populations file\n'
     &   //               ' File: '//mod_hyd_in1)

901      call error(KSTOP,'lee_hydro','Error reading H populations file\n'
     &   //               ' File: '//mod_hyd_in1)

         return
         end

c-----------------------------------------------------------------------------

         subroutine calculate_hydro(mod_hyd_in1,ntau,tau,t,pe,hydro1)

         implicit real*4 (a-h,o-z)
         include 'PARAMETER'
         real*4 hydro1(6,kt),popH(6),tau(*),t(*),pe(*)
         character*(*) mod_hyd_in1

         do i=1,ntau
            call hpopulations(t(i),pe(i),popH)
            do j=1,6
               hydro1(j,i)=popH(j)
            end do
         end do

c        11/11/20 epm: No es necesario las poblaciones de H en disco.
c        call write_hydro(mod_hyd_in1,ntau,tau,hydro1)

         return
         end

c-----------------------------------------------------------------------------

         subroutine write_hydro(filename,ntau,tau,hydro1)

         implicit real*4 (a-h,o-z)
         include 'PARAMETER'
         real*4 hydro1(6,kt),tau(*)
         character*(*) filename
         
         open(8,file=filename)
         write(8,100)'*   Index','nH(1)','nH(2)','nH(3)','nH(4)','nH(5)','nH+'
         do i=1,ntau
            write(8,101)tau(i),(hydro1(j,i),j=1,6)
         end do
         close(8)

100      format(a,6(11x,a))
101      format(1x,1pe15.8,6(1x,1pe15.8))

         return
         end

c-----------------------------------------------------------------------------

c rutina que copia las poblaciones de hidrogeno en otra variable 
         subroutine copia_hydro(ntau,hidro1in,hidro1out,imodel2,hidro2in,hidro2out)
         
         implicit real*4 (a-h,o-z)
         include 'PARAMETER'
         real*4 hydro1in(6,kt),hidro1out(6,kt)
         real*4 hydro2in(6,kt),hidro2out(6,kt)

         do j=1,6
            do i=1,ntau
               hidro1out(j,i)=hydro1in(j,i)
            end do
         end do
         if(imodel2 .eq. 1)then
            do j=1,6
               do i=1,ntau
                  hidro2out(j,i)=hydro2in(j,i)
               end do
            end do
         end if

         return
         end      

c-----------------------------------------------------------------------------
