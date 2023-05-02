c
c File          _______ : error.f
c Description   _______ : Generic error routines.
c Project       _______ : DeSIRe
c Creation date _______ : 05/05/20
c Author        _______ : epm@iac.es
c


c_____________________________________________________________________________
c
c  Routine: error()
c
c  Generic error routine.
c  @param  level error level.
c  @param  routine routine name.
c  @param  text message to show.
c
c_____________________________________________________________________________

      subroutine error( level, routine, text )

      implicit none
      character*(*) routine, text
      integer*4 level
      integer*4 flagv, flagquiet
      integer*4 KSTOP, KERR, KWARN, KPARA, KLINE, KLITE, KTEXT
      parameter (KSTOP=1, KERR=2, KWARN=3, KPARA=4, KLINE=5, KLITE=6, KTEXT=7)

c     05/05/20 epm: Common to save command line flags.
      common/commandline/flagv,flagquiet

c     Format '$': no print the default new line, '/': print a new line.
      if (level .eq. KTEXT) then
         if (flagquiet .ne. 1) then
            call writing(text)
         end if

      else if (level .eq. KLITE) then
         if (flagquiet .ne. 1) then
            call writing(text)
            write(*, *)
         end if

      else if (level .eq. KLINE) then
         if (flagquiet .ne. 1) then
            write(*,'(a,$)') " "
            call writing(text)
            write(*, *)
         end if

      else if (level .eq. KPARA) then
         if (flagquiet .ne. 1) then
            write(*,'(a,$)') " "
            call writing(text)
            write(*, '(/)')
         end if

      else if (level .eq. KWARN) then
         if (flagquiet .ne. 1) then
            write(*, '(a,a,/,a,$)') "-WARNING in routine ", routine, " "
            call writing(text)
            write(*, '(/)')
         end if

      else if (level .eq. KERR) then
         if (flagquiet .ne. 1) then
            write(*, '(a,a,/,a,$)') "-ERROR in routine ", routine, " "
            call writing(text)
            write(*, '(/,a,/)') " Trying to continue....."
         end if

      else if (level .eq. KSTOP) then
         write(*, '(/,a,a,/,a,$)')"-TERMINATING ERROR in routine ",routine," "
         call writing(text)
         write(*, '(/,a,/)') " Exiting....."
         stop

      else
         write(*, '(/,a,a,$)') "-UNKNOWN ERROR LEVEL in routine ", routine
         write(*, '(a,/,a,$)') " (please, report this error)", " "
         call writing(text)
         write(*, '(/,a,/)') " Exiting....."
         stop

      end if

      return
      end


c_____________________________________________________________________________
c
c  Routine: writing()
c
c  Write the text looking for new lines.
c  @param  text message to show.
c
c_____________________________________________________________________________

      subroutine writing(text)

      implicit none
      character*(*) text
      integer*4 nlen, iof, nl

c     Trailing blanks are all blanks located at the end of a line, so we
c     can go through the string until len_trim() position.
      nlen = len_trim(text)
      iof = 1
      nl = index(text(iof:nlen), "\n")
      do while (nl .gt. 0)
         write(*,'(a)') text(iof:iof+nl-2)
         iof = iof+nl+1
         if (iof .le. nlen) then
            nl = index(text(iof:nlen), "\n")
         else
            nl = 0
         end if
      end do
      write(*,'(a,$)') text(iof:nlen)

      return
      end


c_____________________________________________________________________________
