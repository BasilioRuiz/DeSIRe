c
c File          _______ : etable.f
c Description   _______ : Program for calculating errors in tables.
c Project       _______ : DeSIRe
c Creation date _______ : 12/12/22
c Author        _______ : epm@iac.es
c


c_____________________________________________________________________________
c
c  Program: etable
c
c  Main function.
c
c_____________________________________________________________________________

      program etable

      implicit none
      integer*4 bilineal, opentable, closetable  !C functions
      integer*4 openfile, writefile, closefile   !C functions

      interface
         integer*4 function template( x, y, f )
         real*4 x, y, f
         end function
         integer*4 function Pg_T_Pe( x, y, f )
         real*4 x, y, f
         end function
      end interface

      procedure(template), pointer:: fpointer => NULL()  !function pointer

      character s1*20, s2*20
      character filename*64, file_err*72, file_cal*72, file_int*72
      integer*4 codes(72), ncodes
      integer*4 which, itable, ifile_err, ifile_cal, ifile_int
      integer*4 maxboxes, nboxes, nelem
      integer*4 ninfo, ncols, nrows, nxdiv, nydiv
      integer*4 ier, i, j, k, n
      integer*4 ilog  !not logarithmic table values = 0, otherwise = 1
      real*4    xini, xdelta, yini, ydelta
      real*4    x, y, f, p, e, emax, xmax, ymax

      parameter (maxboxes = 100, ninfo = 6)
      real*4    vheader(2 + (maxboxes * ninfo))

      write(*,*) ' __________________________________________________ '
      write(*,*) '|                                                  |'
      write(*,*) '|        Error of the interpolation tables         |'
      write(*,*) '|__________________________________________________|'
      write(*,*) ''
      write(*,*) '  0) table.demo1 '
      write(*,*) '  1) table.demo2 '
      write(*,*) '  2) table.Pg_T-Pe '
      write(*,*) ''
      write(*,'(a,$)') '  Table number: '
      read(*,'(i2)',iostat=ier) which
      if (ier .ne. 0) then
         write(*,'(/a/)') '  ** Error: Invalid table number'
         stop
      end if
      write(*,'(a,$)') '  Divisions between columns (X axis): '
      read(*,'(i4)',iostat=ier) nxdiv
      if (ier.ne.0 .or. nxdiv.le.0) then
         write(*,'(/a/)') '  ** Error: Invalid sampling'
         stop
      end if
      write(*,'(a,$)') '  Divisions between rows    (Y axis): '
      read(*,'(i4)',iostat=ier) nydiv
      if (ier.ne.0 .or. nydiv.le.0) then
         write(*,'(/a/)') '  ** Error: Invalid sampling'
         stop
      end if

c     Select the table.
      select case (which)
         case (0)
            fpointer => template
            filename = 'table.demo1'
            ilog     = 0

         case (1)
            fpointer => template
            filename = 'table.demo2'
            ilog     = 0

         case (2)
            fpointer => Pg_T_Pe
            filename = 'table.Pg_T-Pe'
            ilog     = 1

         case default
            write(*,'(/a/)') '  ** Error: Invalid table number'
            stop
      end select

c     Open the interpolated table.
      ncodes = len_trim(filename)
      call toascii(trim(filename), ncodes, codes)
      itable = opentable(ncodes, codes, vheader)
      if (itable .eq. -1) then
         write(*,'(/a/)') '  ** Error: Unable to open the table'
         stop
      else if (itable .eq. -2) then
         write(*,'(/a/)') '  ** Error: Unable to read the table'
         stop
      else if (itable . eq. -3) then
         write(*,'(/a/)') '  ** Error: Unable to alloc memory'
         stop
      end if

      nboxes = int(vheader(1))
      nelem = 2 + (nboxes * ninfo)

c     Sample all of the boxes overwriting the header.
      do k = 1, nboxes
         ncols   = int(vheader( (k-1)*ninfo + 3 ))
         xini    =     vheader( (k-1)*ninfo + 4 )
         xdelta  =     vheader( (k-1)*ninfo + 5 )
         nrows   = int(vheader( (k-1)*ninfo + 6 ))
         yini    =     vheader( (k-1)*ninfo + 7 )
         ydelta  =     vheader( (k-1)*ninfo + 8 )

         ncols = (ncols-1) * nxdiv  !open interval in the right: [xini, xfin)
         nrows = (nrows-1) * nydiv  !open interval in the right: [yini, yfin)
         xdelta = xdelta / nxdiv
         ydelta = ydelta / nydiv
         vheader( (k-1)*ninfo + 3 ) = real(ncols)
         vheader( (k-1)*ninfo + 5 ) = xdelta
         vheader( (k-1)*ninfo + 6 ) = real(nrows)
         vheader( (k-1)*ninfo + 8 ) = ydelta
      end do

c     Define the output files.
      file_err = trim(filename)//'_err'  !file with the error values
      file_cal = trim(filename)//'_cal'  !file with the calculated values
      file_int = trim(filename)//'_int'  !file with the interpolated values
      write(*,*) ''
      write(*,*) ' Creating error file       : "', trim(file_err), '"'
      write(*,*) ' Creating calculated file  : "', trim(file_cal), '"'
      write(*,*) ' Creating interpolated file: "', trim(file_int), '"'
      write(*,*) ''

c     Open the output files.
      ncodes = len_trim(file_err)
      call toascii(trim(file_err), ncodes, codes)
      ifile_err = openfile(ncodes, codes)
      if (ifile_err .eq. -1) then
         write(*,'(/a/)') '  ** Error: Unable to create the error file'
         stop
      end if
      ncodes = len_trim(file_cal)
      call toascii(trim(file_cal), ncodes, codes)
      ifile_cal = openfile(ncodes, codes)
      if (ifile_cal .eq. -1) then
         write(*,'(/a/)') '  ** Error: Unable to create the calculated file'
         stop
      end if
      ncodes = len_trim(file_int)
      call toascii(trim(file_int), ncodes, codes)
      ifile_int = openfile(ncodes, codes)
      if (ifile_int .eq. -1) then
         write(*,'(/a/)') '  ** Error: Unable to create the interpolated file'
         stop
      end if

c     Write the new header in the output files.
      if (writefile(ifile_err, nelem, vheader) .ne. nelem) then
         write(*,'(/a/)') '  ** Error: Unable to write the error file'
         stop
      end if
      if (writefile(ifile_cal, nelem, vheader) .ne. nelem) then
         write(*,'(/a/)') '  ** Error: Unable to write the calculated file'
         stop
      end if
      if (writefile(ifile_int, nelem, vheader) .ne. nelem) then
         write(*,'(/a/)') '  ** Error: Unable to write the interpolated file'
         stop
      end if

c     Go across the sampled table.
      do k = 1, nboxes
         ncols   = int(vheader( (k-1)*ninfo + 3 ))
         xini    =     vheader( (k-1)*ninfo + 4 )
         xdelta  =     vheader( (k-1)*ninfo + 5 )
         nrows   = int(vheader( (k-1)*ninfo + 6 ))
         yini    =     vheader( (k-1)*ninfo + 7 )
         ydelta  =     vheader( (k-1)*ninfo + 8 )

         emax = 0
         do i = 1, nrows
            y = yini + ydelta*(i-1)
            do j = 1, ncols
               x = xini + xdelta*(j-1)

c              Calculate the real value.
               ier = fpointer(x, y, f)  !f = real value
               if (ier .eq. -1) then
                  write(*,'(/a/)') '  ** Error: Table with NaN values'
                  stop
               end if

c              Calculate the interpolated value.
               ier = bilineal(itable, x, y, p)  !p = interpolated value
               if (ier .eq. -1) then
                  write(*,'(/a/)') '  ** Error: Table not opened'
                  stop
               else if (ier .eq. -2) then
                  write(s1,*) xini
                  write(s2,*) xini + (xdelta * (ncols-1))
                  write(*,'(/a/)') '  ** Error: X coordinate out of range ['
     &               //trim(adjustl(s1))//', '//trim(adjustl(s2))//')'
                  stop
               else if (ier .eq. -3) then
                  write(s1,*) yini
                  write(s2,*) yini + (ydelta * (nrows-1))
                  write(*,'(/a/)') '  ** Error: Y coordinate out of range ['
     &               //trim(adjustl(s1))//', '//trim(adjustl(s2))//')'
                  stop
               else if (ier .eq. -4) then
                  write(*,'(/a/)') '  ** Error: Unable to read the table'
                  stop
               end if

c              Calculate the relative error.
               if (ilog .ne. 0) then  !logarithmic table values
                  e = abs(1 - 10**(abs(p-f)))
               else                   !no logarithmic table values
                  if (abs(f) .lt. 1e-6) then
                     e = abs(p - f)   !absolute error if f = 0
                  else
                     e = abs(p - f) / abs(f)
                  end if
               end if
               if (e .gt. emax) then
                  emax = e
                  xmax = x
                  ymax = y
               end if

c              Write the real, interpolated and error value.
               n = 1
               if (writefile(ifile_cal, n, f) .ne. n) then
                  write(*,'(/a/)') '  ** Error: Unable to write the real value'
                  stop
               end if
               if (writefile(ifile_int, n, p) .ne. n) then
                  write(*,'(/a/)') '  ** Error: Unable to write the'
     &            //               ' interpolated value'
                  stop
               end if
               if (writefile(ifile_err, n, e) .ne. n) then
                  write(*,'(/a/)') '  ** Error: Unable to write the error value'
                  stop
               end if

            end do  !end do columns
         end do  !end do rows
         write(*,'(2x,a,i3,a,e10.3,$)') 'Box', k, ' - maximum error = ', emax
         write(*,'(2x,a,e10.3,a,e10.3,a)') 'in (', xmax, ',', ymax, ')'
      end do

c     Close the output files.
      if (closefile(ifile_err) .ne. 0) then
         write(*,'(/a/)') '  ** Error: Unable to close the error file'
         stop
      end if
      if (closefile(ifile_cal) .ne. 0) then
         write(*,'(/a/)') '  ** Error: Unable to close the calculated file'
         stop
      end if
      if (closefile(ifile_int) .ne. 0) then
         write(*,'(/a/)') '  ** Error: Unable to close the interpolated file'
         stop
      end if

c     Close the interpolated table.
      if (closetable(itable) .ne. 0) then
         write(*,'(/a/)') '  ** Error: Unable to close the table'
         stop
      end if

c     Happy ending.
      write(*,*) ''

      stop
      end


c_____________________________________________________________________________
