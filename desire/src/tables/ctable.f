c
c File          _______ : ctable.f
c Description   _______ : Program for creating interpolation tables.
c Project       _______ : DeSIRe
c Creation date _______ : 12/12/22
c Author        _______ : epm@iac.es
c


c_____________________________________________________________________________
c
c  Program: ctable
c
c  Main function.
c
c_____________________________________________________________________________

      program ctable

      implicit none
      integer*4 openfile, writefile, closefile  !C functions

      interface
         integer*4 function template( x, y, f )
         real*4 x, y, f
         end function
         integer*4 function Pg_T_Pe( x, y, f )
         real*4 x, y, f
         end function
      end interface

      procedure(template), pointer:: fpointer => NULL()  !function pointer

      character filename*64
      integer*4 codes(64), ncodes
      integer*4 which, it, itable1, itable2
      integer*4 maxboxes, nboxes, xyscan  !1 = scan in X; 2 = scan in Y
      integer*4 ifile, ncols, nrows, ninfo
      integer*4 ier, i, j, k, n
      real*4    xini, xdelta, yini, ydelta
      real*4    x, y, f

      parameter (maxboxes = 100, ninfo = 6)
      real*4    info(maxboxes, ninfo)  !in fortran: matrix(#rows, #cols)

      write(*,*) ' ___________________________________________________ '
      write(*,*) '|                                                   |'
      write(*,*) '|         Creation of interpolation tables          |'
      write(*,*) '|___________________________________________________|'
      write(*,*) ''
      write(*,*) '  0) Fake table => table.demo1 (1 box)               '
      write(*,*) '  1) Fake table => table.demo2 (3 boxes)             '
      write(*,*) '  2) Table lg(Pg) from T & lg(Pe) => table.Pg_T-Pe   '
      write(*,*) '  3) All                                             '
      write(*,*) ''
      write(*,'(a,$)') '  Table number: '
      read(*,'(i2)',iostat=ier) which
      if (ier .ne. 0) then
         write(*,'(/a/)') '  ** Error: Invalid table number'
         stop
      end if
      if (which .eq. 3) then
         itable1 = 0  !only in tests itable1 = 0, otherwise itable1 = 1
         itable2 = 2
      else
         itable1 = which
         itable2 = which
      end if

      do it = itable1, itable2  !foreach table

c        Table information.
c        real*4 constant notation: cc.ddE[+|-]ee  or  cc.dd
c        real*8 constant notation: cc.ddD[+|-]ee
         select case (it)
            case (0)
               fpointer => template
               filename = 'table.demo1'
               nboxes = 1
               xyscan = 1
               call boxinfo(5000., 8000., 1000., -4., 3., 1., info(1,:))

            case (1)
               fpointer => template
               filename = 'table.demo2'
               nboxes = 3
               xyscan = 2
               call boxinfo(4.,  8., 2., -2., -1., 1., info(1,:))
               call boxinfo(8., 14., 2.,  0.,  2., 1., info(2,:))
               call boxinfo(6., 10., 2.,  3.,  7., 1., info(3,:))

            case (2)
               fpointer => Pg_T_Pe
               filename = 'table.Pg_T-Pe'
               nboxes = 17
               xyscan = 2
               call boxinfo( 2490.,  3400.,  4., -4.0, -3.0,  0.0115, info(1,:))
               call boxinfo( 2490.,  6500.,  3., -3.0, -2.5,  0.020,  info(2,:))
               call boxinfo( 2490.,  9500.,  3., -2.5, -2.0,  0.020,  info(3,:))
               call boxinfo( 2490., 10500.,  4., -2.0, -1.5,  0.020,  info(4,:))
               call boxinfo( 2490., 12000.,  5., -1.5, -1.0,  0.025,  info(5,:))
               call boxinfo( 3500., 13000.,  5., -1.0, -0.5,  0.025,  info(6,:))
               call boxinfo( 3700., 10000.,  6., -0.5,  0.0,  0.025,  info(7,:))
               call boxinfo( 4150.,  6450.,  7.,  0.0,  0.5,  0.025,  info(8,:))
               call boxinfo( 4600.,  6150.,  8.,  0.5,  1.0,  0.025,  info(9,:))
               call boxinfo( 5300.,  6400., 17.,  1.0,  1.5,  0.025,  info(10,:))
               call boxinfo( 5950.,  6850., 25.,  1.5,  2.0,  0.100,  info(11,:))
               call boxinfo( 6470.,  7530., 25.,  2.0,  2.5,  0.100,  info(12,:))
               call boxinfo( 7050.,  8350., 30.,  2.5,  3.0,  0.100,  info(13,:))
               call boxinfo( 7750.,  9350., 30.,  3.0,  3.5,  0.050,  info(14,:))
               call boxinfo( 8500., 11350., 30.,  3.5,  4.0,  0.040,  info(15,:))
               call boxinfo( 9400., 11370., 20.,  4.0,  4.5,  0.010,  info(16,:))
               call boxinfo(10825., 11600., 30.,  4.5,  4.67, 0.025,  info(17,:))

            case default
               write(*,'(/a/)') '  ** Error: Invalid table number'
               stop
         end select

         write(*,*) ''
         write(*,*) ' Creating interpolation table: "', trim(filename), '"'

c        Filename in ASCII codes.
         ncodes = len_trim(filename)
         call toascii(trim(filename), ncodes, codes)

c        Create the table.
         ifile = openfile(ncodes, codes)
         if (ifile .eq. -1) then
            write(*,'(/a/)') '  ** Error: Unable to create the table'
            stop
         end if

c        Write the number of boxes and the scan type.
         n = 1
         if (writefile(ifile, n, real(nboxes)) .ne. n) then
            write(*,'(/a/)') '  ** Error: Unable to write the table'
            stop
         end if
         if (writefile(ifile, n, real(xyscan)) .ne. n) then
            write(*,'(/a/)') '  ** Error: Unable to write the table'
            stop
         end if

c        Foreach box, write the box information.
         do k = 1, nboxes
            if (writefile(ifile, ninfo, info(k,:)) .ne. ninfo) then
               write(*,'(/a/)') '  ** Error: Unable to write the table'
               stop
            end if
         end do

c        Foreach box, write the box values.
         n = 1
         do k = 1, nboxes
            ncols  = int(info(k, 1))
            xini   =     info(k, 2)
            xdelta =     info(k, 3)
            nrows  = int(info(k, 4))
            yini   =     info(k, 5)
            ydelta =     info(k, 6)
            do i = 1, nrows
               y = yini + ydelta*(i-1)
               do j = 1, ncols
                  x = xini + xdelta*(j-1)
                  ier = fpointer(x, y, f)
                  if (writefile(ifile, n, f) .ne. n) then
                     write(*,'(/a/)') '  ** Error: Unable to write the table'
                     stop
                  end if
               end do
            end do
         end do

c        Close the table.
         if (closefile(ifile) .ne. 0) then
            write(*,'(/a/)') '  ** Error: Unable to close the table'
            stop
         end if

      end do  !foreach table

c     Happy ending.
      write(*,*) ''
      stop
      end


c_____________________________________________________________________________
c
c  Routine: boxinfo()
c
c  Save the box (subtable) information.
c  @param  xini initial value in x.
c  @param  xfin final value in x.
c  @param  xdelta offset in x.
c  @param  yini initial value in y.
c  @param  yfin final value in y.
c  @param  ydelta offset in y.
c  @param  vinfo vector with the info values.
c_____________________________________________________________________________

      subroutine boxinfo( xini, xfin, xdelta, yini, yfin, ydelta, vinfo )

      implicit none
      integer*4 ncols, nrows
      real*4    xini, xfin, xdelta, yini, yfin, ydelta, vinfo(6), modx, mody

      modx = mod(xfin-xini, xdelta)
      mody = mod(yfin-yini, ydelta)
      ncols = int((xfin-xini) / xdelta) + 1
      nrows = int((yfin-yini) / ydelta) + 1
      if (modx .gt. 0) ncols = ncols + 1
      if (mody .gt. 0) nrows = nrows + 1

      vinfo(1) = real(ncols)
      vinfo(2) = xini
      vinfo(3) = xdelta
      vinfo(4) = real(nrows)
      vinfo(5) = yini
      vinfo(6) = ydelta

      return
      end


c_____________________________________________________________________________
