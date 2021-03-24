## ------- file: ------------------------- Makefile ----------------- ##
#
#      Version:       rh2.0, 1-D plane-parallel
#      Author:        Han Uitenbroek (huitenbroek@nso.edu)
#      Last modified: Mon May 21 14:30:59 2018 --
#
#      Updates (epm): Adaptation to the DeSIRe project.
#                     Interoperability C-Fortran.
#                     Conversion of main() as a function for Fortran.
#                     Addition to the "../librh.a" library of all the
#                     objects we need.
#
##     --------------------------                      ----------RH-- ##

# 10/10/19 epm: CPU and OS defined with system calls.
GEOMETRY = ONE_D_PLANE
CPU := $(shell uname -m)
OS  := $(shell uname -s)

include ../makefiles/makefile.$(CPU).$(OS)


## --- If no FORTRAN compiler is available remove -lrh_f90 in following

      LIBS = -lrh -lrh_f90 $(ARCHLIBS)


## --- Define groups of object files --                -------------- ##

ONE_D_OBJS = rhf1d_main.o

ONE_D_OBRC = conversion_main.o

  RAY_OBJS = solveray_main.o

BACKGROPAC = backgrcontr.o

# 27/05/19 epm: Objects to add to the library ("../librh.a").
ADDRH_OBJS = anglequad.o \
             bezier_1D.o \
             conversion.o \
             feautrier.o \
             formal.o \
             hydrostat.o \
             multiatmos.o \
             piecestokes_1D.o \
             piecewise_1D.o \
             project.o \
             rh_chrono.o \
             rh_close.o \
             rh_glob.o \
             rh_sir.o \
             rhf1d.o \
             riiplane.o \
             solveray.o \
             writeflux_xdr.o \
             writegeom_xdr.o

# 20/06/19 epm: Directory for the executables.
BINDIR = ../../../bin


## --- Rules for the executables --                    -------------- ##

# 20/06/19 epm: DeSIRe doesn't need the executables, only the libraries.
libraries:  librh  ../librh.a

all:  rhf1d  conversion  solveray  backgrcontr

# 27/05/19 epm: 3 dependencies.
# 1) compile the executable's sources.
# 2) compile and update the libraries "librh.a" and "librh_f90.a".
# 3) compile and add objects to the library "librh.a".

rhf1d:  $(ONE_D_OBJS)  librh  ../librh.a
	$(LD) -o $@  $(LDFLAGS) $(ONE_D_OBJS) $(LIBS)
	cp -p $@ $(BINDIR)

conversion:  $(ONE_D_OBRC)  librh  ../librh.a
	$(LD) -o $@  $(LDFLAGS) $(ONE_D_OBRC) $(LIBS)
	cp -p $@ $(BINDIR)

solveray:  $(RAY_OBJS)  librh  ../librh.a
	$(LD) -o $@  $(LDFLAGS) $(RAY_OBJS) $(LIBS)
	cp -p $@ $(BINDIR)

backgrcontr:  $(BACKGROPAC)  librh  ../librh.a
	$(LD) -o $@  $(LDFLAGS) $(BACKGROPAC) $(LIBS)
	cp -p $@ $(BINDIR)


## --- If no FORTRAN compiler is available remove librh_f90.a in following.

librh:
	cd ..; $(MAKE) FC=$(FC) librh.a librh_f90.a


## --- 27/05/19 epm: Addition to the library of all the objects we need.

../librh.a: $(ADDRH_OBJS)
	ar $(ARFLAGS) $@ $?


## --- Clean up --                                     -------------- ##

clean:
	rm -f *.o  rhf1d  conversion  solveray  backgrcontr


## --- Explicit dependencies on include files --       -------------- ##


anglequad.o:            ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h

bezier_1D.o:            ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../bezier.h

backgrcontr.o:          ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../background.h  \
                        ../inputs.h  ../error.h

conversion.o:           ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../statistics.h  ../inputs.h

feautrier.o:            ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../error.h

formal.o:               ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../constant.h  \
                        ../background.h  ../inputs.h  ../error.h  ../xdr.h

hydrostat.o:            ../rh.h  ../atom.h  ../atmos.h  geometry.h \
                        ../background.h  ../constant.h  ../accelerate.h \
                        ../error.h

multiatmos.o:           ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../constant.h  ../error.h  ../inputs.h  \
                        ../statistics.h  ../xdr.h  ../desire.h

piecestokes_1D.o:       ../rh.h  ../error.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h

piecewise_1D.o:         ../rh.h  ../error.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h

project.o:              ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h

rh_close.o:             ../rh.h  ../atom.h  ../atmos.h  ../spectrum.h  \
                        ../statistics.h  ../inputs.h

rh_glob.o:              ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  \
                        ../statistics.h  ../inputs.h  ../desire.h

rh_sir.o:               ../desire.h

rhf1d.o:                ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../statistics.h  ../inputs.h  ../xdr.h

riiplane.o:             ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../constant.h  ../error.h

solveray.o:             ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../statistics.h  ../inputs.h  ../error.h  ../xdr.h

writeflux_xdr.o:        ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../constant.h  \
                        ../error.h  ../xdr.h

writegeom_xdr.o:        ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../error.h  ../inputs.h  ../xdr.h

## ------- end ---------------------------- Makefile ---------------- ##
