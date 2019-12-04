# ****************************************************************************
#
# Fichero   : Makefile
# Proyecto  : Metodos de inversion
# Dpto.     : Departamento de Astrofisica - ULL
# Autor     : Basilio Ruiz Cobo - brc@iac.es
# Proposito : Mantenimiento de los ejecutables de inversion
#
# Uso       : > make [fc=compilador]             ... genera desire
#             > make [fc=compilador] lande       ... genera lande
#             > make [fc=compilador] modelador3  ... genera modelador3
#
#             donde fc puede ser gfortran, ifort o nada (defecto = ifort).
#             Por ejemplo:
#
#             > make desire
#             > make fc=gfortran lande
#             > make fc=ifort modelador3
#
# 27/01/04  : primera version.
# 12/11/15  : el compilador a usar puede entrar por linea de comandos.
# 15/05/19  : ejecutables de RH convertidos en libreria.
# 20/06/19  : reestructuracion de directorios.
#
# ****************************************************************************


# Recordatorio de la ultima compilacion (creo que solo para Solaris).

.KEEP_STATE:

# Compilador, flags compilador, librerias del sistema y flags de librerias.
#
# Comprobacion de arrays : -fcheck=bounds (gfortran) / -check bounds (ifort).
# Rastreo error ejecucion: -fbacktrace               / -traceback.
# Extension de lineas    : -ffixed-line-length-none  / -extend_source.
# 24/06/19 epm: Ahora que los ejecutables de RH son rutinas, necesitamos la
# libreria -lpthread para el compilador gfortran.
# 11/11/19 brc: No capturar ninguna excepcion (por la no convergencia) y
# utilizar la opcion que salva las variables locales en memoria estatica
# (si las capturamos y no ponemos -fno-automatic, gfortran produce un SIGFPE).

ifeq ($(fc),gfortran)
FC       = gfortran
FFLAGS   = -Ofast -ffixed-line-length-none -fno-automatic \
           -ffpe-summary=invalid,zero,overflow
ARCHLIBS = -lpthread
else
FC       = ifort
FFLAGS   = -O2 -extend_source -noautomatic
ARCHLIBS =
endif

CC       = gcc
CFLAGS   = -O2

ARFLAGS  = ruv

# Fichero de configuracion (definicion de dimensiones).

FILECONF = PARAMETER
FILEMAKE = .make.lib

# 20/06/17 brc: Fichero en C que sera llamado desde Fortran (para XDR).
# 20/06/17 brc: Libreria C+Fortran que sera llamada desde Fortran (para XDR).
#
# EXTERNC   = xdrlib_.o
# EXTERNLIB = libfxdr.a

# 15/05/19 epm: Ejecutables de RH convertidos en libreria.
# Las librerias "librhf1d.a" (modulos del directorio rhf1d) y "librh.a"
# (modulos del directorio rh) tienen dependencias circulares. Si queremos
# mantener las dos librerias separadas tenemos que referenciarlas dos veces.
# Ver https://eli.thegreenplace.net/2013/07/09/library-order-in-static-linking
# RHLIB     = ../rhf1d/librhf1d.a -L.. -lrh -lrh_f90
# EXTERNLIB = $(RHLIB) $(RHLIB)
# O bien podemos no crear "librhf1d.a" y anadir a "librh.a" todos los modulos
# del directorio rhf1d (opcion elegida).
#
# If no FORTRAN 90 compiler is available remove -lrh_f90 in following.

EXTERNLIB = -L./rh  -lrh  -lrh_f90

# Ficheros fuentes de los ejecutables y directorio bin.

MAIN1SRC = desire.f
MAIN2SRC = lande.f
MAIN3SRC = modelador3.f
MAIN4SRC = PemufromPgT_i.f
MAIN5SRC = PgmufromPeT_i.f

BINDIR = ../bin

# Libreria de inversion y todos sus modulos.

LIB = libdesire.a

LIBSRC =amp2.f \
	atmdatb.f \
	automatico.f \
	automaticop.f \
	blends2.f \
	blendscon2.f \
	cabecera.f \
	chachi.f \
	comprime2.f \
	conhsra.f \
	deconv.f \
	deconv2.f \
	delta1.f \
	departures.f \
	deriva.f \
	deriva4cuad.f \
	derivacuad.f \
	dsaha_db.f \
	dsaha.f \
	dtplanck.f \
	equi.f \
	equisubmu_cont.f \
	equisubmu_contmag.f \
	equisubmu.f \
	escribeFR.f \
	fint.f \
	fperfil2.f \
	fperfil2err.f \
	gasb_db.f \
	gasb_thermo.f \
	gasb.f \
	gase.f \
	hermite_c.f \
	hermite.f \
	inicia_pefrompgt.f \
	interpolate1.f \
	intmodel.f \
	intmodel4.f \
	kappach_db.f \
	kappach.f \
	leeabun.f \
	leelineasii.f \
	leelineasiii.f \
	leemalla2.f \
	leemallab.f \
	leemodi222.f \
	leepesos.f \
	leeuve2.f \
	leeuveobs.f \
	leeuveobsindic.f \
	lines_kurucz.f \
	locate.f \
	lubksb.f \
	ludcmp.f \
	marqcoef2.f \
	marqcoeferr.f \
	marquarderr.f \
	marquardt2.f \
	matinx.f \
	matinx2.f \
	mensaje.f \
	meves.f \
	model.f \
	molecb.f \
	mreadi2.f \
	mreadi3.f \
	mreadi33.f \
	mvoigt.f \
	mvoigtc.f \
	nodos2.f \
	nodos2aut.f \
	numericalsub.f \
	PemufromPgT_isub.f \
	pemufrompgtaut.f \
	penalty.f \
	penalty2.f \
	petaufrompgzro.f \
	PgmufromPeT_isub.f \
	pgpefromrho.f \
	pgzrofrompetau.f \
	planck.f \
	plank2.f \
	polint.f \
	read_atmos_RH.f \
	read_keyword_input_RH.f \
	read_model_atmos.f \
	reduce2.f \
	refrax.f \
	rnorma.f \
	saha_db.f \
	saha.f \
	splinb.f \
	splines.f \
	splines22.f \
	sublande.f \
	sublande2.f \
	svbksb.f \
	svdcmp2.f \
	svdmatriz2.f \
	taulinea.f \
	taulinea2.f \
	taulinea4.f \
	taulinea4sub.f\
	thermosub.f \
	voigt.f \
	write_abun_RH.f \
	write_atmos_RH.f \
	write_keyword_input_RH.f \
	zeeman_jk.f \
	zeeman.f

# Referencias para el reemplazamiento de sufijos.

MAIN1OBJ = $(MAIN1SRC:.f=.o)
MAIN2OBJ = $(MAIN2SRC:.f=.o)
MAIN3OBJ = $(MAIN3SRC:.f=.o)
MAIN4OBJ = $(MAIN4SRC:.f=.o)
MAIN5OBJ = $(MAIN5SRC:.f=.o)

LIBOBJ = $(LIBSRC:.f=.o)


#
# A) Ejecutable desire ('uso> make' o 'uso> make desire').
#
desire: $(MAIN1OBJ) $(LIB) $(FILEMAKE) rhf1d
	$(LINK.f) -o $@ $(MAIN1OBJ) $(LIB) $(EXTERNLIB) $(ARCHLIBS)
	cp -p $@ $(BINDIR)

#
# B) Ejecutable lande ('uso> make lande').
#
lande: $(MAIN2OBJ) $(LIB)
	$(LINK.f) -o $@ $(MAIN2OBJ) $(LIB)
	cp -p $@ $(BINDIR)

#
# C) Ejecutable modelador3 ('uso> make modelador3').
#
modelador3: $(MAIN3OBJ) $(LIB)
	$(LINK.f) -o $@ $(MAIN3OBJ) $(LIB)
	cp -p $@ $(BINDIR)

#
# D) Ejecutable PemufromPgT_i ('uso> make PemufromPgT_i').
#
PemufromPgT_i: $(MAIN4OBJ) $(LIB)
	$(LINK.f) -o $@ $(MAIN4OBJ) $(LIB)
	cp -p $@ $(BINDIR)

PgmufromPeT_i: $(MAIN5OBJ) $(LIB)
	$(LINK.f) -o $@ $(MAIN5OBJ) $(LIB)
	cp -p $@ $(BINDIR)

#
# E) Actualizacion de la libreria de DeSIRe.
#
$(LIB): $(LIBOBJ)
	$(AR) $(ARFLAGS) $@ $?

#
# F) Reconstruccion de la libreria de DeSIRe.
#
$(FILEMAKE): $(FILECONF)
	$(COMPILE.f) $(LIBSRC)
	$(AR) $(ARFLAGS) $(LIB) $(LIBOBJ)
	@touch $@

#
# G) Reconstruccion de las librerias de RH.
#
rhf1d:
	$(MAKE) -C rh/rhf1d

#
# H) Todos los ejecutables de DeSIRe.
#
all:  desire  lande  modelador3  PemufromPgT_i  PgmufromPeT_i

#
# I) Instalacion completa de los ejecutables de DeSIRe y de RH.
#
install:
#	Limpieza en la distribucion del 10/10/19 (solve1d no funciona).
	rm -f rh/rhf1d/iter_1d.c rh/rhf1d/solve1d.c rh/rhf1d/solve1d
	rm -f ../bin/solve1d
#	Limpieza en la distribucion del 11/11/19 (codigo RH nuevo).
	rm -f rh/rhf1d/bezier.c rh/rhf1d/bezier.h
	rm -f rh/rhf1d/piecestokes.c rh/rhf1d/piecewise.c
	find . -name "*.a" -exec rm {} \;
	find . -name "*.o" -exec rm {} \;
	$(MAKE) fc=$(FC) all
	$(MAKE) -C rh/rhf1d all


# ****************************************************************************
# ============================================================================
#
#                            GFORTRAN OPTIONS
#                            ----------------
#
# -ffpe-trap=list
#     Specify a list of floating point exception traps to enable. On
#     most systems, if a floating point exception occurs and the trap for
#     that exception is enabled, a SIGFPE signal will be sent and the
#     program being aborted, producing a core file useful for debugging.
#     "list" is a (possibly empty) comma-separated list of the following
#     exceptions: invalid (invalid floating point operation, such as
#     "SQRT(-1.0)"), zero (division by zero), overflow (overflow in a
#     floating point operation), underflow (underflow in a floating point
#     operation), inexact (loss of precision during operation), and
#     denormal (operation performed on a denormal value). The first five
#     exceptions correspond to the five IEEE 754 exceptions, whereas the
#     last one (denormal) is not part of the IEEE 754 standard but is
#     available on some common architectures such as x86.
#     (Denormal numbers are numbers where this representation would result
#      in an exponent that is below the smallest representable exponent.)
#
#     The first three exceptions (invalid, zero, and overflow) often
#     indicate serious errors, and unless the program has provisions for
#     dealing with these exceptions, enabling traps for these three
#     exceptions is probably a good idea.
#
#     Many, if not most, floating point operations incur loss of
#     precision due to rounding, and hence the "ffpe-trap=inexact" is
#     likely to be uninteresting in practice.
#
#     By default no exception traps are enabled.
#
#
# -ffpe-summary=list
#     Specify a list of floating-point exceptions, whose flag status is
#     printed to "ERROR_UNIT" when invoking "STOP" and "ERROR STOP".
#     "list" can be either none, all or a comma-separated list of the
#     following exceptions: invalid, zero, overflow, underflow, inexact
#     and denormal. (See -ffpe-trap for a description of the exceptions.)
#
#     By default, a summary for all exceptions but inexact is shown.
#
#
# -fno-automatic
#     Treat each program unit (except those marked as RECURSIVE) as if the
#     SAVE statement were specified for every local variable and array
#     referenced in it. Does not affect common blocks. The default, which is
#     "-fautomatic", uses the stack for local variables smaller than the value
#     given by "-fmax-stack-var-size".  Use the option "-frecursive" to use no
#     static memory.
#
#
# ============================================================================
#
#                            IFORT OPTIONS
#                            -------------
#
# -fpe<n>
#     Specifies floating-point exception handling for the main program at
#     run-time. You can specify one of the following values for <n>:
#
#     0 - Floating-point invalid, divide-by-zero, and overflow exceptions
#     are enabled. If any such exceptions occur, execution is aborted.
#     Underflow results will be set to zero unless you explicitly specify
#     -no-ftz. On Itanium(R)-based systems, underflow behavior is equivalent
#     to specifying option -ftz. On IA-32 or Intel(R) EM64T systems,
#     underflow results from SSE instructions, as well as x87 instructions,
#     will be set to zero. By contrast, option -ftz only sets SSE underflow
#     results to zero.
#
#     To get more detailed location information about where the error
#     occurred, use -traceback.
#
#     1 - All floating-point exceptions are disabled. Underflow results will
#     be set to zero unless you explicitly specify -no-ftz.
#
#     3 - All floating-point exceptions are disabled. Floating-point
#     underflow is gradual, unless you explicitly specify a compiler option
#     that enables flush-to-zero. This is the default; it provides full IEEE
#     support.
#
# -ftz
#     Flushes denormal results to zero when the application is in the gradual
#     underflow mode. It may improve performance if the denormal values aren't
#     critical to the behavior of your application. The default is -no-ftz.
#
#     The following options set the -ftz option: -fpe0, -fpe1, and on
#     Itanium(R)-based systems option -O3. Option -O2 sets the -no-ftz option.
#
#     Note: When SSE instructions are used on IA-32 systems, option -no-ftz is
#     ignored. However, you can enable gradual underflow by calling a function
#     in C in the main program that clears the FTZ and DAZ bits in the MXCSR
#     or by calling the function for_set_fpe in the main program to clear
#     those bits. Be aware that denormal processing can significantly slow
#     down computation.
#
# -noautomatic
#     Places variables, except those declared as AUTOMATIC, in static memory.
#     The default is -auto-scalar. However, if you specify -recursive or
#     -openmp, the default is -automatic.
#
#
# ============================================================================
