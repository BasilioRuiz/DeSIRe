# ****************************************************************************
#
# Fichero   : Makefile
# Proyecto  : DeSIRe
# Dpto.     : Departamento de Astrofisica - ULL
# Autor     : Basilio Ruiz Cobo - brc@iac.es
# Proposito : Mantenimiento de los ejecutables de inversion
#
# Uso       : > make [fc=compilador]             ... genera desire
#             > make [fc=compilador] modelador3  ... genera modelador3
#
#             donde fc puede ser gfortran, ifort o nada (defecto = ifort).
#
# Nota      : Aunque basicamente esta codificado en Fortran 77, hay unas
#             pocas funciones mas avanzadas.
#             get_command_argument() => Fortran 2003
#             iachar()               => Fortran 95
#             trim() y len_trim()    => Fortran 90
#
# 27/01/04  : Primera version.
# 12/11/15  : El compilador a usar puede entrar por linea de comandos.
# 15/05/19  : (epm) Ejecutables de RH convertidos en libreria.
# 29/05/19  : (epm) Incluyo 'make rhf1d' para compilar todo desde aqui.
# 20/06/19  : (epm) Reestructuracion de directorios.
# 04/04/20  : (epm) Ya no distribuimos los ejecutables de RH.
# 05/05/20  : (epm) Borro modulos obsoletos que ya no se usan.
# 06/06/20  : (epm) Anoto como compilar si 'kld' es muy grande.
# 07/07/20  : (epm) Borro el ejecutable 'lande' y modulos no usados.
# 08/08/20  : (epm) Renombro leelineasiii.f por leelineas.f.
#             (epm) Borro write_abund_RH.f y las rutinas de Planck.
# 09/09/20  : (epm) Renombro mreadi<n>.f por mread<n>.f.
# 10/10/20  : (epm) Nuevo modulo lines_wave.f.
#             (epm) Borro write_keyword_input_RH.f.
# 11/11/20  : (epm) Borro read_keyword_input.f y lecturas de atmosferas.
#             (epm) Nuevos modulos para las poblaciones de H:
#                   hpopulations.f, lee_hydro.f, mread4.f.
# 02/02/21  : (epm) Libreria externa para detectar fugas de memoria.
#
# ****************************************************************************


# Recordatorio de la ultima compilacion (creo que solo para Solaris).

.KEEP_STATE:

# Compilador, flags compilador, librerias del sistema y flags de librerias.

# Comprobacion de arrays : -fcheck=bounds (gfortran) / -check bounds (ifort).
# Rastreo error ejecucion: -fbacktrace               / -traceback.
# Extension de lineas    : -ffixed-line-length-none  / -extend_source.
# 24/06/19 epm: Ahora que los ejecutables de RH son rutinas, necesitamos la
# libreria -lpthread para el compilador gfortran.
# 11/11/19 epm: No capturar ninguna excepcion (por la no convergencia) y
# utilizar la opcion que salva las variables locales en memoria estatica
# (si las capturamos y no ponemos -fno-automatic, gfortran produce un SIGFPE).
# 06/06/20 epm: Se puede evitar el error del linker si 'kld' de PARAMETER
# es muy grande (>= 1000). Leer los comentarios del final.

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

# 12/12/20 han: Required to the latest Ubuntu installation (20.10).
ifneq ("$(wildcard /usr/include/tirpc)","")
    ARCHLIBS += -ltirpc
endif

# 02/02/21 epm: Libreria externa para detectar fugas de memoria.
MEMLIB = -L./rh/extern -lcMemDbg

# 15/05/19 epm: Ejecutables de RH convertidos en libreria.
# Las librerias "librhf1d.a" (modulos del directorio rhf1d) y "librh.a"
# (modulos del directorio rh) tienen dependencias circulares. Si queremos
# mantener las dos librerias separadas tenemos que referenciarlas dos veces.
# Ver https://eli.thegreenplace.net/2013/07/09/library-order-in-static-linking
# RHLIB  = ../rhf1d/librhf1d.a -L.. -lrh -lrh_f90
# RHLIB2 = $(RHLIB) $(RHLIB)
# O bien podemos no crear "librhf1d.a" y anadir a "librh.a" todos los modulos
# del directorio rhf1d (opcion elegida).
#
# If no FORTRAN 90 compiler is available remove -lrh_f90 in following.

RHLIB = -L./rh -lrh -lrh_f90

ARFLAGS  = rv

# Fichero de configuracion (definicion de dimensiones).

FILECONF = PARAMETER
FILEMAKE = .make.lib

# Ficheros fuentes de los ejecutables y directorio bin.

MAIN1SRC = desire.f
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
	delo_bezier3_c.f \
	delo_bezier3.f \
	delo_bezier3log.f \
	departures.f \
	deriva.f \
	deriva4cuad.f \
	derivacuad.f \
	dsaha.f \
	dsaha_db.f \
	equi.f \
	equisubmu.f \
	equisubmu_cont.f \
	equisubmu_contmag.f \
	error.f \
	escribeFR.f \
	fint.f \
	fperfil2.f \
	fperfil2err.f \
	gasb.f \
	gasb_db.f \
	gasb_thermo.f \
	gase.f \
	hermite.f \
	hermite_c.f \
	hpopulations.f \
	inicia_pefrompgt.f \
	interpolate1.f \
	intmodel.f \
	intmodel4.f \
	kappach.f \
	kappach_db.f \
	lee_hydro.f \
	leeabun.f \
	leelineas.f \
	leemalla2.f \
	leemallab.f \
	leemodi222.f \
	leepesos.f \
	leeuve2.f \
	leeuveobsindic.f \
	linearsystem.f \
	lines_kurucz.f \
	lines_wave.f \
	locate.f \
	lubksb.f \
	ludcmp.f \
	marqcoef2.f \
	marqcoeferr.f \
	marquarderr.f \
	marquardt2.f \
	matinx.f \
	matinx2.f \
	meves.f \
	molecb.f \
	mread2.f \
	mread3.f \
	mread33.f \
	mread4.f \
	mvoigt.f \
	mvoigtc.f \
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
	polint.f \
	reduce2.f \
	refrax.f \
	rnorma.f \
	saha.f \
	saha_db.f \
	spectra.f \
	splinb.f \
	splines.f \
	splines22.f \
	svbksb.f \
	svdcmp2.f \
	svdmatriz2.f \
	taulinea.f \
	taulinea2.f \
	taulinea4.f \
	taulinea4sub.f\
	thermosub.f \
	voigt.f \
	write_atmos_RH.f \
	zeeman.f \
	zeeman_jk.f

# Referencias para el reemplazamiento de sufijos.

MAIN1OBJ = $(MAIN1SRC:.f=.o)
MAIN3OBJ = $(MAIN3SRC:.f=.o)
MAIN4OBJ = $(MAIN4SRC:.f=.o)
MAIN5OBJ = $(MAIN5SRC:.f=.o)

LIBOBJ = $(LIBSRC:.f=.o)


#
# A) Ejecutable desire ('uso> make' o 'uso> make desire').
#
desire: $(MAIN1OBJ) $(LIB) $(FILEMAKE) rhf1d cmemlib
	$(LINK.f) -o $@ $(MAIN1OBJ) $(LIB) $(RHLIB) $(MEMLIB) $(ARCHLIBS)
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
# H) Reconstruccion de la libreria de deteccion de fugas de memoria.
#
cmemlib:
	$(MAKE) -C rh/extern

#
# I) Todos los ejecutables de DeSIRe.
#
all:  desire  modelador3  PemufromPgT_i  PgmufromPeT_i

#
# J) Instalacion completa de los ejecutables de DeSIRe.
#
install:
#	10/10/19 epm: Limpieza en la distribucion (solve1d no funciona).
	rm -f rh/rhf1d/iter_1d.c rh/rhf1d/solve1d.c rh/rhf1d/solve1d
	rm -f ../bin/solve1d
#	11/11/19 epm: Limpieza en la distribucion (codigo RH nuevo).
	rm -f rh/rhf1d/bezier.c rh/rhf1d/bezier.h
	rm -f rh/rhf1d/piecestokes.c rh/rhf1d/piecewise.c
#	04/04/20 epm: Limpieza en la distribucion (no mas ejecutables de RH).
	rm -f rh/rhf1d/backgrcontr rh/rhf1d/rhf1d rh/rhf1d/solveray
	rm -f ../bin/backgrcontr ../bin/rhf1d ../bin/solveray
#	04/04/20 epm: Directorio no funcional al suprimir las salidas de RH.
	rm -fr rh/idl
#	04/04/20 epm: Fichero renombrado (rh_desire.c -> rh_sir.c).
	rm -f rh/rhf1d/rh_desire.c
#	05/05/20 epm: Suprimo los siguientes modulos porque no son usados.
	rm -f leeuveobs.f mensaje.f model.f nodos2.f util.c
#	07/07/20 epm: Suprimo leelineasii() y el ejecutable lande.
	rm -f leelineasii.f sublande.f sublande2.f lande.f lande
	rm -f ../bin/lande
#	08/08/20 epm: Fichero renombrado (leelineasiii.f -> leelineas.f).
	rm -f leelineasiii.f
#	08/08/20 epm: Suprimo write_abun_RH()y las rutinas de Planck.
	rm -f write_abun_RH.f dtplanck.f planck.f plank2.f
#	09/09/20 epm: Ficheros renombrados (mreadi<n>.f -> mread<n>.f).
	rm -f mreadi2.f mreadi3.f mreadi33.f
#	10/10/20 epm: Suprimo write_keyword_input_RH().
	rm -f write_keyword_input_RH.f
#	11/11/20 epm: Suprimo read_* de 3 modulos que ya no se usan.
	rm -f read_keyword_input_RH.f read_atmos_RH.f read_model_atmos.f
#	11/11/20 epm: Con los modelos por memoria no funciona 'conversion'.
	rm -f rh/rhf1d/conversion ../bin/conversion
#	11/11/20 epm: Por si acaso cqn tiene una version pasada por brc.
	rm -f Hpopulations.f

#	10/10/20 epm: Borro ficheros del directorio run que ya no se usan.
	@find ../run/example ../run/example2 ../run/example_visp \
	       -name "abundance.input" -exec rm {} \; -or \
	       -name "kurucz.input" -exec rm {} \; -or \
	       -name "*.kur" -exec rm {} \; -or \
	       -name "*.wave" -exec rm {} \;

#	Finalmente borramos objetos y recompilamos.
	find . -name "*.a" -exec rm {} \;
	find . -name "*.o" -exec rm {} \;
	$(MAKE) fc=$(FC) all
#	04/04/20 epm: Ya no distribuimos los ejecutables de RH.
#	11/11/20 epm: Con los modelos por memoria no funciona 'conversion'.
#	$(MAKE) -C rh/rhf1d conversion


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
#
#       RELOCATION TRUNCATED TO FIT: ERROR FORTRAN WITH LARGE ARRAYS
#       ------------------------------------------------------------
#
# Setting 'parameter (kld=10000)' in the file "PARAMETER" we get this
# linker error:
# > gfortran -Ofast -ffixed-line-length-none -fno-automatic
#            -ffpe-summary=invalid,zero,overflow 
#            -o desire desire.o libdesire.a -L./rh  -lrh  -lrh_f90 -lpthread
# libdesire.a(fperfil2.o): In function `fperfil2_':
# fperfil2.f:(.text+0xc36): relocation truncated to fit:
#    R_X86_64_PC32 against symbol `numero_lte_' defined in COMMON section in
#    libdesire.a(fperfil2.o)
# fperfil2.f:(.text+0x100d): relocation truncated to fit:
#    R_X86_64_PC32 against symbol `mprimera_' defined in COMMON section in
#    libdesire.a(fperfil2.o)
# fperfil2.f:(.text+0x108f): relocation truncated to fit:
#    R_X86_64_PC32 against symbol `brklm_' defined in COMMON section in
#    libdesire.a(fperfil2.o)
# fperfil2.f:(.text+0x16b0): additional relocation overflows omitted from the
#    output
# collect2: error: ld returned 1 exit status
# Makefile:225: recipe for target 'desire' failed
#
# The error you get is returned by the linker because the size of the
# statically-allocated block exceeds the range of what can be addressed by a
# 32-bit addressing instruction, which is 2 GB. This is irrelevant of whether
# you index your array using 32-bit or 64-bit integers. The problem is related
# to the total size of a statically-allocated array.
# 
# In order to avoid this, you can compile your code with -mcmodel=medium or
# -mcmodel=large. Statically-allocated arrays larger than 2GB are then allowed.
#
# This option tells the compiler to use a specific memory model to generate
# code and store data. It can affect code size and performance. If your
# program has COMMON blocks and local data with a total size smaller than 2GB,
# -mcmodel=small is sufficient. COMMONs larger than 2GB require
# -mcmodel=medium or -mcmodel=large.
#
# The option -mcmodel=<model> is valid for both gfortran and ifort.
#
# (Efectivamente poner el flag -mcmodel=medium empeora un poco la s/n y la
# chi**2 y tarda bastante mas, sobre un 20 por ciento).
#
#
# ============================================================================
