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
# Extension de lineas de codigo: -ffixed-line-length-none o -extend_source.
# 24/06/19 epm: Ahora que los ejecutables de RH son rutinas, necesitamos la
# libreria -lpthread para el compilador gfortran.

ifeq ($(fc),gfortran)
FC       = gfortran
FFLAGS   = -Ofast -ffixed-line-length-none -fno-automatic
ARCHLIBS = -lpthread
else ifeq ($(fc),ifort)
FC       = ifort
FFLAGS   = -O2 -extend_source
ARCHLIBS =
else
FC       = ifort
FFLAGS   = -O2 -extend_source
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
	rm -f rh/rhf1d/iter_1d.c rh/rhf1d/iter_1d.o
	rm -f rh/rhf1d/solve1d.c rh/rhf1d/solve1d.o rh/rhf1d/solve1d
	rm -f ../bin/solve1d
	find . -name "*.a" -exec rm {} \;
	find . -name "*.o" -exec rm {} \;
	$(MAKE) fc=$(FC) all
	$(MAKE) -C rh/rhf1d all


# ****************************************************************************