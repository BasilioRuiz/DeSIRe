# ****************************************************************************
#
# File          : Makefile
# Description   : Compilations and links.
# Project       : DeSIRe
# Creation date : 12/12/22
# Author        : epm@iac.es
#
# ****************************************************************************


ifeq ($(fc),gfortran)
   FFLAGS = -ffixed-line-length-none
   FC     = gfortran
else
   FFLAGS = -extend_source
   FC     = ifort
endif

CC = gcc

LIB = libdesire.a

CTABLE_OBJS = io_fns.o interface_fns.o ctable.o
ETABLE_OBJS = io_fns.o interface_fns.o bilineal_fns.o etable.o
PTABLE_OBJS = ptable.o

#TABLES = $(filter table.%,$(shell ls tables/*))
TABLES = $(filter table.%,$(shell ls *))

tables:  ctable  etable  ptable

all:  libdesire  ctable  etable  ptable

libdesire:
	@$(MAKE) --no-print-directory -C .. fc=$(FC) $(LIB)

ctable: $(CTABLE_OBJS) ../$(LIB)
	$(FC) -o $@ $(CTABLE_OBJS) ../$(LIB)

etable: $(ETABLE_OBJS) ../$(LIB)
	$(FC) -o $@ $(ETABLE_OBJS) ../$(LIB)

ptable: $(PTABLE_OBJS)
	$(CC) -o $@ $(PTABLE_OBJS) -lm

pepe:
ifeq ($(strip $(TABLES)),)
	@echo
	@echo "  ---------------------------------"
	@echo "   INTERPOLATED TABLES NOT CREATED"
	@echo "  ---------------------------------"
	@echo "   Please, go to directory tables ??? and run ctable."
	@echo "   > cd tables"
	@echo "   > ctable"
	@echo
else
	@echo -n
endif


# ****************************************************************************
