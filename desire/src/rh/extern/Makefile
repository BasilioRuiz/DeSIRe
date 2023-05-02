# ****************************************************************************
#
# File          : Makefile
# Description   : Compilation of the memory leak detection library.
# Project       : DeSIRe
# Creation date : 02/02/21
# Author        : epm@iac.es
#
# ****************************************************************************


# To use valgrind we need to compile with flag -g.
#
# > valgrind --leak-check=full [--show-leak-kinds=all] demo

CC     = gcc
CFLAGS =

libcMemDbg.a: cMemDbg.c cMemDbg.h
	$(CC) -c $(CFLAGS) cMemDbg.c
	$(AR) $(ARFLAGS) $@ cMemDbg.o

demo: demo.c libcMemDbg.a
	$(CC) -c $(CFLAGS) demo.c
	$(LINK.c) -o $@  $(LDFLAGS) demo.o libcMemDbg.a

all:  libcMemDbg.a  demo


# ****************************************************************************
