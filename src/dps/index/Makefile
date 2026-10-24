# Makefile for DPS, Ingo Steller
#
# make test: for making al routines and the test program
# make:		 make all
#
#
# Version 0.1
STDFLAGS = -g -O0 
UTIL_INC = ../../../src/dps/utils
LIBDIR = ../../../lib
BIN = ../../../bin
F77 = f77 -c
FFLAGS = -g -O2 -fallow-argument-mismatch
CC = gcc -c
# Pre-C99 DPS sources call undeclared functions; keep that a warning.
CFLAGS = ${STDFLAGS} -I${UTIL_INC} -Wno-error=implicit-function-declaration
AR_FLAGS = ruv

IOBJS = index.o error_exit.o \
	ind_gen_direct.o ind_sort_merge.o ind_choose.o \
	ind_help.o ind_fortran.o
DOBJS = dps_index.o dps_index_read.o
FOBJS = dps_fortran_test.o
LIBS = -ljpeg -lm 
LOCALLIBS = libdps_index.a libdps.a libfftpack.a

DIST-COMMON = Makefile.in

HEADERS     = allocate.h call_blas.h cryst_matrix.h dps_index_read.h \
	fftn.h ind_choose.h ind_gen_direct.h ind_help.h ind_sort_merge.h \
	ind_stat.h ind_transform.h index.h linux_defs.h mat_inversion.h \
	nrutil.h read_frame.h read_om.h read_peaks.h ref_cell_on_xyz.h \
	util.h util_def.h vector_math.h index_def.h

SOURCES     = dps_fortran_test.f  dps_index.c dps_index_read.c error_exit.f \
	ind_choose.c ind_fortran.c ind_gen_direct.c ind_help.c \
	ind_sort_merge.c index.c 

all: index 

dps_fortran: index $(FOBJS)
	$(F77) $(FFLAGS) -o $(BIN)/dps_fortran_test $(FOBJS) -L$(LIBDIR) \
	$(LOCALLIBS) $(LIBS) 

index: $(IOBJS) 
	${AR} $(AR_FLAGS) $(LIBDIR)/libdps_index.a $(IOBJS)
	chmod 644 $(LIBDIR)/libdps_index.a

dps_index: $(LIBDIR)/libdps.a index $(DOBJS)
	${CC} $(CFLAGS) -o $(BIN)/dps_index $(DOBJS) $(LIB)

dps_index_read.o: dps_index_read.c dps_index_read.h

dps_index.o: dps_index.c index.h

index.o: index.c index.h

invert.o: invert.f global.inc include.inc
	${F77} ${FFLAGS} -o invert.o invert.f

error_exit.o: error_exit.f
	${F77} ${FFLAGS} -o error_exit.o error_exit.f

ind_gen_direct.o: ind_gen_direct.c ind_gen_direct.h

ind_sort_merge.o: ind_sort_merge.c ind_sort_merge.h index.h

ind_choose.o: ind_choose.c index.h ind_choose.h

ind_help.o:

ind_fortran.o: ind_fortran.c index.h

dps_fortran_test.o: 

mostlyclean-compile:
	-rm -f *.o core *.core

clean-generic:;  -rm -f *~

distclean-generic:
	-rm -f Makefile

clean-bin:
	-rm -f $(LIBDIR)/libdps_index.a $(BIN)/dps_fortran_test $(BIN)/dps_index

clean: clean-generic mostlyclean-compile

distclean: distclean-generic clean-bin clean

realclean: distclean
	-rm -f ${SOURCES} ${HEADERS} ${DIST-COMMON}

uninstall: clean-bin
