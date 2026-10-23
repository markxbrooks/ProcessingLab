#-----------------------------------------------------------------------------
# Makefile for DPS to MOSFLM, Harry Powell
#
# Version 1.0; Builds DPSLIB.a under OSF on an Alpha; if the DPS object
#              libraries don't exist, we should make them via the "make all"
#              command. Bear in mind that we are using a local version of 
#              BRAVAIS.C (bravais using four arguments rather than one), so we 
#              have to include an extra make for a local version of that 
#              library here as well.
#
#-----------------------------------------------------------------------------
DPS        = ..
BIN        = ${DPS}/bin
AR_FLAGS   = ruv
TIDY       = rm -f sta*
VERBOSE    = @VERBOSE@
CC         = gcc -c
# Pre-C99 DPS sources call undeclared functions; keep that a warning.
CFLAGS     = -g -O0  -Wno-error=implicit-function-declaration
RANLIB     = ranlib

# 	location of the modified DPS libraries and source code
LIB        = ${DPS}/lib
SRC        = ${DPS}/src

# 	object files from local interfaces to DPS routines
HOBJS      = write_ub_matrix.o bravais_fortran.o \
              new_matrix.o

# 	object files from DPS bravais lattice generation/reduction routines
IOBJS      = bravais.o  bravais_goof.o \
        bravais_permutate.o bravais_transform.o

# these are the Cambridge _modified_ library routines from DPS
INDEX      = ${DPS}/index
DPSINDEX   = ${SRC}/dps/index
DPSBRAVAIS = ${SRC}/dps/bravais
DPSPEAK    = ${SRC}/dps/peak_search
DPSUTIL    = ${SRC}/dps/util

# libraries ready for linking with MOSFLM objects
LIBRARIES  = ${LIB}/libdps.a ${LIB}/libdps_index.a \
             ${LIB}/libpeak.a ${LIB}/libdps_bravais.a

DEBUG      = 

SOURCES = bravais.c bravais.h bravais.fortran.c bravais_goof.h \
	bravais_help.h bravais_permutate.h bravais_transform.h \
	new_matrix.c write_ub_matrix.c

DIST-OTHER = Makefile.in

all:	    ${LIB}/DPSLIB.a

${LIB}/DPSLIB.a: ${LIB}/dps_index.a ${LIBRARIES}
	cd ${LIB} ; ${TIDY} ; \
	ar x dps_index.a ; ${TIDY} ; \
	ar x libdps.a ; ${TIDY} ; \
	ar x libfftpack.a ; ${TIDY} ; \
	ar x libdps_index.a ; ${TIDY} ; \
	ar ${AR_FLAGS} DPSLIB.a `ls *.o` ; ${TIDY} ; \
	rm -f *.o ; \
	${RANLIB} DPSLIB.a

check:
# This only checks the syntax of this Makefile!
	cat -v -t -e Makefile|more

# Here are the makes for DPS libraries
${LIB}/libdps.a: 
	cd ${DPSUTIL}; \
	${MAKE} -f Makefile; 

${LIB}/libdps_index.a: 
	cd ${DPSINDEX} ; \
	${MAKE} -f Makefile ; 

# local version of libdps_index.a
${LIB}/dps_index.a: ${IOBJS} ${HOBJS}
	mkdir -p ${LIB}
	ar ru ${LIB}/dps_index.a ${IOBJS} ${HOBJS}

${LIB}/libpeak.a: 

${LIB}/libdps_bravais.a: 

bravais_fortran.o: bravais_fortran.c bravais.o
	${CC} ${DEBUG} ${CFLAGS} -o bravais_fortran.o \
	bravais_fortran.c

bravais.o: bravais.c bravais.h
	${CC} ${DEBUG} ${CFLAGS}  -o bravais.o bravais.c

bravais_goof.o: ${DPSBRAVAIS}/bravais_goof.c
	${CC} ${DEBUG} ${CFLAGS}  -o bravais_goof.o \
	${DPSBRAVAIS}/bravais_goof.c 

bravais_permutate.o: ${DPSBRAVAIS}/bravais_permutate.c
	${CC} ${DEBUG} ${CFLAGS}  -o bravais_permutate.o \
	${DPSBRAVAIS}/bravais_permutate.c 

bravais_transform.o: 
	${CC} ${DEBUG} ${CFLAGS}  -o bravais_transform.o \
	${DPSBRAVAIS}/bravais_transform.c ; 

new_matrix.o: new_matrix.c
	${CC} ${DEBUG} ${CFLAGS}  new_matrix.c

write_ub_matrix.o: write_ub_matrix.c
	${CC} ${DEBUG} ${CFLAGS}  write_ub_matrix.c

mostlyclean-compile:
	-rm -f *.o core *.core
	cd ${DPSBRAVAIS} ; rm -f *.o core *.core
	cd ${LIB} ; rm -f *.o core *.core

clean-generic:
	-rm -f *~
	cd ${DPSBRAVAIS} ; rm -f *~

clean-lib: 
	-test -z ${LIB}/DPSLIB.a || rm -f ${LIB}/DPSLIB.a
 
clean: clean-generic mostlyclean-compile
	-test -z ${LIB}/dps_index.a || rm -f ${LIB}/dps_index.a
	cd ${DPSINDEX}; ${MAKE} clean 
	cd ${DPSUTIL} ; ${MAKE} clean

distclean-generic:
	-rm -f Makefile

distclean: distclean-generic clean clean-lib
	cd ${DPSINDEX} ; ${MAKE} distclean
	cd ${DPSUTIL} ; ${MAKE} distclean

realclean: distclean-generic clean clean-lib
	rm -f ${SOURCES} ${DIST-OTHERS}
	cd ${DPSINDEX} ; ${MAKE} realclean
	cd ${DPSBRAVAIS} ; rm -f ${BRAVAIS-FILES}
	cd ${DPSUTIL} ; ${MAKE} realclean
	cd ${DPSPEAK} ; rm -f ${PEAKSEARCH-FILES}

BRAVAIS-FILES = bravais_help.h dps_bravais_read.c read_om.h \
	Makefile.in bravais_permutate.c dps_bravais_read.h  read_peaks.h \
	allocate.h bravais_permutate.h fftn.h ref_cell_on_xyz.h \
	bravais.c bravais_transform.c ind_stat.h util.h \
	bravais.h bravais_transform.h ind_transform.h util_def.h \
	bravais_goof.c call_blas.h mat_inversion.h vector_math.h \
	bravais_goof.h cryst_matrix.h nrutil.h \
	bravais_help.c dps_bravais.c read_frame.h

PEAKSEARCH-FILES = dps_peaksearch.c

uninstall: clean-lib

