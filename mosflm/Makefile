# Makefile tail for building ipmosflm on any UNIX box...
# This must be "inc"'d by a header Makefile
DPS    = ..
CC     = gcc -c
F77    = gfortran -c
FLINK  = gfortran
FFLAGS = -g -O2 -fno-second-underscore -fallow-argument-mismatch
MOSFLAGS =  -g -O2 -fno-second-underscore -fallow-argument-mismatch -cpp
CFLAGS = -O0 -g -O0 -Wno-implicit-function-declaration -Wno-error=implicit-function-declaration -Wno-return-type -I../jpg -I/opt/homebrew/opt/hdf5/include -I../lz4
LIBS = -L/Applications/ccp4-9/lib -lccp4f -lccp4c -lmmdb2 -L/usr/X11/lib -lSM -lICE \
   -lgfortran -lm -lncurses -L/opt/homebrew/opt/hdf5/lib -lhdf5_hl -lhdf5
LDFLAGS = ../jpg/libjpeg.a
ipmosflm_LDFLAGS = 
DPSBIN = ${DPS}/bin
DPSLIB = ${DPS}/lib
M4     = @M4@

SOURCES =  control.f mosflm_all_ip_inc.f unoptimizable.f pck.c \
	ctimer.c datetime.c chelp.c open_socket.c close_socket.c read_socket.c \
	write_socket.c \
	mosflm_jpeg.c fileaccess.c matzerr.c hklmanage.c torefix.f \
	pilatus_c.c pilatus_pck.c cbfwrap.c corecounter.c byteoffset.c \
	mosflm_all_ip_inc_noX.f torefix_noX.f

DIST-OTHERS = add_refix.com

#OBJS   = control.o unoptimizable.o pck.o ctimer.o datetime.o chelp.o open_socket.o \
	#close_socket.o read_socket.o write_socket.o mosflm_jpeg.o fileaccess.o matzerr.o \
	#hklmanage.o pilatus_c.o pilatus_pck.o cbfwrap.o corecounter.o byteoffset.o \
	#fromhdf5.o lz4.o h5zlz4.o bshuf_h5filter.o bshuf_h5plugin.o bitshuffle.o \
	#mosflm_all_ip_inc.o torefix.o 
OBJS   = control.o unoptimizable.o pck.o ctimer.o datetime.o open_socket.o \
	close_socket.o read_socket.o write_socket.o mosflm_jpeg.o fileaccess.o matzerr.o \
	hklmanage.o pilatus_c.o pilatus_pck.o cbfwrap.o corecounter.o byteoffset.o \
	fromhdf5.o lz4.o h5zlz4.o bshuf_h5filter.o bshuf_h5plugin.o bitshuffle.o \
	mosflm_all_ip_inc_noX.o torefix_noX.o

LOCALLIBS   = ${DPSLIB}/libdps_index.a ${DPSLIB}/DPSLIB.a ${DPSLIB}/libfftpack.a ${DPS}/cbf/lib/libcbf.a

.SUFFIXES:	.o .c .f

.c.o:
	${CC} ${CFLAGS} -o $@ $<

.f.o:
	${F77} ${FFLAGS} -o $@ $<

# Phase D / D6: Fortran amalgamation is opt-in.
#   make MOSFLM_FORTRAN_ORACLE=1 ipmosflm
# Optional Phase C C++ PREDICT splice (also requires oracle flag):
#   cmake -S ../mosflm-cpp -B ../mosflm-cpp/build -DMOSFLM_FORTRAN_ORACLE=ON
#   cmake --build ../mosflm-cpp/build --target mosflm_predict_hook
#   make MOSFLM_FORTRAN_ORACLE=1 MOSFLM_CPP_PREDICT=1 ipmosflm
# Default ProcessingLab binary (no amalgamation):
#   cmake --build ../mosflm-cpp/build --target ProcessingLab
MOSFLM_CPP_DIR = ../mosflm-cpp
MOSFLM_CPP_BUILD = $(MOSFLM_CPP_DIR)/build
# Link flags (must not pass gfortran-only FFLAGS to g++).
LINKFLAGS =
ifeq ($(MOSFLM_CPP_PREDICT),1)
MOSFLAGS += -DMOSFLM_CPP_PREDICT
FLINK = g++
# Homebrew gfortran is not on the default clang++ library search path.
GFORTRAN_LIBDIR := $(shell dirname `gfortran -print-file-name=libgfortran.dylib`)
LIBS += -L$(MOSFLM_CPP_BUILD) -lmosflm_predict_hook -lmosflm_geometry \
	-L$(GFORTRAN_LIBDIR) -lgfortran -lc++
else
LINKFLAGS = ${FFLAGS}
endif

ipmosflm:
ifeq ($(MOSFLM_FORTRAN_ORACLE),1)
	${MAKE} ipmosflm-oracle
else
	@echo "Phase D / D6: Fortran amalgamation is opt-in."
	@echo "  make MOSFLM_FORTRAN_ORACLE=1 ipmosflm"
	@echo "Default C++ CONTROL binary (ProcessingLab):"
	@echo "  cmake --build $(MOSFLM_CPP_BUILD) --target ProcessingLab"
	@echo "  # or: make -C .. cpp"
	@false
endif

ipmosflm-oracle:	${OBJS}
		mkdir -p ${DPSBIN}
		${FLINK} ${ipmosflm_LDFLAGS} ${LINKFLAGS} ${OBJS} ${LOCALLIBS} ${LIBS} \
                -o ${DPSBIN}/ipmosflm ${LDFLAGS}

mosflm_all_ip_inc.o: mosflm_all_ip_inc.f
	${F77} ${MOSFLAGS} mosflm_all_ip_inc.f

#mosflm_all_ip_inc.o: mosflm_all_ip_inc.f
	#/bin/mv -f mosflm_all_ip_inc.f mosflm_all_ip_inc.f.bak
	#sed s/"integer\*1"/"BYTE"/ mosflm_all_ip_inc.f.bak > mosflm_all_ip_inc.f
	#${F77} ${MOSFLAGS} mosflm_all_ip_inc.f
	#/bin/mv -f mosflm_all_ip_inc.f.bak mosflm_all_ip_inc.f

mosflm_all_ip_inc_noX.o: mosflm_all_ip_inc_noX.f
	${F77} ${MOSFLAGS} mosflm_all_ip_inc_noX.f

#mosflm_all_ip_inc_noX.o: mosflm_all_ip_inc_noX.f
	#/bin/mv -f mosflm_all_ip_inc_noX.f mosflm_all_ip_inc_noX.f.bak
	#sed s/"integer\*1"/"BYTE"/ mosflm_all_ip_inc_noX.f.bak > mosflm_all_ip_inc_noX.f
	#${F77} ${MOSFLAGS} mosflm_all_ip_inc_noX.f
	#/bin/mv -f mosflm_all_ip_inc_noX.f.bak mosflm_all_ip_inc_noX.f


control.o:  control.f
		${F77} ${FFLAGS} -O0 control.f

unoptimizable.o: unoptimizable.f
		${F77} ${FFLAGS} -O0 unoptimizable.f

pck.o: pck.c
		${CC} ${CFLAGS} pck.c

chelp.o: chelp.c
		${CC} ${CFLAGS} -I${DPS} chelp.c

mosflm_jpeg.o: mosflm_jpeg.c
		${CC} ${CFLAGS} -I../jpg mosflm_jpeg.c

# for HDF5 files
fromhdf5.o: fromhdf5.c
	${CC} -c  ${CFLAGS} fromhdf5.c

lz4.o:  ../lz4/lz4.c
	${CC} -c ${CFLAGS} -I../lz4 ../lz4/lz4.c

h5zlz4.o:       ../lz4/h5zlz4.c
	${CC} -c -std=c99 ${CFLAGS} -I../lz4 ../lz4/h5zlz4.c

bshuf_h5filter.o:        ../bitshuffle/bshuf_h5filter.c
	${CC} -c ${CFLAGS} -I../lz4 ../bitshuffle/bshuf_h5filter.c

bshuf_h5plugin.o:        ../bitshuffle/bshuf_h5plugin.c
	${CC} -c ${CFLAGS} -I../lz4 ../bitshuffle/bshuf_h5plugin.c

bitshuffle.o:    ../bitshuffle/bitshuffle.c
	${CC} -c -std=c99 ${CFLAGS} -I../lz4 ../bitshuffle/bitshuffle.c

clean-compile:
	rm -f ${OBJS} core *core

clean-bin:
	-test -z ${DPSBIN}/ipmosflm || rm -f ${DPSBIN}/ipmosflm

clean-generic:
	-rm -f *~

clean: clean-generic clean-compile clean-bin

distclean-generic:
	-rm -f Makefile

distclean: distclean-generic clean

realclean: distclean
	-rm -f ${SOURCES} ${DIST-OTHERS}

uninstall:

check:
# This only checks the syntax of this Makefile!
	cat -v -t -e Makefile|more

