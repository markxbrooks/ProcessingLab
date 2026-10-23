CC      = gcc -c
# Pre-C99 DPS sources call undeclared functions; keep that a warning.
CFLAGS  = -g -O0  -Wno-error=implicit-function-declaration
BASEDIR = ../../../..
LIBDIR  = $(BASEDIR)/lib
INCDIR  = .
LIBFILE = $(LIBDIR)/libfilec.a
RANLIB  = ranlib
SHELL   = /bin/sh
ARFLAGS = ruv
RM      = rm -f
CP      = cp

$(LIBFILE) : \
		$(LIBFILE)(wrfile.o) \
		$(LIBFILE)(wrmad.o) \
		$(LIBFILE)(wrfits.o) \
		$(LIBFILE)(wrswap.o) \
		$(LIBFILE)(wrhead.o) \
		$(LIBFILE)(wrcal.o) \
		$(LIBFILE)(rdcal.o) \
		$(LIBFILE)(rdfile.o) \
		$(LIBFILE)(rdfits.o) \
		$(LIBFILE)(rdsmv.o) \
		$(LIBFILE)(rdbinary.o) \
		$(LIBFILE)(rdmar.o) \
		$(LIBFILE)(rdlum.o) \
		$(LIBFILE)(rdmad.o) \
		$(LIBFILE)(rdhead.o) \
		$(LIBFILE)(rdctbl.o) \
		$(LIBFILE)(wrctbl.o) \
		$(LIBFILE)(dskhd.o) \
		$(LIBFILE)(dskio.o) \
		$(LIBFILE)(swpbyt.o) \
		$(LIBFILE)(getbo.o) \
		$(LIBFILE)(filec_debug.o) \
		$(LIBFILE)(jzext.o) \
		$(LIBFILE)(namfil.o)
	${MAKE} $?
	${AR} $(ARFLAGS) $@  $?
	if [ -f $(RANLIB) ]; then $(RANLIB) $@; fi
	$(RM) $?
.c.a:;
.f.a:;

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) $<
 
$(LIBFILE)(wrfile.o)  : wrfile.c $(INCDIR)/filec.h
$(LIBFILE)(wrmad.o)  : wrmad.c $(INCDIR)/filec.h
$(LIBFILE)(wrfits.o)  : wrfits.c $(INCDIR)/filec.h
$(LIBFILE)(wrswap.o)  : wrswap.c $(INCDIR)/filec.h
$(LIBFILE)(wrhead.o)  : wrhead.c $(INCDIR)/filec.h
$(LIBFILE)(wrcal.o)  : wrcal.c $(INCDIR)/filec.h
$(LIBFILE)(rdcal.o)  : rdcal.c $(INCDIR)/filec.h
$(LIBFILE)(rdfits.o)  : rdfits.c $(INCDIR)/filec.h
$(LIBFILE)(rdmar.o)  : rdmar.c $(INCDIR)/filec.h
$(LIBFILE)(rdsmv.o)  : rdsmv.c $(INCDIR)/filec.h
$(LIBFILE)(rdbinary.o)  : rdbinary.c $(INCDIR)/filec.h
$(LIBFILE)(rdlum.o)  : rdlum.c $(INCDIR)/filec.h
$(LIBFILE)(rdmad.o)  : rdmad.c $(INCDIR)/filec.h
$(LIBFILE)(rdhead.o)  : rdhead.c $(INCDIR)/filec.h
$(LIBFILE)(rdctbl.o)  : rdctbl.c $(INCDIR)/filec.h
$(LIBFILE)(wrctbl.o)  : wrctbl.c $(INCDIR)/filec.h
$(LIBFILE)(namfil.o)  : namfil.c $(INCDIR)/filec.h
$(LIBFILE)(dskhd.o)  : dskhd.c $(INCDIR)/filec.h
$(LIBFILE)(dskio.o)  : dskio.c $(INCDIR)/filec.h
$(LIBFILE)(swpbyt.o) : swpbyt.c $(INCDIR)/filec.h
$(LIBFILE)(filec_debug.o) : filec_debug.c $(INCDIR)/filec.h
$(LIBFILE)(getbo.o) : getbo.c $(INCDIR)/filec.h

clean :
	rm -f *~
	rm -f *.o
	rm -f file.tar*

clean-bin:
	-test -z ${LIBDIR}/libfilec.a || rm -f ${LIBDIR}/libfilec.a

distclean: clean-bin clean
	-rm -f Makefile

realclean: distclean
	-rm ${SOURCES} ${HEADERS} ${DIST-COMMON}

uninstall: clean-bin

dskhd : dskhd.c
	$(CC) $(CFLAGS) -O0 -g -o $@ $@.c -DTEST_DSKHD $(EXTLIBS)

file_test : file_test.c
	$(CC) $(CFLAGS) -o $@ $@.c $(EXTLIBS)

test_file : test_file.c
	$(CC) $(CFLAGS) -o $@ $@.c -I$(INCDIR) $(EXTLIBS) -ltermioc

test_dskhd : dskhd.c
	$(CC) $(CFLAGS) -DTEST_DSKHD -o $@ dskhd.c -I$(INCDIR) $(EXTLIBS)

recover: recover.c
	$(CC) $(CFLAGS) -o $@ $@.c $(EXTLIBS)

cuthead : cuthead.c
	$(CC)  $(CFLAGS) -o $@ $@.c -I$(INCDIR) $(EXTLIBS)

DIST-COMMON = Makefile.in documentation.txt

SOURCES = dskhd.c dskio.c file_test.c filec_debug.c getbo.c jzext.c \
	namfil.c rdcal.c rdctbl.c rdfile.c rdfits.c rdhead.c rdlum.c \
	rdmad.c rdmar.c rdsmv.c swpbyt.c wrcal.c wrctbl.c wrfile.c \
	wrfits.c wrhead.c wrmad.c wrswap.c rdbinary.c

HEADERS = filec.h
