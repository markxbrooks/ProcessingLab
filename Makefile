# $Id: Makefile.in,v 1.15 2008/08/06 10:51:49 ccb Exp $
#
# Makefile for DPS to MOSFLM, Harry Powell
#
# This goes in the top level directory of the MOSFLM distribution with DPS
# indexing.
#
#-----------------------------------------------------------------------------
prefix =  /Applications/ccp4-9
exec_prefix = ${prefix}
bindir  = ${exec_prefix}/bin
MOSHOME = .
RANLIB  = ranlib
INSTALL = /opt/homebrew/bin/ginstall -c
INSTALL_PROGRAM = ${INSTALL}

bin_PROGRAMS = ipmosflm

all:	bin/ipmosflm

libs: 
	cd index ; ${MAKE}
	if [ -f $(RANLIB) ]; then $(RANLIB) lib/DPSLIB.a; fi
	cd cbf ; ${MAKE} ; if [ -f $(RANLIB) ]; \
	then $(RANLIB) lib/libcbf.a; fi

force:
	rm bin/ipmosflm
	${MAKE} ipmosflm

bin/ipmosflm:
	${MAKE} ipmosflm

ipmosflm: libs
	cd mosflm ; ${MAKE}

distclean-generic:
	-rm -f ${MAKE}file Makefile
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*
	-rm -f lib/*.a

clean-generic:
	-rm -f *~
	-rm -f lib/*.a

clean: clean-generic
	cd index ; ${MAKE} clean
	cd cbf ; ${MAKE} clean
	cd mosflm ; ${MAKE} clean	

distclean: distclean-generic
	-rm -f config.status config.h
	cd index ; ${MAKE} distclean
	cd cbf ; ${MAKE} distclean
	cd mosflm ; ${MAKE} distclean

realclean: clean

install: all
	mkdir -p "$(DESTDIR)$(bindir)"
	$(INSTALL_PROGRAM) bin/$(bin_PROGRAMS) "$(DESTDIR)$(bindir)"

uninstall: 
	cd index ; ${MAKE} uninstall
	cd cbf ; ${MAKE} uninstall
	cd mosflm ; ${MAKE} uninstall

DIST_COMMON = install-sh README.html README.txt README600.html README600.txt \
	README610.html README610.txt SUMMARY \
	aclocal.m4 config.h.in config.sub configure config.guess

EXTRA_DIST = doc/mosflm_user_guide.doc doc/mosflm_user_guide.html to_big.csh \
	to_small.csh Q314.csh 
