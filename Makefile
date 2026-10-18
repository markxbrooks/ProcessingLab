# ProcessingLab top-level Makefile (Phase D / D6).
#
# Default: C++ core + thin ipmosflm (no Fortran amalgamation).
# Oracle:  make fortran-oracle   (or MOSFLM_FORTRAN_ORACLE=1)
#
# Autotools variables below remain for the optional oracle path.
#-----------------------------------------------------------------------------
prefix =  /Applications/ccp4-9
exec_prefix = ${prefix}
bindir  = ${exec_prefix}/bin
MOSHOME = .
RANLIB  = ranlib
INSTALL = /opt/homebrew/bin/ginstall -c
INSTALL_PROGRAM = ${INSTALL}

CPP_BUILD = mosflm-cpp/build
bin_PROGRAMS = ipmosflm

.PHONY: all cpp cbf-lib dps-lib fortran-oracle libs force install uninstall \
	clean distclean realclean test

all: cpp

# Vendored CBFlib (required by mosflm-cpp CMake).
cbf-lib:
	cd cbf && ${MAKE}
	if [ -f "$(RANLIB)" ] || command -v ranlib >/dev/null 2>&1; then \
		ranlib cbf/lib/libcbf.a 2>/dev/null || true; \
	fi

# DPS archives for BravaisClassifier / AutoIndexer.
dps-lib:
	cd index && ${MAKE}
	cd src/dps/util && ${MAKE}
	cd src/dps/index && ${MAKE}

# Default release artifact: mosflm-cpp CLI + thin CCP4-compatible ipmosflm.
cpp: cbf-lib dps-lib
	cmake -S mosflm-cpp -B $(CPP_BUILD) -DCMAKE_BUILD_TYPE=Release \
		-DMOSFLM_FORTRAN_ORACLE=OFF
	cmake --build $(CPP_BUILD) --target mosflm_cpp_cli ipmosflm
	mkdir -p bin
	cp -f $(CPP_BUILD)/ipmosflm bin/ipmosflm
	cp -f $(CPP_BUILD)/mosflm_cpp_cli bin/mosflm_cpp_cli
	@echo "Installed thin binaries under bin/ (no amalgamation)."

test: cpp
	ctest --test-dir $(CPP_BUILD) --output-on-failure

# Optional Fortran science amalgamation (nightly / differentials).
libs: cbf-lib
	cd index ; ${MAKE}
	if [ -f $(RANLIB) ]; then $(RANLIB) lib/DPSLIB.a; fi

fortran-oracle: libs
	cd mosflm ; ${MAKE} MOSFLM_FORTRAN_ORACLE=1 ipmosflm

force:
	rm -f bin/ipmosflm
	${MAKE} cpp

# Back-compat alias: ``make ipmosflm`` → thin C++ binary.
ipmosflm: cpp

bin/ipmosflm: cpp

distclean-generic:
	-rm -f ${MAKE}file Makefile
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*
	-rm -f lib/*.a

clean-generic:
	-rm -f *~
	-rm -f lib/*.a

clean: clean-generic
	cd index ; ${MAKE} clean || true
	cd cbf ; ${MAKE} clean || true
	cd mosflm ; ${MAKE} clean || true
	rm -rf $(CPP_BUILD)
	rm -f bin/ipmosflm bin/mosflm_cpp_cli

distclean: distclean-generic clean
	-rm -f config.status config.h
	cd index ; ${MAKE} distclean || true
	cd cbf ; ${MAKE} distclean || true
	cd mosflm ; ${MAKE} distclean || true

realclean: clean

install: cpp
	mkdir -p "$(DESTDIR)$(bindir)"
	$(INSTALL_PROGRAM) bin/ipmosflm "$(DESTDIR)$(bindir)"
	$(INSTALL_PROGRAM) bin/mosflm_cpp_cli "$(DESTDIR)$(bindir)"

uninstall:
	rm -f "$(DESTDIR)$(bindir)/ipmosflm" "$(DESTDIR)$(bindir)/mosflm_cpp_cli"
