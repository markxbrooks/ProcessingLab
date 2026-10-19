# ProcessingLab top-level Makefile (Phase D / E).
#
# Default: C++ core + ProcessingLab (no Fortran amalgamation).
#   make / make cpp
# Oracle amalgamation (F77 ipmosflm; not selected by env alone at this level):
#   make dps-lib fortran-oracle
#   # MOSFLM_FORTRAN_ORACLE=1 is only read inside mosflm/Makefile
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

CPP_BUILD = processing-lab/build
# C++ release binary; F77 oracle remains ``ipmosflm``.
ifeq ($(OS),Windows_NT)
  EXEEXT = .exe
else
  EXEEXT =
endif
PROCESSINGLAB = ProcessingLab$(EXEEXT)
bin_PROGRAMS = $(PROCESSINGLAB)

.PHONY: all cpp cbf-lib dps-lib fortran-oracle libs force install uninstall \
	clean distclean realclean test

all: cpp

# Vendored CBFlib (required by processing-lab CMake).
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

# Default release artifact: processing-lab CLI + ProcessingLab (C++ CONTROL).
cpp: cbf-lib dps-lib
	cmake -S processing-lab -B $(CPP_BUILD) -DCMAKE_BUILD_TYPE=Release \
		-DMOSFLM_FORTRAN_ORACLE=OFF
	cmake --build $(CPP_BUILD) --target mosflm_cpp_cli ProcessingLab
	mkdir -p bin
	cp -f $(CPP_BUILD)/$(PROCESSINGLAB) bin/$(PROCESSINGLAB)
	cp -f $(CPP_BUILD)/mosflm_cpp_cli$(EXEEXT) bin/mosflm_cpp_cli$(EXEEXT)
	@echo "Installed thin binaries under bin/ (ProcessingLab; no amalgamation)."

test: cpp
	ctest --test-dir $(CPP_BUILD) --output-on-failure

# Optional Fortran science amalgamation (nightly / differentials).
libs: cbf-lib
	cd index ; ${MAKE}
	if [ -f $(RANLIB) ]; then $(RANLIB) lib/DPSLIB.a; fi

# Needs DPS archives (libdps_index.a); ``libs`` alone is insufficient.
# Installs F77 oracle as bin/ipmosflm (does not overwrite ProcessingLab).
fortran-oracle: libs dps-lib
	cd mosflm ; ${MAKE} MOSFLM_FORTRAN_ORACLE=1 ipmosflm

force:
	rm -f bin/$(PROCESSINGLAB)
	${MAKE} cpp

# Back-compat alias: ``make ipmosflm`` → C++ ProcessingLab.
ipmosflm: cpp
ProcessingLab: cpp

bin/$(PROCESSINGLAB): cpp

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
	rm -f bin/$(PROCESSINGLAB) bin/mosflm_cpp_cli$(EXEEXT) bin/ipmosflm

distclean: distclean-generic clean
	-rm -f config.status config.h
	cd index ; ${MAKE} distclean || true
	cd cbf ; ${MAKE} distclean || true
	cd mosflm ; ${MAKE} distclean || true

realclean: clean

install: cpp
	mkdir -p "$(DESTDIR)$(bindir)"
	$(INSTALL_PROGRAM) bin/$(PROCESSINGLAB) "$(DESTDIR)$(bindir)"
	$(INSTALL_PROGRAM) bin/mosflm_cpp_cli$(EXEEXT) "$(DESTDIR)$(bindir)"

uninstall:
	rm -f "$(DESTDIR)$(bindir)/$(PROCESSINGLAB)" \
		"$(DESTDIR)$(bindir)/mosflm_cpp_cli$(EXEEXT)"
