# ProcessingLab top-level Makefile.in (Phase D / E); ./configure → Makefile.
#
# Default: C++ core + ProcessingLab (no Fortran amalgamation).
#   make / make cpp
# Oracle amalgamation (F77 ipmosflm; not selected by env alone at this level):
#   make dps-lib fortran-oracle
#   # MOSFLM_FORTRAN_ORACLE=1 is only read inside mosflm/Makefile
#
# Autotools variables below remain for the optional oracle path.
#-----------------------------------------------------------------------------
prefix =  /Users/brooks/.local
exec_prefix = ${prefix}
bindir  = ${exec_prefix}/bin
MOSHOME = .
RANLIB  = ranlib
INSTALL = /opt/homebrew/bin/ginstall -c
INSTALL_PROGRAM = ${INSTALL}

# FFTW3 from ./configure (--with-fftw / --without-fftw).
HAVE_FFTW = yes
FFTW_CFLAGS = -I/opt/homebrew/Cellar/fftw/3.3.11/include
FFTW_LIBS = -L/opt/homebrew/Cellar/fftw/3.3.11/lib -lfftw3

CPP_BUILD = processing-lab/build
# C++ release binary; F77 oracle remains ``ipmosflm``.
ifeq ($(OS),Windows_NT)
  EXEEXT = .exe
else
  EXEEXT =
endif
PROCESSINGLAB = ProcessingLab$(EXEEXT)
bin_PROGRAMS = $(PROCESSINGLAB)

# Pass Autotools FFTW detection into the (legacy) processing-lab CMake bridge.
# MOSFLM_FFTW_FROM_CONFIGURE=ON stops CMake from re-searching when configure
# already decided (including --without-fftw).
ifeq ($(HAVE_FFTW),yes)
  FFTW_CMAKE_ARGS = -DMOSFLM_FFTW_FROM_CONFIGURE=ON -DMOSFLM_HAS_FFTW=ON \
	"-DMOSFLM_FFTW_CFLAGS=$(FFTW_CFLAGS)" \
	"-DMOSFLM_FFTW_LIBS=$(FFTW_LIBS)"
else
  FFTW_CMAKE_ARGS = -DMOSFLM_FFTW_FROM_CONFIGURE=ON -DMOSFLM_HAS_FFTW=OFF
endif

.PHONY: all cpp cpp-core cbf-lib dps-lib fortran-oracle libs force install \
	uninstall clean distclean realclean test

all: cpp

# Vendored CBFlib (required by processing-lab CMake).
cbf-lib:
	cd third_party/cbf && ${MAKE}
	if [ -f "$(RANLIB)" ] || command -v ranlib >/dev/null 2>&1; then \
		ranlib third_party/cbf/lib/libcbf.a 2>/dev/null || true; \
	fi

# DPS archives for BravaisClassifier / AutoIndexer.
# Needs gfortran (invert.f, fftpack/*.f, index/error_exit.f).
dps-lib:
	@if ! command -v gfortran >/dev/null 2>&1; then \
		echo "error: gfortran not found (required for DPS util/fftpack)."; \
		echo "  Ubuntu/Debian:  sudo apt install gfortran"; \
		echo "  Fedora/RHEL:    sudo dnf install gcc-gfortran"; \
		echo "  Or skip indexing and build C++ only:  make cpp-core"; \
		exit 1; \
	fi
	cd index && ${MAKE}
	cd src/dps/util && ${MAKE}
	cd src/dps/index && ${MAKE}

# C++ CONTROL binary without DPS archives (no Bravais / AutoIndex tests).
# FFTW is detected by ./configure (HAVE_FFTW=yes) and injected here;
# do not rely on CMake's own FFTW search when building via Autotools.
cpp-core: cbf-lib
	@if [ "$(HAVE_FFTW)" = yes ]; then \
		echo "Using FFTW3 from configure: $(FFTW_LIBS)"; \
	else \
		echo "FFTW3 not enabled (configure --with-fftw / apt install libfftw3-dev)"; \
	fi
	cmake -S processing-lab -B $(CPP_BUILD) -DCMAKE_BUILD_TYPE=Release \
		-DMOSFLM_FORTRAN_ORACLE=OFF \
		$(FFTW_CMAKE_ARGS)
	cmake --build $(CPP_BUILD) --target processinglab_dev ProcessingLab
	mkdir -p bin
	# Copy to a temp name then mv: overwriting a Mach-O in place on macOS
	# leaves a stale code-signature cache and the binary dies with SIGKILL.
	cp -f $(CPP_BUILD)/$(PROCESSINGLAB) bin/$(PROCESSINGLAB).tmp && \
		mv -f bin/$(PROCESSINGLAB).tmp bin/$(PROCESSINGLAB)
	# processinglab_dev (developer diagnostics) stays in $(CPP_BUILD).
	@echo "Installed thin binaries under bin/ (ProcessingLab; no amalgamation)."

# Default release: C++ + DPS when gfortran is available.
cpp: cbf-lib
	@if command -v gfortran >/dev/null 2>&1; then \
		$(MAKE) dps-lib; \
	else \
		echo "warning: gfortran not found — skipping dps-lib (no Bravais/AutoIndex)."; \
		echo "  Install with:  sudo apt install gfortran"; \
		echo "  Then re-run:   make dps-lib"; \
	fi
	$(MAKE) cpp-core

# ``cpp`` builds only the shipped binaries; build the test programs too.
test: cpp
	cmake --build $(CPP_BUILD)
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
	cd third_party/cbf ; ${MAKE} clean || true
	cd mosflm ; ${MAKE} clean || true
	rm -rf $(CPP_BUILD)
	rm -f bin/$(PROCESSINGLAB) bin/mosflm_cpp_cli$(EXEEXT) bin/ipmosflm

distclean: distclean-generic clean
	-rm -f config.status config.h
	cd index ; ${MAKE} distclean || true
	cd third_party/cbf ; ${MAKE} distclean || true
	cd mosflm ; ${MAKE} distclean || true

realclean: clean

# Installs what ``make`` built; never builds, so ``sudo make install`` leaves
# no root-owned files in the tree.
install:
	@test -x bin/$(PROCESSINGLAB) || { \
		echo "error: bin/$(PROCESSINGLAB) not built; run 'make' first (without sudo)."; \
		exit 1; \
	}
	mkdir -p "$(DESTDIR)$(bindir)"
	$(INSTALL_PROGRAM) bin/$(PROCESSINGLAB) "$(DESTDIR)$(bindir)"

# Also removes mosflm_cpp_cli left by older installs.
uninstall:
	rm -f "$(DESTDIR)$(bindir)/$(PROCESSINGLAB)" \
		"$(DESTDIR)$(bindir)/mosflm_cpp_cli$(EXEEXT)"

DIST_COMMON = install-sh README.html README.txt README600.html README600.txt \
	README610.html README610.txt SUMMARY \
	aclocal.m4 config.h.in config.sub configure config.guess

EXTRA_DIST = doc/mosflm_user_guide.doc doc/mosflm_user_guide.html to_big.csh \
	to_small.csh Q314.csh
