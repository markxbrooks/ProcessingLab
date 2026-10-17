# Process this file with autoconf to produce a configure script.
AC_INIT([ipmosflm], [7.3.0], [harry@lmb-mrc.cam.ac.uk])
AC_CONFIG_SRCDIR([mosflm/mosflm_all_ip_inc.f])
AC_CONFIG_HEADER([config.h])
 
AC_ARG_WITH(ccp4-prefix,
  AS_HELP_STRING([--with-ccp4-prefix=PFX], [set ccp4 libraries prefix]),
  [if test "x${withval}" != xyes; then
   ccp4_prefix=$withval
   else
AC_MSG_ERROR([incorrect use of --with-ccp4-prefix])
   fi ],
   ccp4_prefix=""
   )

AC_ARG_VAR(CCP4, [CCP4 top dir])

if test "x$ccp4_prefix" != x; then
  CCP4=$ccp4_prefix
fi

if test "$CCP4"
then
  AC_PREFIX_DEFAULT([$CCP4])
  ac_status=$?
else
  AC_MSG_ERROR([ CCP4 library not defined ])
fi

#CCP4 uses FC for fortran
if test "$FC"; then
  F77=$FC
fi

AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_CXX_LIBRARY_LDFLAGS
AC_PROG_F77
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB

_AC_F77_NAME_MANGLING
dnl The alternatives in the cases here could be left over from a previous
dnl config.cache -- autoconf changed.
echo "$ac_cv_f77_mangling"--
case "$ac_cv_f77_mangling" in
  "upper case, no underscore, no extra underscore" | \
  "upper case, no underscores")
    AC_DEFINE([LINKTYP], [2],[fortran linking]) ;;
  "lower case, no underscore, no extra underscore" | \
  "lower case, no underscores" )
    AC_DEFINE([LINKTYP], [3]) ;;
  "lower case, underscore, no extra underscore" | \
  "lower case, single underscores")
    AC_DEFINE([LINKTYP], [1]) ;;
  "lower case, underscore, extra underscore" | \
  "lower case, double underscores")
    AC_DEFINE([LINKTYP], [1])
    # We aren't set up for f2c-style double underscoring.
    if test "$G77" = yes; then
      AC_MSG_WARN([
*** You need to use option -fno-second-underscore with g77
*** to compile Fortran to link with xdl_view.
])
      FFLAGS="-fno-second-underscore $FFLAGS"
    else
      AC_MSG_WARN([
*** Your Fortran compiler looks like f2c.
*** You can't currently link xdl_view routines with Fortran compiled with it.
*** You could use g77 with option `'-fno-second-underscore' instead.])
    fi ;;
  *)
    AC_MSG_ERROR([Unknown sort of Fortran name mangling]) ;;
esac

MOSFLM_F77_INT1()
AM_CONDITIONAL(INT_BYTE, test "$ac_mosflm_f77_int" = yes)
AM_CONDITIONAL(WITH_X,[test "x$with_x" = xyes])
# Checks for libraries.

AC_CHECK_LIB([tinfo], [tgoto])
AC_CHECK_LIB([ncurses], [initscr])
# test for curses if no ncurses
if test $ac_cv_lib_ncurses_initscr != yes; then
AC_CHECK_LIB([curses], [initscr])
fi

AC_CHECK_LIB([m], [cos])

AC_CHECK_LIB([jpeg],[jpeg_start_compress])

AC_CHECK_LIB([hdf5],[H5Fopen])

AC_CHECK_LIB([hdf5_hl],[H5DSset_scale])

# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h fcntl.h malloc.h sys/malloc.h malloc/malloc.h netdb.h netinet/in.h stddef.h string.h strings.h sys/socket.h unistd.h])

AC_CHECK_HEADERS([hdf5.h])

# Define size_t to a suitable type, if standard headers do not define it.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([bzero gethostbyname gethostname inet_ntoa rint socket])

#system services
AC_PATH_XTRA

AC_MOSFLM_OPTIONS

CCP4_CHECK_LIB([$EXTRA_LIBS $CXXLIBS])

if test "x$with_x" = xyes; then
CCP4_CHECK_XDL_VIEW
X_EXTRA_LIBS="-lXt -lX11 $X_EXTRA_LIBS"
fi

AC_CHECK_HEADER([winsock2.h],
		[AC_DEFINE([HAVE_LIBWS2_32],[1],[Define if libws2_32 exists.])
		LIBS="-lws2_32 $LIBS"],
		[],
		[#include <windows.h>])
AC_SUBST(AR_FLAGS)
AC_SUBST(TIDY)
AC_SUBST(CCP4_LIBS)
AC_SUBST(XDL_LIB)
AC_SUBST(EXTRA_LIBS)
AC_SUBST(ipmosflm_LDFLAGS)

AC_CONFIG_SUBDIRS([cbf])

# Provide some information about the compiler.
echo "$as_me:$LINENO:" \
     "mosflm header $F77 compiler version" >&5
if { ( eval echo "$as_me:$LINENO: \"$F77 --version </dev/null >&5\"") >&5
  (eval $F77 --version </dev/null >&5) 2>&5
  ac_status=$?
  echo "$as_me:$LINENO: \$? = $ac_status" >&5
  (exit $ac_status); } ; then
  ac_fortran_vers=`$F77 --version`
elif { (eval echo "$as_me:$LINENO: \"$F77 -v </dev/null >&5\"") >&5
  (eval $F77 -v </dev/null >&5) 2>&5
  ac_status=$?
  echo "$as_me:$LINENO: \$? = $ac_status" >&5
  (exit $ac_status); } ; then
  ac_fortran_vers="$ac_status"
elif { (eval echo "$as_me:$LINENO: \"$F77 -V </dev/null >&5\"") >&5
  (eval $F77 -V </dev/null >&5) 2>&5
  ac_status=$?
  echo "$as_me:$LINENO: \$? = $ac_status" >&5
  (exit $ac_status); } ; then
  ac_fortran_vers="$ac_status"
elif { (eval echo "$as_me:$LINENO: \"$F77 -qversion </dev/null >&5\"") >&5
  (eval $F77 -qversion </dev/null >&5) 2>&5
  ac_status=$?
  echo "$as_me:$LINENO: \$? = $ac_status" >&5
  (exit $ac_status); } ; then
  ac_fortran_vers="$ac_status"
else
  ac_fortran_vers="$unknown"
fi
FORTRAN_VERSION=`echo $ac_fortran_vers | sed '1!d'`

cat >conftest.f <<_ACEOF
      PROGRAM CONFTEST
      REAL A
      A = 1.0
      STOP
      END
_ACEOF
 echo
 test_compile="$F77 $FFLAGS -o conftest.o conftest.f"
 if { (eval echo "$as_me:$LINENO: mosflm header ABI using: \"$F77 $FFLAGS\"" ) >&5
  (eval $test_compile ) 2>&1
  ac_status=$?
  (exit $ac_status); }; then
      FMOSFLM_EXEC=`file conftest.o | sed 's/conftest.o://'`
  else
      FMOSFLM_EXEC="unknown"
  fi
  rm -rf conftest*
echo "$as_me:$LINENO: mosflm header ABI: $FMOSFLM_EXEC" >&5

BUILD_DATE=`date +"%A, %B %d %Y at %H:%M"`
echo "$as_me:$LINENO: mosflm header build date: $BUILD_DATE" >&5

BUILD_USER=`echo $USER`
echo "$as_me:$LINENO: mosflm header builder: $BUILD_USER" >&5

AC_LANG_PUSH(Fortran 77)
AX_CHECK_COMPILE_FLAG([-fallow-argument-mismatch],[FFLAGS="$FFLAGS -fallow-argument-mismatch"],, [-Werror])
AC_LANG_POP(Fortran 77)

AC_DEFINE_UNQUOTED([FORTRAN_COMPILER],["$F77"],[fortran compiler])
AC_DEFINE_UNQUOTED([FORTRAN_VERSION],["$FORTRAN_VERSION"],[fortran compiler version])
AC_DEFINE_UNQUOTED([FMOSFLM_EXEC],["$FMOSFLM_EXEC"],[executable file information])
AC_DEFINE_UNQUOTED([BUILD_DATE],["$BUILD_DATE"],[build date])
AC_DEFINE_UNQUOTED([BUILD_USER],["$BUILD_USER"],[build user])

AC_CONFIG_FILES([Makefile
                 index/Makefile
                 mosflm/Makefile
                 src/dps/index/Makefile
                 src/dps/util/Makefile
                 src/dps/util/filec/Makefile
                 src/dps/util/fftpack/Makefile])
AC_OUTPUT
