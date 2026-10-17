AC_DEFUN([AC_MOSFLM_OPTIONS],
[
#our specifics
#module specific flags
MOSFLM_EXTRA_FFLAGS=""

case "$host" in
  *osf4* )
    if test "`basename $CC`" = cc; then
      if ! (echo "$CFLAGS") | grep std1 >/dev/null ; then
        CFLAGS="$CFLAGS -std1"
      fi
      DEFS_SYS="-Dalpha"
    elif test "`basename $CC`" = c89; then
      DEFS_SYS="-Dalpha"
    else
      DEFS_SYS="-D__OSF1__"
    fi
    if test "`basename $F77`" = f77; then
      FFLAGS="$FFLAGS -fpe0 -assume dummy_aliases -nowarn"
    elif test "`basename $F77`" = g77; then
      FFLAGS="$FFLAGS -fno-globals -fno-automatic"
    fi
    X_EXTRA_FLAGS="-dnet_stub $X_EXTRA_FLAGS"
    AR_FLAGS="ruv"
    TIDY=${TIDY:-"rm -f ________64ELEL_"}
    ;;
  *osf5* )
    if test "`basename $CC`" = cc; then
      CFLAGS="$CFLAGS -w -DPROTOTYPE"
      DEFS_SYS="-Dalpha"
    else
      CFLAGS="$CFLAGS -fpe0 -assume dummy_aliases -assume no2underscores -nowarn"
      DEFS_SYS="-D__OSF1__"
    fi 
    FFLAGS="$FFLAGS -fpe0 -assume dummy_aliases -nowarn"
    X_EXTRA_FLAGS="-dnet_stub $X_EXTRA_FLAGS"
    AR_FLAGS="ruv"
    TIDY=${TIDY:-"rm -f ________64ELEL_"}
    ;;
  *sun* )
   AR_FLAGS="ruv"
   TIDY=${TIDY:-":"}
  ;;
  *irix6* )
   FFLAGS="$FFLAGS -w -Nn12000 -Nq12000"
   AR_FLAGS="vrus"
   TIDY=${TIDY:-"rm -f ________64ELEL_"}
   VERBOSE=""
   DEFS_SYS="-Dsgi"
    # Probably consistent with libccp4 (should have a hairy test for this):
    if test "`basename $CC`" = cc; then
      case "$CFLAGS" in
      *-64* ) ;;   # assume they meant it
      *-o32* | *-32* ) ;;   # presumably OK
      *-n32* ) ;;
      * ) 
      CFLAGS="$CFLAGS -n32"
      FFLAGS="$FFLAGS -n32"
      AC_MSG_WARN([

*** Using compiler flag -n32 for probable compatibility with
*** the CCP4 library.  You should check this, and re-run configure
*** specifying an appropriate value of CFLAGS if necessary.
*** See the INSTALL file for how to do so.
])
        ;;
    esac
    elif test "`basename $ac_cv_prog_CC`" = "gcc"; then
      AC_MSG_WARN([

*** Using gcc configured for Irix 6 won't be compatible with
*** the CCP4 library if that was compiled with `cc -o32'.  You should
*** check this and re-run configure specifying appropriate values of CC
*** and CFLAGS if necessary.  The `--without-gcc' option may be useful.
*** See the INSTALL file for more info.
])
    fi ;;
  *linux* )
    case $CC in
    *ccc )
      CFLAGS="$CFLAGS -O -c -tune generic" ;;
    *icc | *ecc)
      CFLAGS="$CFLAGS -c -w"
       case "$CFLAGS" in
         *-O* | *-O1* | *-O2* | *-O3* )
          CFLAGS=`echo $CFLAGS | sed s%-O[[\ 123]]%-O0\ %g`
          ;;
         *-O0* ) ;;
         *) CFLAGS="$CFLAGS -O0"
         esac ;;
    *gcc)
      case "$CFLAGS" in
         *-O* | *-O1* | *-O2* | *-O3* )
          CFLAGS=`echo $CFLAGS | sed s%-O[[\ 123]]%-O0\ %g`
          ;;
         *-O0* ) ;;
         *) CFLAGS="$CFLAGS -O0"
         esac
      esac

    ac_ext=f    
    ac_compile='$F77 $ac_cv_prog_f77_v -c $FFLAGS conftest.$ac_ext 2>&1'
cat >conftest.$ac_ext <<_ACEOF
      program main

      end
_ACEOF
    rm -f conftest.$ac_objext 

    case $F77 in
    *g77* )
      EXTRA_LIBS="$EXTRA_LIBS -lg2c"
      case "$FFLAGS" in
        *-fno-globals* ) ;;
        *) FFLAGS="$FFLAGS -fno-globals"
        esac
      case "$FFLAGS" in
        *fno-automatic* ) ;;
        *) FFLAGS="$FFLAGS -fno-automatic"
        esac
      case "$FFLAGS" in
        *fno-second-underscore* ) ;;
        *) FFLAGS="$FFLAGS -fno-second-underscore"
        esac
      case `$F77 -v 2>&1` in
      *2.96*) 
       case "$FFLAGS" in
         *-O* | *-O1* | *-O2* | *-O3* )
          FFLAGS=`echo $FFLAGS | sed s%-O[[\ 123]]%-O0\ %g`
          ;;
         *-O0* ) ;;
         *) FFLAGS="$FFLAGS -O0"
         esac
         ;;
      *) ;;
      esac ;;
    *gfortran* )
      EXTRA_LIBS="$EXTRA_LIBS -lgfortran"
      case "$FFLAGS" in
        *fno-second-underscore* ) ;;
        *) FFLAGS="$FFLAGS -fno-second-underscore"
        esac ;;
    *fort* ) 
      if (eval $ac_compile) | grep 'INTEL_COMPILER' >/dev/null ; then
         if ! (echo "$F77 $FFLAGS $LDFLAGS") | grep Vaxlib >/dev/null ; then
          LDFLAGS="$LDFLAGS -Vaxlib"
         fi
      else
         FFLAGS="$FFLAGS -O -fpe0 -assume dummy_aliases -assume no2underscores -c -tune generic" 
      fi
      ;;
    *ifc | *efc) 
      if ! (echo "$F77 $FFLAGS $LDFLAGS") | grep Vaxlib >/dev/null ; then
        LDFLAGS="$LDFLAGS -Vaxlib"
      fi
    esac 
    DEFS_SYS="-D__linux__"
    AR_FLAGS="ruv"
    TIDY=${TIDY:-"rm -f sta*"}
    VERBOSE=""
    ;;
  *darwin* )
    case $CC in
       *icc )
         case "$CFLAGS" in
           *-O* | *-O1* | *-O2* | *-O3* )
             CFLAGS=`echo $CFLAGS | sed s%-O[[\ 123]]%-O0\ %g`
             ;;
           *-O0* ) ;;
           *) CFLAGS="$CFLAGS -O0"
         esac ;;
       *xlc* )
         case "$CFLAGS" in
           *-O* | *-O1* | *-O2* | *-O3* )
             CFLAGS=`echo $CFLAGS | sed s%-O[[\ 123]]%-O0\ %g`
             ;;
           *-O0* ) ;;
           *) CFLAGS="$CFLAGS -O0"
         esac ;;
      *gcc* )
         case "$CFLAGS" in
           *-O* | *-O1* | *-O2* | *-O3* )
             CFLAGS=`echo $CFLAGS | sed s%-O[[\ 123]]%-O0\ %g`
             ;;
           *-O0* ) ;;
           *) CFLAGS="$CFLAGS -O0"
         esac
      esac
    case $F77 in
    *gfortran* )
      EXTRA_LIBS="$EXTRA_LIBS -lgfortran"
      case "$FFLAGS" in
        *fno-second-underscore* ) ;;
        *) FFLAGS="$FFLAGS -fno-second-underscore"
        esac ;;
    *g77 )
      EXTRA_LIBS="$EXTRA_LIBS -lg2c"
      case "$FFLAGS" in
        *-fno-globals* ) ;;
        *) FFLAGS="$FFLAGS -fno-globals"
        esac
      case "$FFLAGS" in
        *fno-automatic* ) ;;
        *) FFLAGS="$FFLAGS -fno-automatic"
        esac
      case "$FFLAGS" in 
        *fno-second-underscore* ) ;;
        *) FFLAGS="$FFLAGS -fno-second-underscore"
        esac
    case `$F77 -v 2>&1` in
     *3.1*) 
       case "$FFLAGS" in
         *-O* | *-O1* | *-O2* | *-O3* )
          FFLAGS=`echo $FFLAGS | sed s%-O[[\ 123]]%-O0\ %g`
          ;;
         *-O0* ) ;;
         *) FFLAGS="$FFLAGS -O0"
         esac
         ;;
     *3.4*)
       case "$FFLAGS" in
         *funroll-loops* )
          FFLAGS=`echo $FFLAGS | sed s%funroll-loops%fno-unroll-loops%g`
          ;;
         *) FFLAGS="$FFLAGS -fno-unroll-loops"
         esac
         ;;
     *)
     esac ;;
    *xlf* )
      case "$FFLAGS" in
        *-qextname* ) ;;
        *) FFLAGS="$FFLAGS -qextname"
        esac ;;
    *ifc )
      if ! (echo "$F77 $FFLAGS $LDFLAGS") | grep Vaxlib >/dev/null ; then
        LDFLAGS="$LDFLAGS -Vaxlib"
      fi
      ;;

    esac 
    DEFS_SYS="-D__APPLE__"
    AR_FLAGS="ruv"
    TIDY=${TIDY:-"rm -f sta*"}
    VERBOSE=""
    ;;
  *freebsd*)
    EXTRA_LIBS="$EXTRA_LIBS -lg2c"
    case $F77 in
    *gfortran* )
      EXTRA_LIBS="$EXTRA_LIBS -lgfortran"
      case "$FFLAGS" in
        *fno-second-underscore* ) ;;
        *) FFLAGS="$FFLAGS -fno-second-underscore"
        esac ;;
    *g77* )
      case "$FFLAGS" in
        *-fno-globals* ) ;;
        *) FFLAGS="$FFLAGS -fno-globals"
        esac
      case "$FFLAGS" in 
        *fno-automatic* ) ;;
        *) FFLAGS="$FFLAGS -fno-automatic"
        esac
      case "$FFLAGS" in 
        *fno-second-underscore* ) ;;
        *) FFLAGS="$FFLAGS -fno-second-underscore"
        esac
      ;;
      *)
    esac
    DEFS_SYS="-D__FreeBSD__"
    AR_FLAGS="ruv"
    TIDY=${TIDY:-"rm -f sta*"}
    VERBOSE=""
    ;;
  *)
esac

AC_SUBST(MOSFLM_EXTRA_FFLAGS)

])

AC_DEFUN([AC_F77_LIBRARY_LDFLAGS],
[AC_LANG_PUSH(Fortran 77)dnl    
_AC_PROG_F77_V
AC_CACHE_CHECK([for Fortran 77 libraries], ac_cv_flibs,
[if test "x$FLIBS" != "x"; then
  ac_cv_flibs="$FLIBS" # Let the user override the test.
else

_AC_PROG_F77_V_OUTPUT

ac_cv_flibs=

# Save positional arguments (if any)
ac_save_positional="$[@]"

set X $ac_f77_v_output
while test $[@%:@] != 1; do
  shift
  ac_arg=$[1]
  case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_flibs, ,
              ac_cv_flibs="$ac_cv_flibs $ac_arg")
          ;;
        -bI:*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_flibs, ,
             [_AC_LINKER_OPTION([$ac_arg], ac_cv_flibs)])
          ;;
          # Ignore these flags.
        -lang* | -lcrt0.o | -lc | -lgcc | -libmil | -LANG:=*)
          ;;
        -lcrt1.o | -lcrtbegin.o | -lcrt2.o | -lcrt1.10.[[5-6]].o )
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ac_cv_flibs="$ac_cv_flibs $ac_arg"
          ;;
        -[[LRuY]])
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
          case $[2] in
             "" | -*);;
             *)
                ac_arg="$ac_arg$[2]"
                shift; shift
                set X $ac_arg "$[@]"
                ;;
          esac
          ;;
        -YP,*)
          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
            _AC_LIST_MEMBER_IF($ac_j, $ac_cv_flibs, ,
                               [ac_arg="$ac_arg $ac_j"
                               ac_cv_flibs="$ac_cv_flibs $ac_j"])
          done
          ;;
        -[[lLR]]*)       
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_flibs, ,
                             ac_cv_flibs="$ac_cv_flibs $ac_arg")
          ;;
          # Ignore everything else.      
  esac
done
# restore positional arguments    
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`echo $ac_f77_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&
        _AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_flibs)
      ;;
esac
fi # test "x$FLIBS" = "x"

#workaround for -lompstubs ordering under solaris
if test "$GCC" != yes; then
  case $host_os in
    solaris* ) ac_cv_flibs="-xopenmp=stubs $ac_cv_flibs" ;;
    esac
fi

])
FLIBS="$ac_cv_flibs"       
AC_SUBST(FLIBS)
AC_LANG_POP(Fortran 77)dnl
])# AC_F77_LIBRARY_LDFLAGS


# _AC_PROG_F77_V_OUTPUT([FLAG = $ac_cv_prog_f77_v])
# -------------------------------------------------
# Link a trivial Fortran program, compiling with a verbose output FLAG
# (which default value, $ac_cv_prog_f77_v, is computed by
# _AC_PROG_F77_V), and return the output in $ac_f77_v_output.  This
# output is processed in the way expected by AC_F77_LIBRARY_LDFLAGS,
# so that any link flags that are echoed by the compiler appear as
# space-separated items.
AC_DEFUN([_AC_PROG_F77_V_OUTPUT],
[AC_REQUIRE([AC_PROG_F77])dnl
AC_LANG_PUSH(Fortran 77)dnl 
                
AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
        
# Compile and link our simple test program by passing a flag (argument
# 1 to this macro) to the Fortran 77 compiler in order to get
# "verbose" output that we can then parse for the Fortran 77 linker
# flags.
ac_save_FFLAGS=$FFLAGS 
FFLAGS="$FFLAGS m4_default([$1], [$ac_cv_prog_f77_v])"
(eval echo $as_me:__oline__: \"$ac_link\") >&AS_MESSAGE_LOG_FD
ac_f77_v_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1 | grep -v 'Driving:'`
echo "$ac_f77_v_output" >&AS_MESSAGE_LOG_FD
FFLAGS=$ac_save_FFLAGS
  
rm -f conftest*
AC_LANG_POP(Fortran 77)dnl
          
# If we are using xlf then replace all the commas with spaces.
if echo $ac_f77_v_output | grep xlfentry >/dev/null 2>&1; then
  ac_f77_v_output=`echo $ac_f77_v_output | sed 's/,/ /g'`
fi 

# If we are using Cray Fortran then delete quotes. 
# Use "\"" instead of '"' for font-lock-mode. 
# FIXME: a more general fix for quoted arguments with spaces?  
if echo $ac_f77_v_output | grep cft90 >/dev/null 2>&1; then
  ac_f77_v_output=`echo $ac_f77_v_output | sed "s/\"//g"`
fi

if echo $ac_f77_v_output | grep ifort >/dev/null 2>&1; then
  ac_f77_v_output=`echo $ac_f77_v_output | sed "s/\"//g"`
fi[]dnl
])# _AC_PROG_F77_V_OUTPUT
                        
      
# _AC_PROG_F77_V
# --------------
#
# Determine the flag that causes the Fortran 77 compiler to print
# information of library and object files (normally -v)
# Needed for AC_F77_LIBRARY_FLAGS
# Some compilers don't accept -v (Lahey: -verbose, xlf: -V, Fujitsu: -###)
AC_DEFUN([_AC_PROG_F77_V],
[AC_CACHE_CHECK([how to get verbose linking output from $F77],
                [ac_cv_prog_f77_v],
[AC_LANG_ASSERT(Fortran 77)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_f77_v=
# Try some options frequently used verbose output
for ac_verb in -v -verbose --verbose -V -\#\#\#; do
  _AC_PROG_F77_V_OUTPUT($ac_verb)
  # look for -l* and *.a constructs in the output
  for ac_arg in $ac_f77_v_output; do
     case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a | -[[lLRu]]*)
          ac_cv_prog_f77_v=$ac_verb
          break 2 ;;
     esac
  done
done
if test -z "$ac_cv_prog_f77_v"; then
   AC_MSG_WARN([cannot determine how to obtain linking information from $F77])
fi],
                  [AC_MSG_WARN([compilation failed])])
])])# _AC_PROG_F77_V

