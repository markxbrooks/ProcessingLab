dnl Ensure that the CCP4 library is linkable, else guess.
AC_DEFUN([CCP4_CHECK_LIB],
 [ ccp4_check_lib_save_LIBS=$LIBS
   ccp4_check_lib_save_LDFLAGS=$LDFLAGS
   AC_LANG_PUSH([Fortran 77])
   LIBS=""
   if test "$CLIB"
   then
     LDFLAGS="-L$CLIB"
   else
     LDFLAGS="-L$CCP4/lib"
   fi
# library tests
  AC_CHECK_LIB(ccp4f, qisnan, ,[AC_MSG_WARN([
*** Can't find the CCP4 library against which to link.
*** You can set the environment variable CLIB to point to the right
    directory.  Currently this is "$CLIB".  If necessary, check the
    file `config.log' for other problems.
*** Setting to $CCP4/lib ])
    LDFLAGS="-L$CCP4/lib $LDFLAGS"],[-lccp4c -lmmdb2])
  CCP4_LIBS="$LDFLAGS -lccp4f -lccp4c -lmmdb2"
  LIBS="$ccp4_check_lib_save_LIBS"
  LDFLAGS="$ccp4_check_lib_save_LDFLAGS"
  AC_LANG_POP([Fortran 77])
])

