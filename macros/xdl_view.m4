dnl Ensure that the xdl_view library is linkable, else die.
AC_DEFUN([CCP4_CHECK_XDL_VIEW],
  [ ccp4_check_xdl_view_save_LIBS=$LIBS
    ccp4_check_xdl_view_save_LDFLAGS=$LDFLAGS
    LIBS=""
   if test "$CLIB"
   then
     LDFLAGS="-L$CLIB"
   else
     LDFLAGS="-L$CCP4/lib"
   fi
    LDFLAGS="$LDFLAGS $X_LIBS"
AC_CHECK_LIB(xdl_view, xdl_open_view, ,[ AC_MSG_WARN([
*** Can't find the xdl_view library against which to link.
*** You can set the environment variable CLIB to point to the right
    directory.  Currently this is "$CLIB".  If necessary, check the
    file `config.log' for other problems.
*** Setting to $CCP4/lib ])
    LDFLAGS="-L$CCP4/lib $X_LIBS"
    LIBS="-lxdl_view" ], -lXt -lX11) 
    XDL_LIB="$LDFLAGS $LIBS"
    LIBS="$ccp4_check_xdl_view_save_LIBS"
    LDFLAGS="$ccp4_check_xdl_view_save_LDFLAGS"    
])
