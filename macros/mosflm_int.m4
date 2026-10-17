dnl Ensure that the integer*1 is supported
AC_DEFUN([MOSFLM_F77_INT1],
 [ AC_LANG_PUSH(Fortran 77)dnl
   AC_LINK_IFELSE([AC_LANG_PROGRAM([],[
      integer*1 tmp
      tmp = 1
])],
      [ac_mosflm_f77_int=yes],[ac_mosflm_f77_int=no])
   AC_LANG_POP(Fortran 77)
])
