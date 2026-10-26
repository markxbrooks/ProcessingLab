# Process this file with autoconf to produce a configure script.
AC_INIT([cbf], [0.7], [yaya@bernstein-plus-sons.com])
AC_CONFIG_SRCDIR([src/cbf.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([alloca.h limits.h malloc.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memmove memset pow strchr strspn strstr strtoul])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
