AC_INIT([jpeg],[6.2])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE
AC_PROG_LIBTOOL
AC_PROG_CC
AC_HEADER_STDC
AC_PROG_INSTALL

AM_C_PROTOTYPES
AM_CONFIG_HEADER(jconfig.h)
AC_CHECK_HEADERS(stddef.h)
AC_TYPE_SIZE_T
AC_C_CONST
AC_C_INLINE
 
echo $ac_n "checking for type unsigned char""... $ac_c" 1>&6
echo "configure:$LINENO: checking for type unsigned char" >&5
cat > conftest.$ac_ext <<EOF
#line $LINENO "configure"
#include "confdefs.h"

int main() {
 unsigned char un_char; 
; return 0; }
EOF
if { (eval echo configure:$LINENO: \"$ac_compile\") 1>&5; (eval $ac_compile) 2>&5; }; then
  rm -rf conftest*
  echo "$ac_t""yes" 1>&6
cat >> confdefs.h <<\EOF
#define HAVE_UNSIGNED_CHAR 
EOF

else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  echo "$ac_t""no" 1>&6
fi
rm -f conftest*
echo $ac_n "checking for type unsigned short""... $ac_c" 1>&6
echo "configure:$LINENO: checking for type unsigned short" >&5
cat > conftest.$ac_ext <<EOF
#line $LINENO "configure"
#include "confdefs.h"

int main() {
 unsigned short un_short; 
; return 0; }
EOF
if { (eval echo configure:$LINENO: \"$ac_compile\") 1>&5; (eval $ac_compile) 2>&5; }; then
  rm -rf conftest*
  echo "$ac_t""yes" 1>&6
cat >> confdefs.h <<\EOF
#define HAVE_UNSIGNED_SHORT 
EOF

else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  echo "$ac_t""no" 1>&6
fi
rm -f conftest*

echo $ac_n "checking for type void""... $ac_c" 1>&6
echo "configure:$LINENO: checking for type void" >&5
cat > conftest.$ac_ext <<EOF
#line $LINENO "configure"
#include "confdefs.h"

/* Caution: a C++ compiler will insist on valid prototypes */
typedef void * void_ptr;	/* check void * */
#ifdef PROTOTYPES		/* check ptr to function returning void */
typedef void (*void_func) (int a, int b);
#else
typedef void (*void_func) ();
#endif

#ifdef PROTOTYPES		/* check void function result */
void test3function (void_ptr arg1, void_func arg2)
#else
void test3function (arg1, arg2)
     void_ptr arg1;
     void_func arg2;
#endif
{
  char * locptr = (char *) arg1; /* check casting to and from void * */
  arg1 = (void *) locptr;
  (*arg2) (1, 2);		/* check call of fcn returning void */
}

int main() {
 
; return 0; }
EOF
if { (eval echo configure:$LINENO: \"$ac_compile\") 1>&5; (eval $ac_compile) 2>&5; }; then
  rm -rf conftest*
  echo "$ac_t""yes" 1>&6
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  echo "$ac_t""no" 1>&6
cat >> confdefs.h <<\EOF
#define void char
EOF

fi
rm -f conftest*

echo $ac_n "checking for broken incomplete types""... $ac_c" 1>&6
echo "configure:$LINENO: checking for broken incomplete types" >&5
cat > conftest.$ac_ext <<EOF
#line $LINENO "configure"
#include "confdefs.h"
 typedef struct undefined_structure * undef_struct_ptr; 
int main() {

; return 0; }
EOF
if { (eval echo configure:$LINENO: \"$ac_compile\") 1>&5; (eval $ac_compile) 2>&5; }; then
  rm -rf conftest*
  echo "$ac_t""ok" 1>&6
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  echo "$ac_t""broken" 1>&6
cat >> confdefs.h <<\EOF
#define INCOMPLETE_TYPES_BROKEN 
EOF

fi
rm -f conftest*

echo $ac_n "checking for short external names""... $ac_c" 1>&6
echo "configure:$LINENO: checking for short external names" >&5
cat > conftest.$ac_ext <<EOF
#line $LINENO "configure"
#include "confdefs.h"

int possibly_duplicate_function () { return 0; }
int possibly_dupli_function () { return 1; }

int main() {
 
; return 0; }
EOF
if { (eval echo configure:$LINENO: \"$ac_link\") 1>&5; (eval $ac_link) 2>&5; } && test -s conftest; then
  rm -rf conftest*
  echo "$ac_t""ok" 1>&6
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  echo "$ac_t""short" 1>&6
cat >> confdefs.h <<\EOF
#define NEED_SHORT_EXTERNAL_NAMES 
EOF

fi
rm -f conftest*
echo $ac_n "checking to see if char is signed""... $ac_c" 1>&6
echo "configure:$LINENO: checking to see if char is signed" >&5
if test "$cross_compiling" = yes; then
  echo Assuming that char is signed on target machine.
echo If it is unsigned, this will be a little bit inefficient.

else
  cat > conftest.$ac_ext <<EOF
#line $LINENO "configure"
#include "confdefs.h"

#ifdef PROTOTYPES
int is_char_signed (int arg)
#else
int is_char_signed (arg)
     int arg;
#endif
{
  if (arg == 189) {		/* expected result for unsigned char */
    return 0;			/* type char is unsigned */
  }
  else if (arg != -67) {	/* expected result for signed char */
    printf("Hmm, it seems 'char' is not eight bits wide on your machine.\n");
    printf("I fear the JPEG software will not work at all.\n\n");
  }
  return 1;			/* assume char is signed otherwise */
}
char signed_char_check = (char) (-67);
main() {
  exit(is_char_signed((int) signed_char_check));
}
EOF
if { (eval echo configure:$LINENO: \"$ac_link\") 1>&5; (eval $ac_link) 2>&5; } && test -s conftest && (./conftest; exit) 2>/dev/null
then
  echo "$ac_t""no" 1>&6
cat >> confdefs.h <<\EOF
#define CHAR_IS_UNSIGNED 
EOF

else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -fr conftest*
  echo "$ac_t""yes" 1>&6
fi
rm -fr conftest*
fi

echo $ac_n "checking to see if right shift is signed""... $ac_c" 1>&6
echo "configure:$LINENO: checking to see if right shift is signed" >&5
if test "$cross_compiling" = yes; then
  echo "$ac_t""Assuming that right shift is signed on target machine." 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line $LINENO "configure"
#include "confdefs.h"

#ifdef PROTOTYPES
int is_shifting_signed (long arg)
#else
int is_shifting_signed (arg)
     long arg;
#endif
/* See whether right-shift on a long is signed or not. */
{
  long res = arg >> 4;

  if (res == -0x7F7E80CL) {	/* expected result for signed shift */
    return 1;			/* right shift is signed */
  }
  /* see if unsigned-shift hack will fix it. */
  /* we can't just test exact value since it depends on width of long... */
  res |= (~0L) << (32-4);
  if (res == -0x7F7E80CL) {	/* expected result now? */
    return 0;			/* right shift is unsigned */
  }
  printf("Right shift isn't acting as I expect it to.\n");
  printf("I fear the JPEG software will not work at all.\n\n");
  return 0;			/* try it with unsigned anyway */
}
main() {
  exit(is_shifting_signed(-0x7F7E80B1L));
}
EOF
if { (eval echo configure:$LINENO: \"$ac_link\") 1>&5; (eval $ac_link) 2>&5; } && test -s conftest && (./conftest; exit) 2>/dev/null
then
  echo "$ac_t""no" 1>&6
cat >> confdefs.h <<\EOF
#define RIGHT_SHIFT_IS_UNSIGNED 
EOF

else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -fr conftest*
  echo "$ac_t""yes" 1>&6
fi
rm -fr conftest*
fi

echo $ac_n "checking to see if fopen accepts b spec""... $ac_c" 1>&6
echo "configure:$LINENO: checking to see if fopen accepts b spec" >&5
if test "$cross_compiling" = yes; then
  echo "$ac_t""Assuming that it does." 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line $LINENO "configure"
#include "confdefs.h"

#include <stdio.h>
main() {
  if (fopen("conftestdata", "wb") != NULL)
    exit(0);
  exit(1);
}
EOF
if { (eval echo configure:$LINENO: \"$ac_link\") 1>&5; (eval $ac_link) 2>&5; } && test -s conftest && (./conftest; exit) 2>/dev/null
then
  echo "$ac_t""yes" 1>&6
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -fr conftest*
  echo "$ac_t""no" 1>&6
cat >> confdefs.h <<\EOF
#define DONT_USE_B_MODE 
EOF

fi
rm -fr conftest*
fi

# Select memory manager depending on user input.
# If no "-enable-maxmem", use jmemnobs
MEMORYMGR='jmemnobs'
MAXMEM="no"
# Check whether --enable-maxmem or --disable-maxmem was given.
if test "${enable_maxmem+set}" = set; then
  enableval="$enable_maxmem"
  MAXMEM="$enableval"
fi

# support --with-maxmem for backwards compatibility with IJG V5.
# Check whether --with-maxmem or --without-maxmem was given.
if test "${with_maxmem+set}" = set; then
  withval="$with_maxmem"
  MAXMEM="$withval"
fi

if test "x$MAXMEM" = xyes; then
  MAXMEM=1
fi
if test "x$MAXMEM" != xno; then
  if test -n "`echo $MAXMEM | sed 's/@<:@0-9@:>@//g'`"; then
    { echo "configure: error: non-numeric argument to --enable-maxmem" 1>&2; exit 1; }
  fi
  DEFAULTMAXMEM=`expr $MAXMEM \* 1048576`
cat >> confdefs.h <<EOF
#define DEFAULT_MAX_MEM ${DEFAULTMAXMEM}
EOF

echo $ac_n "checking for 'tmpfile()'""... $ac_c" 1>&6
echo "configure:$LINENO: checking for 'tmpfile()'" >&5
cat > conftest.$ac_ext <<EOF
#line $LINEO "configure"
#include "confdefs.h"
#include <stdio.h>
int main() {
 FILE * tfile = tmpfile(); 
; return 0; }
EOF
if { (eval echo configure:$LINENO: \"$ac_link\") 1>&5; (eval $ac_link) 2>&5; } && test -s conftest; then
  rm -rf conftest*
  echo "$ac_t""yes" 1>&6
MEMORYMGR='jmemansi'
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  echo "$ac_t""no" 1>&6
MEMORYMGR='jmemname'
cat >> confdefs.h <<\EOF
#define NEED_SIGNAL_CATCHER 
EOF

echo $ac_n "checking for 'mktemp()'""... $ac_c" 1>&6
echo "configure:$LINENO: checking for 'mktemp()'" >&5
cat > conftest.$ac_ext <<EOF
#line $LINENO "configure"
#include "confdefs.h"

int main() {
 char fname@<:@80@:>@; mktemp(fname); 
; return 0; }
EOF
if { (eval echo configure:$LINENO: \"$ac_link\") 1>&5; (eval $ac_link) 2>&5; } && test -s conftest; then
  rm -rf conftest*
  echo "$ac_t""yes" 1>&6
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  echo "$ac_t""no" 1>&6
cat >> confdefs.h <<\EOF
#define NO_MKTEMP 
EOF

fi
rm -f conftest*
fi
rm -f conftest*
fi


# Extract the library version ID from jpeglib.h.
echo $ac_n "checking libjpeg version number""... $ac_c" 1>&6
echo "configure:$LINENO: checking libjpeg version number" >&5
JPEG_LIB_VERSION=`sed -e '/^#define JPEG_LIB_VERSION/!d' -e 's/^[^0-9]*\([0-9][0-9]*\).*$/\1/' $srcdir/jpeglib.h`
echo "$ac_t""$JPEG_LIB_VERSION" 1>&6

# ansi2knr needs -DBSD if string.h is missing
if test $ac_cv_header_string_h = no; then
  ANSI2KNRFLAGS="-DBSD"
else
  ANSI2KNRFLAGS=""
fi

AC_SUBST(JPEG_LIB_VERSION)
AC_SUBST(ANSI2KNRFLAGS)
AC_SUBST(INSTALL_LIB)
AC_SUBST(MEMORYMGR)

AC_CONFIG_FILES(Makefile)

AC_OUTPUT
