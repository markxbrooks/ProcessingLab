CFLAGS= @CFLAGS@ @CPPFLAGS@
# Generally, we recommend defining any configuration symbols in jconfig.h,
# NOT via -D switches here.
# However, any special defines for ansi2knr.c may be included here:
ANSI2KNRFLAGS= @ANSI2KNRFLAGS@

AUTOMAKE_OPTIONS = ansi2knr

# Library version ID; libtool uses this for the shared library version number.
# Note: we suggest this match the macro of the same name in jpeglib.h.
JPEG_LIB_VERSION = @JPEG_LIB_VERSION@

# Put here the object file name for the correct system-dependent memory
# manager file.  For Unix this is usually jmemnobs.o, but you may want
# to use jmemansi.o or jmemname.o if you have limited swap space.
SYSDEPMEM= @MEMORYMGR@

# source files: JPEG library proper
libjpeg_la_SOURCES = jcapimin.c jcapistd.c jccoefct.c jccolor.c jcdctmgr.c jchuff.c \
        jcinit.c jcmainct.c jcmarker.c jcmaster.c jcomapi.c jcparam.c \
        jcphuff.c jcprepct.c jcsample.c jctrans.c jdapimin.c jdapistd.c \
        jdatadst.c jdatasrc.c jdcoefct.c jdcolor.c jddctmgr.c jdhuff.c \
        jdinput.c jdmainct.c jdmarker.c jdmaster.c jdmerge.c jdphuff.c \
        jdpostct.c jdsample.c jdtrans.c jerror.c jfdctflt.c jfdctfst.c \
        jfdctint.c jidctflt.c jidctfst.c jidctint.c jidctred.c jquant1.c \
        jquant2.c jutils.c jmemmgr.c

libjpeg_a_SOURCES = $(libjpeg_la_SOURCES)

lib_LTLIBRARIES = libjpeg.la
lib_LIBRARIES = libjpeg.a
 
# memmgr back ends: compile only one of these into a working library
SYSDEPSOURCES= jmemansi.c jmemname.c jmemnobs.c jmemdos.c jmemmac.c

EXTRA_libjpeg_la_SOURCES = $(SYSDEPSOURCES)
EXTRA_libjpeg_a_SOURCES = $(EXTRA_libjpeg_la_SOURCES)

libjpeg_la_LIBADD = $(SYSDEPMEM).lo
libjpeg_la_DEPENDENCIES = $(SYSDEPMEM).lo
libjpeg_a_LIBADD = $(SYSDEPMEM).$(OBJEXT)
libjpeg_a_DEPENDENCIES = $(SYSDEPMEM).$(OBJEXT)

libjpeg_la_LDFLAGS = -version-info $(JPEG_LIB_VERSION):0
# files included by source files
noinst_HEADERS = jchuff.h jdhuff.h jdct.h jinclude.h jmemsys.h \
        jpegint.h jversion.h cdjpeg.h cderror.h transupp.h
include_HEADERS = jconfig.h jpeglib.h jmorecfg.h jerror.h

# documentation, test, and support files
DOCS= install.doc usage.doc \
        wizard.doc example.c libjpeg.doc structure.doc \
        coderules.doc filelist.doc
MKFILES= makefile.cfg makefile.ansi makefile.unix makefile.bcc \
        makefile.mc6 makefile.dj makefile.wat makefile.vc makelib.ds \
        makeapps.ds makeproj.mac makcjpeg.st makdjpeg.st makljpeg.st \
        maktjpeg.st makefile.manx makefile.sas makefile.mms makefile.vms \
        BJECTmakvms.opt
CONFIGFILES= jconfig.cfg jconfig.bcc jconfig.mc6 jconfig.dj jconfig.wat \
        jconfig.vc jconfig.mac jconfig.st jconfig.manx jconfig.sas \
        jconfig.vms
OTHERFILES= jconfig.doc ckconfig.c jmemdosa.asm
TESTFILES= testorig.jpg testimg.ppm testimg.bmp testimg.jpg testprog.jpg \
        testimgp.jpg
EXTRA_DIST = $(DOCS) $(MKFILES) $(CONFIGFILES) \
        $(OTHERFILES) $(TESTFILES)

# library object files common to compression and decompression
COMOBJECTS= jcomapi.$(OBJEXT) jutils.$(OBJEXT) jerror.$(OBJEXT) jmemmgr.$(OBJEXT) $(SYSDEPMEM)
# compression library object files
CLIBOBJECTS= jcapimin.$(OBJEXT) jcapistd.$(OBJEXT) jctrans.$(OBJEXT) jcparam.$(OBJEXT) \
        jdatadst.$(OBJEXT) jcinit.$(OBJEXT) jcmaster.$(OBJEXT) jcmarker.$(OBJEXT) jcmainct.$(OBJEXT) \
        jcprepct.$(OBJEXT) jccoefct.$(OBJEXT) jccolor.$(OBJEXT) jcsample.$(OBJEXT) jchuff.$(OBJEXT) \
        jcphuff.$(OBJEXT) jcdctmgr.$(OBJEXT) jfdctfst.$(OBJEXT) jfdctflt.$(OBJEXT) \
        jfdctint.$(OBJEXT)
# decompression library object files
DLIBOBJECTS= jdapimin.$(OBJEXT) jdapistd.$(OBJEXT) jdtrans.$(OBJEXT) jdatasrc.$(OBJEXT) \
        jdmaster.$(OBJEXT) jdinput.$(OBJEXT) jdmarker.$(OBJEXT) jdhuff.$(OBJEXT) jdphuff.$(OBJEXT) \
        jdmainct.$(OBJEXT) jdcoefct.$(OBJEXT) jdpostct.$(OBJEXT) jddctmgr.$(OBJEXT) \
        jidctfst.$(OBJEXT) jidctflt.$(OBJEXT) jidctint.$(OBJEXT) jidctred.$(OBJEXT) \
        jdsample.$(OBJEXT) jdcolor.$(OBJEXT) jquant1.$(OBJEXT) jquant2.$(OBJEXT) jdmerge.$(OBJEXT)
# These objectfiles are included in libjpeg.a
LIBOBJECTS= $(CLIBOBJECTS) $(DLIBOBJECTS) $(COMOBJECTS)

bin_PROGRAMS = cjpeg djpeg jpegtran rdjpgcom wrjpgcom

cjpeg_SOURCES = cjpeg.c rdppm.c rdgif.c rdtarga.c rdrle.c \
        rdbmp.c rdswitch.c cdjpeg.c
cjpeg_LDFLAGS = -ljpeg
djpeg_SOURCES = djpeg.c wrppm.c wrgif.c wrtarga.c wrrle.c \
        wrbmp.c rdcolmap.c cdjpeg.c
djpeg_LDFLAGS = -ljpeg
jpegtran_SOURCES = jpegtran.c rdswitch.c cdjpeg.c transupp.c
jpegtran_LDFLAGS = -ljpeg
rdjpgcom_SOURCES = rdjpgcom.c
rdjpgcom_LDFLAGS = -ljpeg
wrjpgcom_SOURCES = wrjpgcom.c
wrjpgcom_LDFLAGS = -ljpeg

dist_man_MANS = cjpeg.1 djpeg.1 jpegtran.1 rdjpgcom.1 wrjpgcom.1

test: cjpeg djpeg jpegtran
	$(RM) testout*
	./djpeg -dct int -ppm -outfile testout.ppm  $(srcdir)/testorig.jpg
	./djpeg -dct int -bmp -colors 256 -outfile testout.bmp  $(srcdir)/testorig.jpg
	./cjpeg -dct int -outfile testout.jpg  $(srcdir)/testimg.ppm
	./djpeg -dct int -ppm -outfile testoutp.ppm $(srcdir)/testprog.jpg
	./cjpeg -dct int -progressive -opt -outfile testoutp.jpg $(srcdir)/testimg.ppm
	./jpegtran -outfile testoutt.jpg $(srcdir)/testprog.jpg
	cmp $(srcdir)/testimg.ppm testout.ppm
	cmp $(srcdir)/testimg.bmp testout.bmp
	cmp $(srcdir)/testimg.jpg testout.jpg
	cmp $(srcdir)/testimg.ppm testoutp.ppm
	cmp $(srcdir)/testimgp.jpg testoutp.jpg
	cmp $(srcdir)/testorig.jpg testoutt.jpg

check: test

# Mistake catcher:

jconfig.h: jconfig.doc
	echo You must prepare a system-dependent jconfig.h file.
	echo Please read the installation directions in install.doc.
	exit 1

# GNU Make likes to know which target names are not really files to be made:
.PHONY: all install install-lib install-headers clean distclean test check


jcapimin.$(OBJEXT): jcapimin.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jcapistd.$(OBJEXT): jcapistd.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jccoefct.$(OBJEXT): jccoefct.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jccolor.$(OBJEXT): jccolor.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jcdctmgr.$(OBJEXT): jcdctmgr.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdct.h
jchuff.$(OBJECT): jchuff.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jchuff.h
jcinit.$(OBJEXT): jcinit.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jcmainct.$(OBJEXT): jcmainct.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jcmarker.$(OBJEXT): jcmarker.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jcmaster.$(OBJEXT): jcmaster.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jcomapi.$(OBJEXT): jcomapi.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jcparam.$(OBJEXT): jcparam.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jcphuff.$(OBJEXT): jcphuff.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jchuff.h
jcprepct.$(OBJEXT): jcprepct.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jcsample.$(OBJEXT): jcsample.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jctrans.$(OBJEXT): jctrans.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdapimin.$(OBJEXT): jdapimin.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdapistd.$(OBJEXT): jdapistd.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdatadst.$(OBJEXT): jdatadst.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h
jdatasrc.$(OBJEXT): jdatasrc.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h
jdcoefct.$(OBJEXT): jdcoefct.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdcolor.$(OBJEXT): jdcolor.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jddctmgr.$(OBJEXT): jddctmgr.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdct.h
jdhuff.$(OBJEXT): jdhuff.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdhuff.h
jdinput.$(OBJEXT): jdinput.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdmainct.$(OBJEXT): jdmainct.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdmarker.$(OBJEXT): jdmarker.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdmaster.$(OBJEXT): jdmaster.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdmerge.$(OBJEXT): jdmerge.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdphuff.$(OBJEXT): jdphuff.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdhuff.h
jdpostct.$(OBJEXT): jdpostct.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdsample.$(OBJEXT): jdsample.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jdtrans.$(OBJEXT): jdtrans.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jerror.$(OBJEXT): jerror.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jversion.h jerror.h
jfdctflt.$(OBJEXT): jfdctflt.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdct.h
jfdctfst.$(OBJEXT): jfdctfst.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdct.h
jfdctint.$(OBJEXT): jfdctint.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdct.h
jidctflt.$(OBJEXT): jidctflt.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdct.h
jidctfst.$(OBJEXT): jidctfst.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdct.h
jidctint.$(OBJEXT): jidctint.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdct.h
jidctred.$(OBJEXT): jidctred.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jdct.h
jquant1.$(OBJEXT): jquant1.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jquant2.$(OBJEXT): jquant2.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jutils.$(OBJEXT): jutils.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h
jmemmgr.$(OBJEXT): jmemmgr.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jmemsys.h
jmemansi.$(OBJEXT): jmemansi.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jmemsys.h
jmemname.$(OBJEXT): jmemname.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jmemsys.h
jmemnobs.$(OBJEXT): jmemnobs.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jmemsys.h
jmemdos.$(OBJEXT): jmemdos.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jmemsys.h
jmemmac.$(OBJEXT): jmemmac.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h jmemsys.h
cjpeg.$(OBJEXT): cjpeg.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h jversion.h
djpeg.$(OBJEXT): djpeg.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h jversion.h
jpegtran.$(OBJEXT): jpegtran.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h transupp.h jversion.h
rdjpgcom.$(OBJEXT): rdjpgcom.c jinclude.h jconfig.h
wrjpgcom.$(OBJEXT): wrjpgcom.c jinclude.h jconfig.h
cdjpeg.$(OBJEXT): cdjpeg.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
rdcolmap.$(OBJEXT): rdcolmap.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
rdswitch.$(OBJEXT): rdswitch.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
transupp.$(OBJEXT): transupp.c jinclude.h jconfig.h jpeglib.h jmorecfg.h jpegint.h jerror.h transupp.h
rdppm.$(OBJEXT): rdppm.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
wrppm.$(OBJEXT): wrppm.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
rdgif.$(OBJEXT): rdgif.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
wrgif.$(OBJEXT): wrgif.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
rdtarga.$(OBJEXT): rdtarga.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
wrtarga.$(OBJEXT): wrtarga.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
rdbmp.$(OBJEXT): rdbmp.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
wrbmp.$(OBJEXT): wrbmp.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
rdrle.$(OBJEXT): rdrle.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
wrrle.$(OBJEXT): wrrle.c cdjpeg.h jinclude.h jconfig.h jpeglib.h jmorecfg.h jerror.h cderror.h
