LIBDIR = ../../../../lib
AR_FLAGS = ruv
F77    = gfortran -c
LFLAGS = 
FFLAGS = -g -O2 -fno-second-underscore -fallow-argument-mismatch
#


OBJS1 = \
	 cfftb.o cfftb1.o cfftf.o cfftf1.o cffti.o cffti1.o cosqb.o cosqb1.o\
	 cosqf.o cosqf1.o cosqi.o cost.o costi.o ezfft1.o ezfftb.o ezfftf.o\
	 ezffti.o passb.o passb2.o passb3.o passb4.o passb5.o passf.o\
	 passf2.o passf3.o passf4.o passf5.o pimach.o radb2.o radb3.o\
	 radb4.o radb5.o radbg.o radf2.o radf3.o radf4.o radf5.o radfg.o\
	 rfftb.o rfftb1.o rfftf.o rfftf1.o rffti.o rffti1.o sinqb.o sinqf.o\
	 sinqi.o sint.o sint1.o sinti.o tfftpk.o

FSRC1 = \
	 cfftb.f cfftb1.f cfftf.f cfftf1.f cffti.f cffti1.f cosqb.f cosqb1.f\
	 cosqf.f cosqf1.f cosqi.f cost.f costi.f ezfft1.f ezfftb.f ezfftf.f\
	 ezffti.f passb.f passb2.f passb3.f passb4.f passb5.f passf.f\
	 passf2.f passf3.f passf4.f passf5.f pimach.f radb2.f radb3.f\
	 radb4.f radb5.f radbg.f radf2.f radf3.f radf4.f radf5.f radfg.f\
	 rfftb.f rfftb1.f rfftf.f rfftf1.f rffti.f rffti1.f sinqb.f sinqf.f\
	 sinqi.f sint.f sint1.f sinti.f tfftpk.f

OTHR_SOURCES = Makefile README

TARGET	= $(LIBDIR)/libfftpack.a

.SUFFIXES:      .o .f

.f.o:
	$(F77) $(FFLAGS) -o $@ $<

all: $(TARGET)

$(TARGET):  $(OBJS1)
	${AR} $(AR_FLAGS) $(TARGET)  $(OBJS1)
	chmod 644 $(TARGET)

clean:
	- rm -f $(OBJS1) 

distclean: clean
	-rm -f ${TARGET} Makefile

realclean: distclean
	-rm ${FSRC1} Makefile.in README

uninstall:
	-rm -f ${TARGET}

