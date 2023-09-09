# Set the list of programs to compile

PROGRAMS = amvw.a

# Compile all of the test programs and the library
all	      	             : amvw.a

AMVW_FILES                   = balance.o chasebulge.o crgivens.o diagblock.o  factor.o \
                              init_random_seed.o normalpoly.o rescheck.o zamvw2.o      \
                              buildbulge.o cnormalpoly.o deflation.o dto4.o            \
                              fuse.o modified_quadratic.o  randpolyjt.o throughdiag.o  \
                              zamvw.o


amvw.a                      : $(AMVW_FILES)
	echo $(FCOMP)
	ar cr amvw.a $(AMVW_FILES) 

%.o		: %.f90
	$(FCOMP) -c $(FOPTS)  $<


.PHONY: clean

clean:
	rm -f *.o *.a


