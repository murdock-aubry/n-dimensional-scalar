.phony:  all PROGRAM clean

# Choose the compiler and set compiler options
OMP_NUM_THREADS      := 1024#$(shell nproc)
export OMP_NUM_THREADS

#OMP_STACKSIZE        := 1024M
#export OMP_STACKSIZE

FCOMP     =  gfortran
FOPTS     =  -fopenmp -O3 -fexternal-blas  -march=native -w -fallow-argument-mismatch $(PREC)
#FOPTS    +=  -freal-8-real-10 -fdefault-integer-8
#FOPTS    +=  -fdefault-real-8 -fdefault-integer-8

# specify BLAS and LAPACK library
LDOPT     =  -lblas -llapack -flto
#LDOPT     = OpenBLAS/libopenblas.a
#LDOPT     =  FLAME/libflame.a BLIS/libblis.a

FC        = $(FCOMP)
FFLAGS    = $(FOPTS)

CC        = gcc
COPTS     = -O3  -I./include

export FCOMP FOPTS

# Set the list of programs to compile

PROGRAMS = experiment1 experiment2 experiment5 experiment_n

# Compile all of the test programs and the library
all	      	             : clean $(PROGRAMS) 


# List the dependencies for each module's test program


EXPERIMENT_N_FILES             = utils.o                                                  \
                                 linalg0.o                                                \
                                 legendre.o                                               \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 odesolve_n.o                                             \
                                 odesolve.o                                               \
                                 odetwo.o                                                 \
                                 scalar_n.o                                               \
                                 scalar.o                                                 \
                                 AMVW/amvw.a


EXPERIMENT5_FILES              = utils.o                                                  \
                                 linalg0.o                                                \
                                 legendre.o                                               \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 odesolve.o                                               \
                                 odetwo.o                                                 \
                                 scalar.o                                                 \
                                 AMVW/amvw.a


EXPERIMENT2_FILES              = utils.o                                                  \
                                 linalg0.o                                                \
                                 legendre.o                                               \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 odesolve.o                                               \
                                 odetwo.o                                                 \
                                 scalar.o                                                 \
                                 AMVW/amvw.a

EXPERIMENT1_FILES              = utils.o                                                  \
                                 linalg0.o                                                \
                                 legendre.o                                               \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 odesolve.o                                               \
                                 odetwo.o                                                 \
                                 scalar.o                                                 \
                                 AMVW/amvw.a


AMVW/amvw.a                    : 
	$(MAKE) -C AMVW

experiment_n.o                  : $(EXPERIMENT_N_FILES) experiment_n.f90
experiment_n                    : $(EXPERIMENT_N_FILES) experiment_n.o

experiment5.o                  : $(EXPERIMENT5_FILES) experiment5.f90
experiment5                    : $(EXPERIMENT5_FILES) experiment5.o

experiment2.o                  : $(EXPERIMENT2_FILES) experiment2.f90
experiment2                    : $(EXPERIMENT2_FILES) experiment2.o

experiment1.o                  : $(EXPERIMENT1_FILES) experiment1.f90
experiment1                    : $(EXPERIMENT1_FILES) experiment1.o

# Setup the general compilation rules

%		: %.o
	$(FCOMP) $(FOPTS) -o $@ $^ $(LDOPT)
	@echo  
	@echo 
	@echo "---------[ $@     ]--------------------------------------------"
	@echo 
	@./$@
	@echo 
	@echo "--------------------------------------------------------------------------"
	@echo 


%.o		: %.f90
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.f
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.cpp
	$(CPPCOMP) -c $(CPPOPTS)  $<

%.o		: %.c
	$(CC) -c $(COPTS)  $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~ fort.* a.out *.a
	rm -f $(PROGRAMS)
	rm -f *.dat
	rm -f *.py
	rm -f *.pdf
	rm -f *.tex
	cd AMVW; make clean; cd ..



