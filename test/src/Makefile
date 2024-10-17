FC=mpif90
FFLAG = -O2 -ffree-line-length-none #-coverage
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

.SUFFIXES: .f90
.PHONY: clean

GLOB_INC:=$(shell pkg-config --cflags scifor )

OBJS=ASSERTING.o

all: $(OBJS) spin2 hubbard #spin3
	@echo "all done"


spin2:
	make -C Spin1d_SU2

spin3:
	make -C Spin1d_SU3

hubbard:
	make -C Hubbard1d

.f90.o:
	$(FC) $(FFLAGS) -c $< ${GLOB_INC}

clean:
	rm -rf *.mod *.o 

