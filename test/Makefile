FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}


TSTDIR=src
TSTBIN=bin

all: bin
	make -C $(TSTDIR)

bin:
	mkdir $(TSTBIN)
	mkdir $(TSTBIN)/Spin1d_SU2
#	mkdir $(TSTBIN)/Spin1d_SU3
	mkdir $(TSTBIN)/Hubbard1d
####### INSERT FOLDER TO BE TESTED IN LIST_DIR
	echo Spin1d_SU2 >> $(TSTBIN)/list_dir
#       echo Spin1d_SU3  >> $(TSTBIN)/list_dir
	echo Hubbard1d >> $(TSTBIN)/list_dir


test:
	./test.sh

clean:
	make clean -C $(TSTDIR)
	rm -rf $(TSTBIN)
