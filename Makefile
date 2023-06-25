FC=mpif90


SF_INC:=$(shell pkg-config --cflags scifor)
SF_LIB:=$(shell pkg-config --libs scifor)

#FFLAG=-ffree-line-length-none
#FFLAG = -O2 -funroll-loops -ffree-line-length-none -fPIC -w -fallow-argument-mismatch
FFLAG = -O0 -p -g -Wsurprising -Waliasing -fwhole-file -fcheck=all -fbacktrace -fbounds-check  -ffree-line-length-none -fPIC -w -fallow-argument-mismatch
#FFLAG=-O0 -p -g  -fbacktrace -fwhole-file -fcheck=all -fbounds-check -fsanitize=address -fdebug-aux-vars -Wall -Waliasing -Wsurprising -Wampersand -Warray-bounds -Wc-binding-type -Wcharacter-truncation -Wconversion -Wdo-subscript -Wfunction-elimination -Wimplicit-interface -Wimplicit-procedure -Wintrinsic-shadow -Wintrinsics-std -Wno-align-commons -Wno-overwrite-recursive -Wno-tabs -Wreal-q-constant -Wunderflow -Wunused-parameter -Wrealloc-lhs -Wrealloc-lhs-all -Wfrontend-loop-interchange -Wtarget-lifetime -Wextra -Wimplicit-interface -Wno-unused-function -fPIC -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999999

OBJS     = AUX_FUNCS.o HLOCAL.o MATRIX_SPARSE.o MATRIX_BLOCKS.o  LIST_SECTORS.o LIST_OPERATORS.o  SITES.o BLOCKS.o
##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90


all: $(OBJS)

hlocal: MYOBJ=HLOCAL.o
hlocal: HLOCAL.o
	${FC} $(FFLAG) $(MYOBJ) ./example/testHLOCAL.f90 -o ./test/testHLOCAL ${SF_INC} ${SF_LIB}


matrix_sparse: MYOBJ=AUX_FUNCS.o MATRIX_SPARSE.o
matrix_sparse: AUX_FUNCS.o MATRIX_SPARSE.o
	${FC} $(FFLAG) $(MYOBJ) ./example/testMATRIX_SPARSE.f90 -o ./test/testMATRIX_SPARSE ${SF_INC} ${SF_LIB}



matrix_blocks: MYOBJ=AUX_FUNCS.o MATRIX_BLOCKS.o
matrix_blocks: AUX_FUNCS.o MATRIX_BLOCKS.o
	${FC} $(FFLAG) $(MYOBJ) ./example/testMATRIX_BLOCKS.f90 -o ./test/testMATRIX_BLOCKS ${SF_INC} ${SF_LIB}


sectors: MYOBJ=AUX_FUNCS.o LIST_SECTORS.o
sectors: AUX_FUNCS.o LIST_SECTORS.o
	${FC} $(FFLAG) ${MYOBJ} ./example/testLIST_SECTORS.f90 -o ./test/testLIST_SECTORS ${SF_INC} ${SF_LIB}


operators: MYOBJ=AUX_FUNCS.o MATRIX_SPARSE.o LIST_OPERATORS.o
operators: AUX_FUNCS.o MATRIX_SPARSE.o LIST_OPERATORS.o
	${FC} $(FFLAG) ${MYOBJ} ./example/testLIST_OPERATORS.f90 -o ./test/testLIST_OPERATORS ${SF_INC} ${SF_LIB}



sites: MYOBJ=AUX_FUNCS.o MATRIX_SPARSE.o LIST_OPERATORS.o LIST_SECTORS.o HLOCAL.o SITES.o
sites: AUX_FUNCS.o MATRIX_SPARSE.o LIST_OPERATORS.o LIST_SECTORS.o HLOCAL.o SITES.o
	${FC} $(FFLAG) ${MYOBJ} ./example/testSITES.f90 -o ./test/testSITES ${SF_INC} ${SF_LIB}

blocks: $(OBJS)
	${FC} $(FFLAG) ${OBJS} ./example/testBLOCKS.f90 -o ./test/testBLOCKS ${SF_INC} ${SF_LIB}

dmrg:  $(OBJS)
	@echo "compiling test_iDMRG"
	${FC} $(FFLAG) ${OBJS} ./test_iDMRG.f90 -o ./test/test_iDMRG ${SF_INC} ${SF_LIB}
	${FC} $(FFLAG) ${OBJS} ./test_iDMRGqn.f90 -o ./test/test_iDMRGqn ${SF_INC} ${SF_LIB}


qn:  $(OBJS)
	@echo "compiling testQN"
	${FC} $(FFLAG) ${OBJS} ./testQN.f90 -o ./test/testQN ${SF_INC} ${SF_LIB}

ed:  $(OBJS)
	@echo "compiling testEDkronUD"
	${FC} $(FFLAG) ${OBJS} example/testEDnsites.f90 -o ./test/testEDnsites ${SF_INC} ${SF_LIB}
	${FC} $(FFLAG) ${OBJS} ./testEDkron.f90 -o ./test/testEDkron ${SF_INC} ${SF_LIB}

clean: 
	@echo "Cleaning:"
	@rm -fv *.mod *.o *~ 
	@echo ""

.f90.o:	
	$(FC) $(FFLAG) $(SF_INC) -c $< 
