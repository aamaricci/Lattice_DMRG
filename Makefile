FC=mpif90
EXE=Heisenberg1d_iDMRG

SF_INC:=$(shell pkg-config --cflags scifor)
SF_LIB:=$(shell pkg-config --libs scifor)

OFLAG= -O2 -funroll-loops -ffree-line-length-none -fPIC -w -fallow-argument-mismatch
DFLAG= -cpp -D_DEBUG -O0 -p -g -Wsurprising -Waliasing -fwhole-file -fcheck=all -fbacktrace -fbounds-check  -ffree-line-length-none -fPIC -w -fallow-argument-mismatch
DDFLAG=-cpp -D_DEBUG -O0 -p -g  -fbacktrace -fwhole-file -fcheck=all -fbounds-check -fsanitize=address -fdebug-aux-vars -Wall -Waliasing -Wsurprising -Wampersand -Warray-bounds -Wc-binding-type -Wcharacter-truncation -Wconversion -Wdo-subscript -Wfunction-elimination -Wimplicit-interface -Wimplicit-procedure -Wintrinsic-shadow -Wintrinsics-std -Wno-align-commons -Wno-overwrite-recursive -Wno-tabs -Wreal-q-constant -Wunderflow -Wunused-parameter -Wrealloc-lhs -Wrealloc-lhs-all -Wfrontend-loop-interchange -Wtarget-lifetime -Wextra -Wimplicit-interface -Wno-unused-function -fPIC -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999999

OBJS = VERSION.o INPUT_VARS.o AUX_FUNCS.o HLOCAL.o MATRIX_SPARSE.o MATRIX_BLOCKS.o TUPLE_BASIS.o LIST_SECTORS.o LIST_OPERATORS.o  SITES.o BLOCKS.o SYSTEM.o HUBBARD_DMRG.o

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: git_code_version = "$(REV)"' > git_version.inc

all: FLAG=$(OFLAG)
all: code

debug: FLAG=$(DFLAG)
debug: code

ddebug: FLAG=$(DDFLAG)
ddebug: code

code: version $(OBJS)
	@echo "compiling test_iDMRG"
	${FC} $(FLAG) ${OBJS} drivers/$(EXE).f90 -o $(HOME)/.bin/$(EXE) ${SF_INC} ${SF_LIB}


test: hlocal matrix_sparse matrix_blocks tuples sectors operators sites blocks


hlocal: FLAG=$(DFLAG)
hlocal: MYOBJ=HLOCAL.o
hlocal: HLOCAL.o
	${FC} $(FLAG) $(MYOBJ) drivers/test/testHLOCAL.f90 -o $(HOME)/.bin/testHLOCAL ${SF_INC} ${SF_LIB}

matrix_sparse: FLAG=$(DFLAG)
matrix_sparse: MYOBJ=AUX_FUNCS.o MATRIX_SPARSE.o
matrix_sparse: AUX_FUNCS.o MATRIX_SPARSE.o
	${FC} $(FLAG) $(MYOBJ) drivers/test/testMATRIX_SPARSE.f90 -o $(HOME)/.bin/testMATRIX_SPARSE ${SF_INC} ${SF_LIB}


matrix_blocks: FLAG=$(DFLAG)
matrix_blocks: MYOBJ=AUX_FUNCS.o MATRIX_BLOCKS.o
matrix_blocks: AUX_FUNCS.o MATRIX_BLOCKS.o
	${FC} $(FLAG) $(MYOBJ) drivers/test/testMATRIX_BLOCKS.f90 -o $(HOME)/.bin/testMATRIX_BLOCKS ${SF_INC} ${SF_LIB}

tuples: FLAG=$(DFLAG)
tuples: MYOBJ=AUX_FUNCS.o TUPLE_BASIS.o
tuples: AUX_FUNCS.o TUPLE_BASIS.o
	${FC} $(FLAG) ${MYOBJ} drivers/test/testTUPLE_BASIS.f90 -o $(HOME)/.bin/testTUPLE_BASIS ${SF_INC} ${SF_LIB}

sectors: FLAG=$(DFLAG)
sectors: MYOBJ=AUX_FUNCS.o TUPLE_BASIS.o LIST_SECTORS.o
sectors: AUX_FUNCS.o TUPLE_BASIS.o LIST_SECTORS.o
	${FC} $(FLAG) ${MYOBJ} drivers/test/testLIST_SECTORS.f90 -o $(HOME)/.bin/testLIST_SECTORS ${SF_INC} ${SF_LIB}

operators: FLAG=$(DFLAG)
operators: MYOBJ=AUX_FUNCS.o MATRIX_SPARSE.o LIST_OPERATORS.o
operators: AUX_FUNCS.o MATRIX_SPARSE.o LIST_OPERATORS.o
	${FC} $(FLAG) ${MYOBJ} drivers/test/testLIST_OPERATORS.f90 -o $(HOME)/.bin/testLIST_OPERATORS ${SF_INC} ${SF_LIB}


sites: FLAG=$(DFLAG)
sites: MYOBJ=AUX_FUNCS.o MATRIX_SPARSE.o LIST_OPERATORS.o LIST_SECTORS.o HLOCAL.o SITES.o
sites: AUX_FUNCS.o MATRIX_SPARSE.o LIST_OPERATORS.o LIST_SECTORS.o HLOCAL.o SITES.o
	${FC} $(FLAG) ${MYOBJ} drivers/test/testSITES.f90 -o $(HOME)/.bin/testSITES ${SF_INC} ${SF_LIB}

blocks: FLAG=$(DFLAG)
blocks: $(OBJS)
	${FC} $(FLAG) ${OBJS} drivers/test/testBLOCKS.f90 -o $(HOME)/.bin/testBLOCKS ${SF_INC} ${SF_LIB}


clean: 
	@echo "Cleaning:"
	@rm -fv *.mod *.o *~ git_version.inc
	@echo ""

version:
	@echo "git_version: $(REV)"
	@echo $(VER)

.f90.o:	
	$(FC) $(FLAG) $(SF_INC) -c $< 




# qn:  $(OBJS)
# 	@echo "compiling testQN:"
# 	${FC} $(FLAG) ${OBJS} drivers/test/testQN.f90 -o $(HOME)/.bin/testQN ${SF_INC} ${SF_LIB}

# nsites:  $(OBJS)
# 	@echo "compiling EDnsites:"
# 	${FC} $(FLAG) ${OBJS} drivers/test/testEDnsites.f90 -o $(HOME)/.bin/testEDnsites ${SF_INC} ${SF_LIB}
