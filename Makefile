FC=mpif90
EXE=dmrg_hm_1d
#EXE=kron_hm_1d

SF_INC:=$(shell pkg-config --cflags scifor)
SF_LIB:=$(shell pkg-config --libs scifor)

NFLAG= -cpp -D_ -ffree-line-length-none -fPIC -w -fallow-argument-mismatch
OFLAG= -cpp -D_ -O3 -funroll-loops -ffree-line-length-none -fPIC -w -fallow-argument-mismatch -fopenmp -Ofast -ftree-parallelize-loops=4
DFLAG= -cpp -D_DEBUG -O0 -p -g -Wsurprising -Waliasing -fwhole-file -fcheck=all -fbacktrace -fbounds-check  -ffree-line-length-none -fPIC -w -fallow-argument-mismatch -fopenmp
DDFLAG=-cpp -D_DEBUG -O0 -p -g  -fbacktrace -fwhole-file -fcheck=all -fbounds-check  -fdebug-aux-vars -Wall -Waliasing -Wsurprising -Wampersand -Warray-bounds -Wc-binding-type -Wcharacter-truncation -Wconversion -Wdo-subscript -Wfunction-elimination -Wimplicit-interface -Wimplicit-procedure -Wintrinsic-shadow -Wintrinsics-std -Wno-align-commons -Wno-overwrite-recursive -Wno-tabs -Wreal-q-constant -Wunderflow -Wunused-parameter -Wrealloc-lhs -Wrealloc-lhs-all -Wfrontend-loop-interchange -Wtarget-lifetime -Wextra -Wimplicit-interface -Wno-unused-function -fPIC -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999999 
#-fsanitize=address
OBJS = VERSION.o INPUT_VARS.o AUX_FUNCS.o HLOCAL.o MATRIX_SPARSE.o  TUPLE_BASIS.o LIST_SECTORS.o MATRIX_BLOCKS.o LIST_OPERATORS.o  SITES.o BLOCKS.o VARS_GLOBAL.o DMRG_CONNECT.o DMRG_SUPERBLOCK_COMMON.o DMRG_SUPERBLOCK_SPARSE.o DMRG_SUPERBLOCK_DIRECT.o DMRG_SUPERBLOCK.o DMRG_RDM.o DMRG_MEASURE.o DMRG_SYSTEM.o  DMRG.o

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: git_code_version = "$(REV)"' > git_version.inc

all: FLAG=$(NFLAG)
all: code

opt: FLAG=$(OFLAG)
opt: code

debug: FLAG=$(DFLAG)
debug: code

ddebug: FLAG=$(DDFLAG)
ddebug: code

code: version $(OBJS)
	@echo "compiling ${EXE}"
	${FC} $(FLAG) ${OBJS} drivers/$(EXE).f90 -o $(HOME)/.bin/$(EXE) ${SF_INC} ${SF_LIB}
	@echo " "
	@echo "> created: ${EXE}"	

test: hlocal matrix_sparse matrix_blocks tuples sectors operators sites blocks


hlocal: FLAG=$(DFLAG)
hlocal: MYOBJ=VERSION.o INPUT_VARS.o AUX_FUNCS.o
hlocal: version VERSION.o INPUT_VARS.o AUX_FUNCS.o
	${FC} $(FLAG) -D_TEST $(MYOBJ) HLOCAL.f90 -o $(HOME)/.bin/testHLOCAL ${SF_INC} ${SF_LIB}

matrix_sparse: FLAG=$(DFLAG)
matrix_sparse: MYOBJ=VERSION.o INPUT_VARS.o AUX_FUNCS.o 
matrix_sparse: version VERSION.o INPUT_VARS.o AUX_FUNCS.o
	${FC} $(FLAG) -D_TEST $(MYOBJ) MATRIX_SPARSE.f90 -o $(HOME)/.bin/testMATRIX_SPARSE ${SF_INC} ${SF_LIB}

tuples: FLAG=$(DFLAG)
tuples: MYOBJ=VERSION.o INPUT_VARS.o AUX_FUNCS.o
tuples: version VERSION.o INPUT_VARS.o AUX_FUNCS.o
	${FC} $(FLAG) -D_TEST ${MYOBJ}  TUPLE_BASIS.f90 -o $(HOME)/.bin/testTUPLE_BASIS ${SF_INC} ${SF_LIB}


sectors: FLAG=$(DFLAG)
sectors: MYOBJ=VERSION.o INPUT_VARS.o AUX_FUNCS.o TUPLE_BASIS.o 
sectors: version VERSION.o INPUT_VARS.o AUX_FUNCS.o TUPLE_BASIS.o
	${FC} $(FLAG) -D_TEST ${MYOBJ} LIST_SECTORS.f90 -o $(HOME)/.bin/testLIST_SECTORS ${SF_INC} ${SF_LIB}

matrix_blocks: FLAG=$(DFLAG)
matrix_blocks: MYOBJ=VERSION.o INPUT_VARS.o AUX_FUNCS.o TUPLE_BASIS.o LIST_SECTORS.o MATRIX_SPARSE.o
matrix_blocks: version VERSION.o INPUT_VARS.o AUX_FUNCS.o TUPLE_BASIS.o  LIST_SECTORS.o MATRIX_SPARSE.o 
	${FC} $(FLAG) -D_TEST $(MYOBJ) MATRIX_BLOCKS.f90 -o $(HOME)/.bin/testMATRIX_BLOCKS ${SF_INC} ${SF_LIB}


operators: FLAG=$(DFLAG)
operators: MYOBJ= VERSION.o INPUT_VARS.o AUX_FUNCS.o MATRIX_SPARSE.o
operators: version VERSION.o INPUT_VARS.o AUX_FUNCS.o MATRIX_SPARSE.o
	${FC} $(FLAG) -D_TEST ${MYOBJ} LIST_OPERATORS.f90 -o $(HOME)/.bin/testLIST_OPERATORS ${SF_INC} ${SF_LIB}


sites: FLAG=$(DFLAG)
sites: MYOBJ=VERSION.o INPUT_VARS.o HLOCAL.o AUX_FUNCS.o MATRIX_SPARSE.o TUPLE_BASIS.o LIST_OPERATORS.o LIST_SECTORS.o 
sites: version VERSION.o INPUT_VARS.o HLOCAL.o AUX_FUNCS.o MATRIX_SPARSE.o TUPLE_BASIS.o LIST_OPERATORS.o LIST_SECTORS.o
	${FC} $(FLAG) -D_TEST ${MYOBJ} SITES.f90 -o $(HOME)/.bin/testSITES ${SF_INC} ${SF_LIB}

blocks: FLAG=$(DFLAG)	
blocks: MYOBJ=VERSION.o INPUT_VARS.o AUX_FUNCS.o HLOCAL.o MATRIX_SPARSE.o MATRIX_BLOCKS.o TUPLE_BASIS.o LIST_SECTORS.o LIST_OPERATORS.o  SITES.o
blocks: version VERSION.o INPUT_VARS.o AUX_FUNCS.o HLOCAL.o MATRIX_SPARSE.o MATRIX_BLOCKS.o TUPLE_BASIS.o LIST_SECTORS.o LIST_OPERATORS.o  SITES.o
	${FC} $(FLAG) -D_TEST ${MYOBJ} BLOCKS.f90 -o $(HOME)/.bin/testBLOCKS ${SF_INC} ${SF_LIB}


clean: 
	@echo "Cleaning:"
	@rm -fv *.mod *.o *~ git_version.inc
	@echo ""

version:
	@echo "git_version: $(REV)"
	@echo $(VER)

.f90.o:	
	$(FC) $(FLAG) $(SF_INC) -c $< 
