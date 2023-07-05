MODULE INPUT_VARS
  USE SF_VERSION
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str
  ! USE VERSION
  implicit none

  !input variables
  !=========================================================
  integer              :: Ldmrg               !# of iDMRG steps to take, Ldmrg=max length of the SB.
  integer              :: Mdmrg               !# of states to retain at truncation. If 0 use Edmrg as threshold.  
  real(8)              :: Edmrg               !Threshold energy used to evaluate the number of states to keep. If 0d0 use fixed Mdmrg.
  integer              :: Lsweep              !# of DMRG sweep to take for finite DMRG algorithm.
  !
  integer              :: Norb                !# of orbitals
  integer              :: Nspin               !Nspin=# spin degeneracy (max 2)
  !
  real(8),allocatable  :: Uloc(:)             !local interactions
  real(8)              :: Ust                 !intra-orbitals interactions
  real(8)              :: Jh                  !J_Hund: Hunds' coupling constant 
  real(8)              :: Jx                  !J_X: coupling constant for the spin-eXchange interaction term
  real(8)              :: Jp                  !J_P: coupling constant for the Pair-hopping interaction term
  !
  real(8)              :: xmu                 !chemical potential
  real(8)              :: temp                !temperature
  !
  real(8)              :: eps                 !broadening
  real(8)              :: wini,wfin           !frequency range
  logical              :: HFmode              !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)              :: cutoff              !cutoff for spectral summation
  real(8)              :: gs_threshold        !Energy threshold for ground state degeneracy loop up
  !

  logical,allocatable  :: gf_flag(:)           !evaluate Green's functions for Norb+1
  logical,allocatable  :: chispin_flag(:)      !evaluate spin susceptibility for Norb+1
  !
  !
  integer              :: filling          !Total number of allowed electrons
  logical              :: sparse_H         !flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--
  integer              :: verbose          !
  !
  character(len=12)    :: lanc_method         !select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only) 
  real(8)              :: lanc_tolerance      !Tolerance for the Lanczos iterations as used in Arpack and plain lanczos. 
  integer              :: lanc_niter          !Max number of Lanczos iterations
  integer              :: lanc_ngfiter        !Max number of iteration in resolvant tri-diagonalization
  integer              :: lanc_ncv_factor     !Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_ncv_add        !Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_neigen         !Number of required eigenvalues per sector in the SuperBlock
  integer              :: lanc_dim_threshold  !Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.


  !Some parameters for function dimension:
  !=========================================================
  integer              :: Lmats
  integer              :: Lreal
  integer              :: Ltau

  !LOG
  !=========================================================
  integer,save         :: LOGfile


  !THIS IS JUST A RELOCATED GLOBAL VARIABLE
  character(len=200)   :: input_file=""


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine read_input(INPUTunit)
#ifdef _MPI
    USE MPI
    USE SF_MPI
#endif
    character(len=*) :: INPUTunit
    logical          :: master=.true.
    integer          :: i,rank=0,add,dim
#ifdef _MPI
    if(check_MPI())then
       master=get_Master_MPI(MPI_COMM_WORLD)
       rank  =get_Rank_MPI(MPI_COMM_WORLD)
    endif
#endif
    !
    !Store the name of the input file:
    input_file=str(INPUTunit)
    !
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Ldmrg,"Ldmrg",INPUTunit,default=5,comment="iDMRG steps to take=max length of the SB.")
    call parse_input_variable(Mdmrg,"Mdmrg",INPUTunit,default=20,comment="# of states to retain at truncation. If 0 use Edmrg as threshold.  ")
    call parse_input_variable(Edmrg,"Rdmrg",INPUTunit,default=0d0,comment="Threshold energy used to evaluate the number of states to keep. If 0d0 use fixed Mdmrg.")
    call parse_input_variable(Lsweep,"Lsweep",INPUTunit,default=1,comment="# of DMRG sweep to take for finite DMRG algorithm.")


    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of impurity orbitals.")
    ! call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy (max 2)")
    ! call parse_input_variable(filling,"FILLING",INPUTunit,default=0,comment="Total number of allowed electrons")
    allocate(Uloc(Norb))
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=(/( 2d0,i=1,size(Uloc) )/),comment="Values of the local interaction per orbital")
    ! call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    ! call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    ! call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="S-E coupling")
    ! call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="P-H coupling")
    !
    ! call parse_input_variable(temp,"TEMP",INPUTunit,default=0d0,comment="temperature")
    !
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")
    !
    ! call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
    ! call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
    !
    ! call parse_input_variable(Lmats,"LMATS",INPUTunit,default=4096,comment="Number of Matsubara frequencies.")
    ! call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of real-axis frequencies.")
    ! call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1024,comment="Number of imaginary time points.")
    !
    ! allocate(gf_flag(Norb))
    ! allocate(chispin_flag(Norb))
    ! call parse_input_variable(gf_flag,"GF_FLAG",INPUTunit,default=(/( .false.,i=1,size(gf_flag) )/),comment="Flag to activate Greens functions calculation")
    ! call parse_input_variable(chispin_flag,"CHISPIN_FLAG",INPUTunit,default=(/( .false.,i=1,size(chispin_flag) )/),comment="Flag to activate spin susceptibility calculation.")
    ! !
    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.,comment="Flag to set the Hartree form of the interaction (n-1/2). see xmu.")
    ! call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
    ! call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
    ! call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")
    !
    ! call parse_input_variable(sparse_H,"SPARSE_H",INPUTunit,default=.true.,comment="flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE ")   
    ! call parse_input_variable(verbose,"VERBOSE",INPUTunit,default=3,comment="Verbosity level: 0=almost nothing --> 5:all. Really: all")
    !
    call parse_input_variable(lanc_method,"LANC_METHOD",INPUTunit,default="arpack",comment="select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only)")
    call parse_input_variable(lanc_neigen,"LANC_NEIGEN",INPUTunit,default=2,comment="Number of states per SB sector to be determined.")
    call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",INPUTunit,default=10,comment="Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",INPUTunit,default=0,comment="Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_niter,"LANC_NITER",INPUTunit,default=512,comment="Number of Lanczos iteration in spectrum determination.")
    call parse_input_variable(lanc_ngfiter,"LANC_NGFITER",INPUTunit,default=200,comment="Number of Lanczos iteration in GF determination. Number of momenta.")
    call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",INPUTunit,default=1d-18,comment="Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.")
    call parse_input_variable(lanc_dim_threshold,"LANC_DIM_THRESHOLD",INPUTunit,default=1024,comment="Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.")
    !
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")


#ifdef _MPI
    if(check_MPI())then
       if(.not.master)then
          LOGfile=1000-rank
          open(LOGfile,file="stdOUT.rank"//str(rank)//".ed")
          do i=1,get_Size_MPI(MPI_COMM_WORLD)
             if(i==rank)write(*,"(A,I0,A,I0)")"Rank ",rank," writing to unit: ",LOGfile
          enddo
       endif
    endif
#endif
    !
    !
    if(master)then
       call print_input()
       call save_input(INPUTunit)
       call scifor_version()
       ! call code_version(version)
    endif
  end subroutine read_input




  subroutine substring_delete (s,sub)
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop


END MODULE INPUT_VARS
