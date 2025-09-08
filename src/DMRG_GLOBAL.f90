MODULE DMRG_GLOBAL
  USE SCIFOR
  USE INPUT_VARS
  USE AUX_FUNCS
  USE MATRIX_GRAPH
  USE MATRIX_SPARSE, id=>sp_eye
  USE MATRIX_BLOCKS
  USE TUPLE_BASIS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none


  ! abstract interface
  !    function UserHconnect(left,right) 
  !      USE BLOCKS
  !      implicit none
  !      type(block)                          :: left
  !      type(block)                          :: right
  !      real(8),dimension(:,:),allocatable :: UserHconnect ![Nso,Nso]
  !    end function UserHconnect
  ! end interface
  ! procedure(UserHconnect),pointer,public  :: Hmodel=>null()
  ! !

#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable :: HopH
#else
  real(8),dimension(:,:),allocatable    :: HopH
#endif




  abstract interface
     subroutine sparse_HxV(Nloc,v,Hv)
       integer                          :: Nloc
#ifdef _CMPLX
       complex(8),dimension(Nloc)       :: v
       complex(8),dimension(Nloc)       :: Hv
#else
       real(8),dimension(Nloc)          :: v
       real(8),dimension(Nloc)          :: Hv
#endif
     end subroutine sparse_HxV
  end interface
  procedure(sparse_HxV),pointer,public  :: spHtimesV_p=>null()



  type(sparse_matrix)                   :: spHsb
  !
  real(8)                               :: truncation_error_left,truncation_error_right
  character(len=:),allocatable          :: suffix
  real(8),dimension(:),allocatable      :: target_Qn,current_target_QN
  type(block)                           :: init_left,init_right
  logical                               :: init_called=.false.
#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable :: gs_vector
#else
  real(8),dimension(:,:),allocatable    :: gs_vector
#endif
  real(8),dimension(:),allocatable      :: gs_energy
  real(8),dimension(:),allocatable      :: rho_left_evals
  real(8),dimension(:),allocatable      :: rho_right_evals
  type(blocks_matrix)                   :: psi_left
  type(blocks_matrix)                   :: psi_right
  type(blocks_matrix)                   :: rho_left
  type(blocks_matrix)                   :: rho_right
  !GLOBAL LEFT & RIGHT & DOT 
  type(block)                           :: left,right
  type(site),dimension(:),allocatable   :: dot
  !
  !SUPERBLOCK SHARED THINGS
  integer,dimension(:),allocatable      :: sb_states
  type(sectors_list)                    :: sb_sector
  !
  integer                               :: Mstates
  real(8)                               :: Estates




  !Memory pool for HxV direct product 
  type(sparse_matrix),allocatable,dimension(:)   :: Hleft,Hright
  type(sparse_matrix),allocatable,dimension(:,:) :: A,B
  integer,dimension(:),allocatable               :: Dls,Drs,Offset
  integer,dimension(:,:),allocatable             :: RowOffset,ColOffset,isb2jsb
  logical,dimension(:,:),allocatable             :: IsHconjg


  !This is the internal Mpi Communicator and variables.
  !=========================================================
#ifdef _MPI
  integer                               :: MpiComm_Global=MPI_COMM_NULL
  integer                               :: MpiComm=MPI_COMM_NULL
  integer                               :: MpiGroup_Global=MPI_GROUP_NULL
  integer                               :: MpiGroup=MPI_GROUP_NULL
#else
  integer                               :: MpiComm_Global=0
  integer                               :: MpiComm=0
  integer                               :: MpiGroup_Global=0
  integer                               :: MpiGroup=0
#endif
  logical                               :: MpiStatus=.false.
  logical                               :: MpiMaster=.true.
  integer                               :: MpiRank=0
  integer                               :: MpiSize=1
  integer,allocatable,dimension(:)      :: mpiDls
  integer,allocatable,dimension(:)      :: mpiDrs
  integer,allocatable,dimension(:)      :: mpiDl,mpiDr
  integer,allocatable,dimension(:)      :: mpiOffset
  integer,allocatable,dimension(:,:)    :: mpiRowOffset,mpiColOffset
  integer                               :: mpiL=0,mpiR=0
  !


#ifdef _MPI
  ! public :: scatter_vector_MPI
  ! public :: scatter_basis_MPI
  public :: gather_vector_MPI
  public :: allgather_vector_MPI
#endif





contains



  !##################################################################
  !##################################################################
  !                   MPI AUX FUNCTIONS
  !##################################################################
  !##################################################################

  subroutine dmrg_set_MpiComm()
#ifdef _MPI
    integer :: ierr
    MpiComm_Global = MPI_COMM_WORLD
    MpiComm        = MPI_COMM_WORLD
    call Mpi_Comm_group(MpiComm_Global,MpiGroup_Global,ierr)
    MpiStatus      = .true.
    MpiSize        = get_Size_MPI(MpiComm_Global)
    MpiRank        = get_Rank_MPI(MpiComm_Global)
    MpiMaster      = get_Master_MPI(MpiComm_Global)
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG dmrg_set_MpiComm: setting MPI comm"
#endif
#endif
  end subroutine dmrg_set_MpiComm

  subroutine dmrg_del_MpiComm()
#ifdef _MPI    
    MpiComm_Global = MPI_UNDEFINED
    MpiComm        = MPI_UNDEFINED
    MpiGroup_Global= MPI_GROUP_NULL
    MpiStatus      = .false.
    MpiSize        = 1
    MpiRank        = 0
    MpiMaster      = .true.
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG dmrg_del_MpiComm: deleting MPI comm"
#endif
#endif
  end subroutine dmrg_del_MpiComm



  !##################################################################
  !##################################################################
  !              MPI AUXILIARY FUNCTIONS
  !##################################################################
  !##################################################################
  !
  !####################################################################
  !               ALL-2-ALL-V VECTOR MPI TRANSPOSITION 
  !####################################################################
#ifdef _MPI
  subroutine vector_transpose_MPI(nrow,qcol,a,ncol,qrow,b)
    !
    integer                            :: nrow !Global number of rows 
    integer                            :: ncol !Global number of columns
    integer                            :: qrow !Local number of rows on each thread
    integer                            :: qcol !Local number of columns on each thread
    real(8)                            :: a(nrow,qcol) ! Input vector to be transposed
    real(8)                            :: b(ncol,qrow) ! Output vector :math:`b = v^T`
    real(8),dimension(:),allocatable   :: Vtmp
    integer,allocatable,dimension(:,:) :: send_counts,send_offset
    integer,allocatable,dimension(:,:) :: recv_counts,recv_offset
    integer                            :: counts,Ntot
    integer                            :: i,j,irank,ierr
    !
    counts = Nrow/MpiSize
    Ntot   = Ncol/MpiSize
    if(mod(Ncol,MpiSize)/=0)Ntot=Ntot+1
    !
    allocate(send_counts(0:MpiSize-1,Ntot));send_counts=0
    allocate(send_offset(0:MpiSize-1,Ntot));send_offset=0
    allocate(recv_counts(0:MpiSize-1,Ntot));recv_counts=0
    allocate(recv_offset(0:MpiSize-1,Ntot));recv_offset=0
    !
    do i=1,qcol
       do irank=0,MpiSize-1
          if(irank < mod(Nrow,MpiSize))then
             send_counts(irank,i) = counts+1
          else
             send_counts(irank,i) = counts
          endif
       enddo
    enddo
    !
    do i=1,Ntot
       call MPI_AllToAll(&
            send_counts(0:,i),1,MPI_INTEGER,&
            recv_counts(0:,i),1,MPI_INTEGER,&
            MpiComm,ierr)
    enddo
    !
    do i=1,Ntot
       do irank=1,MpiSize-1
          send_offset(irank,i) = send_counts(irank-1,i) + send_offset(irank-1,i)
       enddo
    enddo
    !
    !Get the irank=0 elements, i.e. first entries:
    recv_offset(0,1) = 0
    do i=2,Ntot
       recv_offset(0,i) = sum(recv_counts(0,:i-1))
    enddo
    !the rest of the entries:
    do i=1,Ntot
       do irank=1,MpiSize-1
          recv_offset(irank,i) = recv_offset(irank-1,i) + sum(recv_counts(irank-1,:))
       enddo
    enddo
    !
    !
    do j=1,Ntot
       !Fix issue with empty columns arising from having few MPI nodes
       if(j<=size(A,2))then
          Vtmp = A(:,j)            !automatic allocation
       else
          allocate(Vtmp(0))
       endif
       call MPI_AllToAllV(& ! A(:,j),send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
            Vtmp,send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
            B(:,:),recv_counts(:,j),recv_offset(:,j),MPI_DOUBLE_PRECISION,&
            MpiComm,ierr)
       deallocate(Vtmp)
    enddo
    !
    call local_transpose(b,ncol,qrow)
    !
    return
  end subroutine vector_transpose_MPI
  !
  subroutine local_transpose(mat,nrow,ncol)
    integer                      :: nrow,ncol
    real(8),dimension(Nrow,Ncol) :: mat
    mat = transpose(reshape(mat,[Ncol,Nrow]))
  end subroutine local_transpose
  !=========================================================








  !! Scatter V into the arrays Vloc on each thread: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine scatter_vector_MPI(MpiComm,v,vloc)
    integer                          :: MpiComm
#ifdef _CMPLX
    complex(8),dimension(:)          :: v    !size[N]
    complex(8),dimension(:)          :: vloc !size[Nloc]
#else
    real(8),dimension(:)             :: v    !size[N]
    real(8),dimension(:)             :: vloc !size[Nloc]
#endif
    integer                          :: i,k,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: pCounts,pOffset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG d_scatter_vector_MPI: scatter v into vloc"
#endif
    !
    if( MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null )return
    ! stop "scatter_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "scatter_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(pCounts(0:MpiSize-1)) ; pCounts=0
    allocate(pOffset(0:MpiSize-1)) ; pOffset=0
    !
    Vloc=0d0
    do k=1,size(sb_sector)
       !
       !Get Counts;
       pCounts=0
       pOffset=0
       call MPI_AllGather(mpiDl(k),1,MPI_INTEGER,pCounts,1,MPI_INTEGER,MpiComm,MpiIerr)
       !
       !Get Offset:
       pOffset(0)=0
       do i=1,MpiSize-1
          pOffset(i) = pOffset(i-1) + pCounts(i-1)
       enddo
       !
       if(MpiMaster)then
          v_start = 1 + Offset(k)
          v_end = Drs(k)*Dls(k)+OffSet(k)
       else
          v_start = 1
          v_end = 1
       endif
       !
       vloc_start = 1 + mpiOffset(k)
       vloc_end   = mpiDl(k)+mpiOffSet(k)
       !
#ifdef _CMPLX
       call MPI_Scatterv(V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_COMPLEX,&
            Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_COMPLEX,0,MpiComm,MpiIerr)
#else
       call MPI_Scatterv(V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_PRECISION,&
            Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_PRECISION,0,MpiComm,MpiIerr)
#endif
    enddo
    !
    return
  end subroutine scatter_vector_MPI



  subroutine scatter_basis_MPI(MpiComm,v,vloc)
    integer                   :: MpiComm
#ifdef _CMPLX
    complex(8),dimension(:,:) :: v    !size[N,N]
    complex(8),dimension(:,:) :: vloc !size[Nloc,Neigen]
#else
    real(8),dimension(:,:)    :: v    !size[N,N]
    real(8),dimension(:,:)    :: vloc !size[Nloc,Neigen]
#endif
    integer                   :: N,Nloc,Neigen,i
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG d_scatter_basis_MPI: scatter many v"
#endif
    N      = size(v,1)
    Nloc   = size(vloc,1)
    Neigen = size(vloc,2)
    if( size(v,2) < Neigen ) stop "error scatter_basis_MPI: size(v,2) < Neigen"
    !
    do i=1,Neigen
       call scatter_vector_MPI(MpiComm,v(:,i),vloc(:,i))
    end do
    !
    return
  end subroutine scatter_basis_MPI





  !! gather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine gather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
#ifdef _CMPLX
    complex(8),dimension(:)          :: vloc !size[Nloc]
    complex(8),dimension(:)          :: v    !size[N]
#else
    real(8),dimension(:)             :: vloc !size[Nloc]
    real(8),dimension(:)             :: v    !size[N]
#endif
    integer                          :: i,k,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: pCounts,pOffset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(verbose>4)write(Logfile,"(A)")"DEBUG d_gather_basis_MPI: gather  v"
#endif
    !
    if(  MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null ) return
    !stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "gather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(pCounts(0:MpiSize-1)) ; pCounts=0
    allocate(pOffset(0:MpiSize-1)) ; pOffset=0
    !
    V = 0d0
    do k=1,size(sb_sector)
       !
       !Get Counts;
       pCounts=0
       pOffset=0
       call MPI_AllGather(mpiDl(k),1,MPI_INTEGER,pCounts,1,MPI_INTEGER,MpiComm,MpiIerr)
       !
       !Get Offset:
       pOffset(0)=0
       do i=1,MpiSize-1
          pOffset(i) = pOffset(i-1) + pCounts(i-1)
       enddo
       !
       vloc_start = 1 + mpiOffset(k)
       vloc_end   = mpiDl(k)+mpiOffSet(k)
       !
       if(MpiMaster)then
          v_start = 1 + Offset(k)
          v_end = Drs(k)*Dls(k)+OffSet(k)
       else
          v_start = 1
          v_end = 1
       endif
       !
#ifdef _CMPLX
       call MPI_Gatherv(Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_COMPLEX,&
            V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_COMPLEX,0,MpiComm,MpiIerr)
#else
       call MPI_Gatherv(Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_PRECISION,&
            V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_PRECISION,0,MpiComm,MpiIerr)
#endif
    enddo
    !
    return
  end subroutine gather_vector_MPI




  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine allgather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
#ifdef _CMPLX
    complex(8),dimension(:)          :: vloc !size[Nloc]
    complex(8),dimension(:)          :: v    !size[N]
#else
    real(8),dimension(:)             :: vloc !size[Nloc]
    real(8),dimension(:)             :: v    !size[N]
#endif
    integer                          :: i,k,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: pCounts,pOffset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(verbose>4)write(Logfile,"(A)")"DEBUG d_allgather_basis_MPI: allgather v"
#endif
    !
    if(  MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null ) return
    ! stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N    = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "allgather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(pCounts(0:MpiSize-1)) ; pCounts=0
    allocate(pOffset(0:MpiSize-1)) ; pOffset=0
    !
    V = 0d0
    do k=1,size(sb_sector)
       !
       !Get Counts;
       pCounts=0
       pOffset=0
       call MPI_AllGather(mpiDl(k),1,MPI_INTEGER,pCounts,1,MPI_INTEGER,MpiComm,MpiIerr)
       !
       !Get Offset:
       pOffset(0)=0
       do i=1,MpiSize-1
          pOffset(i) = pOffset(i-1) + pCounts(i-1)
       enddo
       !
       vloc_start = 1 + mpiOffset(k)
       vloc_end   = mpiDl(k)+mpiOffSet(k)
       !
       v_start = 1 + Offset(k)
       v_end = Drs(k)*Dls(k)+OffSet(k)
       !
#ifdef _CMPLX
       call MPI_AllGatherv(Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_COMPLEX,&
            V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_COMPLEX,MpiComm,MpiIerr)
#else
       call MPI_AllGatherv(Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_PRECISION,&
            V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_PRECISION,MpiComm,MpiIerr)
#endif
    enddo
    !
    return
  end subroutine allgather_vector_MPI


#endif












  function sector_qn_dim(self,qn) result(dim)
    type(sectors_list)   :: self
    real(8),dimension(:) :: qn
    integer              :: dim
    dim = 0
    if(.not.self%has_qn(qn))return
    dim =  size(self%map(qn=qn))
  end function sector_qn_dim





  subroutine dmrg_graphic(label)
    integer                  :: label
    integer                  :: i,N,Mleft,Mright,LMleft,LMright,index,Ltot
    character(:),allocatable :: Ldot,Rdot
    real(8)                  :: eps=1d-6
    integer                  :: M=50
    !
    Ltot = Ldmrg
    Ldot = bold_green('=')
    Rdot = bold_red('-')
    ! if(Ltot>M)then
    !    Ldot = bg_green('=')
    !    Rdot = bg_red('-')
    ! endif
    !
    N = int(Ltot/(M+eps))+1
    !
    ! call execute_command_line("clear")
    do i=1,3
       write(LOGfile,*)""
    enddo
    select case(label)
    case default; stop "dmrg_graphic error: label != 1(L),2(R)"
    case(0)
       Mleft  = int(left%length/(N+eps))+1
       Mright = int(right%length/(N+eps))+1
       LMleft = Ltot/N-Mleft
       LMright= Ltot/N-Mright
       index=nint(mod(dble(left%length),N+eps))
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")&
            "left; right=",left%length,right%length,"|"
       if(LMleft>0)write(LOGfile,"("//str(LMleft)//"A)",advance="no")(" ",i=1,LMleft)
       write(LOGfile,"("//str(Mleft)//"A)",advance="no")(trim(Ldot),i=1,Mleft)
       write(LOGfile,"(A)",advance="no")bold_green("*")//bold("|")//bold_red("*")
       write(LOGfile,"("//str(Mright)//"A)",advance="no")(trim(Rdot),i=1,Mright)
       if(LMright>0)write(LOGfile,"("//str(LMright)//"A)",advance="no")(" ",i=1,LMright)
       if(Ltot<=M)then
          write(LOGfile,"(A1,2x,2I4)",advance='yes')"|",left%length+1,right%length+1
       else
          write(LOGfile,"(A1,2x,2I4,2x,I3,2x,A,1x,A)",advance='yes')"|",left%length+1,right%length+1, &
               index,Ldot//"->"//str(N)//"= ;",Rdot//"->"//str(N)//"-"
       endif
    case(1)
       Mleft  = int(left%length/(N+eps))+1
       Mright = int(right%length/(N+eps))+1
       LMleft = 0
       LMright= 0
       index=nint(mod(dble(left%length),N+eps))
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")&
            "left; right=",left%length,right%length,"|"
       if(LMleft>0)write(LOGfile,"("//str(LMleft)//"A)",advance="no")(" ",i=1,LMleft)
       write(LOGfile,"("//str(Mleft)//"A)",advance="no")(trim(Ldot),i=1,Mleft)
       write(LOGfile,"(A)",advance="no")bg_green(">")//"|"//bold_red("*")
       write(LOGfile,"("//str(Mright)//"A)",advance="no")(trim(Rdot),i=1,Mright)
       if(LMright>0)write(LOGfile,"("//str(LMright)//"A)",advance="no")(" ",i=1,LMright)
       if(Ltot<=M)then
          write(LOGfile,"(A1,2x,2I4)",advance='yes')"|",left%length+1,right%length+1
       else
          write(LOGfile,"(A1,2x,2I4,2x,I3,2x,A,1x,A)",advance='yes')"|",left%length+1,right%length+1, &
               index,Ldot//"->"//str(N)//"= ;",Rdot//"->"//str(N)//"-"
       endif
    case(2)
       Mleft  = int(left%length/(N+eps))+1
       Mright = int(right%length/(N+eps))+1
       LMleft = 0
       LMright= 0
       index=nint(mod(dble(left%length),N+eps))
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")&
            "left; right=",left%length,right%length,"|"
       if(LMleft>0)write(LOGfile,"("//str(LMleft)//"A)",advance="no")(" ",i=1,LMleft)
       write(LOGfile,"("//str(Mleft)//"A)",advance="no")(trim(Ldot),i=1,Mleft)
       write(LOGfile,"(A)",advance="no")bold_green("*")//"|"//bg_red("<")
       write(LOGfile,"("//str(Mright)//"A)",advance="no")(trim(Rdot),i=1,Mright)
       if(LMright>0)write(LOGfile,"("//str(LMright)//"A)",advance="no")(" ",i=1,LMright)
       if(Ltot<=M)then
          write(LOGfile,"(A1,2x,2I4)",advance='yes')"|",left%length+1,right%length+1
       else
          write(LOGfile,"(A1,2x,2I4,2x,I3,2x,A,1x,A)",advance='yes')"|",left%length+1,right%length+1, &
               index,Ldot//"->"//str(N)//"= ;",Rdot//"->"//str(N)//"-"
       endif
    end select
    call wait(250)
  end subroutine dmrg_graphic







END MODULE DMRG_GLOBAL





