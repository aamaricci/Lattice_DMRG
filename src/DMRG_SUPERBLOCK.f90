MODULE DMRG_SUPERBLOCK
  USE DMRG_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK_SETUP
  implicit none
  private


  public :: sb_get_states
  public :: sb_diag

  public :: sb_build_Hv
  public :: sb_vecDim_Hv
  public :: sb_delete_Hv
  
  integer :: i,j
  integer :: ispin
  integer :: iorb,jorb
  integer :: io,jo


contains



  !-----------------------------------------------------------------!
  ! Purpose: build the list of states compatible with the specified
  ! quantum numbers
  !-----------------------------------------------------------------!
  subroutine sb_get_states()
    integer                          :: ileft,iright
    integer                          :: i,j,istate,unit
    real(8),dimension(:),allocatable :: left_qn,right_qn
    integer,dimension(:),allocatable :: left_map,right_map
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock get states"
#endif
    !
    if(MpiMaster)call start_timer("Get SB states")
    !
    if(allocated(sb_states))deallocate(sb_states)
    !
    call sb_sector%free()
    !
#ifdef _DEBUG
    if(MpiMaster.AND.verbose>5)unit = fopen('SB_list_'//str(left%length)//'.dat')
#endif
    do ileft=1,size(left%sectors(1))
       left_qn   = left%sectors(1)%qn(index=ileft)
       right_qn  = current_target_qn - left_qn
       if(.not.right%sectors(1)%has_qn(right_qn))cycle
       !
       left_map  = left%sectors(1)%map(qn=left_qn)
       right_map = right%sectors(1)%map(qn=right_qn)
       !
#ifdef _DEBUG
       if(MpiMaster.AND.verbose>5)then
          write(unit,*)""
          write(unit,*)left_qn,right_qn
          write(unit,*)size(left_map),size(right_map)
       endif
#endif
       do i=1,size(left_map)
          do j=1,size(right_map)
             istate=right_map(j) + (left_map(i)-1)*right%Dim
             call append(sb_states, istate)
             call sb_sector%append(qn=left_qn,istate=size(sb_states))
#ifdef _DEBUG
             if(MpiMaster.AND.verbose>5)write(unit,*)left_map(i),right_map(j),istate
#endif
          enddo
       enddo
    enddo
#ifdef _DEBUG
    if(MpiMaster.AND.verbose>5)close(unit)
#endif
    !
    call stop_timer()
    !
  end subroutine sb_get_states




  !-----------------------------------------------------------------!
  ! Purpose: Diagonalize the SuperBlock problem.
  !-----------------------------------------------------------------!
  subroutine sb_diag()
    integer                               :: m_sb
    integer                               :: Nitermax,Neigen,Nblock
    real(8),dimension(:),allocatable      :: evals
#ifdef _MPI
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: mpiEvec
#else
    real(8),dimension(:,:),allocatable    :: mpiEvec
#endif
#endif       
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hsb
#else
    real(8),dimension(:,:),allocatable    :: Hsb
#endif
    integer                               :: vecDim,Nloc,m_tmp
    logical                               :: exist,lanc_solve
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock diagonalization"
#endif
    !
    m_sb = size(sb_states)
    if(m_sb==0)stop "sb_diag ERROR: size(sb_states)==0"
    !
    !Set Lanczos params
    Neigen   = min(m_sb,lanc_neigen)
    Nitermax = min(m_sb,lanc_niter)
    Nblock   = min(m_sb,lanc_ncv_factor*Lanc_Neigen+lanc_ncv_add)
    !
    !Decide how to operate on H_sb
    lanc_solve  = .true.
    if(Neigen==m_sb)lanc_solve=.false.
    if(m_sb <= lanc_dim_threshold)lanc_solve=.false.
    !
    !Allocate EigPairs
    if(allocated(gs_energy))deallocate(gs_energy)
    allocate(gs_energy(Neigen));gs_energy=0d0
    !
    if(MpiMaster)call start_timer("Diag H_sb")
    if(lanc_solve)then          !Use (P)-Arpack
       !
       call sb_build_Hv()
       !
       if(allocated(gs_vector))deallocate(gs_vector)
       vecDim = sb_vecDim_Hv()
       !
#ifdef _MPI
       if(MpiStatus)then          
          allocate(mpiEvec(vecDim,Neigen));mpiEvec=zero
          call sp_eigh(MpiComm,spHtimesV_p,gs_energy,mpiEvec,&
               Nblock,&
               Nitermax,&
               tol=lanc_tolerance,&
               iverbose=(verbose>4))

          m_tmp= 0
          call AllReduce_MPI(MpiComm,vecDim,m_tmp)
          if(m_tmp/=m_sb)stop "DMRG.sb_diag ERROR: reduced vecDim != m_sb"
          allocate(gs_vector(m_sb,Neigen));gs_vector=zero
          do i=1,Neigen
             call allgather_vector_MPI(MpiComm,mpiEvec(:,i),gs_vector(:,i))
          enddo
          deallocate(mpiEvec)
       else
          allocate(gs_vector(vecDim,Neigen));gs_vector=zero
          call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
               Nblock,&
               Nitermax,&
               tol=lanc_tolerance,&
               iverbose=(verbose>4))
       endif
#else
       allocate(gs_vector(vecDim,Neigen));gs_vector=zero
       call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
            Nblock,&
            Nitermax,&
            tol=lanc_tolerance,&
            iverbose=(verbose>4))
#endif
       !
    else !use LAPACK
       !
       call sb_build_Hv(Hsb)
       !
       if(allocated(gs_vector))deallocate(gs_vector)
       allocate(gs_vector(m_sb,Neigen));gs_vector=zero
       allocate(evals(m_sb))
       !
       if(MpiMaster)call eigh(Hsb,evals)
#ifdef _MPI
       if(MpiStatus)then
          call Bcast_MPI(MpiComm,Hsb)
          call Bcast_MPI(MpiComm,evals)
       endif
#endif
       gs_vector(:,1:Neigen) = Hsb(:,1:Neigen)
       gs_energy(1:Neigen)   = evals(1:Neigen)
       !
       deallocate(Hsb,evals)
       !
    endif
    call stop_timer()
    !    
    !Free Memory
    call sb_delete_Hv()
  end subroutine sb_diag









  !##################################################################
  !         SETUP THE SUPERBLOCK HAMILTONIAN PROBLEM
  ! . if Hmat: returb H^SB as dense matrix there for Lapack use
  ! . if sparse_H = T: build H^SB as sparse matrix
  ! . if sparse_H = F: setup H^SB terms and blocks for H*v procedure
  !##################################################################
  subroutine sb_build_Hv(Hmat)
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable,optional :: Hmat
#else
    real(8),dimension(:,:),allocatable,optional    :: Hmat
#endif
    real(8),dimension(:),allocatable               :: qn
    integer                                        :: q,m_sb,Nsb,irank
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock build H*v"
#endif
    !
    if(.not.allocated(sb_states))stop "build_Hv_superblock ERROR: sb_states not allocated"
    m_sb = size(sb_states)
    Nsb  = size(sb_sector)
    !
#ifdef _MPI
#ifdef _DEBUG
    if(allocated(Dls))deallocate(Dls)
    if(allocated(Drs))deallocate(Drs)
    if(allocated(mpiDls))deallocate(mpiDls)
    allocate(Dls(Nsb),Drs(Nsb),mpiDls(Nsb))
    do q=1,Nsb
       qn  = sb_sector%qn(index=q)
       Dls(q) = sector_qn_dim(left%sectors(1),qn)
       Drs(q) = sector_qn_dim(right%sectors(1),current_target_qn - qn)
       mpiDls(q) = Dls(q)/MpiSize
       if(MpiRank < mod(Dls(q),MpiSize))mpiDls(q) = mpiDls(q)+1
    enddo
    if(MpiStatus.AND.verbose>4.AND.(MpiComm/=Mpi_Comm_Null).AND.MpiSize>=1)then
       if(MpiMaster)write(*,*)"         mpiRank,          mpiDls        -        mpiDl        -      mpiL    -   mpiOffset"
       do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)write(*,*)MpiRank,mpiDls,"-",&
               Drs(:)*mpiDls(:),"-",sum(Drs(:)*mpiDls(:)),"-",&
               (sum(Drs(1:q-1)*mpiDls(1:q-1)),q=1,Nsb)
       enddo
       call Barrier_MPI(MpiComm)
    endif
    deallocate(Dls,Drs,mpiDls)
#endif
#endif
    !
    !IF PRESENT HMAT: get SB_H sparse > dump it to dense Hmat > return
    if(present(Hmat))then
       if(allocated(Hmat))deallocate(Hmat)
       allocate(Hmat(m_sb,m_sb));Hmat=zero
       !Nullify HxV function pointer:
       spHtimesV_p => null()
       !
       !>Build Sparse Hsb:
       if(MpiMaster)call start_timer("get H_sb Dense: LAPACK")
       call Setup_SuperBlock_Sparse() !<- no MPI here
       if(MpiMaster)call stop_timer()
       !
       !Dump Hsb to dense matrix as required:
       call spHsb%dump(Hmat)
       return
    endif
    !
    !Build SuperBLock HxV operation: stored or direct
    select case(sparse_H)
    case(.true.)
       if(MpiMaster)call start_timer("get H_sb Sparse: ARPACK")
       call Setup_SuperBlock_Sparse() !<- no MPI here yet
       if(MpiMaster)call stop_timer()
       !
       !Set HxV function pointer:
       spHtimesV_p => spMatVec_sparse_main
       !
    case(.false.)
       if(MpiMaster)call start_timer("get H_sb Direct: ARPACK")
       call Setup_SuperBlock_Direct() !<- SETUP MPI here
       if(MpiMaster)call stop_timer()
       !
       !Set HxV function pointer:
       spHtimesV_p => spMatVec_direct_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_p => spMatVec_MPI_direct_main
#endif
       !
    end select
  end subroutine sb_build_Hv





  subroutine sb_delete_Hv()
    !
    spHtimesV_p => null()
    !
    call spHsb%free()
    if(allocated(Hleft))then
       do concurrent(i=1:size(Hleft))
          call Hleft(i)%free()
       enddo
       deallocate(Hleft)
    endif
    if(allocated(Hright))then
       do concurrent(i=1:size(Hright))
          call Hright(i)%free()
       enddo
       deallocate(Hright)
    endif
    if(allocated(A))then
       do concurrent(i=1:size(A,1),j=1:size(A,2))
          call A(i,j)%free()
       enddo
       deallocate(A)
    endif
    if(allocated(B))then
       do concurrent(i=1:size(B,1),j=1:size(B,2))
          call B(i,j)%free()
       enddo
       deallocate(B)
    endif
    if(allocated(Dls))deallocate(Dls)
    if(allocated(Drs))deallocate(Drs)
    if(allocated(Offset))deallocate(Offset)
    if(allocated(RowOffset))deallocate(RowOffset)
    if(allocated(ColOffset))deallocate(ColOffset)
    if(allocated(mpiDls))deallocate(mpiDls)
    if(allocated(mpiDrs))deallocate(mpiDrs)
    if(allocated(mpiDl))deallocate(mpiDl)
    if(allocated(mpiDr))deallocate(mpiDr)
    if(allocated(mpiOffset))deallocate(mpiOffset)
    if(allocated(mpiRowOffset))deallocate(mpiRowOffset)
    if(allocated(mpiColOffset))deallocate(mpiColOffset)
  end subroutine sb_delete_Hv






  function sb_vecDim_Hv() result(vecDim)
    integer                          :: vecDim           !vector or vector chunck dimension
    real(8),dimension(:),allocatable :: qn
    integer                          :: q
    !

#ifdef _MPI

    if(MpiStatus)then
       do q=1,size(sb_sector)
          mpiDls(q) = Dls(q)/MpiSize
          if(MpiRank < mod(Dls(q),MpiSize))mpiDls(q) = mpiDls(q)+1
          mpiDl(q)  = Drs(q)*mpiDls(q)
       enddo
    else
       do q=1,size(sb_sector)
          mpiDl(q) = Drs(q)*Dls(q)
       enddo
       if(sum(mpiDl)/=size(sb_states))stop "sb_vecDim_Hv error: no MPI but vecDim != m_sb"
    endif
#else
    do q=1,size(sb_sector)
       mpiDl(q) = Drs(q)*Dls(q)
    enddo
    if(sum(mpiDl)/=size(sb_states))stop "sb_vecDim_Hv error: no MPI but vecDim != m_sb"
#endif
    !
    vecDim = sum(mpiDl)
    !
  end function sb_vecDim_Hv





END MODULE DMRG_SUPERBLOCK
