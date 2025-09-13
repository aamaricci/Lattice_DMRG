MODULE DMRG_SUPERBLOCK
  USE DMRG_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK_SETUP
  implicit none
  private


  public :: sb_get_states
  public :: sb_get_states_
  public :: sb_diag

  public :: sb_build_Hv
  public :: sb_vecDim_Hv
  public :: sb_delete_Hv

  integer                          :: ispin
  integer                          :: iorb,jorb
  integer                          :: io,jo
  integer,dimension(:),allocatable :: sb_states_tmp
  type(sectors_list)               :: sb_sector_tmp


contains



  !-----------------------------------------------------------------!
  ! Purpose: build the list of states compatible with the specified
  ! quantum numbers
  !SUPERBLOCK SHARED THINGS:
  ! integer,dimension(:),allocatable               :: sb_states
  ! type(sectors_list)                             :: sb_sector
  !-----------------------------------------------------------------!
  subroutine sb_get_states()
    integer                          :: ql,iright
    integer                          :: ir,il,istate,unit,k,Nl,Nr
    real(8),dimension(:),allocatable :: left_qn,right_qn
    integer,dimension(:),allocatable :: left_map,right_map
    !
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> d733a73 (Updated code.)
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock get states"
#endif
    !
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
>>>>>>> d733a73 (Updated code.)
    call start_timer("Get SB states")
=======
    if(MpiMaster)call start_timer("Get SB states")
>>>>>>> 4f09a08 (Extended the MPI algorithm to measure operators.)
    !
    !INIT SB STATES OBJECTS:
    if(allocated(sb_states))deallocate(sb_states)
    call sb_sector%free()
    !
#ifdef _DEBUG
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
    if(verbose>5)unit = fopen('SB_list_'//str(left%length)//'.dat')
=======
    unit = fopen('SB_list_'//str(left%length)//'.dat')
>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
    if(verbose>5)unit = fopen('SB_list_'//str(left%length)//'.dat')
>>>>>>> d733a73 (Updated code.)
=======
    if(MpiMaster.AND.verbose>5)unit = fopen('SB_list_'//str(left%length)//'.dat')
>>>>>>> 4f09a08 (Extended the MPI algorithm to measure operators.)
#endif
    do ql=1,size(left%sectors(1))
       left_qn   = left%sectors(1)%qn(index=ql)
       right_qn  = current_target_qn - left_qn
       if(.not.right%sectors(1)%has_qn(right_qn))cycle
       !
       left_map  = left%sectors(1)%map(qn=left_qn)
       right_map = right%sectors(1)%map(qn=right_qn)
       !
#ifdef _DEBUG
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> d733a73 (Updated code.)
       if(verbose>5)then
=======
       if(MpiMaster.AND.verbose>5)then
<<<<<<< HEAD
>>>>>>> 4f09a08 (Extended the MPI algorithm to measure operators.)
          write(unit,*)""
=======
          write(unit,*)" "
<<<<<<< HEAD
>>>>>>> c8aed3d (Updated code.)
          write(unit,*)left_qn,right_qn
          write(unit,*)size(left_map),size(right_map)
=======
          write(unit,*)left_qn,right_qn,size(left_map),size(right_map)
>>>>>>> d4c1240 (Porting sb_get_states to MPI.)
       endif
<<<<<<< HEAD
=======
       write(unit,*)""
       write(unit,*)left_qn,right_qn
       write(unit,*)size(left_map),size(right_map)
>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
>>>>>>> d733a73 (Updated code.)
#endif
       ! do i=1,size(left_map)
       !    do j=1,size(right_map)
       Nl = size(left_map)
       Nr = size(right_map)
       do k=1,Nl*Nr
          ir = mod(k,Nr);if(ir==0)ir=Nr
          il = (k-1)/Nr+1
          istate=right_map(ir) + (left_map(il)-1)*right%Dim
          call append(sb_states, istate)
          call sb_sector%append(qn=left_qn,istate=size(sb_states))
#ifdef _DEBUG
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
             if(verbose>5)write(unit,*)left_map(i),right_map(j),istate
=======
             if(MpiMaster.AND.verbose>5)write(unit,*)left_map(i),right_map(j),istate
>>>>>>> 4f09a08 (Extended the MPI algorithm to measure operators.)
=======
          if(MpiMaster.AND.verbose>5)write(unit,*)left_map(il),right_map(ir),istate,size(sb_states)
>>>>>>> d4c1240 (Porting sb_get_states to MPI.)
#endif
       enddo
    enddo
    ! enddo
#ifdef _DEBUG
<<<<<<< HEAD
    if(verbose>5)close(unit)
=======
             write(unit,*)left_map(i),right_map(j),istate
=======
             if(verbose>5)write(unit,*)left_map(i),right_map(j),istate
>>>>>>> d733a73 (Updated code.)
#endif
          enddo
       enddo
    enddo
#ifdef _DEBUG
<<<<<<< HEAD
    close(unit)
>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
    if(verbose>5)close(unit)
>>>>>>> d733a73 (Updated code.)
=======
    if(MpiMaster.AND.verbose>5)close(unit)
>>>>>>> 4f09a08 (Extended the MPI algorithm to measure operators.)
#endif
    !
    ! if(MpiMaster)print*,"master TRUE sb_sector:"
    if(MpiMaster)then
       do k=1,size(sb_states)
          write(200,*)k,sb_states(k)
       enddo
    endif
    if(MpiMaster)call stop_timer()
    !
  end subroutine sb_get_states









  subroutine sb_get_states_()
    integer                          :: ql,iright
    integer                          :: ir,il,istate,unit,k,Nsl
    real(8),dimension(:),allocatable :: left_qn,right_qn
    integer,dimension(:),allocatable :: left_map,right_map
    integer,dimension(:),allocatable :: Nl,Nr,Offset,Q,R,Qoffset
    integer :: istart,iend
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock get states 2"
#endif
    !
    if(MpiMaster)call start_timer("Get SB states")
    !
    !INIT SB STATES OBJECTS:
    if(allocated(sb_states))deallocate(sb_states)
    call sb_sector%free()
#ifdef _MPI
    if(allocated(sb_states_tmp))deallocate(sb_states_tmp)
    call sb_sector_tmp%free()
#endif
    !
    Nsl = size(left%sectors(1))
    allocate(Nl(Nsl),Nr(Nsl),Offset(Nsl),Q(Nsl),R(Nsl),Qoffset(Nsl))
    !
    do ql=1,Nsl
       left_qn   = left%sectors(1)%qn(index=ql)
       right_qn  = current_target_qn - left_qn
       if(.not.right%sectors(1)%has_qn(right_qn))cycle
       !
       left_map  = left%sectors(1)%map(qn=left_qn)
       right_map = right%sectors(1)%map(qn=right_qn)
       !
       write(100+MpiRank,*)" "
       ! write(100+MpiRank,*)left_qn,right_qn,size(left_map),size(right_map)

       Nl(ql) = size(left_map)
       Nr(ql) = size(right_map)
       Offset(ql)=sum(Nl(1:ql-1)*Nr(1:ql-1))
       Istart = 1
       Iend   = Nr(ql)*Nl(ql)
       !
       Q(ql) = (Nr(ql)*Nl(ql))/MpiSize
       R(ql) = 0
       if(MpiRank == MpiSize-1)R(ql)=mod(Nr(ql)*Nl(ql),MpiSize)
       Qoffset(ql) = sum(Q(1:ql-1)+R(1:ql-1))
       !
       Istart = 1+MpiRank*Q(ql)
       Iend   = (MpiRank+1)*Q(ql) + R(ql)
       
       do k=Istart,Iend
          ir = mod(k,Nr(ql));if(ir==0)ir=Nr(ql)
          il = (k-1)/Nr(ql)+1
          istate=right_map(ir) + (left_map(il)-1)*right%Dim
#ifdef _MPI
          if(MpiStatus)then
             call append(sb_states_tmp, istate)
             call sb_sector_tmp%append(qn=left_qn,istate=k+Offset(ql))
          else
             call append(sb_states, istate)
             call sb_sector%append(qn=left_qn,istate=size(sb_states))
          endif
          write(100+MpiRank,*)k,ql,istart,iend,left_map(il),right_map(ir),istate,k+Offset(ql)
#endif
       enddo
    enddo
    !
    ! call sb_sector_tmp%show()
    if(MpiMaster)call stop_timer()

    !If MPI we need to recollect all the data:
    !
#ifdef _MPI
    if(MpiStatus)then
       call allgather_sb_states()
       if(MpiMaster)then
          do k=1,size(sb_states)
             write(201,*)k,sb_states(k)
          enddo
       endif

    endif
#endif
    !

  contains


    subroutine allgather_sb_states()
      integer                          :: i,Ntot,ql
      integer                          :: itmp_start,itmp_end,i_start,i_end
      integer,dimension(:),allocatable :: pCounts,pOffset
      integer                          :: MpiIerr
      !
      Ntot = sum(Nl*Nr)         !size of the sb_array
      allocate(sb_states(Ntot))
      !
      allocate(pCounts(0:MpiSize-1)) ; pCounts=0
      allocate(pOffset(0:MpiSize-1)) ; pOffset=0
      !
      do ql=1,Nsl
         pCounts=0 
         pOffset=0
         call MPI_AllGather(Q(ql)+R(ql),1,MPI_INTEGER,pCounts,1,MPI_INTEGER,MpiComm,MpiIerr)
         !
         !Get Offset:
         pOffset(0)=0
         do i=1,MpiSize-1
            pOffset(i) = pOffset(i-1) + pCounts(i-1)
         enddo
         !
         itmp_start = 1+Qoffset(ql)
         itmp_end   = Q(ql) + R(ql) +Qoffset(ql)
         !
         i_start = 1 + Offset(ql)
         i_end   = Nl(ql)*Nr(ql) + OffSet(ql)
         !
         write(300+MpiRank,*)itmp_start,itmp_end,i_start,i_end
         !
         call MPI_AllGatherv(&
              sb_states_tmp(itmp_start:itmp_end),Q(ql)+R(ql),MPI_INTEGER,&
              sb_states(i_start:i_end),pCounts,pOffset,MPI_INTEGER,MpiComm,MpiIerr)
      enddo


    end subroutine allgather_sb_states


  end subroutine sb_get_states_















  !-----------------------------------------------------------------!
  ! Purpose: Diagonalize the SuperBlock problem.
  !-----------------------------------------------------------------!
  subroutine sb_diag()
    integer                               :: m_sb
    integer                               :: Nitermax,Neigen,Nblock
    real(8),dimension(:),allocatable      :: evals
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hsb
#else
    real(8),dimension(:,:),allocatable    :: Hsb
#endif
    integer                               :: vecDim,Nloc,m_tmp
    logical                               :: exist,lanc_solve
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> d733a73 (Updated code.)
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock diagonalization"
#endif
    !
<<<<<<< HEAD
=======
    !    
>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
>>>>>>> d733a73 (Updated code.)
    m_sb = size(sb_states)
    if(m_sb==0)stop "sb_diag ERROR: size(sb_states)==0"
    !
    !Set Lanczos params
    Neigen   = min(m_sb,lanc_neigen)
    Nitermax = min(m_sb,lanc_niter)
    Nblock   = min(m_sb,lanc_ncv_factor*Lanc_Neigen+lanc_ncv_add)
    !
    !Decide how to operate on H_sb
    lanc_solve = .true.
    if(Neigen==m_sb)lanc_solve=.false.
    if(m_sb  <= lanc_dim_threshold)lanc_solve=.false.
    !
    !Allocate EigPairs:
    if(allocated(gs_energy))deallocate(gs_energy)
    allocate(gs_energy(Neigen));gs_energy=zero
    !
    if(MpiMaster)call start_timer("Diag H_sb")
    if(lanc_solve)then          !Use (P)-Arpack
       !
       call sb_build_Hv()
       !
       if(allocated(gs_vector))deallocate(gs_vector)
       vecDim = sb_vecDim_Hv()
       allocate(gs_vector(vecDim,Neigen));gs_vector=zero
       !
#ifdef _MPI
       if(MpiStatus)then          
          call sp_eigh(MpiComm,spHtimesV_p,gs_energy,gs_vector,&
               Nblock,&
               Nitermax,&
               tol=lanc_tolerance,&
               iverbose=(verbose>4))
       else
          call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
               Nblock,&
               Nitermax,&
               tol=lanc_tolerance,&
               iverbose=(verbose>4))
       endif
#else
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
       vecDim = sb_vecDim_Hv()
       allocate(gs_vector(vecDim,Neigen));gs_vector=zero
       allocate(evals(m_sb))
<<<<<<< HEAD
<<<<<<< HEAD
       call eigh(Hsb,evals)
<<<<<<< HEAD
<<<<<<< HEAD
       gs_vector(:,1:Neigen) = Hsb(:,1:Neigen)
       gs_energy(1:Neigen)   = evals(1:Neigen)
=======
       gs_vector(:,1:Lanc_Neigen) = Hsb(:,1:Lanc_Neigen)
       gs_energy(1:Lanc_Neigen)   = evals(1:Lanc_Neigen)
>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
       gs_vector(:,1:Neigen) = Hsb(:,1:Neigen)
       gs_energy(1:Neigen)   = evals(1:Neigen)
>>>>>>> d733a73 (Updated code.)
=======
       !
=======

>>>>>>> c8aed3d (Updated code.)
       if(MpiMaster)call eigh(Hsb,evals)
#ifdef _MPI
       if(MpiStatus)then
          call Bcast_MPI(MpiComm,evals)
          call sb_build_dims()
          call scatter_vector_MPI(MpiComm,Hsb(:,1:Neigen),gs_vector)
          gs_energy(1:Neigen)   = evals(1:Neigen)
          call sb_delete_dims()
       else
          gs_vector(:,1:Neigen) = Hsb(:,1:Neigen)
          gs_energy(1:Neigen)   = evals(1:Neigen)
       endif
#else
       gs_vector(:,1:Neigen) = Hsb(:,1:Neigen)
       gs_energy(1:Neigen)   = evals(1:Neigen)          
#endif
       !
>>>>>>> 6deacad (Updating code, implementing MPI for direct Hv)
       deallocate(Hsb,evals)
       !
    endif
    if(MpiMaster)call stop_timer()
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
       !>Build Sparse Hsb: (Root only)
       if(MpiMaster)then
          call start_timer("get H_sb Dense: LAPACK")
          call Setup_SuperBlock_Sparse() !<- no MPI here
          call stop_timer()
          !
          !Dump Hsb to dense matrix as required:
          call spHsb%dump(Hmat)
       endif
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
    integer :: i,j
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
    !
    call sb_delete_dims()
    if(allocated(RowOffset))deallocate(RowOffset)
<<<<<<< HEAD
    if(allocated(ColOffset))deallocate(ColOffset)    
<<<<<<< HEAD
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
    real(8),dimension(:,:),allocatable,optional :: Hmat
#endif
    integer                                        :: m_sb
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> d733a73 (Updated code.)
    !
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: SuperBlock build H*v"
#endif
    !
<<<<<<< HEAD
=======

>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
>>>>>>> d733a73 (Updated code.)
    if(.not.allocated(sb_states))stop "build_Hv_superblock ERROR: sb_states not allocated"
    m_sb = size(sb_states)

    !IF PRESENT HMAT: get SB_H sparse > dump it to dense Hmat > return
    if(present(Hmat))then
       if(allocated(Hmat))deallocate(Hmat)
       allocate(Hmat(m_sb,m_sb));Hmat=zero
       !Nullify HxV function pointer:
       spHtimesV_p => null()
       !
       !>Build Sparse Hsb:
       call start_timer("get H_sb Dense: LAPACK")
       call Setup_SuperBlock_Sparse()
       call stop_timer()
       !
       !Dump Hsb to dense matrix as required:
       call spHsb%dump(Hmat)
       return
    endif
    !
    !Build SuperBLock HxV operation: stored or direct
    select case(sparse_H)
    case(.true.)
       call start_timer("get H_sb Sparse: ARPACK")
       call Setup_SuperBlock_Sparse()
       call stop_timer()
       !
       !Set HxV function pointer:
       spHtimesV_p => spMatVec_sparse_main
       !
    case(.false.)
       call start_timer("get H_sb Direct: ARPACK")
       call Setup_SuperBlock_Direct()
       call stop_timer()
       !
       !Set HxV function pointer:
       spHtimesV_p => spMatVec_direct_main
       !
    end select
  end subroutine sb_build_Hv



=======
=======
    if(allocated(ColOffset))deallocate(ColOffset)
<<<<<<< HEAD
    if(allocated(mpiDls))deallocate(mpiDls)
    if(allocated(mpiDrs))deallocate(mpiDrs)
<<<<<<< HEAD
    if(allocated(mpiQs))deallocate(mpiQs)
>>>>>>> 6deacad (Updating code, implementing MPI for direct Hv)
=======
    if(allocated(mpiDl))deallocate(mpiDl)
    if(allocated(mpiDr))deallocate(mpiDr)
    if(allocated(mpiOffset))deallocate(mpiOffset)
<<<<<<< HEAD
>>>>>>> e9d8bd1 (Updated code.)
=======
    if(allocated(mpiRowOffset))deallocate(mpiRowOffset)
    if(allocated(mpiColOffset))deallocate(mpiColOffset)
>>>>>>> 4eb4ad9 (IT IS FUCKING WORKING!!!)
=======
#ifdef _MPI
    if(allocated(mpiRowOffset))deallocate(mpiRowOffset)
    if(allocated(mpiColOffset))deallocate(mpiColOffset)
#endif
    !
>>>>>>> c8aed3d (Updated code.)
  end subroutine sb_delete_Hv
>>>>>>> 6adc5a4 (Updated code, importing MPI)






  function sb_vecDim_Hv() result(vecDim)
    integer                          :: vecDim           !vector or vector chunck dimension
    real(8),dimension(:),allocatable :: qn
    integer                          :: q,Nsb
    !
    Nsb  = size(sb_sector)
    if(allocated(Dls))deallocate(Dls)
    if(allocated(Drs))deallocate(Drs)
    if(allocated(mpiDls))deallocate(mpiDls)
    if(allocated(mpiDl))deallocate(mpiDl)
    allocate(Dls(Nsb),Drs(Nsb),mpiDls(Nsb),mpiDl(Nsb))
    do q=1,Nsb
       qn     = sb_sector%qn(index=q)
       Dls(q) = sector_qn_dim(left%sectors(1),qn)
       Drs(q) = sector_qn_dim(right%sectors(1),current_target_qn - qn)
    enddo
#ifdef _MPI
    if(MpiStatus)then
       do q=1,Nsb
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
