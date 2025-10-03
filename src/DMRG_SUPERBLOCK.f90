MODULE DMRG_SUPERBLOCK
  USE DMRG_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK_SETUP
  USE MPI
  implicit none
  private


  public :: sb_get_states
  public :: sb_diag

  public :: sb_build_Hv
  public :: sb_vecDim_Hv
  public :: sb_delete_Hv

  integer                          :: ispin
  integer                          :: iorb,jorb
  integer                          :: io,jo



contains




  !-----------------------------------------------------------------!
  ! purpose: build the list of states compatible with the specified
  ! quantum numbers
  !SUPERBLOCK SHARED THINGS:
  ! integer,dimension(:),allocatable               :: sb_states
  ! type(sectors_list)                             :: sb_sector
  !-----------------------------------------------------------------!
  subroutine sb_get_states()
    integer                          :: ql,iright
    integer                          :: ir,il,istate,unit,k,Nsl
    real(8),dimension(:),allocatable :: left_qn,right_qn
    integer,dimension(:),allocatable :: left_map,right_map
    integer,dimension(:),allocatable :: Astates
    integer,dimension(:),allocatable :: Nl,Nr,Offset
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock get states"
#endif
    !
    if(MpiMaster)call start_timer("Build SB states")
    t0=t_start()
    !
    !INIT SB STATES OBJECTS:
    if(allocated(sb_states))deallocate(sb_states)
    call sb_sector%free()
    !
#ifdef _DEBUG
    if(MpiMaster.AND.verbose>5)unit = fopen('SB_list_'//str(left%length)//'.dat')
#endif
    !
    Nsl = size(left%sectors(1))
    allocate(Nl(Nsl),Nr(Nsl),Offset(Nsl))
    !
    do ql=1,size(left%sectors(1))
       left_qn   = left%sectors(1)%qn(index=ql)
       right_qn  = current_target_qn - left_qn
       if(.not.right%sectors(1)%has_qn(right_qn))cycle
       !
       left_map  = left%sectors(1)%map(qn=left_qn)
       right_map = right%sectors(1)%map(qn=right_qn)
       !
#ifdef _DEBUG
       if(MpiMaster.AND.verbose>5)then
          write(unit,*)" "
          write(unit,*)left_qn,right_qn,size(left_map),size(right_map)
       endif
#endif
       !
       Nl(ql) = size(left_map)
       Nr(ql) = size(right_map)
       Offset(ql)=sum(Nl(1:ql-1)*Nr(1:ql-1))
       if(allocated(Astates))deallocate(Astates)
       allocate(Astates(Nr(ql)*Nl(ql)));Astates=0
       !
       do k=1,Nr(ql)*Nl(ql)
          ir = mod(k,Nr(ql));if(ir==0)ir=Nr(ql)
          il = (k-1)/Nr(ql)+1
          istate=right_map(ir) + (left_map(il)-1)*right%Dim
          call append(sb_states, istate)
          Astates(k) = k+Offset(ql) !==size(sb_states)
#ifdef _DEBUG
          if(MpiMaster.AND.verbose>5)&
               write(unit,*)left_map(il),right_map(ir),istate,size(sb_states),k+Offset(ql)
#endif
       enddo
       call sb_sector%appends(qn=left_qn,istates=Astates)
    enddo
    !
#ifdef _DEBUG
    if(MpiMaster.AND.verbose>5)close(unit)
#endif
    !
    if(MpiMaster)call stop_timer("Build SB states")
    t_sb_get_states=t_stop()
  end subroutine sb_get_states











  !-----------------------------------------------------------------!
  ! Purpose: Diagonalize the SuperBlock problem.
  !-----------------------------------------------------------------!
  subroutine sb_diag()
    integer                               :: m_sb
    integer                               :: Nitermax,Neigen,Nblock
    real(8),dimension(:),allocatable      :: evals
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hsb,gs_tmp
#else
    real(8),dimension(:,:),allocatable    :: Hsb,gs_tmp
#endif
    integer                               :: vecDim,Nloc,m_tmp
    logical                               :: exist,lanc_solve,fMpi
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
    lanc_solve = .true.
    if(Neigen==m_sb)lanc_solve=.false.
    if(m_sb  <= lanc_dim_threshold)lanc_solve=.false.
    !
    !Allocate EigPairs:
    if(allocated(gs_energy))deallocate(gs_energy)
    allocate(gs_energy(Neigen));gs_energy=zero
    !
    !
    if(lanc_solve)then          !Use (P)-Arpack
       !
       call sb_build_Hv()
       !
       if(allocated(gs_vector))deallocate(gs_vector)
       vecDim = sb_vecDim_Hv()
       !
       allocate(gs_vector(vecDim,Neigen));gs_vector=zero
       !
       if(MpiMaster)call start_timer("Diag H_sb")
       t0=t_start()
#ifdef _MPI
       if(MpiStatus)then
          !This condition protect against problems too small compared to MpiSize:
          !Solve serial and scatter the result
          call sb_check_Hv(fMPI)
          if(fMpi)then
             allocate(gs_tmp(m_sb,Neigen))
             if(MpiMaster)call sp_eigh(spHtimesV_p,gs_energy,gs_tmp,&
                  Nblock,&
                  Nitermax,&
                  tol=lanc_tolerance,&
                  iverbose=(verbose>4),NumOp=NumOp)
             call Bcast_MPI(MpiComm,gs_energy)
             call sb_build_dims()
             call scatter_vector_MPI(MpiComm,gs_tmp,gs_vector)
             call sb_delete_dims()
             deallocate(gs_tmp)
          else
             call sp_eigh(MpiComm,spHtimesV_p,gs_energy,gs_vector,&
                  Nblock,&
                  Nitermax,&
                  tol=lanc_tolerance,&
                  iverbose=(verbose>4),NumOp=NumOp)
          endif
       else
          call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
               Nblock,&
               Nitermax,&
               tol=lanc_tolerance,&
               iverbose=(verbose>4),NumOp=NumOp)
       endif
#else
       call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
            Nblock,&
            Nitermax,&
            tol=lanc_tolerance,&
            iverbose=(verbose>4),NumOp=NumOp)
#endif
       if(MpiMaster)call stop_timer("Diag H_sb")
       t_sb_diag=t_stop()
       !
    else !use LAPACK
       !
       call sb_build_Hv(Hsb)
       !
       if(allocated(gs_vector))deallocate(gs_vector)
       vecDim = sb_vecDim_Hv()
       allocate(gs_vector(vecDim,Neigen));gs_vector=zero
       allocate(evals(m_sb))
       !
       if(MpiMaster)call start_timer("Diag H_sb")
       t0=t_start()
       if(MpiMaster)call eigh(Hsb,evals)
       if(MpiMaster)call stop_timer("Diag H_sb")
       t_sb_diag=t_stop()
       !
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
       deallocate(Hsb,evals)
       !
    endif
    !
    !Free Memory
    call sb_delete_Hv()
    !
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
    integer                                        :: q,m_sb,Nsb,irank,vecDim
    logical :: fMPI
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
#ifdef _DEBUG
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
#endif
    deallocate(Dls,Drs,mpiDls)
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
          call Setup_SuperBlock_Sparse() !<- no MPI here
          !
          !Dump Hsb to dense matrix as required:
          write(LogFile,"(A)")"LAPACK: spHsb.dump(Hmat)"
          call spHsb%dump(Hmat)
       endif
    else
       !
       !Build SuperBLock HxV operation: stored or direct
       select case(sparse_H)
       case(.true.)
          call Setup_SuperBlock_Sparse() !<- no MPI here yet
          spHtimesV_p => spMatVec_sparse_main
          !
       case(.false.)
          call Setup_SuperBlock_Direct() !<- SETUP MPI here
          spHtimesV_p => spMatVec_direct_main
#ifdef _MPI
          if(MpiStatus)then
             call sb_check_Hv(fMPI)
             if(fMPI)then  !this is true for all nodes at once see sb_vecDim_Hv
                write(LogFile,"(A)")"ARPACK: using SERIAL over MPI as MpiSize > N"
                spHtimesV_p => spMatVec_direct_main
             else
                spHtimesV_p => spMatVec_MPI_direct_main
             endif
          endif
#endif
          !
       end select
    endif
    !
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
    if(allocated(ColOffset))deallocate(ColOffset)
    !
  end subroutine sb_delete_Hv






  function sb_vecDim_Hv() result(vecDim)
    integer                          :: vecDim           !vector or vector chunck dimension
    real(8),dimension(:),allocatable :: qn
    integer                          :: q,Nsb,i
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
  end function sb_vecDim_Hv





  subroutine sb_check_Hv(anyZero)
    integer :: vecDim           !vector or vector chunck dimension
    integer :: ierr
    logical :: hasZero,anyZero
    !
    anyZero=.false.
    !
#ifdef _MPI
    if(MpiStatus)then
       vecDim  = sb_vecDim_Hv()
       hasZero = (vecDim == 0)  !T if a rank has vecDim==0
       call MPI_ALLREDUCE(hasZero, anyZero, 1, MPI_LOGICAL, MPI_LOR, MpiComm, ierr)
       !anyZero=T if at least 1 nodes has vecDim==0
    end if
#endif
  end subroutine sb_check_Hv

END MODULE DMRG_SUPERBLOCK















!   ! !This is the MPI version
!   ! NEEDS TO BE FURTHER TESTES AS NODES GET DIFFERENT SB_STATES IN THE END.
!   ! The AllGather is apparently wrong... needs to study this better.
!   subroutine sb_get_states()
!     integer                          :: ql,iright,i
!     integer                          :: ir,il,istate,unit,k,Nsl
!     real(8),dimension(:),allocatable :: left_qn,right_qn
!     integer,dimension(:),allocatable :: left_map,right_map
!     integer,dimension(:),allocatable :: Astates
!     integer,dimension(:),allocatable :: Nl,Nr,Offset,Q,R,Qoffset
!     integer                          :: istart,iend,ierr,Rdim
!     integer,dimension(:),allocatable :: sb_states_tmp      
!     !
! #ifdef _DEBUG
!     if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock get states"
! #endif
!     !
!     if(MpiMaster)call start_timer("Get SB states")
!     t0=t_start()
!     !
!     !INIT SB STATES OBJECTS:
!     if(allocated(sb_states))deallocate(sb_states)
!     if(allocated(sb_states_tmp))deallocate(sb_states_tmp)
!     call sb_sector%free()
!     !
! #ifdef _DEBUG
!     if(verbose>5)unit = fopen('SB_list_'//str(left%length)//"_n"//str(MpiRank)//".dat")
! #endif
!     !
!     Nsl = size(left%sectors(1))
!     allocate(Nl(Nsl),Nr(Nsl),Offset(Nsl),Q(Nsl),R(Nsl),Qoffset(Nsl))
!     !
!     do ql=1,Nsl
!        left_qn   = left%sectors(1)%qn(index=ql)
!        right_qn  = current_target_qn - left_qn
!        if(.not.right%sectors(1)%has_qn(right_qn))cycle
!        !
!        left_map  = left%sectors(1)%map(qn=left_qn)
!        right_map = right%sectors(1)%map(qn=right_qn)
!        !
! #ifdef _DEBUG
!        if(verbose>5)then
!           write(unit,*)" "
!           write(unit,*)left_qn,right_qn,size(left_map),size(right_map)
!        endif
! #endif
!        !
!        Rdim   = right%dim
!        Nl(ql) = size(left_map)
!        Nr(ql) = size(right_map)
!        Offset(ql)=sum(Nl(1:ql-1)*Nr(1:ql-1))
!        if(allocated(Astates))deallocate(Astates)
!        allocate(Astates(Nr(ql)*Nl(ql)))
!        Astates=0
!        !
!        Istart = 1
!        Iend   = Nr(ql)*Nl(ql)
! #ifdef _MPI
!        if(MpiStatus)then
!           Q(ql) = (Nr(ql)*Nl(ql))/MpiSize
!           R(ql) = 0
!           if(MpiRank == MpiSize-1)R(ql)=mod(Nr(ql)*Nl(ql),MpiSize)
!           Qoffset(ql) = sum(Q(1:ql-1)+R(1:ql-1))
!           Istart = 1+MpiRank*Q(ql)
!           Iend   = (MpiRank+1)*Q(ql) + R(ql)
!           if(MpiMaster)Rdim=right%dim
!           call Bcast_MPI(MpiComm,Rdim)
!        endif
! #endif
!        !
!        do k=Istart,Iend
!           ir = mod(k,Nr(ql));if(ir==0)ir=Nr(ql)
!           il = (k-1)/Nr(ql)+1
!           istate=right_map(ir) + (left_map(il)-1)*Rdim
! #ifdef _MPI
!           if(MpiStatus)then
!              call append(sb_states_tmp, istate)
! #ifdef _DEBUG
!              if(verbose>5)write(unit,*)k,il,ir,left_map(il),right_map(ir),istate,size(sb_states_tmp),k+Offset(ql)
! #endif
!           else
!              call append(sb_states, istate)
! #ifdef _DEBUG
!              if(verbose>5)write(unit,*)k,il,ir,left_map(il),right_map(ir),istate,size(sb_states),k+Offset(ql)
! #endif
!           endif
! #else
!           call append(sb_states, istate)
! #ifdef _DEBUG
!           if(verbose>5)write(unit,*)k,il,ir,left_map(il),right_map(ir),istate,size(sb_states),k+Offset(ql)
! #endif
! #endif
!           Astates(k) = k+Offset(ql) !==size(sb_states)
!        enddo
!        !
!        !       
! #ifdef _MPI
!        if(MpiStatus)then
!           call MPI_Allreduce(MPI_IN_PLACE, Astates, size(Astates), MPI_INTEGER, MPI_SUM, MpiComm, ierr)
!           call sb_sector%appends(qn=left_qn,istates=Astates)
!        else
!           call sb_sector%appends(qn=left_qn,istates=Astates)
!        endif
! #else
!        call sb_sector%appends(qn=left_qn,istates=Astates)
! #endif
!        !      
!     enddo
!     !
!     if(MpiMaster)call stop_timer("Get SB states")
!     t_sb_get_states=t_stop()
!     !
! #ifdef _MPI
!     !
!     if(MpiStatus)then
!        call gather_sb_states()
!        ! call Bcast_MPI(MpiComm,sb_states)
!     endif

!     allocate(sb_states_tmp(sum(Nl*Nr)))
!     sb_states_tmp=0

!     !
!     if(MpiMaster)then
!        do i=1,size(sb_states)
!           write(400,*)sb_states(i)
!        enddo
!        rewind(400)
!     endif
!     call Barrier_MPI()


!     do i=1,size(sb_states)
!        read(400,*)sb_states_tmp(i)
!     enddo
!     rewind(400)

!     do i=1,size(sb_states)
!        write(500+MpiRank,*)sb_states_tmp(i)
!     enddo
!     rewind(500+MpiRank)

!     if(any(sb_states_tmp/=sb_states))then
!        stop "ERROR In SB_states"
!     endif

!     !   
!   contains
!     !
!     subroutine gather_sb_states()
!       integer                          :: i,Ntot,ql
!       integer                          :: itmp_start,itmp_end,i_start,i_end
!       integer,dimension(:),allocatable :: pCounts,pOffset
!       integer                          :: MpiIerr
!       !
!       Ntot = sum(Nl*Nr)
!       allocate(sb_states(Ntot))
!       sb_states=0
!       !
!       allocate(pCounts(0:MpiSize-1))
!       allocate(pOffset(0:MpiSize-1))
!       !
!       do ql=1,Nsl
!          pCounts=0 
!          pOffset=0
!          call MPI_AllGather(Q(ql)+R(ql),1,MPI_INTEGER,pCounts,1,MPI_INTEGER,MpiComm,MpiIerr)
!          kb_agthr_sb_states = kb_agthr_sb_states + sum(pCounts)*DATA_kb_size
!          !
!          !Get Offset:
!          pOffset(0)=0
!          do i=1,MpiSize-1
!             pOffset(i) = pOffset(i-1) + pCounts(i-1)
!          enddo
!          !
!          itmp_start = 1+Qoffset(ql)
!          itmp_end   = Q(ql) + R(ql) + Qoffset(ql)
!          !
!          i_start = 1 + Offset(ql)
!          i_end   = Nl(ql)*Nr(ql) + OffSet(ql)
!          !
!          call MPI_AllGatherv(&
!               sb_states_tmp(itmp_start:itmp_end),Q(ql)+R(ql),MPI_INTEGER,&
!               sb_states(i_start:i_end),pCounts,pOffset,MPI_INTEGER,MpiComm,MpiIerr)
!       enddo
!       deallocate(sb_states_tmp)
!       !
!     end subroutine gather_sb_states
! #endif
!   end subroutine sb_get_states
