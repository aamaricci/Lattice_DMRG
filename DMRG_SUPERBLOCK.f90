MODULE DMRG_SUPERBLOCK
  USE VARS_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK_COMMON
  USE DMRG_SUPERBLOCK_SPARSE
  USE DMRG_SUPERBLOCK_DIRECT
  implicit none
  private


  public :: sb_get_states
  public :: sb_diag


contains



  !-----------------------------------------------------------------!
  ! Purpose: build the list of states compatible with the specified
  ! quantum numbers
  !-----------------------------------------------------------------!
  subroutine sb_get_states()
    integer                          :: ileft,iright
    integer                          :: i,j,istate
    real(8),dimension(:),allocatable :: left_qn,right_qn
    integer,dimension(:),allocatable :: left_map,right_map
    !
    call start_timer()
    !
    if(allocated(sb_states))deallocate(sb_states)
    !
    call sb_sector%free()
    !
#ifdef _DEBUG
    open(102,file='SB_list_'//str(left%length)//'.dat')
#endif
    do ileft=1,size(left%sectors(1))
       left_qn   = left%sectors(1)%qn(index=ileft)
       right_qn  = current_target_qn - left_qn
       if(.not.right%sectors(1)%has_qn(right_qn))cycle
       !
       left_map = left%sectors(1)%map(qn=left_qn)
       right_map = right%sectors(1)%map(qn=right_qn)
       !
#ifdef _DEBUG
       write(102,*)left_qn,right_qn
       write(102,*)size(left_map),size(right_map)
#endif
       do i=1,size(left_map)
          do j=1,size(right_map)
             istate=right_map(j) + (left_map(i)-1)*right%Dim
#ifdef _DEBUG
             write(102,*)left_map(i),right_map(j),istate,right%Dim
#endif
             call append(sb_states, istate)
             call sb_sector%append(qn=left_qn,istate=size(sb_states))
          enddo
       enddo
#ifdef _DEBUG
       write(102,*)""
#endif
       ! call eta(ileft,size(left%sectors(1)))
    enddo
#ifdef _DEBUG
    close(102)
#endif
    !
    call stop_timer("Get SB states")
    !
    call sb_sector%show(file='SB_sector_'//str(left%length)//'.dat')
  end subroutine sb_get_states




  !##################################################################
  !
  !         SETUP THE SUPERBLOCK HAMILTONIAN PROBLEM
  ! . if Hmat: returb H^SB as dense matrix there for Lapack use
  ! . if sparse_H = T: build H^SB as sparse matrix
  ! . if sparse_H = F: setup H^SB terms and blocks for H*v procedure
  !##################################################################
  subroutine sb_build_Hv(Hmat)
    real(8),dimension(:,:),allocatable,optional :: Hmat
    integer                                     :: m_sb

    if(.not.allocated(sb_states))stop "build_Hv_superblock ERROR: sb_states not allocated"
    m_sb = size(sb_states)

    !IF PRESENT HMAT: get SB_H sparse > dump it to dense Hmat > return
    if(present(Hmat))then
       if(allocated(Hmat))deallocate(Hmat)
       allocate(Hmat(m_sb,m_sb));Hmat=0d0
       !Nullify HxV function pointer:
       spHtimesV_p => null()
       !
       !>Build Sparse Hsb:
       call start_timer("get H_sb")
       call Setup_SuperBlock_Sparse()
       call stop_timer("Done H_sb")
       !
       !Dump Hsb to dense matrix as required:
       call spHsb%dump(Hmat)
       return
    endif
    !
    !Build SuperBLock HxV operation: stored or direct
    select case(sparse_H)
    case(.true.)
       call start_timer("get H_sb")
       call Setup_SuperBlock_Sparse()
       call stop_timer("Done H_sb")
       !
       !Set HxV function pointer:
       spHtimesV_p => spMatVec_sparse_main
       !
    case(.false.)
       call start_timer("get H_sb")
       call Setup_SuperBlock_Direct()
       call stop_timer("Done H_sb")
       !
       !Set HxV function pointer:
       spHtimesV_p => spMatVec_direct_main
       !
    end select
  end subroutine sb_build_Hv




  !-----------------------------------------------------------------!
  ! Purpose: Diagonalize the SuperBlock problem.
  !-----------------------------------------------------------------!
  subroutine sb_diag()
    integer                            :: m_sb
    integer                            :: Nitermax,Neigen,Nblock
    real(8),dimension(:),allocatable   :: evals
    real(8),dimension(:,:),allocatable :: Hsb
    logical                            :: exist,lanc_solve
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
    if(Lanc_Neigen==m_sb)lanc_solve=.false.
    if(m_sb <= lanc_dim_threshold)lanc_solve=.false.
    !
    !Allocate EigPairs
    if(allocated(gs_energy))deallocate(gs_energy)
    if(allocated(gs_vector))deallocate(gs_vector)
    allocate(gs_energy(Neigen))     ;gs_energy=0d0
    allocate(gs_vector(m_sb,Neigen));gs_vector=0d0
    !
    call start_timer()
    if(lanc_solve)then
       call sb_build_Hv()
       call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
            Nblock,&
            Nitermax,&
            tol=lanc_tolerance,&
            iverbose=.false.)
    else !use LAPACK
       call sb_build_Hv(Hsb)
       allocate(evals(m_sb))
       call eigh(Hsb,evals)
       gs_vector(:,1:Lanc_Neigen) = Hsb(:,1:Lanc_Neigen)
       gs_energy(1:Lanc_Neigen)   = evals(1:Lanc_Neigen)
       deallocate(Hsb,evals)
    endif
    call stop_timer("Diag H_sb")

    
    !Free Memory
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
  end subroutine sb_diag
















END MODULE DMRG_SUPERBLOCK
