MODULE DMRG_SUPERBLOCK_SPARSE
  USE VARS_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK_COMMON
  implicit none
  private

  public :: Setup_SuperBlock_Sparse
  public :: spMatVec_sparse_main
  
contains



  !##################################################################
  !              SETUP THE SUPERBLOCK HAMILTONIAN
  !                      SPARSE MODE
  !    H^SB = H^L x 1^R  + 1^L x H^R + H^LR
  !    H^LR = sum_p A_p x B_p
  ! 
  ! * sparse: get the sparse global SB Hamiltonian spHsb
  !##################################################################
  subroutine Setup_SuperBlock_Sparse()
    integer                         :: m_left,m_right
    character(len=:),allocatable    :: type
    type(sparse_matrix)             :: H2
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Sparse ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Sparse ERROR: Missing right.H operator in the list"
    !
    type=left%type()
    if(type/=right%type())&
         stop "Setup_SuperBlock_Sparse ERROR: left.Type != right.Type"
    !
    m_left = left%dim
    m_right= right%dim
    !
    select case(type)
    case default;stop "Setup_SuperBlock_Sparse ERROR: wrong left/right.Type"
    case ("spin","s")
       H2 = connect_spin_blocks(left,right,sb_states)
    case ("fermion","f,","electron","e")
       H2 = connect_fermion_blocks(left,right,sb_states)
    end select
    !
    spHsb = H2 + &
         sp_kron(left%operators%op("H"),id(m_right),sb_states) + &
         sp_kron(id(m_left),right%operators%op("H"),sb_states)  
    !
  end subroutine Setup_SuperBlock_Sparse



  !##################################################################
  !              SuperBlock MATRIX-VECTOR PRODUCTS
  !              using shared quantities in GLOBAL
  !##################################################################
  subroutine spMatVec_sparse_main(Nloc,v,Hv)
    integer                 :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                 :: i,j,jcol
    Hv=zero
    do i=1,Nloc
       matmul: do jcol=1, spHsb%row(i)%Size
          val = spHsb%row(i)%vals(jcol)
          j   = spHsb%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end subroutine spMatVec_sparse_main



END MODULE DMRG_SUPERBLOCK_SPARSE
