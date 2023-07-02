MODULE BLOCKS
  USE SCIFOR, only: str,assert_shape
  USE AUX_FUNCS
  USE MATRIX_SPARSE
  USE TUPLE_BASIS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  implicit none
  private


  
  type block
     integer                                     :: length=0
     integer                                     :: Dim=1
     type(sectors_list),dimension(:),allocatable :: sectors
     type(operators_list)                        :: operators
   contains
     procedure,pass     :: free      => free_block
     procedure,pass     :: put       => put_op_block
     procedure,pass     :: load      => load_op_block
     procedure,pass     :: get_basis => get_basis_block
     procedure,pass     :: set_basis => set_basis_block
     procedure,pass     :: show      => show_block
     procedure,pass     :: is_valid  => is_valid_block
     procedure,pass     :: renormalize => rotate_operators_block
  end type block


  !GENERIC CONSTRUCTOR
  interface block
     module procedure :: constructor_from_scrath
     module procedure :: constructor_from_site
     module procedure :: constructor_from_block
  end interface block

  !GENERIC CONSTRUCTOR
  interface as_block
     module procedure :: constructor_from_scrath
     module procedure :: constructor_from_site
     module procedure :: constructor_from_block
  end interface as_block


  public :: block
  public :: as_block


  integer :: i,j

contains


  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Destructor
  !+------------------------------------------------------------------+
  subroutine free_block(self)
    class(block) :: self
    self%length = 0
    self%Dim   = 1
    call self%operators%free()
    if(allocated(self%sectors))then
       call self%sectors%free()
       deallocate(self%sectors)
    endif
  end subroutine free_block



  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function constructor_from_scrath(length,Dim,sectors,operators) result(self)
    integer,intent(in)              :: length
    integer,intent(in)              :: Dim
    type(sectors_list),intent(in)   :: sectors(:)
    type(operators_list),intent(in) :: operators
    type(block)                     :: self
    self%length    = length
    self%Dim       = Dim
    self%operators = operators
    allocate(self%sectors(size(sectors)))
    do i=1,size(self%sectors)
       self%sectors(i)   = sectors(i)
    enddo
  end function constructor_from_scrath


  function constructor_from_site(ssite) result(self)
    type(site),intent(in) :: ssite
    type(block)           :: self
    self%length    = 1
    self%Dim       = ssite%Dim
    self%operators = ssite%operators
    allocate(self%sectors(size(ssite%sectors)))
    do i=1,size(self%sectors)
       self%sectors(i)   = ssite%sectors(i)
    enddo
  end function constructor_from_site

  function constructor_from_block(b) result(self)
    type(block),intent(in) :: b
    type(block)           :: self
    self%length    = b%length
    self%Dim       = b%Dim
    self%operators = b%operators
    allocate(self%sectors(size(b%sectors)))
    do i=1,size(self%sectors)
       self%sectors(i) = b%sectors(i)
    enddo
  end function constructor_from_block




  !##################################################################
  !##################################################################
  !       PUT/LOAD  - GET/DUMP OPERATORS IN/FROM A LIST
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Load a sparse operator in the block dictionary
  !+------------------------------------------------------------------+
  subroutine put_op_block(self,key,op)
    class(block)                   :: self
    character(len=*),intent(in)    :: key
    type(sparse_matrix),intent(in) :: op
    call self%operators%put(str(key),op)
  end subroutine put_op_block


  !+------------------------------------------------------------------+
  !PURPOSE:  Load a dense operator in the block dictionary
  !+------------------------------------------------------------------+
  subroutine load_op_block(self,key,op)
    class(block)                      :: self
    character(len=*),intent(in)       :: key
    real(8),dimension(:,:),intent(in) :: op
    call self%operators%load(str(key),op)
  end subroutine load_op_block



  !+------------------------------------------------------------------+
  !PURPOSE:  Get Basis of the sector
  !+------------------------------------------------------------------+
  subroutine get_basis_block(self,basis,indx)
    class(block)     :: self
    type(tbasis)     :: basis
    integer,optional :: indx
    integer          :: indx_
    indx_=1;if(present(indx))indx_=indx
    if(indx_<1.OR.indx_>size(self%sectors))stop "SET_SECTORS_BLOCK ERROR: indx out of range"
    call basis%free()
    basis  = self%sectors(indx_)%basis()
  end subroutine get_basis_block


  !+------------------------------------------------------------------+
  !PURPOSE:  Put a QN array in the site
  !+------------------------------------------------------------------+
  subroutine set_basis_block(self,basis,indx)
    class(block)     :: self
    type(tbasis)     :: basis
    integer,optional :: indx
    integer          :: indx_
    indx_=1;if(present(indx))indx_=indx
    if(indx_<1.OR.indx_>size(self%sectors))stop "SET_SECTORS_BLOCK ERROR: indx out of range"
    self%sectors(indx_) = sectors_list( basis )
  end subroutine set_basis_block






  !+------------------------------------------------------------------+
  !PURPOSE:  
  !+------------------------------------------------------------------+
  subroutine rotate_operators_block(self,Umat)
    class(block)                 :: self
    real(8),dimension(:,:)       :: Umat   ![N,M]
    integer                      :: i,N,M  !N=self%dim,M=truncated dimension
    type(sparse_matrix)          :: Op
    character(len=:),allocatable :: key
    !
    N = size(Umat,1)
    M = size(Umat,2)
    if(N/=self%dim) stop "rotate_operators_block error: size(Umat,1) != self.dim"
    do i=1,size(self%operators)
       key = self%operators%key(index=i)
       Op  = self%operators%op(index=i)
       call self%put(str(key),rotate_and_truncate(Op,Umat,N,M))
    enddo
    self%dim = M
    !
    call Op%free()
    !
  end subroutine rotate_operators_block
  !
  function rotate_and_truncate(Op,trRho,N,M) result(RotOp)
    type(sparse_matrix),intent(in) :: Op
    real(8),dimension(N,M)         :: trRho  ![Nesys,M]
    integer                        :: N,M
    type(sparse_matrix)            :: RotOp
    real(8),dimension(M,M)         :: Umat
    real(8),dimension(N,N)         :: OpMat
    N = size(trRho,1)
    M = size(trRho,2)
    if( any( shape(Op) /= [N,N] ) ) stop "rotate_and_truncate error: shape(Op) != [N,N] N=size(Rho,1)"
    OpMat= Op%as_matrix()
    Umat = matmul( transpose(trRho), matmul(OpMat,trRho)) ![M,N].[N,N].[N,M]=[M,M]
    call RotOp%load( Umat )
  end function rotate_and_truncate


  !##################################################################
  !##################################################################
  !              OPERATIONS / ASSIGNEMENTS
  !##################################################################
  !##################################################################
  ! !+------------------------------------------------------------------+
  ! !PURPOSE:  Equality between two blocks
  ! !+------------------------------------------------------------------+
  ! subroutine block_equal_block(A,B)
  !   type(block),intent(inout) :: A
  !   type(block),intent(in)    :: B
  !   call A%free
  !   A%length    = B%length
  !   A%Dim       = B%Dim
  !   A%operators = B%operators
  !   allocate(A%sectors(size(B%sectors)))
  !   do i=1,size(A%sectors)
  !      A%sectors(i)   = B%sectors(i)
  !   enddo
  ! end subroutine block_equal_block


  function is_valid_block(self,nobasis) result(bool)
    class(block)                          :: self
    logical,optional                      :: nobasis    
    logical                               :: bool
    logical                               :: nobasis_
    integer,dimension(size(self%sectors)) :: Lvec
    !
    nobasis_=.false.;if(present(nobasis))nobasis_=nobasis
    !
    bool = self%operators%is_valid(self%Dim)
    if(nobasis_)return
    do i=1,size(self%sectors)
       Lvec(i) = len(self%sectors(i))
       write(900,*)Lvec(i)
    enddo
    write(900,*)self%dim,product(Lvec)
    bool=bool.AND.(self%dim==product(Lvec))
  end function is_valid_block




  !##################################################################
  !##################################################################
  !              SHOW 
  !##################################################################
  !##################################################################
  subroutine show_block(self,fmt,wOP)
    class(block)              :: self
    character(len=*),optional :: fmt
    logical,optional          :: wOP
    character(len=32)         :: fmt_
    logical :: wOP_
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    wOP_=.true.;if(present(wOP))wOP_=wOP
    write(*,"(A15,I6)")"Block Length  =",self%length
    write(*,"(A15,I6)")"Block Dim     =",self%Dim
    write(*,"(A15,I6)")"Block Sectors =",size(self%sectors)
    do i=1,size(self%sectors)
       write(*,"(A14,I6)")"Block Sector  =",i
       call self%sectors(i)%show()
    enddo
    if(wOP_)then
       write(*,"(A14)")"Block Ops     ="
       call self%operators%show(fmt=fmt_)
    endif
  end subroutine show_block









END MODULE BLOCKS





