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
     procedure,pass     :: set_basis => set_basis_block
     procedure,pass     :: show      => show_block
     procedure,pass     :: is_valid  => is_valid_block
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
       self%sectors(i)   = b%sectors(i)
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
  !PURPOSE:  Put a QN array in the site
  !+------------------------------------------------------------------+
  subroutine set_basis_block(self,basis,indx)
    class(block)         :: self
    type(tbasis)         :: basis
    integer,optional     :: indx
    integer              :: indx_
    indx_=1;if(present(indx))indx_=indx
    if(indx_<1.OR.indx_>size(self%sectors))stop "SET_SECTORS_BLOCK ERROR: indx out of range"
    self%sectors(indx_) = sectors_list( basis )
  end subroutine set_basis_block



  !##################################################################
  !##################################################################
  !              OPERATIONS / ASSIGNEMENTS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Equality between two blocks
  !+------------------------------------------------------------------+
  subroutine block_equal_block(A,B)
    type(block),intent(inout) :: A
    type(block),intent(in)    :: B
    call A%free
    A%length    = B%length
    A%Dim       = B%Dim
    A%operators = B%operators
    allocate(A%sectors(size(B%sectors)))
    do i=1,size(A%sectors)
       A%sectors(i)   = B%sectors(i)
    enddo
  end subroutine block_equal_block


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





