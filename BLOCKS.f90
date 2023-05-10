MODULE BLOCKS
  USE SCIFOR, only: str,assert_shape
  USE AUX_FUNCS
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  implicit none
  private


  type block
     integer              :: length=0
     integer              :: dim=1
     type(sectors_list)   :: sectors
     type(operators_list) :: operators
   contains
     procedure,pass     :: free     => free_block
     procedure,pass     :: put      => put_op_block
     procedure,pass     :: load     => load_op_block
     procedure,pass     :: set_sectors => set_sectors_block
     procedure,pass     :: show     => show_block
     procedure,pass     :: is_valid => is_valid_block
  end type block


  !GENERIC CONSTRUCTOR
  interface block
     module procedure :: build_block_from_scrath
     module procedure :: build_block_from_site
  end interface block

  !GENERIC CONSTRUCTOR
  interface as_block
     module procedure :: build_block_from_scrath
     module procedure :: build_block_from_site
  end interface as_block

  ! !EQUALITY 
  ! interface assignment(=)
  !    module procedure :: block_equal_block
  ! end interface assignment(=)

  public :: block
  public :: as_block
  ! public :: assignment(=)


contains


  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function build_block_from_scrath(length,dim,sectors,operators) result(self)
    integer,intent(in)              :: length
    integer,intent(in)              :: dim
    type(sectors_list),intent(in)   :: sectors
    type(operators_list),intent(in) :: operators
    type(block)                     :: self
    self%length    = length
    self%dim       = dim
    self%sectors   = sectors
    self%operators = operators
  end function build_block_from_scrath

  function build_block_from_site(ssite) result(self)
    type(site),intent(in) :: ssite
    type(block)           :: self
    self%length    = 1
    self%dim       = ssite%dim
    self%sectors   = ssite%sectors        
    self%operators = ssite%operators
  end function build_block_from_site


  !+------------------------------------------------------------------+
  !PURPOSE:  Destructor
  !+------------------------------------------------------------------+
  subroutine free_block(self)
    class(block) :: self
    self%length = 0
    self%dim    = 1
    call self%operators%free()
    call self%sectors%free()
  end subroutine free_block



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
  subroutine set_sectors_block(self,vec)
    class(block)         :: self
    real(8),dimension(:) :: vec
    self%sectors = sectors_list(vec)
  end subroutine set_sectors_block



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
    A%length = B%length
    A%dim    = B%dim
    A%operators  = B%operators
    A%sectors= B%sectors
  end subroutine block_equal_block


  function is_valid_block(self) result(bool)
    class(block) :: self
    logical      :: bool
    bool = self%operators%is_valid(self%dim).AND.(self%dim==len(self%sectors))
  end function is_valid_block





  !##################################################################
  !##################################################################
  !              SHOW 
  !##################################################################
  !##################################################################
  subroutine show_block(self,dble,fmt)
    class(block)              :: self
    logical,optional          :: dble
    character(len=*),optional :: fmt
    logical                   :: dble_
    character(len=32)         :: fmt_
    dble_=show_dble;if(present(dble))dble_=dble
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    write(*,*)"Block Length  =",self%length
    write(*,*)"Block Dim     =",self%dim
    write(*,*)"Block Sectors :"
    call self%sectors%show()
    write(*,*)"Site Operators:"
    call self%operators%show(dble=dble_,fmt=fmt_)
  end subroutine show_block









END MODULE BLOCKS





