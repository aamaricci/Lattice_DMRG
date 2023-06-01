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
     integer                                     :: length=0
     integer                                     :: dim=1
     type(sectors_list),dimension(:),allocatable :: sectors
     type(operators_list)                        :: operators
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
     module procedure :: constructor_from_scrath
     module procedure :: constructor_from_site
  end interface block

  !GENERIC CONSTRUCTOR
  interface as_block
     module procedure :: constructor_from_scrath
     module procedure :: constructor_from_site
  end interface as_block

  !EQUALITY 
  interface assignment(=)
     module procedure :: block_equal_block
  end interface assignment(=)

  public :: block
  public :: as_block
  public :: assignment(=)


  integer :: i,j

contains


  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function constructor_from_scrath(length,dim,sectors,operators) result(self)
    integer,intent(in)              :: length
    integer,intent(in)              :: dim
    type(sectors_list),intent(in)   :: sectors(:)
    type(operators_list),intent(in) :: operators
    type(block)                     :: self
    self%length    = length
    self%dim       = dim
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
    self%dim       = ssite%dim
    self%operators = ssite%operators
    allocate(self%sectors(size(ssite%sectors)))
    do i=1,size(self%sectors)
       self%sectors(i)   = ssite%sectors(i)
    enddo
  end function constructor_from_site


  !+------------------------------------------------------------------+
  !PURPOSE:  Destructor
  !+------------------------------------------------------------------+
  subroutine free_block(self)
    class(block) :: self
    self%length = 0
    self%dim    = 1
    call self%operators%free()
    if(allocated(self%sectors))then
       call self%sectors%free()
       deallocate(self%sectors)
    endif
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
  subroutine set_sectors_block(self,indx,vec)
    class(block)         :: self
    integer              :: indx
    real(8),dimension(:) :: vec
    if(indx<1.OR.indx>size(self%sectors))stop "SET_SECTORS_BLOCK ERROR: indx out of range"
    self%sectors(indx) = sectors_list(vec)
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
    allocate(A%sectors(size(B%sectors)))
    do i=1,size(A%sectors)
       A%sectors(i)   = B%sectors(i)
    enddo
  end subroutine block_equal_block


  function is_valid_block(self) result(bool)
    class(block) :: self
    logical      :: bool
    integer,dimension(size(self%sectors)) :: Lvec
    bool = self%operators%is_valid(self%dim)
    do i=1,size(self%sectors)
       Lvec = len(self%sectors(i))
    enddo
    bool=bool.AND.(self%dim==product(Lvec))
  end function is_valid_block


  

  !##################################################################
  !##################################################################
  !              SHOW 
  !##################################################################
  !##################################################################
  subroutine show_block(self,fmt)
    class(block)              :: self
    character(len=*),optional :: fmt
    character(len=32)         :: fmt_
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    write(*,*)"Block Length  =",self%length
    write(*,*)"Block Dim     =",self%dim
    do i=1,size(self%sectors)
       write(*,*)"Block Sectors: ",i
       call self%sectors(i)%show()
    enddo
    write(*,*)"Block Operators:"
    call self%operators%show(fmt=fmt_)
  end subroutine show_block









END MODULE BLOCKS





