MODULE MATRIX_BLOCKS
  USE SCIFOR, only: str,free_unit,zero,assert_shape,zeye,eigh,sort_quicksort
  USE AUX_FUNCS
  implicit none
  private

  !BLOCK MATRIX COMPONENT
  type block_type
     integer                            :: index=0
     real(8),dimension(:),allocatable   :: qn
     real(8),dimension(:),allocatable   :: E
     integer,dimension(:),allocatable   :: map
     real(8),dimension(:,:),allocatable :: M
     type(block_type),pointer           :: next=>null()
  end type block_type


  !BLOCK MATRIX STRUCTURE
  type blocks_matrix
     integer                          :: Nblock=0
     integer                          :: Ndim=0
     real(8),dimension(:),allocatable :: evalues
     integer,dimension(:),allocatable :: eorder
     logical                          :: diag=.false.
     type(block_type),pointer         :: root=>null()
   contains
     procedure,pass :: free   => free_blocks_matrix
     procedure,pass :: load   => load_blocks_matrix
     procedure,pass :: append => append_blocks_matrix
     procedure,pass :: push   => push_blocks_matrix
     procedure,pass :: block  => get_block_blocks_matrix
     procedure,pass :: qn     => get_qn_blocks_matrix
     procedure,pass :: map    => get_map_blocks_matrix
     procedure,pass :: index  => get_index_blocks_matrix
     procedure,pass :: evals  => evals_blocks_matrix
     procedure,pass :: eval   => evals_blocks_matrix
     procedure,pass :: evec   => evec_blocks_matrix
     procedure,pass :: has_qn => has_qn_blocks_matrix
     procedure,pass :: dump   => dump_blocks_matrix
     procedure,pass :: dgr    => dgr_blocks_matrix
     procedure,pass :: show   => show_blocks_matrix
     procedure,pass :: shape  => shape_block_blocks_matrix
     procedure,pass :: size   => size_block_blocks_matrix
     procedure,pass :: dims   => dimensions_blocks_matrix
     procedure,pass :: eigh   => eigh_blocks_matrix
     procedure,pass :: find   => find_indices_blocks_matrix
  end type blocks_matrix



  interface blocks_matrix
     module procedure :: construct_blocks_matrix
  end interface blocks_matrix

  interface as_blocks
     module procedure :: construct_blocks_matrix
  end interface as_blocks


  !EQUALITY with scalar and function (A=B)
  interface assignment(=)
     module procedure :: blocks_matrix_equal_blocks_matrix
  end interface assignment(=)



  !INTRINSIC FUNCTION SIZE
  intrinsic :: size
  interface size
     module procedure :: size_blocks_matrix
  end interface size


  !RETURN SHAPE OF THE WHOLE BLOCK MATRIX [Ndim,Ndim]
  intrinsic :: shape
  interface shape
     module procedure :: shape_blocks_matrix
  end interface shape

  intrinsic :: transpose
  interface transpose
     module procedure :: transpose_blocks_matrix
  end interface transpose

  interface hconjg
     module procedure :: transpose_blocks_matrix
  end interface hconjg

  interface as_matrix
     module procedure :: as_matrix_blocks_matrix
  end interface as_matrix

  public :: blocks_matrix
  public :: as_blocks
  public :: size
  public :: shape
  public :: transpose
  public :: hconjg
  public :: assignment(=)



contains       


  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function construct_blocks_matrix(matrix,qn,map) result(self)
    real(8),dimension(:,:),intent(in)           :: matrix
    real(8),dimension(:),allocatable,intent(in) :: qn
    integer,dimension(:)                        :: map
    type(blocks_matrix)                         :: self
    call self%free()
    allocate(self%root)
    call self%append(matrix,qn,map)
  end function construct_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: load a block matrix into a block one
  !+------------------------------------------------------------------+
  subroutine load_blocks_matrix(self,matrix,qn,map)
    class(blocks_matrix),intent(inout) :: self
    real(8),dimension(:),allocatable    :: qn
    integer,dimension(:)               :: map
    real(8),dimension(:,:),intent(in)  :: matrix
    call self%free()
    allocate(self%root)
    call self%append(matrix,qn,map)
  end subroutine load_blocks_matrix



  !+------------------------------------------------------------------+
  !PURPOSE: free an entire block matrix
  !+------------------------------------------------------------------+
  subroutine free_blocks_matrix(self)    
    class(blocks_matrix),intent(inout) :: self
    type(block_type),pointer          :: p,c
    if(.not.associated(self%root))return
    do
       p => self%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next
       c%next => null()
       c%index=  0
       if(allocated(c%qn))deallocate(c%qn)
       if(allocated(c%M))deallocate(c%M)
       if(allocated(c%map))deallocate(c%map)
       deallocate(c)
    enddo
    if(allocated(self%evalues))deallocate(self%evalues)
    if(allocated(self%eorder))deallocate(self%eorder)
    self%Nblock=0
    self%Ndim=0
    self%diag=.false.
    self%root=>null()
    p=>null()
    c=>null()
  end subroutine free_blocks_matrix






  !##################################################################
  !##################################################################
  !       APPEND/PUSH  a block to the list
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Append a block
  !+------------------------------------------------------------------+
  subroutine append_blocks_matrix(self,matrix,qn,map)
    class(blocks_matrix),intent(inout)          :: self
    real(8),dimension(:,:),intent(in)           :: matrix
    real(8),dimension(:),allocatable,intent(in) :: qn
    integer,dimension(:)                        :: map
    integer                                     :: Dim
    type(block_type),pointer                    :: p,c
    !    
    if(.not.associated(self%root))allocate(self%root)
    !
    Dim = size(matrix,1);call assert_shape(matrix,[Dim,Dim],"append_block","matrix")
    !
    p => self%root
    c => p%next
    do
       if(.not.associated(c))exit
       p => c
       c => c%next
    end do
    !
    allocate(p%next)
    allocate(p%next%M, source=matrix)
    ! p%next%M     = matrix
    allocate(p%next%map, source=map)
    ! p%next%map   = map
    allocate(p%next%qn, source=qn)
    ! p%next%qn    = qn
    p%next%index = p%index+1
    !
    if(.not.associated(c))then !end of the list special case (c=>c%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    self%Nblock = self%Nblock+1
    !
    self%Ndim = self%Ndim+size(Matrix,1)
    !
    p=>null()
    c=>null()
  end subroutine append_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE:  PUSH a block
  !+------------------------------------------------------------------+
  subroutine push_blocks_matrix(self,matrix,qn,map)
    class(blocks_matrix),intent(inout)          :: self
    real(8),dimension(:,:),intent(in)           :: matrix
    real(8),dimension(:),allocatable,intent(in) :: qn
    integer,dimension(:)                        :: map
    logical                                     :: iupdate
    type(block_type),pointer                    :: p,c
    integer                                     :: Dim
    !    
    if(.not.associated(self%root))allocate(self%root)
    !
    Dim = size(matrix,1);call assert_shape(matrix,[Dim,Dim],"append_block","matrix")
    !
    if(.not.self%has_qn(qn))then
       call self%append(matrix,qn,map)
       return
    endif
    !
    p => self%root
    c => p%next
    do
       if(.not.associated(c))exit
       if(all(c%qn == qn))then
          exit
       endif
       p => c
       c => c%next
    end do
    !
    self%Ndim = self%Ndim-size(c%M,1)
    if(allocated(c%M))deallocate(c%M)
    allocate(c%M, source=Matrix)
    ! c%M   = Matrix
    if(allocated(c%map))deallocate(c%map)
    allocate(c%map, source=map)
    ! c%map = map
    allocate(c%qn, source=qn)
    ! c%qn  = qn
    self%Ndim = self%Ndim+size(Matrix,1)
    !
    p=>null()
    c=>null()
  end subroutine push_blocks_matrix





  !##################################################################
  !##################################################################
  !              RETRIEVE CONTENT: GET BLOCK, DUMP WHOLE MATRIX
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: get 
  !+------------------------------------------------------------------+
  function get_block_blocks_matrix(self,index,m) result(matrix)
    class(blocks_matrix)               :: self
    integer,optional                   :: index,m
    real(8),dimension(:,:),allocatable :: matrix
    integer                            :: index_,m_
    logical                            :: ifound
    type(block_type),pointer           :: c
    !    
    index_=1
    if(present(index))then
       index_=index
    elseif(present(m))then
       m_ = m
       if(self%diag)m_=self%eorder(m)
       call self%find(m_,index_)
    else
       stop "get_block_blocks_matrix error: !present(index) + !present(m)"
    endif
    if(index_>self%Nblock.OR.index_<=0)stop "get_block_blocks_matrix error: block_index !in [1,self.size]"
    !
    ifound=.false.
    c => self%root%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if(c%index == index_) then
          ifound=.true.
          exit          
       endif
       c => c%next
    end do
    if(.not.ifound)stop "get_blocks_matrix error: not found"
    !
    if(allocated(matrix))deallocate(matrix)
    allocate(matrix, source=c%M)
    ! matrix = c%M
    !
    c=>null()
  end function get_block_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: get 
  !+------------------------------------------------------------------+
  function get_qn_blocks_matrix(self,index,m) result(qn)
    class(blocks_matrix)             :: self
    integer,optional                 :: index
    integer,optional                 :: m
    integer                          :: index_,q,m_
    real(8),dimension(:),allocatable :: qn
    logical                          :: ifound
    type(block_type),pointer         :: c
    !
    if(allocated(qn))deallocate(qn)
    !
    index_=1
    if(present(index))then
       index_=index
    elseif(present(m))then
       m_ = m
       if(self%diag)m_=self%eorder(m)
       call self%find(m_,index_)
    else
       stop "get_qn_blocks_matrix error: !present(index) + !present(m)"
    endif
    if(index_>self%Nblock.OR.index_<=0)stop "get_qn_blocks_matrix error: block_index !in [1,self.size]"
    ifound=.false.
    c => self%root%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if(c%index == index_) then
          ifound=.true.
          exit          
       endif
       c => c%next
    end do
    if(.not.ifound)stop "get_qn_matrix error: not found"
    !
    allocate(qn, source=c%qn)
    ! qn = c%qn
    !
    c=>null()
  end function get_qn_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: get 
  !+------------------------------------------------------------------+
  function get_map_blocks_matrix(self,index,m) result(map)
    class(blocks_matrix)             :: self
    integer,optional                 :: index
    integer,optional                 :: m
    integer                          :: index_,q,m_
    integer,dimension(:),allocatable :: map
    logical                          :: ifound
    type(block_type),pointer         :: c
    !    
    index_=1
    if(present(index))then
       index_=index
    elseif(present(m))then
       m_ = m
       if(self%diag)m_=self%eorder(m)
       call self%find(m_,index_)
    else
       stop "get_qn_blocks_matrix error: !present(index) + !present(m)"
    endif
    if(index_>self%Nblock.OR.index_<=0)stop "get_qn_blocks_matrix error: block_index !in [1,self.size]"
    ifound=.false.
    c => self%root%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if(c%index == index_) then
          ifound=.true.
          exit          
       endif
       c => c%next
    end do
    if(.not.ifound)stop "get_qn_matrix error: not found"
    !
    if(allocated(map))deallocate(map)
    allocate(map, mold=c%map)
    map = c%map
    !
    c=>null()
  end function get_map_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: get the element value at position (i,j) in the block matrix
  !+------------------------------------------------------------------+
  function get_index_blocks_matrix(self,qn) result(index)
    class(blocks_matrix)                        :: self
    real(8),dimension(:),allocatable,intent(in) :: qn
    integer                                     :: index
    logical                                     :: ifound
    type(block_type),pointer                    :: c
    !
    !qn does not exist: return with index=0
    index=0
    if(.not.self%has_qn(qn))return
    !
    !else qn exists so we can find it
    c => self%root%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if(all(c%qn == qn)) then
          exit          
       endif
       c => c%next
    end do
    !
    index = c%index
    !
    c=>null()
  end function get_index_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: dump a block matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  function dump_blocks_matrix(self) result(matrix)
    class(blocks_matrix),intent(inout) :: self
    real(8),dimension(:,:),allocatable :: matrix,mtmp
    integer                            :: Offset
    integer                            :: N,M
    type(block_type),pointer           :: c
    if(allocated(matrix))deallocate(matrix)
    allocate(matrix(self%Ndim,self%Ndim));matrix=0d0
    Offset=0
    c => self%root%next
    do
       if(.not.associated(c))exit
       N = size(c%M,1)
       Matrix(Offset+1:Offset+N,Offset+1:Offset+N) = c%M
       Offset = Offset + N
       c => c%next
    end do
    c=>null()
    !
  end function dump_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  function as_matrix_blocks_matrix(self) result(matrix)
    class(blocks_matrix),intent(inout) :: self
    real(8),dimension(:,:),allocatable :: matrix
    if(allocated(matrix))deallocate(matrix)
    matrix = self%dump()
  end function as_matrix_blocks_matrix




  !##################################################################
  !##################################################################
  !              LINEAR ALGEBRA 
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: performs eigh on each block of the blocks_matrix,
  ! the EigenVectors matrix overwrites the block and the EigenValues
  ! are stored in the array E of the block (so far unused). 
  !+------------------------------------------------------------------+
  subroutine eigh_blocks_matrix(self,sort,reverse,order)
    class(blocks_matrix),intent(inout)        :: self
    logical,optional                          :: sort,reverse
    integer,dimension(:),allocatable,optional :: order
    logical                                   :: sort_,reverse_
    type(block_type),pointer                  :: c
    integer                                   :: i,Nloc,Offset,N
    !
    sort_   =.true.;if(present(sort))sort_=sort
    reverse_=.true.;if(present(reverse))reverse_=reverse
    !
    N = self%Ndim
    if(allocated(self%evalues))deallocate(self%evalues)
    if(allocated(self%eorder))deallocate(self%eorder)
    allocate(self%evalues(N));self%evalues=0d0
    allocate(self%eorder(N));self%eorder=0
    !
    Offset=0
    c => self%root%next
    do                          !loop over all blocks
       if(.not.associated(c))exit
       Nloc = size(c%M,1)
       if(allocated(c%E))deallocate(c%E)
       allocate(c%E(Nloc))
       call eigh(c%M,c%E)  !<- overwrites blocks with eigenvec matrix
       !
       self%evalues(Offset+1:Offset+Nloc) = c%E
       Offset = Offset + Nloc
       !
       c => c%next
    end do
    c=>null()
    !
    !init Eorder vector to 1..N
    self%eorder=(/(i,i=1,N)/)
    if(sort_)then
       call sort_quicksort(self%evalues,self%eorder)
       if(reverse_)then
          self%evalues = self%evalues(N:1:-1)
          self%eorder  = self%eorder(N:1:-1)
       endif
    endif
    if(present(order))then
       if(allocated(order))deallocate(order)
       allocate(order, mold=self%eorder)
       order = self%eorder
    endif
    !
    self%diag=.true.
    !
  end subroutine eigh_blocks_matrix



  !+------------------------------------------------------------------+
  !PURPOSE: returns an array with all evals of the block matrix 
  !+------------------------------------------------------------------+
  function evals_blocks_matrix(self) result(evals)
    class(blocks_matrix)             :: self
    real(8),dimension(:),allocatable :: evals
    type(block_type),pointer         :: c
    !
    if(.not.self%diag)stop "Evals_blocks_matrix error: self.diag=F, call self.eigh() before trying to get an eigenvector"
    !
    if(allocated(evals))deallocate(evals)
    !
    allocate(evals, source=self%evalues)
    !
  end function evals_blocks_matrix



  !+------------------------------------------------------------------+
  !PURPOSE: returns a single eigenvalue with all evals of the block matrix 
  !+------------------------------------------------------------------+
  function eval_blocks_matrix(self,m) result(eval)
    class(blocks_matrix)             :: self
    integer                          :: m
    real(8)                          :: eval
    if(.not.self%diag)stop "Evals_blocks_matrix error: self.diag=F, call self.eigh() before trying to get an eigenvector"
    if(m>self%Ndim.OR.m<=0)stop "evals_block_matrix warning: m !in [1,self.Ndim]"
    eval = self%evalues(m)
  end function eval_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: 
  !+------------------------------------------------------------------+
  function evec_blocks_matrix(self,m) result(vec)
    class(blocks_matrix)             :: self
    integer                          :: m
    real(8),dimension(:),allocatable :: vec
    integer                          :: m_,i,q,pos
    type(block_type),pointer         :: c
    !
    if(.not.self%diag)stop "Evec_blocks_matrix error: self.diag=F, call self.eigh() before trying to get an eigenvector"
    !
    if(m>self%Ndim.OR.m<=0)stop "evals_block_matrix warning: m !in [1,self.Ndim]"
    m_=self%eorder(m)
    call self%find(m_,q,pos)
    !
    c => self%root
    do i=1,q
       c => c%next
    enddo
    !
    if(allocated(vec))deallocate(vec)
    allocate(vec, source=c%M(:,pos))
    !
    c=>null()
  end function evec_blocks_matrix






  !##################################################################
  !##################################################################
  !              OVERLOAD INTRINSIC and OPERATIONS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Returns the size of given qmap
  !+------------------------------------------------------------------+
  function size_blocks_matrix(self) result(size)
    class(blocks_matrix),intent(in) :: self
    integer                        :: size
    size = self%Nblock
  end function size_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE:  Return the shape of a sparse matrix
  !+------------------------------------------------------------------+
  function shape_blocks_matrix(self) result(shape)
    class(blocks_matrix),intent(in) :: self
    integer,dimension(2)           :: shape
    shape = [self%Ndim,self%Ndim]
  end function shape_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: 
  !+------------------------------------------------------------------+
  function shape_block_blocks_matrix(self,index) result(Bshape)
    class(blocks_matrix)                  :: self
    integer,optional                      :: index
    integer                               :: index_
    integer,dimension(2)                  :: Bshape
    logical                               :: ifound
    type(block_type),pointer              :: c
    !
    if(.not.present(index))then
       Bshape = shape(self)
       return
    endif
    !
    index_=1;if(present(index))index_=index
    if(index_>self%Nblock.OR.index_<=0)stop "shape_block_blocks_matrix error: block_index !in [1,self.size]"
    !
    ifound=.false.
    c => self%root%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if(c%index == index_) then
          ifound=.true.
          exit          
       endif
       c => c%next
    end do
    if(.not.ifound)stop "get_blocks_matrix error: not found"
    !
    Bshape = shape(c%M)
    !
    c=>null()
  end function shape_block_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: 
  !+------------------------------------------------------------------+
  function size_block_blocks_matrix(self,index,dim) result(Bsize)
    class(blocks_matrix)     :: self
    integer,optional         :: index
    integer,optional         :: dim
    integer                  :: index_
    integer                  :: dim_
    integer                  :: Bsize
    integer                  :: Mshape(2)
    logical                  :: ifound
    type(block_type),pointer :: c
    !
    dim_=1;if(present(dim))dim_=dim
    if(dim_ < 1 .OR. dim_ > 2)then
       write(0,*)"get warning: dim_ != [1,2]: reset to dim=1"
       dim_=1
    endif
    !
    index_=1;if(present(index))index_=index
    if(index_>self%Nblock.OR.index_<=0)stop "shape_block_blocks_matrix error: block_index !in [1,self.size]"
    !
    if(.not.present(index))then !return the dimension of the whole self.
       Mshape= shape(self)
       Bsize = Mshape(dim_)
       return
    endif
    !
    ifound=.false.
    c => self%root%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if(c%index == index_) then
          ifound=.true.
          exit          
       endif
       c => c%next
    end do
    if(.not.ifound)stop "get_blocks_matrix error: not found"
    !
    Bsize = size(c%M,dim_)
    !
    c=>null()
  end function size_block_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE:  Returns dimensions D_q=1,size(self) of the Block matrix
  !+------------------------------------------------------------------+
  function dimensions_blocks_matrix(self,dim) result(Dvec)
    class(blocks_matrix),intent(inout) :: self
    integer,dimension(:),allocatable   :: Dvec
    integer,optional                   :: dim
    integer                            :: dim_
    type(block_type),pointer           :: c
    integer                            :: q
    dim_=1;if(present(dim))dim_=dim
    if(dim_ < 1 .OR. dim_ > 2)then
       write(0,*)"get warning: dim_ != [1,2]: reset to dim=1"
       dim_=1
    endif
    if(allocated(Dvec))deallocate(Dvec)
    allocate(Dvec(size(self)))
    do q=1,size(self)
       Dvec(q) = self%size(index=q,dim=dim_)
    enddo
  end function dimensions_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE:  Given an integer index m=1,self.Ndim,
  ! returns the corresponding BLock index containing it and the relative
  ! position pos inside the block.
  ! [1...D_1][D_1+1...D_2]...[D_{q-1}...D_q] = D=self.Ndim
  !+------------------------------------------------------------------+  
  subroutine find_indices_blocks_matrix(self,m,index,pos)
    class(blocks_matrix)             :: self
    integer                          :: m
    integer                          :: index
    integer                          :: q
    integer,optional                 :: pos
    integer                          :: pos_
    integer,dimension(:),allocatable :: Dq
    Dq = self%dims(dim=2)
    if(m>sum(Dq).OR.m<=0)stop "find_indices_blocks_matrix error: m !in [1,D==self.Ndim]"
    do q=1,size(Dq)
       if(m <= sum(Dq(1:q)))then
          index = q
          pos_  = m - sum(Dq(1:q-1))          
          exit
       endif
    enddo
    if(present(pos))pos=pos_
  end subroutine find_indices_blocks_matrix

  !+------------------------------------------------------------------+
  !PURPOSE:  Returns True is qn exists, False otherwise
  !+------------------------------------------------------------------+
  function has_qn_blocks_matrix(self, qn) result(bool)
    class(blocks_matrix),intent(inout)          :: self
    real(8),dimension(:),allocatable,intent(in) :: qn
    type(block_type),pointer                    :: c
    logical                                     :: bool
    bool=.false.
    c => self%root%next
    do                            !traverse the list until index is found
       if(.not.associated(c))exit
       if(all(c%qn == qn)) then
          bool=.true.
          exit
       endif
       c => c%next
    end do
  end function has_qn_blocks_matrix




  !##################################################################
  !##################################################################
  !              SHOW  
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: show matrix
  !+------------------------------------------------------------------+
  subroutine show_blocks_matrix(self,dble,fmt,file)
    class(blocks_matrix),intent(in) :: self
    character(len=*),optional       :: fmt
    logical,optional                :: dble
    character(len=*),optional       :: file
    character(len=12)               :: fmt_
    logical                         :: dble_
    integer                         :: i,j,unit_
    character(len=64)               :: format
    real(8)                         :: val
    type(block_type),pointer        :: c
    !
    unit_=6
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    dble_=show_dble;if(present(dble))dble_=dble
    if(present(file))open(free_unit(unit_),file=str(file))
    !
    if(dble_)then
       format='('//str(fmt_)//',1x)'
    else       
       format='(A1,'//str(fmt_)//',A1,'//str(fmt_)//',A1,1x)'
    endif
    !
    write(*,"(A6,I12)")"Size :",size(self)
    write(*,"(A6,2I6)")"Shape:",shape(self)
    write(*,"(A18)")"------------------"
    c => self%root%next
    do
       if(.not.associated(c))exit
       write(unit_,"(A6,I12)")"Index:",c%index
       write(unit_,"(A6,"//str(size(c%map))//"I6)")"map  :",c%map
       write(unit_,"(A6,"//str(size(c%qn))//"F12.4)")"QN   :",c%qn
       if(allocated(c%E))&
            write(unit_,"(A6,"//str(size(c%E))//"F12.4)")"E    :",c%E
       write(unit_,"(A6)")"Block:"
       do i=1,size(c%M,1)
          do j=1,size(c%M,2)
             val = c%M(i,j)
             ! if(dble_)then
             write(unit_,"("//str(self%Ndim)//str(format)//")",advance='no')val
             ! else       
             !    write(unit_,"("//str(self%Ndim)//str(format)//")",advance='no')&
             !         "(",dreal(val),",",dimag(val),")"
             ! endif
          enddo
          write(unit_,*)
       enddo
       write(unit_,*)
       if(present(file))close(unit_)
       c => c%next
    enddo
  end subroutine show_blocks_matrix





  !##################################################################
  !##################################################################
  !               BLOCK MATRIX BASIC ALGEBRA 
  !##################################################################
  !##################################################################
  function dgr_blocks_matrix(a) result(adg)
    class(blocks_matrix), intent(in) :: a
    type(blocks_matrix)              :: adg
    integer                               :: i    
    call adg%free()
    do i=1,size(a)
       call adg%append((transpose(a%block(index=i))), a%qn(index=i), a%map(index=i))
    enddo
  end function dgr_blocks_matrix

  function transpose_blocks_matrix(a) result(adg)
    class(blocks_matrix), intent(in) :: a
    type(blocks_matrix)              :: adg
    integer                               :: i    
    call adg%free()
    do i=1,size(a)
       call adg%append((transpose(a%block(index=i))), a%qn(index=i), a%map(index=i))
    enddo
  end function transpose_blocks_matrix


  !+------------------------------------------------------------------+
  !PURPOSE:  Block matrix equality A = B. Deep copy
  !+------------------------------------------------------------------+
  subroutine blocks_matrix_equal_blocks_matrix(a,b)
    type(blocks_matrix),intent(inout)      :: a
    type(blocks_matrix),intent(in)         :: b
    integer                               :: i    
    call a%free()
    do i=1,size(b)
       call a%append(b%block(index=i), b%qn(index=i), b%map(index=i))
    enddo
  end subroutine blocks_matrix_equal_blocks_matrix














end module MATRIX_BLOCKS






! !+------------------------------------------------------------------+
!  !PURPOSE: get the element value at position (i,j) in the block matrix
!  !+------------------------------------------------------------------+
!  function evec_blocks_matrix(self,index,pos) result(vec)
!    class(blocks_matrix)             :: self
!    integer                          :: index
!    integer                          :: pos
!    real(8),dimension(:),allocatable :: vec
!    integer                          :: index_
!    logical                          :: ifound
!    type(block_type),pointer         :: c
!    !    
!    index_=index
!    if(index_>self%Nblock)then
!       write(0,*)"get warning: index > self.size: reset to index=self.size"
!       index_=self%Nblock
!    endif
!    if(index_<=0)then
!       write(0,*)"get warning: index <= 0: reset to index=1"
!       index_=1
!    endif
!    !
!    if(allocated(vec))deallocate(vec)
!    !
!    ifound=.false.
!    c => self%root%next
!    do                            !traverse the list until QN is found
!       if(.not.associated(c))exit
!       if(c%index == index_) then
!          ifound=.true.
!          exit          
!       endif
!       c => c%next
!    end do
!    if(.not.ifound)stop "get_vec error: not found"
!    !
!    if((pos < 0).OR.(pos > size(c%M,2)))stop "get_vec error: pos<0 or pos>size(M,2) "
!    allocate(vec, source=c%M(:,pos))
!    !
!    c=>null()
!  end function evec_blocks_matrix
