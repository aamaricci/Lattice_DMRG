!
! A tuple is just an array q=[q_1,...,q_n].
! Yet we need to construct and operate with vectors of tuples.
! Here we construct specific objects to operate on these entities.
!
MODULE TUPLE_BASIS
  USE SCIFOR, only:str,sort_quicksort,assert_shape
  USE AUX_FUNCS
  implicit none
  private

  type tuple
     real(8),dimension(:),allocatable     :: qn
  end type tuple

  type tbasis
     type(tuple),dimension(:),allocatable :: basis
     integer                              :: qdim=0
     integer                              :: size=0
   contains
     procedure,pass      :: free     => free_tbasis     !destructor
     procedure,pass      :: dimq     => dimq_tbasis     !return qdim of the QN tuples
     procedure,pass      :: dump     => dump_tbasis     !dump basis to a rank-2 array [Nbasis,Qdim]
     procedure,pass      :: flat     => flat_tbasis     !flat basis to a rank-1 array
     procedure,pass      :: qn       => qn_tbasis     !flat basis to a rank-1 array
     procedure,pass      :: expand   => expand_tbasis   !expand the entire basis adding a new layer of qn
     procedure,pass      :: index    => index_tbasis    !create index from a given qn
     procedure,pass      :: show     => show_tbasis     !show
  end type tbasis


  interface tbasis
     module procedure :: construct_tbasis_I
     module procedure :: construct_tbasis_D
  end interface tbasis

  !RETURN SHAPE OF THE TBASIS
  intrinsic :: shape
  interface shape
     module procedure :: shape_tbasis
  end interface shape


  !ADDITION
  interface operator (+)
     module procedure :: sum_tbasis
  end interface operator (+)

  public :: tbasis
  public :: shape
  public :: operator(+)


contains


  !##################################################################
  !##################################################################
  !       DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Free a Tbasis (destructor) 
  !+------------------------------------------------------------------+
  subroutine free_tbasis(self)
    class(tbasis),intent(inout) :: self
    integer                     :: i
    if(.not.allocated(self%basis))return
    do i=1,size(self%basis)
       if(allocated(self%basis(i)%qn))deallocate(self%basis(i)%qn)
    enddo
    if(allocated(self%basis))deallocate(self%basis)
    self%qdim=0
    self%size=0
  end subroutine free_tbasis




  !##################################################################
  !##################################################################
  !       CONSTRUCT TVECTOR FROM ARRAY 
  !##################################################################
  !##################################################################
  function construct_tbasis_I(array,qdim) result(self)
    integer,dimension(:)               :: array
    integer                            :: Qdim
    type(tbasis),target                :: self
    integer                            :: Nbasis,N,i
    real(8),dimension(:,:),allocatable :: tvec
    call self%free()
    !
    N = size(array)
    if(mod(N,Qdim)/=0)stop "construct_tbasis_I error: size(array)%Qdim!=0"
    Nbasis = N/Qdim
    tvec = dble(transpose(reshape(array, shape=[Qdim,Nbasis])))
    allocate(self%basis(Nbasis))    
    do i=1,Nbasis
       allocate(self%basis(i)%qn, source=tvec(i,:))
    enddo
    self%qdim=Qdim
    self%size=Nbasis
  end function construct_tbasis_I

  function construct_tbasis_D(array,qdim) result(self)
    real(8),dimension(:)               :: array
    integer                            :: Qdim
    type(tbasis),target                :: self
    integer                            :: Nbasis,N,i
    real(8),dimension(:,:),allocatable :: tvec
    call self%free()
    !
    N = size(array)
    if(mod(N,Qdim)/=0)stop "construct_tbasis_I error: size(array)%Qdim!=0"
    Nbasis = N/Qdim
    tvec = dble(transpose(reshape(array, shape=[Qdim,Nbasis])))
    allocate(self%basis(Nbasis))
    do i=1,Nbasis
       allocate(self%basis(i)%qn, source=tvec(i,:))
    enddo
    self%qdim=Qdim
    self%size=Nbasis
  end function construct_tbasis_D


  !##################################################################
  !##################################################################
  !       EXPAND THE BASIS VERTICALLY ADDING A NEW LAYER
  !       OF QN GIVEN A SUITABLE ARRAY OF STATES (ie LIST OF QNs)
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: expand a tbasis set by set, ie add a layer of states for a give
  ! value of one of the quantum numbers. 
  ! if basis is not allocated: create a new instance using the array and qdim=1
  ! else: append to each qn the new component from the array
  !+------------------------------------------------------------------+
  subroutine expand_tbasis(self,array)
    class(tbasis),intent(inout)        :: self
    real(8),dimension(:)               :: array
    integer                            :: Nbasis,N,i,Qdim
    real(8),dimension(:,:),allocatable :: tvec
    !
    N = size(array)
    !
    if(.not.allocated(self%basis))then
       Nbasis = N
       tvec = dble(transpose(reshape(array, shape=[1,Nbasis])))
       allocate(self%basis(Nbasis))
       do i=1,Nbasis
          allocate(self%basis(i)%qn, source=tvec(i,:))
       enddo
       self%qdim=1
       self%size=Nbasis
       return
    endif
    !
    if(N/=size(self%basis))stop "grow_tbasis error: size(array) != size(self.basis)"
    do i=1,N
       call append(self%basis(i)%qn,array(i))
    enddo
    self%qdim=self%qdim+1
  end subroutine expand_tbasis




  !##################################################################
  !##################################################################
  !       DUMP TO A RANK-2 ARRAY/ FLAT to a RANK-1 ARRAY
  !##################################################################
  !##################################################################
  function dump_tbasis(self) result(tvec)
    class(tbasis),intent(inout)        :: self
    real(8),dimension(:,:),allocatable :: tvec
    integer                            :: i,Nbasis,Qdim
    if(allocated(tvec))deallocate(tvec)
    Nbasis = size(self%basis)
    Qdim   = self%qdim
    allocate(tvec(Nbasis,Qdim))
    do i=1,Nbasis
       tvec(i,:) = self%basis(i)%qn
    enddo
  end function dump_tbasis



  function flat_tbasis(self) result(array)
    class(tbasis),intent(inout)        :: self
    real(8),dimension(:),allocatable   :: array
    real(8),dimension(:,:),allocatable :: tvec
    integer                            :: i,Nbasis,Qdim
    Nbasis = size(self%basis)
    Qdim   = self%qdim
    allocate(tvec(Nbasis,Qdim))
    do i=1,Nbasis
       tvec(i,:) = self%basis(i)%qn
    enddo
    if(allocated(array))deallocate(array)
    array = pack(transpose(tvec),.true.)
  end function flat_tbasis


  function qn_tbasis(self,index) result(qn)
    class(tbasis),intent(inout)      :: self
    real(8),dimension(:),allocatable :: qn
    integer                          :: index
    integer                          :: index_
    integer                          :: i,Nbasis,Qdim
    index_=index
    if(index_>self%size.OR.index_<=0)stop "qn_tbasis: index !in [1,self.size]"
    if(allocated(qn))deallocate(qn)
    Nbasis = size(self%basis)
    Qdim   = self%qdim
    allocate(qn, source=self%basis(index_)%qn)
  end function qn_tbasis


  !##################################################################
  !##################################################################
  !       RETURN QDIM
  !##################################################################
  !##################################################################  
  function dimq_tbasis(self) result(qdim)
    class(tbasis),intent(inout)        :: self
    integer                            :: Qdim
    qdim = self%qdim
  end function dimq_tbasis





  !##################################################################
  !##################################################################
  !              INDEX THE BASIS GIVEN A QN TUPLE
  !##################################################################
  !##################################################################
  function index_tbasis(self,qn) result(index)
    class(tbasis)                    :: self
    real(8),dimension(self%qdim)     :: qn
    integer,dimension(:),allocatable :: index
    logical,dimension(:),allocatable  :: mask
    integer                           :: i,N,pos,Nbasis
    !
    if(allocated(index))deallocate(index)
    allocate(mask(self%size))
    forall(i=1:size(mask))mask(i) = all(self%basis(i)%qn==qn)
    N   = count(mask)
    pos = 0
    do i=1,N
       pos = pos+findloc(mask(pos+1:),value=.true.,dim=1)
       call append(index,pos)
    enddo
    call sort_quicksort(index)
  end function index_tbasis






  !+------------------------------------------------------------------+
  !PURPOSE:  Return the shape of a sparse matrix
  !+------------------------------------------------------------------+
  function shape_tbasis(self) result(shape)
    class(tbasis),intent(in) :: self
    integer,dimension(2)     :: shape
    integer                  :: Nbasis,Qdim
    Qdim  = self%qdim
    Nbasis = self%size
    shape = [Nbasis,Qdim]
  end function shape_tbasis


 


  !##################################################################
  !##################################################################
  !              OPERATIONS
  !##################################################################
  !##################################################################
  function sum_tbasis(a,b) result(c)
    type(tbasis), intent(in)  :: a,b
    type(tbasis)              :: c
    integer                   :: N,i,Qdim
    real(8),dimension(a%qdim) :: qn
    if(a%size/=b%size)stop "sum_tbasis error: size(a.basis) != size(b.basis)"
    if(a%qdim/=b%qdim)stop "sum_tbasis error: a.Qdim != b.Qdim"
    call c%free()
    N=a%size
    Qdim=a%qdim
    c%qdim=Qdim
    c%size=N
    allocate(c%basis(N))
    do i=1,N
       qn =  a%basis(i)%qn + b%basis(i)%qn
       allocate(c%basis(i)%qn, source=qn)
    enddo
  end function sum_tbasis





  !##################################################################
  !##################################################################
  !              SHOW TVECTOR
  !##################################################################
  !##################################################################
  subroutine show_tbasis(self)
    class(tbasis) :: self
    integer       :: Nqn,Qdim,i,j
    !
    do i=1,size(self%basis)
       write(*,"(I6,2x)",advance='no')i
       write(*,"(A1)",advance='no')"["
       write(*,"(F6.2,A1)",advance='no')self%basis(i)%qn(1)
       do j=2,self%Qdim
          write(*,"(A1,F6.2)",advance='no')",",self%basis(i)%qn(j)
       enddo
       write(*,"(A1)",advance='no')"]"
       write(*,*)""
    enddo
    write(*,*)""
  end subroutine show_tbasis



END MODULE TUPLE_BASIS
