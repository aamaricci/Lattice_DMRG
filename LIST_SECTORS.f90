! We assume quantum numbers are given in terms of a tuple Q=[q_1,...,q_M],
! e.g. M=2 with q_1=n_up, q_2=n_dw
! Then we employ a TUPLE_BASIS, ie a vector of tuples as V=[Q_1,...,Q_N]^T.
! For each i=1,...,N we have a tuple Q_i.
! Here we create and index, ie MAP, the sectors based on such tuples of QN.
! ie we return the state indices (as positions in the tuple_basis) with
! a given value of Q.
MODULE LIST_SECTORS
  USE SCIFOR, only:str,sort_quicksort,assert_shape
  USE AUX_FUNCS
  USE TUPLE_BASIS
  implicit none
  private


  type qtype
     integer                          :: index=0
     real(8),dimension(:),allocatable :: qn           !tuple of quantum numbers
     integer,dimension(:),allocatable :: map          !map of the tuples, the states with a given tuple of qn.
     type(qtype),pointer              :: next=>null()  !link to next box (chain)
  end type qtype


  type sectors_list
     integer             :: qdim=0
     integer             :: size=0
     type(qtype),pointer :: root=>null()
   contains
     procedure,pass      :: free     => free_sectors_list     !destructor
     procedure,pass      :: put      => put_sectors_list      !put sectors_list qn array
     procedure,pass      :: load     => load_sectors_list     !load=sequential put=constructor
     procedure,pass      :: append   => append_sectors_list   !append map state given a QN
     procedure,pass      :: get      => get_sectors_list      !get qn and map for a given index
     procedure,pass      :: map      => get_map_sectors_list  !get map for a given qn/index
     procedure,pass      :: qn       => get_qn_sectors_list   !return qn for a given index
     procedure,pass      :: basis    => basis_sectors_list    !return the basis of the sector list
     procedure,pass      :: has_qn   => has_qn_sectors_list   !True if qn exists
     procedure,pass      :: show     => show_sectors_list     !show sectors_list to screen
  end type sectors_list


  !GENERIC CONSTRUCTOR
  interface sectors_list
     module procedure :: construct_sectors_tbasis
  end interface sectors_list

  !GENERIC CONSTRUCTOR
  interface sectors
     module procedure :: construct_sectors_tbasis
  end interface sectors




  !intrinsic FUNCTION SIZE(OTUPLE)
  intrinsic :: size
  interface size
     module procedure :: size_sectors_list
  end interface size


  intrinsic :: len
  interface len
     module procedure :: len_sectors_list
  end interface len

  public :: sectors_list
  public :: sectors
  public :: size
  public :: len


contains





  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Free an sectors_list (destructor) 
  !+------------------------------------------------------------------+
  elemental subroutine free_sectors_list(self)
    class(sectors_list),intent(inout) :: self
    type(qtype),pointer               :: p,c
    if(.not.associated(self%root))return
    do
       p => self%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next
       c%next => null()
       c%index=  0
       if(allocated(c%qn))deallocate(c%qn)
       if(allocated(c%map))deallocate(c%map)
       deallocate(c)
    enddo
    self%size=0
    self%root=>null()
    p=>null()
    c=>null()
  end subroutine free_sectors_list

  

  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor: 
  !+------------------------------------------------------------------+
  function construct_sectors_tbasis(basis) result(self)
    type(sectors_list),target        :: self
    type(tbasis)                     :: basis ![Nqn,Qdim]
    real(8),dimension(:),allocatable :: qn   ![Qdim]
    integer                          :: i,Nbasis
    !
    call self%free()
    allocate(self%root)
    !
    Nbasis = basis%size
    !
    do i=1,Nbasis
       qn = basis%qn(i)
       call self%put(qn,basis)
    enddo
  end function construct_sectors_tbasis











  !##################################################################
  !##################################################################
  !       PUT/LOAD  QNs list
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Put 
  !+------------------------------------------------------------------+
  subroutine put_sectors_list(self,qn,basis)
    class(sectors_list),intent(inout) :: self
    real(8),dimension(:),intent(in)   :: qn
    type(tbasis)                      :: basis
    integer                           :: i,pos,N,Qdim
    type(qtype),pointer               :: p,c
    logical                           :: iadd
    !
    if(.not.associated(self%root))allocate(self%root)
    Qdim = size(qn)
    if(self%qdim==0)self%qdim=qdim
    if(self%qdim/=qdim)stop "put_sectors_list error: size(qn) != self.qdim"
    !
    iadd = .false.
    p => self%root
    c => p%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if (all(c%qn == qn)) then
          iadd = .true.
          exit
       endif
       p => c
       c => c%next
    end do
    !
    if(iadd)then                !QN exists: create a new map
       if(allocated(c%map))deallocate(c%map)
       c%map = basis%index(qn)
    else                        !QN does not exist: create a new element
       allocate(p%next)
       p%next%qn    = qn
       p%next%index = p%index+1
       p%next%map   =  basis%index(qn)
       if(.not.associated(c))then !end of the list special case (c=>c%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
       self%size = self%size+1
    endif
    p=>null()
    c=>null()
  end subroutine put_sectors_list


  !+------------------------------------------------------------------+
  !PURPOSE:  Load 
  !+------------------------------------------------------------------+
  subroutine load_sectors_list(self,basis)
    class(sectors_list),intent(inout) :: self
    type(tbasis)                      :: basis
    real(8),dimension(:),allocatable  :: qn
    integer                           :: i
    !
    if(.not.associated(self%root))allocate(self%root)
    !
    do i=1,basis%size
       qn = basis%qn(i)
       call self%put(qn,basis)
    enddo
    !
  end subroutine load_sectors_list



  
  !+------------------------------------------------------------------+
  !PURPOSE:  Append a state in the map of a given QN if existing.
  ! If not, create it and append there.
  !+------------------------------------------------------------------+
  subroutine append_sectors_list(self,qn,istate)
    class(sectors_list),intent(inout) :: self
    integer,intent(in)                :: istate
    real(8),dimension(:),intent(in)   :: qn
    type(qtype),pointer               :: p,c
    logical                           :: iadd
    iadd = .false.
    if(.not.associated(self%root))allocate(self%root)
    p => self%root
    c => p%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if (all(c%qn == qn)) then
          iadd = .true.
          exit
       endif
       p => c
       c => c%next
    end do
    !
    if(iadd)then                !QN exists: create a new map
       call append(c%map,istate)
       call sort_quicksort(c%map)
    else                        !QN does not exist: create a new element
       allocate(p%next)
       p%next%qn    = qn
       p%next%index = p%index+1
       call append(p%next%map,istate)
       call sort_quicksort(p%next%map)
       if(.not.associated(c))then !end of the list special case (c=>c%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
       self%size = self%size+1
    endif
    p=>null()
    c=>null()
  end subroutine append_sectors_list







  !##################################################################
  !##################################################################
  !              RETRIEVE CONTENT: QN, MAP
  !##################################################################
  !##################################################################
  subroutine get_sectors_list(self,index,qn,map)
    class(sectors_list),intent(inout) :: self
    integer                           :: index
    real(8),dimension(:),allocatable  :: qn
    integer,dimension(:),allocatable  :: map
    integer                           :: index_
    type(qtype),pointer               :: c
    logical                           :: ifound
    !
    index_=index
    if(index_>self%size.OR.index_<=0)stop "get_sectors_list: index !in [1,self.size]"
    !
    if(allocated(map))deallocate(map)
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
    if(.not.ifound)stop "get error: not found"
    !
    allocate(qn, source=c%qn)
    allocate(map, source=c%map)
    !
    c=>null()
  end subroutine get_sectors_list




  !+------------------------------------------------------------------+
  !PURPOSE: Return map of the sectors_list given QN
  !+------------------------------------------------------------------+
  function get_map_sectors_list(self,index,qn) result(map)
    class(sectors_list)              :: self
    integer,optional                 :: index
    real(8),dimension(:),optional    :: qn
    integer,dimension(:),allocatable :: map
    integer                          :: index_
    type(qtype),pointer              :: c
    logical                          :: ifound
    !
    index_=self%size;if(present(index))index_=index
    if(index_>self%size.OR.index_<=0)stop "get_map_sectors_list: index !in [1,self.size]"
    if(.not.present(index).AND..not.present(qn))stop "get_map_sectors_list: no input given: use index=i OR qn=QN"
    !
    if(allocated(map))deallocate(map)
    !
    ifound=.false.
    c => self%root%next
    loop:do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if(present(qn))then
          if(all(c%qn == qn)) then
             ifound=.true.
             exit loop
          endif
       elseif(c%index == index_)then
          ifound=.true.
          exit
       endif
       c => c%next
    end do loop
    if(.not.ifound)stop "get_map error: not found"
    !
    allocate(map, source=c%map)
    !
    c=>null()
  end function get_map_sectors_list



  !+------------------------------------------------------------------+
  !PURPOSE: Return key of the sectors_list  corresponding to a given index
  !+------------------------------------------------------------------+  
  function get_qn_sectors_list(self,index) result(qn)
    class(sectors_list),intent(inout) :: self
    integer                           :: index
    integer                           :: index_
    real(8),dimension(:),allocatable  :: qn
    type(qtype),pointer               :: c
    logical                           :: ifound
    !
    index_=index
    if(index_>self%size.OR.index_<=0)stop "get_qn_sectors_list: index !in [1,self.size]"
    !
    c => self%root%next
    do                            !traverse the list until index is found
       if(.not.associated(c))exit
       if(c%index == index_) then
          ifound=.true.
          exit
       endif
       c => c%next
    end do
    if(.not.ifound)stop "get_qn error: not found"
    !
    allocate(qn, source = c%qn)
    !
    c=>null()
  end function get_qn_sectors_list



  !+------------------------------------------------------------------+
  !PURPOSE: Return the basis as a Tuple_basis
  !+------------------------------------------------------------------+  
  function basis_sectors_list(self) result(basis)
    class(sectors_list),intent(inout)  :: self
    type(tbasis)                       :: basis
    real(8),dimension(:,:),allocatable :: tvec
    real(8),dimension(:),allocatable   :: qn
    integer,dimension(:),allocatable   :: map
    integer                            :: i,j,Nbasis,Qdim
    !
    call basis%free()
    !
    Nbasis=len(self)
    Qdim  =self%qdim
    !
    allocate(tvec(Nbasis,Qdim))
    do i=1,size(self)
       qn  = self%qn(index=i)
       map = self%map(index=i)
       do j=1,size(map)
          tvec(map(j),:)=qn(:)
       enddo
    enddo
    !
    basis = tbasis(pack(transpose(tvec),.true.),qdim=qdim)
    !
  end function basis_sectors_list

  ! !+------------------------------------------------------------------+
  ! !PURPOSE: Return all the keys in the sectors_list
  ! !+------------------------------------------------------------------+  
  ! function basis_sectors_list(self) result(basis)
  !   class(sectors_list),intent(inout)  :: self
  !   real(8),dimension(:,:),allocatable :: basis
  !   real(8),dimension(:),allocatable   :: qn
  !   integer,dimension(:),allocatable   :: map
  !   integer                            :: i,j,io,Nbasis,Qdim
  !   if(allocated(basis))deallocate(basis)
  !   Nbasis=len(self)
  !   Qdim   = self%qdim
  !   allocate(basis(Nbasis,Qdim))
  !   io = 0
  !   do i=1,size(self)
  !      qn  = self%qn(index=i)
  !      map = self%map(index=i)
  !      do j=1,size(map)
  !         basis(map(j),:)=qn(:)
  !      enddo
  !   enddo
  ! end function basis_sectors_list




  !##################################################################
  !##################################################################
  !              ENUMERATOR & ITERATORS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Returns the size of given sectors_list
  !+------------------------------------------------------------------+
  function size_sectors_list(self) result(size)
    class(sectors_list),intent(in) :: self
    integer                        :: size
    size = self%size
  end function size_sectors_list


  !+------------------------------------------------------------------+
  !PURPOSE:  Returns the length of the basis states of the sectors_list
  !+------------------------------------------------------------------+
  function len_sectors_list(self) result(length)
    class(sectors_list),intent(in)   :: self
    integer                          :: length,i
    integer,dimension(:),allocatable :: map
    length=0
    do i=1,size(self)
       map = self%map(index=i)
       length = length + size(map)
    enddo
  end function len_sectors_list

  
  !+------------------------------------------------------------------+
  !PURPOSE:  Returns True is qn exists, False otherwise
  !+------------------------------------------------------------------+
  function has_qn_sectors_list(self, qn) result(bool)
    class(sectors_list),intent(inout) :: self
    real(8),dimension(:),intent(in)   :: qn
    type(qtype),pointer               :: c
    logical                           :: bool
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
    c=>null()
  end function has_qn_sectors_list








  !##################################################################
  !##################################################################
  !              PRINT  / SHOW / SPY
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Pretty print an sectors_list
  !+------------------------------------------------------------------+
  subroutine show_sectors_list(self)
    class(sectors_list),intent(inout) :: self
    integer                           :: i,count=0
    type(qtype),pointer               :: c
    !
    write(*,"(A6,I12)")"Size :",self%size
    write(*,"(A6,I12)")"Len  :",len(self)
    write(*,"(A18)")"------------------"
    c => self%root%next
    do
       if(.not.associated(c))exit
       count=count+1
       write(*,"(A6,I12)")"Index:",c%index
       write(*,"(A6)",advance='no')"Qn   :"
       call show_qn(c%qn)
       write(*,"(A6,"//str(size(c%map))//"I6)")&
            "Map  :",(c%map(i),i=1,size(c%map))
       write(*,*)""
       c => c%next
    end do
  contains
    subroutine show_qn(qn)
      real(8),dimension(:) :: qn
      integer              :: j
      write(*,"(A1)",advance='no')"["
      write(*,"(F6.2,A1)",advance='no')qn(1)
      do j=2,size(qn)
         write(*,"(A1,F6.2)",advance='no')",",qn(j)
      enddo
      write(*,"(A1)",advance='no')"]"
      write(*,*)""
    end subroutine show_qn
  end subroutine show_sectors_list



















END MODULE LIST_SECTORS
