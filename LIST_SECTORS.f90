MODULE LIST_SECTORS
  USE SCIFOR, only:str,sort_quicksort
  USE AUX_FUNCS
  implicit none
  private


  type qtype
     integer                          :: index=0
     real(8)                          :: qn=huge(1d0)
     integer,dimension(:),allocatable :: map
     type(qtype),pointer              :: next=>null()  !link to next box (chain)
  end type qtype


  type sectors_list
     integer                          :: size=0
     type(qtype),pointer              :: root=>null()
   contains
     procedure,pass :: free     => free_sectors_list     !destructor
     procedure,pass :: put      => put_sectors_list      !put sectors_list qn array
     procedure,pass :: load     => load_sectors_list     !load=sequential put=constructor
     procedure,pass :: append   => append_sectors_list   !append map state given a QN
     procedure,pass :: get      => get_sectors_list      !get qn and map for a given index
     procedure,pass :: map      => get_map_sectors_list  !get map for a given qn/index
     procedure,pass :: qn       => get_qn_sectors_list   !return qn for a given index
     procedure,pass :: basis    => basis_sectors_list    !return the basis of the sector list
     procedure,pass :: has_qn   => has_qn_sectors_list   !True if qn exists
     procedure,pass :: show     => show_sectors_list     !show sectors_list to screen
  end type sectors_list


  !GENERIC CONSTRUCTOR
  interface sectors_list
     module procedure :: construct_sectors_list_I
     module procedure :: construct_sectors_list_D
  end interface sectors_list

  !GENERIC CONSTRUCTOR
  interface sectors
     module procedure :: construct_sectors_list_I
     module procedure :: construct_sectors_list_D
  end interface sectors


  !INTRINSIC FUNCTION SIZE(OTUPLE)
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
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function construct_sectors_list_I(vec) result(self)
    type(sectors_list),target :: self
    integer,dimension(:)      :: vec
    real(8)                   :: qn
    integer                   :: iqn
    call self%free()
    allocate(self%root)
    do iqn=1,size(vec)
       qn = dble(vec(iqn))
       call self%put(qn,dble(vec))
    enddo
  end function construct_sectors_list_I

  function construct_sectors_list_D(vec) result(self)
    type(sectors_list),target :: self
    real(8),dimension(:)      :: vec
    real(8)                   :: qn
    integer                   :: iqn
    call self%free()
    allocate(self%root)
    do iqn=1,size(vec)
       qn = vec(iqn)
       call self%put(qn,vec)
    enddo
  end function construct_sectors_list_D







  !+------------------------------------------------------------------+
  !PURPOSE:  Free an sectors_list (destructor) 
  !+------------------------------------------------------------------+
  subroutine free_sectors_list(self)
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
       c%qn   = huge(1d0)
       if(allocated(c%map))deallocate(c%map)
       deallocate(c)
    enddo
    self%size=0
    self%root=>null()
    p=>null()
    c=>null()
  end subroutine free_sectors_list



  !##################################################################
  !##################################################################
  !       PUT/LOAD  QNs list
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Put a sparse matrix as operator in the sectors_list
  !+------------------------------------------------------------------+
  subroutine put_sectors_list(self,qn,vec)
    class(sectors_list),intent(inout) :: self
    real(8),dimension(:),intent(in)   :: vec
    real(8),intent(in)                :: qn
    integer                           :: i,pos,N
    type(qtype),pointer               :: p,c
    logical                           :: iadd
    !
    if(.not.associated(self%root))allocate(self%root)
    !
    iadd = .false.
    p => self%root
    c => p%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if (c%qn == qn) then
          iadd = .true.
          exit
       endif
       p => c
       c => c%next
    end do
    !
    if(iadd)then                !QN exists: create a new map
       if(allocated(c%map))deallocate(c%map)
       N   = count(vec==qn)
       pos = 0
       do i=1,N
          pos = pos+findloc(vec(pos+1:),value=qn,dim=1)
          call append(c%map,pos)
       enddo
       call sort_quicksort(c%map)
    else                        !QN does not exist: create a new element
       allocate(p%next)
       p%next%qn = qn
       p%next%index = p%index+1
       N   = count(vec==qn)
       pos = 0
       do i=1,N
          pos = pos+findloc(vec(pos+1:),value=qn,dim=1)
          call append(p%next%map,pos)
       enddo
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
  end subroutine put_sectors_list


  !+------------------------------------------------------------------+
  !PURPOSE:  Load 
  !+------------------------------------------------------------------+
  subroutine load_sectors_list(self,vec)
    class(sectors_list),intent(inout) :: self
    real(8),dimension(:)              :: vec
    integer                           :: i
    ! call self%free()
    if(.not.associated(self%root))allocate(self%root)
    do i=1,size(vec)
       call self%put(vec(i),vec)
    enddo
  end subroutine load_sectors_list



  !+------------------------------------------------------------------+
  !PURPOSE:  Append a state in the map of a given QN if existing.
  ! If not, create it and append there.
  !+------------------------------------------------------------------+
  subroutine append_sectors_list(self,qn,istate)
    class(sectors_list),intent(inout) :: self
    integer,intent(in)                :: istate
    real(8),intent(in)                :: qn
    type(qtype),pointer               :: p,c
    logical                           :: iadd
    iadd = .false.
    if(.not.associated(self%root))allocate(self%root)
    p => self%root
    c => p%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if (c%qn == qn) then
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
    integer                           :: index_
    integer,dimension(:),allocatable  :: map
    real(8)                           :: qn
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
    qn = c%qn
    allocate(map, source=c%map)
    !
    c=>null()
  end subroutine get_sectors_list


  !+------------------------------------------------------------------+
  !PURPOSE: Return map of the sectors_list given QN
  !+------------------------------------------------------------------+
  function get_map_sectors_list(self,index,qn) result(map)
    class(sectors_list)               :: self
    integer,dimension(:),allocatable  :: map
    integer,optional                  :: index
    real(8),optional                  :: qn
    integer                           :: index_
    type(qtype),pointer               :: c
    logical                           :: ifound
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
          if (c%qn == qn) then
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
    real(8)                           :: qn
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
    qn = c%qn
    !
    c=>null()
  end function get_qn_sectors_list



  !+------------------------------------------------------------------+
  !PURPOSE: Return all the keys in the sectors_list
  !+------------------------------------------------------------------+  
  function basis_sectors_list(self) result(basis)
    class(sectors_list),intent(inout) :: self
    real(8),dimension(:),allocatable  :: basis
    real(8)                           :: qn
    integer,dimension(:),allocatable  :: map
    integer                           :: i,j,io,Nbasis
    if(allocated(basis))deallocate(basis)
    Nbasis=len(self)
    allocate(basis(Nbasis))
    io = 0
    do i=1,size(self)
       qn  = self%qn(index=i)
       map = self%map(index=i)
       do j=1,size(map)
          basis(map(j))=qn
       enddo
    enddo
  end function basis_sectors_list






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
    real(8),intent(in)                :: qn
    type(qtype),pointer               :: c
    logical                           :: bool
    bool=.false.
    c => self%root%next
    do                            !traverse the list until index is found
       if(.not.associated(c))exit
       if(c%qn == qn) then
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
       write(*,"(A6,F12.6)")"Qn   :",c%qn
       write(*,"(A6,"//str(size(c%map))//"I6)")&
            "Map  :",(c%map(i),i=1,size(c%map))
       write(*,*)""
       c => c%next
    end do
  end subroutine show_sectors_list















END MODULE LIST_SECTORS
