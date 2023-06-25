!
! We assume quantum numbers are given in terms of a tuple Q=[n_1,...,n_qn]^T
! Then we define a Tvec, ie a vector of tuples as, V=[Q_1,...,Q_N]
! as an array rank-2 with dimensions [N,N_qn]. For each i=1,...,N we have a tuple
! Q_i.
! The goal here is to create and index, or a MAP, of these tuples in a givem Tvec.
! ie return the state indices (as positions in Tvec) with a given value of Q.
!
MODULE LIST_SECTORS
  USE SCIFOR, only:str,sort_quicksort,assert_shape
  USE AUX_FUNCS
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
     module procedure :: construct_sectors_list_I
     module procedure :: construct_sectors_list_D
  end interface sectors_list

  !GENERIC CONSTRUCTOR
  interface sectors
     module procedure :: construct_sectors_list_I
     module procedure :: construct_sectors_list_D
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



  interface tvec
     module procedure :: construct_tvec_I
     module procedure :: construct_tvec_D
  end interface tvec

  interface tflat
     module procedure :: pack_tvec_I
     module procedure :: pack_tvec_D
  end interface tflat

  public :: tvec
  public :: tflat
  public :: show_tvec

contains





  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor: 
  !+------------------------------------------------------------------+
  function construct_sectors_list_I(tvec) result(self)
    type(sectors_list),target       :: self
    integer,dimension(:,:)          :: tvec ![Nqn,Qdim]
    real(8),dimension(size(tvec,2)) :: qn   ![Qdim]
    integer                         :: iqn,Nqn
    call self%free()
    allocate(self%root)
    !
    Nqn = size(tvec,1)           !Number of tuples in this tvec
    !
    do iqn=1,Nqn
       qn = dble(tvec(iqn,:))
       call self%put(qn,dble(tvec))
    enddo
  end function construct_sectors_list_I

  function construct_sectors_list_D(tvec) result(self)
    type(sectors_list),target       :: self
    real(8),dimension(:,:)          :: tvec ![Nqn,Qdim]
    real(8),dimension(size(tvec,2)) :: qn   ![Qdim]
    integer                         :: iqn,Nqn
    call self%free()
    allocate(self%root)
    !
    Nqn = size(tvec,1)           !Number of tuples in this tvec
    !
    do iqn=1,Nqn
       qn = tvec(iqn,:)
       call self%put(qn,tvec)
    enddo
  end function construct_sectors_list_D


  

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







  !##################################################################
  !##################################################################
  !       PUT/LOAD  QNs list
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Put 
  !+------------------------------------------------------------------+
  subroutine put_sectors_list(self,qn,tvec)
    class(sectors_list),intent(inout) :: self
    real(8),dimension(:),intent(in)   :: qn  ![Qdim]
    real(8),dimension(:,:),intent(in) :: tvec ![Nqn,Qdim]
    integer                           :: i,pos,N,Qdim,Nqn
    type(qtype),pointer               :: p,c
    logical                           :: iadd
    !
    if(.not.associated(self%root))allocate(self%root)
    Qdim = size(qn)
    Nqn  = size(tvec,1)
    if(self%qdim==0)self%qdim=qdim
    if(self%qdim/=qdim)stop "put_sectors_list error: size(qn) != self.qdim"
    call assert_shape(tvec,[Nqn,Qdim],"put_sectors_list","tvec")
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
       c%map = index_tvec(tvec,qn)
    else                        !QN does not exist: create a new element
       allocate(p%next)
       p%next%qn    = qn
       p%next%index = p%index+1
       p%next%map   = index_tvec(tvec,qn)
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
  subroutine load_sectors_list(self,tvec)
    class(sectors_list),intent(inout) :: self
    real(8),dimension(:,:)            :: tvec ![Nqn,Qdim]
    real(8),dimension(size(tvec,2))   :: qn
    integer                           :: i
    ! call self%free()
    if(.not.associated(self%root))allocate(self%root)
    do i=1,size(tvec,1)
       qn = tvec(i,:)
       call self%put(qn,tvec)
    enddo
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
  !PURPOSE: Return all the keys in the sectors_list
  !+------------------------------------------------------------------+  
  function basis_sectors_list(self) result(basis)
    class(sectors_list),intent(inout)  :: self
    real(8),dimension(:,:),allocatable :: basis
    real(8),dimension(:),allocatable   :: qn
    integer,dimension(:),allocatable   :: map
    integer                            :: i,j,io,Nbasis,Qdim
    if(allocated(basis))deallocate(basis)
    Nbasis=len(self)
    Qdim   = self%qdim
    allocate(basis(Nbasis,Qdim))
    io = 0
    do i=1,size(self)
       qn  = self%qn(index=i)
       map = self%map(index=i)
       do j=1,size(map)
          basis(map(j),:)=qn(:)
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












  !##################################################################
  !##################################################################
  !       CONSTRUCT TVECTOR FROM ARRAY 
  !##################################################################
  !##################################################################
  function construct_tvec_I(array,qdim) result(tvec)
    integer,dimension(:)               :: array
    integer                            :: Qdim,Nqn,N
    real(8),dimension(:,:),allocatable :: tvec
    if(allocated(tvec))deallocate(tvec)
    N = size(array)
    Nqn = N/Qdim;if(mod(N,Qdim)/=0)stop "construct_tvec_I error: size(array)%Qdim!=0"
    tvec = dble(transpose(reshape(array, shape=[Qdim,Nqn])))
  end function construct_tvec_I

  function construct_tvec_D(array,qdim) result(tvec)
    real(8),dimension(:)               :: array
    integer                            :: Qdim,Nqn,N
    real(8),dimension(:,:),allocatable :: tvec
    if(allocated(tvec))deallocate(tvec)
    N = size(array)
    Nqn = N/Qdim;if(mod(N,Qdim)/=0)stop "construct_tvec_I error: size(array)%Qdim!=0"
    tvec = transpose(reshape(array, shape=[Qdim,Nqn]))
  end function construct_tvec_D



  !##################################################################
  !##################################################################
  !              PACK TVECTOR TO ARRAY
  !##################################################################
  !##################################################################
  function pack_tvec_I(tvec) result(array)
    integer,dimension(:,:)           :: tvec
    integer,dimension(:),allocatable :: array
    if(allocated(array))deallocate(array)
    array = pack(transpose(tvec),.true.)
  end function pack_tvec_I

  function pack_tvec_D(tvec) result(array)
    real(8),dimension(:,:)           :: tvec
    real(8),dimension(:),allocatable :: array
    if(allocated(array))deallocate(array)
    array = pack(transpose(tvec),.true.)
  end function pack_tvec_D



  !##################################################################
  !##################################################################
  !              INDEX FROM TVECTOR
  !##################################################################
  !##################################################################
  function index_tvec(tvec,qn) result(index)
    real(8),dimension(:,:)           :: tvec
    real(8),dimension(size(tvec,2))  :: qn
    integer,dimension(:),allocatable :: index
    logical,dimension(size(tvec,1))  :: mask
    integer                          :: i,N,pos
    !
    if(allocated(index))deallocate(index)
    forall(i=1:size(mask))mask(i) = all(tvec(i,:)==qn)
    N   = count(mask)
    pos = 0
    do i=1,N
       pos = pos+findloc(mask(pos+1:),value=.true.,dim=1)
       call append(index,pos)
    enddo
    call sort_quicksort(index)
  end function index_tvec



  !##################################################################
  !##################################################################
  !              SHOW TVECTOR
  !##################################################################
  !##################################################################
  subroutine show_tvec(tvec)
    real(8),dimension(:,:) :: tvec
    integer                :: Nqn,Qdim,i,j
    Nqn = size(tvec,1)
    Qdim= size(tvec,2)
    do i=1,Nqn
       write(*,"(A1)",advance='no')"["
       write(*,"(F6.2,A1)",advance='no')tvec(i,1)
       do j=2,Qdim
          write(*,"(A1,F6.2)",advance='no')",",tvec(i,j)
       enddo
       write(*,"(A1)",advance='no')"]"
       write(*,*)""
    enddo
    write(*,*)""
  end subroutine show_tvec








END MODULE LIST_SECTORS
