! We assume quantum numbers are given in terms of a tuple Q=[q_1,...,q_M],
! e.g. M=2 with q_1=n_up, q_2=n_dw
! Then we employ a TUPLE_BASIS, ie a vector of tuples as V=[Q_1,...,Q_N]^T.
! For each i=1,...,N we have a tuple Q_i.
! Here we create and index, ie MAP, the sectors based on such tuples of QN.
! ie we return the state indices (as positions in the tuple_basis) with
! a given value of Q.
MODULE LIST_SECTORS
  USE SCIFOR, only:str,sort_quicksort,assert_shape,free_unit
  USE AUX_FUNCS
  USE TUPLE_BASIS
  implicit none
  private


  type qtype
     integer                          :: index=0
     real(8),dimension(:),allocatable :: qn   !tuple of quantum numbers
     integer,dimension(:),allocatable :: map  !map of the tuples, the states with a given tuple of qn.
     type(qtype),pointer              :: next=>null()
  end type qtype


  type sectors_list
     integer             :: qdim=0
     integer             :: size=0
     type(qtype),pointer :: root=>null()
   contains
<<<<<<< HEAD
<<<<<<< HEAD:src/LIST/LIST_SECTORS.f90
<<<<<<< HEAD:src/LIST/LIST_SECTORS.f90
=======
>>>>>>> d3539b5 (2.1.0 UPDATED STABLE VERSION):LIST_SECTORS.f90
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
     procedure,pass :: free   => free_sectors_list     !destructor
     procedure,pass :: put    => put_sectors_list      !put sectors_list qn array
     procedure,pass :: load   => load_sectors_list     !load=sequential put=constructor
     procedure,pass :: append => append_sectors_list   !append map state given a QN
     procedure,pass :: appends=> appends_sectors_list   !append map state given a QN
     procedure,pass :: get    => get_sectors_list      !get qn and map for a given index
     procedure,pass :: map    => get_map_sectors_list  !get map for a given qn/index
     procedure,pass :: qn     => get_qn_sectors_list   !return qn for a given index
     procedure,pass :: index  => get_index_sectors_list!return index for a given qn
     procedure,pass :: basis  => basis_sectors_list    !return the basis of the sector list
     procedure,pass :: has_qn => has_qn_sectors_list   !True if qn exists
     procedure,pass :: show   => show_sectors_list     !show sectors_list to screen
<<<<<<< HEAD
<<<<<<< HEAD:src/LIST/LIST_SECTORS.f90
=======
     procedure,pass      :: free     => free_sectors_list     !destructor
     procedure,pass      :: put      => put_sectors_list      !put sectors_list qn array
     procedure,pass      :: load     => load_sectors_list     !load=sequential put=constructor
     procedure,pass      :: append   => append_sectors_list   !append map state given a QN
     procedure,pass      :: get      => get_sectors_list      !get qn and map for a given index
     procedure,pass      :: map      => get_map_sectors_list  !get map for a given qn/index
     procedure,pass      :: qn       => get_qn_sectors_list   !return qn for a given index
     procedure,pass      :: index    => get_index_sectors_list   !return index for a given qn
     procedure,pass      :: basis    => basis_sectors_list    !return the basis of the sector list
     procedure,pass      :: has_qn   => has_qn_sectors_list   !True if qn exists
     procedure,pass      :: show     => show_sectors_list     !show sectors_list to screen
>>>>>>> 4174253 (Intermediate commit, working on Kron_hm_1d_2bands):LIST_SECTORS.f90
=======
>>>>>>> d3539b5 (2.1.0 UPDATED STABLE VERSION):LIST_SECTORS.f90
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
  end type sectors_list


  !GENERIC CONSTRUCTOR
  interface sectors_list
     module procedure :: construct_sectors_tbasis
  end interface sectors_list

  !GENERIC CONSTRUCTOR
  interface sectors
     module procedure :: construct_sectors_tbasis
  end interface sectors


  !EQUALITY operator A=B  [deep copy]
  interface assignment(=)
     module procedure :: equality_sector_list
  end interface assignment(=)


  !INTRINSIC FUNCTION SIZE:
  intrinsic :: size
  interface size
     module procedure :: size_sectors_list
  end interface size

  !FUNCTION DIM:
  interface dim
     module procedure :: len_sectors_list
  end interface dim

  public :: sectors_list
  public :: sectors
  public :: size
  public :: dim
  public :: assignment(=)


<<<<<<< HEAD
<<<<<<< HEAD:src/LIST/LIST_SECTORS.f90
<<<<<<< HEAD:src/LIST/LIST_SECTORS.f90


=======
  
  
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):LIST_SECTORS.f90
=======


>>>>>>> 4174253 (Intermediate commit, working on Kron_hm_1d_2bands):LIST_SECTORS.f90
=======


>>>>>>> 7e90d6a (Updating Cmake library construction)
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


  !+------------------------------------------------------------------+
  ! !PURPOSE:  Intrinsic constructor: 
  ! !+------------------------------------------------------------------+
  ! function construct_from_sectors(b) result(self)
  !   type(sectors_list),intent(in)    :: b
  !   type(sectors_list)               :: self
  !   type(tbasis)                     :: basis ![Nqn,Qdim]
  !   real(8),dimension(:),allocatable :: qn   ![Qdim]
  !   integer                          :: i,Nbasis
  !   !
  !   call self%free()
  !   do i=1,Nbasis
  !      qn = basis%qn(i)
  !      call self%put(qn,basis)
  !   enddo
  ! end function construct_from_sectors









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
       p%next%map   = basis%index(qn)
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
  pure subroutine append_sectors_list(self,qn,istate)
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
       ! call sort_quicksort(c%map)
    else                        !QN does not exist: create a new element
       allocate(p%next)
       p%next%qn    = qn
       p%next%index = p%index+1
       call append(p%next%map,istate)
       ! call sort_quicksort(p%next%map)
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


  pure subroutine appends_sectors_list(self,qn,istates)
    class(sectors_list),intent(inout) :: self
    integer,dimension(:),intent(in)   :: istates
    real(8),dimension(:),intent(in)   :: qn
    type(qtype),pointer               :: p,c
    logical                           :: iadd
    integer :: i
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
       do i=1,size(istates)
          call append(c%map,istates(i))
       enddo
       ! call sort_quicksort(c%map)
    else                        !QN does not exist: create a new element
       allocate(p%next)
       p%next%qn    = qn
       p%next%index = p%index+1
       do i=1,size(istates)
          call append(p%next%map,istates(i))
       enddo
       ! call sort_quicksort(p%next%map)
       if(.not.associated(c))then !end of the list special case (c=>c%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
       self%size = self%size+1
    endif
    p=>null()
    c=>null()
  end subroutine appends_sectors_list





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
    class(sectors_list),intent(in)   :: self
    integer                          :: index
    integer                          :: index_
    real(8),dimension(:),allocatable :: qn
    type(qtype),pointer              :: c
    logical                          :: ifound
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
  !PURPOSE: Return key of the sectors_list  corresponding to a given index
  !+------------------------------------------------------------------+  
  function get_index_sectors_list(self,qn) result(index)
    class(sectors_list),intent(inout) :: self
    real(8),dimension(:),intent(in)   :: qn
    integer                           :: index
    real(8),dimension(:),allocatable  :: qn_
    integer                           :: i
    type(qtype),pointer               :: c
    logical                           :: ifound=.false.
    !
    index=0
    if(.not. self%has_qn(qn) )return
    !
    do i=1,size(self)
       qn_ = self%qn(index=i)
       if(all(qn == qn_))then
          ifound=.true.
          index = i
          exit
       endif
    end do
    if(.not.ifound)stop "get_index error: QN not found"
  end function get_index_sectors_list



  !+------------------------------------------------------------------+
  !PURPOSE: Return the basis as a Tuple_basis
  !+------------------------------------------------------------------+  
  function basis_sectors_list(self) result(basis)
    class(sectors_list),intent(in)     :: self
    type(tbasis)                       :: basis
    real(8),dimension(:,:),allocatable :: tvec
    real(8),dimension(:),allocatable   :: qn
    integer,dimension(:),allocatable   :: map
    integer                            :: i,j,Nbasis,Qdim
    !
    call basis%free()
    !
    Nbasis=dim(self)
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





  !+------------------------------------------------------------------+
  !PURPOSE:  Tbasis equality A = B. Deep copy
  !+------------------------------------------------------------------+
  subroutine equality_sector_list(a,b)
    type(sectors_list),intent(inout) :: a
    type(sectors_list),intent(in)    :: b
    type(tbasis)                     :: basis ![Nqn,Qdim]
    real(8),dimension(:),allocatable :: qn    ![Qdim]
    integer                          :: i,Nbasis
    call a%free()
    basis = b%basis()
    do i=1,basis%size
       qn = basis%qn(i)
       call a%put(qn,basis)
    enddo
  end subroutine equality_sector_list



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
    class(sectors_list),intent(in) :: self
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
  subroutine show_sectors_list(self,unit,file)
    class(sectors_list),intent(inout) :: self
    integer,optional :: unit
    integer                           :: i,count=0
    type(qtype),pointer               :: c
    character(len=*),optional       :: file
    integer                         :: unit_
    unit_=6
    if(present(unit))unit_=unit
    if(present(file))open(free_unit(unit_),file=str(file))
    !
    write(unit_,"(A6,I12)")"Size :",self%size
    write(unit_,"(A6,I12)")"Dim  :",dim(self)
    write(unit_,"(A18)")"------------------"
    c => self%root%next
    do
       if(.not.associated(c))exit
       count=count+1
       write(unit_,"(A6,I12)")"Index:",c%index
       write(unit_,"(A6)",advance='no')"Qn   :"
       call show_qn(c%qn)
       write(unit_,"(A6,"//str(size(c%map))//"I6)")&
            "Map  :",(c%map(i),i=1,size(c%map))
       write(unit_,*)""
       c => c%next
    end do
    if(present(file))close(unit_)
  contains
    subroutine show_qn(qn)
      real(8),dimension(:) :: qn
      integer              :: j
      write(unit_,"(A1)",advance='no')"["
      write(unit_,"(F6.2,A1)",advance='no')qn(1)
      do j=2,size(qn)
         write(unit_,"(A1,F6.2)",advance='no')",",qn(j)
      enddo
      write(unit_,"(A1)",advance='no')"]"
      write(unit_,*)""
    end subroutine show_qn
  end subroutine show_sectors_list



















END MODULE LIST_SECTORS









!##################################################################
!##################################################################
!##################################################################
!##################################################################
!                          /_  __/ ____/ ___/_  __/
!                           / / / __/  \__ \ / /   
!                          / / / /___ ___/ // /    
!                         /_/ /_____//____//_/     
!##################################################################
!##################################################################
!##################################################################
!##################################################################
#ifdef _TEST
program testLIST_SECTORS
  USE SCIFOR
  USE AUX_FUNCS
  USE TUPLE_BASIS
  USE LIST_SECTORS
  implicit none

  type(sectors_list)               :: a,dot,c
  type(sectors_list)               :: b(2)
  integer,dimension(:),allocatable :: map,qn
  type(tbasis)                     :: vec1,vec2,vec3,basis,a_basis,dot_basis
  logical,dimension(:),allocatable :: mask
  integer                          :: i,j

  vec1 = tbasis([0,0, 1,0, 0,1, 1,1, 0,0, 0,1],Qdim=2)
  vec3 = tbasis([0,0, 5,1, 5,1, 7,2, 0,0, 3,2],2)
  vec2 = tbasis([0,0, 1,0, 0,1, 1,1],Qdim=2)
  print*,shape(vec1)
  print*,shape(vec2)
  print*,shape(vec3)

  print*,"TEST CONSTRUCTOR 1: vec2"
  call vec2%show
  a = sectors_list( vec2 )
  call a%show()
  print*,"size(a)=",size(a)
  call a%free()
  print*,""

  print*,"TEST CONSTRUCTOR 2: vec1"
  call vec1%show
  call a%load( vec1)
  call a%show()
  call a%free()
  print*,""


  print*,size(b)
  b(1) = sectors_list(vec1)
  b(2) = sectors_list(vec2)
  print*,size(b)
  print*,size(b(1)),size(b(2))
  call b%free()

  print*,"TEST APPEND:"
  call vec1%show()
  call a%load( vec1 )
  call a%show()
  print*,"append qn=[0d0,0d0] at state=10"
  call a%append(qn=[0d0,0d0],istate=10)
  print*,"append qn=[2d0,0d0] at state=8"
  call a%append(qn=[2d0,0d0],istate=8)
  call a%show()
  call a%free()

  print*,"TEST GET MAP"
  a = sectors( vec1 )
  call a%show
  map= a%map(qn=[0d0,0d0])
  print*,"QN: [",[0d0,0d0],"]",allocated(map),map
  print*,""
  map= a%map(qn=[1d0,0d0])
  print*,"QN: [",[1d0,0d0],"]",allocated(map),map
  print*,""
  map= a%map(index=3)
  print*,"Indx:",3,allocated(map),map
  print*,""



  print*,"TEST SIZE/LEN"
  print*,size(a)
  print*,dim(a)


  print*,"TEST RETURN BASIS INTO A TUPLE_BASIS"  
  basis = a%basis()
  call basis%show()



  print*,"HAS_QN [1,0]:T"
  print*,a%has_qn([1d0,0d0])

  print*,"HAS_QN [5,23]:F"
  print*,a%has_qn([5d0,23d0])


  print*,"TEST RETURN BASIS INTO A TUPLE_BASIS"    
  a   = sectors_list( tbasis([0,0, 1,0, 0,1, 1,1],Qdim=2) )
  dot = sectors_list( tbasis([0,0, 1,0, 0,1, 1,1],Qdim=2) )

  basis = a%basis().o.dot%basis()
  call basis%show()
  c = sectors_list( basis)
  call c%show()


  print*,"Check DEEP COPY ="
  a   = sectors_list( tbasis([0,0, 1,0, 0,1, 1,1],Qdim=2) )
  dot = sectors_list( tbasis([0,1, 0,0, 0,1, 1,1],Qdim=2) )  

  b(1) = a
  b(2) = dot
  call a%show()
  call dot%show()

  call b(1)%show()
  call b(2)%show()

end program testLIST_SECTORS
#endif
