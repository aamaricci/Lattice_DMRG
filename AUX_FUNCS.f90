MODULE AUX_FUNCS
  USE SCIFOR, only: free_unit
  implicit none
  private


  !AUX:
  interface append
     module procedure :: append_I
     module procedure :: append_D
     module procedure :: append_C
  end interface append

  interface add_to
     module procedure :: append_I
     module procedure :: append_D
     module procedure :: append_C
  end interface add_to

  public :: outsum
  public :: append
  public :: add_to
  public :: binary_search
  public :: KId
  public :: KSz
  public :: fopen
  
  logical,parameter,public           :: show_dble=.true.
  character(len=12),parameter,public :: show_fmt='F9.3'

contains



  function KId(n) result(A)
    integer, intent(in) :: n
    real(8)             :: A(2**n, 2**n)
    integer             :: i
    A = 0d0
    forall(i=1:2**n)A(i,i) = 1d0
  end function KId

  recursive function KSz(n) result(A)
    integer, intent(in) :: n
    real(8)             :: A(2**n, 2**n)
    integer             :: d(2**n)
    integer             :: i
    d = szvec(n)
    A = 0d0
    forall(i=1:2**n)A(i,i) = dble(d(i))
  end function KSz

  recursive function szvec(n) result(vec)
    integer,intent(in)      :: n
    integer,dimension(2**n) :: vec
    if(n==1)then
       vec = [1,-1]
    else
       vec = [szvec(n-1),-szvec(n-1)]
    endif
  end function szvec



  !##################################################################
  !##################################################################
  !              AUXILIARY COMPUTATIONAL ROUTINES
  !##################################################################
  !##################################################################
  function outsum(a,b) result(c)
    real(8),dimension(:)               :: a
    real(8),dimension(:)               :: b
    real(8),dimension(size(a)*size(b)) :: c
    integer                            :: i,j,io
    c=0d0
    do i=1,size(a)
       do j=1,size(b)
          io = j + (i-1)*size(b)
          c(io) = a(i) + b(j)
       enddo
    enddo
  end function outsum


  pure subroutine append_I(vec,val)
    integer,dimension(:),allocatable,intent(inout) :: vec
    integer,intent(in)                             :: val  
    integer,dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine append_I

  pure subroutine append_D(vec,val)
    real(8),dimension(:),allocatable,intent(inout) :: vec
    real(8),intent(in)                             :: val  
    real(8),dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine append_D

  pure subroutine append_C(vec,val)
    complex(8),dimension(:),allocatable,intent(inout) :: vec
    complex(8),intent(in)                             :: val  
    complex(8),dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine append_C






  !+------------------------------------------------------------------+
  !PURPOSE  : Sort an array, gives the new ordering of the label.
  !+------------------------------------------------------------------+
  recursive function binary_search(Ain,value) result(bsresult)
    integer,intent(in)           :: Ain(:), value
    integer                      :: bsresult, mid
    integer,dimension(size(Ain)) :: A,Order
    !
    a = ain
    call sort_array(a,Order)
    !
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
    !
    bsresult = Order(bsresult)
    !
  end function binary_search



  !+------------------------------------------------------------------+
  !PURPOSE  : Sort an array, gives the new ordering of the label.
  !+------------------------------------------------------------------+
  subroutine sort_array(array,order)
    implicit none
    integer,dimension(:)                    :: array
    integer,dimension(size(array))          :: order
    integer,dimension(size(array))          :: backup
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort(array, order,1, size(array))
    do i=1,size(array)
       backup(i)=array(order(i))
    enddo
    array=backup
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:) :: array
      integer, dimension(:) :: order
      integer               :: left
      integer               :: right
      integer               :: i
      integer               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:) :: order
      integer               :: first, second
      integer               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    integer function qsort_rand( lower, upper )
      integer               :: lower, upper
      real(8)               :: r
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    !---------------------------------------------!
    function compare(f,g)
      implicit none
      integer               :: f,g
      integer               :: compare
      if(f<g) then
         compare=-1
      else
         compare=1
      endif
    end function compare
  end subroutine sort_array





  function fopen(fname,append) result(unit)
    character(len=*) :: fname
    logical,optional :: append
    integer          :: unit
    logical          :: append_,bool
    append_=.true.;if(present(append))append_=append
    select case(append_)
    case (.true.)    
       inquire(file=trim(fname), exist=bool)
       unit = free_unit()
       if (bool) then
          open(unit,file=trim(fname),status="old",position="append",action="write")
       else
          open(unit,file=trim(fname),status="new",action="write")
       end if
    case(.false.)
       open(unit,file=trim(fname),status="new",action="write")
    end select
  end function fopen

END MODULE AUX_FUNCS
