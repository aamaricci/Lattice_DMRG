program testLIST_SECTORS
  USE SCIFOR
  USE LIST_SECTORS
  implicit none

  type(sectors_list)               :: a
  type(sectors_list) :: b(2)
  integer,dimension(:),allocatable :: map,basis
  integer                          :: i
  real(8)                          :: qn


  print*,"TEST CONSTRUCTOR 1"
  print*,[3, 5, 5, 7, 3]
  a = sectors_list([3, 5, 5, 7, 3])
  call a%show()
  print*,"size(a)=",size(a)
  call a%free()


  print*,"TEST CONSTRUCTOR 2"
  call a%load([-1d0,0d0,0d0,1d0])
  call a%show()
  call a%free()


  print*,size(b)
  b(1) = sectors_list([3, 5, 5, 7, 3])
  b(2) = sectors_list([3, 5, 5, 7, 3,10])
  print*,size(b)
  print*,size(b(1)),size(b(2))
  call b%free()

  print*,"TEST APPEND:"
  call a%load([-1d0,0d0,0d0,1d0])
  call a%show()
  call a%append(qn=1d0,istate=10)
  call a%append(qn=2d0,istate=8)
  call a%show()
  call a%free()

  print*,"TEST GET MAP"
  a = sectors(dble([1, 3, 5, 5, 7, 3, 0, 3, 1]))
  call a%show
  map= a%map(qn=3d0)
  print*,"QN:",3,allocated(map),map
  print*,""
  map= a%map(qn=1d0)
  print*,"QN:",1,allocated(map),map
  print*,""
  map= a%map(index=5)
  print*,"Indx:",5,allocated(map),map
  print*,""
  call a%free()

  call a%load([-1d0,0d0,0d0,1d0])
  call a%show()
  ! call a%append(qn=1d0,istate=10)
  ! call a%append(qn=2d0,istate=8)
  call a%show()

  print*,"TEST SIZE/LEN"
  print*,size(a)
  print*,len(a)

  basis = a%basis()
  print*,basis


  deallocate(map)
  map= a%map(index=size(a))
  print*,"Curent:",allocated(map),map
  print*,""

  basis = a%basis()
  print*,basis



  print*,a%has_qn(1d0)


  print*,""
  print*,"TEST ITERATION 2"
  do i=1,size(a)
     qn  = a%qn(index=i)
     map = a%map(index=i)
     ! call a%get(index=i,qn=qn,map=map)
     print*,i,qn,":",map
  enddo

end program testLIST_SECTORS
