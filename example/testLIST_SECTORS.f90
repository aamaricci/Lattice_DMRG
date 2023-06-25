program testLIST_SECTORS
  USE SCIFOR
  USE AUX_FUNCS
  USE LIST_SECTORS
  implicit none

  type(sectors_list)                 :: a
  type(sectors_list)                 :: b(2)
  integer,dimension(:),allocatable   :: map,qn
  integer,dimension(:,:),allocatable :: vec1,vec2,vec3,basis
  logical,dimension(:),allocatable   :: mask
  integer                            :: i,j

  qn  = [0,0]
  vec1 = tvec([0,0, 1,0, 0,1, 1,1, 0,0, 0,1],Qdim=2)
  vec2 = tvec([0,0, 5,1, 5,1, 7,2, 0,0, 3,2],2)
  print*,shape(vec1)
  call show_tvec(dble(vec1))


 
  
  print*,"map [0,0]"
  map = index_tvec(vec1,[0,0])
  print*,map

  print*,"map [1,0]"
  map = index_tvec(vec1,[1,0])
  print*,map

  print*,"map [0,1]"
  map = index_tvec(vec1,[0,1])
  print*,map

  print*,"map [1,1]"
  map = index_tvec(vec1,[1,1])
  print*,map
  print*,""

  
  print*,[0,0, 1,0, 0,1, 1,1, 0,0, 0,1]
  print*,pack(transpose(vec1),.true.)
  print*,tflat(vec1)

  vec3 = tvec([tflat(vec1),0,0],Qdim=2)
  call show_tvec(dble(vec3))
  
  print*,"TEST CONSTRUCTOR 1"
  a = sectors_list( vec2 )
  call a%show()
  print*,"size(a)=",size(a)
  call a%free()
  print*,""

  print*,"TEST CONSTRUCTOR 2"
  call a%load( dble(vec1))
  call a%show()
  call a%free()
  print*,""

  
  print*,size(b)
  b(1) = sectors_list(vec1)
  b(2) = sectors_list(vec3)
  print*,size(b)
  print*,size(b(1)),size(b(2))
  call b%free()

  print*,"TEST APPEND:"
  call show_tvec(dble(vec1))
  call a%load( dble(vec1) )
  call a%show()
  call a%append(qn=[0d0,0d0],istate=10)
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
  print*,len(a)

  print*,vec1
  call show_tvec(dble(vec1))
  basis = a%basis()
  print*,basis
  call show_tvec(dble(basis))

  ! deallocate(map)
  ! map= a%map(index=size(a))
  ! print*,"Curent:",allocated(map),map
  ! print*,""

  ! basis = a%basis()
  ! print*,basis


  print*,"HAS_QN"
  print*,a%has_qn([1d0,0d0])

  print*,"HAS_QN"
  print*,a%has_qn([5d0,23d0])

  ! print*,""
  ! print*,"TEST ITERATION 2"
  ! do i=1,size(a)
  !    qn  = a%qn(index=i)
  !    map = a%map(index=i)
  !    ! call a%get(index=i,qn=qn,map=map)
  !    print*,i,qn,":",map
  ! enddo

contains

  function index_tvec(tvec,qn) result(index)
    integer,dimension(:,:)           :: tvec
    integer,dimension(size(tvec,2))  :: qn
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


end program testLIST_SECTORS
