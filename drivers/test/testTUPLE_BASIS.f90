program testTUPLE_BASIS
  USE SCIFOR
  USE AUX_FUNCS
  USE TUPLE_BASIS
  implicit none

  integer,dimension(:),allocatable   :: qn
  integer                            :: i
  type(tbasis)                       :: basis,a,b,v1,v2
  integer,dimension(:),allocatable   :: map
  real(8),dimension(:,:),allocatable :: tvec
  real(8),dimension(:),allocatable   :: array

  ! vec1 = tvec([0,0, 1,0, 0,1, 1,1, 0,0, 0,1],Qdim=2)
  ! vec2 = tvec([0,0, 5,1, 5,1, 7,2, 0,0, 3,2],2)

  print*,"TEST CONSTRUCTOR"
  qn  = [0,0]
  basis = tbasis([0,0, 1,0, 0,1, 1,1, 0,0, 0,1],Qdim=2)
  call basis%show
  print*,"shape:",shape(basis)
  print*,"size :",basis%size
  print*,""


  print*,"TEST GET QN"
  do i=1,basis%size
     qn = basis%qn(i)
     print*,"QN:",qn
  enddo
  print*,""



  print*,"TEST INDEX"
  print*,"map [0,0]"
  map = basis%index([0d0,0d0])
  print*,map

  print*,"map [1,0]"
  map = basis%index([1d0,0d0])
  print*,map

  print*,"map [0,1]"
  map = basis%index([0d0,1d0])
  print*,map

  print*,"map [1,1]"
  map = basis%index([1d0,1d0])
  print*,map


  call basis%free()

  print*,"TEST EXPAND 1"
  basis = tbasis([0d0, 1d0, 0d0, 1d0, 0d0, 0d0],Qdim=1)
  call basis%show
  print*,shape(basis)
  call basis%expand([0d0, 0d0, 1d0, 1d0, 0d0, 1d0])
  call basis%show
  print*,shape(basis)


  print*,"TEST EXPAND 2"
  call basis%free()
  call basis%expand([0d0, 1d0, 0d0, 1d0, 0d0, 0d0])
  call basis%show
  call basis%expand([0d0, 0d0, 1d0, 1d0, 0d0, 1d0])
  call basis%show
  print*,shape(basis)


  print*,"TEST DUMP"
  tvec = basis%dump()
  call show_tvec(tvec)

  print*,"TEST FLAT"
  array = basis%flat()
  print*,[0,0, 1,0, 0,1, 1,1, 0,0, 0,1]
  print*,int(array)


  print*,"TEST SUM"
  a = tbasis([0,0, 1,0, 0,1, 1,1, 0,0, 0,1],Qdim=2)
  b = tbasis([0,1, 1,1, 0,1, 0,0, 0,1, 0,0],Qdim=2)
  basis=a+b
  call basis%show
  print*,[0,0, 1,0, 0,1, 1,1, 0,0, 0,1]+[0,1, 1,1, 0,1, 0,0, 0,1, 0,0]
  print*,int(basis%flat())




  print*,"TEST OUTSUM"
  a = tbasis([0,0, 1,0, 0,1, 1,1],Qdim=2)
  b = tbasis([0,0, 1,0, 0,1, 1,1],Qdim=2)
  basis=a.o.b
  call basis%show
  print*,int(basis%flat())
  print*,""

  call basis%free()
  basis = tsum(a,b)
  call basis%show


  call basis%free()
  call basis%append([0d0,0d0])
  call basis%show
  call basis%append([1d0,0d0])
  call basis%show
  call basis%append([0d0,1d0])
  call basis%show
  call basis%append([1d0,1d0])
  call basis%show

  
contains


  subroutine show_tvec(tvec)
    real(8),dimension(:,:) :: tvec
    integer                :: Nqn,Qdim,i,j
    Nqn = size(tvec,1)
    Qdim= size(tvec,2)
    write(*,*)"Show Tvec:"
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

end program testTUPLE_BASIS
