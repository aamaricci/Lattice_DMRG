program testOPERATORS_TUPLE
  USE SCIFOR
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  implicit none

  type(operators_list)               :: my_list
  type(operators_list)               :: copy_list
  type(sparse_matrix)                :: spSz,spSp,spH,spK,a,b,c
  real(8),dimension(:,:),allocatable :: mat
  integer                            :: i,j,n
  logical                            :: bool
  real(8),dimension(2,2),parameter   :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
  real(8),dimension(2,2),parameter   :: Sz=dble(pauli_z)
  real(8),dimension(2,2),parameter   :: Sx=dble(pauli_x)
  real(8),dimension(2,2),parameter   :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
  real(8),dimension(4,4)             :: Gamma13,Gamma03
  character(len=10) :: key
  character(len=10),allocatable :: keys(:)

  Gamma13=kron(Sx,Sz)
  Gamma03=kron(eye(2),Sz)


  print*,"TEST CONSTRUCTOR, PUT, SHOW, FREE"
  my_list = operators_list(&
       ['H0','Sz','Sp'],&
       [sparse(Hzero),sparse(Sz),sparse(Splus)])
  call my_list%show()
  call my_list%free()


  print*,"TEST LOAD matrices"
  call my_list%load("H0",Hzero)
  call my_list%load("Sz",Sz)
  call my_list%load("Sp",Splus)
  print*,"TEST SHOW"
  call my_list%show()
  call my_list%free()

  print*,"TEST APPEND + LOAD matrices"
  call my_list%append("H0",as_sparse(Hzero))
  call my_list%append("Sz",as_sparse(Sz))
  call my_list%load("Sp",Splus)  
  print*,"TEST SHOW"
  call my_list%show()



  print*,"TEST RETRIEVE FUNCTIONALITIES"
  print*,"TEST .DUMP"
  print*,"Mat.allocated:",allocated(mat)
  print*,"dump Sp -> Mat"
  mat = my_list%dump("Sp")
  print*,"Mat.allocated:",allocated(mat)
  do i=1,size(mat,1)
     write(*,*)(mat(i,j),j=1,size(mat,2))
  enddo
  deallocate(mat)

  print*,"TEST .GET"
  do i=1,size(my_list)
     call my_list%get(index=i,key=key,op=a)
     print*,i
     print*,key
     call a%show()
  enddo

  print*,"TEST .KEY + .OP + ITERATION over index"
  do i=1,size(my_list)
     a = my_list%op(index=i)
     print*,i
     call a%show()
  enddo
  print*,""

  do i=1,size(my_list)
     key = my_list%key(index=i)
     a = my_list%op(key=key)
     print*,i,key
     call a%show
  enddo


  print*,"TEST HAS_KEY"
  print*,"list has key Sz",my_list%has_key("Sz")
  print*,"list has key SZ",my_list%has_key("SZ")
  print*,""



  print*,"TEST IS_VALID "
  print*,my_list%is_valid()
  print*,"is valid with dim=2"
  print*,my_list%is_valid(dim=2)
  print*,"is not valid with dim=3"
  print*,my_list%is_valid(dim=3)
  print*,"is not valid once appended s_0.x.s_3"
  call my_list%append("W",as_sparse(Gamma03))
  print*,my_list%is_valid()
  call my_list%free
  print*,""



  call my_list%append("H0",as_sparse(Hzero))
  call my_list%append("Sz",as_sparse(Sz))
  call my_list%load("Sp",Splus)



  print*,"TEST ="
  copy_list = my_list
  call copy_list%show
  print*,""  



  print*,"TEST my_list.o('key')"
  print*,"before a=empty"
  call a%free
  call a%show
  print*,"a = my_list%op('Sz')"
  a = my_list%op("Sz")
  print*,"a.print"
  call a%show

  print*,""

  print*,"TEST ITERATION SIZE:"
  do i=1,size(my_list)
     a = my_list%op(index=i)
     print*,i,my_list%key(i)
     call a%show()
  enddo


  print*,"TEST ITERATION KEYS:"
  keys = my_list%keys(len(keys))
  do i=1,size(keys)
     a = my_list%op(key=str(keys(i)))
     print*,i,str(keys(i))
     call a%show()
  enddo
end program testOPERATORS_TUPLE
