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
  print*,len(a)


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

end program testLIST_SECTORS
