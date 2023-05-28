program testSPARSE_MATRICES
  USE MATRIX_SPARSE
  USE SCIFOR
  implicit none


  integer                            :: i,j
  type(sparse_matrix)                :: spH,spK,a,b,c,avec(2)
  real(8),dimension(4,4)             :: GammaX
  real(8),dimension(:,:),allocatable :: Amat,Bmat,Cmat

  real(8),dimension(2,2),parameter   :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
  real(8),dimension(2,2),parameter   :: Sz=dble(pauli_z)
  real(8),dimension(2,2),parameter   :: Sx=dble(pauli_x)
  real(8),dimension(2,2),parameter   :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
  real(8),dimension(4,4)             :: Gamma13,Gamma03

  Gamma13=kron(Sx,Sz)
  Gamma03=kron(eye(2),Sz)

  print*,"test INIT"
  call spH%init(2,2)
  print*,shape(spH)
  print*,all(shape(spH) == [2,2])
  print*,all(shape(spH) == [3,3])
  print*,""


  print*,"test CONSTRUCTOR 1: sparse_matrix(matrix)"
  a = sparse_matrix(Sz)
  call a%show()
  print*,shape(a)
  print*,all(shape(a) == [2,2])
  print*,all(shape(a) == [3,3])
  print*,"a.NNZ=",a%nnz()
  call a%free()
  print*,""


  print*,"test CONSTRUCTOR 2: as_sparse(pauli_x*pauli_z)"
  a = as_sparse(Gamma13)
  call a%show()
  print*,shape(a)
  print*,all(shape(a) == [2,2])
  print*,all(shape(a) == [4,4])
  print*,"a.NNZ=",a%nnz()
  print*,""
  call a%free()



  print*,"test CONSTRUCTOR 3: avec(1:2)= [as_sparse(pauli_x),as_sparse(pauli_z)]"  
  avec = [sparse(Sx),as_sparse(Sz)]
  call avec(1)%show()
  call avec(2)%show()
  print*,shape(avec(1))
  call avec%free()
  print*,""


  print*,"test FREE"
  call spH%free()
  print*,""



  print*,"test LOAD and PRINT"
  call spH%load(dble(kron(pauli_0,pauli_z)))
  call spH%show()
  print*,"spH.NNZ=",spH%nnz()
  print*,""

  print*,"test GET ELEMENT"
  write(*,*)"spH(2,2)=",spH%get(2,2)
  write(*,*)"spH(3,4)=",spH%get(3,4)
  print*,""


  print*,"test INSERT ELEMENT"
  call spH%insert(1d0,1,4)
  call spH%insert(-1d0,4,4)
  call spH%show()  
  print*,""


  print*,"test SPY"
  call spH%spy("spH")
  print*,""


  print*,"test DUMP"
  gammaX=0d0
  do i=1,4
     write(*,*)(gammaX(i,j),j=1,4)
  enddo
  print*,""
  call spH%dump(gammaX)
  do i=1,4
     write(*,*)(gammaX(i,j),j=1,4)
  enddo

  gammaX=0d0
  do i=1,4
     write(*,*)(gammaX(i,j),j=1,4)
  enddo
  print*,""
  gammaX = spH%as_matrix()
  call spH%show()
  do i=1,4
     write(*,*)(gammaX(i,j),j=1,4)
  enddo
  print*,""


  print*,"test spK=spH"  
  spK=spH
  call spK%show()


  print*,"test spH=zero"  
  spH=0d0
  call spH%show()
  spH=spK



  print*,"test ADDITION a+b=c"
  print*,"a=sigma_0"
  call a%init(2,2)
  call a%load(eye(2))
  call a%show()

  print*,"b=sigma_X"  
  call b%init(2,2)
  call b%load(Sx)
  call b%show()

  print*,"c=sigma_0 + sigma_X"
  c = a+b
  call c%show()

  call a%free()
  call b%free()
  call c%free()



  print*,"test SUBTRACTION a-b=c"
  print*,"a=sigma_0"
  call a%init(2,2)
  call a%load(eye(2))
  call a%show()

  print*,"b=sigma_Z"  
  call b%init(2,2)
  call b%load(Sz)
  call b%show()


  print*,"c=sigma_0 - sigma_Z"  
  c = a-b
  call c%show()


  call a%free()
  call b%free()
  call c%free()



  print*,"test LEFT SCALAR PRODUCT b=a*const"
  print*,"a=sigma_0"
  call a%init(2,2)
  call a%load(eye(2))
  call a%show()

  print*,"b=2*a"  
  b = 2*a
  call b%show()

  print*,"b=2d0*a"  
  b = 2d0*a
  call b%show()

  print*,"test RIGHT SCALAR PRODUCT b=const*a"
  print*,"b=a*2"  
  b = a*2
  call b%show()

  print*,"b=a*2d0"  
  b = a*2d0
  call b%show()


  print*,"test RIGHT SCALAR DIVISDION b=a/const"
  print*,"b=a/2"  
  b = a/2
  call b%show()

  print*,"b=a/2d0"  
  b = a/2d0
  call b%show()




  print*,"test KRON PRODUCT 1"
  call a%free()
  call b%free()
  call c%free()
  call spH%free()
  !
  call spH%load(gamma13)
  call a%load(Sz)
  call b%load(Sx)

  print*,"a=sigma_0"  
  call a%show()  
  print*,"b=sigma_Z"  
  call b%show()  
  print*,"c=a.x.b"  
  c = a.x.b
  call c%show()
  print*,"spH=sigma_0xsigma_Z"  
  call spH%show()
  print*,""


  print*,""
  print*,"test KRON PRODUCT 2"
  allocate(Amat(2,2),Bmat(2,2))
  allocate(Cmat(4,4))
  Amat = dble(transpose(reshape([1,2,3,4],[2,2])))
  Bmat = dble(transpose(reshape([0,5,6,7],[2,2])))
  Cmat = dble(transpose(reshape([0,5,0,10,5,7,12,14,0,15,0,20,18,21,24,28],[4,4])))
  call a%load(Amat)
  call b%load(Bmat)
  print*,"A = 1 2    B = 0 5"
  print*,"    3 4        6 7"
  call a%show()  
  call b%show()  

  print*,"c=a.x.b"  
  c = a.x.b
  call c%show()
  print*,"C = 0  5  0  10"
  print*,"    6  7  12 14"
  print*,"    0  15 0  20"
  print*,"    18 21 24 28"

  call a%free()
  call b%free()
  call c%free()
  call spH%free()
  deallocate(Amat,Bmat,Cmat)
  print*,""




  print*,""
  print*,"test KRON PRODUCT 3"
  allocate(Amat(3,2),Bmat(2,3))
  Amat = dble(transpose(reshape([1,2,3,4,1,0],[2,3])))
  Bmat = dble(transpose(reshape([0,5,2,6,7,3],[3,2])))
  call a%load(Amat)
  call b%load(Bmat)
  !
  print*," A = 1 2    B = 0 5 2"
  print*,"     3 4        6 7 3"
  print*,"     1 0             "
  call a%show()  
  call b%show()
  !

  allocate(Cmat(6,6))
  Cmat = dble(transpose(reshape([&
       0,5,2,0,10,4, &
       6,7,3,12,14,6,&
       0,15,6,0,20,8,&
       18,21,9,24,28,12,&
       0,5,2,0,0,0,&
       6,7,3,0,0,0],&
       [4,4])))

  print*,"c=a.x.b"  
  c = a.x.b
  call c%show()
  print*,"C = 0      5    2    0     10    4"
  print*,"    6      7    3   12     14    6"
  print*,"    0     15    6    0     20    8"    
  print*,"    18     21    9   24     28   12"    
  print*,"    0      5    2    0      0    0"    
  print*,"    6      7    3    0      0    0"




  call a%free()
  call b%free()
  call c%free()
  call spH%free()

  print*, "TEST TRANSPOSE CONJUGATE"
  call a%load(Sx+Sz)
  call a%show()
  b = transpose(a)
  call b%show()

  if(any( a%as_matrix()-b%as_matrix() /= 0d0) )then
     write(*,*)"Wrong TRANSPOSE"
  else
     write(*,*)"Good TRANSPOSE"
  endif





  print*, "TEST TRANSPOSE CONJUGATE"
  call a%load(Splus)
  call a%show()
  b = a%t()
  call b%show()

  if(any( a%as_matrix()-b%as_matrix() /= 0d0) )then
     write(*,*)"Wrong TRANSPOSE"
  else
     write(*,*)"Good TRANSPOSE"
  endif


  deallocate(Amat,Bmat,Cmat)

  print*,""
  print*,"test KRON PRODUCT 3"

  allocate(Amat(5,5));Amat=0d0
  Amat(1,2) = 1d0
  do i=2,5-1
     Amat(i,i-1) = 1d0
     Amat(i,i+1) = 1d0    
  enddo
  Amat(5,5-1) = 1d0

  allocate(Bmat(5,5))
  Bmat = dble((reshape([1,0,1,0,1,1,0,1,0,1,1,0,1,0,1,1,0,1,0,1,1,0,1,0,1],[5,5])))

  call a%load(Amat)
  call b%load(Bmat)
  !
  print*,"A"
  call a%show()
  print*,"B"
  call b%show()
  print*,""


  allocate(Cmat(5,5))
  Cmat = matmul(Amat,Bmat)
  do i=1,5
     write(*,"(5F9.3,1x)")(Cmat(i,j),j=1,5)
  enddo

  print*,""
  c = a.m.b
  call c%show()
  print*,c%nnz()
end program testSPARSE_MATRICES
