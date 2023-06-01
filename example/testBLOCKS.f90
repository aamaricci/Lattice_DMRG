program testBLOCKS
  USE SCIFOR,  id => zeye
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
  implicit none


  type(block)          :: my_block,a
  type(operators_list) :: op
  type(sectors_list)   :: sect
  integer              :: i
  real(8),dimension(2,2),parameter   :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
  real(8),dimension(2,2),parameter   :: Sz=dble(pauli_z)
  real(8),dimension(2,2),parameter   :: Sx=dble(pauli_x)
  real(8),dimension(2,2),parameter   :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
  real(8),dimension(4,4)             :: Gamma13,Gamma03

  Gamma13=kron(Sx,Sz)
  Gamma03=kron(eye(2),Sz)



  print*,"TEST Constructor 1: from_scratch"
  my_block=block(&
       length=1, &       
       dim=2,&
       sectors=[sectors_list([-0.5d0,0.5d0])],&
       operators=operators_list(['H0','Sz','Sp'],&
       [sparse(Hzero),sparse(Sz),sparse(Splus)]))
  print*,"Showing the operator list:"
  call my_block%show()
  print*,""


  print*,"Check my_block is valid"
  print*,my_block%is_valid()
  print*,""

  print*,"Free the block"
  call my_block%free()
  print*,""





  print*,"TEST Constructor 2: from_site"
  my_block=block(spin_onehalf_site())
  print*,"Showing the operator list:"
  call my_block%show()
  print*,""

  print*,"Check my_block is valid"
  print*,my_block%is_valid()
  print*,""

  print*,"Test equality; a.show()"
  a = my_block



  print*,"a is valid:",a%is_valid()
  call a%show()
  print*,""


  ! print*,"Free the block"
  ! call my_block%free()
  ! print*,""


  print*,"Test retrieve Sector: sect=a.sectors"
  sect = a%sectors(1)
  call sect%show()
  print*,""

  print*,"Get operator list to op:"
  op = a%operators
  print*,"Showing it:"
  call op%show()
  print*,""

  call a%load("G",Gamma13)
  print*,"Show Block a:"
  call a%show()
  print*,""

  print*,"Show op= a.operators:"
  call op%show



contains


  subroutine i_random(A)
    integer,dimension(:) :: A
    integer :: i1,i2
    do i1=1,size(A,1)
       A(i1)=mt_uniform(1,10)
    enddo
  end subroutine i_random


  subroutine print_vec(M)
    integer,dimension(:) :: M
    integer :: i,j
    do i=1,size(M,1)
       write(*,"(I3)")M(i)
    enddo
  end subroutine print_vec

  subroutine print_mat(M)
    integer,dimension(:,:) :: M
    integer :: i,j
    do i=1,size(M,1)
       write(*,"("//str(size(M,2))//"(I3,1x))")(M(i,j),j=1,size(M,2))
    enddo
  end subroutine print_mat


end program testBLOCKS
