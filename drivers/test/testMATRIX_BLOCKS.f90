program testBLOCK_MATRICES
  USE MATRIX_BLOCKS
  USE SCIFOR
  implicit none


  type(blocks_matrix)                :: blH,a,adg
  real(8),dimension(:,:),allocatable :: Matrix
  real(8),dimension(:),allocatable   :: Vec
  real(8),dimension(:),allocatable   :: evals
  integer,dimension(:),allocatable   :: eorder,Dq
  integer                            :: i,j,q,N,count
  real(8),dimension(2,2),parameter   :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
  real(8),dimension(2,2),parameter   :: Sz=dble(pauli_z)
  real(8),dimension(2,2),parameter   :: Sx=dble(pauli_x)
  real(8),dimension(2,2),parameter   :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
  real(8),dimension(4,4)             :: Gamma13,Gamma03

  Gamma13=kron(Sx,Sz)
  Gamma03=kron(eye(2),Sz)

  

  print*,"test CONSTRUCTOR 1: block_matrix(matrix)"
  a = blocks_matrix(Sz,qn=0d0)
  call a%show(dble=.true.)
  print*,shape(a)
  call a%free()
  print*,""


  print*,"test CONSTRUCTOR 2: load(matrix)"
  call blH%load(Sx,qn=1d0)
  call blH%show(dble=.true.)
  call blH%free()
  print*,""

  print*,"test CONSTRUCTOR 3: append"
  call a%append(eye(2),qn=1d0)
  call a%append(kron(eye(2),Sz),qn=2d0)
  call a%show(dble=.true.)
  call a%free()
  print*,""

  print*,"test APPEND:"
  call a%load(eye(2),qn=0d0)
  call a%append(kron(eye(2),Sz),qn=1d0)
  call a%show(dble=.true.)
  call a%free()
  print*,""




  print*,"test GET BLOCK MATRIX:"
  call a%load(eye(2),qn=0d0)
  call a%append(kron(Sz,Sx),qn=1d0)
  call a%show(dble=.true.)
  matrix = a%get(index=2)
  do i=1,4
     write(*,"(100F5.2)")(matrix(i,j),j=1,4)
  enddo
  print*,""


  print*,"test DUMP WHOLE MATRIX:"
  call a%append(Sz,qn=0.5d0)
  matrix = a%dump()
  print*,shape(matrix),shape(a)
  do i=1,size(matrix,1)
     write(*,"(100F5.2)")(matrix(i,j),j=1,size(matrix,2))
  enddo
  call a%free()
  print*,""



  print*,"test EQUALITY:"
  call a%load(Sx,qn=0.5d0)
  call a%append(kron(eye(2),Sz),qn=1d0)
  call a%show(dble=.true.)
  blH = a
  call blH%show(dble=.true.)
  call a%free()
  call blH%free()
  print*,""





  print*,"test PUSH:"
  a = as_blocks(eye(2),qn=0.5d0)
  call a%show()  
  call a%push(kron(Sx,Sz),qn=1d0)
  call a%push(Sx,qn=2d0)
  call a%show()
  print*,"test PUSH as UPDATE:"
  call a%push(Sz,qn=1d0)
  call a%show()
  call a%free()
  print*,""

  print*,"test DGR:"
  a = blocks_matrix(reshape([0d0,0d0,1d0,0d0],[2,2]),qn=1d0)
  call a%append(kron(Splus,Sz),qn=0.5d0)
  call a%show()
  adg = a%dgr()
  ! adg = transpose(a)
  call adg%show()
  call a%free()
  call adg%free()
  print*,""


  print*,"test EIGH:"
  call a%append(mersenne()*eye(2),qn=0d0)
  call a%append(mersenne()*Sx,qn=1d0)
  call a%append(mersenne()*Sz,qn=2d0)
  call a%append(mersenne()*kron(Sz,Sx),qn=3d0)
  call a%show()
  print*,""
  call a%eigh(reverse=.false.)
  print*,"E,V:"
  call a%show(dble=.true.)
  print*,""
  evals = a%evals()
  write(*,"(A3,100I12)")"O:",arange(1,size(evals))
  write(*,"(A3,100F12.5)")"E:",evals
  print*,""



  ! allocate(eorder(size(evals)))
  ! call sort_array(evals,eorder)
  ! write(*,"(A3,100I12)")"O:",Eorder
  ! write(*,"(A3,100F12.5)")"E:",evals
  ! deallocate(evals,eorder)

  evals = a%evals(sort=.true.,reverse=.true.,order=eorder)
  write(*,"(A3,100I12)")"O:",Eorder
  write(*,"(A3,100F12.5)")"E:",evals
  print*,""
  print*,""

  ! vec= a%evec(index=2,pos=2)
  N = size(evals)
  do i=1,N
     call block_find_indices(a,eorder(i),q,j)
     vec = a%evec(index=q,pos=j)
     write(*,"(I3,A2,I3,A10,I4,A10,I4,A1,100F12.5)")i,"->",eorder(i),"Block:",q,"Position:",j,":",dble(vec)
  enddo




  print*,""
  call a%free()
  call a%append(mersenne()*eye(2),qn=0d0)
  call a%append(mersenne()*Sx,qn=1d0)
  call a%append(mersenne()*Sz,qn=2d0)
  call a%append(mersenne()*kron(Sz,Sx),qn=3d0)
  call a%show(dble=.true.)  
  ! print*,a%dims()
  Dq = a%dims()
  N = size(Dq)
  ![1....,D1][D1+1,....,D1+D2][...][sum(D_1:D_{q-1})+1,...,sum(D_1:D_{q-1})+Dq]
  do q=1,N
     do j=1,Dq(q)
        i = j
        if(q>1)i = j + sum(Dq(1:q-1))
        print*,i,q,j
     enddo
  enddo
  print*,""

  do i=1,sum(Dq)
     do q=1,N
        if(i <= sum(Dq(1:q)))then
           j = i - sum(Dq(1:q-1))
           exit
        endif
     enddo
     print*,i,q,j
  enddo



contains

  subroutine block_find_indices(self,in,q,pos)
    type(blocks_matrix)              :: self
    integer                          :: in
    integer                          :: q
    integer                          :: pos
    integer,dimension(:),allocatable :: Dq
    Dq = self%dims()            !retrieve dimensions of all blocks
    do q=1,size(Dq)
       if(in <= sum(Dq(1:q)))then
          pos = in - sum(Dq(1:q-1))
          exit
       endif
    enddo
  end subroutine block_find_indices


end program testBLOCK_MATRICES






! ![1....,D1][D1+1,....,D1+D2][...][sum(D_1:D_{q-1})+1,...,sum(D_1:D_{q-1})+Dq]
! do q=1,N
!    do j=1,Dq(q)
!       i = j
!       if(q>1)i = j + sum(Dq(1:q-1))
!       print*,i,q,j
!    enddo
! enddo
! print*,""

! do i=1,sum(Dq)
!    do q=1,N
!       if(i <= sum(Dq(1:q)))then
!          j = i - sum(Dq(1:q-1))
!          exit
!       endif
!    enddo
!    print*,i,q,j
! enddo
