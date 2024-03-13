MODULE BLOCKS
  USE SCIFOR, only: str,assert_shape,zeye,eye
  USE AUX_FUNCS
  USE MATRIX_SPARSE
  USE TUPLE_BASIS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  implicit none
  private


  type block
     integer                                     :: length=0
     integer                                     :: Dim=1
     type(sectors_list),dimension(:),allocatable :: sectors
     type(operators_list)                        :: operators
     type(operators_list)                        :: omatrices
   contains
     procedure,pass :: free        => free_block
     procedure,pass :: put_op      => put_op_block
     procedure,pass :: put_omat    => put_omat_block
     procedure,pass :: get_basis   => get_basis_block
     procedure,pass :: set_basis   => set_basis_block
     procedure,pass :: show        => show_block
     procedure,pass :: is_valid    => is_valid_block
     procedure,pass :: renormalize => rotate_operators_block
     procedure,pass :: okey        => okey_block
  end type block


  !GENERIC CONSTRUCTOR
  interface block
     module procedure :: constructor_from_scrath
     module procedure :: constructor_from_site
  end interface block

  !GENERIC CONSTRUCTOR
  interface as_block
     module procedure :: constructor_from_scrath
     module procedure :: constructor_from_site
  end interface as_block


  !EQUALITY 
  interface assignment(=)
     module procedure :: equality_block
  end interface assignment(=)


  public :: block
  public :: as_block
  public :: assignment(=)

  integer :: i,j

contains


  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Destructor
  !+------------------------------------------------------------------+
  subroutine free_block(self)
    class(block) :: self
    self%length = 0
    self%Dim   = 1
    call self%operators%free()
    if(allocated(self%sectors))then
       call self%sectors%free()
       deallocate(self%sectors)
    endif
  end subroutine free_block



  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function constructor_from_scrath(length,Dim,sectors,operators,omatrices) result(self)
    integer,intent(in)              :: length
    integer,intent(in)              :: Dim
    type(sectors_list),intent(in)   :: sectors(:)
    type(operators_list),intent(in) :: operators
    type(operators_list),intent(in) :: omatrices
    type(block)                     :: self
    self%length    = length
    self%Dim       = Dim
    self%operators = operators
    self%omatrices = omatrices
    allocate(self%sectors(size(sectors)))
    do i=1,size(self%sectors)
       self%sectors(i)   = sectors(i)
    enddo
  end function constructor_from_scrath


  function constructor_from_site(ssite) result(self)
    type(site),intent(in) :: ssite
    type(block)           :: self
    self%length    = 1
    self%Dim       = ssite%Dim
    self%operators = ssite%operators
    call self%omatrices%put("1",sparse(eye(self%Dim)))
    allocate(self%sectors(size(ssite%sectors)))
    do i=1,size(self%sectors)
       self%sectors(i)   = ssite%sectors(i)
    enddo
  end function constructor_from_site




  !##################################################################
  !##################################################################
  !       PUT  - GET/DUMP OPERATORS IN/FROM A LIST
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Load a sparse operator in the block dictionary
  !+------------------------------------------------------------------+
  subroutine put_op_block(self,key,op)
    class(block)                   :: self
    character(len=*),intent(in)    :: key
    type(sparse_matrix),intent(in) :: op
    call self%operators%put(str(key),op)
  end subroutine put_op_block


  !+------------------------------------------------------------------+
  !PURPOSE:  Load a sparse matrix in the block dictionary
  !+------------------------------------------------------------------+
  subroutine put_omat_block(self,key,op)
    class(block)                   :: self
    character(len=*),intent(in)    :: key
    type(sparse_matrix),intent(in) :: op
    call self%omatrices%put(str(key),op)
  end subroutine put_omat_block


  !+------------------------------------------------------------------+
  !PURPOSE:  Get Basis of the sector
  !+------------------------------------------------------------------+
  subroutine get_basis_block(self,basis,indx)
    class(block)     :: self
    type(tbasis)     :: basis
    integer,optional :: indx
    integer          :: indx_
    indx_=1;if(present(indx))indx_=indx
    if(indx_<1.OR.indx_>size(self%sectors))stop "SET_SECTORS_BLOCK ERROR: indx out of range"
    call basis%free()
    basis  = self%sectors(indx_)%basis()
  end subroutine get_basis_block


  !+------------------------------------------------------------------+
  !PURPOSE:  Put a QN array in the site
  !+------------------------------------------------------------------+
  subroutine set_basis_block(self,basis,indx)
    class(block)     :: self
    type(tbasis)     :: basis
    integer,optional :: indx
    integer          :: indx_
    indx_=1;if(present(indx))indx_=indx
    if(indx_<1.OR.indx_>size(self%sectors))stop "SET_SECTORS_BLOCK ERROR: indx out of range"
    self%sectors(indx_) = sectors_list( basis )
  end subroutine set_basis_block






  !+------------------------------------------------------------------+
  !PURPOSE:  
  !+------------------------------------------------------------------+
  subroutine rotate_operators_block(self,Umat)
    class(block)                 :: self
    real(8),dimension(:,:)    :: Umat   ![N,M]
    integer                      :: i,N,M  !N=self%dim,M=truncated dimension
    type(sparse_matrix)          :: Op
    character(len=:),allocatable :: key
    !
    N = size(Umat,1)
    M = size(Umat,2)
    if(N/=self%dim) stop "rotate_operators_block error: size(Umat,1) != self.dim"
    do i=1,size(self%operators)
       key = self%operators%key(index=i)
       Op  = self%operators%op(index=i)
       call self%put_op(str(key),rotate_and_truncate(Op,Umat,N,M))
    enddo
    self%dim = M
    !
    call Op%free()
    !
  end subroutine rotate_operators_block
  !
  function rotate_and_truncate(Op,trRho,N,M) result(RotOp)
    type(sparse_matrix),intent(in) :: Op
    real(8),dimension(N,M)      :: trRho  ![Nesys,M]
    integer                        :: N,M
    type(sparse_matrix)            :: RotOp
    real(8),dimension(M,M)      :: Umat
    real(8),dimension(N,N)      :: OpMat
    N = size(trRho,1)
    M = size(trRho,2)
    if( any( [Op%Nrow,Op%Ncol] /= [N,N] ) ) stop "rotate_and_truncate error: shape(Op) != [N,N] N=size(Rho,1)"
    OpMat= Op%as_matrix()
    ! Umat = matmul( conjg(transpose(trRho)), matmul(OpMat,trRho)) ![M,N].[N,N].[N,M]=[M,M]
    Umat = matmul( (transpose(trRho)), matmul(OpMat,trRho)) ![M,N].[N,N].[N,M]=[M,M]
    call RotOp%load( Umat )
  end function rotate_and_truncate


  !##################################################################
  !##################################################################
  !              OPERATIONS / ASSIGNEMENTS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Equality between two blocks
  !+------------------------------------------------------------------+
  subroutine equality_block(A,B)
    type(block),intent(inout) :: A
    type(block),intent(in)    :: B
    call A%free
    A%length    = B%length
    A%Dim       = B%Dim
    A%operators = B%operators
    A%omatrices = B%omatrices
    allocate(A%sectors(size(B%sectors)))
    do i=1,size(A%sectors)
       A%sectors(i) = B%sectors(i)
    enddo
  end subroutine equality_block


  function is_valid_block(self,nobasis) result(bool)
    class(block)                          :: self
    logical,optional                      :: nobasis    
    logical                               :: bool
    logical                               :: nobasis_
    integer,dimension(size(self%sectors)) :: Dims
    !
    integer                               :: i
    type(sparse_matrix)                   :: ope
    type(operators_list)                  :: op

    nobasis_=.false.;if(present(nobasis))nobasis_=nobasis
    !
    bool = self%operators%is_valid(self%Dim)
    op = self%operators
    do i=1,size(op)
       ope = op%op(index=i)
    enddo
    if(nobasis_)return
    do i=1,size(self%sectors)
       Dims(i) = dim(self%sectors(i))
    enddo
    bool=bool.AND.(self%dim==product(Dims))
  end function is_valid_block


  function okey_block(self,iorb,ispin,isite) result(string)
    class(block)                 :: self
    integer                      :: iorb
    integer,optional             :: ispin,isite
    integer                      :: isite_,ispin_
    character(len=:),allocatable :: string
    ispin_=0;if(present(ispin))ispin_=ispin
    isite_=0;if(present(isite))isite_=isite
    !
    string = okey(iorb,ispin_,isite_)
    !
  end function okey_block




  !##################################################################
  !##################################################################
  !              SHOW 
  !##################################################################
  !##################################################################
  subroutine show_block(self,fmt,wOP,wOMAT)
    class(block)              :: self
    character(len=*),optional :: fmt
    logical,optional          :: wOP,wOMAT
    character(len=32)         :: fmt_
    logical :: wOP_,wOMAT_
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    wOP_  =.false.;if(present(wOP))  wOP_  =wOP
    wOMAT_=.false.;if(present(wOMAT))wOMAT_=wOMAT
    write(*,"(A15,I6)")"Block Length  =",self%length
    write(*,"(A15,I6)")"Block Dim     =",self%Dim
    write(*,"(A15,I6)")"Block Sectors =",size(self%sectors)
    do i=1,size(self%sectors)
       write(*,"(A14,I6)")"Block Sector  =",i
       call self%sectors(i)%show()
    enddo
    if(wOP_)then
       write(*,"(A14)")"Block Ops     ="
       call self%operators%show(fmt=fmt_)
    endif
    if(wOMAT_)then
       write(*,"(A14)")"Block Omats   ="
       call self%omatrices%show(fmt=fmt_)
    endif
  end subroutine show_block







END MODULE BLOCKS






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
program testBLOCKS
  USE SCIFOR,  id => eye
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  USE TUPLE_BASIS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
  implicit none


  type(block)                      :: my_block,a
  type(block),allocatable          :: my_blocks(:)
  type(operators_list)             :: op
  type(sectors_list)               :: sect
  type(tbasis)                     :: sz_basis
  integer                          :: i
  real(8),dimension(2,2),parameter   :: Hzero=reshape([zero,zero,zero,zero],[2,2])
  real(8),dimension(2,2),parameter   :: Sz=pauli_z
  real(8),dimension(2,2),parameter   :: Sx=pauli_x
  real(8),dimension(2,2),parameter   :: Splus=reshape([zero,zero,one,zero],[2,2])
  real(8),dimension(4,4)             :: Gamma13,Gamma03

  
  Gamma13=kron(Sx,Sz)
  Gamma03=kron(id(2),Sz)


  sz_basis = tbasis([0.5d0,-0.5d0],Qdim=1)

  print*,"TEST Constructor 1: from_scratch"
  my_block=block(&
       length=1, &       
       dim=2,&
       sectors=[sectors_list(sz_basis)],&
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


  print*,"Show op= a.operators:"
  call op%show



  print*,"Check allocatable array of blocks:"
  allocate(my_blocks(5))
  print*,"Copy spin 1/2 into BlockArray(1:5)"
  my_blocks(1)=block(spin_onehalf_site())
  my_blocks(2)=block(spin_onehalf_site())
  my_blocks(3)=block(spin_one_site())
  my_blocks(4)=block(spin_one_site())
  my_blocks(5)=block(spin_onehalf_site())
  !
  print*,"Check each of the block is correct" 
  do i=2,5
     print*,i-1,my_blocks(i-1)%is_valid()
     print*,i,my_blocks(i)%is_valid()
     print*,""
  enddo
  a = my_blocks(5);print*,a%is_valid()
  a = my_blocks(4);print*,a%is_valid()
  a = my_blocks(3);print*,a%is_valid()
  a = my_blocks(2);print*,a%is_valid()
  a = my_blocks(1);print*,a%is_valid()


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
#endif





