MODULE BLOCKS
  USE SCIFOR, only: str,assert_shape,zeye,eye,to_lower,free_unit,file_gzip,file_gunzip,set_store_size,t_start,t_stop,wait,check_MPI,get_Master_MPI
  USE AUX_FUNCS
  USE INPUT_VARS, only: block_file,umat_file,LOGfile
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
     character(len=:),allocatable                :: Opname
     character(len=:),allocatable                :: SiteType
   contains
     procedure,pass :: free        => free_block
     procedure,pass :: put_op      => put_op_block
     procedure,pass :: put_omat    => put_omat_block
     procedure,pass :: write_omat  => write_omat_block
     procedure,pass :: get_basis   => get_basis_block
     procedure,pass :: set_basis   => set_basis_block
     procedure,pass :: show        => show_block
     procedure,pass :: is_valid    => is_valid_block
     procedure,pass :: renormalize => rotate_operators_block
     procedure,pass :: okey        => okey_block
     procedure,pass :: name        => Opname_block
     procedure,pass :: type        => SiteType_block
     procedure,pass :: write       => write_block
     procedure,pass :: save        => save_block
     procedure,pass :: read        => read_block
     procedure,pass :: load        => load_block
     procedure,pass :: load_umat   => load_umat_block
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
    self%Dim    = 1
    call self%operators%free()
    call self%omatrices%free()
    if(allocated(self%sectors))then
       call self%sectors%free()
       deallocate(self%sectors)
    endif
    if(allocated(self%Opname))deallocate(self%Opname)
    if(allocated(self%SiteType))deallocate(self%SiteType)
  end subroutine free_block



  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function constructor_from_scrath(length,Dim,sectors,operators,omatrices,opname,SiteType) result(self)
    integer,intent(in)                       :: length
    integer,intent(in)                       :: Dim
    type(sectors_list),intent(in)            :: sectors(:)
    type(operators_list),intent(in)          :: operators
    type(operators_list),intent(in)  :: omatrices
    character(len=:),allocatable             :: OpName
    character(len=:),allocatable             :: SiteType
    type(block)                              :: self
    self%length    = length
    self%Dim       = Dim
    self%operators = operators
    self%omatrices = omatrices
    allocate(self%sectors(size(sectors)))
    do i=1,size(self%sectors)
       self%sectors(i) = sectors(i)
    enddo
    allocate(self%OpName, source=OpName)
    allocate(self%SiteType, source=SiteType)
  end function constructor_from_scrath


  function constructor_from_site(ssite) result(self)
    type(site),intent(in) :: ssite
    type(block)           :: self
    self%length    = 1
    self%Dim       = ssite%Dim
    self%operators = ssite%operators
#ifdef _CMPLX
    call self%omatrices%put("1",sparse(zeye(self%Dim)))
#else
    call self%omatrices%put("1",sparse(eye(self%Dim)))
#endif
    allocate(self%sectors(size(ssite%sectors)))
    do i=1,size(self%sectors)
       self%sectors(i)   = ssite%sectors(i)
    enddo
    allocate(self%OpName, source=ssite%OpName)
    allocate(self%SiteType, source=ssite%SiteType)
  end function constructor_from_site




  !##################################################################
  !##################################################################
  !       PUT  - GET/DUMP OPERATORS IN/FROM A LIST
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Load a sparse operator in the block dictionary
  !+------------------------------------------------------------------+
  subroutine put_op_block(self,key,op,type)
    class(block)                   :: self
    character(len=*),intent(in)    :: key
    type(sparse_matrix),intent(in) :: op
    character(len=*),intent(in)    :: type
    call self%operators%put(str(key),op,type)
  end subroutine put_op_block


  !+------------------------------------------------------------------+
  !PURPOSE:  Load a sparse matrix in the block dictionary
  !+------------------------------------------------------------------+
  subroutine put_omat_block(self,key,op,type)
    class(block)                   :: self
    character(len=*),intent(in)    :: key
    type(sparse_matrix),intent(in) :: op
    character(len=*),intent(in)    :: type
    call self%omatrices%put(str(key),op,type)
  end subroutine put_omat_block



  subroutine write_omat_block(self,key,op,type,suffix)
    class(block)                   :: self
    character(len=*),intent(in)    :: key
    type(sparse_matrix),intent(in) :: op
    character(len=*),intent(in)    :: type
    character(len=*),intent(in)    :: suffix
    integer :: unit
    unit = fopen(str(umat_file)//str(suffix),append=.true.)
    write(unit,*)str(key)
    if(str(type)=="")then
       write(unit,*)"none"
    else
       write(unit,*)str(type)
    endif
    call op%write(unit=unit)
    close(unit)
  end subroutine write_omat_block





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
  subroutine rotate_operators_block(self,Urho)
    class(block)                     :: self
    type(sparse_matrix)              :: Urho
    integer                          :: i,N,M  !N=self%dim,M=truncated dimension
    type(sparse_matrix)              :: Op
    character(len=:),allocatable     :: key,type
    !
    N = Urho%Nrow
    M = Urho%Ncol
    if(N/=self%dim) stop "self.renormalize error: size(Umat,1) != self.dim"
    do i=1,size(self%operators)
       key  = self%operators%key(index=i)
       type = self%operators%type(index=i)
       Op   = self%operators%op(index=i)
       call self%put_op(str(key),rotate_and_truncate(Op), type)
    enddo
    self%dim = M
    !
    call Op%free()
    !
  contains
    !
    !Udgr.O.U: [M,N].[N,N].[N,M]=[M,M]
    function rotate_and_truncate(Op) result(RotOp)
      type(sparse_matrix),intent(in) :: Op
      type(sparse_matrix)            :: RotOp
      if( any( [Op%Nrow,Op%Ncol] /= [N,N] ) ) &
           stop "self.renormalize error: shape(Op) != [N,N] N=size(Rho,1)"
      RotOp = matmul( matmul(Urho%dgr(),Op), Urho)
    end function rotate_and_truncate
    !
  end subroutine rotate_operators_block
  !



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
    allocate(A%OpName, source=B%OpName)
    allocate(A%SiteType, source=B%SiteType)
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
    integer,optional             :: iorb,isite,ispin
    integer                      :: iorb_,isite_,ispin_
    character(len=:),allocatable :: string
    iorb_ =0;if(present(iorb))iorb_=iorb
    ispin_=0;if(present(ispin))ispin_=ispin
    isite_=0;if(present(isite))isite_=isite
    !
    if(iorb_==0.AND.ispin_==0)stop "Okey_Block ERROR: iorb == ispin == 0"
    string = okey(iorb_,ispin_,isite_)
    !
    !
  end function okey_block


  function OpName_block(self) result(string)
    class(block)                  :: self
    character(len=:),allocatable :: string
    allocate(string, source=self%OpName)
  end function OpName_block


  function SiteType_block(self) result(string)
    class(block)                  :: self
    character(len=:),allocatable :: string
    string = to_lower(str(self%SiteType))
  end function SiteType_block



  !##################################################################
  !##################################################################
  !              SHOW 
  !##################################################################
  !##################################################################
  subroutine show_block(self,fmt,wOP,wOMAT,file)
    class(block)              :: self
    character(len=*),optional :: fmt
    logical,optional          :: wOP,wOMAT
    character(len=*),optional :: file
    character(len=32)         :: fmt_
    logical                   :: wOP_,wOMAT_
    integer                   :: unit_
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    wOP_  =.false.;if(present(wOP))  wOP_  =wOP
    wOMAT_=.false.;if(present(wOMAT))wOMAT_=wOMAT
    unit_=6;if(present(file))open(free_unit(unit_),file=str(file))
    !
    write(unit_,"(A,I6)")"Block Length  =",self%length
    write(unit_,"(A,I6)")"Block Dim     =",self%Dim
    write(unit_,"(A,A)") "Block Type    = ",self%SiteType
    write(unit_,"(A,I6)")"Block Sectors =",size(self%sectors)
    write(unit_,*)""
    do i=1,size(self%sectors)
       write(unit_,"(A15,I6)")"Block Sector  =",i
       call self%sectors(i)%show(unit=unit_)
    enddo
    if(wOP_)then
       write(unit_,"(A)")"Block Operators:"
       write(unit_,"(A,A)")"Op Name: ",self%OpName
       call self%operators%show(fmt=fmt_,unit=unit_)

    endif
    if(wOMAT_)then
       write(unit_,"(A)")"Block Omats:"
       call self%omatrices%show(fmt=fmt_,unit=unit_)
    endif
    if(present(file))close(unit_)
  end subroutine show_block





  subroutine write_block(self,suffix)
    class(block)              :: self
    character(len=*)          :: suffix
    integer                   :: Bunit,Uunit
    !
    open(free_unit(Bunit),file=str(block_file)//str(suffix) )
    !
    !General info:
    write(Bunit,*)self%length
    write(Bunit,*)self%Dim
    write(Bunit,*)self%SiteType
    write(Bunit,*)self%OpName
    !
    !Sectors:
    write(Bunit,*)size(self%sectors)
    !
    do i=1,size(self%sectors)
       call self%sectors(i)%write(unit=Bunit)
    enddo
    !
    !Operators
    call self%operators%write(unit=Bunit)
    !
    !Omatrices
    call self%omatrices%write(unit=Bunit)
    !
    close(Bunit)
    !
  end subroutine write_block




  subroutine save_block(self,suffix,gzip)
    class(block)     :: self
    character(len=*) :: suffix
    integer          :: unit_
    logical          :: gzip,fbool
    call self%write(str(suffix))
    if(gzip)then
       call file_gzip(str(block_file)//str(suffix))
       inquire(file=str(umat_file)//str(suffix), exist=fbool)
       if(fbool)call file_gzip(str(umat_file)//str(suffix))
    endif
  end subroutine save_block







  subroutine read_block(self,suffix)
    class(block)                   :: self
    character(len=*)               :: suffix
    integer                        :: Bunit,Uunit
    integer                        :: SectorSize
    integer                        :: length,il
    integer                        :: Dim
    type(sectors_list),allocatable :: sectors(:)
    type(operators_list)           :: operators
    type(operators_list)           :: omatrices
    type(sparse_matrix)            :: umat
    character(len=32)              :: OpName,key,type
    character(len=32)              :: SiteType
    logical                        :: fbool
    !
    open(free_unit(Bunit),file=str(block_file)//str(suffix) )
    !
    !General info:
    read(Bunit,*)length
    read(Bunit,*)Dim
    read(Bunit,*)SiteType
    read(Bunit,*)OpName
    !
    !Sectors:
    read(Bunit,*)SectorSize
    allocate(sectors(SectorSize))
    do i=1,SectorSize
       call sectors(i)%read(unit=Bunit)
    enddo
    !
    !Operators:
    call operators%read(unit=Bunit)
    !
    !Omatrices
    call omatrices%read(unit=Bunit)
    !if no Umat are included then it reads only the '1'.
    !try to read the rest from file. however this is enough
    !to restart basic calculation, not to measure.
    inquire(file=str(umat_file)//str(suffix), exist=fbool)
    if(fbool)then
       open(free_unit(Uunit),file=str(umat_file)//str(suffix))
       do il=2,length
          read(Uunit,*)key
          read(Uunit,*)type
          call umat%read(unit=Uunit)
          call omatrices%put(str(key),umat,str(type))
       enddo
       call umat%free()
       close(Uunit)
    endif
    !
    if(size(omatrices)/=length)then
       write(LOGfile,*)"read_block WARNING: could not read self.omatrices. No measurements will be possible."
    endif
    !
    !
    self = block(length,Dim,sectors,operators,omatrices,str(opname),str(SiteType))
    !
    close(Bunit)
    do i=1,SectorSize
       call sectors(i)%free()
    enddo
    call operators%free()
    call omatrices%free()
    !
  end subroutine read_block




  subroutine load_block(self,suffix)
    class(block)     :: self
    character(len=*) :: suffix
    logical          :: bool,gzbool
    !Check if block_file exists:
    inquire(file=str(block_file)//str(suffix), exist=bool)
    if(.not.bool)return
    write(LOGfile,*)"Loading from: "//str(block_file)//str(suffix)
    !
    !Check if umat_file exists:
    inquire(file=str(umat_file)//str(suffix), exist=bool)
    if(bool)write(LOGfile,*)"Loading from: "//str(umat_file)//str(suffix)
    !
    call self%read(str(suffix))
    !
  end subroutine load_block


  subroutine load_umat_block(self,suffix,length)
    class(block)                   :: self
    character(len=*)               :: suffix
    integer                        :: length
    integer                        :: il,Uunit
    type(sparse_matrix)            :: umat
    character(len=32)              :: OpName,key,type
    character(len=32)              :: SiteType
    logical                        :: fbool
    !to restart basic calculation, not to measure.
    inquire(file=str(umat_file)//str(suffix), exist=fbool)
    if(fbool)then
       open(free_unit(Uunit),file=str(umat_file)//str(suffix))
       do il=2,length
          read(Uunit,*)key
          read(Uunit,*)type
          call umat%read(unit=Uunit)
          call self%omatrices%put(str(key),umat,str(type))
       enddo
       call umat%free()
       close(Uunit)
    endif
  end subroutine load_umat_block


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
  USE SCIFOR
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  USE TUPLE_BASIS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
  implicit none
  type(block)                         :: my_block,a
  type(block),allocatable             :: my_blocks(:)
  type(operators_list)                :: op
  type(sectors_list)                  :: sect
  type(tbasis)                        :: sz_basis
  integer                             :: i
#ifdef _CMPLX
  complex(8),dimension(2,2),parameter :: Hzero=reshape([zero,zero,zero,zero],[2,2])
  complex(8),dimension(2,2),parameter :: S0=pauli_0
  complex(8),dimension(2,2),parameter :: Sz=pauli_z
  complex(8),dimension(2,2),parameter :: Sx=pauli_x
  complex(8),dimension(2,2),parameter :: Splus=reshape([zero,zero,one,zero],[2,2])
  complex(8),dimension(4,4)           :: Gamma13,Gamma03
#else
  real(8),dimension(2,2),parameter    :: Hzero=reshape([zero,zero,zero,zero],[2,2])
  real(8),dimension(2,2),parameter    :: S0=pauli_0
  real(8),dimension(2,2),parameter    :: Sz=pauli_z
  real(8),dimension(2,2),parameter    :: Sx=pauli_x
  real(8),dimension(2,2),parameter    :: Splus=reshape([zero,zero,one,zero],[2,2])
  real(8),dimension(4,4)              :: Gamma13,Gamma03
#endif

  Gamma13=kron(Sx,Sz)
  Gamma03=kron(S0,Sz)


  sz_basis = tbasis([0.5d0,-0.5d0],Qdim=1)

  print*,"TEST Constructor 1: from_scratch"
  my_block=block(&
       length=1, &       
       dim=2,&
       sectors=[sectors_list(sz_basis)],&
       operators=operators_list(['H0','Sz','Sp'],&
       [sparse(Hzero),sparse(Sz),sparse(Splus)],['b','s','b']),&
       opname='S',&
       sitetype='spin')
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
  my_block=block(spin_site(2))
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
  my_blocks(1)=block(spin_site(2))
  my_blocks(2)=block(spin_site(2))
  my_blocks(3)=block(spin_site(3))
  my_blocks(4)=block(spin_site(3))
  my_blocks(5)=block(spin_site(2))
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




  print*,"Check Write/Read blocks:"
  call a%free()
  a = block(spin_site(3))

  call a%write(file="block_a.dat")
  call a%show()

  print*,"Free memory..."
  call a%free()


  print*,"Read"
  call a%read(file="block_a.dat")
  print*,"Show"
  call a%show(wOP=.true.,wOMAT=.true.)



  print*,"Save"
  call set_store_size(0)
  call a%save(file="a_block.dat")



  print*,"Free memory..."
  call a%free()


  print*,"Load from file"
  call a%load('a_block.dat')
  print*,"show:"
  call a%show()

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










!   subroutine rotate_operators_block(self,Umat)
!     class(block)                     :: self
! #ifdef _CMPLX
!     complex(8),dimension(:,:)        :: Umat   ![N,M]
! #else
!     real(8),dimension(:,:)           :: Umat   ![N,M]
! #endif
!     integer                          :: i,N,M  !N=self%dim,M=truncated dimension
!     type(sparse_matrix)              :: Op
!     character(len=:),allocatable     :: key,type
!     !
!     N = size(Umat,1)
!     M = size(Umat,2)
!     if(N/=self%dim) stop "self.renormalize error: size(Umat,1) != self.dim"
!     do i=1,size(self%operators)
!        key  = self%operators%key(index=i)
!        type = self%operators%type(index=i)
!        Op   = self%operators%op(index=i)
!        call self%put_op(str(key),rotate_and_truncate(Op,Umat,N,M), type)
!     enddo
!     self%dim = M
!     !
!     call Op%free()
!     !
!   contains
!     !
!     !Udgr.rho.U [M,N].[N,N].[N,M]=[M,M]
!     function rotate_and_truncate(Op,trRho,N,M) result(RotOp)
!       type(sparse_matrix),intent(in) :: Op
!       integer                        :: N,M
!       type(sparse_matrix)            :: RotOp
! #ifdef _CMPLX
!       complex(8),dimension(N,M)      :: trRho
!       complex(8),dimension(M,M)      :: Umat
!       complex(8),dimension(N,N)      :: OpMat
! #else
!       real(8),dimension(N,M)         :: trRho
!       real(8),dimension(M,M)         :: Umat
!       real(8),dimension(N,N)         :: OpMat
! #endif
!       N = size(trRho,1)
!       M = size(trRho,2)
!       if( any( [Op%Nrow,Op%Ncol] /= [N,N] ) ) &
!            stop "self.renormalize error: shape(Op) != [N,N] N=size(Rho,1)"
!       OpMat= Op%as_matrix()
! #ifdef _CMPLX
!       Umat = matmul( matmul( conjg(transpose(trRho)),OpMat),trRho) 
! #else
!       Umat = matmul( matmul(transpose(trRho),OpMat),trRho) 
! #endif
!       call RotOp%load( Umat )
!     end function rotate_and_truncate
!     !
!   end subroutine rotate_operators_block
!   !
