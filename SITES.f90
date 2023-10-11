MODULE SITES
  USE SCIFOR, only: str,diag,zeros,eye,kron,pauli_x,pauli_y,pauli_z,xi,operator(.kx.)
  USE INPUT_VARS
  USE AUX_FUNCS
  USE MATRIX_SPARSE
  USE TUPLE_BASIS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE HLOCAL
  implicit none
  private


  type site
     integer                                     :: Dim=1
     type(sectors_list),dimension(:),allocatable :: sectors
     type(operators_list)                        :: operators
   contains
     procedure,pass     :: free        => free_site
     procedure,pass     :: put         => put_op_site
     procedure,pass     :: load        => load_op_site
     procedure,pass     :: get_basis   => get_basis_site
     procedure,pass     :: set_basis   => set_basis_site
     procedure,pass     :: show        => show_site
     procedure,pass     :: is_valid    => is_valid_site
  end type site


  !GENERIC CONSTRUCTOR
  interface site
     module procedure :: constructor_site
  end interface site

  !EQUALITY 
  interface assignment(=)
     module procedure :: equality_site
  end interface assignment(=)


  public :: site
  public :: pauli_site
  public :: spin_onehalf_site
  public :: spin_one_site
  public :: hubbard_site
  public :: assignment(=)

  integer :: i

contains


  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Destructor
  !+------------------------------------------------------------------+
  subroutine free_site(self)
    class(site) :: self
    self%dim   = 1
    call self%operators%free()
    if(allocated(self%sectors))then
       call self%sectors%free()
       deallocate(self%sectors)
    endif
  end subroutine free_site



  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function constructor_site(Dim,sectors,operators) result(self)
    integer,intent(in)              :: Dim
    type(sectors_list),intent(in)   :: sectors(:)
    type(operators_list),intent(in) :: operators
    type(site)                      :: self
    self%Dim       = Dim
    self%operators = operators
    allocate(self%sectors(size(sectors)))
    do i=1,size(self%sectors)
       self%sectors(i)   = sectors(i)
    enddo
  end function constructor_site








  !##################################################################
  !##################################################################
  !       PUT/LOAD
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Put a sparse operator in the site dictionary
  !+------------------------------------------------------------------+
  subroutine put_op_site(self,key,op)
    class(site)                    :: self
    character(len=*),intent(in)    :: key
    type(sparse_matrix),intent(in) :: op
    call self%operators%put(str(key),op)
  end subroutine put_op_site



  !+------------------------------------------------------------------+
  !PURPOSE:  Load a dense operator in the site dictionary
  !+------------------------------------------------------------------+
  subroutine load_op_site(self,key,op)
    class(site)                       :: self
    character(len=*),intent(in)       :: key
    real(8),dimension(:,:),intent(in) :: op
    call self%operators%load(str(key),op)
  end subroutine load_op_site


  !+------------------------------------------------------------------+
  !PURPOSE:  Get Basis of the sector
  !+------------------------------------------------------------------+
  subroutine get_basis_site(self,basis,indx)
    class(site)      :: self
    type(tbasis)     :: basis
    integer,optional :: indx
    integer          :: indx_
    indx_=1;if(present(indx))indx_=indx
    if(indx_<1.OR.indx_>size(self%sectors))stop "SET_SECTORS_BLOCK ERROR: indx out of range"
    call basis%free()
    basis  = self%sectors(indx_)%basis()
  end subroutine get_basis_site


  !+------------------------------------------------------------------+
  !PURPOSE:  Put a QN array in the site
  !+------------------------------------------------------------------+
  subroutine set_basis_site(self,basis,indx)
    class(site)      :: self
    type(tbasis)     :: basis
    integer,optional :: indx
    integer          :: indx_
    indx_=1;if(present(indx))indx_=indx
    if(indx_<1.OR.indx_>size(self%sectors))stop "SET_SECTORS_SITE ERROR: indx out of range"
    self%sectors(indx) = sectors_list( basis )
  end subroutine set_basis_site



  !##################################################################
  !##################################################################
  !              OPERATIONS / ASSIGNEMENTS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Equality between two sites
  !+------------------------------------------------------------------+
  subroutine equality_site(A,B)
    type(site),intent(inout) :: A
    type(site),intent(in)    :: B
    call A%free
    A%Dim       = B%Dim
    A%operators = B%operators
    allocate(A%sectors(size(B%sectors)))
    do i=1,size(A%sectors)
       A%sectors(i)   = B%sectors(i)
    enddo
  end subroutine equality_site


  function is_valid_site(self) result(bool)
    class(site)                           :: self
    logical                               :: bool
    integer,dimension(size(self%sectors)) :: Dims
    bool = self%operators%is_valid(self%Dim)
    do i=1,size(self%sectors)
       Dims = dim(self%sectors(i))
    enddo
    bool=bool.AND.(self%dim==product(Dims))
  end function is_valid_site





  !##################################################################
  !##################################################################
  !              SHOW 
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Pretty print the site structure
  !+------------------------------------------------------------------+
  subroutine show_site(self,fmt)
    class(site)               :: self
    character(len=*),optional :: fmt
    character(len=32)         :: fmt_
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    write(*,"(A14,I6)")"Site Dim     =",self%Dim
    do i=1,size(self%sectors)
       write(*,"(A14,I6)")"Site Sectors =",i
       call self%sectors(i)%show()
    enddo
    write(*,"(A14)")"Site Ops     ="
    call self%operators%show(fmt=fmt_)
  end subroutine show_site





  !##################################################################
  !##################################################################
  !              PRESET
  !##################################################################
  !##################################################################
  function pauli_site(hvec) result(self)
    real(8),dimension(3),optional :: hvec
    type(site)                    :: self
    real(8),dimension(2,2)        :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
    real(8),dimension(2,2)        :: Szeta=reshape([0.5d0,0d0,0d0,0.5d0],[2,2])
    real(8),dimension(2,2)        :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
    call self%free()
    self%Dim = 2
    if(present(hvec))then
       Hzero = Hzero+hvec(1)*pauli_x+hvec(2)*pauli_z
    endif
    call self%put("H",sparse(Hzero))
    call self%put("Sz",sparse(Szeta))
    call self%put("Sp",sparse(Splus))
    allocate(self%sectors(1))
    self%sectors(1) = sectors_list( tbasis([0.5d0,-0.5d0],qdim=1) )
  end function pauli_site


  function spin_onehalf_site(hvec) result(self)
    real(8),dimension(3),optional :: hvec
    type(site)                    :: self
    real(8),dimension(2,2)        :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
    real(8),dimension(2,2)        :: Szeta=reshape([0.5d0,0d0,0d0,0.5d0],[2,2])
    real(8),dimension(2,2)        :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
    call self%free()
    self%Dim = 2
    if(present(hvec))then
       Hzero = Hzero+hvec(1)*pauli_x+hvec(2)*pauli_z
    endif
    call self%put("H",sparse(Hzero))
    call self%put("Sz",sparse(Szeta))
    call self%put("Sp",sparse(Splus))
    allocate(self%sectors(1))
    self%sectors(1) = sectors_list( tbasis([0.5d0,-0.5d0],qdim=1) )
  end function spin_onehalf_site


  function spin_one_site(hvec) result(self)
    real(8),dimension(3),optional :: hvec
    type(site)                        :: self
    real(8),dimension(3,3)            :: Hzero
    real(8),dimension(3,3)            :: Szeta
    real(8),dimension(3,3)            :: Splus
    real(8),dimension(3,3)            :: Sx=&
         reshape([0d0,1d0,0d0,1d0,0d0,1d0,0d0,1d0,0d0],[3,3])/sqrt(2d0)
    Hzero = 0d0
    Szeta = diag([1d0,0d0,-1d0])
    if(present(hvec))then
       Hzero = Hzero+hvec(1)*Sx+hvec(2)*Szeta
    endif
    !
    Splus = 0d0
    Splus(1,2) = sqrt(2d0)
    Splus(2,3) = sqrt(2d0)
    !
    call self%free()
    self%Dim = 3
    call self%put("H",sparse(Hzero))
    call self%put("Sz",sparse(Szeta))
    call self%put("Sp",sparse(Splus))
    allocate(self%sectors(1))
    self%sectors(1) = sectors_list( tbasis([1d0,0d0,-1d0],qdim=1) )
  end function spin_one_site


  !Fock = [|0,0>,|1,0>,|0,1>,|1,1>] <- |up,dw> <- cycle UP first 1_dw x 1_up
  function hubbard_site(hloc) result(self)
    real(8),dimension(Nspin,Norb,Norb),optional :: hloc  
    real(8),dimension(Nspin,Norb,Norb)          :: hloc_
    type(site)                                  :: self
    real(8)                                     :: uloc
    real(8)                                     :: xmu
    real(8),dimension(:,:),allocatable          :: H
    real(8),dimension(:,:),allocatable          :: Cup,Cdw,Cp
    !
    print*,Nspin,Norb,size(hloc_)
    hloc_ = zeros(Nspin,Norb,Norb);if(present(hloc))hloc_=hloc
    !
    call Init_LocalFock_Sectors()
    !
    H   = build_Hlocal_operator(hloc_)
    Cp  = build_C_operator(1,1);Cup=KId(1).kx.Cp
    Cp  = build_C_operator(1,2);Cdw=Cp.kx.KSz(1)
    !
    call self%free()
    self%Dim = 4
    call self%put("H",sparse(H))
    call self%put("Cup",sparse(Cup))
    call self%put("Cdw",sparse(Cdw))
    call self%put("P",sparse(kSz(2)))
    allocate(self%sectors(1))
    self%sectors(1) = sectors_list( tbasis([0,0, 1,0, 0,1, 1,1],Qdim=2) )    
  end function hubbard_site



END MODULE SITES








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
program testSITES
  USE SCIFOR,  id => zeye
  USE INPUT_VARS
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  USE TUPLE_BASIS
  USE LIST_SECTORS
  USE SITES
  implicit none


  type(site)                       :: my_site,a,b
  type(tbasis)                     :: sz_basis
  real(8),dimension(2,2),parameter :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
  real(8),dimension(2,2),parameter :: Sz=dble(pauli_z)
  real(8),dimension(2,2),parameter :: Sx=dble(pauli_x)
  real(8),dimension(2,2),parameter :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
  real(8),dimension(4,4)           :: Gamma13,Gamma03

  Gamma13=kron(Sx,Sz)
  Gamma03=kron(eye(2),Sz)

  sz_basis = tbasis([0.5d0,-0.5d0],Qdim=1)

  my_site = site(&
       dim      = 2, &
       sectors  = [sectors_list(sz_basis)],&
       operators= operators_list(['H0','Sz','Sp'],&
       [sparse(Hzero),sparse(Sz),sparse(Splus)]))
  print*,"Is site valid:",my_site%is_valid()
  call my_site%show()



  a = my_site

  print*,"Test =: a=my_site"
  call a%show()
  print*,a%is_valid()
  print*,"modify a, check a and my_site"
  call a%put("G5",sparse(Gamma03))
  print*,"Is A site valid:",a%is_valid()
  print*,"Is My_site site valid:",my_site%is_valid()


  print*,"Test PAULI SITE"
  b = pauli_site()
  call b%show()
  call b%free()


  call read_input("DMRG.conf")
  b = hubbard_site()
  call b%show()
  print*,"Is site valid:",b%is_valid()
  call b%free()
end program testSITES
#endif
