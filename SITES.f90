MODULE SITES
  USE SCIFOR, only: str,diag,zeros,eye,kron,pauli_x,pauli_y,xi,operator(.kx.)
  USE AUX_FUNCS
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE HLOCAL
  implicit none
  private


  type site
     integer                                     :: Dim=1
     integer                                     :: DimUp=1
     integer                                     :: DimDw=1
     type(sectors_list),dimension(:),allocatable :: sectors
     type(operators_list)                        :: operators
   contains
     procedure,pass     :: free        => free_site
     procedure,pass     :: put         => put_op_site
     procedure,pass     :: load        => load_op_site
     procedure,pass     :: set_sectors => set_sectors_site
     procedure,pass     :: show        => show_site
     procedure,pass     :: is_valid    => is_valid_site
  end type site


  !GENERIC CONSTRUCTOR
  interface site
     module procedure :: constructor_site
  end interface site

  !EQUALITY 
  interface assignment(=)
     module procedure :: site_equal_site
  end interface assignment(=)


  public :: site
  public :: pauli_site
  public :: spin_onehalf_site
  public :: spin_one_site
  public :: hubbard_site
  public :: hubbard_site_ud
  public :: assignment(=)

  integer :: i

contains


  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function constructor_site(Dims,sectors,operators) result(self)
    integer,intent(in)              :: Dims(2)
    type(sectors_list),intent(in)   :: sectors(:)
    type(operators_list),intent(in) :: operators
    type(site)                      :: self
    self%DimUp     = Dims(1)
    self%DimDw     = Dims(2)
    self%Dim       = product(Dims)
    self%operators = operators
    allocate(self%sectors(size(sectors)))
    do i=1,size(self%sectors)
       self%sectors(i)   = sectors(i)
    enddo
  end function constructor_site



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
  !PURPOSE:  Put a QN array in the site
  !+------------------------------------------------------------------+
  subroutine set_sectors_site(self,indx,vec)
    class(site)          :: self
    integer              :: indx
    real(8),dimension(:) :: vec
    if(indx<1.OR.indx>size(self%sectors))stop "SET_SECTORS_SITE ERROR: indx out of range"
    self%sectors(indx) = sectors_list(vec)
  end subroutine set_sectors_site



  !##################################################################
  !##################################################################
  !              OPERATIONS / ASSIGNEMENTS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Equality between two sites
  !+------------------------------------------------------------------+
  subroutine site_equal_site(A,B)
    type(site),intent(inout) :: A
    type(site),intent(in)    :: B
    call A%free
    A%DimUp  = B%DimUp
    A%DimDw  = B%DimDw
    A%Dim    = B%Dim
    A%operators = B%operators
    allocate(A%sectors(size(B%sectors)))
    do i=1,size(A%sectors)
       A%sectors(i)   = B%sectors(i)
    enddo
  end subroutine site_equal_site


  function is_valid_site(self) result(bool)
    class(site)                           :: self
    logical                               :: bool
    integer,dimension(size(self%sectors)) :: Lvec
    bool = self%operators%is_valid([self%DimUp,self%DimDw])
    do i=1,size(self%sectors)
       Lvec = len(self%sectors(i))
    enddo
    bool=bool.AND.(self%dim==product(Lvec))
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
    write(*,*)"Site DimUp    =",self%DimUp
    write(*,*)"Site DimDw    =",self%DimDw
    write(*,*)"Site Dim      =",self%Dim
    do i=1,size(self%sectors)
       write(*,*)"Site Sectors: ",i
       call self%sectors(i)%show()
    enddo
    write(*,*)"Site Operators:"
    call self%operators%show(fmt=fmt_)
  end subroutine show_site





  !##################################################################
  !##################################################################
  !              PRESET
  !##################################################################
  !##################################################################
  function pauli_site() result(self)
    type(site)                       :: self
    real(8),dimension(2,2),parameter :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
    real(8),dimension(2,2),parameter :: Szeta=reshape([1d0,0d0,0d0,-1d0],[2,2])
    real(8),dimension(2,2),parameter :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
    call self%free()
    self%DimUp=2
    self%DimDw=1
    self%Dim=2
    call self%put("H",sparse(Hzero))
    call self%put("Sz",sparse(0.5d0*Szeta))
    call self%put("Sp",sparse(Splus))
    allocate(self%sectors(1))
    self%sectors(1) = sectors_list([0.5d0,-0.5d0])
  end function pauli_site


  function spin_onehalf_site() result(self)
    type(site)                       :: self
    real(8),dimension(2,2),parameter :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
    real(8),dimension(2,2),parameter :: Szeta=reshape([1d0,0d0,0d0,-1d0],[2,2])
    real(8),dimension(2,2),parameter :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
    call self%free()
    self%DimUp=2
    self%DimDw=1
    self%Dim=2
    call self%put("H",sparse(Hzero))
    call self%put("Sz",sparse(0.5d0*Szeta))
    call self%put("Sp",sparse(Splus))
    allocate(self%sectors(1))
    self%sectors(1) = sectors_list([0.5d0,-0.5d0])
  end function spin_onehalf_site


  function spin_one_site() result(self)
    type(site) :: self
    real(8),dimension(3,3) :: Hzero
    real(8),dimension(3,3) :: Szeta
    real(8),dimension(3,3) :: Splus
    !
    Hzero = 0d0
    Szeta = diag([-1d0,0d0,1d0])
    Splus = 0d0
    Splus(1,2) = sqrt(2d0)
    Splus(2,3) = sqrt(2d0)
    !
    call self%free()
    self%DimUp=3
    self%DimDw=1
    self%Dim=3
    call self%put("H",sparse(Hzero))
    call self%put("Sz",sparse(Szeta))
    call self%put("Sp",sparse(Splus))
    allocate(self%sectors(1))
    self%sectors(1) = sectors_list([-1d0,0d0,1d0])
  end function spin_one_site


  !Fock = [|0,0>,|0,1>,|1,0>,|1,1>] <- |up,dw> <- fill DW first = 1_dw x 1_up
  function hubbard_site(uloc,xmu) result(self)
    type(site)                         :: self
    real(8)                            :: uloc
    real(8)                            :: xmu
    real(8),dimension(:,:),allocatable :: Hloc
    real(8),dimension(:,:),allocatable :: Cup,Cdw,Cp
    !
    call Init_LocalFock_Sectors(1,1)
    !
    Hloc = build_Hlocal_operator(hloc=dreal(zeros(1,1,1)),xmu=xmu,uloc=uloc)
    Cp  = build_C_operator(1,1);Cup=KId(1).kx.Cp
    Cp  = build_C_operator(1,2);Cdw=Cp.kx.KSz(1)
    !
    call self%free()
    self%DimUp=2
    self%DimDw=2
    self%Dim=4
    call self%put("H",sparse(Hloc))
    call self%put("Cup",sparse(Cup))
    call self%put("Cdw",sparse(Cdw))
    call self%put("P",sparse(kSz(2)))
    allocate(self%sectors(2))
    self%sectors(1) = sectors_list([0d0,1d0])
    self%sectors(2) = sectors_list([0d0,1d0])
  end function hubbard_site



  !Fock =  [|0>,|1>]_dw x [|0>,|1>]_up reproduces the same as above
  function hubbard_site_ud(uloc,xmu) result(self)
    type(site)                         :: self
    real(8)                            :: uloc
    real(8)                            :: xmu
    real(8),dimension(:,:),allocatable :: Hd
    real(8),dimension(:,:),allocatable :: Cup,Cdw,Cp
    !
    call Init_LocalFock_Sectors(1,1)
    !
    Hd = build_Hlocal_operator(hloc=dreal(zeros(1,1,1)),xmu=xmu,uloc=uloc)
    Cp = build_C_operator(1,1)
    !== identical for Up and Dw as it acts in reduced Fock space of a given spin
    !
    call self%free()
    self%Dim=4
    self%DimUp=2
    self%DimDw=2
    call self%put("Hd",sparse(Hd)) !this is diagonal at the beginning
    call self%put("Hup",sparse(dble(zeros(2,2))))
    call self%put("Hdw",sparse(dble(zeros(2,2)))) 
    call self%put("Cp",sparse(Cp))   !2x2 
    call self%put("P",sparse(kSz(1))) !2x2
    allocate(self%sectors(2))
    self%sectors(1) = sectors_list([0d0,1d0])
    self%sectors(2) = sectors_list([0d0,1d0])
  end function hubbard_site_ud


END MODULE SITES
