MODULE SITES
  USE SCIFOR, only: str,diag,zeros,eye,kron
  USE AUX_FUNCS
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE HLOCAL
  implicit none
  private


  type site
     integer              :: dim=1
     type(sectors_list)   :: sectors
     type(operators_list) :: operators
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
  public :: assignment(=)


contains


  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function constructor_site(dim,sectors,operators) result(self)
    integer,intent(in)              :: dim
    type(sectors_list),intent(in)   :: sectors
    type(operators_list),intent(in) :: operators
    type(site)                      :: self
    self%dim       = dim
    self%sectors   = sectors
    self%operators = operators
  end function constructor_site



  !+------------------------------------------------------------------+
  !PURPOSE:  Destructor
  !+------------------------------------------------------------------+
  subroutine free_site(self)
    class(site)               :: self
    self%dim   = 1
    call self%operators%free()
    call self%sectors%free()
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
  subroutine set_sectors_site(self,vec)
    class(site)                 :: self
    real(8),dimension(self%dim) :: vec
    self%sectors = sectors_list(vec)
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
    A%dim       = B%dim
    A%sectors   = B%sectors
    A%operators = B%operators
  end subroutine site_equal_site


  function is_valid_site(self) result(bool)
    class(site) :: self
    logical     :: bool
    bool = self%operators%is_valid(self%dim).AND.(self%dim==len(self%sectors))
  end function is_valid_site





  !##################################################################
  !##################################################################
  !              SHOW 
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Pretty print the site structure
  !+------------------------------------------------------------------+
  subroutine show_site(self,dble,fmt)
    class(site)               :: self
    logical,optional          :: dble
    character(len=*),optional :: fmt
    logical                   :: dble_
    character(len=32)         :: fmt_
    dble_=show_dble;if(present(dble))dble_=dble
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    write(*,*)"Site Dim      =",self%dim
    write(*,*)"Site Sectors  ="
    call self%sectors%show()
    write(*,*)"Site Operators="
    call self%operators%show(dble=dble_,fmt=fmt_)
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
    self%dim=2
    call self%put("H",sparse(Hzero))
    call self%put("Sz",sparse(0.5d0*Szeta))
    call self%put("Sp",sparse(Splus))
    self%sectors = sectors_list([0.5d0,-0.5d0])
  end function pauli_site


  function spin_onehalf_site() result(self)
    type(site)                       :: self
    real(8),dimension(2,2),parameter :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
    real(8),dimension(2,2),parameter :: Szeta=reshape([1d0,0d0,0d0,-1d0],[2,2])
    real(8),dimension(2,2),parameter :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
    call self%free()
    self%dim=2
    call self%put("H",sparse(Hzero))
    call self%put("Sz",sparse(0.5d0*Szeta))
    call self%put("Sp",sparse(Splus))
    self%sectors = sectors_list([0.5d0,-0.5d0])
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
    self%dim=3
    call self%put("H",sparse(Hzero))
    call self%put("Sz",sparse(Szeta))
    call self%put("Sp",sparse(Splus))
    self%sectors = sectors_list([-1d0,0d0,1d0])
  end function spin_one_site


  !up : |0>, |1>
  !dw : |0>, |1>
  !Fock = up x dw
  function hubbard_site(uloc,xmu) result(self)
    type(site)                         :: self
    real(8)                            :: uloc
    real(8)                            :: xmu
    real(8),dimension(:,:),allocatable :: Hloc
    real(8),dimension(:,:),allocatable :: Cup,Cdw
    !
    call Init_LocalFock_Sectors(1,1)
    !
    Hloc = build_Hlocal_operator(hloc=dreal(zeros(1,1,1)),xmu=xmu,uloc=uloc)
    Cup  = build_C_operator(ispin=1,iorb=1)
    Cdw  = build_C_operator(ispin=2,iorb=1)
    !
    call self%free()
    self%dim=4
    call self%put("H",sparse(Hloc))
    !One should avoid this (which makes useless the simplifications in HLOCAL)
    !However this requires operator_list to accept different dimensions (2 rather than 4 here).
    call self%put("Cup",sparse(kron(eye(2),Cup)))
    call self%put("Cdw",sparse(kron(Cdw,eye(2))))
    self%sectors = sectors_list([0d0,1d0])
  end function hubbard_site



END MODULE SITES
