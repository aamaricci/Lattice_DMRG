program testEDkron
  USE SCIFOR
  USE AUX_FUNCS
  USE MATRIX_SPARSE, id=>sp_eye
  USE MATRIX_BLOCKS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
  implicit none

  character(len=64)                  :: finput
  real(8)                            :: ts,uloc
  integer                            :: Lmax,i,m,Neigen,current_L
  integer                            :: lanc_ncv_factor
  integer                            :: lanc_ncv_add
  integer                            :: lanc_niter
  real(8)                            :: lanc_tolerance
  integer                            :: lanc_threshold
  type(block)                        :: my_block,dimer,trimer,sys,env
  type(sparse_matrix)                :: spHsb,spH
  type(site)                         :: dot
  real(8)                            :: gs_energy,target_Sz
  integer                            :: unit
  integer                            :: model_d=4

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"ts",finput,default=-1d0,comment="Hopping amplitude")
  call parse_input_variable(uloc,"uloc",finput,default=0d0,comment="Hubbard U")
  call parse_input_variable(Lmax,"LMAX",finput,default=1,comment="Final chain length ")
  call parse_input_variable(m,"M",finput,default=20,&
       comment="Number of states retained at truncation of \rho")
  call parse_input_variable(Neigen,"NEIGEN",finput,default=1,&
       comment="Number of eigenpairs required to SB Hamiltonian")
  call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",finput,default=10,&
       comment="Arpack Block size parameters")
  call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",finput,default=0,&
       comment="Arpack Block size parameters")
  call parse_input_variable(lanc_niter,"LANC_NITER",finput,default=512,&
       comment="Number of Lanczos iteration in spectrum determination.")
  call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",finput,default=1d-12,&
       comment="Lanczos tolerance ")
  call parse_input_variable(lanc_threshold,"LANC_THRESHOLD",finput,default=2,&
       comment="Lapack threshold for Arpack ")
  target_Sz=0d0




  !
  !Init the single dot structure:
  dot = hubbard_site(xmu=0d0,Uloc=uloc)
  !
  !Init block from single dot
  sys=block(dot)
  env=sys


  !Run DMRG algorithm
  open(free_unit(unit),file="energyVSlength_L"//str(Lmax)//"_m"//str(m)//".dmrg")
  do i=1,Lmax
     current_L = 2*sys%length + 2
     write(*,*)"SuperBlock Length  =",current_L
     call start_timer()
     call single_dmrg_step(m,sys,env,energy=gs_energy)
     call stop_timer("dmrg_step")
     write(*,*)"E/L, E=",gs_energy/current_L,gs_energy
     write(unit,*)current_L,gs_energy/current_L,gs_energy
     write(*,*)"-------"
     write(*,*)""
  enddo
  close(unit)




contains






  function enlarge_block(self,dot,grow) result(enl_self)
    type(block),intent(inout)        :: self
    type(site),intent(inout)         :: dot
    character(len=*),optional        :: grow
    character(len=16)                :: grow_
    type(block)                      :: enl_self
    integer                          :: mblock,len
    real(8),dimension(:),allocatable :: self_basis, dot_basis
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    mblock =  self%dim
    len    =  self%length
    !
    enl_self%length = len + 1
    enl_self%dim    = mblock*model_d
    !
    select case(str(grow_))
    case ("left","l")
       call enl_self%put("H", (self%operators%op("H").x.id(model_d)) + (id(mblock).x.dot%operators%op("H")) + &
            H2model(self,as_block(dot)))
       call enl_self%put("Cup", Id(mblock).x.dot%operators%op("Cup"))
       call enl_self%put("Cdw", Id(mblock).x.dot%operators%op("Cdw"))
       call enl_self%put("P"  , Id(mblock).x.dot%operators%op("P"))
    case ("right","r")
       call enl_self%put("H", (id(model_d).x.self%operators%op("H")) +  (dot%operators%op("H").x.id(mblock)) + &
            H2model(as_block(dot),self))
       call enl_self%put("Cup", dot%operators%op("Cup").x.Id(mblock))
       call enl_self%put("Cdw", dot%operators%op("Cdw").x.Id(mblock))
       call enl_self%put("P"  , dot%operators%op("P").x.Id(mblock))
    end select
    !
    ! self_basis = self%sectors%basis()
    ! dot_basis  = dot%sectors%basis()
    ! call enl_self%set_sectors( outsum(self_basis,dot_basis) )
    ! !
    ! deallocate(self_basis,dot_basis)
  end function enlarge_block



  !H_lr = -t \sum_sigma[ (C^+_{l,sigma}@P_l) x C_{r,sigma}]  + H.c.
  function H2model(left,right) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    type(sparse_matrix)           :: CupL,CdwL
    type(sparse_matrix)           :: CupR,CdwR
    type(sparse_matrix)           :: H2,P
    P    = left%operators%op("P") 
    CupL = left%operators%op("Cup")
    CdwL = left%operators%op("Cdw")
    CupR = right%operators%op("Cup")
    CdwR = right%operators%op("Cdw")
    H2   = ts*(matmul(CupL%t(),P).x.CupR)  + ts*(matmul(CdwL%t(),P).x.CdwR)
    H2   = H2 + H2%dgr()
    call P%free
    call CupL%free
    call CdwL%free
    call CupR%free
    call CdwR%free
  end function H2model



  subroutine single_dmrg_step(m,sys,env,energy)
    integer                            :: m
    type(block)                        :: sys,env
    type(block)                        :: esys,eenv
    integer                            :: m_esys,m_eenv,m_sb,m_,i,j,Ncv,isys,Nsys,Nenv,im
    real(8),dimension(:),allocatable   :: eig_values
    real(8),dimension(:,:),allocatable :: eig_basis
    real(8),dimension(:),allocatable   :: evals,evalsL,evalsR
    real(8),dimension(:,:),allocatable :: Hsb,rhoL,rhoR,truncation_rhoL,truncation_rhoR
    real(8)                            :: truncation_errorL,truncation_errorR,energy
    !
    !Check if blocks are valid ones
    if(.not.sys%is_valid())stop "single_dmrg_step error: sys is not a valid block"
    !
    !Enlarge blocks
    print*,"Enlarge blocks:"
    esys = enlarge_block(sys,dot,grow='left')
    eenv = enlarge_block(env,dot,grow='right')
    if(.not.esys%is_valid())stop "single_dmrg_step error: enlarged_sys is not a valid block"
    if(.not.eenv%is_valid())stop "single_dmrg_step error: enlarged_env is not a valid block"
    print*,"Done"
    !
    !Get Enlarged Sys/Env actual dimension
    print*,"Enlarged Blocks dimensions:",esys%dim,eenv%dim    
    m_esys = esys%dim
    m_eenv = eenv%dim
    !
    !Build SuperBLock Sector
    call start_timer()
    m_sb = m_esys*m_eenv
    print*,"Super Block dimensions:", m_sb
    !BUild SuperBlock matrix:
    spHsb =  (esys%operators%op("H").x.id(m_eenv)) + (id(m_esys).x.eenv%operators%op("H")) + H2model(esys,eenv)
    call stop_timer("Build H_sb")
    !
    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    allocate(eig_values(Neigen))    ;eig_values=0d0
    allocate(eig_basis(m_sb,Neigen));eig_basis =0d0
    call start_timer()
    if(m_sb < lanc_threshold)then
       print*,"diag SB w/ Lapack"
       allocate(Hsb(m_sb,m_sb));Hsb=0d0
       allocate(evals(m_sb))
       call spHsb%dump(Hsb)
       call eigh(Hsb,evals)
       eig_basis(:,1:Neigen) = Hsb(:,1:Neigen)
       eig_values(1:Neigen)  = evals(1:Neigen)
       deallocate(Hsb,evals)
    else
       print*,"diag SB w/ Arpack"
       Ncv = min(lanc_ncv_factor*Neigen+lanc_ncv_add,m_sb)
       call sp_eigh(sb_HxV,eig_values,eig_basis,&
            Ncv,&
            lanc_niter,&
            tol=lanc_tolerance,&
            iverbose=.false.)
    end if
    energy = eig_values(1)
    call stop_timer("Diag H_sb")
    !
    call start_timer()
    allocate(rhoL(m_esys,m_esys));rhoL=0d0
    allocate(rhoR(m_eenv,m_eenv));rhoR=0d0
    rhoL = build_density_matrix(m_esys,m_eenv,eig_basis(:,1),'left')
    rhoR = build_density_matrix(m_esys,m_eenv,eig_basis(:,1),'right')
    !
    allocate(evalsL(m_esys))
    allocate(evalsR(m_eenv))
    call eigh(rhoL,evalsL)
    call eigh(rhoR,evalsR)
    evalsL    = evalsL(m_esys:1:-1)
    evalsR    = evalsR(m_esys:1:-1)
    rhoL(:,:) = rhoL(:,m_esys:1:-1)
    rhoR(:,:) = rhoR(:,m_esys:1:-1)
    call stop_timer("Get Rho_Left-dot / Rho_dot-Right")
    !
    !
    call start_timer()
    m_ = min(m,m_esys);
    allocate(truncation_rhoL(m_esys,m_))
    allocate(truncation_rhoR(m_eenv,m_))
    truncation_rhoL(:,1:m_) = rhoL(:,1:m_)
    truncation_rhoR(:,1:m_) = rhoR(:,1:m_)
    truncation_errorL = 1d0 - sum(evalsL(1:m_))
    truncation_errorR = 1d0 - sum(evalsR(1:m_))
    call stop_timer("Build U")
    write(*,*)"Truncating      :",m<=m_esys
    write(*,*)"Truncation dim  :",m_
    write(*,*)"Truncation error:",truncation_errorL,truncation_errorR
    !
    !find a clever way to iterate here:
    !Return updated Block:
    call start_timer()
    sys%length = esys%length
    sys%dim    = m_
    call sys%put("H",rotate_and_truncate(esys%operators%op("H"),truncation_rhoL,m_esys,m_))
    call sys%put("Cup",rotate_and_truncate(esys%operators%op("Cup"),truncation_rhoL,m_esys,m_))
    call sys%put("Cdw",rotate_and_truncate(esys%operators%op("Cdw"),truncation_rhoL,m_esys,m_))
    call sys%put("P",rotate_and_truncate(esys%operators%op("P"),truncation_rhoL,m_esys,m_))


    env%length = eenv%length
    env%dim    = m_
    call env%put("H",rotate_and_truncate(eenv%operators%op("H"),truncation_rhoR,m_eenv,m_))
    call env%put("Cup",rotate_and_truncate(eenv%operators%op("Cup"),truncation_rhoR,m_eenv,m_))
    call env%put("Cdw",rotate_and_truncate(eenv%operators%op("Cdw"),truncation_rhoR,m_eenv,m_))
    call env%put("P",rotate_and_truncate(eenv%operators%op("P"),truncation_rhoR,m_eenv,m_))

    call stop_timer("Update+Truncate Block")
    !
    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    if(allocated(evalsL))deallocate(evalsL)
    if(allocated(evalsR))deallocate(evalsR)
    if(allocated(rhoL))deallocate(rhoL)
    if(allocated(rhoR))deallocate(rhoR)
    if(allocated(truncation_rhoL))deallocate(truncation_rhoL)
    if(allocated(truncation_rhoR))deallocate(truncation_rhoR)
    call esys%free()
    call eenv%free()
    call spHsb%free()
  end subroutine single_dmrg_step



  function build_density_matrix(nsys,nenv,psi,direction) result(rho)
    integer                      :: nsys
    integer                      :: nenv
    real(8),dimension(nsys*nenv) :: psi
    character(len=*)             :: direction
    real(8),dimension(nsys,nenv) :: rho_restricted
    real(8),dimension(nsys,nsys) :: rho
    rho_restricted = transpose(reshape(psi, [nenv,nsys]))
    select case(to_lower(str(direction)))
    case ('left','l')
       rho  = matmul(rho_restricted,  transpose(rho_restricted) )
    case ('right','r')
       rho  = matmul(transpose(rho_restricted), rho_restricted  )
    end select
  end function build_density_matrix



  function rotate_and_truncate(Op,trRho,N,M) result(RotOp)
    type(sparse_matrix),intent(in) :: Op
    real(8),dimension(N,M)         :: trRho  ![Nesys,M]
    integer                        :: N,M
    type(sparse_matrix)            :: RotOp
    real(8),dimension(M,M)         :: Umat
    real(8),dimension(N,N)         :: OpMat
    N = size(trRho,1)
    M = size(trRho,2)
    if( any( shape(Op) /= [N,N] ) ) stop "rotate_and_truncate error: shape(Op) != [N,N] N=size(Rho,1)"
    OpMat= Op%as_matrix()
    Umat = matmul( transpose(trRho), matmul(OpMat,trRho)) ![M,N].[N,N].[N,M]=[M,M]
    call RotOp%load( Umat )
    !
    ! RotOp = Op
  end function rotate_and_truncate


  subroutine sb_HxV(Nloc,v,Hv)
    integer                    :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                    :: i,j,jcol
    Hv=0d0
    do i=1,Nloc
       matmul: do jcol=1, spHsb%row(i)%Size
          val = spHsb%row(i)%vals(jcol)
          j   = spHsb%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end subroutine sb_HxV


  function to_upper(StrIn) result(StrOut)
    character(len=*), intent(in) :: strIn
    character(len=len(strIn))    :: strOut
    integer :: i
    do i = 1,len(StrIn)
       select case(StrIn(i:i))
       case("a":"z")
          StrOut(i:i) = achar(iachar(StrIn(i:i))-32)
       case default
          StrOut(i:i) = StrIn(i:i)
       end select
    end do
  end function to_upper

  function to_lower(StrIn) result(StrOut)
    character(len=*), intent(in) :: strIn
    character(len=len(strIn))    :: strOut
    integer :: i
    do i = 1,len(StrIn)
       select case(StrIn(i:i))
       case("A":"Z")
          StrOut(i:i) = achar(iachar(StrIn(i:i))+32)
       case default
          StrOut(i:i) = StrIn(i:i)
       end select
    end do
  end function to_lower

  subroutine print_mat(M,unit)
    real(8),dimension(:,:) :: M
    integer :: i,j,stride,unit,n
    stride=1
    write(unit,*)"Matrix: "
    do i=1,size(M,1)
       do j=1,size(M,2)
          write(unit,"(("//str(stride)//"I2))",advance="no")int(M(i,j))
          if(mod(j,stride)==0)write(unit,"(A1)",advance="no")""
       enddo
       write(unit,*)
       if(mod(i,stride)==0)write(unit,*)
    enddo
    write(unit,*)""
  end subroutine print_mat

end program testEDkron
