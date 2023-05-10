program testDMRGinfty
  USE SCIFOR
  USE AUX_FUNCS
  USE MATRIX_SPARSE, id=>sp_eye
  USE MATRIX_BLOCKS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
  implicit none

  character(len=64)   :: finput
  real(8)             :: Jx,Jz
  integer             :: Lmax,i,m,Neigen,current_L
  integer             :: lanc_ncv_factor
  integer             :: lanc_ncv_add
  integer             :: lanc_niter
  real(8)             :: lanc_tolerance
  integer             :: lanc_threshold
  type(block)         :: my_block
  type(sparse_matrix) :: Szeta,Splus,Hzero,Sminus
  type(sparse_matrix) :: sb_hamiltonian,spHsb
  type(site)          :: dot
  real(8)             :: gs_energy,esys_qn,eenv_qn,target_Sz,sb_qn,val
  integer,allocatable :: esys_map(:),eenv_map(:),sb_states(:)
  integer             :: unit,j,ioffset,current_index
  integer,parameter   :: model_d=2

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(Jx,"Jx",finput,default=1d0,comment="Jx term in the XXZ model")
  call parse_input_variable(Jz,"Jz",finput,default=1d0,comment="Jz term in the XXZ model")
  call parse_input_variable(Lmax,"LMAX",finput,default=20,comment="Final chain length ")
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
  call parse_input_variable(lanc_threshold,"LANC_THRESHOLD",finput,default=256,&
       comment="Lapack threshold for Arpack ")


  !
  !Init the single dot structure:
  dot = spin_onehalf_site()
  !
  !Init block from single dot
  my_block=block(dot)
  call my_block%show()


  !RUN DMRG:
  open(free_unit(unit),file="energyVSlength_L"//str(Lmax)//"_m"//str(m)//".dmrg")
  do i=1,Lmax
     current_L = 2*my_block%length + 2
     write(*,*)"SuperBlock Length  =",current_L
     call start_timer()
     call single_dmrg_step(M,my_block,energy=gs_energy)
     call stop_timer("dmrg_step")
     write(*,*)"E/L=",gs_energy/current_L
     write(unit,*)current_L,gs_energy/current_L,gs_energy
     write(*,*)"-------"
     write(*,*)""
  enddo
  close(unit)


contains

  function H2model(left,right) result(H2)
    type(block)         :: left
    type(block)         :: right
    type(sparse_matrix) :: Sz1,Sp1
    type(sparse_matrix) :: Sz2,Sp2
    type(sparse_matrix) :: H2
    Sz1 = left%operators%op("Sz")
    Sp1 = left%operators%op("Sp")
    Sz2 = right%operators%op("Sz")
    Sp2 = right%operators%op("Sp")
    H2 = Jx/2d0*(Sp1.x.Sp2%dgr()) +  Jx/2d0*(Sp1%dgr().x.Sp2)  + Jz*(Sz1.x.Sz2)
    call Sz1%free()
    call Sp1%free()
    call Sz2%free()
    call Sp2%free()
  end function H2model



  function enlarge_block(self,dot) result(enl_self)
    class(block),intent(inout) :: self
    type(site),intent(inout)   :: dot
    type(block)                :: enl_self
    integer                    :: mblock,len
    !
    mblock =  self%dim
    len    =  self%length
    !
    call enl_self%put("H", &
         (self%operators%op("H").x.id(model_d)) + &
         (id(mblock).x.dot%operators%op("H")) + &
         H2model(self, as_block(dot)) )
    call enl_self%put("Sz", id(mblock).x.dot%operators%op("Sz"))
    call enl_self%put("Sp", id(mblock).x.dot%operators%op("Sp"))

    enl_self%length = self%length + 1
    enl_self%dim    = mblock*model_d
    call enl_self%set_sectors( dble(arange(1,enl_self%dim)) )
  end function enlarge_block




  subroutine single_dmrg_step(m,sys,energy)
    integer                            :: m
    class(block)                       :: sys
    !
    type(block)                        :: sys_enl
    integer                            :: m_esys,m_sb,m_,i,Ncv
    real(8),dimension(:),allocatable   :: eig_values,evals
    real(8),dimension(:,:),allocatable :: eig_basis,evecs
    real(8),dimension(:,:),allocatable :: rho,truncation_rho
    real(8)                            :: truncation_error,energy

    !Check if blocks are valid ones
    if(.not.sys%is_valid())stop "single_dmrg_step error: sys is not a valid block"
    !
    !Enlarge blocks
    sys_enl = enlarge_block(sys,dot)
    if(.not.sys_enl%is_valid())stop "single_dmrg_step error: enlarged_sys is not a valid block"
    !
    !Get Enlarged Sys/Env actual dimension
    m_esys = sys_enl%dim


    !Build SuperBlock Hamiltonian
    m_sb = m_esys*m_esys
    print*,"Super Block dimension:",m_sb
    call start_timer()
    sb_hamiltonian = (&
         sys_enl%operators%op("H").x.id(m_esys))  + &
         (id(m_esys).x.sys_enl%operators%op("H")) + &
         H2model(sys_enl,sys_enl)
    call stop_timer("Build H_sb")


    allocate(eig_values(Neigen))    ;eig_values=0d0
    allocate(eig_basis(m_sb,Neigen));eig_basis =0d0
    call start_timer()
    if(m_sb < lanc_threshold)then
       print*,"diag SB w/ Lapack"
       allocate(rho(m_sb,m_sb));rho=zero
       allocate(evals(m_sb))
       call sb_hamiltonian%dump(rho)
       call eigh(rho,evals)
       eig_basis(:,1:Neigen) = rho(:,1:Neigen)
       eig_values(1:Neigen)  = evals(1:Neigen)
       print*,evals(1:min(5,size(evals)))
       deallocate(rho,evals)
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


    call start_timer()
    allocate(rho(m_esys,m_esys));rho=0d0
    allocate(evals(m_esys))
    rho = build_density_matrix(m_esys,m_esys,eig_basis(:,1))
    call eigh(rho,evals)
    evals    = evals(m_esys:1:-1)
    rho(:,:) = rho(:,m_esys:1:-1)
    call stop_timer("Get Rho")


    call start_timer()
    m_ = min(m,m_esys);
    allocate(truncation_rho(m_esys,m_))
    truncation_rho(:,1:m_) = rho(:,1:m_)
    truncation_error = 1d0 - sum(evals(1:m_))
    call stop_timer("Build U")
    write(*,*)"Truncating      :",m<=m_esys
    write(*,*)"Truncation dim  :",m_
    write(*,*)"Truncation error:",truncation_error


    !Return updated Block:
    call start_timer()
    sys%length = sys_enl%length
    sys%dim    = m_
    !find a clever way to iterate here:
    call sys%put("H",&
         rotate_and_truncate(sys_enl%operators%op("H"),truncation_rho,m_esys,m_))
    call sys%put("Sz",&
         rotate_and_truncate(sys_enl%operators%op("Sz"),truncation_rho,m_esys,m_))
    call sys%put("Sp",&
         rotate_and_truncate(sys_enl%operators%op("Sp"),truncation_rho,m_esys,m_))
    call sys%set_sectors( dble(arange(1,sys%dim)) )
    call stop_timer("Update+Truncate Block")



    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    if(allocated(evals))deallocate(evals)
    if(allocated(rho))deallocate(rho)
    call sb_hamiltonian%free()
    call sys_enl%free()
  end subroutine single_dmrg_step



  function build_density_matrix(nsys,nenv,psi) result(rho)
    integer                      :: nsys
    integer                      :: nenv
    real(8),dimension(nsys,nenv) :: psi
    real(8),dimension(nsys,nsys) :: rho
    rho  = matmul(psi,  (transpose(psi)) )
  end function build_density_matrix


  function rotate_and_truncate(Op,trRho,N,M) result(RotOp)
    type(sparse_matrix),intent(in)        :: Op
    real(8),dimension(N,M)             :: trRho  ![Nesys,M]
    integer                               :: N,M
    type(sparse_matrix)                   :: RotOp
    real(8),dimension(M,M)             :: Umat
    N = size(trRho,1)
    M = size(trRho,2)
    if( any( shape(Op) /= [N,N] ) ) stop "rotate_and_truncate error: shape(Op) != [N,N] N=size(Rho,1)"
    Umat = matmul( (transpose(trRho)), matmul(Op%as_matrix(),trRho)) ![M,N].[N,N].[N,M]=[M,M]
    call RotOp%load( Umat )
  end function rotate_and_truncate


  subroutine sb_HxV(Nloc,v,Hv)
    integer                    :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                    :: i,j,jcol
    Hv=zero
    do i=1,Nloc
       matmul: do jcol=1, sb_hamiltonian%row(i)%Size
          val = sb_hamiltonian%row(i)%vals(jcol)
          j   = sb_hamiltonian%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end subroutine sb_HxV



  subroutine dmrg_graphic(sys,label)
    class(block)     :: sys
    character(len=1) :: label
    integer :: i
    select case(label)
    case default; stop "dmrg_graphic error: label != l,r"
    case("l","L")
       write(*,"("//str(sys%length)//"A1)",advance="no")("=",i=1,sys%length)
       write(*,"(A2)",advance="no")"**"
       write(*,"("//str(sys%length)//"A1)",advance="yes")("-",i=1,sys%length)
    case("r","R")
       write(*,"("//str(sys%length)//"A1)",advance="no")("=",i=1,sys%length)
       write(*,"(A2)",advance="no")"**"
       write(*,"("//str(sys%length)//"A1)",advance="yes")("-",i=1,sys%length)
    end select
  end subroutine dmrg_graphic





  subroutine print_mat(M)
    real(8),dimension(:,:) :: M
    integer :: i,j
    do i=1,size(M,1)
       write(*,"("//str(size(M,2))//"(ES12.3,1x))")(M(i,j),j=1,size(M,2))
    enddo
  end subroutine print_mat


end program testDMRGinfty






! function truncate(Op,Umat,N,M) result(Rmat)
!   complex(8),dimension(N,N)      :: Op
!   integer                        :: N,M
!   complex(8),dimension(N,M)      :: Umat
!   complex(8),dimension(M,M)      :: Rmat
!   Rmat =  matmul( conjg(transpose(Umat)), matmul(Op,Umat))
! end function truncate
