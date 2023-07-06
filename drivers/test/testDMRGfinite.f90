program testDMRGinfty
  USE SCIFOR
  USE SPARSE_MATRICES, id=>sp_eye
  USE OPERATORS_TUPLE
  USE SITES
  USE BLOCKS
  implicit none

  character(len=64)                      :: finput
  real(8)                                :: Jx,Jz
  integer                                :: Lmax,i,m,Neigen
  integer                                :: lanc_ncv_factor
  integer                                :: lanc_ncv_add
  integer                                :: lanc_niter
  real(8)                                :: lanc_tolerance
  type(block)                            :: my_block,sys,env
  type(sparse_matrix)                    :: Szeta,Splus,Hzero,Sminus
  type(sparse_matrix)                    :: sb_hamiltonian
  real(8)                                :: gs_energy
  integer                                :: unit
  integer,parameter                      :: model_d=2
  type(block),dimension(:,:),allocatable :: block_list
  integer                                :: sys_label,env_label

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(Jx,"Jx",finput,default=1d0,comment="Jx term in the XXZ model")
  call parse_input_variable(Jz,"Jz",finput,default=1d0,comment="Jz term in the XXZ model")
  call parse_input_variable(Lmax,"LMAX",finput,default=20,comment="Final chain length ")
  call parse_input_variable(m,"M",finput,default=10,&
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



  !RUN DMRG:
  open(free_unit(unit),file="energyVSlength_L"//str(Lmax)//"_m"//str(m)//".dmrg")
  my_block=block(spin_onehalf_site())
  Hzero=my_block%odict%op("H")
  Szeta=my_block%odict%op("Sz")
  Splus=my_block%odict%op("Sp")


  allocate(block_list(2,Lmax))
  block_list(1,my_block%length) = my_block
  block_list(2,my_block%length) = my_block
  do while (2*my_block%length < Lmax)
     write(*,*)"Block  Length  =",my_block%length
     write(*,*)"System Length  =",2*my_block%length + 2
     call single_dmrg_step(M,my_block,energy=gs_energy)
     write(*,*)"E/L=",gs_energy/(my_block%length*2)
     write(unit,*)my_block%length*2,gs_energy/(my_block%length*2),gs_energy
     !Save blocks
     block_list(1,my_block%length) = my_block
     block_list(2,my_block%length) = my_block
     write(*,*)""
  enddo
  close(unit)

  call sleep(2)
  print*,"Starting Finite Algorithm"
  sys_label=1
  env_label=2
  sys = my_block; call my_block%free
  do m=10,40,10
     print*,"M:",m
     call sleep(1)
     sweep: do while(.true.)
        env = block_list(env_label,Lmax-sys%length-2)
        if(env%length==1)then
           sys_label=3-sys_label !exchange 1<-->2
           env_label=3-env_label !exchange 1<-->2
           !Sys <--> Env
           my_block = sys
           sys      = env
           env      = my_block
           call my_block%free()
        endif

        print*,sys%length,env%length
        call dmrg_graphic(sys,env,sys_label)
        call single_dmrg_step(M,sys,env,energy=gs_energy)
        write(*,*)"E/L=",gs_energy/(Lmax),gs_energy
        write(100,*)gs_energy/(Lmax),gs_energy

        block_list(sys_label,sys%length) = sys
        if(sys_label==1.AND.2*sys%length==Lmax)exit sweep        
     enddo sweep

  enddo

contains

  function H2model(Sz1,Sp1,Sz2,Sp2) result(H2)
    type(sparse_matrix),intent(in) :: Sz1,Sp1
    type(sparse_matrix),intent(in) :: Sz2,Sp2
    type(sparse_matrix)            :: H2
    H2 = Jx/2d0*(Sp1.x.Sp2%dgr()) +  Jx/2d0*(Sp1%dgr().x.Sp2)  + Jz*(Sz1.x.Sz2)
  end function H2model


  function enlarge_block(self) result(enl_self)
    class(block),intent(inout) :: self
    type(block)                :: enl_self
    integer                    :: mblock,len
    enl_self = self
    mblock =  self%dim
    len    =  self%length
    call enl_self%put("H", (self%odict%op("H").x.id(model_d)) + (id(mblock).x.Hzero) + &
         H2model(self%odict%op("Sz"), self%odict%op("Sp"), Szeta, Splus) )
    call enl_self%put("Sz", id(mblock).x.Szeta)
    call enl_self%put("Sp", id(mblock).x.Splus)
    enl_self%length = self%length + 1
    enl_self%dim    = mblock*model_d
  end function enlarge_block




  subroutine single_dmrg_step(m,sys,env,energy)
    integer                               :: m
    class(block)                          :: sys
    class(block),optional                 :: env
    !
    type(block)                           :: sys_enl
    type(block)                           :: env_enl
    !
    integer                               :: m_esys
    integer                               :: m_eenv
    integer                               :: m_sb
    integer                               :: m_
    integer                               :: i,Ncv
    !
    real(8),dimension(:),allocatable      :: eig_values,evals
    complex(8),dimension(:,:),allocatable :: eig_basis,evecs
    complex(8),dimension(:,:),allocatable :: rho,truncation_rho
    real(8)                               :: truncation_error,energy

    !Check if blocks are valid ones
    if(.not.sys%is_valid())stop "single_dmrg_step error: sys is not a valid block"
    if(present(env))then
       if(.not.env%is_valid())stop "single_dmrg_step error: env is not a valid block"
    endif
    !
    !Enlarge blocks
    sys_enl = enlarge_block(sys)
    if(present(env))then
       env_enl = enlarge_block(env)
    else
       env_enl = sys_enl
    endif
    !
    if(.not.sys_enl%is_valid())stop "single_dmrg_step error: enlarged_sys is not a valid block"
    if(.not.env_enl%is_valid())stop "single_dmrg_step error: enlarged_env is not a valid block"
    !Get Enlarged Sys/Env actual dimension
    m_esys = sys_enl%dim
    m_eenv = env_enl%dim
    print*,"Enlarged Block Dimension:",m_esys
    !
    !Build SuperBlock Hamiltonian
    m_sb = m_esys*m_eenv
    print*,"Super Block dimension:",m_sb
    sb_hamiltonian = (&
         sys_enl%odict%op("H").x.id(m_eenv))  + (id(m_esys).x.env_enl%odict%op("H")) + &
         H2model(&
         sys_enl%odict%op("Sz"),sys_enl%odict%op("Sp"),&
         env_enl%odict%op("Sz"),env_enl%odict%op("Sp"))


    allocate(eig_values(Neigen))    ;eig_values=0d0
    allocate(eig_basis(m_sb,Neigen));eig_basis =zero
    if(m_sb < 1024)then
       print*,"diag SB w/ Lapack"
       allocate(rho(m_sb,m_sb));rho=zero
       allocate(evals(m_sb))
       call sb_hamiltonian%dump(rho)
       call eigh(rho,evals)
       eig_basis(:,1:Neigen) = rho(:,1:Neigen)
       eig_values(1:Neigen)  = evals(1:Neigen)
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


    print*,"DM dimension:",m_esys
    allocate(rho(m_esys,m_esys));rho=zero

    allocate(evals(m_esys))
    rho = build_density_matrix(m_esys,m_eenv,eig_basis(:,1))
    call eigh(rho,evals)

    m_ = min(m,m_esys);
    allocate(truncation_rho(m_esys,m_))
    truncation_rho(:,1:m_) = rho(:,m_esys:m_esys-m_+1:-1)
    truncation_error = 1d0 - sum(evals(m_esys:m_esys-m_+1:-1))
    write(*,*)"Truncation dim  :",m_
    write(*,*)"Truncation error:",truncation_error
    print*,"E0=",eig_values(:)

    !Return updated Block:
    sys%length = sys_enl%length
    sys%dim    = m_
    !find a clever way to iterate here:
    call sys%put("H",rotate_and_truncate(sys_enl%odict%op("H"),truncation_rho,m_esys,m_))
    call sys%put("Sz",rotate_and_truncate(sys_enl%odict%op("Sz"),truncation_rho,m_esys,m_))
    call sys%put("Sp",rotate_and_truncate(sys_enl%odict%op("Sp"),truncation_rho,m_esys,m_))


    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    if(allocated(evals))deallocate(evals)
    if(allocated(rho))deallocate(rho)
    call sb_hamiltonian%free()
    call sys_enl%free()
    call env_enl%free()
  end subroutine single_dmrg_step



  function build_density_matrix(nsys,nenv,psi) result(rho)
    integer                         :: nsys
    integer                         :: nenv
    complex(8),dimension(nsys,nenv) :: psi
    complex(8),dimension(nsys,nsys) :: rho
    rho  = matmul(psi,  conjg(transpose(psi)) )
  end function build_density_matrix


  function rotate_and_truncate(Op,trRho,N,M) result(RotOp)
    type(sparse_matrix),intent(in)        :: Op
    complex(8),dimension(N,M)             :: trRho  ![Nesys,M]
    integer                               :: N,M
    type(sparse_matrix)                   :: RotOp
    complex(8),dimension(M,M)             :: Umat
    N = size(trRho,1)
    M = size(trRho,2)
    if( any( shape(Op) /= [N,N] ) ) stop "rotate_and_truncate error: shape(Op) != [N,N] N=size(Rho,1)"
    Umat = matmul( conjg(transpose(trRho)), matmul(Op%as_matrix(),trRho)) ![M,N].[N,N].[N,M]=[M,M]
    call RotOp%load( Umat )
  end function rotate_and_truncate


  subroutine sb_HxV(Nloc,v,Hv)
    integer                    :: Nloc
    complex(8),dimension(Nloc) :: v
    complex(8),dimension(Nloc) :: Hv
    complex(8)                 :: val
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



  subroutine dmrg_graphic(sys,env,label)
    type(block) :: sys,env
    integer     :: label
    integer     :: i
    select case(label)
    case default; stop "dmrg_graphic error: label != 1(L),2(R)"
    case(1)
       write(*,"("//str(sys%length)//"A1)",advance="no")("=",i=1,sys%length)
       write(*,"(A2)",advance="no")"**"
       write(*,"("//str(env%length)//"A1)",advance="yes")("-",i=1,env%length)
    case(2)
       write(*,"("//str(env%length)//"A1)",advance="no")("=",i=1,env%length)
       write(*,"(A2)",advance="no")"**"
       write(*,"("//str(sys%length)//"A1)",advance="yes")("-",i=1,sys%length)
    end select
  end subroutine dmrg_graphic




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


end program testDMRGinfty






! function truncate(Op,Umat,N,M) result(Rmat)
!   complex(8),dimension(N,N)      :: Op
!   integer                        :: N,M
!   complex(8),dimension(N,M)      :: Umat
!   complex(8),dimension(M,M)      :: Rmat
!   Rmat =  matmul( conjg(transpose(Umat)), matmul(Op,Umat))
! end function truncate
