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
  real(8)                            :: ts
  integer                            :: Lmax,i,m,Neigen,current_L
  integer                            :: lanc_ncv_factor
  integer                            :: lanc_ncv_add
  integer                            :: lanc_niter
  real(8)                            :: lanc_tolerance
  integer                            :: lanc_threshold
  type(block)                        :: my_block,dimer,trimer
  type(sparse_matrix)                :: spHsb,spH
  type(site)                         :: dot
  real(8)                            :: gs_energy,target_Sz
  integer                            :: unit
  integer                            :: model_d=4

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"ts",finput,default=-1d0,comment="Hopping amplitude")
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
  dot = hubbard_site(0d0,0d0)
  !
  !Init block from single dot
  my_block=block(dot)
  call my_block%show()



  !Run DMRG algorithm
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






  function enlarge_block(self,dot) result(enl_self)
    type(block),intent(inout)        :: self
    type(site),intent(inout)         :: dot
    type(block)                      :: enl_self
    integer                          :: mblock,len
    real(8),dimension(:),allocatable :: self_basis, dot_basis
    real(8),dimension(4,4)           :: Cup,Cdw,P
    !
    P = kSz(2)
    Cup =  as_matrix(dot%operators%op("Cup"))
    Cdw =  as_matrix(dot%operators%op("Cdw"))
    !
    mblock =  self%dim
    len    =  self%length
    !
    enl_self%length = len + 1
    enl_self%dim    = mblock*model_d
    !
    call enl_self%put("H", (self%operators%op("H").x.id(model_d)) + (id(mblock).x.dot%operators%op("H")) + &
         H2model(self,as_block(dot)))
    call enl_self%put("Cup", Id(mblock).x.as_sparse(matmul(P,Cup)) ) !Cup = Id(Mblock).x.(P_dot.Cup_dot)
    call enl_self%put("Cdw", Id(mblock).x.as_sparse(matmul(P,Cdw)) ) !Cdw = Id(Mblock).x.(P_dot.Cdw_dot)
    !
    self_basis = self%sectors%basis()
    dot_basis  = dot%sectors%basis()
    call enl_self%set_sectors( outsum(self_basis,dot_basis) )
    !
    deallocate(self_basis,dot_basis)
  end function enlarge_block


  function H2model(left,right) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    type(sparse_matrix)           :: CupL,CdwL
    type(sparse_matrix)           :: CupR,CdwR
    type(sparse_matrix)           :: H2
    CupL = left%operators%op("Cup")
    CdwL = left%operators%op("Cdw")
    CupR = right%operators%op("Cup")
    CdwR = right%operators%op("Cdw")
    !H_lr        = H_lr(up) + H_lr(dw)
    !H_lr(sigma) = -t( C^+_l,sigma . C_r,sigma)  + H.c.
    H2 = ts*(CupL%dgr().x.CupR) + ts*(CdwL%dgr().x.CdwR)  
    H2 = H2 + H2%dgr()
    call CupL%free
    call CdwL%free
    call CupR%free
    call CdwR%free
  end function H2model

  
  subroutine print_mat(M,name,n)
    real(8),dimension(:,:) :: M
    character(len=*) :: name
    integer :: i,j,stride,unit,n
    stride=2**n
    open(free_unit(unit),file=str(name)//".dat")
    write(unit,*)"Matrix: "//str(name)
    do i=1,size(M,1)
       do j=1,size(M,2)
          write(unit,"(("//str(stride)//"I2))",advance="no")int(M(i,j))
          if(mod(j,stride)==0)write(unit,"(A1)",advance="no")""
       enddo
       write(unit,*)
       if(mod(i,stride)==0)write(unit,*)
    enddo
    write(unit,*)""
    close(unit)
  end subroutine print_mat




  subroutine single_dmrg_step(m,sys,energy)
    integer                            :: m
    class(block)                       :: sys
    type(block)                        :: esys,eenv
    type(sectors_list)                 :: sb_sector
    real(8)                            :: esys_qn,eenv_qn,sb_qn
    integer,dimension(:),allocatable   :: esys_map,eenv_map,sb_states,sb_map
    integer                            :: m_esys,m_sb,m_,i,j,Ncv,isys,Nsys,Nenv,im
    real(8),dimension(:),allocatable   :: eig_values,evals,sys_basis
    real(8),dimension(:,:),allocatable :: eig_basis
    real(8),dimension(:),allocatable   :: rho_vec
    type(blocks_matrix)                :: rho_sector
    real(8),dimension(:,:),allocatable :: rho,truncation_rho
    real(8)                            :: truncation_error,energy

    !Check if blocks are valid ones
    if(.not.sys%is_valid())stop "single_dmrg_step error: sys is not a valid block"
    !
    !Enlarge blocks
    esys = enlarge_block(sys,dot)
    !
    if(.not.esys%is_valid())stop "single_dmrg_step error: enlarged_sys is not a valid block"
    eenv = esys
    !
    !Get Enlarged Sys/Env actual dimension
    m_esys = esys%dim

    !Build SuperBLock Sector
    call start_timer()
    m_sb = m_esys*m_esys
    print*,"Super Block dimensions:", m_sb
    !BUild SuperBlock matrix:
    spHsb =         sp_kron(esys%operators%op("H"),id(m_esys))
    spHsb = spHsb + sp_kron(id(m_esys),esys%operators%op("H"))
    spHsb = spHsb + H2model(esys,esys)
    call stop_timer("Build H_sb")


    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    allocate(eig_values(Neigen))    ;eig_values=0d0
    allocate(eig_basis(m_sb,Neigen));eig_basis =0d0
    call start_timer()
    if(m_sb < lanc_threshold)then
       print*,"diag SB w/ Lapack"
       allocate(rho(m_sb,m_sb));rho=zero
       allocate(evals(m_sb))
       call spHsb%dump(rho)
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


    !find a clever way to iterate here:
    !Return updated Block:
    call start_timer()
    sys%length = esys%length
    sys%dim    = m_
    call sys%put("H",rotate_and_truncate(esys%operators%op("H"),truncation_rho,m_esys,m_))
    call sys%put("Sz",rotate_and_truncate(esys%operators%op("Sz"),truncation_rho,m_esys,m_))
    call sys%put("Sp",rotate_and_truncate(esys%operators%op("Sp"),truncation_rho,m_esys,m_))
    call sys%set_sectors( sys_basis )
    call stop_timer("Update+Truncate Block")


    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    if(allocated(evals))deallocate(evals)
    if(allocated(rho))deallocate(rho)
    call esys%free()
    call eenv%free()
    call spHsb%free()

  end subroutine single_dmrg_step





  function build_density_matrix(nsys,nenv,psi) result(rho)
    integer                      :: nsys
    integer                      :: nenv
    real(8),dimension(nsys,nenv) :: psi
    real(8),dimension(nsys,nsys) :: rho
    rho  = matmul(psi,  (transpose(psi)) )
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


end program
