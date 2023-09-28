MODULE SYSTEM
  USE SCIFOR
  USE INPUT_VARS
  USE AUX_FUNCS
  USE MATRIX_SPARSE, id=>sp_eye
  USE MATRIX_BLOCKS
  USE TUPLE_BASIS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
  implicit none
  private


  abstract interface
     function Hconnect(left,right,states)
       USE BLOCKS
       USE MATRIX_SPARSE
       type(block)                   :: left
       type(block)                   :: right
       integer,dimension(:),optional :: states
       type(sparse_matrix)           :: Hconnect
     end function Hconnect
  end interface
  procedure(Hconnect),pointer,public :: H2model=>null()
  type(sparse_matrix)                :: spHsb
  type(site)                         :: dot
  real(8),dimension(:),allocatable   :: target_Qn,current_target_QN
  type(block)                        :: sys,env

  public :: init_dmrg
  public :: finalize_dmrg
  public :: idmrg_step

contains


  subroutine init_dmrg(h2,QN,ModelDot)
    procedure(Hconnect)  :: h2
    real(8),dimension(:) :: QN
    type(site)           :: ModelDot
    if(associated(H2model))nullify(H2model); H2model => h2
    allocate(target_qn, source=qn)
    dot       = ModelDot
    !Init block from single dot
    sys=block(dot)
    env=block(dot)
  end subroutine init_dmrg


  subroutine finalize_dmrg()
    if(associated(H2model))nullify(H2model)
    call dot%free()
    call sys%free()
    call env%free()
  end subroutine finalize_dmrg


  subroutine idmrg_step()
    integer                            :: m_sb,m_s,m_e
    integer                            :: m_sys,m_env,Nsys,Nenv
    integer                            :: isys,ienv,isb
    integer                            :: i,j,iqn,Ncv,im,unit,current_L    
    integer,dimension(:),allocatable   :: sb_states
    integer,dimension(:),allocatable   :: sys_map,env_map,sb_map
    real(8),dimension(:),allocatable   :: sb_qn,qn
    real(8),dimension(:),allocatable   :: eig_values
    real(8),dimension(:),allocatable   :: rho_vec
    real(8),dimension(:),allocatable   :: evals,evals_sys,evals_env
    real(8),dimension(:,:),allocatable :: eig_basis,Hsb
    type(blocks_matrix)                :: rho_sys,rho_env
    type(blocks_matrix)                :: Mpsi
    type(tbasis)                       :: sys_basis,env_basis
    type(sparse_matrix)                :: trRho_sys,trRho_env,psiM
    type(sectors_list)                 :: sb_sector
    real(8)                            :: truncation_error_sys,truncation_error_env,corr
    integer                            :: Eunit
    logical                            :: exist
    !
    !Clean screen
    call system('clear')
    !
    !Start DMRG step timer
    call start_timer()
    !
    current_L         = 2*sys%length + 2
    current_target_QN = int(target_qn*current_L/2d0)
    write(LOGfile,*)"SuperBlock Length  =",current_L
    write(LOGfile,*)"Target_QN=",current_target_QN,"filling=",sum(current_target_QN)
    !
    !Check if blocks are valid ones
    if(.not.sys%is_valid(.true.))stop "single_dmrg_step error: sys is not a valid block"
    if(.not.env%is_valid(.true.))stop "single_dmrg_step error: env is not a valid block"
    !
    !Enlarge the Blocks:
    call start_timer()
    call enlarge_block(sys,dot,grow='left')
    call enlarge_block(env,dot,grow='right')
    if(.not.sys%is_valid())stop "single_dmrg_step error: enlarged_sys is not a valid block"
    if(.not.env%is_valid())stop "single_dmrg_step error: enlarged_env is not a valid block"
    call stop_timer("Enlarge blocks")
    !
    !Get Enlarged Sys/Env actual dimension
    m_sys = sys%dim
    m_env = env%dim
    !
    !Build SuperBLock Sector
    call start_timer()
    call get_sb_states(sys,env,sb_states,sb_sector)
    call stop_timer("Get SB states")
    m_sb = size(sb_states)
    !
    write(LOGfile,"(A,I12,I12))")&
         "Enlarged Blocks dimensions           :", sys%dim,env%dim  
    write(LOGfile,"(A,I12)")&
         "SuperBlock Total Dimension           :", m_sys*m_env
    write(LOGfile,"(A,I12)")&
         "SuperBlock Sector Dimension          :", m_sb
    call start_timer()
    write(LOGfile,"(A)")"Building H_sb"
    spHsb = sp_kron(sys%operators%op("H"),id(m_env),sb_states) + &
         sp_kron(id(m_sys),env%operators%op("H"),sb_states)  + &
         H2model(sys,env,sb_states)
    call stop_timer("Done H_sb")
    !
    !
    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    allocate(eig_values(Lanc_Neigen))    ;eig_values=0d0
    allocate(eig_basis(m_sb,Lanc_Neigen));eig_basis =0d0
    call start_timer()
    if(m_sb < lanc_dim_threshold)then
       allocate(Hsb(m_sb,m_sb));Hsb=0d0
       allocate(evals(m_sb))
       call spHsb%dump(Hsb)
       call eigh(Hsb,evals)
       eig_basis(:,1:Lanc_Neigen) = Hsb(:,1:Lanc_Neigen)
       eig_values(1:Lanc_Neigen)  = evals(1:Lanc_Neigen)
       deallocate(Hsb,evals)
    else
       Ncv = min(lanc_ncv_factor*Lanc_Neigen+lanc_ncv_add,m_sb)
       call sp_eigh(sb_HxV,eig_values,eig_basis,&
            Ncv,&
            lanc_niter,&
            tol=lanc_tolerance,&
            iverbose=.false.)
    end if
    call stop_timer("Diag H_sb")
    !
    !
    !BUILD RHO:
    call start_timer()
    do isb=1,size(sb_sector)
       sb_qn  = sb_sector%qn(index=isb)
       sb_map = sb_sector%map(index=isb)
       Nsys   = size(sys%sectors(1)%map(qn=sb_qn))
       Nenv   = size(env%sectors(1)%map(qn=(current_target_qn - sb_qn)))
       if(Nsys*Nenv==0)cycle
       !build PsiMat
       qn = sb_qn
       call Mpsi%append(&
            build_PsiMat(Nsys,Nenv,eig_basis(:,1),sb_map),&
            qn=qn,map=sys%sectors(1)%map(qn=qn))
       !SYSTEM
       qn = sb_qn
       call rho_sys%append(&
            build_density_matrix(Nsys,Nenv,eig_basis(:,1),sb_map,'left'),&
            qn=qn,map=sys%sectors(1)%map(qn=qn))
       !ENVIRONMENT
       qn = current_target_qn-sb_qn
       call rho_env%append(&
            build_density_matrix(Nsys,Nenv,eig_basis(:,1),sb_map,'right'),&
            qn=qn,map=env%sectors(1)%map(qn=qn))
    enddo
    !
    psiM = as_sparse(rho_sys%dump())
    !
    !Build Truncated Density Matrices:
    m_s = min(Mdmrg,m_sys)
    m_e = min(Mdmrg,m_env)
    call rho_sys%eigh(sort=.true.,reverse=.true.)
    call rho_env%eigh(sort=.true.,reverse=.true.)
    !>truncation errors
    evals_sys = rho_sys%evals()
    evals_env = rho_env%evals()
    truncation_error_sys = 1d0 - sum(evals_sys(1:m_s))
    truncation_error_env = 1d0 - sum(evals_env(1:m_e))
    !>truncation-rotation matrices 
    call trRho_sys%init(m_sys,m_s)
    call trRho_env%init(m_env,m_e)
    do im=1,m_s
       rho_vec = rho_sys%evec(m=im)
       sys_map = rho_sys%map(m=im)
       do i=1,size(rho_vec)
          call trRho_sys%insert(rho_vec(i),sys_map(i),im)
       enddo
    enddo
    do im=1,m_e
       rho_vec = rho_env%evec(m=im)
       env_map = rho_env%map(m=im)
       do i=1,size(rho_vec)
          call trRho_env%insert(rho_vec(i),env_map(i),im)
       enddo
    enddo
    !Store all the rotation/truncation matrices:
    call sys%put_omat(str(sys%length),trRho_sys)
    call env%put_omat(str(env%length),trRho_env)
    call stop_timer("Get Rho + O")
    !
    !Renormalize Blocks:
    call start_timer()
    call sys%renormalize(as_matrix(trRho_sys))
    call env%renormalize(as_matrix(trRho_env))
    call stop_timer("Renormalize Blocks")
    !
    !Prepare output and update basis state
    call start_timer()
    do im=1,m_s
       call sys_basis%append( qn=rho_sys%qn(m=im) )
    enddo
    do im=1,m_e
       call env_basis%append( qn=rho_env%qn(m=im) )
    enddo
    call sys%set_basis(basis=sys_basis)
    call env%set_basis(basis=env_basis)
    call stop_timer("SetUp New Basis")
    !
    call start_timer()
    Eunit = fopen(str("energyVSlength_L"//str(Ldmrg)//"_m"//str(Mdmrg)//".dmrg"),append=.true.)
    write(Eunit,*)current_L,eig_values/current_L
    close(Eunit)
    !
    ! print*,shape(SiSj),shape(psiM),shape(eig_basis(:,1))
    ! corr = trace( as_matrix(psiM.m.SiSj) )
    ! print*, corr
    ! write(250,*)current_L,corr
    call stop_timer("Get Observables")
    !
    !
    call stop_timer("dmrg_step")
    write(LOGfile,"(A,L12,12X,L12)")&
         "Truncating                           :",Mdmrg<=m_sys,Mdmrg<=m_env
    write(LOGfile,"(A,I12,12X,I12)")&
         "Truncation Dim                       :",m_s,m_e
    write(LOGfile,"(A,2ES24.15)")&
         "Truncation Errors                    :",truncation_error_sys,truncation_error_env
    write(LOGfile,"(A,"//str(Lanc_Neigen)//"F24.15)")&
         "Energies/L                           :",eig_values/current_L
    write(LOGfile,*)""
    write(LOGfile,*)""
    call wait(500)
    !
    !
    !Clean memory:
    if(allocated(sb_states))deallocate(sb_states)
    if(allocated(sys_map))deallocate(sys_map)
    if(allocated(env_map))deallocate(env_map)
    if(allocated(sb_map))deallocate(sb_map)
    if(allocated(sb_qn))deallocate(sb_qn)
    if(allocated(qn))deallocate(qn)
    if(allocated(rho_vec))deallocate(rho_vec)
    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    if(allocated(evals_sys))deallocate(evals_sys)
    if(allocated(evals_env))deallocate(evals_env)
    call rho_sys%free()
    call rho_env%free()
    call sys_basis%free()
    call env_basis%free()
    call trRho_sys%free()
    call trRho_env%free()
    call sb_sector%free()  
    call spHsb%free()
    !
    return
  end subroutine idmrg_step




  subroutine enlarge_block(self,dot,grow)
    type(block)                  :: self
    type(site)                   :: dot
    character(len=*),optional    :: grow
    character(len=16)            :: grow_
    character(len=:),allocatable :: key
    type(tbasis)                 :: self_basis,dot_basis,enl_basis
    type(sparse_matrix)          :: Hb,Hd,H2
    integer                      :: i
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    if(.not.self%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: Missing self.H operator in the list"
    if(.not.dot%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: Missing dot.H operator in the list"
    !
    !> Update Hamiltonian:
    select case(str(grow_))
    case ("left","l")
       Hb = self%operators%op("H").x.id(dot%dim)
       Hd = id(self%dim).x.dot%operators%op("H")
       H2 = H2model(self,as_block(dot))
    case ("right","r")
       Hb = id(dot%dim).x.self%operators%op("H")
       Hd = dot%operators%op("H").x.id(self%dim)
       H2 = H2model(as_block(dot),self)
    end select
    call self%put_op("H", Hb +  Hd + H2)
    !
    !> Update all the other operators in the list:
    do i=1,size(self%operators)
       key = self%operators%key(index=i)
       if(str(key)=="H")cycle
       select case(str(grow_))
       case ("left","l")
          call self%put_op(str(key), Id(self%dim).x.dot%operators%op(str(key)))
       case ("right","r")
          call self%put_op(str(key), dot%operators%op(str(key)).x.Id(self%dim))
       end select
    enddo
    !
    !> Enlarge dimensions
    self%length = self%length + 1
    self%Dim    = self%Dim*dot%Dim
    !
    !> Enlarge the basis states
    call self%get_basis(self_basis)
    call dot%get_basis(dot_basis)
    !
    select case(str(grow_))
    case ("left","l")
       enl_basis = (self_basis.o.dot_basis)
       call self%set_basis( basis=enl_basis )
    case ("right","r")
       enl_basis = (dot_basis.o.self_basis)
       call self%set_basis( basis=enl_basis )
    end select
    !
    !Free the memory:
    call Hb%free()
    call Hd%free()
    call H2%free()
    call self_basis%free()
    call dot_basis%free()
    call enl_basis%free()
    !
  end subroutine enlarge_block









  subroutine get_sb_states(sys,env,sb_states,sb_sector)
    type(block)                      :: sys,env
    integer,dimension(:),allocatable :: sb_states
    type(sectors_list)               :: sb_sector
    integer                          :: isys,ienv
    integer                          :: i,j,istate
    real(8),dimension(:),allocatable :: sys_qn,env_qn
    integer,dimension(:),allocatable :: sys_map,env_map
    !
    if(allocated(sb_states))deallocate(sb_states)
    !
    do isys=1,size(sys%sectors(1))
       sys_qn  = sys%sectors(1)%qn(index=isys)
       env_qn = current_target_qn - sys_qn
       if(.not.env%sectors(1)%has_qn(env_qn))cycle
       !
       sys_map = sys%sectors(1)%map(qn=sys_qn)
       env_map= env%sectors(1)%map(qn=env_qn)
       !
       do i=1,size(sys_map)
          do j=1,size(env_map)
             istate=env_map(j) + (sys_map(i)-1)*env%Dim
             call append(sb_states, istate)
             call sb_sector%append(qn=sys_qn,istate=size(sb_states))
          enddo
       enddo
    enddo
    !
  end subroutine get_sb_states


  function build_PsiMat(Nsys,Nenv,psi,map) result(psi_mat)
    integer                            :: Nsys,Nenv
    real(8),dimension(:)               :: psi
    integer,dimension(nsys*nenv)       :: map
    real(8),dimension(:,:),allocatable :: psi_mat
    if(allocated(psi_mat))deallocate(psi_mat)
    allocate(psi_mat(nsys,nenv));psi_mat=0d0
    psi_mat = transpose(reshape(psi(map), [nenv,nsys]))
  end function build_psimat


  function build_density_matrix(Nsys,Nenv,psi,map,direction) result(rho)
    integer                            :: Nsys,Nenv
    real(8),dimension(:)               :: psi
    integer,dimension(nsys*nenv)       :: map
    character(len=*)                   :: direction
    real(8),dimension(:,:),allocatable :: rho
    real(8),dimension(nsys,nenv)       :: psi_tmp
    !
    if(allocated(rho))deallocate(rho)
    !
    psi_tmp = transpose(reshape(psi(map), [nenv,nsys]))
    !
    select case(to_lower(str(direction)))
    case ('left','l')
       allocate(rho(nsys,nsys));rho=0d0
       rho  = matmul(psi_tmp,  transpose(psi_tmp) )
    case ('right','r')
       allocate(rho(nenv,nenv));rho=0d0
       rho  = matmul(transpose(psi_tmp), psi_tmp  )
    end select
  end function build_density_matrix



  ! function get_SziSzj(left,right) result(SiSj)
  !   type(operators_list) :: left
  !   type(operators_list) :: right
  !   type(sparse_matrix)  :: Sz1,Sp1
  !   type(sparse_matrix)  :: Sz2,Sp2
  !   type(sparse_matrix)  :: SiSj
  !   Sz1 = left%op("Sz")
  !   Sp1 = left%op("Sp")
  !   Sz2 = right%op("Sz")
  !   Sp2 = right%op("Sp")
  !   print*,shape(Sz1),shape(Sz2)
  !   ! if(present(states))then
  !   !    SiSj = Jx/2d0*sp_kron(Sp1,Sp2%dgr(),states) +  Jx/2d0*sp_kron(Sp1%dgr(),Sp2,states)  + Jp*sp_kron(Sz1,Sz2,states)
  !   ! else
  !   SiSj = Jx/2d0*(Sp1.x.Sp2%dgr()) +  Jx/2d0*(Sp1%dgr().x.Sp2)  + Jp*(Sz1.x.Sz2)
  !   ! endif
  !   call Sz1%free()
  !   call Sp1%free()
  !   call Sz2%free()
  !   call Sp2%free()
  ! end function get_SziSzj


  
  subroutine sb_HxV(Nloc,v,Hv)
    integer                 :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                 :: i,j,jcol
    Hv=0d0
    do i=1,Nloc
       matmul: do jcol=1, spHsb%row(i)%Size
          val = spHsb%row(i)%vals(jcol)
          j   = spHsb%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end subroutine sb_HxV

END MODULE SYSTEM





