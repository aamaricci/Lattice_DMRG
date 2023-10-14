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
  real(8),dimension(:),allocatable   :: eig_values
  real(8),dimension(:,:),allocatable :: eig_basis
  type(block)                        :: init_sys,init_env
  logical                            :: init_called=.false.
  integer,dimension(:),allocatable   :: sb_states

  public :: init_dmrg
  public :: finalize_dmrg
  public :: step_dmrg
  !
  public :: infinite_DMRG
  public :: finite_DMRG

contains

  !j  = L+1,.....L+R
  !j-L= 1,...,R
  !R+1-(j-L)= R+1-(1),....R+1-(R)
  function measure_dmrg(left,right,op,pos) result(avOp)
    type(block),intent(in)             :: left,right
    type(sparse_matrix)                :: op
    integer                            :: pos,j
    !
    character(len=1)                   :: label
    real(8)                            :: avOp
    type(sparse_matrix)                :: Oj,Otmp
    type(sparse_matrix)                :: U,IR,IL
    integer                            :: it,dims(2),dim,L,R,Nsys,Nenv
    real(8),dimension(:),allocatable   :: Vpsi
    real(8),dimension(:,:),allocatable   :: VMat
    !
    L = left%length
    R = right%length
    if(pos<1.OR.pos>L+R)stop "measure_dmrg error: Pos not in [1,Lchain]"
    !
    label='l';if(pos>L)label='r'
    !
    j=pos
    if(pos>L)j=R+1-(pos-L)            !R,...,1
    !
    !Retrieve the dimensions of the initial Dot
    select case(label)
    case ("l","L");dims = shape(left%omatrices%op(index=1)) ;dim=dims(1)
    case ("r","R");dims = shape(right%omatrices%op(index=1));dim=dims(1)
    end select
    !
    !Step 1: build the operator at the site j
    Oj  = Op                    !.x.Id(dim)
    if(j>1)then
       select case(label)
       case ("l","L")
          U    = left%omatrices%op(index=j-1)
          Otmp = matmul(U%t(),U)
          Oj   = Otmp.x.Op
       case ("r","R")
          U    = right%omatrices%op(index=j-1)
          Otmp = matmul(U%t(),U)
          Oj   = Op.x.Otmp
       end select
    end if
    !
    !Step 2: bring the operator to the actual basis
    select case(label)
    case ("l","L")
       do it=j,L
          U  = left%omatrices%op(index=it)
          Oj = matmul(U%t(),matmul(Oj,U))
          Oj = Oj.x.Id(dim)
       enddo
    case ("r","R")
       do it=j,R
          U  = right%omatrices%op(index=it)
          Oj = matmul(U%t(),matmul(Oj,U))
          Oj = Id(dim).x.Oj
       enddo
    end select
    !
    !Step 3 evaluate the average:
    Otmp = Oj
    select case(label)
    case ("l","L")
       U  = right%omatrices%op(index=R)
       IR = Id(dim).x.matmul(U%t(),U)
       Oj = Oj.x.IR
    case ("r","R")
       U  = left%omatrices%op(index=L)
       IL = matmul(U%t(),U).x.Id(dim)
       Oj = IL.x.Oj
    end select
    !
    print*,shape(Oj)


    allocate(Vpsi(Oj%Ncol));Vpsi=0d0
    Vpsi(sb_states) = eig_basis(:,1)
    AvOp = dot_product(Vpsi,Oj%dot(Vpsi))


    print*,shape(Otmp);dims=shape(Otmp)
    Nsys=dims(1)
    Nenv=dims(2)
    allocate(Vmat(Nsys,Nenv));Vmat=0d0
    Vmat= transpose(reshape(Vpsi, [Nenv,Nsys]))
    print*,shape(Vmat)
    print*,trace(matmul(transpose(Vmat),matmul(as_matrix(Otmp),Vmat)))


    deallocate(Vpsi,Vmat)
    !
  end function measure_dmrg



  subroutine infinite_DMRG()
    type(block) :: sys,env
    type(sparse_matrix) :: Op
    real(8) :: val
    integer :: i
    if(.not.init_called)stop "infinite_DMRG ERROR: DMRG not initialized. Call init_dmrg first."
    sys=init_sys
    env=init_env
    do while (2*sys%length < Ldmrg)
       call step_dmrg(sys,env)
       call get_observables(sys,env,suffix="L"//str(Ldmrg)//"_iDMRG")
    enddo
    Op = dot%operators%op(key='Sz')
    do i=1,Ldmrg/2
       val= measure_dmrg(sys,env,Op,i)
       print*,i,val
       write(100,*)i,val
    enddo
  end subroutine infinite_DMRG



  subroutine finite_DMRG()
    type(block)                            :: sys,env
    integer                                :: i,im,env_label,sys_label,current_L
    type(block),dimension(:,:),allocatable :: blocks_list
    type(block)                            :: tmp
    !
    if(mod(Ldmrg,2)/=0)stop "finite_DMRG ERROR: Ldmrg%2 != 0. Ldmrg input must be an even number."
    !
    if(.not.init_called)stop "finite_DMRG ERROR: DMRG not initialized. Call init_dmrg first."
    !
    sys=init_sys
    env=init_env
    !
    sys_label=1
    env_label=2
    !Infinite DMRG
    allocate(blocks_list(2,Ldmrg))
    blocks_list(sys_label,1)=sys
    blocks_list(env_label,1)=env
    do while (2*sys%length < Ldmrg)
       call step_dmrg(sys,env)
       call get_observables(sys,env,suffix="L"//str(Ldmrg)//"_iDMRG")
       blocks_list(sys_label,sys%length)=sys
       blocks_list(env_label,env%length)=env
    enddo
    print*,""
    print*,""
    print*,"START FINITE DMRG:"
    print*,""
    print*,""
    !Finite DMRG: start sweep forth&back:
    do im=1,size(Msweep)
       write(*,"(A,I3,I6)")"Sweep, M:",im,Msweep(im)
       sweep: do while(.true.)
          env = blocks_list(env_label,Ldmrg - sys%length)
          if(env%length==1)then          !end of the chain for env, switch roles sys<-->env
             env_label= 3-env_label
             sys_label= 3-sys_label
             tmp      = sys
             sys      = env
             env      = tmp
             call tmp%free()
          endif
          !
          call step_dmrg(sys,env,sys_label,im)
          call get_observables(sys,env,suffix="L"//str(Ldmrg)//"_sweep"//str(im))
          !
          blocks_list(sys_label,sys%length) = sys
          print*,""
          print*,""
          if(sys_label==1.AND.sys%length==Ldmrg/2)exit sweep
       enddo sweep
    enddo
    !
  end subroutine finite_DMRG





  subroutine init_dmrg(h2,QN,ModelDot)
    procedure(Hconnect)  :: h2
    real(8),dimension(:) :: QN
    type(site)           :: ModelDot
    if(associated(H2model))nullify(H2model); H2model => h2
    allocate(target_qn, source=qn)
    dot        = ModelDot
    init_sys   = block(dot)
    init_env   = block(dot)
    init_called=.true.
  end subroutine init_dmrg


  subroutine finalize_dmrg()
    if(associated(H2model))nullify(H2model)
    call dot%free()
    call init_sys%free()
    call init_env%free()
    init_called=.false.
  end subroutine finalize_dmrg


  subroutine step_dmrg(sys,env,label,isweep)
    type(block),intent(inout)          :: sys,env
    integer,optional                   :: label,isweep
    integer                            :: iLabel
    integer                            :: Mstates
    real(8)                            :: Estates
    integer                            :: m_sb,m_s,m_e
    integer                            :: m_sys,m_env,Nsys,Nenv
    integer                            :: isys,ienv,isb,m_err(1),m_threshold
    integer                            :: i,j,iqn,Ncv,im,unit,current_L    

    integer,dimension(:),allocatable   :: sys_map,env_map,sb_map
    real(8),dimension(:),allocatable   :: sb_qn,qn
    real(8),dimension(:),allocatable   :: rho_vec
    real(8),dimension(:),allocatable   :: evals,evals_sys,evals_env
    real(8),dimension(:,:),allocatable :: Hsb
    type(blocks_matrix)                :: rho_sys,rho_env
    type(tbasis)                       :: sys_basis,env_basis
    type(sparse_matrix)                :: trRho_sys,trRho_env,psiM
    type(sectors_list)                 :: sb_sector
    real(8)                            :: truncation_error_sys,truncation_error_env
    logical                            :: exist
    !
    iLabel=1;if(present(label))iLabel=label
    Mstates=Mdmrg
    Estates=Edmrg
    if(present(isweep))then
       if(isweep>Nsweep)stop "step_dmrg ERROR: isweep > Nsweep"
       if(isweep<1)stop "step_dmrg ERROR: isweep < 1"
       Mstates = Msweep(isweep)
       Estates = Esweep(isweep)
    endif
    !
    call dmrg_graphic(sys,env,iLabel)
    !
    !Start DMRG step timer
    call start_timer()
    !
    current_L         = sys%length + env%length + 2
    current_target_QN = int(target_qn*current_L/2d0)
    write(LOGfile,"(A,I12)")"SuperBlock Length  =",current_L
    write(LOGfile,"(A,"//str(size(current_target_QN))//"F12.7,2x,A,F12.7)")&
         "Target_QN=",current_target_QN,"filling=",sum(current_target_QN)
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
    write(LOGfile,"(A,I12,A1,I12,A1,F10.5,A1)")&
         "SuperBlock Dimension  (tot)          :", m_sb,"(",m_sys*m_env,")",dble(m_sb)/m_sys/m_env,"%"
    call start_timer()
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
    !BUILD RHO:
    call start_timer()
    do isb=1,size(sb_sector)
       sb_qn  = sb_sector%qn(index=isb)
       sb_map = sb_sector%map(index=isb)
       Nsys   = size(sys%sectors(1)%map(qn=sb_qn))
       Nenv   = size(env%sectors(1)%map(qn=(current_target_qn - sb_qn)))
       if(Nsys*Nenv==0)cycle
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
    !Build Truncated Density Matrices:
    call rho_sys%eigh(sort=.true.,reverse=.true.)
    call rho_env%eigh(sort=.true.,reverse=.true.)
    !>truncation errors
    evals_sys = rho_sys%evals()
    evals_env = rho_env%evals()
    m_s = min(Mstates,m_sys,size(evals_sys))
    m_e = min(Mstates,m_env,size(evals_env))
    if(Estates/=0d0)then!< ACTHUNG: treatment of degenerate rho-eigenvalues is recommended 
       !find the value of m such that 1-sum_i lambda_i < Edmrg
       m_err = minloc(abs(1d0-cumulate(evals_sys)-Estates))
       m_s   = m_err(1)
       m_err = minloc(abs(1d0-cumulate(evals_env)-Estates))
       m_e   = m_err(1)
    endif
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
    !
    !Renormalize Blocks:
    call sys%renormalize(as_matrix(trRho_sys))
    call env%renormalize(as_matrix(trRho_env))
    !
    !Prepare output and update basis state
    do im=1,m_s
       call sys_basis%append( qn=rho_sys%qn(m=im) )
    enddo
    do im=1,m_e
       call env_basis%append( qn=rho_env%qn(m=im) )
    enddo
    call sys%set_basis(basis=sys_basis)
    call env%set_basis(basis=env_basis)
    call stop_timer("Get Rho, U->Renormalize, Setup New Basis")
    !
    !
    !
    write(LOGfile,"(A,L12,12X,L12)")&
         "Truncating                           :",Mstates<=m_sys,Mstates<=m_env
    write(LOGfile,"(A,I12,12X,I12)")&
         "Truncation Dim                       :",m_s,m_e
    write(LOGfile,"(A,2ES24.15)")&
         "Truncation Errors                    :",truncation_error_sys,truncation_error_env
    write(LOGfile,"(A,"//str(Lanc_Neigen)//"F24.15)")&
         "Energies/L                           :",eig_values/current_L
    call stop_timer("dmrg_step")
    !
    !
    !Clean memory:
    ! if(allocated(sb_states))deallocate(sb_states)
    if(allocated(sys_map))deallocate(sys_map)
    if(allocated(env_map))deallocate(env_map)
    if(allocated(sb_map))deallocate(sb_map)
    if(allocated(sb_qn))deallocate(sb_qn)
    if(allocated(qn))deallocate(qn)
    if(allocated(rho_vec))deallocate(rho_vec)
    ! if(allocated(eig_values))deallocate(eig_values)
    ! if(allocated(eig_basis))deallocate(eig_basis)
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
  end subroutine step_dmrg



  subroutine dmrg_graphic(sys,env,label)
    type(block) :: sys,env
    integer     :: label,M
    integer     :: i
    write(LOGfile,*)""
    write(LOGfile,*)""
    select case(label)
    case default; stop "dmrg_graphic error: label != 1(L),2(R)"
    case(1)
       write(LOGfile,"(A,2I4,A6,1x)",advance="no")"sys.L + env.L:",sys%length,env%length,"|"
       write(LOGfile,"("//str(sys%length)//"A1)",advance="no")("=",i=1,sys%length)
       write(LOGfile,"(A2)",advance="no")"**"
       write(LOGfile,"("//str(env%length)//"A1)",advance="no")("-",i=1,env%length)
    case(2)
       write(LOGfile,"(A,2I4,A6,1x)",advance="no")"sys.L + env.L:",env%length,sys%length,"|"
       write(LOGfile,"("//str(env%length)//"A1)",advance="no")("=",i=1,env%length)
       write(LOGfile,"(A2)",advance="no")"**"
       write(LOGfile,"("//str(sys%length)//"A1)",advance="no")("-",i=1,sys%length)
    end select
    write(LOGfile,"(1x,A1,2I6)",advance='yes')"|",Ldmrg
  end subroutine dmrg_graphic



  subroutine get_observables(sys,env,suffix)
    type(block),intent(inout) :: sys,env
    character(len=*)          :: suffix
    integer                   :: current_L
    integer                   :: Eunit
    !
    current_L         = sys%length + env%length
    !
    call start_timer()
    Eunit = fopen("energyVSsys.length_"//str(suffix)//".dmrg",append=.true.)
    write(Eunit,*)sys%length,eig_values/current_L
    close(Eunit)
    call stop_timer("Get Observables")
  end subroutine get_observables


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
       env_qn  = current_target_qn - sys_qn
       if(.not.env%sectors(1)%has_qn(env_qn))cycle
       !
       sys_map = sys%sectors(1)%map(qn=sys_qn)
       env_map = env%sectors(1)%map(qn=env_qn)
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





