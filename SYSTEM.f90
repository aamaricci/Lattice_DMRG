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
  procedure(Hconnect),pointer,public    :: H2model=>null()
  !
  type(sparse_matrix)                   :: spHsb
  !
  type(site)                            :: dot
  character(len=:),allocatable          :: suffix
  real(8),dimension(:),allocatable      :: target_Qn,current_target_QN
  type(block)                           :: init_sys,init_env
  logical                               :: init_called=.false.
  real(8),dimension(:),allocatable      :: energies
  real(8),dimension(:),allocatable      :: rho_sys_evals
  real(8),dimension(:),allocatable      :: rho_env_evals
  type(blocks_matrix)                   :: psi_sys,rho_sys
  type(blocks_matrix)                   :: psi_env,rho_env
  !GLOBAL SYS & ENV 
  type(block)                           :: sys,env


  public :: init_DMRG
  public :: finalize_DMRG
  public :: step_DMRG
  !
  public :: infinite_DMRG
  public :: finite_DMRG

  !Measuring:
  public :: Measure_Op_DMRG
  public :: Build_Op_DMRG
  public :: Advance_Op_DMRG
  public :: Advance_Corr_DMRG
  public :: Average_Op_DMRG



contains




  subroutine init_dmrg(h2,QN,ModelDot)
    procedure(Hconnect)          :: h2
    real(8),dimension(:)         :: QN
    type(site)                   :: ModelDot
    if(associated(H2model))nullify(H2model)
    H2model => h2
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
    call sys%free()
    call env%free()
    call psi_sys%free()
    call psi_env%free()
    init_called  = .false.
  end subroutine finalize_dmrg




  subroutine infinite_DMRG()
    !
    if(.not.init_called)stop "infinite_DMRG ERROR: DMRG not initialized. Call init_dmrg first."
    sys=init_sys
    env=init_env
    !
    suffix = label_DMRG('i',1)
    !
    do while (2*sys%length <= Ldmrg)
       call step_dmrg()
       call write_energy()
       call write_entanglement()
    enddo
  end subroutine infinite_DMRG





  subroutine finite_DMRG()
    integer                                :: i,im,env_label,sys_label,current_L
    type(block),dimension(:,:),allocatable :: blocks_list
    type(block)                            :: tmp
    real(8),allocatable                    :: val(:)
    integer                                :: j,Ncorr
    type(sparse_matrix),allocatable        :: Ocorr(:)
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
    !
    suffix = label_DMRG('i',1)
    !
    allocate(blocks_list(2,Ldmrg))
    blocks_list(sys_label,1)=sys
    blocks_list(env_label,1)=env
    do while (2*sys%length <= Ldmrg)
       call step_dmrg()
       call write_energy()
       call write_entanglement()
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
       if(Esweep(im)==0)then
          write(*,"(A,I3,I6)")"Sweep, M:",im,Msweep(im)
       else
          write(*,"(A,I3,F8.4)")"Sweep, E:",im,Esweep(im)
       endif
       suffix = label_DMRG('f',im)
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
          call step_dmrg(sys_label,im)
          call write_energy()
          call write_entanglement()
          !
          blocks_list(sys_label,sys%length) = sys
          print*,""
          print*,""
          if(sys_label==1.AND.sys%length==Ldmrg/2)exit sweep
       enddo sweep
    enddo
    !
  end subroutine finite_DMRG








  subroutine step_dmrg(label,isweep)
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
    real(8),dimension(:),allocatable   :: evals
    complex(8),dimension(:,:),allocatable :: Hsb
    real(8),dimension(:),allocatable   :: eig_values
    complex(8),dimension(:,:),allocatable :: eig_basis
    integer,dimension(:),allocatable   :: sb_states
    type(tbasis)                       :: sys_basis,env_basis
    type(sparse_matrix)                :: trRho_sys,trRho_env
    type(sectors_list)                 :: sb_sector
    real(8)                            :: truncation_error_sys,truncation_error_env
    logical                            :: exist
    !
    iLabel=0;if(present(label))iLabel=label
    Mstates=Mdmrg
    Estates=Edmrg
    if(present(isweep))then
       if(isweep>Nsweep)stop "step_dmrg ERROR: isweep > Nsweep"
       if(isweep<1)stop "step_dmrg ERROR: isweep < 1"
       Mstates = Msweep(isweep)
       Estates = Esweep(isweep)
    endif
    !
    call dmrg_graphic(iLabel)
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
    call get_sb_states(sb_states,sb_sector)
    call stop_timer("Get SB states")
    m_sb = size(sb_states)
    !
    write(LOGfile,"(A,I12,I12))")&
         "Enlarged Blocks dimensions           :", sys%dim,env%dim  
    write(LOGfile,"(A,I12,A1,I12,A1,F10.5,A1)")&
         "SuperBlock Dimension  (tot)          :", m_sb,"(",m_sys*m_env,")",dble(m_sb)/m_sys/m_env,"%"


    !Build SuperBLock Hamiltonian
    call start_timer()
    spHsb = sp_kron(sys%operators%op("H"),id(m_env),sb_states) + &
         sp_kron(id(m_sys),env%operators%op("H"),sb_states)  + &
         H2model(sys,env,sb_states)
    call stop_timer("Done H_sb")
    !
    !Diagonalize SuperBLock Hamiltonian
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
    if(allocated(energies))deallocate(energies)
    allocate(energies, source=eig_values)
    !
    !BUILD RHO and PSI:
    call start_timer()
    call rho_sys%free()
    call rho_env%free()
    call psi_sys%free()
    call psi_env%free()
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
       call psi_sys%append(&
            build_PsiMat(Nsys,Nenv,eig_basis(:,1),sb_map,'left'),&
            qn=qn,map=sys%sectors(1)%map(qn=qn))
       !ENVIRONMENT
       qn = current_target_qn-sb_qn
       call rho_env%append(&
            build_density_matrix(Nsys,Nenv,eig_basis(:,1),sb_map,'right'),&
            qn=qn,map=env%sectors(1)%map(qn=qn))
       call psi_env%append(&
            build_PsiMat(Nsys,Nenv,eig_basis(:,1),sb_map,'right'),&
            qn=qn,map=env%sectors(1)%map(qn=qn))
    enddo
    !
    !Build Truncated Density Matrices:
    call rho_sys%eigh(sort=.true.,reverse=.true.)
    call rho_env%eigh(sort=.true.,reverse=.true.)
    !>truncation errors
    rho_sys_evals = rho_sys%evals()
    rho_env_evals = rho_env%evals()
    m_s = min(Mstates,m_sys,size(rho_sys_evals))
    m_e = min(Mstates,m_env,size(rho_env_evals))
    !< ACTHUNG: treatment of degenerate rho-eigenvalues is recommended
    if(Estates/=0d0)then
       !find the value of m such that 1-sum_i lambda_i < Edmrg
       m_err = minloc(abs(1d0-cumulate(rho_sys_evals)-Estates))
       m_s   = m_err(1)
       m_err = minloc(abs(1d0-cumulate(rho_env_evals)-Estates))
       m_e   = m_err(1)
    endif
    truncation_error_sys = 1d0 - sum(rho_sys_evals(1:m_s))
    truncation_error_env = 1d0 - sum(rho_env_evals(1:m_e))
    !
    !>truncation-rotation matrices:
    trRho_sys = rho_sys%sparse(m_sys,m_s)
    trRho_env = rho_env%sparse(m_env,m_e)
    !
    !>Store all the rotation/truncation matrices:
    call sys%put_omat(str(sys%length),trRho_sys)
    call env%put_omat(str(env%length),trRho_env)
    !
    !>Renormalize Blocks:
    call sys%renormalize(as_matrix(trRho_sys))
    call env%renormalize(as_matrix(trRho_env))
    !
    !>Prepare output and update basis state
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
         "Energies/L                           :",energies/current_L
    call stop_timer("dmrg_step")
    call wait(100)
    !
    !Clean memory:
    if(allocated(sb_states))deallocate(sb_states)
    if(allocated(sys_map))deallocate(sys_map)
    if(allocated(env_map))deallocate(env_map)
    if(allocated(sb_map))deallocate(sb_map)
    if(allocated(sb_qn))deallocate(sb_qn)
    if(allocated(qn))deallocate(qn)
    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    ! if(allocated(evals_sys))deallocate(evals_sys)
    ! if(allocated(evals_env))deallocate(evals_env)
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









  subroutine get_sb_states(sb_states,sb_sector)
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



  subroutine write_energy()
    integer                   :: current_L
    integer                   :: Eunit
    current_L = sys%length + env%length
    Eunit     = fopen("energyVSsys.length_"//str(suffix),append=.true.)
    write(Eunit,*)sys%length,energies/current_L
    close(Eunit)
  end subroutine write_energy


  subroutine write_entanglement()
    integer :: Eunit,i
    real(8) :: entropy
    !
    entropy=0d0
    do i=1,size(rho_sys_evals)
       if(rho_sys_evals(i)<0d0)cycle       
       entropy = entropy-rho_sys_evals(i)*log(rho_sys_evals(i))
    enddo
    Eunit     = fopen("SentropyVSsys.length_"//str(suffix),append=.true.)
    write(Eunit,*)sys%length,entropy
    close(Eunit)
  end subroutine write_entanglement

  subroutine write_user(file,vals,x)
    character(len=*) :: file
    real(8),optional :: x
    real(8)          :: x_
    real(8)          :: vals(:)
    integer          :: Eunit
    x_ = sys%length;if(present(x))x_=x
    Eunit     = fopen(str(file)//"_"//str(suffix),append=.true.)
    write(Eunit,*)x_,vals
    close(Eunit)
  end subroutine write_user









  subroutine Measure_Op_dmrg(Op,ivec,file,avOp)
    integer,dimension(:),intent(in)        :: ivec
    type(sparse_matrix),intent(in)         :: Op
    character(len=*)                       :: file
    character(len=1)                       :: label
    real(8),dimension(size(ivec)),optional :: avOp
    real(8)                                :: val
    integer                                :: it,i
    type(sparse_matrix)                    :: Oi
    !
    suffix=label_DMRG('u')
    do it=1,size(ivec)
       i   = ivec(it)
       label='l';if(i>sys%length)label='r'
       Oi  = Build_Op_dmrg(Op,i)
       Oi  = Advance_Op_dmrg(Oi,i)
       val = Average_Op_dmrg(Oi,label)   
       call write_user(trim(file),[val],x=dble(i))
       call Oi%free()
       if(present(avOp))avOp(it)=val
    enddo
  end subroutine Measure_Op_dmrg


  !Return the Operator Oi at a site I of the chain given the operator O in the dot basis:
  function Build_Op_dmrg(Op,i) result(Oi)
    integer                        :: i
    type(sparse_matrix),intent(in) :: Op
    type(sparse_matrix)            :: Oi
    !
    character(len=1)               :: label
    type(sparse_matrix)            :: U
    integer                        :: L,R
    integer                        :: pos
    !
    L = sys%length
    R = env%length
    if(i<1.OR.i>L+R)stop "construct_op_dmrg error: I not in [1,Lchain]"
    !
    label='l';if(i>L)label='r'
    !
    pos=i;if(i>L)pos=R+1-(i-L)
    !
    !Step 1: build the operator at the site pos
    Oi  = Op
    if(pos>1)then
       select case(to_lower(str(label)))
       case ("l")
          U   = sys%omatrices%op(index=pos-1) !get Id of the block_{I-1} 
          Oi  = matmul(U%t(),U).x.Op          !build right-most Op of the block_I as Id(I-1)xOp_I
          U   = sys%omatrices%op(index=pos)   !retrieve O_I 
          Oi  = matmul(U%t(),matmul(Oi,U))    !rotate+truncate Oi at the basis of Block_I 
       case ("r")
          U   = env%omatrices%op(index=pos-1) !get Id of the block_{I-1} 
          Oi  = Op.x.matmul(U%t(),U)          !build right-most Op of the block_I as Id(I-1)xOp_I
          U   = env%omatrices%op(index=pos)   !retrieve O_I 
          Oi  = matmul(U%t(),matmul(Oi,U))    !rotate+truncate Oi at the basis of Block_I 
       end select
       call U%free()
    end if
  end function Build_Op_dmrg


  function Advance_Op_dmrg(Op,i,nstep) result(Oi)
    type(sparse_matrix),intent(in)   :: Op
    integer                          :: i
    integer,optional                 :: nstep
    type(sparse_matrix)              :: Oi,U
    character(len=1)                 :: label
    integer                          :: L,R
    integer                          :: istart,iend,it
    !
    L = sys%length
    R = env%length
    if(i<1.OR.i>L+R)stop "Advance_Op_DMRG error: I not in [1,Lchain]"
    !
    label ='l';if(i>L) label ='r'
    select case(label)
    case ("l")
       istart = i         ; iend   = L ; if(present(nstep))iend=istart+nstep
       if(iend>L)stop "Advance_Op_DMRG ERROR: iend > L"
    case ("r") 
       istart = R+1-(i-L) ; iend   = R ; if(present(nstep))iend=istart+nstep
       if(iend>R)stop "Advance_Op_DMRG ERROR: iend > R"
    end select
    !
    !Bring the operator to the L/R basis of the target states
    !Here we try to reproduce the same strategy of the dmrg_step: increase and rotate 
    Oi = Op
    select case(label)
    case ("l")
       Oi = Oi.x.Id(dot%dim)
       do it=istart+1,iend
          U  = sys%omatrices%op(index=it)
          Oi = matmul(U%t(),matmul(Oi,U))
          Oi = Oi.x.Id(dot%dim)
       enddo
    case ("r")
       Oi = Id(dot%dim).x.Oi
       do it=istart+1,iend
          U  = env%omatrices%op(index=it)
          Oi = matmul(U%t(),matmul(Oi,U))
          Oi = Id(dot%dim).x.Oi
       enddo
    end select
    call U%free()
  end function Advance_Op_dmrg


  function Advance_Corr_dmrg(Op,i,nstep) result(Oi)
    type(sparse_matrix),intent(in)   :: Op
    integer                          :: i
    integer,optional                 :: nstep
    type(sparse_matrix)              :: Oi,U
    character(len=1)                 :: label
    integer                          :: L,R
    integer                          :: istart,iend,it
    !
    L = sys%length
    R = env%length
    if(i<1.OR.i>L+R)stop "Advance_Op_DMRG error: I not in [1,Lchain]"
    !
    label ='l';if(i>L) label ='r'
    select case(label)
    case ("l")
       istart = i         ; iend   = L ; if(present(nstep))iend=istart+nstep
       if(iend>L)stop "Advance_Op_DMRG ERROR: iend > L"
    case ("r") 
       istart = R+1-(i-L) ; iend   = R ; if(present(nstep))iend=istart+nstep
       if(iend>R)stop "Advance_Op_DMRG ERROR: iend > R"
    end select
    !
    !Bring the operator to the L/R basis of the target states
    !Here we try to reproduce the same strategy of the dmrg_step: increase and rotate 
    Oi = Op
    select case(label)
    case ("l")
       do it=istart+1,iend
          U  = sys%omatrices%op(index=it)
          Oi = matmul(U%t(),matmul(Oi,U))
          Oi = Oi.x.Id(dot%dim)
       enddo
    case ("r")
       do it=istart+1,iend
          U  = env%omatrices%op(index=it)
          Oi = matmul(U%t(),matmul(Oi,U))
          Oi = Id(dot%dim).x.Oi
       enddo
    end select
    call U%free()
  end function Advance_Corr_dmrg

  function Average_Op_dmrg(Op,label) result(avOp)
    type(sparse_matrix),intent(in)   :: Op
    character(len=1)                 :: label
    real(8)                          :: avOp
    type(sparse_matrix)              :: Psi
    integer                          :: L,R,i,pos,dims(2)
    real(8),dimension(:),allocatable :: psi_vec
    integer,dimension(:),allocatable :: psi_map
    !
    !Measure using PSI matrix:
    select case(label)
    case ("l")
       if(any(shape(psi_sys)/=shape(Op)))stop "average_op_dmrg ERROR: shape(psi) != shape(Op)"
       Psi  = as_sparse(psi_sys)
       AvOp = trace(as_matrix(matmul(Psi%t(),matmul(Op,Psi))))
    case ("r")
       if(any(shape(psi_env)/=shape(Op)))stop "average_op_dmrg ERROR: shape(psi) != shape(Op)"
       Psi  = as_sparse(psi_env)
       AvOp = trace(as_matrix(matmul(Psi%t(),matmul(Op,Psi))))
    end select
    call Psi%free()
  end function Average_Op_dmrg











  function build_PsiMat(Nsys,Nenv,psi,map,direction) result(psi_mat)
    integer                               :: Nsys,Nenv
    complex(8),dimension(:)               :: psi
    integer,dimension(nsys*nenv)          :: map
    character(len=*)                      :: direction
    complex(8),dimension(:,:),allocatable :: psi_mat
    if(allocated(psi_mat))deallocate(psi_mat)
    select case(to_lower(str(direction)))
    case ('left','l')
       allocate(psi_mat(nsys,nenv));psi_mat=0d0
       psi_mat = transpose(reshape(psi(map), [nenv,nsys]))
    case ('right','r')
       allocate(psi_mat(nenv,nsys));psi_mat=0d0
       psi_mat = reshape(psi(map), [nenv,nsys])
    end select
  end function build_psimat


  function build_density_matrix(Nsys,Nenv,psi,map,direction) result(rho)
    integer                               :: Nsys,Nenv
    complex(8),dimension(:)               :: psi
    integer,dimension(nsys*nenv)          :: map
    character(len=*)                      :: direction
    complex(8),dimension(:,:),allocatable :: rho
    complex(8),dimension(nsys,nenv)       :: psi_tmp
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



  subroutine dmrg_graphic(label)
    integer                  :: label
    integer                  :: i,N,Msys,Menv,LMsys,LMenv,index,Ltot
    character(:),allocatable :: Ldot,Rdot
    real(8)                  :: eps=1d-6
    integer                  :: M=50

    call wait(50)
    !call system("clear")
    call execute_command_line("clear")
    Ltot = Ldmrg/2
    Ldot = bold_green('=')
    Rdot = bold_red('-')
    ! if(Ltot>M)then
    !    Ldot = bg_green('=')
    !    Rdot = bg_red('-')
    ! endif
    !
    N = int(Ltot/(M+eps))+1
    !
    select case(label)
    case default; stop "dmrg_graphic error: label != 1(L),2(R)"
    case(0)
       Msys = int(sys%length/(N+eps))+1
       Menv = int(env%length/(N+eps))+1
       LMsys= Ltot/N-Msys
       LMenv= Ltot/N-Menv
       index=nint(mod(dble(sys%length),N+eps))
       write(LOGfile,*)""
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")"sys; env=",sys%length,env%length,"|"
       if(LMsys>0)write(LOGfile,"("//str(LMsys)//"A)",advance="no")(" ",i=1,LMsys)
       write(LOGfile,"("//str(Msys)//"A)",advance="no")(trim(Ldot),i=1,Msys)
       write(LOGfile,"(A)",advance="no")bold_green("*")//bold("|")//bold_red("*")
       write(LOGfile,"("//str(Menv)//"A)",advance="no")(trim(Rdot),i=1,Menv)
       if(LMenv>0)write(LOGfile,"("//str(LMenv)//"A)",advance="no")(" ",i=1,LMenv)
    case(1)
       Msys = int(sys%length/(N+eps))+1
       Menv = int(env%length/(N+eps))+1
       LMsys= 0
       LMenv= 0
       index=nint(mod(dble(sys%length),N+eps))
       write(LOGfile,*)""
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")"sys; env=",sys%length,env%length,"|"
       if(LMsys>0)write(LOGfile,"("//str(LMsys)//"A)",advance="no")(" ",i=1,LMsys)
       write(LOGfile,"("//str(Msys)//"A)",advance="no")(trim(Ldot),i=1,Msys)
       write(LOGfile,"(A)",advance="no")bg_green(">")//"|"//bold_red("*")
       write(LOGfile,"("//str(Menv)//"A)",advance="no")(trim(Rdot),i=1,Menv)
       if(LMenv>0)write(LOGfile,"("//str(LMenv)//"A)",advance="no")(" ",i=1,LMenv)
    case(2)
       Msys = int(env%length/(N+eps))+1
       Menv = int(sys%length/(N+eps))+1
       LMsys= 0
       LMenv= 0
       index=nint(mod(dble(env%length),N+eps))
       write(LOGfile,*)""
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")"sys; env=",sys%length,env%length,"|"
       if(LMsys>0)write(LOGfile,"("//str(LMsys)//"A)",advance="no")(" ",i=1,LMsys)
       write(LOGfile,"("//str(Msys)//"A)",advance="no")(trim(Ldot),i=1,Msys)
       write(LOGfile,"(A)",advance="no")bold_green("*")//"|"//bg_red("<")
       write(LOGfile,"("//str(Menv)//"A)",advance="no")(trim(Rdot),i=1,Menv)
       if(LMenv>0)write(LOGfile,"("//str(LMenv)//"A)",advance="no")(" ",i=1,LMenv)
    end select
    if(Ltot<=M)then
       write(LOGfile,"(A1)",advance='yes')"|"
    else
       write(LOGfile,"(A1,I3,2x,A,1x,A)",advance='yes')"|",index,Ldot//"->"//str(N)//"= ;",Rdot//"->"//str(N)//"-"
    endif
    call wait(150)
  end subroutine dmrg_graphic



  subroutine sb_HxV(Nloc,v,Hv)
    integer                    :: Nloc
    complex(8),dimension(Nloc) :: v
    complex(8),dimension(Nloc) :: Hv
    complex(8)                 :: val
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

END MODULE SYSTEM





