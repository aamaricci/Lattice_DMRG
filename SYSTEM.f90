MODULE SYSTEM
  USE GLOBAL
  implicit none
  private

  public :: init_DMRG
  public :: finalize_DMRG
  public :: step_DMRG
  !
  public :: infinite_DMRG
  public :: finite_DMRG

contains



  !##################################################################
  !              INIT/FINALIZE DMRG ALGORITHM
  !##################################################################
  subroutine init_dmrg(h2,ModelDot)
    procedure(Hconnect) :: h2
    type(site)          :: ModelDot
    if(associated(H2model))nullify(H2model)
    H2model => h2
    allocate(target_qn, source=DMRG_QN)
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



  !##################################################################
  !              FINITE/INFINITE DMRG ALGORITHM
  !##################################################################
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
       call write_truncation()
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







  !##################################################################
  !              WORKHORSE: STEP DMRG 
  !##################################################################
  subroutine step_dmrg(label,isweep)
    integer,optional                      :: label,isweep
    integer                               :: iLabel
    integer                               :: Mstates,j_
    real(8)                               :: Estates,e_,err
    integer                               :: m_sb,m_s,m_e
    integer                               :: m_sys,m_env,Nsys,Nenv
    integer                               :: isys,ienv,isb,m_err(1),m_threshold
    integer                               :: i,j,r,l,iqn,Ncv,im,unit,current_L,istate
    integer,dimension(:),allocatable      :: sys_map,env_map,sb_map
    real(8),dimension(:),allocatable      :: sb_qn,qn
    real(8),dimension(:),allocatable      :: evals
    real(8),dimension(:,:),allocatable :: Hsb
    real(8),dimension(:),allocatable      :: eig_values
    real(8),dimension(:,:),allocatable :: eig_basis
    integer,dimension(:),allocatable      :: sb_states
    type(tbasis)                          :: sys_basis,env_basis
    type(sparse_matrix)                   :: trRho_sys,trRho_env
    type(sectors_list)                    :: sb_sector
    logical                               :: exist
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
    current_target_QN = int(target_qn*current_L*Norb)
    write(LOGfile,"(A20,I12)")&
         "SuperBlock Length = ",current_L
    write(LOGfile,"(A20,"//str(size(current_target_QN))//"G12.7)")&
         "Target_QN         = ",current_target_QN
    write(LOGfile,"(A20,G12.7)")&
         "Total Occupation  = ",sum(current_target_QN)
    write(LOGfile,"(A20,2G12.7)")&
         "Filling (1/Norb)  = ",sum(current_target_QN)/current_L,sum(current_target_QN)/current_L/Norb
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
    call start_timer("Get SB states")
    call get_sb_states(sb_states,sb_sector)
    call stop_timer("")
    m_sb = size(sb_states)
    !    
    write(LOGfile,"(A,I12,I12))")&
         "Enlarged Blocks dimensions           :", sys%dim,env%dim  
    write(LOGfile,"(A,I12,A1,I12,A1,F10.5,A1)")&
         "SuperBlock Dimension  (tot)          :", m_sb,"(",m_sys*m_env,")",100*dble(m_sb)/m_sys/m_env,"%"
    !
    !Build SuperBLock Hamiltonian
    call start_timer("get H_sb")
    spHsb = sp_kron(sys%operators%op("H"),id(m_env),sb_states) + &
         sp_kron(id(m_sys),env%operators%op("H"),sb_states)    + &
         H2model(sys,env,sb_states)
    call stop_timer("Done H_sb")
    !
    !Diagonalize SuperBLock Hamiltonian
    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    allocate(eig_values(Lanc_Neigen))    ;eig_values=0d0
    allocate(eig_basis(m_sb,Lanc_Neigen));eig_basis =zero
    call start_timer()
    if(m_sb < lanc_dim_threshold)then
       allocate(Hsb(m_sb,m_sb));Hsb=zero
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
       ! print*,sb_qn,size(sb_map)
       ! print*,sb_map
       ! print*,"=="
       ! print*,sb_states(sb_map)
       ! do i=1,size(sb_map)
       !    istate = sb_states(sb_map(i))
       !    l = (istate-1)/env%Dim+1
       !    r = mod(istate,env%Dim);if(r==0)r=env%Dim
       !    print*,l,r,istate
       ! enddo
       ! print*,""
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


    !Build Truncated Density Matrices:
    call rho_sys%eigh(sort=.true.,reverse=.true.)
    call rho_env%eigh(sort=.true.,reverse=.true.)    
    rho_sys_evals = rho_sys%evals()
    rho_env_evals = rho_env%evals()
    if(Mstates/=0)then
       m_s = min(Mstates,m_sys,size(rho_sys_evals))
       m_e = min(Mstates,m_env,size(rho_env_evals))       
    elseif(Estates/=0d0)then
       m_err = minloc(abs(1d0-cumulate(rho_sys_evals)-Estates))
       m_s   = m_err(1)
       m_err = minloc(abs(1d0-cumulate(rho_env_evals)-Estates))
       m_e   = m_err(1)
    else
       stop "Mdmrg=Edmrg=0 can not fix a threshold for the RDM"
    endif
    e_=rho_sys_evals(m_s)
    j_=m_s
    do i=j_+1,size(rho_sys_evals)
       err=abs(rho_sys_evals(i)-e_)/e_
       if(err<=1d-1)m_s=m_s+1
    enddo
    e_=rho_env_evals(m_e)
    j_=m_e
    do i=j_+1,size(rho_env_evals)
       err=abs(rho_env_evals(i)-e_)/e_
       if(err<=1d-1)m_e=m_e+1
    enddo
    truncation_error_sys = 1d0 - sum(rho_sys_evals(1:m_s))
    truncation_error_env = 1d0 - sum(rho_env_evals(1:m_e))
    !
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
         "Energies/L                           :",energies/sum(current_target_QN)
    call stop_timer("dmrg_step")
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



  !-----------------------------------------------------------------!
  ! Purpose: enlarge a given BLOCK "SELF" growing it left/right with
  ! a SITE "DOT" (specified in the init)
  !-----------------------------------------------------------------!
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


  !-----------------------------------------------------------------!
  ! Purpose: build the list of states compatible with the specified
  ! quantum numbers
  !-----------------------------------------------------------------!
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
    call sb_sector%free()
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






  !##################################################################
  !              BUILD REDUCED DENSITY MATRIX 
  !##################################################################
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
       allocate(rho(nsys,nsys));rho=zero
       rho  = matmul(psi_tmp,  (transpose(psi_tmp)) )
    case ('right','r')
       allocate(rho(nenv,nenv));rho=zero
       rho  = matmul((transpose(psi_tmp)), psi_tmp  )
    end select
  end function build_density_matrix




  !##################################################################
  !              RESHAPE GROUND STATE AS MATRIX 
  !##################################################################
  function build_PsiMat(Nsys,Nenv,psi,map,direction) result(psi_mat)
    integer                            :: Nsys,Nenv
    real(8),dimension(:)               :: psi
    integer,dimension(nsys*nenv)       :: map
    character(len=*)                   :: direction
    real(8),dimension(:,:),allocatable :: psi_mat
    if(allocated(psi_mat))deallocate(psi_mat)
    select case(to_lower(str(direction)))
    case ('left','l')
       allocate(psi_mat(nsys,nenv));psi_mat=zero
       psi_mat = transpose(reshape(psi(map), [nenv,nsys]))
    case ('right','r')
       allocate(psi_mat(nenv,nsys));psi_mat=zero
       psi_mat = reshape(psi(map), [nenv,nsys])
    end select
  end function build_psimat





  !##################################################################
  !              SuperBlock MATRIX-VECTOR PRODUCT 
  !##################################################################
  subroutine sb_HxV(Nloc,v,Hv)
    integer                 :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                 :: i,j,jcol
    Hv=zero
    do i=1,Nloc
       matmul: do jcol=1, spHsb%row(i)%Size
          val = spHsb%row(i)%vals(jcol)
          j   = spHsb%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end subroutine sb_HxV





  !##################################################################
  !                 WRITE TO FILE
  !##################################################################
  !-----------------------------------------------------------------!
  ! Purpose: write energy to file
  !-----------------------------------------------------------------!
  subroutine write_energy()
    integer                   :: current_L
    integer                   :: Eunit
    current_L = sys%length + env%length
    Eunit     = fopen("energyVSsys.length_"//str(suffix),append=.true.)
    write(Eunit,*)sys%length,energies/current_L/Norb
    close(Eunit)
  end subroutine write_energy


  !-----------------------------------------------------------------!
  ! Purpose: write entanglement entropy to file
  !-----------------------------------------------------------------!
  subroutine write_truncation()
    integer                   :: current_L
    integer                   :: Eunit
    current_L = sys%length + env%length
    Eunit     = fopen("truncationVSsys.length_"//str(suffix),append=.true.)
    write(Eunit,*)sys%length,truncation_error_sys/current_L/Norb,truncation_error_env/current_L/Norb
    close(Eunit)
  end subroutine write_truncation


  !-----------------------------------------------------------------!
  ! Purpose: write entanglement entropy to file
  !-----------------------------------------------------------------!
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

END MODULE SYSTEM





