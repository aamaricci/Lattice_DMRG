program test_iDMRG
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

  character(len=64)                  :: finput
  real(8)                            :: ts
  integer                            :: i,m,current_L
  type(block)                        :: sys,env
  type(sparse_matrix)                :: spHsb
  type(site)                         :: dot
  real(8)                            :: gs_energy,target_Qn(2),current_target_QN(2)
  integer                            :: unit

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"ts",finput,default=-1d0,comment="Hopping amplitude")
  call read_input(finput)

  !
  !Init the single dot structure:
  dot = hubbard_site(xmu=xmu,Uloc=uloc(1))
  !
  !Init block from single dot
  sys=block(dot)
  env=block(dot)

  target_qn=[1d0,1d0]

  !Run DMRG algorithm
  open(free_unit(unit),file="iDMRGqn_energyVSlength_L"//str(Ldmrg)//"_m"//str(Mdmrg)//".dmrg")
  do i=1,Ldmrg
     current_L = 2*sys%length + 2
     write(*,*)"SuperBlock Length  =",current_L
     current_target_QN = int(target_qn*current_L/2d0)
     print*,"Target_QN=",current_target_QN,"filling=",sum(current_target_QN)
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
    type(block)               :: self
    type(site)                :: dot
    character(len=*),optional :: grow
    character(len=16)         :: grow_
    type(block)               :: enl_self
    type(tbasis)              :: self_basis,dot_basis,enl_basis
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    enl_self%length = self%length + 1
    enl_self%Dim    = self%Dim*dot%Dim
    !
    allocate(enl_self%sectors(size(self%sectors)))
    !
    call self%get_basis(self_basis)
    call dot%get_basis(dot_basis)
    !
    select case(str(grow_))
    case ("left","l")
       call enl_self%put("H", &
            (self%operators%op("H").x.id(dot%dim)) +  (id(self%dim).x.dot%operators%op("H")) + &
            H2model(self,as_block(dot)))
       call enl_self%put("Cup", Id(self%dim).x.dot%operators%op("Cup"))
       call enl_self%put("Cdw", Id(self%dim).x.dot%operators%op("Cdw"))
       call enl_self%put("P"  , Id(self%dim).x.dot%operators%op("P"))
       enl_basis = (self_basis.o.dot_basis)
       call enl_self%set_basis( basis=enl_basis )
       !
    case ("right","r")
       call enl_self%put("H", &
            (id(dot%dim).x.self%operators%op("H")) +  (dot%operators%op("H").x.id(self%dim)) + &
            H2model(as_block(dot),self))
       call enl_self%put("Cup", dot%operators%op("Cup").x.Id(self%dim))
       call enl_self%put("Cdw", dot%operators%op("Cdw").x.Id(self%dim))
       call enl_self%put("P"  , dot%operators%op("P").x.Id(self%dim))
       enl_basis = (dot_basis.o.self_basis)
       call enl_self%set_basis( basis=enl_basis )       
       !
    end select
    !
    call self_basis%free()
    call dot_basis%free()
    call enl_basis%free()
    !
  end function enlarge_block



  !H_lr = -t \sum_sigma[ (C^+_{l,sigma}@P_l) x C_{r,sigma}]  + H.c.
  function H2model(left,right,states) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    integer,dimension(:),optional :: states
    type(sparse_matrix)           :: CupL,CdwL
    type(sparse_matrix)           :: CupR,CdwR
    type(sparse_matrix)           :: H2,P
    P    = left%operators%op("P")
    CupL = left%operators%op("Cup")
    CdwL = left%operators%op("Cdw")
    CupR = right%operators%op("Cup")
    CdwR = right%operators%op("Cdw")
    if(present(states))then
       H2   = ts*sp_kron(matmul(CupL%t(),P),CupR,states)  + ts*sp_kron(matmul(CdwL%t(),P),CdwR,states)
       H2   = H2 + H2%dgr()
    else
       H2   = ts*(matmul(CupL%t(),P).x.CupR)  + ts*(matmul(CdwL%t(),P).x.CdwR)
       H2   = H2 + H2%dgr()
    endif
    call P%free
    call CupL%free
    call CdwL%free
    call CupR%free
    call CdwR%free
  end function H2model



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



  subroutine single_dmrg_step(m,sys,env,energy)
    integer                            :: m
    type(block)                        :: sys,env
    type(block)                        :: esys,eenv
    type(sectors_list)                 :: sb_sector
    !
    type(blocks_matrix)                :: rho_sys,rho_env
    type(tbasis)                       :: sys_basis,env_basis

    integer                            :: m_sb,m_s,m_e
    integer                            :: m_esys,m_eenv,Nsys,Nenv
    integer                            :: isys,ienv,isb
    integer                            :: i,j,iqn,Ncv,im

    integer,dimension(:),allocatable   :: sb_states
    integer,dimension(:),allocatable   :: esys_map,eenv_map,sb_map
    real(8),dimension(:),allocatable   :: sb_qn,qn

    real(8),dimension(:),allocatable   :: eig_values
    real(8),dimension(:,:),allocatable :: eig_basis,Hmat
    real(8),dimension(:),allocatable   :: rho_vec

    real(8),dimension(:),allocatable   :: evals,evalsL,evalsR
    real(8),dimension(:,:),allocatable :: Hsb,rhoL,rhoR,truncation_rhoL,truncation_rhoR
    real(8)                            :: truncation_errorL,truncation_errorR,energy

    !Check if blocks are valid ones
    if(.not.sys%is_valid(.true.))stop "single_dmrg_step error: sys is not a valid block"
    if(.not.env%is_valid(.true.))stop "single_dmrg_step error: env is not a valid block"
    !
    !Enlarge blocks
    call start_timer()
    esys = enlarge_block(sys,dot,grow='left')
    eenv = enlarge_block(env,dot,grow='right')
    if(.not.esys%is_valid())stop "single_dmrg_step error: enlarged_sys is not a valid block"
    if(.not.eenv%is_valid())stop "single_dmrg_step error: enlarged_env is not a valid block"
    call stop_timer("Enlarge blocks")
    !
    !Get Enlarged Sys/Env actual dimension
    print*,"Enlarged Blocks dimensions:",esys%dim,eenv%dim    
    m_esys = esys%dim
    m_eenv = eenv%dim
    !
    !Build SuperBLock Sector
    call start_timer()
    call get_sb_states(esys,eenv,sb_states,sb_sector)
    call stop_timer("Get SB states")
    m_sb = size(sb_states)


    print*,"Super Block dimensions:", m_esys*m_eenv, m_sb
    call start_timer()
    spHsb = sp_kron(esys%operators%op("H"),id(m_eenv),sb_states) + &
         sp_kron(id(m_esys),eenv%operators%op("H"),sb_states)  + &
         H2model(esys,eenv,sb_states)
    call stop_timer("Build H_sb")



    !
    if(allocated(eig_values))deallocate(eig_values)
    if(allocated(eig_basis))deallocate(eig_basis)
    allocate(eig_values(Lanc_Neigen))    ;eig_values=0d0
    allocate(eig_basis(m_sb,Lanc_Neigen));eig_basis =0d0
    call start_timer()
    if(m_sb < lanc_dim_threshold)then
       print*,"diag SB w/ Lapack"
       allocate(Hsb(m_sb,m_sb));Hsb=0d0
       allocate(evals(m_sb))
       call spHsb%dump(Hsb)
       call eigh(Hsb,evals)
       eig_basis(:,1:Lanc_Neigen) = Hsb(:,1:Lanc_Neigen)
       eig_values(1:Lanc_Neigen)  = evals(1:Lanc_Neigen)
       deallocate(Hsb,evals)
    else
       print*,"diag SB w/ Arpack"
       Ncv = min(lanc_ncv_factor*Lanc_Neigen+lanc_ncv_add,m_sb)
       call sp_eigh(sb_HxV,eig_values,eig_basis,&
            Ncv,&
            lanc_niter,&
            tol=lanc_tolerance,&
            iverbose=.false.)
    end if
    energy = eig_values(1)
    call stop_timer("Diag H_sb")



    !BUILD RHO:
    call start_timer()
    do isb=1,size(sb_sector)
       sb_qn  = sb_sector%qn(index=isb)
       sb_map = sb_sector%map(index=isb)
       Nsys   = size(esys%sectors(1)%map(qn=sb_qn))
       Nenv   = size(eenv%sectors(1)%map(qn=(current_target_qn - sb_qn)))
       if(Nsys*Nenv==0)cycle
       !LEFT :
       qn = sb_qn
       call rho_sys%append(&
            build_density_matrix(Nsys,Nenv,eig_basis(:,1),sb_map,'left'),&
            qn=qn,map=esys%sectors(1)%map(qn=qn))
       !RIGHT:
       qn = current_target_qn-sb_qn
       call rho_env%append(&
            build_density_matrix(Nsys,Nenv,eig_basis(:,1),sb_map,'right'),&
            qn=qn,map=eenv%sectors(1)%map(qn=qn))
    enddo
    !
    !SYSTEM:
    m_s = min(Mdmrg,m_esys)
    m_e = min(Mdmrg,m_eenv)
    allocate(truncation_rhoL(m_esys,m_s));truncation_rhoL=0d0
    allocate(truncation_rhoR(m_eenv,m_e));truncation_rhoR=0d0
    call rho_sys%eigh(sort=.true.,reverse=.true.)
    call rho_env%eigh(sort=.true.,reverse=.true.)
    do im=1,m_s
       rho_vec  = rho_sys%evec(m=im)
       esys_map = rho_sys%map(m=im)
       truncation_rhoL(esys_map,im) = rho_vec
       call sys_basis%append( qn=rho_sys%qn(m=im) )
    enddo
    do im=1,m_e
       rho_vec  = rho_env%evec(m=im)
       eenv_map = rho_env%map(m=im)
       truncation_rhoR(eenv_map,im) = rho_vec
       call env_basis%append( qn=rho_env%qn(m=im) )
    enddo
    evalsL = rho_sys%evals()
    truncation_errorL = 1d0 - sum(evalsL(1:m_s))
    evalsR = rho_env%evals()
    truncation_errorR = 1d0 - sum(evalsR(1:m_e))
    call stop_timer("Get Rho + U")
    write(*,"(A,L12,12X,L12)")"Truncating      :",m<=m_esys,m<=m_eenv
    write(*,"(A,I12,12X,I12)")"Truncation dim  :",m_s,m_e
    write(*,"(A,2ES24.15)")"Truncation error:",truncation_errorL,truncation_errorR
    write(*,"(A,"//str(Lanc_Neigen)//"F24.15)")"Energies        :",eig_values
    write(*,*)""


    !Renormalize Blocks:
    call start_timer()
    call esys%renormalize(truncation_rhoL)
    call eenv%renormalize(truncation_rhoR)
    call stop_timer("Renormalize Blocks")
    !
    !Prepare output and update basis state
    call start_timer()
    sys     = block(esys)
    env     = block(eenv)
    call sys%set_basis(basis=sys_basis)
    call env%set_basis(basis=env_basis)
    call stop_timer("Set Basis Output")

    !Clean memory:
    if(allocated(esys_map))deallocate(esys_map)
    if(allocated(eenv_map))deallocate(eenv_map)
    if(allocated(sb_map))deallocate(sb_map)
    if(allocated(sb_states))deallocate(sb_states)
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
    call sys_basis%free()
    call env_basis%free()
    call sb_sector%free()
    call rho_sys%free()
    call rho_env%free()
    !
  end subroutine single_dmrg_step



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

end program test_iDMRG
