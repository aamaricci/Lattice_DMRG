program test_iDMRG
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
  real(8)                            :: gs_energy,target_Qn(2),current_target_QN(2)
  integer                            :: unit

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"ts",finput,default=-1d0,comment="Hopping amplitude")
  call parse_input_variable(uloc,"uloc",finput,default=0d0,comment="Hubbard U")
  call parse_input_variable(Lmax,"LMAX",finput,default=1,comment="Final chain length ")
  call parse_input_variable(m,"M",finput,default=20,&
       comment="Number of states retained at truncation of \rho")
  call parse_input_variable(Neigen,"NEIGEN",finput,default=2,&
       comment="Number of eigenpairs required to SB Hamiltonian")
  call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",finput,default=10,&
       comment="Arpack Block size parameters")
  call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",finput,default=0,&
       comment="Arpack Block size parameters")
  call parse_input_variable(lanc_niter,"LANC_NITER",finput,default=512,&
       comment="Number of Lanczos iteration in spectrum determination.")
  call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",finput,default=1d-12,&
       comment="Lanczos tolerance ")
  call parse_input_variable(lanc_threshold,"LANC_THRESHOLD",finput,default=512,&
       comment="Lapack threshold for Arpack ")





  !
  !Init the single dot structure:
  dot = hubbard_site(xmu=0d0,Uloc=uloc)
  !
  !Init block from single dot
  sys=block(dot)
  env=sys

  target_qn=[1d0,1d0]
  !Run DMRG algorithm
  open(free_unit(unit),file="energyVSlength_L"//str(Lmax)//"_m"//str(m)//".dmrg")
  do i=1,Lmax
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
    type(block),intent(inout)        :: self
    type(site)                       :: dot
    character(len=*),optional        :: grow
    character(len=16)                :: grow_
    type(block)                      :: enl_self
    real(8),dimension(:),allocatable :: self_basis,dot_basis
    integer :: iqn
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    enl_self%length = self%length + 1
    enl_self%Dim    = self%Dim*dot%Dim
    enl_self%DimUp  = self%DimUp*dot%DimUp
    enl_self%DimDw  = self%DimDw*dot%DimDw   !
    !
    select case(str(grow_))
    case ("left","l")
       call enl_self%put("H", (self%operators%op("H").x.id(dot%dim)) + (id(self%dim).x.dot%operators%op("H")) + &
            H2model(self,as_block(dot)))
       call enl_self%put("Cup", Id(self%dim).x.dot%operators%op("Cup"))
       call enl_self%put("Cdw", Id(self%dim).x.dot%operators%op("Cdw"))
       call enl_self%put("P"  , Id(self%dim).x.dot%operators%op("P"))
    case ("right","r")
       call enl_self%put("H", &
            (id(dot%dim).x.self%operators%op("H")) +  (dot%operators%op("H").x.id(self%dim)) + &
            H2model(as_block(dot),self))
       call enl_self%put("Cup", dot%operators%op("Cup").x.Id(self%dim))
       call enl_self%put("Cdw", dot%operators%op("Cdw").x.Id(self%dim))
       call enl_self%put("P"  , dot%operators%op("P").x.Id(self%dim))
    end select
    !
    allocate(enl_self%sectors(size(self%sectors)))
    do iqn=1,size(self%sectors)       
       self_basis = self%sectors(iqn)%basis()
       dot_basis  = dot%sectors(iqn)%basis()
       print*,outsum(self_basis,dot_basis)
       call enl_self%set_sectors( indx=iqn, vec=outsum(self_basis,dot_basis) )       
       deallocate(self_basis,dot_basis)
    enddo
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



  subroutine single_dmrg_step(m,sys,env,energy)
    integer                            :: m
    type(block)                        :: sys,env
    type(block)                        :: esys,eenv
    type(sectors_list)                 :: sb_sys_sector(2)
    type(sectors_list)                 :: sb_env_sector(2)
    integer                            :: Dims(2)
    real(8)                            :: esys_qn,eenv_qn,sb_qn
    real(8)                            :: esys_qn_dw,eenv_qn_dw
    real(8)                            :: esys_qn_up,eenv_qn_up
    integer,dimension(:),allocatable   :: esys_map,eenv_map,sb_map,states
    integer,dimension(:),allocatable   :: esys_map_dw,eenv_map_dw
    integer,dimension(:),allocatable   :: esys_map_up,eenv_map_up
    integer,dimension(:),allocatable   :: sb_states_up,sb_states_dw
    integer,dimension(:),allocatable   :: sb_states
    integer                            :: m_sb
    integer                            :: m_esys,m_eenv,Nsys,Nenv
    integer                            :: isys,isys_up,isys_dw
    integer                            :: ienv,ienv_up,ienv_dw
    integer                            :: iup,idw,jup,jdw
    integer                            :: i,j,iqn,Ncv,im
    real(8),dimension(:),allocatable   :: eig_values
    real(8),dimension(:,:),allocatable :: eig_basis,Hmat
    real(8),dimension(:),allocatable   :: rho_vec
    type(blocks_matrix)                :: rho_sector
    real(8),dimension(:),allocatable   :: evals,evalsL,evalsR
    real(8),dimension(:,:),allocatable :: Hsb,rhoL,rhoR,truncation_rhoL,truncation_rhoR
    real(8)                            :: truncation_errorL,truncation_errorR,energy

    !Check if blocks are valid ones
    if(.not.sys%is_valid())stop "single_dmrg_step error: sys is not a valid block"
    if(.not.env%is_valid())stop "single_dmrg_step error: env is not a valid block"
    !
    !Enlarge blocks
    print*,"Enlarge blocks:"
    esys = enlarge_block(sys,dot,grow='left')
    if(.not.esys%is_valid())stop "single_dmrg_step error: enlarged_sys is not a valid block"
    eenv = enlarge_block(env,dot,grow='right')
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
    ! !SYSTEM:
    ! !DW electrons
    ! iqn=1
    ! do isys=1,size(esys%sectors(1))
    !    esys_qn = esys%sectors(1)%qn(index=isys)
    !    eenv_qn = current_target_qn(1) - esys_qn
    !    if(.not.eenv%sectors(1)%has_qn(eenv_qn))cycle
    !    esys_map= esys%sectors(1)%map(qn=esys_qn)
    !    eenv_map= eenv%sectors(1)%map(qn=eenv_qn)
    !    do i=1,size(esys_map)
    !       do j=1,size(eenv_map)
    !          call append(sb_states_dw, eenv_map(j)-1 + (esys_map(i)-1)*eenv%DimDw)
    !          call sb_sys_sector(1)%append(qn=esys_qn,istate=size(sb_states_dw))
    !       enddo
    !    enddo
    ! enddo
    ! !UP electrons
    ! do isys=1,size(esys%sectors(2))
    !    esys_qn = esys%sectors(2)%qn(index=isys)
    !    eenv_qn = current_target_qn(2) - esys_qn
    !    if(.not.eenv%sectors(2)%has_qn(eenv_qn))cycle
    !    esys_map= esys%sectors(2)%map(qn=esys_qn)
    !    eenv_map= eenv%sectors(2)%map(qn=eenv_qn)
    !    do i=1,size(esys_map)
    !       do j=1,size(eenv_map)
    !          call append(sb_states_up, eenv_map(j)-1 + (esys_map(i)-1)*eenv%DimUp)
    !          call sb_sys_sector(2)%append(qn=esys_qn,istate=size(sb_states_up))
    !       enddo
    !    enddo
    ! enddo
    ! !
    ! m_sb = size(sb_states_up)*size(sb_states_dw)
    ! allocate(sb_states(m_sb))
    ! do idw=1,size(sb_states_dw)
    !    do iup=1,size(sb_states_up)
    !       i = iup+(idw-1)*size(sb_states_up)
    !       sb_states(i) = sb_states_up(iup)+sb_states_up(idw)*esys%Dim
    !       print*,i,sb_states_up(iup),sb_states_dw(idw),sb_states(i)
    !    enddo
    ! enddo
    ! call stop_timer()
    ! deallocate(sb_states_up,sb_states_dw)
    ! print*,sb_states


    do isys_up=1,size(esys%sectors(1))
       esys_qn_up = esys%sectors(1)%qn(index=isys_up)
       eenv_qn_up = current_target_qn(1) - esys_qn_up
       if(.not.eenv%sectors(1)%has_qn(eenv_qn_up))cycle
       do isys_dw=1,size(esys%sectors(2))
          esys_qn_dw = esys%sectors(2)%qn(index=isys_dw)
          eenv_qn_dw = current_target_qn(2) - esys_qn_dw
          if(.not.eenv%sectors(2)%has_qn(eenv_qn_dw))cycle
          !
          !
          esys_map_up= esys%sectors(1)%map(qn=esys_qn_up)
          esys_map_dw= esys%sectors(2)%map(qn=esys_qn_dw)
          eenv_map_up= eenv%sectors(1)%map(qn=eenv_qn_up)
          eenv_map_dw= eenv%sectors(2)%map(qn=eenv_qn_dw)
          !
          print*,esys_qn_up,esys_qn_dw,"|",eenv_qn_up,eenv_qn_dw
          do idw=1,size(esys_map_dw)
             do iup=1,size(esys_map_up)
                isys = esys_map_up(iup) + (esys_map_dw(idw)-1)*esys%DimUp
                !
                do jdw=1,size(eenv_map_dw)
                   do jup=1,size(eenv_map_up)
                      ienv = eenv_map_up(jup) + (eenv_map_dw(jdw)-1)*eenv%DimUp
                      i=ienv + (isys-1)*eenv%Dim
                      print*,isys,ienv,i
                      call append(sb_states, i)
                      ! call sb_sys_sector%append(qn=esys_qn,istate=size(sb_states))
                   enddo
                enddo
             enddo
          enddo
          print*,""
       enddo
    enddo

    m_sb = size(sb_states)
    print*,sb_states
    ! stop

    !m_sb=m_esys*m_eenv
    print*,"Super Block dimensions:", m_esys*m_eenv, m_sb
    !BUild SuperBlock matrix:
    !The two are indeed equivalent: so the error is not here!!! WTFF
    spH =  (esys%operators%op("H").x.id(m_eenv)) + (id(m_esys).x.eenv%operators%op("H")) + H2model(esys,eenv)
    allocate(Hmat(spH%Nrow,spH%Ncol))
    call spH%dump(Hmat)
    call spHsb%load(Hmat(sb_states,sb_states))
    call spHsb%show(file="Hsb_1.dat")
    print*,shape(spHsb)
    ! call spHsb%free
    ! spHsb = sp_kron(esys%operators%op("H"),id(m_eenv),sb_states) + &
    !      sp_kron(id(m_esys),eenv%operators%op("H"),sb_states)  + &
    !      H2model(esys,eenv,sb_states)
    call stop_timer("Build H_sb")
    ! call spHsb%display()

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
    print*,eig_values
    stop    
    ! !
    ! call start_timer()
    ! !Left
    ! rhoL = build_density_matrix(esys,eenv,eig_basis(:,1),'left')
    ! allocate(evalsL(m_esys))
    ! call eigh(rhoL,evalsL)
    ! evalsL    = evalsL(m_esys:1:-1)
    ! rhoL(:,:) = rhoL(:,m_esys:1:-1)
    ! !
    ! !Right
    ! rhoR = build_density_matrix(esys,eenv,eig_basis(:,1),'right')
    ! allocate(evalsR(m_eenv))
    ! call eigh(rhoR,evalsR)
    ! evalsR    = evalsR(m_eenv:1:-1)
    ! rhoR(:,:) = rhoR(:,m_eenv:1:-1)
    ! call stop_timer("Get Rho_Left-dot / Rho_dot-Right")
    ! !
    ! !
    ! call start_timer()
    ! m_s = min(m,m_esys);
    ! allocate(truncation_rhoL(m_esys,m_s))
    ! truncation_rhoL(:,1:m_s) = rhoL(:,1:m_s)
    ! truncation_errorL = 1d0 - sum(evalsL(1:m_s))
    ! !
    ! m_e = min(m,m_eenv);
    ! allocate(truncation_rhoR(m_eenv,m_e))
    ! truncation_rhoR(:,1:m_e) = rhoR(:,1:m_e)
    ! truncation_errorR = 1d0 - sum(evalsR(1:m_e))
    ! call stop_timer("Build U")
    ! write(*,*)"Truncating      :",m<=m_esys,m<=m_eenv
    ! write(*,*)"Truncation dim  :",m_s,m_e
    ! write(*,*)"Truncation error:",truncation_errorL,truncation_errorR
    ! !
    ! !find a clever way to iterate here:
    ! !Return updated Block:
    ! call start_timer()
    ! sys = esys
    ! sys%dim    = m_s
    ! call sys%put("H",rotate_and_truncate(esys%operators%op("H"),truncation_rhoL,m_esys,m_s))
    ! call sys%put("Cup",rotate_and_truncate(esys%operators%op("Cup"),truncation_rhoL,m_esys,m_s))
    ! call sys%put("Cdw",rotate_and_truncate(esys%operators%op("Cdw"),truncation_rhoL,m_esys,m_s))
    ! call sys%put("P",rotate_and_truncate(esys%operators%op("P"),truncation_rhoL,m_esys,m_s))
    ! do i=1,size(sys%sectors)
    !    call sys%set_sectors(indx=i,vec=sys%sectors(i)%basis())
    ! enddo


    ! env      = esys
    ! env%dim  = m_e
    ! call env%put("H",rotate_and_truncate(eenv%operators%op("H"),truncation_rhoR,m_eenv,m_e))
    ! call env%put("Cup",rotate_and_truncate(eenv%operators%op("Cup"),truncation_rhoR,m_eenv,m_e))
    ! call env%put("Cdw",rotate_and_truncate(eenv%operators%op("Cdw"),truncation_rhoR,m_eenv,m_e))
    ! call env%put("P",rotate_and_truncate(eenv%operators%op("P"),truncation_rhoR,m_eenv,m_e))
    ! do i=1,size(env%sectors)
    !    call env%set_sectors(indx=i,vec=env%sectors(i)%basis())
    ! enddo
    ! call stop_timer("Update+Truncate Block")
    ! !
    if(allocated(esys_map))deallocate(esys_map)
    if(allocated(eenv_map))deallocate(eenv_map)
    if(allocated(sb_map))deallocate(sb_map)

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
    call sb_sys_sector%free()
    call sb_env_sector%free()
    call rho_sector%free()
  end subroutine single_dmrg_step



  function build_density_matrix(sys,env,psi,direction) result(rho)
    type(block)                        :: sys,env
    real(8),dimension(sys%Dim*env%Dim) :: psi
    character(len=*)                   :: direction
    real(8),dimension(:,:),allocatable :: rho
    real(8),dimension(sys%Dim,env%Dim) :: rho_tmp
    !
    if(allocated(rho))deallocate(rho)
    !
    rho_tmp = transpose(reshape(psi, [env%Dim,sys%Dim]))
    !
    select case(to_lower(str(direction)))
    case ('left','l')
       allocate(rho(sys%Dim,sys%Dim));rho=0d0
       rho  = matmul(rho_tmp,  transpose(rho_tmp) )
    case ('right','r')
       allocate(rho(env%Dim,env%Dim));rho=0d0
       rho  = matmul(transpose(rho_tmp), rho_tmp  )
    end select
    !
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

end program test_iDMRG
