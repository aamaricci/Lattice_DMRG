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
  dot = hubbard_site_ud(xmu=0d0,Uloc=uloc)
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






  !We want to arrange states such that the basis is {_n .... _1}_dw {_n ... _1}_up
  ! H = h_d + 1_dw x h_up + h_dw x 1_up
  ! h_d  = H_d(block) x 1(dot) + 1(block) x H_d(dot) : contains the local parts
  ! h_up = h_up(block) x 1(dot) + H2(up)
  function enlarge_block(self,dot,grow) result(enl_self)
    type(block),intent(inout)        :: self
    type(site),intent(inout)         :: dot
    character(len=*),optional        :: grow
    character(len=16)                :: grow_
    type(block)                      :: enl_self
    integer                          :: DsBlock,DsDot
    real(8),dimension(:),allocatable :: self_basis, dot_basis
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    DsBlock = 2**self%length
    DsDot   = 2
    !
    enl_self%length = self%length + 1
    enl_self%DimUp  = self%DimUp*dot%DimUp
    enl_self%DimDw  = self%DimDw*dot%DimDw
    enl_self%dim    = self%dim*dot%dim
    !
    select case(str(grow_))
    case ("left","l")
       call enl_self%put("Hd",  (self%operators%op("Hd") .x.Id(dot%dim)) + (Id(self%dim).x.dot%operators%op("Hd")))
       call enl_self%put("Hup", (self%operators%op("Hup").x.Id(DsDot))   +  H2model(self,as_block(dot)))
       call enl_self%put("Hdw", (self%operators%op("Hdw").x.Id(DsDot))   +  H2model(self,as_block(dot)))
       call enl_self%put("Cp",  Id(DsBlock).x.dot%operators%op("Cp"))
       call enl_self%put("P"  , Id(DsBlock).x.dot%operators%op("P"))
    case ("right","r")
       call enl_self%put("Hd",  (Id(dot%dim).x.self%operators%op("Hd") ) + (dot%operators%op("Hd").x.Id(self%dim)) )
       call enl_self%put("Hup", (Id(DsDot).x.self%operators%op("Hup"))   +  H2model(as_block(dot),self) )
       call enl_self%put("Hdw", (Id(DsDot).x.self%operators%op("Hdw"))   +  H2model(as_block(dot),self) )
       call enl_self%put("Cp",  dot%operators%op("Cp").x.Id(DsBlock))
       call enl_self%put("P"  , dot%operators%op("P").x.Id(DsBlock))
    end select
    !
    allocate(enl_self%sectors(2))
    self_basis = self%sectors(1)%basis()
    dot_basis  = dot%sectors(1)%basis()
    call enl_self%set_sectors( indx=1, vec=outsum(self_basis,dot_basis) )
    self_basis = self%sectors(2)%basis()
    dot_basis  = dot%sectors(2)%basis()
    call enl_self%set_sectors( indx=2, vec=outsum(self_basis,dot_basis) )
    deallocate(self_basis,dot_basis)
  end function enlarge_block



  !H_lr = -t (C^+_{l}@P_l) x C_{r}  + H.c.
  function H2model(left,right) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    type(sparse_matrix)           :: CpL
    type(sparse_matrix)           :: CpR
    type(sparse_matrix)           :: H2,P
    P   = left%operators%op("P") 
    CpL = left%operators%op("Cp")
    CpR = right%operators%op("Cp")
    H2  = ts*(matmul(CpL%t(),P).x.CpR)
    H2  = H2 + H2%dgr()
    call CpL%free
    call CpR%free
    call P%free
  end function H2model


  function build_Hblock(self) result(spH)
    type(block),intent(in) :: self
    type(sparse_matrix)    :: spH
    spH = (self%operators%op("Hd")) + (Id(self%DimDw).x.self%operators%op("Hup")) + (self%operators%op("Hdw").x.Id(self%DimUp))
  end function build_Hblock


  function build_Hsb(left,right) result(spH)
    type(block),intent(in) :: left,right
    type(sparse_matrix)    :: spH,spL,spR,Hup,Hdw
    integer                :: DimDw,DimUp
    DimDw  = left%DimDw*right%DimDw
    DimUp  = left%DimUp*right%DimUp
    Hup = (left%operators%op("Hup").x.Id(right%DimUp)) + (Id(left%DimUp).x.right%operators%op("Hup")) + H2model(left,right)
    Hdw = (left%operators%op("Hdw").x.Id(right%DimDw)) + (Id(left%DimDw).x.right%operators%op("Hdw")) + H2model(left,right)
    spH = (left%operators%op("Hd").x.Id(right%dim)) + (Id(left%dim).x.right%operators%op("Hd")) + (Id(DimDw).x.Hup) + (Hdw.x.Id(DimUp))
  end function build_Hsb




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

    m_sb = m_esys*m_eenv
    print*,"Super Block dimensions:", m_sb
    !BUild SuperBlock matrix:
    call start_timer()
    spHsb  = build_Hsb(esys,eenv)
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
    ! allocate(rhoL(m_esys,m_esys));rhoL=0d0
    ! allocate(rhoR(m_eenv,m_eenv));rhoR=0d0
    ! rhoL = build_density_matrix(m_esys,m_eenv,eig_basis(:,1),'left')
    ! rhoR = build_density_matrix(m_esys,m_eenv,eig_basis(:,1),'right')
    rhoL = build_density_matrix(esys,eenv,eig_basis(:,1),'left')
    rhoR = build_density_matrix(esys,eenv,eig_basis(:,1),'right')
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
    call sys%put("Hd",rotate_and_truncate(esys%operators%op("Hd"),truncation_rhoL,m_esys,m_))
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



  function build_density_matrix(sys,env,psi,direction) result(rho)
    type(block),intent(in)                 :: sys,env
    integer                                :: DsysUp,DsysDw
    integer                                :: DenvUp,DenvDw
    real(8),dimension(sys%Dim*env%Dim)     :: psi
    character(len=*)                       :: direction
    real(8),dimension(:,:),allocatable     :: rho  ![DsysUp*DsysDw,DenvUp*DenvDw] < output
    real(8),dimension(:,:),allocatable     :: rho2 ![DsysUp*DenvUp,DsysDw*DenvDw]
    real(8),dimension(:,:,:,:),allocatable :: rho4 ![DsysUp,DenvUp,DsysDw,DenvDw]
    integer                                :: iup,idw,i2(2),j2(2)
    integer                                :: s,sp,e,ep
    integer                                :: sup,sdw,spup,spdw
    integer                                :: eup,edw,epup,epdw
    ! real(8),dimension(sys%DimUp*env%DimUp,sys%DimDw*env%DimDw) :: rho2
    ! real(8),dimension(sys%DimUp,env%DimUp,sys%DimDw,env%DimDw) :: rho4
    ! real(8),dimension(sys%DimUp*sys%DimDw,env%DimUp*env%DimDw) :: rho
    !
    if(allocated(rho))deallocate(rho)
    !
    DsysUp=sys%DimUp
    DsysDw=sys%DimDw
    DenvUp=env%DimUp
    DenvDw=env%DimDw
    !
    allocate(rho2(DsysUp*DenvUp,DsysDw*DenvDw));rho2=0d0    
    rho2 = transpose(reshape(psi, [DenvDw*DsysDw,DenvUp*DsysUp]))
    !
    allocate(rho4(DsysUp,DenvUp,DsysDw,DenvDw));rho4=0d0    
    do concurrent (iup=1:DsysUp*DenvUp, idw=1:DsysDw*DenvDw)
       i2 = state2indices(iup,[DsysUp,DenvUp])
       j2 = state2indices(idw,[DsysDw,DenvDw])
       rho4(i2(1),i2(2),j2(1),j2(2)) = rho2(iup,idw)
    enddo
    !
    allocate(rho(DsysUp*DsysDw,DenvUp*DenvDw));rho=0d0
    select case(to_lower(str(direction)))
    case ('left','l')
       do concurrent(s=1:DsysUp*DsysDw,sp=1:DsysUp*DsysDw)
          i2 = state2indices(s,[DsysUp,DsysDw])
          j2 = state2indices(sp,[DsysUp,DsysDw])
          sup = i2(1);spup=j2(1)
          sdw = i2(2);spdw=j2(2)       
          do concurrent(eup=1:DenvUp,edw=1:DenvDw)
             rho(s,sp) = rho(s,sp) + rho4(sup,eup,sdw,edw)*rho4(spup,eup,spdw,edw)
          enddo
       enddo
    case ('right','r')
       do concurrent(e=1:DenvUp*DenvDw,ep=1:DenvUp*DenvDw)
          i2 = state2indices(e,[DenvUp,DenvDw])
          j2 = state2indices(ep,[DenvUp,DenvDw])
          eup = i2(1);epup=j2(1)
          edw = i2(2);epdw=j2(2)       
          do concurrent(sup=1:DsysUp,sdw=1:DsysDw)
             rho(e,ep) = rho(e,ep) + rho4(sup,eup,sdw,edw)*rho4(sup,epup,sdw,epdw)
          enddo
       enddo
    end select
    !
  end function build_density_matrix

  function build_spin_density_matrix(self,rho,spin) result(rhoS)
    type(block),intent(in)                                         :: self
    real(8),dimension(self%DimUp*self%DimDw,self%DimUp*self%DimDw) :: rho
    character(len=*)                                               :: spin
    integer                                                        :: DimUp,DimDw
    real(8),dimension(:,:),allocatable                             :: rhoS
    integer                                                        :: iup,idw,i2(2),j2(2)
    integer                                                        :: s,sp,sup,sdw,spup,spdw
    !
    if(allocated(rhoS))deallocate(rhoS)
    !
    DimUp=self%DimUp
    DimDw=self%DimDw
    !
    select case(to_lower(str(spin)))
    case("up","u")
       allocate(rhoS(DimUp,DimUp))
       rhoS=0d0
       do concurrent(sup=1:DimUp,spup=1:DimUp)
          do concurrent(sdw=1:DimDw)
             s  = indices2state([sup,sdw],[DimUp,DimDw])
             sp = indices2state([spup,sdw],[DimUp,DimDw])
             rhoS(sup,spup) = rhoS(sup,spup) + rho(s,sp)
          enddo
       enddo
    case("dw","down","d")
       allocate(rhoS(DimDw,DimDw))
       rhoS = 0d0
       do concurrent(sdw=1:DimDw,spdw=1:DimDw)
          do concurrent(sup=1:DimUp)
             s  = indices2state([sup,sdw],[DimUp,DimDw])
             sp = indices2state([sup,spdw],[DimUp,DimDw])
             rhoS(sdw,spdw) = rhoS(sdw,spdw) + rho(s,sp)
          enddo
       enddo
    case default
       stop "build_spin_density_matrix error: spin is not defined"
    end select
  end function build_spin_density_matrix


  ! function build_density_matrix(nsys,nenv,psi,direction) result(rho)
  !   integer                      :: nsys
  !   integer                      :: nenv
  !   real(8),dimension(nsys*nenv) :: psi
  !   character(len=*)             :: direction
  !   real(8),dimension(nsys,nenv) :: rho_restricted
  !   real(8),dimension(nsys,nsys) :: rho
  !   rho_restricted = transpose(reshape(psi, [nenv,nsys]))
  !   select case(to_lower(str(direction)))
  !   case ('left','l')
  !      rho = matmul(rho_restricted,  transpose(rho_restricted) )
  !   case ('right','r')
  !      rho = matmul(transpose(rho_restricted), rho_restricted  )
  !   end select
  ! end function build_density_matrix



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




  pure function indices2state(ivec,Nvec) result(istate)
    integer,dimension(:),intent(in)          :: ivec
    integer,dimension(size(ivec)),intent(in) :: Nvec
    integer                                  :: istate,i
    istate=ivec(1)
    do i=2,size(ivec)
       istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
    enddo
  end function indices2state

  pure function state2indices(istate,Nvec) result(ivec)
    integer,intent(in)              :: istate
    integer,dimension(:),intent(in) :: Nvec
    integer,dimension(size(Nvec))   :: Ivec
    integer                         :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end function state2indices




end program test_iDMRG
