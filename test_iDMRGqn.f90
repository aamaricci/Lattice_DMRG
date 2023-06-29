program test_iDMRG
  USE SCIFOR
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
    type(block),intent(inout) :: self
    type(site)                :: dot
    character(len=*),optional :: grow
    character(len=16)         :: grow_
    type(block)               :: enl_self
    type(tbasis)              :: self_basis,dot_basis,enl_basis
    integer                   :: iqn
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    enl_self%length = self%length + 1
    enl_self%Dim    = self%Dim*dot%Dim
    !
    allocate(enl_self%sectors(size(self%sectors)))
    !
    select case(str(grow_))
    case ("left","l")
       call enl_self%put("H", &
            (self%operators%op("H").x.id(dot%dim)) +  (id(self%dim).x.dot%operators%op("H")) + &
            H2model(self,as_block(dot)))
       call enl_self%put("Cup", Id(self%dim).x.dot%operators%op("Cup"))
       call enl_self%put("Cdw", Id(self%dim).x.dot%operators%op("Cdw"))
       call enl_self%put("P"  , Id(self%dim).x.dot%operators%op("P"))
       do iqn=1,size(self%sectors)
          enl_basis = self%sectors(iqn)%basis().o.dot%sectors(iqn)%basis()
          call enl_self%set_basis( indx=iqn,  basis=enl_basis )       
       enddo
       !
    case ("right","r")
       call enl_self%put("H", &
            (id(dot%dim).x.self%operators%op("H")) +  (dot%operators%op("H").x.id(self%dim)) + &
            H2model(as_block(dot),self))
       call enl_self%put("Cup", dot%operators%op("Cup").x.Id(self%dim))
       call enl_self%put("Cdw", dot%operators%op("Cdw").x.Id(self%dim))
       call enl_self%put("P"  , dot%operators%op("P").x.Id(self%dim))
       do iqn=1,size(self%sectors)       
          enl_basis = dot%sectors(iqn)%basis().o.self%sectors(iqn)%basis()
          call enl_self%set_basis( indx=iqn,  basis=enl_basis )       
       enddo
       !
    end select
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



  subroutine get_sb_states(left,right,target_qn,sb_states,sb_sector)
    type(block)                      :: left,right
    real(8),dimension(2)             :: target_qn
    integer,dimension(:),allocatable :: sb_states
    type(sectors_list)               :: sb_sector
    integer                          :: ileft,iright
    integer                          :: i,j,istate
    real(8),dimension(:),allocatable :: left_qn,right_qn
    integer,dimension(:),allocatable :: left_map,right_map
    !
    if(allocated(sb_states))deallocate(sb_states)
    !
    do ileft=1,size(left%sectors(1))
       left_qn  = left%sectors(1)%qn(index=ileft)
       right_qn = target_qn - left_qn
       if(.not.right%sectors(1)%has_qn(right_qn))cycle
       !
       left_map = left%sectors(1)%map(qn=left_qn)
       right_map= right%sectors(1)%map(qn=right_qn)
       !
       ! print*,left_qn,"|",right_qn
       do i=1,size(left_map)
          do j=1,size(right_map)
             istate=right_map(j) + (left_map(i)-1)*right%Dim
             ! print*,left_map(i),":",right_map(j),"=",istate
             call append(sb_states, istate)
             call sb_sector%append(qn=left_qn,istate=size(sb_states))
          enddo
          ! print*,""
       enddo
    enddo
    !
  end subroutine get_sb_states



  subroutine single_dmrg_step(m,sys,env,energy)
    integer                            :: m
    type(block)                        :: sys,env
    type(block)                        :: esys,eenv
    type(sectors_list)                 :: sb_sector
    integer,dimension(:),allocatable   :: sb_states
    integer,dimension(:),allocatable   :: esys_map,eenv_map,sb_map
    real(8),dimension(:),allocatable   :: sb_qn,qn
    integer                            :: m_sb,m_s,m_e
    integer                            :: m_esys,m_eenv,Nsys,Nenv
    integer                            :: isys
    integer                            :: ienv
    integer                            :: i,j,iqn,Ncv,im
    type(tbasis)                       :: sys_basisL,sys_basisR
    real(8),dimension(:),allocatable   :: eig_values
    real(8),dimension(:,:),allocatable :: eig_basis,Hmat
    real(8),dimension(:),allocatable   :: rho_vec
    type(blocks_matrix)                :: rho_sectorL,rho_sectorR
    real(8),dimension(:),allocatable   :: evals,evalsL,evalsR
    real(8),dimension(:,:),allocatable :: Hsb,rhoL,rhoR,truncation_rhoL,truncation_rhoR
    real(8)                            :: truncation_errorL,truncation_errorR,energy

    !Check if blocks are valid ones
    if(.not.sys%is_valid())stop "single_dmrg_step error: sys is not a valid block"
    if(.not.env%is_valid())stop "single_dmrg_step error: env is not a valid block"
    !
    !Enlarge blocks
    print*,"Enlarge blocks:"
    call start_timer()
    esys = enlarge_block(sys,dot,grow='left')
    if(.not.esys%is_valid())stop "single_dmrg_step error: enlarged_sys is not a valid block"
    eenv = enlarge_block(env,dot,grow='right')
    if(.not.eenv%is_valid())stop "single_dmrg_step error: enlarged_env is not a valid block"
    print*,"Done"
    call stop_timer("Enlarge blocks")
    !
    !Get Enlarged Sys/Env actual dimension
    print*,"Enlarged Blocks dimensions:",esys%dim,eenv%dim    
    m_esys = esys%dim
    m_eenv = eenv%dim
    !
    !Build SuperBLock Sector
    call start_timer()
    call get_sb_states(esys,eenv,current_target_qn,sb_states,sb_sector)
    call stop_timer("Get SB states")

    m_sb = size(sb_states)
    ! print*,sb_states

    print*,"Super Block dimensions:", m_esys*m_eenv, m_sb
    !BUild SuperBlock matrix:
    spHsb = sp_kron(esys%operators%op("H"),id(m_eenv),sb_states) + &
         sp_kron(id(m_esys),eenv%operators%op("H"),sb_states)  + &
         H2model(esys,eenv,sb_states)
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
    print*,eig_values


    !
    call start_timer()
    do isys=1,size(sb_sector)
       sb_qn  = sb_sector%qn(index=isys)
       sb_map = sb_sector%map(index=isys)
       Nsys   = size(esys%sectors(1)%map(qn=sb_qn))
       Nenv   = size(eenv%sectors(1)%map(qn=(current_target_qn - sb_qn)))
       if(Nsys*Nenv==0)cycle
       !LEFT :
       qn = sb_qn
       call rho_sectorL%append(&
            build_density_matrix(Nsys,Nenv,eig_basis(:,1),sb_map,'left'),&
            qn=qn,map=esys%sectors(1)%map(qn=qn))
       !RIGHT:
       qn = current_target_qn-sb_qn
       call rho_sectorR%append(&
            build_density_matrix(Nsys,Nenv,eig_basis(:,1),sb_map,'right'),&
            qn=qn,map=eenv%sectors(1)%map(qn=qn))
    enddo
    call stop_timer("Build Rho")



    call start_timer()
    !LEFT:
    m_s = min(m,m_esys)
    allocate(truncation_rhoL(m_esys,m_s));truncation_rhoL=0d0
    call rho_sectorL%eigh(sort=.true.,reverse=.true.)
    do im=1,m_s
       rho_vec  = rho_sectorL%evec(m=im)
       esys_map = rho_sectorL%map(m=im)
       truncation_rhoL(esys_map,im) = rho_vec
       call sys_basisL%append( qn=rho_sectorL%qn(m=im) )
    enddo
    evalsL = rho_sectorL%evals()
    truncation_errorL = 1d0 - sum(evalsL(1:m_s))
    !
    m_e = min(m,m_eenv)
    allocate(truncation_rhoR(m_eenv,m_e));truncation_rhoR=0d0
    call rho_sectorR%eigh(sort=.true.,reverse=.true.)
    do im=1,m_e
       rho_vec  = rho_sectorR%evec(m=im)
       eenv_map = rho_sectorR%map(m=im)
       truncation_rhoR(eenv_map,im) = rho_vec
       call sys_basisR%append( qn=rho_sectorR%qn(m=im) )
    enddo
    evalsR = rho_sectorR%evals()
    truncation_errorR = 1d0 - sum(evalsR(1:m_e))
    call stop_timer("Get Rho_Left-Dot / Rho_Dot-Right")
    write(*,*)"Truncating      :",m<=m_esys,m<=m_eenv
    write(*,*)"Truncation dim  :",m_s,m_e
    write(*,*)"Truncation error:",truncation_errorL,truncation_errorR


    !find a clever way to iterate here:
    !Return updated Block:
    call start_timer()
    sys%length = esys%length
    sys%dim    = m_s
    call sys%put("H",rotate_and_truncate(esys%operators%op("H"),truncation_rhoL,m_esys,m_s))
    call sys%put("Cup",rotate_and_truncate(esys%operators%op("Cup"),truncation_rhoL,m_esys,m_s))
    call sys%put("Cdw",rotate_and_truncate(esys%operators%op("Cdw"),truncation_rhoL,m_esys,m_s))
    call sys%put("P",rotate_and_truncate(esys%operators%op("P"),truncation_rhoL,m_esys,m_s))
    do i=1,size(sys%sectors)
       call sys%set_basis(indx=i, basis=sys_basisL)
    enddo

    env%length = eenv%length
    env%dim    = m_e
    call env%put("H",rotate_and_truncate(eenv%operators%op("H"),truncation_rhoR,m_eenv,m_e))
    call env%put("Cup",rotate_and_truncate(eenv%operators%op("Cup"),truncation_rhoR,m_eenv,m_e))
    call env%put("Cdw",rotate_and_truncate(eenv%operators%op("Cdw"),truncation_rhoR,m_eenv,m_e))
    call env%put("P",rotate_and_truncate(eenv%operators%op("P"),truncation_rhoR,m_eenv,m_e))
    do i=1,size(env%sectors)
       call env%set_basis(indx=i,basis=sys_basisR)
    enddo
    call stop_timer("Update+Truncate Block")
    !

    !Clean memory:
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
    call sys_basisL%free()
    call sys_basisR%free()
    call sb_sector%free()
    call rho_sectorL%free()
    call rho_sectorR%free()
  end subroutine single_dmrg_step



  function build_density_matrix(Nsys,Nenv,psi,map,direction) result(rho)
    ! type(block)                        :: sys,env
    integer :: Nsys,Nenv
    real(8),dimension(:)               :: psi
    integer,dimension(nsys*nenv)       :: map
    character(len=*)                   :: direction
    real(8),dimension(:,:),allocatable :: rho
    real(8),dimension(sys%Dim,env%Dim) :: rho_tmp
    !
    if(allocated(rho))deallocate(rho)
    !
    rho_tmp = transpose(reshape(psi(map), [nenv,nsys]))
    !
    select case(to_lower(str(direction)))
    case ('left','l')
       allocate(rho(nsys,nsys));rho=0d0
       rho  = matmul(rho_tmp,  transpose(rho_tmp) )
    case ('right','r')
       allocate(rho(nenv,nenv));rho=0d0
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

end program
