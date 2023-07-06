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
  type(block)                        :: my_block,dimerL,dimerR,trimerL,trimerR,tetramerL,tetramerR
  type(sparse_matrix)                :: spHsb,spH,spH_Ld,spH_dR,spRho
  type(site)                         :: dot
  real(8)                            :: gs_energy,target_Sz
  integer                            :: unit
  real(8),dimension(:,:),allocatable :: Hmatrix,Evecs,rho,rhoUP,rhoDW,rho2
  real(8),dimension(:),allocatable   :: Evals,ev,evUP,evDW
  type(sparse_matrix)       :: op

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"ts",finput,default=-1d0,comment="Hopping amplitude")
  target_Sz=0d0



  !
  !Init the single dot structure:
  dot = hubbard_site_ud(0d0,0d0)
  call dot%show()
  print*,dot%is_valid()
  !
  !Init block from single dot
  my_block=block(dot)
  call my_block%show()
  print*,my_block%is_valid()


  print*,"o->o"
  dimerL = enlarge_block(my_block,dot,grow='left')
  Hmatrix = as_matrix(build_Hblock(dimerL))  !(dimerL%operators%op("Hd")) + (Id(dimerL%DimDw).x.dimerL%operators%op("Hup")) + (dimerL%operators%op("Hdw").x.Id(dimerL%DimUp))
  print*,dimerL%is_valid()
  allocate(Evals(dimerL%dim))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""


  print*,"o-o->o"
  trimerL = enlarge_block(dimerL,dot,grow='left')
  print*,trimerL%is_valid()
  Hmatrix = as_matrix(build_Hblock(trimerL))
  allocate(Evals(trimerL%dim))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""

  print*,"o-o-o->o"
  tetramerL = enlarge_block(trimerL,dot,grow='left')
  print*,tetramerL%is_valid()
  Hmatrix = as_matrix(build_Hblock(tetramerL))
  allocate(Evals(tetramerL%dim))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""



  print*,"o<-o"
  dimerR = enlarge_block(my_block,dot,grow='right')
  print*,dimerR%is_valid()
  !(dimerR%operators%op("Hd")) + (Id(dimerR%DimDw).x.dimerR%operators%op("Hup")) + (dimerR%operators%op("Hdw").x.Id(dimerR%DimUp))  
  Hmatrix = as_matrix(build_Hblock(dimerR))
  allocate(Evals(dimerR%dim))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""



  print*,"o<-o-o"
  trimerR = enlarge_block(dimerR,dot,grow='right')
  print*,trimerR%is_valid()
  Hmatrix = as_matrix(build_Hblock(trimerR))
  allocate(Evals(trimerR%Dim))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""

  print*,"o<-o-o-o"
  tetramerR = enlarge_block(trimerR,dot,grow='right')
  print*,tetramerR%is_valid()
  Hmatrix = as_matrix(build_Hblock(tetramerR))
  allocate(Evals(tetramerR%Dim))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""





  !Init block from single dot
  print*,"o->o++o<-o"
  dimerL = enlarge_block(my_block,dot,grow='left')
  dimerR = enlarge_block(my_block,dot,grow='right')
  spH_Ld = build_Hblock(dimerL)!
  spH_dR = build_Hblock(dimerR)
  spHsb  = (spH_Ld.x.id(dimerR%dim)) + (id(dimerL%dim).x.spH_dR) + (Id(2**4).x.H2model(dimerL,dimerR)) + (H2model(dimerL,dimerR).x.Id(2**4))
  allocate(Evals(2))
  allocate(Evecs(4**4,2))
  call sp_eigh(sb_HxV,evals,evecs,&
       10,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,2
     print*,i,Evals(i)
  enddo
  deallocate(evals,evecs)
  print*,""


  !Init block from single dot
  print*,"o->o++o<-o"
  dimerL = enlarge_block(my_block,dot,grow='left')
  dimerR = enlarge_block(my_block,dot,grow='right')
  spHsb  = build_Hsb(dimerL,dimerR)
  allocate(Evals(2))
  allocate(Evecs(4**4,2))
  call sp_eigh(sb_HxV,evals,evecs,&
       10,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,2
     print*,i,Evals(i)
  enddo

  rho = build_density_matrix(dimerL,dimerR,evecs(:,1),'left')
  where(abs(rho)<1d-10)rho=0d0
  rhoUP = build_spin_density_matrix(dimerL,rho,'up')
  rhoDW = build_spin_density_matrix(dimerL,rho,'dw')
  call spRho%load(rho)
  print*,shape(spRho)
  call spRho%spy("dimer_rho")
  call spRho%free()

  call spRho%load(rhoUP)
  print*,shape(spRho)
  call spRho%spy("dimer_rhoUP")
  call spRho%free()

  call spRho%load(rhoDW)
  print*,shape(spRho)
  call spRho%spy("dimer_rhoDW")
  call spRho%free()

  allocate(rho2(dimerL%dim,dimerL%dim));rho2=0d0
  rho2 = kron(rhoDW,rhoUP)
  call spRho%load(rho2)
  print*,shape(spRho)
  call spRho%spy("dimer_rho2")
  call spRho%free()
  print*,all(abs(rho-rho2)<1d-10)
  deallocate(rho2)

  allocate(ev(dimerL%Dim))
  allocate(evUP(dimerL%DimUp))
  allocate(evDW(dimerL%DimDw))
  call eigh(rho,ev)
  call eigh(rhoUP,evUP)
  call eigh(rhoDW,evDW)

  call spRho%load(rho)
  call spRho%spy("eigh_rho")
  call spRho%free()

  call spRho%load(rhoUP)
  call spRho%spy("eigh_rhoUP")
  call spRho%free()

  call spRho%load(rhoDW)
  call spRho%spy("eigh_rhoDW")
  call spRho%free()


  allocate(rho2(dimerL%dim,dimerL%dim));rho2=0d0
  rho2 = kron(rhoDW,rhoUP)
  call spRho%load(rho2)
  call spRho%spy("eigh_rho2")
  call spRho%free()
  print*,all(abs(rho-rho2)<1d-10)
  deallocate(rho2)

  deallocate(rho,rhoUP,rhoDW)

  stop

  deallocate(evals,evecs)
  print*,""


  print*,"o--o->o++o<-o--o"
  trimerL = enlarge_block(dimerL,dot,grow='left')
  trimerR = enlarge_block(dimerR,dot,grow='right')
  spH_Ld = build_Hblock(trimerL)
  spH_dR = build_Hblock(trimerR)
  spHsb  = (spH_Ld.x.id(trimerR%dim)) + (id(trimerL%dim).x.spH_dR) + (Id(2**6).x.H2model(trimerL,trimerR)) + (H2model(trimerL,trimerR).x.Id(2**6))
  allocate(Evals(2))
  allocate(Evecs(4**6,2))
  call sp_eigh(sb_HxV,evals,evecs,&
       10,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,2
     print*,i,Evals(i)
  enddo

  allocate(rho(trimerL%dim,trimerL%dim));rho=0d0
  rho = build_density_matrix(trimerL,trimerR,evecs(:,1),'left')
  where(abs(rho)<1d-10)rho=0d0
  call spRho%load(rho)
  print*,shape(spRho)
  call spRho%spy("trimer_rho")
  call spRho%free()

  allocate(rhoUP(trimerL%DimUp,trimerL%DimUp));rhoUP=0d0
  rhoUP = build_density_matrix_UP(trimerL,rho)
  call spRho%load(rhoUP)
  print*,shape(spRho)
  call spRho%spy("trimer_rhoUP")
  call spRho%free()

  allocate(rhoDW(trimerL%DimDw,trimerL%DimDw));rhoDW=0d0
  rhoDW = build_density_matrix_DW(trimerL,rho)
  call spRho%load(rhoDW)
  print*,shape(spRho)
  call spRho%spy("trimer_rhoDW")
  call spRho%free()

  allocate(rho2(trimerL%dim,trimerL%dim));rho2=0d0
  rho2 = kron(rhoDW,rhoUP)
  call spRho%load(rho2)
  print*,shape(spRho)
  call spRho%spy("trimer_rho2")
  call spRho%free()

  print*,all(abs(rho-rho2)<1d-10)

  deallocate(rho,rho2,rhoUP,rhoDW)


  deallocate(evals,evecs)
  print*,""


  print*,"o--o--o->o++o<-o--o--o"
  tetramerL = enlarge_block(trimerL,dot,grow='left')
  tetramerR = enlarge_block(trimerR,dot,grow='right')
  spH_Ld = build_Hblock(tetramerL)
  spH_dR = build_Hblock(tetramerR)
  spHsb  = (spH_Ld.x.id(tetramerR%dim)) + (id(tetramerL%dim).x.spH_dR) + (Id(2**8).x.H2model(tetramerL,tetramerR)) + (H2model(tetramerL,tetramerR).x.Id(2**8))
  allocate(Evals(2))
  allocate(Evecs(4**8,2))
  call sp_eigh(sb_HxV,evals,evecs,&
       10,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,2
     print*,i,Evals(i)
  enddo
  deallocate(evals,evecs)
  print*,""



  print*,"o--o--o->o++o<-o--o--o"
  tetramerL = enlarge_block(trimerL,dot,grow='left')
  tetramerR = enlarge_block(trimerR,dot,grow='right')
  spHsb  = build_Hsb(tetramerL,tetramerR)
  allocate(Evals(2))
  allocate(Evecs(4**8,2))
  call sp_eigh(sb_HxV,evals,evecs,&
       10,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,2
     print*,i,Evals(i)
  enddo


  rho = build_density_matrix(tetramerL,tetramerR,evecs(:,1),'left')
  where(abs(rho)<1d-10)rho=0d0
  call spRho%load(rho)
  print*,shape(spRho)
  call spRho%spy("tetramer_rho")
  call spRho%free()

  rhoUP = build_spin_density_matrix(tetramerL,rho,'up')
  rhoDW = build_spin_density_matrix(tetramerL,rho,'dw')
  call spRho%load(rhoUP)
  print*,shape(spRho)
  call spRho%spy("tetramer_rhoUP")
  call spRho%free()

  call spRho%load(rhoDW)
  print*,shape(spRho)
  call spRho%spy("tetramer_rhoDW")
  call spRho%free()

  allocate(rho2(tetramerL%dim,tetramerL%dim));rho2=0d0
  rho2 = kron(rhoDW,rhoUP)
  call spRho%load(rho2)
  print*,shape(spRho)
  call spRho%spy("tetramer_rho2")
  call spRho%free()

  print*,all(abs(rho-rho2)<1d-10)

  deallocate(rho,rho2,rhoUP,rhoDW)


  ! allocate(rho(tetramerL%dim,tetramerL%dim));rho=0d0
  ! rho = build_density_matrix(tetramerL,tetramerR,evecs(:,1),'left')
  ! where(abs(rho)<1d-10)rho=0d0
  ! call spRho%load(rho)
  ! print*,shape(spRho)
  ! call spRho%spy("tetramer_rho")
  ! call spRho%free()

  ! allocate(rhoUP(tetramerL%DimUp,tetramerL%DimUp));rhoUP=0d0
  ! rhoUP = build_density_matrix_UP(tetramerL,rho)
  ! call spRho%load(rhoUP)
  ! print*,shape(spRho)
  ! call spRho%spy("tetramer_rhoUP")
  ! call spRho%free()

  ! allocate(rhoDW(tetramerL%DimDw,tetramerL%DimDw));rhoDW=0d0
  ! rhoDW = build_density_matrix_DW(tetramerL,rho)
  ! call spRho%load(rhoDW)
  ! print*,shape(spRho)
  ! call spRho%spy("tetramer_rhoDW")
  ! call spRho%free()

  ! allocate(rho2(tetramerL%dim,tetramerL%dim));rho2=0d0
  ! rho2 = kron(rhoDW,rhoUP)
  ! call spRho%load(rho2)
  ! print*,shape(spRho)
  ! call spRho%spy("tetramer_rho2")
  ! call spRho%free()

  ! print*,all(abs(rho-rho2)<1d-10)

  ! deallocate(rho,rho2,rhoUP,rhoDW)


  deallocate(evals,evecs)
  print*,""





contains


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


  function build_density_matrix_UP(self,rho) result(rhoUP)
    type(block),intent(in)                                         :: self
    integer                                                        :: DimUp,DimDw
    real(8),dimension(self%DimUp*self%DimDw,self%DimUp*self%DimDw) :: rho
    real(8),dimension(:,:),allocatable                             :: rhoUP ![DimUp,DimUp]
    integer                                                        :: iup,idw,i2(2),j2(2)
    integer                                                        :: s,sp,sup,sdw,spup,spdw
    !
    if(allocated(rhoUP))deallocate(rhoUP)
    !
    DimUp=self%DimUp
    DimDw=self%DimDw
    !
    allocate(rhoUP(DimUp,DimUp))
    rhoUP=0d0
    do concurrent(sup=1:DimUp,spup=1:DimUp)
       do concurrent(sdw=1:DimDw)
          s  = indices2state([sup,sdw],[DimUp,DimDw])
          sp = indices2state([spup,sdw],[DimUp,DimDw])
          rhoUP(sup,spup) = rhoUP(sup,spup) + rho(s,sp)
       enddo
    enddo
  end function build_density_matrix_UP


  function build_density_matrix_DW(self,rho) result(rhoDW)
    type(block),intent(in)                                         :: self
    integer                                                        :: DimUp,DimDw
    real(8),dimension(self%DimUp*self%DimDw,self%DimUp*self%DimDw) :: rho
    real(8),dimension(:,:),allocatable                             :: rhoDW ![DimDw,DimDw]
    integer                                                        :: iup,idw,i2(2),j2(2)
    integer                                                        :: s,sp,sup,sdw,spup,spdw
    !
    if(allocated(rhoDW))deallocate(rhoDW)
    !
    DimUp=self%DimUp
    DimDw=self%DimDw
    !
    allocate(rhoDW(DimDw,DimDw))
    rhoDW = 0d0
    do concurrent(sdw=1:DimDw,spdw=1:DimDw)
       do concurrent(sup=1:DimUp)
          s  = indices2state([sup,sdw],[DimUp,DimDw])
          sp = indices2state([sup,spdw],[DimUp,DimDw])
          rhoDW(sdw,spdw) = rhoDW(sdw,spdw) + rho(s,sp)
       enddo
    enddo
  end function build_density_matrix_DW




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
    !
    ! spL    = build_Hblock(left)
    ! spR    = build_Hblock(right)
    ! spH    = (spL.x.id(right%dim)) + (id(left%dim).x.spR) + (Id(DimDw).x.H2model(left,right)) + (H2model(left,right).x.Id(DimUp))
  end function build_Hsb




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
    enl_self%Dim    = self%Dim*dot%Dim
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






  subroutine print_mat(M,name,n)
    real(8),dimension(:,:) :: M
    character(len=*) :: name
    integer :: i,j,stride,unit,n
    stride=n!2**n
    open(free_unit(unit),file=str(name)//".dat")
    write(unit,*)"Matrix: "//str(name)
    do i=1,size(M,1)
       do j=1,size(M,2)
          write(unit,"(("//str(stride)//"F6.1))",advance="no")M(i,j)
          if(mod(j,stride)==0)write(unit,"(A1)",advance="no")""
       enddo
       write(unit,*)
       if(mod(i,stride)==0)write(unit,*)
    enddo
    write(unit,*)""
    close(unit)
  end subroutine print_mat






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



end program testEDkron
