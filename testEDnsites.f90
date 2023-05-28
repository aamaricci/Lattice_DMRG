MODULE Hsimple
  USE SCIFOR
  implicit none
  private

  integer,save,public              :: Ns       !Number of levels per spin
  integer,save                     :: Nspin,Norb
  integer,save,public              :: Nfock    !Dim of the local Fock space
  integer,save,public              :: Nsectors !Number of sectors
  integer,allocatable,dimension(:) :: getDim             ! [Nsectors]

  !---------------- SECTOR-TO-LOCAL-FOCK SPACE STRUCTURE -------------------!
  type sector_map
     integer,dimension(:),allocatable          :: map
     logical                                   :: status=.false.
  end type sector_map

  type local_fock_sector
     integer                                   :: index       !
     type(sector_map),dimension(:),allocatable :: H
     integer,dimension(:),allocatable          :: DimUps
     integer,dimension(:),allocatable          :: DimDws
     integer                                   :: DimUp
     integer                                   :: DimDw
     integer                                   :: Dim
     integer,dimension(:),allocatable          :: Nups
     integer,dimension(:),allocatable          :: Ndws
     integer                                   :: Nup
     integer                                   :: Ndw
     logical                                   :: status=.false.
  end type local_fock_sector

  interface map_allocate
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate

  interface map_deallocate
     module procedure :: map_deallocate_scalar
     module procedure :: map_deallocate_vector
  end interface map_deallocate


  public :: Init_LocalFock

  public :: Build_H_operator
  public :: Build_C_Operator
  public :: Build_Cdg_Operator
  public :: Build_Dens_Operator
  public :: Build_Docc_Operator


  ! interface operator(.k.)
  !    module procedure :: d_kronecker_product
  !    module procedure :: c_kronecker_product
  ! end interface operator(.k.)
  ! public :: operator(.k.)

contains


  function d_kronecker_product(A,B) result(AxB)
    real(8),intent(in) :: A(:,:), B(:,:)
    integer            :: i,j
    integer            :: rowA,colA
    integer            :: rowB,colB
    real(8)            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))
    AxB = 0
    rowA=size(A,1) ; colA=size(A,2)
    rowB=size(B,1) ; colB=size(B,2)
    forall(i=1:rowA,j=1:colA)
       AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
    end forall
  end function d_kronecker_product
  !
  function c_kronecker_product(A,B) result(AxB)
    complex(8),intent(in) :: A(:,:), B(:,:)
    integer               :: i,j
    integer               :: rowA,colA
    integer               :: rowB,colB
    complex(8)            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))
    AxB = 0
    rowA=size(A,1) ; colA=size(A,2)
    rowB=size(B,1) ; colB=size(B,2)
    forall(i=1:rowA,j=1:colA)
       AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
    end forall
  end function c_kronecker_product


  !##################################################################
  !##################################################################
  !               CREATE THE LOCAL FOCK SPACE
  !##################################################################
  !##################################################################
  subroutine Init_LocalFock(Nsites)
    integer                 :: Nsites
    type(local_fock_sector) :: self
    integer                 :: DimUp,DimDw
    integer                 :: DimUps(1),DimDws(1)
    integer                 :: Nups(1),Ndws(1)
    integer                 :: isector,i
    !
    Nspin    = 1
    Norb     = 1
    Ns       = Nsites
    Nfock    = (2**Ns)*(2**Ns)
    Nsectors = (Ns+1)*(Ns+1)
    !
    write(*,"(A)")"Init Local Fock space:"
    write(*,"(A,I15)") '# of levels           = ',Ns
    write(*,"(A,I15)") 'Fock space dimension  = ',Nfock
    write(*,"(A,I15)") 'Number of sectors     = ',Nsectors
    write(*,"(A)")"--------------------------------------------"
    !
    !Allocate some indexing arrays
    allocate(getDim(Nsectors))          ;getDim=0
    !
    !
    do isector=1,Nsectors
       call get_DimUp(isector,DimUps)
       call get_DimDw(isector,DimDws)
       DimUp = product(DimUps)
       DimDw = product(DimDws)  
       getDim(isector)  = DimUp*DimDw
    enddo
    !
    do i=1,Nfock
       call print_conf(i,Ns,.true.)
    end do
    print*,""
    !
    return
  end subroutine Init_LocalFock



  !##################################################################
  !##################################################################
  !               BUILD FOCK OPERATORS
  !##################################################################
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE: 
  !+-------------------------------------------------------------------+
  function build_H_operator(Hij,xmu,uloc) result(Hmat)
    real(8),dimension(:,:)             :: Hij ![Ns,Ns]
    real(8)                            :: xmu,uloc
    type(local_fock_sector)            :: sectorI
    real(8),dimension(:,:),allocatable :: Hmat
    real(8)                            :: htmp,sg1,sg2
    integer                            :: isector,i,m,k,k1,k2,io,jo
    integer                            :: iup,idw,mup,mdw
    integer                            :: nup(Ns),ndw(Ns),nvec(2*Ns)
    logical                            :: Jcondition
    !
    if(allocated(Hmat))deallocate(Hmat)
    allocate(Hmat(Nfock,Nfock))
    !
    call assert_shape(Hij,[Ns,Ns],"build_hlocal_operator","Hloc")
    !
    !
    Hmat=0d0
    do m=1,Nfock
       nvec = bdecomp(m,2*Ns)
       nup  = nvec(1:Ns)
       ndw  = nvec(Ns+1:2*Ns)
       !
       !LOCAL HAMILTONIAN PART:
       htmp = 0d0
       do io=1,Ns
          htmp = htmp + Hij(io,io)*(nup(io)+ndw(io)) - xmu*(nup(io)+ndw(io))
          htmp = htmp + Uloc*nup(io)*ndw(io) - 0.5d0*Uloc*(nup(io)+ndw(io))
       enddo
       Hmat(m,m)=Hmat(m,m) + htmp
       !
       !Non local part
       do io=1,Ns
          do jo=1,Ns
             !UP electrons
             Jcondition = (Hij(io,jo)/=zero) .AND. (Nup(jo)==1) .AND. (Nup(io)==0)
             if (Jcondition) then
                call c(jo,m,k1,sg1)
                call cdg(io,k1,k,sg2)
                htmp = Hij(io,jo)*sg1*sg2
                !
                Hmat(k,m) = Hmat(k,m)+htmp
             endif
             !DW electrons
             Jcondition = (Hij(io,jo)/=zero) .AND. (Ndw(jo)==1) .AND. (Ndw(io)==0)
             if (Jcondition) then
                call c(jo+Ns,m,k1,sg1)
                call cdg(io+Ns,k1,k,sg2)
                htmp = Hij(io,jo)*sg1*sg2
                !
                Hmat(k,m) = Hmat(k,m)+htmp
             endif
          enddo
       enddo
       !
    enddo
  end function build_H_operator





  function build_C_operator(iorb,ispin) result(Cmat)
    integer                            :: ispin,iorb
    type(local_fock_sector)            :: sectorI
    real(8),dimension(:,:),allocatable :: Cmat
    real(8)                            :: c_
    integer                            :: imp,isector,i,min,mout
    integer                            :: nvec(2*Ns)
    !
    if(allocated(Cmat))deallocate(Cmat)
    allocate(Cmat(Nfock,Nfock))
    if(iorb<1.OR.iorb>Ns)stop "build_C_operator error: iorb<1 OR iorb>Ns. check."
    !
    imp = iorb+(ispin-1)*Ns
    !
    Cmat=0d0
    do min=1,Nfock
       nvec = bdecomp(min,2*Ns)
       if(nvec(imp) == 0) cycle
       call c(imp,min,mout,c_)          
       Cmat(mout,min) = c_
    enddo
  end function build_c_operator


  function build_Cdg_operator(iorb,ispin) result(CDGmat)
    integer                            :: ispin,iorb
    type(local_fock_sector)            :: sectorI
    real(8),dimension(:,:),allocatable :: CDGmat
    real(8)                            :: cdg_
    integer                            :: imp,isector,i,min,mout
    integer                            :: nvec(2*Ns)
    !
    if(allocated(CDGmat))deallocate(CDGmat)
    allocate(CDGmat(Nfock,Nfock))
    if(iorb<1.OR.iorb>Ns)stop "build_CDG_operator error: iorb<1 OR iorb>Ns. check."
    !
    imp = iorb+(ispin-1)*Norb
    !
    CDGmat=0d0
    do min=1,Nfock
       nvec = bdecomp(min,2*Ns)
       if(nvec(imp) == 1) cycle
       call cdg(imp,min,mout,cdg_)          
       CDGmat(mout,min) = cdg_
    enddo
  end function build_cdg_operator


  function build_Dens_operator(iorb,ispin) result(Dmat)
    integer                            :: ispin,iorb
    type(local_fock_sector)            :: sectorI
    real(8),dimension(:,:),allocatable :: Dmat
    real(8)                            :: dens
    integer                            :: imp,isector,i,m
    integer                            :: nvec(2*Ns)
    !
    if(allocated(Dmat))deallocate(Dmat)
    allocate(Dmat(Nfock,Nfock))
    if(iorb<1.OR.iorb>Ns)stop "build_Dens_operator error: iorb<1 OR iorb>Ns. check."
    !
    imp = iorb+(ispin-1)*Norb
    !
    Dmat=0d0
    do m=1,Nfock
       nvec = bdecomp(m,2*Ns)
       dens = dble(nvec(imp))
       Dmat(m,m) = dens
    enddo
  end function build_dens_operator


  function build_Docc_operator(iorb) result(Dmat)
    integer                            :: iorb
    type(local_fock_sector)            :: sectorI
    real(8),dimension(:,:),allocatable :: Dmat
    real(8)                            :: n_up,n_dw
    integer                            :: imp_up,imp_dw,isector,i,m
    integer                            :: nvec(2*Ns)
    !
    if(allocated(Dmat))deallocate(Dmat)
    allocate(Dmat(Nfock,Nfock))
    if(iorb<1.OR.iorb>Ns)stop "build_docc_operator error: iorb<1 OR iorb>Ns. check."
    !Put IF condition on iorb value: iorb=1,...,Ns
    imp_up = iorb
    imp_dw = iorb + Ns
    !
    Dmat=0d0
    do m=1,Nfock
       nvec = bdecomp(m,2*Ns)
       n_up = dble(nvec(imp_up))
       n_dw = dble(nvec(imp_dw))
       Dmat(m,m) = n_up*n_dw
    enddo
  end function build_docc_operator

















  !##################################################################
  !##################################################################
  !            BUILD LOCAL FOCK SPACE SECTORS (Nup,Ndw)
  !##################################################################
  !##################################################################
  subroutine Build_LocalFock_Sector(isector,self)
    integer,intent(in)      :: isector
    type(local_fock_sector) :: self
    integer                 :: iup,idw
    integer                 :: nup_,ndw_
    integer                 :: dim,iud
    !
    if(self%status)call Delete_LocalFock_Sector(self)
    !
    self%index = isector
    !
    allocate(self%H(1))
    allocate(self%DimUps(1))
    allocate(self%DimDws(1))
    allocate(self%Nups(1))
    allocate(self%Ndws(1))
    !
    call get_Nup(isector,self%Nups);self%Nup=sum(self%Nups);
    call get_Ndw(isector,self%Ndws);self%Ndw=sum(self%Ndws);
    call get_DimUp(isector,self%DimUps);self%DimUp=product(self%DimUps);
    call get_DimDw(isector,self%DimDws);self%DimDw=product(self%DimDws);
    self%Dim=self%DimUp*self%DimDw
    !
    call map_allocate(self%H,[self%Dim])
    !
    dim=0
    do idw=0,2**Ns-1
       ndw_= popcnt(idw)
       if(ndw_ /= self%Ndws(1))cycle
       do iup=0,2**Ns-1
          nup_ = popcnt(iup)
          if(nup_ /= self%Nups(1))cycle
          dim      = dim+1
          self%H(1)%map(dim) = iup + idw*2**Ns          
       enddo
    enddo
    !
    self%status=.true.
    !
  end subroutine Build_LocalFock_Sector




  subroutine Delete_LocalFock_Sector(self)
    type(local_fock_sector) :: self
    call map_deallocate(self%H)
    if(allocated(self%H))deallocate(self%H)
    if(allocated(self%DimUps))deallocate(self%DimUps)
    if(allocated(self%DimDws))deallocate(self%DimDws)
    if(allocated(self%Nups))deallocate(self%Nups)
    if(allocated(self%Ndws))deallocate(self%Ndws)
    self%index=0
    self%DimUp=0
    self%DimDw=0
    self%Dim=0
    self%Nup=0
    self%Ndw=0
    self%status=.false.
  end subroutine Delete_LocalFock_Sector







  !##################################################################
  !##################################################################
  !               RETRIEVE INFO PROCEDURES: Nup,Ndw,DimUp,DimDw
  !##################################################################
  !##################################################################
  subroutine get_Nup(isector,Nup)
    integer              :: isector,Nup(1)
    integer              :: i,count
    integer,dimension(2) :: indices_
    count=isector-1
    do i=1,2
       indices_(i) = mod(count,Ns+1)
       count      = count/(Ns+1)
    enddo
    Nup = indices_(2:2:-1)
  end subroutine get_Nup


  subroutine get_Ndw(isector,Ndw)
    integer              :: isector,Ndw(1)
    integer              :: i,count
    integer,dimension(2) :: indices_
    count=isector-1
    do i=1,2
       indices_(i) = mod(count,Ns+1)
       count      = count/(Ns+1)
    enddo
    Ndw = indices_(1:1:-1)
  end subroutine get_Ndw


  subroutine  get_DimUp(isector,DimUps)
    integer                :: isector,DimUps(1)
    integer                :: Nups(1)
    call get_Nup(isector,Nups)
    DimUps(1) = binomial(Ns,Nups(1))
  end subroutine get_DimUp


  subroutine get_DimDw(isector,DimDws)
    integer                :: isector,DimDws(1)
    integer                :: Ndws(1)
    call get_Ndw(isector,Ndws)
    DimDws(1) = binomial(Ns,Ndws(1))
  end subroutine get_DimDw


  subroutine indices2state(ivec,Nvec,istate)
    integer,dimension(:)          :: ivec
    integer,dimension(size(ivec)) :: Nvec
    integer                       :: istate,i
    istate=ivec(1)
    do i=2,size(ivec)
       istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
    enddo
  end subroutine indices2state

  subroutine state2indices(istate,Nvec,ivec)
    integer                       :: istate
    integer,dimension(:)          :: Nvec
    integer,dimension(size(Nvec)) :: Ivec
    integer                       :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end subroutine state2indices


  function iup_index(i,DimUp) result(iup)
    integer :: i
    integer :: DimUp
    integer :: iup
    iup = mod(i,DimUp);if(iup==0)iup=DimUp
  end function iup_index


  function idw_index(i,DimUp) result(idw)
    integer :: i
    integer :: DimUp
    integer :: idw
    idw = (i-1)/DimUp+1
  end function idw_index




  !##################################################################
  !##################################################################
  !               CREATION / DESTRUCTION OPERATORS
  !##################################################################
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=C_pos|in>  OR  |out>=C^+_pos|in> ; 
  !   the sign of |out> has the phase convention, pos labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(.not.btest(in-1,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in-1,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in-1,pos-1) + 1
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in-1,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in-1,l-1))fsgn=-fsgn
    enddo
    out = ibset(in-1,pos-1) + 1
  end subroutine cdg



  !##################################################################
  !##################################################################
  !               BITWISE OPERATIONS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i-1,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp






  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  elemental function binomial(n1,n2) result(nchoos)
    integer,intent(in) :: n1,n2
    real(8)            :: xh
    integer            :: i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial


  subroutine print_conf(i,Ns,advance)
    integer :: dim,i,j,unit_,Ntot,Ns,iup,idw
    logical :: advance
    integer :: ivec(2*Ns)
    unit_=6
    iup  = iup_index(i,2**Ns)
    idw  = idw_index(i,2**Ns)
    Ntot = 2*Ns
    ivec = bdecomp(i,Ntot)
    write(unit_,"(I3,1x,2I2)",advance="no")i,iup,idw
    write(unit_,"(A1)",advance="no")"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    if(advance)then
       write(unit_,"(A1)",advance="yes")">"
    else
       write(unit_,"(A1)",advance="no")">"
    endif
  end subroutine print_conf


  !##################################################################
  !##################################################################
  !               MAP ALLOCATION/DESTRUCTION
  !##################################################################
  !##################################################################
  subroutine map_allocate_scalar(H,N)
    type(sector_map) :: H
    integer          :: N
    if(H%status) call map_deallocate_scalar(H)
    allocate(H%map(N))
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:)       :: H
    integer,dimension(size(H))          :: N
    integer                             :: i
    do i=1,size(H)
       call map_allocate_scalar(H(i),N(i))
    enddo
  end subroutine map_allocate_vector


  subroutine map_deallocate_scalar(H)
    type(sector_map) :: H
    if(.not.H%status)then
       write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    H%status=.false.
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer                       :: i
    do i=1,size(H)
       call map_deallocate_scalar(H(i))
    enddo
  end subroutine map_deallocate_vector




END MODULE HSIMPLE






program testEDnsites
  USE SCIFOR
  USE AUX_FUNCS
  USE HSIMPLE
  implicit none


  integer                            :: Nsites
  real(8),dimension(:,:),allocatable :: hij,H
  real(8),dimension(:),allocatable   :: Evals
  real(8),dimension(:,:),allocatable :: Cup,Cdw,Cdgup,Cdgdw,Hp,Hp2,P,A,B,C,D,HH
  real(8),dimension(2,2) :: Cp!,Sz,Id
  real(8)                            :: ts
  integer                            :: i




  call parse_cmd_variable(nsites,"NSITES",default=1)
  ts    = -1d0
  call Init_LocalFock(Nsites)



  ! Id = pauli_0
  ! Sz = pauli_z
  Cp = (pauli_x+xi*pauli_y)/2d0





  allocate(Hij(Nsites,Nsites));Hij=0d0
  if(Nsites>1)then
     Hij(1,2) = ts
     do i=2,Nsites-1
        Hij(i,i-1) = ts
        Hij(i,i+1) = ts     
     enddo
     Hij(Nsites,Nsites-1) = ts
  endif
  call print_mat(Hij,"Hij")



  H =  build_H_operator(hij=Hij,xmu=0d0,uloc=0d0)
  call print_mat(H,"H")
  allocate(evals(size(H,1)))
  call eigh(H,Evals)
  do i=1,min(10,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)




  ! !C_up(i) = (Id(n) x ... x Id(1))_dw x (Id(n) x ... x C_p(i) x Sz(i-1) x...xSz(1))_up
  ! !        = Id^n x (Id^{n-i} x C_p(i) x Sz^{i-1})
  ! !C_dw(i) = (Id(n) x ... x C_p(i) x Sz(i-1) x ... x Sz(1))_dw x (Id(n) x ... x Id(1))_up
  ! !        = (Id^{n-i} x C_p(i) x Sz^{i-1}) x Sz^n 
  ! do i=1,Nsites
  !    Cup = build_C_operator(i,ispin=1)  
  !    call print_mat(Cup,"C_up_"//str(i))
  ! enddo

  ! if(Nsites==1)then
  !    !C_up(1) = Id^1 x C_p(1)
  !    Cup = KId(1).kx.Cp
  !    call print_mat(Cup,"newC_up_1")
  ! elseif(Nsites==2)then
  !    !Cup(1) = Id^2 x (Id^{2-1} x C_p x Sz^{1-1})
  !    !Cup(2) = Id^2 x (Id^{2-2} x C_p x Sz^{2-1})
  !    Cup = KId(2).kx.KId(1).kx.Cp
  !    call print_mat(Cup,"newC_up_1")
  !    Cup = KId(2).kx.Cp.kx.KSz(1)
  !    call print_mat(Cup,"newC_up_2")
  ! elseif(Nsites==3)then
  !    !Cup(1) = Id^3 x (Id^{3-1} x C_p x Sz^{1-1})
  !    !Cup(2) = Id^3 x (Id^{2-1} x C_p x Sz^{2-1})
  !    !Cup(3) = Id^3 x (Id^{1-1} x C_p x Sz^{3-1})
  !    Cup = KId(3).kx.(KId(2).kx.Cp)
  !    call print_mat(Cup,"newC_up_1")
  !    Cup = KId(3).kx.(KId(1).kx.Cp.kx.KSz(1))
  !    call print_mat(Cup,"newC_up_2")
  !    Cup = KId(3).kx.(Cp.kx.KSz(2))
  !    call print_mat(Cup,"newC_up_3")
  ! endif




  ! do i=1,Nsites
  !    Cdw = build_C_operator(i,ispin=2)  
  !    call print_mat(Cdw,"C_dw_"//str(i))
  ! enddo
  ! if(Nsites==1)then
  !    !C_dw(1) = C_p(1) x Sz
  !    Cdw = Cp.kx.KSz(1)
  !    call print_mat(Cdw,"newC_dw_1")
  ! elseif(Nsites==2)then
  !    !Cdw(1) = (Id x C_p) x Sz^2
  !    !Cdw(2) = (C_p x Sz) x Sz^2
  !    Cdw = (KId(1).kx.Cp).kx.KSz(2)
  !    call print_mat(Cdw,"newC_dw_1")
  !    Cdw = (Cp.kx.KSz(1)).kx.KSz(2)
  !    call print_mat(Cdw,"newC_dw_2")
  ! elseif(Nsites==3)then
  !    !Cdw(1) = (Id^{2} x C_p)          x Sz^3
  !    !Cdw(2) = (Id^{1} x C_p x Sz^{1}) x Sz^3
  !    !Cdw(3) = (C_p x Sz^{2})          x Sz^3
  !    Cdw = KId(2).kx.Cp.kx.KSz(3)
  !    call print_mat(Cdw,"newC_dw_1")
  !    Cdw = KId(1).kx.Cp.kx.KSz(1).kx.KSz(3)
  !    call print_mat(Cdw,"newC_dw_2")
  !    Cdw = Cp.kx.KSz(2).kx.KSz(3)
  !    call print_mat(Cdw,"newC_dw_3")
  ! endif








contains


  subroutine print_mat(M,name,file)
    real(8),dimension(:,:) :: M
    character(len=*)       :: name
    logical,optional       :: file
    logical                :: file_
    integer                :: i,j,stride,unit
    stride=2**Nsites
    file_=.true.;if(present(file))file_=file
    unit=6
    if(file_)open(free_unit(unit),file=str(name)//".dat")
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
    if(file_)close(unit)
  end subroutine print_mat


end program testEDnsites
