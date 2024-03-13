MODULE HLOCAL
  USE SCIFOR
  USE INPUT_VARS
  USE AUX_FUNCS
  implicit none
  private

  integer,save                          :: Ns       !Number of levels per spin == Norb
  integer,save                          :: Ns_Ud
  integer,save                          :: Ns_Orb
  integer,dimension(:),allocatable,save :: Nfocks    !Dim of the local Fock space per channel
  integer,save                          :: Nfock    !Dim of the local Fock space
  integer,save                          :: Nsectors !Number of sectors
  integer,allocatable,dimension(:)      :: getDim             ! [Nsectors]

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


  public :: Init_LocalFock_Space
  public :: Delete_LocalFock_Space
  public :: Build_Hlocal_operator
  public :: Build_C_Operator
  public :: Build_Cdg_Operator
  public :: Build_Dens_Operator

  public :: Build_BasisStates
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
  public :: Build_FermionicSign
=======

>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
contains




  !##################################################################
  !##################################################################
  !               CREATE THE LOCAL FOCK SPACE
  !##################################################################
  !##################################################################
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
  subroutine Init_LocalFock_Space()
=======
  subroutine Init_LocalFock_Space(Norb)
    integer :: Norb
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
    integer :: DimUp,DimDw
    integer :: DimUps(1),DimDws(1)
    integer :: Nups(1),Ndws(1)
    integer :: isector
    !
    Ns       = Norb
    Ns_Orb   = Ns
    Ns_Ud    = 1
    allocate(Nfocks(2*Ns_Ud))
    Nfocks(1)= 2**Ns
    Nfocks(2)= 2**Ns
    Nfock    = product(Nfocks)!(2**Ns)*(2**Ns)
    Nsectors = ((Ns_Orb+1)*(Ns_Orb+1))**Ns_Ud
    !
    write(*,"(A)")"Init Local Fock space:"
    write(*,"(A,I15)") '# of levels           = ',Ns
    write(*,"(A,2I15)")'Hilbert dim per spin  = ',Nfocks
    write(*,"(A,I15)") 'Fock space dimension  = ',Nfock
    write(*,"(A,I15)") 'Number of sectors     = ',Nsectors
    write(*,"(A)")"--------------------------------------------"
    write(LOGfile,"(A)") bg_red('Warning: DIAGONAL HAMILTONIAN OPERATORS only')
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
    !
    call Write_FockStates()
    !
    return
  end subroutine Init_LocalFock_Space


  !##################################################################
  !##################################################################
  !               DELETE THE LOCAL FOCK SPACE
  !##################################################################
  !##################################################################
  subroutine Delete_LocalFock_Space()
    integer :: DimUp,DimDw
    integer :: DimUps(1),DimDws(1)
    integer :: Nups(1),Ndws(1)
    integer :: isector
    !
    Ns       = 0
    Ns_Orb   = Ns
    Ns_Ud    = 1
    deallocate(Nfocks)
    Nfock    = 0
    Nsectors = 0
    deallocate(getDim)
    return
  end subroutine Delete_LocalFock_Space




  !##################################################################
  !##################################################################
  !            BUILD LOCAL FOCK SPACE SECTORS (Nup,Ndw)
  !##################################################################
  !##################################################################
  subroutine Build_LocalFock_Sector(isector,self)
    integer,intent(in)                  :: isector
    type(local_fock_sector)             :: self
    integer                             :: iup,idw
    integer                             :: nup_,ndw_
    integer                             :: dim,iud
    !
    if(self%status)call Delete_LocalFock_Sector(self)
    !
    self%index = isector
    !
    allocate(self%H(Ns_Ud))
    allocate(self%DimUps(Ns_Ud))
    allocate(self%DimDws(Ns_Ud))
    allocate(self%Nups(Ns_Ud))
    allocate(self%Ndws(Ns_Ud))
    !
    call get_Nup(isector,self%Nups);self%Nup=sum(self%Nups)
    call get_Ndw(isector,self%Ndws);self%Ndw=sum(self%Ndws)
    call get_DimUp(isector,self%DimUps);self%DimUp=product(self%DimUps)
    call get_DimDw(isector,self%DimDws);self%DimDw=product(self%DimDws)
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
  !               BUILD LOCAL FOCK OPERATORS
  !##################################################################
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE: 
  !+-------------------------------------------------------------------+
  function build_Hlocal_operator(H) result(Hmat)
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
#ifdef _CMPLX
    complex(8),dimension(:,:)             :: H    ![Nspin*Norb,Nspin*Norb]
    complex(8),dimension(:,:),allocatable :: Hmat
    complex(8),dimension(Nspin,Ns,Ns)     :: Hloc
    complex(8)                            :: htmp
#else
    real(8),dimension(:,:)                :: H    ![Nspin*Norb,Nspin*Norb]
    real(8),dimension(:,:),allocatable    :: Hmat
    real(8),dimension(Nspin,Ns,Ns)        :: Hloc    
    real(8)                               :: htmp
#endif
    type(local_fock_sector)               :: sectorI
    integer                               :: isector,i,m,io,jo,iorb,ispin,jorb,ii,jj
    integer                               :: iup,idw,mup,mdw,k1,k2
    real(8)                               :: sg1,sg2
    integer                               :: nup(Ns),ndw(Ns),nvec(2*Ns)
    logical                               :: Jcondition
=======
    complex(8),dimension(:,:)             :: H    ![Nspin*Norb,Nspin*Norb]
    complex(8),dimension(:,:),allocatable :: Hmat
    complex(8),dimension(Nspin,Ns,Ns)     :: Hloc
=======
    real(8),dimension(:,:)             :: H    ![Nspin*Norb,Nspin*Norb]
    real(8),dimension(:,:),allocatable :: Hmat
    real(8),dimension(Nspin,Ns,Ns)     :: Hloc
>>>>>>> f63915b (Testing the code.):HLOCAL.f90
    type(local_fock_sector)               :: sectorI
    real(8)                               :: htmp
    integer                               :: isector,i,m,io,jo,iorb,ispin,jorb
    integer                               :: iup,idw,mup,mdw
    integer                               :: nup(Ns),ndw(Ns),nvec(2*Ns)
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
    !

#ifdef _CMPLX
    write(*,*)"Using CMPLX code"
    call wait(1000)
#endif
    
    
    if(allocated(Hmat))deallocate(Hmat)
    allocate(Hmat(Nfock,Nfock))
    !
    call assert_shape(H,[Nspin*Ns,Nspin*Ns],"build_hlocal_operator","Hloc")
    !
    do ispin=1,Nspin
       do iorb=1,Norb
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
          io = iorb+(ispin-1)*Norb      
          do jorb=1,Norb
             jo = jorb + (ispin-1)*Norb
=======
          io = iorb+(ispin-1)*Nspin          
          do jorb=1,Norb
             jo = jorb + (ispin-1)*Nspin
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
=======
          io = iorb+(ispin-1)*Norb      
          do jorb=1,Norb
             jo = jorb + (ispin-1)*Norb
>>>>>>> f63915b (Testing the code.):HLOCAL.f90
             Hloc(ispin,iorb,jorb) = H(io,jo)
          enddo
       enddo
    enddo
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
=======
    !
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
    !
    !
    !
    Hmat=zero
    do isector=1,Nsectors
       call Build_LocalFock_Sector(isector,SectorI)
       do i=1,sectorI%Dim
          m    = sectorI%H(1)%map(i)
          nvec = bdecomp(m,2*Ns)
          nup  = nvec(1:Ns)
          ndw  = nvec(Ns+1:2*Ns)
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
          ii   = m+1
=======
          m    = m+1
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
          !
          htmp = zero
          !LOCAL HAMILTONIAN PART:
          ! htmp = htmp - xmu*(sum(nup)+sum(ndw))
          do io=1,Ns
             htmp = htmp + Hloc(1,io,io)*nup(io)
             htmp = htmp + Hloc(Nspin,io,io)*ndw(io)
          enddo
          Hmat(ii,ii)=Hmat(ii,ii) + htmp
          !
          !
          !Density-density interaction: same orbital, opposite spins
          !Euloc=\sum=i U_i*(n_u*n_d)_i
          htmp = 0d0
          do io=1,Ns
             htmp = htmp + Uloc(io)*nup(io)*ndw(io)
          enddo
          !
          if(Norb>1)then
             !density-density interaction: different orbitals, opposite spins:
             ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
             do io=1,Ns
                do jo=io+1,Ns
                   htmp = htmp + Ust*(nup(io)*ndw(jo) + nup(jo)*ndw(io))
                enddo
             enddo
             !density-density interaction: different orbitals, parallel spins
             ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
             do io=1,Ns
                do jo=io+1,Ns
                   htmp = htmp + (Ust-Jh)*(nup(io)*nup(jo) + ndw(io)*ndw(jo))
                enddo
             enddo
          endif
          !if using the Hartree-shifted chemical potential: mu=0 for half-filling
          !sum up the contributions of hartree terms:
          if(hfmode)then
             do io=1,Ns
                htmp = htmp - 0.5d0*Uloc(io)*(nup(io)+ndw(io))
             enddo
             if(Norb>1)then
                do io=1,Ns
                   do jo=io+1,Ns
                      htmp=htmp-0.5d0*Ust*(nup(io)+ndw(io)+nup(jo)+ndw(jo))
                      htmp=htmp-0.5d0*(Ust-Jh)*(nup(io)+ndw(io)+nup(jo)+ndw(jo))
                   enddo
                enddo
             endif
          endif
          !
          Hmat(ii,ii)=Hmat(ii,ii) + htmp
          !
          !UP electrons
          do io=1,Ns
             do jo=1,Ns
                Jcondition = (Hloc(1,io,jo)/=zero) .AND. (Nup(jo)==1) .AND. (Nup(io)==0)
                if (Jcondition) then
                   call c(jo,m,k1,sg1)
                   call cdg(io,k1,k2,sg2)
                   jj    = binary_search(sectorI%H(1)%map,k2)
                   htmp = Hloc(1,io,jo)*sg1*sg2
                   Hmat(ii,jj)=Hmat(ii,jj) + htmp
                endif
             enddo
          enddo
          !
          !DW electrons
          do io=1,Ns
             do jo=1,Ns
                Jcondition = (Hloc(Nspin,io,jo)/=zero) .AND. (Ndw(jo)==1) .AND. (Ndw(io)==0)
                if (Jcondition) then
                   call c(jo+Ns,m,k1,sg1)
                   call cdg(io+Ns,k1,k2,sg2)
                   jj    = binary_search(sectorI%H(1)%map,k2)
                   htmp = Hloc(Nspin,io,jo)*sg1*sg2
                   Hmat(ii,jj)=Hmat(ii,jj) + htmp
                endif
             enddo
          enddo
          !
          !
       enddo
       call Delete_LocalFock_Sector(sectorI)
    enddo
  end function build_hlocal_operator





  function build_C_operator(iorb,ispin) result(Cmat)
    integer                               :: iorb,ispin
    type(local_fock_sector)               :: sectorI
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Cmat
#else
    real(8),dimension(:,:),allocatable    :: Cmat
#endif
=======
    complex(8),dimension(:,:),allocatable :: Cmat
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
=======
    real(8),dimension(:,:),allocatable :: Cmat
>>>>>>> f63915b (Testing the code.):HLOCAL.f90
    real(8)                               :: c_
    integer                               :: isector,i,m,l,Dim
    integer                               :: nvec(2*Ns),alfa
    !
    if(iorb<1.OR.iorb>Ns)stop "build_C_operator error: iorb<1 OR iorb>Ns. check."
    if(ispin<1.OR.ispin>2)stop "build_C_operator error: ispin<1 OR ispin>2. check."
    if(allocated(Cmat))deallocate(Cmat)
    allocate(Cmat(Nfock,Nfock))
    !
    alfa=iorb;if(ispin==2)alfa=iorb+Ns
    !
    Cmat=zero
    do isector=1,Nsectors
       call Build_LocalFock_Sector(isector,SectorI)
       Dim = SectorI%Dim
       do i=1,Dim
          m   = sectorI%H(1)%map(i)
          nvec= bdecomp(m,2*Ns)
          if(nvec(alfa) == 0) cycle
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
          call c(alfa,m,l,c_)
=======
          call c(alfa,m,l,c_)          
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
          Cmat(l+1,m+1) = one*c_
       enddo
       call Delete_LocalFock_Sector(sectorI)
    enddo
  end function build_C_operator

  function build_Cdg_operator(iorb,ispin) result(CDGmat)
    integer                               :: iorb,ispin
    type(local_fock_sector)               :: sectorI
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: CDGmat
#else
    real(8),dimension(:,:),allocatable    :: CDGmat
#endif
=======
    complex(8),dimension(:,:),allocatable :: CDGmat
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
=======
    real(8),dimension(:,:),allocatable :: CDGmat
>>>>>>> f63915b (Testing the code.):HLOCAL.f90
    real(8)                               :: cdg_
    integer                               :: isector,i,m,l,Dim
    integer                               :: nvec(2*Ns),alfa
    !
    if(iorb<1.OR.iorb>Ns)stop "build_CDG_operator error: iorb<1 OR iorb>Ns. check."
    if(ispin<1.OR.ispin>2)stop "build_CDG_operator error: ispin<1 OR ispin>2. check."
    if(allocated(CDGmat))deallocate(CDGmat)
    allocate(CDGmat(Nfock,Nfock))
    !
    alfa=iorb;if(ispin==2)alfa=iorb+Ns
    !
    CDGmat=0d0
    do isector=1,Nsectors
       call Build_LocalFock_Sector(isector,SectorI)
       Dim = SectorI%Dim
       do i=1,Dim
          m    = sectorI%H(1)%map(i)
          nvec = bdecomp(m,2*Ns)
          if(nvec(alfa) == 1) cycle
          call cdg(alfa,m,l,cdg_)          
          CDGmat(l+1,m+1) = one*cdg_
       enddo
       call Delete_LocalFock_Sector(sectorI)
    enddo
  end function build_cdg_operator



  function build_Dens_operator(iorb,ispin) result(Dmat)
    integer                               :: ispin,iorb
    type(local_fock_sector)               :: sectorI
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Dmat
#else
    real(8),dimension(:,:),allocatable :: Dmat
#endif
=======
    complex(8),dimension(:,:),allocatable :: Dmat
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
=======
    real(8),dimension(:,:),allocatable :: Dmat
>>>>>>> f63915b (Testing the code.):HLOCAL.f90
    real(8)                               :: dens
    integer                               :: imp,isector,i,m,Dim
    integer                               :: nvec(2*Ns),alfa
    !
    if(iorb<1.OR.iorb>Ns)stop "build_Dens_operator error: iorb<1 OR iorb>Ns. check."
    if(ispin<1.OR.ispin>2)stop "build_Dens_operator error: ispin<1 OR ispin>2. check."
    if(allocated(Dmat))deallocate(Dmat)
    allocate(Dmat(Nfock,Nfock))
    !
    alfa=iorb;if(ispin==2)alfa=iorb+Ns
    !
    Dmat=zero
    do isector=1,Nsectors
       call Build_LocalFock_Sector(isector,SectorI)
       Dim = SectorI%Dim
       do i=1,Dim
          m    = sectorI%H(1)%map(i)
          nvec = bdecomp(m,2*Ns)
          dens = dble(nvec(alfa))
          Dmat(m+1,m+1) = one*dens
       enddo
       call Delete_LocalFock_Sector(sectorI)
    enddo
  end function build_dens_operator




  !##################################################################
  !##################################################################
  !         WRITE OUT THE STATES OF THE LOCAL FOCK SPACE
  !##################################################################
  !##################################################################
  subroutine Write_FockStates()
    integer,dimension(:),allocatable :: Hvec
    integer                          :: i,iup,idw,Nup,Ndw,NN
    !
    if(allocated(Hvec))deallocate(Hvec)
    !
    NN = 2**Ns
    do idw = 0,NN-1
       Ndw = popcnt(idw)
       do iup = 0,NN-1
          Nup = popcnt(iup)
          i      = iup + idw*NN
          call print_conf(iup,Ns,.false.)
          call print_conf(idw,Ns,.true.)
       enddo
    enddo
  end subroutine Write_FockStates







  !##################################################################
  !##################################################################
  !         BUILD and RETURN THE LOCAL FOCK BASIS QN
  !##################################################################
  !##################################################################
  function Build_BasisStates() result(Hvec)
    integer,dimension(:),allocatable :: Hvec
    integer                          :: i,iup,idw,Nup,Ndw,NN
    !
    if(allocated(Hvec))deallocate(Hvec)
    !
    NN = 2**Ns
    do idw = 0,NN-1
       Ndw = popcnt(idw)
       do iup = 0,NN-1
          Nup = popcnt(iup)
          i      = iup + idw*NN
          call append(Hvec,nup)
          call append(Hvec,ndw)
       enddo
    enddo
  end function Build_BasisStates




<<<<<<< HEAD:src/FOCK/HLOCAL.f90
  !##################################################################
  !##################################################################
  !         BUILD and RETURN THE LOCAL FOCK BASIS QN
  !##################################################################
  !##################################################################
  function Build_BasisStates() result(Hvec)
    integer,dimension(:),allocatable :: Hvec
    integer                          :: i,iup,idw,Nup,Ndw,NN
    !
    if(allocated(Hvec))deallocate(Hvec)
    !
    NN = 2**Ns
    do idw = 0,NN-1
       Ndw = popcnt(idw)
       do iup = 0,NN-1
          Nup = popcnt(iup)
          i   = iup + idw*NN
          call append(Hvec,nup)
          call append(Hvec,ndw)
       enddo
    enddo
  end function Build_BasisStates




  !##################################################################
  !##################################################################
  !         BUILD and RETURN FERMIONIC SIGN MATRIX
  !##################################################################
  !##################################################################
  function Build_FermionicSign() result(P)
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: P
#else
    real(8),dimension(:,:),allocatable :: P
#endif
    integer,dimension(:),allocatable   :: Hvec
    integer                            :: i,iup,idw,Nup,Ndw,NN
    !
    if(allocated(Hvec))deallocate(Hvec)
    !
    NN = 2**Ns
    do idw = 0,NN-1
       Ndw = popcnt(idw)
       do iup = 0,NN-1
          Nup = popcnt(iup)
          i = (-1)**(Nup+Ndw)
          call append(Hvec,i)
       enddo
    enddo
    P = diag(one*Hvec)
  end function Build_FermionicSign



=======
  
>>>>>>> f63915b (Testing the code.):HLOCAL.f90

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
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
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
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp


  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Nlevels) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bjoin(ib,Ntot) result(i)
    integer                 :: Ntot
    integer,dimension(Ntot) :: ib
    integer                 :: i,j
    i=0
    do j=0,Ntot-1
       i=i+ib(j+1)*2**j
    enddo
  end function bjoin





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


<<<<<<< HEAD:src/FOCK/HLOCAL.f90



  !##################################################################
  !##################################################################
  !               PRINT LOCAL STATE
  !##################################################################
  !##################################################################  
  subroutine print_conf(i,Ntot,advance)
    integer :: dim,i,j,Ntot
    logical :: advance
    integer :: ivec(Ntot)
    ivec = bdecomp(i,Ntot)
    write(LOGfile,"(A1)",advance="no")"|"
    write(LOGfile,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    if(advance)then
       write(LOGfile,"(A1)",advance="yes")">"
    else
       write(LOGfile,"(A1)",advance="no")">"
    endif
  end subroutine print_conf


=======
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
END MODULE HLOCAL



!##################################################################
!##################################################################
!##################################################################
!##################################################################
!                          /_  __/ ____/ ___/_  __/
!                           / / / __/  \__ \ / /   
!                          / / / /___ ___/ // /    
!                         /_/ /_____//____//_/     
!##################################################################
!##################################################################
!##################################################################
!##################################################################
#ifdef _TEST
program testHLOCAL
  USE SCIFOR
  USE INPUT_VARS
  USE AUX_FUNCS
  USE HLOCAL
  implicit none
  character(len=64)                     :: finput
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable :: Docc,Cup,Cdw,CDGup,CDGdw,Dens,Hlocal,P
  complex(8),dimension(:,:),allocatable :: hloc
#else
  real(8),dimension(:,:),allocatable :: Docc,Cup,Cdw,CDGup,CDGdw,Dens,Hlocal,P
  real(8),dimension(:,:),allocatable :: hloc
#endif
=======
  complex(8),dimension(:,:),allocatable :: Docc,Cup,Cdw,CDGup,CDGdw,Dens,Hlocal
  complex(8),dimension(:,:),allocatable :: hloc
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
=======
  real(8),dimension(:,:),allocatable :: Docc,Cup,Cdw,CDGup,CDGdw,Dens,Hlocal
  real(8),dimension(:,:),allocatable :: hloc
>>>>>>> f63915b (Testing the code.):HLOCAL.f90
  integer,dimension(:),allocatable      :: Hvec

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')  
  call read_input(finput)


<<<<<<< HEAD:src/FOCK/HLOCAL.f90



  print*,"Imp Fock space"  
  call Init_LocalFock_Space()

  allocate(hloc(Nspin*Norb,Nspin*Norb))
  hloc=0d0;if(Norb==2)hloc=0.5d0*kron(pauli_0,pauli_z)
  call print_matrix(Hloc)

  print*,""
  print*,"H:"
  Hlocal = build_Hlocal_operator(hloc)
  call print_matrix(Hlocal)
=======
  print*,"Imp Fock space"  
  call Init_LocalFock_Space(Norb)
  allocate(hloc(nspin*norb,nspin*norb))
  hloc=0d0
  if(Norb==2)then
     hloc=0.5d0*pauli_z
  ! elseif(Norb==3)then
  !    hloc=0.5d0*spin1_z
  endif
  print*,""
  print*,"H:"
  Hlocal = build_Hlocal_operator(hloc)
  call print_mat(Hlocal)
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90

  print*,""
  print*,"C_up:"
  Cup = build_C_operator(iorb=1,ispin=1)
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
  call print_matrix(Cup)
  ! call print_matrix(kron(eye(2),Cop))
=======
  call print_mat(Cup)
  ! call print_mat(kron(eye(2),Cop))
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90

  print*,""
  print*,"C_dw:"
  Cdw = build_C_operator(iorb=1,ispin=2)
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
  call print_matrix(Cdw)
  ! call print_matrix(kron(Cop,dble(pauli_z)))
=======
  call print_mat(Cdw)
  ! call print_mat(kron(Cop,dble(pauli_z)))
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90

  print*,""
  print*,"CDG_up:"
  CDGup = build_CDG_operator(iorb=1,ispin=1)
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
  call print_matrix(CDGup)
  ! call print_matrix(transpose(kron(eye(2),Cop)))
=======
  call print_mat(CDGup)
  ! call print_mat(transpose(kron(eye(2),Cop)))
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90

  print*,""
  print*,"CDG_dw:"
  CDGdw = build_CDG_operator(iorb=1,ispin=2)
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
  call print_matrix(CDGdw)
  ! call print_matrix(transpose(kron(Cop,dble(pauli_z))))
=======
  call print_mat(CDGdw)
  ! call print_mat(transpose(kron(Cop,dble(pauli_z))))
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90


  print*,""
  print*,"Dens_up:"
  Dens = build_Dens_operator(iorb=1,ispin=1)
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
  call print_matrix(dens)
  ! call print_matrix(kron(eye(2),Dens))
=======
  call print_mat(dens)
  ! call print_mat(kron(eye(2),Dens))
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90

  print*,""
  print*,"Dens_dw:"  
  Dens = build_Dens_operator(iorb=1,ispin=2)
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
  call print_matrix(dens)
  ! call print_matrix(kron(Dens,eye(2)))
=======
  call print_mat(dens)
  ! call print_mat(kron(Dens,eye(2)))
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90


  print*,""
  print*,"Docc = CDG_up.C_up.CDG_dw.C_dw"  
  allocate(Docc, mold=Dens)
  Docc = CDGup.x.Cup.x.CDGdw.x.Cdw
<<<<<<< HEAD:src/FOCK/HLOCAL.f90
  call print_matrix(Docc)
=======
  call print_mat(Docc)
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90


  Hvec = Build_BasisStates()
  print*,size(Hvec)
  print*,Hvec

<<<<<<< HEAD:src/FOCK/HLOCAL.f90



  P = Build_FermionicSign()
  call print_matrix(P)
  print*,""
  call print_matrix(kSz(2*Norb))

  print*,diagonal(P)-diagonal(kSz(2*Norb))
=======
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90
  call Delete_LocalFock_space()




contains


<<<<<<< HEAD:src/FOCK/HLOCAL.f90

  recursive function KSz(n) result(A)
    integer, intent(in) :: n
    real(8)          :: A(2**n, 2**n)
    integer             :: d(2**n)
    integer             :: i
    d = szvec(n)
    A = zero
    forall(i=1:2**n)A(i,i) = one*d(i)
  end function KSz

  recursive function szvec(n) result(vec)
    integer,intent(in)      :: n
    integer,dimension(2**n) :: vec
    if(n==1)then
       vec = [1,-1]
    else
       vec = [szvec(n-1),-szvec(n-1)]
    endif
  end function szvec

=======
  subroutine print_mat(M)
    real(8),dimension(:,:) :: M
    integer                   :: i,j
    do i=1,size(M,1)
       write(*,"("//str(size(M,2))//"(F4.1,1x))")((M(i,j)),j=1,size(M,2))
    enddo
  end subroutine print_mat
>>>>>>> cc4f705 (Major Update: code entirely moved from DBLE to CMPLX.):HLOCAL.f90

end program testHLOCAL
#endif






