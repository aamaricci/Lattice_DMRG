module MATRIX_GRAPH
  USE SCIFOR, only: free_unit,reg,file_length,assert_shape
  USE AUX_FUNCS, only: show_fmt,append
  implicit none
  private


#ifdef _CMPLX
  complex(8) :: zero=cmplx(0d0,0d0)
#else
  real(8)    :: zero=0d0
#endif


  type link
     integer            :: siteI,siteJ
     integer            :: orbI,orbJ
     integer            :: spin
#ifdef _CMPLX
     complex(8)         :: Hval
#else
     real(8)            :: Hval
#endif
     type(link),pointer :: next !link to next box (chain)
  end type link


  type hij_matrix
     type(link),pointer                      :: root     !head/root of the list
     character(len=:),allocatable            :: file     !Name of the W90 file 
     integer                                 :: Size=0   !Number of hopping elements
     integer,allocatable,dimension(:)        :: Nsites   !Number of sites per orbital
     integer                                 :: Ns=0     !Number of levels per spin 
     integer                                 :: Norb=0   !Number of orbitals
     integer                                 :: Nspin=0  !Number of spin channels (1 or 2)
#ifdef _CMPLX
     complex(8),allocatable,dimension(:,:,:) :: Hij      !The H(Ri,Rj)_ab^ss Hamiltonian
     complex(8),allocatable,dimension(:,:,:) :: Hloc     !The local part of H(Ri,Ri)_ab^ss
#else
     real(8),allocatable,dimension(:,:,:)    :: Hij      !The H(Ri,Rj)_ab^ss Hamiltonian
     real(8),allocatable,dimension(:,:,:)    :: Hloc     !The local part of H(Ri,Ri)_ab^ss
#endif
     logical                                 :: built  =.false. !Hij is built
     logical                                 :: status =.false. !Allocated
   contains
     procedure,pass :: init       => hij_init_matrix
     procedure,pass :: free       => hij_free_matrix
     procedure,pass :: info       => hij_info_matrix
     procedure,pass :: write      => hij_write_matrix
     procedure,pass :: read       => hij_read_matrix
     procedure,pass :: add_link   => hij_add_link_matrix
     procedure,pass :: build      => hij_build_matrix
     generic        :: get_Hij    => Hij_get_Hij_main_matrix,Hij_get_Hij_spin_matrix
     generic        :: get_Hloc   => Hij_get_Hloc_main_matrix,Hij_get_Hloc_spin_matrix
     procedure,pass :: pack_Hij   => Hij_get_Hij_pack_matrix
     procedure,pass :: pack_Hloc  => Hij_get_Hloc_pack_matrix
     procedure,pass :: &
          Hij_get_Hij_main_matrix ,Hij_get_Hij_spin_matrix,&
          Hij_get_Hloc_main_matrix,Hij_get_Hloc_spin_matrix

  end type hij_matrix


  interface hij_matrix
     module procedure :: hij_constructor_matrix
  end interface hij_matrix


  public :: hij_matrix

  
contains



  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the Hij structure: constructor
  !+------------------------------------------------------------------+
  subroutine hij_init_matrix(self,Nsites,Nspin)
    class(hij_matrix),intent(inout) :: self
    integer,dimension(:),intent(in) :: Nsites
    integer,intent(in),optional     :: Nspin
    integer                         :: Nspin_
    integer                         :: Ns
    integer                         :: Norb
    !
    Nspin_ = 1; if(present(Nspin))Nspin_=Nspin
    !
    if(self%status)call self%free()
    !
    Norb = size(Nsites)
    Ns   = sum(Nsites)
    !
    allocate(Self%root)
    Self%root%next => null()
    !
    allocate(Self%Nsites(Norb))
    Self%Nsites  = Nsites
    Self%file   = ""
    Self%size   = 0
    Self%Ns     = Ns
    Self%Norb   = Norb
    Self%Nspin  = Nspin_
    !
    allocate( Self%Hij(Nspin_,Ns,Ns) )
    Self%Hij    = zero
    allocate( Self%Hloc(Nspin_,Ns,Ns) )
    Self%Hloc   = zero
    Self%built  =.false.
    Self%status =.true.
    return
  end subroutine hij_init_matrix
  !

  function hij_constructor_matrix(Nsites,Nspin) result(self)
    integer,dimension(:),intent(in) :: Nsites
    integer,intent(in),optional     :: Nspin
    integer                         :: Nspin_
    type(hij_matrix)                :: self
    !
    Nspin_ = 1; if(present(Nspin))Nspin_=Nspin
    call self%init(Nsites,Nspin_)
  end function hij_constructor_matrix
  !



  !+------------------------------------------------------------------+
  !PURPOSE:  free the graph matrix object
  !+------------------------------------------------------------------+
  subroutine hij_free_matrix(self)
    class(hij_matrix),intent(inout) :: self
    type(link),pointer                 :: p,c
    !
    if(.not.Self%status)return
    !
    do
       p => Self%root
       c => p%next                 !current is the first node (root's next)
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next => null()
       call del_link(c)
       deallocate(c)
    end do
    deallocate(Self%root)
    deallocate(Self%Nsites)
    deallocate(Self%file)
    deallocate(Self%Hij)
    deallocate(Self%Hloc)
    Self%size   = 0
    Self%Ns     = 0
    Self%Norb   = 0
    Self%Nspin  = 0
    Self%built  =.false.
    Self%status =.false.
  end subroutine hij_free_matrix
  !








  !+------------------------------------------------------------------+
  !PURPOSE:  print info for the graph matrix 
  !+------------------------------------------------------------------+
  subroutine hij_info_matrix(self)
    class(hij_matrix),intent(inout) :: self
    !
    if(.not.Self%status)then
       write(*,"(A)")"Hij_info WARNING: Self not allocated"
       return
    endif
    !
    write(*,"(A,10I8)")"Nsites  =",Self%Nsites
    write(*,"(A,I8)")  "Ns      =",Self%Ns
    write(*,"(A,I8)")  "Nspin   =",Self%Nspin
    write(*,"(A,I8)")  "Norb    =",Self%Norb
    write(*,"(A,I8)")  "# link  =",Self%size
    if(Self%file /= "")&
         write(*,"(A,1x,A)")"file    =",Self%file
    write(*,"(A,L8)")  "built H =",Self%built
    write(*,"(A,A)")     "Order   = ","Sites -> Orbs -> Spin"
    return
  end subroutine hij_info_matrix
  !




  !+------------------------------------------------------------------+
  !PURPOSE:  write the graph matrix to stdout or unit or file 
  !+------------------------------------------------------------------+  
  subroutine hij_write_matrix(self,unit,file)
    class(hij_matrix),intent(inout) :: self
    character(len=*),optional       :: file
    integer,optional                :: unit
    integer                         :: unit_
    type(link),pointer              :: c           
    !
    if(.not.Self%status)stop "Hij_write ERROR: Self not allocated"
    !
    !
    unit_ = 6 ; if(present(unit))unit_=unit
    if(present(file))open(free_unit(unit_),file=reg(file))
    c => Self%root%next
    do
       if(.not.associated(c))exit
       call print_link(c,unit_)
       c => c%next
    enddo
    if(present(file))close(unit_)
    return
  end subroutine hij_write_matrix
  !

  !+------------------------------------------------------------------+
  !PURPOSE:  Read  graph matrix from file
  !+------------------------------------------------------------------+  
  subroutine Hij_read_matrix(self,file)
    class(hij_matrix),intent(inout) :: self
    character(len=*)                :: file
    integer                         :: Nsize
    integer                         :: unitIO
    integer                         :: ih,i,j,a,b,spin
    real(8)                         :: re,im
    !
    if(.not.Self%status)stop "Hij_read ERROR: Self not allocated"
    !
    !Read from file initial info OR reconstruct them is Header is missing
    open(free_unit(unitIO),file=reg(file),status="old",action="read")
    Nsize = file_length(reg(file))
    Self%file = reg(file)
    Self%Size = Nsize
    !
    do ih=1,Self%Size
#ifdef _CMPLX
       read(unitIO,*)i,j,a,b,spin,re,im
       call self%add_link(i,j,a,b,spin,dcmplx(re,im))
#else
       read(unitIO,*)i,j,a,b,spin,re
       call self%add_link(i,j,a,b,spin,re)
#endif
    enddo
    close(unitIO)
    return
  end subroutine Hij_read_matrix
  !










  !+------------------------------------------------------------------+
  !PURPOSE:  add link to graph matrix to build the lattice structure
  !+------------------------------------------------------------------+  
  subroutine Hij_add_link_matrix(self,siteI,siteJ,orbI,orbJ,spin,Hval)
    class(hij_matrix),intent(inout) :: self
    integer                         :: siteI,siteJ
    integer                         :: orbI,orbJ
    integer                         :: spin
#ifdef _CMPLX
    complex(8) ,intent(in)          :: Hval
#else
    real(8) ,intent(in)             :: Hval
#endif
    type(link),pointer              :: p,c
    integer                         :: k
    !
    !
    if(.not.Self%status)stop "Hij_add_link_matrix ERROR: Self not allocated"
    !
    !
    if(spin>Self%Nspin .OR. spin<=0)&
         stop "Hij_add_link_matrix ERROR: spin index either > Nspin OR < 0"
    if(orbI>Self%Norb  .OR. orbJ>Self%Norb  .OR. orbI<=0 .OR. orbJ<=0)&
         stop "Hij_add_link_matrix ERROR: orbI or orbJ either > Norb OR < 0"
    if(siteI>Self%Nsites(orbI) .OR. siteJ>Self%Nsites(orbJ) .OR. siteI<=0 .OR. siteI<=0 )&
         stop "Hij_add_link_matrix ERROR: siteI or siteJ either > Nlat OR < 0"
    !
    p => Self%root
    c => p%next
    do k=1,Self%size
       if(.not.associated(c))exit
       p => c
       c => c%next
    end do
    allocate(p%next)
    !
    p%next%siteI = siteI
    p%next%siteJ = siteJ
    p%next%orbI  = orbI
    p%next%orbJ  = orbJ
    p%next%spin  = spin
    p%next%Hval  = Hval
    Self%size=Self%size+1
    !
    if(.not.associated(c))then !end of the Self special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
  end subroutine Hij_add_link_matrix
  !










  !+------------------------------------------------------------------+
  !PURPOSE:  Build the actual Hij from graph matrix 
  !+------------------------------------------------------------------+  
  subroutine Hij_build_matrix(self)
    class(hij_matrix),intent(inout) :: self
    type(link),pointer              :: c           
    integer                         :: i,j,a,b,spin
    integer                         :: io,jo
#ifdef _CMPLX
    complex(8)                      :: val
#else
    real(8)                         :: val
#endif
    !
    if(.not.Self%status)stop "Hij_build ERROR: Self not allocated"
    !
    c => Self%root%next
    do
       if(.not.associated(c))exit
       call get_link(c,i,j,a,b,spin,val)
       io = indices_link(i,a,spin)
       jo = indices_link(j,b,spin)
       Self%Hij(spin,io,jo) = val
       if(i==j)Self%Hloc(spin,io,jo) = val
       c => c%next
    enddo
    do spin=1,Self%Nspin
       if(.not.herm_check(Self%Hij(spin,:,:),1d-6))stop "Hij_build ERROR: Hij is not Hermitian"
    enddo
    Self%built=.true.
    return
  contains
    function indices_link(siteI,orbI,spin) result(Indx)
      integer,intent(in)          :: siteI,orbI
      integer,intent(in),optional :: spin
      integer                     :: Indx
      select case(orbI)
      case(1)
         Indx = siteI
      case default
         Indx = siteI + (orbI-1)*Self%Nsites(orbI-1)
      end select
      if(present(spin))Indx = Indx + (spin-1)*Self%Ns
    end function indices_link
  end subroutine Hij_build_matrix
  !












  !+------------------------------------------------------------------+
  !PURPOSE:  Retrieve the  actual Hij matrix from graph matrix type
  !+------------------------------------------------------------------+  
  function Hij_get_Hij_main_matrix(self) result(Hij)
    class(hij_matrix),intent(inout)           :: self
#ifdef _CMPLX
    complex(8),dimension(:,:,:),allocatable :: Hij
#else
    real(8),dimension(:,:,:),allocatable    :: Hij
#endif
    integer                                   :: Nlso,Nspin
    !
    if(.not.Self%status)stop "Hij_get_Hij_matrix ERROR: Self not allocated"
    !
    Nlso = Self%Ns
    Nspin= Self%Nspin
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij(Nspin,Nlso,Nlso))
    if(.not.Self%built)call self%build()
    Hij = Self%Hij
  end function Hij_get_Hij_main_matrix

  function Hij_get_Hij_spin_matrix(self,spin) result(Hij)
    class(hij_matrix),intent(inout)       :: self
    integer,intent(in)                    :: spin
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hij
#else
    real(8),dimension(:,:),allocatable    :: Hij
#endif
    integer                               :: Nlso
    !
    if(.not.Self%status)stop "Hij_get_Hij_ispin_matrix ERROR: Self not allocated"
    !
    Nlso = Self%Ns
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij(Nlso,Nlso))
    if(.not.Self%built)call self%build()
    Hij = Self%Hij(spin,:,:)
  end function Hij_get_Hij_spin_matrix











  !+------------------------------------------------------------------+
  !PURPOSE:  Retrieve the actual Hloc matrix from graph matrix type
  !+------------------------------------------------------------------+  
  function Hij_get_Hloc_main_matrix(self) result(Hloc)
    class(hij_matrix),intent(inout)           :: self
#ifdef _CMPLX
    complex(8),dimension(:,:,:),allocatable :: Hloc
#else
    real(8),dimension(:,:,:),allocatable    :: Hloc
#endif
    integer                                   :: Nlso,Nspin
    !
    if(.not.Self%status)stop "Hij_get_Hloc_matrix ERROR: Self not allocated"
    !
    Nlso = Self%Ns
    Nspin= Self%Nspin
    if(allocated(Hloc))deallocate(Hloc)
    allocate(Hloc(Nspin,Nlso,Nlso))
    if(.not.Self%built)call self%build()
    Hloc = Self%Hloc
  end function Hij_get_Hloc_main_matrix

  function Hij_get_Hloc_spin_matrix(self,spin) result(Hloc)
    class(hij_matrix),intent(inout)       :: self
    integer,intent(in)                    :: spin
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hloc
#else
    real(8),dimension(:,:),allocatable    :: Hloc
#endif
    integer                               :: Nlso
    !
    if(.not.Self%status)stop "Hij_get_Hloc_ispin_matrix ERROR: Self not allocated"
    !
    Nlso = Self%Ns
    if(allocated(Hloc))deallocate(Hloc)
    allocate(Hloc(Nlso,Nlso))
    if(.not.Self%built)call self%build()
    Hloc = Self%Hloc(spin,:,:)
  end function Hij_get_Hloc_spin_matrix






  !+------------------------------------------------------------------+
  !PURPOSE:  Pack graph matrix Hij component into a Hij matrix
  !+------------------------------------------------------------------+  
  function Hij_get_Hij_pack_matrix(self) result(Hij)
    class(hij_matrix),intent(inout)       :: self
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hij
#else
    real(8),dimension(:,:),allocatable    :: Hij
#endif
    integer                               :: Nlso
    !
    if(.not.Self%status)stop "Hij_get_Hij_ispin_matrix ERROR: Self not allocated"
    !
    Nlso = self%Nspin*Self%Ns
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij(Nlso,Nlso))
    if(.not.Self%built)call self%build()
    Hij = pack_matrix(Self%Hij,Self%Ns,Self%Nspin)
  end function Hij_get_Hij_pack_matrix





  !+------------------------------------------------------------------+
  !PURPOSE:  Pack graph matrix Hloc component into a Hloc matrix
  !+------------------------------------------------------------------+  
  function Hij_get_Hloc_pack_matrix(self) result(Hloc)
    class(hij_matrix),intent(inout)       :: self
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hloc
#else
    real(8),dimension(:,:),allocatable    :: Hloc
#endif
    integer                               :: Nlso
    !
    if(.not.Self%status)stop "Hij_get_Hloc_ispin_matrix ERROR: Self not allocated"
    !
    Nlso = self%Nspin*Self%Ns
    if(allocated(Hloc))deallocate(Hloc)
    allocate(Hloc(Nlso,Nlso))
    if(.not.Self%built)call self%build()
    Hloc = pack_matrix(Self%Hloc,Self%Ns,Self%Nspin)
  end function Hij_get_Hloc_pack_matrix

























  !##################################################################
  !                   COMPUTATIOAL ROUTINES
  !##################################################################

  !Set/Get/Delete link: INTERNAL USE
  subroutine set_link(self,siteI,siteJ,orbI,orbJ,spin,Hval)
    type(link),intent(inout) :: self
    integer,intent(in)       :: siteI,siteJ
    integer,intent(in)       :: orbI,orbJ
    integer,intent(in)       :: spin
#ifdef _CMPLX
    complex(8),intent(in)    :: Hval
#else
    real(8),intent(in)       :: Hval
#endif
    self%siteI = siteI
    self%siteJ = siteJ
    self%orbI  = orbI
    self%orbJ  = orbJ
    self%spin  = spin
    self%Hval  = Hval
  end subroutine set_link
  !

  subroutine get_link(self,siteI,siteJ,orbI,orbJ,spin,Hval)
    type(link),intent(inout) :: self
    integer,intent(inout)    :: siteI,siteJ
    integer,intent(inout)    :: orbI,orbJ
    integer,intent(inout)    :: spin
#ifdef _CMPLX
    complex(8),intent(inout) :: Hval
#else
    real(8),intent(inout)    :: Hval
#endif
    siteI = self%siteI
    siteJ = self%siteJ
    orbI  = self%orbI
    orbJ  = self%orbJ
    spin  = self%spin
    Hval  = self%Hval
  end subroutine get_link
  !
  subroutine del_link(self)
    type(link),intent(inout) :: self
    self%siteI = 0
    self%siteJ = 0
    self%orbI  = 0
    self%orbJ  = 0
    self%spin  = 0
    self%Hval  = zero
    self%next  => null()
  end subroutine del_link
  !
  subroutine print_link(self,unit)
    type(link),intent(in) :: self
    integer               :: unit
#ifdef _CMPLX
    write(unit,"(5I6,2F12.4)")&
         self%siteI,self%siteJ,&
         self%orbI,self%orbJ,&
         self%spin,dreal(self%Hval),dimag(self%Hval)
#else
    write(unit,"(5I6,F12.4)")&
         self%siteI,self%siteJ,&
         self%orbI,self%orbJ,&
         self%spin,self%Hval
#endif
  end subroutine print_link
  !















  function herm_check(A,err) result(bool)
#ifdef _CMPLX
    complex(8),intent(in) :: A(:,:)
#else
    real(8),intent(in) :: A(:,:)
#endif
    real(8),optional      :: err
    real(8)               :: err_
    logical               :: bool
    err_ = 0d0;if(present(err))err_=err
    bool = .true.
#ifdef _CMPLX
    if( any( abs(A-conjg(transpose(A))) > err_ ) )bool=.false.
#else
    if( any( abs(A-transpose(A)) > err_ ) )bool=.false.
#endif
  end function herm_check


  function pack_matrix(Hin,Ns,Nspin) result(Hout)
#ifdef _CMPLX
    complex(8),dimension(Nspin,Ns,Ns)       :: Hin
    complex(8),dimension(Nspin*Ns,Nspin*Ns) :: Hout
#else
    real(8),dimension(Nspin,Ns,Ns)          :: Hin
    real(8),dimension(Nspin*Ns,Nspin*Ns)    :: Hout
#endif
    integer                                 :: Ns,Nspin
    integer                                 :: ispin,is,js
    integer                                 :: io,jo
    do ispin=1,Nspin
       do is=1,Ns
          do js=1,Ns
             io = is + (ispin-1)*Ns
             jo = js + (ispin-1)*Ns
             Hout(io,jo) = Hin(ispin,is,js)
          enddo
       enddo
    enddo
  end function pack_matrix
  !
  function unpack_matrix(Hin,Ns,Nspin) result(Hout)
#ifdef _CMPLX
    complex(8),dimension(Nspin*Ns,Nspin*Ns) :: Hin
    complex(8),dimension(Nspin,Ns,Ns)       :: Hout
#else
    real(8),dimension(Nspin*Ns,Nspin*Ns)    :: Hin
    real(8),dimension(Nspin,Ns,Ns)          :: Hout
#endif
    integer                                 :: Ns,Nspin
    integer                                 :: ispin,is,js
    integer                                 :: io,jo
    do ispin=1,Nspin
       do is=1,Ns
          do js=1,Ns
             io = is + (ispin-1)*Ns
             jo = js + (ispin-1)*Ns
             Hout(ispin,is,js) = Hin(io,jo)
          enddo
       enddo
    enddo
  end function unpack_matrix






END MODULE MATRIX_GRAPH




















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
program testGRAPH_MATRIX
  USE MATRIX_GRAPH
  USE SCIFOR
  implicit none

  type(hij_matrix)                        :: H
#ifdef _CMPLX
  complex(8)                              :: eloc=dcmplx(0.2d0,0d0)
  complex(8)                              :: v=dcmplx(0d0,0.1d0)
  complex(8)                              :: ts=dcmplx(-1d0,0d0)
  complex(8),dimension(:,:,:),allocatable :: Hij,Hloc
  complex(8),dimension(:,:),allocatable   :: Hs
#else
  real(8)                                 :: eloc=0.2d0
  real(8)                                 :: v=0.1d0
  real(8)                                 :: ts=-1d0
  real(8),dimension(:,:,:),allocatable    :: Hij,Hloc
  real(8),dimension(:,:),allocatable      :: Hs
#endif
  integer                                 :: iorb

  H = hij_matrix([4])


  !> ionic potential
  call H%add_link(1,1,1,1,1, eloc)
  call H%add_link(2,2,1,1,1,-eloc)
  call H%add_link(3,3,1,1,1, eloc)
  call H%add_link(4,4,1,1,1,-eloc)
  !> hoppings
  call H%add_link(1,2,1,1,1,ts)
  call H%add_link(1,4,1,1,1,ts)
  call H%add_link(2,1,1,1,1,ts)
  call H%add_link(2,3,1,1,1,ts)
  call H%add_link(3,2,1,1,1,ts)
  call H%add_link(3,4,1,1,1,ts)
  call H%add_link(4,1,1,1,1,ts)
  call H%add_link(4,3,1,1,1,ts)
  !
  print*,"H.info"
  call H%info()
  print*,"H.write"
  call H%write()

  print*,"H.build (silent)"
  call H%build()


  print*,"H.get_Hij"
  Hij = H%get_hij()
  call print_matrix(Hij(1,:,:))
  print*,""

  print*,"H.get_Hloc"
  Hloc = H%get_hloc()
  call print_matrix(Hloc(1,:,:))
  print*,""

  print*,"H.get_Hspin=1"  
  Hs = H%get_hij(1)
  call print_matrix(Hs)
  print*,""



  print*,"H.pack_Hij"  
  Hs = H%pack_hij()
  call print_matrix(Hs)
  print*,""


  call H%free()
  print*,H%status









  H = hij_matrix([2,2])


  !> crystal field potential
  call H%add_link(1,1,1,1,1, eloc)
  call H%add_link(1,1,2,2,1,-eloc)
  call H%add_link(2,2,1,1,1, eloc)
  call H%add_link(2,2,2,2,1,-eloc)

  call H%add_link(1,1,1,2,1, v)
  call H%add_link(1,1,2,1,1, conjg(v))
  call H%add_link(2,2,1,2,1, v)
  call H%add_link(2,2,2,1,1, conjg(v))

  !> hoppings
  call H%add_link(1,2,1,1,1,ts)
  call H%add_link(2,1,1,1,1,ts)
  !
  call H%add_link(1,2,2,2,1,-ts)
  call H%add_link(2,1,2,2,1,-ts)
  !
  print*,"H.info"
  call H%info()
  print*,"H.write"
  call H%write()




  print*,"H.build (silent)"
  call H%build()


  print*,"H.get_Hij"
  Hij = H%get_hij()
  call print_matrix(Hij(1,:,:))
  print*,""

  print*,"H.get_Hloc"
  Hloc = H%get_hloc()
  call print_matrix(Hloc(1,:,:))
  print*,""

  print*,"H.get_Hspin=1"  
  Hs = H%get_hij(1)
  call print_matrix(Hs)
  print*,""



  print*,"H.pack_Hij"  
  Hs = H%pack_hij()
  call print_matrix(Hs)
  print*,""


  call H%free()
  print*,H%status




end program testGRAPH_MATRIX
#endif
