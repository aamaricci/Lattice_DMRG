module MATRIX_GRAPH
  USE SCIFOR, only: free_unit,reg,file_length,assert_shape,print_matrix,eye
  USE AUX_FUNCS, only: show_fmt,append
  implicit none
  private


#ifdef _CMPLX
  complex(8) :: zero=cmplx(0d0,0d0)
#else
  real(8)    :: zero=0d0
#endif



  !This represents the Operator T_\delta between two sites Ri and Rj, possibly the same. 
  type link
     integer                :: Nso=0           !Norb*max(Npsin,2)
     integer                :: siteI=0,siteJ=0 !Ri and Rj <Ri|T_ij|Rj>
#ifdef _CMPLX
     complex(8),allocatable :: Hval(:,:)   !The operator CMPLX
#else
     real(8),allocatable    :: Hval(:,:)   !The operator DBLE
#endif
     type(link),pointer     :: next=>null() !link to next box (chain)
  end type link


  type hij_matrix
     type(link),pointer                        :: root=>null() !head/root of the list
     type(link)                                :: homo     !homogeneous model
     character(len=:),allocatable              :: file     !Name of the output file 
     integer                                   :: Size=0   !Number of hopping elements
     integer                                   :: Nsites=0 !Number of sites per dimesion
     integer                                   :: Ndim=0   !Number of lattice dimensions
     integer                                   :: Nlat=0   !Number of sites
     integer                                   :: Norb=0   !Number of orbitals
     integer                                   :: Nspin=0  !Number of spin channels (1 or 2)
     integer                                   :: Nso=0    !Number of degrees of freedom
#ifdef _CMPLX
     complex(8),allocatable,dimension(:,:,:,:) :: Hij      !The non-local Hamiltonian H(Ri,Rj)-H(Ri,Ri)
     complex(8),allocatable,dimension(:,:,:)   :: Hloc     !The local Hamiltonian     H(Ri,Ri) 
#else
     real(8),allocatable,dimension(:,:,:,:)    :: Hij      !The H(Ri,Rj)_ab^ss Hamiltonian
     real(8),allocatable,dimension(:,:,:)      :: Hloc     !The local part of H(Ri,Ri)_ab^ss
#endif
     logical                                   :: ihomo  =.false. !Homo model is set
     logical                                   :: built  =.false. !Hij is built
     logical                                   :: status =.false. !Allocated
   contains
     procedure,pass :: init       => hij_init_matrix
     procedure,pass :: free       => hij_free_matrix
     procedure,pass :: info       => hij_info_matrix
     procedure,pass :: write      => hij_write_matrix
     generic        :: add_link   => hij_add_link_scalar_matrix,hij_add_link_array_matrix
     generic        :: model1d    => Hij_1d_model_scalar_matrix,Hij_1d_model_array_matrix
     procedure,pass :: build      => hij_build_matrix
     procedure,pass :: get        => Hij_get_matrix
     procedure,pass :: get_Hij    => Hij_get_Hij_matrix
     procedure,pass :: get_Hloc   => Hij_get_Hloc_matrix
     procedure,pass :: pack       => Hij_pack_matrix
     procedure,pass :: pack_Hij   => Hij_pack_Hij_matrix
     procedure,pass :: pack_Hloc  => Hij_pack_Hloc_matrix
     procedure,pass :: hij_add_link_scalar_matrix,hij_add_link_array_matrix
     procedure,pass :: Hij_1d_model_scalar_matrix,Hij_1d_model_array_matrix
  end type hij_matrix


  interface hij_matrix
     module procedure :: hij_constructor_matrix
  end interface hij_matrix


  public :: hij_matrix


contains



  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the Hij structure: constructor
  !+------------------------------------------------------------------+
  subroutine hij_init_matrix(self,Nsites,Norb,Nspin)
    class(hij_matrix),intent(inout) :: self
    integer,intent(in)              :: Nsites
    integer,intent(in),optional     :: Norb,Nspin
    integer                         :: Norb_,Nspin_,Nso,Nlat,Ndim
    !
    Norb_  = 1; if(present(Norb))Norb_=Norb
    Nspin_ = 2; if(present(Nspin))Nspin_=Nspin
    Nso    = Norb_*Nspin_
    Nlat   = Nsites!product(Nsites)
    Ndim   = 0!size(Nsites)
    !
    if(self%status)call self%free()
    !
    allocate(Self%root)
    Self%root%next => null()
    !
    Self%file   = ""
    Self%size   = 0
    !allocate(Self%Nsites, source=Nsites)
    Self%Nsites = Nsites
    Self%Nlat   = Nlat
    Self%Ndim   = Ndim
    Self%Norb   = Norb_
    Self%Nspin  = Nspin_
    Self%Nso    = Nso
    !
    allocate( Self%Hij(Nlat,Nlat,Nso,Nso) )
    Self%Hij    = zero
    allocate( Self%Hloc(Nlat,Nso,Nso) ) !local part
    Self%Hloc   = zero
    Self%built  =.false.
    Self%status =.true.
    return
  end subroutine hij_init_matrix
  !

  function hij_constructor_matrix(Nsites,Norb,Nspin) result(self)
    type(hij_matrix)            :: self
    integer,intent(in)          :: Nsites
    integer,intent(in),optional :: Norb,Nspin
    integer                     :: Norb_,Nspin_,Nso
    !
    Norb_  = 1; if(present(Norb))Norb_=Norb
    Nspin_ = 2; if(present(Nspin))Nspin_=Nspin
    call self%init(Nsites,Norb_,Nspin_)
  end function hij_constructor_matrix
  !



  !+------------------------------------------------------------------+
  !PURPOSE:  free the graph matrix object
  !+------------------------------------------------------------------+
  subroutine hij_free_matrix(self)
    class(hij_matrix),intent(inout) :: self
    type(link),pointer              :: p,c
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
    !deallocate(Self%Nsites)
    deallocate(Self%file)
    deallocate(Self%Hij)
    deallocate(Self%Hloc)
    Self%size   = 0
    Self%Nso    = 0
    Self%Nsites = 0
    Self%Nlat   = 0
    Self%Ndim   = 0
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
    write(*,"(A,I8)")  "Nlat    =",Self%Nlat
    write(*,"(A,I8)")  "Dim     =",Self%Ndim
    write(*,"(A,I8)")  "Nspin   =",Self%Nspin
    write(*,"(A,I8)")  "Norb    =",Self%Norb
    write(*,"(A,I8)")  "Nso     =",Self%Nso
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


  !   !+------------------------------------------------------------------+
  !   !PURPOSE:  Read  graph matrix from file
  !   !+------------------------------------------------------------------+  
  !   subroutine Hij_read_matrix(self,file)
  !     class(hij_matrix),intent(inout) :: self
  !     character(len=*)                :: file
  !     integer                         :: Nsize
  !     integer                         :: unitIO
  !     integer                         :: ih,i,j,a,b,spin
  !     real(8)                         :: re,im
  !     !
  !     if(.not.Self%status)stop "Hij_read ERROR: Self not allocated"
  !     !
  !     !Read from file initial info OR reconstruct them is Header is missing
  !     open(free_unit(unitIO),file=reg(file),status="old",action="read")
  !     Nsize = file_length(reg(file))
  !     Self%file = reg(file)
  !     Self%Size = Nsize
  !     !
  !     do ih=1,Self%Size
  ! #ifdef _CMPLX
  !        read(unitIO,*)i,j,a,b,spin,re,im
  !        call self%add_link(i,j,a,b,spin,dcmplx(re,im))
  ! #else
  !        read(unitIO,*)i,j,a,b,spin,re
  !        call self%add_link(i,j,a,b,spin,re)
  ! #endif
  !     enddo
  !     close(unitIO)
  !     return
  !   end subroutine Hij_read_matrix
  !   !










  !+------------------------------------------------------------------+
  !PURPOSE:  add link to graph matrix to build the lattice structure
  !+------------------------------------------------------------------+  
  subroutine Hij_add_link_scalar_matrix(self,siteI,siteJ,t)
    class(hij_matrix),intent(inout) :: self
    integer,intent(in)              :: siteI,siteJ
#ifdef _CMPLX
    complex(8),intent(in)           :: t
    complex(8),allocatable          :: Hval(:,:)
#else
    real(8),intent(in)              :: t
    real(8),allocatable             :: Hval(:,:)
#endif
    type(link),pointer              :: p,c
    integer                         :: k
    !
    !
    if(.not.Self%status)stop "Hij_add_link_matrix ERROR: Self not allocated"
    if( siteI < 0  .OR. siteJ < 0 )&
         stop "Hij_add_link_matrix ERROR: any(siteI) or any(siteJ)  < 0"
    if( siteI>Self%Nsites .OR. siteJ>Self%Nsites )&
         stop "Hij_add_link_matrix ERROR: any(siteI) or any(siteJ)  > Nsites"
    !
    !
    allocate(Hval(self%Nso,self%Nso))
    Hval = t*eye(self%Nso)
    !
    p => Self%root
    c => p%next
    do k=1,Self%size
       if(.not.associated(c))exit
       p => c
       c => c%next
    end do
    !   
    allocate(p%next)
    call set_link(p%next, siteI, siteJ, Hval)
    c=>p%next
    !
    if(siteI/=siteJ)then
       allocate(p%next%next)
#ifdef _CMPLX
       call set_link(p%next%next, siteJ, siteI, conjg(transpose((Hval))) )
#else
       call set_link(p%next%next, siteJ, siteI, transpose((Hval)) )
#endif
       c=>p%next%next
    endif
    !
    c%next  => null()
    !
    Self%size=Self%size+2
  end subroutine Hij_add_link_scalar_matrix


  subroutine Hij_add_link_array_matrix(self,siteI,siteJ,Hval)
    class(hij_matrix),intent(inout) :: self
    integer,intent(in)              :: siteI,siteJ
#ifdef _CMPLX
    complex(8),intent(in)           :: Hval(:,:)
#else
    real(8),intent(in)              :: Hval(:,:)
#endif
    type(link),pointer              :: p,c
    integer                         :: k
    !
    !
    if(.not.Self%status)stop "Hij_add_link_matrix ERROR: Self not allocated"
    if( siteI < 0  .OR. siteJ < 0 )&
         stop "Hij_add_link_matrix ERROR: any(siteI) or any(siteJ)  < 0"
    if( siteI>Self%Nsites .OR. siteJ>Self%Nsites )&
         stop "Hij_add_link_matrix ERROR: any(siteI) or any(siteJ)  > Nsites"
    !
    call assert_shape(Hval,[Self%Nso,Self%Nso],"Hij_add_link_matrices","Hij")
    if(.not.herm_check(Hval,1d-6))stop "Hij_add_link_matrix ERROR: Hval is not hermitian"
    !
    p => Self%root
    c => p%next
    do k=1,Self%size
       if(.not.associated(c))exit
       p => c
       c => c%next
    end do
    !   
    allocate(p%next)
    call set_link(p%next, siteI, siteJ, Hval)
    c=>p%next
    !
    if(siteI/=siteJ)then
       allocate(p%next%next)
#ifdef _CMPLX
       call set_link(p%next%next, siteJ, siteI, conjg(transpose((Hval))) )
#else
       call set_link(p%next%next, siteJ, siteI, transpose((Hval)) )
#endif
       c=>p%next%next
    endif
    !
    c%next  => null()
    !
    Self%size=Self%size+2
  end subroutine Hij_add_link_array_matrix









  subroutine Hij_1d_model_scalar_matrix(self,Tx,T0,pbc)
    class(hij_matrix),intent(inout) :: self
#ifdef _CMPLX
    complex(8),intent(in)           :: Tx
    complex(8),intent(in),optional  :: T0
#else
    real(8),intent(in)              :: Tx
    real(8),intent(in),optional     :: T0
#endif
    logical,intent(in),optional     :: pbc
    logical                         :: pbc_
    integer                         :: ilat,jlat
    !
    pbc_=.false. ; if(present(pbc))pbc_=pbc
    !
    if(.not.Self%status)stop "Hij_set_model_matrix ERROR: Self not allocated"
    !
    do ilat=1,self%Nsites-1
       call self%add_link(ilat,ilat+1,Tx*eye(self%Nso))
    enddo
    if(pbc_)call self%add_link(self%Nsites,1,Tx*eye(self%Nso))
    !
    if(present(T0))then
       do ilat=1,self%Nsites
          call self%add_link(ilat,ilat,T0*eye(self%Nso))
       enddo
    endif
    !
  end subroutine Hij_1d_model_scalar_matrix


  subroutine Hij_1d_model_array_matrix(self,Tx,T0,pbc)
    class(hij_matrix),intent(inout) :: self
#ifdef _CMPLX
    complex(8),intent(in)           :: Tx(:,:)
    complex(8),intent(in),optional  :: T0(:,:)
#else
    real(8),intent(in)              :: Tx(:,:)
    real(8),intent(in),optional     :: T0(:,:)
#endif
    logical,intent(in),optional     :: pbc
    logical                         :: pbc_
    integer                         :: ilat,jlat
    !
    pbc_=.false. ; if(present(pbc))pbc_=pbc
    !
    if(.not.Self%status)stop "Hij_set_model_matrix ERROR: Self not allocated"
    !
    call assert_shape(Tx,[Self%Nso,Self%Nso],"Hij_set_model_matrices","Tx")
    if(.not.herm_check(Tx,1d-6))stop "Hij_set_model_matrix ERROR: Tx is not hermitian"
    !
    do ilat=1,self%Nsites-1
       call self%add_link(ilat,ilat+1,Tx)
    enddo
    if(pbc_)call self%add_link(self%Nsites,1,Tx)
    !
    if(present(T0))then
       call assert_shape(T0,[Self%Nso,Self%Nso],"Hij_set_model_matrices","T0")
       if(.not.herm_check(T0,1d-6))stop "Hij_set_model_matrix ERROR: T0 is not hermitian"
       do ilat=1,self%Nsites
          call self%add_link(ilat,ilat,T0)
       enddo
    endif
    !
  end subroutine Hij_1d_model_array_matrix











  !+------------------------------------------------------------------+
  !PURPOSE:  Build the actual Hij from graph matrix 
  !+------------------------------------------------------------------+  
  subroutine Hij_build_matrix(self)
    class(hij_matrix),intent(inout) :: self
    type(link),pointer              :: c           
    integer                         :: i,j,a,b,spin
    integer                         :: io,jo
#ifdef _CMPLX
    complex(8),allocatable          :: val(:,:)
#else
    real(8),allocatable             :: val(:,:)
#endif
    !
    if(.not.Self%status)stop "Hij_build ERROR: Self not allocated"
    !
    c => Self%root%next
    do
       if(.not.associated(c))exit
       call get_link(c,i,j,val)
       if(i==j)then
          Self%Hloc(i,:,:) = val
       else
          Self%Hij(i,j,:,:) = val
       endif
       if(allocated(val))deallocate(val)
       c => c%next
    enddo

    if(.not.herm_check( pack_matrix(self%Hij,Self%Nlat,Self%Nso) ,1d-6))stop "Hij_build ERROR: Hij is not Hermitian"    
    Self%built=.true.
    return
  end subroutine Hij_build_matrix
  !












  !+------------------------------------------------------------------+
  !PURPOSE:  Retrieve the  actual Hij matrix from graph matrix type
  !+------------------------------------------------------------------+  
  function Hij_get_Hij_matrix(self) result(Hij)
    class(hij_matrix),intent(inout)           :: self
#ifdef _CMPLX
    complex(8),dimension(:,:,:,:),allocatable :: Hij
#else
    real(8),dimension(:,:,:,:),allocatable    :: Hij
#endif
    !
    if(.not.Self%status)stop "Hij_get_Hij_matrix ERROR: Self not allocated"
    !
    if(.not.Self%built)call self%build()
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij,  source=Self%Hij)
  end function Hij_get_Hij_matrix



  !+------------------------------------------------------------------+
  !PURPOSE:  Retrieve the actual Hloc matrix from graph matrix type
  !+------------------------------------------------------------------+  
  function Hij_get_Hloc_matrix(self) result(Hloc)
    class(hij_matrix),intent(inout)           :: self
#ifdef _CMPLX
    complex(8),dimension(:,:,:),allocatable :: Hloc
#else
    real(8),dimension(:,:,:),allocatable    :: Hloc
#endif
    !
    if(.not.Self%status)stop "Hij_get_Hloc_matrix ERROR: Self not allocated"
    !
    if(.not.Self%built)call self%build()
    if(allocated(Hloc))deallocate(Hloc)
    allocate(Hloc, source=Self%Hloc)
  end function Hij_get_Hloc_matrix





  !+------------------------------------------------------------------+
  !PURPOSE:  Retrieve the  actual Hij matrix from graph matrix type
  !+------------------------------------------------------------------+  
  function Hij_get_matrix(self) result(H)
    class(hij_matrix),intent(inout)           :: self
#ifdef _CMPLX
    complex(8),dimension(:,:,:,:),allocatable :: H
#else
    real(8),dimension(:,:,:,:),allocatable    :: H
#endif
    integer                                   :: ilat,is,js,io,jo
    !
    if(.not.Self%status)stop "Hij_get_Hij_matrix ERROR: Self not allocated"
    !
    if(.not.Self%built)call self%build()
    if(allocated(H))deallocate(H)
    allocate(H(self%Nlat,self%Nlat,self%Nso,self%Nso))

    H = self%Hij
    do ilat=1,Self%Nlat
       H(ilat,ilat,:,:) = H(ilat,ilat,:,:) +  Self%Hloc(ilat,:,:)
    enddo
  end function Hij_get_matrix










  !+------------------------------------------------------------------+
  !PURPOSE:  Pack graph matrix Hij component into a Hij matrix
  !+------------------------------------------------------------------+  
  function Hij_pack_Hij_matrix(self) result(Hij)
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
    Nlso = self%Nlat*Self%Nso
    if(.not.Self%built)call self%build()
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij(Nlso,Nlso))
    Hij = pack_matrix(Self%Hij,Self%Nlat,Self%Nso)
  end function Hij_pack_Hij_matrix



  !+------------------------------------------------------------------+
  !PURPOSE:  Pack graph matrix Hloc component into a Hloc matrix
  !+------------------------------------------------------------------+  
  function Hij_pack_Hloc_matrix(self) result(Hloc)
    class(hij_matrix),intent(inout)       :: self
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hloc
#else
    real(8),dimension(:,:),allocatable    :: Hloc
#endif
    integer                               :: Nlso,ilat,is,js,io,jo
    !
    if(.not.Self%status)stop "Hij_get_Hloc_ispin_matrix ERROR: Self not allocated"
    !
    Nlso = self%Nlat*Self%Nso
    if(.not.Self%built)call self%build()
    if(allocated(Hloc))deallocate(Hloc)
    allocate(Hloc(Nlso,Nlso))
    do concurrent(ilat=1:Self%Nlat,is=1:Self%Nso,js=1:Self%Nso)
       io = is + (ilat-1)*Self%Nso
       jo = js + (ilat-1)*Self%Nso
       Hloc(io,jo) = Self%Hloc(ilat,is,js)
    enddo
  end function Hij_pack_Hloc_matrix



















  !+------------------------------------------------------------------+
  !PURPOSE:  Pack graph matrix Hij component into a Hij matrix
  !+------------------------------------------------------------------+  
  function Hij_pack_matrix(self) result(H)
    class(hij_matrix),intent(inout)       :: self
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: H
#else
    real(8),dimension(:,:),allocatable    :: H
#endif
    integer                               :: Nlso
    !
    if(.not.Self%status)stop "Hij_get_Hij_ispin_matrix ERROR: Self not allocated"
    !      
    Nlso = self%Nlat*Self%Nso
    if(.not.Self%built)call self%build()
    if(allocated(H))deallocate(H)
    allocate(H(Nlso,Nlso))
    H = self%pack_Hij() + self%pack_Hloc()
  end function Hij_pack_matrix










  !   !##################################################################
  !   !                   COMPUTATIOAL ROUTINES
  !   !##################################################################

  !Set/Get/Delete link: INTERNAL USE
  subroutine set_link(self,siteI,siteJ,Hval)
    type(link),intent(inout) :: self
    integer,intent(in)       :: siteI,siteJ
#ifdef _CMPLX
    complex(8),intent(in)    :: Hval(:,:)
#else
    real(8),intent(in)       :: Hval(:,:)
#endif    
    self%siteI=siteI
    self%siteJ=siteJ
    allocate(self%Hval,  source=Hval)
  end subroutine set_link
  !   !

  subroutine get_link(self,siteI,siteJ,Hval)
    type(link),intent(inout)             :: self
    integer,intent(inout)                :: siteI,siteJ
#ifdef _CMPLX
    complex(8),allocatable,intent(inout) :: Hval(:,:)
#else
    real(8),allocatable,intent(inout)    :: Hval(:,:)
#endif
    siteI = self%siteI
    siteJ = self%siteJ
    if(allocated(Hval))deallocate(Hval)
    allocate(Hval, source=self%Hval)
  end subroutine get_link
  !
  subroutine del_link(self)
    type(link),intent(inout) :: self
    if(allocated(self%Hval))deallocate(self%Hval)
    self%siteI=0
    self%siteJ=0
    self%next  => null()
  end subroutine del_link
  !
  subroutine print_link(self,unit)
    type(link),intent(in) :: self
    integer               :: unit
    write(unit,"(2I6)")self%siteI,self%siteJ
    call print_matrix(self%Hval)
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


  function pack_matrix(Hin,Nlat,Nso) result(Hout)
#ifdef _CMPLX
    complex(8),dimension(Nlat,Nlat,Nso,Nso) :: Hin
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Hout
#else
    real(8),dimension(Nlat,Nlat,Nso,Nso)    :: Hin
    real(8),dimension(Nlat*Nso,Nlat*Nso)    :: Hout
#endif
    integer                                 :: Nlat,Nso
    integer                                 :: ilat,jlat,is,js
    integer                                 :: io,jo
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,is=1:Nso,js=1:Nso)
       io = is + (ilat-1)*Nso
       jo = js + (jlat-1)*Nso
       Hout(io,jo) = Hin(ilat,jlat,is,js)
    enddo
  end function pack_matrix
  !
  function unpack_matrix(Hin,Nlat,Nso) result(Hout)
#ifdef _CMPLX
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Hin
    complex(8),dimension(Nlat,Nlat,Nso,Nso) :: Hout
#else
    real(8),dimension(Nlat*Nso,Nlat*Nso)    :: Hin
    real(8),dimension(Nlat,Nlat,Nso,Nso)    :: Hout
#endif
    integer                                 :: Nlat,Nso
    integer                                 :: ilat,jlat,is,js
    integer                                 :: io,jo
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,is=1:Nso,js=1:Nso)
       io = is + (ilat-1)*Nso
       jo = js + (jlat-1)*Nso
       Hout(ilat,jlat,is,js) = Hin(io,jo)
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

  type(hij_matrix)                          :: H
#ifdef _CMPLX
  complex(8)                                :: eloc=dcmplx(0.2d0,0d0)
  complex(8)                                :: v=dcmplx(0.1d0,0d0),Jxy=dcmplx(0.25d0,0d0)
  complex(8)                                :: ts=dcmplx(-1d0,0d0),Jz=dcmplx(1d0,0d0)
  complex(8),dimension(:,:),allocatable     :: e,t,Hvec,J
  complex(8),dimension(:,:,:,:),allocatable :: Hij
  complex(8),dimension(:,:,:),allocatable   :: Hloc
  complex(8),dimension(:,:),allocatable     :: Hs
#else
  real(8)                                   :: eloc=0.2d0
  real(8)                                   :: v=0.1d0,Jxy=0.25d0
  real(8)                                   :: ts=-1d0,Jz=1d0
  real(8),dimension(:,:),allocatable        :: e,t,Hvec,J
  real(8),dimension(:,:,:,:),allocatable    :: Hij
  real(8),dimension(:,:,:),allocatable      :: Hloc
  real(8),dimension(:,:),allocatable        :: Hs
#endif

  integer                                   :: iorb,Nspin

  !chain of 4 sites, spin=1
  print*,""
  print*," ---------------------- "
  print*," 4 sites chain, spin=1 "
  print*," ---------------------- "
  print*,""

  H = hij_matrix(4,Nspin=1)

  ! !> ionic potential
  call H%add_link(1,1, eloc)
  call H%add_link(2,2,-eloc)
  call H%add_link(3,3, eloc)
  call H%add_link(4,4,-eloc)

  !> hoppings
  call H%add_link(1,2,ts)
  call H%add_link(2,3,ts)
  call H%add_link(3,4,ts)
  call H%add_link(4,1,ts)

  print*,"H.info"
  call H%info()
  print*,"H.write"
  call H%write()

  print*,"H.build (silent)"
  call H%build()


  print*,"H.get_Hij (print 1,2)"
  Hij = H%get_hij()
  call print_matrix(Hij(1,2,:,:))
  print*,""

  print*,"H.get_Hloc (print 2,2)"
  Hloc = H%get_hloc()
  call print_matrix(Hloc(2,:,:))
  print*,""

  print*,"H.get "
  Hij = H%get()
  call print_matrix(Hij(1,2,:,:))
  call print_matrix(Hij(2,2,:,:))
  print*,""



  print*,"H.pack_Hij"  
  Hs = H%pack_hij()
  call print_matrix(Hs)
  print*,""

  print*,"H.pack_Hloc"  
  Hs = H%pack_hloc()
  call print_matrix(Hs)
  print*,""

  print*,"H.pack"  
  Hs = H%pack()
  call print_matrix(Hs)
  print*,""


  call H%free()
  print*,H%status
  print*,""
  print*,""
  print*,""





  print*,""
  print*," ---------------------- "
  print*," 4 sites chain, spin=2 "
  print*," ---------------------- "
  print*,""

  !Op are 2x2 by spin:
  allocate(e(2,2),t(2,2))
  e = diag([eloc,eloc])
  t = diag([ts,ts])

  H = hij_matrix(4,Nspin=2)

  ! !> ionic potential
  call H%add_link(1,1, e)
  call H%add_link(2,2,-eloc)
  call H%add_link(3,3, eloc)
  call H%add_link(4,4,-e)

  !> hoppings
  call H%add_link(1,2,ts)
  call H%add_link(2,3,t)
  call H%add_link(3,4,t)
  call H%add_link(4,1,ts)

  print*,"H.info"
  call H%info()
  print*,"H.write"
  call H%write()

  print*,"H.build (silent)"
  call H%build()


  print*,"H.get_Hij (print 1,2)"
  Hij = H%get_hij()
  call print_matrix(Hij(1,2,:,:))
  print*,""

  print*,"H.get_Hloc (print 2,2)"
  Hloc = H%get_hloc()
  call print_matrix(Hloc(2,:,:))
  print*,""

  print*,"H.get "
  Hij = H%get()
  call print_matrix(Hij(1,2,:,:))
  call print_matrix(Hij(2,2,:,:))
  print*,""



  print*,"H.pack_Hij"  
  Hs = H%pack_hij()
  call print_matrix(Hs)
  print*,""

  print*,"H.pack_Hloc"  
  Hs = H%pack_hloc()
  call print_matrix(Hs)
  print*,""

  print*,"H.pack"  
  Hs = H%pack()
  call print_matrix(Hs)
  print*,""


  call H%free()
  print*,H%status
  print*,""
  print*,""
  print*,""












  print*,""
  print*," ---------------------- "
  print*," 3 sites chain OBC, spin=1 "
  print*," ---------------------- "
  print*,""


  e = diag([eloc,-eloc])
#ifdef _CMPLX
  t = ts*pauli_tau_z + v*pauli_tau_y
#else
  t = ts*pauli_tau_z + v*pauli_tau_x
#endif
  H = hij_matrix(3,Norb=2,Nspin=1)

  ! !> crystal field potential
  call H%add_link(1,1, e)
  call H%add_link(2,2, e)
  call H%add_link(3,3, e)
  !> hoppings
  call H%add_link(1,2,t)
  call H%add_link(2,3,t)

  !
  print*,"H.info"
  call H%info()
  print*,"H.write"
  call H%write()




  print*,"H.build (silent)"
  call H%build()
  print*,"H.pack_Hij"  
  Hs = H%pack_hij()
  call print_matrix(Hs)
  print*,""

  print*,"H.pack_Hloc"  
  Hs = H%pack_hloc()
  call print_matrix(Hs)
  print*,""

  print*,"H.pack"  
  Hs = H%pack()
  call print_matrix(Hs)
  print*,""


  call H%free()
  print*,H%status














  print*,""
  print*," ---------------------- "
  print*," 4 sites chain, spin=1 "
  print*," ---------------------- "
  print*,""

  H = hij_matrix(6,Norb=1,Nspin=1)

  call H%model1d(Tx=ts,pbc=.true.)
  ! !> ionic potential
  call H%add_link(1,1, eloc)
  call H%add_link(2,2,-eloc)
  call H%add_link(3,3, eloc)
  call H%add_link(4,4,-eloc)


  print*,"H.info"
  call H%info()
  print*,"H.write"
  call H%write()

  print*,"H.build (silent)"
  call H%build()

  print*,"H.pack"  
  Hs = H%pack()
  call print_matrix(Hs)
  print*,""


  call H%free()
  print*,H%status
  print*,""
  print*,""
  print*,""









  print*,""
  print*," ------------------------- "
  print*," 4sites spin-chain Nspin=2 "
  print*," two ways: NORMAL "
  print*," ------------------------- "
  print*,""

  !Op are 2x2 by spin:
  Nspin=2
  allocate(Hvec(Nspin,Nspin))
  allocate(J(Nspin,Nspin))

  Hvec = diag([0d0,0.01d0])
  J    = diag([Jxy,Jz/2d0])


  H = hij_matrix(4,Nspin=2)

  ! !> ionic potential
  call H%add_link(1,1, Hvec)
  call H%add_link(2,2, Hvec)
  call H%add_link(3,3, Hvec)
  call H%add_link(4,4, Hvec)

  !> hoppings
  call H%add_link(1,2,J)
  call H%add_link(2,3,J)
  call H%add_link(3,4,J)

  print*,"H.info"
  call H%info()
  print*,"H.write"
  call H%write()

  print*,"H.build (silent)"
  call H%build()


  print*,"H.get_Hij (print 1,2)"
  Hij = H%get_hij()
  call print_matrix(Hij(1,2,:,:))
  print*,""

  print*,"H.get_Hloc (print 2,2)"
  Hloc = H%get_hloc()
  call print_matrix(Hloc(2,:,:))
  print*,""


  print*,"H.pack_Hij"  
  Hs = H%pack_hij()
  call print_matrix(Hs)
  print*,""

  print*,"H.pack_Hloc"  
  Hs = H%pack_hloc()
  call print_matrix(Hs)
  print*,""

  print*,"H.pack"  
  Hs = H%pack()
  call print_matrix(Hs)
  print*,""


  call H%free()
  print*,H%status
  print*,""
  print*,""
  print*,""




  print*,""
  print*," ------------------------- "
  print*," 4sites spin-chain Nspin=2 "
  print*," two ways: MODEL1d "
  print*," ------------------------- "
  print*,""

  H = hij_matrix(4,Nspin=2) 
  call H%model1d(Tx=J,T0=Hvec,pbc=.false.)
  print*,"H.pack"  
  Hs = H%pack()
  call print_matrix(Hs)
  print*,""


  call H%free()
  print*,H%status
  print*,""
  print*,""
  print*,""



end program testGRAPH_MATRIX
#endif
