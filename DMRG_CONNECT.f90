MODULE DMRG_CONNECT
  USE VARS_GLOBAL
  implicit none
  private


  public :: enlarge_block
  public :: connect_fermion_blocks
  public :: connect_spin_blocks


contains

  
  !-----------------------------------------------------------------!
  ! Purpose: enlarge a given BLOCK "SELF" growing it left/right with
  ! a SITE "DOT" (specified in the init)
  !-----------------------------------------------------------------!
  subroutine enlarge_block(self,dot,grow)
    type(block)                  :: self
    type(site)                   :: dot
    character(len=*),optional    :: grow
    character(len=16)            :: grow_
    character(len=:),allocatable :: key,dtype,otype
    type(tbasis)                 :: self_basis,dot_basis,enl_basis
    type(sparse_matrix)          :: Hb,Hd,H2,eO
    integer                      :: i
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    call start_timer("Enlarge blocks "//str(grow_))
    !
    if(.not.self%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: Missing self.H operator in the list"
    if(.not.dot%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: Missing dot.H operator in the list"
    !
    dtype=dot%type()
    if(dtype/=self%type())&
         stop "Enlarge_Block ERROR: Dot.Type != Self.Type"
    !    
    !> Update Hamiltonian:
    select case(str(grow_))
    case ("left","l")
       Hb = self%operators%op("H").x.id(dot%dim)
       Hd = id(self%dim).x.dot%operators%op("H")
       select case(dtype)
       case default;stop "Enlarge_Block ERROR: wrong dot.Type"
       case ("spin","s")
          H2 = connect_spin_blocks(self,as_block(dot))
       case ("fermion","f,","electron","e")
          H2 = connect_fermion_blocks(self,as_block(dot))
       end select
    case ("right","r")
       Hb = id(dot%dim).x.self%operators%op("H")
       Hd = dot%operators%op("H").x.id(self%dim)
       select case(dtype)
       case default;stop "Enlarge_Block ERROR: wrong dot.Type"
       case ("spin","s")
          H2 = connect_spin_blocks(as_block(dot),self)
       case ("fermion","f,","electron","e")
          H2 = connect_fermion_blocks(as_block(dot),self)
       end select
    end select
    call self%put_op("H", Hb +  Hd + H2, type="bosonic")
    !
    !> Update all the other operators in the list:
    do i=1,size(self%operators)
       key   = str(self%operators%key(index=i))
       otype = str(self%operators%type(index=i))
       if(key=="H")cycle
       !
       !Bosonic operators:
       !O_L -> I_L.x.O_d  | O_R -> O_d.x.I_R
       !Fermionic operators:
       !O_L -> P_L.x.O_d  | O_R -> O_d.x.I_R
       !Sign operator:
       !P_L -> P_L.x.P_d  | P_R -> P_d.x.P_R
       select case(str(grow_))
       case ("left","l")
          select case(to_lower(str(otype)))
          case default;stop "Enlarge_BLock ERROR: wrong self.operators.type L: !\in['Bosonic','Fermionic','Sign']"
          case('b','bose','bosonic')
             eO = Id(self%dim)
          case('f','fermi','fermionic')
             eO = self%operators%op(key="P")
          case('s','sign')
             eO = self%operators%op(key="P")
          end select
          call self%put_op(str(key), eO.x.dot%operators%op(str(key)), type=otype)
       case ("right","r")
          select case(to_lower(str(otype)))
          case default;stop "Enlarge_BLock ERROR: wrong self.operators.type R: !\in['Bosonic','Fermionic','Sign']"
          case('b','bose','bosonic')
             eO = Id(self%dim)
          case('f','fermi','fermionic')
             eO = Id(self%dim)
          case('s','sign')
             eO = self%operators%op(key="P")
          end select
          call self%put_op(str(key), dot%operators%op(str(key)).x.eO, type=otype )
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
    select case(str(grow_))
    case ("left","l")
       enl_basis = (self_basis.o.dot_basis)
       call self%set_basis( basis=enl_basis )
    case ("right","r")
       enl_basis = (dot_basis.o.self_basis)
       call self%set_basis( basis=enl_basis )
    end select
    !
#ifdef _DEBUG
    call self%show(file="Enl"//str(grow_)//"_"//str(self%length)//".dat")
#endif
    !
    if(.not.self%is_valid())then
       write(LOGfile,*)"dmrg_step error: enlarged_block "//str(grow_)// "is not a valid block"
       stop
    endif
    !
    !Free the memory:
    call Hb%free()
    call Hd%free()
    call H2%free()
    call self_basis%free()
    call dot_basis%free()
    call enl_basis%free()
    !
    call stop_timer()
    !
  end subroutine enlarge_block


  !H_lr = \sum_{a}h_aa*(C^+_{left,a}@P_left) x C_{right,a}] + H.c.
  function connect_fermion_blocks(left,right,states) result(H2)
    type(block)                               :: left
    type(block)                               :: right
    integer,dimension(:),optional             :: states
    type(sparse_matrix),dimension(Nspin*Norb) :: Cl,Cr
    type(sparse_matrix)                       :: P,A
    type(sparse_matrix)                       :: H2
    integer                                   :: ispin,iorb,jorb,io,jo
    integer,dimension(2)                      :: Hdims,dleft,dright
    character(len=:),allocatable              :: key
    complex(8),dimension(:,:),allocatable     :: Hij
    !
    !Hij is shared:
    !Hij = Hmodel(left,right)
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij, source=HopH)
    !
    !> Get H2 dimensions:
    dleft = shape(left%operators)
    dright= shape(right%operators)
    Hdims = dleft*dright
    if(present(states))Hdims = [size(states),size(states)]
    call H2%init(Hdims(1),Hdims(2))
    !
    !FERMION SPECIFIC:
    !>Retrieve operators:
    do ispin=1,Nspin
       do iorb=1,Norb
          key = "C"//left%okey(iorb,ispin)
          io = iorb + (ispin-1)*Norb
          Cl(io) = left%operators%op(key)
          Cr(io) = right%operators%op(key)
       enddo
    enddo
    !
    !
    !>Build H2:
    P = left%operators%op("P")  !always acts on the Left Block
    do io=1,Nspin*Norb
       do jo=1,Nspin*Norb
          if(Hij(io,jo)==0d0)cycle
          if(present(states))then
             H2 = H2 &
                  + Hij(io,jo)*sp_kron(matmul(Cl(io)%dgr(),P),Cr(jo),states) &
                  + Hij(io,jo)*sp_kron(matmul(P,Cl(io)),Cr(jo)%dgr(),states)
          else
             H2 = H2 &
                  + Hij(io,jo)*(matmul(Cl(io)%dgr(),P).x.Cr(jo)) &
                  + Hij(io,jo)*(matmul(P,Cl(io)).x.Cr(jo)%dgr())
          endif
       enddo
    enddo
    !
    !> free memory
    call P%free
    do io=1,Nspin*Norb
       call Cl(io)%free
       call Cr(io)%free
    enddo
  end function connect_fermion_blocks




  function connect_spin_blocks(left,right,states) result(H2)
    type(block)                           :: left
    type(block)                           :: right
    integer,dimension(:),optional         :: states
    type(sparse_matrix)                   :: Sl(Nspin)![Sz,Sp]
    type(sparse_matrix)                   :: Sr(Nspin)![Sz,Sp]
    type(sparse_matrix)                   :: H2
    integer,dimension(2)                  :: Hdims
    integer                               :: ispin
    complex(8),dimension(:,:),allocatable :: Hij
    !
    !Hij is shared:
    !Hij = Hmodel(left,right)
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij, source=HopH)    
    !
    !> Get H2 dimensions:
    Hdims = shape(left%operators)*shape(right%operators)
    if(present(states))Hdims = [size(states),size(states)]
    call H2%init(Hdims(1),Hdims(2))
    !
    !>Retrieve operators:
    do ispin=1,Nspin
       Sl(ispin) = left%operators%op("S"//left%okey(0,ispin))
       Sr(ispin) = right%operators%op("S"//right%okey(0,ispin))
    enddo
    !
    !>Build H2:
    if(present(states))then
       H2 = H2 &
            + Hij(1,1)*sp_kron(Sl(1),Sr(1),states)       &
            + Hij(2,2)*sp_kron(Sl(2),Sr(2)%dgr(),states) &
            + Hij(2,2)*sp_kron(Sl(2)%dgr(),Sr(2),states)
    else
       H2 = H2 &
            + Hij(1,1)*(Sl(1).x.Sr(1))       &
            + Hij(2,2)*(Sl(2).x.Sr(2)%dgr()) &
            + Hij(2,2)*(Sl(2)%dgr().x.Sr(2))
    endif
    !
    !> Free memory
    do ispin=1,Nspin
       call Sl(ispin)%free
       call Sr(ispin)%free
    enddo
  end function connect_spin_blocks




END MODULE DMRG_CONNECT
