MODULE DMRG_CONNECT
  USE VARS_GLOBAL
  implicit none
  private


  public :: enlarge_block
  ! public :: enlarge_block_old
  public :: connect_fermion_blocks
  public :: connect_spin_blocks


contains


  !-----------------------------------------------------------------!
  ! Purpose: enlarge a given BLOCK "SELF" growing it left/right with
  ! a SITE "DOT" (specified in the init)
  !-----------------------------------------------------------------!
  subroutine enlarge_block(self)
    type(block)                  :: self
    type(site)                   :: dot    
    character(len=:),allocatable :: tag,key,dtype,otype,type_error_msg,tag_error_msg
    type(tbasis)                 :: self_basis,dot_basis,enl_basis
    type(sparse_matrix)          :: Hb,Hd,H2,eO
    integer                      :: i,dot_index,L,n
    !
    tag= to_lower(self%tag(1))  !tag is a single character = l,s,r,e....
    L  = self%length
    !
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: ENLARGE block "//tag
#endif
    !
    tag_error_msg="Enlarge_block ERROR: self.tag !=[{left,l,sys,s},{right,r,env,e}"
    !
    call start_timer("Enlarge blocks "//tag)
    !1. get the dot at the correct site
    !Get index of the dot to pick
    !OBC: = = = = =  *   x  - - - - - -
    !                |   |            
    !     1 . . . L L+1 L+2  . . .   2L+2
    !                |                | 
    !PBC: = = = = =  *   -  - - - - - x
    if(PBCdmrg)then
       select case(tag)
       case default;stop 
       case ("l","s"); dot_index = L+1
       case ("r","e"); dot_index = 2*L+2
       end select
    else
       select case(tag)
       case default;stop tag_error_msg
       case ("l","s"); dot_index = L+1
       case ("r","e"); dot_index = L+2
       end select
    endif
    !
    dot  = dots(dot_index)
    dtype= dot%type()
    !
    if(.not.self%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: Missing self.H operator in the list"
    if(.not.dot%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: Missing dot.H operator in the list"
    if(dtype/=self%type())&
         stop "Enlarge_Block ERROR: Dot.Type != Self.Type"
    !
    !> Update Hamiltonian:
    if(PBCdmrg)then
       ! ===== * OR ----- x
       select case(tag)
       case default;stop tag_error_msg
       case ("l","s","r","e")
          Hb = self%operators%op("H").x.id(dot%dim) !H_b x 1_d
          Hd = id(self%dim).x.dot%operators%op("H") !1_b x H_d
          select case(to_lower(dtype(1:1)))         !H_bd
          case default;stop "Enlarge_Block ERROR: wrong dot.Type"
          case ("s"); H2 = connect_spin_blocks(self,as_block(dot))
          case ("f"); H2 = connect_fermion_blocks(self,as_block(dot))
             ! case ("b"); H2 = connect_boson_blocks(self,as_block(dot))
          end select
       end select
    else
       ! ===== * OR x -----
       select case(tag)
       case default;stop tag_error_msg
       case ("l","s")
          Hb = self%operators%op("H").x.id(dot%dim) !H_b x 1_d
          Hd = id(self%dim).x.dot%operators%op("H") !1_b x H_d
          select case(to_lower(dtype(1:1)))         !H_bd
          case default;stop "Enlarge_Block ERROR: wrong dot.Type"
          case ("s"); H2 = connect_spin_blocks(self,as_block(dot))
          case ("f"); H2 = connect_fermion_blocks(self,as_block(dot))
             ! case ("b"); H2 = connect_boson_blocks(self,as_block(dot))
          end select
       case ("r","e")
          Hb = id(dot%dim).x.self%operators%op("H") !1_d x H_b
          Hd = dot%operators%op("H").x.id(self%dim) !H_d x 1_b
          select case(to_lower(dtype(1:1)))         !H_bd
          case default;stop "Enlarge_Block ERROR: wrong dot.Type"
          case ("s"); H2 = connect_spin_blocks(as_block(dot),self)
          case ("f"); H2 = connect_fermion_blocks(as_block(dot),self)
             ! case ("b"); H2 = connect_boson_blocks(as_block(dot),self)
          end select
       end select
    endif
    !
    call self%put_op("H", Hb +  Hd + H2, type="bosonic")
    !
    !
    !> Update all the other operators in the list:
    do i=1,size(self%operators)
       key   = str(self%operators%key(index=i))
       otype = to_lower(str(self%operators%type(index=i)))
       type_error_msg  = "Enlarge_BLock ERROR: self.operators.type("//&
            str(i)//") !\in['Bosonic','Fermionic','Sign']"
       !
       if(key=="H")cycle
       !
       !Note that indices have different meaning for OBC and PBC
       if(PBCdmrg)then
          !                       ANY BLOCK
          !b_Left/Right    [_L]====*   ;          ====[_R]*
          !Spin  : O_bL <- O_bL.x.I_d  ; O_bR <- I_b .x.O_d
          !Bose  : O_bL <- O_bL.x.I_d  ; O_bR <- I_b .x.O_d
          !Fermi : O_bL <- O_bL.x.I_d  ; O_bR <- P_bR.x.O_d
          !Psign : P_bL <- P_bL.x.P_d  ; P_bR <- P_bR.x.P_d
          select case(to_lower(reverse(key,1))) !Op key have the form "O_{orb,spin,site}_{L,R}"
          case default;stop "Enlarge_Block ERROR: PBC and Op.key does not end in _L or _R"
          case("l")
             select case(tag)
             case default;stop tag_error_msg
             case ("l","s","r""e") !Block left == right
                select case(otype(1:1))
                case default;stop type_error_msg
                case('s'); eO = Id(dot%dim)
                case('f'); eO = Id(dot%dim)
                case('p'); eO = dot%operators%op(str(key))!"P_L"
                end select
                call self%put_op(str(key), self%operators%op(str(key)).x.eO, type=otype)
             end select
          case("r")
             select case(tag)
             case default;stop tag_error_msg
             case ("l","s","r""e") !Block left == right
                select case(otype(1:1))
                case default;stop type_error_msg
                case('s'); eO = Id(self%dim)
                case('f'); eO = self%operators%op(str(key))!"P_R"
                case('p'); eO = self%operators%op(str(key))!"P_R"
                end select
                call self%put_op(str(key), eO.x.dot%operators%op(str(key)), type=otype)
             end select
          end select
       else
          !           LEFT BLOCK     |   RIGHT BLOCK  
          !Bose  : O_L <- I_L.x.O_d  | O_R <- O_d.x.I_R
          !Fermi : O_L <- P_L.x.O_d  | O_R <- O_d.x.I_R
          !Sign  : P_L <- P_L.x.P_d  | P_R <- P_d.x.P_R
          !                ===== *            x -----
          select case(tag)
          case default;stop tag_error_msg
          case ("l","s")      
             select case(otype(1:1))
             case default;stop type_error_msg
             case('b'); eO = Id(self%dim)
             case('f'); eO = self%operators%op(key="P")
             case('s'); eO = self%operators%op(key="P")
             end select
             call self%put_op(str(key), eO.x.dot%operators%op(str(key)), type=otype)
          case ("r","e") 
             select case(otype(1:1))
             case default;stop type_error_msg
             case('b'); eO = Id(self%dim)
             case('f'); eO = Id(self%dim)
             case('s'); eO = self%operators%op(key="P")
             end select
             call self%put_op(str(key), dot%operators%op(str(key)).x.eO, type=otype )
          end select
       endif
       !
    enddo
    !
    !> Enlarge dimensions
    self%length = self%length + 1
    self%Dim    = self%Dim*dot%Dim
    !
    !> Enlarge the basis states
    call self%get_basis(self_basis)
    call dot%get_basis(dot_basis)
    !
    if(PBCdmrg)then
       select case(tag)
       case default;stop tag_error_msg
       case ("l","s")
          enl_basis = (self_basis.o.dot_basis)
          call self%set_basis( basis=enl_basis )
       case ("r","e")
          enl_basis = (self_basis.o.dot_basis)
          call self%set_basis( basis=enl_basis )
       end select
    else
       select case(tag)
       case default;stop tag_error_msg
       case ("l","s")
          enl_basis = (self_basis.o.dot_basis)
          call self%set_basis( basis=enl_basis )
       case ("r","e")
          enl_basis = (dot_basis.o.self_basis)
          call self%set_basis( basis=enl_basis )
       end select
    endif
    !
    if(.not.self%is_valid())then
       write(LOGfile,*)"dmrg_step error: enlarged_block "//tag// "is not a valid block"
       stop
    endif
    !
    !Free the memory:
    call Hb%free()
    call Hd%free()
    call H2%free()
    call eO%free()
    call self_basis%free()
    call dot_basis%free()
    call enl_basis%free()
    !
    call stop_timer()
    !
  end subroutine enlarge_block

















  function connect_spin_blocks(left,right,states) result(H2)
    type(block)                           :: left
    type(block)                           :: right
    integer,dimension(:),optional         :: states
    ! type(sparse_matrix),allocatable       :: Sl(:,:)![1-2,Nspin];Nspin=[1.Sz,2.Sp]
    ! type(sparse_matrix),allocatable       :: Sr(:,:)![1-2,Nspin];Nspin=[1.Sz,2.Sp]
    type(sparse_matrix),allocatable       :: SLz(:),SLp(:)
    type(sparse_matrix),allocatable       :: SRz(:),SRp(:)
    type(sparse_matrix)                   :: H2
    integer,dimension(2)                  :: Hdims
    integer                               :: ibc,ispin,bcdim
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hij
    complex(8)                            :: Jz,Jxy
#else
    real(8),dimension(:,:),allocatable    :: Hij
    real(8)                               :: Jzz,Jxy
#endif
    !
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: connect Spin blocks"
#endif
    !
    !Hij is shared:
    !Hij = Hmodel(left,right)
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij, source=HopH)
    Jzz = Hij(1,1)
    Jxy = Hij(2,2)
    !
    !> Get H2 dimensions:
    Hdims = shape(left%operators)*shape(right%operators)
    if(present(states))Hdims = [size(states),size(states)]
    call H2%init(Hdims(1),Hdims(2))
    !
    !>Retrieve operators:
    bcdim=1
    if(PBCdmrg)bcdim=2
    allocate(SLz(bcdim),SLp(bcdim))
    allocate(SRz(bcdim),SRp(bcdim))
    !
    if(PBCdmrg)then
       SLz(1) = left%operators%op("S"//left%okey(0,1,ilink="L"))   !left.S_z_L
       SLz(2) = left%operators%op("S"//left%okey(0,1,ilink="R"))   !left.S_z_R
       SLp(1) = left%operators%op("S"//left%okey(0,2,ilink="L"))   !left.S_p_L
       SLp(2) = left%operators%op("S"//left%okey(0,2,ilink="R"))   !left.S_p_R
       !          
       SRz(1) = right%operators%op("S"//right%okey(0,1,ilink="L")) !right.S_z_L
       SRz(2) = right%operators%op("S"//right%okey(0,1,ilink="R")) !right.S_z_R
       SRp(1) = right%operators%op("S"//right%okey(0,2,ilink="L")) !right.S_p_L
       SRp(2) = right%operators%op("S"//right%okey(0,2,ilink="R")) !right.S_p_R
    else
       SLz(1) = left%operators%op("S"//left%okey(0,1))
       SLp(1) = left%operators%op("S"//left%okey(0,2))
       SRz(1) = right%operators%op("S"//right%okey(0,1))
       SRp(1) = right%operators%op("S"//right%okey(0,2))
    endif
    !
    !
    !>Build H2:
    select case(PBCdmrg)
    case (.true.)
       !===[L,r]-[R,l]---
       !S^z_L,r x S^z_R,l + S^+_L,r x S^-_R,l + S^-_L,r  x S^+_R,l 
       !SLz(2) x SRz(1)   + SLp(2) x SRp(1)*  + SLp(2)*  x SRp(1)
       !
       !|[L,l]===---[R,r]|
       !\---------------/
       !S^z_L,l x S^z_R,r + S^+_L,l x S^-_R,r + S^-_L,l x S^+_R,r 
       !SLz(1)  x SRz(2)  + SLp(1) x SRp(2)*  + SLp(1)* x SRp(2)
       if(present(states))then
          H2 = H2 + &
               Jzz*sp_kron(SLz(2),SRz(1),states)      + &
               Jxy*sp_kron(SLp(2),SRp(1)%dgr(),states)+ &
               Jxy*sp_kron(SLp(2)%dgr(),SRp(1),states)+ &
               Jzz*sp_kron(SLz(1),SRz(2),states)      + &
               Jxy*sp_kron(SLp(1),SRp(2)%dgr(),states)+ &
               Jxy*sp_kron(SLp(1)%dgr(),SRp(2),states)
       else
          H2 = H2 + &
               Jzz*(SLz(2).x.SRz(1))      + &
               Jxy*(SLp(2).x.SRp(1)%dgr())+ &
               Jxy*(SLp(2)%dgr().x.SRp(1))+ &
               Jzz*(SLz(1).x.SRz(2))      + &
               Jxy*(SLp(1).x.SRp(2)%dgr())+ &
               Jxy*(SLp(1)%dgr().x.SRp(2))
       endif
    case (.false.)
       !===[L,r]-[R,l]---
       !S^z_L,r x S^z_R,l + S^+_L,r x S^-_R,l + S^-_L,r  x S^+_R,l 
       !SLz(1) x SRz(1)   + SLp(1)  x SRp(1)* + SLp(1)* x SRp(1)
       if(present(states))then
          H2 = H2 + &
               Jzz*sp_kron(SLz(1),SRz(1),states)      + &
               Jxy*sp_kron(SLp(1),SRp(1)%dgr(),states)+ &
               Jxy*sp_kron(SLp(1)%dgr(),SRp(1),states)
       else
          H2 = H2 + &
               Jzz*(SLz(1).x.SRz(1))       + &
               Jxy*(SLp(1).x.SRp(1)%dgr()) + &
               Jxy*(SLp(1)%dgr().x.SRp(1))
       endif
    end select
    !
    !
    !> Free memory
    do ibc=1,bcdim
       call SLz(ibc)%free
       call SLp(ibc)%free
       call SRz(ibc)%free
       call SRp(ibc)%free
    enddo
  end function connect_spin_blocks












  !H_lr = \sum_{a}h_aa*(C^+_{left,a}@P_left) x C_{right,a}] + H.c.
  function connect_fermion_blocks(left,right,states) result(H2)
    type(block)                           :: left
    type(block)                           :: right
    integer,dimension(:),optional         :: states
    type(sparse_matrix),allocatable       :: CL(:,:),CR(:,:),P(:)
    type(sparse_matrix)                   :: H2
    integer                               :: ispin,iorb,jorb,io,jo,ibc,bcdim
    integer,dimension(2)                  :: Hdims,dleft,dright
    character(len=:),allocatable          :: key
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hij
    complex(8)                            :: Tr,Tl
#else
    real(8),dimension(:,:),allocatable    :: Hij
    real(8)                               :: Tr,Tl
#endif

#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: connect Fermion blocks"
#endif
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
    if(present(states))Hdims = size(states)
    call H2%init(Hdims(1),Hdims(2))
    !

    bcdim=1 ; if(PBCdmrg)bcdim=2
    allocate(Cl(bcdim,Nspin*Norb))
    allocate(Cr(bcdim,Nspin*Norb))
    allocate(P(bcdim))
    !
    !>Retrieve operators:
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          if(PBCdmrg)then
             CL(1,io) = left%operators%op( "C"//left%okey(iorb,ispin,ilink="L"))
             CL(2,io) = left%operators%op( "C"//left%okey(iorb,ispin,ilink="R"))
             CR(1,io) = right%operators%op("C"//left%okey(iorb,ispin,ilink="L"))
             CR(2,io) = right%operators%op("C"//left%okey(iorb,ispin,ilink="R"))
          else
             CL(1,io) = left%operators%op( "C"//left%okey(iorb,ispin))
             CR(1,io) = right%operators%op("C"//left%okey(iorb,ispin))
          end if
       enddo
    enddo
    !
    !
    !>Build H2:
    if(PBCdmrg)then
       !===[L,r]-[R,l]---
       ! [Cdg(a)_L,r]@[P_L,r] x C(b)_R,l + [P_L,r]@C(a)_L,r x Cdg(b)_R,l
       ! [CL(2,io)^+.P(2)]    x CR(1,jo) + [P(2).CL(2,io)] x CR^+(1,jo)
       !
       !-[L,l]===---[R,r]-
       !|----------------|
       ! [Cdg(a)_L,l]@[P_L,l] x C(b)_R,r + [P_L,l]@C(a)_L,l x Cdg(b)_R,r
       ! [CL(1,io)^+.P(1)]    x CR(2,jo) + [P(1).Cl(1,io)]  x CR^+(2,jo)
       !
       P(1) = left%operators%op("P"//left%okey(0,0,ilink="L"))
       P(2) = left%operators%op("P"//left%okey(0,0,ilink="R"))
       do io=1,Nspin*Norb
          do jo=1,Nspin*Norb
             if(Hij(io,jo)==0d0)cycle
#ifdef _CMPLX
             Tr = Hij(io,jo)
             Tl = conjg(Tr)
#else
             Tr = Hij(io,jo)
             Tl = Tr
#endif
             if(present(states))then
                H2 = H2 + &
                     Tr*sp_kron( matmul(Cl(2,io)%dgr(),P(2)), Cr(1,jo), states) + &
                     Tl*sp_kron( matmul(P(2),Cl(2,io)), Cr(1,jo)%dgr(), states) + &
                     Tr*sp_kron( matmul(Cl(1,io)%dgr(),P(1)), Cr(2,jo), states) + &
                     Tl*sp_kron( matmul(P(1),Cl(1,io)), Cr(2,jo)%dgr(), states) 
             else
                H2 = H2 + &
                     Tr*( matmul(Cl(2,io)%dgr(),P(2)).x.Cr(1,jo)) + &
                     Tl*( matmul(P(2),Cl(2,io)).x.Cr(1,jo)%dgr()) + &
                     Tr*( matmul(Cl(1,io)%dgr(),P(1)).x.Cr(2,jo)) + &
                     Tl*( matmul(P(1),Cl(1,io)).x.Cr(2,jo)%dgr()) 
             endif
          enddo
       enddo
    else
       !===[L,r]-[R,l]---
       ! [Cdg(a)_L,r.P] x C(b)_R,l + [P.C(a)_L,r x Cdg(b)_R,l
       ! [Cl(1,io)^+.P] x Cr(1,jo) + [P.Cl(1,io)] x Cr^+(1,jo)
       P(1) = left%operators%op("P")
       do io=1,Nspin*Norb
          do jo=1,Nspin*Norb
             if(Hij(io,jo)==0d0)cycle
#ifdef _CMPLX
             Tr = Hij(io,jo)
             Tl = conjg(Tr)
#else
             Tr = Hij(io,jo)
             Tl = Tr
#endif
             if(present(states))then
                H2 = H2 + &
                     Tr*sp_kron(matmul(Cl(1,io)%dgr(),P(1)),Cr(1,jo),states) + &
                     Tl*sp_kron(matmul(P(1),Cl(1,io)),Cr(1,jo)%dgr(),states)
             else
                H2 = H2 + &
                     Tr*(matmul(Cl(1,io)%dgr(),P(1)).x.Cr(1,jo)) + &
                     Tl*(matmul(P(1),Cl(1,io)).x.Cr(1,jo)%dgr())
             endif
          enddo
       enddo
    end if
    !
    !> free memory
    do ibc=1,bcdim
       call P(ibc)%free
       do io=1,Nspin*Norb
          call Cl(ibc,io)%free
          call Cr(ibc,io)%free
       enddo
    enddo
  end function connect_fermion_blocks






END MODULE DMRG_CONNECT








  
!   subroutine enlarge_block_old(self,dot,grow)
!     type(block)                  :: self
!     type(site)                   :: dot
!     character(len=*),optional    :: grow
!     character(len=16)            :: grow_
!     character(len=:),allocatable :: key,dtype,otype
!     type(tbasis)                 :: self_basis,dot_basis,enl_basis
!     type(sparse_matrix)          :: Hb,Hd,H2,eO
!     integer                      :: i
!     !
!     grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
!     !
! #ifdef _DEBUG
!     write(LOGfile,*)"DEBUG: ENLARGE block"//str(grow_)
! #endif
!     !
!     call start_timer("Enlarge blocks "//str(grow_))
!     !
!     if(.not.self%operators%has_key("H"))&
!          stop "Enlarge_Block ERROR: Missing self.H operator in the list"
!     if(.not.dot%operators%has_key("H"))&
!          stop "Enlarge_Block ERROR: Missing dot.H operator in the list"
!     !
!     dtype=dot%type()
!     if(dtype/=self%type())&
!          stop "Enlarge_Block ERROR: Dot.Type != Self.Type"
!     !    
!     !> Update Hamiltonian:
! #ifdef _DEBUG
!     write(LOGfile,*)"DEBUG: ENLARGE block: update H"
! #endif
!     select case(str(grow_))
!     case ("left","l")
!        Hb = self%operators%op("H").x.id(dot%dim)
!        Hd = id(self%dim).x.dot%operators%op("H")
!        select case(dtype)
!        case default;stop "Enlarge_Block ERROR: wrong dot.Type"
!        case ("spin","s")
!           H2 = connect_spin_blocks(self,as_block(dot))
!        case ("fermion","f,","electron","e")
!           H2 = connect_fermion_blocks(self,as_block(dot))
!        end select
!     case ("right","r")
!        Hb = id(dot%dim).x.self%operators%op("H")
!        Hd = dot%operators%op("H").x.id(self%dim)
!        select case(dtype)
!        case default;stop "Enlarge_Block ERROR: wrong dot.Type"
!        case ("spin","s")
!           H2 = connect_spin_blocks(as_block(dot),self)
!        case ("fermion","f,","electron","e")
!           H2 = connect_fermion_blocks(as_block(dot),self)
!        end select
!     end select
!     call self%put_op("H", Hb +  Hd + H2, type="bosonic")
!     !
!     !> Update all the other operators in the list:
! #ifdef _DEBUG
!     write(LOGfile,*)"DEBUG: ENLARGE block: update Op list"
! #endif
!     do i=1,size(self%operators)
!        key   = str(self%operators%key(index=i))
!        otype = str(self%operators%type(index=i))
!        if(key=="H")cycle
!        !
!        !Bosonic operators:
!        !O_L -> I_L.x.O_d  | O_R -> O_d.x.I_R
!        !Fermionic operators:
!        !O_L -> P_L.x.O_d  | O_R -> O_d.x.I_R
!        !Sign operator:
!        !P_L -> P_L.x.P_d  | P_R -> P_d.x.P_R
!        select case(str(grow_))
!        case ("left","l")
!           select case(to_lower(str(otype)))
!           case default;stop "Enlarge_BLock ERROR: wrong self.operators.type L: !\in['Bosonic','Fermionic','Sign']"
!           case('b','bose','bosonic')
!              eO = Id(self%dim)
!           case('f','fermi','fermionic')
!              eO = self%operators%op(key="P")
!           case('s','sign')
!              eO = self%operators%op(key="P")
!           end select
!           call self%put_op(str(key), eO.x.dot%operators%op(str(key)), type=otype)
!        case ("right","r")
!           select case(to_lower(str(otype)))
!           case default;stop "Enlarge_BLock ERROR: wrong self.operators.type R: !\in['Bosonic','Fermionic','Sign']"
!           case('b','bose','bosonic')
!              eO = Id(self%dim)
!           case('f','fermi','fermionic')
!              eO = Id(self%dim)
!           case('s','sign')
!              eO = self%operators%op(key="P")
!           end select
!           call self%put_op(str(key), dot%operators%op(str(key)).x.eO, type=otype )
!        end select
!     enddo
!     !
!     !> Enlarge dimensions
!     self%length = self%length + 1
!     self%Dim    = self%Dim*dot%Dim
!     !
!     !> Enlarge the basis states
!     call self%get_basis(self_basis)
!     call dot%get_basis(dot_basis)
!     select case(str(grow_))
!     case ("left","l")
!        enl_basis = (self_basis.o.dot_basis)
!        call self%set_basis( basis=enl_basis )
!     case ("right","r")
!        enl_basis = (dot_basis.o.self_basis)
!        call self%set_basis( basis=enl_basis )
!     end select
!     !
! #ifdef _DEBUG
!     if(verbose>5)call self%show(file="Enl"//str(grow_)//"_"//str(self%length)//".dat")
! #endif
!     !
!     if(.not.self%is_valid())then
!        write(LOGfile,*)"dmrg_step error: enlarged_block "//str(grow_)// "is not a valid block"
!        stop
!     endif
!     !
!     !Free the memory:
!     call Hb%free()
!     call Hd%free()
!     call H2%free()
!     call self_basis%free()
!     call dot_basis%free()
!     call enl_basis%free()
!     !
!     call stop_timer()
!     !
!   end subroutine enlarge_block_old
