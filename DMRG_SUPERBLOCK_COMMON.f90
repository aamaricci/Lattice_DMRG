MODULE DMRG_SUPERBLOCK_COMMON
  USE VARS_GLOBAL
  USE DMRG_CONNECT
  implicit none


  !TO BE MOVED TO SUPERBLOCK
  !
  integer                                        :: tNso
  integer                                        :: isb,jsb
  type(sparse_matrix),allocatable,dimension(:)   :: Hleft,Hright
  type(sparse_matrix),allocatable,dimension(:,:) :: A,B
  integer,dimension(:),allocatable               :: Dls,Drs,Offset
  integer,dimension(:,:),allocatable             :: RowOffset,ColOffset
  !
  integer                                        :: i,j
  integer                                        :: ispin
  integer                                        :: iorb,jorb
  integer                                        :: io,jo



contains


  !##################################################################
  !              RETURN LEFT.states or RIGHT.states 
  !              contributing to the SUPERBLOCK with
  !              a given QN
  !. sb_sector: inherited from VARS_GLOBAL
  !##################################################################
  function sb2block_states(q,label) result(states)
    real(8),dimension(:)             :: q
    character(len=*)                 :: label
    integer,dimension(:),allocatable :: tmp,states,sb_map
    integer                          :: i,istate,l,r,isb
    !
    if(.not.associated(sb_sector%root))&
         stop "sb2block_states error: sb_sector is not allocated"
    !
    if(allocated(states))deallocate(states)
    !
    !> get the map from the QN of the sector:
    sb_map = sb_sector%map(qn=q)
    !
    !> left,right, sb_sector and sb_states have to be known at this time:
    ! add a check
    select case(to_lower(str(label)))
    case("left","l","sys","s")
       do i=1,size(sb_map)
          istate = sb_states(sb_map(i))
          l = (istate-1)/right%Dim+1
          call append(tmp,l)
       enddo
    case("right","r","env","e")
       do i=1,size(sb_map)
          istate = sb_states(sb_map(i))
          r = mod(istate,right%Dim);if(r==0)r=right%Dim
          call append(tmp,r)
       enddo
    end select
    allocate(states, source=uniq(tmp))
  end function sb2block_states


END MODULE DMRG_SUPERBLOCK_COMMON









! subroutine Setup_SuperBlock_Direct()
!   integer                                      :: Nso,Nsb
!   integer                                      :: it,isb,jsb
!   real(8),dimension(:),allocatable             :: qn,qm
!   type(tstates),dimension(:),allocatable       :: Ai,Aj,Bi,Bj
!   real(8),dimension(:),allocatable             :: dq
!   real(8),dimension(2),parameter               :: qnup=[1d0,0d0],qndw=[0d0,1d0]
!   integer,dimension(:,:,:),allocatable         :: tMap
!   type(sparse_matrix)                          :: P
!   type(sparse_matrix),allocatable,dimension(:) :: Cleft,Cright
!   type(sparse_matrix),allocatable,dimension(:) :: Sleft,Sright
!   real(8),dimension(:,:),allocatable           :: Hij
!   integer                                      :: m_left,m_right
!   character(len=:),allocatable                 :: type
!   !
!   if(.not.left%operators%has_key("H"))&
!        stop "Setup_SuperBlock_Direct ERROR: Missing left.H operator in the list"
!   if(.not.right%operators%has_key("H"))&
!        stop "Setup_SuperBlock_Direct ERROR: Missing right.H operator in the list"
!   !
!   type=left%type()
!   if(type/=right%type())&
!        stop "Setup_SuperBlock_Direct ERROR: left.Type != right.Type"
!   !
!   m_left = left%dim
!   m_right= right%dim
!   !
!   !> GET THE USER DEFINED MODEL HAMILTONIAN PARAMETERS:
!   Hij = Hmodel(left,right)
!   !
!   !
!   select case(type)
!   case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.Type"
!   case ("spin","s")
!      Nso = Nspin
!   case ("fermion","f,","electron","e")
!      Nso = Nspin*Norb
!   end select
!   !
!   Nsb = size(sb_sector)
!   !    
!   !Creating the sequence of operators A*_q, B*_q
!   ! which decompose the term H^LR of the
!   ! super-block Hamiltonian.
!   select case(type)
!   case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.Type"
!   case ("spin","s")
!      tNso = 3
!      allocate(tMap(3,1,1))
!      it = 0
!      do i=1,tNso
!         it = it+1
!         tMap(i,1,1)=it
!      enddo
!   case ("fermion","f,","electron","e")
!      tNso = 2*count(Hij/=0d0)
!      allocate(tMap(2,Nso,Nso))
!      it = 0
!      do i=1,2
!         do io=1,Nso
!            do jo=1,Nso
!               if(Hij(io,jo)==0d0)cycle
!               it = it+1
!               tMap(i,io,jo)=it
!            enddo
!         enddo
!      enddo
!   end select
!   !
!   allocate(Dls(Nsb),Drs(Nsb),Offset(Nsb))
!   Offset=0
!   do isb=1,Nsb
!      qn   = sb_sector%qn(index=isb)
!      Dls(isb)= sector_qn_dim(left%sectors(1),qn)
!      Drs(isb)= sector_qn_dim(right%sectors(1),current_target_qn - qn)
!      if(isb>1)Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
!   enddo
!   !
!   print*,"DIMS:",sum(Dls),sum(Drs)
!   allocate(AI(Nsb),BI(Nsb))
!   allocate(AJ(Nsb),BJ(Nsb))
!   select case(type)
!   case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.Type"
!   case ("spin","s")
!      dq=[1d0]
!      do isb=1,size(sb_sector)
!         qn   = sb_sector%qn(index=isb)
!         AI(isb)%states = sb2block_states(qn,'left')
!         BI(isb)%states = sb2block_states(qn,'right')
!         print*,size(BI(isb)%states),Drs(isb),Dls(isb)
!         qm   = qn - dq
!         if(.not.sb_sector%has_qn(qm))cycle
!         jsb = sb_sector%index(qn=qm)
!         AJ(jsb)%states = sb2block_states(qm,'left')
!         BJ(jsb)%states = sb2block_states(qm,'right')
!      enddo
!   case ("fermion","f,","electron","e")
!      do isb=1,size(sb_sector)
!         qn   = sb_sector%qn(index=isb)
!         AI(isb)%states = sb2block_states(qn,'left')
!         BI(isb)%states = sb2block_states(qn,'right')
!         do ispin=1,Nspin
!            dq = qnup ; if(ispin==2)dq=qndw
!            qm = qn - dq
!            if(.not.sb_sector%has_qn(qm))cycle
!            jsb = sb_sector%index(qn=qm)
!            AJ(jsb)%states = sb2block_states(qm,'left')
!            BJ(jsb)%states = sb2block_states(qm,'right')
!         enddo
!      enddo
!   end select
!   !
!   !
!   allocate(A(tNso,Nsb),B(tNso,Nsb))
!   allocate(Hleft(Nsb),Hright(Nsb))
!   allocate(RowOffset(tNso,Nsb),ColOffset(tNso,Nsb))
!   select case(type)
!   case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.Type"
!   case ("spin","s")
!      allocate(Sleft(Nspin),Sright(Nspin))
!      do ispin=1,Nspin
!         Sleft(ispin)  = left%operators%op("S"//left%okey(0,ispin))
!         Sright(ispin) = right%operators%op("S"//right%okey(0,ispin))
!      enddo
!      !
!      print*,"Constructing H^{L,R}*: the filtered block Hamiltonians"
!      do isb=1,size(sb_sector)
!         qn = sb_sector%qn(index=isb)
!         !
!         !> get: H*^L  and H*^R
!         Hleft(isb) = sp_filter(left%operators%op("H"),AI(isb)%states)
!         Hright(isb)= sp_filter(right%operators%op("H"),BI(isb)%states)
!         print*,"L",shape(Hleft(isb)),Dls(isb)
!         print*,"R",shape(Hright(isb)),Drs(isb)

!         !> get: A = Jp*S_lz .x. B = S_rz + Row/Col Offsets       
!         it=tMap(1,1,1)
!         A(it,isb) = Hij(1,1)*sp_filter(Sleft(1),AI(isb)%states,AI(isb)%states)
!         B(it,isb) = sp_filter(Sright(1),BI(isb)%states,BI(isb)%states)
!         RowOffset(it,isb)=Offset(isb)           
!         ColOffset(it,isb)=Offset(isb)
!         !
!         dq = [1d0]
!         qm = qn - dq
!         if(.not.sb_sector%has_qn(qm))cycle
!         jsb = sb_sector%index(qn=qm)
!         !> get: A = Jp*S_l- .x. B = [S_r-]^+=S_r+ + Row/Col Offsets 
!         it=tMap(2,1,1)
!         A(it,isb) = Hij(2,2)*sp_filter(Sleft(2),AI(isb)%states,AJ(jsb)%states)
!         B(it,isb) = sp_filter(hconjg(Sright(2)),BI(isb)%states,BJ(jsb)%states)
!         RowOffset(it,isb)=Offset(isb)
!         ColOffset(it,isb)=Offset(jsb)!right%sectors(1)%index(qn=qm))
!         !> get H.c. + Row/Col Offsets
!         it=tMap(3,1,1)
!         A(it,isb) = hconjg(A(tMap(2,1,1),isb))
!         B(it,isb) = hconjg(B(tMap(2,1,1),isb))
!         RowOffset(it,isb)=Offset(jsb)!right%sectors(1)%index(qn=qm))
!         ColOffset(it,isb)=Offset(isb)
!      enddo


!   case ("fermion","f,","electron","e")
!      allocate(Cleft(Nso),Cright(Nso))
!      P = left%operators%op("P")
!      do ispin=1,Nspin
!         do iorb=1,Norb
!            io = iorb+(ispin-1)*Norb
!            Cleft(io) = left%operators%op("C"//left%okey(iorb,ispin))
!            Cright(io) = right%operators%op("C"//right%okey(iorb,ispin))
!         enddo
!      enddo

!      print*,"Constructing H^{L,R}*: the filtered block Hamiltonians"
!      do isb=1,size(sb_sector)
!         qn   = sb_sector%qn(index=isb)
!         !
!         !> get: H*^L  and H*^R
!         Hleft(isb) = sp_filter(left%operators%op("H"),AI(isb)%states)
!         Hright(isb)= sp_filter(right%operators%op("H"),BI(isb)%states)
!         !
!         !> get A.x.B terms:
!         do ispin=1,Nspin
!            dq = qnup ; if(ispin==2)dq=qndw
!            qm = qn - dq
!            if(.not.sb_sector%has_qn(qm))cycle
!            jsb = sb_sector%index(qn=qm)
!            !
!            do iorb=1,Norb
!               do jorb=1,Norb
!                  io = iorb+(ispin-1)*Norb
!                  jo = jorb+(ispin-1)*Norb
!                  if(Hij(io,jo)==0d0)cycle
!                  !
!                  !> get: A = H(a,b)*[Cl(a,s)^+@P] .x. B = Cr(b,s)  + Row/Col Offsets
!                  it=tMap(1,io,jo)
!                  A(it,isb) = Hij(io,jo)*sp_filter(matmul(Cleft(io)%dgr(),P),AI(isb)%states,AJ(jsb)%states)
!                  B(it,isb) = sp_filter(Cright(jo),BI(isb)%states,BJ(jsb)%states)
!                  RowOffset(it,isb)=Offset(isb)           
!                  ColOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
!                  !
!                  !> get H.c. + Row/Col Offsets
!                  it=tMap(2,io,jo)
!                  A(it,isb) = hconjg(A(tMap(1,io,jo),isb))
!                  B(it,isb) = hconjg(B(tMap(1,io,jo),isb))
!                  RowOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
!                  ColOffset(it,isb)=Offset(isb)
!               enddo
!            enddo
!            !
!         enddo
!      enddo
!   end select
!   !
!   print*,"Done:"
!   print*,"######################################"
!   print*,""
! end subroutine Setup_SuperBlock_Direct
