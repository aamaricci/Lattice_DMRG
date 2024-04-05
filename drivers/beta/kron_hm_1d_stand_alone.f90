program testEDkron
  USE SCIFOR
  USE DMRG, id=>sp_eye
  USE OMP_LIB, only: omp_get_wtime
  implicit none

  integer                                        :: Nso
  character(len=64)                              :: finput
  integer                                        :: i,j,Nsb,k,l,r,m,iorb,jorb,ispin,it,ip,io,jo
  character(len=1)                               :: DMRGtype
  real(8)                                        :: ts(2),Mh,target_qn(2),suml,lambda,vh
  type(site)                                     :: Dot
  type(sparse_matrix)                            :: C,N,Cl,Cr,P
  type(sparse_matrix),allocatable,dimension(:)   :: Hsb
  real(8),dimension(:,:),allocatable             :: Hloc,Hij
  real(8),dimension(:),allocatable               :: sb_qn
  type(block)                                    :: left,right
  type(sparse_matrix)                            :: spHsb,spH
  real(8)                                        :: gs_energy
  integer                                        :: m_sb
  real(8),dimension(:,:),allocatable             :: Evecs,Rho
  real(8),dimension(:),allocatable               :: Evals
  real(8),dimension(:),allocatable               :: Vec,Hvec,Hvec_
  integer,dimension(:),allocatable               :: sb_states,sb_map
  type(sectors_list)                             :: sb_sector
  integer                                        :: Neigen=2
  integer                                        :: current_L
  real(8),dimension(:),allocatable               :: current_target_QN





  !TO BE MOVED TO SUPERBLOCK
  !
  integer                                        :: tNso
  integer                                        :: isb,jsb
  type(sparse_matrix),allocatable,dimension(:)   :: Hleft,Hright
  type(sparse_matrix),allocatable,dimension(:,:) :: A,B
  integer,dimension(:),allocatable               :: Dls,Drs,Offset
  integer,dimension(:,:),allocatable             :: RowOffset,ColOffset
  real                                           :: t0


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=(/( -0.5d0,i=1,2 )/),&
       comment="Hopping amplitudes")
  call parse_input_variable(lambda,"lambda",finput,default=0d0,&
       comment="Hybridization")
  call parse_input_variable(Mh,"MH",finput,default=0d0,&
       comment="Crystal field splittings")
  call parse_input_variable(vh,"vh",finput,default=0d0,&
       comment="local hybridization amplitude")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)


  ! if(Norb/=2)stop "This code is for Norb=2. STOP"
  Nso = Nspin*Norb

  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = 0d0
  if(Norb==2)Hloc = Mh*kron(pauli_0,pauli_z) + vh*kron(pauli_0,pauli_x)

  Dot  = electron_site(Hloc)
  call dot%show()


  call init_dmrg(hm_1d_model,ModelDot=Dot)
  target_qn = DMRG_qn
  !







  print*,""
  print*,""
  print*,"######################################"
  print*,"   o->o + o<-o"
  print*,"######################################"
  print*,""
  print*,""

  t0 = omp_get_wtime()
  left  = block(dot)
  call enlarge_block_(left,dot,grow='left')
  right = block(dot)
  call enlarge_block_(right,dot,grow='right')
  print*,omp_get_wtime() - t0
  !
  print*,"Get SB states"
  t0 = omp_get_wtime()
  call get_sb_states_(left,right,sb_states,sb_sector)
  print*,omp_get_wtime() - t0


  call build_superblock(left,right,sb_states,sb_sector)



  ! !> ONE STEP TEST:
  ! allocate(Vec(size(sb_states)))
  ! call random_number(vec)
  ! allocate(Hvec, source=vec)
  ! allocate(Hvec_, source=vec)
  ! !
  ! print*,"######################################"
  ! print*,"Test H.v OLD vs NEW method random vector v:"
  ! ! !Set spHsb to H^Lx1^R to be used in the default spHv (which only uses spHsb)
  ! t0 = omp_get_wtime()
  ! spHsb  =  h2_model(left,right,sb_states) + &
  !      sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
  !      sp_kron(id(left%dim),right%operators%op("H"),sb_states)
  ! call sb_HxV(size(sb_states),vec,Hvec_)
  ! print*,omp_get_wtime() - t0
  ! !
  ! !Apply H to v new style
  ! t0 = omp_get_wtime()
  ! call sb_HxV_new(size(sb_states),vec,Hvec)
  ! print*,omp_get_wtime() - t0
  ! !
  ! do isb=1,size(sb_sector)
  !    print*,"Sector Dim:",Dls(isb)*Drs(isb)
  !    print*,"       Hv_old             ","       Hv_new             ","       Error              "
  !    do i=1,Dls(isb)*Drs(isb)
  !       k = i + Offset(isb)
  !       print*,Hvec_(k),Hvec(k),abs(Hvec(k)-Hvec_(k))
  !    enddo
  ! enddo
  ! print*,"Any error:",any(abs(Hvec-Hvec_)>1d-12)
  ! print*,"######################################"
  ! print*,""
  ! print*,""




  !Old Style solution:
  print*,"Old style solution: spH --> H.v --> \psi"
  print*,"######################################"
  t0 = omp_get_wtime()
  spHsb  =  h2_model(left,right,sb_states)   + &
       sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
       sp_kron(id(left%dim),right%operators%op("H"),sb_states)
  allocate(Evals(Neigen))
  allocate(Evecs(size(sb_states),Neigen))
  call sp_eigh(sb_HxV,evals,evecs,&
       3*Neigen,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,Neigen
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""
  !So it won't work anymore:
  call spHsb%free()
  print*,omp_get_wtime() - t0


  !
  !Using _new to get energy:
  print*,"Solving the same problem using new Method:"
  print*,"######################################"
  t0 = omp_get_wtime()
  allocate(Evals(Neigen))
  allocate(Evecs(size(sb_states),Neigen))
  call sp_eigh(sb_HxV_new,evals,evecs,&
       3*Neigen,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,Neigen
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""
  print*,omp_get_wtime() - t0
















contains


  !This is user defined Function to be passed to the SYSTEM
  function hm_1d_model(left,right) result(Hlr)
    type(block)                        :: left
    type(block)                        :: right
    real(8),dimension(:,:),allocatable :: Hlr
    !
    if(allocated(Hlr))deallocate(Hlr)
    allocate(Hlr(Nspin*Norb,Nspin*Norb))
    !
    !workout local part, like random local field
    !if(left%Dim==1 AND right%Dim==1) then operate over local H
    !if(left%Dim==1 OR right%Dim==1) then operate over local H
    Hlr = diag([ts(1:Norb),ts(1:Norb)])
    if(Norb==2)Hlr = Hlr + lambda*kron(pauli_0,pauli_x)
  end function hm_1d_model



  subroutine build_superblock(left,right,sb_states,sb_sector)
    type(block)                                  :: left,right
    integer,dimension(:),intent(in)              :: sb_states
    type(sectors_list)                           :: sb_sector
    integer                                      :: it,isb,jsb
    real(8),dimension(:),allocatable             :: qn,qm
    type(tstates),dimension(:),allocatable       :: Ai,Aj,Bi,Bj
    real(8),dimension(2)                         :: dq
    type(sparse_matrix),allocatable,dimension(:) :: Cleft,Cright
    real(8),dimension(2),parameter               :: qnup=[1d0,0d0],qndw=[0d0,1d0]
    real                                         :: t0
    integer,dimension(:,:,:),allocatable         :: tMap



    print*,"######################################"
    print*,"Set up new solution method: [H*].v --> \psi"
    !Build filterd H*_L, H*_R
    Nsb = size(sb_sector)
    print*,"There are Nsb_states:",Nsb

    Hij = hm_1d_model(left,right)


    !Creating the sequence of operators A*_q, B*_q which decompose the term H^LR of the
    ! super-block Hamiltonian.
    print*,"Constructing \sum_a=1,tNso A*.B*"
    tNso = 2*count(Hij/=0d0)
    allocate(tMap(2,Nso,Nso))
    it = 0
    do i=1,2
       do io=1,Nso
          do jo=1,Nso
             if(Hij(io,jo)==0d0)cycle
             it = it+1
             tMap(i,io,jo)=it
          enddo
       enddo
    enddo
    print*,"There are tNso non-zero elements ",tNso


    print*,"Constructing Offset and D_l,r"
    allocate(Dls(Nsb),Drs(Nsb),Offset(Nsb))
    Offset=0
    do isb=1,Nsb
       sb_qn   = sb_sector%qn(index=isb)
       Dls(isb)= dim_sector(left%sectors(1),sb_qn)
       Drs(isb)= dim_sector(right%sectors(1),current_target_qn - sb_qn)
       if(isb>1)Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
    enddo


    !
    print*,"Constructing filtered SB States"  
    allocate(AI(Nsb),BI(Nsb))
    allocate(AJ(Nsb),BJ(Nsb))
    do isb=1,size(sb_sector)
       qn   = sb_sector%qn(index=isb)
       AI(isb)%states = sb2block_states(qn,'left')
       BI(isb)%states = sb2block_states(qn,'right')
       do ispin=1,Nspin
          dq = qnup ; if(ispin==2)dq=qndw
          qm = qn - dq
          if(.not.sb_sector%has_qn(qm))cycle
          jsb = sb_sector%index(qn=qm)
          AJ(jsb)%states = sb2block_states(qm,'left')
          BJ(jsb)%states = sb2block_states(qm,'right')
       enddo
    enddo


    !> Retrieve Fermions operators to speed up the process below:
    allocate(Hleft(Nsb),Hright(Nsb))
    allocate(A(tNso,Nsb),B(tNso,Nsb))
    allocate(RowOffset(tNso,Nsb),ColOffset(tNso,Nsb))
    allocate(Cleft(Nso),Cright(Nso))

    P = left%operators%op("P")
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb+(ispin-1)*Norb
          Cleft(io) = left%operators%op("C"//left%okey(iorb,ispin))
          Cright(io) = right%operators%op("C"//right%okey(iorb,ispin))
       enddo
    enddo
    !
    print*,"Constructing H^{L,R}*: the filtered block Hamiltonians"
    do isb=1,size(sb_sector)
       qn   = sb_sector%qn(index=isb)
       !
       !> get: H*^L  and H*^R
       Hleft(isb) = sp_filter(left%operators%op("H"),AI(isb)%states)
       Hright(isb)= sp_filter(right%operators%op("H"),BI(isb)%states)
       !
       !> get A.x.B terms:
       do ispin=1,Nspin
          dq = qnup ; if(ispin==2)dq=qndw
          qm   = qn - dq
          if(.not.sb_sector%has_qn(qm))cycle
          jsb = sb_sector%index(qn=qm)
          !
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb+(ispin-1)*Norb
                jo = jorb+(ispin-1)*Norb
                if(Hij(io,jo)==0d0)cycle
                !
                !> get: A = H(a,b)*[Cl(a,s)^+@P] .x. B = Cr(b,s)  + Row/Col Offsets
                it=tMap(1,io,jo)
                A(it,isb) = Hij(io,jo)*sp_filter(matmul(Cleft(io)%dgr(),P),AI(isb)%states,AJ(jsb)%states)
                B(it,isb) = sp_filter(Cright(jo),BI(isb)%states,BJ(jsb)%states)
                RowOffset(it,isb)=Offset(isb)           
                ColOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
                !
                !> get H.c. + Row/Col Offsets
                it=tMap(2,io,jo)
                A(it,isb) = hconjg(A(tMap(1,io,jo),isb))
                B(it,isb) = hconjg(B(tMap(1,io,jo),isb))
                RowOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
                ColOffset(it,isb)=Offset(isb)
             enddo
          enddo
          !
       enddo
    enddo
    print*,"Done:"
    print*,"######################################"
    print*,""
  end subroutine build_superblock



  function dim_sector(self,qn) result(dim)
    type(sectors_list)   :: self
    real(8),dimension(:) :: qn
    integer              :: dim
    dim = 0
    if(.not.self%has_qn(qn))return
    dim =  size(self%map(qn=qn))
  end function dim_sector










  !-----------------------------------------------------------------!
  ! Purpose: enlarge a given BLOCK "SELF" growing it left/right with
  ! a SITE "DOT" (specified in the init)
  !-----------------------------------------------------------------!
  subroutine enlarge_block_(self,dot,grow)
    type(block)                  :: self
    type(site)                   :: dot
    character(len=*),optional    :: grow
    character(len=16)            :: grow_
    character(len=:),allocatable :: key
    type(tbasis)                 :: self_basis,dot_basis,enl_basis
    type(sparse_matrix)          :: Hb,Hd,H2
    integer                      :: i
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    if(.not.self%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: self.H not in the list"
    if(.not.dot%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: dot.H no in the list"
    !
    !> Update Hamiltonian:
    select case(str(grow_))
    case ("left","l")
       Hb = self%operators%op("H").x.id(dot%dim)
       Hd = id(self%dim).x.dot%operators%op("H")
       H2 = h2_model(self,as_block(dot))       
    case ("right","r")
       Hb = id(dot%dim).x.self%operators%op("H")
       Hd = dot%operators%op("H").x.id(self%dim)
       H2 = h2_model(as_block(dot),self)
    end select
    call self%put_op("H", Hb +  Hd + H2)
    !
    !> Update all the other operators in the list: 
    do i=1,size(self%operators)
       key = self%operators%key(index=i)
       if(str(key)=="H")cycle
       select case(str(grow_))
       case ("left","l")
          call self%put_op(str(key),Id(self%dim).x.dot%operators%op(str(key)))
       case ("right","r")
          call self%put_op(str(key),dot%operators%op(str(key)).x.Id(self%dim))
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
    !
    select case(str(grow_))
    case ("left","l")
       enl_basis = (self_basis.o.dot_basis)
       call self%set_basis( basis=enl_basis )
    case ("right","r")
       enl_basis = (dot_basis.o.self_basis)
       call self%set_basis( basis=enl_basis )
    end select
    !
    !Free the memory:
    call Hb%free()
    call Hd%free()
    call H2%free()
    call self_basis%free()
    call dot_basis%free()
    call enl_basis%free()
    !
  end subroutine enlarge_block_



  !H_lr = \sum_{a}h_aa*(C^+_{left,a}@P_left) x C_{right,a}] + H.c.
  function h2_model(left,right,states) result(H2)
    type(block)                           :: left
    type(block)                           :: right
    integer,dimension(:),optional         :: states
    type(sparse_matrix),dimension(Nspin*Norb) :: Cl,Cr
    type(sparse_matrix)                   :: P,A
    type(sparse_matrix)                   :: H2
    integer                               :: ispin,iorb,jorb,io,jo
    integer,dimension(2)                  :: Hdims
    character(len=:),allocatable          :: key
    !
    !Hij is shared:
    print*,"Hij:"
    Hij = hm_1d_model(left,right)
    print*,shape(Hij)
    !
    !> Get H2 dimensions:
    Hdims = shape(left%operators)*shape(right%operators)
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
             H2 = H2 + Hij(io,jo)*sp_kron(matmul(Cl(io)%dgr(),P),Cr(jo),states)
             H2 = H2 + Hij(io,jo)*sp_kron(matmul(P,Cl(io)),Cr(jo)%dgr(),states)
          else
             H2 = H2 + Hij(io,jo)*(matmul(Cl(io)%dgr(),P).x.Cr(jo))
             H2 = H2 + Hij(io,jo)*(matmul(P,Cl(io)).x.Cr(jo)%dgr())
          endif
       enddo
    enddo
    !
    !
    !> free memory
    call P%free
    do io=1,Nspin*Norb
       call Cl(io)%free
       call Cr(io)%free
    enddo
  end function h2_model



  !-----------------------------------------------------------------!
  ! Purpose: build the list of states compatible with the specified
  ! quantum numbers
  !-----------------------------------------------------------------!
  subroutine get_sb_states_(left,right,sb_states,sb_sector)
    type(block)                      :: left,right    
    integer,dimension(:),allocatable :: sb_states
    type(sectors_list)               :: sb_sector
    integer                          :: ileft,iright
    integer                          :: i,j,istate
    real(8),dimension(:),allocatable :: left_qn,right_qn
    integer,dimension(:),allocatable :: left_map,right_map

    !
    current_L         = left%length + right%length
    current_target_QN = int(target_qn*current_L*Norb)
    write(LOGfile,"(A22,I12)")"SuperBlock Length = ",current_L
    write(LOGfile,"(A22,"//str(size(current_target_QN))//"F12.7)")"Target_QN = ",current_target_QN
    write(LOGfile,"(A22,G12.7)")"Total         = ",sum(current_target_QN)
    write(LOGfile,"(A22,G12.7)")"Filling       = ",sum(current_target_QN)/current_L
    write(LOGfile,"(A22,G12.7)")"Filling/Norb  = ",sum(current_target_QN)/current_L/Norb
    if(allocated(sb_states))deallocate(sb_states)
    !
    call sb_sector%free()
    do ileft=1,size(left%sectors(1))
       left_qn  = left%sectors(1)%qn(index=ileft)
       right_qn  = current_target_qn - left_qn
       if(.not.right%sectors(1)%has_qn(right_qn))cycle
       !
       left_map = left%sectors(1)%map(qn=left_qn)
       right_map = right%sectors(1)%map(qn=right_qn)
       !
       do i=1,size(left_map)
          do j=1,size(right_map)
             istate=right_map(j) + (left_map(i)-1)*right%Dim
             call append(sb_states, istate)
             call sb_sector%append(qn=left_qn,istate=size(sb_states))
          enddo
       enddo
    enddo
    !
  end subroutine get_sb_states_





  ! function sb2block_states(left,right,sb_states,sb_map,label) result(states)
  function sb2block_states(q,label) result(states)
    ! type(block)                      :: left,right    
    ! integer,dimension(:)             :: sb_states,sb_map
    real(8),dimension(:)             :: q
    character(len=*)                 :: label
    integer,dimension(:),allocatable :: tmp,states,sb_map
    integer                          :: i,istate,l,r,isb
    !
    if(allocated(states))deallocate(states)
    !
    !> get the map from the QN of the sector:
    sb_map = sb_sector%map(qn=q)
    !> left,right, sb_sector and sb_states have to be known at this time:
    ! add a check
    !
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















  subroutine sb_HxV(Nloc,v,Hv)
    integer                 :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                 :: i,j,jcol
    Hv=zero
    do i=1,Nloc
       matmul: do jcol=1, spHsb%row(i)%Size
          val = spHsb%row(i)%vals(jcol)
          j   = spHsb%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end subroutine sb_HxV



  subroutine sb_HxV_new(Nloc,v,Hv)
    integer                 :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                 :: i,j,jcol,k,ir,il,jr,jl,n,arow,brow,acol,bcol,ia,ib,ic,ja,jb,jc
    real(8) :: aval,bval
    real(8),dimension(:,:),allocatable :: psi,Hpsi
    Hv=zero

    do k=1,size(sb_sector)

       !> apply the H^L x 1^r: need to T v and Hv
       do concurrent(ir=1:Drs(k))!< fix the column: iterate the row:
          do il=1,Dls(k)
             i = ir + (il-1)*Drs(k) + offset(k)
             do jcol=1,Hleft(k)%row(il)%Size
                val = Hleft(k)%row(il)%vals(jcol)
                jl  = Hleft(k)%row(il)%cols(jcol)
                j   = ir + (jl-1)*Drs(k) + offset(k)
                Hv(i) = Hv(i) + val*v(j)
             end do
          enddo
       enddo

       !
       !> apply the 1^L x H^r
       do concurrent( il=1:Dls(k))!< fix teh row: iterate the column:
          do ir=1,Drs(k)
             i = il + (ir-1)*Dls(k) + offset(k)           
             do jcol=1,Hright(k)%row(il)%Size
                val = Hright(k)%row(il)%vals(jcol)
                jl  = Hright(k)%row(il)%cols(jcol)
                j   = jl + (ir-1)*Dls(k) + offset(k)
                Hv(i) = Hv(i) + val*v(j)
             end do
          enddo
       enddo


       do it=1,tNso
          if(.not.A(it,k)%status.OR..not.B(it,k)%status)cycle
          do concurrent (ic=1:A(it,k)%Nrow*B(it,k)%Nrow)
             arow = (ic-1)/B(it,k)%Nrow+1
             brow = mod(ic,B(it,k)%Nrow);if(brow==0)brow=B(it,k)%Nrow
             if(A(it,k)%row(arow)%Size==0.OR.B(it,k)%row(brow)%Size==0)cycle
             i = ic + RowOffset(it,k)
             do ja=1,A(it,k)%row(arow)%Size
                acol = A(it,k)%row(arow)%cols(ja)
                aval = A(it,k)%row(arow)%vals(ja)
                do jb=1,B(it,k)%row(brow)%Size
                   bcol = B(it,k)%row(brow)%cols(jb)
                   bval = B(it,k)%row(brow)%vals(jb)
                   j = bcol+(acol-1)*B(it,k)%Ncol + ColOffset(it,k)
                   Hv(i) = Hv(i) + aval*bval*v(j)
                enddo
             enddo
          enddo
       enddo

    enddo

  end subroutine sb_HxV_new



end program testEDkron










!NEW

! print*,"######################################"
! print*,"Set up new solution method: [H*].v --> \psi"
! !Build filterd H*_L, H*_R
! Nsb = size(sb_sector)
! print*,"There are Nsb_states:",Nsb

! Hij = hm_1d_model(left,right)


! !Creating the sequence of operators A*_q, B*_q which decompose the term H^LR of the
! ! super-block Hamiltonian.
! print*,"Constructing \sum_a=1,tNso A*.B*"
! tNso = 2*count(Hij/=0d0)
! allocate(tMap(2,Nso,Nso))
! it = 0
! do i=1,2
!    do io=1,Nso
!       do jo=1,Nso
!          if(Hij(io,jo)==0d0)cycle
!          it = it+1
!          tMap(i,io,jo)=it
!       enddo
!    enddo
! enddo
! print*,"There are tNso non-zero elements ",tNso


! print*,"Constructing Offset and D_l,r"
! allocate(Dls(Nsb),Drs(Nsb),Offset(Nsb))
! Offset=0
! do isb=1,Nsb
!    sb_qn   = sb_sector%qn(index=isb)
!    !
!    Dls(isb)= dim_sector(left%sectors(1),sb_qn)
!    Drs(isb)= dim_sector(right%sectors(1),current_target_qn - sb_qn)
!    !
!    if(isb>1)Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
! enddo


! allocate(A(tNso,Nsb),B(tNso,Nsb))
! allocate(Hleft(Nsb),Hright(Nsb))
! allocate(RowOffset(tNso,Nsb),ColOffset(tNso,Nsb))
! !
! !> Retrieve Fermions operators to speed up the process below:
! P = left%operators%op("P")
! allocate(Cleft(Nso),Cright(Nso))
! do ispin=1,Nspin
!    do iorb=1,Norb
!       io = iorb+(ispin-1)*Norb
!       Cleft(io) = left%operators%op("C"//left%okey(iorb,ispin))
!       Cright(io) = right%operators%op("C"//right%okey(iorb,ispin))
!    enddo
! enddo


! print*,"Constructing filtered SB States"  
! allocate(AI(Nsb),BI(Nsb))
! allocate(AJ(Nsb),BJ(Nsb))
! do isb=1,size(sb_sector)
!    qn   = sb_sector%qn(index=isb)
!    AI(isb)%states = sb2block_states(left,right,sb_states,sb_sector%map(qn=qn),'left')
!    BI(isb)%states = sb2block_states(left,right,sb_states,sb_sector%map(qn=qn),'right')
!    do ispin=1,Nspin
!       dq = qnup ; if(ispin==2)dq=qndw
!       qm = qn - dq
!       if(.not.sb_sector%has_qn(qm))cycle
!       jsb = sb_sector%index(qn=qm)
!       AJ(jsb)%states = sb2block_states(left,right,sb_states,sb_sector%map(qn=qm),'left')
!       BJ(jsb)%states = sb2block_states(left,right,sb_states,sb_sector%map(qn=qm),'right')
!    enddo
! enddo

! print*,"Constructing H^{L,R}*: the filtered block Hamiltonians"
! t0 = omp_get_wtime()
! do isb=1,size(sb_sector)
!    qn   = sb_sector%qn(index=isb)
!    !
!    !> get: H*^L  and H*^R
!    Hleft(isb) = sp_filter(left%operators%op("H"),AI(isb)%states)
!    Hright(isb)= sp_filter(right%operators%op("H"),BI(isb)%states)
!    !
!    !> get A.x.B terms:
!    do ispin=1,Nspin
!       dq = qnup ; if(ispin==2)dq=qndw
!       qm   = qn - dq
!       if(.not.sb_sector%has_qn(qm))cycle
!       jsb = sb_sector%index(qn=qm)
!       !
!       do iorb=1,Norb
!          do jorb=1,Norb
!             io = iorb+(ispin-1)*Norb
!             jo = jorb+(ispin-1)*Norb
!             if(Hij(io,jo)==0d0)cycle
!             !
!             !> get: A = H(a,b)*[Cl(a,s)^+@P] .x. B = Cr(b,s)  + Row/Col Offsets
!             it=tMap(1,io,jo)
!             A(it,isb) = Hij(io,jo)*sp_filter(matmul(Cleft(io)%dgr(),P),AI(isb)%states,AJ(jsb)%states)
!             B(it,isb) = sp_filter(Cright(jo),BI(isb)%states,BJ(jsb)%states)
!             RowOffset(it,isb)=Offset(isb)           
!             ColOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
!             !
!             !> get H.c. + Row/Col Offsets
!             it=tMap(2,io,jo)
!             A(it,isb) = hconjg(A(tMap(1,io,jo),isb))
!             B(it,isb) = hconjg(B(tMap(1,io,jo),isb))
!             RowOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
!             ColOffset(it,isb)=Offset(isb)
!          enddo
!       enddo
!       !
!    enddo
! enddo
! print*,"Done:",omp_get_wtime() - t0
! print*,"######################################"
! print*,""



!OLD
! do isb=1,size(sb_sector)
!    qn   = sb_sector%qn(index=isb)
!    AIstates = sb2block_states(left,right,sb_states,sb_sector%map(qn=qn),'left')
!    BIstates = sb2block_states(left,right,sb_states,sb_sector%map(qn=qn),'right')
!    !
!    !> get: H*^L  and H*^R
!    Hleft(isb) = sp_filter(left%operators%op("H"),AIstates)
!    Hright(isb)= sp_filter(right%operators%op("H"),BIstates)
!    !
!    !
!    !> get A.x.B terms:
!    do ispin=1,Nspin
!       dq = qnup ; if(ispin==2)dq=qndw
!       qm   = qn - dq
!       if(.not.sb_sector%has_qn(qm))cycle
!       AJstates = sb2block_states(left,right,sb_states,sb_sector%map(qn=qm),'left')
!       BJstates = sb2block_states(left,right,sb_states,sb_sector%map(qn=qm),'right')
!       !
!       do iorb=1,Norb
!          io = iorb+(ispin-1)*Norb
!          do jorb=1,Norb
!             jo = jorb+(ispin-1)*Norb
!             if(Hij(io,jo)==0d0)cycle
!             !> get: A = H(a,b)*[Cl(a,s)^+@P] .x. B = Cr(b,s)  + Row/Col Offsets
!             it=tMap(1,io,jo)
!             A(it,isb) = Hij(io,jo)*sp_filter(matmul(Cleft(io)%dgr(),P),AIstates,AJstates)
!             B(it,isb) = sp_filter(Cright(jo),BIstates,BJstates)
!             RowOffset(it,isb)=Offset(isb)           
!             ColOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
!             !
!             !> get H.c. + Row/Col Offsets
!             it=tMap(2,io,jo)
!             A(it,isb) = hconjg(A(tMap(1,io,jo),isb))
!             B(it,isb) = hconjg(B(tMap(1,io,jo),isb))
!             RowOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
!             ColOffset(it,isb)=Offset(isb)
!          enddo
!       enddo
!       !
!    enddo
! enddo
! print*,"Done:",omp_get_wtime() - t0
! print*,"######################################"
! print*,""







! print*,""
! print*,""
! print*,"######################################"
! print*,"   o + o"
! print*,"######################################"
! print*,""
! print*,""

! print*,"_o+o_: with QN"
! left = block(dot)
! right = block(dot)
! call get_sb_states_(left,right,sb_states,sb_sector)
! spHsb  =  h2_model(left,right,sb_states)  + &
!      sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
!      sp_kron(id(left%dim),right%operators%op("H"),sb_states)

! allocate(Evals(Neigen))
! allocate(Evecs(size(sb_states),Neigen))
! call sp_eigh(sb_HxV,evals,evecs,&
!      5*Neigen,&
!      500,&
!      tol=1d-12,&
!      iverbose=.false.)
! do i=1,2
!    print*,i,Evals(i)/2/left%length/Norb
! enddo
! deallocate(evals,evecs)
! print*,""
