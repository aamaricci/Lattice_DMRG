program testEDkron
  USE SCIFOR
  USE DMRG, id=>sp_eye
  USE OMP_LIB, only: omp_get_wtime
  implicit none

  integer                                        :: Nso,tNso
  character(len=64)                              :: finput
  integer                                        :: i,j,unit,iorb,Nsb,k,l,r,m,arow,acol,brow,bcol,jorb,ispin,it,ip,SUN,jsb
  character(len=1)                               :: DMRGtype
  real(8)                                        :: target_qn(1)
  real(8)                                        :: aval,bval
  type(site)                                     :: Dot
  type(sparse_matrix)                            :: C,N,Cl,Cr,P
  type(sparse_matrix),allocatable,dimension(:)   :: Hleft,Hright,Hsb
  type(sparse_matrix),allocatable,dimension(:,:) :: A,B ![nz{Nspin*Norb*Norb},Nsb]
  integer,dimension(:,:),allocatable             :: RowOffset,ColOffset
  real(8),dimension(:,:),allocatable             :: Hloc,Hij
  real(8),dimension(:),allocatable               :: sb_qn,qn,qm
  type(block)                                    :: my_block,dimer,trimer,dimerL,dimerR,trimerL,trimerR,left,right,sys
  type(sparse_matrix)                            :: spHsb,spH
  real(8)                                        :: gs_energy
  integer                                        :: m_sb,isb,Ncv
  integer                                        :: model_d=4
  real(8),dimension(:,:),allocatable             :: Hmatrix,Evecs,Rho
  real(8),dimension(:),allocatable               :: Evals,Vec,Hvec,Hvec_
  integer,dimension(:),allocatable               :: sb_states,sb_map,Dls,Drs,Offset
  type(tstates),dimension(:),allocatable         :: Ai,Aj,Bi,Bj
  type(sectors_list)                             :: sb_sector
  integer                                        :: Neigen=2
  integer                                        :: current_L
  real(8),dimension(:),allocatable               :: current_target_QN
  real(8),dimension(1)                           :: dq=1d0
  real                                           :: t0
  integer,dimension(:,:,:),allocatable           :: tMap

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(SUN,"SUN",finput,default=2,&
       comment="Spin SU(N) value. 2=> spin 1/2, 3=> spin 1")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)



  dot = spin_site(sun=SUN)
  call dot%show()


  !Init DMRG
  call init_dmrg(heisenberg_1d_model,ModelDot=Dot)
  target_qn = DMRG_qn
  !


  print*,""
  print*,""
  print*,"######################################"
  print*,"   o + o"
  print*,"######################################"
  print*,""
  print*,""

  print*,"_o+o_: with QN"
  left = block(dot)
  right = block(dot)
  call get_sb_states_(left,right,sb_states,sb_sector)
  print*,size(sb_states)
  spHsb  =  heisenberg_1d_model(left,right,sb_states)  + &
       sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
       sp_kron(id(left%dim),right%operators%op("H"),sb_states)

  Neigen=1
  Ncv = min(size(sb_states),5*Neigen)
  allocate(Evals(Neigen))
  allocate(Evecs(size(sb_states),Neigen))
  call sp_eigh(sb_HxV,evals,evecs,&
       Ncv,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,Neigen
     print*,i,Evals(i)/2/left%length
  enddo
  deallocate(evals,evecs)
  print*,""






  print*,""
  print*,""
  print*,"######################################"
  print*,"   o->o + o<-o"
  print*,"######################################"
  print*,""
  print*,""

  print*,"_o->o++o<-o_: with QN"
  left  = block(dot)
  call enlarge_block_(left,dot,grow='left')
  right = block(dot)
  call enlarge_block_(right,dot,grow='right')

  print*,"Get SB states"
  t0 = omp_get_wtime()
  call get_sb_states_(left,right,sb_states,sb_sector)
  print*,omp_get_wtime() - t0






  print*,"######################################"
  print*,"Set up new solution method: [H*].v --> \psi"
  Nsb = size(sb_sector)
  print*,"There are Nsb_states:",Nsb

  Hij = hmodel_user(left,right)
  dq = 1d0


  print*,"Constructing \sum_a=1,tNso A*.B*"
  tNso = 3
  allocate(tMap(3,1,1))
  it = 0
  do i=1,tNso
     it = it+1
     tMap(i,1,1)=it
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


  print*,"Constructing filtered SB States"  
  allocate(AI(Nsb),BI(Nsb))
  allocate(AJ(Nsb),BJ(Nsb))
  do isb=1,size(sb_sector)
     qn   = sb_sector%qn(index=isb)
     AI(isb)%states = sb2block_states(left,right,sb_states,sb_sector%map(qn=qn),'left')
     BI(isb)%states = sb2block_states(left,right,sb_states,sb_sector%map(qn=qn),'right')
     qm   = qn - dq
     if(.not.sb_sector%has_qn(qm))cycle
     jsb = sb_sector%index(qn=qm)
     AJ(jsb)%states = sb2block_states(left,right,sb_states,sb_sector%map(qn=qm),'left')
     BJ(jsb)%states = sb2block_states(left,right,sb_states,sb_sector%map(qn=qm),'right')
  enddo

  !


  allocate(A(tNso,Nsb),B(tNso,Nsb))
  allocate(Hleft(Nsb),Hright(Nsb))
  allocate(RowOffset(tNso,Nsb),ColOffset(tNso,Nsb))
  print*,"Constructing H^{L,R}*: the filtered block Hamiltonians"
  t0 = omp_get_wtime()
  do isb=1,size(sb_sector)
     qn   = sb_sector%qn(index=isb)
     !
     !> get: H*^L  and H*^R
     Hleft(isb) = sp_filter(left%operators%op("H"),AI(isb)%states)
     Hright(isb)= sp_filter(right%operators%op("H"),BI(isb)%states)
     !
     !> get: A = Jp*S_lz .x. B = S_rz + Row/Col Offsets
     it=tMap(1,1,1)
     A(it,isb) = Hij(1,1)*sp_filter(left%operators%op("S_z"),AI(isb)%states,AI(isb)%states)
     B(it,isb) = sp_filter(right%operators%op("S_z"),BI(isb)%states,BI(isb)%states)
     RowOffset(it,isb)=Offset(isb)           
     ColOffset(it,isb)=Offset(isb)
     !
     qm   = qn - dq
     if(.not.sb_sector%has_qn(qm))cycle
     jsb = sb_sector%index(qn=qm)
     !> get: A = Jp*S_l- .x. B = [S_r-]^+=S_r+ + Row/Col Offsets 
     it=tMap(2,1,1)
     A(it,isb) = Hij(2,2)*sp_filter(left%operators%op("S_p"),AI(isb)%states,AJ(jsb)%states)
     B(it,isb) = sp_filter(hconjg(right%operators%op("S_p")),BI(isb)%states,BJ(jsb)%states)
     RowOffset(it,isb)=Offset(isb)
     ColOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
     !> get H.c. + Row/Col Offsets
     it=tMap(3,1,1)
     A(it,isb) = hconjg(A(tMap(2,1,1),isb))
     B(it,isb) = hconjg(B(tMap(2,1,1),isb))
     RowOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
     ColOffset(it,isb)=Offset(isb)
  enddo
  print*,"Done:",omp_get_wtime() - t0
  print*,"######################################"
  print*,""





  !> ONE STEP TEST:
  allocate(Vec(size(sb_states)))
  !call random_number(vec)
  vec = 1d0
  allocate(Hvec, source=vec)
  allocate(Hvec_, source=vec)
  !
  print*,"######################################"
  print*,"Test H.v OLD vs NEW method random vector v:"
  ! !Set spHsb to H^Lx1^R to be used in the default spHv (which only uses spHsb)
  t0 = omp_get_wtime()
  spHsb  =  Heisenberg_1d_model(left,right,sb_states) + &
       sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
       sp_kron(id(left%dim),right%operators%op("H"),sb_states)
  call sb_HxV(size(sb_states),vec,Hvec_)
  print*,omp_get_wtime() - t0
  !
  !Apply H to v new style
  t0 = omp_get_wtime()
  call sb_HxV_new(size(sb_states),vec,Hvec)
  print*,omp_get_wtime() - t0
  !
  do isb=1,size(sb_sector)
     print*,"Sector Dim:",Dls(isb)*Drs(isb)
     print*,"       Hv_old             ","       Hv_new             ","       Error              "
     do i=1,Dls(isb)*Drs(isb)
        k = i + Offset(isb)
        print*,Hvec_(k),Hvec(k),abs(Hvec(k)-Hvec_(k))
     enddo
  enddo
  print*,"Any error:",any(abs(Hvec-Hvec_)>1d-12)
  print*,"######################################"
  print*,""
  print*,""




  !Old Style solution:
  print*,"Old style solution: spH --> H.v --> \psi"
  print*,"######################################"
  t0 = omp_get_wtime()
  spHsb  =  Heisenberg_1d_model(left,right,sb_states)   + &
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
     print*,i,Evals(i)/2/left%length
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
     print*,i,Evals(i)/2/left%length
  enddo
  deallocate(evals,evecs)
  print*,""
  print*,omp_get_wtime() - t0
















contains



  function dim_sector(self,qn) result(dim)
    type(sectors_list)   :: self
    real(8),dimension(:) :: qn
    integer              :: dim
    dim = 0
    if(.not.self%has_qn(qn))return
    dim =  size(self%map(qn=qn))
  end function dim_sector



  !This is user defined Function to be passed to the SYSTEM
  function hmodel_user(left,right) result(Hlr)
    type(block)                        :: left,right
    real(8),dimension(:,:),allocatable :: Hlr
    if(allocated(Hlr))deallocate(Hlr)
    allocate(Hlr(Nspin*Norb,Nspin*Norb))
    !
    !workout local part, like random local field
    !if(left%Dim==1 AND right%Dim==1) then operate over local H
    !if(left%Dim==1 OR right%Dim==1) then operate over local H
    Hlr(1,1) = Jp
    Hlr(2,2) = Jx/2d0
  end function hmodel_user



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


  function h2_model(left,right) result(H2)
    type(block)              :: left
    type(block)              :: right
    type(sparse_matrix)      :: Sl(Nspin)![Sz,Sp]
    type(sparse_matrix)      :: Sr(Nspin)![Sz,Sp]
    type(sparse_matrix)      :: H2
    integer,dimension(2)     :: Hdims
    real(8),dimension(Nspin) :: g
    !
    !Hij is shared:
    Hij = hmodel_user(left,right)
    !
    !> Get H2 dimensions:
    Hdims = shape(left%operators)*shape(right%operators)
    call H2%init(Hdims(1),Hdims(2))

    !>Retrieve operators:
    do ispin=1,Nspin
       Sl(ispin) = left%operators%op("S"//left%okey(0,ispin))
       Sr(ispin) = right%operators%op("S"//right%okey(0,ispin))
    enddo
    !
    !
    !>Build H2:
    H2 = H2 + Hij(1,1)*(Sl(1).x.Sr(1)) + &
         Hij(2,2)*(Sl(2).x.Sr(2)%dgr())+ &
         Hij(2,2)*(Sl(2)%dgr().x.Sr(2))
    !
    !> Free memory
    do ispin=1,Nspin
       call Sl(ispin)%free
       call Sr(ispin)%free
    enddo
  end function h2_model





  !-----------------------------------------------------------------!
  ! Purpose: build the list of states compatible with the specified
  ! quantum numbers
  !-----------------------------------------------------------------!
  subroutine get_sb_states_(sys,env,sb_states,sb_sector)
    type(block)                      :: sys,env    
    integer,dimension(:),allocatable :: sb_states
    type(sectors_list)               :: sb_sector
    integer                          :: isys,ienv
    integer                          :: i,j,istate
    real(8),dimension(:),allocatable :: sys_qn,env_qn
    integer,dimension(:),allocatable :: sys_map,env_map

    !
    current_L         = sys%length + env%length
    current_target_QN = int(target_qn*current_L*Norb)
    write(LOGfile,"(A22,I12)")"SuperBlock Length = ",current_L
    write(LOGfile,"(A22,"//str(size(current_target_QN))//"F12.7)")"Target_QN = ",current_target_QN
    write(LOGfile,"(A22,G12.7)")"Total         = ",sum(current_target_QN)
    write(LOGfile,"(A22,G12.7)")"Filling       = ",sum(current_target_QN)/current_L
    write(LOGfile,"(A22,G12.7)")"Filling/Norb  = ",sum(current_target_QN)/current_L/Norb
    if(allocated(sb_states))deallocate(sb_states)
    !
    call sb_sector%free()
    do isys=1,size(sys%sectors(1))
       sys_qn  = sys%sectors(1)%qn(index=isys)
       env_qn  = current_target_qn - sys_qn
       if(.not.env%sectors(1)%has_qn(env_qn))cycle
       !
       sys_map = sys%sectors(1)%map(qn=sys_qn)!;print*,size(sys_map)
       env_map = env%sectors(1)%map(qn=env_qn)!;print*,size(env_map)
       !
       do i=1,size(sys_map)
          do j=1,size(env_map)
             istate=env_map(j) + (sys_map(i)-1)*env%Dim
             call append(sb_states, istate)
             call sb_sector%append(qn=sys_qn,istate=size(sb_states))
          enddo
       enddo
    enddo
    !
  end subroutine get_sb_states_





  function sb2block_states(sys,env,sb_states,sb_map,label) result(states)
    type(block)                      :: sys,env    
    integer,dimension(:)             :: sb_states,sb_map
    character(len=*)                 :: label
    integer,dimension(:),allocatable :: tmp,states
    integer :: i,istate,l,r,isb
    !
    if(allocated(states))deallocate(states)
    !
    select case(to_lower(str(label)))
    case("left","l","sys","s")
       do i=1,size(sb_map)
          istate = sb_states(sb_map(i))
          l = (istate-1)/env%Dim+1
          call append(tmp,l)
       enddo
    case("right","r","env","e")
       do i=1,size(sb_map)
          istate = sb_states(sb_map(i))
          r = mod(istate,env%Dim);if(r==0)r=env%Dim
          call append(tmp,r)
       enddo
    end select
    allocate(states, source=uniq(tmp))
  end function sb2block_states



  function heisenberg_1d_model(left,right,states) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    integer,dimension(:),optional :: states
    type(sparse_matrix)           :: Sz1,Sp1
    type(sparse_matrix)           :: Sz2,Sp2
    type(sparse_matrix)           :: Hp,Hz,H2
    integer,dimension(2)                  :: Hdims


    Sz1 = left%operators%op("S_z")
    Sp1 = left%operators%op("S_p")
    Sz2 = right%operators%op("S_z")
    Sp2 = right%operators%op("S_p")
    !
    !> Get H2 dimensions:
    Hdims = shape(left%operators)*shape(right%operators)
    if(present(states))Hdims = [size(states),size(states)]
    !
    !>Build H2:
    call H2%init(Hdims(1),Hdims(2))
    if(present(states))then
       H2 = H2 + Jp*sp_kron(Sz1,Sz2,states) + Jx/2d0*sp_kron(Sp1,Sp2%dgr(),states) + Jp/2d0*sp_kron(Sp1%dgr(),Sp2,states)
    else
       H2 = H2 + Jp*(Sz1.x.Sz2) + Jp/2d0*(Sp1.x.Sp2%dgr()) + Jx/2d0*(Sp1%dgr().x.Sp2)
    endif
    call Sz1%free()
    call Sp1%free()
    call Sz2%free()
    call Sp2%free()
  end function Heisenberg_1d_Model



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







  ! function sp_restricted_kron(A,B,states) result(AxB)
  !   type(sparse_matrix), intent(in) :: A,B
  !   integer,dimension(:),intent(in) :: states
  !   type(sparse_matrix)             :: AxB,Ap,Bp
  !   integer                         :: i,icol,j,k,kcol,l,istate,jstate
  !   integer                         :: indx_row,indx_col
  !   real(8)                         :: val,Aval,Bval
  !   !
  !   call AxB%free()
  !   call AxB%init(size(states),size(states))
  !   !
  !   print*,shape(A)
  !   print*,shape(B)
  !   do istate = 1,size(states)
  !      indx_row=states(istate)
  !      i = (indx_row-1)/B%Nrow+1
  !      k = mod(indx_row,B%Nrow);if(k==0)k=B%Nrow
  !      !
  !      do jstate=1,size(states)
  !         indx_col=states(jstate)
  !         j = (indx_col-1)/B%Ncol+1
  !         l = mod(indx_col,B%Ncol);if(l==0)l=B%Ncol
  !         if(&
  !              (.not.any(A%row(i)%cols==j))&
  !              .OR. &
  !              (.not.any(B%row(k)%cols==l)) )cycle
  !         Aval = A%get(i,j)
  !         Bval = B%get(k,l)
  !         val  = Aval*Bval
  !         !
  !         print*,indx_row," > ", indx_col,"  -  A row,col:",i,j,int(Aval)," - B row,col:",k,l,int(Bval)," < ", istate, jstate 
  !         !
  !         call append(AxB%row(istate)%vals,val)
  !         call append(AxB%row(istate)%cols,jstate)
  !         AxB%row(istate)%Size = AxB%row(istate)%Size + 1
  !      enddo
  !      !
  !   enddo
  !   print*,""
  ! end function sp_restricted_kron


end program testEDkron








