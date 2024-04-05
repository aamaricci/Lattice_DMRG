MODULE SUPERBLOCK
  USE GLOBAL
  USE CONNECT
  implicit none
  private


  public :: sb_get_states
  public :: sb_diag
  public :: sb_build_Hv
  
  public :: spMatVec_sparse_main
  public :: spMatVec_direct_main
  
  integer :: i,j
  integer :: ispin
  integer :: iorb,jorb
  integer :: io,jo



contains



  !-----------------------------------------------------------------!
  ! Purpose: build the list of states compatible with the specified
  ! quantum numbers
  !-----------------------------------------------------------------!
  subroutine sb_get_states()!left,right,sb_states,sb_sector)
    ! type(block),intent(in)           :: left,right
    ! integer,dimension(:),allocatable :: sb_states
    ! type(sectors_list)               :: sb_sector
    integer                          :: ileft,iright
    integer                          :: i,j,istate
    real(8),dimension(:),allocatable :: left_qn,right_qn
    integer,dimension(:),allocatable :: left_map,right_map
    !
    call start_timer("Get SB states")
    !
    if(allocated(sb_states))deallocate(sb_states)
    !
    call sb_sector%free()
    !
    do ileft=1,size(left%sectors(1))
       left_qn   = left%sectors(1)%qn(index=ileft)
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
       call eta(ileft,size(left%sectors(1)))
    enddo
    !
    call stop_timer("")
    !
  end subroutine sb_get_states





  subroutine sb_diag()
    integer                            :: m_sb
    integer                            :: Nitermax,Neigen,Nblock
    real(8),dimension(:),allocatable   :: evals
    real(8),dimension(:,:),allocatable :: Hsb
    logical                            :: exist,lanc_solve

    
    m_sb = size(sb_states)
    if(m_sb==0)stop "sb_diag ERROR: size(sb_states)==0"

    !Set Lanczos params
    Neigen   = min(m_sb,lanc_neigen)
    Nitermax = min(m_sb,lanc_niter)
    Nblock   = min(m_sb,lanc_ncv_factor*Lanc_Neigen+lanc_ncv_add)
    !
    !Decide how to operate on H_sb
    lanc_solve  = .true.
    if(Lanc_Neigen==m_sb)lanc_solve=.false.
    if(m_sb <= lanc_dim_threshold)lanc_solve=.false.
    !
    !Allocate EigPairs
    if(allocated(gs_energy))deallocate(gs_energy)
    if(allocated(gs_vector))deallocate(gs_vector)
    allocate(gs_energy(Neigen))     ;gs_energy=0d0
    allocate(gs_vector(m_sb,Neigen));gs_vector=0d0
    !
    call start_timer()
    if(lanc_solve)then
       call sb_build_Hv()
       call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
            Nblock,&
            Nitermax,&
            tol=lanc_tolerance,&
            iverbose=.false.)
    else !use LAPACK
       call sb_build_Hv(Hsb)
       allocate(evals(m_sb))
       call eigh(Hsb,evals)
       gs_vector(:,1:Lanc_Neigen) = Hsb(:,1:Lanc_Neigen)
       gs_energy(1:Lanc_Neigen)   = evals(1:Lanc_Neigen)
       deallocate(Hsb,evals)
    endif
    call stop_timer("Diag H_sb")
  end subroutine sb_diag

  


  subroutine sb_build_Hv(Hmat)
    real(8),dimension(:,:),allocatable,optional :: Hmat
    integer                                     :: m_sb

    if(.not.allocated(sb_states))stop "build_Hv_superblock ERROR: sb_states not allocated"
    m_sb = size(sb_states)

    !IF PRESENT HMAT: get SB_H sparse > dump it to dense Hmat > return
    if(present(Hmat))then
       if(allocated(Hmat))deallocate(Hmat)
       allocate(Hmat(m_sb,m_sb));Hmat=0d0
       !Nullify HxV function pointer:
       spHtimesV_p => null()
       !
       !>Build Sparse Hsb:
       call start_timer("get H_sb")
       call Setup_SuperBlock_Sparse()
       call stop_timer("Done H_sb")
       !
       !Dump Hsb to dense matrix as required:
       call spHsb%dump(Hmat)
       return
    endif

    !Build SuperBLock HxV operation: stored or direct
    select case(sparse_H)
    case(.true.)
       call start_timer("get H_sb")
       call Setup_SuperBlock_Sparse()
       call stop_timer("Done H_sb")
       !
       !Set HxV function pointer:
       spHtimesV_p => spMatVec_sparse_main
       !
    case(.false.)
       call start_timer("get H_sb")
       call Setup_SuperBlock_Direct()
       call stop_timer("Done H_sb")
       !
       !Set HxV function pointer:
       spHtimesV_p => spMatVec_direct_main
       !
    end select

  end subroutine sb_build_Hv






  
  !##################################################################
  !              SETUP THE SUPERBLOCK HAMILTONIAN
  !    H^SB = H^L x 1^R  + 1^L x H^R + H^LR
  !    H^LR = sum_p A_p x B_p
  ! 
  ! * sparse: get the sparse global SB Hamiltonian spHsb
  ! * direct: get vec{A},vec{B},vec{Hleft},vec{Hright}
  !   the sparse matrices which reconstruct H^SB above in terms of
  !   conserved QN blocks:
  !   + H^a=L,R = [H^a[q1], H^a[q2], ... , [H^a[qN]]
  !     with each H^a[qi] a given block of conserved QN qi
  !   + A or B  = [A[q1,q1'], ... , A[qN,qN']]
  !     with each A[q,q'] connectin the blocks of differend QN q and q'
  !##################################################################
  subroutine Setup_SuperBlock_Sparse()
    integer                         :: m_left,m_right
    character(len=:),allocatable    :: type
    type(sparse_matrix)             :: H2
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Sparse ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Sparse ERROR: Missing right.H operator in the list"
    !
    type=left%type()
    if(type/=right%type())&
         stop "Setup_SuperBlock_Sparse ERROR: left.Type != right.Type"
    !
    m_left = left%dim
    m_right= right%dim
    !
    select case(type)
    case default;stop "Setup_SuperBlock_Sparse ERROR: wrong left/right.Type"
    case ("spin","s")
       H2 = connect_spin_blocks(left,right,sb_states)
    case ("fermion","f,","electron","e")
       H2 = connect_fermion_blocks(left,right,sb_states)
    end select
    !
    spHsb = H2 + &
         sp_kron(left%operators%op("H"),id(m_right),sb_states) + &
         sp_kron(id(m_left),right%operators%op("H"),sb_states)  
    !
  end subroutine Setup_SuperBlock_Sparse




  subroutine Setup_SuperBlock_Direct()
    integer                                      :: Nso,Nsb
    integer                                      :: it,isb,jsb
    real(8),dimension(:),allocatable             :: qn,qm
    type(tstates),dimension(:),allocatable       :: Ai,Aj,Bi,Bj
    real(8),dimension(:),allocatable             :: dq
    real(8),dimension(2),parameter               :: qnup=[1d0,0d0],qndw=[0d0,1d0]
    integer,dimension(:,:,:),allocatable         :: tMap
    type(sparse_matrix)                          :: P
    type(sparse_matrix),allocatable,dimension(:) :: Cleft,Cright
    type(sparse_matrix),allocatable,dimension(:) :: Sleft,Sright
    real(8),dimension(:,:),allocatable           :: Hij
    integer                                      :: m_left,m_right
    character(len=:),allocatable                 :: type
    !
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing right.H operator in the list"
    !
    type=left%type()
    if(type/=right%type())&
         stop "Setup_SuperBlock_Direct ERROR: left.Type != right.Type"
    !
    m_left = left%dim
    m_right= right%dim

    !Hij is shared:
    print*,"Hij:"
    Hij = Hmodel(left,right)
    print*,shape(Hij)
    !
    !
    select case(type)
    case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.Type"
    case ("spin","s")
       Nso = Nspin
    case ("fermion","f,","electron","e")
       Nso = Nspin*Norb
    end select
    !
    print*,"######################################"
    print*,"Set up new solution method: [H*].v --> \psi"
    Nsb = size(sb_sector)
    print*,"There are Nsb_states:",Nsb


    !    
    !Creating the sequence of operators A*_q, B*_q which decompose the term H^LR of the
    ! super-block Hamiltonian.
    print*,"Constructing \sum_a=1,tNso A*.B*"
    select case(type)
    case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.Type"
    case ("spin","s")
       tNso = 3
       allocate(tMap(3,1,1))
       it = 0
       do i=1,tNso
          it = it+1
          tMap(i,1,1)=it
       enddo
    case ("fermion","f,","electron","e")
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
    end select
    print*,"There are tNso non-zero elements ",tNso



    print*,"Constructing Offset and D_l,r"
    allocate(Dls(Nsb),Drs(Nsb),Offset(Nsb))
    Offset=0
    do isb=1,Nsb
       qn   = sb_sector%qn(index=isb)
       Dls(isb)= sector_qn_dim(left%sectors(1),qn)
       Drs(isb)= sector_qn_dim(right%sectors(1),current_target_qn - qn)
       if(isb>1)Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
    enddo


    print*,"Constructing filtered SB States"
    allocate(AI(Nsb),BI(Nsb))
    allocate(AJ(Nsb),BJ(Nsb))
    select case(type)
    case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.Type"
    case ("spin","s")
       dq=[1d0]
       do isb=1,size(sb_sector)
          qn   = sb_sector%qn(index=isb)
          AI(isb)%states = sb2block_states(qn,'left')
          BI(isb)%states = sb2block_states(qn,'right')
          qm   = qn - dq
          if(.not.sb_sector%has_qn(qm))cycle
          jsb = sb_sector%index(qn=qm)
          AJ(jsb)%states = sb2block_states(qm,'left')
          BJ(jsb)%states = sb2block_states(qm,'right')
       enddo
    case ("fermion","f,","electron","e")
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
    end select



    allocate(A(tNso,Nsb),B(tNso,Nsb))
    allocate(Hleft(Nsb),Hright(Nsb))
    allocate(RowOffset(tNso,Nsb),ColOffset(tNso,Nsb))
    select case(type)
    case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.Type"
    case ("spin","s")
       allocate(Sleft(Nspin),Sright(Nspin))
       do ispin=1,Nspin
          Sleft(ispin)  = left%operators%op("S"//left%okey(0,ispin))
          Sright(ispin) = right%operators%op("S"//right%okey(0,ispin))
       enddo
       !
       print*,"Constructing H^{L,R}*: the filtered block Hamiltonians"
       do isb=1,size(sb_sector)
          qn = sb_sector%qn(index=isb)
          !
          !> get: H*^L  and H*^R
          Hleft(isb) = sp_filter(left%operators%op("H"),AI(isb)%states)
          Hright(isb)= sp_filter(right%operators%op("H"),BI(isb)%states)
          !
          !> get: A = Jp*S_lz .x. B = S_rz + Row/Col Offsets       
          it=tMap(1,1,1)
          A(it,isb) = Hij(1,1)*sp_filter(Sleft(1),AI(isb)%states,AI(isb)%states)
          B(it,isb) = sp_filter(Sright(1),BI(isb)%states,BI(isb)%states)
          RowOffset(it,isb)=Offset(isb)           
          ColOffset(it,isb)=Offset(isb)
          !
          dq = [1d0]
          qm = qn - dq
          if(.not.sb_sector%has_qn(qm))cycle
          jsb = sb_sector%index(qn=qm)
          !> get: A = Jp*S_l- .x. B = [S_r-]^+=S_r+ + Row/Col Offsets 
          it=tMap(2,1,1)
          A(it,isb) = Hij(2,2)*sp_filter(Sleft(2),AI(isb)%states,AJ(jsb)%states)
          B(it,isb) = sp_filter(hconjg(Sright(2)),BI(isb)%states,BJ(jsb)%states)
          RowOffset(it,isb)=Offset(isb)
          ColOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
          !> get H.c. + Row/Col Offsets
          it=tMap(3,1,1)
          A(it,isb) = hconjg(A(tMap(2,1,1),isb))
          B(it,isb) = hconjg(B(tMap(2,1,1),isb))
          RowOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
          ColOffset(it,isb)=Offset(isb)
       enddo


    case ("fermion","f,","electron","e")
       allocate(Cleft(Nso),Cright(Nso))
       P = left%operators%op("P")
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb+(ispin-1)*Norb
             Cleft(io) = left%operators%op("C"//left%okey(iorb,ispin))
             Cright(io) = right%operators%op("C"//right%okey(iorb,ispin))
          enddo
       enddo

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
             qm = qn - dq
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
    end select
    !
    print*,"Done:"
    print*,"######################################"
    print*,""
  end subroutine Setup_SuperBlock_Direct



  !##################################################################
  !              SuperBlock MATRIX-VECTOR PRODUCTS
  !              using shared quantities in GLOBAL
  !##################################################################
  subroutine spMatVec_sparse_main(Nloc,v,Hv)
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
  end subroutine spMatVec_sparse_main



  subroutine spMatVec_direct_main(Nloc,v,Hv)
    integer                            :: Nloc
    real(8),dimension(Nloc)            :: v
    real(8),dimension(Nloc)            :: Hv
    real(8)                            :: val
    integer                            :: i,j,k,n
    integer                            :: ir,il,jr,jl,it
    integer                            :: arow,brow,acol,bcol,jcol
    integer                            :: ia,ib,ic,ja,jb,jc
    real(8)                            :: aval,bval
    !
    Hv=zero
    !> loop over all the SB sectors:
    sector: do k=1,size(sb_sector)
       !
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
       !
       !> apply the term sum_k sum_it A_it(k).x.B_it(k)
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
    enddo sector
  end subroutine spMatVec_direct_main







  !##################################################################
  !              RETURN LEFT.states or RIGHT.states 
  !              contributing to the SUPERBLOCK with
  !              a given QN 
  !##################################################################
  function sb2block_states(q,label) result(states)
    real(8),dimension(:)             :: q
    character(len=*)                 :: label
    integer,dimension(:),allocatable :: tmp,states,sb_map
    integer                          :: i,istate,l,r,isb
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

END MODULE SUPERBLOCK
