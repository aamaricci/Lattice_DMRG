MODULE DMRG_SUPERBLOCK_SETUP
  USE VARS_GLOBAL
  USE DMRG_CONNECT
  implicit none
  private

  !Memory pool:
  type(sparse_matrix),allocatable,dimension(:)   :: Hleft,Hright
  type(sparse_matrix),allocatable,dimension(:,:) :: A,B
  integer,dimension(:),allocatable               :: Dls,Drs,Offset
  integer,dimension(:,:),allocatable             :: RowOffset,ColOffset,isb2jsb
  logical,dimension(:,:),allocatable             :: IsHconjg
  !
  !Local variables
  integer                                        :: tNso
  integer                                        :: isb,jsb
  integer                                        :: i,j
  integer                                        :: ispin
  integer                                        :: iorb,jorb
  integer                                        :: io,jo


  public :: Setup_SuperBlock_Sparse
  public :: Setup_SuperBlock_Direct
  public :: spMatVec_sparse_main
  public :: spMatVec_direct_main
  !
  !-> used in DMRG_MEASURE to perform H|gs>
  public :: sb2block_states
  !
  !-> used in DMRG_SUPERBLOCK to deallocate
  public :: Hleft
  public :: Hright
  public :: A,B
  public :: Dls,Drs,Offset
  public :: RowOffset,ColOffset


contains



  !##################################################################
  !              SETUP THE SUPERBLOCK HAMILTONIAN
  !                      SPARSE MODE
  !    H^SB = H^L x 1^R  + 1^L x H^R + H^LR
  !    H^LR = sum_p A_p x B_p
  ! 
  ! * sparse: get the sparse global SB Hamiltonian spHsb
  !##################################################################
  !POSSIBLY INCLUDE MPI HERE... this is probably a less efficient version
  !note that only the SB Hamiltonian needs to be constructed in parallel form
  !the other operators (small) are stored by each cpu.
  !In principle one could store any sparse matrix in parallel and build H^SB as MPI too.
  subroutine Setup_SuperBlock_Sparse()
    integer                         :: m_left,m_right
    character(len=:),allocatable    :: type
    type(sparse_matrix)             :: H2
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> d733a73 (Updated code.)
    !
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Setup SB Sparse"
#endif
    !
<<<<<<< HEAD
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
>>>>>>> d733a73 (Updated code.)
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
    select case(to_lower(type))
    case default;stop "Setup_SuperBlock_Sparse ERROR: wrong left/right.Type"
    case ("spin","s")
       H2 = connect_spin_blocks(left,right,sb_states)
    case ("fermion","f,","electron","e")
       H2 = connect_fermion_blocks(left,right,sb_states)
    end select
    !
    spHsb = H2 & 
         + sp_kron(left%operators%op("H"),id(m_right),sb_states) &
         + sp_kron(id(m_left),right%operators%op("H"),sb_states)
    !
  end subroutine Setup_SuperBlock_Sparse





  

  !##################################################################
  !              SETUP THE SUPERBLOCK HAMILTONIAN
  !                       DIRECT MODE
  !    H^SB = H^L x 1^R  + 1^L x H^R + H^LR
  !    H^LR = sum_k A_k x B_k; k={q,p}, Q.N. + elements in H^LR
  ! 
  ! * direct: get vec{A},vec{B},vec{Hleft},vec{Hright}
  !   the sparse matrices which reconstruct H^SB above in terms of
  !   conserved QN blocks:
  !   + H^a=L,R = [H^a[q1], H^a[q2], ... , [H^a[qN]]
  !     with each H^a[qi] a given block of conserved QN qi
  !   + A or B  = [A[q1,q1'], ... , A[qN,qN']]
  !     with each A[q,q'] connectin the blocks of differend QN q and q'
  !##################################################################
  subroutine Setup_SuperBlock_Direct()
    character(len=:),allocatable                 :: type
    !
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> d733a73 (Updated code.)
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Setup SB Direct"
#endif
    !
<<<<<<< HEAD
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
>>>>>>> d733a73 (Updated code.)
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing right.H operator in the list"
    !
    type=left%type()
    if(type/=right%type())&
         stop "Setup_SuperBlock_Direct ERROR: left.Type != right.Type"
    !
    !
    select case(to_lower(type))
    case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.Type"
    case ("spin","s")
       call Setup_SuperBlock_Spin_Direct()
    case ("fermion","f,","electron","e")
       call Setup_SuperBlock_Fermion_Direct()
    end select
  end subroutine Setup_SuperBlock_Direct




  
  !##################################################################
  !                          SPIN CASE
  !##################################################################
  subroutine Setup_SuperBlock_Spin_Direct()
    integer                                      :: Nso,Nsb
    integer                                      :: it,isb,jsb
    real(8),dimension(:),allocatable             :: qn,qm
    type(tstates),dimension(:),allocatable       :: Ai,Aj,Bi,Bj
    real(8),dimension(:),allocatable             :: dq
    integer,dimension(:,:,:),allocatable         :: tMap
    type(sparse_matrix),allocatable,dimension(:) :: Sleft,Sright
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable        :: Hij
#else
    real(8),dimension(:,:),allocatable           :: Hij
#endif
    !
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> d733a73 (Updated code.)
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Setup SB Direct - spin"
#endif
    !
<<<<<<< HEAD
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
>>>>>>> d733a73 (Updated code.)
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing right.H operator in the list"
    !
    !
    !> GET THE USER DEFINED MODEL HAMILTONIAN PARAMETERS:
    ! Hij = Hmodel(left,right)
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij, source=HopH)
    !
    !
    Nso  = Nspin
    tNso = 3
    Nsb  = size(sb_sector)
    !
    !Massive allocation
    if(allocated(tMap))deallocate(tMap)
    if(allocated(Offset))deallocate(Offset)
    if(allocated(AI))deallocate(AI)
    if(allocated(BI))deallocate(BI)
    if(allocated(AJ))deallocate(AJ)
    if(allocated(BJ))deallocate(BJ)
    if(allocated(A))deallocate(A)
    if(allocated(B))deallocate(B)
    if(allocated(Hleft))deallocate(Hleft)
    if(allocated(Hright))deallocate(Hright)
    if(allocated(RowOffset))deallocate(RowOffset)
    if(allocated(ColOffset))deallocate(ColOffset)
    if(allocated(isb2jsb))deallocate(isb2jsb)
    if(allocated(IsHconjg))deallocate(IsHconjg)
    if(allocated(Sleft))deallocate(Sleft)
    if(allocated(Sright))deallocate(Sright)
    allocate(tMap(tNso,1,1))
    allocate(Offset(Nsb))
    allocate(AI(Nsb),BI(Nsb))
    allocate(AJ(Nsb),BJ(Nsb))
    allocate(A(tNso,Nsb),B(tNso,Nsb))
    allocate(Hleft(Nsb),Hright(Nsb))
    allocate(RowOffset(tNso,Nsb),ColOffset(tNso,Nsb))
    allocate(isb2jsb(tNso,Nsb))
    allocate(IsHconjg(tNso,Nsb))
    allocate(Sleft(Nspin),Sright(Nspin))

    !
    !Creating the sequence of operators A*_q, B*_q
    ! which decompose the term H^LR of the
    ! super-block Hamiltonian.
    it = 0
    do i=1,tNso
       it = it+1
       tMap(i,1,1)=it
    enddo
    !
    Offset=0
    do isb=1,Nsb
       qn      = sb_sector%qn(index=isb)
       Dls(isb)= sector_qn_dim(left%sectors(1),qn)
       Drs(isb)= sector_qn_dim(right%sectors(1),current_target_qn - qn)
       if(isb>1)Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
       if(MpiMaster)print*,int(qn), Drs(isb),Dls(isb),Dls(isb)*Drs(isb),Offset(isb)
    enddo
    if(MpiMaster)print*,""
    !
    dq=[1d0]
    do isb=1,Nsb
       qn             = sb_sector%qn(index=isb)
       AI(isb)%states = sb2block_states(qn,'left')
       BI(isb)%states = sb2block_states(qn,'right')
       qm             = qn - dq
       if(.not.sb_sector%has_qn(qm))cycle
       jsb            = sb_sector%index(qn=qm)
       AJ(jsb)%states = sb2block_states(qm,'left')
       BJ(jsb)%states = sb2block_states(qm,'right')
    enddo
    !
    do ispin=1,Nspin
       Sleft(ispin)   = left%operators%op("S"//left%okey(0,ispin))
       Sright(ispin)  = right%operators%op("S"//right%okey(0,ispin))
    enddo
    !
    isb2jsb=0
    do isb=1,Nsb
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
       qm  = qn
       jsb = isb
       RowOffset(it,isb)=Offset(isb)           
       ColOffset(it,isb)=Offset(isb)
       Isb2Jsb(it,isb)=jsb
       IsHconjg(it,isb)=.false.
       !
       if(MpiMaster)print*,"it,k,q,qn,qm:",it,isb,jsb,int(qn),int(qm),"A, B:",&
            shape(A(it,isb)),shape(B(it,isb)),"-",Dls(isb),Dls(jsb),Drs(isb),Drs(jsb)
       !Drs(isb),Dls(isb),Drs(jsb),Dls(jsb)
       !,B(it,isb)%Nrow,A(it,isb)%Nrow
       !
       !
       dq = [1d0]
       qm = qn - dq
       if(.not.sb_sector%has_qn(qm))cycle
       jsb = sb_sector%index(qn=qm)
       !
       !> get: A = Jxy*S_l- .x. B = [S_r-]^+=S_r+ + Row/Col Offsets 
       it=tMap(2,1,1)
       A(it,isb) = Hij(2,2)*sp_filter(Sleft(2),AI(isb)%states,AJ(jsb)%states)
       B(it,isb) = sp_filter(hconjg(Sright(2)),BI(isb)%states,BJ(jsb)%states)
       RowOffset(it,isb)=Offset(isb)
       ColOffset(it,isb)=Offset(jsb)
       Isb2Jsb(it,isb)=jsb
       IsHconjg(it,isb)=.false.
       !
       if(MpiMaster)print*,"it,k,q,qn,qm:",it,isb,jsb,int(qn),int(qm),"A, B:",&
            shape(A(it,isb)),shape(B(it,isb)),"-",Dls(isb),Dls(jsb),Drs(isb),Drs(jsb)
       !,B(it,isb)%Nrow,A(it,isb)%Nrow
       !       
       !> get H.c. + Row/Col Offsets
       it=tMap(3,1,1)
       A(it,isb) = hconjg(A(tMap(2,1,1),isb))
       B(it,isb) = hconjg(B(tMap(2,1,1),isb))
       RowOffset(it,isb)=Offset(jsb)
       ColOffset(it,isb)=Offset(isb)
       Isb2Jsb(it,isb)=jsb
       IsHconjg(it,isb)=.true.  !exchange jsb and isb
       !
       if(MpiMaster)print*,"it,k,q,qn,qm:",it,isb,jsb,int(qn),int(qm),"A, B:",&
            shape(A(it,isb)),shape(B(it,isb)),"-",Dls(jsb),Dls(isb),Drs(jsb),Drs(isb)
       !,B(it,isb)%Nrow,A(it,isb)%Nrow
       if(MpiMaster)print*,""
    enddo
  end subroutine Setup_SuperBlock_Spin_Direct











  !##################################################################
  !                        FERMION CASE
  !##################################################################
  subroutine Setup_SuperBlock_Fermion_Direct()
    integer                                      :: Nso,Nsb
    integer                                      :: it,isb,jsb
    real(8),dimension(:),allocatable             :: qn,qm
    type(tstates),dimension(:),allocatable       :: Ai,Aj,Bi,Bj
    real(8),dimension(:),allocatable             :: dq
    real(8),dimension(2),parameter               :: qnup=[1d0,0d0]
    real(8),dimension(2),parameter               :: qndw=[0d0,1d0]
    integer,dimension(:,:,:),allocatable         :: tMap
    type(sparse_matrix)                          :: P
    type(sparse_matrix),allocatable,dimension(:) :: Cleft,Cright
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable        :: Hij
#else
    real(8),dimension(:,:),allocatable           :: Hij
#endif
    !
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> d733a73 (Updated code.)
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Setup SB Direct - fermion"
#endif
    !
<<<<<<< HEAD
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
=======
>>>>>>> d733a73 (Updated code.)
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing right.H operator in the list"
    !
    !
    !> GET THE USER DEFINED MODEL HAMILTONIAN PARAMETERS:
    ! Hij = Hmodel(left,right)
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij, source=HopH)
    !
    !
    Nso  = Nspin*Norb
    tNso = 2*count(Hij/=zero)
    Nsb  = size(sb_sector)
    !
    !Massive allocation
    allocate(tMap(2,Nso,Nso))
    allocate(Offset(Nsb))
    allocate(AI(Nsb),BI(Nsb))
    allocate(AJ(Nsb),BJ(Nsb))
    allocate(A(tNso,Nsb),B(tNso,Nsb))
    allocate(Hleft(Nsb),Hright(Nsb))
    allocate(RowOffset(tNso,Nsb),ColOffset(tNso,Nsb))
    allocate(Cleft(Nso),Cright(Nso))
    !
    !Creating the sequence of operators A*_q, B*_q
    ! which decompose the term H^LR of the
    ! super-block Hamiltonian.
    it = 0
    do i=1,2
       do io=1,Nso
          do jo=1,Nso
             if(Hij(io,jo)==zero)cycle
             it = it+1
             tMap(i,io,jo)=it
          enddo
       enddo
    enddo
    !
    !
    Offset=0
    do isb=1,Nsb
       qn      = sb_sector%qn(index=isb)
       Dls(isb)= sector_qn_dim(left%sectors(1),qn)
       Drs(isb)= sector_qn_dim(right%sectors(1),current_target_qn - qn)
       if(isb>1)Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
    enddo
    !
    do isb=1,size(sb_sector)
       qn             = sb_sector%qn(index=isb)
       AI(isb)%states = sb2block_states(qn,'left')
       BI(isb)%states = sb2block_states(qn,'right')
       do ispin=1,Nspin
          dq = qnup   ; if(ispin==2)dq=qndw
          qm = qn - dq
          if(.not.sb_sector%has_qn(qm))cycle
          jsb            = sb_sector%index(qn=qm)
          AJ(jsb)%states = sb2block_states(qm,'left')
          BJ(jsb)%states = sb2block_states(qm,'right')
       enddo
    enddo
    !
    !
    P = left%operators%op("P")
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb+(ispin-1)*Norb
          Cleft(io)  = left%operators%op("C"//left%okey(iorb,ispin))
          Cright(io) = right%operators%op("C"//right%okey(iorb,ispin))
       enddo
    enddo
    !
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
                ColOffset(it,isb)=Offset(jsb)
                !
                !> get H.c. + Row/Col Offsets
                it=tMap(2,io,jo)
                A(it,isb) = hconjg(A(tMap(1,io,jo),isb))
                B(it,isb) = hconjg(B(tMap(1,io,jo),isb))
                RowOffset(it,isb)=Offset(jsb)
                ColOffset(it,isb)=Offset(isb)
             enddo
          enddo
          !
       enddo
    enddo
  end subroutine Setup_SuperBlock_Fermion_Direct













  !##################################################################
  !              SuperBlock MATRIX-VECTOR PRODUCTS
  !              using shared quantities in GLOBAL
  !##################################################################
  subroutine spMatVec_sparse_main(Nloc,v,Hv)
    integer                    :: Nloc
    integer                    :: i,j,jcol
#ifdef _CMPLX
    complex(8),dimension(Nloc) :: v
    complex(8),dimension(Nloc) :: Hv
    complex(8)                 :: val
#else
    real(8),dimension(Nloc)    :: v
    real(8),dimension(Nloc)    :: Hv
    real(8)                    :: val
#endif
    Hv=zero
    do i=1,Nloc
       matmul: do jcol=1, spHsb%row(i)%Size
          val = spHsb%row(i)%vals(jcol)
          j   = spHsb%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end subroutine spMatVec_sparse_main









  !##################################################################
  !              SuperBlock MATRIX-VECTOR PRODUCTS
  !              using shared quantities in GLOBAL
  !##################################################################
  subroutine spMatVec_direct_main(Nloc,v,Hv)
    integer                               :: Nloc
#ifdef _CMPLX
    complex(8),dimension(Nloc)            :: v
    complex(8),dimension(Nloc)            :: Hv
    complex(8)                            :: val
    complex(8)                            :: aval,bval
    complex(8),dimension(:,:),allocatable :: C,Ct
#else
    real(8),dimension(Nloc)               :: v
    real(8),dimension(Nloc)               :: Hv
    real(8)                               :: val
    real(8)                               :: aval,bval
    real(8),dimension(:,:),allocatable    :: C,Ct
#endif
    integer                               :: i,j,k,q,n
    integer                               :: ir,il,jr,jl,it
    integer                               :: ai,aj,bi,bj
    integer                               :: ia,ib,ic,ja,jb,jc
    logical                               :: isDiag
    !
    Hv=zero
    !> loop over all the SB sectors:
    sector: do k=1,size(sb_sector)
       !
       !> apply the 1^L x H^r
       !> this loop can be made parallel:
       do il=1,Dls(k)           !Fix the column il: v_il 
          !
          do ir=1,Drs(k)        !H^r.v_il
             i = ir + (il-1)*Drs(k) + offset(k)
             do jcol=1,Hright(k)%row(ir)%Size
                val = Hright(k)%row(ir)%vals(jcol)
                jr  = Hright(k)%row(ir)%cols(jcol)
                j   = jr + (il-1)*Drs(k) + offset(k)
                Hv(i) = Hv(i) + val*v(j)
             end do
          enddo
          !
       enddo
       !
       !> apply the H^L x 1^r
       !> this loop can be made parallel, provided an MPI_trasposition
       do ir=1,Drs(k)           !Fix the row ir: v_ir
          !
          do il=1,Dls(k)        !H^l.v_ir
             i = ir + (il-1)*Drs(k) + offset(k)
             do jcol=1,Hleft(k)%row(il)%Size
                val = Hleft(k)%row(il)%vals(jcol)
                jl  = Hleft(k)%row(il)%cols(jcol)
                j   = ir + (jl-1)*Drs(k) + offset(k)
                Hv(i) = Hv(i) + val*v(j)
             end do
          enddo
          !
       enddo
       !
       !> apply the term sum_k sum_it A_it(k).x.B_it(k)
       !Hv = (A.x.B)vec(V) --> (A.x.B).V  -> vec(B.V.A^T)
       !
       !  B.V.A^T : [B.Nrow,B.Ncol].[B.Ncol,A.Ncol].[A.Ncol,A.Nrow]
       !            [Dr(k),Dr(k')].[Dr(k'),Dl(k')].[Dl(k'),Dl(k)]
       !   C.A^T  : [B.Nrow,A.Ncol].[A.Ncol,A.Nrow]
       !(A.C^T)^T : [ [A.Nrow,A.Ncol].[A.Ncol,B.Nrow] ]^T
       !              [B.Nrow,A.Nrow] = vec(Hv)
       !
       !A(it,k).Nrow = DLk ;if(Hconjg) DLq
       !A(it,k).Ncol = DLq ;if(Hconjg) DLk
       !B(it,k).Nrow = DRk ;if(Hconjg) DRq
       !B(it,k).Ncol = DRq ;if(Hconjg) DRk
       !

       !
       do it=1,tNso
          q = isb2jsb(it,k)
          if(.not.A(it,k)%status.OR..not.B(it,k)%status)cycle
          !
          isDiag=.false. ; if(q==k)isDiag=.true.
          !
          allocate(C(B(it,k)%Nrow,A(it,k)%Ncol));C=0d0
          !
          !1. evaluate MMP: C = B.vec(V)
          !   \sum_bcol B(bi,bj)V_q(bj,aj)=C(bi,aj)
          !   j = bj+(aj-1)B.Ncol + ColOffset_q
          !   \sum_bcol B(bi,bj)v_q(j)=C(bi,aj)
          !
          !> this loop can be made parallel: A(it,k).Ncol == DLq or DLk
          do aj=1,A(it,k)%Ncol !DLq ;if(isHconjg(k))DLk
             !
             do bi=1,B(it,k)%Nrow
                if(B(it,k)%row(bi)%Size==0)cycle
                !
                do jb=1,B(it,k)%row(bi)%Size
                   bj   = B(it,k)%row(bi)%cols(jb)
                   val  = B(it,k)%row(bi)%vals(jb)
                   jc   = bj + (aj-1)*B(it,k)%Ncol
                   j    = jc + ColOffset(it,k)
                   C(bi,aj) = C(bi,aj) + val*v(j)
                enddo
                !
             enddo
          enddo
          !2. evaluate MMP: C.A^t
          !   \sum_aj C(bi,aj)A^t(aj,ai)
          !  =\sum_aj [A(ai,aj)C^t(aj,bi)]^T
          !
          !> this loop can be made parallel: B(it,k).Nrow == DRq or DRk
          !Ct = transpose(C)
          !
          do bi=1,B(it,k)%Nrow
             !
             do ai=1,A(it,k)%Nrow
                if(A(it,k)%row(ai)%Size==0)cycle
                ic =  bi + (ai-1)*B(it,k)%Nrow
                i  =  ic + RowOffset(it,k)
                !
                do ja=1,A(it,k)%row(ai)%Size
                   aj  = A(it,k)%row(ai)%cols(ja)
                   val = A(it,k)%row(ai)%vals(ja)
                   Hv(i) = Hv(i) + val*C(bi,aj)!Ct(aj,bi)
                enddo
                !
             enddo
          enddo
          deallocate(C)
       enddo
    enddo sector
  end subroutine spMatVec_direct_main










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




END MODULE DMRG_SUPERBLOCK_SETUP













!   subroutine spMatVec_direct_main_(Nloc,v,Hv)
!     integer                    :: Nloc
! #ifdef _CMPLX
!     complex(8),dimension(Nloc) :: v
!     complex(8),dimension(Nloc) :: Hv
!     complex(8)                 :: val
!     complex(8)                 :: aval,bval
! #else
!     real(8),dimension(Nloc)    :: v
!     real(8),dimension(Nloc)    :: Hv
!     real(8)                    :: val
!     real(8)                    :: aval,bval
! #endif
!     integer                    :: i,j,k,n
!     integer                    :: ir,il,jr,jl,it
!     integer                    :: arow,brow,acol,bcol,jcol
!     integer                    :: ia,ib,ic,ja,jb,jc,sum
!     !
!     Hv=zero
!     !> loop over all the SB sectors:
!     sector: do k=1,size(sb_sector)
!        !
!        !> apply the 1^L x H^r
!        do il=1,Dls(k)           !Fix the column il: v_il 
!           !
!           do ir=1,Drs(k)        !H^r.v_il
!              i = ir + (il-1)*Drs(k) + offset(k)
!              do jcol=1,Hright(k)%row(ir)%Size
!                 val = Hright(k)%row(ir)%vals(jcol)
!                 jr  = Hright(k)%row(ir)%cols(jcol)
!                 j   = jr + (il-1)*Drs(k) + offset(k)
!                 Hv(i) = Hv(i) + val*v(j)
!              end do
!           enddo
!           !
!        enddo
!        !
!        !> apply the H^L x 1^r
!        do ir=1,Drs(k)           !Fix the row ir: v_ir
!           !
!           do il=1,Dls(k)        !H^l.v_ir
!              i = ir + (il-1)*Drs(k) + offset(k)
!              do jcol=1,Hleft(k)%row(il)%Size
!                 val = Hleft(k)%row(il)%vals(jcol)
!                 jl  = Hleft(k)%row(il)%cols(jcol)
!                 j   = ir + (jl-1)*Drs(k) + offset(k)
!                 Hv(i) = Hv(i) + val*v(j)
!              end do
!           enddo
!           !
!        enddo
!        !
!        !> apply the term sum_k sum_it A_it(k).x.B_it(k)
!        do it=1,tNso          
!           if(.not.A(it,k)%status.OR..not.B(it,k)%status)cycle
!           do ic=1,A(it,k)%Nrow*B(it,k)%Nrow
!              arow = (ic-1)/B(it,k)%Nrow+1
!              brow = mod(ic,B(it,k)%Nrow);if(brow==0)brow=B(it,k)%Nrow
!              if(A(it,k)%row(arow)%Size==0.OR.B(it,k)%row(brow)%Size==0)cycle
!              i = ic + RowOffset(it,k)
!              do ja=1,A(it,k)%row(arow)%Size
!                 acol = A(it,k)%row(arow)%cols(ja)
!                 aval = A(it,k)%row(arow)%vals(ja)
!                 do jb=1,B(it,k)%row(brow)%Size
!                    bcol = B(it,k)%row(brow)%cols(jb)
!                    bval = B(it,k)%row(brow)%vals(jb)
!                    j = bcol+(acol-1)*B(it,k)%Ncol + ColOffset(it,k)
!                    Hv(i) = Hv(i) + aval*bval*v(j)
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo sector
!   end subroutine spMatVec_direct_main_







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
