MODULE DMRG_MEASURE
  USE SCIFOR, only: to_lower
  USE VARS_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK_COMMON
  implicit none
  private


  !Measuring:
  public :: Init_Measure_DMRG
  public :: Measure_Op_DMRG
  public :: NuMeasure_Op_DMRG
  public :: Build_Op_DMRG
  public :: Advance_Op_DMRG
  public :: Advance_Corr_DMRG
  public :: Average_Op_DMRG
  public :: NuAverage_Op_DMRG
  public :: dmrg_write


  integer                                      :: Nso,Nsb
  real(8),dimension(:),allocatable             :: qn,qm
  real(8),dimension(:),allocatable             :: dq
  type(sparse_matrix),allocatable,dimension(:) :: Oleft,Oright,Olist
  type(tstates),dimension(:),allocatable       :: Ai,Bi

  logical ::     measure_status=.false.

contains

  subroutine Init_Measure_dmrg
    !
    print*,"Initialize DMRG measure:"
    !
    Nsb  = size(sb_sector)
    allocate(Dls(Nsb),Drs(Nsb),Offset(Nsb))
    allocate(AI(Nsb),BI(Nsb))
    Offset=0
    do isb=1,Nsb
       qn   = sb_sector%qn(index=isb)
       Dls(isb)= sector_qn_dim(left%sectors(1),qn)
       Drs(isb)= sector_qn_dim(right%sectors(1),current_target_qn - qn)
       if(isb>1)Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
       AI(isb)%states = sb2block_states(qn,'left')
       BI(isb)%states = sb2block_states(qn,'right')
    enddo
    !
    measure_status=.true.
  end subroutine Init_Measure_dmrg



  subroutine End_measure_DMRG
    if(allocated(Dls))deallocate(Dls)
    if(allocated(Drs))deallocate(Drs)
    if(allocated(Offset))deallocate(Offset)
    if(allocated(Olist))deallocate(Olist)
    if(allocated(Ai))deallocate(Ai)
    if(allocated(Bi))deallocate(Bi)
    measure_status=.false.
  end subroutine End_measure_DMRG



  !##################################################################
  !              MEASURE LOCAL OPERATOR
  !##################################################################
  subroutine NuMeasure_Op_dmrg(Op,pos,file,avOp)
    type(sparse_matrix),intent(in)            :: Op
    integer,dimension(:)                      :: pos
    character(len=*)                          :: file
    real(8),dimension(:),allocatable,optional :: avOp
    !
    real(8)                                   :: val
    character(len=1)                          :: label
    integer                                   :: i,ipos,L,R,N,Np
    type(sparse_matrix)                       :: Oi
    !
    !
    suffix=label_DMRG('u')
    !
    Np= size(pos)
    !
    L = left%length-1           !the len of the last block used to create the SB->\psi
    R = right%length-1
    N = L+R
    !
    if(present(avOp))then
       if(allocated(avOp))deallocate(avOp)
       allocate(avOp(Np))
    endif

    !
    call Init_measure_dmrg()
    call start_timer()
    do i=1,Np
       ipos=pos(i)
       Oi = Build_Op_dmrg(Op,ipos)
       Oi = Advance_Op_dmrg(Oi,ipos)
       val= NuAverage_Op_dmrg(Oi,ipos)
       print*,ipos,val
       if(present(avOp))avOp(i)=val
       call write_user(trim(file),[val],x=ipos)
       ! call progress(i,Np)
       call Oi%free()
    enddo
    call stop_timer("Done "//str(file))
    call End_measure_dmrg()
  end subroutine NuMeasure_Op_dmrg


  subroutine Measure_Op_dmrg(Op,file,ref,avOp)
    type(sparse_matrix),intent(in)            :: Op
    character(len=*)                          :: file
    character(len=1)                          :: label
    real(8) :: ref
    real(8),dimension(:),allocatable,optional :: avOp
    real(8)                                   :: val
    integer                                   :: it,i,L,R,N,j,pos,dims(2)
    type(sparse_matrix)                       :: Oi,U,Psi,Ok,I_R,I_L
    !

    suffix=label_DMRG('u')
    !
    L = left%length-1           !the len of the last block used to create the SB->\psi
    R = right%length-1
    N = L+R
    !
    if(present(avOp))then
       if(allocated(avOp))deallocate(avOp)
       allocate(avOp(N))
    endif
    !
    call start_timer()
    do pos=1,N
       Oi = Build_Op_dmrg(Op,pos)
       Oi = Advance_Op_dmrg(Oi,pos)
       val= Average_Op_dmrg(Oi,pos)
       if(present(avOp))avOp(pos)=val
       call write_user(trim(file),[val],x=pos)
       call Oi%free()
       call progress(pos,N)
    enddo
    call stop_timer("Done "//str(file))




    U =  right%omatrices%op(index=R)
    dims=shape(U)
    ! I_R = Id(dot%dim).x.Id(dims(2))!(matmul(U%t(),U))
    I_R = Id(dot%dim*dims(2))
    U =  left%omatrices%op(index=R)
    dims=shape(U)
    ! I_L = Id(dims(2)).x.Id(dot%dim)!(matmul(U%t(),U)).x.Id(dot%dim)
    I_L = Id(dims(2)*dot%dim)!(matmul(U%t(),U)).x.Id(dot%dim)


    print*,"size(GSpsi)",size(gs_vector(:,1))



    print*,""
    print*,"- - - - - - - - - - - - - - - - -"
    print*," METHOD 1: O.x.I - I.x.O full"
    print*,"- - - - - - - - - - - - - - - - -"
    print*,""


    print*,""
    print*,"pos=1"
    print*,""
    print*,"Method <psi|O.x.I_R|psi>"
    print*,shape(I_R)
    pos=1
    Oi = Build_Op_Dmrg(Op,pos)
    do it=1,L
       U  = left%omatrices%op(index=it)
       Oi = (matmul(matmul(U%t(),Oi),U)).x.Id(dot%dim)
    enddo
    Oi = sp_kron(Oi,I_R,sb_states)
    print*,shape(Oi)
    val = dot_product(gs_vector(:,1),Oi%dot(gs_vector(:,1)))
    print*,pos,val    


    print*,""
    print*,"pos=N"
    print*,""
    print*,"Method <psi|I_L.x.O|psi>"
    pos=N
    print*,shape(I_L)
    Oi = Build_Op_Dmrg(Op,pos)
    do it=1,R
       U  = right%omatrices%op(index=it)
       Oi = Id(dot%dim).x.(matmul(matmul(U%t(),Oi),U))
    enddo
    Oi = sp_kron(I_L,Oi,sb_states)
    print*,shape(Oi)
    val = dot_product(gs_vector(:,1),OI%dot(gs_vector(:,1)))
    print*,pos,val    



    print*,""
    print*,"- - - - - - - - - - - - - - - - -"
    print*," METHOD 3: direct O_L, O_R"
    print*,"- - - - - - - - - - - - - - - - -"
    print*,""


    Nsb  = size(sb_sector)
    allocate(Dls(Nsb),Drs(Nsb),Offset(Nsb),Oleft(Nsb),Oright(Nsb))
    allocate(AI(Nsb),BI(Nsb))
    Offset=0
    do isb=1,Nsb
       qn   = sb_sector%qn(index=isb)
       Dls(isb)= sector_qn_dim(left%sectors(1),qn)
       Drs(isb)= sector_qn_dim(right%sectors(1),current_target_qn - qn)
       if(isb>1)Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
       AI(isb)%states = sb2block_states(qn,'left')
       BI(isb)%states = sb2block_states(qn,'right')
    enddo

    pos=1
    Oi = Build_Op_Dmrg(Op,pos)
    do it=1,L
       U  = left%omatrices%op(index=it)
       Oi = (matmul(matmul(U%t(),Oi),U)).x.Id(dot%dim)
    enddo
    do isb=1,Nsb
       !> get: Oi*^L 
       Oleft(isb) = sp_filter(Oi,AI(isb)%states)
    enddo
    val = dot_product(gs_vector(:,1), OdotV_direct(Oleft,gs_vector(:,1),'left'))
    print*,pos,val


    pos=N
    Oi = Build_Op_Dmrg(Op,pos)
    do it=1,R
       U  = right%omatrices%op(index=it)
       Oi = Id(dot%dim).x.(matmul(matmul(U%t(),Oi),U))
    enddo
    do isb=1,Nsb
       !> get: Oi*^R
       Oright(isb) = sp_filter(Oi,BI(isb)%states)
    enddo
    val = dot_product(gs_vector(:,1), OdotV_direct(Oright,gs_vector(:,1),'right'))
    print*,pos,val


    if(allocated(Dls))deallocate(Dls)
    if(allocated(Drs))deallocate(Drs)
    if(allocated(Offset))deallocate(Offset)
    if(allocated(Oleft))deallocate(Oleft)
    if(allocated(Oright))deallocate(Oright)
    if(allocated(Ai))deallocate(Ai)
    if(allocated(Bi))deallocate(Bi)


  contains



    function OdotV_direct(Op,v,direction) result(Ov)
      integer                          :: Nsb,Nloc
      type(sparse_matrix),dimension(:) :: Op
      real(8),dimension(:)             :: v
      character(len=*)                 :: direction
      real(8),dimension(size(v))       :: Ov
      real(8)                          :: val
      integer                          :: i,j,k,n
      integer                          :: ir,il,jr,jl,it
      integer                          :: ia,ib,ic,ja,jb,jc,jcol
      real(8)                          :: aval,bval
      !
      Ov=zero
      !> loop over all the SB sectors:
      select case(to_lower(direction))
      case("l","left","sys","s")
         do k=1,size(sb_sector)
            !> apply the H^L x 1^r: need to T v and Ov
            do ir=1,Drs(k)
               do il=1,Dls(k)
                  i = ir + (il-1)*Drs(k) + offset(k)
                  do jcol=1,Op(k)%row(il)%Size
                     val = Op(k)%row(il)%vals(jcol)
                     jl  = Op(k)%row(il)%cols(jcol)
                     j   = ir + (jl-1)*Drs(k) + offset(k)
                     Ov(i) = Ov(i) + val*v(j)
                  end do
               enddo
            enddo
         enddo
      case("r","right","env","e")
         do k=1,size(sb_sector)
            !> apply the 1^L x H^r
            do il=1,Drs(k)
               do ir=1,Dls(k)
                  i = il + (ir-1)*Drs(k) + offset(k)           
                  do jcol=1,Op(k)%row(il)%Size
                     val = Op(k)%row(il)%vals(jcol)
                     jl  = Op(k)%row(il)%cols(jcol)
                     j   = jl + (ir-1)*Drs(k) + offset(k)
                     Ov(i) = Ov(i) + val*v(j)
                  end do
               enddo
            enddo
         enddo
      end select
      !
    end function OdotV_direct

  end subroutine Measure_Op_dmrg












  !##################################################################
  !              BUILD LOCAL OPERATOR 
  !Purpose: return the O(i) at a site I of the chain given an 
  !         operator O in the local dot basis:
  !##################################################################
  function Build_Op_dmrg(Op,pos,set_basis) result(Oi)
    type(sparse_matrix),intent(in) :: Op
    integer                        :: pos
    type(sparse_matrix)            :: Oi
    logical,optional               :: set_basis   
    !
    character(len=1)               :: label
    type(sparse_matrix)            :: U
    integer                        :: L,R,N
    integer                        :: i
    logical                        :: set_basis_
    !
    set_basis_ = .false. ;if(present(set_basis))set_basis_=set_basis
    !
    !The lenght of the last block contributing to the SB construction-> \psi
    L = left%length-1
    R = right%length-1
    N = L+R
    !
    !Check:
    if(pos<1.OR.pos>N)stop "Build_op_dmrg error: Pos not in [1,Ldmrg]"
    !
    !Get label of the block holding the site at position pos:
    label='l'; if(pos>L)label='r'
    !
    !Get index in the block from the position pos in the chain:
    i=pos    ; if(pos>L)i=N+1-pos
    !
    !Build Operator on the chain at position pos:
    select case(label)
    case('l')
       if(i==1)then
          Oi = Op
       else
          U  = left%omatrices%op(index=i-1)
          Oi = matmul(U%dgr(),U).x.Op
          if(set_basis_)then
             U  = left%omatrices%op(index=i)
             Oi = matmul(U%dgr(),matmul(Oi,U))
          endif
       endif
    case('r')
       if(i==1)then
          Oi = Op
       else
          U  = right%omatrices%op(index=i-1)
          Oi = Op.x.matmul(U%dgr(),U)
          if(set_basis_)then
             U  = right%omatrices%op(index=i)
             Oi = matmul(U%dgr(),matmul(Oi,U))
          endif
       endif
    end select
    call U%free()
  end function Build_Op_dmrg





  !##################################################################
  !                   ADVANCE OPERATOR 
  !Purpose: advance the operator O(i) Nstep from site I 
  !##################################################################
  function Advance_Op_dmrg(Op,pos,nstep) result(Oi)
    type(sparse_matrix),intent(in)   :: Op
    integer                          :: pos
    integer,optional                 :: nstep
    type(sparse_matrix)              :: Oi,U
    character(len=1)                 :: label
    integer                          :: L,R,N
    integer                          :: i,istart,iend,it
    !
    !The lenght of the last block contributing to the SB construction-> \psi
    L = left%length-1
    R = right%length-1
    N = L+R
    !
    !Check:
    if(pos<1.OR.pos>N)stop "Advance_op_dmrg error: Pos not in [1,Ldmrg]"
    !
    !Get label of the block holding the site at position pos:
    label='l'; if(pos>L)label='r'
    !
    !Get index in the block from the position pos in the chain:
    i=pos    ; if(pos>L)i=N+1-pos
    !
    istart  = i
    select case(label)
    case ("l")
       istart = i ; iend   = L ; if(present(nstep))iend=istart+nstep
       if(iend>L)stop "Advance_Op_DMRG ERROR: iend > L"
    case ("r") 
       istart = i ; iend   = R ; if(present(nstep))iend=istart+nstep
       if(iend>R)stop "Advance_Op_DMRG ERROR: iend > R"
    end select
    !
    Oi = Op
    !Evolve to SB basis
    select case(label)
    case ("l")
       do it=istart,iend
          U  = left%omatrices%op(index=it)
          Oi = (matmul(matmul(U%dgr(),Oi),U)).x.Id(dot%dim)
       enddo
    case ("r")
       do it=istart,iend
          U  = right%omatrices%op(index=it)
          Oi = Id(dot%dim).x.(matmul(matmul(U%dgr(),Oi),U))
       enddo
    end select
    call U%free()
  end function Advance_Op_dmrg



  ! !##################################################################
  ! !                   ADVANCE CORRELATION FUNCTION 
  ! !Purpose: advance the correlation O(i) Nstep from site I 
  ! !##################################################################
  function Advance_Corr_dmrg(Op,pos,nstep) result(Oi)
    type(sparse_matrix),intent(in)   :: Op
    integer                          :: pos
    integer,optional                 :: nstep
    type(sparse_matrix)              :: Oi,U
    character(len=1)                 :: label
    integer                          :: L,R,N
    integer                          :: i,istart,iend,it
    !
    !The lenght of the last block contributing to the SB construction-> \psi
    L = left%length-1             !
    R = right%length-1            !
    N = L+R                       !== Ldmrg
    !
    !Check:
    if(pos<1.OR.pos>N)stop "Advance_op_dmrg error: Pos not in [1,Ldmrg]"
    !
    !Get label of the block holding the site at position pos:
    label='l'; if(pos>L)label='r'
    !
    !Get index in the block from the position pos in the chain:
    i=pos    ; if(pos>L)i=N+1-pos
    !
    istart  = i
    select case(label)
    case ("l")
       istart = i ; iend   = L ; if(present(nstep))iend=istart+nstep
       if(iend>L)stop "Advance_Op_DMRG ERROR: iend > L"
    case ("r") 
       istart = i ; iend   = R ; if(present(nstep))iend=istart+nstep
       if(iend>R)stop "Advance_Op_DMRG ERROR: iend > R"
    end select
    !
    !
    Oi = Op
    select case(label)
    case ("l")
       do it=istart+1,iend
          U  = left%omatrices%op(index=it)
          Oi = (matmul(matmul(U%dgr(),Oi),U)).x.Id(dot%dim)
       enddo
    case ("r")
       do it=istart+1,iend
          U  = right%omatrices%op(index=it)
          Oi = Id(dot%dim).x.(matmul(matmul(U%dgr(),Oi),U))
       enddo
    end select
    call U%free()
  end function Advance_Corr_dmrg







  !##################################################################
  !                   AVERAGE OPERATOR 
  !Purpose: take the average of an operator O on the last step basis 
  !##################################################################
  function Average_Op_dmrg(Oi,pos) result(Oval)
    type(sparse_matrix),intent(in)   :: Oi
    integer                          :: pos
    character(len=1)                 :: label
    real(8)                          :: Oval
    type(sparse_matrix)              :: Psi
    integer                          :: L,R,N
    !
    !The lenght of the last block contributing to the SB construction-> \psi
    L = left%length-1
    R = right%length-1
    N = L+R
    !
    !Check:
    if(pos<1.OR.pos>N)stop "Average_op_dmrg error: Pos not in [1,Ldmrg]"
    !
    !Get label of the block holding the site at position pos:
    label='l'; if(pos>L)label='r'
    !
    !Measure using PSI matrix:
    select case(label)
    case ("l")
       if(any(shape(psi_left)/=shape(Oi)))stop "average_op_dmrg ERROR: shape(psi_left) != shape(Oi)"
       Psi  = as_sparse(psi_left)
       Oval = trace(as_matrix(matmul(matmul(Psi%dgr(),Oi),Psi)))
    case ("r")
       if(any(shape(psi_right)/=shape(Oi)))stop "average_op_dmrg ERROR: shape(psi_right) != shape(Oi)"
       Psi  = as_sparse(psi_right)
       Oval = trace(as_matrix(matmul(matmul(Psi%dgr(),Oi),Psi)))
    end select
    call Psi%free()
  end function Average_Op_dmrg



  !##################################################################
  !                   AVERAGE OPERATOR 
  !Purpose: take the average of an operator O on the last step basis 
  !##################################################################
  function NuAverage_Op_dmrg(Oi,pos) result(Oval)
    type(sparse_matrix),intent(in)   :: Oi
    integer                          :: pos
    character(len=1)                 :: label
    real(8)                          :: Oval
    type(sparse_matrix)              :: Psi
    integer                          :: L,R,N
    !
    !The lenght of the last block contributing to the SB construction-> \psi
    L = left%length-1
    R = right%length-1
    N = L+R
    !
    !Check:
    if(pos<1.OR.pos>N)stop "Average_op_dmrg error: Pos not in [1,Ldmrg]"
    !
    !Get label of the block holding the site at position pos:
    label='l'; if(pos>L)label='r'
    !
    allocate(Olist(Nsb))
    !
    !Measure using PSI matrix:
    select case(label)
    case ("l")
       !> get: Oi*^L 
       do isb=1,Nsb
          Olist(isb) = sp_filter(Oi,AI(isb)%states)
       enddo
    case ("r")
       !> get: Oi*^R
       do isb=1,Nsb
          Olist(isb) = sp_filter(Oi,BI(isb)%states)
       enddo
    end select
    !
    Oval = dot_product(gs_vector(:,1), OdotV_direct(Olist,gs_vector(:,1),label))
    !
    do isb=1,Nsb
       call Olist(isb)%free()
    enddo
    deallocate(Olist)
    !
  end function NuAverage_Op_dmrg


  !#################################
  !#################################




  function OdotV_direct(Op,v,direction) result(Ov)
    integer                          :: Nsb,Nloc
    type(sparse_matrix),dimension(:) :: Op
    real(8),dimension(:)             :: v
    character(len=*)                 :: direction
    real(8),dimension(size(v))       :: Ov
    real(8)                          :: val
    integer                          :: i,j,k,n
    integer                          :: ir,il,jr,jl,it
    integer                          :: ia,ib,ic,ja,jb,jc,jcol
    real(8)                          :: aval,bval
    !
    Ov=zero
    !> loop over all the SB sectors:
    select case(to_lower(direction))
    case("l","left","sys","s")
       do k=1,size(sb_sector)
          !> apply the H^L x 1^r: need to T v and Ov
          do ir=1,Drs(k)
             do il=1,Dls(k)
                i = ir + (il-1)*Drs(k) + offset(k)
                do jcol=1,Op(k)%row(il)%Size
                   val = Op(k)%row(il)%vals(jcol)
                   jl  = Op(k)%row(il)%cols(jcol)
                   j   = ir + (jl-1)*Drs(k) + offset(k)
                   Ov(i) = Ov(i) + val*v(j)
                end do
             enddo
          enddo
       enddo
    case("r","right","env","e")
       do k=1,size(sb_sector)
          !> apply the 1^L x H^r
          do il=1,Drs(k)
             do ir=1,Dls(k)
                i = il + (ir-1)*Drs(k) + offset(k)           
                do jcol=1,Op(k)%row(il)%Size
                   val = Op(k)%row(il)%vals(jcol)
                   jl  = Op(k)%row(il)%cols(jcol)
                   j   = jl + (ir-1)*Drs(k) + offset(k)
                   Ov(i) = Ov(i) + val*v(j)
                end do
             enddo
          enddo
       enddo
    end select
    !
  end function OdotV_direct



  !#################################
  !#################################



  subroutine write_user(file,vals,x)
    character(len=*) :: file
    integer,optional :: x
    integer          :: x_
    real(8)          :: vals(:)
    integer          :: Eunit
    x_ = left%length;if(present(x))x_=x
    Eunit     = fopen(str(file)//"_"//str(suffix),append=.true.)
    write(Eunit,*)x_,vals
    close(Eunit)
  end subroutine write_user



  subroutine dmrg_write(file,vals,x)
    character(len=*) :: file
    integer,optional :: x
    integer          :: x_
    real(8)          :: vals(:)
    integer          :: Eunit
    x_ = left%length;if(present(x))x_=x
    Eunit     = fopen(str(file)//"_"//str(suffix),append=.true.)
    write(Eunit,*)x_,vals
    close(Eunit)
  end subroutine dmrg_write








END MODULE DMRG_MEASURE









  ! subroutine Measure_Corr_dmrg(file)
  !   character(len=*)                          :: file
  !   character(len=1)                          :: label
  !   real(8)                                   :: val
  !   integer                                   :: it,i,L,R,N,j,pos
  !   type(sparse_matrix)                       :: Oi,U,Psi,PsiL,PsiR,Sz,Sp,Sz0,Sp0,SzL,SpL,SzR,SpR,U1,U2
  !   !
  !   suffix=label_DMRG('u')
  !   !
  !   L = left%length-1           !the len of the last block used to create the SB->\psi
  !   R = right%length-1
  !   N = L+R
  !   !
  !   Sz0  = dot%operators%op("S_z")
  !   Sp0  = dot%operators%op("S_p")
  !   !
  !   !
  !   !Get index in the block from the position pos in the chain:
  !   do pos=1,N
  !      if(pos==L.OR.pos==L+1)cycle
  !      !
  !      !Get label of the block holding the site at position pos:
  !      label='l'; if(pos>L)label='r'
  !      !
  !      Sz = Build_Op_dmrg(Sz0,pos,set_basis=.true.)
  !      Sp = Build_Op_dmrg(Sp0,pos,set_basis=.true.)
  !      !
  !      select case(label)
  !      case("l")
  !         Oi = Jp*(Sz.x.Sz0) + Jx/2d0*(Sp.x.Sp0%dgr())+ Jx/2d0*(Sp%dgr().x.Sp0)
  !      case("r")
  !         Oi = Jp*(Sz0.x.Sz) + Jx/2d0*(Sp0.x.Sp%dgr())+ Jx/2d0*(Sp0%dgr().x.Sp)
  !      end select
  !      !
  !      Oi = Advance_Corr_DMRG(Oi,pos)
  !      !
  !      val= Average_Op_dmrg(Oi,pos)
  !      call write_user(trim(file),[val],x=pos)
  !   enddo
  ! end subroutine Measure_Corr_dmrg


!#################################
!#################################




! label='l'
! do i=1,L
!    print*,i,label
!    if(i==1)then
!       Oi = Op
!    else                     !i>=2
!       U  = left%omatrices%op(index=i-1)
!       Oi = matmul(U%dgr(),U).x.Op
!    endif
!    !Evolve to SB basis
!    do it=i,L
!       U  = left%omatrices%op(index=it)
!       Oi = (matmul(matmul(U%dgr(),Oi),U)).x.Id(dot%dim)
!    enddo
!    !Measure:
!    Psi = as_sparse(psi_left)
!    val = trace(as_matrix(matmul(matmul(Psi%t(),Oi),Psi)))
!    write(100,*)i,val
! enddo
! !
! label='r'
! do j=1,R
!    i=R+1-j
!    print*,i,label
!    if(i==1)then
!       Oi = Op
!    else                     !i>=2
!       U  = right%omatrices%op(index=i-1)
!       Oi = Op.x.matmul(U%dgr(),U)
!    endif
!    !Evolve to SB basis
!    do it=i,R
!       U  = right%omatrices%op(index=it)
!       Oi = Id(dot%dim).x.(matmul(matmul(U%dgr(),Oi),U))
!    enddo
!    Psi = as_sparse(psi_right)
!    val = trace(as_matrix(matmul(matmul(Psi%dgr(),Oi),Psi)))
!    ! print*,N+1-i,val
!    write(100,*)N+1-i,val
! enddo








! function Build_Op_dmrg(Op,i) result(Oi)
!   integer                        :: i
!   type(sparse_matrix),intent(in) :: Op
!   type(sparse_matrix)            :: Oi
!   !
!   character(len=1)               :: label
!   type(sparse_matrix)            :: U
!   integer                        :: L,R
!   integer                        :: pos
!   !
!   L = left%length
!   R = right%length
!   if(i<1.OR.i>L+R)stop "construct_op_dmrg error: I not in [1,Lchain]"
!   !
!   label='l';if(i>L)label='r'
!   !
!   pos=i;if(i>L)pos=R+1-(i-L)
!   !
!   !Step 1: build the operator at the site pos
!   Oi  = Op
!   if(pos>1)then
!      select case(to_lower(str(label)))
!      case ("l")
!         U   = left%omatrices%op(index=pos-1) !get Id of the block_{I-1} 
!         Oi  = matmul(U%dgr(),U).x.Op          !build right-most Op of the block_I as Id(I-1)xOp_I
!         U   = left%omatrices%op(index=pos)   !retrieve O_I 
!         Oi  = matmul(U%dgr(),matmul(Oi,U))    !rotate+truncate Oi at the basis of Block_I 
!      case ("r")
!         U   = right%omatrices%op(index=pos-1) !get Id of the block_{I-1} 
!         Oi  = Op.x.matmul(U%dgr(),U)          !build right-most Op of the block_I as Id(I-1)xOp_I
!         U   = right%omatrices%op(index=pos)   !retrieve O_I 
!         Oi  = matmul(U%dgr(),matmul(Oi,U))    !rotate+truncate Oi at the basis of Block_I 
!      end select
!      call U%free()
!   end if
! end function Build_Op_dmrg

!WORKING VERSION:
! !Build Operator:
! select case(label)
! case('l')
!    if(i==1)then
!       Oi = Op
!    else
!       U  = left%omatrices%op(index=i-1)
!       Oi = matmul(U%dgr(),U).x.Op
!    endif
! case('r')
!    if(i==1)then
!       Oi = Op
!    else
!       U  = right%omatrices%op(index=i-1)
!       Oi = Op.x.matmul(U%dgr(),U)
!    endif
! end select








! function Advance_Op_dmrg(Op,i,nstep) result(Oi)
!   type(sparse_matrix),intent(in)   :: Op
!   integer                          :: i
!   integer,optional                 :: nstep
!   type(sparse_matrix)              :: Oi,U
!   character(len=1)                 :: label
!   integer                          :: L,R
!   integer                          :: istart,iend,it
!   !
!   L = left%length-1
!   R = right%length-1
!   if(i<1.OR.i>L+R)stop "Advance_Op_DMRG error: I not in [1,Lchain]"
!   !
!   label ='l';if(i>L) label ='r'
!   select case(label)
!   case ("l")
!      istart = i         ; iend   = L ; if(present(nstep))iend=istart+nstep
!      if(iend>L)stop "Advance_Op_DMRG ERROR: iend > L"
!   case ("r") 
!      istart = R+1-(i-L) ; iend   = R ; if(present(nstep))iend=istart+nstep
!      if(iend>R)stop "Advance_Op_DMRG ERROR: iend > R"
!   end select
!   !
!   !Bring the operator to the L/R basis of the target states
!   !Here we try to reproduce the same strategy of the dmrg_step: increase and rotate 
!   Oi = Op
!   select case(label)
!   case ("l")
!      Oi = Oi.x.Id(dot%dim)
!      do it=istart,iend
!         U  = left%omatrices%op(index=it)
!         Oi = matmul(U%dgr(),matmul(Oi,U))
!         Oi = Oi.x.Id(dot%dim)
!      enddo
!   case ("r")
!      Oi = Id(dot%dim).x.Oi
!      do it=istart,iend
!         U  = right%omatrices%op(index=it)
!         Oi = matmul(U%dgr(),matmul(Oi,U))
!         Oi = Id(dot%dim).x.Oi
!      enddo
!   end select
!   call U%free()
! end function Advance_Op_dmrg

!WORKING VERSION:
! select case(label)
! case('l')
!    !Evolve to SB basis
!    do it=i,L
!       U  = left%omatrices%op(index=it)
!       Oi = (matmul(matmul(U%dgr(),Oi),U)).x.Id(dot%dim)
!    enddo
! case('r')
!    do it=i,R
!       U  = right%omatrices%op(index=it)
!       Oi = Id(dot%dim).x.(matmul(matmul(U%dgr(),Oi),U))
!    enddo
! end select






! function Average_Op_dmrg(Op,label) result(avOp)
!   type(sparse_matrix),intent(in)   :: Op
!   character(len=1)                 :: label
!   real(8)                          :: avOp
!   type(sparse_matrix)              :: Psi
!   integer                          :: L,R,i,pos,dims(2)
!   real(8),dimension(:),allocatable :: psi_vec
!   integer,dimension(:),allocatable :: psi_map
!   !
!   !Measure using PSI matrix:
!   select case(label)
!   case ("l")
!      if(any(shape(psi_left)/=shape(Op)))stop "average_op_dmrg ERROR: shape(psi) != shape(Op)"
!      Psi  = as_sparse(psi_left)
!      AvOp = trace(as_matrix(matmul(Psi%dgr(),matmul(Op,Psi))))
!   case ("r")
!      if(any(shape(psi_right)/=shape(Op)))stop "average_op_dmrg ERROR: shape(psi) != shape(Op)"
!      Psi  = as_sparse(psi_right)
!      AvOp = trace(as_matrix(matmul(Psi%dgr(),matmul(Op,Psi))))
!   end select
!   call Psi%free()
! end function Average_Op_dmrg

!WORKING VERSION:
! select case(label)
! case('l')
!    Psi = as_sparse(psi_left)
!    val = trace(as_matrix(matmul(matmul(Psi%t(),Oi),Psi)))
! case('r')
!    Psi = as_sparse(psi_right)
!    val = trace(as_matrix(matmul(matmul(Psi%dgr(),Oi),Psi)))
! end select







  ! subroutine Measure_Op_dmrg(Op,file,ref,avOp)
  !   type(sparse_matrix),intent(in)            :: Op
  !   character(len=*)                          :: file
  !   character(len=1)                          :: label
  !   real(8) :: ref
  !   real(8),dimension(:),allocatable,optional :: avOp
  !   real(8)                                   :: val
  !   integer                                   :: it,i,L,R,N,j,pos,dims(2)
  !   type(sparse_matrix)                       :: Oi,U,Psi,Ok,I_R,I_L
  !   !

  !   suffix=label_DMRG('u')
  !   !
  !   L = left%length-1           !the len of the last block used to create the SB->\psi
  !   R = right%length-1
  !   N = L+R
  !   !
  !   if(present(avOp))then
  !      if(allocated(avOp))deallocate(avOp)
  !      allocate(avOp(N))
  !   endif
  !   !
  !   ! call start_timer()
  !   ! do pos=1,N
  !   !    Oi = Build_Op_dmrg(Op,pos)
  !   !    Oi = Advance_Op_dmrg(Oi,pos)
  !   !    val= Average_Op_dmrg(Oi,pos)
  !   !    if(present(avOp))avOp(pos)=val
  !   !    call write_user(trim(file),[val],x=pos)
  !   !    call Oi%free()
  !   !    call progress(pos,N)
  !   ! enddo
  !   ! call stop_timer("Done "//str(file))




  !   ! U =  right%omatrices%op(index=R)
  !   ! dims=shape(U)
  !   ! ! I_R = Id(dot%dim).x.Id(dims(2))!(matmul(U%t(),U))
  !   ! I_R = Id(dot%dim*dims(2))
  !   ! U =  left%omatrices%op(index=R)
  !   ! dims=shape(U)
  !   ! ! I_L = Id(dims(2)).x.Id(dot%dim)!(matmul(U%t(),U)).x.Id(dot%dim)
  !   ! I_L = Id(dims(2)*dot%dim)!(matmul(U%t(),U)).x.Id(dot%dim)


  !   ! print*,"size(GSpsi)",size(gs_vector(:,1))



  !   ! print*,""
  !   ! print*,"- - - - - - - - - - - - - - - - -"
  !   ! print*," METHOD 1: O.x.I - I.x.O full"
  !   ! print*,"- - - - - - - - - - - - - - - - -"
  !   ! print*,""


  !   ! print*,""
  !   ! print*,"pos=1"
  !   ! print*,""
  !   ! print*,"Method <psi|O.x.I_R|psi>"
  !   ! print*,shape(I_R)
  !   ! pos=1
  !   ! Oi = Build_Op_Dmrg(Op,pos)
  !   ! do it=1,L
  !   !    U  = left%omatrices%op(index=it)
  !   !    Oi = (matmul(matmul(U%t(),Oi),U)).x.Id(dot%dim)
  !   ! enddo
  !   ! Oi = sp_kron(Oi,I_R,sb_states)
  !   ! print*,shape(Oi)
  !   ! val = dot_product(gs_vector(:,1),Oi%dot(gs_vector(:,1)))
  !   ! print*,pos,val    


  !   ! print*,""
  !   ! print*,"pos=N"
  !   ! print*,""
  !   ! print*,"Method <psi|I_L.x.O|psi>"
  !   ! pos=N
  !   ! print*,shape(I_L)
  !   ! Oi = Build_Op_Dmrg(Op,pos)
  !   ! do it=1,R
  !   !    U  = right%omatrices%op(index=it)
  !   !    Oi = Id(dot%dim).x.(matmul(matmul(U%t(),Oi),U))
  !   ! enddo
  !   ! Oi = sp_kron(I_L,Oi,sb_states)
  !   ! print*,shape(Oi)
  !   ! val = dot_product(gs_vector(:,1),OI%dot(gs_vector(:,1)))
  !   ! print*,pos,val    



  !   ! print*,""
  !   ! print*,"- - - - - - - - - - - - - - - - -"
  !   ! print*," METHOD 2: O_L.psi_L, O_R.psi_R  "
  !   ! print*,"- - - - - - - - - - - - - - - - -"
  !   ! print*,""



  !   ! print*,""
  !   ! print*,"Method <psiL|O_L|psiL>"
  !   ! pos=1
  !   ! Oi = Build_Op_Dmrg(Op,pos)
  !   ! do it=1,L
  !   !    U  = left%omatrices%op(index=it)
  !   !    Oi = (matmul(matmul(U%t(),Oi),U)).x.Id(dot%dim)
  !   ! enddo
  !   ! print*,"psiL_blocks:",shape(psi_left)
  !   ! Psi  = as_sparse(psi_left)
  !   ! print*,"psiL_blocks:",shape(psi_left)
  !   ! print*,shape(Psi%t()),shape(Oi),shape(Psi)
  !   ! val  = trace(as_matrix(matmul(matmul(Psi%t(),Oi),Psi)))
  !   ! print*,pos,val



  !   ! print*,""
  !   ! print*,"Method <psiR|O_R|psiR>"
  !   ! pos=N
  !   ! Oi = Build_Op_Dmrg(Op,pos)
  !   ! do it=1,R
  !   !    U  = right%omatrices%op(index=it)
  !   !    Oi = Id(dot%dim).x.(matmul(matmul(U%t(),Oi),U))
  !   ! enddo
  !   ! print*,"psiR_blocks:",shape(psi_right)
  !   ! Psi  = as_sparse(psi_right)
  !   ! print*,shape(Psi%t()),shape(Oi),shape(Psi)
  !   ! val  = trace(as_matrix(matmul(matmul(Psi%t(),Oi),Psi)))
  !   ! print*,pos,val





  !   print*,""
  !   print*,"- - - - - - - - - - - - - - - - -"
  !   print*," METHOD 3: direct O_L, O_R"
  !   print*,"- - - - - - - - - - - - - - - - -"
  !   print*,""


  !   Nsb  = size(sb_sector)
  !   allocate(Dls(Nsb),Drs(Nsb),Offset(Nsb),Oleft(Nsb),Oright(Nsb))
  !   allocate(AI(Nsb),BI(Nsb))
  !   Offset=0
  !   do isb=1,Nsb
  !      qn   = sb_sector%qn(index=isb)
  !      Dls(isb)= sector_qn_dim(left%sectors(1),qn)
  !      Drs(isb)= sector_qn_dim(right%sectors(1),current_target_qn - qn)
  !      if(isb>1)Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
  !      AI(isb)%states = sb2block_states(qn,'left')
  !      BI(isb)%states = sb2block_states(qn,'right')
  !   enddo

  !   pos=1
  !   Oi = Build_Op_Dmrg(Op,pos)
  !   do it=1,L
  !      U  = left%omatrices%op(index=it)
  !      Oi = (matmul(matmul(U%t(),Oi),U)).x.Id(dot%dim)
  !   enddo
  !   do isb=1,Nsb
  !      !> get: Oi*^L 
  !      Oleft(isb) = sp_filter(Oi,AI(isb)%states)
  !   enddo
  !   val = dot_product(gs_vector(:,1), OdotV_direct(Oleft,gs_vector(:,1),'left'))
  !   print*,pos,val


  !   pos=N
  !   Oi = Build_Op_Dmrg(Op,pos)
  !   do it=1,R
  !      U  = right%omatrices%op(index=it)
  !      Oi = Id(dot%dim).x.(matmul(matmul(U%t(),Oi),U))
  !   enddo
  !   do isb=1,Nsb
  !      !> get: Oi*^R
  !      Oright(isb) = sp_filter(Oi,BI(isb)%states)
  !   enddo
  !   val = dot_product(gs_vector(:,1), OdotV_direct(Oright,gs_vector(:,1),'right'))
  !   print*,pos,val


  !   if(allocated(Dls))deallocate(Dls)
  !   if(allocated(Drs))deallocate(Drs)
  !   if(allocated(Offset))deallocate(Offset)
  !   if(allocated(Oleft))deallocate(Oleft)
  !   if(allocated(Oright))deallocate(Oright)
  !   if(allocated(Ai))deallocate(Ai)
  !   if(allocated(Bi))deallocate(Bi)


  ! contains



  !   function OdotV_direct(Op,v,direction) result(Ov)
  !     integer                          :: Nsb,Nloc
  !     type(sparse_matrix),dimension(:) :: Op
  !     real(8),dimension(:)             :: v
  !     character(len=*)                 :: direction
  !     real(8),dimension(size(v))       :: Ov
  !     real(8)                          :: val
  !     integer                          :: i,j,k,n
  !     integer                          :: ir,il,jr,jl,it
  !     integer                          :: ia,ib,ic,ja,jb,jc,jcol
  !     real(8)                          :: aval,bval
  !     !
  !     Ov=zero
  !     !> loop over all the SB sectors:
  !     select case(to_lower(direction))
  !     case("l","left","sys","s")
  !        do k=1,size(sb_sector)
  !           !> apply the H^L x 1^r: need to T v and Ov
  !           do ir=1,Drs(k)
  !              do il=1,Dls(k)
  !                 i = ir + (il-1)*Drs(k) + offset(k)
  !                 do jcol=1,Op(k)%row(il)%Size
  !                    val = Op(k)%row(il)%vals(jcol)
  !                    jl  = Op(k)%row(il)%cols(jcol)
  !                    j   = ir + (jl-1)*Drs(k) + offset(k)
  !                    Ov(i) = Ov(i) + val*v(j)
  !                 end do
  !              enddo
  !           enddo
  !        enddo
  !     case("r","right","env","e")
  !        do k=1,size(sb_sector)
  !           !> apply the 1^L x H^r
  !           do il=1,Drs(k)
  !              do ir=1,Dls(k)
  !                 i = il + (ir-1)*Drs(k) + offset(k)           
  !                 do jcol=1,Op(k)%row(il)%Size
  !                    val = Op(k)%row(il)%vals(jcol)
  !                    jl  = Op(k)%row(il)%cols(jcol)
  !                    j   = jl + (ir-1)*Drs(k) + offset(k)
  !                    Ov(i) = Ov(i) + val*v(j)
  !                 end do
  !              enddo
  !           enddo
  !        enddo
  !     end select
  !     !
  !   end function OdotV_direct

  ! end subroutine Measure_Op_dmrg
