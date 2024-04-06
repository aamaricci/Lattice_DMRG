MODULE DMRG_MEASURE
  USE VARS_GLOBAL
  implicit none
  private


  !Measuring:
  public :: Measure_Op_DMRG
  public :: Build_Op_DMRG
  public :: Advance_Op_DMRG
  public :: Advance_Corr_DMRG
  public :: Average_Op_DMRG



contains


  !##################################################################
  !              MEASURE LOCAL OPERATOR
  !##################################################################
  subroutine Measure_Op_dmrg(Op,ivec,file,avOp)
    integer,dimension(:),intent(in)        :: ivec
    type(sparse_matrix),intent(in)         :: Op
    character(len=*)                       :: file
    character(len=1)                       :: label
    real(8),dimension(size(ivec)),optional :: avOp
    real(8)                                :: val
    integer                                :: it,i
    type(sparse_matrix)                    :: Oi
    !
    suffix=label_DMRG('u')
    do it=1,size(ivec)
       i   = ivec(it)
       label='l';if(i>left%length)label='r'
       Oi  = Build_Op_dmrg(Op,i)
       Oi  = Advance_Op_dmrg(Oi,i)
       val = Average_Op_dmrg(Oi,label)   
       call write_user(trim(file),[val],x=dble(i))
       call Oi%free()
       if(present(avOp))avOp(it)=val
    enddo
  end subroutine Measure_Op_dmrg


  !##################################################################
  !              BUILD LOCAL OPERATOR 
  !Purpose: return the O(i) at a site I of the chain given an 
  !         operator O in the local dot basis:
  !##################################################################
  function Build_Op_dmrg(Op,i) result(Oi)
    integer                        :: i
    type(sparse_matrix),intent(in) :: Op
    type(sparse_matrix)            :: Oi
    !
    character(len=1)               :: label
    type(sparse_matrix)            :: U
    integer                        :: L,R
    integer                        :: pos
    !
    L = left%length
    R = right%length
    if(i<1.OR.i>L+R)stop "construct_op_dmrg error: I not in [1,Lchain]"
    !
    label='l';if(i>L)label='r'
    !
    pos=i;if(i>L)pos=R+1-(i-L)
    !
    !Step 1: build the operator at the site pos
    Oi  = Op
    if(pos>1)then
       select case(to_lower(str(label)))
       case ("l")
          U   = left%omatrices%op(index=pos-1) !get Id of the block_{I-1} 
          Oi  = matmul(U%dgr(),U).x.Op          !build right-most Op of the block_I as Id(I-1)xOp_I
          U   = left%omatrices%op(index=pos)   !retrieve O_I 
          Oi  = matmul(U%dgr(),matmul(Oi,U))    !rotate+truncate Oi at the basis of Block_I 
       case ("r")
          U   = right%omatrices%op(index=pos-1) !get Id of the block_{I-1} 
          Oi  = Op.x.matmul(U%dgr(),U)          !build right-most Op of the block_I as Id(I-1)xOp_I
          U   = right%omatrices%op(index=pos)   !retrieve O_I 
          Oi  = matmul(U%dgr(),matmul(Oi,U))    !rotate+truncate Oi at the basis of Block_I 
       end select
       call U%free()
    end if
  end function Build_Op_dmrg



  !##################################################################
  !                   ADVANCE OPERATOR 
  !Purpose: advance the operator O(i) Nstep from site I 
  !##################################################################
  function Advance_Op_dmrg(Op,i,nstep) result(Oi)
    type(sparse_matrix),intent(in)   :: Op
    integer                          :: i
    integer,optional                 :: nstep
    type(sparse_matrix)              :: Oi,U
    character(len=1)                 :: label
    integer                          :: L,R
    integer                          :: istart,iend,it
    !
    L = left%length
    R = right%length
    if(i<1.OR.i>L+R)stop "Advance_Op_DMRG error: I not in [1,Lchain]"
    !
    label ='l';if(i>L) label ='r'
    select case(label)
    case ("l")
       istart = i         ; iend   = L ; if(present(nstep))iend=istart+nstep
       if(iend>L)stop "Advance_Op_DMRG ERROR: iend > L"
    case ("r") 
       istart = R+1-(i-L) ; iend   = R ; if(present(nstep))iend=istart+nstep
       if(iend>R)stop "Advance_Op_DMRG ERROR: iend > R"
    end select
    !
    !Bring the operator to the L/R basis of the target states
    !Here we try to reproduce the same strategy of the dmrg_step: increase and rotate 
    Oi = Op
    select case(label)
    case ("l")
       Oi = Oi.x.Id(dot%dim)
       do it=istart+1,iend
          U  = left%omatrices%op(index=it)
          Oi = matmul(U%dgr(),matmul(Oi,U))
          Oi = Oi.x.Id(dot%dim)
       enddo
    case ("r")
       Oi = Id(dot%dim).x.Oi
       do it=istart+1,iend
          U  = right%omatrices%op(index=it)
          Oi = matmul(U%dgr(),matmul(Oi,U))
          Oi = Id(dot%dim).x.Oi
       enddo
    end select
    call U%free()
  end function Advance_Op_dmrg


  !##################################################################
  !                   ADVANCE CORRELATION FUNCTION 
  !Purpose: advance the correlation O(i) Nstep from site I 
  !##################################################################
  function Advance_Corr_dmrg(Op,i,nstep) result(Oi)
    type(sparse_matrix),intent(in)   :: Op
    integer                          :: i
    integer,optional                 :: nstep
    type(sparse_matrix)              :: Oi,U
    character(len=1)                 :: label
    integer                          :: L,R
    integer                          :: istart,iend,it
    !
    L = left%length
    R = right%length
    if(i<1.OR.i>L+R)stop "Advance_Op_DMRG error: I not in [1,Lchain]"
    !
    label ='l';if(i>L) label ='r'
    select case(label)
    case ("l")
       istart = i         ; iend   = L ; if(present(nstep))iend=istart+nstep
       if(iend>L)stop "Advance_Op_DMRG ERROR: iend > L"
    case ("r") 
       istart = R+1-(i-L) ; iend   = R ; if(present(nstep))iend=istart+nstep
       if(iend>R)stop "Advance_Op_DMRG ERROR: iend > R"
    end select
    !
    !Bring the operator to the L/R basis of the target states
    !Here we try to reproduce the same strategy of the dmrg_step: increase and rotate 
    Oi = Op
    select case(label)
    case ("l")
       do it=istart+1,iend
          U  = left%omatrices%op(index=it)
          Oi = matmul(U%dgr(),matmul(Oi,U))
          Oi = Oi.x.Id(dot%dim)
       enddo
    case ("r")
       do it=istart+1,iend
          U  = right%omatrices%op(index=it)
          Oi = matmul(U%dgr(),matmul(Oi,U))
          Oi = Id(dot%dim).x.Oi
       enddo
    end select
    call U%free()
  end function Advance_Corr_dmrg



  !##################################################################
  !                   AVERAGE OPERATOR 
  !Purpose: take the average of an operator O on the last step basis 
  !##################################################################
  function Average_Op_dmrg(Op,label) result(avOp)
    type(sparse_matrix),intent(in)   :: Op
    character(len=1)                 :: label
    real(8)                          :: avOp
    type(sparse_matrix)              :: Psi
    integer                          :: L,R,i,pos,dims(2)
    real(8),dimension(:),allocatable :: psi_vec
    integer,dimension(:),allocatable :: psi_map
    !
    !Measure using PSI matrix:
    select case(label)
    case ("l")
       if(any(shape(psi_left)/=shape(Op)))stop "average_op_dmrg ERROR: shape(psi) != shape(Op)"
       Psi  = as_sparse(psi_left)
       AvOp = trace(as_matrix(matmul(Psi%dgr(),matmul(Op,Psi))))
    case ("r")
       if(any(shape(psi_right)/=shape(Op)))stop "average_op_dmrg ERROR: shape(psi) != shape(Op)"
       Psi  = as_sparse(psi_right)
       AvOp = trace(as_matrix(matmul(Psi%dgr(),matmul(Op,Psi))))
    end select
    call Psi%free()
  end function Average_Op_dmrg






  subroutine write_user(file,vals,x)
    character(len=*) :: file
    real(8),optional :: x
    real(8)          :: x_
    real(8)          :: vals(:)
    integer          :: Eunit
    x_ = left%length;if(present(x))x_=x
    Eunit     = fopen(str(file)//"_"//str(suffix),append=.true.)
    write(Eunit,*)x_,vals
    close(Eunit)
  end subroutine write_user











END MODULE DMRG_MEASURE





