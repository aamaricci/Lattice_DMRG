MODULE DMRG_MEASURE
  USE VARS_GLOBAL
  implicit none
  private


  !Measuring:
  public :: Measure_Op_DMRG
  public :: Measure_Corr_dmrg
  public :: Build_Op_DMRG
  public :: Advance_Op_DMRG
  public :: Advance_Corr_DMRG
  public :: Average_Op_DMRG



contains


  !##################################################################
  !              MEASURE LOCAL OPERATOR
  !##################################################################
  subroutine Measure_Op_dmrg(Op,file,avOp)
    type(sparse_matrix),intent(in)            :: Op
    character(len=*)                          :: file
    character(len=1)                          :: label
    real(8),dimension(:),allocatable,optional :: avOp
    real(8)                                   :: val
    integer                                   :: it,i,L,R,N,j,pos
    type(sparse_matrix)                       :: Oi,U,Psi,Sz,Sp,Sz0,Sp0,SiSj
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
  end subroutine Measure_Op_dmrg






  subroutine Measure_Corr_dmrg(file)
    character(len=*)                          :: file
    character(len=1)                          :: label
    real(8)                                   :: val
    integer                                   :: it,i,L,R,N,j,pos
    type(sparse_matrix)                       :: Oi,U,Psi,PsiL,PsiR,Sz,Sp,Sz0,Sp0,SzL,SpL,SzR,SpR,U1,U2
    !
    suffix=label_DMRG('u')
    !
    L = left%length-1           !the len of the last block used to create the SB->\psi
    R = right%length-1
    N = L+R
    !
    Sz0  = dot%operators%op("S_z")
    Sp0  = dot%operators%op("S_p")
    !
    !
    !Get index in the block from the position pos in the chain:
    do pos=1,L-1
       Sz = Build_Op_dmrg(Sz0,pos,set_basis=.true.)
       Sp = Build_Op_dmrg(Sp0,pos,set_basis=.true.)
       !
       Oi = Jp*(Sz.x.Sz0) + Jx/2d0*(Sp.x.Sp0%dgr())+ Jx/2d0*(Sp%dgr().x.Sp0)
       !
       Oi = Advance_Corr_DMRG(Oi,pos)
       !
       val= Average_Op_dmrg(Oi,pos)
       write(100,*)pos,val
    enddo



    do pos=L+2,N
       i = N+1-pos
       Sz = Build_Op_dmrg(Sz0,pos,set_basis=.true.)
       Sp = Build_Op_dmrg(Sp0,pos,set_basis=.true.)
       !
       Oi   = Jp*(Sz0.x.Sz) + Jx/2d0*(Sp0.x.Sp%dgr())+ Jx/2d0*(Sp0%dgr().x.Sp)
       !
       Oi = Advance_Corr_DMRG(Oi,pos)
       !
       val= Average_Op_dmrg(Oi,pos)
       write(200,*)pos,val
    enddo





    !Get index in the block from the position pos in the chain:
    do pos=1,N
       if(pos==L.OR.pos==L+1)cycle
       !
       !Get label of the block holding the site at position pos:
       label='l'; if(pos>L)label='r'

       !
       Sz = Build_Op_dmrg(Sz0,pos,set_basis=.true.)
       Sp = Build_Op_dmrg(Sp0,pos,set_basis=.true.)
       !
       select case(label)
       case("l")
          Oi = Jp*(Sz.x.Sz0) + Jx/2d0*(Sp.x.Sp0%dgr())+ Jx/2d0*(Sp%dgr().x.Sp0)
       case("r")
          Oi = Jp*(Sz0.x.Sz) + Jx/2d0*(Sp0.x.Sp%dgr())+ Jx/2d0*(Sp0%dgr().x.Sp)
       end select
       !
       Oi = Advance_Corr_DMRG(Oi,pos)
       !
       val= Average_Op_dmrg(Oi,pos)

       print*,pos,val
       write(300,*)pos,val
    enddo

  end subroutine Measure_Corr_dmrg






  

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











END MODULE DMRG_MEASURE










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
