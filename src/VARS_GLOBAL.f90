MODULE VARS_GLOBAL
  USE SCIFOR
  USE INPUT_VARS
  USE AUX_FUNCS
  USE MATRIX_GRAPH
  USE MATRIX_SPARSE, id=>sp_eye
  USE MATRIX_BLOCKS
  USE TUPLE_BASIS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
  implicit none


  ! abstract interface
  !    function UserHconnect(left,right) 
  !      USE BLOCKS
  !      implicit none
  !      type(block)                          :: left
  !      type(block)                          :: right
  !      real(8),dimension(:,:),allocatable :: UserHconnect ![Nso,Nso]
  !    end function UserHconnect
  ! end interface
  ! procedure(UserHconnect),pointer,public  :: Hmodel=>null()
  ! !

#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable   :: HopH
#else
  real(8),dimension(:,:),allocatable      :: HopH
#endif




  abstract interface
     subroutine sparse_HxV(Nloc,v,Hv)
       integer                          :: Nloc
#ifdef _CMPLX
       complex(8),dimension(Nloc)       :: v
       complex(8),dimension(Nloc)       :: Hv
#else
       real(8),dimension(Nloc)          :: v
       real(8),dimension(Nloc)          :: Hv
#endif
     end subroutine sparse_HxV
  end interface
  procedure(sparse_HxV),pointer,public  :: spHtimesV_p=>null()



  type(sparse_matrix)                   :: spHsb
  !
  real(8)                               :: truncation_error_left,truncation_error_right
  character(len=:),allocatable          :: suffix
  real(8),dimension(:),allocatable      :: target_Qn,current_target_QN
  type(block)                           :: init_left,init_right
  logical                               :: init_called=.false.
#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable :: gs_vector
#else
  real(8),dimension(:,:),allocatable    :: gs_vector
#endif
  real(8),dimension(:),allocatable      :: gs_energy
  real(8),dimension(:),allocatable      :: rho_left_evals
  real(8),dimension(:),allocatable      :: rho_right_evals
  type(blocks_matrix)                   :: psi_left
  type(blocks_matrix)                   :: psi_right
  type(blocks_matrix)                   :: rho_left
  type(blocks_matrix)                   :: rho_right
  !GLOBAL LEFT & RIGHT & DOT 
  type(block)                           :: left,right
  type(site),dimension(:),allocatable   :: dot
  !
  !SUPERBLOCK SHARED THINGS
  integer,dimension(:),allocatable      :: sb_states
  type(sectors_list)                    :: sb_sector
  !
  integer                               :: Mstates
  real(8)                               :: Estates


contains



  function sector_qn_dim(self,qn) result(dim)
    type(sectors_list)   :: self
    real(8),dimension(:) :: qn
    integer              :: dim
    dim = 0
    if(.not.self%has_qn(qn))return
    dim =  size(self%map(qn=qn))
  end function sector_qn_dim



  subroutine dmrg_graphic(label)
    integer                  :: label
    integer                  :: i,N,Mleft,Mright,LMleft,LMright,index,Ltot
    character(:),allocatable :: Ldot,Rdot
    real(8)                  :: eps=1d-6
    integer                  :: M=50
    !
    Ltot = Ldmrg
    Ldot = bold_green('=')
    Rdot = bold_red('-')
    ! if(Ltot>M)then
    !    Ldot = bg_green('=')
    !    Rdot = bg_red('-')
    ! endif
    !
    N = int(Ltot/(M+eps))+1
    !
    ! call execute_command_line("clear")
    do i=1,3
       write(LOGfile,*)""
    enddo
    select case(label)
    case default; stop "dmrg_graphic error: label != 1(L),2(R)"
    case(0)
       Mleft  = int(left%length/(N+eps))+1
       Mright = int(right%length/(N+eps))+1
       LMleft = Ltot/N-Mleft
       LMright= Ltot/N-Mright
       index=nint(mod(dble(left%length),N+eps))
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")&
            "left; right=",left%length,right%length,"|"
       if(LMleft>0)write(LOGfile,"("//str(LMleft)//"A)",advance="no")(" ",i=1,LMleft)
       write(LOGfile,"("//str(Mleft)//"A)",advance="no")(trim(Ldot),i=1,Mleft)
       write(LOGfile,"(A)",advance="no")bold_green("*")//bold("|")//bold_red("*")
       write(LOGfile,"("//str(Mright)//"A)",advance="no")(trim(Rdot),i=1,Mright)
       if(LMright>0)write(LOGfile,"("//str(LMright)//"A)",advance="no")(" ",i=1,LMright)
       if(Ltot<=M)then
          write(LOGfile,"(A1,2x,2I4)",advance='yes')"|",left%length+1,right%length+1
       else
          write(LOGfile,"(A1,2x,2I4,2x,I3,2x,A,1x,A)",advance='yes')"|",left%length+1,right%length+1, &
               index,Ldot//"->"//str(N)//"= ;",Rdot//"->"//str(N)//"-"
       endif
    case(1)
       Mleft  = int(left%length/(N+eps))+1
       Mright = int(right%length/(N+eps))+1
       LMleft = 0
       LMright= 0
       index=nint(mod(dble(left%length),N+eps))
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")&
            "left; right=",left%length,right%length,"|"
       if(LMleft>0)write(LOGfile,"("//str(LMleft)//"A)",advance="no")(" ",i=1,LMleft)
       write(LOGfile,"("//str(Mleft)//"A)",advance="no")(trim(Ldot),i=1,Mleft)
       write(LOGfile,"(A)",advance="no")bg_green(">")//"|"//bold_red("*")
       write(LOGfile,"("//str(Mright)//"A)",advance="no")(trim(Rdot),i=1,Mright)
       if(LMright>0)write(LOGfile,"("//str(LMright)//"A)",advance="no")(" ",i=1,LMright)
       if(Ltot<=M)then
          write(LOGfile,"(A1,2x,2I4)",advance='yes')"|",left%length+1,right%length+1
       else
          write(LOGfile,"(A1,2x,2I4,2x,I3,2x,A,1x,A)",advance='yes')"|",left%length+1,right%length+1, &
               index,Ldot//"->"//str(N)//"= ;",Rdot//"->"//str(N)//"-"
       endif
    case(2)
       Mleft  = int(left%length/(N+eps))+1
       Mright = int(right%length/(N+eps))+1
       LMleft = 0
       LMright= 0
       index=nint(mod(dble(left%length),N+eps))
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")&
            "left; right=",left%length,right%length,"|"
       if(LMleft>0)write(LOGfile,"("//str(LMleft)//"A)",advance="no")(" ",i=1,LMleft)
       write(LOGfile,"("//str(Mleft)//"A)",advance="no")(trim(Ldot),i=1,Mleft)
       write(LOGfile,"(A)",advance="no")bold_green("*")//"|"//bg_red("<")
       write(LOGfile,"("//str(Mright)//"A)",advance="no")(trim(Rdot),i=1,Mright)
       if(LMright>0)write(LOGfile,"("//str(LMright)//"A)",advance="no")(" ",i=1,LMright)
       if(Ltot<=M)then
          write(LOGfile,"(A1,2x,2I4)",advance='yes')"|",left%length+1,right%length+1
       else
          write(LOGfile,"(A1,2x,2I4,2x,I3,2x,A,1x,A)",advance='yes')"|",left%length+1,right%length+1, &
               index,Ldot//"->"//str(N)//"= ;",Rdot//"->"//str(N)//"-"
       endif
    end select
    call wait(250)
  end subroutine dmrg_graphic


END MODULE VARS_GLOBAL





