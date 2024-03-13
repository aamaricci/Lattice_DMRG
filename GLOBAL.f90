MODULE GLOBAL
  USE SCIFOR
  USE INPUT_VARS
  USE AUX_FUNCS
  USE MATRIX_SPARSE, id=>sp_eye
  USE MATRIX_BLOCKS
  USE TUPLE_BASIS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
  implicit none


  abstract interface
     function Hconnect(left,right,states)
       USE BLOCKS
       USE MATRIX_SPARSE
       type(block)                   :: left
       type(block)                   :: right
       integer,dimension(:),optional :: states
       type(sparse_matrix)           :: Hconnect
     end function Hconnect
  end interface
  procedure(Hconnect),pointer,public    :: H2model=>null()
  !
  type(sparse_matrix)                   :: spHsb
  !
  real(8)                               :: truncation_error_sys,truncation_error_env
  type(site)                            :: dot
  character(len=:),allocatable          :: suffix
  real(8),dimension(:),allocatable      :: target_Qn,current_target_QN
  type(block)                           :: init_sys,init_env
  logical                               :: init_called=.false.
  real(8),dimension(:),allocatable      :: energies
  real(8),dimension(:),allocatable      :: rho_sys_evals
  real(8),dimension(:),allocatable      :: rho_env_evals
  type(blocks_matrix)                   :: psi_sys,rho_sys
  type(blocks_matrix)                   :: psi_env,rho_env
  !GLOBAL SYS & ENV 
  type(block)                           :: sys,env



contains



  subroutine dmrg_graphic(label)
    integer                  :: label
    integer                  :: i,N,Msys,Menv,LMsys,LMenv,index,Ltot
    character(:),allocatable :: Ldot,Rdot
    real(8)                  :: eps=1d-6
    integer                  :: M=50
    !
    call wait(50)
    !call system("clear")
    call execute_command_line("clear")
    Ltot = Ldmrg/2
    Ldot = bold_green('=')
    Rdot = bold_red('-')
    ! if(Ltot>M)then
    !    Ldot = bg_green('=')
    !    Rdot = bg_red('-')
    ! endif
    !
    N = int(Ltot/(M+eps))+1
    !
    select case(label)
    case default; stop "dmrg_graphic error: label != 1(L),2(R)"
    case(0)
       Msys = int(sys%length/(N+eps))+1
       Menv = int(env%length/(N+eps))+1
       LMsys= Ltot/N-Msys
       LMenv= Ltot/N-Menv
       index=nint(mod(dble(sys%length),N+eps))
       write(LOGfile,*)""
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")"sys; env=",sys%length,env%length,"|"
       if(LMsys>0)write(LOGfile,"("//str(LMsys)//"A)",advance="no")(" ",i=1,LMsys)
       write(LOGfile,"("//str(Msys)//"A)",advance="no")(trim(Ldot),i=1,Msys)
       write(LOGfile,"(A)",advance="no")bold_green("*")//bold("|")//bold_red("*")
       write(LOGfile,"("//str(Menv)//"A)",advance="no")(trim(Rdot),i=1,Menv)
       if(LMenv>0)write(LOGfile,"("//str(LMenv)//"A)",advance="no")(" ",i=1,LMenv)
    case(1)
       Msys = int(sys%length/(N+eps))+1
       Menv = int(env%length/(N+eps))+1
       LMsys= 0
       LMenv= 0
       index=nint(mod(dble(sys%length),N+eps))
       write(LOGfile,*)""
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")"sys; env=",sys%length,env%length,"|"
       if(LMsys>0)write(LOGfile,"("//str(LMsys)//"A)",advance="no")(" ",i=1,LMsys)
       write(LOGfile,"("//str(Msys)//"A)",advance="no")(trim(Ldot),i=1,Msys)
       write(LOGfile,"(A)",advance="no")bg_green(">")//"|"//bold_red("*")
       write(LOGfile,"("//str(Menv)//"A)",advance="no")(trim(Rdot),i=1,Menv)
       if(LMenv>0)write(LOGfile,"("//str(LMenv)//"A)",advance="no")(" ",i=1,LMenv)
    case(2)
       Msys = int(env%length/(N+eps))+1
       Menv = int(sys%length/(N+eps))+1
       LMsys= 0
       LMenv= 0
       index=nint(mod(dble(env%length),N+eps))
       write(LOGfile,*)""
       write(LOGfile,"(A,2I4,2x,A1)",advance="no")"sys; env=",sys%length,env%length,"|"
       if(LMsys>0)write(LOGfile,"("//str(LMsys)//"A)",advance="no")(" ",i=1,LMsys)
       write(LOGfile,"("//str(Msys)//"A)",advance="no")(trim(Ldot),i=1,Msys)
       write(LOGfile,"(A)",advance="no")bold_green("*")//"|"//bg_red("<")
       write(LOGfile,"("//str(Menv)//"A)",advance="no")(trim(Rdot),i=1,Menv)
       if(LMenv>0)write(LOGfile,"("//str(LMenv)//"A)",advance="no")(" ",i=1,LMenv)
    end select
    if(Ltot<=M)then
       write(LOGfile,"(A1)",advance='yes')"|"
    else
       write(LOGfile,"(A1,I3,2x,A,1x,A)",advance='yes')"|",index,Ldot//"->"//str(N)//"= ;",Rdot//"->"//str(N)//"-"
    endif
    call wait(150)
  end subroutine dmrg_graphic


END MODULE GLOBAL





