program hubbard_1d
  USE SCIFOR
  USE DMRG
  implicit none

  character(len=64) :: finput
  integer           :: i
  character(len=1)  :: DMRGtype
  real(8)           :: target_Qn(2)

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(target_qn,"target_qn",finput,default=[1d0,1d0],&
       comment="Target Sector Occupation [Nup,Ndw] in units [0:1]")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)


  !Init DMRG
  call init_dmrg(hubbard_1d_model,target_qn,ModelDot=hubbard_site())
  !Run DMRG algorithm
  select case(DMRGtype)
  case default;stop "DMRGtype != [Infinite,Finite]"
  case("i","I")
     call infinite_DMRG()
  case("f","F")
     call finite_DMRG()
  end select
  !Finalize DMRG
  call finalize_dmrg()


contains


  !H_lr = -t \sum_sigma[ (C^+_{l,sigma}@P_l) x C_{r,sigma}]  + H.c.
  function hubbard_1d_model(left,right,states) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    integer,dimension(:),optional :: states
    type(sparse_matrix)           :: CupL,CdwL
    type(sparse_matrix)           :: CupR,CdwR
    type(sparse_matrix)           :: P
    type(sparse_matrix)           :: H2
    P    = left%operators%op("P")
    CupL = left%operators%op("Cup")
    CdwL = left%operators%op("Cdw")
    CupR = right%operators%op("Cup")
    CdwR = right%operators%op("Cdw")
    if(present(states))then
       H2   = ts(1)*sp_kron(matmul(CupL%t(),P),CupR,states)  + ts(1)*sp_kron(matmul(CdwL%t(),P),CdwR,states)
       H2   = H2 + H2%dgr()
    else
       H2   = ts(1)*(matmul(CupL%t(),P).x.CupR)  + ts(1)*(matmul(CdwL%t(),P).x.CdwR)
       H2   = H2 + H2%dgr()
    endif
    call P%free
    call CupL%free
    call CdwL%free
    call CupR%free
    call CdwR%free
  end function Hubbard_1d_Model




end program hubbard_1d



