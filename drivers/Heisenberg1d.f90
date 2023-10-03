program heisenberg_1d
  USE SCIFOR
  USE DMRG
  implicit none

  character(len=64) :: finput
  character(len=1)  :: DMRGtype
  integer           :: i
  real(8)           :: target_Sz(1)

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(target_Sz,"target_Sz",finput,default=[0d0],&
       comment="Target Sector Magnetizatio Sz in units [-1:1]")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)


  !Init DMRG
  call init_dmrg(heisenberg_1d_model,target_Sz,ModelDot=spin_onehalf_site())
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


  function heisenberg_1d_model(left,right,states) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    integer,dimension(:),optional :: states
    type(sparse_matrix)           :: Sz1,Sp1
    type(sparse_matrix)           :: Sz2,Sp2
    type(sparse_matrix)           :: H2
    Sz1 = left%operators%op("Sz")
    Sp1 = left%operators%op("Sp")
    Sz2 = right%operators%op("Sz")
    Sp2 = right%operators%op("Sp")
    if(present(states))then
       H2 = Jx/2d0*sp_kron(Sp1,Sp2%dgr(),states) +  Jx/2d0*sp_kron(Sp1%dgr(),Sp2,states)  + Jp*sp_kron(Sz1,Sz2,states)
    else
       H2 = Jx/2d0*(Sp1.x.Sp2%dgr()) +  Jx/2d0*(Sp1%dgr().x.Sp2)  + Jp*(Sz1.x.Sz2)
    endif
    call Sz1%free()
    call Sp1%free()
    call Sz2%free()
    call Sp2%free()
  end function Heisenberg_1d_Model



end program heisenberg_1d






!
!Here we need to construct two things:
!1. the indices of the reduced SB hamiltonian, ie the indices of the row/col
!   corresponding to the sectors contributing to target_Sz. How?
!   One loops over the enlarged blocksys  QNs. For each get the enlarged env QNs such
!   that the sum is the target QN: esys.qn + eenv.qn = targetQN. {This can probably be
!   simplified building a double list [esys.qn, eenv.qn].}
!   For any couple of QNs use the sector map to get the list of states corresponding to
!   each QN of the esys and eenv, say istate and jstate.
!   For any couple of istate, jstate form the index "jstate + (istate-1)*eenv.dim"
!   corresponding to the state in the SB hamiltonian.
!2. the map of the states for the new updated block. How?
!   Labelling the states found at 1. in terms of the QNs of the system. This map
!   will be used later on to block decompose the density matrix, the rotation matrix
!   and so on.
!Store here the states in 1.
! allocate(sb_states(0))
