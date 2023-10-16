program heisenberg_1d
  USE SCIFOR
  USE DMRG
  implicit none

  character(len=64)   :: finput
  character(len=1)    :: DMRGtype
  integer             :: i,SUN,Unit
  real(8)             :: target_Sz(1),Hvec(3)
  !
  type(site)          :: Dot
  type(sparse_matrix) :: Op,bSz,bSp
  real(8),allocatable :: vals(:)
  real(8)             :: val


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(SUN,"SUN",finput,default=2,&
       comment="Spin SU(N) value. 2=> spin 1/2, 3=> spin 1")
  call parse_input_variable(target_Sz,"target_Sz",finput,default=[0d0],&
       comment="Target Sector Magnetizatio Sz in units [-1:1]")
  call parse_input_variable(Hvec,"Hvec",finput,default=[0d0,0d0,0d0],&
       comment="Target Sector Magnetizatio Sz in units [-1:1]")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)


  !Init DMRG
  select case(SUN)
  case default;stop "SU(N) Spin. Allowed values N=2,3 => Spin 1/2, Spin 1"
  case(2)
     dot = spin_onehalf_site(hvec)
     call init_dmrg(heisenberg_1d_model,target_Sz,ModelDot=spin_onehalf_site(hvec))
  case(3)
     dot = spin_one_site(Hvec)
     call init_dmrg(heisenberg_1d_model,target_Sz,ModelDot=spin_one_site(Hvec))
  end select

  !Run DMRG algorithm
  select case(DMRGtype)
  case default;stop "DMRGtype != [Infinite,Finite]"
  case("i","I")
     call infinite_DMRG()
  case("f","F")
     call finite_DMRG()
  end select

  !Measure Sz
  unit=fopen("SzVSj.dmrg")
  do i=1,Ldmrg/2
     vals = measure_local_dmrg([dot%operators%op(key='Sz')],[i])
     write(unit,*)i,vals
  enddo
  close(unit)



  ! do i=1,2
  !    bSz = construct_op_DMRG(dot%operators%op(key='Sz'),i)
  !    bSp = construct_op_DMRG(dot%operators%op(key='Sp'),i)
  !    !Build SiSj
  !    Op = heisenberg_1d_SiSj(bSz,bSp,dot%operators%op(key='Sz'),dot%operators%op(key='Sp'))
  !    Op = actualize_op_DMRG(Op,i+1)
  !    val= measure_op_DMRG(Op,i+1)
  !    write(*,*)i,val
  !    write(200,*)i,val
  ! enddo

  
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



  function heisenberg_1d_SiSj(lSz,lSp,rSz,rSp) result(SiSj)
    type(sparse_matrix)           :: lSz,lSp
    type(sparse_matrix)           :: rSz,rSp
    type(sparse_matrix)           :: SiSj
    ! if(present(states))then
    !    H2 = Jx/2d0*sp_kron(Sp1,Sp2%dgr(),states) +  Jx/2d0*sp_kron(Sp1%dgr(),Sp2,states)  + Jp*sp_kron(Sz1,Sz2,states)
    ! else
    SiSj = 0.5d0*(lSp.x.rSp%dgr()) + 0.5d0*(lSp%dgr().x.rSp)  + (lSz.x.rSz)
    ! endif
  end function heisenberg_1d_SiSj

end program heisenberg_1d





