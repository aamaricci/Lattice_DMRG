program heisenberg_1d
  USE SCIFOR
  USE DMRG
  implicit none

  character(len=64)   :: finput
  character(len=1)    :: DMRGtype
  integer             :: i,SUN,Unit
  real(8)             :: target_Sz(1),Hvec(3)
  type(site)          :: Dot
  type(sparse_matrix) :: bSz,bSp,SiSj


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
     dot = spin_onehalf_site(one*hvec)
     call init_dmrg(heisenberg_1d_model,target_Sz,ModelDot=Dot)
  case(3)
     dot = spin_one_site(one*Hvec)
     call init_dmrg(heisenberg_1d_model,target_Sz,ModelDot=Dot)
  end select

  !Run DMRG algorithm
  select case(DMRGtype)
  case default;stop "DMRGtype != [Infinite,Finite]"
  case("i","I")
     call infinite_DMRG()
  case("f","F")
     call finite_DMRG()
  end select

  !Post-processing and measure quantities:
  
  !Measure <Sz(i)>
  call Measure_Op_DMRG(dot%operators%op(key='Sz'),arange(1,Ldmrg/2),file="SzVSj")

  
  unit=fopen("SiSjVSj"//str(label_DMRG('u')),append=.true.)
  !Measure <S(i).S(i+1)>
  do i=1,Ldmrg/2-1
     !Get Block Operators at the site I 
     bSz = Build_Op_DMRG(dot%operators%op("Sz"),i)
     bSp = Build_Op_DMRG(dot%operators%op("Sp"),i)
     !Build the correlation:
     SiSj= get_SiSj(bSz,bSp,dot%operators%op("Sz"),dot%operators%op("Sp"))
     !Advance correlation to Psi + measure
     SiSj= Advance_Corr_DMRG(SiSj,i)
     write(unit,*)i,Average_Op_DMRG(SiSj,'l')
  enddo
  close(unit)

  !Finalize DMRG
  call finalize_dmrg()

contains


  function heisenberg_1d_model(left,right,states) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    integer,dimension(:),optional :: states
    type(sparse_matrix)           :: Sz1,Sp1
    type(sparse_matrix)           :: Sz2,Sp2
    type(sparse_matrix)           :: Hp,Hz,H2
    Sz1 = left%operators%op("Sz")
    Sp1 = left%operators%op("Sp")
    Sz2 = right%operators%op("Sz")
    Sp2 = right%operators%op("Sp")
    if(present(states))then
       Hz = Jp*sp_kron(Sz1,Sz2,states)
       Hp = Jx/2d0*sp_kron(Sp1,Sp2%dgr(),states)
    else
       Hz = Jp*(Sz1.x.Sz2)
       Hp = Jx/2d0*(Sp1.x.Sp2%dgr())
    endif
    H2 = Hz + Hp + Hp%dgr()
    call Sz1%free()
    call Sp1%free()
    call Sz2%free()
    call Sp2%free()
  end function Heisenberg_1d_Model


  function get_SiSj(Sz1,Sp1,Sz2,Sp2) result(sisj)
    type(sparse_matrix) :: sisj
    type(sparse_matrix) :: Sz1,Sp1,Sz2,Sp2
    SiSj = 0.5d0*(Sp1.x.Sp2%dgr()) +  0.5d0*(Sp1%dgr().x.Sp2)  + (Sz1.x.Sz2)
  end function get_SiSj

end program





