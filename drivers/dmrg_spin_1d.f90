program dmrg_spin_1d
  USE SCIFOR
  USE DMRG
  implicit none
  character(len=64)   :: finput
  character(len=1)    :: DMRGtype
  integer             :: i,SUN,Unit,pos
  real(8)             :: Hvec(3)
  type(site)          :: Dot
  type(sparse_matrix) :: bSz,bSp,SiSj
  real(8),dimension(:,:),allocatable :: Hlr


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(SUN,"SUN",finput,default=2,&
       comment="Spin SU(N) value. 2=> spin 1/2, 3=> spin 1")
  call parse_input_variable(Hvec,"Hvec",finput,default=[0d0,0d0,0d0],&
       comment="Magnetization Vector field")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)


  !Init DMRG
  print*,hvec
  Dot = spin_site(sun=SUN,Hvec=Hvec)
  call Dot%show()


  if(allocated(Hlr))deallocate(Hlr)
  allocate(Hlr(Nspin*Norb,Nspin*Norb))
  Hlr(1,1) = Jp
  Hlr(2,2) = Jx/2d0
  call init_dmrg(Hlr,ModelDot=Dot)


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
  call Measure_DMRG(dot%operators%op(key="S_z"),file="SzVSj")

  !Measure <S(i).S(i+1)>
  unit=fopen("SiSjVSj"//str(label_DMRG('u')),append=.true.)
  call Init_measure_dmrg()
  do pos=1,Ldmrg-1
     bSz = Build_Op_DMRG(dot%operators%op("S_z"),pos,set_basis=.true.)
     bSp = Build_Op_DMRG(dot%operators%op("S_p"),pos,set_basis=.true.)
     SiSj= get_SiSj(bSz,bSp,dot%operators%op("S_z"),dot%operators%op("S_p"))
     SiSj= Advance_Corr_DMRG(SiSj,pos)
     write(unit,*)pos,Average_Op_DMRG(SiSj,pos)
  enddo
  call End_measure_dmrg()
  close(unit)


  !Finalize DMRG
  call finalize_dmrg()

contains


  function get_SiSj(Sz1,Sp1,Sz2,Sp2) result(sisj)
    type(sparse_matrix) :: sisj
    type(sparse_matrix) :: Sz1,Sp1,Sz2,Sp2
    SiSj = 0.5d0*(Sp1.x.Sp2%dgr()) +  0.5d0*(Sp1%dgr().x.Sp2)  + (Sz1.x.Sz2)
  end function get_SiSj


end program dmrg_spin_1d





