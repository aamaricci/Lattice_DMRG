program dmrg_spin_1d
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
       comment="Target Sector MagnetizationSz in units [-1:1]")
  call parse_input_variable(Hvec,"Hvec",finput,default=[0d0,0d0,0d0],&
       comment="Magnetization Vector field")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)


  !Init DMRG
  Dot = spin_site(sun=SUN,Hvec=Hvec)
  call Dot%show()

  call init_dmrg(heisenberg_1d_model,ModelDot=Dot)


  !Run DMRG algorithm
  select case(DMRGtype)
  case default;stop "DMRGtype != [Infinite,Finite]"
  case("i","I")
     call infinite_DMRG()
  case("f","F")
     call finite_DMRG()
  end select

  ! !Post-processing and measure quantities:
  ! !Measure <Sz(i)>
  ! call Measure_Op_DMRG(dot%operators%op(key='Sz'),arange(1,Ldmrg/2),file="SzVSj")


  ! unit=fopen("SiSjVSj"//str(label_DMRG('u')),append=.true.)
  ! !Measure <S(i).S(i+1)>
  ! do i=1,Ldmrg/2-1
  !    !Get Block Operators at the site I 
  !    bSz = Build_Op_DMRG(dot%operators%op("Sz"),i)
  !    bSp = Build_Op_DMRG(dot%operators%op("Sp"),i)
  !    !Build the correlation:
  !    SiSj= get_SiSj(bSz,bSp,dot%operators%op("Sz"),dot%operators%op("Sp"))
  !    !Advance correlation to Psi + measure
  !    SiSj= Advance_Corr_DMRG(SiSj,i)
  !    write(unit,*)i,Average_Op_DMRG(SiSj,'l')
  ! enddo
  ! close(unit)

  !Finalize DMRG
  call finalize_dmrg()

contains


  function heisenberg_1d_model(left,right) result(Hlr)
    type(block)                   :: left
    type(block)                   :: right
    real(8),dimension(:,:),allocatable :: Hlr
    if(allocated(Hlr))deallocate(Hlr)
    allocate(Hlr(Nspin*Norb,Nspin*Norb))
    !
    !workout local part, like random local field
    !if(left%Dim==1 AND right%Dim==1) then operate over local H
    !if(left%Dim==1 OR right%Dim==1) then operate over local H
    Hlr(1,1) = Jp
    Hlr(2,2) = Jx/2d0
  end function Heisenberg_1d_Model


  function get_SiSj(Sz1,Sp1,Sz2,Sp2) result(sisj)
    type(sparse_matrix) :: sisj
    type(sparse_matrix) :: Sz1,Sp1,Sz2,Sp2
    SiSj = 0.5d0*(Sp1.x.Sp2%dgr()) +  0.5d0*(Sp1%dgr().x.Sp2)  + (Sz1.x.Sz2)
  end function get_SiSj


end program dmrg_spin_1d





