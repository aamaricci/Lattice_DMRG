program testEDkron
  USE SCIFOR
  USE DMRG, id=>sp_eye
  USE OMP_LIB, only: omp_get_wtime
  implicit none

  integer                                      :: Nso
  character(len=64)                            :: finput
  integer                                      :: i,j,Nsb,SUN
  character(len=1)                             :: DMRGtype
  integer                                      :: current_L,unit,pos,N
  type(site)                                   :: my_dot
  type(sparse_matrix)                          :: C,Cl,Cr,P
  integer                                      :: m_sb
  type(sparse_matrix) :: bSz,bSp,SiSj,OpSz
  real(8)             :: Hvec(3)


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(SUN,"SUN",finput,default=2,&
       comment="Spin SU(N) value. 2=> spin 1/2, 3=> spin 1")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call parse_input_variable(Hvec,"Hvec",finput,default=[0d0,0d0,0d0],&
       comment="Magnetization Vector field")
  call read_input(finput)

  print*,hvec
  my_dot = spin_site(sun=SUN,hvec=hvec)
  call my_dot%show()

  OpSz = my_dot%operators%op(key="S_z")
  OpSz = matmul(OpSz,OpSz)
  
  !Init DMRG
  call init_dmrg(spin_1d_hmodel,ModelDot=my_dot)
  target_qn = DMRG_qn
  !

  left=init_left
  right=init_right
  suffix = label_DMRG('i',1)


  do i=1,Ldmrg
     call step_dmrg()
     call write_energy()
     ! call write_truncation()
     call Measure_Op_DMRG(Op=OpSz,pos=arange(1,left%length+right%length))
     print*,""
     print*,""
  enddo


  !Post-processing and measure quantities:
  !Measure <Sz(i)>
  call Measure_Op_DMRG(file="sz2VSj",Op=OpSz)
  call Measure_Op_DMRG(file="szVSj",Op=dot%operators%op(key="S_z"))

  ! !Measure <S(i).S(i+1)>
  ! call Init_measure_dmrg()
  ! unit=fopen("SiSjVSj"//str(label_DMRG('u')),append=.true.)
  ! do pos=1,N/2-1
  !    bSz = Build_Op_DMRG(dot%operators%op("S_z"),pos,set_basis=.true.)
  !    bSp = Build_Op_DMRG(dot%operators%op("S_p"),pos,set_basis=.true.)
  !    SiSj= get_SiSj(bSz,bSp,dot%operators%op("S_z"),dot%operators%op("S_p"))
  !    SiSj= Advance_Corr_DMRG(SiSj,pos)
  !    write(unit,*)pos,Average_Op_DMRG(SiSj,pos)
  ! enddo
  ! close(unit)
  ! call End_measure_dmrg()

  call finalize_dmrg()


contains



  !This is user defined Function to be passed to the SYSTEM
  function spin_1d_hmodel(left,right) result(Hlr)
    type(block)                        :: left,right
    real(8),dimension(:,:),allocatable :: Hlr
    if(allocated(Hlr))deallocate(Hlr)
    allocate(Hlr(Nspin*Norb,Nspin*Norb))
    !
    !workout local part, like random local field
    !if(left%Dim==1 AND right%Dim==1) then operate over local H
    !if(left%Dim==1 OR right%Dim==1) then operate over local H
    Hlr(1,1) = Jp
    Hlr(2,2) = Jx/2d0
  end function spin_1d_hmodel



  function get_SiSj(Sz1,Sp1,Sz2,Sp2) result(sisj)
    type(sparse_matrix) :: sisj
    type(sparse_matrix) :: Sz1,Sp1,Sz2,Sp2
    SiSj = 0.5d0*(Sp1.x.Sp2%dgr()) +  0.5d0*(Sp1%dgr().x.Sp2)  + (Sz1.x.Sz2)
  end function get_SiSj


end program testEDkron
