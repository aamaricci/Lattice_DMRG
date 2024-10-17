program dmrg_spin_1d
  USE SCIFOR
  USE DMRG
  USE ASSERTING
  implicit none
  character(len=64)                  :: finput
  character(len=1)                   :: DMRGtype
  integer                            :: i,Unit,L
  type(site)                         :: Dot
  type(sparse_matrix)                :: bSz,bSp,SiSj
  real(8),dimension(:,:),allocatable :: Hlr
  real(8),dimension(:),allocatable   :: avSz,x,data,data_
  real(8),dimension(:),allocatable   :: e,s

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call read_input(finput)


  !Init DMRG
  Dot = spin_site(sun=2)
  allocate(Hlr(Nspin*Norb,Nspin*Norb))
  Hlr(1,1) = Jp
  Hlr(2,2) = Jx/2d0
  call init_dmrg(Hlr,ModelDot=Dot)
  call infinite_DMRG()

  !Post-processing and measure quantities:
  call Measure_DMRG(dot%operators%op(key="S_z"),pos=arange(1,Ldmrg),avOp=avSz)

<<<<<<< HEAD

=======
>>>>>>> 7fd23fe (Test testing update)
  !Check energy:
  L = file_length("energyVSleft.length_L40_M20_iDMRG.dmrg")
  allocate(x(L),data(L))
  call sread("energyVSleft.length_L40_M20_iDMRG.dmrg",x,data)
  L = file_length("energy.check")
  deallocate(x)
  allocate(x(L),data_(L))
  call sread("energy.check",x,data_)
  if(size(data)/=size(data_))stop "Energy files have different sizes"
  call assert(data,data_,"E")
  deallocate(x,data,data_)
  !
  !
  L = file_length("SentropyVSleft.length_L40_M20_iDMRG.dmrg")
  allocate(x(L),data(L))
  call sread("SentropyVSleft.length_L40_M20_iDMRG.dmrg",x,data)
  L = file_length("Sentropy.check")
  deallocate(x)
  allocate(x(L),data_(L))
  call sread("Sentropy.check",x,data_)
  if(size(data)/=size(data_))stop "Sentropy files have different sizes"
  call assert(data,data_,"S")
  deallocate(x,data,data_)
  !
  !
<<<<<<< HEAD
  L = file_length("sz.check")
  allocate(data_(L))
  call read_array("sz.check",data_)
  if(size(avSz)/=size(data_))stop "Sz files have different sizes"
  call assert(avSz,data_,"Sz",tol=1d-8)
  deallocate(avSz,data_)
=======
  L = file_length("SzVSj.check")
  allocate(x(L),data_(L))
  call sread("SzVSj.check",x,data_)
  if(size(avSz)/=size(data_))stop "SzVSj files have different sizes"
  call assert(avSz,data_,"Sz")
  deallocate(x,avSz,data_)
>>>>>>> 7fd23fe (Test testing update)



  !Finalize DMRG
  call finalize_dmrg()



end program dmrg_spin_1d





