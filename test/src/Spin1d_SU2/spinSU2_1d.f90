program dmrg_spin_1d
  USE SCIFOR
  USE DMRG
  USE ASSERTING
  implicit none
  character(len=64)                   :: finput
  character(len=1)                    :: DMRGtype
  integer                             :: i,Unit,L
  type(site),dimension(:),allocatable :: Dot
  type(sparse_matrix)                 :: bSz,bSp,SiSj
  real(8),dimension(:,:),allocatable  :: Hlr
  real(8),dimension(:),allocatable    :: avSz,x,data,data_
  real(8),dimension(:),allocatable    :: e,s

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call read_input(finput)


  !Init DMRG
  allocate(Dot(1))              !Homogeneous system
  Dot = spin_site(sun=2)        !spin_site handles array!
  Hlr = diag([Jp,Jx/2d0])       !implicit fortran allocation
  call init_dmrg(Hlr,ModelDot=Dot)
  call infinite_DMRG()

  !Post-processing and measure quantities:
<<<<<<< HEAD
  call Measure_DMRG(dot%operators%op(key="S_z"),pos=arange(1,Ldmrg),avOp=avSz)

<<<<<<< HEAD
<<<<<<< HEAD
=======
  call Measure_DMRG(dot(1)%operators%op(key="S_z"),pos=arange(1,Ldmrg),avOp=avSz)
>>>>>>> 8cabf7d (Although we have included MATRIX_GRAPH the code now)

=======
>>>>>>> 7fd23fe (Test testing update)
=======

>>>>>>> 56ab977 (test testing 2)
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
<<<<<<< HEAD
=======
>>>>>>> 56ab977 (test testing 2)
  L = file_length("sz.check")
  allocate(data_(L))
  call read_array("sz.check",data_)
  if(size(avSz)/=size(data_))stop "Sz files have different sizes"
  call assert(avSz,data_,"Sz",tol=1d-8)
  deallocate(avSz,data_)
<<<<<<< HEAD
=======
  L = file_length("SzVSj.check")
  allocate(x(L),data_(L))
  call sread("SzVSj.check",x,data_)
  if(size(avSz)/=size(data_))stop "SzVSj files have different sizes"
  call assert(avSz,data_,"Sz")
  deallocate(x,avSz,data_)
>>>>>>> 7fd23fe (Test testing update)
=======
>>>>>>> 56ab977 (test testing 2)



  !Finalize DMRG
  call finalize_dmrg()



end program dmrg_spin_1d





