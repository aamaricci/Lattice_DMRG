program hubbard_1d
  USE SCIFOR
  USE DMRG
  implicit none

  integer                            :: Nso
  character(len=64)                  :: finput
  integer                            :: i,unit,iorb
  character(len=1)                   :: DMRGtype
  real(8)                            :: ts(2),Mh(2),lambda
  type(site)                         :: Dot
  real(8),dimension(:,:),allocatable :: Hloc

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=(/( -0.5d0,i=1,2 )/),&
       comment="Hopping amplitudes")
  call parse_input_variable(Mh,"MH",finput,default=(/(0d0,i=1,2 )/),&
       comment="Crystal field splittings")
  call parse_input_variable(lambda,"LAMBDA",finput,default=0d0,&
       comment="off-diagonal amplitude")  
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)


  if(Norb>2)stop "This code is for Norb<=2. STOP"
  Nso = Nspin*Norb


  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = diag([Mh(1:Norb),Mh(1:Norb)])
  Dot  = electron_site(Hloc)


  !Init DMRG
  call init_dmrg(hm_1d_model,ModelDot=Dot)

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



  !This is user defined Function to be passed to the SYSTEM
  function hm_1d_model(left,right) result(Hlr)
    type(block)                        :: left
    type(block)                        :: right
    real(8),dimension(:,:),allocatable :: Hlr
    !
    if(allocated(Hlr))deallocate(Hlr)
    allocate(Hlr(Nspin*Norb,Nspin*Norb))
    !
    !workout local part, like random local field
    !if(left%Dim==1 AND right%Dim==1) then operate over local H
    !if(left%Dim==1 OR right%Dim==1) then operate over local H
    Hlr = diag([ts(1:Norb),ts(1:Norb)])
    if(Norb==2)Hlr = Hlr + lambda*kron(pauli_0,pauli_x)
  end function hm_1d_model


end program hubbard_1d







