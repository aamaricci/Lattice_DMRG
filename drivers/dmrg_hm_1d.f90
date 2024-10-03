program hubbard_1d
  USE SCIFOR
  USE DMRG
  implicit none

  integer                                        :: Nso
  character(len=64)                              :: finput
  integer                                        :: i,unit,iorb,ispin
  character(len=1)                               :: DMRGtype
  real(8)                                        :: ts(2),Mh(2),lambda
  type(site)                                     :: Dot
  real(8),dimension(:,:),allocatable             :: Hloc,Hlr
  type(sparse_matrix),dimension(:,:),allocatable :: N,C
  type(sparse_matrix),dimension(:),allocatable   :: dens,docc,sz,s2z,Mvec

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

  if(allocated(Hlr))deallocate(Hlr)
  allocate(Hlr(Nso,Nso))
  Hlr = diag([ts(1:Norb),ts(1:Norb)])
  if(Norb==2)Hlr = Hlr + lambda*kron(pauli_0,pauli_x)


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
  allocate(C(Norb,Nspin),N(Norb,Nspin))
  do ispin=1,Nspin
     do iorb=1,Norb
        C(iorb,ispin) = dot%operators%op(key="C"//dot%okey(iorb,ispin))
        N(iorb,ispin) = matmul(C(iorb,ispin)%dgr(),C(iorb,ispin))
     enddo
  enddo
  allocate(Mvec(3*Norb),sz(Norb))
  do iorb=1,Norb
     sz(iorb)          = n(iorb,1)-n(iorb,2)
     Mvec(iorb)        = n(iorb,1)+n(iorb,2)
     Mvec(iorb+Norb)   = matmul(n(iorb,1),n(iorb,2))
     Mvec(iorb+2*Norb) = matmul(sz(iorb),sz(iorb))
  enddo


  call Measure_DMRG(Mvec,file="n_d_s2z_l1VSj", pos=arange(1,Ldmrg))



  !Finalize DMRG
  call finalize_dmrg()


end program hubbard_1d







