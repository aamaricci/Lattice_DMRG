program BHZ_1d
  USE SCIFOR
  USE DMRG
  implicit none

  integer                                        :: Nso
  character(len=64)                              :: finput
  integer                                        :: i,unit,iorb,ispin
  character(len=1)                               :: DMRGtype
  real(8)                                        :: eh,mh,lambda
  type(site)                                     :: Dot
  complex(8),dimension(4,4)                      :: Gamma5,Gamma1
  complex(8),dimension(:,:),allocatable          :: Hloc,Hlr
  type(sparse_matrix),dimension(:,:),allocatable :: N,C
  type(sparse_matrix),dimension(:),allocatable   :: dens,docc,sz,s2z,Mvec

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(eh,"EH",finput,default=1.d0)
  call parse_input_variable(mh,"MH",finput,default=0.5d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)



  if(Nspin/=2.OR.Norb/=2)&
       stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb

  gamma1=kron( pauli_sigma_z, pauli_tau_x)
  gamma5=kron( pauli_sigma_0, pauli_tau_z)

  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = Mh*Gamma5
  Dot  = electron_site(Hloc)

  !>Hopping Hamiltonian (i->i+1, right hop direction)
  if(allocated(Hlr))deallocate(Hlr)
  allocate(Hlr(Nso,Nso))
  Hlr = -eh/2d0*Gamma5 + xi*lambda/2d0*Gamma1


  !Init DMRG
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




end program BHZ_1d



