program BHZ_1d
  USE SCIFOR
  USE DMRG
  implicit none

  integer                               :: Nso
  character(len=64)                     :: finput
  integer                               :: i
  character(len=1)                      :: DMRGtype
  real(8)                               :: target_Qn(2)
  real(8)                               :: eh,mh,lambda
  complex(8),dimension(4,4)             :: Gamma5,Gamma1
  complex(8),dimension(:,:),allocatable :: Hloc,Tx


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(eh,"EH",finput,default=1.d0)
  call parse_input_variable(mh,"MH",finput,default=0.5d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(target_qn,"target_qn",finput,default=[1d0,1d0],&
       comment="Target Sector Occupation [Nup,Ndw] in units [0:2]")
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

  !>Hopping Hamiltonian (i->i+1, right hop direction)
  allocate(Tx(Nso,Nso))
  Tx = -eh/2d0*Gamma5 + xi*lambda/2d0*Gamma1


  !Init DMRG
  call init_dmrg(bhz_1d_model,target_qn,ModelDot=hubbard_site(Hloc))

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


  !H_LR = \sum_{s,ab}[ H_{s,ab} (C^+_{s,a}@P)_L x (C_{s,b})_R]  + H.c.
  !
  !where H_{s,ab} = T_x for the 1d-BHZ model:
  !T_x = -e/2 \Gamma_5 + i\lambda/2 \Gamma_x
  function bhz_1d_model(left,right,states) result(H2)
    type(block)                            :: left
    type(block)                            :: right
    integer,dimension(:),optional          :: states
    type(sparse_matrix),dimension(Nso)     :: Cl,Cr
    type(sparse_matrix)                    :: P
    type(sparse_matrix)                    :: H2
    integer                                :: ispin,iorb,jorb,io,jo
    character(len=:),allocatable           :: lkey,rkey
    !>Retrieve operators:
    P = left%operators%op("P")
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          lkey = "C"//left%okey(iorb,ispin)
          rkey = "C"//right%okey(iorb,ispin)
          Cl(io) = left%operators%op(lkey)
          Cr(io) = right%operators%op(rkey)
       enddo
    enddo
    !>Build connection Hamiltonian:
    H2=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          do jorb=1,Norb
             jo = jorb + (ispin-1)*Norb
             if(Tx(io,jo)==0d0)cycle
             if(present(states))then
                H2 =  H2 + Tx(io,jo)*sp_kron(matmul(Cl(io)%dgr(),P),Cr(jo),states)
             else
                H2 =  H2 + Tx(io,jo)*(matmul(Cl(io)%dgr(),P).x.Cr(jo))
             endif
          enddo
       enddo
    enddo
    H2 = H2 + H2%dgr()
    !> free memory
    call P%free
    do io=1,Nso
       call Cl(io)%free
       call Cr(io)%free
    enddo
  end function bhz_1d_model




end program



