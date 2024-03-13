program hubbard_1d
  USE SCIFOR
  USE DMRG
  implicit none

  integer                               :: Nso
  character(len=64)                     :: finput
  integer                               :: i,unit,iorb
  character(len=1)                      :: DMRGtype
  real(8)                               :: ts(2),Mh(2)
  type(site)                            :: Dot
  type(sparse_matrix)                   :: C,N
  real(8),dimension(:,:),allocatable :: Hloc,Hij

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=(/( -1d0,i=1,2 )/),&
       comment="Hopping amplitudes")
  call parse_input_variable(Mh,"MH",finput,default=(/(0d0,i=1,2 )/),&
       comment="Crystal field splittings")  
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)

  

  if(Norb>2)stop "This code is for Norb<=2. STOP"
  Nso = Nspin*Norb


  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = one*diag([Mh,Mh])
  Dot  = hubbard_site(Hloc)


  !Init DMRG
  call init_dmrg(hubbard_1d_model,ModelDot=Dot)

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



  !H_lr = \sum_{a}h_aa*(C^+_{left,a}@P_left) x C_{right,a}] + H.c.
  function hubbard_1d_model(left,right,states) result(H2)
    type(block)                           :: left
    type(block)                           :: right
    integer,dimension(:),optional         :: states
    type(sparse_matrix),dimension(Norb,2) :: Cl,Cr
    type(sparse_matrix)                   :: P
    type(sparse_matrix)                   :: H2
    integer                               :: ispin,iorb,jorb,io,jo
    integer,dimension(2)                  :: Hdims
    character(len=:),allocatable          :: key
    !
    !>Retrieve operators:
    P = left%operators%op("P")
    do ispin=1,Nspin
       do iorb=1,Norb
          key = "C"//left%okey(iorb,ispin)
          Cl(iorb,ispin) = left%operators%op(key)
          key = "C"//right%okey(iorb,ispin)
          Cr(iorb,ispin) = right%operators%op(key)
       enddo
    enddo
    !
    !> Get H2 dimensions:
    Hdims = shape(left%operators)*shape(right%operators)
    if(present(states))Hdims = [size(states),size(states)]
    !
    !>Build H2:
    call H2%init(Hdims(1),Hdims(2))
    do iorb=1,Norb
       if(present(states))then
          H2 = H2 + ts(iorb)*sp_kron(matmul(Cl(iorb,1)%dgr(),P),Cr(iorb,1),states)
          H2 = H2 + ts(iorb)*sp_kron(matmul(Cl(iorb,2)%dgr(),P),Cr(iorb,2),states)
       else
          H2 = H2 + ts(iorb)*(matmul(Cl(iorb,1)%dgr(),P).x.Cr(iorb,1))
          H2 = H2 + ts(iorb)*(matmul(Cl(iorb,2)%dgr(),P).x.Cr(iorb,2))
       endif
    enddo
    H2 = H2 + H2%dgr()
    !
    !
    !> free memory
    call P%free
    do ispin=1,Nspin
       do iorb=1,Norb
          call Cl(iorb,ispin)%free
          call Cr(iorb,ispin)%free
       enddo
    enddo
  end function hubbard_1d_model


  ! !H_lr = -t \sum_sigma[ (C^+_{l,sigma}@P_l) x C_{r,sigma}]  + H.c.
  ! function hubbard_1d_model_(left,right,states) result(H2)
  !   type(block)                           :: left
  !   type(block)                           :: right
  !   integer,dimension(:),optional         :: states
  !   type(sparse_matrix)                   :: CupL,CdwL
  !   type(sparse_matrix)                   :: CupR,CdwR
  !   type(sparse_matrix)                   :: P
  !   type(sparse_matrix)                   :: H2
  !   P    = left%operators%op("P")
  !   CupL = left%operators%op("C_l1_s1")
  !   CdwL = left%operators%op("C_l1_s2")
  !   CupR = right%operators%op("C_l1_s1")
  !   CdwR = right%operators%op("C_l1_s2")
  !   if(present(states))then
  !      H2   = ts(1)*sp_kron(matmul(CupL%dgr(),P),CupR,states)  + ts(1)*sp_kron(matmul(CdwL%dgr(),P),CdwR,states)
  !      H2   = H2 + H2%dgr()
  !   else
  !      H2   = ts(1)*(matmul(CupL%dgr(),P).x.CupR)  + ts(1)*(matmul(CdwL%dgr(),P).x.CdwR)
  !      H2   = H2 + H2%dgr()
  !   endif
  !   call P%free
  !   call CupL%free
  !   call CdwL%free
  !   call CupR%free
  !   call CdwR%free
  ! end function hubbard_1d_model_


end program hubbard_1d







