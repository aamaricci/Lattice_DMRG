program hubbard_1d
  USE SCIFOR
  USE DMRG
  implicit none

  integer                               :: Nso
  character(len=64)                     :: finput
  integer                               :: i,unit,iorb
  character(len=1)                      :: DMRGtype
  real(8)                               :: target_Qn(2)
  type(site)                            :: Dot
  type(sparse_matrix)                   :: C,N
  real(8)                               :: ts
  complex(8),dimension(:,:),allocatable :: Hloc,Hij

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(target_qn,"target_qn",finput,default=[1d0,1d0],&
       comment="Target Sector Occupation [Nup,Ndw] in units [0:1]")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call parse_input_variable(ts,"TS",finput,default=-1d0,&
       comment="Hopping amplitude")
  call read_input(finput)


  if(Norb>1)stop "This code is for Norb=1. Stop."
  Nso = Nspin*Norb


  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = zero
  Dot  = hubbard_site(Hloc)



  !Init DMRG
  call init_dmrg(hubbard_1d_model,target_qn,ModelDot=Dot)

  !Run DMRG algorithm
  select case(DMRGtype)
  case default;stop "DMRGtype != [Infinite,Finite]"
  case("i","I")
     call infinite_DMRG()
  case("f","F")
     call finite_DMRG()
  end select



  ! !Post-processing and measure quantities:
  ! !Measure <n(a)>
  ! do iorb=1,Norb
  !    call Measure_Op_DMRG(get_N(Dot,iorb),arange(1,Ldmrg/2),file="densVSj_l"//str(iorb))
  ! enddo


  !Finalize DMRG
  call finalize_dmrg()


contains



  !H_lr = \sum_{a}h_aa*(C^+_{left,a}@P_left) x C_{right,a}] + H.c.
  function hubbard_1d_model(left,right,states) result(H2)
    type(block)                           :: left
    type(block)                           :: right
    integer,dimension(:),optional         :: states
    type(sparse_matrix),dimension(2) :: Cl,Cr
    type(sparse_matrix)                   :: P
    type(sparse_matrix)                   :: H2
    integer                               :: ispin,iorb,jorb,io,jo
    integer,dimension(2)                  :: Hdims
    character(len=:),allocatable          :: key
    !
    !>Retrieve operators:
    P     = left%operators%op("P")
    Cl(1) = left%operators%op("C"//left%okey(1,1))
    Cl(2) = left%operators%op("C"//left%okey(1,2))
    Cr(1) = right%operators%op("C"//right%okey(1,1))
    Cr(2) = right%operators%op("C"//right%okey(1,2))
    !
    !>Get H2 dimensions:
    Hdims = shape(left%operators)*shape(right%operators)
    if(present(states))Hdims = [size(states),size(states)]
    !
    !>Build H2
    call H2%init(Hdims(1),Hdims(2))
    if(present(states))then
       H2 = H2 + ts*sp_kron(matmul(Cl(1)%dgr(),P),Cr(1),states)
       H2 = H2 + ts*sp_kron(matmul(Cl(2)%dgr(),P),Cr(2),states)
    else
       H2 = H2 + ts*(matmul(Cl(1)%dgr(),P).x.Cr(1))
       H2 = H2 + ts*(matmul(Cl(2)%dgr(),P).x.Cr(2))
    endif
    H2 = H2 + H2%dgr()
    !
    !> free memory
    call P%free
    call Cl(1)%free
    call Cl(2)%free
    call Cr(1)%free
    call Cr(2)%free
  end function hubbard_1d_model

  
  function get_N(dot,iorb) result(dens)
    type(site)                            :: dot
    integer                               :: iorb,ispin
    type(sparse_matrix)                   :: dens
    type(sparse_matrix)                   :: Cup,Cdw,P
    P    = dot%operators%op("P")
    Cup  = dot%operators%op("C"//dot%okey(iorb,1))
    Cdw  = dot%operators%op("C"//dot%okey(iorb,2))
    dens = matmul(Cup%dgr(),Cup) + matmul(Cdw%dgr(),Cdw)
  end function get_N



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







