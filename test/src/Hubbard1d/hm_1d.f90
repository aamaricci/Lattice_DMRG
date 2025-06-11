program hubbard_1d
  USE SCIFOR
  USE DMRG
  USE ASSERTING
  implicit none

  integer                                        :: Nso
  character(len=64)                              :: finput
  integer                                        :: i,unit,iorb,ispin,L
  character(len=1)                               :: DMRGtype
  real(8)                                        :: ts(2),Mh(2)
  type(site),dimension(:),allocatable            :: Dot
  real(8),dimension(:,:),allocatable             :: Hloc,Hlr
  type(sparse_matrix),dimension(:,:),allocatable :: Nop,Cop
  type(sparse_matrix),dimension(:),allocatable   :: dens,docc,s2z
  real(8),dimension(:,:),allocatable             :: avO
  real(8),dimension(:),allocatable               :: x,data,data_
  real(8),dimension(:),allocatable               :: n,d,m,e,s


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=(/( -1d0,i=1,2 )/),&
       comment="Hopping amplitudes")
  call parse_input_variable(Mh,"MH",finput,default=(/(0d0,i=1,2 )/),&
       comment="Crystal field splittings")
  call read_input(finput)


  Nso = Nspin*Norb

  ! allocate(Hloc(Nso,Nso))
  Hloc = diag([Mh(1:Norb),Mh(1:Norb)])


  allocate(Dot(1))
  Dot = electron_site()

  if(allocated(Hlr))deallocate(Hlr)
  allocate(Hlr(Nso,Nso))
  Hlr = diag([ts(1:Norb),ts(1:Norb)])

  call init_dmrg(Hlr,ModelDot=Dot)

  !Run DMRG algorithm
  call infinite_DMRG()


  !Post-processing and measure quantities:
  !Measure <Sz(i)>
  allocate(Cop(Norb,Nspin),Nop(Norb,Nspin))
  do ispin=1,Nspin
     do iorb=1,Norb
        Cop(iorb,ispin) = dot(1)%operators%op(key="C"//dot(1)%okey(iorb,ispin))
        Nop(iorb,ispin) = matmul(Cop(iorb,ispin)%dgr(),Cop(iorb,ispin))
     enddo
  enddo
  allocate(dens(Norb),docc(Norb),s2z(Norb))
  do iorb=1,Norb
     dens(iorb) = Nop(iorb,1)+Nop(iorb,2)
     docc(iorb) = matmul(Nop(iorb,1),Nop(iorb,2))
     s2z(iorb)  = matmul((Nop(iorb,1)-Nop(iorb,2)),(Nop(iorb,1)-Nop(iorb,2)))
  enddo


  call Measure_DMRG([dens,docc,s2z],pos=arange(1,Ldmrg),avOp=avO)


  call save_array("n.out",avO(1,:))
  call save_array("d.out",avO(2,:))
  call save_array("s2z.out",avO(3,:))

  
  !Check energy:
  L = file_length("energyVSleft.length_L40_M40_iDMRG.dmrg")
  allocate(x(L),data(L))
  call sread("energyVSleft.length_L40_M40_iDMRG.dmrg",x,data)
  L = file_length("energy.check")
  allocate(data_(L))
  call read_array("energy.check",data_)
  if(size(data)/=size(data_))stop "Energy files have different sizes"
  call assert(data,data_,"E",1d-6)
  deallocate(x,data,data_)

  !
  L = file_length("SentropyVSleft.length_L40_M40_iDMRG.dmrg")
  allocate(x(L),data(L))
  call sread("SentropyVSleft.length_L40_M40_iDMRG.dmrg",x,data)
  L = file_length("Sentropy.check")
  allocate(data_(L))
  call read_array("Sentropy.check",data_)
  if(size(data)/=size(data_))stop "Sentropy files have different sizes"
  call assert(data,data_,"S",1d-6)
  deallocate(x,data,data_)
  !
  i = 1
  L = file_length("n.check")
  allocate(data_(L))
  call read_array("n.check",data_)
  if(size(avO,2)/=size(data_))stop "N files have different sizes"
  call assert(avO(i,:),data_,"N",1d-6)
  deallocate(data_)
  !
  i = 2
  L = file_length("d.check")
  allocate(data_(L))
  call read_array("d.check",data_)
  if(size(avO,2)/=size(data_))stop "N files have different sizes"
  call assert(avO(i,:),data_,"D",1d-6)
  deallocate(data_)
  !
  i = 3
  L = file_length("s2z.check")
  allocate(data_(L))
  call read_array("s2z.check",data_)
  if(size(avO,2)/=size(data_))stop "S2z files have different sizes"
  call assert(avO(i,:),data_,"S2z",1d-6)
  deallocate(data_)




  !Finalize DMRG
  call finalize_dmrg()


end program hubbard_1d







