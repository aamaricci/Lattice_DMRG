program testEDkron
  USE SCIFOR
  USE DMRG, id=>sp_eye
  implicit none

  integer                                      :: Nso,current_L
  character(len=64)                            :: finput
  integer                                      :: i,j,Nsb,k,l,r,m,iorb,jorb,ispin,it,ip,io,jo
  character(len=1)                             :: DMRGtype
  real(8)                                      :: ts(2),Mh,suml,lambda,vh
  type(site)                                   :: my_dot
  type(sparse_matrix),allocatable,dimension(:) :: Hsb
  real(8)                                      :: t0
  integer                                      :: m_sb
  real(8),dimension(:,:),allocatable           :: Evecs,Rho,Hloc
  real(8),dimension(:),allocatable             :: Evals
  integer                                      :: Neigen=2

  type(sparse_matrix),dimension(:,:),allocatable :: N,C
  type(sparse_matrix) :: dens,docc,P,sz,m2



  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=(/( -0.5d0,i=1,2 )/),&
       comment="Hopping amplitudes")
  call parse_input_variable(lambda,"lambda",finput,default=0d0,&
       comment="Hybridization")
  call parse_input_variable(Mh,"MH",finput,default=0d0,&
       comment="Crystal field splittings")
  call parse_input_variable(vh,"vh",finput,default=0d0,&
       comment="local hybridization amplitude")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)


  ! if(Norb/=2)stop "This code is for Norb=2. STOP"
  Nso = Nspin*Norb

  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = 0d0
  if(Norb==2)Hloc = Mh*kron(pauli_0,pauli_z) + vh*kron(pauli_0,pauli_x)

  my_Dot  = electron_site(Hloc)
  call my_dot%show()


  call init_dmrg(hm_1d_model,ModelDot=my_Dot)
  target_qn = DMRG_qn
  !


  !Post-processing and measure quantities:
  !Measure <Sz(i)>
  P = dot%operators%op(key="P")
  allocate(C(Norb,Nspin),N(Norb,Nspin))
  do ispin=1,Nspin
     do iorb=1,Norb
        C(iorb,ispin) = dot%operators%op(key="C"//dot%okey(iorb,ispin))
        n(iorb,ispin) = matmul(C(iorb,ispin)%dgr(),C(iorb,ispin))
     enddo
  enddo

  dens = n(1,1)!+n(1,2)
  ! print*,"DENS:"
  ! call dens%show()

  docc = matmul(n(1,1),n(1,2))
  ! print*,"DOCC:"
  ! call docc%show()
  m2 = matmul((n(1,2)-n(1,1)),(n(1,2)-n(1,1)))

  left=init_left
  right=init_right
  suffix = label_DMRG('i',1)

  do i=1,Ldmrg

     call step_dmrg()
     call write_energy()
     L = left%length+right%length
     call Measure_Op_DMRG(Op=dens,pos=[1,2,L-1,L])
     print*,""
     print*,""
     print*,""
  enddo



  call Measure_Op_DMRG(file='densVSj',Op=dens)
  ! call Measure_Op_DMRG(file='doccVSj',Op=docc)
  ! call Measure_Op_DMRG(file='mVSj',Op=m2)
  ! suffix = label_DMRG('i',1)
  ! print*,""
  ! print*,""
  ! print*,""
  ! print*,""


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







end program testEDkron
