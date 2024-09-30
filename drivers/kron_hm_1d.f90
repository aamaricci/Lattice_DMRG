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
  type(sparse_matrix)                          :: C,N,Cl,Cr,P
  type(sparse_matrix),allocatable,dimension(:) :: Hsb
  real(8)                                      :: t0
  integer                                      :: m_sb
  real(8),dimension(:,:),allocatable           :: Evecs,Rho,Hloc,Hlr
  real(8),dimension(:),allocatable             :: Evals
  ! integer                                      :: Neigen=2





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




  if(allocated(Hlr))deallocate(Hlr)
  allocate(Hlr(Nso,Nso))
  Hlr = diag([ts(1:Norb),ts(1:Norb)])
  if(Norb==2)Hlr = Hlr + lambda*kron(pauli_0,pauli_x)


  call init_dmrg(Hlr,ModelDot=my_Dot)
  target_qn = DMRG_qn
  !






  print*,""
  print*,""
  print*,"######################################"
  print*,"   o->o + o<-o"
  print*,"######################################"
  print*,""
  print*,""


  write(LOGfile,"(A22,2I12)")"Blocks Length (L-R) = ",left%length,right%length
  left  = block(my_dot)
  right = block(my_dot)
  call enlarge_block(left,my_dot,grow='left')
  call enlarge_block(right,my_dot,grow='right')
  !


  !#################################
  !    Build SUPER-BLOCK Sector
  !#################################
  current_L         = left%length + right%length
  current_target_QN = int(target_qn*current_L*Norb)
  write(LOGfile,"(A22,I12)")"SuperBlock Length = ",current_L
  write(LOGfile,"(A22,"//str(size(current_target_QN))//"F12.7)")"Target_QN = ",current_target_QN
  write(LOGfile,"(A22,G12.7)")"Total         = ",sum(current_target_QN)
  write(LOGfile,"(A22,G12.7)")"Filling       = ",sum(current_target_QN)/current_L
  write(LOGfile,"(A22,G12.7)")"Filling/Norb  = ",sum(current_target_QN)/current_L/Norb
  call sb_get_states()
  print*,size(sb_states)


  !Old Style solution:
  print*,"Old style solution: spH --> H.v --> \psi"
  print*,"######################################"
  sparse_H = .true.
  m_sb = size(sb_states)
  call sb_diag()
  do i=1,size(gs_energy)
     print*,i,gs_energy(i)/2/left%length/Norb
  enddo
  print*,""
  !So it won't work anymore:
  call spHsb%free()




  
  ! !
  ! !Using _new to get energy:
  ! print*,"Solving the same problem using new Method:"
  ! print*,"######################################"
  ! sparse_H = .false.
  ! call sb_build_Hv()
  ! allocate(Evals(Neigen))
  ! allocate(Evecs(size(sb_states),Neigen))
  ! call sp_eigh(spHtimesV_p,evals,evecs,&
  !      3*Neigen,&
  !      500,&
  !      tol=1d-12,&
  !      iverbose=.false.)
  ! do i=1,Neigen
  !    print*,i,Evals(i)/2/left%length/Norb
  ! enddo
  ! deallocate(evals,evecs)
  ! print*,""


end program testEDkron
