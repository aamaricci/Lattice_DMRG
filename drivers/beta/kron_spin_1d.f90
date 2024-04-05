program testEDkron
  USE SCIFOR
  USE DMRG, id=>sp_eye
  USE OMP_LIB, only: omp_get_wtime
  implicit none

  integer                                      :: Nso
  character(len=64)                            :: finput
  integer                                      :: i,j,Nsb,SUN
  character(len=1)                             :: DMRGtype
  real(8),dimension(:),allocatable             :: Vec,Hvec,Hvec_
  integer                                      :: current_L
  type(site)                                   :: my_dot
  type(sparse_matrix)                          :: C,N,Cl,Cr,P
  type(sparse_matrix),allocatable,dimension(:) :: Hsb
  real(8)                                      :: t0
  integer                                      :: m_sb
  real(8),dimension(:,:),allocatable           :: Evecs,Rho,Hloc
  real(8),dimension(:),allocatable             :: Evals
  integer                                      :: Neigen=2

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(SUN,"SUN",finput,default=2,&
       comment="Spin SU(N) value. 2=> spin 1/2, 3=> spin 1")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)



  my_dot = spin_site(sun=SUN)
  call my_dot%show()


  !Init DMRG
  call init_dmrg(spin_1d_hmodel,ModelDot=my_dot)
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
  t0 = omp_get_wtime()
  sparse_H = .true.
  call sb_build_Hv()
  allocate(Evals(Neigen))
  allocate(Evecs(size(sb_states),Neigen))
  call sp_eigh(spMatVec_sparse_main,evals,evecs,&
       3*Neigen,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,Neigen
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""
  !So it won't work anymore:
  call spHsb%free()
  print*,omp_get_wtime() - t0


  !
  !Using _new to get energy:
  print*,"Solving the same problem using new Method:"
  print*,"######################################"
  t0 = omp_get_wtime()
  sparse_H = .false.
  call sb_build_Hv()
  allocate(Evals(Neigen))
  allocate(Evecs(size(sb_states),Neigen))
  call sp_eigh(spMatVec_direct_main,evals,evecs,&
       3*Neigen,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,Neigen
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""
  print*,omp_get_wtime() - t0
















contains



  !This is user defined Function to be passed to the SYSTEM
  function spin_1d_hmodel(left,right) result(Hlr)
    type(block)                        :: left,right
    real(8),dimension(:,:),allocatable :: Hlr
    if(allocated(Hlr))deallocate(Hlr)
    allocate(Hlr(Nspin*Norb,Nspin*Norb))
    !
    !workout local part, like random local field
    !if(left%Dim==1 AND right%Dim==1) then operate over local H
    !if(left%Dim==1 OR right%Dim==1) then operate over local H
    Hlr(1,1) = Jp
    Hlr(2,2) = Jx/2d0
  end function spin_1d_hmodel


end program testEDkron
