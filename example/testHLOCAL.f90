program testHLOCAL
  USE SCIFOR
  USE HLOCAL
  implicit none


  integer :: Norb,Nspin
  real(8),dimension(:,:),allocatable   :: Docc,Cop,CDGop,Dens,Hlocal
  real(8),dimension(:,:,:),allocatable :: h0
  real(8),dimension(:,:),allocatable :: CdgC

  Nspin=1
  Norb =1
  call Init_LocalFock_Sectors(Norb,Nspin)

  allocate(h0(nspin,norb,norb));h0=0d0
  print*,""
  print*,"H:"
  Hlocal = build_Hlocal_operator(h0,xmu=0d0,uloc=5d0)
  call print_mat(Hlocal)


  print*,""
  print*,"C:"
  Cop = build_C_operator(ispin=1,iorb=1)
  call print_mat(Cop)


  print*,""
  print*,"CDG:"
  CDGop = build_CDG_operator(ispin=1,iorb=1)
  call print_mat(CDGop)

  print*,""
  print*,"CDG:"
  CDGop = build_CDG_operator(ispin=2,iorb=1)
  call print_mat(CDGop)


  print*,""
  print*,"Dens_up:"
  Dens = build_Dens_operator(ispin=1,iorb=1)
  call print_mat(dens)
  print*,""
  print*,"Dens_dw:"  
  Dens = build_Dens_operator(ispin=2,iorb=1)
  call print_mat(dens)

  print*,""
  print*,"Docc:"
  docc = build_Docc_operator(iorb=1)
  call print_mat(docc)

  print*,""
  CdgC = build_CdgC_operator(ispin=1,iorb=1,jorb=2)
  call print_mat(CdgC)

contains


  subroutine print_mat(M)
    real(8),dimension(:,:) :: M
    integer :: i,j
    do i=1,size(M,1)
       write(*,"("//str(size(M,2))//"(F4.1,1x))")(M(i,j),j=1,size(M,2))
    enddo
  end subroutine print_mat



end program testHLOCAL
