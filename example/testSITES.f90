program testSITES
  USE SCIFOR,  id => zeye
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  implicit none


  type(site) :: my_site,a,b
  real(8),dimension(2,2),parameter   :: Hzero=reshape([0d0,0d0,0d0,0d0],[2,2])
  real(8),dimension(2,2),parameter   :: Sz=dble(pauli_z)
  real(8),dimension(2,2),parameter   :: Sx=dble(pauli_x)
  real(8),dimension(2,2),parameter   :: Splus=reshape([0d0,0d0,1d0,0d0],[2,2])
  real(8),dimension(4,4)             :: Gamma13,Gamma03

  Gamma13=kron(Sx,Sz)
  Gamma03=kron(eye(2),Sz)

  
  my_site = site(&
       dim      = 2, &
       sectors  = [sectors_list([-0.5d0,0.5d0])],&
       operators= operators_list(['H0','Sz','Sp'],&
       [sparse(Hzero),sparse(Sz),sparse(Splus)]))
  print*,"Is site valid:",my_site%is_valid()
  call my_site%show()



  a = my_site

  print*,"Test ="
  call a%show()
  call a%put("G5",sparse(Gamma03))
  print*,"Is site valid:",my_site%is_valid()


  print*,"Test PAULI SITE"
  b = pauli_site()
  call b%show()
  call b%free()

  b = hubbard_site(0d0,0d0)
  call b%show()
  print*,"Is site valid:",b%is_valid()
  call b%free()
end program testSITES
