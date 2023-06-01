program testEDkron
  USE SCIFOR
  USE AUX_FUNCS
  USE MATRIX_SPARSE, id=>sp_eye
  USE MATRIX_BLOCKS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
  implicit none

  character(len=64)                  :: finput
  real(8)                            :: ts
  integer                            :: Lmax,i,m,Neigen,current_L
  integer                            :: lanc_ncv_factor
  integer                            :: lanc_ncv_add
  integer                            :: lanc_niter
  real(8)                            :: lanc_tolerance
  integer                            :: lanc_threshold
  type(block)                        :: my_block,dimer,trimer,dimerL,dimerR,trimerL,trimerR
  type(sparse_matrix)                :: spHsb,spH
  type(site)                         :: dot
  real(8)                            :: gs_energy,target_Sz
  integer                            :: unit
  integer                            :: model_d=4
  real(8),dimension(:,:),allocatable :: Hmatrix,Evecs
  real(8),dimension(:),allocatable   :: Evals
  type(sparse_matrix)       :: op


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"ts",finput,default=-1d0,comment="Hopping amplitude")
  target_Sz=0d0




  !
  !Init the single dot structure:
  dot = hubbard_site_ud(0d0,0d0)
  !
  !Init block from single dot
  my_block=block(dot)
  call my_block%show()



  print*,"o->o"
  dimer = enlarge_block(my_block,dot,grow='left')
  spH = (dimer%operators%op("Hd")) + (Id(2**2).x.dimer%operators%op("Hup")) + (dimer%operators%op("Hdw").x.Id(2**2))
  Hmatrix = as_matrix(spH)!dimer%operators%op("H"))
  allocate(Evals(4**2))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""



  print*,"o-o->o"
  trimer = enlarge_block(dimer,dot,grow='left')
  spH = (trimer%operators%op("Hd")) + (Id(2**3).x.trimer%operators%op("Hup")) + (trimer%operators%op("Hdw").x.Id(2**3))
  Hmatrix = as_matrix(spH)!dimer%operators%op("H"))
  allocate(Evals(4**3))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""

  print*,"o-o-o->o"
  dimer = enlarge_block(trimer,dot,grow='left')  
  spH =   (Id(2**4).x.dimer%operators%op("Hup")) + (dimer%operators%op("Hdw").x.Id(2**4))  + (dimer%operators%op("Hd"))
  Hmatrix = as_matrix(spH)!dimer%operators%op("H"))
  allocate(Evals(4**4))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""


contains



  !We want to arrange states such that the basis is {_n .... _1}_dw {_n ... _1}_up
  ! H = H_d + 1_dw x h_up + h_dw x 1_up
  ! we need to grow accordingly each element of the dictionary:
  ! H_d  = H_d(block) x 1(dot) + 1(block) x H_d(dot) : contains the local parts
  ! h_up = h_up(block) x 1(dot) + H2(up)

  function enlarge_block(self,dot,grow) result(enl_self)
    type(block),intent(inout)        :: self
    type(site),intent(inout)         :: dot
    character(len=*),optional        :: grow
    character(len=16)                :: grow_
    type(block)                      :: enl_self
    integer                          :: mblock,len
    real(8),dimension(:),allocatable :: self_basis, dot_basis
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    mblock =  self%dim
    len    =  self%length
    !
    enl_self%length = len + 1
    enl_self%dim    = mblock*model_d
    !
    ! select case(str(grow_))
    ! case ("left","l")
    call enl_self%put("Hd", (self%operators%op("Hd").x.id(model_d)) + (id(mblock).x.dot%operators%op("Hd")))
    call enl_self%put("Hup", (self%operators%op("Hup").x.id(2)) + H2model(self,as_block(dot)))
    call enl_self%put("Hdw", (self%operators%op("Hdw").x.id(2)) + H2model(self,as_block(dot)))
    call enl_self%put("Cp", Id(2**len).x.dot%operators%op("Cp"))
    call enl_self%put("P"  , Id(2**len).x.dot%operators%op("P"))
    !   case ("right","r")
    !      call enl_self%put("H", (id(model_d).x.self%operators%op("H")) +  (dot%operators%op("H").x.id(mblock)) + &
    !           H2model(as_block(dot),self))
    !      call enl_self%put("Cup", dot%operators%op("Cup").x.Id(mblock))
    !      call enl_self%put("Cdw", dot%operators%op("Cdw").x.Id(mblock))
    !      call enl_self%put("P"  , dot%operators%op("P").x.Id(mblock))
    !   end select
    !   !
    allocate(enl_self%sectors(2))
    self_basis = self%sectors(1)%basis()
    dot_basis  = dot%sectors(1)%basis()
    call enl_self%set_sectors( indx=1, vec=outsum(self_basis,dot_basis) )
    self_basis = self%sectors(2)%basis()
    dot_basis  = dot%sectors(2)%basis()
    call enl_self%set_sectors( indx=2, vec=outsum(self_basis,dot_basis) )
    !   ! !
    deallocate(self_basis,dot_basis)
  end function enlarge_block



  !H_lr = -t (C^+_{l}@P_l) x C_{r}  + H.c.
  function H2model(left,right) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    type(sparse_matrix)           :: CpL
    type(sparse_matrix)           :: CpR
    type(sparse_matrix)           :: H2,P
    P    = left%operators%op("P") 
    CpL = left%operators%op("Cp")
    CpR = right%operators%op("Cp")
    H2 = ts*(matmul(CpL%t(),P).x.CpR)
    H2 = H2 + H2%dgr()
    call CpL%free
    call CpR%free
    call P%free
  end function H2model

  subroutine print_mat(M,name,n)
    real(8),dimension(:,:) :: M
    character(len=*) :: name
    integer :: i,j,stride,unit,n
    stride=2**n
    open(free_unit(unit),file=str(name)//".dat")
    write(unit,*)"Matrix: "//str(name)
    do i=1,size(M,1)
       do j=1,size(M,2)
          write(unit,"(("//str(stride)//"I2))",advance="no")int(M(i,j))
          if(mod(j,stride)==0)write(unit,"(A1)",advance="no")""
       enddo
       write(unit,*)
       if(mod(i,stride)==0)write(unit,*)
    enddo
    write(unit,*)""
    close(unit)
  end subroutine print_mat






  subroutine sb_HxV(Nloc,v,Hv)
    integer                 :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                 :: i,j,jcol
    Hv=0d0
    do i=1,Nloc
       matmul: do jcol=1, spHsb%row(i)%Size
          val = spHsb%row(i)%vals(jcol)
          j   = spHsb%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end subroutine sb_HxV

end program testEDkron
