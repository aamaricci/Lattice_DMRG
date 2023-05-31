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

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"ts",finput,default=-1d0,comment="Hopping amplitude")
  call parse_input_variable(Lmax,"LMAX",finput,default=1,comment="Final chain length ")
  call parse_input_variable(m,"M",finput,default=20,&
       comment="Number of states retained at truncation of \rho")
  call parse_input_variable(Neigen,"NEIGEN",finput,default=1,&
       comment="Number of eigenpairs required to SB Hamiltonian")
  call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",finput,default=10,&
       comment="Arpack Block size parameters")
  call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",finput,default=0,&
       comment="Arpack Block size parameters")
  call parse_input_variable(lanc_niter,"LANC_NITER",finput,default=512,&
       comment="Number of Lanczos iteration in spectrum determination.")
  call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",finput,default=1d-12,&
       comment="Lanczos tolerance ")
  call parse_input_variable(lanc_threshold,"LANC_THRESHOLD",finput,default=2,&
       comment="Lapack threshold for Arpack ")
  target_Sz=0d0




  !
  !Init the single dot structure:
  dot = hubbard_site(0d0,0d0)
  !
  !Init block from single dot
  my_block=block(dot)
  call my_block%show()


  print*,"o->o"
  dimer = enlarge_block(my_block,dot,grow='left')
  Hmatrix = as_matrix(dimer%operators%op("H"))
  allocate(Evals(4**2))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""


  print*,"o--o->o"
  trimer= enlarge_block(dimer,dot,grow='left')
  Hmatrix = as_matrix(trimer%operators%op("H"))
  allocate(Evals(4**3))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""



  print*,"o--o--o->o"
  my_block= enlarge_block(trimer,dot,grow='left')
  Hmatrix = as_matrix(my_block%operators%op("H"))
  allocate(Evals(4**4))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""


  print*,"o--o--o--o->o"
  dimer= enlarge_block(my_block,dot,grow='left')
  Hmatrix = as_matrix(dimer%operators%op("H"))
  allocate(Evals(4**5))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""


  print*,"o--o--o--o--o->o"
  my_block= enlarge_block(dimer,dot,grow='left')
  spHsb = my_block%operators%op("H")
  allocate(Evals(2))
  allocate(Evecs(4**6,2))
  call sp_eigh(sb_HxV,evals,evecs,&
       10,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,2
     print*,i,Evals(i)
  enddo
  deallocate(evals,evecs)
  print*,""




  !Init block from single dot
  my_block=block(dot)

  print*,"o->o++o<-o"
  dimerL = enlarge_block(my_block,dot,grow='left')
  dimerR = enlarge_block(my_block,dot,grow='right')
  spHsb  = (dimerL%operators%op("H").x.id(dimerR%dim)) + (id(dimerL%dim).x.dimerR%operators%op("H")) + H2model(dimerL,dimerR)
  allocate(Evals(2))
  allocate(Evecs(4**4,2))
  call sp_eigh(sb_HxV,evals,evecs,&
       10,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,2
     print*,i,Evals(i)
  enddo
  deallocate(evals,evecs)
  print*,""




  print*,"o--o->o++o<-o--o"
  trimerL = enlarge_block(dimerL,dot,grow='left')
  trimerR = enlarge_block(dimerR,dot,grow='right')
  spHsb  = (trimerL%operators%op("H").x.id(trimerR%dim)) + (id(trimerL%dim).x.trimerR%operators%op("H")) + H2model(trimerL,trimerR)
  allocate(Evals(2))
  allocate(Evecs(4**6,2))
  call sp_eigh(sb_HxV,evals,evecs,&
       10,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,2
     print*,i,Evals(i)
  enddo
  deallocate(evals,evecs)
  print*,""





  stop

  !Init block from single dot
  my_block=block(dot)

  print*,"o<-o"
  dimer = enlarge_block(my_block,dot,grow='right')
  Hmatrix = as_matrix(dimer%operators%op("H"))
  allocate(Evals(4**2))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""

  print*,"o<-o--o"
  trimer= enlarge_block(dimer,dot,grow='right')
  Hmatrix = as_matrix(trimer%operators%op("H"))
  allocate(Evals(4**3))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""



  print*,"o<-o--o--o"
  my_block= enlarge_block(trimer,dot,grow='right')
  Hmatrix = as_matrix(my_block%operators%op("H"))
  allocate(Evals(4**4))
  call eigh(Hmatrix,Evals)
  do i=1,min(3,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""







contains



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
    select case(str(grow_))
    case ("left","l")
       call enl_self%put("H", (self%operators%op("H").x.id(model_d)) + (id(mblock).x.dot%operators%op("H")) + &
            H2model(self,as_block(dot)))
       call enl_self%put("Cup", Id(mblock).x.dot%operators%op("Cup"))
       call enl_self%put("Cdw", Id(mblock).x.dot%operators%op("Cdw"))
       call enl_self%put("P"  , Id(mblock).x.dot%operators%op("P"))
    case ("right","r")
       call enl_self%put("H", (id(model_d).x.self%operators%op("H")) +  (dot%operators%op("H").x.id(mblock)) + &
            H2model(as_block(dot),self))
       call enl_self%put("Cup", dot%operators%op("Cup").x.Id(mblock))
       call enl_self%put("Cdw", dot%operators%op("Cdw").x.Id(mblock))
       call enl_self%put("P"  , dot%operators%op("P").x.Id(mblock))
    end select
    !
    self_basis = self%sectors%basis()
    dot_basis  = dot%sectors%basis()
    call enl_self%set_sectors( outsum(self_basis,dot_basis) )
    !
    deallocate(self_basis,dot_basis)
  end function enlarge_block



  !H_lr = -t \sum_sigma[ (C^+_{l,sigma}@P_l) x C_{r,sigma}]  + H.c.
  function H2model(left,right) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    type(sparse_matrix)           :: CupL,CdwL
    type(sparse_matrix)           :: CupR,CdwR
    type(sparse_matrix)           :: H2,P
    P    = left%operators%op("P") 
    CupL = left%operators%op("Cup")
    CdwL = left%operators%op("Cdw")
    CupR = right%operators%op("Cup")
    CdwR = right%operators%op("Cdw")
    H2 = ts*(matmul(CupL%t(),P).x.CupR) + ts*(matmul(CdwL%t(),P).x.CdwR)  
    H2 = H2 + H2%dgr()
    call CupL%free
    call CdwL%free
    call CupR%free
    call CdwR%free
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








  function to_upper(StrIn) result(StrOut)
    character(len=*), intent(in) :: strIn
    character(len=len(strIn))    :: strOut
    integer :: i
    do i = 1,len(StrIn)
       select case(StrIn(i:i))
       case("a":"z")
          StrOut(i:i) = achar(iachar(StrIn(i:i))-32)
       case default
          StrOut(i:i) = StrIn(i:i)
       end select
    end do
  end function to_upper

  function to_lower(StrIn) result(StrOut)
    character(len=*), intent(in) :: strIn
    character(len=len(strIn))    :: strOut
    integer :: i
    do i = 1,len(StrIn)
       select case(StrIn(i:i))
       case("A":"Z")
          StrOut(i:i) = achar(iachar(StrIn(i:i))+32)
       case default
          StrOut(i:i) = StrIn(i:i)
       end select
    end do
  end function to_lower

end program testEDkron
