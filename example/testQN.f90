program testEDkron
  USE SCIFOR
  USE AUX_FUNCS
  USE MATRIX_SPARSE, id=>sp_eye
  USE MATRIX_BLOCKS
  USE TUPLE_BASIS
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
  type(block)                        :: my_block,dimer,trimer,dimerL,dimerR,trimerL,trimerR,left,right,tetramer
  type(sparse_matrix)                :: spHsb,spH
  type(site)                         :: dot
  real(8)                            :: gs_energy,target_Sz
  integer                            :: unit,m_sb
  integer                            :: model_d=4
  real(8),dimension(:,:),allocatable :: Hmatrix,Evecs
  real(8),dimension(:),allocatable   :: Evals
  integer,dimension(:),allocatable :: sb_states

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
  print*,"MY BLOCK = DOT"
  call my_block%show()
  print*,""


  print*,"o->o"
  dimer    = enlarge_block(my_block,dot,grow='left')
  print*,"DIMER"
  call dimer%show()
  Hmatrix  = as_matrix(dimer%operators%op("H"))
  call print_mat(Hmatrix,"dimer_growH",n=2)
  allocate(Evals(4**2))
  call eigh(Hmatrix,Evals)
  do i=1,min(100,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""



  print*,"o+o"
  left = block(dot)
  right= left
  !
  spHsb  = (left%operators%op("H").x.id(right%dim)) + (id(left%dim).x.right%operators%op("H")) + H2model(left,right)
  Hmatrix = as_matrix(spHsb)
  !
  call print_mat(Hmatrix,"dimer_Hsb",n=2)
  allocate(Evals(4**2))
  call eigh(Hmatrix,Evals)
  do i=1,min(16,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""

  sb_states = get_sb_states(left,right,[1,1])
  Hmatrix = as_matrix(spHsb)  
  call print_mat(Hmatrix(sb_states,sb_states),"dimer_sectorHsb",n=1)
  allocate(Evals(4),Evecs(4,4))
  Evecs = Hmatrix(sb_states,sb_states)
  call print_mat(Evecs,"Evecs",n=1)
  call eigh(Evecs,Evals)
  do i=1,min(4,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evecs,evals)



  print*,""
  print*,""
  print*,""




  print*,"o--o--o->o"
  trimer   = enlarge_block(dimer,dot,grow='left')
  tetramer = enlarge_block(trimer,dot,grow='left')
  call tetramer%show(wOP=.false.)
  Hmatrix  = as_matrix(tetramer%operators%op("H"))
  call print_mat(Hmatrix,"tetramer_H",n=4)
  allocate(Evals(4**4))
  call eigh(Hmatrix,Evals)
  do i=1,min(16,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""


  print*,"o->o++o<-o"
  my_block=block(dot)
  left  = enlarge_block(my_block,dot,grow='left')
  right = enlarge_block(my_block,dot,grow='right')
  !
  spHsb  = (left%operators%op("H").x.id(right%dim))
  spHsb  = spHsb + (id(left%dim).x.right%operators%op("H")) 
  spHsb  = spHsb + H2model(left,right)
  Hmatrix  = as_matrix(spHsb)
  call print_mat(Hmatrix,"tetramer_Hsb",n=4)
  allocate(Evals(4**4))
  call eigh(Hmatrix,Evals)
  do i=1,min(16,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evals)
  print*,""


  sb_states = get_sb_states(left,right,[2,2])
  m_sb      = size(sb_states)
  Hmatrix = as_matrix(spHsb)
  call print_mat(Hmatrix(sb_states,sb_states),"tetramer_sectorHsb",n=0)
  allocate(Evals(m_sb),Evecs(m_sb,m_sb))
  Evecs = Hmatrix(sb_states,sb_states)
  call print_mat(Evecs,"Evecs",n=1)
  call eigh(Evecs,Evals)
  do i=1,min(16,size(evals))
     print*,i,Evals(i)
  enddo
  deallocate(evecs,evals)
  print*,""









contains


  function enlarge_block(self,dot,grow) result(enl_self)
    type(block),intent(inout) :: self
    type(site)                :: dot
    character(len=*),optional :: grow
    character(len=16)         :: grow_
    type(block)               :: enl_self
    type(tbasis)              :: self_basis,dot_basis,enl_basis
    integer                   :: iqn
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    allocate(enl_self%Dims(2))
    !
    enl_self%length = self%length + 1
    enl_self%Dim    = self%Dim*dot%Dim
    enl_self%Dims   = self%Dims*dot%Dims
    !
    allocate(enl_self%sectors(size(self%sectors)))
    !
    select case(str(grow_))
    case ("left","l")
       call enl_self%put("H", &
            (self%operators%op("H").x.id(dot%dim)) +  (id(self%dim).x.dot%operators%op("H")) + &
            H2model(self,as_block(dot)))
       call enl_self%put("Cup", Id(self%dim).x.dot%operators%op("Cup"))
       call enl_self%put("Cdw", Id(self%dim).x.dot%operators%op("Cdw"))
       call enl_self%put("P"  , Id(self%dim).x.dot%operators%op("P"))
       do iqn=1,size(self%sectors)
          enl_basis = self%sectors(iqn)%basis().o.dot%sectors(iqn)%basis()
          call enl_self%set_basis( indx=iqn,  basis=enl_basis )       
       enddo
       !
    case ("right","r")
       call enl_self%put("H", &
            (id(dot%dim).x.self%operators%op("H")) +  (dot%operators%op("H").x.id(self%dim)) + &
            H2model(as_block(dot),self))
       call enl_self%put("Cup", dot%operators%op("Cup").x.Id(self%dim))
       call enl_self%put("Cdw", dot%operators%op("Cdw").x.Id(self%dim))
       call enl_self%put("P"  , dot%operators%op("P").x.Id(self%dim))
       do iqn=1,size(self%sectors)       
          enl_basis = dot%sectors(iqn)%basis().o.self%sectors(iqn)%basis()
          call enl_self%set_basis( indx=iqn,  basis=enl_basis )       
       enddo
       !
    end select
    !
  end function enlarge_block



  !H_lr = -t \sum_sigma[ (C^+_{l,sigma}@P_l) x C_{r,sigma}]  + H.c.
  function H2model(left,right) result(H2)
    type(block)                   :: left
    type(block)                   :: right
    type(sparse_matrix)           :: CupL,CdwL
    type(sparse_matrix)           :: CupR,CdwR
    type(sparse_matrix)           :: H2,P
    P   = left%operators%op("P") 
    CupL = left%operators%op("Cup")
    CdwL = left%operators%op("Cdw")
    CupR = right%operators%op("Cup")
    CdwR = right%operators%op("Cdw")
    H2 =   ts*(matmul(CupL%t(),P).x.CupR) + ts*(matmul(CdwL%t(),P).x.CdwR)  &
         - ts*(matmul(CupL,P).x.CupR%t()) - ts*(matmul(CdwL,P).x.CdwR%t())
    ! H2 = H2 + H2%dgr()
    call CupL%free
    call CdwL%free
    call CupR%free
    call CdwR%free
  end function H2model



  function get_sb_states(left,right,target_qn) result(sb_states)
    type(block)                      :: left,right
    integer,dimension(2)             :: target_qn
    integer,dimension(:),allocatable :: sb_states

    integer :: ileft,iright
    real(8),dimension(:),allocatable :: left_qn,right_qn
    integer,dimension(:),allocatable :: left_map,right_map,sb_map,states

    integer,dimension(:),allocatable :: sb_states_up,sb_states_dw
    integer                          :: m_sb
    integer                          :: m_left,m_right,Nleft,Nright
    integer                          :: iup,idw,jup,jdw
    integer                          :: i,j,istate,Ncv,im

    if(allocated(sb_states))deallocate(sb_states)
    
    do ileft=1,size(left%sectors(1))
       left_qn  = left%sectors(1)%qn(index=ileft)
       right_qn = target_qn - left_qn
       if(.not.right%sectors(1)%has_qn(right_qn))cycle
       !
       left_map = left%sectors(1)%map(qn=left_qn)
       right_map= right%sectors(1)%map(qn=right_qn)
       !
       print*,left_qn,"|",right_qn
       do i=1,size(left_map)
          do j=1,size(right_map)
             istate=right_map(j) + (left_map(i)-1)*right%Dim
             print*,left_map(i),":",right_map(j),"=",istate
             call append(sb_states, istate)
             ! call sb_left_sector%append(qn=left_qn,istate=size(sb_states))
          enddo
          print*,""
       enddo
    enddo
    !
    m_sb = size(sb_states)
    print*,sb_states
    !
  end function get_sb_states





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
          if(n>0.and.mod(j,stride)==0)write(unit,"(A1)",advance="no")""
       enddo
       write(unit,*)
       if(n>0.and.mod(i,stride)==0)write(unit,*)
    enddo
    write(unit,*)""
    close(unit)
  end subroutine print_mat




end program testEDkron





! function enlarge_block(self,dot,grow) result(enl_self)
!   type(block),intent(inout)        :: self
!   type(site),intent(inout)         :: dot
!   character(len=*),optional        :: grow
!   character(len=16)                :: grow_
!   type(block)                      :: enl_self
!   integer                          :: mblock,len,iqn
!   real(8),dimension(:),allocatable :: self_basis, dot_basis
!   !
!   grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
!   !
!   mblock =  self%dim
!   len    =  self%length
!   !
!   enl_self%length = len + 1
!   enl_self%dim    = mblock*model_d
!   !
!   select case(str(grow_))
!   case ("left","l")
!      call enl_self%put("H", (self%operators%op("H").x.id(model_d)) + (id(mblock).x.dot%operators%op("H")) + &
!           H2model(self,as_block(dot)))
!      call enl_self%put("Cup", Id(mblock).x.dot%operators%op("Cup"))
!      call enl_self%put("Cdw", Id(mblock).x.dot%operators%op("Cdw"))
!      call enl_self%put("P"  , Id(mblock).x.dot%operators%op("P"))
!   case ("right","r")
!      call enl_self%put("H", (id(model_d).x.self%operators%op("H")) +  (dot%operators%op("H").x.id(mblock)) + &
!           H2model(as_block(dot),self))
!      call enl_self%put("Cup", dot%operators%op("Cup").x.Id(mblock))
!      call enl_self%put("Cdw", dot%operators%op("Cdw").x.Id(mblock))
!      call enl_self%put("P"  , dot%operators%op("P").x.Id(mblock))
!   end select
!   !
!   allocate(enl_self%sectors(size(self%sectors)))
!   do iqn=1,size(self%sectors)       
!      self_basis = self%sectors(iqn)%basis()
!      dot_basis  = dot%sectors(iqn)%basis()
!      print*,outsum(self_basis,dot_basis)
!      call enl_self%set_sectors( indx=iqn, vec=outsum(self_basis,dot_basis) )       
!      deallocate(self_basis,dot_basis)
!   enddo
! end function enlarge_block
