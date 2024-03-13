program testEDkron
  USE SCIFOR
  USE DMRG, id=>sp_eye
  implicit none

  integer                               :: Nso
  character(len=64)                     :: finput
  integer                               :: i,unit,iorb
  character(len=1)                      :: DMRGtype
  real(8)                               :: target_Qn(2)
  real(8)                               :: ts(2),Mh(2)
  type(site)                            :: Dot
  type(sparse_matrix)                   :: C,N
  complex(8),dimension(:,:),allocatable :: Hloc,Hij

  type(block)                        :: my_block,dimer,trimer,dimerL,dimerR,trimerL,trimerR,left,right,sys
  type(sparse_matrix)                :: spHsb,spH
  real(8)                            :: gs_energy,target_Sz
  integer                            :: m_sb
  integer                            :: model_d=4
  real(8),dimension(:,:),allocatable :: Hmatrix,Evecs
  real(8),dimension(:),allocatable   :: Evals
  integer,dimension(:),allocatable      :: sb_states
  type(sectors_list)                    :: sb_sector


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=(/( -1d0,i=1,2 )/),&
       comment="Hopping amplitudes")
  call parse_input_variable(Mh,"MH",finput,default=(/(0d0,i=1,2 )/),&
       comment="Crystal field splittings")
  call parse_input_variable(target_qn,"target_qn",finput,default=[1d0,1d0],&
       comment="Target Sector Occupation [Nup,Ndw] in units [0:1]")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call read_input(finput)


  ! if(Norb/=2)stop "This code is for Norb=2. STOP"
  Nso = Nspin*Norb


  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = one*diag([Mh,Mh])
  Dot  = hubbard_site(Hloc)


  call init_dmrg(hubbard_1d_model,target_qn,ModelDot=Dot)

  !
  ! !Init block from single dot
  ! my_block=block(dot)
  ! call my_block%show()

  print*,""
  print*,""
  print*,"######################################"
  print*,"   o + o"
  print*,"######################################"
  print*,""
  print*,""



  print*,"o->o"
  dimer = block(dot)
  call enlarge_block(dimer,dot)
  Hmatrix = as_matrix(dimer%operators%op("H"))
  allocate(Evals(4**(2*Nso)))
  call eigh(Hmatrix,Evals)
  do i=1,min(4,size(evals))
     print*,i,Evals(i)/dimer%length/Norb
  enddo
  deallocate(evals)
  print*,""






  print*,"o+o: Full Diagonalization:"
  left = block(dot)
  right= left
  !
  spHsb  = (left%operators%op("H").x.id(right%dim)) + (id(left%dim).x.right%operators%op("H")) + hubbard_1d_model(left,right)
  Hmatrix = as_matrix(spHsb)
  !
  allocate(Evals(4**Nso))
  call eigh(Hmatrix,Evals)
  do i=1,size(evals)
     print*,i,Evals(i)/dimer%length/Norb
  enddo
  deallocate(evals)
  call spHsb%free()
  print*,""


  print*,"o+o: Lanc Diagonalization:"
  left = block(dot)
  right = block(dot)
  spHsb  = (left%operators%op("H").x.id(right%dim)) + (id(left%dim).x.right%operators%op("H")) +  hubbard_1d_model(left,right)
  allocate(Evals(2))
  allocate(Evecs(4**(2*Norb),2))
  call sp_eigh(sb_HxV,evals,evecs,&
       10,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,2
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""



  print*,"_o+o_ with SB states"
  left = block(dot)
  right = block(dot)
  ! call left%show()
  call get_sb_states(left,right,sb_states,sb_sector)
  ! call sb_sector%show()
  ! print*,size(sb_states)
  ! print*,sb_states
  spHsb  = sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
       sp_kron(id(left%dim),right%operators%op("H"),sb_states)      + &
       hubbard_1d_model(left,right,sb_states)

  print*,"Full ED"
  allocate(Evals(size(sb_states)))
  Hmatrix = as_matrix(spHsb)
  call eigh(Hmatrix,Evals)
  do i=1,size(evals)
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals)

  print*,"Lanc ED"
  allocate(Evals(2))
  allocate(Evecs(size(sb_states),2))
  call sp_eigh(sb_HxV,evals,evecs,&
       10,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,2
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""










  print*,""
  print*,""
  print*,"######################################"
  print*,"   o->o + o<-o"
  print*,"######################################"
  print*,""
  print*,""


  print*,"o->o++o<-o: SB States"
  left = block(dot)
  right = block(dot)
  ! print*,"Left.show: dot"
  ! call left%show()
  call enlarge_block(left,dot,grow='left')
  call enlarge_block(right,dot,grow='right')
  ! print*,"Left.show: enlarged"
  ! call left%show()
  call get_sb_states(left,right,sb_states,sb_sector)
  ! print*,"SB.show:"
  ! call sb_sector%show()
  ! print*,size(sb_states)

  spHsb  = sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
       sp_kron(id(left%dim),right%operators%op("H"),sb_states)      + &
       hubbard_1d_model(left,right,sb_states)


  allocate(Evals(2))
  allocate(Evecs(spHsb%Nrow,2))
  call sp_eigh(sb_HxV,evals,evecs,&
       2*3,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,5
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""






  print*,"o->o++o<-o: Full ED (compare)"
  left = block(dot)
  right = block(dot)
  call enlarge_block(left,dot,grow='left')
  call enlarge_block(right,dot,grow='right')
  spHsb  = (left%operators%op("H").x.id(right%dim)) + (id(left%dim).x.right%operators%op("H")) +  hubbard_1d_model(left,right)
  allocate(Evals(5))
  allocate(Evecs(4**(4*Norb),5))
  call sp_eigh(sb_HxV,evals,evecs,&
       5*3,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,10
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""

  print*,"Lanc + Filtering SB_states (2nd check):"
  call sp_reduce(spHsb,sb_states)
  allocate(Evals(2))
  allocate(Evecs(spHsb%Nrow,2))
  call sp_eigh(sb_HxV,evals,evecs,&
       2*3,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,5
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""


contains

  subroutine sp_reduce(self,list)
    type(sparse_matrix)                   :: self
    integer,dimension(:)                  :: list
    complex(8),dimension(:,:),allocatable :: M
    M = as_matrix(self)
    self = as_sparse(M(list,list))
  end subroutine sp_reduce

  !-----------------------------------------------------------------!
  ! Purpose: enlarge a given BLOCK "SELF" growing it left/right with
  ! a SITE "DOT" (specified in the init)
  !-----------------------------------------------------------------!
  subroutine enlarge_block(self,dot,grow)
    type(block)                  :: self
    type(site)                   :: dot
    character(len=*),optional    :: grow
    character(len=16)            :: grow_
    character(len=:),allocatable :: key
    type(tbasis)                 :: self_basis,dot_basis,enl_basis
    type(sparse_matrix)          :: Hb,Hd,H2
    integer                      :: i
    !
    grow_=str('left');if(present(grow))grow_=to_lower(str(grow))
    !
    if(.not.self%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: Missing self.H operator in the list"
    if(.not.dot%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: Missing dot.H operator in the list"
    !
    !> Update Hamiltonian:
    select case(str(grow_))
    case ("left","l")
       Hb = self%operators%op("H").x.id(dot%dim)
       Hd = id(self%dim).x.dot%operators%op("H")
       H2 = hubbard_1d_model(self,as_block(dot))       
    case ("right","r")
       Hb = id(dot%dim).x.self%operators%op("H")
       Hd = dot%operators%op("H").x.id(self%dim)
       H2 = hubbard_1d_model(as_block(dot),self)
    end select
    call self%put_op("H", Hb +  Hd + H2)
    !
    !> Update all the other operators in the list: 
    do i=1,size(self%operators)
       key = self%operators%key(index=i)
       if(str(key)=="H")cycle
       select case(str(grow_))
       case ("left","l")
          call self%put_op(str(key), Id(self%dim).x.dot%operators%op(str(key)))
       case ("right","r")
          call self%put_op(str(key), dot%operators%op(str(key)).x.Id(self%dim))
       end select
    enddo
    !
    !> Enlarge dimensions
    self%length = self%length + 1
    self%Dim    = self%Dim*dot%Dim
    !
    !> Enlarge the basis states
    call self%get_basis(self_basis)
    call dot%get_basis(dot_basis)
    !
    select case(str(grow_))
    case ("left","l")
       enl_basis = (self_basis.o.dot_basis)
       call self%set_basis( basis=enl_basis )
    case ("right","r")
       enl_basis = (dot_basis.o.self_basis)
       call self%set_basis( basis=enl_basis )
    end select
    !
    !Free the memory:
    call Hb%free()
    call Hd%free()
    call H2%free()
    call self_basis%free()
    call dot_basis%free()
    call enl_basis%free()
    !
  end subroutine enlarge_block

  !-----------------------------------------------------------------!
  ! Purpose: build the list of states compatible with the specified
  ! quantum numbers
  !-----------------------------------------------------------------!
  subroutine get_sb_states(sys,env,sb_states,sb_sector)
    type(block)                      :: sys,env    
    integer,dimension(:),allocatable :: sb_states
    type(sectors_list)               :: sb_sector
    integer                          :: isys,ienv
    integer                          :: i,j,istate
    real(8),dimension(:),allocatable :: sys_qn,env_qn
    integer,dimension(:),allocatable :: sys_map,env_map
    integer :: current_L
    real(8),dimension(:),allocatable      :: current_target_QN
    !
    current_L         = sys%length + env%length
    current_target_QN = int(target_qn*current_L*Norb)
    write(LOGfile,"(A22,I12)")"SuperBlock Length = ",current_L
    write(LOGfile,"(A22,"//str(size(current_target_QN))//"F12.7)")"Target_QN = ",current_target_QN
    write(LOGfile,"(A22,G12.7)")"Total         = ",sum(current_target_QN)
    write(LOGfile,"(A22,G12.7)")"Filling       = ",sum(current_target_QN)/current_L
    write(LOGfile,"(A22,G12.7)")"Filling/Norb  = ",sum(current_target_QN)/current_L/Norb

    if(allocated(sb_states))deallocate(sb_states)
    !
    call sb_sector%free()
    do isys=1,size(sys%sectors(1))
       sys_qn  = sys%sectors(1)%qn(index=isys)
       env_qn  = current_target_qn - sys_qn
       if(.not.env%sectors(1)%has_qn(env_qn))cycle
       !
       sys_map = sys%sectors(1)%map(qn=sys_qn)
       env_map = env%sectors(1)%map(qn=env_qn)
       !
       do i=1,size(sys_map)
          do j=1,size(env_map)
             istate=env_map(j) + (sys_map(i)-1)*env%Dim
             call append(sb_states, istate)
             call sb_sector%append(qn=sys_qn,istate=size(sb_states))
          enddo
       enddo
    enddo
    !
  end subroutine get_sb_states



  !H_lr = \sum_{a}h_aa*(C^+_{left,a}@P_left) x C_{right,a}] + H.c.
  function hubbard_1d_model(left,right,states) result(H2)
    type(block)                           :: left
    type(block)                           :: right
    integer,dimension(:),optional         :: states
    type(sparse_matrix),dimension(Norb,2) :: Cl,Cr
    type(sparse_matrix)                   :: P
    type(sparse_matrix)                   :: H2
    integer                               :: ispin,iorb,jorb,io,jo
    integer,dimension(2)                  :: Hdims
    character(len=:),allocatable          :: key
    !
    !>Retrieve operators:
    P = left%operators%op("P")
    do ispin=1,Nspin
       do iorb=1,Norb
          key = "C"//left%okey(iorb,ispin)
          Cl(iorb,ispin) = left%operators%op(key)
          Cr(iorb,ispin) = right%operators%op(key)
       enddo
    enddo
    !
    !> Get H2 dimensions:
    Hdims = shape(left%operators)*shape(right%operators)
    if(present(states))Hdims = [size(states),size(states)]
    !
    !>Build H2:
    call H2%init(Hdims(1),Hdims(2))
    do iorb=1,Norb
       if(present(states))then
          H2 = H2 + ts(iorb)*sp_kron(matmul(Cl(iorb,1)%dgr(),P),Cr(iorb,1),states)
          H2 = H2 + ts(iorb)*sp_kron(matmul(Cl(iorb,2)%dgr(),P),Cr(iorb,2),states)
       else
          H2 = H2 + ts(iorb)*(matmul(Cl(iorb,1)%dgr(),P).x.Cr(iorb,1))
          H2 = H2 + ts(iorb)*(matmul(Cl(iorb,2)%dgr(),P).x.Cr(iorb,2))
       endif
    enddo
    H2 = H2 + H2%dgr()
    !
    !
    !> free memory
    call P%free
    do ispin=1,Nspin
       do iorb=1,Norb
          call Cl(iorb,ispin)%free
          call Cr(iorb,ispin)%free
       enddo
    enddo
  end function hubbard_1d_model




  subroutine sb_HxV(Nloc,v,Hv)
    integer                 :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                 :: i,j,jcol
    Hv=zero
    do i=1,Nloc
       matmul: do jcol=1, spHsb%row(i)%Size
          val = dreal(spHsb%row(i)%vals(jcol))
          j   = spHsb%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end subroutine sb_HxV

end program testEDkron









  ! print*,"o--o->o"
  ! trimer = dimer
  ! call enlarge_block(trimer,dot)
  ! spHsb =  trimer%operators%op("H")
  ! print*,shape(spHsb)
  ! allocate(Evals(2))
  ! allocate(Evecs(4**(3*Norb),2))
  ! call sp_eigh(sb_HxV,evals,evecs,&
  !      10,&
  !      500,&
  !      tol=1d-12,&
  !      iverbose=.false.)
  ! do i=1,2
  !    print*,i,Evals(i)/trimer%length/Norb
  ! enddo
  ! deallocate(evals,evecs)
  ! call spHsb%free()
  ! print*,""



  ! print*,"o--o--o->o"
  ! sys = trimer
  ! print*,"go to enlarge:"
  ! call enlarge_block(sys,dot)
  ! print*,"done:"
  ! spHsb =  sys%operators%op("H")  
  ! print*,shape(spHsb),4**(4*Norb)
  ! allocate(Evals(2))
  ! allocate(Evecs(4**(4*Norb),2))
  ! call sp_eigh(sb_HxV,evals,evecs,&
  !      10,&
  !      500,&
  !      tol=1d-12,&
  !      iverbose=.false.)
  ! do i=1,2
  !    print*,i,Evals(i)/sys%length/Norb
  ! enddo
  ! deallocate(evals,evecs)
  ! print*,""
  ! print*,""
