program testEDkron
  USE SCIFOR
  USE DMRG, id=>sp_eye
  implicit none

  integer                                        :: Nso,tNso
  character(len=64)                              :: finput
  integer                                        :: i,j,unit,iorb,Nsb,k,l,r,m,arow,acol,brow,bcol,jorb,ispin,it,ip
  character(len=1)                               :: DMRGtype
  real(8)                                        :: ts(2),Mh,target_qn(2),suml,lambda,vh,aval,bval
  type(site)                                     :: Dot
  type(sparse_matrix)                            :: C,N,Cl,Cr,P
  type(sparse_matrix),allocatable,dimension(:)   :: Hleft,Hright,Hsb
  type(sparse_matrix),allocatable,dimension(:,:) :: A,B ![nz{Nspin*Norb*Norb},Nsb]
  integer,dimension(:,:),allocatable               :: bOffset,aOffset
  integer,dimension(:,:),allocatable               :: RowOffset,ColOffset
  real(8),dimension(:,:),allocatable             :: Hloc,Hij
  real(8),dimension(:),allocatable               :: sb_qn,qn,qm
  type(block)                                    :: my_block,dimer,trimer,dimerL,dimerR,trimerL,trimerR,left,right,sys
  type(sparse_matrix)                            :: spHsb,spH
  real(8)                                        :: gs_energy
  integer                                        :: m_sb,isb
  integer                                        :: model_d=4
  real(8),dimension(:,:),allocatable             :: Hmatrix,Evecs,Rho
  real(8),dimension(:),allocatable               :: Evals,Vec,Hvec,Hvec_
  integer,dimension(:),allocatable               :: sb_states,sb_map,lstates,rstates,Dls,Drs,Istates,Jstates,Offset
  type(sectors_list)                             :: sb_sector
  integer                                        :: Neigen=2
  integer                                        :: current_L
  real(8),dimension(:),allocatable               :: current_target_QN
  real(8),dimension(2),parameter                 :: qnup=[1d0,0d0],qndw=[0d0,1d0]
  real(8),dimension(2)                 :: qk
  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=(/( -1d0,i=1,2 )/),&
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
  print*,Norb

  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = 0d0
  if(Norb==2)Hloc = Mh*kron(pauli_0,pauli_z) + vh*kron(pauli_0,pauli_x)
  Dot  = hubbard_site(Hloc)
  call dot%show()

  allocate(Hij(Norb,Norb))
  Hij = diag(ts(1:Norb))
  if(Norb==2)Hij = Hij + lambda*pauli_x




  call init_dmrg(hubbard_1d_model,ModelDot=Dot)
  target_qn = DMRG_qn
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



  ! print*,"o->o"
  ! dimer = block(dot)
  ! call enlarge_block_(dimer,dot)
  ! Hmatrix = as_matrix(dimer%operators%op("H"))
  ! allocate(Evals(4**(2*Nso)))
  ! call eigh(Hmatrix,Evals)
  ! do i=1,min(4,size(evals))
  !    print*,i,Evals(i)/dimer%length/Norb
  ! enddo
  ! deallocate(evals)
  ! print*,""


  ! print*,"o+o: no QN"
  ! left = block(dot)
  ! right = block(dot)
  ! spHsb  = (left%operators%op("H").x.id(right%dim)) + (id(left%dim).x.right%operators%op("H")) +  hubbard_1d_model(left,right)
  ! allocate(Evals(Neigen))
  ! allocate(Evecs(4**(2*Norb),Neigen))
  ! call sp_eigh(sb_HxV,evals,evecs,&
  !      5*Neigen,&
  !      500,&
  !      tol=1d-12,&
  !      iverbose=.false.)
  ! print*,"Arpack Evals :"
  ! do i=1,Neigen
  !    print*,i,Evals(i)/2/left%length/Norb
  ! enddo
  ! deallocate(Evals)
  ! print*,""
  ! rho = build_density_matrix_(left%dim,right%dim,evecs(:,1),'left')
  ! allocate(Evals(size(rho,1)))
  ! call eigh(rho,Evals)
  ! print*,"\rho Lambdas:"
  ! suml=0d0
  ! do i=size(evals),1,-1
  !    suml=suml+evals(i)
  !    write(*,*)size(evals)-i+1,evals(i),suml
  !    write(300,*)size(evals)-i+1,evals(i),suml
  ! enddo
  ! deallocate(evals,evecs)
  ! print*,""




  print*,"_o+o_: with QN"
  left = block(dot)
  right = block(dot)
  call get_sb_states_(left,right,sb_states,sb_sector)

  spHsb  =  hubbard_1d_model(left,right,sb_states)  + &
       sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
       sp_kron(id(left%dim),right%operators%op("H"),sb_states)


  allocate(Evals(Neigen))
  allocate(Evecs(size(sb_states),Neigen))
  call sp_eigh(sb_HxV,evals,evecs,&
       5*Neigen,&
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





  ! print*,"o->o++o<-o: no QN"
  ! left = block(dot)
  ! right = block(dot)
  ! call enlarge_block_(left,dot,grow='left')
  ! call enlarge_block_(right,dot,grow='right')
  ! spHsb  = (left%operators%op("H").x.id(right%dim)) + (id(left%dim).x.right%operators%op("H")) +  hubbard_1d_model(left,right)
  ! allocate(Evals(Neigen))
  ! allocate(Evecs(4**(4*Norb),Neigen))
  ! call sp_eigh(sb_HxV,evals,evecs,&
  !      Neigen*3,&
  !      500,&
  !      tol=1d-12,&
  !      iverbose=.false.)
  ! print*,"Arpack Evals:"
  ! do i=1,Neigen
  !    print*,i,Evals(i)/2/left%length/Norb
  ! enddo
  ! deallocate(Evals)
  ! rho = build_density_matrix_(left%dim,right%dim,evecs(:,1),'left')
  ! allocate(Evals(size(rho,1)))
  ! print*,"\rho Lambdas:"
  ! call eigh(rho,Evals)
  ! suml=0d0
  ! do i=size(evals),1,-1
  !    suml=suml+evals(i)
  !    write(*,*)size(evals)-i+1,evals(i),suml
  !    write(400,*)size(evals)-i+1,evals(i),suml
  ! enddo
  ! deallocate(evals,evecs)
  ! print*,""








  print*,"_o->o++o<-o_: with QN"
  left  = block(dot)
  right = block(dot)
  call enlarge_block_(left,dot,grow='left')
  call enlarge_block_(right,dot,grow='right')
  call get_sb_states_(left,right,sb_states,sb_sector)


  !Old Style solution:
  spHsb  =  hubbard_1d_model(left,right,sb_states)   + &
       sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
       sp_kron(id(left%dim),right%operators%op("H"),sb_states)

  allocate(Evals(Neigen))
  allocate(Evecs(size(sb_states),Neigen))
  call sp_eigh(sb_HxV,evals,evecs,&
       3*Neigen,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,Neigen
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""





  !Build filterd H*_L, H*_R 
  Nsb = size(sb_sector)
  allocate(Hleft(Nsb),Hright(Nsb))
  allocate(Dls(Nsb),Drs(Nsb))
  do isb=1,size(sb_sector)
     sb_qn   = sb_sector%qn(index=isb)
     !
     Dls(isb)= dim_sector(left%sectors(1),sb_qn)
     Drs(isb)= dim_sector(right%sectors(1),current_target_qn - sb_qn)
     !     
     lstates= sb2lr_states(left,right,sb_states,sb_sector%map(qn=sb_qn),'left')
     rstates= sb2lr_states(left,right,sb_states,sb_sector%map(qn=sb_qn),'right')
     Hleft(isb) = sp_filter(left%operators%op("H"),lstates)
     Hright(isb)= sp_filter(right%operators%op("H"),rstates)
  enddo

  !Build Offsets
  allocate(Offset(Nsb))
  Offset=0
  do isb=2,Nsb
     Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
  enddo




  !This gives us a snapshot of the actual contributions A and B to C=A.x*.B
  !where .x*. is the restricted kron product.
  !It shows that a given set of Aij and Bkl, named A*_ij and B*_kl, contribute
  !to the determination of C*: 
  ! print*,"Cl^+@P x Cr sb_states full:"
  Cl = matmul(hconjg(left%operators%op("C_l1_s1")),left%operators%op("P"))
  Cr = right%operators%op("C_l1_s1")
  P = sp_restricted_kron(Cl,Cr,sb_states)
  call P%display()
  print*,"shape:",shape(P)
  print*,"nnz  :",P%nnz()
  print*,""
  print*,""
  print*,""
  print*,""

  print*,"P@Cl x Cr^+ sb_states full:"
  Cl = matmul(left%operators%op("P"),left%operators%op("C_l1_s1"))
  Cr = right%operators%op("C_l1_s1")
  P = sp_restricted_kron(Cl,Cr%dgr(),sb_states)
  call P%display()
  print*,"shape:",shape(P)
  print*,"nnz  :",P%nnz()
  print*,""
  print*,""
  print*,""
  print*,""


  !Workout the A*_L(\a,k), B*_R(\a,k) for each \a=1,tNso and k=1,Nsb




  tNso = 2*Nspin*count(Hij/=0d0)

  allocate(A(tNso,Nsb),B(tNso,Nsb))
  allocate(aOffset(tNso,Nsb),bOffset(tNso,Nsb))
  allocate(RowOffset(tNso,Nsb),ColOffset(tNso,Nsb))
  P = left%operators%op("P")


  it=0
  do ispin=1,2
     do iorb=1,Norb
        Cl = left%operators%op("C"//left%okey(iorb,ispin))
        do jorb=1,Norb
           if(Hij(iorb,jorb)==0d0)cycle
           Cr = right%operators%op("C"//left%okey(jorb,ispin))
           !
           qk = qnup
           if(ispin==2)qk=qndw

           ! Direct Hopping: A.x.B
           !A = Hij(iorb,jorb)*[Cl(a,s)^+.@.P] .x. B = Cr(b,s) 
           it=it+1
           do isb=1,size(sb_sector)
              qn   = sb_sector%qn(index=isb)
              qm   = qn - qk
              if(.not.sb_sector%has_qn(qm))cycle
              !Get A(\a,k)
              Istates   = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qn),'left')
              Jstates   = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qm),'left')
              A(it,isb) = Hij(iorb,jorb)*sp_filter_2(matmul(Cl%dgr(),P),Istates,Jstates)
              !Get B(\a,k)
              Istates   = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qn),'right')
              Jstates   = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qm),'right')
              B(it,isb) = sp_filter_2(Cr,Istates,Jstates)
              !
              RowOffset(it,isb)=Offset(isb)           
              ColOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
           enddo

           ! Hermitian conjugate: A.x.B
           !A = Hij(iorb,jorb)*[P.@.Cl(a,s)] .x. B = Cr(b,s)^+
           it=it+1
           do isb=1,size(sb_sector)
              qn   = sb_sector%qn(index=isb)
              qm   = qn - qk
              if(.not.sb_sector%has_qn(qm))cycle
              !Get A(\a,k)
              Istates   = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qm),'left')
              Jstates   = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qn),'left')
              A(it,isb) = Hij(iorb,jorb)*sp_filter_2(matmul(P,Cl),Istates,Jstates)
              !Get B(\a,k)
              Istates = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qm),'right')
              Jstates = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qn),'right')
              B(it,isb) = sp_filter_2(Cr%dgr(),Istates,Jstates)
              !
              RowOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
              ColOffset(it,isb)=Offset(isb)
           enddo
        enddo
     enddo
  enddo

  ! k=0
  ! do ip=1,tNso
  !    do isb=1,size(sb_sector)
  !       if(A(ip,isb)%status.AND.B(ip,isb)%status)then
  !          print*,"Index:",isb!,offset(isb),offset(isb-1)
  !          P=A(ip,isb).x.B(ip,isb)
  !          k = k + P%nnz()
  !          call P%display()
  !          do m=1,A(ip,isb)%Nrow*B(ip,isb)%Nrow
  !             arow = (m-1)/B(ip,isb)%Nrow+1
  !             brow = mod(m,B(ip,isb)%Nrow);if(brow==0)brow=B(ip,isb)%Nrow
  !             if(A(ip,isb)%row(arow)%Size==0)cycle
  !             if(B(ip,isb)%row(brow)%Size==0)cycle
  !             do i=1,A(ip,isb)%row(arow)%Size
  !                acol = A(ip,isb)%row(arow)%cols(i)
  !                aval = A(ip,isb)%row(arow)%vals(i)
  !                do j=1,B(ip,isb)%row(brow)%Size
  !                   bcol = B(ip,isb)%row(brow)%cols(j)
  !                   bval = B(ip,isb)%row(brow)%vals(j)
  !                   !print*,arow,acol," - ",brow,bcol," -- ",aval*bval
  !                   !ip is not just that isb-->isb-1 one needs to get the connected sub-sector for the column.
  !                   !for instance isb=5 (2,0) goes to isb=1 (1,0). 
  !                   ! print*,brow + (arow-1)*B(isb)%Nrow + Offset(isb), bcol+(acol-1)*B(isb)%Ncol+Offset(isb-1)
  !                   ! print*,m+Offset(isb),brow + (arow-1)*B(ip,isb)%Nrow + Offset(isb), bcol+(acol-1)*B(ip,isb)%Ncol + bOffset(ip,isb),aval*bval
  !                   print*,m+RowOffset(ip,isb),&
  !                        brow + (arow-1)*B(ip,isb)%Nrow + RowOffset(ip,isb), &
  !                        bcol + (acol-1)*B(ip,isb)%Ncol + ColOffset(ip,isb), &
  !                        aval*bval
  !                enddo
  !             enddo
  !          enddo
  !       endif
  !    enddo
  ! enddo
  ! print*,k


  ! !Continue counting:
  ! !A = Hij(iorb,jorb)*[P.@.Cl(a,s)]
  ! !B = Cr(b,s)^+
  ! do ispin=1,1
  !    do iorb=1,Norb
  !       Cl = left%operators%op("C"//left%okey(iorb,ispin))
  !       do jorb=1,Norb
  !          if(Hij(iorb,jorb)==0d0)cycle
  !          Cr = right%operators%op("C"//left%okey(jorb,ispin))
  !          !
  !          qk = qnup
  !          if(ispin==2)qk=qndw
  !          !
  !          it=it+1
  !          do isb=1,size(sb_sector)
  !             qn   = sb_sector%qn(index=isb)
  !             qm   = qn - qk
  !             if(.not.sb_sector%has_qn(qm))cycle
  !             !Get A(\a,k)
  !             Istates   = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qm),'left')
  !             Jstates   = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qn),'left')
  !             A(it,isb) = Hij(iorb,jorb)*sp_filter_2(matmul(P,Cl),Istates,Jstates)
  !             ! print*,"Size Istates[qn-1up],Jstates[qn]",size(Istates),size(Jstates)
  !             ! print*,"Istates:",Istates
  !             ! print*,"Jstates:",Jstates
  !             ! print*,'Shape(Cl^+)',shape(A(it,isb))
  !             ! call A(it,isb)%show()
  !             !Get B(\a,k)
  !             Istates = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qm),'right')
  !             Jstates = sb2lr_states(left,right,sb_states,sb_sector%map(qn=qn),'right')
  !             B(it,isb) = sp_filter_2(Cr%dgr(),Istates,Jstates)
  !             ! print*,"Size Istates[qn-1up],Jstates[qn]",size(Istates),size(Jstates)
  !             ! print*,"Istates:",Istates
  !             ! print*,"Jstates:",Jstates
  !             ! print*,'Shape(Cl^+)',shape(B(it,isb))
  !             ! call B(it,isb)%show()
  !             !
  !             ! aOffset(it,isb)=Offset(left%sectors(1)%index(qn=qm))
  !             ! bOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
  !             RowOffset(it,isb)=Offset(right%sectors(1)%index(qn=qm))
  !             ColOffset(it,isb)=Offset(isb)
  !          enddo
  !       enddo
  !    enddo
  ! enddo


  ! k=0
  ! do ip=1,tNso
  !    do isb=1,size(sb_sector)
  !       if(A(ip,isb)%status.AND.B(ip,isb)%status)then
  !          print*,"Index:",isb!,offset(isb),offset(isb-1)
  !          P=A(ip,isb).x.B(ip,isb)
  !          k = k + P%nnz()
  !          call P%display()
  !          do m=1,A(ip,isb)%Nrow*B(ip,isb)%Nrow
  !             arow = (m-1)/B(ip,isb)%Nrow+1
  !             brow = mod(m,B(ip,isb)%Nrow);if(brow==0)brow=B(ip,isb)%Nrow
  !             if(A(ip,isb)%row(arow)%Size==0)cycle
  !             if(B(ip,isb)%row(brow)%Size==0)cycle
  !             do i=1,A(ip,isb)%row(arow)%Size
  !                acol = A(ip,isb)%row(arow)%cols(i)
  !                aval = A(ip,isb)%row(arow)%vals(i)
  !                do j=1,B(ip,isb)%row(brow)%Size
  !                   bcol = B(ip,isb)%row(brow)%cols(j)
  !                   bval = B(ip,isb)%row(brow)%vals(j)
  !                   ! print*,m+bOffset(ip,isb),&
  !                   !      brow + (arow-1)*B(ip,isb)%Nrow + bOffset(ip,isb), &
  !                   !      bcol + (acol-1)*B(ip,isb)%Ncol + Offset(isb),&
  !                   !      aval*bval
  !                   print*,m+RowOffset(ip,isb),&
  !                        brow + (arow-1)*B(ip,isb)%Nrow + RowOffset(ip,isb), &
  !                        bcol + (acol-1)*B(ip,isb)%Ncol + ColOffset(ip,isb),&
  !                        aval*bval
  !                enddo
  !             enddo

  !          enddo
  !       endif
  !    enddo
  ! enddo
  ! print*,k


  ! stop





  !Old Style solution:
  spHsb  =  hubbard_1d_model(left,right,sb_states)   + &
       sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
       sp_kron(id(left%dim),right%operators%op("H"),sb_states)

  allocate(Evals(Neigen))
  allocate(Evecs(size(sb_states),Neigen))
  call sp_eigh(sb_HxV,evals,evecs,&
       3*Neigen,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,Neigen
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""


  !
  !Using _new to get energy:
  spHsb  =  hubbard_1d_model(left,right,sb_states) 
  allocate(Evals(Neigen))
  allocate(Evecs(size(sb_states),Neigen))
  call sp_eigh(sb_HxV_new,evals,evecs,&
       3*Neigen,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  do i=1,Neigen
     print*,i,Evals(i)/2/left%length/Norb
  enddo
  deallocate(evals,evecs)
  print*,""





  !test
  allocate(Vec(size(sb_states)))

  call random_number(vec)
  allocate(Hvec, source=vec)
  allocate(Hvec_, source=vec)


  print*,"Test HxV OLD"

  ! !Set spHsb to H^Lx1^R to be used in the default spHv (which only uses spHsb)
  spHsb  =  hubbard_1d_model_(left,right,sb_states) + &
       sp_kron(left%operators%op("H"),id(right%dim),sb_states) + &
       sp_kron(id(left%dim),right%operators%op("H"),sb_states)
  call sb_HxV(size(sb_states),vec,Hvec_)


  print*,"Test HxV NEW"

  !Apply H to v new style
  ! spHsb  =  hubbard_1d_model_(left,right,sb_states)
  call sb_HxV_new(size(sb_states),vec,Hvec)


  j=0
  do isb=1,size(sb_sector)
     do i=1,Dls(isb)*Drs(isb)
        print*,Hvec(i+j),Hvec_(i+j),abs(Hvec(i+j)-Hvec_(i+j))
     enddo
     print*,"  - -  - - - ",Dls(isb)*Drs(isb)
     j=j+Dls(isb)*Drs(isb)
  enddo
  !












contains



  function dim_sector(self,qn) result(dim)
    type(sectors_list)   :: self
    real(8),dimension(:) :: qn
    integer              :: dim
    dim = 0
    if(.not.self%has_qn(qn))return
    dim =  size(self%map(qn=qn))
  end function dim_sector





  subroutine sp_reduce(self,list)
    type(sparse_matrix)                   :: self
    integer,dimension(:)                  :: list
    real(8),dimension(:,:),allocatable :: M
    M = as_matrix(self)
    self = as_sparse(M(list,list))
  end subroutine sp_reduce

  !-----------------------------------------------------------------!
  ! Purpose: enlarge a given BLOCK "SELF" growing it left/right with
  ! a SITE "DOT" (specified in the init)
  !-----------------------------------------------------------------!
  subroutine enlarge_block_(self,dot,grow)
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
  end subroutine enlarge_block_

  !-----------------------------------------------------------------!
  ! Purpose: build the list of states compatible with the specified
  ! quantum numbers
  !-----------------------------------------------------------------!
  subroutine get_sb_states_(sys,env,sb_states,sb_sector)
    type(block)                      :: sys,env    
    integer,dimension(:),allocatable :: sb_states
    type(sectors_list)               :: sb_sector
    integer                          :: isys,ienv
    integer                          :: i,j,istate
    real(8),dimension(:),allocatable :: sys_qn,env_qn
    integer,dimension(:),allocatable :: sys_map,env_map

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
       print*,sys_qn,"|",env_qn
       do i=1,size(sys_map)
          do j=1,size(env_map)
             istate=env_map(j) + (sys_map(i)-1)*env%Dim
             print*,sys_map(i),env_map(j),istate
             call append(sb_states, istate)
             call sb_sector%append(qn=sys_qn,istate=size(sb_states))
          enddo
       enddo
    enddo
    !
    print*,sb_states
  end subroutine get_sb_states_





  function sb2lr_states(sys,env,sb_states,sb_map,label) result(states)
    type(block)                      :: sys,env    
    integer,dimension(:)             :: sb_states,sb_map
    character(len=*)                 :: label
    integer,dimension(:),allocatable :: tmp,states
    integer :: i,istate,l,r,isb
    !
    if(allocated(states))deallocate(states)
    !
    select case(to_lower(str(label)))
    case("left","l","sys","s")
       do i=1,size(sb_map)
          istate = sb_states(sb_map(i))
          l = (istate-1)/env%Dim+1
          call append(tmp,l)
       enddo
    case("right","r","env","e")
       do i=1,size(sb_map)
          istate = sb_states(sb_map(i))
          r = mod(istate,env%Dim);if(r==0)r=env%Dim
          call append(tmp,r)
       enddo
    end select
    allocate(states, source=uniq(tmp))
  end function sb2lr_states



  !H_lr = \sum_{a}h_aa*(C^+_{left,a}@P_left) x C_{right,a}] + H.c.
  function hubbard_1d_model(left,right,states) result(H2)
    type(block)                           :: left
    type(block)                           :: right
    integer,dimension(:),optional         :: states
    type(sparse_matrix),dimension(Norb,2) :: Cl,Cr
    type(sparse_matrix)                   :: P,A
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
    do ispin=1,Nspin
       do iorb=1,Norb
          do jorb=1,Norb
             if(Hij(iorb,jorb)==0d0)cycle
             if(present(states))then
                H2 = H2 + Hij(iorb,jorb)*sp_kron(matmul(Cl(iorb,ispin)%dgr(),P),Cr(jorb,ispin),states)
                H2 = H2 + Hij(iorb,jorb)*sp_kron(matmul(P,Cl(iorb,ispin)),Cr(jorb,ispin)%dgr(),states)
             else
                H2 = H2 + Hij(iorb,jorb)*(matmul(Cl(iorb,ispin)%dgr(),P).x.Cr(jorb,ispin))
                H2 = H2 + Hij(iorb,jorb)*(matmul(P,Cl(iorb,ispin)).x.Cr(jorb,ispin)%dgr())
             endif
          enddo
       enddo
    enddo
    ! H2 = H2 + H2%dgr()
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





  !H_lr = \sum_{a}h_aa*(C^+_{left,a}@P_left) x C_{right,a}] + H.c.
  function hubbard_1d_model_(left,right,states) result(H2)
    type(block)                           :: left
    type(block)                           :: right
    integer,dimension(:),optional         :: states
    type(sparse_matrix),dimension(Norb,2) :: Cl,Cr
    type(sparse_matrix)                   :: P,A
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
    do ispin=1,Nspin
       do iorb=1,Norb
          do jorb=1,Norb
             if(Hij(iorb,jorb)==0d0)cycle
             if(present(states))then
                H2 = H2 + Hij(iorb,jorb)*sp_kron(matmul(Cl(iorb,ispin)%dgr(),P),Cr(jorb,ispin),states)
                H2 = H2 + Hij(iorb,jorb)*sp_kron(matmul(P,Cl(iorb,ispin)),Cr(jorb,ispin)%dgr(),states)
             else
                H2 = H2 + Hij(iorb,jorb)*(matmul(Cl(iorb,ispin)%dgr(),P).x.Cr(jorb,ispin))
                H2 = H2 + Hij(iorb,jorb)*(matmul(P,Cl(iorb,ispin)).x.Cr(jorb,ispin)%dgr())
             endif
          enddo
       enddo
    enddo
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
  end function hubbard_1d_model_







  subroutine sb_HxV(Nloc,v,Hv)
    integer                 :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                 :: i,j,jcol
    Hv=zero
    do i=1,Nloc
       matmul: do jcol=1, spHsb%row(i)%Size
          val = spHsb%row(i)%vals(jcol)
          j   = spHsb%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end subroutine sb_HxV



  subroutine sb_HxV_new(Nloc,v,Hv)
    integer                 :: Nloc
    real(8),dimension(Nloc) :: v
    real(8),dimension(Nloc) :: Hv
    real(8)                 :: val
    integer                 :: i,j,jcol,k,ir,il,jr,jl,n,arow,brow,acol,bcol,ia,ib,ic,ja,jb,jc
    real(8),dimension(:,:),allocatable :: psi,Hpsi
    Hv=zero

    do k=1,size(sb_sector)

       !> apply the H^L x 1^r: need to T v and Hv
       do ir=1,Drs(k)!< fix the column: iterate the row:
          do il=1,Dls(k)
             i = ir + (il-1)*Drs(k) + offset(k)
             do jcol=1,Hleft(k)%row(il)%Size
                val = Hleft(k)%row(il)%vals(jcol)
                jl  = Hleft(k)%row(il)%cols(jcol)
                j   = ir + (jl-1)*Drs(k) + offset(k)
                Hv(i) = Hv(i) + val*v(j)
             end do
          enddo
       enddo

       !
       !> apply the 1^L x H^r
       do il=1,Dls(k)!< fix teh row: iterate the column:
          do ir=1,Drs(k)
             i = il + (ir-1)*Dls(k) + offset(k)           
             do jcol=1,Hright(k)%row(il)%Size
                val = Hright(k)%row(il)%vals(jcol)
                jl  = Hright(k)%row(il)%cols(jcol)
                j   = jl + (ir-1)*Dls(k) + offset(k)
                Hv(i) = Hv(i) + val*v(j)
             end do
          enddo
       enddo


       do it=1,tNso
          if(.not.A(it,k)%status.OR..not.B(it,k)%status)cycle
          do ic=1,A(it,k)%Nrow*B(it,k)%Nrow
             arow = (ic-1)/B(it,k)%Nrow+1
             brow = mod(ic,B(it,k)%Nrow);if(brow==0)brow=B(it,k)%Nrow
             if(A(it,k)%row(arow)%Size==0.OR.B(it,k)%row(brow)%Size==0)cycle
             i = ic + RowOffset(it,k)
             do ja=1,A(it,k)%row(arow)%Size
                acol = A(it,k)%row(arow)%cols(ja)
                aval = A(it,k)%row(arow)%vals(ja)
                do jb=1,B(it,k)%row(brow)%Size
                   bcol = B(it,k)%row(brow)%cols(jb)
                   bval = B(it,k)%row(brow)%vals(jb)
                   j = bcol+(acol-1)*B(it,k)%Ncol + ColOffset(it,k)
                   Hv(i) = Hv(i) + aval*bval*v(j)
                enddo
             enddo
          enddo
       enddo

    enddo



    ! do i=1,Nloc
    !    matmul: do jcol=1, spHsb%row(i)%Size
    !       val = spHsb%row(i)%vals(jcol)
    !       j   = spHsb%row(i)%cols(jcol)
    !       Hv(i) = Hv(i) + val*v(j)
    !    end do matmul
    ! end do


  end subroutine sb_HxV_new




  function build_density_matrix_(Nsys,Nenv,psi,direction) result(rho)
    integer                            :: Nsys,Nenv
    real(8),dimension(nsys*nenv)       :: psi
    character(len=*)                   :: direction
    real(8),dimension(:,:),allocatable :: rho
    real(8),dimension(nsys,nenv)       :: psi_tmp
    !
    if(allocated(rho))deallocate(rho)
    !
    psi_tmp = transpose(reshape(psi, [nenv,nsys]))
    !
    select case(to_lower(str(direction)))
    case ('left','l')
       allocate(rho(nsys,nsys));rho=zero
       rho  = matmul(psi_tmp,  transpose(psi_tmp)  )
    case ('right','r')
       allocate(rho(nenv,nenv));rho=zero
       rho  = matmul(transpose(psi_tmp), psi_tmp  )
    end select
  end function build_density_matrix_





  function sp_restricted_kron(A,B,states) result(AxB)
    type(sparse_matrix), intent(in) :: A,B
    integer,dimension(:),intent(in) :: states
    type(sparse_matrix)             :: AxB,Ap,Bp
    integer                         :: i,icol,j,k,kcol,l,istate,jstate
    integer                         :: indx_row,indx_col
    real(8)                         :: val,Aval,Bval
    !
    call AxB%free()
    call AxB%init(size(states),size(states))
    !
    print*,shape(A)
    print*,shape(B)
    do istate = 1,size(states)
       indx_row=states(istate)
       i = (indx_row-1)/B%Nrow+1
       k = mod(indx_row,B%Nrow);if(k==0)k=B%Nrow
       !
       do jstate=1,size(states)
          indx_col=states(jstate)
          j = (indx_col-1)/B%Ncol+1
          l = mod(indx_col,B%Ncol);if(l==0)l=B%Ncol
          if(&
               (.not.any(A%row(i)%cols==j))&
               .OR. &
               (.not.any(B%row(k)%cols==l)) )cycle
          Aval = A%get(i,j)
          Bval = B%get(k,l)
          val  = Aval*Bval
          !
          print*,indx_row," > ", indx_col,"  -  A row,col:",i,j,int(Aval)," - B row,col:",k,l,int(Bval)," < ", istate, jstate 
          !
          call append(AxB%row(istate)%vals,val)
          call append(AxB%row(istate)%cols,jstate)
          AxB%row(istate)%Size = AxB%row(istate)%Size + 1
       enddo
       !
    enddo
    print*,""
  end function sp_restricted_kron


  function sp_filter_2(A,Istates,Jstates) result(Ak)
    class(sparse_matrix), intent(in) :: A
    integer,dimension(:),intent(in)     :: Istates,Jstates
    type(sparse_matrix)                 :: Ak
    integer                             :: i,j,istate,jstate
    real(8)                             :: val
    !
    call Ak%free()
    call Ak%init(size(Istates),size(Jstates))
    !
    do istate = 1,size(Istates)
       i=Istates(istate)
       do jstate=1,size(Jstates)
          j=Jstates(jstate)
          val = A%get(i,j)
          if(val==0d0)cycle
          print*,istate,jstate,i,j,val
          call append(Ak%row(istate)%vals,val)
          call append(Ak%row(istate)%cols,jstate)
          Ak%row(istate)%Size = Ak%row(istate)%Size + 1
          !
       enddo
    enddo
  end function sp_filter_2



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
