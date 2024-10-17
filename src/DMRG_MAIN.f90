MODULE DMRG_MAIN
  USE VARS_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK
  USE DMRG_RDM
  implicit none
  private

  public :: init_DMRG
  public :: finalize_DMRG
  public :: step_DMRG
  !
  public :: infinite_DMRG
  public :: finite_DMRG
  !
  public :: write_energy
  public :: write_truncation

contains



  !##################################################################
<<<<<<< HEAD
  !              INIT DMRG ALGORITHM
=======
  !              INIT/FINALIZE DMRG ALGORITHM
>>>>>>> 7e90d6a (Updating Cmake library construction)
  !##################################################################
  subroutine init_dmrg(Hij,ModelDot)
#ifdef _CMPLX
    complex(8),dimension(:,:) :: Hij
#else
    real(8),dimension(:,:)    :: Hij
#endif
    type(site)                :: ModelDot
<<<<<<< HEAD
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: init DMRG"
#endif
    !
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
    call assert_shape(Hij,[Nspin*Norb,Nspin*Norb],"init_dmrg","Hij")
    if(allocated(HopH))deallocate(HopH)
    allocate(HopH, source=Hij)
    !
    allocate(target_qn, source=DMRG_QN)
    dot         = ModelDot
    init_left   = block(dot)
    init_right  = block(dot)
    init_called =.true.
  end subroutine init_dmrg

<<<<<<< HEAD
  
  !##################################################################
  !              FINALIZE DMRG ALGORITHM
  !##################################################################
  subroutine finalize_dmrg()
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Finalize DMRG"
#endif
    !
=======

  subroutine finalize_dmrg()
>>>>>>> 7e90d6a (Updating Cmake library construction)
    if(allocated(HopH))deallocate(HopH)    
    call dot%free()
    call init_left%free()
    call init_right%free()
    call left%free()
    call right%free()
    call psi_left%free()
    call psi_right%free()
    !
    call sb_sector%free()
    if(allocated(gs_vector))deallocate(gs_vector)
    if(allocated(sb_states))deallocate(sb_states)
    init_called  = .false.
  end subroutine finalize_dmrg



  !##################################################################
  !              FINITE/INFINITE DMRG ALGORITHM
  !##################################################################
  subroutine infinite_DMRG()
    !
<<<<<<< HEAD
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Infinite Algorithm"
#endif
    !
    if(.not.init_called)&
         stop "infinite_DMRG ERROR: DMRG not initialized. Call init_dmrg first."
    left =init_left
    right=init_right
    !
    do while (left%length < Ldmrg)
       call step_dmrg('i')
=======
    if(.not.init_called)&
         stop "infinite_DMRG ERROR: DMRG not initialized. Call init_dmrg first."
    left=init_left
    right=init_right
    !
    suffix = label_DMRG('i',1)
    !
    do while (left%length < Ldmrg)
       call step_dmrg()
       call write_energy()
       call write_truncation()
       call write_entanglement()
>>>>>>> 7e90d6a (Updating Cmake library construction)
    enddo
  end subroutine infinite_DMRG



  subroutine finite_DMRG()
    integer                                :: i,im,right_label,left_label,current_L
    type(block),dimension(:,:),allocatable :: blocks_list
    type(block)                            :: tmp
    integer                                :: j
<<<<<<< HEAD
    logical                                :: ExitSweep
    integer :: m_rleft,m_rright
    !
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Finite Algorithm"
#endif
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
    !
    if(mod(Ldmrg,2)/=0)&
         stop "finite_DMRG ERROR: Ldmrg%2 != 0. Ldmrg input must be an even number."
    !
    if(.not.init_called)&
         stop "finite_DMRG ERROR: DMRG not initialized. Call init_dmrg first."
    !
<<<<<<< HEAD
    left =init_left
    right=init_right
    !
    left_label =1 
    right_label=2
    !
    !
    !Infinite DMRG
    allocate(blocks_list(2,2*Ldmrg))
    blocks_list(left_label,1)=left
    blocks_list(right_label,1)=right
    do while (left%length < Ldmrg)
       call step_dmrg('i')
=======
    left=init_left
    right=init_right
    !
    left_label=1
    right_label=2
    !Infinite DMRG
    !
    suffix = label_DMRG('i',1)
    !
    allocate(blocks_list(2,Ldmrg))
    blocks_list(left_label,1)=left
    blocks_list(right_label,1)=right
    do while (left%length <= Ldmrg)
       call step_dmrg()
       call write_energy()
       call write_entanglement()
>>>>>>> 7e90d6a (Updating Cmake library construction)
       blocks_list(left_label,left%length)=left
       blocks_list(right_label,right%length)=right
    enddo
    print*,""
    print*,""
    print*,"START FINITE DMRG:"
    print*,""
    print*,""
    !
    !Finite DMRG: start sweep forth&back:
    do im=1,size(Msweep)
       if(Esweep(im)==0)then
          write(*,"(A,I3,I6)")"Sweep, M:",im,Msweep(im)
       else
          write(*,"(A,I3,F8.4)")"Sweep, E:",im,Esweep(im)
       endif
<<<<<<< HEAD
       !
       ExitSweep=.false.
       sweep: do while(.true.)
          right = blocks_list(right_label,2*Ldmrg - left%length)
=======
       suffix = label_DMRG('f',im)
       sweep: do while(.true.)
          right = blocks_list(right_label,Ldmrg - left%length)
>>>>>>> 7e90d6a (Updating Cmake library construction)
          if(right%length==1)then
             right_label= 3-right_label
             left_label = 3-left_label
             tmp        = left
             left       = right
             right      = tmp
             call tmp%free()
<<<<<<< HEAD
             ExitSweep  = .true.
          endif
          !
          call step_dmrg('f',left_label,im)
=======
          endif
          !
          call step_dmrg(left_label,im)
          call write_energy()
          call write_entanglement()
>>>>>>> 7e90d6a (Updating Cmake library construction)
          !
          blocks_list(left_label,left%length) = left
          print*,""
          print*,""
<<<<<<< HEAD
          if(ExitSweep.AND.left_label==1.AND.left%length==Ldmrg+1)exit sweep
=======
          if(left_label==1.AND.left%length==Ldmrg/2)exit sweep
>>>>>>> 7e90d6a (Updating Cmake library construction)
       enddo sweep
    enddo
    !
  end subroutine finite_DMRG












  !##################################################################
  !              WORKHORSE: STEP DMRG 
  !##################################################################
<<<<<<< HEAD
  subroutine step_dmrg(type,label,sweep)
    character(len=1) :: type
    integer,optional :: label,sweep
    integer          :: iLabel
    integer          :: m_sb
    integer          :: m_rleft,m_rright
    integer          :: m_es,m_ee
    integer          :: m_left,m_right
    integer          :: m_eleft,m_eright
    integer          :: current_L
    integer          :: Lleft,Lright
    logical          :: renormalize
    !
    !just 4 DMRG_graphic
    iLabel=0;if(present(label))iLabel=label
    !
    select case(to_lower(type))
    case('f')
       if(.not.present(sweep))stop "step_dmrg ERROR: type=F but no sweep-index provided"
       if(sweep>Nsweep.OR.sweep<1)stop "step_dmrg ERROR: isweep < 1 OR isweep > Nsweep"
       Mstates = Msweep(sweep)
       Estates = Esweep(sweep)
       suffix  = label_DMRG(type,sweep)
    case('i')
       Mstates=Mdmrg
       Estates=Edmrg
       suffix = label_DMRG(type)
    case default
       stop "step_dmrg ERROR: unsupported type. Pick Finite or Infinite"
    end select
    !
=======
  subroutine step_dmrg(label,isweep)
    integer,optional                   :: label,isweep
    integer                            :: iLabel
    integer                            :: m_sb
    integer                            :: m_rleft,m_rright
    integer                            :: m_es,m_ee
    integer                            :: m_left,m_right
    integer                            :: m_eleft,m_eright
    integer                            :: current_L
    !
    iLabel=0;if(present(label))iLabel=label
    !2Bremoved
    suffix = label_DMRG('i',1)
    !
    Mstates=Mdmrg
    Estates=Edmrg
    if(present(isweep))then
       if(isweep>Nsweep)stop "step_dmrg ERROR: isweep > Nsweep"
       if(isweep<1)stop "step_dmrg ERROR: isweep < 1"
       Mstates = Msweep(isweep)
       Estates = Esweep(isweep)
    endif
>>>>>>> 7e90d6a (Updating Cmake library construction)
    !
    if(.not.left%is_valid(.true.))stop "single_dmrg_step error: left is not a valid block"
    if(.not.right%is_valid(.true.))stop "single_dmrg_step error: right is not a valid block"
    !
    !> dimension of the incoming L/R BLOCKS
    m_left  = left%dim
    m_right = right%dim
    !
    !> START DMRG STEP:
    call dmrg_graphic(iLabel)    
    call start_timer()
    !
<<<<<<< HEAD
    ! !#################################
    ! !    Renormalize BLOCKS
    ! !#################################
    ! call renormalize_block('left',m_rleft)
    ! call renormalize_block('right',m_rright)
=======
    !#################################
    !    Renormalize BLOCKS
    !#################################
    call renormalize_block('left',m_rleft)
    call renormalize_block('right',m_rright)
>>>>>>> 7e90d6a (Updating Cmake library construction)
    !
    !
    !#################################
    !    Enlarge L/R BLOCKS: +1 DOT
    !#################################
<<<<<<< HEAD
    Lleft =left%length
    Lright=right%length
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
    call enlarge_block(left,dot,grow='left')
    call enlarge_block(right,dot,grow='right')
    !
    !
    !#################################
    !    Build SUPER-BLOCK Sector
    !#################################
    m_eleft           = left%dim
    m_eright          = right%dim
    current_L         = left%length + right%length
    current_target_QN = int(target_qn*current_L*Norb)
    !
    call sb_get_states()
    m_sb = size(sb_states)
    !
<<<<<<< HEAD
    !
    !#################################
    !       DIAG SUPER-BLOCK
    !#################################
    call sb_diag()
    !
    !#################################
    !      BUILD RDM
    !#################################
    call sb_get_rdm()
    !
    !#################################
    !    Renormalize BLOCKS
    !#################################
    !What are the conditions to skip Renormalization:
    ! 1. Infinite algorithm:
    !   - only at the last loop which is when *enlarged* blocks have length=N, i.e. L.len=N
    !      if(to_lower(type)=='i'.AND.left%length==Ldmrg)
    ! 2. Finite algorithm:
    !   - sweep is made of
    !      [label=1]: L(N/2->N-1)>R,
    !      [lable=2]: R(1->N)->L,
    !      [label=1]: L(1->N/2)->R
    !     so if label=1 and L.len=N:
    !      if(left_label==1.AND.left%length==Ldmrg)
    renormalize=.not.((iLabel==1.OR.to_lower(type)=='i').AND.Lleft==Ldmrg)
    if(renormalize)then
       call renormalize_block('left',m_rleft)
       call renormalize_block('right',m_rright)
    endif
    !
    !> STOP DMRG STEP:
    call stop_timer("dmrg_step")
    !
    !#################################
    !      WRITE AND EXIT
    !#################################
    write(LOGfile,"(A,I12,12X,I12)")&
         "         Blocks Length               :",Lleft,Lright
=======
    write(LOGfile,"(A,I12,12X,I12)")&
         "Blocks Length                        :",left%length-1,right%length-1
    write(LOGfile,"(A,I12,12X,I12,3X,A3,I6,4X,I6,A1)")&
         "Renormalized Blocks Dim              :",m_rleft,m_rright,"< (",m_left,m_right,")"
    write(LOGfile,"(A,L12,12X,L12)")&
         "Truncating                           :",Mstates<=m_left,Mstates<=m_right
    write(LOGfile,"(A,2ES24.15)")&
         "Truncation Errors                    :",truncation_error_left,truncation_error_right
>>>>>>> 7e90d6a (Updating Cmake library construction)
    write(LOGfile,"(A,I12,12X,I12)")&
         "Enlarged Blocks Length               :",left%length,right%length
    write(LOGfile,"(A,I12,12X,I12)")&
         "Enlarged Blocks Dim                  :",m_eleft,m_eright  
    write(LOGfile,"(A,I12)")&
         "SuperBlock Length                    :",current_L
    write(LOGfile,"(A,I12,A2,I12,A1,F10.5,A1)")&
         "SuperBlock Dimension  (tot)          :", &
         m_sb," (",m_eleft*m_eright,")",100*dble(m_sb)/m_eleft/m_eright,"%"
    write(LOGfile,"(A,"//str(size(current_target_QN))//"F24.15)")&
         "Target_QN                            :",current_target_QN
    write(LOGfile,"(A,3x,G24.15)")&
         "Total                                :",sum(current_target_QN)
    write(LOGfile,"(A,3x,G24.15)")&
         "Filling                              :",sum(current_target_QN)/current_L
    write(LOGfile,"(A,3x,G24.15)")&
         "Filling/Norb                         :",sum(current_target_QN)/current_L/Norb
<<<<<<< HEAD
    if(renormalize)then
       write(LOGfile,"(A,I12,12X,I12,3X,A3,I6,4X,I6,A1)")&
            "Renormalized Blocks Dim              :",m_rleft,m_rright,"< (",m_left,m_right,")"
       write(LOGfile,"(A,L12,12X,L12)")&
            "Truncating                           :",Mstates<=m_left,Mstates<=m_right
       write(LOGfile,"(A,2ES24.15)")&
            "Truncation Errors                    :",truncation_error_left,truncation_error_right
    endif
    !
=======
    !
    !
    !#################################
    !       DIAG SUPER-BLOCK
    !#################################
    call sb_diag()
    !
    !#################################
    !      BUILD RDM
    !#################################
    call sb_get_rdm()
    !
    !
    !#################################
    !      WRITE AND EXIT
    !#################################
>>>>>>> 7e90d6a (Updating Cmake library construction)
    select case(left%type())
    case ("fermion","f")
       write(LOGfile,"(A,"//str(Lanc_Neigen)//"F24.15)")&
            "Energies/N                           :",gs_energy/sum(current_target_QN)
    case ("spin","s")
       write(LOGfile,"(A,"//str(Lanc_Neigen)//"F24.15)")&
            "Energies/L                           :",gs_energy/current_L
    end select
<<<<<<< HEAD
    !
    call write_energy()
    call write_truncation()
    call write_entanglement()
    !
=======
    call stop_timer("dmrg_step")
>>>>>>> 7e90d6a (Updating Cmake library construction)
    !
    !Clean memory:
    call spHsb%free()
    !
    return
  end subroutine step_dmrg





















  !##################################################################
  !                 WRITE TO FILE
  !##################################################################
  !-----------------------------------------------------------------!
  ! Purpose: write energy to file
  !-----------------------------------------------------------------!
  subroutine write_energy()
    integer                   :: current_L
    integer                   :: Eunit
    current_L = left%length + right%length
    Eunit     = fopen("energyVSleft.length_"//str(suffix),append=.true.)
    write(Eunit,*)left%length,gs_energy/current_L/Norb
    close(Eunit)
  end subroutine write_energy


  !-----------------------------------------------------------------!
  ! Purpose: write entanglement entropy to file
  !-----------------------------------------------------------------!
  subroutine write_truncation()
    integer                   :: current_L
    integer                   :: Eunit
    current_L = left%length + right%length
    Eunit     = fopen("truncationVSleft.length_"//str(suffix),append=.true.)
    write(Eunit,*)left%length,&
         truncation_error_left/current_L/Norb,&
         truncation_error_right/current_L/Norb
    close(Eunit)
  end subroutine write_truncation


  !-----------------------------------------------------------------!
  ! Purpose: write entanglement entropy to file
  !-----------------------------------------------------------------!
  subroutine write_entanglement()
    integer :: Eunit,i
    real(8) :: entropy
    !
    entropy=0d0
    if(.not.allocated(rho_left_evals))return
    do i=1,size(rho_left_evals)
       if(rho_left_evals(i)<=0d0)cycle
       entropy = entropy-rho_left_evals(i)*log(rho_left_evals(i))
    enddo
    Eunit     = fopen("SentropyVSleft.length_"//str(suffix),append=.true.)
    write(Eunit,*)left%length,entropy
    close(Eunit)
  end subroutine write_entanglement

END MODULE DMRG_MAIN












! !##################################################################
! !              WORKHORSE: STEP DMRG 
! !##################################################################
! subroutine step_dmrg(label,isweep)
!   integer,optional                   :: label,isweep
!   integer                            :: iLabel
!   integer                            :: Mstates,j_
!   real(8)                            :: Estates,e_,err
!   integer                            :: m_sb,m_s,m_e
!   integer                            :: m_left,m_right,Nleft,Nright
!   integer                            :: isb,m_err(1),m_threshold
!   integer                            :: i,j,r,l,iqn,im,unit,current_L,istate,pos
!   integer,dimension(:),allocatable   :: left_map,right_map,sb_map
!   real(8),dimension(:),allocatable   :: sb_qn,qn
!   type(tbasis)                       :: left_basis,right_basis
!   type(sparse_matrix)                :: trRho_left,trRho_right


!   real(8),dimension(:),allocatable   :: El,Er
!   real(8),dimension(:,:),allocatable   :: Ml,Mr

!   !
!   iLabel=0;if(present(label))iLabel=label
!   Mstates=Mdmrg
!   Estates=Edmrg
!   if(present(isweep))then
!      if(isweep>Nsweep)stop "step_dmrg ERROR: isweep > Nsweep"
!      if(isweep<1)stop "step_dmrg ERROR: isweep < 1"
!      Mstates = Msweep(isweep)
!      Estates = Esweep(isweep)
!   endif
!   !
!   !Start DMRG step timer
!   call start_timer()
!   !
!   !
!   !Check if blocks are valid ones
!   if(.not.left%is_valid(.true.))stop "single_dmrg_step error: left is not a valid block"
!   if(.not.right%is_valid(.true.))stop "single_dmrg_step error: right is not a valid block"
!   !
!   !Enlarge the Blocks:
!   call start_timer()
!   write(LOGfile,"(A22,2I12)")"Blocks Length (L-R) = ",left%length,right%length
!   call enlarge_block(left,dot,grow='left')
!   call enlarge_block(right,dot,grow='right')
!   call stop_timer("Enlarge blocks")
!   !
!   print*,"Enl LEFT:"
!   call left%show(file="EnlLEFT_"//str(left%length)//".dat")
!   print*,"Enl RIGHT:"
!   call right%show(file="EnlRIGHT_"//str(right%length)//".dat")

!   m_left  = left%dim
!   m_right = right%dim
!   if(.not.left%is_valid())stop "dmrg_step error: enlarged_left is not a valid block"
!   if(.not.right%is_valid())stop "dmrg_step error: enlarged_right is not a valid block"


!   !#################################
!   !    Build SUPER-BLOCK Sector
!   !#################################
!   current_L         = left%length + right%length
!   current_target_QN = int(target_qn*current_L*Norb)
!   write(LOGfile,"(A22,I12)")"SuperBlock Length = ",current_L
!   write(LOGfile,"(A22,"//str(size(current_target_QN))//"F12.7)")"Target_QN = ",current_target_QN
!   write(LOGfile,"(A22,G12.7)")"Total         = ",sum(current_target_QN)
!   write(LOGfile,"(A22,G12.7)")"Filling       = ",sum(current_target_QN)/current_L
!   write(LOGfile,"(A22,G12.7)")"Filling/Norb  = ",sum(current_target_QN)/current_L/Norb
!   !
!   call sb_get_states()
!   !
!   m_sb = size(sb_states)
!   write(LOGfile,"(A,I12,I12))")&
!        "Enlarged Blocks dimensions           :",m_left,m_right  
!   write(LOGfile,"(A,I12,A1,I12,A1,F10.5,A1)")&
!        "SuperBlock Dimension  (tot)          :", &
!        m_sb,"(",m_left*m_right,")",100*dble(m_sb)/m_left/m_right,"%"



!   !#################################
!   !       DIAG SUPER-BLOCK
!   !#################################
!   call sb_diag()





!   !#################################
!   !      BUILD RHO and PSI
!   !#################################
!   call start_timer()
!   call rho_left%free()
!   call rho_right%free()
!   call psi_left%free()
!   call psi_right%free()
!   do isb=1,size(sb_sector)
!      sb_qn  = sb_sector%qn(index=isb) !this is really the Left qn
!      sb_map = sb_sector%map(index=isb)
!      Nleft  = size(left%sectors(1)%map(qn=sb_qn))
!      Nright = size(right%sectors(1)%map(qn=(current_target_qn - sb_qn)))
!      if(Nleft*Nright==0)cycle
!      !
!      print*,sb_qn,current_target_qn - sb_qn,Nleft,Nright,Nleft*Nright,size(sb_map)
!      qn = sb_qn
!      call rho_left%append(&
!           build_density_matrix(Nleft,Nright,gs_vector(:,1),sb_map,'left'),&
!           qn=qn,map=left%sectors(1)%map(qn=qn))
!      call psi_left%append(&
!           build_PsiMat(Nleft,Nright,gs_vector(:,1),sb_map,'left'),&
!           qn=qn,map=left%sectors(1)%map(qn=qn))
!      !
!      qn = current_target_qn-sb_qn
!      call rho_right%append(&
!           build_density_matrix(Nleft,Nright,gs_vector(:,1),sb_map,'right'),&
!           qn=qn,map=right%sectors(1)%map(qn=qn))
!      call psi_right%append(&
!           build_PsiMat(Nleft,Nright,gs_vector(:,1),sb_map,'right'),&
!           qn=qn,map=right%sectors(1)%map(qn=qn))
!   enddo
!   !
!   print*,"rho_Left:"
!   call rho_left%show(file="Rho_L_"//str(left%length)//".dat")

!   print*,"rho_Right"
!   call rho_right%show(file="Rho_R_"//str(right%length)//".dat")
!   !


!   !Build Truncated Density Matrices:
!   call rho_left%eigh(sort=.true.,reverse=.true.)
!   call rho_right%eigh(sort=.true.,reverse=.true.)
!   rho_left_evals  = rho_left%evals()
!   rho_right_evals = rho_right%evals()
!   !
!   if(Mstates/=0)then
!      m_s   = min(Mstates,m_left,size(rho_left_evals))
!      m_e   = min(Mstates,m_right,size(rho_right_evals))
!   elseif(Estates/=0d0)then
!      m_err = minloc(abs(1d0-cumulate(rho_left_evals)-Estates))
!      m_s   = m_err(1)
!      m_err = minloc(abs(1d0-cumulate(rho_right_evals)-Estates))
!      m_e   = m_err(1)
!   else
!      stop "Mdmrg=Edmrg=0 can not fix a threshold for the RDM"
!   endif
!   !
!   e_=rho_left_evals(m_s)
!   j_=m_s
!   do i=j_+1,size(rho_left_evals)
!      err=abs(rho_left_evals(i)-e_)/e_
!      if(err<=1d-1)m_s=m_s+1
!   enddo
!   e_=rho_right_evals(m_e)
!   j_=m_e
!   do i=j_+1,size(rho_right_evals)
!      err=abs(rho_right_evals(i)-e_)/e_
!      if(err<=1d-1)m_e=m_e+1
!   enddo
!   !
!   !
!   truncation_error_left  = 1d0 - sum(rho_left_evals(1:m_s))
!   truncation_error_right = 1d0 - sum(rho_right_evals(1:m_e))
!   !
!   !
!   print*,"Left_RhoEig"
!   call rho_left%show(file="EighRho_L_"//str(left%length)//".dat")
!   print*,"RIght_RhoEig"
!   call rho_right%show(file="EighRho_R_"//str(right%length)//".dat")

!   !>truncation-rotation matrices:
!   trRho_left  = rho_left%sparse(m_left,m_s)
!   trRho_right = rho_right%sparse(m_right,m_e)

!   print*,"Left_RhoTr"
!   call trRho_left%show(file="TrRho_L_"//str(left%length)//".dat")
!   print*,"Right_RhoTr"
!   call trRho_right%show(file="TrRho_R_"//str(right%length)//".dat")
!   !
!   !>Store all the rotation/truncation matrices:
!   call left%put_omat(str(left%length),trRho_left)
!   call right%put_omat(str(right%length),trRho_right)
!   !
!   !>Renormalize Blocks:
!   call left%renormalize(as_matrix(trRho_left))
!   call right%renormalize(as_matrix(trRho_right))
!   !
!   !>Prepare output and update basis state
!   do im=1,m_s
!      call left_basis%append( qn=rho_left%qn(m=im) )
!   enddo
!   do im=1,m_e
!      call right_basis%append( qn=rho_right%qn(m=im) )
!   enddo
!   call left%set_basis(basis=left_basis)
!   call right%set_basis(basis=right_basis)


!   print*,"Left_NewBasis"
!   call rho_left%show(file="NewBasis_L_"//str(left%length)//".dat")
!   print*,"RIght_NewBasis"
!   call rho_right%show(file="NewBasis_R_"//str(right%length)//".dat")

!   call stop_timer("Get Rho, U->Renormalize, Setup New Basis")
!   !
!   !
!   !
!   write(LOGfile,"(A,L12,12X,L12)")&
!        "Truncating                           :",Mstates<=m_left,Mstates<=m_right
!   write(LOGfile,"(A,I12,12X,I12)")&
!        "Truncation Dim                       :",m_s,m_e
!   write(LOGfile,"(A,2ES24.15)")&
!        "Truncation Errors                    :",truncation_error_left,truncation_error_right
!   select case(left%type())
!   case ("fermion","f")
!      write(LOGfile,"(A,"//str(Lanc_Neigen)//"F24.15)")&
!           "Energies/N                           :",gs_energy/sum(current_target_QN)
!   case ("spin","s")
!      write(LOGfile,"(A,"//str(Lanc_Neigen)//"F24.15)")&
!           "Energies/L                           :",gs_energy/current_L
!   end select
!   call stop_timer("dmrg_step")
!   call dmrg_graphic(iLabel)
!   !
!   !Clean memory:
!   if(allocated(left_map))deallocate(left_map)
!   if(allocated(right_map))deallocate(right_map)
!   if(allocated(sb_map))deallocate(sb_map)
!   if(allocated(sb_qn))deallocate(sb_qn)
!   if(allocated(qn))deallocate(qn)
!   call rho_left%free()
!   call rho_right%free()
!   call left_basis%free()
!   call right_basis%free()
!   call trRho_left%free()
!   call trRho_right%free()
!   call spHsb%free()
!   !
!   return
! end subroutine step_dmrg






! !##################################################################
! !              RESHAPE GROUND STATE AS MATRIX 
! !##################################################################
! function build_PsiMat(Nleft,Nright,psi,map,direction) result(psi_mat)
!   integer                            :: Nleft,Nright
!   real(8),dimension(:)               :: psi
!   integer,dimension(nleft*nright)       :: map
!   character(len=*)                   :: direction
!   real(8),dimension(:,:),allocatable :: psi_mat
!   if(allocated(psi_mat))deallocate(psi_mat)
!   select case(to_lower(str(direction)))
!   case ('left','l')
!      allocate(psi_mat(nleft,nright));psi_mat=zero
!      psi_mat = transpose(reshape(psi(map), [nright,nleft]))
!   case ('right','r')
!      allocate(psi_mat(nright,nleft));psi_mat=zero
!      psi_mat = reshape(psi(map), [nright,nleft])
!   end select
! end function build_psimat





!     call start_timer("Get Rho")
!     call rho_left%free()
!     call rho_right%free()
!     do isb=1,size(sb_sector)
!        sb_qn  = sb_sector%qn(index=isb)
!        sb_map = sb_sector%map(index=isb)
!        Nleft   = size(left%sectors(1)%map(qn=sb_qn))
!        Nright  = size(right%sectors(1)%map(qn=(current_target_qn - sb_qn)))
!        if(Nleft*Nright==0)cycle
!        !
!        qn = sb_qn
!        call rho_left%append(&
!             build_density_matrix(Nleft,Nright,gs_vector(:,1),sb_map,'left'),&
!             qn=qn,map=left%sectors(1)%map(qn=qn))
!        !
!        qn = current_target_qn-sb_qn
!        call rho_right%append(&
!             build_density_matrix(Nleft,Nright,gs_vector(:,1),sb_map,'right'),&
!             qn=qn,map=right%sectors(1)%map(qn=qn))
!     enddo
!     call stop_timer()
! #ifdef _DEBUG
!     write(LOGfile,*)"Show Rho_L/R:"
!     call rho_left%show(file="Rho_L_"//str(left%length)//".dat")
!     call rho_right%show(file="Rho_R_"//str(right%length)//".dat")
! #endif






! !##################################################################
! !              BUILD REDUCED DENSITY MATRIX 
! !##################################################################
! function build_density_matrix(Nleft,Nright,psi,map,direction) result(rho)
!   integer                            :: Nleft,Nright
!   real(8),dimension(:)               :: psi
!   integer,dimension(nleft*nright)    :: map
!   character(len=*)                   :: direction
!   real(8),dimension(:,:),allocatable :: rho
!   real(8),dimension(nleft,nright)    :: psi_tmp
!   integer                            :: il,ir,i
!   !
!   if(allocated(rho))deallocate(rho)
!   !
!   !psi_tmp = transpose(reshape(psi(map), [nright,nleft]))
!   do concurrent(il=1:Nleft,ir=1:Nright)
!      i = map(ir + (il-1)*Nright)
!      psi_tmp(il,ir) = psi(i)
!   enddo
!   !
!   select case(to_lower(str(direction)))
!   case ('left','l')
!      allocate(rho(nleft,nleft));rho=zero
!      rho  = matmul(psi_tmp,  (transpose(psi_tmp)) )
!   case ('right','r')
!      allocate(rho(nright,nright));rho=zero
!      rho  = matmul((transpose(psi_tmp)), psi_tmp  )
!   end select
! end function build_density_matrix








!        call start_timer("Diag Rho")
!        call rho_left%eigh(sort=.true.,reverse=.true.)
!        call rho_right%eigh(sort=.true.,reverse=.true.)
!        rho_left_evals  = rho_left%evals()
!        rho_right_evals = rho_right%evals()
!        call stop_timer()
!        !
!        call start_timer("Renormalize Blocks + Setup Basis")
!        !Build Truncated Density Matrices:
!        if(Mstates/=0)then
!           m_s = min(Mstates,m_left,size(rho_left_evals))
!           m_e = min(Mstates,m_right,size(rho_right_evals))       
!        elseif(Estates/=0d0)then
!           m_err = minloc(abs(1d0-cumulate(rho_left_evals)-Estates))
!           m_s   = m_err(1)
!           m_err = minloc(abs(1d0-cumulate(rho_right_evals)-Estates))
!           m_e   = m_err(1)
!        else
!           stop "Mdmrg=Edmrg=0 can not fix a threshold for the RDM"
!        endif
!        !
!        e_=rho_left_evals(m_s)
!        j_=m_s
!        do i=j_+1,size(rho_left_evals)
!           err=abs(rho_left_evals(i)-e_)/e_
!           if(err<=1d-1)m_s=m_s+1
!        enddo
!        e_=rho_right_evals(m_e)
!        j_=m_e
!        do i=j_+1,size(rho_right_evals)          
!           err=abs(rho_right_evals(i)-e_)/e_
!           if(err<=1d-1)m_e=m_e+1
!        enddo
!        !
!        truncation_error_left  = 1d0 - sum(rho_left_evals(1:m_s))
!        truncation_error_right = 1d0 - sum(rho_right_evals(1:m_e))
!        !
!        !>truncation-rotation matrices:
!        trRho_left  = rho_left%sparse(m_left,m_s)
!        trRho_right = rho_right%sparse(m_right,m_e)
!        !
!        !
!        !>Store all the rotation/truncation matrices:
!        call left%put_omat(str(left%length),trRho_left)
!        call right%put_omat(str(right%length),trRho_right)
!        !
!        !>Renormalize Blocks:
!        call left%renormalize(as_matrix(trRho_left))
!        call right%renormalize(as_matrix(trRho_right))
!        !
!        !>Prepare output and update basis state
!        do im=1,m_s
!           call left_basis%append( qn=rho_left%qn(m=im) )
!        enddo
!        do im=1,m_e
!           call right_basis%append( qn=rho_right%qn(m=im) )
!        enddo
!        call left%set_basis(basis=left_basis)
!        call right%set_basis(basis=right_basis)
!        !
! #ifdef _DEBUG
!        unit     = fopen("lambdas_L_"//str(left%length)//".dat")       
!        do i=1,m_s
!           write(unit,*)i,rho_left_evals(i),floor(log10(abs(rho_left_evals(i))))
!        enddo
!        write(unit,*)" "
!        do i=m_s+1,size(rho_left_evals)
!           write(unit,*)i,rho_left_evals(i),floor(log10(abs(rho_left_evals(i))))
!        enddo
!        close(unit)
!        !
!        unit     = fopen("lambdas_R_"//str(left%length)//".dat")       
!        do i=1,m_e
!           write(unit,*)i,rho_right_evals(i),floor(log10(abs(rho_right_evals(i))))
!        enddo
!        write(unit,*)" "
!        do i=m_s+1,size(rho_right_evals)
!           write(unit,*)i,rho_right_evals(i),floor(log10(abs(rho_right_evals(i))))
!        enddo
!        close(unit)
!        !
!        call trRho_left%show(file="TrRho_L_"//str(left%length)//".dat")
!        call trRho_right%show(file="TrRho_R_"//str(right%length)//".dat")
!        !
!        call rho_left%show(file="NewBasis_L_"//str(left%length)//".dat")
!        call rho_right%show(file="NewBasis_R_"//str(right%length)//".dat")
! #endif
!        !
!        call stop_timer()
!        write(LOGfile,"(A,2ES24.15)")"Truncation Errors                    :",&
!             truncation_error_left,truncation_error_right
!        call left_basis%free()
!        call right_basis%free()
!        call trRho_left%free()
!        call trRho_right%free()
!        call rho_left%free()
!        call rho_right%free()
