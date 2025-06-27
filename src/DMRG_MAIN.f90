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


  type(block),dimension(:,:),allocatable :: blocks_list

contains



  !##################################################################
  !              INIT DMRG ALGORITHM
  !##################################################################
  subroutine init_dmrg(Hij,ModelDot)
#ifdef _CMPLX
    complex(8),dimension(:,:)   :: Hij
#else
    real(8),dimension(:,:)      :: Hij
#endif
    type(site),dimension(:)     :: ModelDot
    integer                     :: ilat
    !
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: init DMRG"
#endif
    !
    !SETUP the Hopping Hamiltonian term:
    call assert_shape(Hij,[Nspin*Norb,Nspin*Norb],"Init_DMRG","Hij")
    if(allocated(HopH))deallocate(HopH)
    allocate(HopH, source=Hij)
    !
    !SETUP the Dots
    if(size(ModelDot)/=2*Ldmrg+2.OR.size(ModelDot)/=1)stop "Init_DMRG ERROR: size(ModelDot) != 2*Ldmrg+2 or 1"
    allocate(dots(2*Ldmrg+2))
    do ilat=1,size(dots)
       select case(size(ModelDot))
       case (1)    ;dots(ilat) = ModelDot(1)
       case default;dots(ilat) = ModelDot(ilat)
       end select
    enddo
    !
    !SETUP the initial DMRG structure
    allocate(target_qn, source=DMRG_QN)
    init_left   = block(dots(1),BlockTag='left')
    init_right  = block(dots(1),BlockTag='right')
    !
    init_called =.true.
  end subroutine init_dmrg


  !##################################################################
  !              FINALIZE DMRG ALGORITHM
  !##################################################################
  subroutine finalize_dmrg()
    integer :: ilat
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Finalize DMRG"
#endif
    if(allocated(HopH))deallocate(HopH)    
    do ilat=1,size(dots)
       call dots(ilat)%free()
    enddo
    deallocate(dots)
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
  !              RUN DMRG ALGORITHM
  !##################################################################
  subroutine DMRG()
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Launching DMRG"
#endif
    if(.not.init_called)&
         stop "DMRG ERROR: DMRG not initialized. Call init_dmrg first."
    !
    select case(to_lower(DMRGtype))
    case('i');call infinite_DMRG()
    case('f');call finite_DMRG()
    case default;stop "ERROR DMRG: unsupported DMRGtype: DMRGtype !=['i','f']"
    end select
  end subroutine DMRG








  !##################################################################
  !              INFINITE/FINITE DMRG ALGORITHM
  !##################################################################
  subroutine infinite_DMRG()
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Infinite Algorithm - ",allocated(blocks_list)
#endif
    !
    left =init_left
    right=init_right
    !
    do while (left%length < Ldmrg)
       call step_dmrg('i')
       if(allocated(blocks_list))then
          blocks_list(1, left%length)=left
          blocks_list(2,right%length)=right
       endif
    enddo
  end subroutine infinite_DMRG



  subroutine finite_DMRG()
    integer :: i,im,current_L
    integer :: j
    logical :: ExitSweep
    !
#ifdef _DEBUG
    write(LOGfile,*)"DEBUG: Finite Algorithm"
#endif
    !
    if(PBCdmrg)stop "ERROR DMRG: PBC+Finite DMRG algorithm is not supported"
    !
    allocate(blocks_list(2,2*Ldmrg))
    !
    blocks_list(1,1)=init_left
    blocks_list(2,1)=init_left
    !
    print*,""
    print*,"RUN IN-FINITE DMRG:"
    print*,""
    call infinite_DMRG()
    !
    print*,""
    print*,"START FINITE DMRG:"
    print*,""
    !Finite DMRG: start sweep forth&back:
    do im=1,size(Msweep)
       if(Esweep(im)==0)then
          write(*,"(A,I3,I6)")"Sweep, M:",im,Msweep(im)
       else
          write(*,"(A,I3,F8.4)")"Sweep, E:",im,Esweep(im)
       endif
       !
       do while(left%length < 2*Ldmrg)
          right = blocks_list(2,2*Ldmrg - left%length)
          call step_dmrg('f',1,im)
          blocks_list(1,left%length) = left
       enddo
       do while(right%length < 2*Ldmrg)
          left = blocks_list(1,2*Ldmrg - right%length)
          call step_dmrg('f',2,im)
          blocks_list(2,right%length)= right
       enddo
       do while(left%length <= Ldmrg)
          right = blocks_list(2,2*Ldmrg - left%length)
          call step_dmrg('f',1,im)
          blocks_list(1,left%length) = left
       enddo
    enddo
    !
  end subroutine finite_DMRG












  !##################################################################
  !              WORKHORSE: STEP DMRG 
  !##################################################################
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
    logical          :: bool1,bool2,renormalize
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
    !
    if(.not.left%is_valid(.true.))stop "single_dmrg_step error: left is not a valid block"
    if(.not.right%is_valid(.true.))stop "single_dmrg_step error: right is not a valid block"
    !
    !> dimension of the incoming L/R BLOCKS
    m_left  = left%dim
    Lleft   = left%length
    m_right = right%dim
    Lright  = right%length
    !
    !> START DMRG STEP:
    call dmrg_graphic(iLabel)    
    call start_timer()
    !
    !
    !#################################
    !    Enlarge L/R BLOCKS: +1 DOT
    !#################################
    ! call enlarge_block(left,dots(left%length+1),grow='left')
    ! call enlarge_block(right,dots(left%length+2),grow='right')
    call enlarge_block(left)
    call enlarge_block(right)
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
    !#################################
    !      WRITE AND EXIT
    !#################################
    write(LOGfile,"(A,I12,12X,I12)")&
         "         Blocks Length               :",Lleft,Lright
    write(LOGfile,"(A,I12,12X,I12)")&
         "Enlarged Blocks Length               :",left%length,right%length
    write(LOGfile,"(A,I12,12X,I12)")&
         "Enlarged Blocks Dim                  :",m_eleft,m_eright
    write(LOGfile,"(A,"//str(size(current_target_QN))//"F24.15)")&
         "Target_QN                            :",current_target_QN
    write(LOGfile,"(A,3x,G24.15)")&
         "Total                                :",sum(current_target_QN)
    write(LOGfile,"(A,3x,G24.15)")&
         "Filling                              :",sum(current_target_QN)/current_L
    write(LOGfile,"(A,3x,G24.15)")&
         "Filling/Norb                         :",sum(current_target_QN)/current_L/Norb
    write(LOGfile,"(A,I12)")&
         "SuperBlock Length                    :",current_L
    write(LOGfile,"(A,I12,A2,I12,A1,F10.5,A1)")&
         "SuperBlock Dimension  (tot)          :", &
         m_sb," (",m_eleft*m_eright,")",100*dble(m_sb)/m_eleft/m_eright,"%"
    !
    !#################################
    !       DIAG SUPER-BLOCK
    !#################################
    call sb_diag()
    write(LOGfile,*)"- - - - - - - - - - - - - - - - - - - - -"
    select case(left%type())
    case ("fermion","f")
       write(LOGfile,"(A,"//str(Lanc_Neigen)//"F24.15)")&
            "Energies/N                           :",gs_energy/sum(current_target_QN)
    case ("spin","s")
       write(LOGfile,"(A,"//str(Lanc_Neigen)//"F24.15)")&
            "Energies/L                           :",gs_energy/current_L
    end select
    write(LOGfile,*)"- - - - - - - - - - - - - - - - - - - - -"
    !
    !#################################
    !      BUILD RDM
    !#################################
    call sb_get_rdm()
    !
    !#################################
    !    Renormalize BLOCKS
    !#################################
    if(to_lower(type)=='i')then
       renormalize = left%length/=Ldmrg
    else
       renormalize=.not.((left%length==Ldmrg+1).AND.(iLabel==1))
    endif
    !
    if(renormalize )then
       call renormalize_block('left',m_rleft)
       call renormalize_block('right',m_rright)
       write(LOGfile,"(A,I12,12X,I12,3X,A3,I6,4X,I6,A1)")&
            "Renormalized Blocks Dim              :",m_rleft,m_rright,"< (",m_left,m_right,")"
       write(LOGfile,"(A,L12,12X,L12)")&
            "Truncating                           :",Mstates<=m_left,Mstates<=m_right
       write(LOGfile,"(A,2ES24.15)")&
            "Truncation Errors                    :",truncation_error_left,truncation_error_right
    endif
    !
    !> STOP DMRG STEP:
    call stop_timer("dmrg_step")
    !
    call write_energy()
    call write_truncation()
    call write_entanglement()
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
    write(Eunit,*)left%length,gs_energy/current_L/Norb,right%length
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
         truncation_error_right/current_L/Norb,right%length
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
    write(Eunit,*)left%length,entropy,right%length
    close(Eunit)
  end subroutine write_entanglement

END MODULE DMRG_MAIN











