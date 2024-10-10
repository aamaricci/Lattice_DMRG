program hm_1d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter                           :: Nspin=1
  integer :: Norb,Nso
  integer                                     :: Nk,Nktot,Nkpath,Nkx,Npts,L
  integer                                     :: Nlat,Nx
  integer                                     :: i,j,k,ik,iorb,jorb,ispin,io,jo
  integer                                     :: ilat,jlat
  integer                                     :: ix
  real(8)                                     :: kx
  real(8),dimension(:,:),allocatable          :: kgrid,kpath,ktrims,Rgrid
  integer,dimension(:,:),allocatable          :: Links
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:,:,:),allocatable   :: Hlat
  real(8)                                     :: mh,lambda
  real(8)                                     :: xmu,beta,eps,Eout(2)
  real(8),allocatable                                     :: dens(:)
  complex(8),allocatable                                  :: Hloc(:,:)
  complex(8) :: arg
  complex(8),dimension(:,:,:,:,:),allocatable :: Gmats,Greal
  complex(8),dimension(:,:,:,:,:,:),allocatable :: iGmats,iGreal
  character(len=20)                           :: file
  logical                                     :: iexist,fftflag
  complex(8),dimension(:,:),allocatable               :: GammaX,Gamma0,GammaZ
  complex(8),dimension(:,:,:),allocatable     :: ftHk
  complex(8),dimension(:,:,:,:),allocatable   :: ftHlat
  real(8),dimension(1)                        :: vecK,vecRi,vecRj

  call parse_input_variable(norb,"NORB","inputHM.conf",default=1)
  call parse_input_variable(nkx,"NKX","inputHM.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputHM.conf",default=500)
  call parse_input_variable(L,"L","inputHM.conf",default=2048)
  call parse_input_variable(mh,"MH","inputHM.conf",default=1d0)
  call parse_input_variable(lambda,"LAMBDA","inputHM.conf",default=0.3d0)
  call parse_input_variable(xmu,"XMU","inputHM.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputHM.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputHM.conf",default=1000.d0)
  call parse_input_variable(fftflag,"FFTFLAG","inputHM.conf",default=.false.)
  call parse_input_variable(file,"FILE","inputHM.conf",default="hkfile_hm.in")
  call save_input_file("inputHM.conf")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")

  Nso = Nspin*Norb
  if(Norb>2)stop "This code is for Norb<=2. STOP"


  Nktot= Nkx
  !
  Nx   = Nkx
  Nlat = Nx



  !SETUP THE GAMMA MATRICES:
  allocate(gammaX(Nso,Nso))
  allocate(gammaZ(Nso,Nso))
  allocate(gamma0(Nso,Nso))
  select case(Norb)
  case(1)
     gammaX=zero
     gamma0=one
     gammaZ=zero
  case(2)
     gammaX=pauli_tau_x
     gamma0=pauli_tau_0
     gammaZ=pauli_tau_z
  end select


  call TB_set_ei([1d0,0d0])
  call TB_set_bk([pi2,0d0])


  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:  
  write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nso,Nso,Nktot))
  call TB_build_model(Hk,hk_model,Nso,[Nkx])


  !GET LOCAL PART OF THE HAMILTONIAN
  Hloc=sum(Hk,dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc)


  !SOLVE ALONG A PATH IN THE BZ.
  Npts=3
  allocate(kpath(Npts,3))
  kpath(1,:)=-kpoint_X1
  kpath(2,:)= kpoint_Gamma
  kpath(3,:)= kpoint_X1
  call TB_Solve_model(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[red1,blue1,red1,blue1],&
       points_name=[character(len=20) :: '-X', 'G', 'X'],&
       file="Eigenband.nint")




  !Build the local GF:
  allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
  allocate(Greal(Nspin,Nspin,Norb,Norb,L))
  call get_gloc(Hk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L),axis='mats')
  call get_gloc(Hk,Greal,zeros(Nspin,Nspin,Norb,Norb,L),axis='real')
  call write_gf(Gmats,"Gloc",axis='mats',iprint=1)
  call write_gf(Greal,"Gloc",axis='real',iprint=1)



  !Occupation:
  allocate(dens(Nso))
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb+(ispin-1)*Norb
        dens(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
     enddo
  enddo
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(dens(iorb),iorb=1,Nso),sum(dens)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(dens(iorb),iorb=1,Nso),sum(dens)



  ! !Kinetic Energy
  ! call dmft_kinetic_energy(Hk,zeros(Nspin,Nspin,Norb,Norb,L))





  !##################################################################

  deallocate(Hk)

  allocate(Hlat(Nso,Nso,Nlat,Nlat))
  allocate(Hk(Nso,Nso,Nktot))


  !>Build direct Lattice Hamiltonian
  allocate(Links(2,1))          !Links: right,left
  Links(1,:) = [1]
  Links(2,:) =-[1]
  call TB_build_model(Hlat,ts_model,Nso,[Nx],Links,pbc=.true.)
  !
  !>Build reciprocal Lattice Hamiltonian
  call TB_build_model(Hk,hk_model,Nso,[Nkx])



  if(fftflag)then
     allocate(iGmats(Nlat,Nspin,Nspin,Norb,Norb,L))
     allocate(iGreal(Nlat,Nspin,Nspin,Norb,Norb,L))
     allocate(ftHk(Nlat*Nso,Nlat*Nso,1))

     do ilat=1,Nlat
        do jlat=1,Nlat
           do io=1,Nso
              do jo=1,Nso
                 i = io + (ilat-1)*Nso
                 j = jo + (jlat-1)*Nso
                 ftHk(i,j,1) = Hlat(io,jo,ilat,jlat)
              enddo
           enddo
        enddo
     enddo

     ! call get_gloc(ftHk,iGmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,L),axis='mats')
     ! call write_gf(iGmats,"ftGij",axis='mats',iprint=4)

     ! call get_gloc(ftHk,iGreal,zeros(Nlat,Nspin,Nspin,Norb,Norb,L),axis='real')
     ! call write_gf(iGreal,"ftGij",axis='real',iprint=4)
     deallocate(iGmats,iGreal,ftHk)
     print*,""
     print*,""
     print*,""


     allocate(Kgrid(Nktot,1))
     call TB_build_kgrid([Nkx],Kgrid)
     allocate(Rgrid(Nlat,1))
     call TB_build_Rgrid([Nx],Rgrid)



     !Build up the K-space Hamiltonian thru FT
     print*,"RECIPROCAL SPACE: H_bhz(kx,ky)=FT[H_bhz(i,j)]"
     allocate(ftHk(Nso,Nso,Nktot))
     allocate(ftHlat(Nso,Nso,Nlat,Nlat))


     !> from Hlat --> Hk = FT(Hlat)
     ftHk=zero
     call start_timer
     do ik=1,Nktot
        vecK = Kgrid(ik,:)
        do ilat=1,Nlat
           vecRi = Rgrid(ilat,:)
           do jlat=1,Nlat
              vecRj = Rgrid(jlat,:)
              !
              arg=dot_product(vecK,vecRj-vecRi)
              !
              ftHk(:,:,ik)=ftHk(:,:,ik) + exp(-xi*arg)*Hlat(:,:,ilat,jlat)/Nlat
              !
           enddo
        enddo
        call eta(ik,Nktot)
     enddo
     where(abs(ftHk)<1.d-6)ftHk=zero
     call stop_timer


     ftHk = ftHk-Hk
     print*,"Identity test for Hk:",sum(abs(ftHk))
     if(sum(abs(ftHk))>1d-6)then
        print*,"Failed! :-("
     else
        print*,"Passed! :-)"
     endif

     print*,""
     print*,""
     print*,""




     !Build up the real-space Hamiltonian thru FT:"
     print*,"REAL SPACE:       H_bhz(i,j)=FT^-1(H_bhz(kx,ky))"
     ftHlat=zero
     call start_timer
     do ik=1,Nktot
        vecK = Kgrid(ik,:)
        !
        do ilat=1,Nlat
           vecRi = Rgrid(ilat,:)
           do jlat=1,Nlat
              vecRj = Rgrid(jlat,:)
              !
              arg=dot_product(vecK,vecRj-vecRi)
              !
              ftHlat(:,:,ilat,jlat)= ftHlat(:,:,ilat,jlat) + exp(xi*arg)*Hk(:,:,ik)/Nktot
           enddo
        enddo
        call eta(ik,Nktot)
     enddo
     where(abs(Hlat)<1.d-6)Hlat=zero
     call stop_timer



     ftHlat = ftHlat-Hlat
     print*,"Identity test for Hlat:"
     if(sum(abs(ftHlat))>1d-6)then
        print*,"Failed! :-("
     else
        print*,"Passed! :-)"
     endif

  endif

contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx
    complex(8),dimension(N,N) :: hk
    if(N/=2)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    Hk = -cos(kx)*Gamma0 + Mh*GammaZ + lambda*cos(kx)*GammaX
  end function hk_model



  function ts_model(link,Nso) result(Hts)
    integer                       :: link
    integer                       :: Nso
    complex(8),dimension(Nso,Nso) :: Hts
    select case(link)
    case (0) !LOCAL PART
       Hts =  Mh*GammaZ
    case (1) !RIGHT HOPPING
       Hts = -0.5d0*Gamma0 + 0.5d0*lambda*GammaX
    case (2) !LEFT HOPPING
       Hts = -0.5d0*Gamma0 + 0.5d0*lambda*GammaX
    case default 
       stop "ts_model ERROR: link != [0:2]"
    end select
  end function ts_model














  subroutine print_Hlat(Hij,file)
    complex(8),dimension(Nso,Nso,Nlat,Nlat) :: Hij
    integer                                 :: ilat,jlat,iso,jso,unit
    character(len=*)                        :: file
    open(free_unit(unit),file=trim(file))
    do ilat=1,Nlat
       do iso=1,Nso
          do jlat=1,Nlat
             do jso=1,Nso
                write(unit,"(F5.2,1x)",advance="no")dreal(Hij(iso,jso,ilat,jlat))
             enddo
             write(unit,"(A2)",advance="no")"  "
          enddo
          write(unit,*)
       enddo
       write(unit,*)""
    enddo
    write(unit,*)
    close(unit)
  end subroutine print_Hlat


  subroutine print_Hk(Hk,file)
    complex(8),dimension(:,:,:)     :: Hk ![Nso][Nso][Nk]
    integer                         :: Nso,Nk,unit
    integer                         :: ik,iso,jso
    character(len=*)                :: file
    open(free_unit(unit),file=trim(file))
    Nk = size(Hk,3)
    Nso= size(Hk,1)
    call assert_shape(Hk,[Nso,Nso,Nk],"print_Hk","Hk")
    do ik=1,Nk
       do iso=1,Nso
          write(unit,"(1000(A1,F6.3,A1,F6.3,A5))")&
               ("(",dreal(Hk(iso,jso,ik)),",",dimag(Hk(iso,jso,ik)),")    ",jso=1,Nso)
       enddo
       write(unit,*)
    enddo
    write(unit,*)
    close(unit)
  end subroutine print_Hk





end program hm_1d


