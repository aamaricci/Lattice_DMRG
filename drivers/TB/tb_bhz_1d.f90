program bhz_1d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                   :: Norb=2,Nspin=2,Nso
  integer                                   :: Nk,Nkpath,Nkx
  integer                                   :: Nlat,Nx,Nup,Ndw
  integer                                   :: i,j,k,ik,iorb,jorb,ispin,io,jo
  integer                                   :: ilat,jlat
  integer                                   :: ix
  real(8)                                   :: kx
  real(8),dimension(:,:),allocatable        :: kgrid,kpath,ktrims,Rgrid
  integer,dimension(:,:),allocatable        :: Links
  complex(8),dimension(:,:,:),allocatable   :: Hk
  complex(8),dimension(:,:),allocatable     :: Hij,Hloc
  complex(8),dimension(:,:,:,:),allocatable :: Hlat
  real(8),dimension(:),allocatable          :: Eij
  real(8),dimension(:),allocatable          :: rhoDiag,kdens,dens
  complex(8),dimension(:,:),allocatable     :: rhoH
  real(8)                                   :: mh,eh,lambda,E0
  real(8)                                   :: xmu,beta
  character(len=20)                         :: file
  logical                                   :: pbc
  complex(8),dimension(4,4)                 :: Gamma1,Gamma5

  call parse_input_variable(nkx,"NKX","input.conf",default=100)
  call parse_input_variable(mh,"MH","input.conf",default=0d0)
  call parse_input_variable(lambda,"LAMBDA","input.conf",default=0d0)
  call parse_input_variable(eh,"EH","input.conf",default=1.d0)
  call parse_input_variable(pbc,"pbc","input.conf",default=.false.)
  call parse_input_variable(beta,"beta","input.conf",default=1000d0)
  call save_input_file("input.conf")


  if(Nspin/=2.OR.Norb/=2)&
       stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb

  gamma1=kron( pauli_sigma_z, pauli_tau_x)
  gamma5=kron( pauli_sigma_0, pauli_tau_z)


  call TB_set_ei([1d0,0d0])
  call TB_set_bk([pi2,0d0])

  allocate(Links(2,1))          !Links: right,left
  Links(1,:) = [1]
  Links(2,:) =-[1]

  allocate(kdens(Nso),dens(Nso))
  do Nx=2,Nkx,2
     Nlat = Nx

     !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:  
     write(*,*) "Using Nk_total="//txtfy(Nx)
     if(allocated(Hk))deallocate(Hk)
     allocate(Hk(Nso,Nso,Nx))
     call TB_build_model(Hk,hk_model,Nso,[Nx])

     !GET LOCAL PART OF THE HAMILTONIAN
     Hloc=sum(Hk,dim=3)/Nx
     where(abs(Hloc)<1d-6)Hloc=zero
     call TB_write_Hloc(Hloc)

     kdens = get_dens_Hk()
     ! write(200,*)Nx,kdens

     ! !Kinetic Energy
     ! call dmft_kinetic_energy(Hk,zeros(Nspin,Nspin,Norb,Norb,L))



     !##################################################################
     if(allocated(Hlat))deallocate(Hlat)
     allocate(Hlat(Nso,Nso,Nlat,Nlat))
     Hlat=zero

     !>Build direct Lattice Hamiltonian

     call TB_build_model(Hlat,ts_model,Nso,[Nlat],Links,pbc=pbc)


     if(allocated(Hij))deallocate(Hij)
     if(allocated(Eij))deallocate(Eij)
     if(allocated(rhoDiag))deallocate(rhoDiag)
     if(allocated(rhoH))deallocate(rhoH)
     allocate(Hij(Nlat*Nso,Nlat*Nso))
     allocate(rhoH(Nlat*Nso,Nlat*Nso))
     allocate(Eij(Nlat*Nso))
     allocate(rhoDiag(Nlat*Nso))

     Hij = zero
     do ilat=1,Nlat
        do io=1,Nso
           i = io + (ilat-1)*Nso
           do jlat=1,Nlat
              do jo=1,Nso
                 j = jo + (jlat-1)*Nso
                 Hij(i,j) = Hlat(io,jo,ilat,jlat)
              enddo
           enddo
           !
        enddo
     enddo

     if(Nlat<10)call print_matrix(Hij)

     call eigh(Hij,Eij)

     rhoDiag = fermi(Eij,beta)
     rhoH    = matmul(Hij , matmul(diag(rhoDiag), conjg(transpose(Hij))) ) 

     dens = 0d0
     do ilat=1,Nlat
        do io=1,Nso
           i = io + (ilat-1)*Nso
           dens(io) = dens(io) + rhoH(i,i)/Nlat
        enddo
     enddo


     Nup = Nlat*Nso/2
     Ndw = Nlat*Nso-Nup

     E0 = sum(Eij(:Nup)) + sum(Eij(:Ndw))
     E0 = E0/Nlat/Nso
     write(*,*)Nlat/2,E0,dens(1),dens(2),kdens(1),kdens(2)!,sum(Eij(:Nup))/Nlat/Nso,sum(Eij(:Ndw))/Nlat/Nso

     write(100,*)Nlat/2,E0,dens(1),dens(2),kdens(1),kdens(2)
  enddo



contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8) :: ek
    real(8)                   :: kx
    complex(8),dimension(N,N) :: hk
    ! if(N/=2)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ek = -eh*cos(kx)
    Hk = (ek+Mh)*Gamma5 + xi*lambda*cos(kx)*Gamma1
  end function hk_model



  function get_dens_hk() result(dens)
    real(8),dimension(Nso)        :: dens
    real(8),dimension(Nso)        :: E,rhoDiag
    complex(8),dimension(Nso,Nso) :: H,RhoH
    integer                       :: ik,io
    rhoH = zero
    do ik=1,Nx
       H = Hk(:,:,ik)
       call eigh(H,E)
       !
       rhoDiag = one*fermi(E,beta)
       rhoH    = rhoH + matmul(H , matmul(diag(rhoDiag), conjg(transpose(H))) )/Nx
       !       
    enddo
    dens = diagonal(rhoH)
  end function get_dens_hk

  function ts_model(link,Nso) result(Hts)
    integer                       :: link
    integer                       :: Nso
    complex(8),dimension(Nso,Nso) :: Hts
    select case(link)
    case (0) !LOCAL PART
       Hts = Mh*Gamma5
    case (1) !RIGHT HOPPING
       Hts = -eh/2d0*Gamma5 + xi*lambda/2d0*Gamma1
    case (2) !LEFT HOPPING
       Hts = -eh/2d0*Gamma5 + xi*lambda/2d0*Gamma1
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

end program bhz_1d


