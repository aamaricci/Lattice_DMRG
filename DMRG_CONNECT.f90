MODULE DMRG_CONNECT
  USE VARS_GLOBAL
  implicit none
  private


  public :: connect_fermion_blocks
  public :: connect_spin_blocks


contains



  !H_lr = \sum_{a}h_aa*(C^+_{left,a}@P_left) x C_{right,a}] + H.c.
  function connect_fermion_blocks(left,right,states) result(H2)
    type(block)                               :: left
    type(block)                               :: right
    integer,dimension(:),optional             :: states
    type(sparse_matrix),dimension(Nspin*Norb) :: Cl,Cr
    type(sparse_matrix)                       :: P,A
    type(sparse_matrix)                       :: H2
    integer                                   :: ispin,iorb,jorb,io,jo
    integer,dimension(2)                      :: Hdims,dleft,dright
    character(len=:),allocatable              :: key
    real(8),dimension(:,:),allocatable        :: Hij
    !
    !Hij is shared:
    Hij = Hmodel(left,right)
    !
    !> Get H2 dimensions:
    dleft = shape(left%operators)
    dright= shape(right%operators)
    Hdims = dleft*dright
    if(present(states))Hdims = [size(states),size(states)]
    call H2%init(Hdims(1),Hdims(2))
    !
    !FERMION SPECIFIC:
    !>Retrieve operators:
    do ispin=1,Nspin
       do iorb=1,Norb
          key = "C"//left%okey(iorb,ispin)
          io = iorb + (ispin-1)*Norb
          Cl(io) = left%operators%op(key)
          Cr(io) = right%operators%op(key)
       enddo
    enddo
    !
    !
    !>Build H2:
    P = left%operators%op("P")  !always acts on the Left Block
    do io=1,Nspin*Norb
       do jo=1,Nspin*Norb
          if(Hij(io,jo)==0d0)cycle
          if(present(states))then
             H2 = H2 + Hij(io,jo)*sp_kron(matmul(Cl(io)%dgr(),P),Cr(jo),states)
             H2 = H2 + Hij(io,jo)*sp_kron(matmul(P,Cl(io)),Cr(jo)%dgr(),states)
          else
             H2 = H2 + Hij(io,jo)*(matmul(Cl(io)%dgr(),P).x.Cr(jo))
             H2 = H2 + Hij(io,jo)*(matmul(P,Cl(io)).x.Cr(jo)%dgr())
          endif

       enddo
    enddo
    !
    !> free memory
    call P%free
    do io=1,Nspin*Norb
       call Cl(io)%free
       call Cr(io)%free
    enddo
  end function connect_fermion_blocks




  function connect_spin_blocks(left,right,states) result(H2)
    type(block)                        :: left
    type(block)                        :: right
    integer,dimension(:),optional      :: states
    type(sparse_matrix)                :: Sl(Nspin)![Sz,Sp]
    type(sparse_matrix)                :: Sr(Nspin)![Sz,Sp]
    type(sparse_matrix)                :: H2
    integer,dimension(2)               :: Hdims
    integer                            :: ispin
    real(8),dimension(:,:),allocatable :: Hij
    !
    !Hij is shared:
    Hij = Hmodel(left,right)
    !
    !> Get H2 dimensions:
    Hdims = shape(left%operators)*shape(right%operators)
    if(present(states))Hdims = [size(states),size(states)]
    call H2%init(Hdims(1),Hdims(2))
    !
    !>Retrieve operators:
    do ispin=1,Nspin
       Sl(ispin) = left%operators%op("S"//left%okey(0,ispin))
       Sr(ispin) = right%operators%op("S"//right%okey(0,ispin))
    enddo
    !
    !
    !>Build H2:
    if(present(states))then
       H2 = H2 + Hij(1,1)*sp_kron(Sl(1),Sr(1),states) + &
            Hij(2,2)*sp_kron(Sl(2),Sr(2)%dgr(),states)+ &
            Hij(2,2)*sp_kron(Sl(2)%dgr(),Sr(2),states)
    else
       H2 = H2 + Hij(1,1)*(Sl(1).x.Sr(1)) + &
            Hij(2,2)*(Sl(2).x.Sr(2)%dgr())+ &
            Hij(2,2)*(Sl(2)%dgr().x.Sr(2))
    endif
    !
    !> Free memory
    do ispin=1,Nspin
       call Sl(ispin)%free
       call Sr(ispin)%free
    enddo
  end function connect_spin_blocks




END MODULE DMRG_CONNECT
