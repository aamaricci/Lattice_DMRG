MODULE MATRIX_SPARSE  
  USE SCIFOR, only: str,free_unit,zero,assert_shape,eye
  USE AUX_FUNCS
  implicit none
  private


  !SPARSE ROW OF THE SPARSE MATRIX: note this is dynamic array
  type sparse_row
     integer                          :: size
     real(8),dimension(:),allocatable :: vals
     integer,dimension(:),allocatable :: cols
  end type sparse_row

  !SPARSE MATRIX STRUCTURE
  type sparse_matrix
     type(sparse_row),dimension(:),allocatable :: row
     integer                                   :: Nrow=0
     integer                                   :: Ncol=0
   contains
     procedure,pass :: init      => sp_init_matrix
     procedure,pass :: free      => sp_free_matrix
     procedure,pass :: load      => sp_load_matrix
     procedure,pass :: dump      => sp_dump_matrix
     procedure,pass :: fast_insert=> sp_fast_insert_element
     procedure,pass :: insert    => sp_insert_element
     procedure,pass :: get       => sp_get_element
     procedure,pass :: show      => sp_show_matrix
     procedure,pass :: print     => sp_print_matrix
     procedure,pass :: spy       => sp_spy_matrix
     procedure,pass :: as_matrix => sp_as_matrix
     procedure,pass :: dgr       => sp_dgr_matrix
     procedure,pass :: nnz       => sp_nnz_matrix
  end type sparse_matrix


  interface sparse_matrix
     module procedure :: sp_construct_matrix
  end interface sparse_matrix

  interface as_sparse
     module procedure :: sp_construct_matrix
  end interface as_sparse

  interface sparse
     module procedure :: sp_construct_matrix
  end interface sparse

  !EQUALITY with scalar and function (A=B, A=cmplx)
  interface assignment(=)
     module procedure :: sp_matrix_equal_scalar
     module procedure :: sp_matrix_equal_matrix
  end interface assignment(=)

  !ADDITION
  interface operator (+)
     module procedure :: sp_plus_matrix
  end interface operator (+)

  !SUBTRACTION
  interface operator (-)
     module procedure :: sp_minus_matrix
  end interface operator (-)

  !SCALAR PRODUCT
  interface operator(*)
     module procedure :: sp_left_product_matrix_i
     module procedure :: sp_left_product_matrix_d
     !
     module procedure :: sp_right_product_matrix_i
     module procedure :: sp_right_product_matrix_d
  end interface operator(*)

  !SCALAR DIVISION
  interface operator(/)
     module procedure :: sp_right_division_matrix_i
     module procedure :: sp_right_division_matrix_d
  end interface operator(/)


  !KRONECKER PRODUCT
  interface operator(.x.)
     module procedure :: sp_kron_matrix
  end interface operator(.x.)

  interface sp_kron
     module procedure :: sp_restricted_kron_matrix
  end interface sp_kron


  !RETURN SHAPE OF THE SPARSE MATRIX [Nrow,Ncol]
  intrinsic :: shape
  interface shape
     module procedure :: sp_shape_matrix
  end interface shape

  intrinsic :: transpose
  interface transpose
     module procedure :: sp_transpose_matrix
  end interface transpose

  interface hconjg
     module procedure :: sp_dgr_matrix
  end interface hconjg

  public :: sparse_matrix
  public :: as_sparse
  public :: sparse
  public :: assignment(=)
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: operator(.x.)
  public :: sp_kron
  public :: shape
  public :: transpose
  public :: hconjg
  public :: sp_eye



contains       







  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  elemental subroutine sp_init_matrix(sparse,Nrow,Ncol)
    class(sparse_matrix),intent(inout) :: sparse
    integer,intent(in)                 :: Nrow,Ncol
    integer                            :: i
    !
    call sparse%free()
    ! 
    sparse%Nrow=Nrow
    sparse%Ncol=Ncol
    !
    if(allocated(sparse%row))deallocate(sparse%row)
    allocate(sparse%row(Nrow))
    do i=1,Nrow
       sparse%row(i)%size=0
       allocate(sparse%row(i)%vals(0)) !empty array
       allocate(sparse%row(i)%cols(0)) !empty array
    end do
    !
  end subroutine sp_init_matrix

  function sp_construct_matrix(matrix) result(self)
    real(8),dimension(:,:),intent(in) :: matrix
    type(sparse_matrix)               :: self
    call self%load(matrix)
  end function sp_construct_matrix



  !+------------------------------------------------------------------+
  !PURPOSE: free an entire sparse matrix
  !+------------------------------------------------------------------+
  elemental subroutine sp_free_matrix(sparse)    
    class(sparse_matrix),intent(inout) :: sparse
    integer                            :: i
    !
    do i=1,sparse%Nrow
       sparse%row(i)%Size  = 0
       if(allocated(sparse%row(i)%vals))deallocate(sparse%row(i)%vals)
       if(allocated(sparse%row(i)%cols))deallocate(sparse%row(i)%cols)
    enddo
    if(allocated(sparse%row))deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
  end subroutine sp_free_matrix




  !+------------------------------------------------------------------+
  !PURPOSE: load a dense matrix into a sparse one
  !+------------------------------------------------------------------+
  subroutine sp_load_matrix(sparse,matrix)
    class(sparse_matrix),intent(inout) :: sparse
    real(8),dimension(:,:),intent(in)  :: matrix
    integer                            :: i,j,Ndim1,Ndim2
    !
    call sparse%free()
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    call sparse%init(Ndim1,Ndim2)
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=0d0)call sp_insert_element(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix(sparse,matrix)
    class(sparse_matrix),intent(in)      :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    integer                              :: i,j,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    call assert_shape(matrix,[sparse%Nrow,sparse%Ncol],"sp_dump_matrix","Matrix")
    Matrix = 0d0
    do i=1,Ndim1
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = matrix(i,sparse%row(i)%cols(j)) + sparse%row(i)%vals(j)
       enddo
    enddo
  end subroutine sp_dump_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  function sp_as_matrix(sparse) result(matrix)
    class(sparse_matrix),intent(in)            :: sparse
    real(8),dimension(sparse%Nrow,sparse%Ncol) :: matrix
    integer                                    :: i,j
    matrix = 0d0
    do i=1,sparse%Nrow
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = sparse%row(i)%vals(j)
       enddo
    enddo
  end function sp_as_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element(sparse,value,i,j)
    class(sparse_matrix),intent(inout) :: sparse
    real(8),intent(in)                 :: value
    integer,intent(in)                 :: i,j
    integer                            :: column,pos
    logical                            :: iadd
    !
    column = j
    !
    iadd = .false.                          !check if column already exist
    if(any(sparse%row(i)%cols == column))then         !
       pos = binary_search(sparse%row(i)%cols,column) !find the position  column in %cols
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then                            !this column exists so just sum it up       
       sparse%row(i)%vals(pos)=sparse%row(i)%vals(pos) + value  !add up value to the current one in %vals
    else                                    !this column is new. increase counter and store it 
       call append(sparse%row(i)%vals,value)
       call append(sparse%row(i)%cols,column)
       sparse%row(i)%Size = sparse%row(i)%Size + 1
    endif
    !
    if(sparse%row(i)%Size > sparse%Ncol)stop "sp_insert_element ERROR: row%Size > sparse%Ncol"
    !
  end subroutine sp_insert_element

  !no addition
  subroutine sp_fast_insert_element(sparse,value,i,j)
    class(sparse_matrix),intent(inout) :: sparse
    real(8),intent(in)                 :: value
    integer,intent(in)                 :: i,j
    integer                            :: column,pos
    !
    column = j
    !
    call append(sparse%row(i)%vals,value)
    call append(sparse%row(i)%cols,column)
    sparse%row(i)%Size = sparse%row(i)%Size + 1
    !
    if(sparse%row(i)%Size > sparse%Ncol)stop "sp_insert_element ERROR: row%Size > sparse%Ncol"
    !
  end subroutine sp_fast_insert_element


  !+------------------------------------------------------------------+
  !PURPOSE: get the element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  function sp_get_element(sparse,i,j) result(value)
    class(sparse_matrix),intent(inout) :: sparse    
    integer,intent(in)                 :: i,j
    real(8)                            :: value
    integer                            :: pos
    value=0d0
    do pos=1,sparse%row(i)%size
       if(j==sparse%row(i)%cols(pos))value=sparse%row(i)%vals(pos)
    enddo
  end function sp_get_element


  !+------------------------------------------------------------------+
  !PURPOSE: return the number of non-zero elements
  !+------------------------------------------------------------------+
  function sp_nnz_matrix(sparse) result(nnz)
    class(sparse_matrix),intent(in) :: sparse    
    integer :: nnz
    integer :: i
    nnz=0
    if(.not.allocated(sparse%row))return
    do i=1,sparse%Nrow
       nnz=nnz + sparse%row(i)%size
    enddo
  end function sp_nnz_matrix





  !##################################################################
  !##################################################################
  !              SHOW  / SPY
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: show matrix
  !+------------------------------------------------------------------+
  subroutine sp_show_matrix(sparse,dble,fmt)
    class(sparse_matrix)      :: sparse
    logical,optional          :: dble
    character(len=*),optional :: fmt
    character(len=12)         :: fmt_
    logical                   :: dble_
    integer                   :: i,j,unit_,Ns,NN
    character(len=64)         :: format
    real(8)                   :: val
    unit_=6
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt) !ES10.3
    dble_=show_dble;if(present(dble))dble_=dble
    ! if(dble_)then
    format='('//str(fmt_)//',1x)'
    Ns=sparse%Nrow
    do i=1,sparse%Nrow
       do j=1,sparse%Ncol
          val = sp_get_element(sparse,i,j)
          write(unit_,"("//str(sparse%Ncol)//str(format)//")",advance='no')val!dreal(val)
       enddo
       write(unit_,*)
    enddo
    !    write(unit_,*)
    ! else       
    !    format='(A1,'//str(fmt_)//',A1,'//str(fmt_)//',A1,1x)'
    !    Ns=sparse%Nrow
    !    do i=1,sparse%Nrow
    !       do j=1,sparse%Ncol
    !          val = sp_get_element(sparse,i,j)
    !          write(unit_,"("//str(sparse%Ncol)//str(format)//")",advance='no')&
    !               "(",dreal(val),",",dimag(val),")"
    !       enddo
    !       write(unit_,*)
    !    enddo
    !    write(unit_,*)
    ! endif

  end subroutine sp_show_matrix

  subroutine sp_print_matrix(sparse,file)
    class(sparse_matrix) :: sparse
    character(len=*)     :: file
    integer              :: unit_
    integer              :: i,j,Ns
    character(len=64)    :: fmt
    real(8)              :: val
    open(free_unit(unit_),file=str(file))
    fmt='(ES10.3,1x)'
    Ns=sparse%Nrow
    do i=1,sparse%Nrow
       do j=1,sparse%Ncol
          val = sp_get_element(sparse,i,j)
          write(unit_,"("//str(sparse%Ncol)//str(fmt)//")",advance='no')val
       enddo
       write(unit_,*)
    enddo
    close(unit_)
  end subroutine sp_print_matrix



  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+  
  subroutine sp_spy_matrix(sparse,header)
    class(sparse_matrix)    :: sparse
    character ( len = * )   :: header
    integer                 :: N1,N2
    character ( len = 255 ) :: command_filename
    integer                 :: command_unit
    character ( len = 255 ) :: data_filename
    integer                 :: data_unit
    integer                 :: i, j
    integer                 :: nz_num
    !
    !  Create data file.
    !
    !
    N1 = sparse%Nrow
    N2 = sparse%Ncol
    data_filename = trim ( header ) // '_data.dat'
    open (unit=free_unit(data_unit), file = data_filename, status = 'replace' )
    nz_num = 0
    do i=1,N1
       do j=1,sparse%row(i)%size
          write(data_unit,'(2x,i6,2x,i6)') sparse%row(i)%cols(j),i
          nz_num = nz_num + 1
       enddo
    enddo
    close(data_unit)
    !
    !  Create command file.
    !
    command_filename = "plot_"//str(header)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)')'set title "',nz_num,' nonzeros for '//str(header)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename)//'" w p pt 5 ps 0.4 lc rgb "red"'
    close ( unit = command_unit )
    return
  end subroutine sp_spy_matrix





  !##################################################################
  !##################################################################
  !               SPARSE MATRIX KRON PRODUCT 
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Perform simple Kroenecker product of two sparse matrices
  !AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
  !+------------------------------------------------------------------+
  function sp_kron_matrix(A,B) result(AxB)
    type(sparse_matrix), intent(in) :: A,B
    type(sparse_matrix)             :: AxB
    integer                         :: i,icol,j,k,kcol,l
    integer                         :: indx_row,indx_col
    real(8)                         :: value
    call AxB%free()
    call AxB%init(a%Nrow*b%Nrow,a%Ncol*b%Ncol)
    do indx_row = 1,A%Nrow*B%Nrow
       k = mod(indx_row,B%Nrow);if(k==0)k=B%Nrow
       i = (indx_row-1)/B%Nrow+1
       do icol=1,A%row(i)%size
          j = A%row(i)%cols(icol)
          do kcol=1,B%row(k)%size
             l = B%row(k)%cols(kcol)
             indx_col = l + (j-1)*B%Ncol
             value    = A%row(i)%vals(icol)*B%row(k)%vals(kcol)
             !
             call append(AxB%row(indx_row)%vals,value)
             call append(AxB%row(indx_row)%cols,indx_col)
             AxB%row(indx_row)%Size = AxB%row(indx_row)%Size + 1
             !
          enddo
       enddo
    enddo
  end function sp_kron_matrix


  function sp_restricted_kron_matrix(A,B,states) result(AxB)
    type(sparse_matrix), intent(in)  :: A,B
    integer,dimension(:),intent(in)  :: states
    type(sparse_matrix)              :: AxB
    integer                          :: i,icol,j,k,kcol,l,istate,jstate
    integer                          :: indx_row,indx_col
    real(8)                          :: value
    integer,dimension(:),allocatable :: inv_states
    call AxB%free()
    call AxB%init(size(states),size(states))
    allocate(inv_states(A%Ncol*B%Ncol))
    inv_states=0
    do i=1,size(states)
       inv_states(states(i)) = i
    enddo
    do istate = 1,size(states)
       indx_row=states(istate)
       k = mod(indx_row,B%Nrow);if(k==0)k=B%Nrow
       i = (indx_row-1)/B%Nrow+1
       do icol=1,A%row(i)%size
          j = A%row(i)%cols(icol)
          do kcol=1,B%row(k)%size
             l = B%row(k)%cols(kcol)
             indx_col = l + (j-1)*B%Ncol
             jstate   = inv_states(indx_col)
             value    = A%row(i)%vals(icol)*B%row(k)%vals(kcol)
             !
             call append(AxB%row(istate)%vals,value)
             call append(AxB%row(istate)%cols,jstate)
             AxB%row(istate)%Size = AxB%row(istate)%Size + 1
             !
          enddo
       enddo
    enddo
  end function sp_restricted_kron_matrix



  !##################################################################
  !##################################################################
  !               SPARSE MATRIX BASIC ALGEBRA 
  !##################################################################
  !##################################################################
  function sp_dgr_matrix(a) result(c)
    class(sparse_matrix), intent(in) :: a
    type(sparse_matrix)              :: c
    integer                          :: col
    real(8)                          :: val
    integer                          :: i,j    
    call c%init(a%Ncol,a%Nrow)      !tranpose
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)
          call c%insert(val,col,i)
       enddo
    enddo
  end function sp_dgr_matrix


  function sp_transpose_matrix(a) result(c)
    class(sparse_matrix), intent(in) :: a
    type(sparse_matrix)              :: c
    integer                          :: col
    real(8)                          :: val
    integer                          :: i,j    
    call c%init(a%Ncol,a%Nrow)      !tranpose
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)
          call c%insert(val,col,i)
       enddo
    enddo
  end function sp_transpose_matrix


  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix equality spA = spB. Deep copy
  !+------------------------------------------------------------------+
  subroutine sp_matrix_equal_matrix(a,b)
    type(sparse_matrix),intent(inout) :: a
    type(sparse_matrix),intent(in)    :: b
    integer                           :: col
    real(8)                           :: val
    integer                           :: i,j    
    call a%free()
    call a%init(b%Nrow,b%Ncol)
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val = b%row(i)%vals(j)
          call a%insert(val,i,col)
       enddo
    enddo
  end subroutine sp_matrix_equal_matrix


  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix scalar equality spA = const. 
  !+------------------------------------------------------------------+
  subroutine sp_matrix_equal_scalar(a,c)
    type(sparse_matrix),intent(inout) :: a
    real(8),intent(in)                :: c
    integer                           :: i,j    
    ! if(.not.a%status)stop "sp_matrix_equal_scalar error: a is not allocated"
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          a%row(i)%vals(j) = c
       enddo
    enddo
  end subroutine sp_matrix_equal_scalar



  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix addition spA + spB = spC
  !+------------------------------------------------------------------+
  function sp_plus_matrix(a,b) result(c)
    type(sparse_matrix), intent(in) :: a,b
    type(sparse_matrix)             :: c
    integer                         :: col
    real(8)                         :: val
    integer                         :: i,j    
    ! if(.not.a%status)stop "sp_plus_matrix error: a is not allocated"
    ! if(.not.b%status)stop "sp_plus_matrix error: b is not allocated"
    if(a%Nrow/=b%Nrow)stop "sp_plus_matrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "sp_plus_matrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val = b%row(i)%vals(j)
          call c%insert(val,i,col)
       enddo
    enddo
  end function sp_plus_matrix



  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix difference spA - spB = spC
  !+------------------------------------------------------------------+
  function sp_minus_matrix(a,b) result(c)
    type(sparse_matrix), intent(in) :: a,b
    type(sparse_matrix)             :: c
    integer                         :: col
    real(8)                         :: val
    integer                         :: i,j    
    if(a%Nrow/=b%Nrow)stop "sp_plus_matrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "sp_plus_matrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val =-b%row(i)%vals(j)
          call c%insert(val,i,col)
       enddo
    enddo
  end function sp_minus_matrix




  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix left scalar product const*spA = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function sp_left_product_matrix_i(C,A) result(B)
    integer,intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function sp_left_product_matrix_i

  function sp_left_product_matrix_d(C,A) result(B)
    real(8),intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function sp_left_product_matrix_d





  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar product spA*const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function sp_right_product_matrix_i(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function sp_right_product_matrix_i

  function sp_right_product_matrix_d(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function sp_right_product_matrix_d





  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar division spA/const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function sp_right_division_matrix_i(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)/C
          call b%insert(val,i,col)
       enddo
    enddo
  end function sp_right_division_matrix_i

  function sp_right_division_matrix_d(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)/C
          call b%insert(val,i,col)
       enddo
    enddo
  end function sp_right_division_matrix_d




  !+------------------------------------------------------------------+
  !PURPOSE:  Return the shape of a sparse matrix
  !+------------------------------------------------------------------+
  function sp_shape_matrix(self) result(shape)
    class(sparse_matrix),intent(in) :: self
    integer,dimension(2)            :: shape
    shape = [self%Nrow,self%Ncol]
  end function sp_shape_matrix



  !+------------------------------------------------------------------+
  !PURPOSE:  Return the identiy sparse matrix of given dimension
  !+------------------------------------------------------------------+
  function sp_eye(ndim) result(self)
    type(sparse_matrix) :: self
    integer             :: ndim
    call self%load(eye(ndim))
  end function sp_eye





end module MATRIX_SPARSE










! !+------------------------------------------------------------------+
! !PURPOSE:  Perform simple Kroenecker product of two sparse matrices
! !AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
! !+------------------------------------------------------------------+
! function sp_kron_matrix(A,B) result(AxB)
!   type(sparse_matrix), intent(in) :: A,B
!   type(sparse_matrix)             :: AxB
!   integer                         :: i,icol,j,k,kcol,l
!   integer                         :: indx_row,indx_col
!   complex(8)                      :: value
!   call AxB%free()
!   call AxB%init(a%Nrow*b%Nrow,a%Ncol*b%Ncol)
!   do i=1,A%Nrow
!      do k=1,B%Nrow

!         do icol=1,A%row(i)%size
!            j = A%row(i)%cols(icol)
!            do kcol=1,B%row(k)%size
!               l = B%row(k)%cols(kcol)
!               !
!               indx_row = k + (i-1)*B%Nrow
!               indx_col = l + (j-1)*B%Ncol
!               value    = A%row(i)%vals(icol)*B%row(k)%vals(kcol)
!               call AxB%insert(value,indx_row,indx_col)
!               !
!            enddo
!         enddo
!      enddo
!   enddo
! end function sp_kron_matrix
