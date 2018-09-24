#include "intrinsic_sizes.h"
#include "redef.h"
#undef DEBUG
#ifdef DEBUG
#define DEBSTART(string) print*,"== START ",string," =="
#define DEBEND(string) print*,"== END   ",string," =="
#else
#define DEBSTART(string) 
#define DEBEND(string)
#endif

MODULE StoreBandedMatrixModule
  USE ProcessGridModule
  USE StoreFullMatrixModule
  use StoreVectorModule
  use MatrixEntryListModule
  use MatrixEntryNodeModule
  use IntegerListModule
  use IntegerNodeModule
  USE mpi
  IMPLICIT NONE

  PUBLIC :: initialize, finalize
  PUBLIC :: isInitialized, ALLOCATE, attach, set_value, add_value, commit_values, autotune,dot_multiply
  PUBLIC :: add_matrix, subtract_matrix, multiply_matrix_with_scalar,  my_sum
  PUBLIC :: get_global_matrix_locally, mat_get_value, get_local_abs_square_sum, show
  public :: mat_get_row_pointer
  PUBLIC :: static_size_of_StoreBandedMatrixObject, set_zero, get_number_of_bands
  PUBLIC :: square_matrix, row_axpy, LU_factor, LU_factor_ok, LU_solve
  public :: print_storage_details, transpose_and_conjugate, transpose_storage, assignment(=)
  complex,parameter :: complex_kind=(1.0D0,1.0D0)
  PRIVATE

  TYPE,public :: StoreBandedMatrixObject 
     TYPE(ProcessGrid) :: PG
     !> the per processor local matrix part of the global matrix
     COMPLEX, DIMENSION(:,:), POINTER :: localMatrix 
     !> data needed for LU factorization
     COMPLEX, DIMENSION(:,:), POINTER :: factoredBand
     COMPLEX, DIMENSION(:),   POINTER :: auxData
     INTEGER, DIMENSION(:),   POINTER :: ipiv

     INTEGER :: NRows, NCols ! number of rows and cols globally
     INTEGER :: LowerBandwidth, UpperBandwidth  ! bandwidths without main diagonal
     INTEGER :: ColsPerBlock ! block size of distributed matrix
     INTEGER :: NumberOfStoredRows
     LOGICAL :: isInitialized=.FALSE.
     LOGICAL :: isTransposed = .FALSE. !< determines storage scheme. true means lda is a row of the matrix
     LOGICAL :: isFactored = .FALSE.   !< determines if matrix is LU factored
     INTEGER,DIMENSION(7) :: Desc      !< matrix descriptor for BLACS, ScaLapack etc. Currently unused!
     !> A temporary list for storing the entries of the matrix, before a commit_values is called.
     TYPE(MatrixEntryList) :: storage_list
  END TYPE StoreBandedMatrixObject

 TYPE,public :: StoreBandedMatrixObjectReal 
     TYPE(ProcessGrid) :: PG
     !> the per processor local matrix part of the global matrix
     REAL, DIMENSION(:,:), POINTER :: localMatrix 
     !> data needed for LU factorization
     REAL, DIMENSION(:,:), POINTER :: factoredBand
     REAL, DIMENSION(:),   POINTER :: auxData
     INTEGER, DIMENSION(:),   POINTER :: ipiv

     INTEGER :: NRows, NCols ! number of rows and cols globally
     INTEGER :: LowerBandwidth, UpperBandwidth  ! bandwidths without main diagonal
     INTEGER :: ColsPerBlock ! block size of distributed matrix
     INTEGER :: NumberOfStoredRows
     LOGICAL :: isInitialized=.FALSE.
     LOGICAL :: isTransposed = .FALSE. !< determines storage scheme. true means lda is a row of the matrix
     LOGICAL :: isFactored = .FALSE.   !< determines if matrix is LU factored
     INTEGER,DIMENSION(7) :: Desc      !< matrix descriptor for BLACS, ScaLapack etc. Currently unused!
     !> A temporary list for storing the entries of the matrix, before a commit_values is called.
     TYPE(MatrixEntryList) :: storage_list
  END TYPE StoreBandedMatrixObjectReal

  ! PG_cols is a process grid of size 1 X NProcs needed for LU Decomposition.
  ! It is needed because pcgbtrs doesn't accept the standerd processgrid PG
  ! which is NProcs X 1 (although the documentation states it differently).
  ! It should eventually go into TYPE StoreBandedMatrixObject ...
  TYPE(ProcessGrid), SAVE :: PG_cols

#ifdef FORTRAN2003_ENUM
  enum,bind(c)
     enumerator :: USE_BLAS=1, F_INTRINSIC, BY_INDEX, BY_ROW
  endenum
#else
  INTEGER, parameter :: USE_BLAS=1, F_INTRINSIC=2, BY_INDEX=3, BY_ROW=4
#endif

  integer :: MatVecMethod=F_INTRINSIC
  !> initializes a banded matrix
  INTERFACE initialize
     MODULE PROCEDURE sbm_initialize_matrix, sbm_initialize_matrix_real
  END INTERFACE

  INTERFACE finalize
     MODULE PROCEDURE sbm_finalize_matrix, sbm_finalize_matrix_real
  END INTERFACE
  
  INTERFACE isInitialized
     module procedure sbm_isInitialized, sbm_isInitialized_real
  END INTERFACE

  INTERFACE allocate
     module procedure sbm_allocate, sbm_allocate_real
  END INTERFACE

  INTERFACE set_value
     MODULE PROCEDURE sbm_set_value, sbm_set_value_real
  END INTERFACE

  INTERFACE add_value
     MODULE PROCEDURE sbm_add_value, sbm_add_value_real
  END INTERFACE

  INTERFACE commit_values
     MODULE PROCEDURE sbm_commit_values, sbm_commit_values_real
  END INTERFACE

  INTERFACE set_zero
     MODULE PROCEDURE sbm_set_zero, sbm_set_zero_real
  END INTERFACE

  INTERFACE get_global_matrix_locally
     module procedure sbm_get_global_matrix_locally, sbm_get_global_matrix_locally_real
  END INTERFACE

  INTERFACE mat_get_value
     MODULE PROCEDURE sbm_get_value, sbm_get_value_real
  END INTERFACE

  interface mat_get_row_pointer
     module procedure sbm_get_row_pointer, sbm_get_row_pointer_real
  end interface

  INTERFACE get_local_abs_square_sum
     module procedure sbm_get_local_abs_square_sum, sbm_get_local_abs_square_sum_real
  END INTERFACE

  interface autotune
     module procedure sbm_autotune, sbm_autotune_real
  end interface

  INTERFACE dot_multiply
     MODULE PROCEDURE sbm_dot_multiply, sbm_matmat_dot_multiply, sbm_dot_multiply_real, sbm_matmat_dot_multiply_real
     module procedure sbm_dot_multiply_Banded_with_Full, sbm_dot_multiply_Banded_with_Full_real
  END INTERFACE

  INTERFACE square_matrix
     MODULE PROCEDURE sbm_square_matrix, sbm_square_matrix_real
  END INTERFACE

  INTERFACE show
     MODULE PROCEDURE sbm_show_on_screen, sbm_show_in_file, sbm_show_on_screen_real, sbm_show_in_file_real
  END INTERFACE

  INTERFACE add_matrix
     MODULE PROCEDURE sbm_add_to_matrix, sbm_add_matrix, sbm_add_matrix_real,  sbm_add_to_matrix_real
  END INTERFACE

  INTERFACE subtract_matrix
     MODULE PROCEDURE sbm_subtract_matrix, sbm_subtract_from_matrix, sbm_subtract_matrix_real, sbm_subtract_from_matrix_real
  END INTERFACE

  INTERFACE multiply_matrix_with_scalar
     MODULE PROCEDURE sbm_multiply_matrix_with_real, sbm_scale_matrix_by_real
     MODULE PROCEDURE sbm_multiply_matrix_with_real_real, sbm_scale_matrix_by_real_real
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE sbm_assign_matrix, sbm_assign_matrix_real
  END INTERFACE

  INTERFACE get_number_of_bands
     module procedure sbm_get_number_of_bands, sbm_get_number_of_bands_real
  END INTERFACE

  INTERFACE my_sum
     MODULE PROCEDURE sbm_matrix_sum, sbm_matrix_sum_real
  END INTERFACE

  INTERFACE row_axpy
     module procedure sbm_row_axpy, sbm_row_axpy_real
  END INTERFACE

  interface transpose_and_conjugate
     module procedure sbm_transpose_and_conjugate, sbm_transpose_and_conjugate_real
  end interface

  interface transpose_storage
     module procedure sbm_transpose_storage, sbm_transpose_storage_real
  end interface

  INTERFACE LU_factor
     module procedure sbm_LU_factor, sbm_LU_factor_real
  END INTERFACE

  INTERFACE LU_factor_ok
     module procedure sbm_LU_factor_ok, sbm_LU_factor_ok_real
  END INTERFACE

  INTERFACE LU_solve
     module procedure sbm_LU_solve, sbm_LU_solve_real
  END INTERFACE

  interface exit_if_factored
     module procedure sbm_exit_if_factored, sbm_exit_if_factored_real
  end interface

  interface print_storage_details
     module procedure sbm_print_storage_details, sbm_print_storage_details_real
  end interface

#ifdef ALL_ROUTINES
  INTERFACE invert
     module procedure sbm_invert_matrix, sbm_invert_matrix_real
  END INTERFACE
#endif

CONTAINS
  FUNCTION static_size_of_StoreBandedMatrixObject() RESULT(memory_need)
    INTEGER :: memory_need

    memory_need = size_of_processgrid()
  END FUNCTION static_size_of_StoreBandedMatrixobject
  
  SUBROUTINE sbm_initialize_matrix(mat,n_rows,n_cols,PG,transp)
    TYPE(StoreBandedMatrixObject), intent(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    type(ProcessGrid) :: PG
    LOGICAL, INTENT(IN),OPTIONAL :: transp

    ! Local variables
  
    mat%NRows = n_rows
    mat%NCols = n_cols
    mat%PG = PG

    ! Distribute matrix on the grid
    ! block-column matrix
    !mat%ColsPerBlock = mat%NCols / PG%NPRows
    mat%ColsPerBlock = mat%NCols / PG%NProcs

    !mat%Desc = (/ 501, PG%Context, mat%NCols, mat%ColsPerBlock, 0, mat%LowerBandwidth+mat%UpperBandwidth+1,0 /)

    mat%LowerBandwidth = 0
    mat%UpperBandwidth = 0
    mat%NumberOfStoredRows = 1

    CALL initialize(mat%storage_list)

    mat%isInitialized = .TRUE.
    IF (PRESENT(transp)) THEN
       mat%isTransposed = transp
    ELSE
       mat%isTransposed = .false.
    END IF

    ! Get pointers into a defined state
    NULLIFY(mat%localMatrix)
    NULLIFY(mat%factoredBand)
    NULLIFY(mat%auxData)
    NULLIFY(mat%ipiv)

    CALL blacs_barrier(PG%Context,'A')
  END SUBROUTINE sbm_initialize_matrix

  LOGICAL FUNCTION sbm_isInitialized(mat) 
    type(StoreBandedMatrixObject) :: mat

    sbm_isInitialized = mat%isInitialized
  END FUNCTION sbm_isInitialized


  SUBROUTINE sbm_allocate(mat,lbw_in,ubw_in)
    TYPE(StoreBandedMatrixObject) :: mat
    integer, optional,intent(IN) :: lbw_in, ubw_in

    if (present(lbw_in)) then
       if (lbw_in.le.mat%NCols) then
          mat%LowerBandwidth = lbw_in
       else
          print*,"Lower Bandwidth can not be larger than dim of matrix!"
          !CALL tracebackqq()

          stop
       end if
    end if

    if (present(ubw_in)) then
       if (ubw_in.le.mat%NCols) then
          mat%UpperBandwidth = ubw_in
       else
          print*,"Upper Bandwidth can not be larger than dim of matrix!"
          stop
       end if
    end if

    if (associated(mat%localMatrix)) deallocate(mat%localMatrix)

    IF (mat%isTransposed) THEN
       ALLOCATE(mat%localMatrix(-mat%LowerBandwidth:mat%UpperBandwidth,mat%ColsPerBlock))
    ELSE
       ALLOCATE(mat%localMatrix(-mat%UpperBandwidth:mat%LowerBandwidth,mat%ColsPerBlock))
    END IF
    mat%NumberOfStoredRows = mat%LowerBandwidth+mat%UpperBandwidth+1
    mat%localMatrix = CMPLX(0.0,0.0,kind(complex_kind))
  END SUBROUTINE sbm_allocate

  SUBROUTINE sbm_finalize_matrix(mat)
    type(StoreBandedMatrixObject) :: mat

    if (ASSOCIATED(mat%localMatrix)) THEN
       DEALLOCATE(mat%localMatrix)
    END IF
    IF (ASSOCIATED(mat%factoredBand)) THEN
       DEALLOCATE(mat%factoredBand)
    END IF
    IF (ASSOCIATED(mat%auxData)) THEN
       DEALLOCATE(mat%auxData)
    END IF
    IF (ASSOCIATED(mat%ipiv)) THEN
       DEALLOCATE(mat%ipiv)
    END IF
    CALL finalize(mat%storage_list)
    mat%isInitialized=.false.
  END SUBROUTINE sbm_finalize_matrix

  SUBROUTINE sbm_set_value(mat,irow,icol,value)
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc, row_index, col_index
    INTEGER :: distance!, coord(2)
    TYPE(MatrixEntryNodeData) :: nodedata
    real :: tolerance

    CALL exit_if_factored(mat)

    DEBSTART("set_value(SBMO)")
    !tolerance = REAL(0.0,KIND(tolerance))
    ! only consider value if value is larger than tolerance (numeric zero)
    IF ( (REAL(value).NE.REAL(0.0,KIND(tolerance))) .OR. (AIMAG(value).NE.REAL(0.0,KIND(tolerance))) ) THEN
    !IF ((ABS(REAL(value)).GT.tolerance).OR.(ABS(AIMAG(value)).GT.tolerance)) THEN
       ! calculate the processor coordinates, just determined by the global column index
       ! or global row index, depending on the storage scheme
       IF (mat%isTransposed) THEN
          proc  = (irow-1)/mat%ColsPerBlock
       ELSE
          proc  = (icol-1)/mat%ColsPerBlock
       END IF

       ! only process the entry if we are on the right processor
       IF (proc.EQ.mat%PG%Rank) THEN
          ! check, if the value is inside the band 
          distance = icol-irow
          IF ((distance>mat%UpperBandwidth).OR.(-distance>mat%LowerBandwidth)) THEN
             ! out of bands
             ! put value and indices into list for later insertion
             nodedata%coord(1)=irow
             nodedata%coord(2)=icol
             nodedata%value = value
             CALL push(mat%storage_list,nodedata)
          ELSE
             IF (mat%isTransposed) THEN
                row_index = MOD(irow-1,mat%ColsPerBlock)+1
                col_index = distance

                mat%localMatrix(col_index, row_index) = value
             ELSE
                col_index = MOD(icol-1,mat%ColsPerBlock)+1
                row_index = -distance

                mat%localMatrix(row_index, col_index) = value
             END IF
          END IF
       END IF
    end if
    DEBEND("set_value(SBMO)")
  END SUBROUTINE sbm_set_value

  SUBROUTINE sbm_add_value(mat,global_row,global_col,value)
    TYPE(StoreBandedMatrixObject) :: mat
    INTEGER, INTENT(IN) :: global_row,global_col
    COMPLEX, INTENT(IN) :: value

    ! Local variables
    COMPLEX :: old_value
    integer :: proc

    DEBSTART("add_value(SFMO)")

    !CALL exit_if_factored(mat)

    IF (mat%isTransposed) THEN
       proc = (global_row-1)/mat%ColsPerBlock
    ELSE
       proc = (global_col-1)/mat%ColsPerBlock
    END IF

    IF (proc.EQ.mat%PG%Rank) THEN
       old_value = mat_get_value(mat,global_row,global_col)
       CALL set_value(mat,global_row,global_col,old_value+value)
    END IF
    DEBEND("add_value(SFMO)")
    
  END SUBROUTINE sbm_add_value

  FUNCTION sbm_get_value(mat,global_row,global_col,trans_in) RESULT(value)
    TYPE(StoreBandedMatrixObject) :: mat
    INTEGER :: global_row, global_col
    complex :: value
    character(len=1),optional,intent(IN) :: trans_in

    ! local variables
    INTEGER :: proc, row_index, col_index, distance
    character(len=1) :: trans

    !CALL exit_if_factored(mat)

    if (present(trans_in)) then
       trans=trans_in
    else
       trans='N'
    end if

    ! calculate the processor coordinates, just determined by the global column index
    IF (mat%isTransposed) THEN
       if (trans.eq.'N') then
          proc = (global_row-1)/mat%ColsPerBlock
       else
          proc = (global_col-1)/mat%ColsPerBlock
       end if
    ELSE
       if (trans.eq.'N') then
          proc = (global_col-1)/mat%ColsPerBlock
       else
          proc = (global_row-1)/mat%ColsPerBlock
       end if
    END IF

    ! only process the entry if we are on the right processor
    IF (proc.EQ.mat%PG%Rank) THEN
       ! check, if the value is inside the band 
       if (trans.eq.'N') then
          distance = global_col-global_row
       else
          ! cols and rows are exchanged
          distance = global_row-global_col
       end if
       IF ((distance.le.mat%UpperBandwidth).and.(distance.ge.-mat%LowerBandwidth)) THEN
          IF (mat%isTransposed) THEN
             if (trans.eq.'N') then
                col_index = MOD(global_row-1,mat%ColsPerBlock)+1
             else
                col_index = MOD(global_col-1,mat%ColsPerBlock)+1
             end if
             row_index = distance
          ELSE
             if (trans.eq.'N') then
                col_index = MOD(global_col-1,mat%ColsPerBlock)+1
             else
                col_index = MOD(global_row-1,mat%ColsPerBlock)+1
             end if
             row_index = -distance
          END IF
          
          if (trans.ne.'C') then
             value = mat%localMatrix(row_index, col_index)
          else
             value = conjg(mat%localMatrix(row_index, col_index))
          end if
       ELSE
          ! out of bands
          ! return zero
          value = CMPLX(0.,0.,kind(value))
       END IF
    ELSE
       PRINT*,"This line should never be reached in mat_get_value of storebandedmatrix.F90."
       !CALL tracebackqq()
       value = CMPLX(1000, 70000)
       ! this value can be arbitrary, as it should never be used in further calculations
    END IF
  END FUNCTION sbm_get_value

  subroutine sbm_get_row_pointer(mat,global_row,ptr_row,start_col,end_col)
    TYPE(StoreBandedMatrixObject) :: mat
    INTEGER,intent(IN) :: global_row
    complex,dimension(:),pointer,intent(OUT) :: ptr_row
    integer, INTENT(OUT) :: start_col, end_col

    ! local variables
    INTEGER :: proc, local_row_index

    !CALL exit_if_factored(mat)

    if (.not.mat%isTransposed) then
       print*,"get_row_pointer only works for transposed matrices in storebandedmatrix"
       stop
    end if

    ! mat is transposed, stored in row-major order
    ! calculate the processor coordinates, just determined by the global column index
    proc = (global_row-1)/mat%ColsPerBlock

    ! only process the entry if we are on the right processor
    IF (proc.EQ.mat%PG%Rank) THEN
       local_row_index = mod(global_row-1,mat%ColsPerBlock)+1
       start_col = max(global_row - mat%LowerBandwidth,1)
       end_col   = min(global_row + mat%UpperBandwidth,mat%NCols)
       
       ptr_row  => mat%localMatrix(start_col-global_row:end_col-global_row,local_row_index)
    else
       nullify(ptr_row)
    end IF
    
  END subroutine sbm_get_row_pointer

  ! Here the cached values, which lie outside the bands are inserted.
  ! This is done, by extending the bandwidth to contain these entries. 
  SUBROUTINE sbm_commit_values(mat)
    TYPE(StoreBandedMatrixObject) :: mat

    ! Local variables
    TYPE(MatrixEntryNodeData), pointer :: nodedata_ptr
    TYPE(MatrixEntryNodeData) :: mynode
    COMPLEX, DIMENSION(:,:), POINTER :: newMatrix
    INTEGER, DIMENSION(2) :: max_bw, global_max_bw, coord
    INTEGER :: irow, icol, ierr, dist, row_index, col_index, old_NumberOfStoredRows

    CALL exit_if_factored(mat)
    
    ! 1. find maximal bandwidths in local list
    max_bw(1) = mat%UpperBandwidth
    max_bw(2) = mat%LowerBandwidth
    !PRINT*, "max_bw = ", max_bw

    call init_Iterator(mat%storage_list)
    do while (.not.isAtEnd(mat%storage_list))
       nodedata_ptr => get_next(mat%storage_list)
       coord = nodedata_ptr%coord
       dist = coord(2)-coord(1)
       IF (dist.GT.0) THEN
          max_bw(1) = MAX(dist, max_bw(1))
       ELSEIF (dist.LT.0) THEN
          max_bw(2) = MAX(-dist,max_bw(2))
       END IF
    END DO

    ! 2. we do a global reduction to get the global maximum
    CALL mpi_allreduce(max_bw, global_max_bw, 2, MPI_INTEGER, MPI_MAX, mat%PG%Communicator, ierr)

    ! 3. allocate the new local matrix if necessary
    IF ((global_max_bw(1).GT.mat%UpperBandwidth).OR.(global_max_bw(2).GT.mat%LowerBandwidth)) THEN
       ! Attention, newMatrix will not be deallocated, as the mat%localMatrix will point to it!
       old_NumberOfStoredRows = mat%NumberOfStoredRows
       mat%NumberOfStoredRows = global_max_bw(1)+global_max_bw(2)+1
       !PRINT*, "Allocating newMatrix(",mat%NumberOfStoredRows,mat%ColsPerBlock,")"
       IF (mat%isTransposed) THEN
          ALLOCATE(newMatrix(-global_max_bw(2):global_max_bw(1),mat%ColsPerBlock))
       ELSE
          ALLOCATE(newMatrix(-global_max_bw(1):global_max_bw(2),mat%ColsPerBlock))
       END IF
       newMatrix = CMPLX(0,0,kind(newMatrix))

       ! 4. copy old matrix into new one
       IF (mat%isTransposed) THEN
          DO icol=1,mat%ColsPerBlock
             DO irow=-mat%LowerBandwidth,mat%UpperBandwidth
                !PRINT*,irow,icol," row_index = ", global_max_bw(1)+1-mat%LowerBandwidth-1+irow
                newMatrix(irow,icol) = mat%localMatrix(irow,icol)
             END DO
          END DO
       ELSE
          DO icol=1,mat%ColsPerBlock
             DO irow=-mat%UpperBandwidth,mat%LowerBandwidth
                !PRINT*,irow,icol," row_index = ", global_max_bw(1)+1-mat%LowerBandwidth-1+irow
                newMatrix(irow,icol) = mat%localMatrix(irow,icol)
             END DO
          END DO
       END IF
       DEALLOCATE(mat%localMatrix)
       mat%localMatrix => newMatrix
       mat%UpperBandwidth = global_max_bw(1)
       mat%LowerBandwidth = global_max_bw(2)

       ! 5. add the entries of the storage_list
       DO WHILE (.NOT.(IsEmpty(mat%storage_list))) 
          CALL pop(mat%storage_list, mynode)

          IF (mat%isTransposed) THEN
             col_index = MOD(mynode%coord(1)-1,mat%ColsPerBlock)+1
             row_index = (mynode%coord(2)-mynode%coord(1))
          ELSE
             col_index = MOD(mynode%coord(2)-1,mat%ColsPerBlock)+1
             row_index = -(mynode%coord(2)-mynode%coord(1))
          END IF
          
          mat%localMatrix(row_index, col_index) = mynode%value
       END DO
    END IF

  END SUBROUTINE sbm_commit_values

  SUBROUTINE sbm_set_zero(this)
    type(StoreBandedMatrixObject) :: this
    
    ! deallocate all local storage and set the diagonal to zero

    IF (ASSOCIATED(this%localMatrix))  DEALLOCATE(this%localMatrix)
    IF (ASSOCIATED(this%factoredBand)) DEALLOCATE(this%factoredBand)
    IF (ASSOCIATED(this%auxData))      DEALLOCATE(this%auxData)
    IF (ASSOCIATED(this%ipiv))         DEALLOCATE(this%ipiv)

    this%UpperBandwidth=0
    this%LowerBandwidth=0
    this%NumberOfStoredRows=1
    ! allocate only the diagonal
    ALLOCATE(this%localMatrix(0:0,this%ColsPerBlock))
    ! set diagonal to zero
    this%localMatrix = CMPLX(0,0, kind(this%localMatrix))
  END SUBROUTINE sbm_set_zero

  !> Tune the Matrix-Vector multiplication
  subroutine sbm_autotune(mat,vec,res)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    TYPE(StoreVectorObject),INTENT(IN) :: vec
    TYPE(StoreVectorObject),INTENT(INOUT) :: res

    ! Local variables
    integer :: ierr, min_index(1), iter,max_iter=1
    real(8) :: start_time,end_time
    real(8) :: measured_time(4)

    ! the first call is always slower, therefore we measure 
    ! the second call
    measured_time = 0.0D0
    do iter=1,max_iter
       MatVecMethod = USE_BLAS
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       start_time = MPI_Wtime()
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       end_time = MPI_Wtime()
       measured_time(1) = measured_time(1) + end_time-start_time
       
       MatVecMethod = F_INTRINSIC
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       start_time = MPI_Wtime()
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       end_time = MPI_Wtime()
       measured_time(2) = measured_time(2) + end_time-start_time
       
       MatVecMethod = BY_INDEX
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       start_time = MPI_Wtime()
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       end_time = MPI_Wtime()
       measured_time(3) = measured_time(3) + end_time-start_time
       
       MatVecMethod = BY_ROW
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       start_time = MPI_Wtime()
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       end_time = MPI_Wtime()
       measured_time(4) = measured_time(4) + end_time-start_time
    end do
    write(*,"(A)",advance='no') "For Sparse Matrix Vector Multiplication, we use "
    min_index = minloc(measured_time)
    select case (min_index(1))
    case(1) 
       MatVecMethod = USE_BLAS
       write(*,"(A)",advance='no') "USE_BLAS"
    case(2) 
       MatVecMethod = F_INTRINSIC
       write(*,"(A)",advance='no') "F_INTRINSIC"
    case(3) 
       MatVecMethod = BY_INDEX
       write(*,"(A)",advance='no') "BY_INDEX"
    case(4) 
       MatVecMethod = BY_ROW
       write(*,"(A)",advance='no') "BY_ROW"
    end select
    write(*,"(A)") " method."

  end subroutine sbm_autotune

  !> Transpose the storage format of matrix mat
  !! The mathematical structure of the matrix is not changed, only
  !! the storage format (row-major or column-major) is changed.
  !! \param mat the input matrix
  !! \param mat_t the output matrix in the other storage format
  subroutine sbm_transpose_storage(mat,mat_t)
    type(StoreBandedMatrixObject), intent(IN) ::mat
    type(StoreBandedMatrixObject), intent(INOUT) :: mat_t

    ! Local variables
    integer :: rank, icol, ierr,irow
    integer :: global_col, global_row
    complex,allocatable,dimension(:,:) :: local_total_mat, local_total_mat_t

    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)

    if (mat%isTransposed) then
       print*,"transpose_storage only works for normal to transposed change"
    else
       if (.not.mat_t%isTransposed) then
          print*,"the target matrix of transpose_storage must have transposed format."
          stop
       end if
       ! mat is not transposed (that means it is in column-major order)
       
       ! simple algorithm (memory intensive)
       call allocate(mat_t,mat%LowerBandwidth,mat%UpperBandwidth)
       allocate(local_total_mat(-mat%UpperBandwidth:mat%LowerBandwidth,mat%NCols))
       allocate(local_total_mat_t(-mat_t%LowerBandwidth:mat_t%UpperBandwidth,mat_t%NRows))
       ! first gather the full matrix on all processes

       call mpi_gather(mat%localmatrix,size(mat%localmatrix),MPI_COMPLEX_TYPE,&
            &local_total_mat,size(mat%localmatrix),MPI_COMPLEX_TYPE,0,mat%PG%communicator,ierr)
       if (rank.eq.0) then
          do icol=1,mat_t%NRows
             do irow=-mat_t%LowerBandwidth,mat_t%UpperBandwidth
                global_row = icol
                global_col = icol+irow
                if ((global_col.ge.1).and.(global_col.le.mat_t%NRows)) then
                   local_total_mat_t(irow,icol) = local_total_mat(global_row-global_col,global_col)
                end if
             end do
          end do
       end if
       call mpi_scatter(local_total_mat_t,size(mat_t%localmatrix),MPI_COMPLEX_TYPE,&
            &mat_t%localmatrix,size(mat_t%localmatrix),MPI_COMPLEX_TYPE,0,mat%PG%communicator,ierr)
       deallocate(local_total_mat,local_total_mat_t)
    end if
  end subroutine sbm_transpose_storage

  !> Multiply BandedMatrix (in column-major storage) with Vector
  !! A StoreBandedMatrixObject in column-major (standard) storage format
  !! is multiplied with a StoreVectorObject. The result is a StoreVectorObject.
  !! \param mat the matrix
  !! \param vec the vector
  !! \param res the result vector, res = mat . vec
  SUBROUTINE sbm_dot_multiply(mat, vec, res)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    TYPE(StoreVectorObject),INTENT(IN) :: vec
    TYPE(StoreVectorObject),INTENT(INOUT) :: res

    !Local variables
    COMPLEX,PARAMETER :: complex_zero=(0,0)
    COMPLEX, DIMENSION(1:mat%NumberOfStoredRows+mat%ColsPerBlock-1) :: temp_res
    COMPLEX, DIMENSION(1:mat%UpperBandwidth) :: from_next
    COMPLEX, DIMENSION(1:mat%LowerBandwidth) :: from_previous
    INTEGER :: rank, ierr, pe_dest, pe_source, tag, t_innerlast,t_innerfirst, innerlast, icol
   ! INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: n_procs, npe_upper_sum, npe_lower_sum
    INTEGER, DIMENSION(1) :: dims, coords
    LOGICAL, DIMENSION(1) :: periods

    CALL exit_if_factored(mat)

    IF (mat%isTransposed) THEN
       ! A similar routine handles the case with transposed matrix
       CALL sbm_dot_multiply_transposed(mat,vec,res)
    ELSE
       CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)
       CALL mpi_cart_get(mat%PG%communicator, 1, dims, periods, coords,ierr)
       CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)

       temp_res = CMPLX(0.0,0.0,KIND(complex_kind))

       DO icol=1,mat%ColsPerBlock
          temp_res(icol:icol+mat%NumberOfStoredRows-1) = &
               & temp_res(icol:icol+mat%NumberOfStoredRows-1)&
               & +vec%localVector(icol)*mat%localMatrix(:,icol)
       END DO

       ! now we need a reduction over the processes

       ! the result vector has locally the following parts
       ! 1:UpperBandwidth, parts of the summation, which have to be sent to the previous process
       ! UpperBandwidth+1:UpperBandwidth+RowsPerBlock, the partial sums, which remain on the process
       ! UpperBandwidth+RowsPerBlock+1:UpperBandwidth+RowsPerBlock+LowerBandwidth, parts of the
       !  summation, which have to be sent to the next process

       ! innerlast is the index of the last inner point in the local res%localMatrix array
       innerlast = res%RowsPerBlock
       ! t_innerlast is the index of the last inner point in the local temp_res array
       t_innerfirst = mat%UpperBandwidth+1
       t_innerlast = mat%UpperBandwidth+res%RowsPerBlock

       IF (n_procs.GT.1) THEN
          ! exchange of the upper bandwidth
          IF (.NOT.periods(1)) THEN
             ! direction is not periodic
             ! calculate how many processes are involved in the sum-exchange
             npe_upper_sum = mat%UpperBandwidth/res%RowsPerBlock
             IF (MODULO(mat%UpperBandwidth,res%RowsPerBlock).NE.0) npe_upper_sum = npe_upper_sum + 1
             npe_lower_sum = mat%LowerBandwidth/res%RowsPerBlock
             IF (MODULO(mat%LowerBandwidth,res%RowsPerBlock).NE.0) npe_lower_sum = npe_lower_sum + 1

             IF (mat%UpperBandwidth.GT.0) THEN
                IF (npe_upper_sum.EQ.1) THEN
                   ! just one process involed in communication, we can do all sends and receives
                   ! in parallel
                   res%localVector = complex_zero
                   tag = 396
                   CALL mpi_cart_shift(mat%PG%communicator,0,-1,pe_source,pe_dest,ierr)
                   from_next = CMPLX(0,0,kind(from_next))
                   CALL mpi_sendrecv(temp_res,mat%UpperBandwidth,MPI_COMPLEX_TYPE,pe_dest,tag,&
                        & from_next,mat%UpperBandwidth, &
                        & MPI_COMPLEX_TYPE, pe_source, tag, &
                        & mat%PG%communicator, MPI_STATUS_IGNORE, ierr)
                   temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) = &
                        &temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) &
                        &+ from_next
                ELSE
                   ! we start sending at the last process and send always to the previous, then
                   ! add the result and send it to the previous, ...
                   tag = 339
                   IF (rank.EQ.n_procs-1) THEN
                      ! last process, only send
                      CALL mpi_send(temp_res,mat%UpperBandwidth,MPI_COMPLEX_TYPE, &
                           & rank-1, tag, mat%PG%communicator,ierr)
                   ELSEIF (rank.eq.0) then
                      ! first process, only receive, add
                      CALL mpi_recv(from_next,mat%UpperBandwidth,MPI_COMPLEX_TYPE, &
                           & rank+1, tag, mat%PG%communicator,MPI_STATUS_IGNORE,ierr)
                      temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) = &
                           &temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) &
                           &+ from_next
                   ELSE
                      ! all other processes, receive, add, send
                      call mpi_recv(from_next,mat%UpperBandwidth,MPI_COMPLEX_TYPE, &
                           & rank+1, tag, mat%PG%communicator,MPI_STATUS_IGNORE,ierr)
                      temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) = &
                           &temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) &
                           &+ from_next
                      call mpi_send(temp_res,mat%UpperBandwidth,MPI_COMPLEX_TYPE, &
                           & rank-1, tag, mat%PG%communicator,ierr)
                   END IF
                END IF
             END IF

             ! LowerBandwidth
             IF (mat%LowerBandwidth.GT.0) THEN
                IF (npe_lower_sum.EQ.1) THEN
                   tag = 407
                   CALL mpi_cart_shift(mat%PG%communicator,0,1,pe_source,pe_dest,ierr)
                   from_previous = CMPLX(0,0,kind(from_previous))
                   CALL mpi_sendrecv(temp_res(t_innerlast+1),mat%LowerBandwidth,MPI_COMPLEX_TYPE,pe_dest,tag,&
                        & from_previous,mat%LowerBandwidth,MPI_COMPLEX_TYPE, pe_source, tag, &
                        & mat%PG%communicator, MPI_STATUS_IGNORE, ierr)
                   temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) = &
                        &temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) &
                        &+ from_previous
                ELSE
                   ! we start sending at the first process and send always to the next, then
                   ! add the result and send it to the next, ...
                   tag = 366
                   IF (rank.EQ.0) THEN
                      ! first process, only send
                      CALL mpi_send(temp_res(t_innerlast+1),&
                           & mat%LowerBandwidth,MPI_COMPLEX_TYPE, &
                           & rank+1, tag, mat%PG%communicator,ierr)
                   ELSEIF (rank.eq.n_procs-1) then
                      ! last process, only receive, add
                      CALL mpi_recv(from_previous,mat%LowerBandwidth,MPI_COMPLEX_TYPE, &
                           & rank-1, tag, mat%PG%communicator,MPI_STATUS_IGNORE,ierr)
                      temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) = &
                           &temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) &
                           &+ from_previous
                   ELSE
                      ! all other processes, receive, add, send
                      call mpi_recv(from_previous,mat%LowerBandwidth,MPI_COMPLEX_TYPE, &
                           & rank-1, tag, mat%PG%communicator,MPI_STATUS_IGNORE,ierr)
                      temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) = &
                           &temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) &
                           &+ from_previous
                      CALL mpi_send(temp_res(t_innerlast+1),&
                           & mat%LowerBandwidth,MPI_COMPLEX_TYPE, &
                           & rank+1, tag, mat%PG%communicator,ierr)
                   END IF
                END IF
             END IF
          ELSE
             ! direction is periodic
             ! we cannot start at the last or first process immediately but we have to 
             ! sum up the entries at the start process
             STOP "periodic boundary condition for banded matrices not yet implemented"
          END IF
       END IF
       ! now temp_res contains the correct sum on all processes. The inner points have
       ! to be copied to the res%localMatrix
       res%localVector = temp_res(t_innerfirst:t_innerlast)
    END IF

  END SUBROUTINE sbm_dot_multiply

  !> Multiply BandedMatrix (in row-major storage) with Vector
  !! A StoreBandedMatrixObject in row-major (transposed) storage format
  !! is multiplied with a StoreVectorObject. The result is a StoreVectorObject.
  !! \param mat the matrix
  !! \param vec the vector
  !! \param res the result vector, res = mat . vec
  SUBROUTINE sbm_dot_multiply_transposed(mat,vec,res)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    TYPE(StoreVectorObject),INTENT(IN) :: vec
    TYPE(StoreVectorObject),INTENT(INOUT) :: res

    !Local variables
    INTEGER :: rank, n_procs, ierr,irow,row_start, mat_col_start, mat_col_end, global_row, i
    integer :: s_global_col, e_global_col, start_row, end_row, LBW, UBW
    !INTEGER, DIMENSION(1) :: dims, coords
    !LOGICAL, DIMENSION(1) :: periods
    COMPLEX, DIMENSION(1:vec%NRows),target :: localfullvec
    COMPLEX, DIMENSION(:), pointer :: ptr_vec, row_ptr
    complex :: prod, cdotu

    CALL exit_if_factored(mat)

    if (mat%LowerBandwidth+mat%UpperBandwidth+1.gt.mat%NRows) then
       MatVecMethod  = F_INTRINSIC
    else
       MatVecMethod = BY_INDEX
    end if

    LBW = mat%LowerBandwidth
    UBW = mat%UpperBandwidth
    ! We know that the matrix is transposed. 
    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)
    !CALL mpi_cart_get(mat%PG%communicator, 1, dims, periods, coords,ierr)
    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)
    
    ! First we gather the vector vec on all ranks.
    IF (n_procs.EQ.1) THEN
       ptr_vec => vec%localVector
    ELSE
       CALL get_global_vector_locally(vec,localfullvec)
       ptr_vec => localfullvec
    END IF
    ! localfullvec contains the full vector vec. This is now multiplied with
    ! the matrix

    !PERFON('sbm_dm')
    select case (MatVecMethod) 
    case(USE_BLAS)
       DO irow=1,mat%ColsPerBlock
          global_row = rank*mat%ColsPerBlock+irow
          call mat_get_row_pointer(mat,global_row,row_ptr,s_global_col,e_global_col)
          prod = cdotu(e_global_col-s_global_col+1,row_ptr,1,ptr_vec(s_global_col),1)
          
          CALL set_value(res,global_row, prod)
       end DO

    case(F_INTRINSIC)
       DO irow=1,mat%ColsPerBlock
          global_row = rank*mat%ColsPerBlock+irow
          call mat_get_row_pointer(mat,global_row,row_ptr,s_global_col,e_global_col)
          prod = dot_product(conjg(row_ptr),ptr_vec(s_global_col:e_global_col))
          CALL set_value(res,global_row, prod)
       end DO

    case(BY_INDEX)
       start_row = rank*mat%ColsPerBlock+1
       end_row   = (rank+1)*mat%ColsPerBlock
       
       if (start_row.le.mat%LowerBandwidth) then
          if (end_row.le.vec%NRows-mat%UpperBandwidth) then
             ! case 1
             do global_row = start_row, min(mat%LowerBandwidth,end_row)
                prod = mat%localMatrix(1-global_row,global_row-(start_row-1))*ptr_vec(1)
                DO i=1,mat%UpperBandwidth-1+global_row
                   prod = prod + mat%localMatrix(1-global_row+i,global_row-(start_row-1))*ptr_vec(1+i)
                END DO
                CALL set_value(res,global_row, prod)
             end do
             do global_row = min(mat%LowerBandwidth,end_row)+1,end_row
                mat_col_start = -mat%LowerBandwidth
                row_start = global_row - mat%LowerBandwidth
                mat_col_end = min(vec%NRows-global_row,mat%UpperBandwidth)
                
                prod = mat%localMatrix(-mat%LowerBandwidth,global_row-(start_row-1))*ptr_vec(row_start)
                DO i=1,mat_col_end-mat_col_start
                   prod = prod + mat%localMatrix(-mat%LowerBandwidth+i,global_row-(start_row-1))*ptr_vec(row_start+i)
                END DO
                CALL set_value(res,global_row, prod)
             end do
          else
             ! case 2
             ! This is mainly the case for n_procs = 1
             DO global_row = start_row,LBW
                !mat_col_start = 1-global_row
                !row_start = 1
                !mat_col_end = mat%UpperBandwidth
                
                prod = mat%localMatrix(1-global_row,global_row-start_row+1)*ptr_vec(1)
                DO i=1,UBW-(1-global_row)
                   prod = prod + mat%localMatrix(1-global_row+i,global_row-start_row+1)*ptr_vec(1+i)
                END DO
                CALL set_value(res,global_row, prod)
             end DO

             DO global_row = LBW+1,mat%NRows-UBW
                !mat_col_start = -mat%LowerBandwidth
                !row_start = global_row-mat%LowerBandwidth
                !mat_col_end = mat%UpperBandwidth

                !res%localvector(global_row) = dot_product(&
                !res%localvector(global_row-(start_row-1)) = dot_product(&
                !     & conjg(mat%localMatrix(-LBW:UBW,global_row-(start_row-1))),&
                !     & ptr_vec(global_row-LBW:global_row+UBW))

                prod = mat%localMatrix(-mat%LowerBandwidth,global_row-(start_row-1)) &
                     & *ptr_vec(global_row-mat%LowerBandwidth)
                DO i=1,mat%UpperBandwidth+mat%LowerBandwidth
                   prod = prod + mat%localMatrix(-mat%LowerBandwidth+i,global_row-(start_row-1)) &
                        & *ptr_vec(global_row-mat%LowerBandwidth+i)
                END DO
                res%localvector(global_row-(start_row-1)) = prod
                !CALL set_value(res,global_row, prod)
             end DO

             DO global_row = mat%NRows-UBW+1,end_row
                !mat_col_start = -mat%LowerBandwidth
                row_start = global_row-LBW
                !mat_col_end = vec%NRows-global_row
                
                prod = mat%localMatrix(-LBW,global_row-(start_row-1))*ptr_vec(row_start)
                DO i=1,vec%NRows-global_row+LBW
                   prod = prod + mat%localMatrix(-LBW+i,global_row-(start_row-1))*ptr_vec(row_start+i)
                END DO
                CALL set_value(res,global_row, prod)
             end DO
          end if
       else
          if (end_row.le.vec%NRows-mat%UpperBandwidth) then
             ! case 3
             DO global_row = start_row,end_row
                mat_col_start = max(1-global_row,-mat%LowerBandwidth)
                row_start = global_row+mat_col_start
                
                mat_col_end = min(vec%NRows-global_row,mat%UpperBandwidth)
                
                prod = mat%localMatrix(mat_col_start,global_row-(start_row-1))*ptr_vec(row_start)
                DO i=1,mat_col_end-mat_col_start
                   prod = prod + mat%localMatrix(mat_col_start+i,global_row-(start_row-1))*ptr_vec(row_start+i)
                END DO
                CALL set_value(res,global_row, prod)
             end DO
          else
             ! case 4
             DO global_row = start_row,min(mat%NRows-mat%UpperBandwidth,end_row)
                mat_col_start = max(1-global_row,-mat%LowerBandwidth)
                row_start = global_row+mat_col_start
                
                mat_col_end = min(vec%NRows-global_row,mat%UpperBandwidth)
                
                prod = mat%localMatrix(mat_col_start,global_row-(start_row-1))*ptr_vec(row_start)
                DO i=1,mat_col_end-mat_col_start
                   prod = prod + mat%localMatrix(mat_col_start+i,global_row-(start_row-1))*ptr_vec(row_start+i)
                END DO
                CALL set_value(res,global_row, prod)
             end DO
             DO global_row = max(start_row,mat%NRows-mat%UpperBandwidth),end_row
                mat_col_start = max(1-global_row,-mat%LowerBandwidth)
                row_start = global_row+mat_col_start
                
                mat_col_end = min(vec%NRows-global_row,mat%UpperBandwidth)
                
                prod = mat%localMatrix(mat_col_start,global_row-(start_row-1))*ptr_vec(row_start)
                DO i=1,mat_col_end-mat_col_start
                   prod = prod + mat%localMatrix(mat_col_start+i,global_row-(start_row-1))*ptr_vec(row_start+i)
                END DO
                CALL set_value(res,global_row, prod)
             end DO
          end if
       end if

    case(BY_ROW)
       DO irow=1,mat%ColsPerBlock
          global_row = rank*mat%ColsPerBlock+irow
          call mat_get_row_pointer(mat,global_row,row_ptr,s_global_col,e_global_col)
          prod = cmplx(0.0,0.0)
          do i=s_global_col,e_global_col
             prod = prod + row_ptr(i-s_global_col+1)*ptr_vec(i)
          end do
          CALL set_value(res,global_row, prod)
       end DO
    end select
    !PERFOFF
  END SUBROUTINE sbm_dot_multiply_transposed

  subroutine sbm_matmat_dot_multiply(mat,mat2,res,transA_in, transB_in)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat, mat2
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in     
    CHARACTER, intent(IN), optional :: transB_in     

    !Local variables
    CHARACTER :: transA='N', transB='N'

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif
    if (present(transB_in)) then
       transB = transB_in
    else
       transB = 'N'
    endif

    CALL exit_if_factored(mat)
    CALL exit_if_factored(mat2)
    
    ! We assume at the moment, that the result is always in row-major order.
    !if (.not.res%isTransposed) then
    !   print*, "In all dot_multiply functions of storebandedmatrix, the result must be in transposed storage."
    !   stop
    !end if

    if ( res%isTransposed ) then
       if ( ((mat%isTransposed).and.(transA.eq.'N')) .or. &
            & ( (.not.mat%isTransposed).and.((transA.eq.'C').or.(transA.eq.'T')) ) ) then
          ! we use the transposed structure of the storage of mat
          if ( ((mat2%isTransposed).and.(transB.eq.'N')) .or. &
               & ( (.not.mat2%isTransposed).and.((transB.eq.'C').or.(transB.eq.'T')) ) ) then
             ! do transp*transp -> transp
             call sbm_matmat_TTT_dot_multiply(mat,mat2,res,transA,transB)
          else
             ! do transp*normal -> transp
             call sbm_matmat_TNT_dot_multiply(mat,mat2,res,transA,transB)
          end if
       else
          ! here the storage is used in column-major order either with
          ! .not.mat%isTransposed (and transA.eq.'N') or with a combination of 
          ! mat%isTransposed and transA.eq.'C' or 'T'
          if ( ((mat2%isTransposed).and.(transB.eq.'N')) .or. &
               & ( (.not.mat2%isTransposed).and.((transB.eq.'C').or.(transB.eq.'T')) ) ) then
             ! we use the transposed structure of the storage of mat2
             ! do normal*transp -> transp
             call sbm_matmat_NTT_dot_multiply(mat,mat2,res,transA,transB)
          else
             ! here the storage of mat2 is used in column-major order either with
             ! .not.mat2%isTransposed (and transB.eq.'N') or with a combination of 
             ! mat2%isTransposed and transB.eq.'C' or 'T'
             ! do normal*normal -> transp
             call sbm_matmat_NNT_dot_multiply(mat,mat2,res,transA,transB)
          end if
       end if
    else
       if ( ((mat%isTransposed).and.(transA.eq.'N')) .or. &
            & ( (.not.mat%isTransposed).and.((transA.eq.'C').or.(transA.eq.'T')) ) ) then
          ! we use the transposed structure of the storage of mat
          if ( ((mat2%isTransposed).and.(transB.eq.'N')) .or. &
               & ( (.not.mat2%isTransposed).and.((transB.eq.'C').or.(transB.eq.'T')) ) ) then
             ! do transp*transp -> normal
             call sbm_matmat_TTN_dot_multiply(mat,mat2,res,transA,transB)
          else
             ! do transp*normal -> normal
             call sbm_matmat_TNN_dot_multiply(mat,mat2,res,transA,transB)
          end if
       else
          ! here the storage is used in column-major order either with
          ! .not.mat%isTransposed (and transA.eq.'N') or with a combination of 
          ! mat%isTransposed and transA.eq.'C' or 'T'
          if ( ((mat2%isTransposed).and.(transB.eq.'N')) .or. &
               & ( (.not.mat2%isTransposed).and.((transB.eq.'C').or.(transB.eq.'T')) ) ) then
             ! we use the transposed structure of the storage of mat2
             ! do normal*transp -> normal
             call sbm_matmat_NTN_dot_multiply(mat,mat2,res,transA,transB)
          else
             ! here the storage of mat2 is used in column-major order either with
             ! .not.mat2%isTransposed (and transB.eq.'N') or with a combination of 
             ! mat2%isTransposed and transB.eq.'C' or 'T'
             ! do normal*normal -> normal
             call sbm_matmat_NNN_dot_multiply(mat,mat2,res,transA,transB)
          end if
       end if
    end if
  end subroutine sbm_matmat_dot_multiply

  SUBROUTINE sbm_matmat_TNT_dot_multiply(mat, mat2, res,transA, transB)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat, mat2
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    ! Local variables
    INTEGER :: rank, ierr, n_procs
    !INTEGER, DIMENSION(1) :: dims, coords
    !LOGICAL, dimension(1) :: periods
    complex :: res_value
    INTEGER :: res_col_start, res_col_end, k, icol, irow, res_upperbandwidth, res_lowerbandwidth
    INTEGER :: global_row, k_start, k_end, res_col
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: localfullmat

    if (transA.ne.'N') then
       print*,"TNT: First matrix with 'C' or 'T' and transp not implemented yet."
       stop
    end if

    if (transB.ne.'N') then
       print*,"TNT: Second matrix with 'C' or 'T' and not_transp not implemented yet."
       stop
    end if

    ! first delete the result matrix. If this is not done, and the matrix has been used before
    ! for example in a loop with larger bandwidth than the result in the actual routine will
    ! use, then we get wrong results.
    call set_zero(res)
    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)
    !    CALL mpi_cart_get(mat%PG%communicator, 1, dims, periods, coords,ierr)
    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)

    ! 1. calculate the upper and lower bandwidth of the result matrix
    res_upperbandwidth = mat%UpperBandwidth+mat2%UpperBandwidth
    res_lowerbandwidth = mat%LowerBandwidth+mat2%LowerBandwidth

    IF (n_procs.GT.1) THEN
       
       ! simplified transfer: copy whole mat2 to all ranks
       ALLOCATE(localfullmat(-mat2%UpperBandwidth:mat2%LowerBandwidth,mat2%NCols))
       CALL mpi_allgather(mat2%localMatrix,SIZE(mat2%localMatrix),MPI_COMPLEX_TYPE,&
            & localfullmat,SIZE(mat2%localMatrix),MPI_COMPLEX_TYPE,&
            & mat2%PG%communicator,ierr)

       ! now multiply the two matrices
       DO irow=1,mat%ColsPerBlock
          global_row = rank*mat%ColsPerBlock+irow

          res_col_start = MAX(1,global_row-res_lowerbandwidth)
          res_col_end   = MIN(res%NCols,global_row+res_upperbandwidth)

          DO icol=res_col_start,res_col_end
             k_start = MAX(global_row-mat%LowerBandwidth,icol-mat2%UpperBandwidth,1)
             k_end   = MIN(global_row+mat%UpperBandwidth,icol+mat2%LowerBandwidth,mat%NCols)
             res_value = CMPLX(0,0,KIND(res_value))
             DO k=k_start, k_end
                !res_value = res_value + mat_get_value(mat,global_row,k)*localfullmat(k-icol,icol)
                res_value = res_value + mat%localMatrix(k-global_row,irow)*localfullmat(k-icol,icol)
             END DO
             CALL set_value(res,global_row,icol,res_value)
          END DO
       END DO
       CALL commit_values(res)
       DEALLOCATE(localfullmat)

       ! 2. determine which columns of mat2 are to be transferred
       ! last row on a rank
       !global_row = (rank+1)*mat%ColsPerBlock-1
       ! last column not equal zero in the result matrix is
       !max_col_index = global_row+res_upperbandwidth
       ! on which process do we have this column
       !col_on_proc = max_col_index/mat%ColsPerBlock
       ! so we need all columns from rank+1 to col_on_procs
       ! allocate the recv buffer
       !n_recv_cols = col_on_proc-(rank+1)+1)*mat%ColsPerBlock
       !ALLOCATE(recvbuf(mat2%UpperBandwidth:mat2%LowerBandwidth,mat2%ColsPerBlock))
       !DO iproc=1,col_on_proc
       !   CALL mpi_recv(recvbuf(:,)
          
       ! 3. start the transfer with non blocking send and receive

       ! 4. In parallel, compute the local entries of the result matrix

       ! 5. wait until transfer ready

       ! 6. then compute the missing entries of the result matrix
       
    ELSE
       ! matrices are not parallelized
       DO global_row = 1, res%NRows
          res_col_start = MAX(1,global_row-res_lowerbandwidth)
          res_col_end   = MIN(res%NCols,global_row+res_upperbandwidth)

          DO res_col = res_col_start,res_col_end
             k_start = MAX(global_row-mat%LowerBandwidth,res_col-mat2%UpperBandwidth,1)
             k_end   = MIN(global_row+mat%UpperBandwidth,res_col+mat2%LowerBandwidth,mat%NCols)

             res_value = CMPLX(0,0,kind(res_value))
             DO k=k_start,k_end
                !res_value = res_value + mat_get_value(mat,global_row,k)*mat_get_value(mat2,k,res_col)
                res_value = res_value + mat%localMatrix(k-global_row,global_row)*mat2%localMatrix(k-res_col,res_col)
                ! with direct access to localMatrix this could be accelerated.
                !res_value = res_value + mat_get_value(mat,global_row,k)*mat_get_value(mat2,k,res_col)
             END DO
             CALL set_value(res,global_row,res_col,res_value)
          END DO
       END DO
       call commit_values(res)
    END IF

  END SUBROUTINE sbm_matmat_TNT_dot_multiply

  !>Multiplication of a BandedMatrix with a BandedMatrix
  !!
  !! This is the routine where the real work is done. A StoreBandedMatrixObject
  !! is multiplied with another StoreBandedMatrixObject to give a StoreBandedMatrixObject.
  !! All BandedMatrix object have to be in transposed storage format (row-major order).
  !! \param mat the banded matrix
  !! \param mat2 the second banded matrix
  !! \param res  the banded result matrix
  !! \param transA_in operation on mat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  !! \param transB_in operation on mat2, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  subroutine sbm_matmat_TTT_dot_multiply(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    type(StoreBandedMatrixObject),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    !Local variables
    complex, dimension(:,:), pointer :: ptr_next_recv, ptr_next_compute, temp_ptr
    !complex, dimension(:), pointer :: ptr_fast_first, ptr_fast_second
    complex, dimension(:,:), allocatable, target :: temp1, temp2
    complex :: sum_value, cdotu
    integer :: ierr,send_request, recv_request, tag, n_procs, rank
    !integer :: status(MPI_STATUS_SIZE)
    integer :: i_pe,i_block,i_block_recv,i_block_send
    integer :: global_res_row, global_res_col, global_index_start,global_index_end
    integer :: s_colindexset_mat, e_colindexset_mat,s_rowindexset_mat2,e_rowindexset_mat2
    integer :: s_intersection, e_intersection
    integer :: s_global_col, e_global_col, loop_length
    integer :: RowsPerBlock, res_lbw, res_ubw, local_row,local_col
    complex, dimension(:),pointer :: first_mat_row
    !complex, dimension(1:res%NCols) :: result_row

    if ((transA.ne.'N').or.(transB.ne.'N')) then
       print*,"At the moment, sbm_matmat_TTT_dot_multiply only works with 'N' for both matrices."
       !call tracebackqq()
       stop
    end if

    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)
    call mpi_comm_rank(mat%PG%communicator, rank, ierr)


    ! allocate the two temporary arrays
    allocate(temp1(-mat2%LowerBandwidth:mat2%UpperBandwidth,mat2%ColsPerBlock))
    allocate(temp2(-mat2%LowerBandwidth:mat2%UpperBandwidth,mat2%ColsPerBlock))

    call set_zero(res)
    RowsPerBlock = mat%NRows/mat%PG%NProcs
    res_lbw = mat%LowerBandwidth+mat2%LowerBandwidth
    res_ubw = mat%UpperBandwidth+mat2%UpperBandwidth
    ! some arbitrary tag, just for debug purposes
    tag = 9784
    ! loop over all processes
    ! we start on all processes with the actual process
    ptr_next_compute => mat2%localMatrix
    ptr_next_recv => temp1
    do i_pe=0,n_procs-1
       ! calculate the running index, which starts at the rank
       ! and wraps around at n_procs
       ! i_block runs from 0 to n_procs-1, with wrapping at n_procs
       i_block = modulo(rank+i_pe,n_procs)
       
       i_block_recv = modulo(i_block+1,n_procs)
       i_block_send = modulo(rank-i_pe-1+n_procs,n_procs)

       !write(*,"(I3,A,4I3)") rank," : ",i_pe,i_block,i_block_recv,i_block_send

       ! do not wait for the first iteration, as nothing has been sent
       if (i_block.ne.rank) then
          ! wait for the received block
          !print*,rank,": Waiting for recv to finalize"
          call mpi_wait(recv_request,MPI_STATUS_IGNORE,ierr)
          ! swap the pointers
          temp_ptr => ptr_next_recv
          ptr_next_recv => ptr_next_compute
          ptr_next_compute => temp_ptr
          ! wait for the sended block
          !print*,rank,": Waiting for send to finalize"
          call mpi_wait(send_request,MPI_STATUS_IGNORE,ierr)
       end if
       ! initiate non-blocking receive for the next block from some following pe
       if (i_block_recv.ne.rank) then
          !print*,rank,": Starting irecv from ",i_block_recv
          call mpi_irecv(ptr_next_recv,mat2%NumberOfStoredRows*mat2%ColsPerBlock,MPI_COMPLEX_TYPE,&
               &i_block_recv,tag,mat2%PG%communicator,recv_request,ierr)
       end if
       ! initiate non-blocking send of the local array to some previous pe
       if (i_block_send.ne.rank) then
          !print*,rank,": Starting isend to ",i_block_send
          call mpi_isend(mat2%localMatrix,mat2%NumberOfStoredRows*mat2%ColsPerBlock,MPI_COMPLEX_TYPE,&
               &i_block_send,tag,mat2%PG%communicator,send_request,ierr)
       end if
       ! do computation
       ! multiply mat%localMatrix with ptr_next_compute
       ! mat is a transposed matrix, also the ptr_next_compute
       ! ptr_next_compute is a subblock of the localMatrix of mat2

       global_index_start = i_block*RowsPerBlock+1
       global_index_end   = (i_block+1)*RowsPerBlock
       !print*,rank,i_block,": global_index = ",global_index_start,global_index_end
       do global_res_row=rank*RowsPerBlock+1,(rank+1)*RowsPerBlock
          call mat_get_row_pointer(mat,global_res_row,first_mat_row,s_global_col,e_global_col)
          !print*,first_mat_row
          do global_res_col=max(global_res_row-res_lbw,1),min(global_res_row+res_ubw,res%NCols)
             !write(*,"(A,2I4)") "global_res = ", global_res_row,global_res_col

             ! calculating the global index set
             ! [s_global,e_global] is already intersection of colindexset_mat and I_block
             s_colindexset_mat = global_res_row-mat%LowerBandwidth
             e_colindexset_mat = global_res_row+mat%UpperBandwidth
             s_rowindexset_mat2= global_res_col-mat2%UpperBandwidth
             e_rowindexset_mat2= global_res_col+mat2%LowerBandwidth
             !print*,"colindexset = ",s_colindexset_mat,e_colindexset_mat
             !print*,"rowindexset = ",s_rowindexset_mat2,e_rowindexset_mat2

             ! now we need the intersection of I_block colindexset and rowindexset
             s_intersection = max(global_index_start,s_colindexset_mat,s_rowindexset_mat2)
             e_intersection = min(global_index_end, e_colindexset_mat,e_rowindexset_mat2)
             !write(*,"(4X,A,2I3)") "intersection = ", s_intersection,e_intersection
             loop_length = max(e_intersection-s_intersection+1,0)
             if (loop_length.gt.0) then

                ! the global index runs from s_intersection to e_intersection
                ! now convert to local indices
                !print*,first_mat_row(1+s_intersection-s_global_col:1+s_intersection-s_global_col+loop_length-1)
                local_row = lbound(ptr_next_compute,1)+global_res_col-s_intersection+mat2%LowerBandwidth
                local_col = mod(s_intersection-1,mat2%ColsPerBlock)+1
                !print*,"local_row,col = ",local_row,local_col
                !print*,ptr_next_compute(local_row,local_col)
                !loop_length = min(e_global_col,global_index_end)-max(s_global_col,global_index_start)+1
                !result_row(global_res_col) = cdotu(loop_length,first_mat_row(1+s_intersection-s_global),1,&
                sum_value = cdotu(loop_length,first_mat_row(1+s_intersection-s_global_col),1,&
                     & ptr_next_compute(local_row,local_col),mat2%UpperBandwidth+mat2%LowerBandwidth)
                call add_value(res,global_res_row,global_res_col,sum_value)
             end if

             !ptr_local_column => ptr_next_compute((global_mat_col-1)*RowsPerBlock+1:global_mat_col*RowsPerBlock)
             !if (s_global_col.ge.global_index_start) then
             !   ptr_fast_first => first_mat_row(1:loop_length)
             !   ptr_fast_second => ptr_local_column(s_global_col-global_index_start+1:s_global_col-global_index_start+loop_length)
             !else
             !   ptr_fast_first => first_mat_row(global_index_start-s_global_col+1:global_index_start-s_global_col+loop_length)
             !   ptr_fast_second => ptr_local_column(1:loop_length)
             !end if

             ! to save memory, we can write directly into the matrix with the following line
             !sum_value = cdotu(loop_length,ptr_fast_first,1,ptr_fast_second,1)
             !call add_value(res,global_bmat_row,global_mat_col,sum_value)
             ! for better performance, we collect a whole row of the result matrix
             !result_row(global_mat_col) = cdotu(loop_length,ptr_fast_first,1,ptr_fast_second,1)
          end do
          ! for better performance, now add the whole result_row to the result matrix
          !call add_to_row(res,global_bmat_row, result_row)
       end do
       call commit_values(res)
       ! ---------------------- Ende des komplizierten Teils --------------

       if (i_block.eq.rank) then
          ! the first iteration is special, as ptr_next_compute points
          ! to localMatrix, but to do the pointer swapping, we have to 
          ! redirect ptr_next_compute to the other temporary storage.
          ptr_next_compute => temp2
       end if
    end do
    deallocate(temp1,temp2)
  end subroutine sbm_matmat_TTT_dot_multiply

  !>Multiplication of a BandedMatrix with a BandedMatrix
  !!
  !! This is the routine where the real work is done. A StoreBandedMatrixObject
  !! is multiplied with another StoreBandedMatrixObject to give a StoreBandedMatrixObject.
  !! The first matrix is now in column-major order, the two other matrices are
  !! stored in row-major order.
  !! \param mat the banded matrix
  !! \param mat2 the second banded matrix
  !! \param res  the banded result matrix
  !! \param transA_in operation on mat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  !! \param transB_in operation on mat2, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  subroutine sbm_matmat_NTT_dot_multiply(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    type(StoreBandedMatrixObject),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    !Local variables
    complex, dimension(:), pointer :: res_row_ptr, row_ptr
    complex, dimension(:,:),pointer :: res_localmatrix_temp
    integer, dimension(:), allocatable :: recvcounts

    integer :: ierr, n_procs, rank
    integer :: first_LBW, first_UBW, sec_LBW, sec_UBW
    integer :: global_res_row, s_global_index,e_global_index
    integer :: s_colindexset_mat, e_colindexset_mat,s_colindexset_mat2,e_colindexset_mat2
    integer :: s_intersection, e_intersection, loop_length
    integer :: s_intersection2, e_intersection2, loop_length2
    integer :: s_global_col, e_global_col,s_res_col, e_res_col
    integer :: res_lbw, res_ubw, inner_index
    logical :: conjugate, conjugate2

    !integer :: zaxpy_calls

    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)
    call mpi_comm_rank(mat%PG%communicator, rank, ierr)

    ! test if we have to use a conjugated matrix or not
    if (transA.eq.'N') then
       ! first matrix is not Transposed
       conjugate=.false.
       first_LBW = mat%LowerBandwidth
       first_UBW = mat%UpperBandwidth
    else
       ! first matrix isTransposed, that means Lower and Upper
       ! Bandwidths are exchanged
       first_LBW = mat%UpperBandwidth
       first_UBW = mat%LowerBandwidth
       if (transA.eq.'C') then
          conjugate=.true.
       else
          conjugate=.false.
       end if
    end if

    if (transB.eq.'N') then
       ! second matrix isTransposed
       sec_LBW=mat2%LowerBandwidth
       sec_UBW=mat2%UpperBandwidth
       conjugate2 = .false.
    else
       ! second matrix .not.isTransposed
       sec_LBW = mat2%UpperBandwidth
       sec_UBW = mat2%LowerBandwidth
       if (transB.eq.'C') then
          conjugate2=.true.
       else
          conjugate2=.false.
       end if
    end if

    ! allocate the two temporary arrays

    res_lbw = min(first_LBW+sec_LBW,res%NCols)
    res_ubw = min(first_UBW+sec_UBW,res%NCols)
    !call set_zero(res)
    ! allocate the storage in res with the calculated bandwidths
    call allocate(res,res_lbw,res_ubw)
    if (n_procs.ne.1) then
       ! for communication we need a temporary global array
       ! if memory is an issue, this can be changed
       !write(*,"(A,I10,A,F8.2,A)") "Allocating res_localmatrix_temp with ",&
       !     &(res_ubw+res_lbw+1)*res%NRows," elements (",&
       !     &(res_ubw+res_lbw+1)*res%NRows*16/1024./1024.," MB)."
       !DEC$ ATTRIBUTES ALIGN: 128:: res_localmatrix_temp
       allocate(res_localmatrix_temp(-res_lbw:res_ubw,1:res%NRows))
       res_localmatrix_temp = cmplx(0.0,0.0)
    else
       res_localmatrix_temp => res%localMatrix
       !write(*,"(A,I10,A,F8.2,A)") "Using res_localmatrix_temp with ",size(res%localMatrix),&
       !     &" elements (",size(res%localMatrix)*16/1024./1024.," MB)."
    end if
    
    !zaxpy_calls = 0
    ! calculate the product on each processor with the localmatrices
    s_global_index = rank*mat%ColsPerBlock+1
    e_global_index = (rank+1)*mat%ColsPerBlock
    do global_res_row=1,res%NRows
       !write(*,"(A,I4)") "global_res_row = ",global_res_row

       ! calculating the global index set
       ! [s_global,e_global] is already intersection of colindexset_mat and I_block
       s_colindexset_mat = global_res_row-first_LBW
       e_colindexset_mat = global_res_row+first_UBW
       !print*,"colindexset = ",s_colindexset_mat,e_colindexset_mat
       !print*,"colindexset2 = ",s_colindexset_mat2,e_colindexset_mat2
       
       ! now we need the intersection of I_block colindexset and rowindexset
       s_intersection = max(s_global_index,s_colindexset_mat)
       e_intersection = min(e_global_index, e_colindexset_mat)
       !print*,"intersection = ",s_intersection,e_intersection
       !print*,"intersection2 = ",s_intersection2,e_intersection2

       loop_length = max(e_intersection-s_intersection+1,0)
       
       !if (loop_length2.gt.0) then
       ! ========= mat_get_row_pointer not working on res but on res_localmatrix_temp
       s_res_col = max(global_res_row - res_LBW,1)
       e_res_col = min(global_res_row + res_UBW,res%NCols)
       
       res_row_ptr  => res_localmatrix_temp(s_res_col-global_res_row:e_res_col-global_res_row,global_res_row)
       !write(*,"(A,2I3,A)") "res_row_ptr(",lbound(res_row_ptr,1),ubound(res_row_ptr,1),")"
       !print*,s_res_col-global_res_row,e_res_col-global_res_row
       !==========
       !call mat_get_row_pointer(res,global_res_row,res_row_ptr,s_res_col,e_res_col)
       !print*,"res_col = ",s_res_col,e_res_col,", res_row_ptr = ",res_row_ptr

       !zaxpy_calls = zaxpy_calls + loop_length

       do inner_index=s_intersection,e_intersection
          !write(*,"(5X,A,I4)") "inner_index = ",inner_index
          call mat_get_row_pointer(mat2,inner_index,row_ptr,s_global_col,e_global_col)
          !print*,"global_col = ",s_global_col,e_global_col,", row_ptr = ",row_ptr
          
          s_colindexset_mat2= inner_index-sec_LBW
          e_colindexset_mat2= inner_index+sec_UBW
          s_intersection2 = max(1,s_colindexset_mat2)
          e_intersection2 = min(mat2%NCols,e_colindexset_mat2)
          loop_length2 = max(e_intersection2-s_intersection2+1,0)

          !write(*,"(A,I4,A)") "zaxpy(",loop_length2,")"
          call caxpy(loop_length2,mat_get_value(mat,global_res_row,inner_index,transA),&
               &row_ptr,1,res_row_ptr(1+s_global_col-s_res_col),1)
       end do
    end do

    !write(*,"(A,I10,A)") "zaxpy was called ",zaxpy_calls," times."
    if (n_procs.ne.1) then
       ! Here now follows the communication of the result matrices. They have all to be reduced
       ! but different parts of the matrix to different processes.
       allocate(recvcounts(n_procs))
       recvcounts=res%NumberOfStoredRows*res%ColsPerBlock
       call mpi_reduce_scatter(res_localMatrix_temp,res%localMatrix,recvcounts,&
            & MPI_COMPLEX_TYPE,MPI_SUM,res%PG%Communicator,ierr)
       deallocate(recvcounts)
       deallocate(res_localmatrix_temp)
    end if
  end subroutine sbm_matmat_NTT_dot_multiply

  subroutine sbm_matmat_NNT_dot_multiply(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    type(StoreBandedMatrixObject),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    print*,"NNT_dot_multiply is not yet implemented."
    stop
  end subroutine sbm_matmat_NNT_dot_multiply

  SUBROUTINE sbm_matmat_TNN_dot_multiply(mat, mat2, res,transA, transB)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat, mat2
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    ! Local variables
    INTEGER :: rank, ierr, n_procs
    !INTEGER, DIMENSION(1) :: dims, coords
    !LOGICAL, dimension(1) :: periods
    complex :: res_value
    INTEGER :: res_col_start, res_col_end, k, icol, irow, res_upperbandwidth, res_lowerbandwidth
    integer :: res_row_start, res_row_end
    INTEGER :: global_row, global_col, k_start, k_end, res_col
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: localfullmat

    if (transA.ne.'N') then
       print*,"TNT: First matrix with 'C' or 'T' and transp not implemented yet."
       stop
    end if

    if (transB.ne.'N') then
       print*,"TNT: Second matrix with 'C' or 'T' and not_transp not implemented yet."
       stop
    end if

    ! first delete the result matrix. If this is not done, and the matrix has been used before
    ! for example in a loop with larger bandwidth than the result in the actual routine will
    ! use, then we get wrong results.
    call set_zero(res)
    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)
    !    CALL mpi_cart_get(mat%PG%communicator, 1, dims, periods, coords,ierr)
    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)

    ! 1. calculate the upper and lower bandwidth of the result matrix
    res_upperbandwidth = mat%UpperBandwidth+mat2%UpperBandwidth
    res_lowerbandwidth = mat%LowerBandwidth+mat2%LowerBandwidth

    IF (n_procs.GT.1) THEN
       
       ! simplified transfer: copy whole mat to all ranks
       ALLOCATE(localfullmat(-mat%LowerBandwidth:mat%UpperBandwidth,mat%NCols))
       CALL mpi_allgather(mat%localMatrix,SIZE(mat%localMatrix),MPI_COMPLEX_TYPE,&
            & localfullmat,SIZE(mat%localMatrix),MPI_COMPLEX_TYPE,&
            & mat%PG%communicator,ierr)


       ! each processor works on a column block of the whole matrix
       !do icol=1,res%ColsPerBlock
       !   global_col = rank*res%ColsPerBlock+icol

       !do global_row = 1,res%NRows
       !      res_value = CMPLX(0,0,KIND(res_value))
       !      do k=1,mat%NRows
                ! mat is stored in tranposed format
       !         res_value += mat(k,global_row)*mat2(k,icol)
       !      end do
       !      CALL set_value(res,global_row,global_col,res_value)
       !   end do
       !end do

       do icol=1,res%ColsPerBlock
          global_col = rank*res%ColsPerBlock+icol

          res_row_start = MAX(1,global_col-res_upperbandwidth)
          res_row_end   = MIN(res%NCols,global_col+res_lowerbandwidth)

          do global_row = res_row_start,res_row_end
             ! now determine the k range
             k_start = max(1,global_row-mat%LowerBandwidth,global_col-mat2%UpperBandwidth)
             k_end   = min(mat%NCols,global_row+mat%UpperBandwidth,global_col+mat2%LowerBandwidth)
             !write(*,"(4(A,I3))") "Calculating r(",global_row,",",global_col,&
             !     & "), k running from ",k_start," to ",k_end
             res_value = CMPLX(0,0,KIND(res_value))
             do k=k_start,k_end
                ! mat is stored in tranposed format
                res_value = res_value + &
                     &localfullmat(k-global_row,global_row)*mat_get_value(mat2,k,global_col)
             end do
             CALL set_value(res,global_row,global_col,res_value)
          end do
       end do

       CALL commit_values(res)
       DEALLOCATE(localfullmat)

       ! 2. determine which columns of mat2 are to be transferred
       ! last row on a rank
       !global_row = (rank+1)*mat%ColsPerBlock-1
       ! last column not equal zero in the result matrix is
       !max_col_index = global_row+res_upperbandwidth
       ! on which process do we have this column
       !col_on_proc = max_col_index/mat%ColsPerBlock
       ! so we need all columns from rank+1 to col_on_procs
       ! allocate the recv buffer
       !n_recv_cols = col_on_proc-(rank+1)+1)*mat%ColsPerBlock
       !ALLOCATE(recvbuf(mat2%UpperBandwidth:mat2%LowerBandwidth,mat2%ColsPerBlock))
       !DO iproc=1,col_on_proc
       !   CALL mpi_recv(recvbuf(:,)
          
       ! 3. start the transfer with non blocking send and receive

       ! 4. In parallel, compute the local entries of the result matrix

       ! 5. wait until transfer ready

       ! 6. then compute the missing entries of the result matrix
       
    ELSE
       ! matrices are not parallelized
       DO global_row = 1, res%NRows
          res_col_start = MAX(1,global_row-res_lowerbandwidth)
          res_col_end   = MIN(res%NCols,global_row+res_upperbandwidth)

          DO res_col = res_col_start,res_col_end
             k_start = MAX(global_row-mat%LowerBandwidth,res_col-mat2%UpperBandwidth,1)
             k_end   = MIN(global_row+mat%UpperBandwidth,res_col+mat2%LowerBandwidth,mat%NCols)

             res_value = CMPLX(0,0,kind(res_value))
             DO k=k_start,k_end
                !res_value = res_value + mat_get_value(mat,global_row,k)*mat_get_value(mat2,k,res_col)
                res_value = res_value + mat%localMatrix(k-global_row,global_row)*mat2%localMatrix(k-res_col,res_col)
                ! with direct access to localMatrix this could be accelerated.
                !res_value = res_value + mat_get_value(mat,global_row,k)*mat_get_value(mat2,k,res_col)
             END DO
             CALL set_value(res,global_row,res_col,res_value)
          END DO
       END DO
       call commit_values(res)
    END IF

  END SUBROUTINE sbm_matmat_TNN_dot_multiply

  !>Multiplication of a BandedMatrix with a BandedMatrix
  !!
  !! This is the routine where the real work is done. A StoreBandedMatrixObject
  !! is multiplied with another StoreBandedMatrixObject to give a StoreBandedMatrixObject.
  !! All BandedMatrix object have to be in transposed storage format (row-major order).
  !! \param mat the banded matrix
  !! \param mat2 the second banded matrix
  !! \param res  the banded result matrix
  !! \param transA_in operation on mat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  !! \param transB_in operation on mat2, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  subroutine sbm_matmat_TTN_dot_multiply(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    type(StoreBandedMatrixObject),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    print*,"TTN_dot_multiply is not yet implemented."
    stop

  end subroutine sbm_matmat_TTN_dot_multiply

  !>Multiplication of a BandedMatrix with a BandedMatrix
  !!
  !! This is the routine where the real work is done. A StoreBandedMatrixObject
  !! is multiplied with another StoreBandedMatrixObject to give a StoreBandedMatrixObject.
  !! The first matrix is now in column-major order, the two other matrices are
  !! stored in row-major order.
  !! \param mat the banded matrix
  !! \param mat2 the second banded matrix
  !! \param res  the banded result matrix
  !! \param transA_in operation on mat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  !! \param transB_in operation on mat2, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  subroutine sbm_matmat_NTN_dot_multiply(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    type(StoreBandedMatrixObject),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    print*,"NTN_dot_multiply is not yet implemented."
    stop

  end subroutine sbm_matmat_NTN_dot_multiply

  subroutine sbm_matmat_NNN_dot_multiply(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    type(StoreBandedMatrixObject),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    print*,"NNN_dot_multiply is not yet implemented."
    stop
  end subroutine sbm_matmat_NNN_dot_multiply

  !>Multiplication of a BandedMatrix with a Matrix
  !!
  !! This is the routine where the real work is done. A StoreBandedMatrixObject
  !! is multiplied with a StoreFullMatrixObject to give a StoreFullMatrixObject.
  !! The BandedMatrix has to be in transposed storage format, the full matrix
  !! in normal storage format (this is the only possibility at the moment). 
  !! The result will be in normal storage format.
  !! \param bmat the banded matrix
  !! \param mat  the full matrix
  !! \param res  the full result matrix
  !! \param transA_in operation on bmat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  !! \param transB_in operation on mat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  subroutine sbm_dot_multiply_Banded_with_Full(bmat,mat,res,transA_in,transB_in)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: bmat
    type(StoreFullMatrixObject),intent(IN) :: mat
    TYPE(StoreFullMatrixObject),INTENT(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in     
    CHARACTER, intent(IN), optional :: transB_in     

    !Local variables
    CHARACTER :: transA='N', transB='N'
    complex, dimension(:), pointer :: ptr_next_recv, ptr_next_compute, temp_ptr, ptr_local_column
    complex, dimension(:), pointer :: ptr_fast_first, ptr_fast_second
    complex, dimension(:), allocatable, target :: temp1, temp2
    complex, dimension(:),pointer :: first_mat_row
    complex, dimension(1:res%NCols) :: result_row
    complex :: cdotu
    integer :: ierr,send_request, recv_request, tag, n_procs, rank
    !integer :: status(MPI_STATUS_SIZE)
    integer :: i_pe,i_block,i_block_recv,i_block_send,global_bmat_row
    integer :: global_mat_col,global_index_start,global_index_end
    integer :: s_global_col, e_global_col, loop_length

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif
    if (present(transB_in)) then
       transB = transB_in
    else
       transB = 'N'
    endif

    if ((transA.ne.'N').or.(transB.ne.'N')) then
       print*,"At the moment, dot_multiply_Banded_with_Full only works with 'N' for both matrices."
       !call tracebackqq()
       stop
    end if

    ! if one of the matrices is LU decomposed, we skip, as the normal
    ! dot multiply does not work with LU decomposed matrices
    CALL exit_if_factored(bmat)
    CALL exit_if_factored(mat)

    IF (.NOT.(bmat%isTransposed)) then
       PRINT*,"dot_multiply for banded with full is only allowed with &
            &the storage schemes transp*not_transp ->not_transp"
       STOP
    END IF

    !call print_storage_details(mat)

    CALL mpi_comm_size(bmat%PG%communicator, n_procs, ierr)
    call mpi_comm_rank(bmat%PG%communicator, rank, ierr)


    if (mat%BlockrowsPerProcess.eq.1) then
       ! we have a distribution over the rows of the full matrix
       ! or over the columns of the full matrix, exactly one block
       ! per process.
       ! allocate the two temporary arrays
       allocate(temp1(mat%NBlocks*mat%Blocksize))
       allocate(temp2(mat%NBlocks*mat%Blocksize))
    else
       ! scalapack block distributed storage of the full matrix
       ! not implemented!
       print*,"dot_multiply for banded with full for a block distributed &
            &full matrix with NPCols.ne.1 .and. NPRows.ne.1 is not implemented."
       stop
    end if

    ! for test purposes, set res matrix to zero
    ! has to be changed for production, just for testing
    res%localMatrix=cmplx(0.0,0.0)

    ! some arbitrary tag, just for debug purposes
    tag = 9783
    ! loop over all processes
    ! we start on all processes with the actual process
    ptr_next_compute => mat%localMatrix
    ptr_next_recv => temp1
    do i_pe=0,n_procs-1
       ! calculate the running index, which starts at the rank
       ! and wraps around at n_procs
       ! i_block runs from 0 to n_procs-1, with wrapping at n_procs
       i_block = modulo(rank+i_pe,n_procs)
       
       i_block_recv = modulo(i_block+1,n_procs)
       i_block_send = modulo(rank-i_pe-1+n_procs,n_procs)

       !write(*,"(I3,A,4I3)") rank," : ",i_pe,i_block,i_block_recv,i_block_send

       ! do not wait for the first iteration, as nothing has been sent
       if (i_block.ne.rank) then
          ! wait for the received block
          !print*,rank,": Waiting for recv to finalize"
          call mpi_wait(recv_request,MPI_STATUS_IGNORE,ierr)
          ! swap the pointers
          temp_ptr => ptr_next_recv
          ptr_next_recv => ptr_next_compute
          ptr_next_compute => temp_ptr
          ! wait for the sended block
          !print*,rank,": Waiting for send to finalize"
          call mpi_wait(send_request,MPI_STATUS_IGNORE,ierr)
       end if
       ! initiate non-blocking receive for the next block from some following pe
       if (i_block_recv.ne.rank) then
          !print*,rank,": Starting irecv from ",i_block_recv
          call mpi_irecv(ptr_next_recv,mat%NBlocks*mat%Blocksize,MPI_COMPLEX_TYPE,&
               &i_block_recv,tag,bmat%PG%communicator,recv_request,ierr)
       end if
       ! initiate non-blocking send of the local array to some previous pe
       if (i_block_send.ne.rank) then
          !print*,rank,": Starting isend to ",i_block_send
          call mpi_isend(mat%localMatrix,mat%NBlocks*mat%Blocksize,MPI_COMPLEX_TYPE,&
               &i_block_send,tag,bmat%PG%communicator,send_request,ierr)
       end if
       ! do computation
       ! multiply bmat%localMatrix with ptr_next_compute
       ! bmat is a transposed matrix
       ! ptr_next_compute is a block of the full matrix
       global_index_start = i_block*mat%RowsPerBlock+1
       global_index_end   = (i_block+1)*mat%RowsPerBlock
       do global_bmat_row=rank*mat%RowsPerBlock+1,(rank+1)*mat%RowsPerBlock
          call mat_get_row_pointer(bmat,global_bmat_row,first_mat_row,s_global_col,e_global_col)
          do global_mat_col=1,mat%NCols
             ptr_local_column => ptr_next_compute((global_mat_col-1)*mat%RowsPerBlock+1:global_mat_col*mat%RowsPerBlock)
             loop_length = min(e_global_col,global_index_end)-max(s_global_col,global_index_start)+1
             if (s_global_col.ge.global_index_start) then
                ptr_fast_first => first_mat_row(1:loop_length)
                ptr_fast_second => ptr_local_column(s_global_col-global_index_start+1:s_global_col-global_index_start+loop_length)
             else
                ptr_fast_first => first_mat_row(global_index_start-s_global_col+1:global_index_start-s_global_col+loop_length)
                ptr_fast_second => ptr_local_column(1:loop_length)
             end if

             ! to save memory, we can write directly into the matrix with the following line
             !sum_value = cdotu(loop_length,ptr_fast_first,1,ptr_fast_second,1)
             !call add_value(res,global_bmat_row,global_mat_col,sum_value)
             ! for better performance, we collect a whole row of the result matrix
             if (loop_length.gt.0) then
                result_row(global_mat_col) = cdotu(loop_length,ptr_fast_first,1,ptr_fast_second,1)
             else
                result_row(global_mat_col) = cmplx(0.0D0,0.0D0)
             end if
          end do
          ! for better performance, now add the whole result_row to the result matrix
          call add_to_row(res,global_bmat_row, result_row)
       end do
       call commit_values(res)

       if (i_block.eq.rank) then
          ! the first iteration is special, as ptr_next_compute points
          ! to localMatrix, but to do the pointer swapping, we have to 
          ! redirect ptr_next_compute to the other temporary storage.
          ptr_next_compute => temp2
       end if
    end do
    deallocate(temp1,temp2)
  end subroutine sbm_dot_multiply_Banded_with_Full

  !> Calculate the minimal k and maximal k index on process proc
  SUBROUTINE get_minmax_k(mat,rank,min_k,max_k)
    TYPE(StoreBandedMatrixObject), intent(IN) :: mat
    INTEGER, intent(IN) :: rank
    INTEGER, INTENT(OUT) :: min_k, max_k

    ! Local variables
    INTEGER :: icol, irow, this_min_k, this_max_k, n_procs,ierr
    INTEGER :: n_local_rows,local_rows_start,local_rows_end,local_cols_start,local_cols_end

    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)
    
    IF (mat%isTransposed) THEN
       n_local_rows = mat%NRows/n_procs
       local_rows_start = rank*n_local_rows+1
       local_rows_end   = (rank+1)*n_local_rows
       local_cols_start = 1
       local_cols_end   = mat%NCols

       min_k = mat%NRows
       max_k = 1
       DO icol=local_cols_start,local_cols_end
          DO irow = local_rows_start,local_rows_end
             this_min_k = MAX(irow-mat%LowerBandwidth,icol-mat%UpperBandwidth,1)
             IF (this_min_k.LT.min_k) min_k = this_min_k
             this_max_k = MIN(irow+mat%UpperBandwidth,icol+mat%LowerBandwidth,mat%NRows)
             IF (this_max_k.GT.max_k) max_k = this_max_k
          END DO
       END DO
    ELSE
    END IF
    
  END SUBROUTINE get_minmax_k

  !> Calculate the square of the given Banded Matrix. The square of a banded matrix is again
  !! a banded matrix with double bandwidths.
  SUBROUTINE sbm_square_matrix(mat,res)
    TYPE(StoreBandedMatrixObject), INTENT(IN) :: mat
    TYPE(StoreBandedMatrixObject), intent(INOUT) :: res

    !Local variables
    INTEGER :: rank, n_procs, ierr, min_k, max_k, min_proc, max_proc, psource,iproc, pdest, tag
    TYPE(IntegerList) :: recv_list, send_list
    type(IntegerNodeData) :: nodedata
    !INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: send_request(MAX_SEND_REQUESTS), recv_request, irequest,n_requests,&
         &calculate_rank
    COMPLEX, DIMENSION(mat%NumberOfStoredRows,mat%ColsPerBlock),TARGET :: rbuf1,rbuf2
    COMPLEX, DIMENSION(:,:),POINTER :: ptr_calculate, ptr_recv
    LOGICAL :: wait_for_recv

    CALL exit_if_factored(mat)

    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)
    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)

    IF (mat%isTransposed) THEN
       call initialize(recv_list)
       call initialize(send_list)
       tag = 1001
       CALL get_minmax_k(mat,rank,min_k,max_k)
       !PRINT*,rank,": min_k = ", min_k,", max_k = ", max_k
       ! min_k and max_k contains the minimal used row or column on the local processor
       ! now we calculate the ranks of the processes from which we need data
       min_proc = (min_k-1)/mat%ColsPerBlock
       max_proc = (max_k-1)/mat%ColsPerBlock
       !PRINT*,"min_proc = ",min_proc,", max_proc = ", max_proc

       DO psource=min_proc,max_proc
          IF (psource.NE.rank) THEN
             nodedata%value = psource
             call push(recv_list,nodedata)
          END IF
       END DO
       !min/max_proc are the ranks of the processors from where we need data. 
       ! Further we need to know, to which ranks we have to send data
       DO iproc=0,n_procs-1
          CALL get_minmax_k(mat,iproc,min_k, max_k)
          min_proc = (min_k-1)/mat%ColsPerBlock
          max_proc = (max_k-1)/mat%ColsPerBlock
          !PRINT*,"for iproc = ", iproc,", min_k = ", min_k, ", max_k = ", max_k,", min_proc = ", min_proc, ", max_proc = ", max_proc
          IF ((rank.NE.iproc) .AND. (rank.LE.max_proc) .AND. (rank.GE.min_proc)) THEN
             nodedata%value=iproc
             CALL push(send_list,nodedata)
          END IF
       END DO

       ! We have now the full receicve list and send list 
       ! Now do all non blocking sends
       !IF (.NOT.isEmpty(send_list)) THEN
       !   WRITE(*,"(A)",advance="no") "send_list = "
       !   CALL output(send_list)
       !   WRITE(*,"(A)") ""
       !ELSE
       !   PRINT*,"send_list is empty"
       !END IF
       irequest = 1
       DO WHILE (.NOT.isEmpty(send_list))
          CALL pop(send_list,nodedata)
          pdest = nodedata%value
          PRINT*,"isend to ", pdest
          IF (irequest.GT.MAX_SEND_REQUESTS) THEN
             PRINT*,"You have to increase MAX_SEND_REQUESTS (actual set to ",MAX_SEND_REQUESTS,")!"
             stop
          END IF
          CALL mpi_isend(mat%localmatrix,mat%NumberOfStoredRows*mat%ColsPerBlock,MPI_COMPLEX_TYPE,&
               &pdest,tag,mat%PG%communicator,send_request(irequest),ierr)
          irequest = irequest + 1
       END DO
       n_requests = irequest-1

       ! receiving the data is overlapped with calculation
       ptr_calculate => mat%LocalMatrix
       ptr_recv => rbuf1
       !IF (.NOT.isEmpty(recv_list)) THEN
       !   WRITE(*,"(A)",advance="no") "recv_list = "
       !   CALL output(recv_list)
       !   WRITE(*,"(A)") ""
       !ELSE
       !   PRINT*,"recv_list is empty"
       !END IF

       !PRINT*,"length(recv_list) = ", length(recv_list),isEmpty(recv_list)
       calculate_rank = rank ! start with local multiplication
       DO 
          IF (.NOT.isEmpty(recv_list)) THEN
             CALL pop(recv_list,nodedata)
             psource = nodedata%value
             PRINT*,"irecv from ",psource
             CALL mpi_irecv(ptr_recv,mat%NumberOfStoredRows*mat%ColsPerBlock,MPI_COMPLEX_TYPE,&
                  &psource,tag,mat%PG%communicator,recv_request,ierr)
             wait_for_recv = .TRUE.
          ELSE
             wait_for_recv = .FALSE.
          END IF
          ! do the calculation LocalMatrix dot ptr_calculate
          PRINT*,"========================== local calculation ==============="
          CALL local_block_multiplication(mat,ptr_calculate,calculate_rank,res)
          IF (wait_for_recv) THEN
             CALL mpi_wait(recv_request,MPI_STATUS_IGNORE,ierr)
             PRINT*,"waiting for receive"
          
             IF (ASSOCIATED(ptr_recv,rbuf1)) THEN
                ptr_calculate => rbuf1
                ptr_recv => rbuf2
             ELSE
                ptr_calculate => rbuf2
                ptr_recv => rbuf1
             END IF
             calculate_rank = psource
          ELSE
             EXIT
          END IF
       END DO

       ! wait for finishing all non blocking sends
       DO irequest=1,n_requests
          !PRINT*,"Now waiting for send to complete for irequest = ",irequest," of ",n_requests
          CALL mpi_wait(send_request(irequest),MPI_STATUS_IGNORE,ierr)
       END DO

    ELSE
    END IF

  CONTAINS

    SUBROUTINE local_block_multiplication(mat,ptr_calculate,calculate_rank,res)
      TYPE(StoreBandedMatrixObject), INTENT(IN) :: mat
      COMPLEX, DIMENSION(:,:), INTENT(IN) :: ptr_calculate
      INTEGER, INTENT(IN) :: calculate_rank
      TYPE(StoreBandedMatrixObject), INTENT(INOUT) :: res
      
    END SUBROUTINE local_block_multiplication
  END SUBROUTINE sbm_square_matrix


  !*****************************************************************
  !* The global matrix, which is distributed over all processes
  !* is gathered and stored in a local array. This routine is
  !* mainly used for debugging and backward compatibility. In the
  !* final production code, it should not be called as it does not
  !* scale in memory.
  !*****************************************************************
  SUBROUTINE sbm_get_global_matrix_locally(mat,localfullmat)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    COMPLEX, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    ! local variables
    INTEGER :: ierr,iCol,iRow, global_col, global_row, rank

    CALL exit_if_factored(mat)

    ! otherwise, we have to take into account the cyclic block matrix
    ! structure, which complicates the operation.
    ! First we allocate a matrix to take the transposed local matrix
    !ALLOCATE(localfullmat(mat%NRows,mat%NCols))
    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)

    !PRINT*,"NumberOfStoredRows = ",mat%NumberOfStoredRows, mat%ColsPerBlock
    localfullmat = CMPLX(0.,0.,kind(localfullmat))
    ! and we fill this local matrix
    IF (mat%isTransposed) THEN
       DO icol=1,mat%ColsPerBlock
          global_col = rank*mat%ColsPerBlock+icol
          DO irow = -mat%LowerBandwidth,mat%UpperBandwidth
             global_row = global_col+irow
             IF (global_row.GE.1 .AND. global_row.LE.mat%NRows) THEN
                localfullmat(global_row,global_col) = mat%localMatrix(irow,icol)
             END IF
          END DO
       END DO
       
       !PRINT*, "Communicator in MPI_Allgather is ", mat%PG%Communicator
       CALL mpi_allgather(MPI_IN_PLACE,mat%NRows*mat%ColsPerBlock,MPI_COMPLEX_TYPE,&
            &localfullmat,mat%NRows*mat%ColsPerBlock,MPI_COMPLEX_TYPE,&
            &mat%PG%Communicator,ierr)
       localfullmat = transpose(localfullmat)
    ELSE
       DO icol=1,mat%ColsPerBlock
          global_col = rank*mat%ColsPerBlock+icol
          DO irow = -mat%UpperBandwidth,mat%LowerBandwidth
             global_row = global_col+irow
             IF (global_row.GE.1 .AND. global_row.LE.mat%NRows) THEN
                localfullmat(global_row,global_col) = mat%localMatrix(irow,icol)
             END IF
          END DO
       END DO
       
       !PRINT*, "Communicator in MPI_Allgather is ", mat%PG%Communicator
       CALL mpi_allgather(MPI_IN_PLACE,mat%NRows*mat%ColsPerBlock,MPI_COMPLEX_TYPE,&
            &localfullmat,mat%NRows*mat%ColsPerBlock,MPI_COMPLEX_TYPE,&
            &mat%PG%Communicator,ierr)
    END IF

  END SUBROUTINE sbm_get_global_matrix_locally

  FUNCTION sbm_get_local_abs_square_sum(mat) RESULT(value)
    type(StoreBandedMatrixObject) :: mat
    REAL :: value

    INTEGER :: icol, global_col, irow, global_row, rank, ierr

    CALL exit_if_factored(mat)

    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)

    value = 0.0D0
    IF (mat%isTransposed) THEN
       DO icol=1,mat%ColsPerBlock
          !global_col = rank*mat%ColsPerBlock+icol
          global_row = rank*mat%ColsPerBlock+icol

          DO irow = -mat%LowerBandwidth,mat%UpperBandwidth
             !global_row = global_col+irow
             global_col = global_row+irow
             !IF (global_row.GE.1 .AND. global_row.LE.mat%NRows) THEN
             IF (global_col.GE.1 .AND. global_col.LE.mat%NCols) THEN
                value = value +  REAL(mat%localMatrix(irow,icol)*CONJG(mat%localMatrix(irow,icol)))
             END IF
          END DO
       END DO
    ELSE
       DO icol=1,mat%ColsPerBlock
          global_col = rank*mat%ColsPerBlock+icol

          DO irow = -mat%UpperBandwidth,mat%LowerBandwidth
             global_row = global_col+irow
             IF (global_row.GE.1 .AND. global_row.LE.mat%NRows) THEN
                value = value +  REAL(mat%localMatrix(irow,icol)*CONJG(mat%localMatrix(irow,icol)))
             END IF
          END DO
       END DO
    END IF
  END FUNCTION sbm_get_local_abs_square_sum

  SUBROUTINE sbm_show_on_screen(mat)
    type(StoreBandedMatrixObject) :: mat

    ! Local variables
    INTEGER :: iProc,row_start,row_end,iCol, iRow

    CALL exit_if_factored(mat)

    ! do some debugging output
    !call mpi_barrier(mat%PG%Communicator,ierr)
    CALL blacs_barrier(mat%PG%Context,'A')
    DO iProc=0,mat%PG%NProcs-1
       !call mpi_barrier(mat%PG%Communicator,ierr)
       CALL blacs_barrier(mat%PG%Context,'A')
       IF (iProc.EQ.mat%PG%Rank) THEN
          PRINT*,"---- Start output Rank = ",mat%PG%Rank," ----"
          IF (mat%isTransposed) THEN
             row_start = -mat%LowerBandwidth
             row_end   = mat%UpperBandwidth
          ELSE
             row_start = -mat%UpperBandwidth
             row_end   = mat%LowerBandwidth
          END IF
          DO iRow=row_start,row_end
             DO iCol = 1,mat%ColsPerBlock
                WRITE(*,"(2ES9.2,1X)",advance='NO') mat%localMatrix(irow,icol)
             END DO
             WRITE(*,"(A)") ""
          END DO
          PRINT*,"---- End   output Rank = ",mat%PG%Rank," ----"
       END IF
    END DO
    !call mpi_barrier(mat%PG%Communicator,ierr)
    CALL blacs_barrier(mat%PG%Context,'A')
  END SUBROUTINE sbm_show_on_screen

  SUBROUTINE sbm_show_in_file(mat,filename)
    type(StoreBandedMatrixObject) :: mat
    CHARACTER(len=*) :: filename

    ! Local variables
    INTEGER :: iCol, iRow, row_start, row_end
    character(len=FILENAME_MAX) :: full_filename

    CALL exit_if_factored(mat)

    ! open a file on each process
    WRITE(full_filename,"(A,I3.3,A)") trim(filename),mat%PG%Rank,".dat"
    OPEN(77,file=TRIM(full_filename))

    ! do some debugging output
    CALL blacs_barrier(mat%PG%Context,'A')
    IF (mat%isTransposed) THEN
       row_start = -mat%LowerBandwidth
       row_end   = mat%UpperBandwidth
    ELSE
       row_start = -mat%UpperBandwidth
       row_end   = mat%LowerBandwidth
    END IF
    DO iRow=row_start,row_end
       DO iCol = 1,mat%ColsPerBlock
          WRITE(77,"(2ES20.10,1X)",advance='NO') mat%localMatrix(irow,icol)
       END DO
       WRITE(77,"(A)") ""
    END DO
    CLOSE(77)
    CALL blacs_barrier(mat%PG%Context,'A')
  END SUBROUTINE sbm_show_in_file

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE sbm_add_to_matrix(this,mat)
    TYPE(StoreBandedMatrixObject), INTENT(INOUT) :: this
    TYPE(StoreBandedMatrixObject), INTENT(IN) :: mat

    INTEGER :: rank, ierr, irow, icol, max_UpperBandwidth, max_LowerBandwidth

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)

    ! nrows and ncols have to be the same for both matrices,
    ! this is NOT checked!
    CALL mpi_comm_rank(this%PG%communicator, rank,ierr)
#if 0
    IF (this%UpperBandwidth.EQ.mat%UpperBandwidth) THEN
       ! 1. BWUs are equal
       DO icol=1,this%ColsPerBlock
          this%localMatrix(-this%UpperBandwidth:0,icol) = &
               & this%localMatrix(1:this%UpperBandwidth+1,icol) &
               & + mat%localMatrix(1:mat%UpperBandwidth+1,icol)
       END DO
    ELSEIF (this%UpperBandwidth.gt.mat%UpperBandwidth) THEN
       BWU_diff = this%UpperBandwidth-mat%UpperBandwidth
       ! 2. BWU of target is larger 
       DO icol=1,this%ColsPerBlock
          this%localMatrix(BWU_diff+1:BWU_diff+mat%UpperBandwidth+1,icol) = &
               & this%localMatrix(BWU_diff+1:BWU_diff+mat%UpperBandwidth+1,icol) + &
               & mat%localMatrix(1:mat%UpperBandwidth+1,icol)
       END DO
    ELSE
       ! 3. BWU of target is smaller than BWU of source
       DO icol=rank*mat%ColsPerBlock+1,(rank+1)*mat%ColsPerBlock
          DO irow=icol-mat%UpperBandwidth,icol
             IF ((irow.ge.1).AND.(irow.le.mat%NRows)) THEN
                CALL add_value(this,irow,icol,mat_get_value(mat,irow,icol))
             END IF
          END DO
       END DO
       call commit_values(this)
    END IF
#endif
    ! Get the values of mat and add them to this
    IF (mat%isTransposed.AND.this%isTransposed) THEN
       max_UpperBandwidth = MAX(mat%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(mat%LowerBandwidth, this%LowerBandwidth)
       DO irow = rank*mat%ColsPerBlock+1,(rank+1)*mat%ColsPerBlock
          DO icol=irow-max_LowerBandwidth,irow+max_UpperBandwidth
             IF ((icol.GE.1).AND.(icol.LE.mat%NCols)) THEN
                CALL add_value(this,irow,icol,mat_get_value(mat,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSEIF (.NOT.mat%isTransposed .AND. .NOT.this%isTransposed) THEN
       max_UpperBandwidth = MAX(mat%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(mat%LowerBandwidth, this%LowerBandwidth)
       DO icol=rank*mat%ColsPerBlock+1,(rank+1)*mat%ColsPerBlock
          DO irow=icol-max_UpperBandwidth,icol+max_LowerBandwidth
             IF ((irow.GE.1).AND.(irow.LE.mat%NRows)) THEN
                CALL add_value(this,irow,icol,mat_get_value(mat,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSE
       PRINT*,"You can add only two transposed matrices or two not transposed matrices, but not mixed."
       STOP
    END IF
  END SUBROUTINE sbm_add_to_matrix

  SUBROUTINE sbm_add_matrix(this,mat,res)
    TYPE(StoreBandedMatrixObject), INTENT(IN) :: this,mat
    TYPE(StoreBandedMatrixObject), INTENT(INOUT) :: res

    res = this

    CALL add_matrix(res,mat)
  END SUBROUTINE sbm_add_matrix

  SUBROUTINE sbm_subtract_matrix(this,mat,res)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: this,mat
    TYPE(StoreBandedMatrixObject),intent(INOUT) :: res
    ! calculate res = this - mat

    res = this

    CALL subtract_matrix(res,mat)

  END SUBROUTINE sbm_subtract_matrix

  SUBROUTINE sbm_subtract_from_matrix(this,mat)
    TYPE(StoreBandedMatrixObject),INTENT(INOUT) :: this
    TYPE(StoreBandedMatrixObject),intent(IN) :: mat
    ! calculate this = this - mat

    INTEGER :: rank, ierr, irow, icol, max_LowerBandwidth, max_UpperBandwidth

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)

    ! nrows and ncols have to be the same for both matrices,
    ! this is NOT checked!
    CALL mpi_comm_rank(this%PG%communicator, rank,ierr)

    ! Get the values of mat and add them to this
    IF (mat%isTransposed.AND.this%isTransposed) THEN
       max_UpperBandwidth = MAX(mat%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(mat%LowerBandwidth, this%LowerBandwidth)
       DO irow = rank*mat%ColsPerBlock+1,(rank+1)*mat%ColsPerBlock
          DO icol=irow-max_LowerBandwidth,irow+max_UpperBandwidth
             IF ((icol.GE.1).AND.(icol.LE.mat%NCols)) THEN
                CALL add_value(this,irow,icol,-mat_get_value(mat,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSEIF (.NOT.mat%isTransposed .AND. .NOT.this%isTransposed) THEN
       max_UpperBandwidth = MAX(mat%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(mat%LowerBandwidth, this%LowerBandwidth)
       DO icol=rank*mat%ColsPerBlock+1,(rank+1)*mat%ColsPerBlock
          DO irow=icol-max_UpperBandwidth,icol+max_LowerBandwidth
             IF ((irow.GE.1).AND.(irow.LE.mat%NRows)) THEN
                CALL add_value(this,irow,icol,-mat_get_value(mat,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSE
       PRINT*,"You can subtract only two transposed matrices or two not transposed matrices, but not mixed."
       STOP
    END IF

  END SUBROUTINE sbm_subtract_from_matrix

  SUBROUTINE sbm_multiply_matrix_with_real(this, scalar, res)
    TYPE(StoreBandedMatrixObject), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(StoreBandedMatrixObject),intent(INOUT) :: res

    INTEGER :: irow,icol,rank,ierr, max_LowerBandwidth, max_UpperBandwidth

    CALL exit_if_factored(this)
    CALL exit_if_factored(res)

    CALL mpi_comm_rank(this%PG%communicator,rank,ierr)

    ! Get the values of mat and add them to this
    IF (res%isTransposed.AND.this%isTransposed) THEN
       max_UpperBandwidth = MAX(res%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(res%LowerBandwidth, this%LowerBandwidth)
       DO irow = rank*res%ColsPerBlock+1,(rank+1)*res%ColsPerBlock
          DO icol=irow-max_LowerBandwidth,irow+max_UpperBandwidth
             IF ((icol.GE.1).AND.(icol.LE.res%NCols)) THEN
                CALL set_value(res,irow,icol,scalar*mat_get_value(this,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSEIF (.NOT.res%isTransposed .AND. .NOT.this%isTransposed) THEN
       max_UpperBandwidth = MAX(res%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(res%LowerBandwidth, this%LowerBandwidth)
       DO icol=rank*res%ColsPerBlock+1,(rank+1)*res%ColsPerBlock
          DO irow=icol-max_UpperBandwidth,icol+max_LowerBandwidth
             IF ((irow.GE.1).AND.(irow.LE.res%NRows)) THEN
                CALL set_value(res,irow,icol,scalar*mat_get_value(this,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSE
       PRINT*,"You can subtract only two transposed matrices or two not transposed matrices, but not mixed."
       STOP
    END IF

  END SUBROUTINE sbm_multiply_matrix_with_real

  SUBROUTINE sbm_scale_matrix_by_real(this, scalar)
    TYPE(StoreBandedMatrixObject), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL exit_if_factored(this)

    this%localMatrix = this%localMatrix*scalar
  END SUBROUTINE sbm_scale_matrix_by_real

  SUBROUTINE sbm_assign_matrix(lmat,rmat)
    TYPE(StoreBandedMatrixObject), intent(INOUT) :: lmat
    TYPE(StoreBandedMatrixObject), intent(IN) :: rmat

    INTEGER :: irow,icol,rank,ierr

    CALL exit_if_factored(lmat)
    CALL exit_if_factored(rmat)

    IF (lmat%PG.EQ.rmat%PG) THEN
       CALL mpi_comm_rank(lmat%PG%communicator,rank,ierr)
       IF (lmat%isTransposed.AND.rmat%isTransposed) THEN
          ! both matrices are stored in row-major order
          DO irow=rank*rmat%ColsPerBlock+1,(rank+1)*rmat%ColsPerBlock
             DO icol=irow-rmat%LowerBandwidth,irow+rmat%UpperBandwidth
                IF ((icol.GE.1).AND.(icol.LE.rmat%NCols)) THEN
                   CALL set_value(lmat,irow,icol,mat_get_value(rmat,irow,icol))
                END IF
             END DO
          END DO
          CALL commit_values(lmat)
       ELSEIF (.not.lmat%isTransposed.and..not.rmat%isTransposed) then
          DO icol=rank*rmat%ColsPerBlock+1,(rank+1)*rmat%ColsPerBlock
             DO irow=icol-rmat%UpperBandwidth,icol+rmat%LowerBandwidth
                IF ((irow.GE.1).AND.(irow.LE.rmat%NRows)) THEN
                   CALL set_value(lmat,irow,icol,mat_get_value(rmat,irow,icol))
                END IF
             END DO
          END DO
          CALL commit_values(lmat)
       ELSE
          PRINT*,"Assignment is only possible if both matrices have the same storage scheme (isTransposed)."
          STOP
       END IF
    ELSE
       PRINT*,"Assignment of matrices is only possible within the same context."
       STOP
    END IF
  END SUBROUTINE sbm_assign_matrix
  
  FUNCTION sbm_get_number_of_bands(mat)
    type(StoreBandedMatrixObject) :: mat
    INTEGER :: sbm_get_number_of_bands

    sbm_get_number_of_bands=mat%LowerBandwidth+mat%UpperBandwidth+1
  END FUNCTION sbm_get_number_of_bands

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE sbm_matrix_sum(mat,communicator)
    TYPE(StoreBandedMatrixObject), INTENT(INOUT) :: mat
    integer :: communicator

    INTEGER :: max_lower_bw, max_upper_bw, ierr, icol,irow
    COMPLEX, DIMENSION(:,:), pointer :: newMatrix

    CALL exit_if_factored(mat)

    ! For a banded matrix, we have to 
    ! 1) get the maximal bandwidths over the communicator
    CALL mpi_allreduce(mat%LowerBandwidth,max_lower_bw,1,MPI_INTEGER,&
         &MPI_MAX,communicator,ierr)
    CALL mpi_allreduce(mat%UpperBandwidth,max_upper_bw,1,MPI_INTEGER,&
         &MPI_MAX,communicator,ierr)
    
    ! 2) extend all matrices to this maximum bandwidth
    IF ((max_lower_bw.GT.mat%LowerBandwidth).OR.(max_upper_bw.GT.mat%UpperBandwidth)) THEN
       ! Attention, newMatrix will not be deallocated, as the mat%localMatrix will point to it!
       IF (mat%isTransposed) THEN
          ALLOCATE(newMatrix(-max_lower_bw:max_upper_bw,mat%ColsPerBlock))
       ELSE
          ALLOCATE(newMatrix(-max_upper_bw:max_lower_bw,mat%ColsPerBlock))
       END IF
       newMatrix = CMPLX(0,0,kind(newMatrix))

       ! 4. copy old matrix into new one
       IF (mat%isTransposed) THEN
          DO icol=1,mat%ColsPerBlock
             DO irow=-mat%LowerBandwidth,mat%UpperBandwidth
                newMatrix(irow,icol) = mat%localMatrix(irow,icol)
             END DO
          END DO
       ELSE
          DO icol=1,mat%ColsPerBlock
             DO irow=-mat%UpperBandwidth,mat%LowerBandwidth
                newMatrix(irow,icol) = mat%localMatrix(irow,icol)
             END DO
          END DO
       END IF
       DEALLOCATE(mat%localMatrix)
       mat%localMatrix => newMatrix
       mat%UpperBandwidth = max_upper_bw
       mat%LowerBandwidth = max_lower_bw
       mat%NumberOfStoredRows = mat%UpperBandwidth+mat%LowerBandwidth+1
    END IF

    ! all matrices have now the same shape, so we can
    ! 3) allreduce over all matrices

    CALL mpi_allreduce(MPI_IN_PLACE,mat%localMatrix,mat%NumberOfStoredRows*mat%ColsPerBlock, &
         & MPI_COMPLEX_TYPE, MPI_SUM, communicator, ierr)
    
  END SUBROUTINE sbm_matrix_sum

  SUBROUTINE sbm_row_axpy(this,this_row,mat,scalar)
    TYPE(StoreBandedMatrixObject),intent(INOUT) :: this
    INTEGER, INTENT(IN) :: this_row
    TYPE(StoreBandedMatrixObject),intent(IN) :: mat
    COMPLEX, intent(IN) :: scalar

    ! Local variables
    INTEGER :: rank, rank_of_global_row,irow, global_col,ierr, min_row,max_row, stored_column

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)

    CALL mpi_comm_rank(this%PG%communicator,rank,ierr)
    
    IF (this%isTransposed) THEN
       IF (.NOT.mat%isTransposed) THEN
          PRINT*,"sbm_row_axpy: both matrices must have the same storage mode"
          STOP
       END IF
       ! both matrices are transposed
       
       ! calculate the local index out of the global row index
       rank_of_global_row = (this_row-1)/this%ColsPerBlock
       ! do only the calculation, if the global row is stored locally on the actual processor
       IF (rank.EQ.rank_of_global_row) THEN
          min_row = MIN(this%LowerBandwidth,mat%LowerBandwidth)
          max_row = MIN(this%UpperBandwidth,mat%UpperBandwidth)
          stored_column = MOD(this_row-1,this%ColsPerBlock)+1
          this%localMatrix(-min_row:max_row,stored_column) = &
               & this%localMatrix(-min_row:max_row,stored_column)&
               & + mat%localMatrix(-min_row:max_row,stored_column)*scalar
          IF (mat%LowerBandwidth.GT.this%LowerBandwidth) THEN
             DO irow=-mat%LowerBandwidth,-min_row-1
                ! convert to global column index
                global_col = this_row+irow
                CALL set_value(this,this_row,global_col,mat%localMatrix(irow,stored_column)*scalar)
             END DO
          END IF
          IF (mat%UpperBandwidth.GT.this%UpperBandwidth) THEN
             DO irow=max_row+1,mat%UpperBandwidth
                ! convert to global column index
                global_col = this_row+irow
                CALL set_value(this,this_row,global_col,mat%localMatrix(irow,stored_column)*scalar)
             END DO
          END IF
       END IF
       CALL commit_values(this)
    ELSE
       IF (mat%isTransposed) THEN
          PRINT*,"sbm_row_axpy: both matrices must have the same storage mode"
          STOP
       END IF
       ! both matrices are not transposed
    END IF

  END SUBROUTINE sbm_row_axpy

  subroutine sbm_transpose_and_conjugate(mat,dagger_mat)
    type(StoreBandedMatrixObject) :: mat
    type(StoreBandedMatrixObject) :: dagger_mat

    ! To transpose, we set the isTransposed flag of dagger_mat
    ! to the inverse of the isTransposed flag of mat
    dagger_mat%isTransposed = .not.mat%isTransposed
    ! Also exchange NRows and NCols, and the Bandwidths
    dagger_mat%NRows = mat%NCols
    dagger_mat%NCols = mat%NRows
    dagger_mat%LowerBandwidth = mat%UpperBandwidth
    dagger_mat%UpperBandwidth = mat%LowerBandwidth
    dagger_mat%NumberOfStoredRows = mat%NumberOfStoredRows
    dagger_mat%ColsPerBlock = mat%ColsPerBlock

    ! now copy the localmatrices and conjugate
    call allocate(dagger_mat)
    !if (associated(dagger_mat%localmatrix)) deallocate(dagger_mat%localmatrix)
    
    !allocate(dagger_mat%localmatrix(-dagger_mat%UpperBandwidth:dagger_mat%LowerBandwidth,&
    !     &1:dagger_mat%ColsPerBlock))
    dagger_mat%localmatrix = conjg(mat%localmatrix)
  end subroutine sbm_transpose_and_conjugate

  SUBROUTINE sbm_LU_factor(mat)
    type(StoreBandedMatrixObject) :: mat
    
    ! Local variables
    COMPLEX :: work(1)
    INTEGER :: i, n, bwu, bwl, bw_max, laf, lwork, info, desca(7)
    INTEGER :: my_rank, n_procs, mpierr
    !INTEGER :: mpi_status(mpi_status_size)
    INTEGER :: np_send, np_recv
    COMPLEX, ALLOCATABLE :: sbuf(:), rbuf(:)

    CALL exit_if_factored(mat)

    ! if PG_cols is not yet initialized we have to do it here

    if(.NOT. PG_cols%isInitialized) CALL initialize(PG_cols,MAT_BLOCKS_OF_COLS,mat%PG%Communicator)

    call mpi_comm_rank(mat%PG%Communicator,my_rank,mpierr)
    call mpi_comm_size(mat%PG%Communicator,n_procs,mpierr)

    IF (mat%Nrows /= mat%Ncols) THEN
      PRINT *,'Matrix must be square for sbm_factor_matrix'
      CALL mpi_abort(MPI_COMM_WORLD, 1, mpierr)
    ENDIF

    bwu = mat%UpperBandwidth
    bwl = mat%LowerBandwidth

    bw_max = max(bwu,bwl)

    ! pcgbtrf/pcgbtrs work only if bwu+bwl+1 <= mat%ColsPerBlock.
    ! NB: The transposition below works only if bw_max <= mat%ColsPerBlock,
    ! but this is a weaker restriction

    IF(bwu+bwl+1 > mat%ColsPerBlock) THEN
      print *,"ERROR sbm_factor: bwu+bwl+1, ColsPerBlock ",bwu+bwl+1,mat%ColsPerBlock
      CALL mpi_abort(MPI_COMM_WORLD, 1, mpierr)
    ENDIF

    ALLOCATE(mat%factoredBand(-2*bwu-bwl:bwl,mat%ColsPerBlock))
    mat%factoredBand(:,:) = 0

    ! fill in mat%factoredBand

    IF(mat%isTransposed) THEN

      ! Please note: pcgbtrf/pcgbtrs are NOT able to deal with transposed matrices
      ! although the first argument of pcgbtrs pretends that.
      ! In reality the first argument of pcgbtrs is completely ignored.
      ! Thus we have to transpose the matrix here ...

      allocate(sbuf(bw_max*(bw_max+1)/2))
      allocate(rbuf(bw_max*(bw_max+1)/2))

      ! Super diagonals

      ! Send pieces at end of localMatrix to next processor

      ! Please note: It is not necessary to do a cyclic exchange
      ! but it is convenient to do so when using sendrecv
      np_send = MOD(my_rank+1,n_procs)
      np_recv = MOD(my_rank+n_procs-1,n_procs)

      n = 0
      DO i = 1, bwu
        sbuf(n+1:n+i) = mat%localMatrix(i,mat%ColsPerBlock-i+1:mat%ColsPerBlock)
        n = n+i
      ENDDO
      if(my_rank==n_procs-1) sbuf(:) = 0 ! Safety only

      call mpi_sendrecv(sbuf, n, MPI_COMPLEX_TYPE, np_send, 111, &
                        rbuf, n, MPI_COMPLEX_TYPE, np_recv, 111, &
                        mat%PG%Communicator, MPI_STATUS_IGNORE, mpierr)

      n = 0
      DO i = 1, bwu
        mat%factoredBand(-i,1+i:mat%ColsPerBlock) = mat%localMatrix(i,1:mat%ColsPerBlock-i)
        mat%factoredBand(-i,1:i) = rbuf(n+1:n+i)
        n = n+i
      ENDDO

      ! Diagonal

      mat%factoredBand(0,1:mat%ColsPerBlock) = mat%localMatrix(0,1:mat%ColsPerBlock)

      ! Sub diagonals

      ! Send pieces at start of localMatrix to previous processor

      np_recv = MOD(my_rank+1,n_procs)
      np_send = MOD(my_rank+n_procs-1,n_procs)

      n = 0
      DO i = 1, bwl
        sbuf(n+1:n+i) = mat%localMatrix(-i,1:i)
        n = n+i
      ENDDO

      if(my_rank==0) sbuf(:) = 0 ! Safety only

      call mpi_sendrecv(sbuf,n,MPI_COMPLEX_TYPE,np_send,111, &
                        rbuf,n,MPI_COMPLEX_TYPE,np_recv,111, &
                        mat%PG%Communicator, MPI_STATUS_IGNORE, mpierr)

      n = 0
      DO i = 1, bwl
        mat%factoredBand(i,1:mat%ColsPerBlock-i)  = mat%localMatrix(-i,1+i:mat%ColsPerBlock)
        mat%factoredBand(i,mat%ColsPerBlock-i+1:mat%ColsPerBlock) = rbuf(n+1:n+i)
        n = n+i
      ENDDO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)

    ELSE

      mat%factoredBand(-bwu:bwl,:) = mat%localMatrix(-bwu:bwl,:)

    ENDIF

    ! The doc in the header of PDGBTRF isn't clear about IPIV,
    ! I hope the following is enough ...
    ALLOCATE(mat%ipiv(2*mat%ColsPerBlock+bwu+bwl+1))

    ! allocate additional storage
    laf = (mat%ColsPerBlock+bwu)*(bwl+bwu)+6*(bwl+bwu)*(bwl+2*bwu)
    ALLOCATE(mat%auxData(laf))

    ! Size of work is not stated in header of pcgbtrf, obviously it is 1
    lwork = SIZE(work)

    desca = (/ 501, PG_cols%Context, mat%Nrows, mat%ColsPerBlock, 0, 2*bwl+2*bwu+1, 0 /)

    CALL pcgbtrf(mat%Nrows, bwl, bwu, mat%factoredBand, 1, desca, &
                 mat%ipiv, mat%auxData, laf, work, lwork, info)

    IF (info.NE.0) THEN
       PRINT*,"pcgbtrf info = ",info
       stop
    END IF

    mat%isFactored = .true.

    ! Now we can deallocate the original matrix data

    deallocate(mat%localMatrix)
    nullify(mat%localMatrix)

  END SUBROUTINE sbm_LU_factor

  LOGICAL FUNCTION sbm_LU_factor_ok(mat)
    type(StoreBandedMatrixObject) :: mat

    ! Checks if the bandwidth requirement for pcgbtrf/pcgbtrs is ok.

    IF(mat%isFactored) THEN
      ! In this case it is of course not ok to factor ...
      sbm_LU_factor_ok = .FALSE.
      RETURN
    ENDIF

    IF(mat%ColsPerBlock >= mat%LowerBandwidth+mat%UpperBandwidth+1) THEN
      sbm_LU_factor_ok = .TRUE.
    ELSE
      sbm_LU_factor_ok = .FALSE.
    ENDIF
  END FUNCTION sbm_LU_factor_ok

  SUBROUTINE sbm_LU_solve(mat, vec, res)
    TYPE(StoreBandedMatrixObject),INTENT(IN) :: mat
    TYPE(StoreVectorObject),INTENT(IN) :: vec
    TYPE(StoreVectorObject),INTENT(INOUT) :: res

    INTEGER :: bwu, bwl, lwork, info, desca(7), descb(7)
    COMPLEX, ALLOCATABLE :: work(:)

    ! checks if vec and res fit to mat

    if(vec%RowsPerBlock /= mat%ColsPerBlock) STOP 'sbm_LU_solve: mat and vec not conforming'
    if(res%RowsPerBlock /= mat%ColsPerBlock) STOP 'sbm_LU_solve: mat and res not conforming'

    if(.not. mat%isFactored) STOP 'solve_LU called with un-factored matrix'

    bwu = mat%UpperBandwidth
    bwl = mat%LowerBandwidth

    ! set up despriptors

    desca = (/ 501, PG_cols%Context, mat%Nrows, mat%ColsPerBlock, 0, 2*bwl+2*bwu+1, 0 /)
    descb = (/ 502, PG_cols%Context, mat%Nrows, mat%ColsPerBlock, 0, mat%ColsPerBlock, 0 /)

    ! Size of work is not stated in header of pcgbtrs, however in the source code:

    lwork = mat%ColsPerBlock + 2*bwl + 4*bwu
    ALLOCATE(work(lwork))

    res%localVector = vec%localVector

    CALL pcgbtrs('N', mat%Nrows, bwl, bwu, 1, mat%factoredBand, 1, desca, &
                 mat%ipiv, res%localVector, 1, descb, mat%auxData, size(mat%auxData), &
                 work, lwork, info)

    DEALLOCATE(work)

    IF (info.NE.0) THEN
       PRINT*,"pcgbtrs info = ",info
       stop
    END IF

  END SUBROUTINE sbm_LU_solve

  SUBROUTINE sbm_exit_if_factored(mat)
    type(StoreBandedMatrixObject) :: mat

    INTEGER :: ierr

    ! This routines is called by all routines which use mat%localMatrix.
    ! It does exit the program if the matrix is LU-factored.

    IF(mat%isFactored) THEN
      PRINT *,'ERROR: Attempt use original matrix data of a LU-factored Matrix'
      CALL mpi_abort(MPI_COMM_WORLD, 1, ierr)
    ENDIF

  END SUBROUTINE sbm_exit_if_factored

  subroutine sbm_print_storage_details(mat)
    type(StoreBandedMatrixObject) :: mat

    print*,"Storage details: "
    print*,"ColsPerBlock,NumberOfStoredRows = ",mat%ColsPerBlock,mat%NumberOfStoredRows
  end subroutine sbm_print_storage_details

#ifdef ALL_ROUTINES
  ! ===================================================
  ! == Define some mathematical routines on matrices ==
  ! ===================================================
  SUBROUTINE sbm_invert_matrix(mat)
    type(StoreBandedMatrixObject) :: mat
    
    ! Local variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv, iwork
    COMPLEX, DIMENSION(:), allocatable :: work
    INTEGER :: info,lwork,liwork,ub_locr

    ! calculate upper bound of ipiv storage needed
    ! LOCr(M_A)+MB_A (ScaLapack User's guide)
    ub_locr = ( (mat%NRows+mat%RowsPerBlock-1)/mat%RowsPerBlock &
         & + mat%PG%NProcs-1 ) / mat%PG%NProcs * mat%RowsPerBlock

    !PRINT*,"ub_locr = ",ub_locr
    ALLOCATE(ipiv(ub_locr+mat%RowsPerBlock))
    !WRITE(*,"(16ES10.2)") mat%localMatrix
    ! First do the LU decomposition in place
    CALL pcgetrf(mat%NRows,mat%NCols,&
         &mat%localMatrix, 1, 1,&
         &mat%Desc,ipiv,info)

    IF (info.NE.0) THEN
       PRINT*,"pcgetrf produced an error in StoreBandedMatrix/sbm_invert_matrix"
       stop
    END IF
    !PRINT*,"after LU decomposition"
    !WRITE(*,"(16ES10.2)") mat%localMatrix
    

    ! Second calculate the inverse of the matrix
    ! determine the work array size
    ALLOCATE(work(1),iwork(1))
    CALL pcgetri(mat%Nrows,mat%localMatrix,1,1,mat%Desc,ipiv,&
         &work,-1,iwork,-1,info)
    lwork = INT(work(1))
    liwork = iwork(1)
    !PRINT*,"lwork = ",lwork,", liwork = ",liwork
    DEALLOCATE(work,iwork)

    !PRINT*,"after first pcgetri call"
    !WRITE(*,"(16ES10.2)") mat%localMatrix

    
    ALLOCATE(work(lwork),iwork(liwork))
    CALL pcgetri(mat%Nrows,mat%localMatrix,1,1,mat%Desc,ipiv,&
         &work,lwork,iwork,liwork,info)

    IF (info.NE.0) THEN
       PRINT*,"pcgetri produced an error in StoreBandedMatrix/sbm_invert_matrix"
       stop
    END IF
    !PRINT*,"after second pcgetri call"
    !WRITE(*,"(16ES10.2)") mat%localMatrix

    DEALLOCATE(work,iwork)
    DEALLOCATE(ipiv)
  END SUBROUTINE sbm_invert_matrix
#endif

  !*****************************************************************
  !* Repeat everything, but now for REAL datatypes
  !*****************************************************************

  SUBROUTINE sbm_initialize_matrix_real(mat,n_rows,n_cols,PG,transp)
    TYPE(StoreBandedMatrixObjectReal), intent(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    type(ProcessGrid) :: PG
    LOGICAL, INTENT(IN),OPTIONAL :: transp

    ! Local variables
  
    mat%NRows = n_rows
    mat%NCols = n_cols
    mat%PG = PG

    ! Distribute matrix on the grid
    ! block-column matrix
    !mat%ColsPerBlock = mat%NCols / PG%NPRows
    mat%ColsPerBlock = mat%NCols / PG%NProcs

    !mat%Desc = (/ 501, PG%Context, mat%NCols, mat%ColsPerBlock, 0, mat%LowerBandwidth+mat%UpperBandwidth+1,0 /)

    mat%LowerBandwidth = 0
    mat%UpperBandwidth = 0
    mat%NumberOfStoredRows = 1

    CALL initialize(mat%storage_list)

    mat%isInitialized = .TRUE.
    IF (PRESENT(transp)) THEN
       mat%isTransposed = transp
    ELSE
       mat%isTransposed = .false.
    END IF

    ! Get pointers into a defined state
    NULLIFY(mat%localMatrix)
    NULLIFY(mat%factoredBand)
    NULLIFY(mat%auxData)
    NULLIFY(mat%ipiv)

    CALL blacs_barrier(PG%Context,'A')
  END SUBROUTINE sbm_initialize_matrix_real

  LOGICAL FUNCTION sbm_isInitialized_real(mat) 
    type(StoreBandedMatrixObjectReal) :: mat

    sbm_isInitialized_real = mat%isInitialized
  END FUNCTION sbm_isInitialized_real


  SUBROUTINE sbm_allocate_real(mat,lbw_in,ubw_in)
    TYPE(StoreBandedMatrixObjectReal) :: mat
    integer, optional,intent(IN) :: lbw_in, ubw_in

    if (present(lbw_in)) then
       if (lbw_in.le.mat%NCols) then
          mat%LowerBandwidth = lbw_in
       else
          print*,"Lower Bandwidth can not be larger than dim of matrix!"
          !CALL tracebackqq()

          stop
       end if
    end if

    if (present(ubw_in)) then
       if (ubw_in.le.mat%NCols) then
          mat%UpperBandwidth = ubw_in
       else
          print*,"Upper Bandwidth can not be larger than dim of matrix!"
          stop
       end if
    end if

    if (associated(mat%localMatrix)) deallocate(mat%localMatrix)

    IF (mat%isTransposed) THEN
       ALLOCATE(mat%localMatrix(-mat%LowerBandwidth:mat%UpperBandwidth,mat%ColsPerBlock))
    ELSE
       ALLOCATE(mat%localMatrix(-mat%UpperBandwidth:mat%LowerBandwidth,mat%ColsPerBlock))
    END IF
    mat%NumberOfStoredRows = mat%LowerBandwidth+mat%UpperBandwidth+1
    mat%localMatrix = 0.0
  END SUBROUTINE sbm_allocate_real

  SUBROUTINE sbm_finalize_matrix_real(mat)
    type(StoreBandedMatrixObjectReal) :: mat

    if (ASSOCIATED(mat%localMatrix)) THEN
       DEALLOCATE(mat%localMatrix)
    END IF
    IF (ASSOCIATED(mat%factoredBand)) THEN
       DEALLOCATE(mat%factoredBand)
    END IF
    IF (ASSOCIATED(mat%auxData)) THEN
       DEALLOCATE(mat%auxData)
    END IF
    IF (ASSOCIATED(mat%ipiv)) THEN
       DEALLOCATE(mat%ipiv)
    END IF
    CALL finalize(mat%storage_list)
    mat%isInitialized=.false.
  END SUBROUTINE sbm_finalize_matrix_real

  SUBROUTINE sbm_set_value_real(mat,irow,icol,value)
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc, row_index, col_index
    INTEGER :: distance!, coord(2)
    TYPE(MatrixEntryNodeData) :: nodedata
    real :: tolerance

    CALL exit_if_factored(mat)

    DEBSTART("set_value(SBMO)")
    !tolerance = REAL(0.0,KIND(tolerance))
    ! only consider value if value is larger than tolerance (numeric zero)
    IF ( (REAL(value).NE.REAL(0.0,KIND(tolerance))) ) THEN
    !IF ((ABS(REAL(value)).GT.tolerance).OR.(ABS(AIMAG(value)).GT.tolerance)) THEN
       ! calculate the processor coordinates, just determined by the global column index
       ! or global row index, depending on the storage scheme
       IF (mat%isTransposed) THEN
          proc  = (irow-1)/mat%ColsPerBlock
       ELSE
          proc  = (icol-1)/mat%ColsPerBlock
       END IF

       ! only process the entry if we are on the right processor
       IF (proc.EQ.mat%PG%Rank) THEN
          ! check, if the value is inside the band 
          distance = icol-irow
          IF ((distance>mat%UpperBandwidth).OR.(-distance>mat%LowerBandwidth)) THEN
             ! out of bands
             ! put value and indices into list for later insertion
             nodedata%coord(1)=irow
             nodedata%coord(2)=icol
             nodedata%value = value
             CALL push(mat%storage_list,nodedata)
          ELSE
             IF (mat%isTransposed) THEN
                row_index = MOD(irow-1,mat%ColsPerBlock)+1
                col_index = distance

                mat%localMatrix(col_index, row_index) = value
             ELSE
                col_index = MOD(icol-1,mat%ColsPerBlock)+1
                row_index = -distance

                mat%localMatrix(row_index, col_index) = value
             END IF
          END IF
       END IF
    end if
    DEBEND("set_value(SBMO)")
  END SUBROUTINE sbm_set_value_real

  SUBROUTINE sbm_add_value_real(mat,global_row,global_col,value)
    TYPE(StoreBandedMatrixObjectReal) :: mat
    INTEGER, INTENT(IN) :: global_row,global_col
    REAL, INTENT(IN) :: value

    ! Local variables
    REAL :: old_value
    integer :: proc

    DEBSTART("add_value(SFMO)")

    !CALL exit_if_factored(mat)

    IF (mat%isTransposed) THEN
       proc = (global_row-1)/mat%ColsPerBlock
    ELSE
       proc = (global_col-1)/mat%ColsPerBlock
    END IF

    IF (proc.EQ.mat%PG%Rank) THEN
       old_value = mat_get_value(mat,global_row,global_col)
       CALL set_value(mat,global_row,global_col,old_value+value)
    END IF
    DEBEND("add_value(SFMO)")
    
  END SUBROUTINE sbm_add_value_real

  FUNCTION sbm_get_value_real(mat,global_row,global_col,trans_in) RESULT(value)
    TYPE(StoreBandedMatrixObjectReal) :: mat
    INTEGER :: global_row, global_col
    real :: value
    character(len=1),optional,intent(IN) :: trans_in

    ! local variables
    INTEGER :: proc, row_index, col_index, distance
    character(len=1) :: trans

    !CALL exit_if_factored(mat)

    if (present(trans_in)) then
       trans=trans_in
    else
       trans='N'
    end if

    ! calculate the processor coordinates, just determined by the global column index
    IF (mat%isTransposed) THEN
       if (trans.eq.'N') then
          proc = (global_row-1)/mat%ColsPerBlock
       else
          proc = (global_col-1)/mat%ColsPerBlock
       end if
    ELSE
       if (trans.eq.'N') then
          proc = (global_col-1)/mat%ColsPerBlock
       else
          proc = (global_row-1)/mat%ColsPerBlock
       end if
    END IF

    ! only process the entry if we are on the right processor
    IF (proc.EQ.mat%PG%Rank) THEN
       ! check, if the value is inside the band 
       if (trans.eq.'N') then
          distance = global_col-global_row
       else
          ! cols and rows are exchanged
          distance = global_row-global_col
       end if
       IF ((distance.le.mat%UpperBandwidth).and.(distance.ge.-mat%LowerBandwidth)) THEN
          IF (mat%isTransposed) THEN
             if (trans.eq.'N') then
                col_index = MOD(global_row-1,mat%ColsPerBlock)+1
             else
                col_index = MOD(global_col-1,mat%ColsPerBlock)+1
             end if
             row_index = distance
          ELSE
             if (trans.eq.'N') then
                col_index = MOD(global_col-1,mat%ColsPerBlock)+1
             else
                col_index = MOD(global_row-1,mat%ColsPerBlock)+1
             end if
             row_index = -distance
          END IF
          
          if (trans.ne.'C') then
             value = mat%localMatrix(row_index, col_index)
          else
             value = mat%localMatrix(row_index, col_index)
          end if
       ELSE
          ! out of bands
          ! return zero
          value = 0.0
       END IF
    ELSE
       PRINT*,"This line should never be reached in mat_get_value of storebandedmatrix.F90."
       !CALL tracebackqq()
       value = 1000
       ! this value can be arbitrary, as it should never be used in further calculations
    END IF
  END FUNCTION sbm_get_value_real

  subroutine sbm_get_row_pointer_real(mat,global_row,ptr_row,start_col,end_col)
    TYPE(StoreBandedMatrixObjectReal) :: mat
    INTEGER,intent(IN) :: global_row
    REAL,dimension(:),pointer,intent(OUT) :: ptr_row
    integer, INTENT(OUT) :: start_col, end_col

    ! local variables
    INTEGER :: proc, local_row_index

    !CALL exit_if_factored(mat)

    if (.not.mat%isTransposed) then
       print*,"get_row_pointer only works for transposed matrices in storebandedmatrix"
       stop
    end if

    ! mat is transposed, stored in row-major order
    ! calculate the processor coordinates, just determined by the global column index
    proc = (global_row-1)/mat%ColsPerBlock

    ! only process the entry if we are on the right processor
    IF (proc.EQ.mat%PG%Rank) THEN
       local_row_index = mod(global_row-1,mat%ColsPerBlock)+1
       start_col = max(global_row - mat%LowerBandwidth,1)
       end_col   = min(global_row + mat%UpperBandwidth,mat%NCols)
       
       ptr_row  => mat%localMatrix(start_col-global_row:end_col-global_row,local_row_index)
    else
       nullify(ptr_row)
    end IF
    
  END subroutine sbm_get_row_pointer_real

  ! Here the cached values, which lie outside the bands are inserted.
  ! This is done, by extending the bandwidth to contain these entries. 
  SUBROUTINE sbm_commit_values_real(mat)
    TYPE(StoreBandedMatrixObjectReal) :: mat

    ! Local variables
    TYPE(MatrixEntryNodeData), pointer :: nodedata_ptr
    TYPE(MatrixEntryNodeData) :: mynode
    REAL, DIMENSION(:,:), POINTER :: newMatrix
    INTEGER, DIMENSION(2) :: max_bw, global_max_bw, coord
    INTEGER :: irow, icol, ierr, dist, row_index, col_index, old_NumberOfStoredRows

    CALL exit_if_factored(mat)
    
    ! 1. find maximal bandwidths in local list
    max_bw(1) = mat%UpperBandwidth
    max_bw(2) = mat%LowerBandwidth
    !PRINT*, "max_bw = ", max_bw

    call init_Iterator(mat%storage_list)
    do while (.not.isAtEnd(mat%storage_list))
       nodedata_ptr => get_next(mat%storage_list)
       coord = nodedata_ptr%coord
       dist = coord(2)-coord(1)
       IF (dist.GT.0) THEN
          max_bw(1) = MAX(dist, max_bw(1))
       ELSEIF (dist.LT.0) THEN
          max_bw(2) = MAX(-dist,max_bw(2))
       END IF
    END DO

    ! 2. we do a global reduction to get the global maximum
    CALL mpi_allreduce(max_bw, global_max_bw, 2, MPI_INTEGER, MPI_MAX, mat%PG%Communicator, ierr)

    ! 3. allocate the new local matrix if necessary
    IF ((global_max_bw(1).GT.mat%UpperBandwidth).OR.(global_max_bw(2).GT.mat%LowerBandwidth)) THEN
       ! Attention, newMatrix will not be deallocated, as the mat%localMatrix will point to it!
       old_NumberOfStoredRows = mat%NumberOfStoredRows
       mat%NumberOfStoredRows = global_max_bw(1)+global_max_bw(2)+1
       !PRINT*, "Allocating newMatrix(",mat%NumberOfStoredRows,mat%ColsPerBlock,")"
       IF (mat%isTransposed) THEN
          ALLOCATE(newMatrix(-global_max_bw(2):global_max_bw(1),mat%ColsPerBlock))
       ELSE
          ALLOCATE(newMatrix(-global_max_bw(1):global_max_bw(2),mat%ColsPerBlock))
       END IF
       newMatrix = CMPLX(0,0,kind(newMatrix))

       ! 4. copy old matrix into new one
       IF (mat%isTransposed) THEN
          DO icol=1,mat%ColsPerBlock
             DO irow=-mat%LowerBandwidth,mat%UpperBandwidth
                !PRINT*,irow,icol," row_index = ", global_max_bw(1)+1-mat%LowerBandwidth-1+irow
                newMatrix(irow,icol) = mat%localMatrix(irow,icol)
             END DO
          END DO
       ELSE
          DO icol=1,mat%ColsPerBlock
             DO irow=-mat%UpperBandwidth,mat%LowerBandwidth
                !PRINT*,irow,icol," row_index = ", global_max_bw(1)+1-mat%LowerBandwidth-1+irow
                newMatrix(irow,icol) = mat%localMatrix(irow,icol)
             END DO
          END DO
       END IF
       DEALLOCATE(mat%localMatrix)
       mat%localMatrix => newMatrix
       mat%UpperBandwidth = global_max_bw(1)
       mat%LowerBandwidth = global_max_bw(2)

       ! 5. add the entries of the storage_list
       DO WHILE (.NOT.(IsEmpty(mat%storage_list))) 
          CALL pop(mat%storage_list, mynode)

          IF (mat%isTransposed) THEN
             col_index = MOD(mynode%coord(1)-1,mat%ColsPerBlock)+1
             row_index = (mynode%coord(2)-mynode%coord(1))
          ELSE
             col_index = MOD(mynode%coord(2)-1,mat%ColsPerBlock)+1
             row_index = -(mynode%coord(2)-mynode%coord(1))
          END IF
          
          mat%localMatrix(row_index, col_index) = mynode%value
       END DO
    END IF

  END SUBROUTINE sbm_commit_values_real

  SUBROUTINE sbm_set_zero_real(this)
    type(StoreBandedMatrixObjectReal) :: this
    
    ! deallocate all local storage and set the diagonal to zero

    IF (ASSOCIATED(this%localMatrix))  DEALLOCATE(this%localMatrix)
    IF (ASSOCIATED(this%factoredBand)) DEALLOCATE(this%factoredBand)
    IF (ASSOCIATED(this%auxData))      DEALLOCATE(this%auxData)
    IF (ASSOCIATED(this%ipiv))         DEALLOCATE(this%ipiv)

    this%UpperBandwidth=0
    this%LowerBandwidth=0
    this%NumberOfStoredRows=1
    ! allocate only the diagonal
    ALLOCATE(this%localMatrix(0:0,this%ColsPerBlock))
    ! set diagonal to zero
    this%localMatrix = 0.0
  END SUBROUTINE sbm_set_zero_real

  !> Tune the Matrix-Vector multiplication
  subroutine sbm_autotune_real(mat,vec,res)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    TYPE(StoreVectorObjectReal),INTENT(IN) :: vec
    TYPE(StoreVectorObjectReal),INTENT(INOUT) :: res

    ! Local variables
    integer :: ierr, min_index(1), iter,max_iter=1
    real(8) :: start_time,end_time
    real(8) :: measured_time(4)

    ! the first call is always slower, therefore we measure 
    ! the second call
    measured_time = 0.0D0
    do iter=1,max_iter
       MatVecMethod = USE_BLAS
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       start_time = MPI_Wtime()
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       end_time = MPI_Wtime()
       measured_time(1) = measured_time(1) + end_time-start_time
       
       MatVecMethod = F_INTRINSIC
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       start_time = MPI_Wtime()
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       end_time = MPI_Wtime()
       measured_time(2) = measured_time(2) + end_time-start_time
       
       MatVecMethod = BY_INDEX
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       start_time = MPI_Wtime()
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       end_time = MPI_Wtime()
       measured_time(3) = measured_time(3) + end_time-start_time
       
       MatVecMethod = BY_ROW
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       start_time = MPI_Wtime()
       call dot_multiply(mat,vec,res)
       call mpi_barrier(mat%PG%communicator,ierr)
       end_time = MPI_Wtime()
       measured_time(4) = measured_time(4) + end_time-start_time
    end do
    write(*,"(A)",advance='no') "For Sparse Matrix Vector Multiplication, we use "
    min_index = minloc(measured_time)
    select case (min_index(1))
    case(1) 
       MatVecMethod = USE_BLAS
       write(*,"(A)",advance='no') "USE_BLAS"
    case(2) 
       MatVecMethod = F_INTRINSIC
       write(*,"(A)",advance='no') "F_INTRINSIC"
    case(3) 
       MatVecMethod = BY_INDEX
       write(*,"(A)",advance='no') "BY_INDEX"
    case(4) 
       MatVecMethod = BY_ROW
       write(*,"(A)",advance='no') "BY_ROW"
    end select
    write(*,"(A)") " method."

  end subroutine sbm_autotune_real

  !> Transpose the storage format of matrix mat
  !! The mathematical structure of the matrix is not changed, only
  !! the storage format (row-major or column-major) is changed.
  !! \param mat the input matrix
  !! \param mat_t the output matrix in the other storage format
  subroutine sbm_transpose_storage_real(mat,mat_t)
    type(StoreBandedMatrixObjectReal), intent(IN) ::mat
    type(StoreBandedMatrixObjectReal), intent(INOUT) :: mat_t

    ! Local variables
    integer :: rank, icol, ierr,irow
    integer :: global_col, global_row
    REAL,allocatable,dimension(:,:) :: local_total_mat, local_total_mat_t

    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)

    if (mat%isTransposed) then
       print*,"transpose_storage only works for normal to transposed change"
    else
       if (.not.mat_t%isTransposed) then
          print*,"the target matrix of transpose_storage must have transposed format."
          stop
       end if
       ! mat is not transposed (that means it is in column-major order)
       
       ! simple algorithm (memory intensive)
       call allocate(mat_t,mat%LowerBandwidth,mat%UpperBandwidth)
       allocate(local_total_mat(-mat%UpperBandwidth:mat%LowerBandwidth,mat%NCols))
       allocate(local_total_mat_t(-mat_t%LowerBandwidth:mat_t%UpperBandwidth,mat_t%NRows))
       ! first gather the full matrix on all processes

       call mpi_gather(mat%localmatrix,size(mat%localmatrix),MPI_REAL_TYPE,&
            &local_total_mat,size(mat%localmatrix),MPI_REAL_TYPE,0,mat%PG%communicator,ierr)
       if (rank.eq.0) then
          do icol=1,mat_t%NRows
             do irow=-mat_t%LowerBandwidth,mat_t%UpperBandwidth
                global_row = icol
                global_col = icol+irow
                if ((global_col.ge.1).and.(global_col.le.mat_t%NRows)) then
                   local_total_mat_t(irow,icol) = local_total_mat(global_row-global_col,global_col)
                end if
             end do
          end do
       end if
       call mpi_scatter(local_total_mat_t,size(mat_t%localmatrix),MPI_REAL_TYPE,&
            &mat_t%localmatrix,size(mat_t%localmatrix),MPI_REAL_TYPE,0,mat%PG%communicator,ierr)
       deallocate(local_total_mat,local_total_mat_t)
    end if
  end subroutine sbm_transpose_storage_real

  !> Multiply BandedMatrix (in column-major storage) with Vector
  !! A StoreBandedMatrixObject in column-major (standard) storage format
  !! is multiplied with a StoreVectorObject. The result is a StoreVectorObject.
  !! \param mat the matrix
  !! \param vec the vector
  !! \param res the result vector, res = mat . vec
  SUBROUTINE sbm_dot_multiply_real(mat, vec, res)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    TYPE(StoreVectorObjectReal),INTENT(IN) :: vec
    TYPE(StoreVectorObjectReal),INTENT(INOUT) :: res

    !Local variables
    COMPLEX,PARAMETER :: complex_zero=(0,0)
    REAL, DIMENSION(1:mat%NumberOfStoredRows+mat%ColsPerBlock-1) :: temp_res
    REAL, DIMENSION(1:mat%UpperBandwidth) :: from_next
    REAL, DIMENSION(1:mat%LowerBandwidth) :: from_previous
    INTEGER :: rank, ierr, pe_dest, pe_source, tag, t_innerlast,t_innerfirst, innerlast, icol
    !INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: n_procs, npe_upper_sum, npe_lower_sum
    INTEGER, DIMENSION(1) :: dims, coords
    LOGICAL, DIMENSION(1) :: periods

    CALL exit_if_factored(mat)

    IF (mat%isTransposed) THEN
       ! A similar routine handles the case with transposed matrix
       CALL sbm_dot_multiply_transposed_real(mat,vec,res)
    ELSE
       CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)
       CALL mpi_cart_get(mat%PG%communicator, 1, dims, periods, coords,ierr)
       CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)

!       temp_res = CMPLX(0.0,0.0,KIND(complex_kind))
       temp_res = 0.0

       DO icol=1,mat%ColsPerBlock
          temp_res(icol:icol+mat%NumberOfStoredRows-1) = &
               & temp_res(icol:icol+mat%NumberOfStoredRows-1)&
               & +vec%localVector(icol)*mat%localMatrix(:,icol)
       END DO

       ! now we need a reduction over the processes

       ! the result vector has locally the following parts
       ! 1:UpperBandwidth, parts of the summation, which have to be sent to the previous process
       ! UpperBandwidth+1:UpperBandwidth+RowsPerBlock, the partial sums, which remain on the process
       ! UpperBandwidth+RowsPerBlock+1:UpperBandwidth+RowsPerBlock+LowerBandwidth, parts of the
       !  summation, which have to be sent to the next process

       ! innerlast is the index of the last inner point in the local res%localMatrix array
       innerlast = res%RowsPerBlock
       ! t_innerlast is the index of the last inner point in the local temp_res array
       t_innerfirst = mat%UpperBandwidth+1
       t_innerlast = mat%UpperBandwidth+res%RowsPerBlock

       IF (n_procs.GT.1) THEN
          ! exchange of the upper bandwidth
          IF (.NOT.periods(1)) THEN
             ! direction is not periodic
             ! calculate how many processes are involved in the sum-exchange
             npe_upper_sum = mat%UpperBandwidth/res%RowsPerBlock
             IF (MODULO(mat%UpperBandwidth,res%RowsPerBlock).NE.0) npe_upper_sum = npe_upper_sum + 1
             npe_lower_sum = mat%LowerBandwidth/res%RowsPerBlock
             IF (MODULO(mat%LowerBandwidth,res%RowsPerBlock).NE.0) npe_lower_sum = npe_lower_sum + 1

             IF (mat%UpperBandwidth.GT.0) THEN
                IF (npe_upper_sum.EQ.1) THEN
                   ! just one process involed in communication, we can do all sends and receives
                   ! in parallel
                   res%localVector = 0.0
                   tag = 396
                   CALL mpi_cart_shift(mat%PG%communicator,0,-1,pe_source,pe_dest,ierr)
                  ! from_next = CMPLX(0,0,kind(from_next))
                   from_next = 0.0
                   CALL mpi_sendrecv(temp_res,mat%UpperBandwidth,MPI_REAL_TYPE,pe_dest,tag,&
                        & from_next,mat%UpperBandwidth, &
                        & MPI_REAL_TYPE, pe_source, tag, &
                        & mat%PG%communicator, MPI_STATUS_IGNORE, ierr)
                   temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) = &
                        &temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) &
                        &+ from_next
                ELSE
                   ! we start sending at the last process and send always to the previous, then
                   ! add the result and send it to the previous, ...
                   tag = 339
                   IF (rank.EQ.n_procs-1) THEN
                      ! last process, only send
                      CALL mpi_send(temp_res,mat%UpperBandwidth,MPI_REAL_TYPE, &
                           & rank-1, tag, mat%PG%communicator,ierr)
                   ELSEIF (rank.eq.0) then
                      ! first process, only receive, add
                      CALL mpi_recv(from_next,mat%UpperBandwidth,MPI_REAL_TYPE, &
                           & rank+1, tag, mat%PG%communicator,MPI_STATUS_IGNORE,ierr)
                      temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) = &
                           &temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) &
                           &+ from_next
                   ELSE
                      ! all other processes, receive, add, send
                      call mpi_recv(from_next,mat%UpperBandwidth,MPI_REAL_TYPE, &
                           & rank+1, tag, mat%PG%communicator,MPI_STATUS_IGNORE,ierr)
                      temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) = &
                           &temp_res(res%RowsPerBlock+1:res%RowsPerBlock+mat%UpperBandwidth) &
                           &+ from_next
                      call mpi_send(temp_res,mat%UpperBandwidth,MPI_REAL_TYPE, &
                           & rank-1, tag, mat%PG%communicator,ierr)
                   END IF
                END IF
             END IF

             ! LowerBandwidth
             IF (mat%LowerBandwidth.GT.0) THEN
                IF (npe_lower_sum.EQ.1) THEN
                   tag = 407
                   CALL mpi_cart_shift(mat%PG%communicator,0,1,pe_source,pe_dest,ierr)
                   from_previous = 0
                   CALL mpi_sendrecv(temp_res(t_innerlast+1),mat%LowerBandwidth,MPI_REAL_TYPE,pe_dest,tag,&
                        & from_previous,mat%LowerBandwidth,MPI_REAL_TYPE, pe_source, tag, &
                        & mat%PG%communicator, MPI_STATUS_IGNORE, ierr)
                   temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) = &
                        &temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) &
                        &+ from_previous
                ELSE
                   ! we start sending at the first process and send always to the next, then
                   ! add the result and send it to the next, ...
                   tag = 366
                   IF (rank.EQ.0) THEN
                      ! first process, only send
                      CALL mpi_send(temp_res(t_innerlast+1),&
                           & mat%LowerBandwidth,MPI_REAL_TYPE, &
                           & rank+1, tag, mat%PG%communicator,ierr)
                   ELSEIF (rank.eq.n_procs-1) then
                      ! last process, only receive, add
                      CALL mpi_recv(from_previous,mat%LowerBandwidth,MPI_REAL_TYPE, &
                           & rank-1, tag, mat%PG%communicator,MPI_STATUS_IGNORE,ierr)
                      temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) = &
                           &temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) &
                           &+ from_previous
                   ELSE
                      ! all other processes, receive, add, send
                      call mpi_recv(from_previous,mat%LowerBandwidth,MPI_REAL_TYPE, &
                           & rank-1, tag, mat%PG%communicator,MPI_STATUS_IGNORE,ierr)
                      temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) = &
                           &temp_res(t_innerfirst:t_innerfirst+mat%LowerBandwidth-1) &
                           &+ from_previous
                      CALL mpi_send(temp_res(t_innerlast+1),&
                           & mat%LowerBandwidth,MPI_REAL_TYPE, &
                           & rank+1, tag, mat%PG%communicator,ierr)
                   END IF
                END IF
             END IF
          ELSE
             ! direction is periodic
             ! we cannot start at the last or first process immediately but we have to 
             ! sum up the entries at the start process
             STOP "periodic boundary condition for banded matrices not yet implemented"
          END IF
       END IF
       ! now temp_res contains the correct sum on all processes. The inner points have
       ! to be copied to the res%localMatrix
       res%localVector = temp_res(t_innerfirst:t_innerlast)
    END IF

  END SUBROUTINE sbm_dot_multiply_real

  !> Multiply BandedMatrix (in row-major storage) with Vector
  !! A StoreBandedMatrixObject in row-major (transposed) storage format
  !! is multiplied with a StoreVectorObject. The result is a StoreVectorObject.
  !! \param mat the matrix
  !! \param vec the vector
  !! \param res the result vector, res = mat . vec
  SUBROUTINE sbm_dot_multiply_transposed_real(mat,vec,res)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    TYPE(StoreVectorObjectReal),INTENT(IN) :: vec
    TYPE(StoreVectorObjectReal),INTENT(INOUT) :: res

    !Local variables
    INTEGER :: rank, n_procs, ierr,irow,row_start, mat_col_start, mat_col_end, global_row, i
    integer :: s_global_col, e_global_col, start_row, end_row, LBW, UBW
    !INTEGER, DIMENSION(1) :: dims, coords
    !LOGICAL, DIMENSION(1) :: periods
    REAL, DIMENSION(1:vec%NRows),target :: localfullvec
    REAL, DIMENSION(:), pointer :: ptr_vec, row_ptr
    REAL :: prod !, cdotu

    CALL exit_if_factored(mat)

    if (mat%LowerBandwidth+mat%UpperBandwidth+1.gt.mat%NRows) then
       MatVecMethod  = F_INTRINSIC
    else
       MatVecMethod = BY_INDEX
    end if

    LBW = mat%LowerBandwidth
    UBW = mat%UpperBandwidth
    ! We know that the matrix is transposed. 
    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)
    !CALL mpi_cart_get(mat%PG%communicator, 1, dims, periods, coords,ierr)
    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)
    
    ! First we gather the vector vec on all ranks.
    IF (n_procs.EQ.1) THEN
       ptr_vec => vec%localVector
    ELSE
       CALL get_global_vector_locally(vec,localfullvec)
       ptr_vec => localfullvec
    END IF
    ! localfullvec contains the full vector vec. This is now multiplied with
    ! the matrix

    !PERFON('sbm_dm')
    select case (MatVecMethod) 
    case(USE_BLAS)
       DO irow=1,mat%ColsPerBlock
          global_row = rank*mat%ColsPerBlock+irow
          call mat_get_row_pointer(mat,global_row,row_ptr,s_global_col,e_global_col)
stop 'fix implementation here as cdotu is not supposed to work with real numbers'
!          prod = cdotu(e_global_col-s_global_col+1,row_ptr,1,ptr_vec(s_global_col),1)
          
          CALL set_value(res,global_row, prod)
       end DO

    case(F_INTRINSIC)
       DO irow=1,mat%ColsPerBlock
          global_row = rank*mat%ColsPerBlock+irow
          call mat_get_row_pointer(mat,global_row,row_ptr,s_global_col,e_global_col)
          prod = dot_product(row_ptr,ptr_vec(s_global_col:e_global_col))
          CALL set_value(res,global_row, prod)
       end DO

    case(BY_INDEX)
       start_row = rank*mat%ColsPerBlock+1
       end_row   = (rank+1)*mat%ColsPerBlock
       
       if (start_row.le.mat%LowerBandwidth) then
          if (end_row.le.vec%NRows-mat%UpperBandwidth) then
             ! case 1
             do global_row = start_row, min(mat%LowerBandwidth,end_row)
                prod = mat%localMatrix(1-global_row,global_row-(start_row-1))*ptr_vec(1)
                DO i=1,mat%UpperBandwidth-1+global_row
                   prod = prod + mat%localMatrix(1-global_row+i,global_row-(start_row-1))*ptr_vec(1+i)
                END DO
                CALL set_value(res,global_row, prod)
             end do
             do global_row = min(mat%LowerBandwidth,end_row)+1,end_row
                mat_col_start = -mat%LowerBandwidth
                row_start = global_row - mat%LowerBandwidth
                mat_col_end = min(vec%NRows-global_row,mat%UpperBandwidth)
                
                prod = mat%localMatrix(-mat%LowerBandwidth,global_row-(start_row-1))*ptr_vec(row_start)
                DO i=1,mat_col_end-mat_col_start
                   prod = prod + mat%localMatrix(-mat%LowerBandwidth+i,global_row-(start_row-1))*ptr_vec(row_start+i)
                END DO
                CALL set_value(res,global_row, prod)
             end do
          else
             ! case 2
             ! This is mainly the case for n_procs = 1
             DO global_row = start_row,LBW
                !mat_col_start = 1-global_row
                !row_start = 1
                !mat_col_end = mat%UpperBandwidth
                
                prod = mat%localMatrix(1-global_row,global_row-start_row+1)*ptr_vec(1)
                DO i=1,UBW-(1-global_row)
                   prod = prod + mat%localMatrix(1-global_row+i,global_row-start_row+1)*ptr_vec(1+i)
                END DO
                CALL set_value(res,global_row, prod)
             end DO

             DO global_row = LBW+1,mat%NRows-UBW
                !mat_col_start = -mat%LowerBandwidth
                !row_start = global_row-mat%LowerBandwidth
                !mat_col_end = mat%UpperBandwidth

                !res%localvector(global_row) = dot_product(&
                !res%localvector(global_row-(start_row-1)) = dot_product(&
                !     & conjg(mat%localMatrix(-LBW:UBW,global_row-(start_row-1))),&
                !     & ptr_vec(global_row-LBW:global_row+UBW))

                prod = mat%localMatrix(-mat%LowerBandwidth,global_row-(start_row-1)) &
                     & *ptr_vec(global_row-mat%LowerBandwidth)
                DO i=1,mat%UpperBandwidth+mat%LowerBandwidth
                   prod = prod + mat%localMatrix(-mat%LowerBandwidth+i,global_row-(start_row-1)) &
                        & *ptr_vec(global_row-mat%LowerBandwidth+i)
                END DO
                res%localvector(global_row-(start_row-1)) = prod
                !CALL set_value(res,global_row, prod)
             end DO

             DO global_row = mat%NRows-UBW+1,end_row
                !mat_col_start = -mat%LowerBandwidth
                row_start = global_row-LBW
                !mat_col_end = vec%NRows-global_row
                
                prod = mat%localMatrix(-LBW,global_row-(start_row-1))*ptr_vec(row_start)
                DO i=1,vec%NRows-global_row+LBW
                   prod = prod + mat%localMatrix(-LBW+i,global_row-(start_row-1))*ptr_vec(row_start+i)
                END DO
                CALL set_value(res,global_row, prod)
             end DO
          end if
       else
          if (end_row.le.vec%NRows-mat%UpperBandwidth) then
             ! case 3
             DO global_row = start_row,end_row
                mat_col_start = max(1-global_row,-mat%LowerBandwidth)
                row_start = global_row+mat_col_start
                
                mat_col_end = min(vec%NRows-global_row,mat%UpperBandwidth)
                
                prod = mat%localMatrix(mat_col_start,global_row-(start_row-1))*ptr_vec(row_start)
                DO i=1,mat_col_end-mat_col_start
                   prod = prod + mat%localMatrix(mat_col_start+i,global_row-(start_row-1))*ptr_vec(row_start+i)
                END DO
                CALL set_value(res,global_row, prod)
             end DO
          else
             ! case 4
             DO global_row = start_row,min(mat%NRows-mat%UpperBandwidth,end_row)
                mat_col_start = max(1-global_row,-mat%LowerBandwidth)
                row_start = global_row+mat_col_start
                
                mat_col_end = min(vec%NRows-global_row,mat%UpperBandwidth)
                
                prod = mat%localMatrix(mat_col_start,global_row-(start_row-1))*ptr_vec(row_start)
                DO i=1,mat_col_end-mat_col_start
                   prod = prod + mat%localMatrix(mat_col_start+i,global_row-(start_row-1))*ptr_vec(row_start+i)
                END DO
                CALL set_value(res,global_row, prod)
             end DO
             DO global_row = max(start_row,mat%NRows-mat%UpperBandwidth),end_row
                mat_col_start = max(1-global_row,-mat%LowerBandwidth)
                row_start = global_row+mat_col_start
                
                mat_col_end = min(vec%NRows-global_row,mat%UpperBandwidth)
                
                prod = mat%localMatrix(mat_col_start,global_row-(start_row-1))*ptr_vec(row_start)
                DO i=1,mat_col_end-mat_col_start
                   prod = prod + mat%localMatrix(mat_col_start+i,global_row-(start_row-1))*ptr_vec(row_start+i)
                END DO
                CALL set_value(res,global_row, prod)
             end DO
          end if
       end if

    case(BY_ROW)
       DO irow=1,mat%ColsPerBlock
          global_row = rank*mat%ColsPerBlock+irow
          call mat_get_row_pointer(mat,global_row,row_ptr,s_global_col,e_global_col)
          prod = cmplx(0.0,0.0)
          do i=s_global_col,e_global_col
             prod = prod + row_ptr(i-s_global_col+1)*ptr_vec(i)
          end do
          CALL set_value(res,global_row, prod)
       end DO
    end select
    !PERFOFF
  END SUBROUTINE sbm_dot_multiply_transposed_real

  subroutine sbm_matmat_dot_multiply_real(mat,mat2,res,transA_in, transB_in)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat, mat2
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in     
    CHARACTER, intent(IN), optional :: transB_in     

    !Local variables
    CHARACTER :: transA='N', transB='N'

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif
    if (present(transB_in)) then
       transB = transB_in
    else
       transB = 'N'
    endif

    CALL exit_if_factored(mat)
    CALL exit_if_factored(mat2)
    
    ! We assume at the moment, that the result is always in row-major order.
    !if (.not.res%isTransposed) then
    !   print*, "In all dot_multiply functions of storebandedmatrix, the result must be in transposed storage."
    !   stop
    !end if

    if ( res%isTransposed ) then
       if ( ((mat%isTransposed).and.(transA.eq.'N')) .or. &
            & ( (.not.mat%isTransposed).and.((transA.eq.'C').or.(transA.eq.'T')) ) ) then
          ! we use the transposed structure of the storage of mat
          if ( ((mat2%isTransposed).and.(transB.eq.'N')) .or. &
               & ( (.not.mat2%isTransposed).and.((transB.eq.'C').or.(transB.eq.'T')) ) ) then
             ! do transp*transp -> transp
             call sbm_matmat_TTT_dot_multiply_real(mat,mat2,res,transA,transB)
          else
             ! do transp*normal -> transp
             call sbm_matmat_TNT_dot_multiply_real(mat,mat2,res,transA,transB)
          end if
       else
          ! here the storage is used in column-major order either with
          ! .not.mat%isTransposed (and transA.eq.'N') or with a combination of 
          ! mat%isTransposed and transA.eq.'C' or 'T'
          if ( ((mat2%isTransposed).and.(transB.eq.'N')) .or. &
               & ( (.not.mat2%isTransposed).and.((transB.eq.'C').or.(transB.eq.'T')) ) ) then
             ! we use the transposed structure of the storage of mat2
             ! do normal*transp -> transp
             call sbm_matmat_NTT_dot_multiply_real(mat,mat2,res,transA,transB)
          else
             ! here the storage of mat2 is used in column-major order either with
             ! .not.mat2%isTransposed (and transB.eq.'N') or with a combination of 
             ! mat2%isTransposed and transB.eq.'C' or 'T'
             ! do normal*normal -> transp
             call sbm_matmat_NNT_dot_multiply_real(mat,mat2,res,transA,transB)
          end if
       end if
    else
       if ( ((mat%isTransposed).and.(transA.eq.'N')) .or. &
            & ( (.not.mat%isTransposed).and.((transA.eq.'C').or.(transA.eq.'T')) ) ) then
          ! we use the transposed structure of the storage of mat
          if ( ((mat2%isTransposed).and.(transB.eq.'N')) .or. &
               & ( (.not.mat2%isTransposed).and.((transB.eq.'C').or.(transB.eq.'T')) ) ) then
             ! do transp*transp -> normal
             call sbm_matmat_TTN_dot_multiply_real(mat,mat2,res,transA,transB)
          else
             ! do transp*normal -> normal
             call sbm_matmat_TNN_dot_multiply_real(mat,mat2,res,transA,transB)
          end if
       else
          ! here the storage is used in column-major order either with
          ! .not.mat%isTransposed (and transA.eq.'N') or with a combination of 
          ! mat%isTransposed and transA.eq.'C' or 'T'
          if ( ((mat2%isTransposed).and.(transB.eq.'N')) .or. &
               & ( (.not.mat2%isTransposed).and.((transB.eq.'C').or.(transB.eq.'T')) ) ) then
             ! we use the transposed structure of the storage of mat2
             ! do normal*transp -> normal
             call sbm_matmat_NTN_dot_multiply_real(mat,mat2,res,transA,transB)
          else
             ! here the storage of mat2 is used in column-major order either with
             ! .not.mat2%isTransposed (and transB.eq.'N') or with a combination of 
             ! mat2%isTransposed and transB.eq.'C' or 'T'
             ! do normal*normal -> normal
             call sbm_matmat_NNN_dot_multiply_real(mat,mat2,res,transA,transB)
          end if
       end if
    end if
  end subroutine sbm_matmat_dot_multiply_real

  SUBROUTINE sbm_matmat_TNT_dot_multiply_real(mat, mat2, res,transA, transB)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat, mat2
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    ! Local variables
    INTEGER :: rank, ierr, n_procs
    !INTEGER, DIMENSION(1) :: dims, coords
    !LOGICAL, dimension(1) :: periods
    REAL :: res_value
    INTEGER :: res_col_start, res_col_end, k, icol, irow, res_upperbandwidth, res_lowerbandwidth
    INTEGER :: global_row, k_start, k_end, res_col
    REAL, DIMENSION(:,:), ALLOCATABLE :: localfullmat

    if (transA.ne.'N') then
       print*,"TNT: First matrix with 'C' or 'T' and transp not implemented yet."
       stop
    end if

    if (transB.ne.'N') then
       print*,"TNT: Second matrix with 'C' or 'T' and not_transp not implemented yet."
       stop
    end if

    ! first delete the result matrix. If this is not done, and the matrix has been used before
    ! for example in a loop with larger bandwidth than the result in the actual routine will
    ! use, then we get wrong results.
    call set_zero(res)
    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)
    !    CALL mpi_cart_get(mat%PG%communicator, 1, dims, periods, coords,ierr)
    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)

    ! 1. calculate the upper and lower bandwidth of the result matrix
    res_upperbandwidth = mat%UpperBandwidth+mat2%UpperBandwidth
    res_lowerbandwidth = mat%LowerBandwidth+mat2%LowerBandwidth

    IF (n_procs.GT.1) THEN
       
       ! simplified transfer: copy whole mat2 to all ranks
       ALLOCATE(localfullmat(-mat2%UpperBandwidth:mat2%LowerBandwidth,mat2%NCols))
       CALL mpi_allgather(mat2%localMatrix,SIZE(mat2%localMatrix),MPI_REAL_TYPE,&
            & localfullmat,SIZE(mat2%localMatrix),MPI_REAL_TYPE,&
            & mat2%PG%communicator,ierr)

       ! now multiply the two matrices
       DO irow=1,mat%ColsPerBlock
          global_row = rank*mat%ColsPerBlock+irow

          res_col_start = MAX(1,global_row-res_lowerbandwidth)
          res_col_end   = MIN(res%NCols,global_row+res_upperbandwidth)

          DO icol=res_col_start,res_col_end
             k_start = MAX(global_row-mat%LowerBandwidth,icol-mat2%UpperBandwidth,1)
             k_end   = MIN(global_row+mat%UpperBandwidth,icol+mat2%LowerBandwidth,mat%NCols)
             res_value = 0
             DO k=k_start, k_end
                !res_value = res_value + mat_get_value(mat,global_row,k)*localfullmat(k-icol,icol)
                res_value = res_value + mat%localMatrix(k-global_row,irow)*localfullmat(k-icol,icol)
             END DO
             CALL set_value(res,global_row,icol,res_value)
          END DO
       END DO
       CALL commit_values(res)
       DEALLOCATE(localfullmat)

       ! 2. determine which columns of mat2 are to be transferred
       ! last row on a rank
       !global_row = (rank+1)*mat%ColsPerBlock-1
       ! last column not equal zero in the result matrix is
       !max_col_index = global_row+res_upperbandwidth
       ! on which process do we have this column
       !col_on_proc = max_col_index/mat%ColsPerBlock
       ! so we need all columns from rank+1 to col_on_procs
       ! allocate the recv buffer
       !n_recv_cols = col_on_proc-(rank+1)+1)*mat%ColsPerBlock
       !ALLOCATE(recvbuf(mat2%UpperBandwidth:mat2%LowerBandwidth,mat2%ColsPerBlock))
       !DO iproc=1,col_on_proc
       !   CALL mpi_recv(recvbuf(:,)
          
       ! 3. start the transfer with non blocking send and receive

       ! 4. In parallel, compute the local entries of the result matrix

       ! 5. wait until transfer ready

       ! 6. then compute the missing entries of the result matrix
       
    ELSE
       ! matrices are not parallelized
       DO global_row = 1, res%NRows
          res_col_start = MAX(1,global_row-res_lowerbandwidth)
          res_col_end   = MIN(res%NCols,global_row+res_upperbandwidth)

          DO res_col = res_col_start,res_col_end
             k_start = MAX(global_row-mat%LowerBandwidth,res_col-mat2%UpperBandwidth,1)
             k_end   = MIN(global_row+mat%UpperBandwidth,res_col+mat2%LowerBandwidth,mat%NCols)

           !  res_value = CMPLX(0,0,kind(res_value))
             res_value = 0.0
             DO k=k_start,k_end
                !res_value = res_value + mat_get_value(mat,global_row,k)*mat_get_value(mat2,k,res_col)
                res_value = res_value + mat%localMatrix(k-global_row,global_row)*mat2%localMatrix(k-res_col,res_col)
                ! with direct access to localMatrix this could be accelerated.
                !res_value = res_value + mat_get_value(mat,global_row,k)*mat_get_value(mat2,k,res_col)
             END DO
             CALL set_value(res,global_row,res_col,res_value)
          END DO
       END DO
       call commit_values(res)
    END IF

  END SUBROUTINE sbm_matmat_TNT_dot_multiply_real

  !>Multiplication of a BandedMatrix with a BandedMatrix
  !!
  !! This is the routine where the real work is done. A StoreBandedMatrixObject
  !! is multiplied with another StoreBandedMatrixObject to give a StoreBandedMatrixObject.
  !! All BandedMatrix object have to be in transposed storage format (row-major order).
  !! \param mat the banded matrix
  !! \param mat2 the second banded matrix
  !! \param res  the banded result matrix
  !! \param transA_in operation on mat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  !! \param transB_in operation on mat2, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  subroutine sbm_matmat_TTT_dot_multiply_real(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    type(StoreBandedMatrixObjectReal),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    !Local variables
    REAL, dimension(:,:), pointer :: ptr_next_recv, ptr_next_compute, temp_ptr
    !complex, dimension(:), pointer :: ptr_fast_first, ptr_fast_second
    REAL, dimension(:,:), allocatable, target :: temp1, temp2
    REAL :: sum_value !, cdotu
    integer :: ierr,send_request, recv_request, tag, n_procs, rank
    !integer :: status(MPI_STATUS_SIZE)
    integer :: i_pe,i_block,i_block_recv,i_block_send
    integer :: global_res_row, global_res_col, global_index_start,global_index_end
    integer :: s_colindexset_mat, e_colindexset_mat,s_rowindexset_mat2,e_rowindexset_mat2
    integer :: s_intersection, e_intersection
    integer :: s_global_col, e_global_col, loop_length
    integer :: RowsPerBlock, res_lbw, res_ubw, local_row,local_col
    REAL, dimension(:),pointer :: first_mat_row
    !complex, dimension(1:res%NCols) :: result_row

    if ((transA.ne.'N').or.(transB.ne.'N')) then
       print*,"At the moment, sbm_matmat_TTT_dot_multiply only works with 'N' for both matrices."
       !call tracebackqq()
       stop
    end if

    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)
    call mpi_comm_rank(mat%PG%communicator, rank, ierr)


    ! allocate the two temporary arrays
    allocate(temp1(-mat2%LowerBandwidth:mat2%UpperBandwidth,mat2%ColsPerBlock))
    allocate(temp2(-mat2%LowerBandwidth:mat2%UpperBandwidth,mat2%ColsPerBlock))

    call set_zero(res)
    RowsPerBlock = mat%NRows/mat%PG%NProcs
    res_lbw = mat%LowerBandwidth+mat2%LowerBandwidth
    res_ubw = mat%UpperBandwidth+mat2%UpperBandwidth
    ! some arbitrary tag, just for debug purposes
    tag = 9784
    ! loop over all processes
    ! we start on all processes with the actual process
    ptr_next_compute => mat2%localMatrix
    ptr_next_recv => temp1
    do i_pe=0,n_procs-1
       ! calculate the running index, which starts at the rank
       ! and wraps around at n_procs
       ! i_block runs from 0 to n_procs-1, with wrapping at n_procs
       i_block = modulo(rank+i_pe,n_procs)
       
       i_block_recv = modulo(i_block+1,n_procs)
       i_block_send = modulo(rank-i_pe-1+n_procs,n_procs)

       !write(*,"(I3,A,4I3)") rank," : ",i_pe,i_block,i_block_recv,i_block_send

       ! do not wait for the first iteration, as nothing has been sent
       if (i_block.ne.rank) then
          ! wait for the received block
          !print*,rank,": Waiting for recv to finalize"
          call mpi_wait(recv_request,MPI_STATUS_IGNORE,ierr)
          ! swap the pointers
          temp_ptr => ptr_next_recv
          ptr_next_recv => ptr_next_compute
          ptr_next_compute => temp_ptr
          ! wait for the sended block
          !print*,rank,": Waiting for send to finalize"
          call mpi_wait(send_request,MPI_STATUS_IGNORE,ierr)
       end if
       ! initiate non-blocking receive for the next block from some following pe
       if (i_block_recv.ne.rank) then
          !print*,rank,": Starting irecv from ",i_block_recv
          call mpi_irecv(ptr_next_recv,mat2%NumberOfStoredRows*mat2%ColsPerBlock,MPI_REAL_TYPE,&
               &i_block_recv,tag,mat2%PG%communicator,recv_request,ierr)
       end if
       ! initiate non-blocking send of the local array to some previous pe
       if (i_block_send.ne.rank) then
          !print*,rank,": Starting isend to ",i_block_send
          call mpi_isend(mat2%localMatrix,mat2%NumberOfStoredRows*mat2%ColsPerBlock,MPI_REAL_TYPE,&
               &i_block_send,tag,mat2%PG%communicator,send_request,ierr)
       end if
       ! do computation
       ! multiply mat%localMatrix with ptr_next_compute
       ! mat is a transposed matrix, also the ptr_next_compute
       ! ptr_next_compute is a subblock of the localMatrix of mat2

       global_index_start = i_block*RowsPerBlock+1
       global_index_end   = (i_block+1)*RowsPerBlock
       !print*,rank,i_block,": global_index = ",global_index_start,global_index_end
       do global_res_row=rank*RowsPerBlock+1,(rank+1)*RowsPerBlock
          call mat_get_row_pointer(mat,global_res_row,first_mat_row,s_global_col,e_global_col)
          !print*,first_mat_row
          do global_res_col=max(global_res_row-res_lbw,1),min(global_res_row+res_ubw,res%NCols)
             !write(*,"(A,2I4)") "global_res = ", global_res_row,global_res_col

             ! calculating the global index set
             ! [s_global,e_global] is already intersection of colindexset_mat and I_block
             s_colindexset_mat = global_res_row-mat%LowerBandwidth
             e_colindexset_mat = global_res_row+mat%UpperBandwidth
             s_rowindexset_mat2= global_res_col-mat2%UpperBandwidth
             e_rowindexset_mat2= global_res_col+mat2%LowerBandwidth
             !print*,"colindexset = ",s_colindexset_mat,e_colindexset_mat
             !print*,"rowindexset = ",s_rowindexset_mat2,e_rowindexset_mat2

             ! now we need the intersection of I_block colindexset and rowindexset
             s_intersection = max(global_index_start,s_colindexset_mat,s_rowindexset_mat2)
             e_intersection = min(global_index_end, e_colindexset_mat,e_rowindexset_mat2)
             !write(*,"(4X,A,2I3)") "intersection = ", s_intersection,e_intersection
             loop_length = max(e_intersection-s_intersection+1,0)
             if (loop_length.gt.0) then

                ! the global index runs from s_intersection to e_intersection
                ! now convert to local indices
                !print*,first_mat_row(1+s_intersection-s_global_col:1+s_intersection-s_global_col+loop_length-1)
                local_row = lbound(ptr_next_compute,1)+global_res_col-s_intersection+mat2%LowerBandwidth
                local_col = mod(s_intersection-1,mat2%ColsPerBlock)+1
                !print*,"local_row,col = ",local_row,local_col
                !print*,ptr_next_compute(local_row,local_col)
                !loop_length = min(e_global_col,global_index_end)-max(s_global_col,global_index_start)+1
                !result_row(global_res_col) = cdotu(loop_length,first_mat_row(1+s_intersection-s_global),1,&
stop 'fix implementation here as cdotu is not supposed to work with real numbers'
!                sum_value = cdotu(loop_length,first_mat_row(1+s_intersection-s_global_col),1,&
!                     & ptr_next_compute(local_row,local_col),mat2%UpperBandwidth+mat2%LowerBandwidth)
                call add_value(res,global_res_row,global_res_col,sum_value)
             end if

             !ptr_local_column => ptr_next_compute((global_mat_col-1)*RowsPerBlock+1:global_mat_col*RowsPerBlock)
             !if (s_global_col.ge.global_index_start) then
             !   ptr_fast_first => first_mat_row(1:loop_length)
             !   ptr_fast_second => ptr_local_column(s_global_col-global_index_start+1:s_global_col-global_index_start+loop_length)
             !else
             !   ptr_fast_first => first_mat_row(global_index_start-s_global_col+1:global_index_start-s_global_col+loop_length)
             !   ptr_fast_second => ptr_local_column(1:loop_length)
             !end if

             ! to save memory, we can write directly into the matrix with the following line
             !sum_value = cdotu(loop_length,ptr_fast_first,1,ptr_fast_second,1)
             !call add_value(res,global_bmat_row,global_mat_col,sum_value)
             ! for better performance, we collect a whole row of the result matrix
             !result_row(global_mat_col) = cdotu(loop_length,ptr_fast_first,1,ptr_fast_second,1)
          end do
          ! for better performance, now add the whole result_row to the result matrix
          !call add_to_row(res,global_bmat_row, result_row)
       end do
       call commit_values(res)
       ! ---------------------- Ende des komplizierten Teils --------------

       if (i_block.eq.rank) then
          ! the first iteration is special, as ptr_next_compute points
          ! to localMatrix, but to do the pointer swapping, we have to 
          ! redirect ptr_next_compute to the other temporary storage.
          ptr_next_compute => temp2
       end if
    end do
    deallocate(temp1,temp2)
  end subroutine sbm_matmat_TTT_dot_multiply_real

  !>Multiplication of a BandedMatrix with a BandedMatrix
  !!
  !! This is the routine where the real work is done. A StoreBandedMatrixObject
  !! is multiplied with another StoreBandedMatrixObject to give a StoreBandedMatrixObject.
  !! The first matrix is now in column-major order, the two other matrices are
  !! stored in row-major order.
  !! \param mat the banded matrix
  !! \param mat2 the second banded matrix
  !! \param res  the banded result matrix
  !! \param transA_in operation on mat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  !! \param transB_in operation on mat2, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  subroutine sbm_matmat_NTT_dot_multiply_real(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    type(StoreBandedMatrixObjectReal),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    !Local variables
    REAL, dimension(:), pointer :: res_row_ptr, row_ptr
    REAL, dimension(:,:),pointer :: res_localmatrix_temp
    integer, dimension(:), allocatable :: recvcounts

    integer :: ierr, n_procs, rank
    integer :: first_LBW, first_UBW, sec_LBW, sec_UBW
    integer :: global_res_row, s_global_index,e_global_index
    integer :: s_colindexset_mat, e_colindexset_mat,s_colindexset_mat2,e_colindexset_mat2
    integer :: s_intersection, e_intersection, loop_length
    integer :: s_intersection2, e_intersection2, loop_length2
    integer :: s_global_col, e_global_col,s_res_col, e_res_col
    integer :: res_lbw, res_ubw, inner_index
    logical :: conjugate, conjugate2

    !integer :: zaxpy_calls

    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)
    call mpi_comm_rank(mat%PG%communicator, rank, ierr)

    ! test if we have to use a conjugated matrix or not
    if (transA.eq.'N') then
       ! first matrix is not Transposed
       conjugate=.false.
       first_LBW = mat%LowerBandwidth
       first_UBW = mat%UpperBandwidth
    else
       ! first matrix isTransposed, that means Lower and Upper
       ! Bandwidths are exchanged
       first_LBW = mat%UpperBandwidth
       first_UBW = mat%LowerBandwidth
       if (transA.eq.'C') then
          conjugate=.true.
       else
          conjugate=.false.
       end if
    end if

    if (transB.eq.'N') then
       ! second matrix isTransposed
       sec_LBW=mat2%LowerBandwidth
       sec_UBW=mat2%UpperBandwidth
       conjugate2 = .false.
    else
       ! second matrix .not.isTransposed
       sec_LBW = mat2%UpperBandwidth
       sec_UBW = mat2%LowerBandwidth
       if (transB.eq.'C') then
          conjugate2=.true.
       else
          conjugate2=.false.
       end if
    end if

    ! allocate the two temporary arrays

    res_lbw = min(first_LBW+sec_LBW,res%NCols)
    res_ubw = min(first_UBW+sec_UBW,res%NCols)
    !call set_zero(res)
    ! allocate the storage in res with the calculated bandwidths
    call allocate(res,res_lbw,res_ubw)
    if (n_procs.ne.1) then
       ! for communication we need a temporary global array
       ! if memory is an issue, this can be changed
       !write(*,"(A,I10,A,F8.2,A)") "Allocating res_localmatrix_temp with ",&
       !     &(res_ubw+res_lbw+1)*res%NRows," elements (",&
       !     &(res_ubw+res_lbw+1)*res%NRows*16/1024./1024.," MB)."
       !DEC$ ATTRIBUTES ALIGN: 128:: res_localmatrix_temp
       allocate(res_localmatrix_temp(-res_lbw:res_ubw,1:res%NRows))
       res_localmatrix_temp = cmplx(0.0,0.0)
    else
       res_localmatrix_temp => res%localMatrix
       !write(*,"(A,I10,A,F8.2,A)") "Using res_localmatrix_temp with ",size(res%localMatrix),&
       !     &" elements (",size(res%localMatrix)*16/1024./1024.," MB)."
    end if
    
    !zaxpy_calls = 0
    ! calculate the product on each processor with the localmatrices
    s_global_index = rank*mat%ColsPerBlock+1
    e_global_index = (rank+1)*mat%ColsPerBlock
    do global_res_row=1,res%NRows
       !write(*,"(A,I4)") "global_res_row = ",global_res_row

       ! calculating the global index set
       ! [s_global,e_global] is already intersection of colindexset_mat and I_block
       s_colindexset_mat = global_res_row-first_LBW
       e_colindexset_mat = global_res_row+first_UBW
       !print*,"colindexset = ",s_colindexset_mat,e_colindexset_mat
       !print*,"colindexset2 = ",s_colindexset_mat2,e_colindexset_mat2
       
       ! now we need the intersection of I_block colindexset and rowindexset
       s_intersection = max(s_global_index,s_colindexset_mat)
       e_intersection = min(e_global_index, e_colindexset_mat)
       !print*,"intersection = ",s_intersection,e_intersection
       !print*,"intersection2 = ",s_intersection2,e_intersection2

       loop_length = max(e_intersection-s_intersection+1,0)
       
       !if (loop_length2.gt.0) then
       ! ========= mat_get_row_pointer not working on res but on res_localmatrix_temp
       s_res_col = max(global_res_row - res_LBW,1)
       e_res_col = min(global_res_row + res_UBW,res%NCols)
       
       res_row_ptr  => res_localmatrix_temp(s_res_col-global_res_row:e_res_col-global_res_row,global_res_row)
       !write(*,"(A,2I3,A)") "res_row_ptr(",lbound(res_row_ptr,1),ubound(res_row_ptr,1),")"
       !print*,s_res_col-global_res_row,e_res_col-global_res_row
       !==========
       !call mat_get_row_pointer(res,global_res_row,res_row_ptr,s_res_col,e_res_col)
       !print*,"res_col = ",s_res_col,e_res_col,", res_row_ptr = ",res_row_ptr

       !zaxpy_calls = zaxpy_calls + loop_length

       do inner_index=s_intersection,e_intersection
          !write(*,"(5X,A,I4)") "inner_index = ",inner_index
          call mat_get_row_pointer(mat2,inner_index,row_ptr,s_global_col,e_global_col)
          !print*,"global_col = ",s_global_col,e_global_col,", row_ptr = ",row_ptr
          
          s_colindexset_mat2= inner_index-sec_LBW
          e_colindexset_mat2= inner_index+sec_UBW
          s_intersection2 = max(1,s_colindexset_mat2)
          e_intersection2 = min(mat2%NCols,e_colindexset_mat2)
          loop_length2 = max(e_intersection2-s_intersection2+1,0)

          !write(*,"(A,I4,A)") "zaxpy(",loop_length2,")"
          call caxpy(loop_length2,mat_get_value(mat,global_res_row,inner_index,transA),&
               &row_ptr,1,res_row_ptr(1+s_global_col-s_res_col),1)
       end do
    end do

    !write(*,"(A,I10,A)") "zaxpy was called ",zaxpy_calls," times."
    if (n_procs.ne.1) then
       ! Here now follows the communication of the result matrices. They have all to be reduced
       ! but different parts of the matrix to different processes.
       allocate(recvcounts(n_procs))
       recvcounts=res%NumberOfStoredRows*res%ColsPerBlock
       call mpi_reduce_scatter(res_localMatrix_temp,res%localMatrix,recvcounts,&
            & MPI_REAL_TYPE,MPI_SUM,res%PG%Communicator,ierr)
       deallocate(recvcounts)
       deallocate(res_localmatrix_temp)
    end if
  end subroutine sbm_matmat_NTT_dot_multiply_real

  subroutine sbm_matmat_NNT_dot_multiply_real(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    type(StoreBandedMatrixObjectReal),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    print*,"NNT_dot_multiply is not yet implemented."
    stop
  end subroutine sbm_matmat_NNT_dot_multiply_real

  SUBROUTINE sbm_matmat_TNN_dot_multiply_real(mat, mat2, res,transA, transB)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat, mat2
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    ! Local variables
    INTEGER :: rank, ierr, n_procs
    !INTEGER, DIMENSION(1) :: dims, coords
    !LOGICAL, dimension(1) :: periods
    REAL :: res_value
    INTEGER :: res_col_start, res_col_end, k, icol, irow, res_upperbandwidth, res_lowerbandwidth
    integer :: res_row_start, res_row_end
    INTEGER :: global_row, global_col, k_start, k_end, res_col
    REAL, DIMENSION(:,:), ALLOCATABLE :: localfullmat

    if (transA.ne.'N') then
       print*,"TNT: First matrix with 'C' or 'T' and transp not implemented yet."
       stop
    end if

    if (transB.ne.'N') then
       print*,"TNT: Second matrix with 'C' or 'T' and not_transp not implemented yet."
       stop
    end if

    ! first delete the result matrix. If this is not done, and the matrix has been used before
    ! for example in a loop with larger bandwidth than the result in the actual routine will
    ! use, then we get wrong results.
    call set_zero(res)
    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)
    !    CALL mpi_cart_get(mat%PG%communicator, 1, dims, periods, coords,ierr)
    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)

    ! 1. calculate the upper and lower bandwidth of the result matrix
    res_upperbandwidth = mat%UpperBandwidth+mat2%UpperBandwidth
    res_lowerbandwidth = mat%LowerBandwidth+mat2%LowerBandwidth

    IF (n_procs.GT.1) THEN
       
       ! simplified transfer: copy whole mat to all ranks
       ALLOCATE(localfullmat(-mat%LowerBandwidth:mat%UpperBandwidth,mat%NCols))
       CALL mpi_allgather(mat%localMatrix,SIZE(mat%localMatrix),MPI_REAL_TYPE,&
            & localfullmat,SIZE(mat%localMatrix),MPI_REAL_TYPE,&
            & mat%PG%communicator,ierr)


       ! each processor works on a column block of the whole matrix
       !do icol=1,res%ColsPerBlock
       !   global_col = rank*res%ColsPerBlock+icol

       !do global_row = 1,res%NRows
       !      res_value = CMPLX(0,0,KIND(res_value))
       !      do k=1,mat%NRows
                ! mat is stored in tranposed format
       !         res_value += mat(k,global_row)*mat2(k,icol)
       !      end do
       !      CALL set_value(res,global_row,global_col,res_value)
       !   end do
       !end do

       do icol=1,res%ColsPerBlock
          global_col = rank*res%ColsPerBlock+icol

          res_row_start = MAX(1,global_col-res_upperbandwidth)
          res_row_end   = MIN(res%NCols,global_col+res_lowerbandwidth)

          do global_row = res_row_start,res_row_end
             ! now determine the k range
             k_start = max(1,global_row-mat%LowerBandwidth,global_col-mat2%UpperBandwidth)
             k_end   = min(mat%NCols,global_row+mat%UpperBandwidth,global_col+mat2%LowerBandwidth)
             !write(*,"(4(A,I3))") "Calculating r(",global_row,",",global_col,&
             !     & "), k running from ",k_start," to ",k_end
             res_value = 0
             do k=k_start,k_end
                ! mat is stored in tranposed format
                res_value = res_value + &
                     &localfullmat(k-global_row,global_row)*mat_get_value(mat2,k,global_col)
             end do
             CALL set_value(res,global_row,global_col,res_value)
          end do
       end do

       CALL commit_values(res)
       DEALLOCATE(localfullmat)

       ! 2. determine which columns of mat2 are to be transferred
       ! last row on a rank
       !global_row = (rank+1)*mat%ColsPerBlock-1
       ! last column not equal zero in the result matrix is
       !max_col_index = global_row+res_upperbandwidth
       ! on which process do we have this column
       !col_on_proc = max_col_index/mat%ColsPerBlock
       ! so we need all columns from rank+1 to col_on_procs
       ! allocate the recv buffer
       !n_recv_cols = col_on_proc-(rank+1)+1)*mat%ColsPerBlock
       !ALLOCATE(recvbuf(mat2%UpperBandwidth:mat2%LowerBandwidth,mat2%ColsPerBlock))
       !DO iproc=1,col_on_proc
       !   CALL mpi_recv(recvbuf(:,)
          
       ! 3. start the transfer with non blocking send and receive

       ! 4. In parallel, compute the local entries of the result matrix

       ! 5. wait until transfer ready

       ! 6. then compute the missing entries of the result matrix
       
    ELSE
       ! matrices are not parallelized
       DO global_row = 1, res%NRows
          res_col_start = MAX(1,global_row-res_lowerbandwidth)
          res_col_end   = MIN(res%NCols,global_row+res_upperbandwidth)

          DO res_col = res_col_start,res_col_end
             k_start = MAX(global_row-mat%LowerBandwidth,res_col-mat2%UpperBandwidth,1)
             k_end   = MIN(global_row+mat%UpperBandwidth,res_col+mat2%LowerBandwidth,mat%NCols)

        !     res_value = CMPLX(0,0,kind(res_value))
             res_value = 0.0
             DO k=k_start,k_end
                !res_value = res_value + mat_get_value(mat,global_row,k)*mat_get_value(mat2,k,res_col)
                res_value = res_value + mat%localMatrix(k-global_row,global_row)*mat2%localMatrix(k-res_col,res_col)
                ! with direct access to localMatrix this could be accelerated.
                !res_value = res_value + mat_get_value(mat,global_row,k)*mat_get_value(mat2,k,res_col)
             END DO
             CALL set_value(res,global_row,res_col,res_value)
          END DO
       END DO
       call commit_values(res)
    END IF

  END SUBROUTINE sbm_matmat_TNN_dot_multiply_real

  !>Multiplication of a BandedMatrix with a BandedMatrix
  !!
  !! This is the routine where the real work is done. A StoreBandedMatrixObject
  !! is multiplied with another StoreBandedMatrixObject to give a StoreBandedMatrixObject.
  !! All BandedMatrix object have to be in transposed storage format (row-major order).
  !! \param mat the banded matrix
  !! \param mat2 the second banded matrix
  !! \param res  the banded result matrix
  !! \param transA_in operation on mat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  !! \param transB_in operation on mat2, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  subroutine sbm_matmat_TTN_dot_multiply_real(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    type(StoreBandedMatrixObjectReal),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    print*,"TTN_dot_multiply is not yet implemented."
    stop

  end subroutine sbm_matmat_TTN_dot_multiply_real

  !>Multiplication of a BandedMatrix with a BandedMatrix
  !!
  !! This is the routine where the real work is done. A StoreBandedMatrixObject
  !! is multiplied with another StoreBandedMatrixObject to give a StoreBandedMatrixObject.
  !! The first matrix is now in column-major order, the two other matrices are
  !! stored in row-major order.
  !! \param mat the banded matrix
  !! \param mat2 the second banded matrix
  !! \param res  the banded result matrix
  !! \param transA_in operation on mat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  !! \param transB_in operation on mat2, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  subroutine sbm_matmat_NTN_dot_multiply_real(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    type(StoreBandedMatrixObjectReal),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    print*,"NTN_dot_multiply is not yet implemented."
    stop

  end subroutine sbm_matmat_NTN_dot_multiply_real

  subroutine sbm_matmat_NNN_dot_multiply_real(mat,mat2,res,transA,transB)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    type(StoreBandedMatrixObjectReal),intent(IN) :: mat2
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN) :: transA
    CHARACTER, intent(IN) :: transB

    print*,"NNN_dot_multiply is not yet implemented."
    stop
  end subroutine sbm_matmat_NNN_dot_multiply_real

  !>Multiplication of a BandedMatrix with a Matrix
  !!
  !! This is the routine where the real work is done. A StoreBandedMatrixObject
  !! is multiplied with a StoreFullMatrixObject to give a StoreFullMatrixObject.
  !! The BandedMatrix has to be in transposed storage format, the full matrix
  !! in normal storage format (this is the only possibility at the moment). 
  !! The result will be in normal storage format.
  !! \param bmat the banded matrix
  !! \param mat  the full matrix
  !! \param res  the full result matrix
  !! \param transA_in operation on bmat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  !! \param transB_in operation on mat, 'C' (conjugate/transpose), 'T' (transpose), 'N' (do nothing),
  !! default 'N'
  subroutine sbm_dot_multiply_Banded_with_Full_real(bmat,mat,res,transA_in,transB_in)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: bmat
    type(StoreFullMatrixObjectReal),intent(IN) :: mat
    TYPE(StoreFullMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in     
    CHARACTER, intent(IN), optional :: transB_in     

    !Local variables
    CHARACTER :: transA='N', transB='N'
    real, dimension(:), pointer :: ptr_next_recv, ptr_next_compute, temp_ptr, ptr_local_column
    real, dimension(:), pointer :: ptr_fast_first, ptr_fast_second
    real, dimension(:), allocatable, target :: temp1, temp2
    real, dimension(:),pointer :: first_mat_row
    real, dimension(1:res%NCols) :: result_row
!    real :: cdotu
    integer :: ierr,send_request, recv_request, tag, n_procs, rank
    !integer :: status(MPI_STATUS_SIZE)
    integer :: i_pe,i_block,i_block_recv,i_block_send,global_bmat_row
    integer :: global_mat_col,global_index_start,global_index_end
    integer :: s_global_col, e_global_col, loop_length

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif
    if (present(transB_in)) then
       transB = transB_in
    else
       transB = 'N'
    endif

    if ((transA.ne.'N').or.(transB.ne.'N')) then
       print*,"At the moment, dot_multiply_Banded_with_Full only works with 'N' for both matrices."
       !call tracebackqq()
       stop
    end if

    ! if one of the matrices is LU decomposed, we skip, as the normal
    ! dot multiply does not work with LU decomposed matrices
    CALL exit_if_factored(bmat)
    CALL exit_if_factored(mat)

    IF (.NOT.(bmat%isTransposed)) then
       PRINT*,"dot_multiply for banded with full is only allowed with &
            &the storage schemes transp*not_transp ->not_transp"
       STOP
    END IF

    !call print_storage_details(mat)

    CALL mpi_comm_size(bmat%PG%communicator, n_procs, ierr)
    call mpi_comm_rank(bmat%PG%communicator, rank, ierr)


    if (mat%BlockrowsPerProcess.eq.1) then
       ! we have a distribution over the rows of the full matrix
       ! or over the columns of the full matrix, exactly one block
       ! per process.
       ! allocate the two temporary arrays
       allocate(temp1(mat%NBlocks*mat%Blocksize))
       allocate(temp2(mat%NBlocks*mat%Blocksize))
    else
       ! scalapack block distributed storage of the full matrix
       ! not implemented!
       print*,"dot_multiply for banded with full for a block distributed &
            &full matrix with NPCols.ne.1 .and. NPRows.ne.1 is not implemented."
       stop
    end if

    ! for test purposes, set res matrix to zero
    ! has to be changed for production, just for testing
    res%localMatrix=00

    ! some arbitrary tag, just for debug purposes
    tag = 9783
    ! loop over all processes
    ! we start on all processes with the actual process
    ptr_next_compute => mat%localMatrix
    ptr_next_recv => temp1
    do i_pe=0,n_procs-1
       ! calculate the running index, which starts at the rank
       ! and wraps around at n_procs
       ! i_block runs from 0 to n_procs-1, with wrapping at n_procs
       i_block = modulo(rank+i_pe,n_procs)
       
       i_block_recv = modulo(i_block+1,n_procs)
       i_block_send = modulo(rank-i_pe-1+n_procs,n_procs)

       !write(*,"(I3,A,4I3)") rank," : ",i_pe,i_block,i_block_recv,i_block_send

       ! do not wait for the first iteration, as nothing has been sent
       if (i_block.ne.rank) then
          ! wait for the received block
          !print*,rank,": Waiting for recv to finalize"
          call mpi_wait(recv_request,MPI_STATUS_IGNORE,ierr)
          ! swap the pointers
          temp_ptr => ptr_next_recv
          ptr_next_recv => ptr_next_compute
          ptr_next_compute => temp_ptr
          ! wait for the sended block
          !print*,rank,": Waiting for send to finalize"
          call mpi_wait(send_request,MPI_STATUS_IGNORE,ierr)
       end if
       ! initiate non-blocking receive for the next block from some following pe
       if (i_block_recv.ne.rank) then
          !print*,rank,": Starting irecv from ",i_block_recv
          call mpi_irecv(ptr_next_recv,mat%NBlocks*mat%Blocksize,MPI_REAL_TYPE,&
               &i_block_recv,tag,bmat%PG%communicator,recv_request,ierr)
       end if
       ! initiate non-blocking send of the local array to some previous pe
       if (i_block_send.ne.rank) then
          !print*,rank,": Starting isend to ",i_block_send
          call mpi_isend(mat%localMatrix,mat%NBlocks*mat%Blocksize,MPI_REAL_TYPE,&
               &i_block_send,tag,bmat%PG%communicator,send_request,ierr)
       end if
       ! do computation
       ! multiply bmat%localMatrix with ptr_next_compute
       ! bmat is a transposed matrix
       ! ptr_next_compute is a block of the full matrix
       global_index_start = i_block*mat%RowsPerBlock+1
       global_index_end   = (i_block+1)*mat%RowsPerBlock
       do global_bmat_row=rank*mat%RowsPerBlock+1,(rank+1)*mat%RowsPerBlock
          call mat_get_row_pointer(bmat,global_bmat_row,first_mat_row,s_global_col,e_global_col)
          do global_mat_col=1,mat%NCols
             ptr_local_column => ptr_next_compute((global_mat_col-1)*mat%RowsPerBlock+1:global_mat_col*mat%RowsPerBlock)
             loop_length = min(e_global_col,global_index_end)-max(s_global_col,global_index_start)+1
             if (s_global_col.ge.global_index_start) then
                ptr_fast_first => first_mat_row(1:loop_length)
                ptr_fast_second => ptr_local_column(s_global_col-global_index_start+1:s_global_col-global_index_start+loop_length)
             else
                ptr_fast_first => first_mat_row(global_index_start-s_global_col+1:global_index_start-s_global_col+loop_length)
                ptr_fast_second => ptr_local_column(1:loop_length)
             end if

             ! to save memory, we can write directly into the matrix with the following line
             !sum_value = cdotu(loop_length,ptr_fast_first,1,ptr_fast_second,1)
             !call add_value(res,global_bmat_row,global_mat_col,sum_value)
             ! for better performance, we collect a whole row of the result matrix
             if (loop_length.gt.0) then
stop 'fix implementation here as cdotu is not supposed to work with real numbers'
!                result_row(global_mat_col) = cdotu(loop_length,ptr_fast_first,1,ptr_fast_second,1)
             else
                result_row(global_mat_col) = cmplx(0.0D0,0.0D0)
             end if
          end do
          ! for better performance, now add the whole result_row to the result matrix
          call add_to_row(res,global_bmat_row, result_row)
       end do
       call commit_values(res)

       if (i_block.eq.rank) then
          ! the first iteration is special, as ptr_next_compute points
          ! to localMatrix, but to do the pointer swapping, we have to 
          ! redirect ptr_next_compute to the other temporary storage.
          ptr_next_compute => temp2
       end if
    end do
    deallocate(temp1,temp2)
  end subroutine sbm_dot_multiply_Banded_with_Full_real

  !> Calculate the minimal k and maximal k index on process proc
  SUBROUTINE get_minmax_k_real(mat,rank,min_k,max_k)
    TYPE(StoreBandedMatrixObjectReal), intent(IN) :: mat
    INTEGER, intent(IN) :: rank
    INTEGER, INTENT(OUT) :: min_k, max_k

    ! Local variables
    INTEGER :: icol, irow, this_min_k, this_max_k, n_procs,ierr
    INTEGER :: n_local_rows,local_rows_start,local_rows_end,local_cols_start,local_cols_end

    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)
    
    IF (mat%isTransposed) THEN
       n_local_rows = mat%NRows/n_procs
       local_rows_start = rank*n_local_rows+1
       local_rows_end   = (rank+1)*n_local_rows
       local_cols_start = 1
       local_cols_end   = mat%NCols

       min_k = mat%NRows
       max_k = 1
       DO icol=local_cols_start,local_cols_end
          DO irow = local_rows_start,local_rows_end
             this_min_k = MAX(irow-mat%LowerBandwidth,icol-mat%UpperBandwidth,1)
             IF (this_min_k.LT.min_k) min_k = this_min_k
             this_max_k = MIN(irow+mat%UpperBandwidth,icol+mat%LowerBandwidth,mat%NRows)
             IF (this_max_k.GT.max_k) max_k = this_max_k
          END DO
       END DO
    ELSE
    END IF
    
  END SUBROUTINE get_minmax_k_real

  !> Calculate the square of the given Banded Matrix. The square of a banded matrix is again
  !! a banded matrix with double bandwidths.
  SUBROUTINE sbm_square_matrix_real(mat,res)
    TYPE(StoreBandedMatrixObjectReal), INTENT(IN) :: mat
    TYPE(StoreBandedMatrixObjectReal), intent(INOUT) :: res

    !Local variables
    INTEGER :: rank, n_procs, ierr, min_k, max_k, min_proc, max_proc, psource,iproc, pdest, tag
    TYPE(IntegerList) :: recv_list, send_list
    type(IntegerNodeData) :: nodedata
    !INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: send_request(MAX_SEND_REQUESTS), recv_request, irequest,n_requests,&
         &calculate_rank
    REAL, DIMENSION(mat%NumberOfStoredRows,mat%ColsPerBlock),TARGET :: rbuf1,rbuf2
    REAL, DIMENSION(:,:),POINTER :: ptr_calculate, ptr_recv
    LOGICAL :: wait_for_recv

    CALL exit_if_factored(mat)

    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)
    CALL mpi_comm_size(mat%PG%communicator, n_procs, ierr)

    IF (mat%isTransposed) THEN
       call initialize(recv_list)
       call initialize(send_list)
       tag = 1001
       CALL get_minmax_k_real(mat,rank,min_k,max_k)
       !PRINT*,rank,": min_k = ", min_k,", max_k = ", max_k
       ! min_k and max_k contains the minimal used row or column on the local processor
       ! now we calculate the ranks of the processes from which we need data
       min_proc = (min_k-1)/mat%ColsPerBlock
       max_proc = (max_k-1)/mat%ColsPerBlock
       !PRINT*,"min_proc = ",min_proc,", max_proc = ", max_proc

       DO psource=min_proc,max_proc
          IF (psource.NE.rank) THEN
             nodedata%value = psource
             call push(recv_list,nodedata)
          END IF
       END DO
       !min/max_proc are the ranks of the processors from where we need data. 
       ! Further we need to know, to which ranks we have to send data
       DO iproc=0,n_procs-1
          CALL get_minmax_k_real(mat,iproc,min_k, max_k)
          min_proc = (min_k-1)/mat%ColsPerBlock
          max_proc = (max_k-1)/mat%ColsPerBlock
          !PRINT*,"for iproc = ", iproc,", min_k = ", min_k, ", max_k = ", max_k,", min_proc = ", min_proc, ", max_proc = ", max_proc
          IF ((rank.NE.iproc) .AND. (rank.LE.max_proc) .AND. (rank.GE.min_proc)) THEN
             nodedata%value=iproc
             CALL push(send_list,nodedata)
          END IF
       END DO

       ! We have now the full receicve list and send list 
       ! Now do all non blocking sends
       !IF (.NOT.isEmpty(send_list)) THEN
       !   WRITE(*,"(A)",advance="no") "send_list = "
       !   CALL output(send_list)
       !   WRITE(*,"(A)") ""
       !ELSE
       !   PRINT*,"send_list is empty"
       !END IF
       irequest = 1
       DO WHILE (.NOT.isEmpty(send_list))
          CALL pop(send_list,nodedata)
          pdest = nodedata%value
          PRINT*,"isend to ", pdest
          IF (irequest.GT.MAX_SEND_REQUESTS) THEN
             PRINT*,"You have to increase MAX_SEND_REQUESTS (actual set to ",MAX_SEND_REQUESTS,")!"
             stop
          END IF
          CALL mpi_isend(mat%localmatrix,mat%NumberOfStoredRows*mat%ColsPerBlock,MPI_REAL_TYPE,&
               &pdest,tag,mat%PG%communicator,send_request(irequest),ierr)
          irequest = irequest + 1
       END DO
       n_requests = irequest-1

       ! receiving the data is overlapped with calculation
       ptr_calculate => mat%LocalMatrix
       ptr_recv => rbuf1
       !IF (.NOT.isEmpty(recv_list)) THEN
       !   WRITE(*,"(A)",advance="no") "recv_list = "
       !   CALL output(recv_list)
       !   WRITE(*,"(A)") ""
       !ELSE
       !   PRINT*,"recv_list is empty"
       !END IF

       !PRINT*,"length(recv_list) = ", length(recv_list),isEmpty(recv_list)
       calculate_rank = rank ! start with local multiplication
       DO 
          IF (.NOT.isEmpty(recv_list)) THEN
             CALL pop(recv_list,nodedata)
             psource = nodedata%value
             PRINT*,"irecv from ",psource
             CALL mpi_irecv(ptr_recv,mat%NumberOfStoredRows*mat%ColsPerBlock,MPI_REAL_TYPE,&
                  &psource,tag,mat%PG%communicator,recv_request,ierr)
             wait_for_recv = .TRUE.
          ELSE
             wait_for_recv = .FALSE.
          END IF
          ! do the calculation LocalMatrix dot ptr_calculate
          PRINT*,"========================== local calculation ==============="
          CALL local_block_multiplication_real(mat,ptr_calculate,calculate_rank,res)
          IF (wait_for_recv) THEN
             CALL mpi_wait(recv_request,MPI_STATUS_IGNORE,ierr)
             PRINT*,"waiting for receive"
          
             IF (ASSOCIATED(ptr_recv,rbuf1)) THEN
                ptr_calculate => rbuf1
                ptr_recv => rbuf2
             ELSE
                ptr_calculate => rbuf2
                ptr_recv => rbuf1
             END IF
             calculate_rank = psource
          ELSE
             EXIT
          END IF
       END DO

       ! wait for finishing all non blocking sends
       DO irequest=1,n_requests
          !PRINT*,"Now waiting for send to complete for irequest = ",irequest," of ",n_requests
          CALL mpi_wait(send_request(irequest),MPI_STATUS_IGNORE,ierr)
       END DO

    ELSE
    END IF

  CONTAINS

    SUBROUTINE local_block_multiplication_real(mat,ptr_calculate,calculate_rank,res)
      TYPE(StoreBandedMatrixObjectReal), INTENT(IN) :: mat
      REAL, DIMENSION(:,:), INTENT(IN) :: ptr_calculate
      INTEGER, INTENT(IN) :: calculate_rank
      TYPE(StoreBandedMatrixObjectReal), INTENT(INOUT) :: res
      
    END SUBROUTINE local_block_multiplication_real
  END SUBROUTINE sbm_square_matrix_real


  !*****************************************************************
  !* The global matrix, which is distributed over all processes
  !* is gathered and stored in a local array. This routine is
  !* mainly used for debugging and backward compatibility. In the
  !* final production code, it should not be called as it does not
  !* scale in memory.
  !*****************************************************************
  SUBROUTINE sbm_get_global_matrix_locally_real(mat,localfullmat)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    REAL, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    ! local variables
    INTEGER :: ierr,iCol,iRow, global_col, global_row, rank

    CALL exit_if_factored(mat)

    ! otherwise, we have to take into account the cyclic block matrix
    ! structure, which complicates the operation.
    ! First we allocate a matrix to take the transposed local matrix
    !ALLOCATE(localfullmat(mat%NRows,mat%NCols))
    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)

    !PRINT*,"NumberOfStoredRows = ",mat%NumberOfStoredRows, mat%ColsPerBlock
    localfullmat = 0
    ! and we fill this local matrix
    IF (mat%isTransposed) THEN
       DO icol=1,mat%ColsPerBlock
          global_col = rank*mat%ColsPerBlock+icol
          DO irow = -mat%LowerBandwidth,mat%UpperBandwidth
             global_row = global_col+irow
             IF (global_row.GE.1 .AND. global_row.LE.mat%NRows) THEN
                localfullmat(global_row,global_col) = mat%localMatrix(irow,icol)
             END IF
          END DO
       END DO
       
       !PRINT*, "Communicator in MPI_Allgather is ", mat%PG%Communicator
       CALL mpi_allgather(MPI_IN_PLACE,mat%NRows*mat%ColsPerBlock,MPI_REAL_TYPE,&
            &localfullmat,mat%NRows*mat%ColsPerBlock,MPI_REAL_TYPE,&
            &mat%PG%Communicator,ierr)
       localfullmat = transpose(localfullmat)
    ELSE
       DO icol=1,mat%ColsPerBlock
          global_col = rank*mat%ColsPerBlock+icol
          DO irow = -mat%UpperBandwidth,mat%LowerBandwidth
             global_row = global_col+irow
             IF (global_row.GE.1 .AND. global_row.LE.mat%NRows) THEN
                localfullmat(global_row,global_col) = mat%localMatrix(irow,icol)
             END IF
          END DO
       END DO
       
       !PRINT*, "Communicator in MPI_Allgather is ", mat%PG%Communicator
       CALL mpi_allgather(MPI_IN_PLACE,mat%NRows*mat%ColsPerBlock,MPI_REAL_TYPE,&
            &localfullmat,mat%NRows*mat%ColsPerBlock,MPI_REAL_TYPE,&
            &mat%PG%Communicator,ierr)
    END IF

  END SUBROUTINE sbm_get_global_matrix_locally_real

  FUNCTION sbm_get_local_abs_square_sum_real(mat) RESULT(value)
    type(StoreBandedMatrixObjectReal) :: mat
    REAL :: value

    INTEGER :: icol, global_col, irow, global_row, rank, ierr

    CALL exit_if_factored(mat)

    CALL mpi_comm_rank(mat%PG%communicator,rank,ierr)

    value = 0.0D0
    IF (mat%isTransposed) THEN
       DO icol=1,mat%ColsPerBlock
          !global_col = rank*mat%ColsPerBlock+icol
          global_row = rank*mat%ColsPerBlock+icol

          DO irow = -mat%LowerBandwidth,mat%UpperBandwidth
             !global_row = global_col+irow
             global_col = global_row+irow
             !IF (global_row.GE.1 .AND. global_row.LE.mat%NRows) THEN
             IF (global_col.GE.1 .AND. global_col.LE.mat%NCols) THEN
                value = value +  REAL(mat%localMatrix(irow,icol)*mat%localMatrix(irow,icol))
             END IF
          END DO
       END DO
    ELSE
       DO icol=1,mat%ColsPerBlock
          global_col = rank*mat%ColsPerBlock+icol

          DO irow = -mat%UpperBandwidth,mat%LowerBandwidth
             global_row = global_col+irow
             IF (global_row.GE.1 .AND. global_row.LE.mat%NRows) THEN
                value = value +  REAL(mat%localMatrix(irow,icol)*mat%localMatrix(irow,icol))
             END IF
          END DO
       END DO
    END IF
  END FUNCTION sbm_get_local_abs_square_sum_real

  SUBROUTINE sbm_show_on_screen_real(mat)
    type(StoreBandedMatrixObjectReal) :: mat

    ! Local variables
    INTEGER :: iProc,row_start,row_end,iCol, iRow

    CALL exit_if_factored(mat)

    ! do some debugging output
    !call mpi_barrier(mat%PG%Communicator,ierr)
    CALL blacs_barrier(mat%PG%Context,'A')
    DO iProc=0,mat%PG%NProcs-1
       !call mpi_barrier(mat%PG%Communicator,ierr)
       CALL blacs_barrier(mat%PG%Context,'A')
       IF (iProc.EQ.mat%PG%Rank) THEN
          PRINT*,"---- Start output Rank = ",mat%PG%Rank," ----"
          IF (mat%isTransposed) THEN
             row_start = -mat%LowerBandwidth
             row_end   = mat%UpperBandwidth
          ELSE
             row_start = -mat%UpperBandwidth
             row_end   = mat%LowerBandwidth
          END IF
          DO iRow=row_start,row_end
             DO iCol = 1,mat%ColsPerBlock
                WRITE(*,"(2ES9.2,1X)",advance='NO') mat%localMatrix(irow,icol)
             END DO
             WRITE(*,"(A)") ""
          END DO
          PRINT*,"---- End   output Rank = ",mat%PG%Rank," ----"
       END IF
    END DO
    !call mpi_barrier(mat%PG%Communicator,ierr)
    CALL blacs_barrier(mat%PG%Context,'A')
  END SUBROUTINE sbm_show_on_screen_real

  SUBROUTINE sbm_show_in_file_real(mat,filename)
    type(StoreBandedMatrixObjectReal) :: mat
    CHARACTER(len=*) :: filename

    ! Local variables
    INTEGER :: iCol, iRow, row_start, row_end
    character(len=FILENAME_MAX) :: full_filename

    CALL exit_if_factored(mat)

    ! open a file on each process
    WRITE(full_filename,"(A,I3.3,A)") trim(filename),mat%PG%Rank,".dat"
    OPEN(77,file=TRIM(full_filename))

    ! do some debugging output
    CALL blacs_barrier(mat%PG%Context,'A')
    IF (mat%isTransposed) THEN
       row_start = -mat%LowerBandwidth
       row_end   = mat%UpperBandwidth
    ELSE
       row_start = -mat%UpperBandwidth
       row_end   = mat%LowerBandwidth
    END IF
    DO iRow=row_start,row_end
       DO iCol = 1,mat%ColsPerBlock
          WRITE(77,"(2ES20.10,1X)",advance='NO') mat%localMatrix(irow,icol)
       END DO
       WRITE(77,"(A)") ""
    END DO
    CLOSE(77)
    CALL blacs_barrier(mat%PG%Context,'A')
  END SUBROUTINE sbm_show_in_file_real

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE sbm_add_to_matrix_real(this,mat)
    TYPE(StoreBandedMatrixObjectReal), INTENT(INOUT) :: this
    TYPE(StoreBandedMatrixObjectReal), INTENT(IN) :: mat

    INTEGER :: rank, ierr, irow, icol, max_UpperBandwidth, max_LowerBandwidth

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)

    ! nrows and ncols have to be the same for both matrices,
    ! this is NOT checked!
    CALL mpi_comm_rank(this%PG%communicator, rank,ierr)
#if 0
    IF (this%UpperBandwidth.EQ.mat%UpperBandwidth) THEN
       ! 1. BWUs are equal
       DO icol=1,this%ColsPerBlock
          this%localMatrix(-this%UpperBandwidth:0,icol) = &
               & this%localMatrix(1:this%UpperBandwidth+1,icol) &
               & + mat%localMatrix(1:mat%UpperBandwidth+1,icol)
       END DO
    ELSEIF (this%UpperBandwidth.gt.mat%UpperBandwidth) THEN
       BWU_diff = this%UpperBandwidth-mat%UpperBandwidth
       ! 2. BWU of target is larger 
       DO icol=1,this%ColsPerBlock
          this%localMatrix(BWU_diff+1:BWU_diff+mat%UpperBandwidth+1,icol) = &
               & this%localMatrix(BWU_diff+1:BWU_diff+mat%UpperBandwidth+1,icol) + &
               & mat%localMatrix(1:mat%UpperBandwidth+1,icol)
       END DO
    ELSE
       ! 3. BWU of target is smaller than BWU of source
       DO icol=rank*mat%ColsPerBlock+1,(rank+1)*mat%ColsPerBlock
          DO irow=icol-mat%UpperBandwidth,icol
             IF ((irow.ge.1).AND.(irow.le.mat%NRows)) THEN
                CALL add_value(this,irow,icol,mat_get_value(mat,irow,icol))
             END IF
          END DO
       END DO
       call commit_values(this)
    END IF
#endif
    ! Get the values of mat and add them to this
    IF (mat%isTransposed.AND.this%isTransposed) THEN
       max_UpperBandwidth = MAX(mat%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(mat%LowerBandwidth, this%LowerBandwidth)
       DO irow = rank*mat%ColsPerBlock+1,(rank+1)*mat%ColsPerBlock
          DO icol=irow-max_LowerBandwidth,irow+max_UpperBandwidth
             IF ((icol.GE.1).AND.(icol.LE.mat%NCols)) THEN
                CALL add_value(this,irow,icol,mat_get_value(mat,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSEIF (.NOT.mat%isTransposed .AND. .NOT.this%isTransposed) THEN
       max_UpperBandwidth = MAX(mat%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(mat%LowerBandwidth, this%LowerBandwidth)
       DO icol=rank*mat%ColsPerBlock+1,(rank+1)*mat%ColsPerBlock
          DO irow=icol-max_UpperBandwidth,icol+max_LowerBandwidth
             IF ((irow.GE.1).AND.(irow.LE.mat%NRows)) THEN
                CALL add_value(this,irow,icol,mat_get_value(mat,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSE
       PRINT*,"You can add only two transposed matrices or two not transposed matrices, but not mixed."
       STOP
    END IF
  END SUBROUTINE sbm_add_to_matrix_real

  SUBROUTINE sbm_add_matrix_real(this,mat,res)
    TYPE(StoreBandedMatrixObjectReal), INTENT(IN) :: this,mat
    TYPE(StoreBandedMatrixObjectReal), INTENT(INOUT) :: res

    res = this

    CALL add_matrix(res,mat)
  END SUBROUTINE sbm_add_matrix_real

  SUBROUTINE sbm_subtract_matrix_real(this,mat,res)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: this,mat
    TYPE(StoreBandedMatrixObjectReal),intent(INOUT) :: res
    ! calculate res = this - mat

    res = this

    CALL subtract_matrix(res,mat)

  END SUBROUTINE sbm_subtract_matrix_real

  SUBROUTINE sbm_subtract_from_matrix_real(this,mat)
    TYPE(StoreBandedMatrixObjectReal),INTENT(INOUT) :: this
    TYPE(StoreBandedMatrixObjectReal),intent(IN) :: mat
    ! calculate this = this - mat

    INTEGER :: rank, ierr, irow, icol, max_LowerBandwidth, max_UpperBandwidth

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)

    ! nrows and ncols have to be the same for both matrices,
    ! this is NOT checked!
    CALL mpi_comm_rank(this%PG%communicator, rank,ierr)

    ! Get the values of mat and add them to this
    IF (mat%isTransposed.AND.this%isTransposed) THEN
       max_UpperBandwidth = MAX(mat%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(mat%LowerBandwidth, this%LowerBandwidth)
       DO irow = rank*mat%ColsPerBlock+1,(rank+1)*mat%ColsPerBlock
          DO icol=irow-max_LowerBandwidth,irow+max_UpperBandwidth
             IF ((icol.GE.1).AND.(icol.LE.mat%NCols)) THEN
                CALL add_value(this,irow,icol,-mat_get_value(mat,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSEIF (.NOT.mat%isTransposed .AND. .NOT.this%isTransposed) THEN
       max_UpperBandwidth = MAX(mat%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(mat%LowerBandwidth, this%LowerBandwidth)
       DO icol=rank*mat%ColsPerBlock+1,(rank+1)*mat%ColsPerBlock
          DO irow=icol-max_UpperBandwidth,icol+max_LowerBandwidth
             IF ((irow.GE.1).AND.(irow.LE.mat%NRows)) THEN
                CALL add_value(this,irow,icol,-mat_get_value(mat,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSE
       PRINT*,"You can subtract only two transposed matrices or two not transposed matrices, but not mixed."
       STOP
    END IF

  END SUBROUTINE sbm_subtract_from_matrix_real

  SUBROUTINE sbm_multiply_matrix_with_real_real(this, scalar, res)
    TYPE(StoreBandedMatrixObjectReal), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(StoreBandedMatrixObjectReal),intent(INOUT) :: res

    INTEGER :: irow,icol,rank,ierr, max_LowerBandwidth, max_UpperBandwidth

    CALL exit_if_factored(this)
    CALL exit_if_factored(res)

    CALL mpi_comm_rank(this%PG%communicator,rank,ierr)

    ! Get the values of mat and add them to this
    IF (res%isTransposed.AND.this%isTransposed) THEN
       max_UpperBandwidth = MAX(res%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(res%LowerBandwidth, this%LowerBandwidth)
       DO irow = rank*res%ColsPerBlock+1,(rank+1)*res%ColsPerBlock
          DO icol=irow-max_LowerBandwidth,irow+max_UpperBandwidth
             IF ((icol.GE.1).AND.(icol.LE.res%NCols)) THEN
                CALL set_value(res,irow,icol,scalar*mat_get_value(this,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)

    ELSEIF (.NOT.res%isTransposed .AND. .NOT.this%isTransposed) THEN
       max_UpperBandwidth = MAX(res%UpperBandwidth,this%UpperBandwidth)
       max_LowerBandwidth = MAX(res%LowerBandwidth, this%LowerBandwidth)
       DO icol=rank*res%ColsPerBlock+1,(rank+1)*res%ColsPerBlock
          DO irow=icol-max_UpperBandwidth,icol+max_LowerBandwidth
             IF ((irow.GE.1).AND.(irow.LE.res%NRows)) THEN
                CALL set_value(res,irow,icol,scalar*mat_get_value(this,irow,icol))
             END IF
          END DO
       END DO
       CALL commit_values(this)
    ELSE
       PRINT*,"You can subtract only two transposed matrices or two not transposed matrices, but not mixed."
       STOP
    END IF

  END SUBROUTINE sbm_multiply_matrix_with_real_real

  SUBROUTINE sbm_scale_matrix_by_real_real(this, scalar)
    TYPE(StoreBandedMatrixObjectReal), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL exit_if_factored(this)

    this%localMatrix = this%localMatrix*scalar
  END SUBROUTINE sbm_scale_matrix_by_real_real

  SUBROUTINE sbm_assign_matrix_real(lmat,rmat)
    TYPE(StoreBandedMatrixObjectReal), intent(INOUT) :: lmat
    TYPE(StoreBandedMatrixObjectReal), intent(IN) :: rmat

    INTEGER :: irow,icol,rank,ierr

    CALL exit_if_factored(lmat)
    CALL exit_if_factored(rmat)

    IF (lmat%PG.EQ.rmat%PG) THEN
       CALL mpi_comm_rank(lmat%PG%communicator,rank,ierr)
       IF (lmat%isTransposed.AND.rmat%isTransposed) THEN
          ! both matrices are stored in row-major order
          DO irow=rank*rmat%ColsPerBlock+1,(rank+1)*rmat%ColsPerBlock
             DO icol=irow-rmat%LowerBandwidth,irow+rmat%UpperBandwidth
                IF ((icol.GE.1).AND.(icol.LE.rmat%NCols)) THEN
                   CALL set_value(lmat,irow,icol,mat_get_value(rmat,irow,icol))
                END IF
             END DO
          END DO
          CALL commit_values(lmat)
       ELSEIF (.not.lmat%isTransposed.and..not.rmat%isTransposed) then
          DO icol=rank*rmat%ColsPerBlock+1,(rank+1)*rmat%ColsPerBlock
             DO irow=icol-rmat%UpperBandwidth,icol+rmat%LowerBandwidth
                IF ((irow.GE.1).AND.(irow.LE.rmat%NRows)) THEN
                   CALL set_value(lmat,irow,icol,mat_get_value(rmat,irow,icol))
                END IF
             END DO
          END DO
          CALL commit_values(lmat)
       ELSE
          PRINT*,"Assignment is only possible if both matrices have the same storage scheme (isTransposed)."
          STOP
       END IF
    ELSE
       PRINT*,"Assignment of matrices is only possible within the same context."
       STOP
    END IF
  END SUBROUTINE sbm_assign_matrix_real
  
  FUNCTION sbm_get_number_of_bands_real(mat)
    type(StoreBandedMatrixObjectReal) :: mat
    INTEGER :: sbm_get_number_of_bands_real

    sbm_get_number_of_bands_real=mat%LowerBandwidth+mat%UpperBandwidth+1
  END FUNCTION sbm_get_number_of_bands_real

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE sbm_matrix_sum_real(mat,communicator)
    TYPE(StoreBandedMatrixObjectReal), INTENT(INOUT) :: mat
    integer :: communicator

    INTEGER :: max_lower_bw, max_upper_bw, ierr, icol,irow
    REAL, DIMENSION(:,:), pointer :: newMatrix

    CALL exit_if_factored(mat)

    ! For a banded matrix, we have to 
    ! 1) get the maximal bandwidths over the communicator
    CALL mpi_allreduce(mat%LowerBandwidth,max_lower_bw,1,MPI_INTEGER,&
         &MPI_MAX,communicator,ierr)
    CALL mpi_allreduce(mat%UpperBandwidth,max_upper_bw,1,MPI_INTEGER,&
         &MPI_MAX,communicator,ierr)
    
    ! 2) extend all matrices to this maximum bandwidth
    IF ((max_lower_bw.GT.mat%LowerBandwidth).OR.(max_upper_bw.GT.mat%UpperBandwidth)) THEN
       ! Attention, newMatrix will not be deallocated, as the mat%localMatrix will point to it!
       IF (mat%isTransposed) THEN
          ALLOCATE(newMatrix(-max_lower_bw:max_upper_bw,mat%ColsPerBlock))
       ELSE
          ALLOCATE(newMatrix(-max_upper_bw:max_lower_bw,mat%ColsPerBlock))
       END IF
       newMatrix = 0

       ! 4. copy old matrix into new one
       IF (mat%isTransposed) THEN
          DO icol=1,mat%ColsPerBlock
             DO irow=-mat%LowerBandwidth,mat%UpperBandwidth
                newMatrix(irow,icol) = mat%localMatrix(irow,icol)
             END DO
          END DO
       ELSE
          DO icol=1,mat%ColsPerBlock
             DO irow=-mat%UpperBandwidth,mat%LowerBandwidth
                newMatrix(irow,icol) = mat%localMatrix(irow,icol)
             END DO
          END DO
       END IF
       DEALLOCATE(mat%localMatrix)
       mat%localMatrix => newMatrix
       mat%UpperBandwidth = max_upper_bw
       mat%LowerBandwidth = max_lower_bw
       mat%NumberOfStoredRows = mat%UpperBandwidth+mat%LowerBandwidth+1
    END IF

    ! all matrices have now the same shape, so we can
    ! 3) allreduce over all matrices

    CALL mpi_allreduce(MPI_IN_PLACE,mat%localMatrix,mat%NumberOfStoredRows*mat%ColsPerBlock, &
         & MPI_REAL_TYPE, MPI_SUM, communicator, ierr)
    
  END SUBROUTINE sbm_matrix_sum_real

  SUBROUTINE sbm_row_axpy_real(this,this_row,mat,scalar)
    TYPE(StoreBandedMatrixObjectReal),intent(INOUT) :: this
    INTEGER, INTENT(IN) :: this_row
    TYPE(StoreBandedMatrixObjectReal),intent(IN) :: mat
    REAL, intent(IN) :: scalar

    ! Local variables
    INTEGER :: rank, rank_of_global_row,irow, global_col,ierr, min_row,max_row, stored_column

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)

    CALL mpi_comm_rank(this%PG%communicator,rank,ierr)
    
    IF (this%isTransposed) THEN
       IF (.NOT.mat%isTransposed) THEN
          PRINT*,"sbm_row_axpy: both matrices must have the same storage mode"
          STOP
       END IF
       ! both matrices are transposed
       
       ! calculate the local index out of the global row index
       rank_of_global_row = (this_row-1)/this%ColsPerBlock
       ! do only the calculation, if the global row is stored locally on the actual processor
       IF (rank.EQ.rank_of_global_row) THEN
          min_row = MIN(this%LowerBandwidth,mat%LowerBandwidth)
          max_row = MIN(this%UpperBandwidth,mat%UpperBandwidth)
          stored_column = MOD(this_row-1,this%ColsPerBlock)+1
          this%localMatrix(-min_row:max_row,stored_column) = &
               & this%localMatrix(-min_row:max_row,stored_column)&
               & + mat%localMatrix(-min_row:max_row,stored_column)*scalar
          IF (mat%LowerBandwidth.GT.this%LowerBandwidth) THEN
             DO irow=-mat%LowerBandwidth,-min_row-1
                ! convert to global column index
                global_col = this_row+irow
                CALL set_value(this,this_row,global_col,mat%localMatrix(irow,stored_column)*scalar)
             END DO
          END IF
          IF (mat%UpperBandwidth.GT.this%UpperBandwidth) THEN
             DO irow=max_row+1,mat%UpperBandwidth
                ! convert to global column index
                global_col = this_row+irow
                CALL set_value(this,this_row,global_col,mat%localMatrix(irow,stored_column)*scalar)
             END DO
          END IF
       END IF
       CALL commit_values(this)
    ELSE
       IF (mat%isTransposed) THEN
          PRINT*,"sbm_row_axpy: both matrices must have the same storage mode"
          STOP
       END IF
       ! both matrices are not transposed
    END IF

  END SUBROUTINE sbm_row_axpy_real

  subroutine sbm_transpose_and_conjugate_real(mat,dagger_mat)
    type(StoreBandedMatrixObjectReal) :: mat
    type(StoreBandedMatrixObjectReal) :: dagger_mat

    ! To transpose, we set the isTransposed flag of dagger_mat
    ! to the inverse of the isTransposed flag of mat
    dagger_mat%isTransposed = .not.mat%isTransposed
    ! Also exchange NRows and NCols, and the Bandwidths
    dagger_mat%NRows = mat%NCols
    dagger_mat%NCols = mat%NRows
    dagger_mat%LowerBandwidth = mat%UpperBandwidth
    dagger_mat%UpperBandwidth = mat%LowerBandwidth
    dagger_mat%NumberOfStoredRows = mat%NumberOfStoredRows
    dagger_mat%ColsPerBlock = mat%ColsPerBlock

    ! now copy the localmatrices and conjugate
    call allocate(dagger_mat)
    !if (associated(dagger_mat%localmatrix)) deallocate(dagger_mat%localmatrix)
    
    !allocate(dagger_mat%localmatrix(-dagger_mat%UpperBandwidth:dagger_mat%LowerBandwidth,&
    !     &1:dagger_mat%ColsPerBlock))
    dagger_mat%localmatrix = mat%localmatrix
  end subroutine sbm_transpose_and_conjugate_real

  SUBROUTINE sbm_LU_factor_real(mat)
    type(StoreBandedMatrixObjectReal) :: mat
    
    ! Local variables
    REAL :: work(1)
    INTEGER :: i, n, bwu, bwl, bw_max, laf, lwork, info, desca(7)
    INTEGER :: my_rank, n_procs, mpierr
    !INTEGER :: mpi_status(mpi_status_size)
    INTEGER :: np_send, np_recv
    REAL, ALLOCATABLE :: sbuf(:), rbuf(:)

    CALL exit_if_factored(mat)

    ! if PG_cols is not yet initialized we have to do it here

    if(.NOT. PG_cols%isInitialized) CALL initialize(PG_cols,MAT_BLOCKS_OF_COLS,mat%PG%Communicator)

    call mpi_comm_rank(mat%PG%Communicator,my_rank,mpierr)
    call mpi_comm_size(mat%PG%Communicator,n_procs,mpierr)

    IF (mat%Nrows /= mat%Ncols) THEN
      PRINT *,'Matrix must be square for sbm_factor_matrix'
      CALL mpi_abort(MPI_COMM_WORLD, 1, mpierr)
    ENDIF

    bwu = mat%UpperBandwidth
    bwl = mat%LowerBandwidth

    bw_max = max(bwu,bwl)

    ! pcgbtrf/pcgbtrs work only if bwu+bwl+1 <= mat%ColsPerBlock.
    ! NB: The transposition below works only if bw_max <= mat%ColsPerBlock,
    ! but this is a weaker restriction

    IF(bwu+bwl+1 > mat%ColsPerBlock) THEN
      print *,"ERROR sbm_factor: bwu+bwl+1, ColsPerBlock ",bwu+bwl+1,mat%ColsPerBlock
      CALL mpi_abort(MPI_COMM_WORLD, 1, mpierr)
    ENDIF

    ALLOCATE(mat%factoredBand(-2*bwu-bwl:bwl,mat%ColsPerBlock))
    mat%factoredBand(:,:) = 0

    ! fill in mat%factoredBand

    IF(mat%isTransposed) THEN

      ! Please note: pcgbtrf/pcgbtrs are NOT able to deal with transposed matrices
      ! although the first argument of pcgbtrs pretends that.
      ! In reality the first argument of pcgbtrs is completely ignored.
      ! Thus we have to transpose the matrix here ...

      allocate(sbuf(bw_max*(bw_max+1)/2))
      allocate(rbuf(bw_max*(bw_max+1)/2))

      ! Super diagonals

      ! Send pieces at end of localMatrix to next processor

      ! Please note: It is not necessary to do a cyclic exchange
      ! but it is convenient to do so when using sendrecv
      np_send = MOD(my_rank+1,n_procs)
      np_recv = MOD(my_rank+n_procs-1,n_procs)

      n = 0
      DO i = 1, bwu
        sbuf(n+1:n+i) = mat%localMatrix(i,mat%ColsPerBlock-i+1:mat%ColsPerBlock)
        n = n+i
      ENDDO
      if(my_rank==n_procs-1) sbuf(:) = 0 ! Safety only

      call mpi_sendrecv(sbuf, n, MPI_REAL_TYPE, np_send, 111, &
                        rbuf, n, MPI_REAL_TYPE, np_recv, 111, &
                        mat%PG%Communicator, MPI_STATUS_IGNORE, mpierr)

      n = 0
      DO i = 1, bwu
        mat%factoredBand(-i,1+i:mat%ColsPerBlock) = mat%localMatrix(i,1:mat%ColsPerBlock-i)
        mat%factoredBand(-i,1:i) = rbuf(n+1:n+i)
        n = n+i
      ENDDO

      ! Diagonal

      mat%factoredBand(0,1:mat%ColsPerBlock) = mat%localMatrix(0,1:mat%ColsPerBlock)

      ! Sub diagonals

      ! Send pieces at start of localMatrix to previous processor

      np_recv = MOD(my_rank+1,n_procs)
      np_send = MOD(my_rank+n_procs-1,n_procs)

      n = 0
      DO i = 1, bwl
        sbuf(n+1:n+i) = mat%localMatrix(-i,1:i)
        n = n+i
      ENDDO

      if(my_rank==0) sbuf(:) = 0 ! Safety only

      call mpi_sendrecv(sbuf,n,MPI_REAL_TYPE,np_send,111, &
                        rbuf,n,MPI_REAL_TYPE,np_recv,111, &
                        mat%PG%Communicator, MPI_STATUS_IGNORE, mpierr)

      n = 0
      DO i = 1, bwl
        mat%factoredBand(i,1:mat%ColsPerBlock-i)  = mat%localMatrix(-i,1+i:mat%ColsPerBlock)
        mat%factoredBand(i,mat%ColsPerBlock-i+1:mat%ColsPerBlock) = rbuf(n+1:n+i)
        n = n+i
      ENDDO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)

    ELSE

      mat%factoredBand(-bwu:bwl,:) = mat%localMatrix(-bwu:bwl,:)

    ENDIF

    ! The doc in the header of PDGBTRF isn't clear about IPIV,
    ! I hope the following is enough ...
    ALLOCATE(mat%ipiv(2*mat%ColsPerBlock+bwu+bwl+1))

    ! allocate additional storage
    laf = (mat%ColsPerBlock+bwu)*(bwl+bwu)+6*(bwl+bwu)*(bwl+2*bwu)
    ALLOCATE(mat%auxData(laf))

    ! Size of work is not stated in header of pcgbtrf, obviously it is 1
    lwork = SIZE(work)

    desca = (/ 501, PG_cols%Context, mat%Nrows, mat%ColsPerBlock, 0, 2*bwl+2*bwu+1, 0 /)

    CALL pcgbtrf(mat%Nrows, bwl, bwu, mat%factoredBand, 1, desca, &
                 mat%ipiv, mat%auxData, laf, work, lwork, info)

    IF (info.NE.0) THEN
       PRINT*,"pcgbtrf info = ",info
       stop
    END IF

    mat%isFactored = .true.

    ! Now we can deallocate the original matrix data

    deallocate(mat%localMatrix)
    nullify(mat%localMatrix)

  END SUBROUTINE sbm_LU_factor_real

  LOGICAL FUNCTION sbm_LU_factor_ok_real(mat)
    type(StoreBandedMatrixObjectReal) :: mat

    ! Checks if the bandwidth requirement for pcgbtrf/pcgbtrs is ok.

    IF(mat%isFactored) THEN
      ! In this case it is of course not ok to factor ...
      sbm_LU_factor_ok_real = .FALSE.
      RETURN
    ENDIF

    IF(mat%ColsPerBlock >= mat%LowerBandwidth+mat%UpperBandwidth+1) THEN
      sbm_LU_factor_ok_real = .TRUE.
    ELSE
      sbm_LU_factor_ok_real = .FALSE.
    ENDIF
  END FUNCTION sbm_LU_factor_ok_real

  SUBROUTINE sbm_LU_solve_real(mat, vec, res)
    TYPE(StoreBandedMatrixObjectReal),INTENT(IN) :: mat
    TYPE(StoreVectorObjectReal),INTENT(IN) :: vec
    TYPE(StoreVectorObjectReal),INTENT(INOUT) :: res

    INTEGER :: bwu, bwl, lwork, info, desca(7), descb(7)
    REAL, ALLOCATABLE :: work(:)

    ! checks if vec and res fit to mat

    if(vec%RowsPerBlock /= mat%ColsPerBlock) STOP 'sbm_LU_solve: mat and vec not conforming'
    if(res%RowsPerBlock /= mat%ColsPerBlock) STOP 'sbm_LU_solve: mat and res not conforming'

    if(.not. mat%isFactored) STOP 'solve_LU called with un-factored matrix'

    bwu = mat%UpperBandwidth
    bwl = mat%LowerBandwidth

    ! set up despriptors

    desca = (/ 501, PG_cols%Context, mat%Nrows, mat%ColsPerBlock, 0, 2*bwl+2*bwu+1, 0 /)
    descb = (/ 502, PG_cols%Context, mat%Nrows, mat%ColsPerBlock, 0, mat%ColsPerBlock, 0 /)

    ! Size of work is not stated in header of pcgbtrs, however in the source code:

    lwork = mat%ColsPerBlock + 2*bwl + 4*bwu
    ALLOCATE(work(lwork))

    res%localVector = vec%localVector

    CALL pcgbtrs('N', mat%Nrows, bwl, bwu, 1, mat%factoredBand, 1, desca, &
                 mat%ipiv, res%localVector, 1, descb, mat%auxData, size(mat%auxData), &
                 work, lwork, info)

    DEALLOCATE(work)

    IF (info.NE.0) THEN
       PRINT*,"pcgbtrs info = ",info
       stop
    END IF

  END SUBROUTINE sbm_LU_solve_real

  SUBROUTINE sbm_exit_if_factored_real(mat)
    type(StoreBandedMatrixObjectReal) :: mat

    INTEGER :: ierr

    ! This routines is called by all routines which use mat%localMatrix.
    ! It does exit the program if the matrix is LU-factored.

    IF(mat%isFactored) THEN
      PRINT *,'ERROR: Attempt use original matrix data of a LU-factored Matrix'
      CALL mpi_abort(MPI_COMM_WORLD, 1, ierr)
    ENDIF

  END SUBROUTINE sbm_exit_if_factored_real

  subroutine sbm_print_storage_details_real(mat)
    type(StoreBandedMatrixObjectReal) :: mat

    print*,"Storage details: "
    print*,"ColsPerBlock,NumberOfStoredRows = ",mat%ColsPerBlock,mat%NumberOfStoredRows
  end subroutine sbm_print_storage_details_real

#ifdef ALL_ROUTINES
  ! ===================================================
  ! == Define some mathematical routines on matrices ==
  ! ===================================================
  SUBROUTINE sbm_invert_matrix_real(mat)
    type(StoreBandedMatrixObjectReal) :: mat
    
    ! Local variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv, iwork
    REAL, DIMENSION(:), allocatable :: work
    INTEGER :: info,lwork,liwork,ub_locr

    ! calculate upper bound of ipiv storage needed
    ! LOCr(M_A)+MB_A (ScaLapack User's guide)
    ub_locr = ( (mat%NRows+mat%RowsPerBlock-1)/mat%RowsPerBlock &
         & + mat%PG%NProcs-1 ) / mat%PG%NProcs * mat%RowsPerBlock

    !PRINT*,"ub_locr = ",ub_locr
    ALLOCATE(ipiv(ub_locr+mat%RowsPerBlock))
    !WRITE(*,"(16ES10.2)") mat%localMatrix
    ! First do the LU decomposition in place
    CALL pcgetrf(mat%NRows,mat%NCols,&
         &mat%localMatrix, 1, 1,&
         &mat%Desc,ipiv,info)

    IF (info.NE.0) THEN
       PRINT*,"pcgetrf produced an error in StoreBandedMatrix/sbm_invert_matrix"
       stop
    END IF
    !PRINT*,"after LU decomposition"
    !WRITE(*,"(16ES10.2)") mat%localMatrix
    

    ! Second calculate the inverse of the matrix
    ! determine the work array size
    ALLOCATE(work(1),iwork(1))
    CALL pcgetri(mat%Nrows,mat%localMatrix,1,1,mat%Desc,ipiv,&
         &work,-1,iwork,-1,info)
    lwork = INT(work(1))
    liwork = iwork(1)
    !PRINT*,"lwork = ",lwork,", liwork = ",liwork
    DEALLOCATE(work,iwork)

    !PRINT*,"after first pcgetri call"
    !WRITE(*,"(16ES10.2)") mat%localMatrix

    
    ALLOCATE(work(lwork),iwork(liwork))
    CALL pcgetri(mat%Nrows,mat%localMatrix,1,1,mat%Desc,ipiv,&
         &work,lwork,iwork,liwork,info)

    IF (info.NE.0) THEN
       PRINT*,"pcgetri produced an error in StoreBandedMatrix/sbm_invert_matrix"
       stop
    END IF
    !PRINT*,"after second pcgetri call"
    !WRITE(*,"(16ES10.2)") mat%localMatrix

    DEALLOCATE(work,iwork)
    DEALLOCATE(ipiv)
  END SUBROUTINE sbm_invert_matrix_real
#endif

END MODULE StoreBandedMatrixModule
