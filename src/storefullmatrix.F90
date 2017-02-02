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
#define OLD_PROCESSGRID

MODULE StoreFullMatrixModule
  USE ProcessGridModule
  use StoreVectorModule
  USE mpi
  IMPLICIT NONE

  !> Type to hold the definition of a full matrix
  !!
  !! In this type, the real storage of a full matrix is defined. It holds
  !! all attributes which are needed to store and use the data of the matrix.
  TYPE,public :: StoreFullMatrixObject 
     TYPE(ProcessGrid) :: PG
     INTEGER :: NRows, NCols !< number of rows and cols globally
     INTEGER :: RowsPerBlock, ColsPerBlock !< block size of distributed matrix, MB, NB
     INTEGER :: NBlocks !< number of blocks per process
     INTEGER :: Blocksize !< size of a block, MB x NB
     INTEGER :: BlockrowsPerProcess !< Number of rows of blocks per process
     integer :: BlockcolsPerProcess !< Number of cols of blocks per process
     LOGICAL :: isVector
     LOGICAL :: isAttached !< switch to show, if the matrix is attached to a field outside of the matrix structure
     LOGICAL :: isInitialized=.FALSE.
     LOGICAL :: isFactored = .FALSE.        !< determines if matrix is LU factored
     INTEGER,DIMENSION(9) :: Desc           !< matrix descriptor for BLACS, ScaLapack etc. 
     ! the per processor local matrix part of the global matrix
     COMPLEX, DIMENSION(:), POINTER :: localMatrix   !< either points to an array outside or to an allocated array

     INTEGER :: bwl, bwu
     COMPLEX, DIMENSION(:,:), POINTER :: factoredBand
     INTEGER, DIMENSION(:), POINTER :: ipiv
     COMPLEX, DIMENSION(:), POINTER :: rhs
  END TYPE StoreFullMatrixObject

TYPE, public :: StoreFullMatrixObjectReal 
     TYPE(ProcessGrid) :: PG
     INTEGER :: NRows, NCols !< number of rows and cols globally
     INTEGER :: RowsPerBlock, ColsPerBlock !< block size of distributed matrix, MB, NB
     INTEGER :: NBlocks !< number of blocks per process
     INTEGER :: Blocksize !< size of a block, MB x NB
     INTEGER :: BlockrowsPerProcess !< Number of rows of blocks per process
     integer :: BlockcolsPerProcess !< Number of cols of blocks per process
     LOGICAL :: isVector
     LOGICAL :: isAttached !< switch to show, if the matrix is attached to a field outside of the matrix structure
     LOGICAL :: isInitialized=.FALSE.
     LOGICAL :: isFactored = .FALSE.        !< determines if matrix is LU factored
     INTEGER,DIMENSION(9) :: Desc           !< matrix descriptor for BLACS, ScaLapack etc. 
     ! the per processor local matrix part of the global matrix
     REAL, DIMENSION(:), POINTER :: localMatrix   !< either points to an array outside or to an allocated array

     INTEGER :: bwl, bwu
     REAL, DIMENSION(:,:), POINTER :: factoredBand
     INTEGER, DIMENSION(:), POINTER :: ipiv
     REAL, DIMENSION(:), POINTER :: rhs
  END TYPE StoreFullMatrixObjectReal

  INTERFACE initialize
     MODULE PROCEDURE mp_initialize_matrix, mp_initialize_matrix_real
  END INTERFACE

  INTERFACE isInitialized
     module procedure mp_isInitialized, mp_isInitialized_real
  END INTERFACE

  INTERFACE allocate
     module procedure mp_allocate, mp_allocate_real
  END INTERFACE

  INTERFACE attach
     MODULE PROCEDURE mp_attach_matrix, mp_attach_vector
     module procedure mp_attach_matrix_real, mp_attach_vector_real
  END INTERFACE

  INTERFACE set_value
     MODULE PROCEDURE mp_set_value, mp_set_value_real
  END INTERFACE

  INTERFACE add_value
     MODULE PROCEDURE mp_add_value, mp_add_value_real
  END INTERFACE

  interface add_to_row
     module procedure sfm_add_to_row, sfm_add_to_row_real
  end interface

  INTERFACE finalize
     MODULE PROCEDURE mp_finalize_matrix, mp_finalize_matrix_real
  END INTERFACE
  
  INTERFACE commit_values
     MODULE PROCEDURE mp_commit_values, mp_commit_values_real
  END INTERFACE

  INTERFACE get_global_matrix_locally
     module procedure mp_get_global_matrix_locally, mp_get_global_matrix_locally_real
  END INTERFACE

  INTERFACE mat_get_value
     MODULE PROCEDURE mp_get_value, mp_get_value_real
  END INTERFACE

  INTERFACE get_local_abs_square_sum
     module procedure mp_get_local_abs_square_sum, mp_get_local_abs_square_sum_real
  END INTERFACE

  INTERFACE dot_multiply
     MODULE PROCEDURE mp_dot_multiply_MatMat, mp_dot_multiply_MatVec
     module procedure mp_dot_multiply_MatMat_real, mp_dot_multiply_MatVec_real
  END INTERFACE

  INTERFACE dot_multiply_MatDagger_Mat
     MODULE PROCEDURE mp_dot_multiply_MatDagger_Mat, mp_dot_multiply_MatDagger_Mat_real
  END INTERFACE

  INTERFACE show
     MODULE PROCEDURE mp_show_on_screen, mp_show_in_file, mp_show_on_screen_real, mp_show_in_file_real
  END INTERFACE

  INTERFACE add_matrix
     MODULE PROCEDURE mp_add_to_matrix, mp_add_matrix, mp_add_to_matrix_real, mp_add_matrix_real
  END INTERFACE

  INTERFACE subtract_matrix
  !INTERFACE OPERATOR(-)
     MODULE PROCEDURE mp_subtract_matrix
     module procedure mp_subtract_from_matrix
     module procedure mp_subtract_matrix_real
     module procedure mp_subtract_from_matrix_real
  END INTERFACE

  interface multiply_matrix_with_scalar
  !INTERFACE OPERATOR(*)
     MODULE PROCEDURE mp_multiply_matrix_with_real
     module procedure mp_scale_matrix_by_real
     module procedure mp_multiply_matrix_with_real_real
     module procedure mp_scale_matrix_by_real_real
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE mp_assign_matrix, mp_assign_matrix_real
  END INTERFACE

  INTERFACE invert
     module procedure mp_invert_matrix, mp_invert_matrix_real
  END INTERFACE

  INTERFACE my_sum
     MODULE PROCEDURE my_matrix_sum, my_matrix_sum_real
  END INTERFACE

  INTERFACE LU_factor
     MODULE PROCEDURE mp_LU_factor, mp_LU_factor_real
  END INTERFACE

  INTERFACE LU_solve
     MODULE PROCEDURE mp_LU_solve, mp_LU_solve_real
  END INTERFACE

  interface exit_if_factored
     module procedure sfm_exit_if_factored, sfm_exit_if_factored_real
  end interface

  interface print_storage_details
     module procedure sfm_print_storage_details, sfm_print_storage_details_real
  end interface

  PRIVATE :: mp_set_value, mp_initialize_matrix, mp_finalize_matrix, mp_commit_values
  PRIVATE :: mp_set_value_real, mp_initialize_matrix_real, mp_finalize_matrix_real, mp_commit_values_real
  PRIVATE :: mp_attach_matrix, mp_attach_vector
  PRIVATE :: mp_attach_matrix_real, mp_attach_vector_real
  private :: mp_attach_matrix_general
  PRIVATE :: mp_attach_matrix_general_real
  PRIVATE :: mp_allocate, mp_dot_multiply_MatVec,mp_dot_multiply_MatMat
  PRIVATE :: mp_allocate_real, mp_dot_multiply_MatVec_real,mp_dot_multiply_MatMat_real
  PRIVATE :: mp_dot_multiply_MatDagger_Mat
  PRIVATE :: mp_dot_multiply_MatDagger_Mat_real
  PRIVATE :: mp_add_value, sfm_add_to_row,mp_show_on_screen, mp_show_in_file
  PRIVATE :: mp_add_value_real, sfm_add_to_row_real, mp_show_on_screen_real, mp_show_in_file_real
  PRIVATE :: mp_get_global_matrix_locally, mp_get_value
  PRIVATE :: mp_get_global_matrix_locally_real, mp_get_value_real

  PRIVATE :: mp_add_to_matrix, mp_add_matrix, mp_subtract_matrix, mp_subtract_from_matrix
  PRIVATE :: mp_add_to_matrix_real, mp_add_matrix_real, mp_subtract_matrix_real, mp_subtract_from_matrix_real
  PRIVATE :: mp_multiply_matrix_with_real, mp_scale_matrix_by_real
  PRIVATE :: mp_multiply_matrix_with_real_real, mp_scale_matrix_by_real_real
  private :: mp_assign_matrix
  private :: mp_assign_matrix_real

  PRIVATE :: mp_get_local_abs_square_sum
  PRIVATE :: mp_get_local_abs_square_sum_real
  private :: mp_invert_matrix
  private :: mp_invert_matrix_real
  private :: mp_LU_factor, mp_LU_solve
  private :: mp_LU_factor_real, mp_LU_solve_real
  private :: get_local_linear_index
  private :: get_local_linear_index_real
  private :: mp_isInitialized, sfm_print_storage_details
  private :: mp_isInitialized_real, sfm_print_storage_details_real
  private :: my_matrix_sum, my_matrix_sum_real

CONTAINS
  FUNCTION static_size_of_storefullmatrixobject() RESULT(memory_need)
    INTEGER :: memory_need

    memory_need = 17*SIZE_OF_INTEGER &
         & + 3*SIZE_OF_LOGICAL &
         & + SIZE_OF_POINTER &
         & + size_of_processgrid()
  END FUNCTION static_size_of_storefullmatrixobject

#ifdef OLD_PROCESSGRID  
  SUBROUTINE mp_initialize_matrix(mat,n_rows,n_cols,PG)
    TYPE(StoreFullMatrixObject), intent(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    type(ProcessGrid) :: PG

    ! Local variables
    INTEGER :: info, minimal_blocksize

  
    mat%NRows = n_rows
    mat%NCols = n_cols
    mat%PG = PG
    mat%isVector = (n_cols.EQ.1)
    ! Distribute matrix on the grid
    mat%RowsPerBlock = mat%NRows / PG%NPRows
    IF (mat%NCols.LE.mat%PG%NPCols) THEN
       mat%ColsPerBlock = 1
    ELSE
       mat%ColsPerBlock = mat%NCols/mat%PG%NPCols
    END IF

    if (.not.mat%isVector) then
       ! a matrix should have square blocks
       minimal_blocksize = MIN(mat%RowsPerBlock,mat%ColsPerBlock)
       mat%RowsPerBlock = minimal_blocksize
       mat%ColsPerBlock = minimal_blocksize
    END IF
    mat%Blocksize = mat%RowsPerBlock*mat%ColsPerBlock
    mat%BlockrowsPerProcess = mat%NRows/(PG%NPRows*mat%RowsPerBlock)
    !PRINT*,"NPCols = ", PG%NPCols,", ColsPerBlock = ", mat%ColsPerBlock
    mat%BlockcolsPerProcess = mat%NCols/(PG%NPCols*mat%ColsPerBlock)
    mat%NBlocks = mat%BlockrowsPerProcess*mat%BlockcolsPerProcess
    !mat%Desc=1000
    CALL descinit( mat%Desc, mat%NRows, mat%NCols, &
         & mat%RowsPerBlock, mat%ColsPerBlock, 0, 0, PG%Context, mat%NRows/PG%NPRows, info )

    !PRINT*,mat%Desc
    IF (info.NE.0) THEN
       PRINT*,'DESC_mat, info = ', info
       STOP
    END IF

    mat%isInitialized = .TRUE.
    mat%isAttached = .FALSE.
    mat%isFactored = .FALSE.

    ! Get pointers into a defined state
    NULLIFY(mat%localMatrix)
    NULLIFY(mat%factoredBand)
    NULLIFY(mat%ipiv)
    NULLIFY(mat%rhs)

    CALL blacs_barrier(PG%Context,'A')
  END SUBROUTINE mp_initialize_matrix
#else
  SUBROUTINE mp_initialize_matrix(mat,n_rows,n_cols,PG)
    TYPE(StoreFullMatrixObject), intent(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    type(ProcessGrid) :: PG

    ! Local variables
    INTEGER :: info, minimal_blocksize

  
    mat%NRows = n_rows
    mat%NCols = n_cols
    mat%PG = PG
    mat%isVector = (n_cols.EQ.1)
    ! Distribute matrix on the grid
    IF (mat%isVector) THEN
       mat%RowsPerBlock = mat%NRows/mat%PG%NProcs
    ELSE
       mat%RowsPerBlock = mat%NRows
    END IF
    mat%ColsPerBlock = MAX(1,mat%NCols/mat%PG%NProcs)

    if (.not.mat%isVector) then
       ! a matrix should have square blocks
       minimal_blocksize = MIN(mat%RowsPerBlock,mat%ColsPerBlock)
       mat%RowsPerBlock = minimal_blocksize
       mat%ColsPerBlock = minimal_blocksize
    END IF
    mat%Blocksize = mat%RowsPerBlock*mat%ColsPerBlock
    IF (mat%isVector) THEN
       mat%BlockrowsPerProcess = 1
    ELSE
       mat%BlockrowsPerProcess = mat%NRows/mat%RowsPerBlock
    END IF
    mat%BlockcolsPerProcess = MAX(1,mat%NCols/(PG%NProcs*mat%ColsPerBlock))
    mat%NBlocks = mat%BlockrowsPerProcess*mat%BlockcolsPerProcess
    CALL descinit( mat%Desc, mat%NRows, mat%NCols, &
         & mat%RowsPerBlock, mat%ColsPerBlock, 0, 0, PG%Context, &
         & mat%RowsPerBlock*mat%BlockrowsPerProcess, info )

    !PRINT*,mat%Desc
    IF (info.NE.0) THEN
       PRINT*,'DESC_mat, info = ', info
       STOP
    END IF

    mat%isInitialized = .TRUE.
    mat%isAttached = .FALSE.
    mat%isFactored = .FALSE.

    ! Get pointers into a defined state
    NULLIFY(mat%localMatrix)
    NULLIFY(mat%factoredBand)
    NULLIFY(mat%ipiv)
    NULLIFY(mat%rhs)

    CALL blacs_barrier(PG%Context,'A')
  END SUBROUTINE mp_initialize_matrix
#endif

  LOGICAL FUNCTION mp_isInitialized(mat) 
    type(StoreFullMatrixObject) :: mat

    mp_isInitialized = mat%isInitialized
  END FUNCTION mp_isInitialized


  SUBROUTINE mp_allocate(mat)
    TYPE(StoreFullMatrixObject) :: mat

    !PRINT*,"Allocating localMatrix with ",mat%NBlocks*mat%Blocksize," complex entries."
    ALLOCATE(mat%localMatrix(mat%NBlocks*mat%Blocksize))
    mat%isAttached = .FALSE.
  END SUBROUTINE mp_allocate

  SUBROUTINE mp_attach_vector(mat,localarray)
    type(StoreFullMatrixObject) :: mat
    COMPLEX, DIMENSION(:),TARGET :: localarray

    CALL exit_if_factored(mat)
    CALL mp_attach_matrix_general(mat,localarray)
  END SUBROUTINE mp_attach_vector

  SUBROUTINE mp_attach_matrix(mat,localarray)
    type(StoreFullMatrixObject) :: mat
    COMPLEX, DIMENSION(:,:),TARGET :: localarray

    CALL exit_if_factored(mat)
    CALL mp_attach_matrix_general(mat,localarray)
  END SUBROUTINE mp_attach_matrix

  SUBROUTINE mp_attach_matrix_general(mat,localarray)
    TYPE(StoreFullMatrixObject) :: mat
    COMPLEX,DIMENSION(*),TARGET :: localarray

    CALL exit_if_factored(mat)
    mat%localMatrix => localarray(1:mat%NBlocks*mat%Blocksize)
    mat%isAttached = .TRUE.

  END SUBROUTINE mp_attach_matrix_general

  SUBROUTINE mp_finalize_matrix(mat)
    type(StoreFullMatrixObject) :: mat

    IF (mat%isAttached) THEN
       ! do nothing
       mat%isAttached=.FALSE.
       IF (ASSOCIATED(mat%localMatrix)) &
            & NULLIFY(mat%localMatrix)
    ELSE
       IF (ASSOCIATED(mat%localMatrix) ) DEALLOCATE(mat%localMatrix)
       IF (ASSOCIATED(mat%factoredBand)) DEALLOCATE(mat%factoredBand)
       IF (ASSOCIATED(mat%ipiv)        ) DEALLOCATE(mat%ipiv)
       IF (ASSOCIATED(mat%rhs)         ) DEALLOCATE(mat%rhs)
    END IF
    mat%isInitialized=.false.
  END SUBROUTINE mp_finalize_matrix

#ifdef OLD_PROCESSGRID
  SUBROUTINE mp_set_value(mat,global_row,global_col,value)
    TYPE(StoreFullMatrixObject),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: global_row,global_col
    COMPLEX, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc_r, proc_c, block_r, block_c, row_index, col_index
    INTEGER :: linear_index !, linear_blockindex

    DEBSTART("set_value(SFMO)")
    CALL exit_if_factored(mat)
    !PRINT*,"In routine set_value of StoreFullMatrixModule on PRow,PCol = ",mat%PRow,mat%PCol
    ! calculate the processor coordinates (starting with 0) of the given global indices
    proc_r = MOD((global_row-1)/mat%RowsPerBlock,mat%PG%NPRows)
    proc_c = MOD((global_col-1)/mat%ColsPerBlock,mat%PG%NPCols)
    !PRINT*,mat%Rank,mat%PRow,mat%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF ((proc_r.EQ.mat%PG%PRow) .AND. (proc_c.EQ.mat%PG%PCol)) THEN
       ! calculate the local block indices (starting from 0) in which the global indices lie
       block_r = (global_row-1)/(mat%PG%NPRows*mat%RowsPerBlock)
       block_c = (global_col-1)/(mat%PG%NPCols*mat%ColsPerBlock)
       ! as a last step, we calculate the indices (starting from 1) in the blocks
       row_index = MOD(global_row-1,mat%RowsPerBlock)+1
       col_index = MOD(global_col-1,mat%ColsPerBlock)+1
       !PRINT*,mat%Rank,"Setting (",row_index,",",col_index,") in block ",block_r,block_c
       !IF (mat%isVector) THEN
       !   mat%localVector(row_index) = value
       !ELSE
       !mat%localMatrix(row_index,col_index) = value
       ! calculate the linear index (starting from 1) in localMatrix
       !linear_blockindex = block_c*mat%BlockrowsPerProcess + block_r
       !linear_index = linear_blockindex*mat%Blocksize &
       !     &+ (col_index-1)*mat%RowsPerBlock+row_index
       linear_index = get_local_linear_index(mat,block_r,block_c,row_index,col_index)
       mat%localMatrix(linear_index) = value
       !PRINT*,"Global = (",global_row,global_col,"), P = (",proc_r,proc_c,"), Bl = (",block_r,block_c,"), ",&
        !    &row_index,col_index,", LI = ",linear_index, &
         !   &get_local_linear_index(mat,block_r,block_c,row_index,col_index)
       !END IF
    END IF
    DEBEND("set_value(SFMO)")
  END SUBROUTINE mp_set_value

  SUBROUTINE mp_add_value(mat,global_row,global_col,value)
    TYPE(StoreFullMatrixObject) :: mat
    INTEGER, INTENT(IN) :: global_row,global_col
    COMPLEX, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc_r, proc_c, block_r, block_c, row_index, col_index, linear_index

    DEBSTART("add_value(SFMO)")

    if (mat%isFactored) stop "Attempt to acces matrix entries of a factored matrix."
    !CALL exit_if_factored(mat)

    ! calculate the processor coordinates of the given global indices
    proc_r = MOD((global_row-1)/mat%RowsPerBlock,mat%PG%NPRows)
    proc_c = MOD((global_col-1)/mat%ColsPerBlock,mat%PG%NPCols)

    ! only process the entry if we are on the right processor
    IF ((proc_r.EQ.mat%PG%PRow) .AND. (proc_c.EQ.mat%PG%PCol)) THEN
       ! calculate the local block indices in which the global indices lie
       block_r = (global_row-1)/(mat%PG%NPRows*mat%RowsPerBlock)
       block_c = (global_col-1)/(mat%PG%NPCols*mat%ColsPerBlock)
       ! as a last step, we calculate the indices in the blocks
       row_index = MOD(global_row-1,mat%RowsPerBlock)+1
       col_index = MOD(global_col-1,mat%ColsPerBlock)+1

       !IF (mat%isVector) THEN
       !   mat%localVector(row_index) = mat%localVector(row_index) + value
       !ELSE
       !mat%localMatrix(row_index,col_index) = mat%localMatrix(row_index,col_index) + value
       linear_index = get_local_linear_index(mat,block_r,block_c,row_index,col_index)
       mat%localMatrix(linear_index) = mat%localMatrix(linear_index) + value
       !mat%localMatrix((col_index-1)*mat%RowsPerBlock+row_index) = &
       !     &mat%localMatrix((col_index-1)*mat%RowsPerBlock+row_index) + value
       !END IF
    END IF
    DEBEND("add_value(SFMO)")
    
  END SUBROUTINE mp_add_value

  !> Add the given array to the indicated row in the matrix
  !!
  !! This routine only works, if the processgrid is 1D and 
  !! PG\%NPCols=1. Only then a whole row of the matrix is local to one process.
  !! \param mat The StoreFullMatrixObject to work on
  !! \param global_row The global index of the row, starting with 1 up to mat%NRows
  !! \param whole_row An array containing the values, which will be added to the already
  !! existing values in row global_row of the matrix mat.
  subroutine sfm_add_to_row(mat,global_row,whole_row)
    TYPE(StoreFullMatrixObject) :: mat
    INTEGER, INTENT(IN) :: global_row
    COMPLEX, dimension(:),INTENT(IN) :: whole_row

    ! Local variables
    complex, parameter :: complex_one=(1.0D0,0.0D0)
    INTEGER :: proc_r, proc_c, row_index

    DEBSTART("add_value(SFMO)")

    if (mat%isFactored) stop "Attempt to acces matrix entries of a factored matrix."
    !CALL exit_if_factored(mat)

    if (mat%PG%NPCols.ne.1) stop "sfm_add_to_row can only be used if PG%NPCols.eq.1"

    ! calculate the processor coordinates of the given global indices
    proc_r = MOD((global_row-1)/mat%RowsPerBlock,mat%PG%NPRows)
    ! as per default PG%NPCols=1, proc_c is always 0, only if the
    ! 2D cyclic block distribution is really 2D, we have to reuse this line
    !proc_c = MOD((global_col-1)/mat%ColsPerBlock,mat%PG%NPCols)
    proc_c = 0

    ! only process the entry if we are on the right processor
    !IF ((proc_r.EQ.mat%PG%PRow) .AND. (proc_c.EQ.mat%PG%PCol)) THEN
    IF (proc_r.EQ.mat%PG%PRow) THEN
       ! calculate the local block indices in which the global indices lie
       !block_r = (global_row-1)/(mat%PG%NPRows*mat%RowsPerBlock)
       !block_c = (global_col-1)/(mat%PG%NPCols*mat%ColsPerBlock)
       ! as a last step, we calculate the indices in the blocks
       row_index = MOD(global_row-1,mat%RowsPerBlock)+1
       !col_index = MOD(global_col-1,mat%ColsPerBlock)+1

       !linear_index = get_local_linear_index(mat,block_r,block_c,row_index,col_index)
       call caxpy(mat%NCols,complex_one,whole_row,1,mat%localMatrix(row_index),mat%RowsPerBlock)

       !mat%localMatrix(linear_index) = mat%localMatrix(linear_index) + value
    END IF
    DEBEND("add_value(SFMO)")
    
  END SUBROUTINE sfm_add_to_row

  FUNCTION mp_get_value(mat,irow,icol) RESULT(value)
    TYPE(StoreFullMatrixObject) :: mat
    INTEGER :: irow, icol
    complex :: value

    ! local variables
    INTEGER :: proc_r, proc_c, block_r, block_c, row_index, col_index, linear_index

    !CALL exit_if_factored(mat)
    proc_r = MOD((irow-1)/mat%RowsPerBlock,mat%PG%NPRows)
    proc_c = MOD((icol-1)/mat%ColsPerBlock,mat%PG%NPCols)
    !PRINT*,mat%Rank,mat%PRow,mat%PCol,proc_r, proc_c
    ! return the entry, if we are on the right processor
    IF ((proc_r.EQ.mat%PG%PRow) .AND. (proc_c.EQ.mat%PG%PCol)) THEN
       ! calculate the local block indices in which the global indices lie
       block_r = (irow-1)/(mat%PG%NPRows*mat%RowsPerBlock)
       block_c = (icol-1)/(mat%PG%NPCols*mat%ColsPerBlock)
       ! as a last step, we calculate the indices in the blocks
       row_index = MOD(irow-1,mat%RowsPerBlock)+1
       col_index = MOD(icol-1,mat%ColsPerBlock)+1
       !PRINT*,mat%Rank,"Setting (",row_index,",",col_index,") in block ",block_r,block_c
       linear_index = get_local_linear_index(mat,block_r,block_c,row_index,col_index)
       value = mat%localMatrix(linear_index)
    ELSE
       ! return nothing
       value = CMPLX(1000,700000)
       ! this can be arbitrary, as it is never used in further calculations
       PRINT*,"This line should never be reached in storefullmatrix.F90."
       !CALL tracebackqq
    END IF
  END FUNCTION mp_get_value

#else
  SUBROUTINE mp_set_value(mat,irow,icol,value)
    TYPE(StoreFullMatrixObject),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc, block_r, block_c, row_index, col_index
    INTEGER :: linear_index !, linear_blockindex

    DEBSTART("set_value(SFMO)")
    CALL exit_if_factored(mat)
    !PRINT*,"In routine set_value of StoreFullMatrixModule on PRow,PCol = ",mat%PRow,mat%PCol
    ! calculate the processor coordinates (starting with 0) of the given global indices
    IF (mat%isVector) THEN
       proc = MOD((irow-1)/mat%RowsPerBlock,mat%PG%NProcs)
    ELSE
       proc = MOD((icol-1)/mat%ColsPerBlock,mat%PG%NProcs)
    END IF
    !PRINT*,mat%Rank,mat%PRow,mat%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF (proc.EQ.mat%PG%Rank) THEN
       ! calculate the local block indices (starting from 0) in which the global indices lie
       IF (mat%isVector) THEN
          !blocknr = (irow-1)/(mat%PG%NProcs*mat%RowsPerBlock) is always zero
          row_index = MOD((irow-1),mat%RowsPerBlock)+1
          !PRINT*,mat%PG%Rank,": row_index = ", row_index,LBOUND(mat%localMatrix),UBOUND(mat%localMatrix)
          mat%localMatrix(row_index) = value
       ELSE
          block_r = (irow-1)/mat%RowsPerBlock
          block_c = (icol-1)/(mat%PG%NProcs*mat%ColsPerBlock)
          row_index = MOD(irow-1,mat%RowsPerBlock)+1
          col_index = MOD(icol-1,mat%ColsPerBlock)+1

          linear_index = get_local_linear_index(mat,block_r,block_c,row_index,col_index)
          mat%localMatrix(linear_index) = value
       END IF
    END IF
    DEBEND("set_value(SFMO)")
  END SUBROUTINE mp_set_value

  SUBROUTINE mp_add_value(mat,irow,icol,value)
    TYPE(StoreFullMatrixObject) :: mat
    INTEGER, INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc, block_r, block_c, row_index, col_index, linear_index

    DEBSTART("add_value(SFMO)")

    CALL exit_if_factored(mat)
    ! calculate the processor coordinates of the given global indices
    IF (mat%isVector) THEN
       proc = MOD((irow-1)/mat%RowsPerBlock,mat%PG%NProcs)
    ELSE
       proc = MOD((icol-1)/mat%ColsPerBlock,mat%PG%NProcs)
    END IF
    !PRINT*,mat%Rank,mat%PRow,mat%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF (proc.EQ.mat%PG%Rank) THEN
       ! calculate the local block indices (starting from 0) in which the global indices lie
       IF (mat%isVector) THEN
          !blocknr = (irow-1)/(mat%PG%NProcs*mat%RowsPerBlock) is always zero
          row_index = MOD((irow-1),mat%RowsPerBlock)+1
          mat%localMatrix(row_index) = mat%localMatrix(row_index) + value
       ELSE
          block_r = (irow-1)/mat%RowsPerBlock
          block_c = (icol-1)/(mat%PG%NProcs*mat%ColsPerBlock)
          row_index = MOD(irow-1,mat%RowsPerBlock)+1
          col_index = MOD(icol-1,mat%ColsPerBlock)+1

          linear_index = get_local_linear_index(mat,block_r,block_c,row_index,col_index)
          mat%localMatrix(linear_index) = mat%localMatrix(linear_index) + value
       END IF
    END IF
    DEBEND("add_value(SFMO)")
    
  END SUBROUTINE mp_add_value

  FUNCTION mp_get_value(mat,irow,icol) RESULT(value)
    TYPE(StoreFullMatrixObject) :: mat
    INTEGER :: irow, icol
    complex :: value

    ! local variables
    INTEGER :: proc, block_r, block_c, row_index, col_index, linear_index

    CALL exit_if_factored(mat)
    IF (mat%isVector) THEN
       proc = MOD((irow-1)/mat%RowsPerBlock,mat%PG%NProcs)
    ELSE
       proc = MOD((icol-1)/mat%ColsPerBlock,mat%PG%NProcs)
    END IF
    !PRINT*,mat%Rank,mat%PRow,mat%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF (proc.EQ.mat%PG%Rank) THEN
       ! calculate the local block indices (starting from 0) in which the global indices lie
       IF (mat%isVector) THEN
          !blocknr = (irow-1)/(mat%PG%NProcs*mat%RowsPerBlock) is always zero
          row_index = MOD((irow-1),mat%RowsPerBlock)+1
          value = mat%localMatrix(row_index)
       ELSE
          block_r = (irow-1)/mat%RowsPerBlock
          block_c = (icol-1)/(mat%PG%NProcs*mat%ColsPerBlock)
          row_index = MOD(irow-1,mat%RowsPerBlock)+1
          col_index = MOD(icol-1,mat%ColsPerBlock)+1

          linear_index = get_local_linear_index(mat,block_r,block_c,row_index,col_index)
          value = mat%localMatrix(linear_index)
       END IF
    ELSE
       ! return nothing
       value = CMPLX(1000,700000)
       ! this can be arbitrary, as it is never used in further calculations
    END IF
  END FUNCTION mp_get_value

#endif

  SUBROUTINE mp_commit_values(mat)
    TYPE(StoreFullMatrixObject) :: mat

    ! nothing happens, as a Full matrix can not be further optimized

  END SUBROUTINE mp_commit_values

  !>Compute Matrix Matrix multiplication
  !! \param mat first matrix
  !! \param vec second matrix or vector
  !! \param res result
  !! \param transA_in optional flag for operator on 1st matrix ('N'=nothing,'T'=transpose, 'C'=transpose/conjugate)
  !! \param transB_in optional flag for operator on 2nd matrix ('N'=nothing,'T'=transpose, 'C'=transpose/conjugate)
  SUBROUTINE mp_dot_multiply_MatMat(mat, vec, res, transA_in, transB_in)
    TYPE(StoreFullMatrixObject),INTENT(IN) :: mat
    TYPE(StoreFullMatrixObject),INTENT(IN) :: vec
    TYPE(StoreFullMatrixObject),INTENT(INOUT) :: res
    CHARACTER, INTENT(IN), OPTIONAL :: transA_in, transB_in

    !Local variables
    CHARACTER :: transA='N', transB='N'
    COMPLEX,PARAMETER :: complex_one=(1,0), complex_zero=(0,0)
    LOGICAL,PARAMETER :: original=.TRUE.

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
    IF ((vec%isVector).AND.(res%isVector)) THEN
       if (transB.ne.'N') then
          stop "transB ne 'N' is not implemented for vectors in storefullmatrix/mp_dot_multiply_MatMat"
       endif
       IF (original) THEN
          res%localMatrix = complex_zero
          CALL pcgemv(transA, mat%Nrows, mat%NCols, &
               & complex_one, mat%localMatrix, 1, 1, mat%Desc, &
               & vec%localMatrix, 1, 1, vec%Desc, 1, &
               & complex_zero, res%localMatrix, 1, 1, res%Desc, 1)
       ELSE
          CALL pcgemv(transA, mat%Nrows, 16, &
               & complex_one, mat%localMatrix, 1, 1, mat%Desc, &
               & vec%localMatrix, 1, 1, vec%Desc, 1, &
               & complex_zero, res%localMatrix, 1, 1, res%Desc, 1)
          PRINT*,mat%PG%Rank,SUM(res%localMatrix*CONJG(res%localMatrix))
       END IF
    ELSE
       CALL pcgemm(transA, transB, mat%NRows, res%NCols, mat%NCols, &
            & complex_one, mat%localMatrix, 1, 1, mat%Desc, &
            & vec%localMatrix, 1, 1, vec%Desc,  &
            & complex_zero, res%localMatrix, 1, 1, res%Desc)
    END IF

  END SUBROUTINE mp_dot_multiply_MatMat

  SUBROUTINE mp_dot_multiply_MatDagger_Mat(mat1, mat2, res)
    TYPE(StoreFullMatrixObject),INTENT(IN) :: mat1, mat2
    TYPE(StoreFullMatrixObject),INTENT(INOUT) :: res

    call mp_dot_multiply_MatMat(mat1,mat2,res,'C')

  END SUBROUTINE mp_dot_multiply_MatDagger_Mat


  SUBROUTINE mp_dot_multiply_MatVec(mat, vec, res, transA_in, transB_in)
    TYPE(StoreFullMatrixObject),INTENT(IN) :: mat
    TYPE(StoreVectorObject),INTENT(IN) :: vec
    TYPE(StoreVectorObject),INTENT(INOUT) :: res
    CHARACTER, INTENT(IN), OPTIONAL :: transA_in, transB_in

    !Local variables
    CHARACTER :: transA='N', transB='N'
    COMPLEX,PARAMETER :: complex_one=(1,0), complex_zero=(0,0)
    LOGICAL,PARAMETER :: original=.TRUE.
    !COMPLEX, DIMENSION(res%NBlocks*res%Blocksize) :: temporary_array

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

    if (transB.ne.'N') then
       stop "transB ne 'N' is not implemented in storefullmatrix/mp_dot_multiply_MatVec"
    endif

    CALL exit_if_factored(mat)
  
    IF (original) THEN
       CALL pcgemv(transA, mat%Nrows, mat%NCols, &
            & complex_one, mat%localMatrix, 1, 1, mat%Desc, &
            & vec%localVector, 1, 1, vec%Desc, 1, &
            & complex_zero, res%localVector, 1, 1, res%Desc, 1)
    ELSE
       CALL pcgemv(transA, mat%Nrows, 16, &
            & complex_one, mat%localMatrix, 1, 1, mat%Desc, &
            & vec%localVector, 1, 1, vec%Desc, 1, &
            & complex_zero, res%localVector, 1, 1, res%Desc, 1)
       PRINT*,mat%PG%Rank,SUM(res%localVector*CONJG(res%localVector))
    END IF
    
  END SUBROUTINE mp_dot_multiply_MatVec

  !*****************************************************************
  !* The global matrix, which is distributed over all processes
  !* is gathered and stored in a local array. This routine is
  !* mainly used for debugging and backward compatibility. In the
  !* final production code, it should not be called as it does not
  !* scale in memory.
  !*****************************************************************
  SUBROUTINE mp_get_global_matrix_locally(mat,localfullmat)
    TYPE(StoreFullMatrixObject),INTENT(IN) :: mat
    COMPLEX, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    ! local variables
    INTEGER :: ierr,iBlockcol, iBlockrow,iCol,iRow,linear_blockindex, linear_index
    COMPLEX, DIMENSION(:,:),ALLOCATABLE :: localtransp,localtranspfullmat, localtemp
    !INTEGER :: desctemp(9)

    CALL exit_if_factored(mat)
    IF (mat%isVector) THEN
       ! If mat is a Vector, just concatenate the local vectors on all processes.
       CALL mpi_allgather(mat%localMatrix,mat%NBlocks*mat%Blocksize,MPI_COMPLEX_TYPE,&
            &localfullmat,mat%NBlocks*mat%Blocksize,&
            &MPI_COMPLEX_TYPE,mat%PG%Communicator,ierr)
    ELSE
       ! otherwise, we have to take into account the cyclic block matrix
       ! structure, which complicates the operation.
       ! First we allocate a matrix to take the transposed local matrix
       ALLOCATE(localtransp(mat%ColsPerBlock*mat%BlockcolsPerProcess,&
            &mat%RowsPerBlock*mat%BlockrowsPerProcess))
       ! and a global matrix to store the transposed global matrix.
       ALLOCATE(localtranspfullmat(mat%NCols,mat%NRows))

       ! additionally we need a local matrix untransposed to store the 
       ! correctly reordered local matrix
       ALLOCATE(localtemp(mat%RowsPerBlock*mat%BlockrowsPerProcess,&
            &mat%ColsPerBlock*mat%BlockcolsPerProcess))

       ! and we fill this local matrix
       DO iBlockcol = 0,mat%BlockcolsPerProcess-1
          DO iCol=1,mat%ColsPerBlock
             DO iBlockrow = 0,mat%BlockrowsPerProcess-1
                DO iRow=1,mat%RowsPerBlock
                   linear_blockindex = iBlockcol*mat%BlockrowsPerProcess + iBlockrow
                   linear_index = get_local_linear_index(mat,iBlockrow,iBlockcol,iRow,iCol)
                   !linear_blockindex*mat%Blocksize &
                   !     &+ (iCol-1)*mat%RowsPerBlock+iRow

                   localtemp(iBlockrow*mat%RowsPerBlock+iRow,&
                        &iBlockcol*mat%ColsPerBlock+iCol) = mat%localMatrix(linear_index)
                END DO
             END DO
          END DO
       END DO

       localtransp = transpose(localtemp)
       !CALL cgetmo(localtemp,mat%RowsPerBlock*mat%BlockrowsPerProcess,&
       !     &mat%RowsPerBlock*mat%BlockrowsPerProcess,&
       !     &mat%ColsPerBlock*mat%BlockcolsPerProcess,&
       !     &localtransp,mat%ColsPerBlock*mat%BlockcolsPerProcess)

       !PRINT*, "Communicator in MPI_Allgather is ", mat%PG%Communicator
       CALL mpi_allgather(localtransp,mat%NBlocks*mat%Blocksize,MPI_COMPLEX_TYPE,&
            &localtranspfullmat,mat%NBlocks*mat%Blocksize,MPI_COMPLEX_TYPE,&
            &mat%PG%Communicator,ierr)
       !CALL mpi_allgather(localtransp,mat%RowsPerBlock*mat%ColsPerBlock,MPI_COMPLEX_TYPE,&
       !     &localtranspfullmat,mat%RowsPerBlock*mat%ColsPerBlock,MPI_COMPLEX_TYPE,&
       !     &mat%PG%Communicator,ierr)
       
       !CALL cgetmo(localtranspfullmat,mat%NCols,mat%NCols,mat%NRows,&
       !     &localfullmat,mat%NRows)
       localfullmat = TRANSPOSE(localtranspfullmat)
       DEALLOCATE(localtransp)
       DEALLOCATE(localtranspfullmat)
       deallocate(localtemp)
    END IF
  END SUBROUTINE mp_get_global_matrix_locally

  FUNCTION mp_get_local_abs_square_sum(mat) RESULT(value)
    type(StoreFullMatrixObject) :: mat
    REAL :: value

    CALL exit_if_factored(mat)
    value = SUM(REAL(mat%localMatrix)**2+AIMAG(mat%localMatrix)**2)

  END FUNCTION mp_get_local_abs_square_sum

  SUBROUTINE mp_show_on_screen(mat)
    type(StoreFullMatrixObject) :: mat

    ! Local variables
    INTEGER :: iProc,iBlockrow,iBlockcol,iCol, iRow

    CALL exit_if_factored(mat)
    ! do some debugging output
    CALL blacs_barrier(mat%PG%Context,'A')
    DO iProc=0,mat%PG%NProcs-1
       CALL blacs_barrier(mat%PG%Context,'A')
       IF (iProc.EQ.mat%PG%Rank) THEN
          PRINT*,"---- Start output Rank = ",mat%PG%Rank," ----"
          DO iBlockrow=0,mat%BlockrowsPerProcess-1
             PRINT*,"  -- Start output Blockrow = ",iBlockrow," --"
             DO iRow=1,mat%RowsPerBlock
                DO iBlockcol=0,mat%BlockcolsPerProcess-1
                   DO iCol = 1,mat%ColsPerBlock
                      WRITE(*,"(2ES9.2,1X)",advance='NO') &
                           &mat%localMatrix(get_local_linear_index(mat,iBlockrow,iBlockcol,iRow,iCol))
                   END DO
                END DO
                WRITE(*,"(A)") ""
             END DO
             PRINT*,"  -- End   output Blockrow = ",iBlockrow," --"
          END DO
          PRINT*,"---- End   output Rank = ",mat%PG%Rank," ----"
       END IF
    END DO
    CALL blacs_barrier(mat%PG%Context,'A')
  END SUBROUTINE mp_show_on_screen

  SUBROUTINE mp_show_in_file(mat,filename)
    type(StoreFullMatrixObject) :: mat
    CHARACTER(len=*) :: filename

    ! Local variables
    INTEGER :: iCol, iRow,iBlockrow,iBlockcol !,iProc
    character(len=FILENAME_MAX) :: full_filename

    CALL exit_if_factored(mat)
    ! open a file on each process
    WRITE(full_filename,"(A,I3.3,A)") filename,mat%PG%Rank,".dat"
    OPEN(77,file=TRIM(full_filename))

    ! do some debugging output
    CALL blacs_barrier(mat%PG%Context,'A')
    !DO iProc=0,mat%PG%NProcs-1
    !   CALL blacs_barrier(mat%PG%Context,'A')
    !   IF (iProc.EQ.mat%PG%Rank) THEN
          !PRINT*,"---- Start output Rank = ",mat%PG%Rank," ----"
          DO iBlockrow=0,mat%BlockrowsPerProcess-1
             !PRINT*,"  -- Start output Blockrow = ",iBlockrow," --"
             DO iRow=1,mat%RowsPerBlock
                !IF (mat%isVector) THEN
                   !WRITE(*,"(2ES10.2,X)") mat%localVector(iRow)
                !   WRITE(*,"(2ES10.2,X)") mat%localMatrix(iRow)
                !ELSE
                DO iBlockcol=0,mat%BlockcolsPerProcess-1
                   DO iCol = 1,mat%ColsPerBlock
                      WRITE(77,"(2ES20.10)",advance='NO') &
                           &mat%localMatrix(get_local_linear_index(mat,iBlockrow,iBlockcol,iRow,iCol))
                   END DO
!                   WRITE(77,"(A)",advance='NO') "|"
                END DO
                WRITE(77,"(A)") ""
                !END IF
             END DO
             !PRINT*,"  -- End   output Blockrow = ",iBlockrow," --"
          END DO
          !PRINT*,"---- End   output Rank = ",mat%PG%Rank," ----"
       !END IF
    !END DO
    CLOSE(77)
    CALL blacs_barrier(mat%PG%Context,'A')
  END SUBROUTINE mp_show_in_file

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix(this,mat)
    TYPE(StoreFullMatrixObject), INTENT(INOUT) :: this
    TYPE(StoreFullMatrixObject), INTENT(IN) :: mat

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)
    this%localMatrix = this%localMatrix + mat%localMatrix
  END SUBROUTINE mp_add_to_matrix

  SUBROUTINE mp_add_matrix(this,mat,res)
    TYPE(StoreFullMatrixObject), INTENT(IN) :: this,mat
    TYPE(StoreFullMatrixObject), INTENT(INOUT) :: res

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)
    CALL exit_if_factored(res)
    res%localMatrix = this%localMatrix+mat%localMatrix
  END SUBROUTINE mp_add_matrix

  SUBROUTINE mp_subtract_matrix(this,mat,res)
    TYPE(StoreFullMatrixObject),INTENT(IN) :: this,mat
    TYPE(StoreFullMatrixObject),intent(INOUT) :: res
    ! calculate res = this - mat

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)
    CALL exit_if_factored(res)
    res%localMatrix = this%localMatrix - mat%localMatrix
  END SUBROUTINE mp_subtract_matrix

  SUBROUTINE mp_subtract_from_matrix(this,mat)
    TYPE(StoreFullMatrixObject),INTENT(INOUT) :: this
    TYPE(StoreFullMatrixObject),intent(IN) :: mat
    ! calculate this = this - mat

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)
    this%localMatrix = this%localMatrix - mat%localMatrix
  END SUBROUTINE mp_subtract_from_matrix

  SUBROUTINE mp_multiply_matrix_with_real(this, scalar, res)
    TYPE(StoreFullMatrixObject), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(StoreFullMatrixObject),intent(INOUT) :: res

    CALL exit_if_factored(this)
    CALL exit_if_factored(res)
    res%localMatrix = this%localMatrix*scalar
  END SUBROUTINE mp_multiply_matrix_with_real

  SUBROUTINE mp_scale_matrix_by_real(this, scalar)
    TYPE(StoreFullMatrixObject), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL exit_if_factored(this)
    this%localMatrix = this%localMatrix*scalar
  END SUBROUTINE mp_scale_matrix_by_real


  SUBROUTINE mp_assign_matrix(lmat,rmat)
    TYPE(StoreFullMatrixObject), intent(INOUT) :: lmat
    TYPE(StoreFullMatrixObject), intent(IN) :: rmat

    CALL exit_if_factored(lmat)
    CALL exit_if_factored(rmat)
    IF (lmat%PG.EQ.rmat%PG) THEN
       lmat%localMatrix = rmat%localMatrix
    ELSE
       PRINT*,"Assignment of matrices is only possible within the same context."
    END IF
  END SUBROUTINE mp_assign_matrix

  ! ===================================================
  ! == Define some mathematical routines on matrices ==
  ! ===================================================
  SUBROUTINE mp_invert_matrix(mat)
    type(StoreFullMatrixObject) :: mat
    
    ! Local variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv, iwork
    COMPLEX, DIMENSION(:), allocatable :: work
    INTEGER :: info,lwork,liwork,ub_locr

    CALL exit_if_factored(mat)
    ! calculate upper bound of ipiv storage needed
    ! LOCr(M_A)+MB_A (ScaLapack User's guide)
    ub_locr = ( (mat%NRows+mat%RowsPerBlock-1)/mat%RowsPerBlock &
         & + mat%PG%NPRows-1 ) / mat%PG%NPRows * mat%RowsPerBlock

    !PRINT*,"ub_locr = ",ub_locr
    ALLOCATE(ipiv(ub_locr+mat%RowsPerBlock))
    !WRITE(*,"(16ES10.2)") mat%localMatrix
    ! First do the LU decomposition in place
    CALL pcgetrf(mat%NRows,mat%NCols,&
         &mat%localMatrix, 1, 1,&
         &mat%Desc,ipiv,info)

    IF (info.NE.0) THEN
       PRINT*,"pcgetrf produced an error in storefullmatrix/mp_invert_matrix"
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
       PRINT*,"pcgetri produced an error in storefullmatrix/mp_invert_matrix"
       stop
    END IF
    !PRINT*,"after second pcgetri call"
    !WRITE(*,"(16ES10.2)") mat%localMatrix

    DEALLOCATE(work,iwork)
    DEALLOCATE(ipiv)
  END SUBROUTINE mp_invert_matrix

  SUBROUTINE mp_LU_factor(mat)
    type(StoreFullMatrixObject) :: mat
    
    ! Local variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: ab
    INTEGER :: i, j, j_off, nbw(2), ierr, nstore

    CALL exit_if_factored(mat)

    ! This routine only works on a grid with distributed rows

    IF(mat%PG%NPCols /= 1 .OR. mat%PG%NPRows /= mat%PG%NProcs) STOP 'Illegal grid for LU_factor'

    ! Get bandwidth

    j_off = mat%PG%PRow * mat%RowsPerBlock

    nbw(1) = 0 ! lower bandwidth
    nbw(2) = 0 ! upper bandwidth
    DO i = 1, mat%NCols
    DO j = 1, mat%RowsPerBlock
       IF(mat%localMatrix((i-1)*mat%RowsPerBlock+j) /= 0.) THEN
          IF(j+j_off > i) THEN
             nbw(1) = MAX(nbw(1),(j+j_off)-i)
          ELSE
             nbw(2) = MAX(nbw(2),i-(j+j_off))
          ENDIF
       ENDIF
    ENDDO
    ENDDO

    CALL MPI_Allreduce(MPI_IN_PLACE, nbw, 2, MPI_INTEGER, MPI_MAX, mat%PG%Communicator, ierr)

    mat%bwl = nbw(1)
    mat%bwu = nbw(2)

    ! Get storage size of LU decomposition (neglecting ipiv)

    nstore = (2*mat%bwl+mat%bwu+1)*mat%RowsPerBlock

    ! The following is just for testing if inversion/LU have the same results
    ! Setting the environment variable GENE_FORCE_INVERT to some value
    ! will force inversion instead of LU decomposition.
    call GET_ENVIRONMENT_VARIABLE('GENE_FORCE_INVERT',status=ierr)
    if(ierr==0) nstore = 1000000000

    ! Check if it pays off to factor the matrix, otherways just use mp_invert_matrix

    IF(mat%NCols*mat%RowsPerBlock <= nstore) THEN

       CALL mp_invert_matrix(mat)

       mat%bwl = -1
       mat%bwu = -1
       mat%isFactored = .TRUE.
       RETURN
    ENDIF

    ! Get the data of the matrix bands on all procs.
    ! Doing it this way (using MPI_Allreduce) is not the most efficient way,
    ! but this routine is called only in the initialization phase.

    ALLOCATE(ab(2*mat%bwl+mat%bwu+1,mat%NCols))
    ab(:,:) = 0

    DO i = 1, mat%NCols
    DO j = 1, mat%RowsPerBlock
       IF((j+j_off)-i <= mat%bwl .AND. i-(j+j_off) <= mat%bwu) THEN
          ab(mat%bwl+mat%bwu+1+(j+j_off)-i,i) = mat%localMatrix((i-1)*mat%RowsPerBlock+j)
       ENDIF
    ENDDO
    ENDDO

    CALL MPI_Allreduce(MPI_IN_PLACE,ab,SIZE(ab),MPI_COMPLEX_TYPE,MPI_SUM,mat%PG%Communicator, ierr)

    ! Use LAPACK routine cgbtrf (on all processors with identical data)
    ! and store the results distributed.
    ! Again, this is not the most efficient way but very simple to do ...

    ALLOCATE(ipiv(mat%NCols))

    CALL cgbtrf(mat%NRows,mat%NCols,mat%bwl,mat%bwu,ab,ubound(ab,1),ipiv,ierr)

    IF (ierr.NE.0) THEN
       PRINT*,"cgbtrf info = ",ierr
       stop
    END IF

    allocate(mat%factoredBand(2*mat%bwl+mat%bwu+1,mat%RowsPerBlock))
    allocate(mat%ipiv(mat%RowsPerBlock))

    mat%factoredBand(:,1:mat%RowsPerBlock) = ab(:,j_off+1:j_off+mat%RowsPerBlock)
    mat%ipiv(1:mat%RowsPerBlock) = ipiv(j_off+1:j_off+mat%RowsPerBlock)

    DEALLOCATE(ipiv)
    DEALLOCATE(ab)

    ! mat%localMatrix is not needed any more
    DEALLOCATE(mat%localMatrix)

    mat%isFactored = .TRUE.

  END SUBROUTINE mp_LU_factor

  SUBROUTINE mp_LU_solve(mat, vec, res, phase)
    TYPE(StoreFullMatrixObject),INTENT(INOUT) :: mat
    TYPE(StoreFullMatrixObject),INTENT(IN) :: vec
    TYPE(StoreFullMatrixObject),INTENT(INOUT) :: res
    INTEGER, INTENT(IN) :: phase

    INTEGER :: j, l, lm, kd, j_off, j_max, j_min, ierr
    INTEGER :: mpi_status(mpi_status_size)
    COMPLEX :: ctmp
    COMPLEX,PARAMETER :: complex_one=(1,0), complex_zero=(0,0)

    ! This routine must be called 3 times (with phase=1,2,3) in order
    ! to solve L*U*x = b completly.
    ! For phase=1, a copy of vec is stored internally, res is not used
    ! For phase=2, the forward substitution is done, vec and res are not used
    ! For phase=3, the backward substitution is done, res is set to the solution, vec is not used
    !
    ! The reason for doing it this complicated way is that a single LU-solution
    ! does not parallelize since it is sweeping once forward and back through
    ! all processors.
    !
    ! For multiple LU solutions, however, the processors may overlap in doing work
    ! (the first one working on solution n, the next one on solution n-1 etc.)
    ! if first all forward substitutions are done and then all backward substitutions.

    ! Safety checks:

    IF(vec%NCols /= 1) stop 'LU_solve: NCols must be 1'
    IF(.NOT. mat%isFactored) stop 'LU_solve: not factored'
    IF(mat%PG%NPCols /= 1 .OR. mat%PG%NPRows /= mat%PG%NProcs) STOP 'Illegal grid for LU_solve'

    ! If it didn't pay off to use the LU factorization of the banded matrix,
    ! the inverse of the full matrix is stored and mat%bwl is < 0

    IF(mat%bwl<0) THEN

       IF(phase == 1) THEN
          ! Allocate mat%rhs and store the complete rhs on every proc
          ALLOCATE(mat%rhs(mat%NRows))
          CALL MPI_Allgather(vec%localMatrix,mat%RowsPerBlock,MPI_COMPLEX_TYPE,&
                             mat%rhs,        mat%RowsPerBlock,MPI_COMPLEX_TYPE,&
                             vec%PG%Communicator,ierr)
       ENDIF
       ! Nothing to do for phase == 2
       IF(phase == 3) THEN
          ! Multiply inverse matrix with stored vector
          CALL cgemv('N', mat%RowsPerBlock, mat%NCols, &
                     complex_one, mat%localMatrix, mat%RowsPerBlock, &
                     mat%rhs, 1, complex_zero, res%localMatrix, 1)
          DEALLOCATE(mat%rhs)
       ENDIF
       RETURN

    ENDIF

    j_off = mat%PG%PRow * mat%RowsPerBlock

    KD = mat%bwl+mat%bwu+1 ! Position of diagonal in mat%factoredBand

    IF(phase == 1) THEN

       ! Store a copy of vec (on proc 0)

       ALLOCATE(mat%rhs(mat%NRows))

       CALL MPI_Gather(vec%localMatrix,mat%RowsPerBlock,MPI_COMPLEX_TYPE,&
                       mat%rhs,        mat%RowsPerBlock,MPI_COMPLEX_TYPE,&
                       0,vec%PG%Communicator,ierr)
    ENDIF

    IF(phase == 2) THEN

       ! Forward substitution

       IF(mat%PG%PRow > 0) THEN
          ! Wait until previous proc sends mat%rhs
          j_min = MAX(j_off + 1 - (mat%bwl+mat%bwu), 1)
          j_max = mat%NRows
          CALL MPI_Recv(mat%rhs(j_min),j_max-j_min+1,MPI_COMPLEX_TYPE,mat%PG%PRow-1,111, &
                        vec%PG%Communicator, mpi_status, ierr)
       ENDIF

       ! My part of forward substitution

       DO j = j_off+1, j_off+mat%RowsPerBlock
          if(j==mat%NRows) exit
          lm = MIN( mat%bwl, mat%NRows-j )
          l = mat%ipiv(j-j_off)
          IF( l.NE.j ) then
             ctmp = mat%rhs(l)
             mat%rhs(l) = mat%rhs(j)
             mat%rhs(j) = ctmp
          endif
          mat%rhs(j+1:j+lm) = mat%rhs(j+1:j+lm) - mat%rhs(j)*mat%factoredBand(KD+1:KD+lm,j-j_off)
       enddo

       IF(mat%PG%PRow < mat%PG%NPRows-1) THEN
          ! Send mat%rhs to next proc
          j_min = MAX(j_off + mat%RowsPerBlock + 1 - (mat%bwl+mat%bwu), 1)
          j_max = mat%NRows
          CALL MPI_Send(mat%rhs(j_min),j_max-j_min+1,MPI_COMPLEX_TYPE,mat%PG%PRow+1,111, &
                        vec%PG%Communicator, ierr)
       ENDIF
    ENDIF

    IF(phase == 3) THEN

       ! Backward substitution

       IF(mat%PG%PRow < mat%PG%NPRows-1) THEN
          ! Wait until next proc sends back mat%rhs
          j_min = MAX(j_off + mat%RowsPerBlock + 1 - (mat%bwl+mat%bwu), 1)
          j_max = j_off + mat%RowsPerBlock
          CALL MPI_Recv(mat%rhs(j_min),j_max-j_min+1,MPI_COMPLEX_TYPE,mat%PG%PRow+1,222, &
                        vec%PG%Communicator, mpi_status, ierr)
       ENDIF

       ! My part of backward substitution

       DO j = j_off+mat%RowsPerBlock, j_off+1, -1
          IF (mat%rhs(j).NE.0.) THEN
             mat%rhs(j) = mat%rhs(j)/mat%factoredBand(KD,j-j_off)
             if(j==1) exit
             lm = MIN(mat%bwl+mat%bwu, j-1)
             mat%rhs(j-lm:j-1) = mat%rhs(j-lm:j-1) - mat%rhs(j)*mat%factoredBand(KD-lm:KD-1,j-j_off)
          END IF
       ENDDO

       IF(mat%PG%PRow > 0) THEN
          ! Send mat%rhs to previous proc
          j_min = MAX(j_off + 1 - (mat%bwl+mat%bwu), 1)
          j_max = j_off
          CALL MPI_Send(mat%rhs(j_min),j_max-j_min+1,MPI_COMPLEX_TYPE,mat%PG%PRow-1,222, &
                        vec%PG%Communicator, ierr)
       ENDIF

       ! Store local part of mat%rhs

       res%localMatrix(1:mat%RowsPerBlock) = mat%rhs(j_off+1:j_off+mat%RowsPerBlock)

       DEALLOCATE(mat%rhs)

    ENDIF

  END SUBROUTINE mp_LU_solve

  SUBROUTINE sfm_exit_if_factored(mat)
    type(StoreFullMatrixObject) :: mat

    INTEGER :: ierr

    ! This routines is called by all routines which use mat%localMatrix.
    ! It does exit the program if the matrix is LU-factored.

    IF(mat%isFactored) THEN
      PRINT *,'ERROR: Attempt use original matrix data of a LU-factored Matrix'
      CALL mpi_abort(MPI_COMM_WORLD, 1, ierr)
    ENDIF

  END SUBROUTINE sfm_exit_if_factored

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_matrix_sum(mat,communicator)
    TYPE(StoreFullMatrixObject), INTENT(INOUT) :: mat
    COMPLEX, DIMENSION(mat%NBlocks*mat%Blocksize) :: totalsum
    integer :: communicator

    integer :: ierr
    CALL mpi_allreduce(mat%localMatrix,totalsum,mat%NBlocks*mat%Blocksize, &
         & MPI_COMPLEX_TYPE, MPI_SUM, communicator, ierr)
    
    mat%localMatrix = totalsum
  END SUBROUTINE my_matrix_sum

  ! ==================================
  ! == only private helper routines ==
  ! ==================================
  FUNCTION get_local_linear_index(mat,blockrow,blockcol,irow,icol) RESULT(index)
    TYPE(StoreFullMatrixObject), intent(IN) :: mat
    INTEGER,intent(IN) :: blockrow,blockcol,irow,icol
    integer :: index
    
    ! calculate the linear index of the element with coordinates irow,icol in the block
    ! with the coordinates blockrow,blockcol of matrix mat
    !linear_blockindex = blockcol*mat%BlockrowsPerProcess + blockrow
    index = (blockcol*mat%BlockrowsPerProcess + blockrow)*mat%Blocksize &
         &+ (icol-1)*mat%RowsPerBlock+irow
  END FUNCTION get_local_linear_index

  subroutine sfm_print_storage_details(mat)
    type(StoreFullMatrixObject) :: mat

    integer :: rank,ierr

    if (.not.mat%isInitialized) then
       print*,"Matrix is NOT initialized."
       return
    end if

    call mpi_comm_rank(mat%PG%communicator,rank,ierr)
    if (rank.eq.0) then
       print*,"NRows,NCols = ",mat%NRows,mat%NCols
       print*,"RowsPerBlock,ColsPerBlock = ",mat%RowsPerBlock,mat%ColsPerBlock
       print*,"NBlocks,Blocksize = ", mat%NBlocks, mat%Blocksize
       if (mat%NBlocks.gt.1) then
          print*,"BlockrowsPerProcess,BlockcolsPerProcess",mat%BlockrowsPerProcess,mat%BlockcolsPerProcess
       end if
       
       if (mat%isAttached) then
          print*,"Matrix is attached."
       else
          print*,"Matrix is NOT attached."
       end if
    end if
  end subroutine sfm_print_storage_details

  ! ===========================================
  ! == Repeat EVERYTHING for REAL data types ==
  ! ===========================================
#ifdef OLD_PROCESSGRID
  SUBROUTINE mp_initialize_matrix_real(mat,n_rows,n_cols,PG)
    TYPE(StoreFullMatrixObjectReal), intent(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    type(ProcessGrid) :: PG

    ! Local variables
    INTEGER :: info, minimal_blocksize

  
    mat%NRows = n_rows
    mat%NCols = n_cols
    mat%PG = PG
    mat%isVector = (n_cols.EQ.1)
    ! Distribute matrix on the grid
    mat%RowsPerBlock = mat%NRows / PG%NPRows
    IF (mat%NCols.LE.mat%PG%NPCols) THEN
       mat%ColsPerBlock = 1
    ELSE
       mat%ColsPerBlock = mat%NCols/mat%PG%NPCols
    END IF

    if (.not.mat%isVector) then
       ! a matrix should have square blocks
       minimal_blocksize = MIN(mat%RowsPerBlock,mat%ColsPerBlock)
       mat%RowsPerBlock = minimal_blocksize
       mat%ColsPerBlock = minimal_blocksize
    END IF
    mat%Blocksize = mat%RowsPerBlock*mat%ColsPerBlock
    mat%BlockrowsPerProcess = mat%NRows/(PG%NPRows*mat%RowsPerBlock)
    !PRINT*,"NPCols = ", PG%NPCols,", ColsPerBlock = ", mat%ColsPerBlock
    mat%BlockcolsPerProcess = mat%NCols/(PG%NPCols*mat%ColsPerBlock)
    mat%NBlocks = mat%BlockrowsPerProcess*mat%BlockcolsPerProcess
    !mat%Desc=1000
    CALL descinit( mat%Desc, mat%NRows, mat%NCols, &
         & mat%RowsPerBlock, mat%ColsPerBlock, 0, 0, PG%Context, mat%NRows/PG%NPRows, info )

    !PRINT*,mat%Desc
    IF (info.NE.0) THEN
       PRINT*,'DESC_mat, info = ', info
       STOP
    END IF

    mat%isInitialized = .TRUE.
    mat%isAttached = .FALSE.
    mat%isFactored = .FALSE.

    ! Get pointers into a defined state
    NULLIFY(mat%localMatrix)
    NULLIFY(mat%factoredBand)
    NULLIFY(mat%ipiv)
    NULLIFY(mat%rhs)

    CALL blacs_barrier(PG%Context,'A')
  END SUBROUTINE mp_initialize_matrix_real
#else
  SUBROUTINE mp_initialize_matrix_real(mat,n_rows,n_cols,PG)
    TYPE(StoreFullMatrixObjectReal), intent(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    type(ProcessGrid) :: PG

    ! Local variables
    INTEGER :: info, minimal_blocksize

  
    mat%NRows = n_rows
    mat%NCols = n_cols
    mat%PG = PG
    mat%isVector = (n_cols.EQ.1)
    ! Distribute matrix on the grid
    IF (mat%isVector) THEN
       mat%RowsPerBlock = mat%NRows/mat%PG%NProcs
    ELSE
       mat%RowsPerBlock = mat%NRows
    END IF
    mat%ColsPerBlock = MAX(1,mat%NCols/mat%PG%NProcs)

    if (.not.mat%isVector) then
       ! a matrix should have square blocks
       minimal_blocksize = MIN(mat%RowsPerBlock,mat%ColsPerBlock)
       mat%RowsPerBlock = minimal_blocksize
       mat%ColsPerBlock = minimal_blocksize
    END IF
    mat%Blocksize = mat%RowsPerBlock*mat%ColsPerBlock
    IF (mat%isVector) THEN
       mat%BlockrowsPerProcess = 1
    ELSE
       mat%BlockrowsPerProcess = mat%NRows/mat%RowsPerBlock
    END IF
    mat%BlockcolsPerProcess = MAX(1,mat%NCols/(PG%NProcs*mat%ColsPerBlock))
    mat%NBlocks = mat%BlockrowsPerProcess*mat%BlockcolsPerProcess
    CALL descinit( mat%Desc, mat%NRows, mat%NCols, &
         & mat%RowsPerBlock, mat%ColsPerBlock, 0, 0, PG%Context, &
         & mat%RowsPerBlock*mat%BlockrowsPerProcess, info )

    !PRINT*,mat%Desc
    IF (info.NE.0) THEN
       PRINT*,'DESC_mat, info = ', info
       STOP
    END IF

    mat%isInitialized = .TRUE.
    mat%isAttached = .FALSE.
    mat%isFactored = .FALSE.

    ! Get pointers into a defined state
    NULLIFY(mat%localMatrix)
    NULLIFY(mat%factoredBand)
    NULLIFY(mat%ipiv)
    NULLIFY(mat%rhs)

    CALL blacs_barrier(PG%Context,'A')
  END SUBROUTINE mp_initialize_matrix_real
#endif

  LOGICAL FUNCTION mp_isInitialized_real(mat) 
    type(StoreFullMatrixObjectReal) :: mat

    mp_isInitialized_real = mat%isInitialized
  END FUNCTION mp_isInitialized_real


  SUBROUTINE mp_allocate_real(mat)
    TYPE(StoreFullMatrixObjectReal) :: mat

    !PRINT*,"Allocating localMatrix with ",mat%NBlocks*mat%Blocksize," complex entries."
    ALLOCATE(mat%localMatrix(mat%NBlocks*mat%Blocksize))
    mat%isAttached = .FALSE.
  END SUBROUTINE mp_allocate_real

  SUBROUTINE mp_attach_vector_real(mat,localarray)
    type(StoreFullMatrixObjectReal) :: mat
    REAL, DIMENSION(:),TARGET :: localarray

    CALL exit_if_factored(mat)
    CALL mp_attach_matrix_general_real(mat,localarray)
  END SUBROUTINE mp_attach_vector_real

  SUBROUTINE mp_attach_matrix_real(mat,localarray)
    type(StoreFullMatrixObjectReal) :: mat
    Real, DIMENSION(:,:),TARGET :: localarray

    CALL exit_if_factored(mat)
    CALL mp_attach_matrix_general_real(mat,localarray)
  END SUBROUTINE mp_attach_matrix_real

  SUBROUTINE mp_attach_matrix_general_real(mat,localarray)
    TYPE(StoreFullMatrixObjectReal) :: mat
    REAL,DIMENSION(*),TARGET :: localarray

    CALL exit_if_factored(mat)
    mat%localMatrix => localarray(1:mat%NBlocks*mat%Blocksize)
    mat%isAttached = .TRUE.

  END SUBROUTINE mp_attach_matrix_general_real

  SUBROUTINE mp_finalize_matrix_real(mat)
    type(StoreFullMatrixObjectReal) :: mat

    IF (mat%isAttached) THEN
       ! do nothing
       mat%isAttached=.FALSE.
       IF (ASSOCIATED(mat%localMatrix)) &
            & NULLIFY(mat%localMatrix)
    ELSE
       IF (ASSOCIATED(mat%localMatrix) ) DEALLOCATE(mat%localMatrix)
       IF (ASSOCIATED(mat%factoredBand)) DEALLOCATE(mat%factoredBand)
       IF (ASSOCIATED(mat%ipiv)        ) DEALLOCATE(mat%ipiv)
       IF (ASSOCIATED(mat%rhs)         ) DEALLOCATE(mat%rhs)
    END IF
    mat%isInitialized=.false.
  END SUBROUTINE mp_finalize_matrix_real

#ifdef OLD_PROCESSGRID
  SUBROUTINE mp_set_value_real(mat,global_row,global_col,value)
    TYPE(StoreFullMatrixObjectReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: global_row,global_col
    REAL, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc_r, proc_c, block_r, block_c, row_index, col_index
    INTEGER :: linear_index !, linear_blockindex

    DEBSTART("set_value(SFMO)")
    CALL exit_if_factored(mat)
    !PRINT*,"In routine set_value of StoreFullMatrixModule on PRow,PCol = ",mat%PRow,mat%PCol
    ! calculate the processor coordinates (starting with 0) of the given global indices
    proc_r = MOD((global_row-1)/mat%RowsPerBlock,mat%PG%NPRows)
    proc_c = MOD((global_col-1)/mat%ColsPerBlock,mat%PG%NPCols)
    !PRINT*,mat%Rank,mat%PRow,mat%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF ((proc_r.EQ.mat%PG%PRow) .AND. (proc_c.EQ.mat%PG%PCol)) THEN
       ! calculate the local block indices (starting from 0) in which the global indices lie
       block_r = (global_row-1)/(mat%PG%NPRows*mat%RowsPerBlock)
       block_c = (global_col-1)/(mat%PG%NPCols*mat%ColsPerBlock)
       ! as a last step, we calculate the indices (starting from 1) in the blocks
       row_index = MOD(global_row-1,mat%RowsPerBlock)+1
       col_index = MOD(global_col-1,mat%ColsPerBlock)+1
       !PRINT*,mat%Rank,"Setting (",row_index,",",col_index,") in block ",block_r,block_c
       !IF (mat%isVector) THEN
       !   mat%localVector(row_index) = value
       !ELSE
       !mat%localMatrix(row_index,col_index) = value
       ! calculate the linear index (starting from 1) in localMatrix
       !linear_blockindex = block_c*mat%BlockrowsPerProcess + block_r
       !linear_index = linear_blockindex*mat%Blocksize &
       !     &+ (col_index-1)*mat%RowsPerBlock+row_index
       linear_index = get_local_linear_index_real(mat,block_r,block_c,row_index,col_index)
       mat%localMatrix(linear_index) = value
       !PRINT*,"Global = (",global_row,global_col,"), P = (",proc_r,proc_c,"), Bl = (",block_r,block_c,"), ",&
        !    &row_index,col_index,", LI = ",linear_index, &
         !   &get_local_linear_index_real(mat,block_r,block_c,row_index,col_index)
       !END IF
    END IF
    DEBEND("set_value(SFMO)")
  END SUBROUTINE mp_set_value_real

  SUBROUTINE mp_add_value_real(mat,global_row,global_col,value)
    TYPE(StoreFullMatrixObjectReal) :: mat
    INTEGER, INTENT(IN) :: global_row,global_col
    REAL, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc_r, proc_c, block_r, block_c, row_index, col_index, linear_index

    DEBSTART("add_value(SFMO)")

    if (mat%isFactored) stop "Attempt to acces matrix entries of a factored matrix."
    !CALL exit_if_factored(mat)

    ! calculate the processor coordinates of the given global indices
    proc_r = MOD((global_row-1)/mat%RowsPerBlock,mat%PG%NPRows)
    proc_c = MOD((global_col-1)/mat%ColsPerBlock,mat%PG%NPCols)

    ! only process the entry if we are on the right processor
    IF ((proc_r.EQ.mat%PG%PRow) .AND. (proc_c.EQ.mat%PG%PCol)) THEN
       ! calculate the local block indices in which the global indices lie
       block_r = (global_row-1)/(mat%PG%NPRows*mat%RowsPerBlock)
       block_c = (global_col-1)/(mat%PG%NPCols*mat%ColsPerBlock)
       ! as a last step, we calculate the indices in the blocks
       row_index = MOD(global_row-1,mat%RowsPerBlock)+1
       col_index = MOD(global_col-1,mat%ColsPerBlock)+1

       !IF (mat%isVector) THEN
       !   mat%localVector(row_index) = mat%localVector(row_index) + value
       !ELSE
       !mat%localMatrix(row_index,col_index) = mat%localMatrix(row_index,col_index) + value
       linear_index = get_local_linear_index_real(mat,block_r,block_c,row_index,col_index)
       mat%localMatrix(linear_index) = mat%localMatrix(linear_index) + value
       !mat%localMatrix((col_index-1)*mat%RowsPerBlock+row_index) = &
       !     &mat%localMatrix((col_index-1)*mat%RowsPerBlock+row_index) + value
       !END IF
    END IF
    DEBEND("add_value(SFMO)")
    
  END SUBROUTINE mp_add_value_real

  !> Add the given array to the indicated row in the matrix
  !!
  !! This routine only works, if the processgrid is 1D and 
  !! PG\%NPCols=1. Only then a whole row of the matrix is local to one process.
  !! \param mat The StoreFullMatrixObject to work on
  !! \param global_row The global index of the row, starting with 1 up to mat%NRows
  !! \param whole_row An array containing the values, which will be added to the already
  !! existing values in row global_row of the matrix mat.
  subroutine sfm_add_to_row_real(mat,global_row,whole_row)
    TYPE(StoreFullMatrixObjectReal) :: mat
    INTEGER, INTENT(IN) :: global_row
    REAL, dimension(:),INTENT(IN) :: whole_row

    ! Local variables
    INTEGER :: proc_r, proc_c, row_index

    DEBSTART("add_value(SFMO)")

    if (mat%isFactored) stop "Attempt to acces matrix entries of a factored matrix."
    !CALL exit_if_factored(mat)

    if (mat%PG%NPCols.ne.1) stop "sfm_add_to_row can only be used if PG%NPCols.eq.1"

    ! calculate the processor coordinates of the given global indices
    proc_r = MOD((global_row-1)/mat%RowsPerBlock,mat%PG%NPRows)
    ! as per default PG%NPCols=1, proc_c is always 0, only if the
    ! 2D cyclic block distribution is really 2D, we have to reuse this line
    !proc_c = MOD((global_col-1)/mat%ColsPerBlock,mat%PG%NPCols)
    proc_c = 0

    ! only process the entry if we are on the right processor
    !IF ((proc_r.EQ.mat%PG%PRow) .AND. (proc_c.EQ.mat%PG%PCol)) THEN
    IF (proc_r.EQ.mat%PG%PRow) THEN
       ! calculate the local block indices in which the global indices lie
       !block_r = (global_row-1)/(mat%PG%NPRows*mat%RowsPerBlock)
       !block_c = (global_col-1)/(mat%PG%NPCols*mat%ColsPerBlock)
       ! as a last step, we calculate the indices in the blocks
       row_index = MOD(global_row-1,mat%RowsPerBlock)+1
       !col_index = MOD(global_col-1,mat%ColsPerBlock)+1

       !linear_index = get_local_linear_index(mat,block_r,block_c,row_index,col_index)
       call caxpy(mat%NCols,1,whole_row,1,mat%localMatrix(row_index),mat%RowsPerBlock)

       !mat%localMatrix(linear_index) = mat%localMatrix(linear_index) + value
    END IF
    DEBEND("add_value(SFMO)")
    
  END SUBROUTINE sfm_add_to_row_real

  FUNCTION mp_get_value_real(mat,irow,icol) RESULT(value)
    TYPE(StoreFullMatrixObjectReal) :: mat
    INTEGER :: irow, icol
    real :: value

    ! local variables
    INTEGER :: proc_r, proc_c, block_r, block_c, row_index, col_index, linear_index

    !CALL exit_if_factored(mat)
    proc_r = MOD((irow-1)/mat%RowsPerBlock,mat%PG%NPRows)
    proc_c = MOD((icol-1)/mat%ColsPerBlock,mat%PG%NPCols)
    !PRINT*,mat%Rank,mat%PRow,mat%PCol,proc_r, proc_c
    ! return the entry, if we are on the right processor
    IF ((proc_r.EQ.mat%PG%PRow) .AND. (proc_c.EQ.mat%PG%PCol)) THEN
       ! calculate the local block indices in which the global indices lie
       block_r = (irow-1)/(mat%PG%NPRows*mat%RowsPerBlock)
       block_c = (icol-1)/(mat%PG%NPCols*mat%ColsPerBlock)
       ! as a last step, we calculate the indices in the blocks
       row_index = MOD(irow-1,mat%RowsPerBlock)+1
       col_index = MOD(icol-1,mat%ColsPerBlock)+1
       !PRINT*,mat%Rank,"Setting (",row_index,",",col_index,") in block ",block_r,block_c
       linear_index = get_local_linear_index_real(mat,block_r,block_c,row_index,col_index)
       value = mat%localMatrix(linear_index)
    ELSE
       ! return nothing
       value = 1000
       ! this can be arbitrary, as it is never used in further calculations
       PRINT*,"This line should never be reached in storefullmatrix.F90."
       !CALL tracebackqq
    END IF
  END FUNCTION mp_get_value_real

#else
  SUBROUTINE mp_set_value_real(mat,irow,icol,value)
    TYPE(StoreFullMatrixObjectReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc, block_r, block_c, row_index, col_index
    INTEGER :: linear_index !, linear_blockindex

    DEBSTART("set_value(SFMO)")
    CALL exit_if_factored(mat)
    !PRINT*,"In routine set_value of StoreFullMatrixModule on PRow,PCol = ",mat%PRow,mat%PCol
    ! calculate the processor coordinates (starting with 0) of the given global indices
    IF (mat%isVector) THEN
       proc = MOD((irow-1)/mat%RowsPerBlock,mat%PG%NProcs)
    ELSE
       proc = MOD((icol-1)/mat%ColsPerBlock,mat%PG%NProcs)
    END IF
    !PRINT*,mat%Rank,mat%PRow,mat%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF (proc.EQ.mat%PG%Rank) THEN
       ! calculate the local block indices (starting from 0) in which the global indices lie
       IF (mat%isVector) THEN
          !blocknr = (irow-1)/(mat%PG%NProcs*mat%RowsPerBlock) is always zero
          row_index = MOD((irow-1),mat%RowsPerBlock)+1
          !PRINT*,mat%PG%Rank,": row_index = ", row_index,LBOUND(mat%localMatrix),UBOUND(mat%localMatrix)
          mat%localMatrix(row_index) = value
       ELSE
          block_r = (irow-1)/mat%RowsPerBlock
          block_c = (icol-1)/(mat%PG%NProcs*mat%ColsPerBlock)
          row_index = MOD(irow-1,mat%RowsPerBlock)+1
          col_index = MOD(icol-1,mat%ColsPerBlock)+1

          linear_index = get_local_linear_index_real(mat,block_r,block_c,row_index,col_index)
          mat%localMatrix(linear_index) = value
       END IF
    END IF
    DEBEND("set_value(SFMO)")
  END SUBROUTINE mp_set_value_real

  SUBROUTINE mp_add_value_real(mat,irow,icol,value)
    TYPE(StoreFullMatrixObjectReal) :: mat
    INTEGER, INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc, block_r, block_c, row_index, col_index, linear_index

    DEBSTART("add_value(SFMO)")

    CALL exit_if_factored(mat)
    ! calculate the processor coordinates of the given global indices
    IF (mat%isVector) THEN
       proc = MOD((irow-1)/mat%RowsPerBlock,mat%PG%NProcs)
    ELSE
       proc = MOD((icol-1)/mat%ColsPerBlock,mat%PG%NProcs)
    END IF
    !PRINT*,mat%Rank,mat%PRow,mat%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF (proc.EQ.mat%PG%Rank) THEN
       ! calculate the local block indices (starting from 0) in which the global indices lie
       IF (mat%isVector) THEN
          !blocknr = (irow-1)/(mat%PG%NProcs*mat%RowsPerBlock) is always zero
          row_index = MOD((irow-1),mat%RowsPerBlock)+1
          mat%localMatrix(row_index) = mat%localMatrix(row_index) + value
       ELSE
          block_r = (irow-1)/mat%RowsPerBlock
          block_c = (icol-1)/(mat%PG%NProcs*mat%ColsPerBlock)
          row_index = MOD(irow-1,mat%RowsPerBlock)+1
          col_index = MOD(icol-1,mat%ColsPerBlock)+1

          linear_index = get_local_linear_index_real(mat,block_r,block_c,row_index,col_index)
          mat%localMatrix(linear_index) = mat%localMatrix(linear_index) + value
       END IF
    END IF
    DEBEND("add_value(SFMO)")
    
  END SUBROUTINE mp_add_value_real

  FUNCTION mp_get_value_real(mat,irow,icol) RESULT(value)
    TYPE(StoreFullMatrixObjectReal) :: mat
    INTEGER :: irow, icol
    real :: value

    ! local variables
    INTEGER :: proc, block_r, block_c, row_index, col_index, linear_index

    CALL exit_if_factored(mat)
    IF (mat%isVector) THEN
       proc = MOD((irow-1)/mat%RowsPerBlock,mat%PG%NProcs)
    ELSE
       proc = MOD((icol-1)/mat%ColsPerBlock,mat%PG%NProcs)
    END IF
    !PRINT*,mat%Rank,mat%PRow,mat%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF (proc.EQ.mat%PG%Rank) THEN
       ! calculate the local block indices (starting from 0) in which the global indices lie
       IF (mat%isVector) THEN
          !blocknr = (irow-1)/(mat%PG%NProcs*mat%RowsPerBlock) is always zero
          row_index = MOD((irow-1),mat%RowsPerBlock)+1
          value = mat%localMatrix(row_index)
       ELSE
          block_r = (irow-1)/mat%RowsPerBlock
          block_c = (icol-1)/(mat%PG%NProcs*mat%ColsPerBlock)
          row_index = MOD(irow-1,mat%RowsPerBlock)+1
          col_index = MOD(icol-1,mat%ColsPerBlock)+1

          linear_index = get_local_linear_index_real(mat,block_r,block_c,row_index,col_index)
          value = mat%localMatrix(linear_index)
       END IF
    ELSE
       ! return nothing
       value = 1000
       ! this can be arbitrary, as it is never used in further calculations
    END IF
  END FUNCTION mp_get_value_real

#endif

  SUBROUTINE mp_commit_values_real(mat)
    TYPE(StoreFullMatrixObjectReal) :: mat

    ! nothing happens, as a Full matrix can not be further optimized

  END SUBROUTINE mp_commit_values_real

  !>Compute Matrix Matrix multiplication
  !! \param mat first matrix
  !! \param vec second matrix or vector
  !! \param res result
  !! \param transA_in optional flag for operator on 1st matrix ('N'=nothing,'T'=transpose, 'C'=transpose/conjugate)
  !! \param transB_in optional flag for operator on 2nd matrix ('N'=nothing,'T'=transpose, 'C'=transpose/conjugate)
  SUBROUTINE mp_dot_multiply_MatMat_real(mat, vec, res, transA_in, transB_in)
    TYPE(StoreFullMatrixObjectReal),INTENT(IN) :: mat
    TYPE(StoreFullMatrixObjectReal),INTENT(IN) :: vec
    TYPE(StoreFullMatrixObjectReal),INTENT(INOUT) :: res
    CHARACTER, INTENT(IN), OPTIONAL :: transA_in, transB_in

    !Local variables
    CHARACTER :: transA='N', transB='N'
    LOGICAL,PARAMETER :: original=.TRUE.

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
    IF ((vec%isVector).AND.(res%isVector)) THEN
       if (transB.ne.'N') then
          stop "transB ne 'N' is not implemented for vectors in storefullmatrix/mp_dot_multiply_MatMat"
       endif
       IF (original) THEN
          res%localMatrix = 0
          CALL psgemv(transA, mat%Nrows, mat%NCols, &
               & 1, mat%localMatrix, 1, 1, mat%Desc, &
               & vec%localMatrix, 1, 1, vec%Desc, 1, &
               & 0, res%localMatrix, 1, 1, res%Desc, 1)
       ELSE
          CALL psgemv(transA, mat%Nrows, 16, &
               & 1, mat%localMatrix, 1, 1, mat%Desc, &
               & vec%localMatrix, 1, 1, vec%Desc, 1, &
               & 0, res%localMatrix, 1, 1, res%Desc, 1)
          PRINT*,mat%PG%Rank,SUM(res%localMatrix*res%localMatrix)
       END IF
    ELSE
       CALL psgemm(transA, transB, mat%NRows, res%NCols, mat%NCols, &
            & 1, mat%localMatrix, 1, 1, mat%Desc, &
            & vec%localMatrix, 1, 1, vec%Desc,  &
            & 0, res%localMatrix, 1, 1, res%Desc)
    END IF

  END SUBROUTINE mp_dot_multiply_MatMat_real

  SUBROUTINE mp_dot_multiply_MatDagger_Mat_real(mat1, mat2, res)
    TYPE(StoreFullMatrixObjectReal),INTENT(IN) :: mat1, mat2
    TYPE(StoreFullMatrixObjectReal),INTENT(INOUT) :: res

    call mp_dot_multiply_MatMat_real(mat1,mat2,res,'C')

  END SUBROUTINE mp_dot_multiply_MatDagger_Mat_real


  SUBROUTINE mp_dot_multiply_MatVec_real(mat, vec, res, transA_in, transB_in)
    TYPE(StoreFullMatrixObjectReal),INTENT(IN) :: mat
    TYPE(StoreVectorObjectReal),INTENT(IN) :: vec
    TYPE(StoreVectorObjectReal),INTENT(INOUT) :: res
    CHARACTER, INTENT(IN), OPTIONAL :: transA_in, transB_in

    !Local variables
    CHARACTER :: transA='N', transB='N'
    LOGICAL,PARAMETER :: original=.TRUE.
    !COMPLEX, DIMENSION(res%NBlocks*res%Blocksize) :: temporary_array

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

    if (transB.ne.'N') then
       stop "transB ne 'N' is not implemented in storefullmatrix/mp_dot_multiply_MatVec"
    endif

    CALL exit_if_factored(mat)
  
    IF (original) THEN
       CALL psgemv(transA, mat%Nrows, mat%NCols, &
            & 1, mat%localMatrix, 1, 1, mat%Desc, &
            & vec%localVector, 1, 1, vec%Desc, 1, &
            & 0, res%localVector, 1, 1, res%Desc, 1)
    ELSE
       CALL psgemv(transA, mat%Nrows, 16, &
            & 1, mat%localMatrix, 1, 1, mat%Desc, &
            & vec%localVector, 1, 1, vec%Desc, 1, &
            & 0, res%localVector, 1, 1, res%Desc, 1)
       PRINT*,mat%PG%Rank,SUM(res%localVector*res%localVector)
    END IF
    
  END SUBROUTINE mp_dot_multiply_MatVec_real

  !*****************************************************************
  !* The global matrix, which is distributed over all processes
  !* is gathered and stored in a local array. This routine is
  !* mainly used for debugging and backward compatibility. In the
  !* final production code, it should not be called as it does not
  !* scale in memory.
  !*****************************************************************
  SUBROUTINE mp_get_global_matrix_locally_real(mat,localfullmat)
    TYPE(StoreFullMatrixObjectReal),INTENT(IN) :: mat
    REAL, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    ! local variables
    INTEGER :: ierr,iBlockcol, iBlockrow,iCol,iRow,linear_blockindex, linear_index
    REAL, DIMENSION(:,:),ALLOCATABLE :: localtransp,localtranspfullmat, localtemp
    !INTEGER :: desctemp(9)

    CALL exit_if_factored(mat)
    IF (mat%isVector) THEN
       ! If mat is a Vector, just concatenate the local vectors on all processes.
       CALL mpi_allgather(mat%localMatrix,mat%NBlocks*mat%Blocksize,MPI_REAL_TYPE,&
            &localfullmat,mat%NBlocks*mat%Blocksize,&
            &MPI_REAL_TYPE,mat%PG%Communicator,ierr)
    ELSE
       ! otherwise, we have to take into account the cyclic block matrix
       ! structure, which complicates the operation.
       ! First we allocate a matrix to take the transposed local matrix
       ALLOCATE(localtransp(mat%ColsPerBlock*mat%BlockcolsPerProcess,&
            &mat%RowsPerBlock*mat%BlockrowsPerProcess))
       ! and a global matrix to store the transposed global matrix.
       ALLOCATE(localtranspfullmat(mat%NCols,mat%NRows))

       ! additionally we need a local matrix untransposed to store the 
       ! correctly reordered local matrix
       ALLOCATE(localtemp(mat%RowsPerBlock*mat%BlockrowsPerProcess,&
            &mat%ColsPerBlock*mat%BlockcolsPerProcess))

       ! and we fill this local matrix
       DO iBlockcol = 0,mat%BlockcolsPerProcess-1
          DO iCol=1,mat%ColsPerBlock
             DO iBlockrow = 0,mat%BlockrowsPerProcess-1
                DO iRow=1,mat%RowsPerBlock
                   linear_blockindex = iBlockcol*mat%BlockrowsPerProcess + iBlockrow
                   linear_index = get_local_linear_index_real(mat,iBlockrow,iBlockcol,iRow,iCol)
                   !linear_blockindex*mat%Blocksize &
                   !     &+ (iCol-1)*mat%RowsPerBlock+iRow

                   localtemp(iBlockrow*mat%RowsPerBlock+iRow,&
                        &iBlockcol*mat%ColsPerBlock+iCol) = mat%localMatrix(linear_index)
                END DO
             END DO
          END DO
       END DO

       localtransp = transpose(localtemp)
       !CALL cgetmo(localtemp,mat%RowsPerBlock*mat%BlockrowsPerProcess,&
       !     &mat%RowsPerBlock*mat%BlockrowsPerProcess,&
       !     &mat%ColsPerBlock*mat%BlockcolsPerProcess,&
       !     &localtransp,mat%ColsPerBlock*mat%BlockcolsPerProcess)

       !PRINT*, "Communicator in MPI_Allgather is ", mat%PG%Communicator
       CALL mpi_allgather(localtransp,mat%NBlocks*mat%Blocksize,MPI_REAL_TYPE,&
            &localtranspfullmat,mat%NBlocks*mat%Blocksize,MPI_REAL_TYPE,&
            &mat%PG%Communicator,ierr)
       !CALL mpi_allgather(localtransp,mat%RowsPerBlock*mat%ColsPerBlock,MPI_REAL_TYPE,&
       !     &localtranspfullmat,mat%RowsPerBlock*mat%ColsPerBlock,MPI_REAL_TYPE,&
       !     &mat%PG%Communicator,ierr)
       
       !CALL cgetmo(localtranspfullmat,mat%NCols,mat%NCols,mat%NRows,&
       !     &localfullmat,mat%NRows)
       localfullmat = TRANSPOSE(localtranspfullmat)
       DEALLOCATE(localtransp)
       DEALLOCATE(localtranspfullmat)
       deallocate(localtemp)
    END IF
  END SUBROUTINE mp_get_global_matrix_locally_real

  FUNCTION mp_get_local_abs_square_sum_real(mat) RESULT(value)
    type(StoreFullMatrixObjectReal) :: mat
    REAL :: value

    CALL exit_if_factored(mat)
    value = SUM(REAL(mat%localMatrix)**2)

  END FUNCTION mp_get_local_abs_square_sum_real

  SUBROUTINE mp_show_on_screen_real(mat)
    type(StoreFullMatrixObjectReal) :: mat

    ! Local variables
    INTEGER :: iProc,iBlockrow,iBlockcol,iCol, iRow

    CALL exit_if_factored(mat)
    ! do some debugging output
    CALL blacs_barrier(mat%PG%Context,'A')
    DO iProc=0,mat%PG%NProcs-1
       CALL blacs_barrier(mat%PG%Context,'A')
       IF (iProc.EQ.mat%PG%Rank) THEN
          PRINT*,"---- Start output Rank = ",mat%PG%Rank," ----"
          DO iBlockrow=0,mat%BlockrowsPerProcess-1
             PRINT*,"  -- Start output Blockrow = ",iBlockrow," --"
             DO iRow=1,mat%RowsPerBlock
                DO iBlockcol=0,mat%BlockcolsPerProcess-1
                   DO iCol = 1,mat%ColsPerBlock
                      WRITE(*,"(2ES9.2,1X)",advance='NO') &
                           &mat%localMatrix(get_local_linear_index_real(mat,iBlockrow,iBlockcol,iRow,iCol))
                   END DO
                END DO
                WRITE(*,"(A)") ""
             END DO
             PRINT*,"  -- End   output Blockrow = ",iBlockrow," --"
          END DO
          PRINT*,"---- End   output Rank = ",mat%PG%Rank," ----"
       END IF
    END DO
    CALL blacs_barrier(mat%PG%Context,'A')
  END SUBROUTINE mp_show_on_screen_real

  SUBROUTINE mp_show_in_file_real(mat,filename)
    type(StoreFullMatrixObjectReal) :: mat
    CHARACTER(len=*) :: filename

    ! Local variables
    INTEGER :: iCol, iRow,iBlockrow,iBlockcol !,iProc
    character(len=FILENAME_MAX) :: full_filename

    CALL exit_if_factored(mat)
    ! open a file on each process
    WRITE(full_filename,"(A,I3.3,A)") filename,mat%PG%Rank,".dat"
    OPEN(77,file=TRIM(full_filename))

    ! do some debugging output
    CALL blacs_barrier(mat%PG%Context,'A')
    !DO iProc=0,mat%PG%NProcs-1
    !   CALL blacs_barrier(mat%PG%Context,'A')
    !   IF (iProc.EQ.mat%PG%Rank) THEN
          !PRINT*,"---- Start output Rank = ",mat%PG%Rank," ----"
          DO iBlockrow=0,mat%BlockrowsPerProcess-1
             !PRINT*,"  -- Start output Blockrow = ",iBlockrow," --"
             DO iRow=1,mat%RowsPerBlock
                !IF (mat%isVector) THEN
                   !WRITE(*,"(2ES10.2,X)") mat%localVector(iRow)
                !   WRITE(*,"(2ES10.2,X)") mat%localMatrix(iRow)
                !ELSE
                DO iBlockcol=0,mat%BlockcolsPerProcess-1
                   DO iCol = 1,mat%ColsPerBlock
                      WRITE(77,"(2ES20.10)",advance='NO') &
                           &mat%localMatrix(get_local_linear_index_real(mat,iBlockrow,iBlockcol,iRow,iCol))
                   END DO
!                   WRITE(77,"(A)",advance='NO') "|"
                END DO
                WRITE(77,"(A)") ""
                !END IF
             END DO
             !PRINT*,"  -- End   output Blockrow = ",iBlockrow," --"
          END DO
          !PRINT*,"---- End   output Rank = ",mat%PG%Rank," ----"
       !END IF
    !END DO
    CLOSE(77)
    CALL blacs_barrier(mat%PG%Context,'A')
  END SUBROUTINE mp_show_in_file_real

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix_real(this,mat)
    TYPE(StoreFullMatrixObjectReal), INTENT(INOUT) :: this
    TYPE(StoreFullMatrixObjectReal), INTENT(IN) :: mat

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)
    this%localMatrix = this%localMatrix + mat%localMatrix
  END SUBROUTINE mp_add_to_matrix_real

  SUBROUTINE mp_add_matrix_real(this,mat,res)
    TYPE(StoreFullMatrixObjectReal), INTENT(IN) :: this,mat
    TYPE(StoreFullMatrixObjectReal), INTENT(INOUT) :: res

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)
    CALL exit_if_factored(res)
    res%localMatrix = this%localMatrix+mat%localMatrix
  END SUBROUTINE mp_add_matrix_real

  SUBROUTINE mp_subtract_matrix_real(this,mat,res)
    TYPE(StoreFullMatrixObjectReal),INTENT(IN) :: this,mat
    TYPE(StoreFullMatrixObjectReal),intent(INOUT) :: res
    ! calculate res = this - mat

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)
    CALL exit_if_factored(res)
    res%localMatrix = this%localMatrix - mat%localMatrix
  END SUBROUTINE mp_subtract_matrix_real

  SUBROUTINE mp_subtract_from_matrix_real(this,mat)
    TYPE(StoreFullMatrixObjectReal),INTENT(INOUT) :: this
    TYPE(StoreFullMatrixObjectReal),intent(IN) :: mat
    ! calculate this = this - mat

    CALL exit_if_factored(this)
    CALL exit_if_factored(mat)
    this%localMatrix = this%localMatrix - mat%localMatrix
  END SUBROUTINE mp_subtract_from_matrix_real

  SUBROUTINE mp_multiply_matrix_with_real_real(this, scalar, res)
    TYPE(StoreFullMatrixObjectReal), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(StoreFullMatrixObjectReal),intent(INOUT) :: res

    CALL exit_if_factored(this)
    CALL exit_if_factored(res)
    res%localMatrix = this%localMatrix*scalar
  END SUBROUTINE mp_multiply_matrix_with_real_real

  SUBROUTINE mp_scale_matrix_by_real_real(this, scalar)
    TYPE(StoreFullMatrixObjectReal), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL exit_if_factored(this)
    this%localMatrix = this%localMatrix*scalar
  END SUBROUTINE mp_scale_matrix_by_real_real


  SUBROUTINE mp_assign_matrix_real(lmat,rmat)
    TYPE(StoreFullMatrixObjectReal), intent(INOUT) :: lmat
    TYPE(StoreFullMatrixObjectReal), intent(IN) :: rmat

    CALL exit_if_factored(lmat)
    CALL exit_if_factored(rmat)
    IF (lmat%PG.EQ.rmat%PG) THEN
       lmat%localMatrix = rmat%localMatrix
    ELSE
       PRINT*,"Assignment of matrices is only possible within the same context."
    END IF
  END SUBROUTINE mp_assign_matrix_real

  ! ===================================================
  ! == Define some mathematical routines on matrices ==
  ! ===================================================
  SUBROUTINE mp_invert_matrix_real(mat)
    type(StoreFullMatrixObjectReal) :: mat
    
    ! Local variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv, iwork
    REAL, DIMENSION(:), allocatable :: work
    INTEGER :: info,lwork,liwork,ub_locr

    CALL exit_if_factored(mat)
    ! calculate upper bound of ipiv storage needed
    ! LOCr(M_A)+MB_A (ScaLapack User's guide)
    ub_locr = ( (mat%NRows+mat%RowsPerBlock-1)/mat%RowsPerBlock &
         & + mat%PG%NPRows-1 ) / mat%PG%NPRows * mat%RowsPerBlock

    !PRINT*,"ub_locr = ",ub_locr
    ALLOCATE(ipiv(ub_locr+mat%RowsPerBlock))
    !WRITE(*,"(16ES10.2)") mat%localMatrix
    ! First do the LU decomposition in place
    CALL psgetrf(mat%NRows,mat%NCols,&
         &mat%localMatrix, 1, 1,&
         &mat%Desc,ipiv,info)

    IF (info.NE.0) THEN
       PRINT*,"pcgetrf produced an error in storefullmatrix/mp_invert_matrix"
       stop
    END IF
    !PRINT*,"after LU decomposition"
    !WRITE(*,"(16ES10.2)") mat%localMatrix
    

    ! Second calculate the inverse of the matrix
    ! determine the work array size
    ALLOCATE(work(1),iwork(1))
    CALL psgetri(mat%Nrows,mat%localMatrix,1,1,mat%Desc,ipiv,&
         &work,-1,iwork,-1,info)
    lwork = INT(work(1))
    liwork = iwork(1)
    !PRINT*,"lwork = ",lwork,", liwork = ",liwork
    DEALLOCATE(work,iwork)

    !PRINT*,"after first pcgetri call"
    !WRITE(*,"(16ES10.2)") mat%localMatrix

    
    ALLOCATE(work(lwork),iwork(liwork))
    CALL psgetri(mat%Nrows,mat%localMatrix,1,1,mat%Desc,ipiv,&
         &work,lwork,iwork,liwork,info)

    IF (info.NE.0) THEN
       PRINT*,"pcgetri produced an error in storefullmatrix/mp_invert_matrix"
       stop
    END IF
    !PRINT*,"after second pcgetri call"
    !WRITE(*,"(16ES10.2)") mat%localMatrix

    DEALLOCATE(work,iwork)
    DEALLOCATE(ipiv)
  END SUBROUTINE mp_invert_matrix_real

  SUBROUTINE mp_LU_factor_real(mat)
    type(StoreFullMatrixObjectReal) :: mat
    
    ! Local variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
    REAL, DIMENSION(:,:), ALLOCATABLE :: ab
    INTEGER :: i, j, j_off, nbw(2), ierr, nstore

    CALL exit_if_factored(mat)

    ! This routine only works on a grid with distributed rows

    IF(mat%PG%NPCols /= 1 .OR. mat%PG%NPRows /= mat%PG%NProcs) STOP 'Illegal grid for LU_factor'

    ! Get bandwidth

    j_off = mat%PG%PRow * mat%RowsPerBlock

    nbw(1) = 0 ! lower bandwidth
    nbw(2) = 0 ! upper bandwidth
    DO i = 1, mat%NCols
    DO j = 1, mat%RowsPerBlock
       IF(mat%localMatrix((i-1)*mat%RowsPerBlock+j) /= 0.) THEN
          IF(j+j_off > i) THEN
             nbw(1) = MAX(nbw(1),(j+j_off)-i)
          ELSE
             nbw(2) = MAX(nbw(2),i-(j+j_off))
          ENDIF
       ENDIF
    ENDDO
    ENDDO

    CALL MPI_Allreduce(MPI_IN_PLACE, nbw, 2, MPI_INTEGER, MPI_MAX, mat%PG%Communicator, ierr)

    mat%bwl = nbw(1)
    mat%bwu = nbw(2)

    ! Get storage size of LU decomposition (neglecting ipiv)

    nstore = (2*mat%bwl+mat%bwu+1)*mat%RowsPerBlock

    ! The following is just for testing if inversion/LU have the same results
    ! Setting the environment variable GENE_FORCE_INVERT to some value
    ! will force inversion instead of LU decomposition.
    call GET_ENVIRONMENT_VARIABLE('GENE_FORCE_INVERT',status=ierr)
    if(ierr==0) nstore = 1000000000

    ! Check if it pays off to factor the matrix, otherways just use mp_invert_matrix

    IF(mat%NCols*mat%RowsPerBlock <= nstore) THEN

       CALL mp_invert_matrix_real(mat)

       mat%bwl = -1
       mat%bwu = -1
       mat%isFactored = .TRUE.
       RETURN
    ENDIF

    ! Get the data of the matrix bands on all procs.
    ! Doing it this way (using MPI_Allreduce) is not the most efficient way,
    ! but this routine is called only in the initialization phase.

    ALLOCATE(ab(2*mat%bwl+mat%bwu+1,mat%NCols))
    ab(:,:) = 0

    DO i = 1, mat%NCols
    DO j = 1, mat%RowsPerBlock
       IF((j+j_off)-i <= mat%bwl .AND. i-(j+j_off) <= mat%bwu) THEN
          ab(mat%bwl+mat%bwu+1+(j+j_off)-i,i) = mat%localMatrix((i-1)*mat%RowsPerBlock+j)
       ENDIF
    ENDDO
    ENDDO

    CALL MPI_Allreduce(MPI_IN_PLACE,ab,SIZE(ab),MPI_REAL_TYPE,MPI_SUM,mat%PG%Communicator, ierr)

    ! Use LAPACK routine cgbtrf (on all processors with identical data)
    ! and store the results distributed.
    ! Again, this is not the most efficient way but very simple to do ...

    ALLOCATE(ipiv(mat%NCols))

    CALL sgbtrf(mat%NRows,mat%NCols,mat%bwl,mat%bwu,ab,ubound(ab,1),ipiv,ierr)

    IF (ierr.NE.0) THEN
       PRINT*,"cgbtrf info = ",ierr
       stop
    END IF

    allocate(mat%factoredBand(2*mat%bwl+mat%bwu+1,mat%RowsPerBlock))
    allocate(mat%ipiv(mat%RowsPerBlock))

    mat%factoredBand(:,1:mat%RowsPerBlock) = ab(:,j_off+1:j_off+mat%RowsPerBlock)
    mat%ipiv(1:mat%RowsPerBlock) = ipiv(j_off+1:j_off+mat%RowsPerBlock)

    DEALLOCATE(ipiv)
    DEALLOCATE(ab)

    ! mat%localMatrix is not needed any more
    DEALLOCATE(mat%localMatrix)

    mat%isFactored = .TRUE.

  END SUBROUTINE mp_LU_factor_real

  SUBROUTINE mp_LU_solve_real(mat, vec, res, phase)
    TYPE(StoreFullMatrixObjectReal),INTENT(INOUT) :: mat
    TYPE(StoreFullMatrixObjectReal),INTENT(IN) :: vec
    TYPE(StoreFullMatrixObjectReal),INTENT(INOUT) :: res
    INTEGER, INTENT(IN) :: phase

    INTEGER :: j, l, lm, kd, j_off, j_max, j_min, ierr
    INTEGER :: mpi_status(mpi_status_size)
    REAL :: ctmp

    ! This routine must be called 3 times (with phase=1,2,3) in order
    ! to solve L*U*x = b completly.
    ! For phase=1, a copy of vec is stored internally, res is not used
    ! For phase=2, the forward substitution is done, vec and res are not used
    ! For phase=3, the backward substitution is done, res is set to the solution, vec is not used
    !
    ! The reason for doing it this complicated way is that a single LU-solution
    ! does not parallelize since it is sweeping once forward and back through
    ! all processors.
    !
    ! For multiple LU solutions, however, the processors may overlap in doing work
    ! (the first one working on solution n, the next one on solution n-1 etc.)
    ! if first all forward substitutions are done and then all backward substitutions.

    ! Safety checks:

    IF(vec%NCols /= 1) stop 'LU_solve: NCols must be 1'
    IF(.NOT. mat%isFactored) stop 'LU_solve: not factored'
    IF(mat%PG%NPCols /= 1 .OR. mat%PG%NPRows /= mat%PG%NProcs) STOP 'Illegal grid for LU_solve'

    ! If it didn't pay off to use the LU factorization of the banded matrix,
    ! the inverse of the full matrix is stored and mat%bwl is < 0

    IF(mat%bwl<0) THEN

       IF(phase == 1) THEN
          ! Allocate mat%rhs and store the complete rhs on every proc
          ALLOCATE(mat%rhs(mat%NRows))
          CALL MPI_Allgather(vec%localMatrix,mat%RowsPerBlock,MPI_REAL_TYPE,&
                             mat%rhs,        mat%RowsPerBlock,MPI_REAL_TYPE,&
                             vec%PG%Communicator,ierr)
       ENDIF
       ! Nothing to do for phase == 2
       IF(phase == 3) THEN
          ! Multiply inverse matrix with stored vector
          CALL cgemv('N', mat%RowsPerBlock, mat%NCols, &
                     1, mat%localMatrix, mat%RowsPerBlock, &
                     mat%rhs, 1, 0, res%localMatrix, 1)
          DEALLOCATE(mat%rhs)
       ENDIF
       RETURN

    ENDIF

    j_off = mat%PG%PRow * mat%RowsPerBlock

    KD = mat%bwl+mat%bwu+1 ! Position of diagonal in mat%factoredBand

    IF(phase == 1) THEN

       ! Store a copy of vec (on proc 0)

       ALLOCATE(mat%rhs(mat%NRows))

       CALL MPI_Gather(vec%localMatrix,mat%RowsPerBlock,MPI_REAL_TYPE,&
                       mat%rhs,        mat%RowsPerBlock,MPI_REAL_TYPE,&
                       0,vec%PG%Communicator,ierr)
    ENDIF

    IF(phase == 2) THEN

       ! Forward substitution

       IF(mat%PG%PRow > 0) THEN
          ! Wait until previous proc sends mat%rhs
          j_min = MAX(j_off + 1 - (mat%bwl+mat%bwu), 1)
          j_max = mat%NRows
          CALL MPI_Recv(mat%rhs(j_min),j_max-j_min+1,MPI_REAL_TYPE,mat%PG%PRow-1,111, &
                        vec%PG%Communicator, mpi_status, ierr)
       ENDIF

       ! My part of forward substitution

       DO j = j_off+1, j_off+mat%RowsPerBlock
          if(j==mat%NRows) exit
          lm = MIN( mat%bwl, mat%NRows-j )
          l = mat%ipiv(j-j_off)
          IF( l.NE.j ) then
             ctmp = mat%rhs(l)
             mat%rhs(l) = mat%rhs(j)
             mat%rhs(j) = ctmp
          endif
          mat%rhs(j+1:j+lm) = mat%rhs(j+1:j+lm) - mat%rhs(j)*mat%factoredBand(KD+1:KD+lm,j-j_off)
       enddo

       IF(mat%PG%PRow < mat%PG%NPRows-1) THEN
          ! Send mat%rhs to next proc
          j_min = MAX(j_off + mat%RowsPerBlock + 1 - (mat%bwl+mat%bwu), 1)
          j_max = mat%NRows
          CALL MPI_Send(mat%rhs(j_min),j_max-j_min+1,MPI_REAL_TYPE,mat%PG%PRow+1,111, &
                        vec%PG%Communicator, ierr)
       ENDIF
    ENDIF

    IF(phase == 3) THEN

       ! Backward substitution

       IF(mat%PG%PRow < mat%PG%NPRows-1) THEN
          ! Wait until next proc sends back mat%rhs
          j_min = MAX(j_off + mat%RowsPerBlock + 1 - (mat%bwl+mat%bwu), 1)
          j_max = j_off + mat%RowsPerBlock
          CALL MPI_Recv(mat%rhs(j_min),j_max-j_min+1,MPI_REAL_TYPE,mat%PG%PRow+1,222, &
                        vec%PG%Communicator, mpi_status, ierr)
       ENDIF

       ! My part of backward substitution

       DO j = j_off+mat%RowsPerBlock, j_off+1, -1
          IF (mat%rhs(j).NE.0.) THEN
             mat%rhs(j) = mat%rhs(j)/mat%factoredBand(KD,j-j_off)
             if(j==1) exit
             lm = MIN(mat%bwl+mat%bwu, j-1)
             mat%rhs(j-lm:j-1) = mat%rhs(j-lm:j-1) - mat%rhs(j)*mat%factoredBand(KD-lm:KD-1,j-j_off)
          END IF
       ENDDO

       IF(mat%PG%PRow > 0) THEN
          ! Send mat%rhs to previous proc
          j_min = MAX(j_off + 1 - (mat%bwl+mat%bwu), 1)
          j_max = j_off
          CALL MPI_Send(mat%rhs(j_min),j_max-j_min+1,MPI_REAL_TYPE,mat%PG%PRow-1,222, &
                        vec%PG%Communicator, ierr)
       ENDIF

       ! Store local part of mat%rhs

       res%localMatrix(1:mat%RowsPerBlock) = mat%rhs(j_off+1:j_off+mat%RowsPerBlock)

       DEALLOCATE(mat%rhs)

    ENDIF

  END SUBROUTINE mp_LU_solve_real

  SUBROUTINE sfm_exit_if_factored_real(mat)
    type(StoreFullMatrixObjectReal) :: mat

    INTEGER :: ierr

    ! This routines is called by all routines which use mat%localMatrix.
    ! It does exit the program if the matrix is LU-factored.

    IF(mat%isFactored) THEN
      PRINT *,'ERROR: Attempt use original matrix data of a LU-factored Matrix'
      CALL mpi_abort(MPI_COMM_WORLD, 1, ierr)
    ENDIF

  END SUBROUTINE sfm_exit_if_factored_real

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_matrix_sum_real(mat,communicator)
    TYPE(StoreFullMatrixObjectReal), INTENT(INOUT) :: mat
    REAL, DIMENSION(mat%NBlocks*mat%Blocksize) :: totalsum
    integer :: communicator

    integer :: ierr
    CALL mpi_allreduce(mat%localMatrix,totalsum,mat%NBlocks*mat%Blocksize, &
         & MPI_REAL_TYPE, MPI_SUM, communicator, ierr)
    
    mat%localMatrix = totalsum
  END SUBROUTINE my_matrix_sum_real

  ! ==================================
  ! == only private helper routines ==
  ! ==================================
  FUNCTION get_local_linear_index_real(mat,blockrow,blockcol,irow,icol) RESULT(index)
    TYPE(StoreFullMatrixObjectReal), intent(IN) :: mat
    INTEGER,intent(IN) :: blockrow,blockcol,irow,icol
    integer :: index
    
    ! calculate the linear index of the element with coordinates irow,icol in the block
    ! with the coordinates blockrow,blockcol of matrix mat
    !linear_blockindex = blockcol*mat%BlockrowsPerProcess + blockrow
    index = (blockcol*mat%BlockrowsPerProcess + blockrow)*mat%Blocksize &
         &+ (icol-1)*mat%RowsPerBlock+irow
  END FUNCTION get_local_linear_index_real

  subroutine sfm_print_storage_details_real(mat)
    type(StoreFullMatrixObjectReal) :: mat

    integer :: rank,ierr

    if (.not.mat%isInitialized) then
       print*,"Matrix is NOT initialized."
       return
    end if

    call mpi_comm_rank(mat%PG%communicator,rank,ierr)
    if (rank.eq.0) then
       print*,"NRows,NCols = ",mat%NRows,mat%NCols
       print*,"RowsPerBlock,ColsPerBlock = ",mat%RowsPerBlock,mat%ColsPerBlock
       print*,"NBlocks,Blocksize = ", mat%NBlocks, mat%Blocksize
       if (mat%NBlocks.gt.1) then
          print*,"BlockrowsPerProcess,BlockcolsPerProcess",mat%BlockrowsPerProcess,mat%BlockcolsPerProcess
       end if
       
       if (mat%isAttached) then
          print*,"Matrix is attached."
       else
          print*,"Matrix is NOT attached."
       end if
    end if
  end subroutine sfm_print_storage_details_real

END MODULE StoreFullMatrixModule



