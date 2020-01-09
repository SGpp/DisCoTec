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

MODULE MatrixModule
  USE ProcessGridModule
  USE StoreFullMatrixModule
  !USE mpi
  IMPLICIT NONE
  ! setting default accessibility to private
  private

  ! explicitly declare the interfaces as public
  PUBLIC :: initialize, initialize_matrix_module, finalize, finalize_matrix_module
  PUBLIC :: allocate, attach, set_value, add_value, add_to_row,commit_values, set_zero
  PUBLIC :: get_global_matrix_locally, get_local_abs_square_sum, mat_get_value
  PUBLIC :: dot_multiply, show, add_matrix, subtract_matrix, multiply_matrix_with_scalar
  PUBLIC :: dot_multiply_MatDagger_Mat
  PUBLIC :: invert, output_Data, isInitialized, my_sum, ASSIGNMENT(=)
  PUBLIC :: static_size_of_matrix, get_processgrid, LU_factor, LU_solve
  ! the following variables are used from ProcessGridModule. To use them outside the 
  ! MatrixModule, we have to publicize them explicitly.
  PUBLIC :: MAT_BLOCKS_OF_ROWS, MAT_BLOCKS_OF_COLS, MAT_BOTH

#include "matrix.h"

  ! Use the same process grid for all Matrices instantiated whithin this module
  TYPE(ProcessGrid),SAVE :: PG

CONTAINS
  FUNCTION static_size_of_matrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = 2*SIZE_OF_INTEGER &
         & + static_size_of_storefullmatrixobject()
    
  END FUNCTION static_size_of_matrix
    
  SUBROUTINE mp_initialize_matrix(mat,n_rows,n_cols)
    TYPE(Matrix), INTENT(INOUT) :: mat
    INTEGER :: n_rows, n_cols

    DEBSTART("mp_initialize_matrix")
    IF (PG%isInitialized) THEN
       mat%NRows = n_rows
       mat%NCols = n_cols
       CALL initialize(mat%Data,n_rows,n_cols, PG)
    ELSE
       PRINT*,"Before initializing a matrix, the process grid has to be initialized."
       STOP
    END IF
    DEBEND("mp_initialize_matrix")
  END SUBROUTINE mp_initialize_matrix

  SUBROUTINE get_processgrid(ao_PG)
    type(ProcessGrid) :: ao_PG
    ao_PG = PG
  END SUBROUTINE get_processgrid

  SUBROUTINE mp_initialize_vector(vec,n_rows)
    type(Matrix) :: vec
    integer :: n_rows

    CALL initialize(vec,n_rows,1)
  END SUBROUTINE mp_initialize_vector

  SUBROUTINE mp_initialize_matrix_module(comm)
    INTEGER :: comm

    DEBSTART("mp_initialize_matrix_module")
    CALL initialize(PG,MAT_BLOCKS_OF_ROWS, comm)
    !CALL initialize(PG_row,MAT_BLOCKS_OF_ROWS, comm)
    DEBEND("mp_initialize_matrix_module")
  END SUBROUTINE mp_initialize_matrix_module

  LOGICAL FUNCTION mp_isInitialized(mat) 
    type(Matrix) :: mat

    mp_isInitialized = isInitialized(mat%Data)
  END FUNCTION mp_isInitialized

  SUBROUTINE mp_allocate(mat)
    TYPE(Matrix), intent(INOUT) :: mat

    CALL allocate(mat%Data)
  END SUBROUTINE mp_allocate

  SUBROUTINE mp_attach_vector(mat,localarray)
    TYPE(Matrix),INTENT(INOUT) :: mat
    COMPLEX, DIMENSION(:),INTENT(IN),TARGET :: localarray

    DEBSTART("mp_attach_vector")
    CALL attach(mat%Data,localarray)
    DEBEND("mp_attach_vector")
  END SUBROUTINE mp_attach_vector

  SUBROUTINE mp_attach_matrix(mat,localarray)
    TYPE(Matrix),INTENT(INOUT) :: mat
    COMPLEX, DIMENSION(:,:),INTENT(IN),TARGET :: localarray

    CALL attach(mat%Data,localarray)
  END SUBROUTINE mp_attach_matrix

  SUBROUTINE mp_finalize_matrix(mat)
    type(Matrix) :: mat
    call finalize(mat%Data)
  END SUBROUTINE mp_finalize_matrix

  SUBROUTINE mp_finalize_module()
    call finalize(PG)
  END SUBROUTINE mp_finalize_module

  SUBROUTINE mp_set_complex_value(mat,irow,icol,value)
    TYPE(Matrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_set_complex_value

  SUBROUTINE mp_set_real_value(mat,irow,icol,value)
    TYPE(Matrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(mat%Data,irow,icol,CMPLX(value,0.0))
  END SUBROUTINE mp_set_real_value

  SUBROUTINE mp_add_complex_value(mat,irow,icol,value)
    TYPE(Matrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    CALL add_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_add_complex_value

  SUBROUTINE mp_add_real_value(mat,irow,icol,value)
    TYPE(Matrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    CALL add_value(mat%Data,irow,icol,CMPLX(value,0.0))
  END SUBROUTINE mp_add_real_value

  SUBROUTINE mat_add_to_row(mat,irow,whole_row)
    TYPE(Matrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow
    COMPLEX, dimension(:),INTENT(IN) :: whole_row

    CALL add_to_row(mat%Data,irow,whole_row)
  END SUBROUTINE mat_add_to_row


  SUBROUTINE mp_commit_values(mat)
    TYPE(Matrix) :: mat
    
    CALL commit_values(mat%Data)
  END SUBROUTINE mp_commit_values

  SUBROUTINE mp_dot_multiply_MatMat(this,mat, resmat, transA_in, transB_in)
    TYPE(Matrix),intent(IN) :: this
    TYPE(Matrix),intent(IN) :: mat
    TYPE(Matrix),intent(INOUT) :: resmat
    CHARACTER, intent(IN), optional :: transA_in, transB_in

    !PRINT*,mat%NCols,mat%NRows,resmat%NCols,resmat%NRows
    IF ((this%NCols.EQ.mat%NRows).AND. &
         &(this%NRows.EQ.resmat%NRows).AND.&
         &(mat%NCols.EQ.resmat%NCols)) THEN
       if (.not.(present(transA_in).or.present(transB_in))) then
          CALL dot_multiply(this%Data,mat%Data,resmat%Data)
       elseif (present(transA_in).and.present(transB_in)) then
          CALL dot_multiply(this%Data,mat%Data,resmat%Data,transA_in,transB_in)
       elseif (present(transA_in)) then
          CALL dot_multiply(this%Data,mat%Data,resmat%Data,transA_in)
       endif
    ELSE
       PRINT*,"Matrix shapes for multiplication with Matrix do not match, ",&
            &this%NCols,this%NRows,mat%NCols,mat%NRows,resmat%NCols,resmat%NRows
       stop
    END IF
  END SUBROUTINE mp_dot_multiply_MatMat

  SUBROUTINE mp_set_zero(mat)
    type(Matrix) :: mat

    ! Local variables
    INTEGER :: icol, irow

    DEBSTART("mp_set_zero")
    DO irow=1,mat%NRows
       DO icol=1,mat%NCols
          CALL set_value(mat,irow,icol,CMPLX(0.0,0.0))
       END DO
    END DO
    CALL commit_values(mat)
    DEBEND("mp_set_zero")
  END SUBROUTINE mp_set_zero

  SUBROUTINE mp_get_global_matrix_locally(mat,localfullmat)
    TYPE(Matrix),INTENT(IN) :: mat
    COMPLEX, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    CALL get_global_matrix_locally(mat%Data,localfullmat)
  END SUBROUTINE mp_get_global_matrix_locally

  FUNCTION mp_get_value(mat,irow,icol) result(value)
    type(Matrix) :: mat
    INTEGER :: irow, icol
    complex :: value

    value = mat_get_value(mat%Data,irow,icol)
  END FUNCTION mp_get_value

  FUNCTION mp_get_local_abs_square_sum(mat) RESULT(value)
    type(Matrix) :: mat
    real :: value

    value = get_local_abs_square_sum(mat%Data)
  END FUNCTION mp_get_local_abs_square_sum

  SUBROUTINE mp_show_on_screen(mat)
    type(Matrix) :: mat
    
    CALL show(mat%Data)
  END SUBROUTINE mp_show_on_screen

  SUBROUTINE mp_show_in_file(mat,filename)
    type(Matrix) :: mat
    CHARACTER(len=*) :: filename

    CALL show(mat%Data,filename)
  END SUBROUTINE mp_show_in_file

  SUBROUTINE mp_output_data_matrix(basename,mat)
    type(Matrix) :: mat
    CHARACTER(len=*) :: basename
    
    character(len=100) :: filename

    WRITE(filename,"(3A)") "./",TRIM(basename),".dat"
    PRINT*,"Writing to file ",TRIM(filename)
    CALL show(mat,filename)
  END SUBROUTINE mp_output_data_matrix

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix(this,mat)
    TYPE(Matrix), intent(INOUT) :: this
    TYPE(Matrix),INTENT(IN) :: mat
    CALL add_matrix(this%Data,mat%Data)
  END SUBROUTINE mp_add_to_matrix

  SUBROUTINE mp_add_matrix(this,mat,res)
    TYPE(Matrix), INTENT(IN) :: this,mat
    TYPE(Matrix), INTENT(INOUT) :: res
    !PRINT*,"in mp_add_matrix of MatrixModule"
    CALL add_matrix(this%Data,mat%Data,res%Data)
  END SUBROUTINE mp_add_matrix

  SUBROUTINE mp_subtract_matrix(this,mat,res)
    TYPE(Matrix),INTENT(IN) :: this,mat
    TYPE(Matrix),intent(INOUT) :: res
    ! calculate res = this - mat

    CALL subtract_matrix(this%Data,mat%Data, res%Data)
  END SUBROUTINE mp_subtract_matrix

  SUBROUTINE mp_subtract_from_matrix(this,mat)
    TYPE(Matrix),INTENT(INOUT) :: this
    TYPE(Matrix), intent(IN) :: mat
    ! calculate this = this - mat

    CALL subtract_matrix(this%Data,mat%Data)
  END SUBROUTINE mp_subtract_from_matrix
    
  SUBROUTINE mp_multiply_matrix_with_real(this, scalar, res)
    TYPE(Matrix), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(Matrix),intent(INOUT) :: res

    CALL multiply_matrix_with_scalar(this%Data,scalar, res%Data)
  END SUBROUTINE mp_multiply_matrix_with_real

  SUBROUTINE mp_scale_matrix_by_real(this, scalar)
    TYPE(Matrix), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL multiply_matrix_with_scalar(this%Data,scalar)
  END SUBROUTINE mp_scale_matrix_by_real

  SUBROUTINE mp_assign_matrix(lmat,rmat)
    TYPE(Matrix), intent(INOUT) :: lmat
    TYPE(Matrix), intent(IN) :: rmat

    IF ((lmat%Nrows.EQ.rmat%Nrows).AND.(lmat%Ncols.EQ.rmat%Ncols)) THEN
       lmat%Data = rmat%Data
    ELSE
       PRINT*, "Left and right matrices in an assignment must have the same dimensions."
    END IF
  END SUBROUTINE mp_assign_matrix

  ! ====================================================
  ! == Define some mathematical routines for matrices ==
  ! ====================================================
  SUBROUTINE mp_invert_matrix(mat)
    TYPE(Matrix) :: mat

    CALL invert(mat%Data)
  END SUBROUTINE mp_invert_matrix

  SUBROUTINE mp_LU_factor(mat)
    TYPE(Matrix) :: mat

    IF (mat%NRows.EQ.mat%NCOLS) THEN
      CALL LU_factor(mat%Data)
    ELSE
      PRINT*,"Matrix for _LU_factor is not square"
      STOP
    ENDIF
  END SUBROUTINE mp_LU_factor

  SUBROUTINE mp_LU_solve(this, mat, resmat, phase)
    TYPE(Matrix),intent(INOUT) :: this
    TYPE(Matrix),intent(IN) :: mat
    TYPE(Matrix),intent(INOUT) :: resmat
    INTEGER,intent(IN) :: phase

    IF ( this%NRows.EQ.this%NCols .AND. &
         this%NRows.EQ.mat%NRows .AND.&
         this%NRows.EQ.resmat%NRows .AND.&
         mat%NCols.EQ.resmat%NCols ) THEN
       CALL LU_solve(this%Data,mat%Data,resmat%Data,phase)
    ELSE
       PRINT*,"Matrix shapes for LU_solve do not match, ",&
            &this%NCols,this%NRows,mat%NCols,mat%NRows,resmat%NCols,resmat%NRows
       STOP
    END IF
  END SUBROUTINE mp_LU_solve

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_matrix_sum_0D(localsum,communicator)
    !INTEGER, intent(IN) :: len
    TYPE(Matrix), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_sum(localsum%Data,communicator)
  END SUBROUTINE my_matrix_sum_0D

  SUBROUTINE my_matrix_sum_2D(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(Matrix), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_matrix_sum_generic(localsum,len,communicator)
  END SUBROUTINE my_matrix_sum_2D

  SUBROUTINE my_matrix_sum_3D(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(Matrix), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_matrix_sum_generic(localsum,len,communicator)
  END SUBROUTINE my_matrix_sum_3D
  
  SUBROUTINE my_matrix_sum_generic(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(Matrix), intent(INOUT) :: localsum(len)
    integer :: communicator

    integer :: iMatrix
    DO iMatrix=1,len
       CALL my_sum(localsum(iMatrix)%Data,communicator)
    END DO
  END SUBROUTINE my_matrix_sum_generic

  ! ====================================================
  ! ==    Everything repeated, for REAL data types    ==
  ! ====================================================

  SUBROUTINE mp_initialize_matrix_real(mat,n_rows,n_cols)
    TYPE(MatrixReal), INTENT(INOUT) :: mat
    INTEGER :: n_rows, n_cols

    DEBSTART("mp_initialize_matrix_real")
    IF (PG%isInitialized) THEN
       mat%NRows = n_rows
       mat%NCols = n_cols
       CALL initialize(mat%Data,n_rows,n_cols, PG)
    ELSE
       PRINT*,"Before initializing a matrix, the process grid has to be initialized."
       STOP
    END IF
    DEBEND("mp_initialize_matrix_real")
  END SUBROUTINE mp_initialize_matrix_real

  SUBROUTINE mp_initialize_vector_real(vec,n_rows)
    type(MatrixReal) :: vec
    integer :: n_rows

    CALL initialize(vec%Data,n_rows,1,PG)
  END SUBROUTINE mp_initialize_vector_real

  LOGICAL FUNCTION mp_isInitialized_real(mat) 
    type(MatrixReal) :: mat

    mp_isInitialized_real = isInitialized(mat%Data)
  END FUNCTION mp_isInitialized_real

  SUBROUTINE mp_allocate_real(mat)
    TYPE(MatrixReal), intent(INOUT) :: mat

    CALL allocate(mat%Data)
  END SUBROUTINE mp_allocate_real

  SUBROUTINE mp_attach_vector_real(mat,localarray)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    REAL, DIMENSION(:),INTENT(IN),TARGET :: localarray

    DEBSTART("mp_attach_vector_real")
    CALL attach(mat%Data,localarray)
    DEBEND("mp_attach_vector_real")
  END SUBROUTINE mp_attach_vector_real

  SUBROUTINE mp_attach_matrix_real(mat,localarray)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    REAL, DIMENSION(:,:),INTENT(IN),TARGET :: localarray

    CALL attach(mat%Data,localarray)
  END SUBROUTINE mp_attach_matrix_real

  SUBROUTINE mp_finalize_matrix_real(mat)
    type(MatrixReal) :: mat
    call finalize(mat%Data)
  END SUBROUTINE mp_finalize_matrix_real

  SUBROUTINE mp_set_complex_value_real(mat,irow,icol,value)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_set_complex_value_real

  SUBROUTINE mp_set_real_value_real(mat,irow,icol,value)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_set_real_value_real

  SUBROUTINE mp_add_complex_value_real(mat,irow,icol,value)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    CALL add_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_add_complex_value_real

  SUBROUTINE mp_add_real_value_real(mat,irow,icol,value)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    CALL add_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_add_real_value_real

  SUBROUTINE mat_add_to_row_real(mat,irow,whole_row)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow
    REAL, dimension(:),INTENT(IN) :: whole_row

    CALL add_to_row(mat%Data,irow,whole_row)
  END SUBROUTINE mat_add_to_row_real


  SUBROUTINE mp_commit_values_real(mat)
    TYPE(MatrixReal) :: mat
    
    CALL commit_values(mat%Data)
  END SUBROUTINE mp_commit_values_real

  SUBROUTINE mp_dot_multiply_MatMat_real(this,mat, resmat, transA_in, transB_in)
    TYPE(MatrixReal),intent(IN) :: this
    TYPE(MatrixReal),intent(IN) :: mat
    TYPE(MatrixReal),intent(INOUT) :: resmat
    CHARACTER, intent(IN), optional :: transA_in, transB_in

    !PRINT*,mat%NCols,mat%NRows,resmat%NCols,resmat%NRows
    IF ((this%NCols.EQ.mat%NRows).AND. &
         &(this%NRows.EQ.resmat%NRows).AND.&
         &(mat%NCols.EQ.resmat%NCols)) THEN
       if (.not.(present(transA_in).or.present(transB_in))) then
          CALL dot_multiply(this%Data,mat%Data,resmat%Data)
       elseif (present(transA_in).and.present(transB_in)) then
          CALL dot_multiply(this%Data,mat%Data,resmat%Data,transA_in,transB_in)
       elseif (present(transA_in)) then
          CALL dot_multiply(this%Data,mat%Data,resmat%Data,transA_in)
       endif
    ELSE
       PRINT*,"Matrix shapes for multiplication with Matrix do not match, ",&
            &this%NCols,this%NRows,mat%NCols,mat%NRows,resmat%NCols,resmat%NRows
       stop
    END IF
  END SUBROUTINE mp_dot_multiply_MatMat_real

  SUBROUTINE mp_set_zero_real(mat)
    type(MatrixReal) :: mat

    ! Local variables
    INTEGER :: icol, irow

    DEBSTART("mp_set_zero_real")
    DO irow=1,mat%NRows
       DO icol=1,mat%NCols
          CALL set_value(mat%Data,irow,icol,0.0)
       END DO
    END DO
    !CALL commit_values(mat)
    DEBEND("mp_set_zero_real")
  END SUBROUTINE mp_set_zero_real

  SUBROUTINE mp_get_global_matrix_locally_real(mat,localfullmat)
    TYPE(MatrixReal),INTENT(IN) :: mat
    REAL, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    CALL get_global_matrix_locally(mat%Data,localfullmat)
  END SUBROUTINE mp_get_global_matrix_locally_real

  FUNCTION mp_get_value_real(mat,irow,icol) result(value)
    type(MatrixReal) :: mat
    INTEGER :: irow, icol
    real :: value

    value = mat_get_value(mat%Data,irow,icol)
  END FUNCTION mp_get_value_real

  FUNCTION mp_get_local_abs_square_sum_real(mat) RESULT(value)
    type(MatrixReal) :: mat
    real :: value

    value = get_local_abs_square_sum(mat%Data)
  END FUNCTION mp_get_local_abs_square_sum_real

  SUBROUTINE mp_show_on_screen_real(mat)
    type(MatrixReal) :: mat
    
    CALL show(mat%Data)
  END SUBROUTINE mp_show_on_screen_real

  SUBROUTINE mp_show_in_file_real(mat,filename)
    type(MatrixReal) :: mat
    CHARACTER(len=*) :: filename

    CALL show(mat%Data,filename)
  END SUBROUTINE mp_show_in_file_real

  SUBROUTINE mp_output_data_matrix_real(basename,mat)
    type(MatrixReal) :: mat
    CHARACTER(len=*) :: basename
    
    character(len=100) :: filename

    WRITE(filename,"(3A)") "./",TRIM(basename),".dat"
    PRINT*,"Writing to file ",TRIM(filename)
    CALL show(mat%Data,filename)
  END SUBROUTINE mp_output_data_matrix_real

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix_real(this,mat)
    TYPE(MatrixReal), intent(INOUT) :: this
    TYPE(MatrixReal),INTENT(IN) :: mat
    CALL add_matrix(this%Data,mat%Data)
  END SUBROUTINE mp_add_to_matrix_real

  SUBROUTINE mp_add_matrix_real(this,mat,res)
    TYPE(MatrixReal), INTENT(IN) :: this,mat
    TYPE(MatrixReal), INTENT(INOUT) :: res
    !PRINT*,"in mp_add_matrix of MatrixModule"
    CALL add_matrix(this%Data,mat%Data,res%Data)
  END SUBROUTINE mp_add_matrix_real

  SUBROUTINE mp_subtract_matrix_real(this,mat,res)
    TYPE(MatrixReal),INTENT(IN) :: this,mat
    TYPE(MatrixReal),intent(INOUT) :: res
    ! calculate res = this - mat

    CALL subtract_matrix(this%Data,mat%Data, res%Data)
  END SUBROUTINE mp_subtract_matrix_real

  SUBROUTINE mp_subtract_from_matrix_real(this,mat)
    TYPE(MatrixReal),INTENT(INOUT) :: this
    TYPE(MatrixReal), intent(IN) :: mat
    ! calculate this = this - mat

    CALL subtract_matrix(this%Data,mat%Data)
  END SUBROUTINE mp_subtract_from_matrix_real
    
  SUBROUTINE mp_multiply_matrix_with_real_real(this, scalar, res)
    TYPE(MatrixReal), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(MatrixReal),intent(INOUT) :: res

    CALL multiply_matrix_with_scalar(this%Data,scalar, res%Data)
  END SUBROUTINE mp_multiply_matrix_with_real_real

  SUBROUTINE mp_scale_matrix_by_real_real(this, scalar)
    TYPE(MatrixReal), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL multiply_matrix_with_scalar(this%Data,scalar)
  END SUBROUTINE mp_scale_matrix_by_real_real

  SUBROUTINE mp_assign_matrix_real(lmat,rmat)
    TYPE(MatrixReal), intent(INOUT) :: lmat
    TYPE(MatrixReal), intent(IN) :: rmat

    IF ((lmat%Nrows.EQ.rmat%Nrows).AND.(lmat%Ncols.EQ.rmat%Ncols)) THEN
       lmat%Data = rmat%Data
    ELSE
       PRINT*, "Left and right matrices in an assignment must have the same dimensions."
    END IF
  END SUBROUTINE mp_assign_matrix_real

  ! ====================================================
  ! == Define some mathematical routines for matrices ==
  ! ====================================================
  SUBROUTINE mp_invert_matrix_real(mat)
    TYPE(MatrixReal) :: mat

    CALL invert(mat%Data)
  END SUBROUTINE mp_invert_matrix_real

  SUBROUTINE mp_LU_factor_real(mat)
    TYPE(MatrixReal) :: mat

    IF (mat%NRows.EQ.mat%NCOLS) THEN
      CALL LU_factor(mat%Data)
    ELSE
      PRINT*,"Matrix for _LU_factor is not square"
      STOP
    ENDIF
  END SUBROUTINE mp_LU_factor_real

  SUBROUTINE mp_LU_solve_real(this, mat, resmat, phase)
    TYPE(MatrixReal),intent(INOUT) :: this
    TYPE(MatrixReal),intent(IN) :: mat
    TYPE(MatrixReal),intent(INOUT) :: resmat
    INTEGER,intent(IN) :: phase

    IF ( this%NRows.EQ.this%NCols .AND. &
         this%NRows.EQ.mat%NRows .AND.&
         this%NRows.EQ.resmat%NRows .AND.&
         mat%NCols.EQ.resmat%NCols ) THEN
       CALL LU_solve(this%Data,mat%Data,resmat%Data,phase)
    ELSE
       PRINT*,"Matrix shapes for LU_solve do not match, ",&
            &this%NCols,this%NRows,mat%NCols,mat%NRows,resmat%NCols,resmat%NRows
       STOP
    END IF
  END SUBROUTINE mp_LU_solve_real

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_matrix_sum_0D_real(localsum,communicator)
    !INTEGER, intent(IN) :: len
    TYPE(MatrixReal), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_sum(localsum%Data,communicator)
  END SUBROUTINE my_matrix_sum_0D_real

  SUBROUTINE my_matrix_sum_2D_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(MatrixReal), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_matrix_sum_generic_real(localsum,len,communicator)
  END SUBROUTINE my_matrix_sum_2D_real

  SUBROUTINE my_matrix_sum_3D_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(MatrixReal), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_matrix_sum_generic_real(localsum,len,communicator)
  END SUBROUTINE my_matrix_sum_3D_real
  
  SUBROUTINE my_matrix_sum_generic_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(MatrixReal), intent(INOUT) :: localsum(len)
    integer :: communicator

    integer :: iMatrix
    DO iMatrix=1,len
       CALL my_sum(localsum(iMatrix)%Data,communicator)
    END DO
  END SUBROUTINE my_matrix_sum_generic_real

END MODULE MatrixModule

