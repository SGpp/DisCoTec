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

!> This is a frontend module for the usage of a BandedMatrix. It more or less has the interfaces to the 
!! underlying module, which does the real storage. At the moment the underlying module is 
!! StoreBandedMatrixModule, which handles the storage.
MODULE BandedMatrixModule
  USE ProcessGridModule
  USE StoreBandedMatrixModule
  USE MatrixModule
  use VectorModule
  !USE mpi
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: BandedMatrix, BandedMatrixReal
  PUBLIC :: initialize, initialize_BandedMatrix_module, finalize, finalize_BandedMatrix_module
  PUBLIC :: allocate, attach, set_value, add_value, commit_values, set_zero
  PUBLIC :: get_global_matrix_locally, convert_Banded_to_Full, convert_Full_to_Banded, &
       &get_local_abs_square_sum, mat_get_value, mat_get_row_pointer
  PUBLIC :: autotune,dot_multiply, show, add_matrix, subtract_matrix, multiply_matrix_with_scalar
  PUBLIC :: invert, output_Data, isInitialized, my_sum, ASSIGNMENT(=)
  PUBLIC :: static_size_of_BandedMatrix, get_number_of_bands, row_axpy
  public :: transpose_and_conjugate, transpose_storage
  PUBLIC :: LU_factor, LU_factor_ok, LU_solve
  ! the following variables are used from ProcessGridModule. To use them outside the 
  ! MatrixModule, we have to publicize them explicitly.
  PUBLIC :: MAT_BLOCKS_OF_ROWS, MAT_BLOCKS_OF_COLS, MAT_BOTH

#include "BandedMatrix.h"

  ! Use the same process grid for all Matrices instantiated whithin this module
  TYPE(ProcessGrid),SAVE :: PG

CONTAINS
  FUNCTION static_size_of_BandedMatrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = 2*SIZE_OF_INTEGER &
         & + static_size_of_storeBandedmatrixobject()
    
  END FUNCTION static_size_of_BandedMatrix
    
  SUBROUTINE mp_initialize_matrix(mat,n_rows,n_cols,transp)
    Type(BandedMatrix), INTENT(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    LOGICAL, OPTIONAL :: transp

    DEBSTART("mp_initialize_matrix")
    IF (PG%isInitialized) THEN
       mat%NRows = n_rows
       mat%NCols = n_cols
       IF (PRESENT(transp)) THEN
          CALL initialize(mat%Data,n_rows,n_cols, PG,transp)
       ELSE
          CALL initialize(mat%Data,n_rows,n_cols, PG)
       END IF
    ELSE
       PRINT*,"Before initializing a matrix, the process grid has to be initialized."
       STOP
    END IF
    DEBEND("mp_initialize_matrix")
  END SUBROUTINE mp_initialize_matrix

  SUBROUTINE mp_initialize_vector(vec,n_rows)
    type(BandedMatrix) :: vec
    integer :: n_rows

    CALL initialize(vec,n_rows,1)
  END SUBROUTINE mp_initialize_vector

  SUBROUTINE mp_initialize_BandedMatrix_module
    CALL get_processgrid(PG)
  END SUBROUTINE mp_initialize_BandedMatrix_module

  LOGICAL FUNCTION mp_isInitialized(mat) 
    type(BandedMatrix) :: mat

    mp_isInitialized = isInitialized(mat%Data)
  END FUNCTION mp_isInitialized

  SUBROUTINE mp_allocate(mat)
    Type(BandedMatrix), intent(INOUT) :: mat

    CALL allocate(mat%Data)
  END SUBROUTINE mp_allocate

  SUBROUTINE mp_finalize_matrix(mat)
    type(BandedMatrix) :: mat
    call finalize(mat%Data)
  END SUBROUTINE mp_finalize_matrix

  SUBROUTINE mp_finalize_BandedMatrix_module()
    !CALL finalize(PG)
  END SUBROUTINE mp_finalize_BandedMatrix_module

  SUBROUTINE mp_set_complex_value(mat,irow,icol,value)
    Type(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_set_complex_value

  SUBROUTINE mp_set_real_value(mat,irow,icol,value)
    Type(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(mat%Data,irow,icol,CMPLX(value,0.0,KIND(value)))
  END SUBROUTINE mp_set_real_value

  SUBROUTINE mp_add_complex_value(mat,irow,icol,value)
    Type(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    CALL add_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_add_complex_value

  SUBROUTINE mp_add_real_value(mat,irow,icol,value)
    Type(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    CALL add_value(mat%Data,irow,icol,CMPLX(value,0.0,KIND(value)))
  END SUBROUTINE mp_add_real_value


  SUBROUTINE mp_commit_values(mat)
    Type(BandedMatrix) :: mat
    
    CALL commit_values(mat%Data)
  END SUBROUTINE mp_commit_values

  subroutine bm_autotune(mat,vec,res)
    Type(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res

    IF ((mat%NCols.EQ.vec%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       call autotune(mat%Data,vec%Data,res%Data)
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Vector do not match, ",&
            &mat%NCols,mat%NRows,vec%NRows,res%NRows
       stop
    END IF

  end subroutine bm_autotune

  SUBROUTINE mp_dot_multiply(mat,vec, res)
    Type(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res

    !PRINT*,vec%NCols,vec%NRows,res%NCols,res%NRows
    IF ((mat%NCols.EQ.vec%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       CALL dot_multiply(mat%Data,vec%Data,res%Data)
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Vector do not match, ",&
            &mat%NCols,mat%NRows,vec%NRows,res%NRows
       stop
    END IF
  END SUBROUTINE mp_dot_multiply

  subroutine bm_matmat_dot_multiply(mat,mat2,res,transA_in,transB_in)
    type(BandedMatrix),intent(IN) :: mat, mat2
    type(BandedMatrix),intent(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in, transB_in

    IF ((mat%NCols.EQ.mat2%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       if (.not.(present(transA_in).or.present(transB_in))) then
          CALL dot_multiply(mat%Data,mat2%Data,res%Data)
       elseif (present(transA_in).and.present(transB_in)) then
          CALL dot_multiply(mat%Data,mat2%Data,res%Data,transA_in,transB_in)
       elseif (present(transA_in)) then
          CALL dot_multiply(mat%Data,mat2%Data,res%Data,transA_in)
       endif
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Matrix do not match, ",&
            &mat%NCols,mat%NRows,mat2%NRows,res%NRows
       stop
    END IF
  end subroutine bm_matmat_dot_multiply

  subroutine bm_dot_multiply_Banded_with_Full(bmat,mat,res,transA_in,transB_in)
    type(BandedMatrix),intent(IN) :: bmat
    type(Matrix), intent(IN) :: mat
    type(Matrix),intent(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in, transB_in

    IF ((bmat%NCols.EQ.mat%NRows).AND. &
         &(bmat%NRows.EQ.res%NRows)) THEN
       if (.not.(present(transA_in).or.present(transB_in))) then
          CALL dot_multiply(bmat%Data,mat%Data,res%Data)
       elseif (present(transA_in).and.present(transB_in)) then
          CALL dot_multiply(bmat%Data,mat%Data,res%Data,transA_in,transB_in)
       elseif (present(transA_in)) then
          CALL dot_multiply(bmat%Data,mat%Data,res%Data,transA_in)
       endif
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Matrix do not match, ",&
            &bmat%NCols,bmat%NRows,mat%NRows,res%NRows
       stop
    END IF
  end subroutine bm_dot_multiply_Banded_with_Full

  SUBROUTINE bm_square_matrix(bmat,res_bmat)
    TYPE(BandedMatrix),intent(IN) :: bmat
    TYPE(BandedMatrix),intent(INOUT) :: res_bmat

    CALL square_matrix(bmat%Data,res_bmat%Data)
  END SUBROUTINE bm_square_matrix

  SUBROUTINE mp_set_zero(mat)
    type(BandedMatrix) :: mat

    DEBSTART("mp_set_zero")
    call set_zero(mat%Data)
    DEBEND("mp_set_zero")
  END SUBROUTINE mp_set_zero

  SUBROUTINE mp_get_global_matrix_locally(mat,localfullmat)
    Type(BandedMatrix),INTENT(IN) :: mat
    COMPLEX, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    CALL get_global_matrix_locally(mat%Data,localfullmat)
  END SUBROUTINE mp_get_global_matrix_locally

  SUBROUTINE mp_convert_Banded_to_Full(bmat, fmat)
    TYPE(BandedMatrix), intent(IN) :: bmat
    TYPE(Matrix), intent(INOUT) :: fmat

    INTEGER :: irow, icol, RowsPerRank

    IF ((fmat%NRows.NE.bmat%NRows).OR.(fmat%NCols.NE.bmat%NCols)) THEN
       PRINT*,"Number of Rows and Number of Cols must be equal for full and banded matrix."
       stop
    END IF

    RowsPerRank = bmat%NRows/PG%NProcs

    IF (bmat%data%isTransposed) THEN
       DO irow=PG%Rank*RowsPerRank+1,(PG%Rank+1)*RowsPerRank
          !DO irow=1,bmat%NRows
          DO icol=1,bmat%NCols
             CALL set_value(fmat,irow,icol,mat_get_value(bmat,irow,icol))
          END DO
       END DO
    ELSE
       PRINT*, "Can only convert banded to full matrix if banded matrix is in transposed storage mode."
       stop
    END IF
  END SUBROUTINE mp_convert_Banded_to_Full

  SUBROUTINE mp_convert_Full_to_Banded(fmat, bmat)
    TYPE(Matrix), intent(IN) :: fmat
    TYPE(BandedMatrix), intent(INOUT) :: bmat

    INTEGER :: irow, icol, RowsPerRank

    IF ((fmat%NRows.NE.bmat%NRows).OR.(fmat%NCols.NE.bmat%NCols)) THEN
       PRINT*,"Number of Rows and Number of Cols must be equal for full and banded matrix."
       stop
    END IF

    call set_zero(bmat)
    RowsPerRank = bmat%NRows/PG%NProcs

    IF (bmat%data%isTransposed) THEN
       DO irow=PG%Rank*RowsPerRank+1,(PG%Rank+1)*RowsPerRank
          DO icol=1,bmat%NCols
             CALL set_value(bmat,irow,icol,mat_get_value(fmat,irow,icol))
          END DO
       END DO
       CALL commit_values(bmat)
    ELSE
       PRINT*, "Can only convert full to banded matrix if banded matrix is in transposed storage mode."
       stop
    END IF
  END SUBROUTINE mp_convert_Full_to_Banded

  FUNCTION mp_get_value(mat,irow,icol,trans_in) result(value)
    type(BandedMatrix) :: mat
    INTEGER :: irow, icol
    complex :: value
    character(len=1),optional,intent(IN) :: trans_in

    if (present(trans_in)) then
       value = mat_get_value(mat%Data,irow,icol,trans_in)
    else
       value = mat_get_value(mat%Data,irow,icol,'N')
    end if
  END FUNCTION mp_get_value

  subroutine bm_get_row_pointer(mat,global_row,ptr_row,s_global_col,e_global_col)
    TYPE(BandedMatrix) :: mat
    INTEGER,intent(IN) :: global_row
    complex,dimension(:),pointer,intent(OUT) :: ptr_row
    integer, INTENT(OUT) :: s_global_col, e_global_col

    call mat_get_row_pointer(mat%Data,global_row,ptr_row,s_global_col,e_global_col)
  end subroutine bm_get_row_pointer

  FUNCTION mp_get_local_abs_square_sum(mat) RESULT(value)
    type(BandedMatrix) :: mat
    real :: value

    value = get_local_abs_square_sum(mat%Data)
  END FUNCTION mp_get_local_abs_square_sum

  SUBROUTINE mp_show_on_screen(mat)
    type(BandedMatrix) :: mat
    
    CALL show(mat%Data)
  END SUBROUTINE mp_show_on_screen

  SUBROUTINE mp_show_in_file(mat,filename)
    type(BandedMatrix) :: mat
    CHARACTER(len=*) :: filename

    CALL show(mat%Data,filename)
  END SUBROUTINE mp_show_in_file

  SUBROUTINE mp_output_data_matrix(basename,mat)
    type(BandedMatrix) :: mat
    CHARACTER(len=*) :: basename
    
    character(len=FILENAME_MAX) :: filename

    WRITE(filename,"(3A)") "./",TRIM(basename),".dat"
    PRINT*,"Writing to file ",TRIM(filename)
    CALL show(mat,filename)
  END SUBROUTINE mp_output_data_matrix

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix(this,mat)
    Type(BandedMatrix), intent(INOUT) :: this
    Type(BandedMatrix),INTENT(IN) :: mat
    CALL add_matrix(this%Data,mat%Data)
  END SUBROUTINE mp_add_to_matrix

  SUBROUTINE mp_add_matrix(this,mat,res)
    Type(BandedMatrix), INTENT(IN) :: this,mat
    Type(BandedMatrix), INTENT(INOUT) :: res
    !PRINT*,"in mp_add_matrix of MatrixModule"
    CALL add_matrix(this%Data,mat%Data,res%Data)
  END SUBROUTINE mp_add_matrix

  SUBROUTINE mp_subtract_matrix(this,mat,res)
    Type(BandedMatrix),INTENT(IN) :: this,mat
    Type(BandedMatrix),intent(INOUT) :: res
    ! calculate res = this - mat

    CALL subtract_matrix(this%Data,mat%Data, res%Data)
  END SUBROUTINE mp_subtract_matrix

  SUBROUTINE mp_subtract_from_matrix(this,mat)
    Type(BandedMatrix),INTENT(INOUT) :: this
    Type(BandedMatrix), intent(IN) :: mat
    ! calculate this = this - mat

    CALL subtract_matrix(this%Data,mat%Data)
  END SUBROUTINE mp_subtract_from_matrix
    
  SUBROUTINE mp_multiply_matrix_with_real(this, scalar, res)
    Type(BandedMatrix), intent(IN) :: this
    REAL, intent(IN) :: scalar
    Type(BandedMatrix),intent(INOUT) :: res

    CALL multiply_matrix_with_scalar(this%Data,scalar, res%Data)
  END SUBROUTINE mp_multiply_matrix_with_real

  SUBROUTINE mp_scale_matrix_by_real(this, scalar)
    Type(BandedMatrix), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL multiply_matrix_with_scalar(this%Data,scalar)
  END SUBROUTINE mp_scale_matrix_by_real

  SUBROUTINE mp_assign_matrix(lmat,rmat)
    Type(BandedMatrix), intent(INOUT) :: lmat
    Type(BandedMatrix), intent(IN) :: rmat

    IF ((lmat%Nrows.EQ.rmat%Nrows).AND.(lmat%Ncols.EQ.rmat%Ncols)) THEN
       lmat%Data = rmat%Data
    ELSE
       PRINT*, "Left and right matrices in an assignment must have the same dimensions."
    END IF
  END SUBROUTINE mp_assign_matrix

  SUBROUTINE bm_print_storage_details(bmat)
    TYPE(BandedMatrix), INTENT(IN) :: bmat

    CALL print_storage_details(bmat%Data)
  END SUBROUTINE bm_print_storage_details

  FUNCTION bm_get_number_of_bands(mat)
    TYPE(BandedMatrix), intent(IN) :: mat
    integer :: bm_get_number_of_bands

    bm_get_number_of_bands = get_number_of_bands(mat%Data)
  END FUNCTION bm_get_number_of_bands

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE bm_sum_0D(localsum,communicator)
!    INTEGER, intent(IN) :: len
    Type(BandedMatrix), INTENT(INOUT) :: localsum
    integer :: communicator

!    IF (len.NE.1) THEN
!       PRINT*,"BandedMatrix:bm_sum_0D: A zero-dimensional array of banded matrices must have a len of 1."
!       stop
!    END IF

    CALL my_sum(localsum%Data,communicator)
  END SUBROUTINE bm_sum_0D

  SUBROUTINE bm_sum_2D(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    Type(BandedMatrix), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL bm_sum_generic(localsum,len,communicator)
  END SUBROUTINE bm_sum_2D

  SUBROUTINE bm_sum_3D(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    Type(BandedMatrix), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL bm_sum_generic(localsum,len,communicator)
  END SUBROUTINE bm_sum_3D
  
  SUBROUTINE bm_sum_generic(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    Type(BandedMatrix), intent(INOUT) :: localsum(len)
    integer :: communicator

    integer :: iMatrix
    DO iMatrix=1,len
       CALL my_sum(localsum(iMatrix)%Data,communicator)
    END DO
  END SUBROUTINE bm_sum_generic

  SUBROUTINE bm_row_axpy(this,this_row,mat,scalar)
    TYPE(BandedMatrix),intent(INOUT) :: this
    INTEGER, INTENT(IN) :: this_row
    TYPE(BandedMatrix), intent(IN) :: mat
    COMPLEX, intent(IN) :: scalar

    IF ((mat%Nrows.NE.this%NRows).OR.(mat%NCols.NE.this%NCols)) THEN
       PRINT*,"row_axpy: mat and this must have the same shape."
       STOP
    END IF

    IF ((this_row.GE.1).AND.(this_row.LE.this%NRows)) THEN
       CALL row_axpy(this%data,this_row,mat%data,scalar)
    ELSE
       PRINT*,"row_axpy: this_row must be between 1 and ",this%NRows," but is ",this_row
       STOP
    END IF
  END SUBROUTINE bm_row_axpy

  subroutine bm_transpose_and_conjugate(mat,dagger_mat)
    type(BandedMatrix) :: mat,dagger_mat

    call transpose_and_conjugate(mat%data,dagger_mat%data)
  end subroutine bm_transpose_and_conjugate

  subroutine bm_transpose_storage(mat,mat_t)
    type(BandedMatrix) :: mat,mat_t

    call transpose_storage(mat%data,mat_t%data)
  end subroutine bm_transpose_storage

  SUBROUTINE bm_LU_factor(mat)
    TYPE(BandedMatrix),intent(INOUT) :: mat

    CALL LU_factor(mat%DATA)

  END SUBROUTINE bm_LU_factor

  LOGICAL FUNCTION bm_LU_factor_ok(mat)
    TYPE(BandedMatrix),intent(IN) :: mat

    bm_LU_factor_ok = LU_factor_ok(mat%DATA)

  END FUNCTION bm_LU_factor_ok

  SUBROUTINE bm_LU_solve(mat,vec,res)
    TYPE(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res

    CALL LU_solve(mat%DATA,vec%data,res%data)

  END SUBROUTINE bm_LU_solve

#ifdef ALL_ROUTINES
  ! ====================================================
  ! == Define some mathematical routines for matrices ==
  ! ====================================================
  SUBROUTINE mp_invert_matrix(mat)
    Type(BandedMatrix) :: mat

    CALL invert(mat%Data)
  END SUBROUTINE mp_invert_matrix
#endif


  SUBROUTINE mp_initialize_matrix_real(mat,n_rows,n_cols,transp)
    Type(BandedMatrixReal), INTENT(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    LOGICAL, OPTIONAL :: transp

    DEBSTART("mp_initialize_matrix_real")
    IF (PG%isInitialized) THEN
       mat%NRows = n_rows
       mat%NCols = n_cols
       IF (PRESENT(transp)) THEN
          CALL initialize(mat%Data,n_rows,n_cols, PG,transp)
       ELSE
          CALL initialize(mat%Data,n_rows,n_cols, PG)
       END IF
    ELSE
       PRINT*,"Before initializing a matrix, the process grid has to be initialized."
       STOP
    END IF
    DEBEND("mp_initialize_matrix_real")
  END SUBROUTINE mp_initialize_matrix_real

  SUBROUTINE mp_initialize_vector_real(vec,n_rows)
    type(BandedMatrixReal) :: vec
    integer :: n_rows

    CALL initialize(vec%DATA,n_rows,1,PG)
  END SUBROUTINE mp_initialize_vector_real

  LOGICAL FUNCTION mp_isInitialized_real(mat) 
    type(BandedMatrixReal) :: mat

    mp_isInitialized_real = isInitialized(mat%Data)
  END FUNCTION mp_isInitialized_real

  SUBROUTINE mp_allocate_real(mat)
    Type(BandedMatrixReal), intent(INOUT) :: mat

    CALL allocate(mat%Data)
  END SUBROUTINE mp_allocate_real

  SUBROUTINE mp_finalize_matrix_real(mat)
    type(BandedMatrixReal) :: mat
    call finalize(mat%Data)
  END SUBROUTINE mp_finalize_matrix_real

  SUBROUTINE mp_set_complex_value_real(mat,irow,icol,value)
    Type(BandedMatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_set_complex_value_real

  SUBROUTINE mp_set_real_value_real(mat,irow,icol,value)
    Type(BandedMatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_set_real_value_real

  SUBROUTINE mp_add_complex_value_real(mat,irow,icol,value)
    Type(BandedMatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    CALL add_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_add_complex_value_real

  SUBROUTINE mp_add_real_value_real(mat,irow,icol,value)
    Type(BandedMatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    CALL add_value(mat%Data,irow,icol,value)
  END SUBROUTINE mp_add_real_value_real


  SUBROUTINE mp_commit_values_real(mat)
    Type(BandedMatrixReal) :: mat
    
    CALL commit_values(mat%Data)
  END SUBROUTINE mp_commit_values_real

  subroutine bm_autotune_real(mat,vec,res)
    Type(BandedMatrixReal),intent(IN) :: mat
    Type(VectorReal),intent(IN) :: vec
    Type(VectorReal),intent(INOUT) :: res

    IF ((mat%NCols.EQ.vec%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       call autotune(mat%Data,vec%Data,res%Data)
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Vector do not match, ",&
            &mat%NCols,mat%NRows,vec%NRows,res%NRows
       stop
    END IF

  end subroutine bm_autotune_real

  SUBROUTINE mp_dot_multiply_real(mat,vec, res)
    Type(BandedMatrixReal),intent(IN) :: mat
    Type(VectorReal),intent(IN) :: vec
    Type(VectorReal),intent(INOUT) :: res

    !PRINT*,vec%NCols,vec%NRows,res%NCols,res%NRows
    IF ((mat%NCols.EQ.vec%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       CALL dot_multiply(mat%Data,vec%Data,res%Data)
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Vector do not match, ",&
            &mat%NCols,mat%NRows,vec%NRows,res%NRows
       stop
    END IF
  END SUBROUTINE mp_dot_multiply_real

  subroutine bm_matmat_dot_multiply_real(mat,mat2,res,transA_in,transB_in)
    type(BandedMatrixReal),intent(IN) :: mat, mat2
    type(BandedMatrixReal),intent(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in, transB_in

    IF ((mat%NCols.EQ.mat2%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       if (.not.(present(transA_in).or.present(transB_in))) then
          CALL dot_multiply(mat%Data,mat2%Data,res%Data)
       elseif (present(transA_in).and.present(transB_in)) then
          CALL dot_multiply(mat%Data,mat2%Data,res%Data,transA_in,transB_in)
       elseif (present(transA_in)) then
          CALL dot_multiply(mat%Data,mat2%Data,res%Data,transA_in)
       endif
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Matrix do not match, ",&
            &mat%NCols,mat%NRows,mat2%NRows,res%NRows
       stop
    END IF
  end subroutine bm_matmat_dot_multiply_real

  subroutine bm_dot_multiply_Banded_with_Full_real(bmat,mat,res,transA_in,transB_in)
    type(BandedMatrixReal),intent(IN) :: bmat
    type(MatrixReal), intent(IN) :: mat
    type(MatrixReal),intent(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in, transB_in

    IF ((bmat%NCols.EQ.mat%NRows).AND. &
         &(bmat%NRows.EQ.res%NRows)) THEN
       if (.not.(present(transA_in).or.present(transB_in))) then
          CALL dot_multiply(bmat%Data,mat%Data,res%Data)
       elseif (present(transA_in).and.present(transB_in)) then
          CALL dot_multiply(bmat%Data,mat%Data,res%Data,transA_in,transB_in)
       elseif (present(transA_in)) then
          CALL dot_multiply(bmat%Data,mat%Data,res%Data,transA_in)
       endif
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Matrix do not match, ",&
            &bmat%NCols,bmat%NRows,mat%NRows,res%NRows
       stop
    END IF
  end subroutine bm_dot_multiply_Banded_with_Full_real

  SUBROUTINE bm_square_matrix_real(bmat,res_bmat)
    TYPE(BandedMatrixReal),intent(IN) :: bmat
    TYPE(BandedMatrixReal),intent(INOUT) :: res_bmat

    CALL square_matrix(bmat%Data,res_bmat%Data)
  END SUBROUTINE bm_square_matrix_real

  SUBROUTINE mp_set_zero_real(mat)
    type(BandedMatrixReal) :: mat

    DEBSTART("mp_set_zero_real")
    call set_zero(mat%Data)
    DEBEND("mp_set_zero_real")
  END SUBROUTINE mp_set_zero_real

  SUBROUTINE mp_get_global_matrix_locally_real(mat,localfullmat)
    Type(BandedMatrixReal),INTENT(IN) :: mat
    REAL, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    CALL get_global_matrix_locally(mat%Data,localfullmat)
  END SUBROUTINE mp_get_global_matrix_locally_real

  SUBROUTINE mp_convert_Banded_to_Full_real(bmat, fmat)
    TYPE(BandedMatrixReal), intent(IN) :: bmat
    TYPE(MatrixReal), intent(INOUT) :: fmat
    REAL :: temp

    INTEGER :: irow, icol, RowsPerRank

    IF ((fmat%NRows.NE.bmat%NRows).OR.(fmat%NCols.NE.bmat%NCols)) THEN
       PRINT*,"Number of Rows and Number of Cols must be equal for full and banded matrix."
       stop
    END IF

    RowsPerRank = bmat%NRows/PG%NProcs

    IF (bmat%data%isTransposed) THEN
       DO irow=PG%Rank*RowsPerRank+1,(PG%Rank+1)*RowsPerRank
          !DO irow=1,bmat%NRows
          DO icol=1,bmat%NCols
             temp=mat_get_value(bmat%DATA,irow,icol)
             CALL set_value(fmat,irow,icol,temp)
          END DO
       END DO
    ELSE
       PRINT*, "Can only convert banded to full matrix if banded matrix is in transposed storage mode."
       stop
    END IF
  END SUBROUTINE mp_convert_Banded_to_Full_real

  SUBROUTINE mp_convert_Full_to_Banded_real(fmat, bmat)
    TYPE(MatrixReal), intent(IN) :: fmat
    TYPE(BandedMatrixReal), intent(INOUT) :: bmat

    INTEGER :: irow, icol, RowsPerRank

    IF ((fmat%NRows.NE.bmat%NRows).OR.(fmat%NCols.NE.bmat%NCols)) THEN
       PRINT*,"Number of Rows and Number of Cols must be equal for full and banded matrix."
       stop
    END IF

    call set_zero(bmat%DATA)
    RowsPerRank = bmat%NRows/PG%NProcs

    IF (bmat%data%isTransposed) THEN
       DO irow=PG%Rank*RowsPerRank+1,(PG%Rank+1)*RowsPerRank
          DO icol=1,bmat%NCols
             CALL set_value(bmat%DATA,irow,icol,mat_get_value(fmat,irow,icol))
          END DO
       END DO
       CALL commit_values(bmat%DATA)
    ELSE
       PRINT*, "Can only convert full to banded matrix if banded matrix is in transposed storage mode."
       stop
    END IF
  END SUBROUTINE mp_convert_Full_to_Banded_real

  FUNCTION mp_get_value_real(mat,irow,icol,trans_in) result(value)
    type(BandedMatrixReal) :: mat
    INTEGER :: irow, icol
    complex :: value
    character(len=1),optional,intent(IN) :: trans_in

    if (present(trans_in)) then
       value = mat_get_value(mat%Data,irow,icol,trans_in)
    else
       value = mat_get_value(mat%Data,irow,icol,'N')
    end if
  END FUNCTION mp_get_value_real

  subroutine bm_get_row_pointer_real(mat,global_row,ptr_row,s_global_col,e_global_col)
    TYPE(BandedMatrixReal) :: mat
    INTEGER,intent(IN) :: global_row
    real,dimension(:),pointer,intent(OUT) :: ptr_row
    integer, INTENT(OUT) :: s_global_col, e_global_col

    call mat_get_row_pointer(mat%Data,global_row,ptr_row,s_global_col,e_global_col)
  end subroutine bm_get_row_pointer_real

  FUNCTION mp_get_local_abs_square_sum_real(mat) RESULT(value)
    type(BandedMatrixReal) :: mat
    real :: value

    value = get_local_abs_square_sum(mat%Data)
  END FUNCTION mp_get_local_abs_square_sum_real

  SUBROUTINE mp_show_on_screen_real(mat)
    type(BandedMatrixReal) :: mat
    
    CALL show(mat%Data)
  END SUBROUTINE mp_show_on_screen_real

  SUBROUTINE mp_show_in_file_real(mat,filename)
    type(BandedMatrixReal) :: mat
    CHARACTER(len=*) :: filename

    CALL show(mat%Data,filename)
  END SUBROUTINE mp_show_in_file_real

  SUBROUTINE mp_output_data_matrix_real(basename,mat)
    type(BandedMatrixReal) :: mat
    CHARACTER(len=*) :: basename
    
    character(len=FILENAME_MAX) :: filename

    WRITE(filename,"(3A)") "./",TRIM(basename),".dat"
    PRINT*,"Writing to file ",TRIM(filename)
    CALL show(mat%DATA,filename)
  END SUBROUTINE mp_output_data_matrix_real

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix_real(this,mat)
    Type(BandedMatrixReal), intent(INOUT) :: this
    Type(BandedMatrixReal),INTENT(IN) :: mat
    CALL add_matrix(this%Data,mat%Data)
  END SUBROUTINE mp_add_to_matrix_real

  SUBROUTINE mp_add_matrix_real(this,mat,res)
    Type(BandedMatrixReal), INTENT(IN) :: this,mat
    Type(BandedMatrixReal), INTENT(INOUT) :: res
    !PRINT*,"in mp_add_matrix of MatrixModule"
    CALL add_matrix(this%Data,mat%Data,res%Data)
  END SUBROUTINE mp_add_matrix_real

  SUBROUTINE mp_subtract_matrix_real(this,mat,res)
    Type(BandedMatrixReal),INTENT(IN) :: this,mat
    Type(BandedMatrixReal),intent(INOUT) :: res
    ! calculate res = this - mat

    CALL subtract_matrix(this%Data,mat%Data, res%Data)
  END SUBROUTINE mp_subtract_matrix_real

  SUBROUTINE mp_subtract_from_matrix_real(this,mat)
    Type(BandedMatrixReal),INTENT(INOUT) :: this
    Type(BandedMatrixReal), intent(IN) :: mat
    ! calculate this = this - mat

    CALL subtract_matrix(this%Data,mat%Data)
  END SUBROUTINE mp_subtract_from_matrix_real
    
  SUBROUTINE mp_multiply_matrix_with_real_real(this, scalar, res)
    Type(BandedMatrixReal), intent(IN) :: this
    REAL, intent(IN) :: scalar
    Type(BandedMatrixReal),intent(INOUT) :: res

    CALL multiply_matrix_with_scalar(this%Data,scalar, res%Data)
  END SUBROUTINE mp_multiply_matrix_with_real_real

  SUBROUTINE mp_scale_matrix_by_real_real(this, scalar)
    Type(BandedMatrixReal), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL multiply_matrix_with_scalar(this%Data,scalar)
  END SUBROUTINE mp_scale_matrix_by_real_real

  SUBROUTINE mp_assign_matrix_real(lmat,rmat)
    Type(BandedMatrixReal), intent(INOUT) :: lmat
    Type(BandedMatrixReal), intent(IN) :: rmat

    IF ((lmat%Nrows.EQ.rmat%Nrows).AND.(lmat%Ncols.EQ.rmat%Ncols)) THEN
       lmat%Data = rmat%Data
    ELSE
       PRINT*, "Left and right matrices in an assignment must have the same dimensions."
    END IF
  END SUBROUTINE mp_assign_matrix_real

  SUBROUTINE bm_print_storage_details_real(bmat)
    TYPE(BandedMatrixReal), INTENT(IN) :: bmat

    CALL print_storage_details(bmat%Data)
  END SUBROUTINE bm_print_storage_details_real

  FUNCTION bm_get_number_of_bands_real(mat)
    TYPE(BandedMatrixReal), intent(IN) :: mat
    integer :: bm_get_number_of_bands_real

    bm_get_number_of_bands_real = get_number_of_bands(mat%Data)
  END FUNCTION bm_get_number_of_bands_real

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE bm_sum_0D_real(localsum,communicator)
!    INTEGER, intent(IN) :: len
    Type(BandedMatrixReal), INTENT(INOUT) :: localsum
    integer :: communicator

!    IF (len.NE.1) THEN
!       PRINT*,"BandedMatrix:bm_sum_0D: A zero-dimensional array of banded matrices must have a len of 1."
!       stop
!    END IF

    CALL my_sum(localsum%Data,communicator)
  END SUBROUTINE bm_sum_0D_real

  SUBROUTINE bm_sum_2D_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    Type(BandedMatrixReal), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL bm_sum_generic_real(localsum,len,communicator)
  END SUBROUTINE bm_sum_2D_real

  SUBROUTINE bm_sum_3D_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    Type(BandedMatrixReal), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL bm_sum_generic_real(localsum,len,communicator)
  END SUBROUTINE bm_sum_3D_real
  
  SUBROUTINE bm_sum_generic_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    Type(BandedMatrixReal), intent(INOUT) :: localsum(len)
    integer :: communicator

    integer :: iMatrix
    DO iMatrix=1,len
       CALL my_sum(localsum(iMatrix)%Data,communicator)
    END DO
  END SUBROUTINE bm_sum_generic_real

  SUBROUTINE bm_row_axpy_real(this,this_row,mat,scalar)
    TYPE(BandedMatrixReal),intent(INOUT) :: this
    INTEGER, INTENT(IN) :: this_row
    TYPE(BandedMatrixReal), intent(IN) :: mat
    REAL, intent(IN) :: scalar

    IF ((mat%Nrows.NE.this%NRows).OR.(mat%NCols.NE.this%NCols)) THEN
       PRINT*,"row_axpy: mat and this must have the same shape."
       STOP
    END IF

    IF ((this_row.GE.1).AND.(this_row.LE.this%NRows)) THEN
       CALL row_axpy(this%data,this_row,mat%data,scalar)
    ELSE
       PRINT*,"row_axpy: this_row must be between 1 and ",this%NRows," but is ",this_row
       STOP
    END IF
  END SUBROUTINE bm_row_axpy_real

  subroutine bm_transpose_and_conjugate_real(mat,dagger_mat)
    type(BandedMatrixReal) :: mat,dagger_mat

    call transpose_and_conjugate(mat%data,dagger_mat%data)
  end subroutine bm_transpose_and_conjugate_real

  subroutine bm_transpose_storage_real(mat,mat_t)
    type(BandedMatrixReal) :: mat,mat_t

    call transpose_storage(mat%data,mat_t%data)
  end subroutine bm_transpose_storage_real

  SUBROUTINE bm_LU_factor_real(mat)
    TYPE(BandedMatrixReal),intent(INOUT) :: mat

    CALL LU_factor(mat%DATA)

  END SUBROUTINE bm_LU_factor_real

  LOGICAL FUNCTION bm_LU_factor_ok_real(mat)
    TYPE(BandedMatrixReal),intent(IN) :: mat

    bm_LU_factor_ok_real = LU_factor_ok(mat%DATA)

  END FUNCTION bm_LU_factor_ok_real

  SUBROUTINE bm_LU_solve_real(mat,vec,res)
    TYPE(BandedMatrixReal),intent(IN) :: mat
    Type(VectorReal),intent(IN) :: vec
    Type(VectorReal),intent(INOUT) :: res

    CALL LU_solve(mat%DATA,vec%data,res%data)

  END SUBROUTINE bm_LU_solve_real

#ifdef ALL_ROUTINES
  ! ====================================================
  ! == Define some mathematical routines for matrices ==
  ! ====================================================
  SUBROUTINE mp_invert_matrix_real(mat)
    Type(BandedMatrixReal) :: mat

    CALL invert(mat%Data)
  END SUBROUTINE mp_invert_matrix_real
#endif

END MODULE BandedMatrixModule
