#include "switches.h"
#include "intrinsic_sizes.h"
#include "redef.h"
!****h* /MatrixModule
! DESCRIPTION
!****
MODULE MatrixModule
  !USE StoreFullMatrixModule
  !USE mpi
  IMPLICIT NONE

  INTEGER, PARAMETER,PUBLIC :: MAT_BLOCKS_OF_ROWS=1000, MAT_BLOCKS_OF_COLS=1001, MAT_BOTH=1002

  TYPE Matrix
     INTEGER :: NCols
     INTEGER :: NRows
  END TYPE matrix

  TYPE MatrixReal
     INTEGER :: NCols
     INTEGER :: NRows       
  END TYPE MatrixReal    

  INTERFACE initialize
     MODULE PROCEDURE mp_initialize_matrix, mp_initialize_vector
     MODULE PROCEDURE mp_initialize_matrix_real, mp_initialize_vector_real
  END INTERFACE

  INTERFACE initialize_matrix_module
     MODULE PROCEDURE mp_initialize_matrix_module
  END INTERFACE

  INTERFACE allocate
     MODULE PROCEDURE mp_allocate, mp_allocate_real
  END INTERFACE

  INTERFACE attach
     MODULE PROCEDURE mp_attach_vector, mp_attach_matrix
     MODULE PROCEDURE mp_attach_vector_real, mp_attach_matrix_real
  END INTERFACE

  INTERFACE finalize
     MODULE PROCEDURE mp_finalize_matrix, mp_finalize_matrix_real
  END INTERFACE

  INTERFACE finalize_matrix_module
     MODULE PROCEDURE mp_finalize_module
  END INTERFACE

  INTERFACE set_value
     MODULE PROCEDURE mp_set_complex_value, mp_set_real_value
     MODULE PROCEDURE mp_set_real_value_real
  END INTERFACE

  INTERFACE add_value
     MODULE PROCEDURE mp_add_complex_value, mp_add_real_value
     MODULE PROCEDURE mp_add_real_value_real
  END INTERFACE

  INTERFACE commit_values
     module procedure mp_commit_values, mp_commit_values_real
  END INTERFACE

  INTERFACE set_zero
     MODULE PROCEDURE mp_set_zero, mp_set_zero_real
  END INTERFACE

  INTERFACE get_global_matrix_locally
     module procedure mp_get_global_matrix_locally
     MODULE PROCEDURE mp_get_global_matrix_locally_real
  END INTERFACE

  INTERFACE get_local_abs_square_sum
     module procedure mp_get_local_abs_square_sum
     MODULE PROCEDURE mp_get_local_abs_square_sum_real
  END INTERFACE

  INTERFACE mat_get_value
     MODULE PROCEDURE mp_get_value, mp_get_value_real
  END INTERFACE

  INTERFACE dot_multiply
     MODULE PROCEDURE mp_dot_multiply, mp_dot_multiply_real
  END INTERFACE

  INTERFACE show
     MODULE PROCEDURE mp_show_on_screen, mp_show_in_file
     MODULE PROCEDURE mp_show_on_screen_real, mp_show_in_file_real
  END INTERFACE

  INTERFACE add_matrix
     MODULE PROCEDURE mp_add_to_matrix, mp_add_matrix
     MODULE PROCEDURE mp_add_to_matrix_real, mp_add_matrix_real
  END INTERFACE

  interface subtract_matrix
  !INTERFACE OPERATOR(-)
     MODULE PROCEDURE mp_subtract_matrix, mp_subtract_from_matrix
     MODULE PROCEDURE mp_subtract_matrix_real, mp_subtract_from_matrix_real
  END INTERFACE

  interface multiply_matrix_with_scalar
  !INTERFACE OPERATOR(*)
     module procedure mp_multiply_matrix_with_real, mp_multiply_matrix_with_real_real
     module procedure mp_scale_matrix_by_real, mp_scale_matrix_by_real_real
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
     module procedure mp_assign_matrix, mp_assign_matrix_real
  END INTERFACE

  INTERFACE invert
     module procedure mp_invert_matrix, mp_invert_matrix_real
  END INTERFACE

  INTERFACE output_data
     module procedure mp_output_data_matrix, mp_output_data_matrix_real
  END INTERFACE

  INTERFACE isInitialized
     module procedure mp_isInitialized, mp_isInitialized_real
  END INTERFACE

  INTERFACE my_sum
     MODULE PROCEDURE my_matrix_sum_generic, my_matrix_sum_2D, &
          &my_matrix_sum_3D, my_matrix_sum_0D
     MODULE PROCEDURE my_matrix_sum_generic_real, my_matrix_sum_2D_real, &
          &my_matrix_sum_3D_real, my_matrix_sum_0D_real
  END INTERFACE

  INTERFACE LU_factor
     MODULE PROCEDURE mp_LU_factor, mp_LU_factor_real
  END INTERFACE

  INTERFACE LU_solve
     MODULE PROCEDURE mp_LU_solve, mp_LU_solve_real
  END INTERFACE

  PRIVATE :: mp_set_complex_value, mp_set_real_value
  PRIVATE :: mp_set_real_value_real
  PRIVATE :: mp_add_complex_value, mp_add_real_value
  PRIVATE :: mp_add_real_value_real
  PRIVATE :: mp_initialize_matrix, mp_initialize_vector, mp_initialize_matrix_module, &
       & mp_finalize_matrix, mp_finalize_module, mp_commit_values, mp_set_zero
  PRIVATE :: mp_initialize_matrix_real, mp_initialize_vector_real, &
       & mp_finalize_matrix_real, mp_commit_values_real, mp_set_zero_real        
  PRIVATE :: mp_dot_multiply, mp_attach_vector, mp_attach_matrix, &
       &mp_show_on_screen, mp_show_in_file, mp_allocate
  PRIVATE :: mp_dot_multiply_real, mp_attach_vector_real, mp_attach_matrix_real, &
       & mp_show_on_screen_real, mp_show_in_file_real, mp_allocate_real     
  PRIVATE :: mp_get_global_matrix_locally, mp_get_value
  PRIVATE :: mp_get_global_matrix_locally_real, mp_get_value_real

  PRIVATE :: mp_add_to_matrix, mp_add_matrix, mp_subtract_matrix, mp_subtract_from_matrix
  PRIVATE :: mp_add_to_matrix_real, mp_add_matrix_real, mp_subtract_matrix_real, mp_subtract_from_matrix_real 
  PRIVATE :: mp_multiply_matrix_with_real, &
       & mp_scale_matrix_by_real
  PRIVATE :: mp_multiply_matrix_with_real_real, &
       &mp_scale_matrix_by_real_real    

  PRIVATE :: mp_get_local_abs_square_sum, mp_get_local_abs_square_sum_real 
  PRIVATE :: mp_invert_matrix, mp_invert_matrix_real 
  private :: mp_output_data_matrix, mp_output_data_matrix_real
  private :: mp_assign_matrix, mp_assign_matrix_real
  private :: my_matrix_sum_generic, my_matrix_sum_2D, my_matrix_sum_3D, my_matrix_sum_0D
  private :: my_matrix_sum_generic_real, my_matrix_sum_2D_real, my_matrix_sum_3D_real, my_matrix_sum_0D_real
  private :: mp_isInitialized, mp_isInitialized_real

CONTAINS
  SUBROUTINE mp_initialize_matrix(mat,n_rows,n_cols)
    TYPE(Matrix), INTENT(INOUT) :: mat
    INTEGER :: n_rows, n_cols

    PRINT*, "You cannot run the global version of gene without the SCALAPACK switch in the makefile switched on."
    stop
  END SUBROUTINE mp_initialize_matrix

  FUNCTION static_size_of_matrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = 2*SIZE_OF_INTEGER
    
  END FUNCTION static_size_of_matrix

  SUBROUTINE mp_initialize_vector(vec,n_rows)
    type(Matrix) :: vec
    integer :: n_rows

  END SUBROUTINE mp_initialize_vector

  SUBROUTINE mp_initialize_matrix_module(comm)
    INTEGER :: comm
  END SUBROUTINE mp_initialize_matrix_module

  LOGICAL FUNCTION mp_isInitialized(mat) 
    type(Matrix) :: mat
    mp_isInitialized=.false.
  END FUNCTION mp_isInitialized

  SUBROUTINE mp_allocate(mat)
    TYPE(Matrix), intent(INOUT) :: mat
  END SUBROUTINE mp_allocate

  SUBROUTINE mp_attach_vector(mat,localarray)
    TYPE(Matrix),INTENT(INOUT) :: mat
    COMPLEX, DIMENSION(:),INTENT(IN),TARGET :: localarray
  END SUBROUTINE mp_attach_vector

  SUBROUTINE mp_attach_matrix(mat,localarray)
    TYPE(Matrix),INTENT(INOUT) :: mat
    COMPLEX, DIMENSION(:,:),INTENT(IN),TARGET :: localarray
  END SUBROUTINE mp_attach_matrix

  SUBROUTINE mp_finalize_matrix(mat)
    type(Matrix) :: mat
  END SUBROUTINE mp_finalize_matrix

  SUBROUTINE mp_finalize_module()
  END SUBROUTINE mp_finalize_module

  SUBROUTINE mp_set_complex_value(mat,irow,icol,value)
    TYPE(Matrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value
  END SUBROUTINE mp_set_complex_value

  SUBROUTINE mp_set_real_value(mat,irow,icol,value)
    TYPE(Matrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value
  END SUBROUTINE mp_set_real_value

  SUBROUTINE mp_add_complex_value(mat,irow,icol,value)
    TYPE(Matrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value
  END SUBROUTINE mp_add_complex_value

  SUBROUTINE mp_add_real_value(mat,irow,icol,value)
    TYPE(Matrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value
  END SUBROUTINE mp_add_real_value


  SUBROUTINE mp_commit_values(mat)
    TYPE(Matrix) :: mat
  END SUBROUTINE mp_commit_values

  SUBROUTINE mp_dot_multiply(mat,vec, res, transA_in, transB_in)
    TYPE(Matrix),intent(IN) :: mat
    TYPE(Matrix),intent(IN) :: vec
    TYPE(Matrix),intent(INOUT) :: res
    CHARACTER, intent(IN),optional :: transA_in, transB_in
  END SUBROUTINE mp_dot_multiply

  SUBROUTINE mp_set_zero(mat)
    type(Matrix) :: mat
  END SUBROUTINE mp_set_zero

  SUBROUTINE mp_get_global_matrix_locally(mat,localfullmat)
    TYPE(Matrix),INTENT(IN) :: mat
    COMPLEX, DIMENSION(:,:) :: localfullmat
  END SUBROUTINE mp_get_global_matrix_locally

  FUNCTION mp_get_value(mat,irow,icol) result(value)
    type(Matrix) :: mat
    INTEGER :: irow, icol
    complex :: value
    value = cmplx(1,0)
  END FUNCTION mp_get_value

  FUNCTION mp_get_local_abs_square_sum(mat) RESULT(value)
    type(Matrix) :: mat
    real :: value
    value = 1.0
  END FUNCTION mp_get_local_abs_square_sum

  SUBROUTINE mp_show_on_screen(mat)
    type(Matrix) :: mat
  END SUBROUTINE mp_show_on_screen

  SUBROUTINE mp_show_in_file(mat,filename)
    type(Matrix) :: mat
    CHARACTER(len=*) :: filename
  END SUBROUTINE mp_show_in_file

  SUBROUTINE mp_output_data_matrix(basename,mat)
    type(Matrix) :: mat
    CHARACTER(len=*) :: basename
  END SUBROUTINE mp_output_data_matrix

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix(this,mat)
    TYPE(Matrix), intent(INOUT) :: this
    TYPE(Matrix),INTENT(IN) :: mat

  END SUBROUTINE mp_add_to_matrix

  SUBROUTINE mp_add_matrix(this,mat,res)
    TYPE(Matrix), INTENT(IN) :: this,mat
    TYPE(Matrix), INTENT(INOUT) :: res
    !PRINT*,"in mp_add_matrix of MatrixModule"

  END SUBROUTINE mp_add_matrix

  SUBROUTINE mp_subtract_matrix(this,mat,res)
    TYPE(Matrix),INTENT(IN) :: this,mat
    TYPE(Matrix),intent(INOUT) :: res
  END SUBROUTINE mp_subtract_matrix

  SUBROUTINE mp_subtract_from_matrix(this,mat)
    TYPE(Matrix),INTENT(INOUT) :: this
    TYPE(Matrix), intent(IN) :: mat
  END SUBROUTINE mp_subtract_from_matrix
    
  SUBROUTINE mp_multiply_matrix_with_real(this, scalar, res)
    TYPE(Matrix), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(Matrix),intent(INOUT) :: res
  END SUBROUTINE mp_multiply_matrix_with_real

  SUBROUTINE mp_scale_matrix_by_real(this, scalar)
    TYPE(Matrix), intent(INOUT) :: this
    REAL, intent(IN) :: scalar
  END SUBROUTINE mp_scale_matrix_by_real

  SUBROUTINE mp_assign_matrix(lmat,rmat)
    TYPE(Matrix), intent(INOUT) :: lmat
    TYPE(Matrix), intent(IN) :: rmat
  END SUBROUTINE mp_assign_matrix

  ! ====================================================
  ! == Define some mathematical routines for matrices ==
  ! ====================================================
  SUBROUTINE mp_invert_matrix(mat)
    TYPE(Matrix) :: mat
  END SUBROUTINE mp_invert_matrix

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_matrix_sum_0D(localsum,communicator)
    TYPE(Matrix), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_0D

  SUBROUTINE my_matrix_sum_2D(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(Matrix), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_2D

  SUBROUTINE my_matrix_sum_3D(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(Matrix), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_3D
  
  SUBROUTINE my_matrix_sum_generic(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(Matrix), intent(INOUT) :: localsum(len)
    integer :: communicator
  END SUBROUTINE my_matrix_sum_generic

  SUBROUTINE blacs_exit(val)
    integer:: val
  END SUBROUTINE blacs_exit

  SUBROUTINE mp_LU_factor(mat)
    TYPE(Matrix) :: mat
  END SUBROUTINE mp_LU_factor

  SUBROUTINE mp_LU_solve(this, mat, resmat, phase)
    TYPE(Matrix),intent(INOUT) :: this
    TYPE(Matrix),intent(IN) :: mat
    TYPE(Matrix),intent(INOUT) :: resmat
    INTEGER,intent(IN) :: phase
  END SUBROUTINE mp_LU_solve


!===================================== Repeat for REAL ==========================

  SUBROUTINE mp_initialize_matrix_real(mat,n_rows,n_cols)
    TYPE(MatrixReal), INTENT(INOUT) :: mat
    INTEGER :: n_rows, n_cols

    PRINT*, "You cannot run the global version of gene without the SCALAPACK switch in the makefile switched on."
    stop
  END SUBROUTINE mp_initialize_matrix_real

  SUBROUTINE mp_initialize_vector_real(vec,n_rows)
    type(MatrixReal) :: vec
    integer :: n_rows

  END SUBROUTINE mp_initialize_vector_real

  LOGICAL FUNCTION mp_isInitialized_real(mat)
    type(MatrixReal) :: mat
    mp_isInitialized_real=.false.
  END FUNCTION mp_isInitialized_real

  SUBROUTINE mp_allocate_real(mat)
    TYPE(MatrixReal), intent(INOUT) :: mat
  END SUBROUTINE mp_allocate_real

  SUBROUTINE mp_attach_vector_real(mat,localarray)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    REAL, DIMENSION(:),INTENT(IN),TARGET :: localarray
  END SUBROUTINE mp_attach_vector_real

  SUBROUTINE mp_attach_matrix_real(mat,localarray)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    REAL, DIMENSION(:,:),INTENT(IN),TARGET :: localarray
  END SUBROUTINE mp_attach_matrix_real

  SUBROUTINE mp_finalize_matrix_real(mat)
    type(MatrixReal) :: mat
  END SUBROUTINE mp_finalize_matrix_real

  SUBROUTINE mp_set_real_value_real(mat,irow,icol,value)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value
  END SUBROUTINE mp_set_real_value_real

  SUBROUTINE mp_add_real_value_real(mat,irow,icol,value)
    TYPE(MatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value
  END SUBROUTINE mp_add_real_value_real


  SUBROUTINE mp_commit_values_real(mat)
    TYPE(MatrixReal) :: mat
  END SUBROUTINE mp_commit_values_real

  SUBROUTINE mp_dot_multiply_real(mat,vec, res, transA_in, transB_in)
    TYPE(MatrixReal),intent(IN) :: mat
    TYPE(MatrixReal),intent(IN) :: vec
    TYPE(MatrixReal),intent(INOUT) :: res
    CHARACTER, intent(IN),optional :: transA_in, transB_in
  END SUBROUTINE mp_dot_multiply_real

  SUBROUTINE mp_set_zero_real(mat)
    type(MatrixReal) :: mat
  END SUBROUTINE mp_set_zero_real

  SUBROUTINE mp_get_global_matrix_locally_real(mat,localfullmat)
    TYPE(MatrixReal),INTENT(IN) :: mat
    REAL, DIMENSION(:,:) :: localfullmat
  END SUBROUTINE mp_get_global_matrix_locally_real

  FUNCTION mp_get_value_real(mat,irow,icol) result(value)
    type(MatrixReal) :: mat
    INTEGER :: irow, icol
    real :: value
    value = 1.0
  END FUNCTION mp_get_value_real

  FUNCTION mp_get_local_abs_square_sum_real(mat) RESULT(value)
    type(MatrixReal) :: mat
    real :: value
    value = 1.0
  END FUNCTION mp_get_local_abs_square_sum_real

  SUBROUTINE mp_show_on_screen_real(mat)
    type(MatrixReal) :: mat
  END SUBROUTINE mp_show_on_screen_real

  SUBROUTINE mp_show_in_file_real(mat,filename)
    type(MatrixReal) :: mat
    CHARACTER(len=*) :: filename
  END SUBROUTINE mp_show_in_file_real

  SUBROUTINE mp_output_data_matrix_real(basename,mat)
    type(MatrixReal) :: mat
    CHARACTER(len=*) :: basename
  END SUBROUTINE mp_output_data_matrix_real

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix_real(this,mat)
    TYPE(MatrixReal), intent(INOUT) :: this
    TYPE(MatrixReal),INTENT(IN) :: mat

  END SUBROUTINE mp_add_to_matrix_real

  SUBROUTINE mp_add_matrix_real(this,mat,res)
    TYPE(MatrixReal), INTENT(IN) :: this,mat
    TYPE(MatrixReal), INTENT(INOUT) :: res
    !PRINT*,"in mp_add_matrix of MatrixModule"

  END SUBROUTINE mp_add_matrix_real

  SUBROUTINE mp_subtract_matrix_real(this,mat,res)
    TYPE(MatrixReal),INTENT(IN) :: this,mat
    TYPE(MatrixReal),intent(INOUT) :: res
  END SUBROUTINE mp_subtract_matrix_real

  SUBROUTINE mp_subtract_from_matrix_real(this,mat)
    TYPE(MatrixReal),INTENT(INOUT) :: this
    TYPE(MatrixReal), intent(IN) :: mat
  END SUBROUTINE mp_subtract_from_matrix_real

  SUBROUTINE mp_multiply_matrix_with_real_real(this, scalar, res)
    TYPE(MatrixReal), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(MatrixReal),intent(INOUT) :: res
  END SUBROUTINE mp_multiply_matrix_with_real_real

  SUBROUTINE mp_scale_matrix_by_real_real(this, scalar)
    TYPE(MatrixReal), intent(INOUT) :: this
    REAL, intent(IN) :: scalar
  END SUBROUTINE mp_scale_matrix_by_real_real

  SUBROUTINE mp_assign_matrix_real(lmat,rmat)
    TYPE(MatrixReal), intent(INOUT) :: lmat
    TYPE(MatrixReal), intent(IN) :: rmat
  END SUBROUTINE mp_assign_matrix_real

  ! ====================================================
  ! == Define some mathematical routines for matrices ==
  ! ====================================================
  SUBROUTINE mp_invert_matrix_real(mat)
    TYPE(MatrixReal) :: mat
  END SUBROUTINE mp_invert_matrix_real

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_matrix_sum_0D_real(localsum,communicator)
    TYPE(MatrixReal), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_0D_real

  SUBROUTINE my_matrix_sum_2D_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(MatrixReal), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_2D_real

  SUBROUTINE my_matrix_sum_3D_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(MatrixReal), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_3D_real

  SUBROUTINE my_matrix_sum_generic_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(MatrixReal), intent(INOUT) :: localsum(len)
    integer :: communicator
  END SUBROUTINE my_matrix_sum_generic_real

  SUBROUTINE mp_LU_factor_real(mat)
    TYPE(MatrixReal) :: mat
  END SUBROUTINE mp_LU_factor_real

  SUBROUTINE mp_LU_solve_real(this, mat, resmat, phase)
    TYPE(MatrixReal),intent(INOUT) :: this
    TYPE(MatrixReal),intent(IN) :: mat
    TYPE(MatrixReal),intent(INOUT) :: resmat
    INTEGER,intent(IN) :: phase
  END SUBROUTINE mp_LU_solve_real

END MODULE MatrixModule

!=========================== VectorModule dummy ===============================

MODULE VectorModule
  implicit none

  private
  PUBLIC :: Vector,initialize,attach, finalize, isInitialized, &
       & initialize_vector_module, finalize_vector_module, &
       & static_size_of_vector, VectorReal

  TYPE Vector
     INTEGER :: NRows
  END TYPE Vector

  TYPE VectorReal
     INTEGER :: NRows
  END TYPE VectorReal       

  INTERFACE initialize
     MODULE PROCEDURE mp_initialize_vector, mp_initialize_vector_real
  END INTERFACE

  INTERFACE attach
     MODULE PROCEDURE mp_attach_vector, mp_attach_vector_real
  END INTERFACE

  INTERFACE finalize
     MODULE PROCEDURE mp_finalize_vector, mp_finalize_vector_real
  END INTERFACE

  INTERFACE isInitialized
     module procedure mp_isInitialized, mp_isInitialized_real
  END INTERFACE

contains
  FUNCTION static_size_of_vector() RESULT(memory_need)
    integer :: memory_need

    memory_need = SIZE_OF_INTEGER
    
  END FUNCTION static_size_of_vector
    
  SUBROUTINE initialize_vector_module
  END SUBROUTINE initialize_vector_module
  
  SUBROUTINE finalize_vector_module
  END SUBROUTINE finalize_vector_module

  SUBROUTINE mp_initialize_vector(vec,n_rows)
    TYPE(Vector), INTENT(INOUT) :: vec
    INTEGER :: n_rows

  END SUBROUTINE mp_initialize_vector

  SUBROUTINE mp_finalize_vector(vec)
    type(Vector) :: vec
  END SUBROUTINE mp_finalize_vector

  SUBROUTINE mp_attach_vector(vec,localarray)
    TYPE(Vector),INTENT(INOUT) :: vec
    COMPLEX, DIMENSION(:),INTENT(IN),TARGET :: localarray
  END SUBROUTINE mp_attach_vector

  LOGICAL FUNCTION mp_isInitialized(vec) 
    type(Vector) :: vec

    mp_isInitialized = .false.
  END FUNCTION mp_isInitialized

!=========================== Repeat for REAL ================================

  SUBROUTINE mp_initialize_vector_real(vec,n_rows)
    TYPE(VectorReal), INTENT(INOUT) :: vec
    INTEGER :: n_rows

  END SUBROUTINE mp_initialize_vector_real

  SUBROUTINE mp_finalize_vector_real(vec)
    type(VectorReal) :: vec
  END SUBROUTINE mp_finalize_vector_real

  SUBROUTINE mp_attach_vector_real(vec,localarray)
    TYPE(VectorReal),INTENT(INOUT) :: vec
    REAL, DIMENSION(:),INTENT(IN),TARGET :: localarray
  END SUBROUTINE mp_attach_vector_real

  LOGICAL FUNCTION mp_isInitialized_real(vec)
    type(VectorReal) :: vec

    mp_isInitialized_real = .false.
  END FUNCTION mp_isInitialized_real

END MODULE VectorModule

!=========================== BandedMatrixModule dummy ========================

MODULE BandedMatrixModule
  use MatrixModule
  use VectorModule
  implicit none

  PRIVATE
  PUBLIC :: BandedMatrix, BandedMatrixReal, initialize,initialize_BandedMatrix_module,&
       & ALLOCATE,finalize,set_value, commit_values, &
       & set_zero, dot_multiply, my_sum, &
       & add_matrix, isInitialized, get_local_abs_square_sum, &
       & get_global_matrix_locally, &
       & static_size_of_bandedmatrix, mat_get_value, add_value, &
       & get_number_of_bands, show, convert_Banded_to_Full,&
       & convert_Full_to_Banded, LU_solve, row_axpy,&
       & transpose_and_conjugate, transpose_storage

  TYPE BandedMatrix
     INTEGER :: NCols
     INTEGER :: NRows
  END TYPE BandedMatrix

  TYPE BandedMatrixReal
     INTEGER :: NCols
     INTEGER :: Nrows
  END TYPE BandedMatrixReal   

  INTERFACE initialize
     MODULE PROCEDURE mp_initialize_matrix, mp_initialize_vector
     MODULE PROCEDURE mp_initialize_matrix_real, mp_initialize_vector_real
  END INTERFACE

  INTERFACE initialize_BandedMatrix_module
     MODULE PROCEDURE mp_initialize_BandedMatrix_module
  END INTERFACE

  INTERFACE allocate
     MODULE PROCEDURE mp_allocate, mp_allocate_real
  END INTERFACE

  INTERFACE finalize
     MODULE PROCEDURE mp_finalize_matrix, mp_finalize_matrix_real
  END INTERFACE

  INTERFACE set_value
     MODULE PROCEDURE mp_set_complex_value, mp_set_real_value
     MODULE PROCEDURE mp_set_real_value_real
  END INTERFACE

  INTERFACE commit_values
     module procedure mp_commit_values
     module procedure mp_commit_values_real
  END INTERFACE

  INTERFACE set_zero
     MODULE PROCEDURE mp_set_zero
     MODULE PROCEDURE mp_set_zero_real
  END INTERFACE

  INTERFACE dot_multiply
     MODULE PROCEDURE mp_dot_multiply, bm_matmat_dot_multiply
     MODULE PROCEDURE bm_dot_multiply_Banded_with_Full
     MODULE PROCEDURE mp_dot_multiply_real, bm_matmat_dot_multiply_real
     MODULE PROCEDURE bm_dot_multiply_Banded_with_Full_real
  END INTERFACE

  INTERFACE add_matrix
    MODULE PROCEDURE mp_add_to_matrix, mp_add_matrix!, mp_add_Matrix_to_BandedMatrix
    MODULE PROCEDURE mp_add_to_matrix_real, mp_add_matrix_real!, mp_add_Matrix_to_BandedMatrix_real
  END INTERFACE

  INTERFACE row_axpy
     module procedure bm_row_axpy, bm_row_axpy_real
  END INTERFACE

  INTERFACE add_value
     MODULE PROCEDURE mp_add_complex_value, mp_add_real_value
     MODULE PROCEDURE mp_add_real_value_real
  END INTERFACE

  INTERFACE isInitialized
     module procedure mp_isInitialized, mp_isInitialized_real
  END INTERFACE

  INTERFACE get_local_abs_square_sum
     module procedure mp_get_local_abs_square_sum, mp_get_local_abs_square_sum_real
  END INTERFACE

  INTERFACE convert_Banded_to_Full
     module procedure mp_convert_Banded_to_Full
     module procedure mp_convert_Banded_to_Full_real
  END INTERFACE

  INTERFACE convert_Full_to_Banded
     module procedure mp_convert_Full_to_Banded
     module procedure mp_convert_Full_to_Banded_real
  END INTERFACE

  INTERFACE get_number_of_bands
     module procedure bm_get_number_of_bands, bm_get_number_of_bands_real
  END INTERFACE

  INTERFACE my_sum
     MODULE PROCEDURE my_matrix_sum_generic, my_matrix_sum_2D, &
          &my_matrix_sum_3D, my_matrix_sum_0D
     MODULE PROCEDURE my_matrix_sum_generic_real, my_matrix_sum_2D_real, &
          &my_matrix_sum_3D_real, my_matrix_sum_0D_real
  END INTERFACE

  INTERFACE mat_get_value
     module procedure mp_mat_get_value, mp_mat_get_value_real
  END INTERFACE

  INTERFACE show
     MODULE PROCEDURE mp_show_on_screen, mp_show_in_file
     MODULE PROCEDURE mp_show_on_screen_real, mp_show_in_file_real
  END INTERFACE

  INTERFACE LU_solve
     module procedure bm_LU_solve, bm_LU_solve_real
  END INTERFACE

  INTERFACE get_global_matrix_locally
     module procedure mp_get_global_matrix_locally, mp_get_global_matrix_locally_real
  END INTERFACE

  INTERFACE transpose_and_conjugate
     module procedure bm_transpose_and_conjugate
     module procedure bm_transpose_and_conjugate_real
  END INTERFACE

  INTERFACE transpose_storage
     module procedure bm_transpose_storage, bm_transpose_storage_real
  END INTERFACE

CONTAINS
  SUBROUTINE mp_initialize_matrix(mat,n_rows,n_cols,transp)
    Type(BandedMatrix), INTENT(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    LOGICAL, OPTIONAL :: transp
  END SUBROUTINE mp_initialize_matrix

  SUBROUTINE mp_initialize_vector(vec,n_rows)
    type(BandedMatrix) :: vec
    integer :: n_rows

  END SUBROUTINE mp_initialize_vector

  SUBROUTINE mp_initialize_BandedMatrix_module
  END SUBROUTINE mp_initialize_BandedMatrix_module

  LOGICAL FUNCTION mp_isInitialized(mat) 
    type(BandedMatrix) :: mat

    mp_isInitialized = .FALSE.
  END FUNCTION mp_isInitialized

  FUNCTION static_size_of_bandedmatrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = 2*SIZE_OF_INTEGER !???
    
  END FUNCTION static_size_of_bandedmatrix

  SUBROUTINE mp_allocate(mat)
    Type(BandedMatrix), intent(INOUT) :: mat

  END SUBROUTINE mp_allocate

  SUBROUTINE mp_finalize_matrix(mat)
    type(BandedMatrix) :: mat
  END SUBROUTINE mp_finalize_matrix

  SUBROUTINE mp_set_complex_value(mat,irow,icol,value)
    Type(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

  END SUBROUTINE mp_set_complex_value

  SUBROUTINE mp_set_real_value(mat,irow,icol,value)
    Type(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

  END SUBROUTINE mp_set_real_value

  SUBROUTINE mp_add_complex_value(mat,irow,icol,value)
    Type(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

  END SUBROUTINE mp_add_complex_value

  SUBROUTINE mp_add_real_value(mat,irow,icol,value)
    Type(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

  END SUBROUTINE mp_add_real_value

  SUBROUTINE mp_commit_values(mat)
    Type(BandedMatrix) :: mat
    
  END SUBROUTINE mp_commit_values

  SUBROUTINE mp_dot_multiply(mat,vec, res)
    Type(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res

  END SUBROUTINE mp_dot_multiply

  SUBROUTINE bm_matmat_dot_multiply(mat,mat2,res,transA_in, transB_in)
    TYPE(BandedMatrix),INTENT(IN) :: mat, mat2
    TYPE(BandedMatrix),INTENT(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in
    CHARACTER, intent(IN), optional :: transB_in
  END SUBROUTINE bm_matmat_dot_multiply
subroutine bm_dot_multiply_Banded_with_Full(bmat,mat,res,transA_in,transB_in)
    type(BandedMatrix),intent(IN) :: bmat
    type(Matrix), intent(IN) :: mat
    type(Matrix),intent(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in, transB_in
  
end subroutine bm_dot_multiply_Banded_with_Full

  SUBROUTINE mp_set_zero(mat)
    type(BandedMatrix) :: mat

  END SUBROUTINE mp_set_zero

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix(this,mat)
    Type(BandedMatrix), intent(INOUT) :: this
    Type(BandedMatrix),INTENT(IN) :: mat
  END SUBROUTINE mp_add_to_matrix

  SUBROUTINE mp_add_matrix(this,mat,res)
    Type(BandedMatrix), INTENT(IN) :: this,mat
    Type(BandedMatrix), INTENT(INOUT) :: res
  END SUBROUTINE mp_add_matrix

  FUNCTION mp_get_local_abs_square_sum(mat) RESULT(value)
    type(BandedMatrix) :: mat
    real :: value

    value = 0.0
  END FUNCTION mp_get_local_abs_square_sum

  FUNCTION mp_mat_get_value(mat,irow,icol) result(value)
    type(BandedMatrix) :: mat
    INTEGER :: irow, icol
    complex :: value
    value = cmplx(1,0)
  END FUNCTION mp_mat_get_value

  FUNCTION bm_get_number_of_bands(mat)
    type(BandedMatrix) :: mat
    INTEGER :: bm_get_number_of_bands

    bm_get_number_of_bands=1 
  END FUNCTION bm_get_number_of_bands

  SUBROUTINE mp_get_global_matrix_locally(bmat,localfullmat)
    TYPE(BandedMatrix),INTENT(IN) :: bmat
    COMPLEX, DIMENSION(:,:) :: localfullmat
  END SUBROUTINE mp_get_global_matrix_locally


  SUBROUTINE mp_convert_Banded_to_Full(bmat, fmat)
    TYPE(BandedMatrix), intent(IN) :: bmat
    TYPE(Matrix), intent(INOUT) :: fmat

  END SUBROUTINE mp_convert_Banded_to_Full

  SUBROUTINE mp_convert_Full_to_Banded(fmat, bmat)
    TYPE(Matrix), intent(IN) :: fmat
    TYPE(BandedMatrix), intent(INOUT) :: bmat

  END SUBROUTINE mp_convert_Full_to_Banded

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_matrix_sum_0D(localsum,communicator)
    TYPE(BandedMatrix), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_0D

  SUBROUTINE my_matrix_sum_2D(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(BandedMatrix), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_2D

  SUBROUTINE my_matrix_sum_3D(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(BandedMatrix), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_3D
  
  SUBROUTINE my_matrix_sum_generic(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(BandedMatrix), intent(INOUT) :: localsum(len)
    integer :: communicator
  END SUBROUTINE my_matrix_sum_generic

  SUBROUTINE mp_show_on_screen(mat)
    type(BandedMatrix) :: mat
    
  END SUBROUTINE mp_show_on_screen

  SUBROUTINE mp_show_in_file(mat,filename)
    type(BandedMatrix) :: mat
    CHARACTER(len=*) :: filename

  END SUBROUTINE mp_show_in_file

  SUBROUTINE bm_LU_solve(mat,vec,res)
    TYPE(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res

  END SUBROUTINE bm_LU_solve

  SUBROUTINE bm_row_axpy(this,this_row,mat,scalar)
    TYPE(BandedMatrix),intent(INOUT) :: this
    INTEGER, INTENT(IN) :: this_row
    TYPE(BandedMatrix), intent(IN) :: mat
    COMPLEX, intent(IN) :: scalar
  END SUBROUTINE bm_row_axpy

  subroutine bm_transpose_and_conjugate(mat,dagger_mat)
    type(BandedMatrix) :: mat,dagger_mat
  end subroutine bm_transpose_and_conjugate

  subroutine bm_transpose_storage(mat,mat_t)
    type(BandedMatrix) :: mat,mat_t
  end subroutine bm_transpose_storage

! ========================== Repeat for REAL ================================

  SUBROUTINE mp_initialize_matrix_real(mat,n_rows,n_cols,transp)
    Type(BandedMatrixReal), INTENT(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    LOGICAL, OPTIONAL :: transp
  END SUBROUTINE mp_initialize_matrix_real

  SUBROUTINE mp_initialize_vector_real(vec,n_rows)
    type(BandedMatrixReal) :: vec
    integer :: n_rows

  END SUBROUTINE mp_initialize_vector_real

  LOGICAL FUNCTION mp_isInitialized_real(mat)
    type(BandedMatrixReal) :: mat

    mp_isInitialized_real = .FALSE.
  END FUNCTION mp_isInitialized_real

  SUBROUTINE mp_allocate_real(mat)
    Type(BandedMatrixReal), intent(INOUT) :: mat

  END SUBROUTINE mp_allocate_real

  SUBROUTINE mp_finalize_matrix_real(mat)
    type(BandedMatrixReal) :: mat
  END SUBROUTINE mp_finalize_matrix_real

  SUBROUTINE mp_set_real_value_real(mat,irow,icol,value)
    Type(BandedMatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

  END SUBROUTINE mp_set_real_value_real

  SUBROUTINE mp_add_real_value_real(mat,irow,icol,value)
    Type(BandedMatrixReal),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

  END SUBROUTINE mp_add_real_value_real

  SUBROUTINE mp_commit_values_real(mat)
    Type(BandedMatrixReal) :: mat

  END SUBROUTINE mp_commit_values_real

  SUBROUTINE mp_dot_multiply_real(mat,vec, res)
    Type(BandedMatrixReal),intent(IN) :: mat
    Type(VectorReal),intent(IN) :: vec
    Type(VectorReal),intent(INOUT) :: res

  END SUBROUTINE mp_dot_multiply_real

  SUBROUTINE bm_matmat_dot_multiply_real(mat,mat2,res,transA_in, transB_in)
    TYPE(BandedMatrixReal),INTENT(IN) :: mat, mat2
    TYPE(BandedMatrixReal),INTENT(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in
    CHARACTER, intent(IN), optional :: transB_in
  END SUBROUTINE bm_matmat_dot_multiply_real
subroutine bm_dot_multiply_Banded_with_Full_real(bmat,mat,res,transA_in,transB_in)
    type(BandedMatrixReal),intent(IN) :: bmat
    type(MatrixReal), intent(IN) :: mat
    type(MatrixReal),intent(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in, transB_in

end subroutine bm_dot_multiply_Banded_with_Full_real

  SUBROUTINE mp_set_zero_real(mat)
    type(BandedMatrixReal) :: mat

  END SUBROUTINE mp_set_zero_real

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix_real(this,mat)
    Type(BandedMatrixReal), intent(INOUT) :: this
    Type(BandedMatrixReal),INTENT(IN) :: mat
  END SUBROUTINE mp_add_to_matrix_real

  SUBROUTINE mp_add_matrix_real(this,mat,res)
    Type(BandedMatrixReal), INTENT(IN) :: this,mat
    Type(BandedMatrixReal), INTENT(INOUT) :: res
  END SUBROUTINE mp_add_matrix_real

  FUNCTION mp_get_local_abs_square_sum_real(mat) RESULT(value)
    type(BandedMatrixReal) :: mat
    real :: value

    value = 0.0
  END FUNCTION mp_get_local_abs_square_sum_real

  FUNCTION mp_mat_get_value_real(mat,irow,icol) result(value)
    type(BandedMatrixReal) :: mat
    INTEGER :: irow, icol
    real :: value
    value = 1.0
  END FUNCTION mp_mat_get_value_real

  FUNCTION bm_get_number_of_bands_real(mat)
    type(BandedMatrixReal) :: mat
    INTEGER :: bm_get_number_of_bands_real

    bm_get_number_of_bands_real=1
  END FUNCTION bm_get_number_of_bands_real

  SUBROUTINE mp_get_global_matrix_locally_real(bmat,localfullmat)
    TYPE(BandedMatrixReal),INTENT(IN) :: bmat
    REAL, DIMENSION(:,:) :: localfullmat
  END SUBROUTINE mp_get_global_matrix_locally_real


  SUBROUTINE mp_convert_Banded_to_Full_real(bmat, fmat)
    TYPE(BandedMatrixReal), intent(IN) :: bmat
    TYPE(MatrixReal), intent(INOUT) :: fmat

  END SUBROUTINE mp_convert_Banded_to_Full_real

  SUBROUTINE mp_convert_Full_to_Banded_real(fmat, bmat)
    TYPE(MatrixReal), intent(IN) :: fmat
    TYPE(BandedMatrixReal), intent(INOUT) :: bmat

  END SUBROUTINE mp_convert_Full_to_Banded_real

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_matrix_sum_0D_real(localsum,communicator)
    TYPE(BandedMatrixReal), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_0D_real

  SUBROUTINE my_matrix_sum_2D_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(BandedMatrixReal), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_2D_real

  SUBROUTINE my_matrix_sum_3D_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(BandedMatrixReal), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_3D_real

  SUBROUTINE my_matrix_sum_generic_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(BandedMatrixReal), intent(INOUT) :: localsum(len)
    integer :: communicator
  END SUBROUTINE my_matrix_sum_generic_real

  SUBROUTINE mp_show_on_screen_real(mat)
    type(BandedMatrixReal) :: mat

  END SUBROUTINE mp_show_on_screen_real

  SUBROUTINE mp_show_in_file_real(mat,filename)
    type(BandedMatrixReal) :: mat
    CHARACTER(len=*) :: filename

  END SUBROUTINE mp_show_in_file_real

  SUBROUTINE bm_LU_solve_real(mat,vec,res)
    TYPE(BandedMatrixReal),intent(IN) :: mat
    Type(VectorReal),intent(IN) :: vec
    Type(VectorReal),intent(INOUT) :: res

  END SUBROUTINE bm_LU_solve_real

  SUBROUTINE bm_row_axpy_real(this,this_row,mat,scalar)
    TYPE(BandedMatrixReal),intent(INOUT) :: this
    INTEGER, INTENT(IN) :: this_row
    TYPE(BandedMatrixReal), intent(IN) :: mat
    REAL, intent(IN) :: scalar
  END SUBROUTINE bm_row_axpy_real

  subroutine bm_transpose_and_conjugate_real(mat,dagger_mat)
    type(BandedMatrixReal) :: mat,dagger_mat
  end subroutine bm_transpose_and_conjugate_real

  subroutine bm_transpose_storage_real(mat,mat_t)
    type(BandedMatrixReal) :: mat,mat_t
  end subroutine bm_transpose_storage_real

END MODULE BandedMatrixModule

! ======================= Derivative Matrix Module Dummy =====================
MODULE DerivativeMatrixModule
  USE Grid1DModule
  USE MatrixModule
  use BandedMatrixModule

  implicit none

  PRIVATE
  PUBLIC :: DerivativeMatrix, initialize,finalize, show, mat_get_value, &
       & static_size_of_DerivativeMatrix, calculate, calculate_real, DerivativeMatrixReal

  TYPE DerivativeMatrix
     type(BandedMatrix) :: Data

     integer :: derivative_order !the order of the finite differences formulas for the derivatives
     logical :: periodic_boundary
     !TYPE(Grid1D), pointer :: grid
  END TYPE DerivativeMatrix

  TYPE DerivativeMatrixReal
     type(BandedMatrixReal) :: Data
     
     integer :: derivative_order
     logical :: periodic_boundary
  END TYPE DerivativeMatrixReal   


  INTERFACE initialize
     module procedure mp_initialize_matrix, mp_initialize_matrix_real
  END INTERFACE

  INTERFACE finalize
     module procedure mp_finalize_matrix, mp_finalize_matrix_real
  END INTERFACE

  INTERFACE show
     MODULE PROCEDURE mp_show_on_screen, mp_show_in_file
     MODULE PROCEDURE mp_show_on_screen_real, mp_show_in_file_real
  END INTERFACE

  INTERFACE mat_get_value
     module procedure mp_mat_get_value, mp_mat_get_value_real
  END INTERFACE


CONTAINS

  FUNCTION static_size_of_DerivativeMatrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = 0
    
  END FUNCTION static_size_of_DerivativeMatrix
    
  SUBROUTINE mp_initialize_module
  END SUBROUTINE mp_initialize_module

  SUBROUTINE mp_initialize_matrix(dmat,grid,p_derivative_order,transposed)
    type(DerivativeMatrix) :: dmat !> derivative matrix
    TYPE(Grid1D),TARGET,intent(in) :: grid     !> grid definition
    INTEGER,intent(in) :: p_derivative_order   !> derivative order
    LOGICAL, OPTIONAL :: transposed
  END SUBROUTINE mp_initialize_matrix

  SUBROUTINE mp_finalize_matrix(dmat)
    TYPE(DerivativeMatrix) :: dmat
  END SUBROUTINE mp_finalize_matrix

  SUBROUTINE calculate(dmat, which_derivative, rad_bc_type)
    TYPE(DerivativeMatrix) :: dmat
    INTEGER, intent(IN) :: which_derivative
    INTEGER,intent(IN)  :: rad_bc_type 
  END SUBROUTINE calculate

  SUBROUTINE mp_show_on_screen(dmat)
    type(DerivativeMatrix) :: dmat
  END SUBROUTINE mp_show_on_screen

  SUBROUTINE mp_show_in_file(dmat,filename)
    type(DerivativeMatrix) :: dmat
    CHARACTER(len=*) :: filename
  END SUBROUTINE mp_show_in_file

  FUNCTION mp_mat_get_value(dmat,irow,icol) RESULT(value)
    type(DerivativeMatrix) :: dmat
    INTEGER :: irow,icol
    complex :: value
    value = CMPLX(1,0)
  END FUNCTION mp_mat_get_value

!============================= Repeat for REAL =================================

  SUBROUTINE mp_initialize_matrix_real(dmat,grid,p_derivative_order,transposed)
    type(DerivativeMatrixReal) :: dmat !> derivative matrix
    TYPE(Grid1D),TARGET,intent(in) :: grid     !> grid definition
    INTEGER,intent(in) :: p_derivative_order   !> derivative order
    LOGICAL, OPTIONAL :: transposed
  END SUBROUTINE mp_initialize_matrix_real

  SUBROUTINE mp_finalize_matrix_real(dmat)
    TYPE(DerivativeMatrixReal) :: dmat
  END SUBROUTINE mp_finalize_matrix_real

  SUBROUTINE calculate_real(dmat, which_derivative, rad_bc_type)
    TYPE(DerivativeMatrixReal) :: dmat
    INTEGER, intent(IN) :: which_derivative
    INTEGER,intent(IN)  :: rad_bc_type
  END SUBROUTINE calculate_real

  SUBROUTINE mp_show_on_screen_real(dmat)
    type(DerivativeMatrixReal) :: dmat
  END SUBROUTINE mp_show_on_screen_real

  SUBROUTINE mp_show_in_file_real(dmat,filename)
    type(DerivativeMatrixReal) :: dmat
    CHARACTER(len=*) :: filename
  END SUBROUTINE mp_show_in_file_real

  FUNCTION mp_mat_get_value_real(dmat,irow,icol) RESULT(value)
    type(DerivativeMatrixReal) :: dmat
    INTEGER :: irow,icol
    real :: value
    value = 1.0
  END FUNCTION mp_mat_get_value_real

END MODULE DerivativeMatrixModule
