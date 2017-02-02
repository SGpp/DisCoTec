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

MODULE VectorModule
  USE ProcessGridModule
  USE MatrixModule, only: get_processgrid
  USE StoreVectorModule
  !USE mpi
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Vector, VectorReal
  PUBLIC :: initialize, initialize_vector_module, finalize, finalize_vector_module
  PUBLIC :: allocate, attach, set_value, add_value, commit_values, set_zero
  PUBLIC :: get_global_vector_locally, get_local_abs_square_sum, vec_get_value
  PUBLIC :: show, add_vector, subtract_vector, multiply_vector_with_scalar
  PUBLIC :: isInitialized, my_sum, ASSIGNMENT(=)
  PUBLIC :: static_size_of_vector
  ! the following variables are used from ProcessGridModule. To use them outside the 
  ! VectorModule, we have to publicize them explicitly.
  PUBLIC :: MAT_BLOCKS_OF_ROWS, MAT_BLOCKS_OF_COLS, MAT_BOTH

  TYPE Vector
     INTEGER :: NRows
     TYPE(StoreVectorObject) :: DATA
  END TYPE Vector

  TYPE VectorReal
     INTEGER :: NRows
     TYPE(StoreVectorObjectReal) :: DATA
  END TYPE VectorReal

  INTERFACE initialize
     MODULE PROCEDURE mp_initialize_vector, mp_initialize_vector_real
  END INTERFACE

  INTERFACE initialize_vector_module
     MODULE PROCEDURE mp_initialize_vector_module
  END INTERFACE

  INTERFACE allocate
     MODULE PROCEDURE mp_allocate, mp_allocate_real
  END INTERFACE

  INTERFACE attach
     MODULE PROCEDURE mp_attach_vector, mp_attach_vector_real
  END INTERFACE

  INTERFACE finalize
     MODULE PROCEDURE mp_finalize_vector, mp_finalize_vector_real
  END INTERFACE

  INTERFACE finalize_vector_module
     MODULE PROCEDURE mp_finalize_vector_module
  END INTERFACE

  INTERFACE set_value
     MODULE PROCEDURE mp_set_complex_value, mp_set_real_value, mp_set_real_value_real
  END INTERFACE

  INTERFACE add_value
     MODULE PROCEDURE mp_add_complex_value, mp_add_real_value, mp_add_real_value_real
  END INTERFACE

  INTERFACE commit_values
     module procedure mp_commit_values, mp_commit_values_real
  END INTERFACE

  INTERFACE set_zero
     MODULE PROCEDURE mp_set_zero, mp_set_zero_real
  END INTERFACE

  INTERFACE get_global_vector_locally
     module procedure mp_get_global_vector_locally, mp_get_global_vector_locally_real
  END INTERFACE

  INTERFACE get_local_abs_square_sum
     module procedure mp_get_local_abs_square_sum, mp_get_local_abs_square_sum_real
  END INTERFACE

  INTERFACE vec_get_value
     MODULE PROCEDURE mp_get_value, mp_get_value_real
  END INTERFACE

  INTERFACE show
     MODULE PROCEDURE mp_show_on_screen, mp_show_in_file, mp_show_on_screen_real, mp_show_in_file_real
  END INTERFACE

  INTERFACE add_vector
     MODULE PROCEDURE mp_add_to_vector, mp_add_vector, mp_add_to_vector_real, mp_add_vector_real
  END INTERFACE

  interface subtract_vector
     MODULE PROCEDURE mp_subtract_vector, mp_subtract_from_vector, mp_subtract_vector_real, mp_subtract_from_vector_real
  END INTERFACE

  interface multiply_vector_with_scalar
     module procedure mp_multiply_vector_with_real, mp_multiply_vector_with_real_real
     module procedure mp_scale_vector_by_real, mp_scale_vector_by_real_real
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
     module procedure mp_assign_vector, mp_assign_vector_real
  END INTERFACE

  INTERFACE isInitialized
     module procedure mp_isInitialized, mp_isInitialized_real
  END INTERFACE

  INTERFACE my_sum
     MODULE PROCEDURE my_vector_sum_generic, my_vector_sum_2D, my_vector_sum_3D, my_vector_sum_0D
     MODULE PROCEDURE my_vector_sum_generic_real, my_vector_sum_2D_real, my_vector_sum_3D_real, my_vector_sum_0D_real
  END INTERFACE

  ! Use the same process grid for all Vectors instantiated whithin this module
  TYPE(ProcessGrid),SAVE :: PG

CONTAINS
  FUNCTION static_size_of_vector() RESULT(memory_need)
    integer :: memory_need

    memory_need = 1*SIZE_OF_INTEGER &
         & + static_size_of_storevectorobject()
    
  END FUNCTION static_size_of_vector
    
  SUBROUTINE mp_initialize_vector(vec,n_rows)
    TYPE(Vector), INTENT(INOUT) :: vec
    INTEGER :: n_rows

    DEBSTART("mp_initialize_vector")
    IF (PG%isInitialized) THEN
       vec%NRows = n_rows
       CALL initialize(vec%Data,n_rows, PG)
    ELSE
       PRINT*,"Before initializing a vector, the process grid has to be initialized."
       STOP
    END IF
    DEBEND("mp_initialize_vector")
  END SUBROUTINE mp_initialize_vector

  SUBROUTINE mp_initialize_vector_module

    DEBSTART("mp_initialize_vector_module")
    ! get the processgrid from a previously initialized matrix module
    call get_processgrid(PG)

    DEBEND("mp_initialize_vector_module")
  END SUBROUTINE mp_initialize_vector_module

  LOGICAL FUNCTION mp_isInitialized(vec) 
    type(Vector) :: vec

    mp_isInitialized = isInitialized(vec%Data)
  END FUNCTION mp_isInitialized

  SUBROUTINE mp_allocate(vec)
    TYPE(Vector), intent(INOUT) :: vec

    CALL allocate(vec%Data)
  END SUBROUTINE mp_allocate

  SUBROUTINE mp_attach_vector(vec,localarray)
    TYPE(Vector),INTENT(INOUT) :: vec
    COMPLEX, DIMENSION(:),INTENT(IN),TARGET :: localarray

    DEBSTART("mp_attach_vector")
    CALL attach(vec%Data,localarray)
    DEBEND("mp_attach_vector")
  END SUBROUTINE mp_attach_vector

  SUBROUTINE mp_finalize_vector(vec)
    type(Vector) :: vec
    call finalize(vec%Data)
  END SUBROUTINE mp_finalize_vector

  SUBROUTINE mp_finalize_vector_module()
    ! do nothing
  END SUBROUTINE mp_finalize_vector_module

  SUBROUTINE mp_set_complex_value(vec,irow,value)
    TYPE(Vector),INTENT(INOUT) :: vec
    INTEGER,INTENT(IN) :: irow
    COMPLEX, INTENT(IN) :: value

    !PRINT*,"In routine set_value of VectorModule."
    CALL set_value(vec%Data,irow,value)
  END SUBROUTINE mp_set_complex_value

  SUBROUTINE mp_set_real_value(vec,irow,value)
    TYPE(Vector),INTENT(INOUT) :: vec
    INTEGER,INTENT(IN) :: irow
    REAL, INTENT(IN) :: value

    !PRINT*,"In routine set_value of VectorModule."
    CALL set_value(vec%Data,irow,CMPLX(value,0.0))
  END SUBROUTINE mp_set_real_value

  SUBROUTINE mp_add_complex_value(vec,irow,value)
    TYPE(Vector),INTENT(INOUT) :: vec
    INTEGER,INTENT(IN) :: irow
    COMPLEX, INTENT(IN) :: value

    CALL add_value(vec%Data,irow,value)
  END SUBROUTINE mp_add_complex_value

  SUBROUTINE mp_add_real_value(vec,irow,value)
    TYPE(Vector),INTENT(INOUT) :: vec
    INTEGER,INTENT(IN) :: irow
    REAL, INTENT(IN) :: value

    CALL add_value(vec%Data,irow,CMPLX(value,0.0))
  END SUBROUTINE mp_add_real_value

  SUBROUTINE mp_commit_values(vec)
    TYPE(Vector) :: vec
    
    CALL commit_values(vec%Data)
  END SUBROUTINE mp_commit_values

  SUBROUTINE mp_set_zero(vec)
    type(Vector) :: vec

    ! Local variables
    INTEGER :: irow

    DEBSTART("mp_set_zero")
    DO irow=1,vec%NRows
       CALL set_value(vec,irow,CMPLX(0.0,0.0))
    END DO
    CALL commit_values(vec)
    DEBEND("mp_set_zero")
  END SUBROUTINE mp_set_zero

  SUBROUTINE mp_get_global_vector_locally(vec,localfullvec)
    TYPE(Vector),INTENT(IN) :: vec
    COMPLEX, DIMENSION(:),INTENT(OUT) :: localfullvec

    CALL get_global_vector_locally(vec%Data,localfullvec)
  END SUBROUTINE mp_get_global_vector_locally

  FUNCTION mp_get_value(vec,irow) result(value)
    type(Vector) :: vec
    INTEGER :: irow
    complex :: value

    value = vec_get_value(vec%Data,irow)
  END FUNCTION mp_get_value

  FUNCTION mp_get_local_abs_square_sum(vec) RESULT(value)
    type(Vector) :: vec
    real :: value

    value = get_local_abs_square_sum(vec%Data)
  END FUNCTION mp_get_local_abs_square_sum

  SUBROUTINE mp_show_on_screen(vec)
    type(Vector) :: vec
    
    CALL show(vec%Data)
  END SUBROUTINE mp_show_on_screen

  SUBROUTINE mp_show_in_file(vec,filename)
    type(Vector) :: vec
    CHARACTER(len=*) :: filename

    CALL show(vec%Data,filename)
  END SUBROUTINE mp_show_in_file

  ! =======================================
  ! == Define some operators on vecrices ==
  ! =======================================
  SUBROUTINE mp_add_to_vector(this,vec)
    TYPE(Vector), intent(INOUT) :: this
    TYPE(Vector),INTENT(IN) :: vec
    CALL add_vector(this%Data,vec%Data)
  END SUBROUTINE mp_add_to_vector

  SUBROUTINE mp_add_vector(this,vec,res)
    TYPE(Vector), INTENT(IN) :: this,vec
    TYPE(Vector), INTENT(INOUT) :: res
    !PRINT*,"in mp_add_vector of VectorModule"
    CALL add_vector(this%Data,vec%Data,res%Data)
  END SUBROUTINE mp_add_vector

  SUBROUTINE mp_subtract_vector(this,vec,res)
    TYPE(Vector),INTENT(IN) :: this,vec
    TYPE(Vector),intent(INOUT) :: res
    ! calculate res = this - vec

    CALL subtract_vector(this%Data,vec%Data, res%Data)
  END SUBROUTINE mp_subtract_vector

  SUBROUTINE mp_subtract_from_vector(this,vec)
    TYPE(Vector),INTENT(INOUT) :: this
    TYPE(Vector), intent(IN) :: vec
    ! calculate this = this - vec

    CALL subtract_vector(this%Data,vec%Data)
  END SUBROUTINE mp_subtract_from_vector
    
  SUBROUTINE mp_multiply_vector_with_real(this, scalar, res)
    TYPE(Vector), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(Vector),intent(INOUT) :: res

    CALL multiply_vector_with_scalar(this%Data,scalar, res%Data)
  END SUBROUTINE mp_multiply_vector_with_real

  SUBROUTINE mp_scale_vector_by_real(this, scalar)
    TYPE(Vector), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL multiply_vector_with_scalar(this%Data,scalar)
  END SUBROUTINE mp_scale_vector_by_real

  SUBROUTINE mp_assign_vector(lvec,rvec)
    TYPE(Vector), intent(INOUT) :: lvec
    TYPE(Vector), intent(IN) :: rvec

    IF (lvec%Nrows.EQ.rvec%Nrows) THEN
       lvec%Data = rvec%Data
    ELSE
       PRINT*, "Left and right vecrices in an assignment must have the same dimensions."
    END IF
  END SUBROUTINE mp_assign_vector

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_vector_sum_0D(localsum,communicator)
    TYPE(Vector), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_sum(localsum%Data,communicator)
  END SUBROUTINE my_vector_sum_0D

  SUBROUTINE my_vector_sum_2D(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(Vector), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_vector_sum_generic(localsum,len,communicator)
  END SUBROUTINE my_vector_sum_2D

  SUBROUTINE my_vector_sum_3D(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(Vector), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_vector_sum_generic(localsum,len,communicator)
  END SUBROUTINE my_vector_sum_3D
  
  SUBROUTINE my_vector_sum_generic(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(Vector), intent(INOUT) :: localsum(len)
    integer :: communicator

    integer :: iVector
    DO iVector=1,len
       CALL my_sum(localsum(iVector)%Data,communicator)
    END DO
  END SUBROUTINE my_vector_sum_generic

  ! ===========================================
  ! == Repeat EVERYTHING for REAL data types ==
  ! ===========================================

  SUBROUTINE mp_initialize_vector_real(vec,n_rows)
    TYPE(VectorReal), INTENT(INOUT) :: vec
    INTEGER :: n_rows

    DEBSTART("mp_initialize_vector_real")
    IF (PG%isInitialized) THEN
       vec%NRows = n_rows
       CALL initialize(vec%Data,n_rows, PG)
    ELSE
       PRINT*,"Before initializing a vector, the process grid has to be initialized."
       STOP
    END IF
    DEBEND("mp_initialize_vector_real")
  END SUBROUTINE mp_initialize_vector_real

  LOGICAL FUNCTION mp_isInitialized_real(vec) 
    type(VectorReal) :: vec

    mp_isInitialized_real = isInitialized(vec%Data)
  END FUNCTION mp_isInitialized_real

  SUBROUTINE mp_allocate_real(vec)
    TYPE(VectorReal), intent(INOUT) :: vec

    CALL allocate(vec%Data)
  END SUBROUTINE mp_allocate_real

  SUBROUTINE mp_attach_vector_real(vec,localarray)
    TYPE(VectorReal),INTENT(INOUT) :: vec
    REAL, DIMENSION(:),INTENT(IN),TARGET :: localarray

    DEBSTART("mp_attach_vector_real")
    CALL attach(vec%Data,localarray)
    DEBEND("mp_attach_vector_real")
  END SUBROUTINE mp_attach_vector_real

  SUBROUTINE mp_finalize_vector_real(vec)
    type(VectorReal) :: vec
    call finalize(vec%Data)
  END SUBROUTINE mp_finalize_vector_real

  SUBROUTINE mp_set_real_value_real(vec,irow,value)
    TYPE(VectorReal),INTENT(INOUT) :: vec
    INTEGER,INTENT(IN) :: irow
    REAL, INTENT(IN) :: value

    !PRINT*,"In routine set_value of VectorModule."
    CALL set_value(vec%Data,irow,value)
  END SUBROUTINE mp_set_real_value_real

  SUBROUTINE mp_add_real_value_real(vec,irow,value)
    TYPE(VectorReal),INTENT(INOUT) :: vec
    INTEGER,INTENT(IN) :: irow
    REAL, INTENT(IN) :: value

    CALL add_value(vec%Data,irow,value)
  END SUBROUTINE mp_add_real_value_real

  SUBROUTINE mp_commit_values_real(vec)
    TYPE(VectorReal) :: vec
    
    CALL commit_values(vec%Data)
  END SUBROUTINE mp_commit_values_real

  SUBROUTINE mp_set_zero_real(vec)
    type(VectorReal) :: vec

    ! Local variables
    INTEGER :: irow

    DEBSTART("mp_set_zero_real")
    DO irow=1,vec%NRows
       CALL set_value(vec,irow,0.0)
    END DO
    CALL commit_values(vec)
    DEBEND("mp_set_zero_real")
  END SUBROUTINE mp_set_zero_real

  SUBROUTINE mp_get_global_vector_locally_real(vec,localfullvec)
    TYPE(VectorReal),INTENT(IN) :: vec
    REAL, DIMENSION(:),INTENT(OUT) :: localfullvec

    CALL get_global_vector_locally(vec%Data,localfullvec)
  END SUBROUTINE mp_get_global_vector_locally_real

  FUNCTION mp_get_value_real(vec,irow) result(value)
    type(VectorReal) :: vec
    INTEGER :: irow
    real :: value

    value = vec_get_value(vec%Data,irow)
  END FUNCTION mp_get_value_real

  FUNCTION mp_get_local_abs_square_sum_real(vec) RESULT(value)
    type(VectorReal) :: vec
    real :: value

    value = get_local_abs_square_sum(vec%Data)
  END FUNCTION mp_get_local_abs_square_sum_real

  SUBROUTINE mp_show_on_screen_real(vec)
    type(VectorReal) :: vec
    
    CALL show(vec%Data)
  END SUBROUTINE mp_show_on_screen_real

  SUBROUTINE mp_show_in_file_real(vec,filename)
    type(VectorReal) :: vec
    CHARACTER(len=*) :: filename

    CALL show(vec%Data,filename)
  END SUBROUTINE mp_show_in_file_real

  ! =======================================
  ! == Define some operators on vecrices ==
  ! =======================================
  SUBROUTINE mp_add_to_vector_real(this,vec)
    TYPE(VectorReal), intent(INOUT) :: this
    TYPE(VectorReal),INTENT(IN) :: vec
    CALL add_vector(this%Data,vec%Data)
  END SUBROUTINE mp_add_to_vector_real

  SUBROUTINE mp_add_vector_real(this,vec,res)
    TYPE(VectorReal), INTENT(IN) :: this,vec
    TYPE(VectorReal), INTENT(INOUT) :: res
    !PRINT*,"in mp_add_vector of VectorModule"
    CALL add_vector(this%Data,vec%Data,res%Data)
  END SUBROUTINE mp_add_vector_real

  SUBROUTINE mp_subtract_vector_real(this,vec,res)
    TYPE(VectorReal),INTENT(IN) :: this,vec
    TYPE(VectorReal),intent(INOUT) :: res
    ! calculate res = this - vec

    CALL subtract_vector(this%Data,vec%Data, res%Data)
  END SUBROUTINE mp_subtract_vector_real

  SUBROUTINE mp_subtract_from_vector_real(this,vec)
    TYPE(VectorReal),INTENT(INOUT) :: this
    TYPE(VectorReal), intent(IN) :: vec
    ! calculate this = this - vec

    CALL subtract_vector(this%Data,vec%Data)
  END SUBROUTINE mp_subtract_from_vector_real
    
  SUBROUTINE mp_multiply_vector_with_real_real(this, scalar, res)
    TYPE(VectorReal), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(VectorReal),intent(INOUT) :: res

    CALL multiply_vector_with_scalar(this%Data,scalar, res%Data)
  END SUBROUTINE mp_multiply_vector_with_real_real

  SUBROUTINE mp_scale_vector_by_real_real(this, scalar)
    TYPE(VectorReal), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL multiply_vector_with_scalar(this%Data,scalar)
  END SUBROUTINE mp_scale_vector_by_real_real

  SUBROUTINE mp_assign_vector_real(lvec,rvec)
    TYPE(VectorReal), intent(INOUT) :: lvec
    TYPE(VectorReal), intent(IN) :: rvec

    IF (lvec%Nrows.EQ.rvec%Nrows) THEN
       lvec%Data = rvec%Data
    ELSE
       PRINT*, "Left and right vecrices in an assignment must have the same dimensions."
    END IF
  END SUBROUTINE mp_assign_vector_real

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_vector_sum_0D_real(localsum,communicator)
    TYPE(VectorReal), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_sum(localsum%Data,communicator)
  END SUBROUTINE my_vector_sum_0D_real

  SUBROUTINE my_vector_sum_2D_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(VectorReal), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_vector_sum_generic_real(localsum,len,communicator)
  END SUBROUTINE my_vector_sum_2D_real

  SUBROUTINE my_vector_sum_3D_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(VectorReal), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_vector_sum_generic_real(localsum,len,communicator)
  END SUBROUTINE my_vector_sum_3D_real
  
  SUBROUTINE my_vector_sum_generic_real(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    TYPE(VectorReal), intent(INOUT) :: localsum(len)
    integer :: communicator

    integer :: iVector
    DO iVector=1,len
       CALL my_sum(localsum(iVector)%Data,communicator)
    END DO
  END SUBROUTINE my_vector_sum_generic_real

END MODULE VectorModule

