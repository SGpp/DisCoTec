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

MODULE StoreVectorModule
  USE ProcessGridModule
  USE mpi
  IMPLICIT NONE
  private

  TYPE,public :: StoreVectorObject 
     TYPE(ProcessGrid) :: PG
     INTEGER :: NRows ! number of rows and cols globally
     INTEGER :: RowsPerBlock ! block size of distributed vector
     LOGICAL :: isAttached
     LOGICAL :: isInitialized=.FALSE.
     INTEGER,DIMENSION(7) :: Desc           ! vector descriptor for BLACS, ScaLapack etc. 
     ! the per processor local vector part of the global vector
     COMPLEX, DIMENSION(:), POINTER :: localVector
  END TYPE StoreVectorObject

 TYPE,public :: StoreVectorObjectReal
     TYPE(ProcessGrid) :: PG
     INTEGER :: NRows ! number of rows and cols globally
     INTEGER :: RowsPerBlock ! block size of distributed vector
     LOGICAL :: isAttached
     LOGICAL :: isInitialized=.FALSE.
     INTEGER,DIMENSION(7) :: Desc           ! vector descriptor for BLACS, ScaLapack etc. 
     ! the per processor local vector part of the global vector
     REAL, DIMENSION(:), POINTER :: localVector
  END TYPE StoreVectorObjectReal

  INTERFACE initialize
     MODULE PROCEDURE svo_initialize_vector, svo_initialize_vector_real
  END INTERFACE

  INTERFACE isInitialized
     module procedure svo_isInitialized, svo_isInitialized_real
  END INTERFACE

  INTERFACE allocate
     module procedure svo_allocate, svo_allocate_real
  END INTERFACE

  INTERFACE attach
     MODULE PROCEDURE svo_attach_vector, svo_attach_vector_real
  END INTERFACE

  INTERFACE set_value
     MODULE PROCEDURE svo_set_value, svo_set_value_real
  END INTERFACE

  INTERFACE add_value
     MODULE PROCEDURE svo_add_value, svo_add_value_real
  END INTERFACE

  INTERFACE finalize
     MODULE PROCEDURE svo_finalize_vector, svo_finalize_vector_real
  END INTERFACE
  
  INTERFACE commit_values
     MODULE PROCEDURE svo_commit_values, svo_commit_values_real
  END INTERFACE

  INTERFACE get_global_vector_locally
     module procedure svo_get_global_vector_locally, svo_get_global_vector_locally_real
  END INTERFACE

  INTERFACE vec_get_value
     MODULE PROCEDURE svo_get_value, svo_get_value_real
  END INTERFACE

  INTERFACE get_local_abs_square_sum
     module procedure svo_get_local_abs_square_sum, svo_get_local_abs_square_sum_real
  END INTERFACE

  INTERFACE show
     MODULE PROCEDURE svo_show_on_screen, svo_show_in_file, svo_show_on_screen_real, svo_show_in_file_real
  END INTERFACE

  INTERFACE add_vector
     MODULE PROCEDURE svo_add_to_vector, svo_add_vector, svo_add_to_vector_real, svo_add_vector_real
  END INTERFACE

  INTERFACE subtract_vector
  !INTERFACE OPERATOR(-)
     MODULE PROCEDURE svo_subtract_vector, svo_subtract_from_vector, svo_subtract_vector_real, svo_subtract_from_vector_real
  END INTERFACE

  interface multiply_vector_with_scalar
  !INTERFACE OPERATOR(*)
     MODULE PROCEDURE svo_multiply_vector_with_real
     module procedure svo_scale_vector_by_real
     module procedure svo_multiply_vector_with_real_real
     module procedure svo_scale_vector_by_real_real
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE svo_assign_vector, svo_assign_vector_real
  END INTERFACE

  INTERFACE my_sum
     MODULE PROCEDURE my_vector_sum, my_vector_sum_real
  END INTERFACE

  PUBLIC :: initialize, finalize, isInitialized, ALLOCATE, attach
  PUBLIC :: set_value, add_value, vec_get_value, commit_values, get_global_vector_locally
  PUBLIC :: get_local_abs_square_sum, show, my_sum
  PUBLIC :: add_vector, subtract_vector, multiply_vector_with_scalar
  public :: static_size_of_storevectorobject, assignment(=)

CONTAINS
  FUNCTION static_size_of_storevectorobject() RESULT(memory_need)
    INTEGER :: memory_need

    memory_need = 9*SIZE_OF_INTEGER &
         & + 2*SIZE_OF_LOGICAL &
         & + SIZE_OF_POINTER &
         & + size_of_processgrid()
  END FUNCTION static_size_of_storevectorobject
  
  SUBROUTINE svo_initialize_vector(vec,n_rows,PG)
    TYPE(StoreVectorObject), intent(INOUT) :: vec
    INTEGER :: n_rows
    type(ProcessGrid) :: PG

    ! Local variables

    vec%NRows = n_rows
    vec%PG = PG
    ! Distribute vector on the grid
    vec%RowsPerBlock = vec%NRows / PG%NProcs

    !PRINT*,"NPCols = ", PG%NPCols,", ColsPerBlock = ", vec%ColsPerBlock
    !vec%Desc=1000
    vec%Desc = (/ 502, PG%Context, vec%NRows, vec%RowsPerBlock, 0, vec%RowsPerBlock,0 /)

    vec%isInitialized = .TRUE.
    CALL blacs_barrier(PG%Context,'A')
  END SUBROUTINE svo_initialize_vector

  LOGICAL FUNCTION svo_isInitialized(vec) 
    type(StoreVectorObject) :: vec

    svo_isInitialized = vec%isInitialized
  END FUNCTION svo_isInitialized


  SUBROUTINE svo_allocate(vec)
    TYPE(StoreVectorObject) :: vec

    ALLOCATE(vec%localVector(vec%RowsPerBlock))

    vec%isAttached = .FALSE.
  END SUBROUTINE svo_allocate

  SUBROUTINE svo_attach_vector(vec,localarray)
    type(StoreVectorObject) :: vec
    COMPLEX, DIMENSION(:),TARGET :: localarray

    vec%localVector => localarray(1:vec%RowsPerBlock)
    vec%isAttached = .TRUE.
  END SUBROUTINE svo_attach_vector

  SUBROUTINE svo_finalize_vector(vec)
    type(StoreVectorObject) :: vec

    IF (vec%isAttached) THEN
       ! do nothing
       vec%isAttached=.FALSE.
    ELSEIF (ASSOCIATED(vec%localVector)) THEN
       DEALLOCATE(vec%localVector)
    END IF
    NULLIFY(vec%localVector)
    vec%isInitialized=.false.
  END SUBROUTINE svo_finalize_vector

  SUBROUTINE svo_set_value(vec,irow,value)
    TYPE(StoreVectorObject),INTENT(INOUT) :: vec
    INTEGER,INTENT(IN) :: irow
    COMPLEX, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc, row_index

    DEBSTART("set_value(SVO)")
    !PRINT*,"In routine set_value of StoreVectorModule on PRow,PCol = ",vec%PRow,vec%PCol
    ! calculate the processor coordinates (starting with 0) of the given global indices
    proc = (irow-1)/vec%RowsPerBlock

    !PRINT*,vec%Rank,vec%PRow,vec%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF (proc.EQ.vec%PG%Rank) THEN
       row_index = MOD(irow-1,vec%RowsPerBlock)+1
       vec%localVector(row_index) = value
    END IF
    DEBEND("set_value(SVO)")
  END SUBROUTINE svo_set_value

  SUBROUTINE svo_add_value(vec,irow,value)
    TYPE(StoreVectorObject) :: vec
    INTEGER, INTENT(IN) :: irow
    COMPLEX, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc

    DEBSTART("add_value(SVO)")

    proc = (irow-1)/vec%RowsPerBlock

    !PRINT*,vec%Rank,vec%PRow,vec%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF (proc.EQ.vec%PG%Rank) THEN
       CALL set_value(vec,irow,vec_get_value(vec,irow)+value)
    END IF
    DEBEND("add_value(SVO)")
    
  END SUBROUTINE svo_add_value

  SUBROUTINE svo_commit_values(vec)
    TYPE(StoreVectorObject) :: vec

    ! nothing happens, as a Vector can not be further optimized

  END SUBROUTINE svo_commit_values

  !*****************************************************************
  !* The global vector, which is distributed over all processes
  !* is gathered an stored in a local array. 
  !*****************************************************************
  SUBROUTINE svo_get_global_vector_locally(vec,localfullvec)
    TYPE(StoreVectorObject),INTENT(IN) :: vec
    COMPLEX, DIMENSION(:),INTENT(OUT) :: localfullvec

    ! local variables
    INTEGER :: ierr

    ! If vec is a Vector, just concatenate the local vectors on all processes.
    CALL mpi_allgather(vec%localVector,vec%RowsPerBlock,MPI_COMPLEX_TYPE,&
         &localfullvec,vec%RowsPerBlock,&
         &MPI_COMPLEX_TYPE,vec%PG%Communicator,ierr)
  END SUBROUTINE svo_get_global_vector_locally

  FUNCTION svo_get_value(vec,irow) RESULT(value)
    TYPE(StoreVectorObject) :: vec
    INTEGER :: irow
    complex :: value

    ! local variables
    INTEGER :: proc, row_index

    proc = (irow-1)/vec%RowsPerBlock
    !PRINT*,vec%Rank,vec%PRow,vec%PCol,proc_r, proc_c
    ! return the entry, if we are on the right processor
    IF (proc.EQ.vec%PG%Rank) THEN
       ! calculate the local block indices in which the global indices lie
       row_index = MOD(irow-1,vec%RowsPerBlock)+1
       value = vec%localVector(row_index)
    ELSE
       ! return nothing
       value = CMPLX(1000,700000)
       ! this can be arbitrary, as it is never used in further calculations
    END IF
  END FUNCTION svo_get_value

  FUNCTION svo_get_local_abs_square_sum(vec) RESULT(value)
    type(StoreVectorObject) :: vec
    REAL :: value

    value = REAL(SUM(vec%localVector*CONJG(vec%localVector)),KIND(value))
  END FUNCTION svo_get_local_abs_square_sum

  SUBROUTINE svo_show_on_screen(vec)
    type(StoreVectorObject) :: vec

    ! Local variables
    INTEGER :: iProc, iRow

    ! do some debugging output
    CALL blacs_barrier(vec%PG%Context,'A')
    DO iProc=0,vec%PG%NProcs-1
       CALL blacs_barrier(vec%PG%Context,'A')
       IF (iProc.EQ.vec%PG%Rank) THEN
          PRINT*,"---- Start output Rank = ",vec%PG%Rank," ----"
          DO iRow=1,vec%RowsPerBlock
             WRITE(*,"(2ES10.2,1X)") vec%localVector(iRow)
          END DO
          PRINT*,"---- End   output Rank = ",vec%PG%Rank," ----"
       END IF
    END DO
    CALL blacs_barrier(vec%PG%Context,'A')
  END SUBROUTINE svo_show_on_screen

  SUBROUTINE svo_show_in_file(vec,filename)
    type(StoreVectorObject) :: vec
    CHARACTER(len=*) :: filename

    ! Local variables
    INTEGER :: iRow
    character(len=FILENAME_MAX) :: full_filename

    ! open a file on each process
    WRITE(full_filename,"(A,I3.3,A)") filename,vec%PG%Rank,".dat"
    OPEN(77,file=TRIM(full_filename))

    ! do some debugging output
    CALL blacs_barrier(vec%PG%Context,'A')
    DO iRow=1,vec%RowsPerBlock
       WRITE(*,"(2ES10.2,1X)") vec%localVector(iRow)
    END DO
    CLOSE(77)
    CALL blacs_barrier(vec%PG%Context,'A')
  END SUBROUTINE svo_show_in_file

  ! =======================================
  ! == Define some operators on vectors  ==
  ! =======================================
  SUBROUTINE svo_add_to_vector(this,vec)
    TYPE(StoreVectorObject), INTENT(INOUT) :: this
    TYPE(StoreVectorObject), INTENT(IN) :: vec

    this%localVector = this%localVector + vec%localVector
  END SUBROUTINE svo_add_to_vector

  SUBROUTINE svo_add_vector(this,vec,res)
    TYPE(StoreVectorObject), INTENT(IN) :: this,vec
    TYPE(StoreVectorObject), INTENT(INOUT) :: res

    res%localVector = this%localVector+vec%localVector
  END SUBROUTINE svo_add_vector

  SUBROUTINE svo_subtract_vector(this,vec,res)
    TYPE(StoreVectorObject),INTENT(IN) :: this,vec
    TYPE(StoreVectorObject),intent(INOUT) :: res
    ! calculate res = this - vec

    res%localVector = this%localVector - vec%localVector
  END SUBROUTINE svo_subtract_vector

  SUBROUTINE svo_subtract_from_vector(this,vec)
    TYPE(StoreVectorObject),INTENT(INOUT) :: this
    TYPE(StoreVectorObject),intent(IN) :: vec
    ! calculate this = this - vec

    this%localVector = this%localVector - vec%localVector
  END SUBROUTINE svo_subtract_from_vector

  SUBROUTINE svo_multiply_vector_with_real(this, scalar, res)
    TYPE(StoreVectorObject), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(StoreVectorObject),intent(INOUT) :: res

    res%localVector = this%localVector*scalar
  END SUBROUTINE svo_multiply_vector_with_real

  SUBROUTINE svo_scale_vector_by_real(this, scalar)
    TYPE(StoreVectorObject), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    this%localVector = this%localVector*scalar
  END SUBROUTINE svo_scale_vector_by_real

  SUBROUTINE svo_assign_vector(lvec,rvec)
    TYPE(StoreVectorObject), intent(INOUT) :: lvec
    TYPE(StoreVectorObject), intent(IN) :: rvec

    IF (lvec%PG.EQ.rvec%PG) THEN
       lvec%localVector = rvec%localVector
    ELSE
       PRINT*,"Assignment of vectors is only possible within the same context."
    END IF
  END SUBROUTINE svo_assign_vector

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_vector_sum(vec,communicator)
    TYPE(StoreVectorObject), INTENT(INOUT) :: vec
    COMPLEX, DIMENSION(vec%RowsPerBlock) :: totalsum
    integer :: communicator

    integer :: ierr
    CALL mpi_allreduce(vec%localVector,totalsum,vec%RowsPerBlock, &
         & MPI_COMPLEX_TYPE, MPI_SUM, communicator, ierr)
    
    vec%localVector = totalsum
  END SUBROUTINE my_vector_sum

  ! ==========================================
  ! == Repeat EVERYTHING for REAL datatypes ==
  ! ==========================================

  SUBROUTINE svo_initialize_vector_real(vec,n_rows,PG)
    TYPE(StoreVectorObjectReal), intent(INOUT) :: vec
    INTEGER :: n_rows
    type(ProcessGrid) :: PG

    ! Local variables

    vec%NRows = n_rows
    vec%PG = PG
    ! Distribute vector on the grid
    vec%RowsPerBlock = vec%NRows / PG%NProcs

    !PRINT*,"NPCols = ", PG%NPCols,", ColsPerBlock = ", vec%ColsPerBlock
    !vec%Desc=1000
    vec%Desc = (/ 502, PG%Context, vec%NRows, vec%RowsPerBlock, 0, vec%RowsPerBlock,0 /)

    vec%isInitialized = .TRUE.
    CALL blacs_barrier(PG%Context,'A')
  END SUBROUTINE svo_initialize_vector_real

  LOGICAL FUNCTION svo_isInitialized_real(vec) 
    type(StoreVectorObjectReal) :: vec

    svo_isInitialized_real = vec%isInitialized
  END FUNCTION svo_isInitialized_real


  SUBROUTINE svo_allocate_real(vec)
    TYPE(StoreVectorObjectReal) :: vec

    ALLOCATE(vec%localVector(vec%RowsPerBlock))

    vec%isAttached = .FALSE.
  END SUBROUTINE svo_allocate_real

  SUBROUTINE svo_attach_vector_real(vec,localarray)
    type(StoreVectorObjectReal) :: vec
    REAL, DIMENSION(:),TARGET :: localarray

    vec%localVector => localarray(1:vec%RowsPerBlock)
    vec%isAttached = .TRUE.
  END SUBROUTINE svo_attach_vector_real

  SUBROUTINE svo_finalize_vector_real(vec)
    type(StoreVectorObjectReal) :: vec

    IF (vec%isAttached) THEN
       ! do nothing
       vec%isAttached=.FALSE.
    ELSEIF (ASSOCIATED(vec%localVector)) THEN
       DEALLOCATE(vec%localVector)
    END IF
    NULLIFY(vec%localVector)
    vec%isInitialized=.false.
  END SUBROUTINE svo_finalize_vector_real

  SUBROUTINE svo_set_value_real(vec,irow,value)
    TYPE(StoreVectorObjectReal),INTENT(INOUT) :: vec
    INTEGER,INTENT(IN) :: irow
    REAL, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc, row_index

    DEBSTART("set_value(SVO)")
    !PRINT*,"In routine set_value of StoreVectorModule on PRow,PCol = ",vec%PRow,vec%PCol
    ! calculate the processor coordinates (starting with 0) of the given global indices
    proc = (irow-1)/vec%RowsPerBlock

    !PRINT*,vec%Rank,vec%PRow,vec%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF (proc.EQ.vec%PG%Rank) THEN
       row_index = MOD(irow-1,vec%RowsPerBlock)+1
       vec%localVector(row_index) = value
    END IF
    DEBEND("set_value(SVO)")
  END SUBROUTINE svo_set_value_real

  SUBROUTINE svo_add_value_real(vec,irow,value)
    TYPE(StoreVectorObjectReal) :: vec
    INTEGER, INTENT(IN) :: irow
    REAL, INTENT(IN) :: value

    ! Local variables
    INTEGER :: proc

    DEBSTART("add_value(SVO)")

    proc = (irow-1)/vec%RowsPerBlock

    !PRINT*,vec%Rank,vec%PRow,vec%PCol,proc_r, proc_c
    ! only process the entry if we are on the right processor
    IF (proc.EQ.vec%PG%Rank) THEN
       CALL set_value(vec,irow,vec_get_value(vec,irow)+value)
    END IF
    DEBEND("add_value(SVO)")
    
  END SUBROUTINE svo_add_value_real

  SUBROUTINE svo_commit_values_real(vec)
    TYPE(StoreVectorObjectReal) :: vec

    ! nothing happens, as a Vector can not be further optimized

  END SUBROUTINE svo_commit_values_real

  !*****************************************************************
  !* The global vector, which is distributed over all processes
  !* is gathered an stored in a local array. 
  !*****************************************************************
  SUBROUTINE svo_get_global_vector_locally_real(vec,localfullvec)
    TYPE(StoreVectorObjectReal),INTENT(IN) :: vec
    REAL, DIMENSION(:),INTENT(OUT) :: localfullvec

    ! local variables
    INTEGER :: ierr

    ! If vec is a Vector, just concatenate the local vectors on all processes.
    CALL mpi_allgather(vec%localVector,vec%RowsPerBlock,MPI_REAL_TYPE,&
         &localfullvec,vec%RowsPerBlock,&
         &MPI_REAL_TYPE,vec%PG%Communicator,ierr)
  END SUBROUTINE svo_get_global_vector_locally_real

  FUNCTION svo_get_value_real(vec,irow) RESULT(value)
    TYPE(StoreVectorObjectReal) :: vec
    INTEGER :: irow
    real :: value

    ! local variables
    INTEGER :: proc, row_index

    proc = (irow-1)/vec%RowsPerBlock
    !PRINT*,vec%Rank,vec%PRow,vec%PCol,proc_r, proc_c
    ! return the entry, if we are on the right processor
    IF (proc.EQ.vec%PG%Rank) THEN
       ! calculate the local block indices in which the global indices lie
       row_index = MOD(irow-1,vec%RowsPerBlock)+1
       value = vec%localVector(row_index)
    ELSE
       ! return nothing
       value = 1000.0
       ! this can be arbitrary, as it is never used in further calculations
    END IF
  END FUNCTION svo_get_value_real

  FUNCTION svo_get_local_abs_square_sum_real(vec) RESULT(value)
    type(StoreVectorObjectReal) :: vec
    REAL :: value

    value = REAL(SUM(vec%localVector*vec%localVector),KIND(value))
  END FUNCTION svo_get_local_abs_square_sum_real

  SUBROUTINE svo_show_on_screen_real(vec)
    type(StoreVectorObjectReal) :: vec

    ! Local variables
    INTEGER :: iProc, iRow

    ! do some debugging output
    CALL blacs_barrier(vec%PG%Context,'A')
    DO iProc=0,vec%PG%NProcs-1
       CALL blacs_barrier(vec%PG%Context,'A')
       IF (iProc.EQ.vec%PG%Rank) THEN
          PRINT*,"---- Start output Rank = ",vec%PG%Rank," ----"
          DO iRow=1,vec%RowsPerBlock
             WRITE(*,"(2ES10.2,1X)") vec%localVector(iRow)
          END DO
          PRINT*,"---- End   output Rank = ",vec%PG%Rank," ----"
       END IF
    END DO
    CALL blacs_barrier(vec%PG%Context,'A')
  END SUBROUTINE svo_show_on_screen_real

  SUBROUTINE svo_show_in_file_real(vec,filename)
    type(StoreVectorObjectReal) :: vec
    CHARACTER(len=*) :: filename

    ! Local variables
    INTEGER :: iRow
    character(len=FILENAME_MAX) :: full_filename

    ! open a file on each process
    WRITE(full_filename,"(A,I3.3,A)") filename,vec%PG%Rank,".dat"
    OPEN(77,file=TRIM(full_filename))

    ! do some debugging output
    CALL blacs_barrier(vec%PG%Context,'A')
    DO iRow=1,vec%RowsPerBlock
       WRITE(*,"(2ES10.2,1X)") vec%localVector(iRow)
    END DO
    CLOSE(77)
    CALL blacs_barrier(vec%PG%Context,'A')
  END SUBROUTINE svo_show_in_file_real

  ! =======================================
  ! == Define some operators on vectors  ==
  ! =======================================
  SUBROUTINE svo_add_to_vector_real(this,vec)
    TYPE(StoreVectorObjectReal), INTENT(INOUT) :: this
    TYPE(StoreVectorObjectReal), INTENT(IN) :: vec

    this%localVector = this%localVector + vec%localVector
  END SUBROUTINE svo_add_to_vector_real

  SUBROUTINE svo_add_vector_real(this,vec,res)
    TYPE(StoreVectorObjectReal), INTENT(IN) :: this,vec
    TYPE(StoreVectorObjectReal), INTENT(INOUT) :: res

    res%localVector = this%localVector+vec%localVector
  END SUBROUTINE svo_add_vector_real

  SUBROUTINE svo_subtract_vector_real(this,vec,res)
    TYPE(StoreVectorObjectReal),INTENT(IN) :: this,vec
    TYPE(StoreVectorObjectReal),intent(INOUT) :: res
    ! calculate res = this - vec

    res%localVector = this%localVector - vec%localVector
  END SUBROUTINE svo_subtract_vector_real

  SUBROUTINE svo_subtract_from_vector_real(this,vec)
    TYPE(StoreVectorObjectReal),INTENT(INOUT) :: this
    TYPE(StoreVectorObjectReal),intent(IN) :: vec
    ! calculate this = this - vec

    this%localVector = this%localVector - vec%localVector
  END SUBROUTINE svo_subtract_from_vector_real

  SUBROUTINE svo_multiply_vector_with_real_real(this, scalar, res)
    TYPE(StoreVectorObjectReal), intent(IN) :: this
    REAL, intent(IN) :: scalar
    TYPE(StoreVectorObjectReal),intent(INOUT) :: res

    res%localVector = this%localVector*scalar
  END SUBROUTINE svo_multiply_vector_with_real_real

  SUBROUTINE svo_scale_vector_by_real_real(this, scalar)
    TYPE(StoreVectorObjectReal), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    this%localVector = this%localVector*scalar
  END SUBROUTINE svo_scale_vector_by_real_real

  SUBROUTINE svo_assign_vector_real(lvec,rvec)
    TYPE(StoreVectorObjectReal), intent(INOUT) :: lvec
    TYPE(StoreVectorObjectReal), intent(IN) :: rvec

    IF (lvec%PG.EQ.rvec%PG) THEN
       lvec%localVector = rvec%localVector
    ELSE
       PRINT*,"Assignment of vecrices is only possible within the same context."
    END IF
  END SUBROUTINE svo_assign_vector_real

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_vector_sum_real(vec,communicator)
    TYPE(StoreVectorObjectReal), INTENT(INOUT) :: vec
    REAL, DIMENSION(vec%RowsPerBlock) :: totalsum
    integer :: communicator

    integer :: ierr
    CALL mpi_allreduce(vec%localVector,totalsum,vec%RowsPerBlock, &
         & MPI_REAL_TYPE, MPI_SUM, communicator, ierr)
    
    vec%localVector = totalsum
  END SUBROUTINE my_vector_sum_real





END MODULE StoreVectorModule



