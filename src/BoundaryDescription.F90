#include "intrinsic_sizes.h"
#include "redef.h"
MODULE BoundaryDescriptionModule
  USE mpi
  IMPLICIT NONE
  
  PUBLIC :: BoundaryDescription, static_size_of_BoundaryDescription
  PUBLIC :: initialize_type, finalize_type, set_mpi_type, set_mpi_type_real
  private

  !> A full description of the boundary, needed for the exchange.
  !!
  TYPE BoundaryDescription
     INTEGER :: lower  !< number of boundary points on lower side
     INTEGER :: upper  !< number of boundary points on upper side
     INTEGER :: count  !< how often exchange the direction, dimension of 2nd direction
     INTEGER :: n_points !< number of points (incl. boundary points) in 1st dimension
     !>index of first inner point, when counting starts with 0 for the most left point
     INTEGER :: innerfirst 
     !> index of last inner point, when counting starts with 0 for the most left point
     INTEGER :: innerlast  
     INTEGER :: subarray_size      !< size of the subarrays of which consists the dimension to exchange
     INTEGER :: n_upper_subarrays  !< number of subarrays in the upper boundary
     INTEGER :: n_lower_subarrays  !< number of subarrays in the lower boundary
     LOGICAL :: mpi_type_allocated !< is true if two mpi datatypes have been set and committed
     INTEGER :: mpi_lower_datatype !< user defined datatypes for the exchange of lower boundary
     INTEGER :: mpi_upper_datatype !< user defined datatypes for the exchange of upper boundary
     !INTEGER :: communicator !< MPI communicator for exchange
     INTEGER :: exchange_direction !< 0:spec, 1:w, 2:v, 3:z, 4:y, 5:x
  END TYPE BoundaryDescription

  INTERFACE initialize_type
     MODULE PROCEDURE mp_initialize_type
  END INTERFACE

  INTERFACE finalize_type
     MODULE PROCEDURE mp_finalize_type
  END INTERFACE

  INTERFACE set_mpi_type
     MODULE PROCEDURE mp_set_mpi_type
  END INTERFACE

  INTERFACE set_mpi_type_real
     MODULE PROCEDURE mp_set_mpi_type_real
  END INTERFACE

CONTAINS

  FUNCTION static_size_of_BoundaryDescription() RESULT(memory_need)
    integer :: memory_need

    memory_need = 12*SIZE_OF_INTEGER&
         & + SIZE_OF_LOGICAL
  END FUNCTION static_size_of_BoundaryDescription

  SUBROUTINE mp_initialize_type(bdesc,a_subarray_size,&
       & a_n_subarrays,&
       & a_n_lower_subarrays,a_n_upper_subarrays,&
       & a_count)
    TYPE(BoundaryDescription),intent(OUT) :: bdesc
    INTEGER, INTENT(IN) :: a_subarray_size, a_n_subarrays, &
         & a_n_lower_subarrays, a_n_upper_subarrays, a_count

    bdesc%lower = a_subarray_size*a_n_lower_subarrays
    bdesc%upper = a_subarray_size*a_n_upper_subarrays
    bdesc%count = a_count
    bdesc%n_points = a_n_subarrays*a_subarray_size
    bdesc%innerfirst = bdesc%lower !array starting with 0, then lower is the first inner index
    bdesc%innerlast = bdesc%n_points-bdesc%upper-1
    bdesc%subarray_size = a_subarray_size
    bdesc%n_upper_subarrays = a_n_upper_subarrays
    bdesc%n_lower_subarrays = a_n_lower_subarrays
    bdesc%mpi_type_allocated = .FALSE.
    !bdesc%exchange_direction = a_exdir
    !bdesc%communicator = a_comm
  END SUBROUTINE mp_initialize_type

  SUBROUTINE mp_set_mpi_type(bdesc)
    type(BoundaryDescription) :: bdesc

    ! Local variables
    integer :: ierr

    CALL mpi_type_vector(bdesc%count,bdesc%lower,bdesc%n_points,&
         & MPI_COMPLEX_TYPE,bdesc%mpi_lower_datatype,ierr)
    CALL mpi_type_commit(bdesc%mpi_lower_datatype, ierr)

    CALL mpi_type_vector(bdesc%count,bdesc%upper,bdesc%n_points,&
         & MPI_COMPLEX_TYPE,bdesc%mpi_upper_datatype,ierr)
    CALL mpi_type_commit(bdesc%mpi_upper_datatype, ierr)

    bdesc%mpi_type_allocated = .TRUE.
  END SUBROUTINE mp_set_mpi_type

  SUBROUTINE mp_set_mpi_type_real(bdesc)
    type(BoundaryDescription) :: bdesc

    ! Local variables
    integer :: ierr

    CALL mpi_type_vector(bdesc%count,bdesc%lower,bdesc%n_points,&
         & MPI_REAL_TYPE,bdesc%mpi_lower_datatype,ierr)
    CALL mpi_type_commit(bdesc%mpi_lower_datatype, ierr)

    CALL mpi_type_vector(bdesc%count,bdesc%upper,bdesc%n_points,&
         & MPI_REAL_TYPE,bdesc%mpi_upper_datatype,ierr)
    CALL mpi_type_commit(bdesc%mpi_upper_datatype, ierr)

    bdesc%mpi_type_allocated = .TRUE.
  END SUBROUTINE mp_set_mpi_type_real


  SUBROUTINE mp_finalize_type(bdesc)
    type(BoundaryDescription) :: bdesc

    integer :: ierr

    IF (bdesc%mpi_type_allocated) THEN
       CALL mpi_type_free(bdesc%mpi_lower_datatype,ierr)
       CALL mpi_type_free(bdesc%mpi_upper_datatype,ierr)
       bdesc%mpi_type_allocated = .FALSE.
    END IF
  END SUBROUTINE mp_finalize_type

END MODULE BoundaryDescriptionModule
