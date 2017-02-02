program test_VectorModule
  USE VectorModule
  use MatrixModule
  use mpi
  IMPLICIT NONE

  logical :: passed_all_tests=.true.
  INTEGER :: ierr, n_procs, mype, comm_cart
  INTEGER :: nrows
  TYPE(Vector) :: avec,bvec
  REAL :: tolerance=1e-10

  CALL mpi_init(ierr)  
  CALL mpi_comm_size(MPI_COMM_WORLD,n_procs, ierr)
  CALL mpi_cart_create(MPI_COMM_WORLD, 1, (/n_procs/), (/.FALSE./), .FALSE., comm_cart,ierr)

  call initialize_matrix_module(comm_cart)
  CALL initialize_Vector_module

  CALL mpi_comm_size(comm_cart,n_procs, ierr)
  CALL mpi_comm_rank(comm_cart,mype, ierr)
  PRINT*,"Using ",n_procs," processors."
  nrows = 8

  !-----------------------------------------
  ! Testing interface initialize
  !-----------------------------------------
  if (test_mp_initialize_vector()) then
     write(*,"(A40,I3,A10)") "test_mp_initialize_vector on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_initialize_vector on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface initialize_vector_module
  !-----------------------------------------
  if (test_mp_initialize_vector_module()) then
     write(*,"(A40,I3,A10)") "test_mp_initialize_vector_module on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_initialize_vector_module on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface allocate
  !-----------------------------------------
  if (test_mp_allocate()) then
     write(*,"(A40,I3,A10)") "test_mp_allocate on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_allocate on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface attach
  !-----------------------------------------
  if (test_mp_attach_vector()) then
     write(*,"(A40,I3,A10)") "test_mp_attach_vector on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_attach_vector on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface finalize
  !-----------------------------------------
  if (test_mp_finalize_vector()) then
     write(*,"(A40,I3,A10)") "test_mp_finalize_vector on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_finalize_vector on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface finalize_vector_module
  !-----------------------------------------
  if (test_mp_finalize_vector_module()) then
     write(*,"(A40,I3,A10)") "test_mp_finalize_vector_module on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_finalize_vector_module on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface set_value
  !-----------------------------------------
  if (test_mp_set_complex_value()) then
     write(*,"(A40,I3,A10)") "test_mp_set_complex_value on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_set_complex_value on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  if (test_mp_set_real_value()) then
     write(*,"(A40,I3,A10)") "test_mp_set_real_value on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_set_real_value on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface add_value
  !-----------------------------------------
  if (test_mp_add_complex_value()) then
     write(*,"(A40,I3,A10)") "test_mp_add_complex_value on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_add_complex_value on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  if (test_mp_add_real_value()) then
     write(*,"(A40,I3,A10)") "test_mp_add_real_value on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_add_real_value on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface commit_values
  !-----------------------------------------
  if (test_mp_commit_values()) then
     write(*,"(A40,I3,A10)") "test_mp_commit_values on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_commit_values on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface set_zero
  !-----------------------------------------
  if (test_mp_set_zero()) then
     write(*,"(A40,I3,A10)") "test_mp_set_zero on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_set_zero on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface get_global_vector_locally
  !-----------------------------------------
  if (test_mp_get_global_vector_locally()) then
     write(*,"(A40,I3,A10)") "test_mp_get_global_vector_locally on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_get_global_vector_locally on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface get_local_abs_square_sum
  !-----------------------------------------
  if (test_mp_get_local_abs_square_sum()) then
     write(*,"(A40,I3,A10)") "test_mp_get_local_abs_square_sum on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_get_local_abs_square_sum on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface vec_get_value
  !-----------------------------------------
  if (test_mp_get_value()) then
     write(*,"(A40,I3,A10)") "test_mp_get_value on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_get_value on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface dot_multiply
  !-----------------------------------------
  if (test_mp_dot_multiply()) then
     write(*,"(A40,I3,A10)") "test_mp_dot_multiply on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_dot_multiply on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface show
  !-----------------------------------------
  if (test_mp_show_on_screen()) then
     write(*,"(A40,I3,A10)") "test_mp_show_on_screen on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_show_on_screen on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  if (test_mp_show_in_file()) then
     write(*,"(A40,I3,A10)") "test_mp_show_in_file on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_show_in_file on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface add_vector
  !-----------------------------------------
  if (test_mp_add_to_vector()) then
     write(*,"(A40,I3,A10)") "test_mp_add_to_vector on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_add_to_vector on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  if (test_mp_add_vector()) then
     write(*,"(A40,I3,A10)") "test_mp_add_vector on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_add_vector on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface subtract_vector
  !-----------------------------------------
  if (test_mp_subtract_vector()) then
     write(*,"(A40,I3,A10)") "test_mp_subtract_vector on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_subtract_vector on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  if (test_mp_subtract_from_vector()) then
     write(*,"(A40,I3,A10)") "test_mp_subtract_from_vector on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_subtract_from_vector on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface multiply_vector_with_scalar
  !-----------------------------------------
  if (test_mp_multiply_vector_with_real()) then
     write(*,"(A40,I3,A10)") "test_mp_multiply_vector_with_real on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_multiply_vector_with_real on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  if (test_mp_scale_vector_by_real()) then
     write(*,"(A40,I3,A10)") "test_mp_scale_vector_by_real on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_scale_vector_by_real on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface ASSIGNMENT_equal
  !-----------------------------------------
  if (test_mp_assign_vector()) then
     write(*,"(A40,I3,A10)") "test_mp_assign_vector on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_assign_vector on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface isInitialized
  !-----------------------------------------
  if (test_mp_isInitialized()) then
     write(*,"(A40,I3,A10)") "test_mp_isInitialized on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_mp_isInitialized on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  !-----------------------------------------
  ! Testing interface my_sum
  !-----------------------------------------
  if (test_my_vector_sum_generic()) then
     write(*,"(A40,I3,A10)") "test_my_vector_sum_generic on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_my_vector_sum_generic on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  if (test_my_vector_sum_2D()) then
     write(*,"(A40,I3,A10)") "test_my_vector_sum_2D on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_my_vector_sum_2D on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  if (test_my_vector_sum_3D()) then
     write(*,"(A40,I3,A10)") "test_my_vector_sum_3D on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_my_vector_sum_3D on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  if (test_my_vector_sum_0D()) then
     write(*,"(A40,I3,A10)") "test_my_vector_sum_0D on process ",mype," passed"
  else
     write(*,"(A40,I3,A10)") "test_my_vector_sum_0D on process ",mype,"  FAILED"
     passed_all_tests=.false.
  end if

  if (passed_all_tests) then
     print*,'All tests on process ',mype,' passed.'
  else
     print*,'Some tests on process ',mype,' FAILED.'
  end if

  CALL finalize_Vector_module
  call finalize_matrix_module
  call mpi_finalize(ierr)

contains

  function test_mp_initialize_vector() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_initialize_vector

  function test_mp_initialize_vector_module() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_initialize_vector_module

  function test_mp_allocate() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_allocate

  function test_mp_attach_vector() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_attach_vector

  function test_mp_finalize_vector() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_finalize_vector

  function test_mp_finalize_vector_module() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_finalize_vector_module

  function test_mp_set_complex_value() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_set_complex_value

  function test_mp_set_real_value() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_set_real_value

  function test_mp_add_complex_value() result(passed)
    logical :: passed

    integer :: irow
    COMPLEX, dimension(nrows) :: localfullvec
    COMPLEX, DIMENSION(nrows/n_procs) :: localvec

    CALL initialize(avec,nrows)
    CALL attach(avec,localvec)
    CALL initialize(bvec,nrows)
    call allocate(bvec)

    DO irow=1,nrows
       CALL set_value(avec,irow,CMPLX(irow,irow**2))
       CALL set_value(bvec, irow, CMPLX(2*irow,3*irow))
    END DO
    call commit_values(avec)
    call commit_values(bvec)
    !CALL show(avec)
    !CALL show(bvec)

    CALL add_vector(avec,bvec)
    !call show(avec)
    
    CALL get_global_vector_locally(avec,localfullvec)
    passed=.true.
    DO irow=1,nrows
       IF (ABS(localfullvec(irow)-CMPLX(3*irow,irow**2+3*irow)).GT.tolerance) passed=.false.
    END DO
    
    CALL finalize(avec)
  end function test_mp_add_complex_value

  function test_mp_add_real_value() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_add_real_value

  function test_mp_commit_values() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_commit_values

  function test_mp_set_zero() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_set_zero

  function test_mp_get_global_vector_locally() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_get_global_vector_locally

  function test_mp_get_local_abs_square_sum() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_get_local_abs_square_sum

  function test_mp_get_value() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_get_value

  function test_mp_dot_multiply() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_dot_multiply

  function test_mp_show_on_screen() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_show_on_screen

  function test_mp_show_in_file() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_show_in_file

  function test_mp_add_to_vector() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_add_to_vector

  function test_mp_add_vector() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_add_vector

  function test_mp_subtract_vector() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_subtract_vector

  function test_mp_subtract_from_vector() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_subtract_from_vector

  function test_mp_multiply_vector_with_real() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_multiply_vector_with_real

  function test_mp_scale_vector_by_real() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_scale_vector_by_real

  function test_mp_assign_vector() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_assign_vector

  function test_mp_isInitialized() result(passed)
    logical :: passed

    passed=.true.
  end function test_mp_isInitialized

  function test_my_vector_sum_generic() result(passed)
    logical :: passed

    passed=.true.
  end function test_my_vector_sum_generic

  function test_my_vector_sum_2D() result(passed)
    logical :: passed

    passed=.true.
  end function test_my_vector_sum_2D

  function test_my_vector_sum_3D() result(passed)
    logical :: passed

    passed=.true.
  end function test_my_vector_sum_3D

  function test_my_vector_sum_0D() result(passed)
    logical :: passed

    passed=.true.
  end function test_my_vector_sum_0D
end program test_VectorModule
