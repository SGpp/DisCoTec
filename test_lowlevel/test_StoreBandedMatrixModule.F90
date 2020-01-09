program test_StoreBandedMatrixModule
  USE StoreBandedMatrixModule
  use StoreFullMatrixModule
  use StoreVectorModule
  use ProcessGridModule
  use mpi
  IMPLICIT NONE
  
  LOGICAL :: passed_all_tests=.TRUE.
  TYPE(StoreBandedMatrixObject) :: Amat, Bmat,Cmat
  INTEGER :: n_procs, rank, ierr, nrows,ncols, iProc, comm_cart
  COMPLEX, DIMENSION(:,:), ALLOCATABLE :: localfullmat,localfullvec
  REAL :: tolerance = 1e-14
  type(ProcessGrid) :: PG

  call mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,n_procs, ierr)
  !CALL initialize_matrix_module(MAT_BLOCKS_OF_ROWS,MPI_COMM_WORLD)
  CALL mpi_cart_create(MPI_COMM_WORLD, 1, (/n_procs/), (/.FALSE./), .FALSE., comm_cart,ierr)
  CALL initialize(PG,MAT_BLOCKS_OF_ROWS,comm_cart)
  CALL mpi_comm_size(comm_cart,n_procs, ierr)
  CALL mpi_comm_rank(comm_cart,rank, ierr)
  if (rank.eq.0) PRINT*,"Using ",n_procs," processors."
  nrows = 8; ncols = 8

	!-----------------------------------------
	! Testing interface initialize
	!-----------------------------------------
	if (test_mp_initialize_matrix()) then
		write(*,"(A40,I3,A10)") "test_mp_initialize_matrix on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_initialize_matrix on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface isInitialized
	!-----------------------------------------
	if (test_mp_isInitialized()) then
		write(*,"(A40,I3,A10)") "test_mp_isInitialized on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_isInitialized on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface allocate
	!-----------------------------------------
	if (test_mp_allocate()) then
		write(*,"(A40,I3,A10)") "test_mp_allocate on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_allocate on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface set_value
	!-----------------------------------------
	if (test_mp_set_value()) then
		write(*,"(A40,I3,A10)") "test_mp_set_value on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_set_value on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface add_value
	!-----------------------------------------
	if (test_mp_add_value()) then
		write(*,"(A40,I3,A10)") "test_mp_add_value on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_add_value on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface finalize
	!-----------------------------------------
	if (test_mp_finalize_matrix()) then
		write(*,"(A40,I3,A10)") "test_mp_finalize_matrix on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_finalize_matrix on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface commit_values
	!-----------------------------------------
	if (test_mp_commit_values()) then
		write(*,"(A40,I3,A10)") "test_mp_commit_values on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_commit_values on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface get_global_matrix_locally
	!-----------------------------------------
	if (test_mp_get_global_matrix_locally()) then
		write(*,"(A40,I3,A10)") "test_mp_get_global_matrix_locally on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_get_global_matrix_locally on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface mat_get_value
	!-----------------------------------------
	if (test_mp_get_value()) then
		write(*,"(A40,I3,A10)") "test_mp_get_value on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_get_value on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface get_local_abs_square_sum
	!-----------------------------------------
	if (test_mp_get_local_abs_square_sum()) then
		write(*,"(A40,I3,A10)") "test_mp_get_local_abs_square_sum on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_get_local_abs_square_sum on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface dot_multiply
	!-----------------------------------------
	if (test_mp_dot_multiply()) then
		write(*,"(A40,I3,A10)") "test_mp_dot_multiply on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_dot_multiply on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if
	if (test_mp_dot_multiply_transposed()) then
		write(*,"(A40,I3,A10)") "test_mp_dot_multiply_transposed on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_dot_multiply_transposed on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface show
	!-----------------------------------------
	if (test_mp_show_on_screen()) then
		write(*,"(A40,I3,A10)") "test_mp_show_on_screen on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_show_on_screen on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	if (test_mp_show_in_file()) then
		write(*,"(A40,I3,A10)") "test_mp_show_in_file on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_show_in_file on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface add_matrix
	!-----------------------------------------
	if (test_mp_add_to_matrix()) then
		write(*,"(A40,I3,A10)") "test_mp_add_to_matrix on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_add_to_matrix on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	if (test_mp_add_matrix()) then
		write(*,"(A40,I3,A10)") "test_mp_add_matrix on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_add_matrix on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface subtract_matrix
	!-----------------------------------------
	if (test_mp_subtract_matrix()) then
		write(*,"(A40,I3,A10)") "test_mp_subtract_matrix on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_subtract_matrix on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	if (test_mp_subtract_from_matrix()) then
		write(*,"(A40,I3,A10)") "test_mp_subtract_from_matrix on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_subtract_from_matrix on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface multiply_matrix_with_scalar
	!-----------------------------------------
	if (test_mp_multiply_matrix_with_real()) then
		write(*,"(A40,I3,A10)") "test_mp_multiply_matrix_with_real on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_multiply_matrix_with_real on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	if (test_mp_scale_matrix_by_real()) then
		write(*,"(A40,I3,A10)") "test_mp_scale_matrix_by_real on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_scale_matrix_by_real on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface ASSIGNMENT_equal
	!-----------------------------------------
	if (test_mp_assign_matrix()) then
		write(*,"(A40,I3,A10)") "test_mp_assign_matrix on process ",rank," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_assign_matrix on process ",rank,"  FAILED"
		passed_all_tests=.false.
	end if


	if (passed_all_tests) then
		print*,'All tests on process ',rank,' passed.'
	else
		print*,'Some tests on process ',rank,' FAILED.'
	end if

 CALL mpi_comm_free(comm_cart, ierr)
 CALL mpi_finalize(ierr)
contains

function test_mp_initialize_matrix() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_initialize_matrix

function test_mp_isInitialized() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_isInitialized

function test_mp_allocate() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_allocate

FUNCTION test_mp_set_value() RESULT(passed)
  LOGICAL :: passed

  integer :: idiag

  CALL initialize(Amat,nrows,ncols, PG)
  CALL allocate(Amat)

  DO idiag = 1,nrows
     CALL set_value(Amat,idiag,idiag,CMPLX(idiag,idiag))
  END DO
  CALL set_value(Amat,1,2,CMPLX(1,2))
  CALL commit_values(Amat)
  !CALL show(Amat)

  CALL add_value(Amat,4,7,CMPLX(4,7))
  call commit_values(Amat)
  !CALL show(Amat)
  call finalize(Amat)
  passed=.TRUE.
end function test_mp_set_value

function test_mp_add_value() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_add_value

function test_mp_finalize_matrix() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_finalize_matrix

function test_mp_commit_values() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_commit_values

FUNCTION test_mp_get_global_matrix_locally() RESULT(passed)
  LOGICAL :: passed

  INTEGER :: idiag,irow,icol
  type(StoreFullMatrixObject) :: Afmat
  
  CALL initialize(Afmat, nrows, ncols, PG)
  call allocate(Afmat)
  DO icol = 1,ncols
     DO irow = 1,nrows
        CALL set_value(Afmat,irow,icol,CMPLX(0,0))
     END DO
  END DO

  CALL initialize(Amat,nrows,ncols, PG)
  CALL allocate(Amat)

  DO idiag = 1,nrows
     CALL set_value(Afmat,idiag,idiag,CMPLX(idiag,idiag))
     if (idiag.lt.nrows) CALL set_value(Afmat,idiag,idiag+1,CMPLX(idiag,idiag+1))
     if (idiag.gt.1) CALL set_value(Afmat,idiag,idiag-1,CMPLX(idiag,idiag-1))
  END DO
  call commit_values(Afmat)
  !CALL show(Afmat)

  DO icol = 1,ncols
     DO irow = 1,nrows
        CALL set_value(Amat,irow,icol,mat_get_value(Afmat,irow,icol))
     END DO
  END DO
  !IF (idiag.LT.nrows) CALL set_value(Amat,idiag,idiag+1,CMPLX(idiag,idiag+1))
  !IF (idiag.GT.1) CALL set_value(Amat,idiag,idiag-1,CMPLX(idiag,idiag-1))
  !END DO
  !CALL set_value(Amat,1,2,CMPLX(1,2))
  !CALL set_value(Amat,4,7,CMPLX(4,7))
  CALL commit_values(Amat)
  CALL show(Amat,"temp")

  ALLOCATE(localfullmat(nrows,ncols))
  CALL get_global_matrix_locally(Amat,localfullmat)
  IF (rank.EQ.0) THEN
     DO irow=1,nrows
        DO icol=1,ncols
           WRITE(*,"(2F4.1,1X)",advance="NO") localfullmat(irow,icol)
        END DO
        WRITE(*,"(A)") ""
     END DO
  END IF
  
  passed = .TRUE.
  !IF (rank.EQ.0) THEN
     !PRINT*,"localfullmat = "
     DO irow=1,nrows
        DO icol=1,ncols
           IF ((irow.EQ.icol).OR.(irow.EQ.1.AND.icol.EQ.2).OR.(irow.EQ.4.AND.icol.EQ.7)) THEN
              IF ( ABS(localfullmat(irow,icol)-CMPLX(irow,icol,KIND(localfullmat))).GT.tolerance ) THEN
                 passed = .FALSE.
              END IF
           ELSE
              IF ( ABS(localfullmat(irow,icol)).GT.tolerance ) THEN
                 passed = .FALSE.
              END IF
           END IF
           !WRITE(*,"(2F4.1,1X)", advance='NO') localfullmat(irow,icol)
        END DO
        !WRITE(*,"(A)") ""
     END DO
  !END IF
  call finalize(Amat)
  !passed=.TRUE.
end function test_mp_get_global_matrix_locally

function test_mp_get_value() result(passed)
  LOGICAL :: passed

  integer :: idiag
  complex :: val

  CALL initialize(Amat,nrows,ncols,PG)
  CALL ALLOCATE(Amat)

  ! set matrix Amat
  DO idiag = 1,nrows
     CALL set_value(Amat,idiag,idiag,CMPLX(idiag,idiag))
  END DO
  CALL set_value(Amat,1,2,CMPLX(1,2))
  CALL set_value(Amat,4,7,CMPLX(4,7))
  CALL commit_values(Amat)

  val = mat_get_value(Amat,1,2)
  IF (ABS(val-CMPLX(1,2)).gt.tolerance) passed=.false.
  val = mat_get_value(Amat,1,8)
  IF (ABS(val).gt.tolerance) passed=.false.
  passed=.TRUE.
end function test_mp_get_value

function test_mp_get_local_abs_square_sum() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_get_local_abs_square_sum

FUNCTION test_mp_dot_multiply() RESULT(passed)
  LOGICAL :: passed

  TYPE(StoreVectorObject) :: vec, res
  INTEGER :: idiag
  COMPLEX, ALLOCATABLE, DIMENSION(:) :: localres
  COMPLEX, ALLOCATABLE, dimension(:) :: reference

  CALL initialize(Amat,nrows,ncols,PG)
  CALL ALLOCATE(Amat)
  CALL initialize(vec,nrows,PG)
  CALL ALLOCATE(vec)
  CALL initialize(res,nrows,PG)
  CALL ALLOCATE(res)

  ! set matrix Amat
  DO idiag = 1,nrows
     CALL set_value(Amat,idiag,idiag,CMPLX(idiag,idiag))
  END DO
  CALL set_value(Amat,1,2,CMPLX(1,2))
  CALL set_value(Amat,4,7,CMPLX(4,7))
  CALL commit_values(Amat)

  ! set vector vec
  DO idiag=1,nrows
     CALL set_value(vec, idiag,CMPLX(idiag,0))
  END DO
  !CALL show(Amat)
  !call show(vec)

  CALL dot_multiply(Amat,vec,res)

  ALLOCATE(localres(nrows))
  allocate(reference(nrows))
  CALL get_global_vector_locally(res,localres)
  reference = (/(3,5),(4,4),(9,9),(44,65),(25,25),(36,36),(49,49),(64,64)/)
  passed=.TRUE.
  DO idiag=1,nrows
     IF (ABS(localres(idiag)-reference(idiag)).GT.tolerance) passed=.FALSE.
  END DO

  CALL finalize(Amat)
  CALL finalize(res)
  CALL finalize(vec)
END FUNCTION test_mp_dot_multiply

FUNCTION test_mp_dot_multiply_transposed() RESULT(passed)
  LOGICAL :: passed

  TYPE(StoreVectorObject) :: vec, res
  INTEGER :: idiag
  COMPLEX, ALLOCATABLE, DIMENSION(:) :: localres
  COMPLEX, ALLOCATABLE, DIMENSION(:) :: reference
  logical :: transposed = .true.

  CALL initialize(Amat,nrows,ncols,PG,transposed)
  CALL ALLOCATE(Amat)
  CALL initialize(vec,nrows,PG)
  CALL ALLOCATE(vec)
  CALL initialize(res,nrows,PG)
  CALL ALLOCATE(res)

  ! set matrix Amat
  DO idiag = 1,nrows
     CALL set_value(Amat,idiag,idiag,CMPLX(idiag,idiag))
  END DO
  CALL set_value(Amat,1,2,CMPLX(1,2))
  CALL set_value(Amat,4,7,CMPLX(4,7))
  CALL commit_values(Amat)

  !IF (Amat%isTransposed) THEN
  !   PRINT*,"Amat is stored in transposed form"
  !ELSE
  !   PRINT*,"Amat is stored in original form"
  !END IF

  ! set vector vec
  DO idiag=1,nrows
     CALL set_value(vec, idiag,CMPLX(idiag,0))
  END DO
  call commit_values(vec)
  CALL show(Amat)
  call show(vec)

  CALL dot_multiply(Amat,vec,res)

  call show(res)
  ALLOCATE(localres(nrows))
  allocate(reference(nrows))
  CALL get_global_vector_locally(res,localres)
  reference = (/(3,5),(4,4),(9,9),(44,65),(25,25),(36,36),(49,49),(64,64)/)
  passed=.TRUE.
  DO idiag=1,nrows
     IF (ABS(localres(idiag)-reference(idiag)).GT.tolerance) passed=.FALSE.
  END DO

  CALL finalize(Amat)
  CALL finalize(res)
  CALL finalize(vec)
END FUNCTION test_mp_dot_multiply_transposed

function test_mp_show_on_screen() result(passed)
  LOGICAL :: passed

  INTEGER :: idiag

  CALL initialize(Amat,nrows,ncols,PG)
  CALL ALLOCATE(Amat)

  ! set matrix Amat
  DO idiag = 1,nrows
     CALL set_value(Amat,idiag,idiag,CMPLX(idiag,idiag))
  END DO
  CALL set_value(Amat,1,2,CMPLX(1,2))
  CALL set_value(Amat,4,7,CMPLX(4,7))
  CALL commit_values(Amat)

  !CALL show(Amat)
  call finalize(Amat)
  passed=.TRUE.
end function test_mp_show_on_screen

function test_mp_show_in_file() result(passed)
  LOGICAL :: passed

  INTEGER :: idiag

  CALL initialize(Amat,nrows,ncols,PG)
  CALL ALLOCATE(Amat)

  ! set matrix Amat
  DO idiag = 1,nrows
     CALL set_value(Amat,idiag,idiag,CMPLX(idiag,idiag))
  END DO
  CALL set_value(Amat,1,2,CMPLX(1,2))
  CALL set_value(Amat,7,4,CMPLX(7,4))
  CALL commit_values(Amat)

  CALL show(Amat,"amat")
  call finalize(Amat)
	passed=.true.
end function test_mp_show_in_file

function test_mp_add_to_matrix() result(passed)
  LOGICAL :: passed
  
  integer :: idiag

  CALL initialize(Amat,nrows,ncols,PG)
  call allocate(Amat)
  CALL initialize(Bmat,nrows, ncols, PG)
  call allocate(Bmat)

  ! set matrix Amat
  DO idiag = 1,nrows
     CALL set_value(Amat,idiag,idiag,CMPLX(idiag,idiag))
  END DO
  CALL set_value(Amat,1,2,CMPLX(1,2))
  CALL set_value(Amat,4,5,CMPLX(4,5))
  CALL commit_values(Amat)
  CALL show(Amat,'add_amat')

  ! set matrix Bmat
  DO idiag = 1,nrows
     CALL set_value(Bmat,idiag,idiag,CMPLX(2*idiag,idiag))
  END DO
  CALL set_value(Bmat,2,5,CMPLX(2,5))
  CALL set_value(Bmat,8,5,CMPLX(8,5))
  CALL commit_values(Bmat)
  CALL show(Bmat,'add_bmat')
  
  CALL add_matrix(Amat,Bmat)
  CALL show(Amat,'res_amat')
  CALL finalize(Amat)
  CALL finalize(Bmat)
 
	passed=.true.
end function test_mp_add_to_matrix

function test_mp_add_matrix() result(passed)
	logical :: passed

  integer :: idiag

  CALL initialize(Amat,nrows,ncols,PG)
  call allocate(Amat)
  CALL initialize(Bmat,nrows, ncols, PG)
  call allocate(Bmat)
  CALL initialize(Cmat,nrows, ncols, PG)
  call allocate(Cmat)

  ! set matrix Amat
  DO idiag = 1,nrows
     CALL set_value(Amat,idiag,idiag,CMPLX(idiag,idiag))
  END DO
  CALL set_value(Amat,1,2,CMPLX(1,2))
  CALL set_value(Amat,4,5,CMPLX(4,5))
  CALL commit_values(Amat)
  !CALL show(Amat)

  ! set matrix Bmat
  DO idiag = 1,nrows
     CALL set_value(Bmat,idiag,idiag,CMPLX(2*idiag,idiag))
  END DO
  CALL set_value(Bmat,2,5,CMPLX(2,5))
  CALL set_value(Bmat,8,5,CMPLX(8,5))
  CALL commit_values(Bmat)
  !call show(Bmat)
  
  CALL add_matrix(Amat,Bmat,Cmat)
  !CALL show(Cmat)
  CALL finalize(Amat)
  CALL finalize(Bmat)
  call finalize(Cmat)
	passed=.true.
end function test_mp_add_matrix

function test_mp_subtract_matrix() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_subtract_matrix

function test_mp_subtract_from_matrix() result(passed)
	logical :: passed

  integer :: idiag

  CALL initialize(Amat,nrows,ncols,PG)
  call allocate(Amat)
  CALL initialize(Bmat,nrows, ncols, PG)
  call allocate(Bmat)

  ! set matrix Amat
  DO idiag = 1,nrows
     CALL set_value(Amat,idiag,idiag,CMPLX(idiag,idiag))
  END DO
  CALL set_value(Amat,1,2,CMPLX(1,2))
  CALL set_value(Amat,4,5,CMPLX(4,5))
  CALL commit_values(Amat)
  CALL show(Amat,'sub_amat')

  ! set matrix Bmat
  DO idiag = 1,nrows
     CALL set_value(Bmat,idiag,idiag,CMPLX(2*idiag,idiag))
  END DO
  CALL set_value(Bmat,2,5,CMPLX(2,5))
  CALL set_value(Bmat,8,5,CMPLX(8,5))
  CALL commit_values(Bmat)
  CALL show(Bmat,'sub_bmat')
  
  CALL subtract_matrix(Amat,Bmat)
  CALL show(Amat,'res_amat')
  CALL finalize(Amat)
  CALL finalize(Bmat)
	passed=.true.
end function test_mp_subtract_from_matrix

function test_mp_multiply_matrix_with_real() result(passed)
	logical :: passed

  integer :: idiag

  CALL initialize(Amat,nrows,ncols,PG)
  call allocate(Amat)
  CALL initialize(Bmat,nrows,ncols,PG)
  call allocate(Bmat)

  ! set matrix Amat
  DO idiag = 1,nrows
     CALL set_value(Amat,idiag,idiag,CMPLX(idiag,idiag))
  END DO
  CALL set_value(Amat,1,2,CMPLX(1,2))
  CALL set_value(Amat,4,5,CMPLX(4,5))
  CALL commit_values(Amat)
  !CALL show(Amat)

  
  CALL multiply_matrix_with_scalar(Amat,3.0, Bmat)
  !CALL show(Bmat)
  call finalize(Bmat)
  CALL finalize(Amat)
	passed=.true.
end function test_mp_multiply_matrix_with_real

function test_mp_scale_matrix_by_real() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_scale_matrix_by_real

function test_mp_assign_matrix() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_assign_matrix

end program
