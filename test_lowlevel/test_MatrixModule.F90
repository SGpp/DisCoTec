PROGRAM test_MatrixModule
  USE MatrixModule
  use VectorModule
  use Matrix_Vector_module
  use mpi
  IMPLICIT NONE

  TYPE(Matrix) :: Amat, Bmat,Cmat
  TYPE(Matrix) :: vec, res
  INTEGER :: n_procs, rank, ierr, nrows,ncols, iProc, comm_cart
  COMPLEX, DIMENSION(:,:), ALLOCATABLE :: localfullmat,localfullvec
  REAL :: tolerance = 1e-14

  call mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,n_procs, ierr)
  CALL mpi_cart_create(MPI_COMM_WORLD, 1, (/n_procs/), (/.FALSE./), .FALSE., comm_cart,ierr)

  call initialize_matrix_module(comm_cart)

  CALL mpi_comm_size(comm_cart,n_procs, ierr)
  CALL mpi_comm_rank(comm_cart,rank, ierr)
  PRINT*,"Using ",n_procs," processors."
  nrows = 8; ncols = 8

  ALLOCATE(localfullmat(nrows,ncols),localfullvec(nrows,1))

  CALL initialize(res, nrows)
  CALL ALLOCATE(res)
  !CALL set_zero(res)

  CALL initialize(vec, nrows)
  CALL ALLOCATE(vec)
  !CALL set_zero(vec)

  CALL initialize(Amat,nrows,ncols)
  call allocate(Amat)
  !call set_zero(Amat)

  CALL initialize(Bmat,nrows,ncols)
  call allocate(Bmat)
  !call set_zero(Bmat)

  CALL initialize(Cmat,nrows,ncols)
  call allocate(Cmat)
  !call set_zero(Cmat)

	if (test_set_value()) then
		write(*,"(I3,A40,A10)") rank,'test_set_value','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_set_value','  FAILED'
	end if

	if (test_mat_get_value()) then
		write(*,"(I3,A40,A10)") rank,'test_mat_get_value','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_mat_get_value','  FAILED'
	end if

	if (test_add_value()) then
		write(*,"(I3,A40,A10)") rank,'test_add_value','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_add_value','  FAILED'
	end if

 if (test_add_to_row()) then
    write(*,"(I3,A40,A10)") rank,'test_add_to_row','passed'
	else
    write(*,"(I3,A40,A10)") rank,'test_add_to_row','  FAILED'
 end if
 
 
	if (test_get_global_matrix_locally()) then
		write(*,"(I3,A40,A10)") rank,'test_get_global_matrix_locally','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_get_global_matrix_locally','  FAILED'
	end if

	if (test_get_local_abs_square_sum()) then
		write(*,"(I3,A40,A10)") rank,'test_get_local_abs_square_sum','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_get_local_abs_square_sum','  FAILED'
	end if

	if (test_dot_multiply()) then
		write(*,"(I3,A40,A10)") rank,'test_dot_multiply','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_dot_multiply','  FAILED'
	end if

	if (test_show()) then
		write(*,"(I3,A40,A10)") rank,'test_show','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_show','  FAILED'
	end if

	if (test_add_matrix()) then
		write(*,"(I3,A40,A10)") rank,'test_add_matrix','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_add_matrix','  FAILED'
	end if

	if (test_subtract_matrix()) then
		WRITE(*,"(I3,A40,A10)") rank,"test_subtract_matrix","passed"
	else
		WRITE(*,"(I3,A40,A10)") rank,"test_subtract_matrix","  FAILED"
	end if

	if (test_multiply_matrix_with_scalar()) then
		WRITE(*,"(I3,A40,A10)") rank,"test_multiply_matrix_with_scalar","passed"
	else
		write(*,"(I3,A40,A10)") rank,"test_multiply_matrix_with_scalar","  FAILED"
	end if

	if (test_ASSIGNMENT_equal()) then
		write(*,"(I3,A40,A10)") rank,"test_ASSIGNMENT_equal","passed"
	else
		write(*,"(I3,A40,A10)") rank,"test_ASSIGNMENT_equal","  FAILED"
	end if
#if 1
        IF (test_invert()) THEN
		write(*,"(I3,A40,A10)") rank,'test_invert','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_invert','  FAILED'
	end if
#endif

	if (test_output_data()) then
		write(*,"(I3,A40,A10)") rank,'test_output_data','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_output_data','  FAILED'
	end if

	if (test_my_sum()) then
		write(*,"(I3,A40,A10)") rank,'test_my_sum','passed'
	else
		write(*,"(I3,A40,A10)") rank,'test_my_sum','  FAILED'
	end if

  CALL finalize(Amat)
  call finalize(Bmat)
  CALL finalize(Cmat)
  CALL finalize(vec)
  CALL finalize(res)
  call finalize_vector_module
  call finalize_matrix_module

  call mpi_finalize(ierr)

contains

function test_initialize() result(passed)
	logical :: passed

	passed=.true.
end function

function test_initialize_matrix_module() result(passed)
	logical :: passed

	passed=.true.
end function

function test_allocate() result(passed)
	logical :: passed

	passed=.true.
end function

function test_attach() result(passed)
	logical :: passed

	passed=.true.
end function

function test_finalize() result(passed)
	logical :: passed

	passed=.true.
end function

function test_set_value() result(passed)
	logical :: passed

        INTEGER :: irow,icol
        DO irow=1,nrows
           DO icol=1,ncols
              CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
           END DO
        END DO

        !call show(Amat)
        
	passed=.true.
end function

function test_add_value() result(passed)
	logical :: passed

        INTEGER :: irow, icol

        call set_zero(Amat)
        DO irow=1,nrows
           DO icol=1,ncols
              CALL add_value(Amat,irow,icol,CMPLX(1,2))
           END DO
        END DO

        CALL get_global_matrix_locally(Amat,localfullmat)
	passed=.true.
        DO irow=1,nrows
           DO icol=1,ncols
              IF (ABS(localfullmat(irow,icol)-CMPLX(1,2)).gt.tolerance) passed=.FALSE.
           END DO
        END DO
end function

function test_add_to_row() result(passed)
	logical :: passed

        INTEGER :: irow, icol, which_row
        complex,dimension(:),allocatable :: whole_row

        call set_zero(Amat)
        DO irow=1,nrows
           DO icol=1,ncols
              CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
           END DO
        END DO
        call commit_values(Amat)
        call show(Amat)

        allocate(whole_row(1:ncols))
        which_row = 7
        whole_row=cmplx(1.0,2.0)
        !print*,whole_row
        call add_to_row(Amat,which_row,whole_row)
        call show(Amat)

        CALL get_global_matrix_locally(Amat,localfullmat)
	passed=.true.
        DO irow=1,nrows
           if (irow.ne.which_row) then
              DO icol=1,ncols
                 IF (ABS(localfullmat(irow,icol)-CMPLX(irow,icol)).gt.tolerance) passed=.FALSE.
              END DO
           else
              do icol=1,ncols
                 IF (ABS(localfullmat(irow,icol)-CMPLX(irow,icol)-cmplx(1.0,2.0)).gt.tolerance) passed=.FALSE.
              END DO
           end if
        END DO
        deallocate(whole_row)
end function

function test_commit_values() result(passed)
	logical :: passed

	passed=.true.
end function

function test_set_zero() result(passed)
	logical :: passed

	passed=.true.
end function

function test_get_global_matrix_locally() result(passed)
	logical :: passed

	passed=.true.
end function

function test_get_local_abs_square_sum() result(passed)
	logical :: passed

	passed=.true.
end function

FUNCTION test_mat_get_value() RESULT(passed)
	logical :: passed

        INTEGER :: RowsPerProcess,irow,icol

        call set_zero(Amat)
        DO irow=1,nrows
           DO icol=1,ncols
              CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
           END DO
        END DO
        !CALL show(Amat)

        RowsPerProcess = nrows/n_procs
        passed=.true.
        do irow=1,nrows
           if (rank.eq.((irow-1)/RowsPerProcess)) then
              DO icol=1,ncols
                 IF (ABS(mat_get_value(Amat,irow,icol)-CMPLX(irow,icol)).GT.tolerance) THEN
                    passed=.FALSE.
                    PRINT*,rank,irow,icol,mat_get_value(Amat,irow,icol)
                 END IF
              END DO
           END IF
        END DO
END FUNCTION test_mat_get_value

function test_dot_multiply() result(passed)
	logical :: passed

        INTEGER :: irow,icol
        COMPLEX :: prefac
        COMPLEX,PARAMETER :: imag=CMPLX(0,1)
        COMPLEX, DIMENSION(1:nrows/n_procs) :: locarr
        type(Matrix) :: locvec

        call set_zero(Amat)
        call set_zero(vec)
        call set_zero(res)

        CALL initialize(locvec,nrows)
        CALL attach(locvec,locarr)
        call set_zero(locvec)
        !CALL show(locvec)

        prefac = CMPLX(1,0)
        DO irow=1,nrows
           DO icol=1,ncols
              CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
           END DO
           CALL set_value(vec,irow,1,prefac)
        END DO
        !call show(Amat)

        CALL dot_multiply(Amat,vec,locvec)
        !CALL show(locvec)
        CALL get_global_matrix_locally(locvec,localfullvec)
	passed=.true.
        DO irow=1,nrows
           IF (ABS(localfullvec(irow,1)&
                &-0.5D0*prefac*(2*irow*ncols+imag*(ncols**2+ncols))).GT.tolerance) then
              !WRITE(*,"(I3,4ES15.4)") irow,localfullvec(irow,1),&
              !     &0.5D0*prefac*(2*irow*ncols+imag*(ncols**2+ncols))
              passed=.FALSE.
           END IF
        END DO
        call finalize(locvec)
end function

function test_show() result(passed)
	logical :: passed

	passed=.true.
end function

FUNCTION test_add_matrix() RESULT(passed)
	logical :: passed

        INTEGER :: irow,icol,ierr
        call set_zero(Amat)
        call set_zero(Bmat)

        DO irow=1,nrows
           DO icol=1,ncols
              CALL set_value(Amat,irow,icol,CMPLX(1,2))
              CALL set_value(Bmat,irow,icol,CMPLX(7,4))
           END DO
        END DO
        CALL add_matrix(Amat,Bmat,Cmat)
        CALL get_global_matrix_locally(Cmat,localfullmat)
	passed=.true.
        DO irow=1,nrows
           DO icol=1,ncols
              IF (localfullmat(irow,icol).NE.CMPLX(8,6)) passed=.false.
           END DO
        END DO

        CALL mpi_barrier(MPI_COMM_WORLD,ierr)
        CALL add_matrix(Amat,Bmat)
        CALL get_global_matrix_locally(Amat,localfullmat)
        DO irow=1,nrows
           DO icol=1,ncols
              IF (localfullmat(irow,icol).NE.CMPLX(8,6)) passed=.FALSE.
           END DO
        END DO
        
END FUNCTION test_add_matrix

function test_subtract_matrix() result(passed)
	logical :: passed

        INTEGER :: irow,icol

        call set_zero(Amat)
        call set_zero(Bmat)
        call set_zero(Cmat)

        DO irow=1,nrows
           DO icol=1,ncols
              CALL set_value(Amat,irow,icol,CMPLX(1,2))
              CALL set_value(Bmat,irow,icol,CMPLX(7,4))
           END DO
        END DO
        CALL subtract_matrix(Amat,Bmat,Cmat)

        CALL get_global_matrix_locally(Cmat,localfullmat)
	passed=.true.
        DO irow=1,nrows
           DO icol=1,ncols
              IF (localfullmat(irow,icol).NE.CMPLX(-6,-2)) passed=.FALSE.
           END DO
        END DO

        CALL subtract_matrix(Amat,Bmat)

        CALL get_global_matrix_locally(Amat,localfullmat)
	passed=.true.
        DO irow=1,nrows
           DO icol=1,ncols
              IF (localfullmat(irow,icol).NE.CMPLX(-6,-2)) passed=.FALSE.
           END DO
        END DO
end function

function test_multiply_matrix_with_scalar() result(passed)
	logical :: passed

        real :: factor
        INTEGER :: icol,irow

        factor = 2.0

        call set_zero(Amat)
        call set_zero(Cmat)

        DO irow=1,nrows
           DO icol=1,ncols
              CALL set_value(Amat,irow,icol,CMPLX(1,2))
           END DO
        END DO
        CALL multiply_matrix_with_scalar(Amat,factor,Cmat)

        CALL get_global_matrix_locally(Cmat,localfullmat)
	passed=.true.
        DO irow=1,nrows
           DO icol=1,ncols
              IF (localfullmat(irow,icol).NE.CMPLX(2,4)) passed=.FALSE.
           END DO
        END DO

        CALL multiply_matrix_with_scalar(Amat,factor)

        CALL get_global_matrix_locally(Amat,localfullmat)
	passed=.true.
        DO irow=1,nrows
           DO icol=1,ncols
              IF (localfullmat(irow,icol).NE.CMPLX(2,4)) passed=.FALSE.
           END DO
        END DO

end function

function test_ASSIGNMENT_equal() result(passed)
	logical :: passed

	passed=.true.
end function

function test_invert() result(passed)
	logical :: passed

        INTEGER :: irow,icol

        call set_zero(Amat)
        call set_zero(Bmat)
        call set_zero(Cmat)

        DO irow=1,nrows
           DO icol=1,ncols
              CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
           END DO
           CALL add_value(Amat,irow,irow,CMPLX(irow*irow,0))
        END DO
        !CALL set_value(Amat,1,1,12.0)
        Bmat = Amat

        call invert(Amat)
        CALL dot_multiply(Amat,Bmat,Cmat)
        
        CALL get_global_matrix_locally(Cmat,localfullmat)
	passed=.true.
        DO irow=1,nrows
           DO icol=1,ncols
              IF (irow.EQ.icol) THEN
                 IF (ABS(localfullmat(irow,icol)-CMPLX(1,0)).GT.1e-14) THEN
                    passed=.FALSE.
                    PRINT*,irow,icol,localfullmat(irow,icol)
                 END IF
              ELSE
                 IF (ABS(localfullmat(irow,icol)).GT.1e-14) THEN
                    passed=.FALSE.
                    PRINT*,irow,icol,localfullmat(irow,icol)
                 END IF
              END IF
           END DO
        END DO
end function

function test_output_data() result(passed)
	logical :: passed

	passed=.true.
end function

function test_my_sum() result(passed)
	logical :: passed

        INTEGER :: ierr, comm_cart, prow,pcol, crank,irow,icol,comm_x,comm_else
        INTEGER :: n_procs_x,my_pex,global_rank

        IF (n_procs.GE.4) THEN

           ! now finalize all matrices on the old processgrid, and finalize the
           ! processgrid, as we need a new one here
           CALL finalize(Amat)
           CALL finalize(Bmat)
           CALL finalize(Cmat)
           CALL finalize(vec)
           CALL finalize(res)
           CALL finalize_Vector_module()
           call finalize_Matrix_module

           ! split the comm_world communicator in at least to directions
           CALL mpi_cart_create(MPI_COMM_WORLD,2,(/2,n_procs/2/),&
                &(/.FALSE.,.FALSE./),.TRUE.,comm_cart,ierr)

           CALL mpi_cart_sub(comm_cart,(/.TRUE.,.FALSE./),comm_x,ierr)
           CALL mpi_cart_sub(comm_cart,(/.FALSE.,.TRUE./),comm_else,ierr)

           !DO prow=0,1
           !   DO pcol=0,(n_procs/2)-1
           !      CALL mpi_cart_rank(comm_cart,(/prow,pcol/),crank,ierr)
           !      PRINT*,"(",prow,pcol,"), crank = ",crank
           !   END DO
           !END DO

           
           CALL mpi_comm_size(comm_x,n_procs_x, ierr)
           CALL mpi_comm_rank(comm_x,my_pex, ierr)
           CALL mpi_comm_rank(MPI_COMM_WORLD,global_rank,ierr)
           !PRINT*,"rank = ",rank,", global_rank = ",global_rank,&
           !     &", my_pex = ",my_pex,", n_procs_x = ",n_procs_x

           ! initialize a new process grid
           CALL initialize_matrix_module(comm_x)

           CALL initialize(Amat,nrows,ncols)
           CALL ALLOCATE(Amat)

           !CALL initialize(Bmat,nrows,ncols)
           !CALL ALLOCATE(Bmat)

           !CALL initialize(Cmat,nrows,ncols)
           !CALL ALLOCATE(Cmat)

           DO irow=1,nrows
              DO icol=1,ncols
                 CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
              END DO
           END DO

           !call show(Amat)
           CALL my_sum(Amat,1,comm_else)

           !call show(Amat)
           CALL get_global_matrix_locally(Amat,localfullmat)
           passed=.TRUE.
           DO irow=1,nrows
              DO icol=1,ncols
                 !WRITE(*,"(2ES9.2,X)",advance='NO') localfullmat(irow,icol)
                 IF (ABS(localfullmat(irow,icol)-(n_procs/2)*CMPLX(irow,icol)).GT.tolerance) THEN
                    !PRINT*,global_rank,irow,icol,localfullmat(irow,icol)
                    passed=.FALSE.
                 END IF
              END DO
              !WRITE(*,"(A)") ""
           END DO
        ELSE
           PRINT*, "For testing the my_sum routine, we need at least 4 processes."
           passed=.FALSE.
        END IF
end function
end program
