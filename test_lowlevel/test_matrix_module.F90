PROGRAM test_matrix
  use ProcessGridModule
  USE MatrixModule
  use mpi
  IMPLICIT NONE
  
  INTEGER :: ierr,n_procs,irow,icol,nrows,ncols, rank, iProc, iCounter
  INTEGER,PARAMETER :: ny0=16
  TYPE(Matrix) :: Amat, vec, res, Bmat,Cmat
  COMPLEX, DIMENSION(:,:), ALLOCATABLE, TARGET :: localarray, localresult
  COMPLEX, DIMENSION(:,:), ALLOCATABLE :: cntrl_A, cntrl_B, cntrl_C, globalMatrixlocally
  COMPLEX :: cmat_val
  LOGICAL :: differences_found, mat_output=.false.

  call mpi_init(ierr)
  call perfinit()
  call perfon('main')
  CALL initialize_matrix_module(MAT_BLOCKS_OF_ROWS,MPI_COMM_WORLD)

  CALL mpi_comm_size(MPI_COMM_WORLD,n_procs, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD,rank, ierr)
  PRINT*,"Using ",n_procs," processors."
  nrows = 4; ncols = 4
  ! allocate the control matrices
  ALLOCATE(cntrl_A(nrows,ncols),cntrl_B(nrows,ncols),cntrl_c(nrows,ncols))
  ALLOCATE(globalMatrixlocally(nrows,ncols))
  cntrl_A = CMPLX(0,0)
  cntrl_B = CMPLX(0,0)
  cntrl_C = CMPLX(0,0)
  globalMatrixlocally = CMPLX(0,0)

  ALLOCATE(localarray(ncols/n_procs,ny0))
  CALL initialize(vec,ncols,1)
  CALL attach(vec,localarray(:,1))
  call set_zero(vec)
  !call show(vec)

  ALLOCATE(localresult(nrows/n_procs,ny0))
  CALL initialize(res,nrows,1)
  CALL attach(res, localresult(:,1))

  CALL initialize(Amat,nrows,ncols)
  call allocate(Amat)
  call set_zero(Amat)

  CALL initialize(Bmat,nrows,ncols)
  call allocate(Bmat)
  call set_zero(Bmat)

  CALL initialize(Cmat,nrows,ncols)
  call allocate(Cmat)
  call set_zero(Cmat)

  call perfon('set_mat')
  DO irow=1,nrows
     DO icol=1,ncols
        cntrl_A(irow,icol) = CMPLX(irow,icol)
        cntrl_B(irow,icol) = CMPLX(irow)/CMPLX(icol)
        CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
        CALL set_value(Bmat,irow,icol,CMPLX(irow)/cmplx(icol))
     END DO
  END DO
  CALL commit_values(Amat)
  CALL commit_values(Bmat)
  call perfoff
  !call show(Amat)
  !CALL show(Bmat)

  call perfon('set_vec')
  DO irow=1,ncols
     CALL set_value(vec,irow,1,CMPLX(1,0))
  END DO
  CALL commit_values(vec)
  call perfoff
  !CALL show(vec)

  call perfon('mult')
  !DO iCounter=1,1000
     CALL dot_multiply(Amat,vec,res)
     CALL dot_multiply(Amat,Bmat,Cmat)
  !END DO
  call perfoff

  cntrl_C=MATMUL(cntrl_A,cntrl_B)
  !CALL show(res)
  CALL show(Cmat)

  CALL get_global_matrix_locally(Cmat,globalMatrixlocally)

  IF (mat_output) THEN
     IF (rank.EQ.0) THEN
        PRINT*, "cntrl_C = "
        DO irow=1,nrows
           WRITE(*,"(I3,A)",advance='NO') irow,": "
           DO icol=1,ncols
              WRITE(*,"(2ES10.2)",advance='NO') cntrl_C(irow,icol)
           END DO
           WRITE(*,"(A)") ""
        END DO
     END IF
     CALL mpi_barrier(MPI_COMM_WORLD,ierr)
     PRINT*, "globalmatrixlocally = "
     DO irow=1,nrows
        WRITE(*,"(I3,A)",advance='NO') irow,": "
        DO icol=1,ncols
           WRITE(*,"(2ES10.2)",advance='NO') globalMatrixlocally(irow,icol)
        END DO
        WRITE(*,"(A)") ""
     END DO
  END IF

#if 1
  DO iProc=0,n_procs-1
     CALL mpi_barrier(MPI_COMM_WORLD,ierr)
     IF (iProc.EQ.rank) THEN
        differences_found = .FALSE.
        DO irow=1,nrows
           DO icol=1,ncols
              !cmat_val = mat_get_value(Cmat,irow,icol)
              cmat_val = globalmatrixlocally(irow,icol)
              IF (cmat_val.NE.cntrl_C(irow,icol)) THEN
                 WRITE(*,"(A,2I3,A,2ES20.10)") "Differences at element ",irow,icol," of ",cmat_val-cntrl_C(irow,icol)
                 differences_found = .TRUE.
              END IF
           END DO
        END DO
        IF (.NOT.differences_found) THEN
           PRINT*,rank,": No differences found, TEST OK!"
        ELSE
           PRINT*,rank,": Differences found, TEST FAILED!"
        END IF
     END IF
  END DO
#endif
  
  CALL finalize(Amat)
  call finalize(Bmat)
  CALL finalize(Cmat)
  CALL finalize(vec)
  CALL finalize(res)
  call finalize()
  call perfoff

  DO iProc=0,n_procs-1
     CALL mpi_barrier(MPI_COMM_WORLD,ierr)
     IF (iProc.EQ.rank) THEN
        PRINT*,"MPI-Rank is ",rank
        CALL perfout('main')
     END IF
  END DO

  call mpi_finalize(ierr)
END PROGRAM test_matrix
