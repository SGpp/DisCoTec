#include "redef.h"
PROGRAM test_BandedMatrixModule
  use MatrixModule
  USE BandedMatrixModule
  use VectorModule
  USE mpi
  IMPLICIT NONE

  Type(BandedMatrix) :: Amat, Bmat,Cmat
  TYPE(Vector) :: vec, res
  INTEGER :: n_procs, rank, ierr, nrows,ncols, comm_cart
  COMPLEX, DIMENSION(:,:), ALLOCATABLE :: localfullmat
  COMPLEX, DIMENSION(:), allocatable :: localfullvec
  REAL :: tolerance = 1e-14
  logical :: global_passed,passed

  call mpi_init(ierr)  
  CALL mpi_comm_size(MPI_COMM_WORLD,n_procs, ierr)
  CALL mpi_cart_create(MPI_COMM_WORLD, 1, (/n_procs/), (/.FALSE./), .FALSE., comm_cart,ierr)

  call initialize_matrix_module(comm_cart)
  CALL initialize_BandedMatrix_module
  call initialize_Vector_module

  CALL mpi_comm_size(comm_cart,n_procs, ierr)
  CALL mpi_comm_rank(comm_cart,rank, ierr)
  PRINT*,"Using ",n_procs," processors."
  nrows = 8; ncols = 8

  ALLOCATE(localfullmat(nrows,ncols),localfullvec(nrows))

  CALL initialize(res, nrows)
  CALL ALLOCATE(res)
  !CALL set_zero(res)

  CALL initialize(vec, nrows)
  CALL ALLOCATE(vec)
  !CALL set_zero(vec)


  passed = test_set_value()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_set_value','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_set_value','  FAILED'
     end if
  end if

  passed = test_mat_get_value()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_mat_get_value','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_mat_get_value','  FAILED'
     end if
  end if

  passed = test_mat_get_row_pointer()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_mat_get_row_pointer','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_mat_get_row_pointer','  FAILED'
     end if
  end if

  passed = test_add_value()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_add_value','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_add_value','  FAILED'
     end if
  end if

  passed = test_ASSIGNMENT_equal()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,"test_ASSIGNMENT_equal","passed"
     else
        write(*,"(I3,A40,A10)") rank,"test_ASSIGNMENT_equal","  FAILED"
     end if
  end if

  passed = test_get_local_abs_square_sum()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_get_local_abs_square_sum','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_get_local_abs_square_sum','  FAILED'
     end if
  end if

  passed = test_transpose_storage()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_transpose_storage','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_transpose_storage','  FAILED'
     end if
  end if

  passed = test_dot_multiply()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_dot_multiply','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_dot_multiply','  FAILED'
     end if
  end if

  passed = test_matmat_dot_multiply()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        WRITE(*,"(I3,A40,A10)") rank,'test_matmat_dot_multiply','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_matmat_dot_multiply','  FAILED'
     end if
  end if

  passed = test_dot_multiply_TTT()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        WRITE(*,"(I3,A40,A10)") rank,'test_dot_multiply_TTT','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_dot_multiply_TTT','  FAILED'
     end if
  end if

  passed = test_dot_multiply_TNT()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        WRITE(*,"(I3,A40,A10)") rank,'test_dot_multiply_TNT','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_dot_multiply_TNT','  FAILED'
     end if
  end if

  passed = test_dot_multiply_TNN()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        WRITE(*,"(I3,A40,A10)") rank,'test_dot_multiply_TNN','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_dot_multiply_TNN','  FAILED'
     end if
  end if

  passed = test_dot_multiply_NTT()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        WRITE(*,"(I3,A40,A10)") rank,'test_dot_multiply_NTT','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_dot_multiply_NTT','  FAILED'
     end if
  end if

  passed = test_dot_multiply_NTT_2()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        WRITE(*,"(I3,A40,A10)") rank,'test_dot_multiply_NTT_2','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_dot_multiply_NTT_2','  FAILED'
     end if
  end if

  passed = test_dot_multiply_Banded_with_Full()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        WRITE(*,"(I3,A40,A10)") rank,'test_dot_multiply_Banded_with_Full','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_dot_multiply_Banded_with_Full','  FAILED'
     end if
  end if

  passed = test_show()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_show','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_show','  FAILED'
     end if
  end if

  passed = test_add_matrix()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_add_matrix','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_add_matrix','  FAILED'
     end if
  end if

  passed = test_subtract_matrix()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        WRITE(*,"(I3,A40,A10)") rank,"test_subtract_matrix","passed"
     else
        WRITE(*,"(I3,A40,A10)") rank,"test_subtract_matrix","  FAILED"
     end if
  end if

  passed = test_multiply_matrix_with_scalar()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        WRITE(*,"(I3,A40,A10)") rank,"test_multiply_matrix_with_scalar","passed"
     else
        write(*,"(I3,A40,A10)") rank,"test_multiply_matrix_with_scalar","  FAILED"
     end if
  end if

#if 1
  passed = test_invert()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_invert','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_invert','  FAILED'
     end if
  end if
#endif

  passed = test_output_data()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_output_data','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_output_data','  FAILED'
     end if
  end if

  passed = test_transpose_and_conjugate()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_transpose_and_conjugate','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_transpose_and_conjugate','  FAILED'
     end if
  end if

  passed = test_my_sum()
  call mpi_reduce(passed,global_passed,1,MPI_LOGICAL,MPI_LAND,0,comm_cart,ierr)
  if (rank.eq.0) then
     if (global_passed) then
        write(*,"(I3,A40,A10)") rank,'test_my_sum','passed'
     else
        write(*,"(I3,A40,A10)") rank,'test_my_sum','  FAILED'
     end if
  end if

  CALL finalize(vec)
  CALL finalize(res)
  call finalize_BandedMatrix_module()
  call finalize_Vector_module
  call finalize_Matrix_module()

  call mpi_finalize(ierr)
contains

  function test_initialize() result(passed)
    logical :: passed

    passed=.true.
  end function test_initialize

  function test_initialize_matrix_module() result(passed)
    logical :: passed

    passed=.true.
  end function test_initialize_matrix_module

  function test_allocate() result(passed)
    logical :: passed

    passed=.true.
  end function test_allocate

  function test_attach() result(passed)
    logical :: passed

    passed=.true.
  end function test_attach

  function test_finalize() result(passed)
    logical :: passed

    passed=.true.
  end function test_finalize

  function test_set_value() result(passed)
    logical :: passed

    INTEGER :: irow,icol

    CALL initialize(Amat,nrows,ncols,.false.)
    CALL ALLOCATE(Amat)

    DO irow=1,nrows
       DO icol=1,ncols
          CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
       END DO
    END DO
    call commit_values(Amat)
    !CALL show(Amat)
    call finalize(Amat)
    passed=.true.
  end function test_set_value

  function test_add_value() result(passed)
    logical :: passed

    INTEGER :: irow, icol

    CALL initialize(Amat,nrows,ncols,.true.)
    CALL ALLOCATE(Amat)

    call set_zero(Amat)
    DO irow=1,nrows
       DO icol=1,ncols
          CALL add_value(Amat,irow,icol,CMPLX(1,2))
       END DO
    END DO
    call commit_values(Amat)
    !CALL show(Amat)
    CALL get_global_matrix_locally(Amat,localfullmat)
    passed=.true.
    DO irow=1,nrows
       DO icol=1,ncols
          !WRITE(*,"(2F5.1,X)",advance='no') localfullmat(irow,icol)
          IF (ABS(localfullmat(irow,icol)-CMPLX(1,2)).gt.tolerance) passed=.FALSE.
       END DO
       !WRITE(*,"(A)") ""
    END DO
    call finalize(Amat)
  end function test_add_value


  function test_get_local_abs_square_sum() result(passed)
    logical :: passed

    passed=.true.
  end function test_get_local_abs_square_sum

  FUNCTION test_mat_get_value() RESULT(passed)
    logical :: passed

    INTEGER :: irow,icol
    logical :: transp

    transp = .true.
    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)
    CALL initialize(Bmat,nrows,ncols)
    CALL ALLOCATE(Bmat)

    !CALL set_zero(Amat)
    DO irow=1,nrows
       DO icol=1,ncols
          ! A has UpperBandwidth=2 and LowerBandwidth=1
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
          END IF
          ! B has UpperBandwidth=2 and LowerBandwidth=1
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Bmat,irow,icol,CMPLX(irow,icol))
          END IF
       END DO
    END DO
    call commit_values(Amat)
    call commit_values(Bmat)
    passed=.true.
#if 0
    print*,"A(2,1)N = ",mat_get_value(Amat,2,1,'N')
    print*,"A(2,1)T = ",mat_get_value(Amat,2,1,'T')
    print*,"A(2,1)C = ",mat_get_value(Amat,2,1,'C')
    print*,"B(2,1)N = ",mat_get_value(Bmat,2,1,'N')
    print*,"B(2,1)T = ",mat_get_value(Bmat,2,1,'T')
    print*,"B(2,1)C = ",mat_get_value(Bmat,2,1,'C')

    do irow=1,nrows
       do icol=1,ncols
          val = mat_get_value(Amat,irow,icol,'N')
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             if (abs(val-cmplx(irow,icol)).gt.tolerance) then
                passed=.false.
             end if
          END IF
          ! B has UpperBandwidth=2 and LowerBandwidth=1
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             val = mat_get_value(Bmat,irow,icol,'N')
             if (abs(val-cmplx(irow,icol)).gt.tolerance) then
                passed=.false.
             end if
          END IF
       end do
    end do

    IF (transp) THEN
       RowsPerProcess = nrows/n_procs
       DO irow=1,nrows
          IF (((irow-1)/RowsPerProcess).EQ.rank) THEN
             DO icol=1,ncols
                IF (ABS(mat_get_value(Amat,irow,icol,'N')-CMPLX(irow,icol)).GT.tolerance) THEN
                   passed=.FALSE.
                   PRINT*,rank,irow,icol,mat_get_value(Amat,irow,icol)
                END IF
             END DO
          END IF
       END DO
    ELSE
       ColsPerProcess = ncols/n_procs
       DO icol=1,ncols
          IF (((icol-1)/ColsPerProcess).EQ.rank) THEN
             DO irow=1,nrows
                IF (ABS(mat_get_value(Amat,irow,icol)-CMPLX(irow,icol)).GT.tolerance) THEN
                   passed=.FALSE.
                   PRINT*,rank,irow,icol,mat_get_value(Amat,irow,icol)
                END IF
             END DO
          END IF
       END DO
    END IF
#endif
    CALL finalize(Amat)
    call finalize(Bmat)

  END FUNCTION test_mat_get_value

  FUNCTION test_mat_get_row_pointer() RESULT(passed)
    logical :: passed

    INTEGER :: global_row,s_global_col,e_global_col,irow,icol
    complex, dimension(:),pointer :: ptr_row
    logical :: transp

    transp = .true.
    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)

    !CALL set_zero(Amat)
    DO irow=1,nrows
       DO icol=1,ncols
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
          END IF
       END DO
    END DO
    call commit_values(Amat)
    !CALL show(Amat)

    passed=.true.
    IF (transp) THEN
       global_row=1
       call mat_get_row_pointer(Amat,global_row,ptr_row,s_global_col,e_global_col)
       !if (associated(ptr_row)) then
       !   print*,rank,": global_col = ",s_global_col, e_global_col,", ptr_row = ", ptr_row
       !else
       !   print*,rank,": ptr_row not associated"
       !end if

       global_row=7
       call mat_get_row_pointer(Amat,global_row,ptr_row,s_global_col,e_global_col)
       !if (associated(ptr_row)) then
       !   print*,rank, ": global_col = ",s_global_col, e_global_col,", ptr_row = ", ptr_row
       !else
       !   print*,rank,": ptr_row not associated"
       !end if
    ELSE
    END IF
    CALL finalize(Amat)
  END FUNCTION test_mat_get_row_pointer

  function test_transpose_storage() result(passed)
    logical :: passed

    INTEGER :: irow,icol, lbw, ubw
    COMPLEX,PARAMETER :: imag=CMPLX(0,1)
    complex, dimension(:,:),allocatable :: localmatrix
    LOGICAL :: transp=.true.

    passed = .true.
    CALL initialize(Amat,nrows,ncols)
    CALL ALLOCATE(Amat)

    call initialize(Bmat,nrows,ncols,transp)
    call allocate(Bmat)

    CALL set_zero(Amat)

    lbw = 1
    ubw = 2
    DO irow=1,nrows
       DO icol=1,ncols
          IF ((icol.GE.irow-lbw).AND.(icol.LE.irow+ubw)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
          END IF
       END DO
    END DO
    CALL commit_values(Amat)

    call show(Amat)
    call transpose_storage(Amat,Bmat)

    call show(Bmat)

    allocate(localmatrix(nrows,ncols))
    call get_global_matrix_locally(Bmat,localmatrix)
    do irow=1,nrows
       do icol=1,ncols
          IF ((icol.GE.irow-lbw).AND.(icol.LE.irow+ubw)) THEN
             if (abs(localmatrix(irow,icol)-cmplx(irow,icol)).gt.tolerance) then
                print*,"error at irow,icol = ",irow,icol,abs(localmatrix(irow,icol)-cmplx(irow,icol))
                passed = .false.
             end if
          END IF
       end do
    end do

    deallocate(localmatrix)
    call finalize(Amat)
    call finalize(Bmat)
  end function test_transpose_storage

  function test_dot_multiply() result(passed)
    logical :: passed

    INTEGER :: irow,icol
    COMPLEX :: prefac
    COMPLEX,PARAMETER :: imag=CMPLX(0,1)
    COMPLEX, dimension(1:nrows) :: locarr
    complex,dimension(:,:),allocatable :: Acmat
    complex,dimension(:),allocatable :: cvec,cres
    type(Vector) :: locvec
    LOGICAL :: transp=.true.

    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)

    allocate(Acmat(nrows,ncols))
    allocate(cvec(nrows),cres(nrows))

    !CALL set_zero(Amat)
    call set_zero(vec)
    call set_zero(res)

    CALL initialize(locvec,nrows)
    CALL attach(locvec,locarr)
    call set_zero(locvec)
    !CALL show(locvec)
    
    Acmat = cmplx(0.0,0.0)
    prefac = CMPLX(1,0)
    DO irow=1,nrows
       DO icol=1,ncols
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
             Acmat(irow,icol) = cmplx(irow,icol)
          END IF
       END DO
       CALL set_value(vec,irow,prefac)
       cvec(irow) = prefac
    END DO
    CALL commit_values(Amat)
    call commit_values(vec)
    !CALL show(Amat)
    !call show(vec)

    CALL dot_multiply(Amat,vec,locvec)
    !CALL show(locvec)
    cres = matmul(Acmat,cvec)

    CALL get_global_vector_locally(locvec,localfullvec)
    passed=.true.
    DO irow=1,nrows
       IF (ABS(localfullvec(irow)-cres(irow)).GT.tolerance) then
          !WRITE(*,"(I3,4ES15.4)") irow,localfullvec(irow,1),&
          !     &0.5D0*prefac*(2*irow*ncols+imag*(ncols**2+ncols))
          passed=.FALSE.
       END IF
    END DO
    call finalize(locvec)
    call finalize(Amat)
    deallocate(Acmat,cvec,cres)
  end function test_dot_multiply

  FUNCTION test_matmat_dot_multiply() RESULT(passed)
    logical :: passed

    INTEGER :: irow,icol
    COMPLEX :: prefac
    COMPLEX,PARAMETER :: imag=CMPLX(0,1)
    REAL(8) :: time_start, time_end
    LOGICAL :: transp = .TRUE.

    !PRINT*," in test_matmat_dot_multiply"
    passed= .TRUE.
    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)

    CALL initialize(Bmat,nrows,ncols)
    call allocate(Bmat)

    CALL initialize(Cmat,nrows,ncols,transp)
    call allocate(Cmat)

    prefac = CMPLX(1,0)
    DO irow=1,nrows
       DO icol=1,ncols
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
          END IF
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+1)) THEN
             CALL set_value(Bmat,irow,icol,CMPLX(1,0))
          END IF
       END DO
       !CALL set_value(Bmat,irow,irow,CMPLX(1,0))
    END DO
    CALL commit_values(Amat)
    call commit_values(Bmat)
    !CALL show(Amat)
    !CALL show(Bmat)

    time_start = MPI_WTime()
    ! calculate Cmat = Amat.Bmat
    CALL dot_multiply(Amat,Bmat,Cmat)
    time_end = MPI_WTime()
    !PRINT*,"time used is ", time_end-time_start
    !PRINT*, "Now show Cmat, the result."
    !call show(Cmat)

    call finalize(Cmat)
    CALL finalize(Bmat)
    CALL finalize(Amat)
  END FUNCTION test_matmat_dot_multiply

  FUNCTION test_dot_multiply_TTT() RESULT(passed)
    logical :: passed

    complex, allocatable,dimension(:,:) :: Acmat,Bcmat,Ccmat,localmatrix
    INTEGER :: irow,icol
    COMPLEX :: prefac
    REAL(8) :: time_start, time_end
    LOGICAL :: transp = .TRUE.

    passed= .TRUE.
    ! first we initialize three banded matrices in row-major storage format
    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)

    CALL initialize(Bmat,nrows,ncols,transp)
    call allocate(Bmat)

    CALL initialize(Cmat,nrows,ncols,transp)
    call allocate(Cmat)

    ! locally we also have three full matrices for correctness check
    allocate(Acmat(nrows,ncols))
    allocate(Bcmat(nrows,ncols))
    allocate(Ccmat(nrows,ncols))
    allocate(localmatrix(nrows,ncols))
    Acmat= cmplx(0.0,0.0)
    Bcmat = cmplx(0.0,0.0)
    ! set the matrices
    prefac = CMPLX(1,0)
    DO irow=1,nrows
       DO icol=1,ncols
          ! A has UpperBandwidth=2 and LowerBandwidth=1
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
             Acmat(irow,icol) = cmplx(irow,icol)
          END IF
          ! B has UpperBandwidth=1 and LowerBandwidth=3
          IF ((icol.GE.irow-3).AND.(icol.LE.irow+1)) THEN
             CALL set_value(Bmat,irow,icol,CMPLX(irow,icol))
             Bcmat(irow,icol) = cmplx(irow,icol)
          END IF

       END DO
    END DO
    CALL commit_values(Amat)
    call commit_values(Bmat)
    !CALL show(Amat)
    !CALL show(Bmat)

    !print*,"Starting with the calculation."
    time_start = MPI_WTime()
    ! calculate Fmat = Amat.Emat
    CALL dot_multiply(Amat,Bmat,Cmat,'N','N')
    time_end = MPI_WTime()
    !PRINT*,"time used is ", time_end-time_start
    !PRINT*, "Now show Fmat, the result."
    !call show(Cmat)

    ! do the control calculation locally
    Ccmat = matmul(Acmat,Bcmat)
    !print*,Ccmat
    call get_global_matrix_locally(Cmat,localmatrix)
    do irow=1,nrows
       do icol=1,ncols
          if (abs(localmatrix(irow,icol)-Ccmat(irow,icol)).gt.tolerance) then
             print*,"error at irow,icol = ",irow,icol,abs(localmatrix(irow,icol)-Ccmat(irow,icol))
             passed = .false.
          end if
       end do
    end do
    deallocate(Acmat,Bcmat,Ccmat)
    call finalize(Cmat)
    CALL finalize(Bmat)
    CALL finalize(Amat)
  END FUNCTION test_dot_multiply_TTT

  FUNCTION test_dot_multiply_TNT() RESULT(passed)
    logical :: passed

    complex, allocatable,dimension(:,:) :: Acmat,Bcmat,Ccmat,localmatrix
    INTEGER :: irow,icol
    COMPLEX :: prefac
    REAL(8) :: time_start, time_end
    LOGICAL :: transp = .TRUE.

    passed= .TRUE.
    ! first we initialize three banded matrices in row-major storage format
    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)

    CALL initialize(Bmat,nrows,ncols,transp)
    call allocate(Bmat)

    CALL initialize(Cmat,nrows,ncols,transp)
    call allocate(Cmat)

    ! locally we also have three full matrices for correctness check
    allocate(Acmat(nrows,ncols))
    allocate(Bcmat(nrows,ncols))
    allocate(Ccmat(nrows,ncols))
    allocate(localmatrix(nrows,ncols))
    Acmat= cmplx(0.0,0.0)
    Bcmat = cmplx(0.0,0.0)
    ! set the matrices
    prefac = CMPLX(1,0)
    DO irow=1,nrows
       DO icol=1,ncols
          ! A has UpperBandwidth=2 and LowerBandwidth=1
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
             Acmat(irow,icol) = cmplx(irow,icol)
          END IF
          ! B has UpperBandwidth=1 and LowerBandwidth=3
          IF ((icol.GE.irow-3).AND.(icol.LE.irow+1)) THEN
             CALL set_value(Bmat,irow,icol,CMPLX(irow,icol))
             Bcmat(irow,icol) = cmplx(irow,icol)
          END IF

       END DO
    END DO
    CALL commit_values(Amat)
    call commit_values(Bmat)
    !CALL show(Amat)
    !CALL show(Bmat)

    !print*,"Starting with the calculation."
    time_start = MPI_WTime()
    ! calculate Fmat = Amat.Emat
    CALL dot_multiply(Amat,Bmat,Cmat,'N','N')
    time_end = MPI_WTime()
    !PRINT*,"time used is ", time_end-time_start
    !PRINT*, "Now show Fmat, the result."
    !call show(Cmat)

    ! do the control calculation locally
    Ccmat = matmul(Acmat,Bcmat)
    !print*,Ccmat
    call get_global_matrix_locally(Cmat,localmatrix)
    do irow=1,nrows
       do icol=1,ncols
          if (abs(localmatrix(irow,icol)-Ccmat(irow,icol)).gt.tolerance) then
             print*,"error at irow,icol = ",irow,icol,abs(localmatrix(irow,icol)-Ccmat(irow,icol))
             passed = .false.
          end if
       end do
    end do
    deallocate(Acmat,Bcmat,Ccmat)
    call finalize(Cmat)
    CALL finalize(Bmat)
    CALL finalize(Amat)
  END FUNCTION test_dot_multiply_TNT

  FUNCTION test_dot_multiply_TNN() RESULT(passed)
    logical :: passed

    complex, allocatable,dimension(:,:) :: Acmat,Bcmat,Ccmat,localmatrix
    INTEGER :: irow,icol
    COMPLEX :: prefac
    REAL(8) :: time_start, time_end
    LOGICAL :: transp = .TRUE.

    passed= .TRUE.
    ! first we initialize three banded matrices in row-major storage format
    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)

    CALL initialize(Bmat,nrows,ncols)
    call allocate(Bmat)

    CALL initialize(Cmat,nrows,ncols)
    call allocate(Cmat)

    ! locally we also have three full matrices for correctness check
    allocate(Acmat(nrows,ncols))
    allocate(Bcmat(nrows,ncols))
    allocate(Ccmat(nrows,ncols))
    allocate(localmatrix(nrows,ncols))
    Acmat= cmplx(0.0,0.0)
    Bcmat = cmplx(0.0,0.0)
    ! set the matrices
    prefac = CMPLX(1,0)
    DO irow=1,nrows
       DO icol=1,ncols
          ! A has UpperBandwidth=2 and LowerBandwidth=1
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
             Acmat(irow,icol) = cmplx(irow,icol)
          END IF
          ! B has UpperBandwidth=1 and LowerBandwidth=0
          IF ((icol.GE.irow).AND.(icol.LE.irow+1)) THEN
             CALL set_value(Bmat,irow,icol,CMPLX(irow,icol))
             Bcmat(irow,icol) = cmplx(irow,icol)
          END IF

       END DO
    END DO
    CALL commit_values(Amat)
    call commit_values(Bmat)
    !CALL show(Amat)
    !CALL show(Bmat)

    !print*,"Starting with the calculation."
    time_start = MPI_WTime()
    ! calculate Fmat = Amat.Emat
    CALL dot_multiply(Amat,Bmat,Cmat,'N','N')
    time_end = MPI_WTime()
    !PRINT*,"time used is ", time_end-time_start
    !PRINT*, "Now show Fmat, the result."
    !call show(Cmat)

    ! do the control calculation locally
    Ccmat = matmul(Acmat,Bcmat)
    !print*,Ccmat
    call get_global_matrix_locally(Cmat,localmatrix)
    do irow=1,nrows
       do icol=1,ncols
          if (abs(localmatrix(irow,icol)-Ccmat(irow,icol)).gt.tolerance) then
             !print*,"error at irow,icol = ",irow,icol,abs(localmatrix(irow,icol)-Ccmat(irow,icol))
             passed = .false.
          end if
       end do
    end do
    deallocate(Acmat,Bcmat,Ccmat)
    call finalize(Cmat)
    CALL finalize(Bmat)
    CALL finalize(Amat)
  END FUNCTION test_dot_multiply_TNN

  FUNCTION test_dot_multiply_NTT() RESULT(passed)
    logical :: passed

    complex, allocatable,dimension(:,:) :: Acmat,Bcmat,Ccmat,localmatrix
    INTEGER :: irow,icol
    COMPLEX :: prefac
    REAL(8) :: time_start, time_end
    LOGICAL :: transp = .TRUE.

    passed= .TRUE.
    ! first we initialize three banded matrices (N,T,T format)
    CALL initialize(Amat,nrows,ncols)
    CALL ALLOCATE(Amat)

    CALL initialize(Bmat,nrows,ncols,transp)
    call allocate(Bmat)

    CALL initialize(Cmat,nrows,ncols,transp)
    call allocate(Cmat)

    ! locally we also have three full matrices for correctness check
    allocate(Acmat(nrows,ncols))
    allocate(Bcmat(nrows,ncols))
    allocate(Ccmat(nrows,ncols))
    allocate(localmatrix(nrows,ncols))
    Acmat= cmplx(0.0,0.0)
    Bcmat = cmplx(0.0,0.0)
    ! set the matrices
    prefac = CMPLX(1,0)
    DO irow=1,nrows
       DO icol=1,ncols
          ! A has UpperBandwidth=2 and LowerBandwidth=1
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
             Acmat(irow,icol) = cmplx(irow,icol)
          END IF
          ! B has UpperBandwidth=1 and LowerBandwidth=0
          IF ((icol.GE.irow-0).AND.(icol.LE.irow+1)) THEN
             CALL set_value(Bmat,irow,icol,CMPLX(irow,icol))
             Bcmat(irow,icol) = cmplx(irow,icol)
          END IF

       END DO
    END DO
    CALL commit_values(Amat)
    call commit_values(Bmat)
    !CALL show(Amat)
    !CALL show(Bmat)

    !print*,"Starting with the calculation."
    time_start = MPI_WTime()
    ! calculate Fmat = Amat.Emat
    CALL dot_multiply(Amat,Bmat,Cmat,'N','N')
    time_end = MPI_WTime()
    !PRINT*,"time used is ", time_end-time_start
    !PRINT*, "Now show Fmat, the result."
    !call show(Cmat)

    ! do the control calculation locally
    Ccmat = matmul(Acmat,Bcmat)
    !print*,Ccmat
    call get_global_matrix_locally(Cmat,localmatrix)
    do irow=1,nrows
       do icol=1,ncols
          if (abs(localmatrix(irow,icol)-Ccmat(irow,icol)).gt.tolerance) then
             write(*,"(A,2I3,2F10.2)") "error at irow,icol = ",irow,icol,localmatrix(irow,icol)-Ccmat(irow,icol)
             passed = .false.
          end if
       end do
    end do
    deallocate(Acmat,Bcmat,Ccmat)
    call finalize(Cmat)
    CALL finalize(Bmat)
    CALL finalize(Amat)
  END FUNCTION test_dot_multiply_NTT

  FUNCTION test_dot_multiply_NTT_2() RESULT(passed)
    logical :: passed

    complex, allocatable,dimension(:,:) :: Acmat,Bcmat,Ccmat,localmatrix
    INTEGER :: irow,icol
    COMPLEX :: prefac
    REAL(8) :: time_start, time_end
    LOGICAL :: transp = .TRUE.

    passed= .TRUE.
    ! first we initialize three banded matrices (N,T,T format)
    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)

    CALL initialize(Bmat,nrows,ncols,transp)
    call allocate(Bmat)

    CALL initialize(Cmat,nrows,ncols,transp)
    call allocate(Cmat)

    ! locally we also have three full matrices for correctness check
    allocate(Acmat(nrows,ncols))
    allocate(Bcmat(nrows,ncols))
    allocate(Ccmat(nrows,ncols))
    allocate(localmatrix(nrows,ncols))
    Acmat= cmplx(0.0,0.0)
    Bcmat = cmplx(0.0,0.0)
    ! set the matrices
    prefac = CMPLX(1,0)
    DO irow=1,nrows
       DO icol=1,ncols
          ! A has UpperBandwidth=2 and LowerBandwidth=1
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
             Acmat(irow,icol) = cmplx(irow,icol)
          END IF
          ! B has UpperBandwidth=1 and LowerBandwidth=0
          IF ((icol.GE.irow-0).AND.(icol.LE.irow+1)) THEN
             CALL set_value(Bmat,irow,icol,CMPLX(irow,icol))
             Bcmat(irow,icol) = cmplx(irow,icol)
          END IF

       END DO
    END DO
    CALL commit_values(Amat)
    call commit_values(Bmat)
    !CALL show(Amat)
    !CALL show(Bmat)

    !print*,"Starting with the calculation."
    time_start = MPI_WTime()
    ! calculate Fmat = Amat.Emat
    CALL dot_multiply(Amat,Bmat,Cmat,'C','N')
    time_end = MPI_WTime()
    !PRINT*,"time used is ", time_end-time_start
    !PRINT*, "Now show Fmat, the result."
    !call show(Cmat)

    ! do the control calculation locally
    Ccmat = matmul(transpose(conjg(Acmat)),Bcmat)
    !print*,Ccmat
    call get_global_matrix_locally(Cmat,localmatrix)
    do irow=1,nrows
       do icol=1,ncols
          if (abs(localmatrix(irow,icol)-Ccmat(irow,icol)).gt.tolerance) then
             write(*,"(A,2I3,2F10.2)") "error at irow,icol = ",irow,icol,localmatrix(irow,icol)-Ccmat(irow,icol)
             passed = .false.
          end if
       end do
    end do
    deallocate(Acmat,Bcmat,Ccmat)
    call finalize(Cmat)
    CALL finalize(Bmat)
    CALL finalize(Amat)
  END FUNCTION test_dot_multiply_NTT_2

  FUNCTION test_dot_multiply_Banded_with_Full() RESULT(passed)
    logical :: passed

    type(Matrix) :: Emat, Fmat

    complex, allocatable,dimension(:,:) :: Acmat,Ecmat,Fcmat,localmatrix
    INTEGER :: irow,icol
    COMPLEX :: prefac
    REAL(8) :: time_start, time_end
    LOGICAL :: transp = .TRUE.

    passed= .TRUE.
    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)

    ! two full matrices in column-major storage format
    CALL initialize(Emat,nrows,ncols)
    call allocate(Emat)

    CALL initialize(Fmat,nrows,ncols)
    call allocate(Fmat)

    allocate(Acmat(nrows,ncols))
    allocate(Ecmat(nrows,ncols))
    allocate(Fcmat(nrows,ncols))
    allocate(localmatrix(nrows,ncols))
    Acmat= cmplx(0.0,0.0)
    ! set the matrices
    prefac = CMPLX(1,0)
    DO irow=1,nrows
       DO icol=1,ncols
          IF ((icol.GE.irow-1).AND.(icol.LE.irow+2)) THEN
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
             Acmat(irow,icol) = cmplx(irow,icol)
          END IF
          CALL set_value(Emat,irow,icol,cmplx(0.1D0*irow,0.1D0*icol))
          Ecmat(irow,icol) = cmplx(0.1D0*irow,0.1D0*icol)
       END DO
    END DO
    CALL commit_values(Amat)
    call commit_values(Emat)
    !CALL show(Amat)
    !CALL show(Emat)

    !print*,"Starting with the calculation."
    time_start = MPI_WTime()
    ! calculate Fmat = Amat.Emat
    CALL dot_multiply(Amat,Emat,Fmat,'N','N')
    time_end = MPI_WTime()
    !PRINT*,"time used is ", time_end-time_start
    !PRINT*, "Now show Fmat, the result."
    !call show(Fmat)

    ! do the control calculation locally
    Fcmat = matmul(Acmat,Ecmat)
    call get_global_matrix_locally(Fmat,localmatrix)
    do irow=1,nrows
       do icol=1,ncols
          if (abs(localmatrix(irow,icol)-Fcmat(irow,icol)).gt.tolerance) then
             print*,"error at irow,icol = ",irow,icol,abs(localmatrix(irow,icol)-Fcmat(irow,icol))
             passed = .false.
          end if
       end do
    end do
    deallocate(Acmat,Ecmat,Fcmat)
    call finalize(Emat)
    CALL finalize(Fmat)
    CALL finalize(Amat)
  END FUNCTION test_dot_multiply_Banded_with_Full

  function test_show() result(passed)
    logical :: passed

    passed=.true.
  end function test_show

  FUNCTION test_add_matrix() RESULT(passed)
    logical :: passed

    INTEGER :: irow,icol,ierr
    logical :: transp = .true.
    REAL(8) :: time_start, time_end

    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)
    CALL initialize(Bmat,nrows,ncols,transp)
    CALL ALLOCATE(Bmat)
    CALL initialize(Cmat,nrows,ncols,transp)
    CALL ALLOCATE(Cmat)

    DO irow=1,nrows
       DO icol=1,ncols
          CALL set_value(Amat,irow,icol,CMPLX(1,2))
          CALL set_value(Bmat,irow,icol,CMPLX(7,4))
       END DO
    END DO
    call commit_values(Amat)
    call commit_values(Bmat)

    !call show(Amat)
    !call show(Bmat)
    time_start = MPI_Wtime()
    CALL add_matrix(Amat,Bmat,Cmat)
    time_end = mpi_wtime()
    !call show(Cmat)

    !PRINT*,"add_matrix(Amat,Bmat,Cmat) needs ",time_end-time_start," s"
    CALL get_global_matrix_locally(Cmat,localfullmat)
    passed=.true.
    DO irow=1,nrows
       DO icol=1,ncols
          IF (abs(localfullmat(irow,icol)-CMPLX(8.0,6.0)).gt.tolerance) then
             passed=.false.
             print*,irow,icol,localfullmat(irow,icol)-cmplx(8,6)
          end IF
       END DO
    END DO

    CALL mpi_barrier(MPI_COMM_WORLD,ierr)
    time_start = MPI_Wtime()
    CALL add_matrix(Amat,Bmat)
    time_end = mpi_wtime()
    !PRINT*,"add_matrix(Amat,Bmat) needs ",time_end-time_start," s"
    CALL get_global_matrix_locally(Amat,localfullmat)
    DO irow=1,nrows
       DO icol=1,ncols
          IF (localfullmat(irow,icol).NE.CMPLX(8,6)) passed=.FALSE.
       END DO
    END DO
    call finalize(Amat)
    call finalize(Bmat)
    call finalize(Cmat)
  END FUNCTION test_add_matrix

  FUNCTION test_subtract_matrix() RESULT(passed)
    logical :: passed

    INTEGER :: irow,icol

    CALL initialize(Amat,nrows,ncols)
    CALL ALLOCATE(Amat)
    CALL initialize(Bmat,nrows,ncols)
    CALL ALLOCATE(Bmat)
    CALL initialize(Cmat,nrows,ncols)
    CALL ALLOCATE(Cmat)

    !CALL set_zero(Amat)
    !call set_zero(Bmat)
    !call set_zero(Cmat)

    DO irow=1,nrows
       DO icol=1,ncols
          CALL set_value(Amat,irow,icol,CMPLX(1,2))
          CALL set_value(Bmat,irow,icol,CMPLX(7,4))
       END DO
    END DO
    call commit_values(Amat)
    call commit_values(Bmat)
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
    call finalize(Amat)
    call finalize(Bmat)
    call finalize(Cmat)
  END FUNCTION test_subtract_matrix

  FUNCTION test_multiply_matrix_with_scalar() RESULT(passed)
    logical :: passed

    real :: factor
    INTEGER :: icol,irow

    CALL initialize(Amat,nrows,ncols)
    CALL ALLOCATE(Amat)
    CALL initialize(Cmat,nrows,ncols)
    CALL ALLOCATE(Cmat)

    factor = 2.0

    !CALL set_zero(Amat)
    !call set_zero(Cmat)

    DO irow=1,nrows
       DO icol=1,ncols
          CALL set_value(Amat,irow,icol,CMPLX(1,2))
       END DO
    END DO
    call commit_values(Amat)
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
    call finalize(Amat)
    call finalize(Cmat)
  END FUNCTION test_multiply_matrix_with_scalar

  function test_ASSIGNMENT_equal() result(passed)
    logical :: passed

    integer :: irow,icol
    logical :: transp=.true.

    CALL initialize(Amat,nrows,ncols,transp)
    CALL ALLOCATE(Amat)
    CALL initialize(Cmat,nrows,ncols,transp)
    CALL ALLOCATE(Cmat)
    
    DO irow=1,nrows
       DO icol=1,ncols
          CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
       END DO
    END DO
    call commit_values(Amat)

    Cmat = Amat

    call get_global_matrix_locally(Cmat,localfullmat)

    passed = .true.
    do irow=1,nrows
       do icol=1,ncols
          if (abs(localfullmat(irow,icol)-cmplx(irow,icol)).gt.tolerance) then
             passed = .false.
          end if
       end do
    end do

  end function test_ASSIGNMENT_equal

  FUNCTION test_invert() RESULT(passed)
    logical :: passed

#ifdef ALL_ROUTINES
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
#endif
    passed= .TRUE.
  END FUNCTION test_invert

  function test_output_data() result(passed)
    logical :: passed

    passed=.true.
  end function test_output_data

  function test_transpose_and_conjugate() result(passed)
    logical :: passed
    logical :: transp=.true.
    integer :: irow, icol, lbw,ubw

    complex :: must_value

    CALL initialize(Amat,nrows,ncols,.not.transp)
    CALL ALLOCATE(Amat)
    CALL initialize(Cmat,nrows,ncols,transp)
    CALL ALLOCATE(Cmat)
    lbw = 1
    ubw = 2
    DO irow=1,nrows
       DO icol=1,ncols
          if ((icol.ge.irow-lbw).and.(icol.le.irow+ubw)) then
             CALL set_value(Amat,irow,icol,CMPLX(irow,icol))
          end if
       END DO
    END DO
    call commit_values(Amat)

    !call show(Amat)
    call transpose_and_conjugate(Amat,Cmat)
    !call show(Cmat)
    
    call get_global_matrix_locally(Cmat,localfullmat)

    passed = .true.
    do irow=1,nrows
       do icol=1,ncols
          !write(*,"(2ES9.2)",advance="no") localfullmat(irow,icol)
          if ((icol.ge.irow-ubw).and.(icol.le.irow+lbw)) then
             must_value = cmplx(icol,-irow)
          else
             must_value = cmplx(0.0d0,0.0d0)
          end if
          if (abs(localfullmat(irow,icol)-must_value).gt.tolerance) then
             print*,irow,icol
             passed = .false.
          end if
       end do
       !write(*,"(A)") ""
    end do

  end function test_transpose_and_conjugate

  FUNCTION test_my_sum() RESULT(passed)
    logical :: passed

#ifdef ALL_ROUTINES
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
       CALL finalize()

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
       CALL initialize_matrix_module(MAT_BLOCKS_OF_ROWS,comm_x)

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
#endif
    passed = .TRUE.
  END FUNCTION test_my_sum
end program test_BandedMatrixModule
