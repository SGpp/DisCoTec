#include "redef.h"
!>Module for efficient computation of Y = A*X + Y operations
module axpy
  use communications
  use discretization, only: mype
  
  IMPLICIT NONE

  logical :: axpy_ij_use_blas = .true.

contains
  
  subroutine initialize_axpy_ij(N,ntests)
    INTEGER :: N, ntests
    INTEGER :: t
    double precision :: time1, time2
    real :: blas_time, loop_time
    COMPLEX, DIMENSION(1:N) :: X, Y, RES
    REAL :: A

    X = (1.0,0.5)
    A = 1.0

    axpy_ij_use_blas = .false.
    Y = (0.1,0.2)
    call my_barrier()
    Call get_systime(time1)    
    !use 1000*ntests to get significant times
    DO t=1, 1000*ntests
       call axpy_ij(N,A,X,Y)
    ENDDO
    Call get_systime(time2) 
    loop_time = time2-time1
    call sync_time(loop_time)
    RES = Y

    axpy_ij_use_blas = .true.
    !need to initialize saxpy with first call?
    call axpy_ij(N,A,X,Y)

    Y = (0.1,0.2)
    call my_barrier()
    Call get_systime(time1)    
    !use 1000*ntests to get significant times
    DO t=1, 1000*ntests
       call axpy_ij(N,A,X,Y)
    ENDDO
    Call get_systime(time2) 
    blas_time = time2-time1
    call sync_time(blas_time)
    
    !if (mype.eq.0) write(*,'(A,2I5,A,2ES20.10)') 'lij0/lb0 ', &
    !     &N, ntests, ' blas/loop time: ', blas_time, loop_time

    IF (abs(sum(Y)-sum(RES)).gt.epsilon(A)) &
         &stop 'different test results in initialize_axpy_ij'
    
    axpy_ij_use_blas = (blas_time.le.loop_time)

  end subroutine initialize_axpy_ij

  subroutine axpy_ij(N, A, X, Y)
    INTEGER, INTENT(IN) :: N
    REAL,INTENT(IN) :: A
    COMPLEX,INTENT(IN)    :: X(1:N)
    COMPLEX,INTENT(INOUT) :: Y(1:N)

    IF (axpy_ij_use_blas) THEN
       !use BLAS routine
       call saxpy(2*N, A, X(1), 1, Y(1), 1)
    ELSE
       !use a simple loop instead of BLAS
       Y(:) = Y(:) + A*X(:)
    ENDIF

  end subroutine axpy_ij

end module axpy
