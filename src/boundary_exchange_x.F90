#include "redef.h"
MODULE boundary_exchange_x
  use discretization, only : lj0 ! this should be changed, to be independent of discretization
  use BoundaryDescriptionModule
  use boundary_exchange_general
  use par_other, only: p_has_0_mode
  IMPLICIT NONE

  PUBLIC :: exchange_x, initialize_boundary_exchange_x, finalize_boundary_exchange_x,&
       & exchange_x_is_initialized

  PRIVATE

  INTERFACE exchange_x
     MODULE PROCEDURE exchange_x_1D,exchange_x_2D,exchange_x_2D_real,exchange_x_3D,&
          & exchange_x_4D,exchange_x_5D, exchange_x_6D
  END INTERFACE

  INTEGER :: rad_bc_type, my_pex, n_procs_x
  LOGICAL :: is_initialized=.FALSE.
  
CONTAINS

  SUBROUTINE initialize_boundary_exchange_x(ai_comm_x,ai_rad_bc_type)
    INTEGER, INTENT(IN) :: ai_comm_x, ai_rad_bc_type

    ! Local variables
    integer :: ierr

    rad_bc_type = ai_rad_bc_type

    CALL mpi_comm_rank(ai_comm_x, my_pex, ierr)
    CALL mpi_comm_size(ai_comm_x, n_procs_x, ierr)

    is_initialized = .TRUE.
  END SUBROUTINE initialize_boundary_exchange_x

  FUNCTION exchange_x_is_initialized()
    logical :: exchange_x_is_initialized
    exchange_x_is_initialized = is_initialized
  END FUNCTION exchange_x_is_initialized

  FUNCTION exchange_x_is_parallelized()
    logical :: exchange_x_is_parallelized
    exchange_x_is_parallelized = (n_procs_x.gt.1)
  END FUNCTION exchange_x_is_parallelized

  SUBROUTINE finalize_boundary_exchange_x
    is_initialized = .FALSE.
  END SUBROUTINE finalize_boundary_exchange_x
    
  SUBROUTINE exchange_x_1D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, dimension(:) :: u

    IF ( (bdesc%n_points.NE.SIZE(u,1))&
         &.OR.bdesc%count.NE.1 ) THEN
       STOP "Wrong BoundaryDescription in call to exchange_x_1D!"
    END IF
    
    CALL exchange_x_general(bdesc,u,1)
  END SUBROUTINE exchange_x_1D

  SUBROUTINE exchange_x_2D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1))&
         &.OR.bdesc%count.NE.SIZE(u,2)) THEN
       PRINT*,"Wrong BoundaryDescription in call to exchange_x_2D!"
       PRINT*,"bdesc%n_points,size(u,1) = ",bdesc%n_points,SIZE(u,1),&
            &", bdesc%count,size(u,2) = ", bdesc%count,SIZE(u,2)
       PRINT*,"bdesc = ",bdesc
       STOP 
    END IF

    CALL exchange_x_general(bdesc, u,SIZE(u,2))
  END SUBROUTINE exchange_x_2D

  SUBROUTINE exchange_x_2D_real(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    REAL, DIMENSION(:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1))&
         &.OR.bdesc%count.NE.SIZE(u,2)) THEN
       PRINT*,"Wrong BoundaryDescription in call to exchange_x_2D!"
       PRINT*,"bdesc%n_points,size(u,1) = ",bdesc%n_points,SIZE(u,1),&
            &", bdesc%count,size(u,2) = ", bdesc%count,SIZE(u,2)
       PRINT*,"bdesc = ",bdesc
       STOP 
    END IF

    CALL exchange_x_general_real(bdesc, u,SIZE(u,2))
  END SUBROUTINE exchange_x_2D_real

  SUBROUTINE exchange_x_3D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1))&
         &.OR.bdesc%count.gt.SIZE(u,2)*SIZE(u,3)) THEN
       STOP "Wrong BoundaryDescription in call to exchange_x_3D!"
    END IF
    
    CALL exchange_x_general(bdesc, u,SIZE(u,2)*SIZE(u,3))
  END SUBROUTINE exchange_x_3D

  SUBROUTINE exchange_x_4D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1))&
         &.OR.bdesc%count.GT.SIZE(u,2)*SIZE(u,3)*SIZE(u,4)) THEN
       STOP "Wrong BoundaryDescription in call to exchange_x_4D!"
    END IF
    
    CALL exchange_x_general(bdesc, u,SIZE(u,2)*SIZE(u,3)*SIZE(u,4))
  END SUBROUTINE exchange_x_4D

  SUBROUTINE exchange_x_5D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1))&
         &.OR.bdesc%count.GT.SIZE(u,2)*SIZE(u,3)*SIZE(u,4)*SIZE(u,5)) THEN
       STOP "Wrong BoundaryDescription in call to exchange_x_5D!"
    END IF
    
    CALL exchange_x_general(bdesc, u,SIZE(u,2)*SIZE(u,3)*SIZE(u,4)*SIZE(u,5))
  END SUBROUTINE exchange_x_5D

  SUBROUTINE exchange_x_6D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1))&
         &.OR.bdesc%count.GT.SIZE(u,2)*SIZE(u,3)*SIZE(u,4)*SIZE(u,5)*SIZE(u,6)) THEN
       STOP "Wrong BoundaryDescription in call to exchange_x_6D!"
    END IF
    
    CALL exchange_x_general(bdesc, u, SIZE(u,2)*SIZE(u,3)*SIZE(u,4)*SIZE(u,5)*SIZE(u,6))
  END SUBROUTINE exchange_x_6D

  SUBROUTINE exchange_x_general(bdesc,u, n_dim_2)
    TYPE(BoundaryDescription) :: bdesc
    integer :: n_dim_2
    COMPLEX, DIMENSION(0:bdesc%n_points-1,1:n_dim_2) :: u

    ! Local variables
    INTEGER :: n_exchanges, i_exchange, i_count,i

    PERFON('ex_x')
    IF (MOD(n_dim_2,bdesc%count).NE.0) THEN
       WRITE(*,"(2(A,I6),A)") "Wrong BoundaryDescription in call to exchange_x_general! We cannot exchange ",n_dim_2,&
         " times in x direction in blocks of ",bdesc%count,". Aborting!"
       STOP
    END IF

    ! Here we apply the physical boundary conditions as pre-process of the exchange
    IF ((rad_bc_type.EQ.1).OR.(rad_bc_type.EQ.2)) THEN
       ! Dirichlet boundary condition
       ! these statements are cause of a slight load imbalance, as only the first and last point
       ! have to do the setting of the boundary points
       IF (my_pex.EQ.0) THEN
          u(0:bdesc%lower-1,:) = CMPLX(0,0,KIND(u))
       END IF
       IF (my_pex.EQ.n_procs_x-1) THEN
          u(bdesc%n_points-bdesc%upper:bdesc%n_points-1,:) = CMPLX(0,0,KIND(u))
       END IF
    END IF

    bdesc%exchange_direction = 5 ! set to x direction

    n_exchanges = n_dim_2/bdesc%count
    DO i_exchange=0,n_exchanges-1
       CALL exchange_general(bdesc,u(:,1+i_exchange*bdesc%count:(i_exchange+1)*bdesc%count))
    END DO

    ! Here we apply the physical boundary conditions as post-process of the exchange
    If ((rad_bc_type.eq.2).and.(p_has_0_mode).and.(my_pex.eq.0)) then
       ! We need the symmetric boundary for all outer indizes
       ! This solution is not nice, as we assume that lj0 is always
       ! the number of y points.
       if (modulo(bdesc%count,lj0).ne.0) then
          print*,"bdesc%count = ",bdesc%count,", lj0 = ",lj0
          stop
       end if
       do i_count=1,bdesc%count/lj0
          Do i=1,bdesc%lower
             u(bdesc%innerfirst-i,1+(i_count-1)*lj0) = u(bdesc%innerfirst+i,1+(i_count-1)*lj0)
          Enddo
       end do
    Endif

    PERFOFF
  END SUBROUTINE exchange_x_general

  SUBROUTINE exchange_x_general_real(bdesc,u, n_dim_2)
    TYPE(BoundaryDescription) :: bdesc
    integer :: n_dim_2
    REAL, DIMENSION(0:bdesc%n_points-1,1:n_dim_2) :: u

    ! Local variables
    INTEGER :: n_exchanges, i_exchange, i, i_count

    !PERFON('ex_x_re')
    IF (MOD(n_dim_2,bdesc%count).NE.0) THEN
       WRITE(*,"(2(A,I6),A)") "Wrong BoundaryDescription in call to exchange_x_general! We cannot exchange ",n_dim_2,&
         " times in x direction in blocks of ",bdesc%count,". Aborting!"
       STOP
    END IF

    ! Here we apply the physical boundary conditions as pre-process of the exchange
    IF ((rad_bc_type.EQ.1).OR.(rad_bc_type.EQ.2)) THEN
       ! Dirichlet boundary condition
       ! these statements are cause of a slight load imbalance, as only the first and last point
       ! have to do the setting of the boundary points
       IF (my_pex.EQ.0) THEN
          u(0:bdesc%lower-1,:) = 0.
       END IF
       IF (my_pex.EQ.n_procs_x-1) THEN
          u(bdesc%n_points-bdesc%upper:bdesc%n_points-1,:) = 0.
       END IF
    END IF

    bdesc%exchange_direction = 5 ! set to x direction

    n_exchanges = n_dim_2/bdesc%count
    DO i_exchange=0,n_exchanges-1
       CALL exchange_general_real(bdesc,u(:,1+i_exchange*bdesc%count:(i_exchange+1)*bdesc%count))
    END DO

    ! Here we apply the physical boundary conditions as post-process of the exchange
    If ((rad_bc_type.eq.2).and.(p_has_0_mode).and.(my_pex.eq.0)) then
       ! We need the symmetric boundary for all outer indizes
       ! This solution is not nice, as we assume that lj0 is always
       ! the number of y points.
       do i_count=1,bdesc%count/lj0
          Do i=1,bdesc%lower
             u(bdesc%innerfirst-i,1+(i_count-1)*lj0) = u(bdesc%innerfirst+i,1+(i_count-1)*lj0)
          Enddo
       end do
    Endif

    !PERFOFF
  END SUBROUTINE exchange_x_general_real


END MODULE boundary_exchange_x
