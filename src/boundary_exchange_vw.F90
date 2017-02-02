#include "redef.h"
MODULE boundary_exchange_vw
  use BoundaryDescriptionModule
  use boundary_exchange_general
  implicit none

  PRIVATE
  PUBLIC :: exchange_v, exchange_mu, initialize_boundary_exchange_vw, finalize_boundary_exchange_vw
  
  INTERFACE exchange_v
     MODULE PROCEDURE exchange_v_4D,exchange_v_5D, exchange_v_6D
  END INTERFACE

  INTERFACE exchange_mu
     MODULE PROCEDURE exchange_mu_5D, exchange_mu_6D
  END INTERFACE

  INTERFACE initialize_boundary_exchange_vw
     MODULE PROCEDURE bevw_initialize_boundary_exchange_vw
  END INTERFACE

  INTEGER :: my_pev, my_pew, n_procs_v, n_procs_w
CONTAINS

  SUBROUTINE bevw_initialize_boundary_exchange_vw(ai_comm_v, ai_comm_w)
    INTEGER, INTENT(IN) :: ai_comm_v, ai_comm_w

    ! Local variables
    integer :: ierr

    CALL mpi_comm_rank(ai_comm_v,my_pev,ierr)
    CALL mpi_comm_rank(ai_comm_w,my_pew,ierr)
    CALL mpi_comm_size(ai_comm_v,n_procs_v, ierr)
    CALL mpi_comm_size(ai_comm_w,n_procs_w, ierr)

  END SUBROUTINE bevw_initialize_boundary_exchange_vw

  SUBROUTINE finalize_boundary_exchange_vw
  END SUBROUTINE finalize_boundary_exchange_vw

  SUBROUTINE exchange_v_4D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u)).OR.bdesc%count.gt.1) THEN
       PRINT*,"Wrong BoundaryDescription in call to exchange_v_4D!"
       PRINT*,"bdesc%n_points = ",bdesc%n_points,", bdesc%count = ", bdesc%count
       STOP
    END IF
    
    CALL exchange_v_general(bdesc, u,1)
  END SUBROUTINE exchange_v_4D

  SUBROUTINE exchange_v_5D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1)*SIZE(u,2)*SIZE(u,3)*SIZE(u,4))&
         &.OR.bdesc%count.GT.SIZE(u,5)) THEN
       STOP "Wrong BoundaryDescription in call to exchange_v_5D!"
    END IF
    
    CALL exchange_v_general(bdesc, u,SIZE(u,5))
  END SUBROUTINE exchange_v_5D

  SUBROUTINE exchange_v_6D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1)*SIZE(u,2)*SIZE(u,3)*SIZE(u,4))&
         &.OR.bdesc%count.GT.SIZE(u,5)*SIZE(u,6)) THEN
       PRINT*,"n_points = ", bdesc%n_points,", count = ", bdesc%count
       STOP "Wrong BoundaryDescription in call to exchange_v_6D!"
    END IF
    
    CALL exchange_v_general(bdesc, u,SIZE(u,5)*SIZE(u,6))
  END SUBROUTINE exchange_v_6D

  SUBROUTINE exchange_mu_5D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u)).OR.bdesc%count.GT.1) THEN
       STOP "Wrong BoundaryDescription in call to exchange_mu_5D!"
    END IF
    
    CALL exchange_mu_general(bdesc, u,1)
  END SUBROUTINE exchange_mu_5D

  SUBROUTINE exchange_mu_6D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1)*SIZE(u,2)*SIZE(u,3)*SIZE(u,4)*SIZE(u,5))&
         &.OR.bdesc%count.GT.SIZE(u,6)) THEN
       STOP "Wrong BoundaryDescription in call to exchange_mu_6D!"
    END IF
    
    CALL exchange_mu_general(bdesc, u,SIZE(u,6))
  END SUBROUTINE exchange_mu_6D

!================= general routines =============

  SUBROUTINE exchange_v_general(bdesc,u, n_dim_2)
    TYPE(BoundaryDescription) :: bdesc
    INTEGER :: n_dim_2
    COMPLEX, DIMENSION(bdesc%n_points,1:n_dim_2) :: u

    ! Local variables
    INTEGER :: n_exchanges, i_exchange
    PERFON('ex_vp')

    IF (MOD(n_dim_2,bdesc%count).NE.0) THEN
       WRITE(*,"(2(A,I6),A)") "Wrong BoundaryDescription in call to exchange_v_general! We cannot exchange ",&
            & n_dim_2," times in v direction in blocks of ",bdesc%count,". Aborting!"
       STOP
    END IF

    ! Here we apply the physical boundary conditions as pre-process of the exchange
    ! Dirichlet boundary conditions in v direction
    !IF (my_pev.EQ.0) THEN
    !   u(1:bdesc%lower,:) = CMPLX(0,0,KIND(u))
    !END IF
    !IF (my_pev.EQ.n_procs_v-1) THEN
    !   u(bdesc%n_points-bdesc%upper+1:bdesc%n_points,:) = CMPLX(0,0,KIND(u))
    !END IF
    ! End of application of physical boundary conditions

    !WRITE(*,"(A,ES20.10)") "-deb- in exchange_vw_general, u = ", DBLE(SUM(u*CONJG(u)))
    !PRINT*, "size = ", SIZE(u,1),SIZE(u,2)
    bdesc%exchange_direction = 2 ! set to v direction
    !bdesc%communicator = my_comm_vw
    n_exchanges = n_dim_2/bdesc%count
    DO i_exchange=0,n_exchanges-1
       !PRINT*,i_exchange,1+i_exchange*bdesc%count,(i_exchange+1)*bdesc%count
       CALL exchange_general(bdesc,u(:,1+i_exchange*bdesc%count:(i_exchange+1)*bdesc%count))
    END DO
    !WRITE(*,"(A,ES20.10)") "-deb- after periodic exchange, u = ", DBLE(SUM(u*CONJG(u)))

    ! Here we apply the physical boundary conditions as post-process of the exchange.
    ! End of application of physical boundary conditions
    PERFOFF
  END SUBROUTINE exchange_v_general

  SUBROUTINE exchange_mu_general(bdesc,u, n_dim_2)
    TYPE(BoundaryDescription) :: bdesc
    INTEGER :: n_dim_2
    COMPLEX, DIMENSION(bdesc%n_points,1:n_dim_2) :: u

    ! Local variables
    INTEGER :: n_exchanges, i_exchange

    PERFON('ex_mu')
    IF (MOD(n_dim_2,bdesc%count).NE.0) THEN
       WRITE(*,"(2(A,I6),A)") "Wrong BoundaryDescription in call to exchange_v_general! We cannot exchange ",&
            & n_dim_2," times in mu direction in blocks of ",bdesc%count,". Aborting!"
       STOP
    END IF

    ! Here we apply the physical boundary conditions as pre-process of the exchange
    !IF (my_pew.EQ.0) THEN
    !   u(1:bdesc%lower,:) = CMPLX(0,0,KIND(u))
    !END IF
    !IF (my_pew.EQ.n_procs_w-1) THEN
    !   u(bdesc%n_points-bdesc%upper+1:bdesc%n_points,:) = CMPLX(0,0,KIND(u))
    !END IF
    ! End of application of physical boundary conditions

    !WRITE(*,"(A,ES20.10)") "-deb- in exchange_vw_general, u = ", DBLE(SUM(u*CONJG(u)))
    !PRINT*, "size = ", SIZE(u,1),SIZE(u,2)
    bdesc%exchange_direction = 1 ! set to mu direction
    !bdesc%communicator = my_comm_vw
    n_exchanges = n_dim_2/bdesc%count
    DO i_exchange=0,n_exchanges-1
       !PRINT*,i_exchange,1+i_exchange*bdesc%count,(i_exchange+1)*bdesc%count
       CALL exchange_general(bdesc,u(:,1+i_exchange*bdesc%count:(i_exchange+1)*bdesc%count))
    END DO
    !WRITE(*,"(A,ES20.10)") "-deb- after periodic exchange, u = ", DBLE(SUM(u*CONJG(u)))

    ! Here we apply the physical boundary conditions as post-process of the exchange.
    ! End of application of physical boundary conditions
    PERFOFF
  END SUBROUTINE exchange_mu_general

END MODULE boundary_exchange_vw
