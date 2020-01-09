#include "redef.h"
MODULE boundary_exchange_general
  use mpi
  USE BoundaryDescriptionModule
  IMPLICIT NONE

  PUBLIC :: exchange_general, initialize_boundary_exchange_general, finalize_boundary_exchange_general
  PUBLIC :: exchange_general_real
  PRIVATE

  INTEGER :: comm_cart
  INTEGER :: cartdim
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: periods
  INTEGER, ALLOCATABLE, DIMENSION(:) :: dims, coords

CONTAINS
  SUBROUTINE initialize_boundary_exchange_general(ai_comm_cart)
    INTEGER, intent(IN) :: ai_comm_cart
    
    ! Local variables
    integer :: ierr

    comm_cart = ai_comm_cart
    CALL mpi_cartdim_get(comm_cart,cartdim, ierr)
    ALLOCATE(periods(0:cartdim-1),dims(0:cartdim-1),coords(0:cartdim-1))
    CALL mpi_cart_get(comm_cart, cartdim, dims, periods, coords, ierr)
  END SUBROUTINE initialize_boundary_exchange_general

  SUBROUTINE finalize_boundary_exchange_general
    DEALLOCATE(periods,dims,coords)
  END SUBROUTINE finalize_boundary_exchange_general

  ! --------------------------------------------
  !> General boundary exchange. Is called for all boundary exchanges.
  !!
  !! This is the main boundary exchange routine. It takes a boundary or
  !! exchange description and a 2D array. It does a boundary exchange along
  !! one direction, which is specified in the BoundaryDescription. No physical
  !! boundary conditions are implemented in this routine, as the physical
  !! boundary conditions vary over the directions. The routine is only called
  !! by the direction specific but general routines.
  !! The routine just copies the content of the inner cells to the boundary (ghost)
  !! cells, either over MPI network or directly (if no parallelization is done
  !! in the specified direction).
  !!
  SUBROUTINE exchange_general(bdesc,u)
    TYPE(BoundaryDescription),INTENT(IN) :: bdesc
    COMPLEX, DIMENSION(0:bdesc%n_points-1,0:bdesc%count-1),INTENT(INOUT) :: u

    ! Local variables
    INTEGER :: lboundary, uboundary, innerfirst, innerlast, n_innerpoints
    INTEGER :: stat(MPI_STATUS_SIZE), ierr, tag, rank_source, rank_dest, n_procs
    INTEGER :: my_proc, my_thread
#ifdef WITHOMP
    integer :: omp_get_thread_num
#endif

    PERFON_I('ex_gen')
#ifdef WITHOMP
    my_thread = omp_get_thread_num()
#else
    my_thread = 0
#endif

    lboundary  = bdesc%lower
    uboundary  = bdesc%upper
    innerfirst = bdesc%innerfirst
    innerlast  = bdesc%innerlast
    n_innerpoints = innerlast-innerfirst+1

    n_procs = dims(bdesc%exchange_direction)
    ! is the direction parallelized over several MPI processes?
    IF (n_procs.GT.1) THEN
       ! check further if there is a boundary at all
       IF ( (lboundary.EQ.0).AND.(uboundary.EQ.0) ) RETURN
       
       my_proc = coords(bdesc%exchange_direction)
       !WRITE(*,"(I3,A,ES20.12)") my_proc,": u(inner) = ", DBLE(SUM(u(innerfirst:innerlast,:)*CONJG(u(innerfirst:innerlast,:))))

       ! there is a boundary, so we do an exchange with MPI

       IF ((n_innerpoints.GE.lboundary).AND.(n_innerpoints.GE.uboundary)) THEN
          !
          ! Send boundary to next pex and receive from previous:
          !
          CALL mpi_cart_shift(comm_cart,bdesc%exchange_direction,1,rank_source,rank_dest,ierr)
          !WRITE(*,"(I3,3(A,I3))") this_rank,": rank_source = ",rank_source,", rank_dest = ", rank_dest,&
          !     &", exdir = ", bdesc%exchange_direction
          !PERFON_I('exgensr')
          tag = 351+my_thread !__LINE__
          CALL MPI_Sendrecv(u(innerlast-lboundary+1,0), 1, &
                  & bdesc%mpi_lower_datatype, rank_dest, tag,&
                  u(0,0), 1, bdesc%mpi_lower_datatype, rank_source, tag,&
                  comm_cart, stat, ierr)
          !PERFOFF_I
          !
          ! Send boundary to previous pex and receive from next:
          !
          CALL mpi_cart_shift(comm_cart,bdesc%exchange_direction,-1,rank_source,rank_dest,ierr)

          !PERFON_I('exgensl')
          tag = 370+my_thread !__LINE__
          CALL MPI_Sendrecv(u(innerfirst,0), 1, &
               &bdesc%mpi_upper_datatype, rank_dest, tag,&
               u(innerlast+1,0), 1, bdesc%mpi_upper_datatype, rank_source, tag,&
               comm_cart, stat, ierr)
          !PERFOFF_I
       ELSE
          ! now we have more boundary points than inner points. To overcome this
          ! problem we do an ordered exchange, and start from the physical boundaries
          my_proc = coords(bdesc%exchange_direction)
          !np_involved = lboundary/n_innerpoints
          !IF (MODULO(lboundary,n_innerpoints).NE.0) np_involved = np_involved+1

          CALL mpi_cart_shift(comm_cart,bdesc%exchange_direction,1,rank_source,rank_dest,ierr)
          tag = 480+my_thread
          IF (my_proc.LT.n_procs-1) THEN
             !PRINT*,my_proc,": lo sending from ",my_proc," to ", rank_dest
             CALL mpi_send(u(innerlast-lboundary+1,0),1,&
                  &bdesc%mpi_lower_datatype,rank_dest,tag,comm_cart,ierr)
             !PRINT*,my_proc,": lo end sending"
          END IF
          IF (my_proc.GT.0) THEN
             !PRINT*,my_proc,": lo receiving on ",my_proc," from ", rank_source
             CALL mpi_recv(u(0,0),1,bdesc%mpi_lower_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
             !PRINT*,my_proc,": lo end receiving"
          END IF

          IF (my_proc.EQ.n_procs-1) THEN
             CALL mpi_send(u(innerlast-lboundary+1,0),1,&
                  &bdesc%mpi_lower_datatype,rank_dest,tag,comm_cart,ierr)
             CALL mpi_recv(u(0,0),1,bdesc%mpi_lower_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
          END IF
          IF (my_proc.LT.n_procs-1) THEN
             CALL mpi_recv(u(0,0),1,bdesc%mpi_lower_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
             CALL mpi_send(u(innerlast-lboundary+1,0),1,&
                  &bdesc%mpi_lower_datatype,rank_dest,tag,comm_cart,ierr)
          END IF

          !PRINT*,mype,"upper boundary exchange"
          ! other direction
          CALL mpi_cart_shift(comm_cart,bdesc%exchange_direction,-1,rank_source,rank_dest,ierr)
          tag = 500+my_thread
          IF (my_proc.GT.0) THEN
             CALL mpi_send(u(innerfirst,0),1,&
                  &bdesc%mpi_upper_datatype,rank_dest,tag,comm_cart,ierr)
          END IF
          IF (my_proc.LT.n_procs-1) THEN
             CALL mpi_recv(u(innerlast+1,0),1,bdesc%mpi_upper_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
          END IF
          
          IF (my_proc.EQ.0) THEN
             CALL mpi_send(u(innerfirst,0),1,&
                  &bdesc%mpi_upper_datatype,rank_dest,tag,comm_cart,ierr)
             CALL mpi_recv(u(innerlast+1,0),1,bdesc%mpi_upper_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
          ELSE
             CALL mpi_recv(u(innerlast+1,0),1,bdesc%mpi_upper_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
             CALL mpi_send(u(innerfirst,0),1,&
                  &bdesc%mpi_upper_datatype,rank_dest,tag,comm_cart,ierr)
          END IF
          
          
          !PRINT*,"more boundary points than innerpoints in exchange_direction ",&
          !     &bdesc%exchange_direction," (0:spec,1:w,2:v,3:z,4:y,5:x)"
          !STOP
       END IF
    ELSE
       !n_procs .eq. 1, just one process in the exchange direction
       ! if the direction is periodic, we exchange the boundaries, if not, we do nothing.
       !CALL mpi_cart_get(comm_cart,6,dims,periods,coords,ierr)
       IF (periods(bdesc%exchange_direction)) THEN
          ! topologically periodic
          u(0:bdesc%lower-1,:) = u(innerlast-bdesc%lower+1:innerlast,:)
          u(innerlast+1:innerlast+bdesc%upper,:)=u(innerfirst:innerfirst+bdesc%upper-1,:)
       ELSE
          ! topologically non-periodic, do nothing
       END IF
    END IF
    PERFOFF_I
  END SUBROUTINE exchange_general

  !> General boundary exchange. Is called for all boundary exchanges.
  !!
  !! This is the main boundary exchange routine. It takes a boundary or
  !! exchange description and a 2D array. It does a boundary exchange along
  !! one direction, which is specified in the BoundaryDescription. No physical
  !! boundary conditions are implemented in this routine, as the physical
  !! boundary conditions vary over the directions. The routine is only called
  !! by the direction specific but general routines.
  !! The routine just copies the content of the inner cells to the boundary (ghost)
  !! cells, either over MPI network or directly (if no parallelization is done
  !! in the specified direction).
  !!
  SUBROUTINE exchange_general_real(bdesc,u)
    TYPE(BoundaryDescription),INTENT(IN) :: bdesc
    REAL, DIMENSION(0:bdesc%n_points-1,0:bdesc%count-1),INTENT(INOUT) :: u

    ! Local variables
    INTEGER :: lboundary, uboundary, innerfirst, innerlast, n_innerpoints
    INTEGER :: stat(MPI_STATUS_SIZE), ierr, tag, rank_source, rank_dest, n_procs
    INTEGER :: my_proc, my_thread
#ifdef WITHOMP
    integer :: omp_get_thread_num
#endif

    PERFON('ex_gen')
#ifdef WITHOMP
    my_thread = omp_get_thread_num()
#else
    my_thread = 0
#endif

    lboundary  = bdesc%lower
    uboundary  = bdesc%upper
    innerfirst = bdesc%innerfirst
    innerlast  = bdesc%innerlast
    n_innerpoints = innerlast-innerfirst+1

    n_procs = dims(bdesc%exchange_direction)
    ! is the direction parallelized over several MPI processes?
    IF (n_procs.GT.1) THEN
       ! check further if there is a boundary at all
       IF ( (lboundary.EQ.0).AND.(uboundary.EQ.0) ) RETURN
       
       my_proc = coords(bdesc%exchange_direction)
       !WRITE(*,"(I3,A,ES20.12)") my_proc,": u(inner) = ", DBLE(SUM(u(innerfirst:innerlast,:)*CONJG(u(innerfirst:innerlast,:))))

       ! there is a boundary, so we do an exchange with MPI

       IF ((n_innerpoints.GE.lboundary).AND.(n_innerpoints.GE.uboundary)) THEN
          !
          ! Send boundary to next pex and receive from previous:
          !
          CALL mpi_cart_shift(comm_cart,bdesc%exchange_direction,1,rank_source,rank_dest,ierr)
          !WRITE(*,"(I3,3(A,I3))") this_rank,": rank_source = ",rank_source,", rank_dest = ", rank_dest,&
          !     &", exdir = ", bdesc%exchange_direction
          tag = 351+my_thread !__LINE__
          CALL MPI_Sendrecv(u(innerlast-lboundary+1,0), 1, &
                  & bdesc%mpi_lower_datatype, rank_dest, tag,&
                  u(0,0), 1, bdesc%mpi_lower_datatype, rank_source, tag,&
                  comm_cart, stat, ierr)

          !
          ! Send boundary to previous pex and receive from next:
          !
          CALL mpi_cart_shift(comm_cart,bdesc%exchange_direction,-1,rank_source,rank_dest,ierr)

          tag = 370+my_thread !__LINE__
          CALL MPI_Sendrecv(u(innerfirst,0), 1, &
               &bdesc%mpi_upper_datatype, rank_dest, tag,&
               u(innerlast+1,0), 1, bdesc%mpi_upper_datatype, rank_source, tag,&
               comm_cart, stat, ierr)
       ELSE
          ! now we have more boundary points than inner points. To overcome this
          ! problem we do an ordered exchange, and start from the physical boundaries
          my_proc = coords(bdesc%exchange_direction)
          !np_involved = lboundary/n_innerpoints
          !IF (MODULO(lboundary,n_innerpoints).NE.0) np_involved = np_involved+1

          CALL mpi_cart_shift(comm_cart,bdesc%exchange_direction,1,rank_source,rank_dest,ierr)
          tag = 480+my_thread
          IF (my_proc.LT.n_procs-1) THEN
             !PRINT*,my_proc,": lo sending from ",my_proc," to ", rank_dest
             CALL mpi_send(u(innerlast-lboundary+1,0),1,&
                  &bdesc%mpi_lower_datatype,rank_dest,tag,comm_cart,ierr)
             !PRINT*,my_proc,": lo end sending"
          END IF
          IF (my_proc.GT.0) THEN
             !PRINT*,my_proc,": lo receiving on ",my_proc," from ", rank_source
             CALL mpi_recv(u(0,0),1,bdesc%mpi_lower_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
             !PRINT*,my_proc,": lo end receiving"
          END IF

          IF (my_proc.EQ.n_procs-1) THEN
             CALL mpi_send(u(innerlast-lboundary+1,0),1,&
                  &bdesc%mpi_lower_datatype,rank_dest,tag,comm_cart,ierr)
             CALL mpi_recv(u(0,0),1,bdesc%mpi_lower_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
          END IF
          IF (my_proc.LT.n_procs-1) THEN
             CALL mpi_recv(u(0,0),1,bdesc%mpi_lower_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
             CALL mpi_send(u(innerlast-lboundary+1,0),1,&
                  &bdesc%mpi_lower_datatype,rank_dest,tag,comm_cart,ierr)
          END IF

          !PRINT*,mype,"upper boundary exchange"
          ! other direction
          CALL mpi_cart_shift(comm_cart,bdesc%exchange_direction,-1,rank_source,rank_dest,ierr)
          tag = 500+my_thread
          IF (my_proc.GT.0) THEN
             CALL mpi_send(u(innerfirst,0),1,&
                  &bdesc%mpi_upper_datatype,rank_dest,tag,comm_cart,ierr)
          END IF
          IF (my_proc.LT.n_procs-1) THEN
             CALL mpi_recv(u(innerlast+1,0),1,bdesc%mpi_upper_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
          END IF
          
          IF (my_proc.EQ.0) THEN
             CALL mpi_send(u(innerfirst,0),1,&
                  &bdesc%mpi_upper_datatype,rank_dest,tag,comm_cart,ierr)
             CALL mpi_recv(u(innerlast+1,0),1,bdesc%mpi_upper_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
          ELSE
             CALL mpi_recv(u(innerlast+1,0),1,bdesc%mpi_upper_datatype,&
                  &rank_source,tag,comm_cart,stat,ierr)
             CALL mpi_send(u(innerfirst,0),1,&
                  &bdesc%mpi_upper_datatype,rank_dest,tag,comm_cart,ierr)
          END IF
          
          
          !PRINT*,"more boundary points than innerpoints in exchange_direction ",&
          !     &bdesc%exchange_direction," (0:spec,1:w,2:v,3:z,4:y,5:x)"
          !STOP
       END IF
    ELSE
       !n_procs .eq. 1, just one process in the exchange direction
       ! if the direction is periodic, we exchange the boundaries, if not, we do nothing.
       !CALL mpi_cart_get(comm_cart,6,dims,periods,coords,ierr)
       IF (periods(bdesc%exchange_direction)) THEN
          ! topologically periodic
          u(0:bdesc%lower-1,:) = u(innerlast-bdesc%lower+1:innerlast,:)
          u(innerlast+1:innerlast+bdesc%upper,:)=u(innerfirst:innerfirst+bdesc%upper-1,:)
       ELSE
          ! topologically non-periodic, do nothing
       END IF
    END IF
    PERFOFF
  END SUBROUTINE exchange_general_real


END MODULE boundary_exchange_general
