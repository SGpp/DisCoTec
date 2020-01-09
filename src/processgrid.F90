#include "intrinsic_sizes.h"

MODULE ProcessGridModule
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: MAT_BLOCKS_OF_ROWS, MAT_BLOCKS_OF_COLS, MAT_BOTH
  PUBLIC :: ProcessGrid, initialize, finalize, operator(.EQ.)
  PUBLIC :: size_of_processgrid

  INTEGER, PARAMETER :: MAT_BLOCKS_OF_ROWS=1000, MAT_BLOCKS_OF_COLS=1001, MAT_BOTH=1002

  TYPE ProcessGrid
     INTEGER :: Context        ! BLACS context for communication
     INTEGER :: Communicator   ! MPI communicator
     INTEGER :: NProcs         ! number of all processes of the context
     INTEGER :: NPRows, NPCols ! number of process block sizes
     INTEGER :: Rank           ! the rank of the process in the given communicator (not the global rank)
     INTEGER :: PRow           ! processor coordinates of the local processor in the grid
     INTEGER :: PCol
     LOGICAL :: isInitialized=.FALSE.
  END TYPE ProcessGrid

  INTERFACE initialize
     MODULE PROCEDURE mp_initialize_PG
  END INTERFACE

  INTERFACE finalize
     MODULE PROCEDURE mp_finalize_PG
  END INTERFACE

  INTERFACE OPERATOR(.EQ.)
     module procedure mp_isequal
  END INTERFACE

CONTAINS
  FUNCTION size_of_processgrid() RESULT(memory_need)
    real :: memory_need
    
    REAL :: mem_integer, mem_logical

    mem_integer = SIZE_OF_INTEGER/1024./1024.
    mem_logical = SIZE_OF_LOGICAL/1024./1024.

    memory_need = 8 * mem_integer + mem_logical
  END FUNCTION size_of_processgrid

  SUBROUTINE mp_initialize_PG(PG,p_dim,comm)
    TYPE(ProcessGrid),intent(INOUT) :: PG
    INTEGER,intent(IN) :: p_dim, comm
    
    ! Local variables
    INTEGER :: np_rows, np_cols, ierr!, blacs_debuglevel

    IF (.NOT.PG%isInitialized) THEN
       CALL mpi_comm_size(comm,PG%NProcs,ierr)
       CALL mpi_comm_rank(comm,PG%Rank,ierr)

       !IF (PG%Rank.EQ.0) THEN
       !   PRINT*,"Initializing ProcessGrid for BLACS and ScaLapack routines with ",&
       !        &PG%NProcs," processes."
       !END IF
       
       SELECT CASE (p_dim)
       CASE (MAT_BLOCKS_OF_ROWS)
          PG%NPRows = PG%NProcs
          PG%NPCols = 1
          
       CASE (MAT_BLOCKS_OF_COLS)
          PG%NPRows = 1
          PG%NPCols = PG%NProcs
          
       CASE (MAT_BOTH)
          PRINT*,"MAT_BOTH is not yet implemented."
       END SELECT
       
       PG%Communicator = comm
       PG%Context = comm
       !CALL blacs_get( PG%Context, 2, blacs_debuglevel )
       !PRINT*,"BLACS debuglevel is ",blacs_debuglevel
       CALL BLACS_GRIDINIT( PG%Context,'C',PG%NPRows,PG%NPCols)
       CALL BLACS_GRIDINFO( PG%Context, np_rows, np_cols, PG%PRow, PG%PCol )
       !PRINT*,PG%Rank, ':',np_rows, np_cols, PG%PRow, PG%PCol

       PG%isInitialized = .TRUE.
    END IF
  END SUBROUTINE mp_initialize_PG

  SUBROUTINE mp_finalize_PG(PG)
    type(ProcessGrid) :: PG

!    INTEGER :: np_rows,np_cols,prow,pcol

    IF (PG%isInitialized) THEN
       !PRINT*,"Freeing the blacs context ",PG%Context
       CALL BLACS_GRIDEXIT(PG%Context)
       !CALL BLACS_GRIDINFO( PG%Context, np_rows, np_cols, prow,pcol )
       !PRINT*,'After gridexit: ',np_rows, np_cols, prow,pcol
       PG%isInitialized=.FALSE.
    END IF
    
  END SUBROUTINE mp_finalize_PG

  LOGICAL FUNCTION mp_isequal(lPG,rPG)
    TYPE(ProcessGrid), intent(IN) :: lPG,rPG

    !PRINT*,"lPG%Context = ",lPG%Context,", rPG%Context = ", rPG%Context
    mp_isequal = (lPG%Context==rPG%Context)
  END FUNCTION mp_isequal

END MODULE ProcessGridModule
