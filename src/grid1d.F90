#include "redef.h"
#include "intrinsic_sizes.h"
MODULE Grid1DModule
  IMPLICIT NONE

  PUBLIC :: Grid1D, initialize, finalize, set_boundaries, set_fourier_boundaries
  PUBLIC :: get_gridspacing, get_number_of_nodes, is_periodic_boundary, get_node, print_grid
  PUBLIC :: size_of_grid1d, set_startindex, get_boundaries, get_startindex
  PRIVATE

  !> This type describes a grid in one dimension. 
  TYPE Grid1D
     integer :: n_nodes !<number of nodes 
     real :: left_boundary !<position of the left boundary
     real :: right_boundary !<position of the right boundary
     logical :: isEquidistant !< set to true if the grid is equidistant
     !> Set to true if the grid is periodic. In the periodic case the rightmost
     !> point is does not lie on the boundary if the grid is equidistant.
     logical :: isPeriodic 
     LOGICAL :: isFourierGrid !< true, if the nodes are k's (Fourier modes)
     integer :: startindex = 1 !< gives the index of the first point (necessary for get_position)
     REAL, DIMENSION(:), ALLOCATABLE :: nodes !<the positions of the grid nodes, starting with index 1
  END TYPE Grid1D

  INTERFACE initialize
     module procedure mp_initialize
  END INTERFACE
  
  INTERFACE finalize
     module procedure mp_finalize
  END INTERFACE

  INTERFACE print_grid
     module procedure grid1d_print_grid
  END INTERFACE

  INTERFACE set_boundaries
     MODULE PROCEDURE mp_set_boundaries
  END INTERFACE

  INTERFACE get_boundaries
     MODULE PROCEDURE grid1d_get_boundaries
  END INTERFACE

  INTERFACE set_fourier_boundaries
     MODULE PROCEDURE mp_set_fourier_boundaries
  END INTERFACE

  INTERFACE get_gridspacing
     MODULE PROCEDURE mp_get_gridspacing
  END INTERFACE

  INTERFACE get_number_of_nodes
     MODULE PROCEDURE mp_get_number_of_nodes
  END INTERFACE

  INTERFACE get_startindex
     module procedure grid1d_get_startindex
  END INTERFACE

  INTERFACE is_periodic_boundary
     module procedure mp_is_periodic_boundary
  END INTERFACE

  INTERFACE get_node
     module procedure grid1d_get_node
  END INTERFACE

  INTERFACE set_startindex
     module procedure grid1d_set_startindex
  END INTERFACE

CONTAINS
  FUNCTION size_of_grid1d(n_points) RESULT(memory_need)
    INTEGER,intent(IN) :: n_points
    INTEGER :: memory_need

    memory_need = 2*SIZE_OF_INTEGER&
         & + 2*SIZE_OF_REAL&
         & + 3*SIZE_OF_LOGICAL&
         & + n_points*SIZE_OF_REAL
  END FUNCTION size_of_grid1d

  !>Initializes the Grid1D datatype. For initialization the number of nodes is necessary. All other
  !>elements are set to default values. Allocates the array for the grid nodes.
  SUBROUTINE mp_initialize(grid,p_n_nodes,p_isFourierGrid)
    TYPE(Grid1D) :: grid !< An instance of a Grid1D.
    INTEGER,intent(IN) :: p_n_nodes !< The number of nodes along the dimension of the grid.
    LOGICAL, INTENT(IN), OPTIONAL :: p_isFourierGrid

    IF (PRESENT(p_isFourierGrid)) THEN
       grid%isFourierGrid = p_isFourierGrid
    ELSE
       ! the default is a grid in real space
       grid%isFourierGrid = .FALSE.
    END IF

    IF (grid%isFourierGrid) THEN
       IF (p_n_nodes.GT.0) THEN
          grid%n_nodes = p_n_nodes
       ELSE
          PRINT*, "Number of Fourier modes must be at least 1!"
          STOP
       END IF
       ALLOCATE(grid%nodes(1:grid%n_nodes))

       ! set default: start with mode k=0.0, interval 1.0, up to n_nodes-1
       CALL set_Fourier_boundaries(grid, 0.0, 1.0)
       
    ELSE
       IF (p_n_nodes.GT.1) THEN
          ! we need at least two nodes for the boundaries
          grid%n_nodes = p_n_nodes
          ALLOCATE(grid%nodes(1:grid%n_nodes))
          
          ! set default interval 0 to 1, equidistant, periodic
          CALL set_boundaries(grid, 0.0, 1.0, .TRUE.)
       ELSE
          PRINT*, "Number of nodes must be at least 2!"
          STOP
       END IF
    END IF
  END SUBROUTINE mp_initialize

  !> Finalizes the instance of the datatype. The grid nodes are deallocated.
  SUBROUTINE mp_finalize(grid)
    type(Grid1D) :: grid
    IF (ALLOCATED(grid%nodes)) DEALLOCATE(grid%nodes)
  END SUBROUTINE mp_finalize

  SUBROUTINE grid1d_print_grid(this)
    type(Grid1D) :: this

    PRINT*, "n_nodes = ", this%n_nodes
    PRINT*, "left/right_boundary = ", this%left_boundary, this%right_boundary
    PRINT*, "isEquidistant = ", this%isEquidistant
    PRINT*, "isPeriodic    = ", this%isPeriodic
    PRINT*, "isFourierGrid = ", this%isFourierGrid
    PRINT*, "startindex = ", this%startindex
    PRINT*, "nodes = ", this%nodes
    
  END SUBROUTINE grid1d_print_grid

  !> This routine sets the boundary positions and the periodic flag of the grid.
  SUBROUTINE mp_set_boundaries(grid,p_left,p_right, p_periodic)
    TYPE(Grid1D) :: grid !<instance of Grid1D
    REAL :: p_left  !<position of the left boundary
    real :: p_right !<position of the right boundary
    LOGICAL :: p_periodic !<flag for periodic boundary

    ! Local variables
    integer :: iNode
    real :: dx

    grid%left_boundary = p_left
    grid%right_boundary = p_right
    grid%isEquidistant = .TRUE.
    grid%isPeriodic = p_periodic

    IF (grid%isPeriodic) THEN
       dx = (grid%right_boundary-grid%left_boundary)/DBLE(grid%n_nodes)
    ELSE
       dx = (grid%right_boundary-grid%left_boundary)/DBLE(grid%n_nodes-1)
    END IF
    DO iNode=1,grid%n_nodes
       grid%nodes(iNode) = grid%left_boundary+(iNode-1)*dx
    END DO
  END SUBROUTINE mp_set_boundaries

  SUBROUTINE grid1d_get_boundaries(grid,left_bnd, right_bnd)
    type(Grid1D) :: grid
    REAL,intent(OUT) :: left_bnd, right_bnd

    left_bnd = grid%left_boundary
    right_bnd = grid%right_boundary
  END SUBROUTINE grid1d_get_boundaries

  SUBROUTINE mp_set_fourier_boundaries(grid,k_start,k_diff)
    type(Grid1D) :: grid
    REAL :: k_start, k_diff

    ! Local variables
    integer :: iNode

    IF (.NOT.grid%isFourierGrid) THEN
       PRINT*," Switching grid to Fourier Grid. Better use initialize."
       grid%isFourierGrid  = .true.
    END IF

    grid%left_boundary  = k_start
    grid%right_boundary = k_start+grid%n_nodes*k_diff
    grid%isEquidistant  = .TRUE.
    grid%isPeriodic     = .TRUE.
    
    DO iNode=1,grid%n_nodes
       grid%nodes(iNode) = grid%left_boundary+(iNode-1)*k_diff
    END DO

  END SUBROUTINE mp_set_fourier_boundaries

  !> Sets the index of the first entry. Default is 1, as usual in Fortran.
  !! 
  !! Can be arbitrary, it is taken into account whenever an index is given
  !! as argument.
  !!
  SUBROUTINE grid1d_set_startindex(grid,startindex)
    TYPE(grid1D) :: grid
    INTEGER :: startindex

    grid%startindex = startindex
  END SUBROUTINE grid1d_set_startindex

  !> Returns the grid spacing, that is the distance between two neighboring grid nodes.
  !> This gives back a value only for an equidistant grid. For non-equidistant grids, -1 is returned.
  FUNCTION mp_get_gridspacing(grid) RESULT(gridspacing)
    TYPE(Grid1D) :: grid
    real :: gridspacing 

    IF (grid%isEquidistant) THEN
       !PRINT*,grid%nodes
       gridspacing = grid%nodes(2)-grid%nodes(1)
    ELSE
       gridspacing = -1.0
    END IF
  END FUNCTION mp_get_gridspacing

  !> Returns the number of nodes of the grid.
  FUNCTION mp_get_number_of_nodes(grid) RESULT(n_nodes)
    type(Grid1D) :: grid
    integer :: n_nodes
    
    n_nodes = grid%n_nodes
  END FUNCTION mp_get_number_of_nodes

  FUNCTION grid1d_get_startindex(grid) RESULT(sind)
    type(Grid1D) :: grid
    integer :: sind

    sind = grid%startindex
  END FUNCTION grid1d_get_startindex

  !> Returns the periodic flag of the grid.
  FUNCTION mp_is_periodic_boundary(grid) RESULT(is_periodic)
    type(Grid1D) :: grid
    logical :: is_periodic

    is_periodic = grid%isPeriodic
  END FUNCTION mp_is_periodic_boundary

  !> Returns a node.
  !!
  !! Returns the node with the index given as argument. Where the indices start
  !! can be set by the set_startindex subroutine.
  !!
  FUNCTION grid1d_get_node(grid,index) RESULT(node)
    TYPE(Grid1D),INTENT(IN) :: grid
    INTEGER,intent(IN) :: index
    real :: node
    
    node = grid%nodes(index-grid%startindex+1)
  END FUNCTION grid1d_get_node
END MODULE Grid1DModule
