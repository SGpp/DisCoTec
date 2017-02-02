#include "redef.h"
#include "intrinsic_sizes.h"
MODULE DerivativeMatrixModule
  use Grid1DModule
  USE BandedMatrixModule

  implicit none

  TYPE DerivativeMatrix
     type(BandedMatrix) :: Data
     integer :: derivative_order !the order of the finite differences formulas for the derivatives
     logical :: periodic_boundary
     TYPE(Grid1D), pointer :: grid
  END TYPE DerivativeMatrix

  TYPE DerivativeMatrixReal
     type(BandedMatrixReal) :: Data
     integer :: derivative_order
     logical :: periodic_boundary
     TYPE(Grid1D), pointer :: grid
  END TYPE DerivativeMatrixReal

  INTERFACE initialize
     module procedure mp_initialize_matrix, mp_initialize_matrix_real
  END INTERFACE

  INTERFACE finalize
     module procedure mp_finalize_matrix, mp_finalize_matrix_real
  END INTERFACE

  INTERFACE show
     MODULE PROCEDURE mp_show_on_screen, mp_show_in_file, mp_show_on_screen_real, mp_show_in_file_real
  END INTERFACE

  INTERFACE mat_get_value
     module procedure mp_mat_get_value, mp_mat_get_value_real
  END INTERFACE

  PRIVATE :: mp_initialize_matrix, mp_finalize_matrix, mp_initialize_matrix_real, mp_finalize_matrix_real
  private :: mp_mat_get_value, mp_mat_get_value_real

CONTAINS

  FUNCTION static_size_of_DerivativeMatrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = static_size_of_BandedMatrix()&
         & + SIZE_OF_INTEGER&
         & + SIZE_OF_LOGICAL
    
  END FUNCTION static_size_of_DerivativeMatrix

  SUBROUTINE mp_initialize_module
  END SUBROUTINE mp_initialize_module

  !> Initialize Derivative matrices
  SUBROUTINE mp_initialize_matrix(dmat,grid,p_derivative_order,transposed)
    type(DerivativeMatrix),intent(out) :: dmat !> derivative matrix
    TYPE(Grid1D),TARGET,intent(in) :: grid     !> grid definition
    INTEGER,intent(in) :: p_derivative_order   !> derivative order
    LOGICAL, OPTIONAL :: transposed

    ! Local variables
    integer :: n_points

    n_points = get_number_of_nodes(grid)

    IF (PRESENT(transposed)) THEN
       CALL initialize(dmat%Data,n_points,n_points,transposed)
    ELSE
       CALL initialize(dmat%Data,n_points,n_points)
    END IF
    
    CALL allocate(dmat%Data)
    dmat%derivative_order = p_derivative_order
    dmat%periodic_boundary = is_periodic_boundary(grid)

    dmat%grid => grid
  END SUBROUTINE mp_initialize_matrix

  SUBROUTINE mp_finalize_matrix(dmat)
    TYPE(DerivativeMatrix) :: dmat
    CALL finalize(dmat%Data)
  END SUBROUTINE mp_finalize_matrix

  !> Calculate derivative matrices
  SUBROUTINE calculate(dmat, which_derivative, rad_bc_type)
    TYPE(DerivativeMatrix) :: dmat             !> matrix to contain the derivative matrix
    INTEGER, intent(IN) :: which_derivative    !> number of the derivative to calculate the matrix for
    !                   can be 1 or 2 for the first or second derivative
    INTEGER,intent(in) :: rad_bc_type          !> radial boundary condition

    !local variables
    REAL :: dx
    REAL, DIMENSION(:),ALLOCATABLE :: sten
    INTEGER :: i, s, n_points

    PERFON('deriv_mat')

    call set_zero(dmat%Data)

    IF (dmat%grid%isEquidistant) THEN

       ALLOCATE(sten(-dmat%derivative_order/2:dmat%derivative_order/2))
       dx = get_gridspacing(dmat%grid)
       
       SELECT CASE (dmat%derivative_order)
       CASE(2)
          !PRINT*,"in calculate, derivative_order=2"
          SELECT CASE (which_derivative)
          CASE(1) ! first derivative
             !PRINT*,"in calculate, derivative_order=2, first derivative"
             sten = (/-0.5,0.0,0.5/) / dx
          CASE(2) ! second derivative
             sten = (/1.0,-2.0,1.0/) / dx**2
          END SELECT
       CASE(4)
          SELECT CASE (which_derivative)
          CASE(1)
             sten = (/1.0,-8.0,0.0,8.0,-1.0/) / dx
          CASE(2)
             sten = (/-1.0,16.0,-30.0,16.0,-1.0/) / dx**2
          END SELECT
          sten=sten/12.0
       CASE(6)
          SELECT CASE (which_derivative)
          CASE(1)
             sten = (/-1.0,9.0,-45.0,0.0,45.0,-9.0,1.0/) / dx
          case default
             PRINT*,"for sixth order derivative, the second derivative is not yet implemented."
             STOP
          END SELECT
          sten = sten/60.0
       CASE default
          PRINT*,"only derivation order 2, 4 and 6 implemented"
          STOP
       END SELECT

       !PRINT*, sten
       n_points = get_number_of_nodes(dmat%grid)
       DO s=-dmat%derivative_order/2,dmat%derivative_order/2 
          DO i=1,n_points
             IF (((i+s).LT.1).OR.((i+s).GT.n_points)) THEN !outside box
                IF (rad_bc_type.EQ.0) THEN !periodic b.c.
                   IF (is_periodic_boundary(dmat%grid)) THEN
                      ! The last and the first point are identical and lie on the
                      ! left or right boundary, respectively. So the derivatives 
                      ! only run over the nintervals columns and rows of the matrix.
                      CALL add_value(dmat%Data,i,MOD(i+s-1+n_points,n_points)+1, sten(s))
                   ELSE
                      STOP 'rad_bc_type is 0, but the grid is nonperiodic'
                   END IF
                ELSEIF (rad_bc_type.EQ.1) THEN !Dirichlet b.c. at both boundaries
                   !do nothing (zero's outside box)
                ELSEIF (rad_bc_type.EQ.2) THEN !v. Neumann at inner and Dirichlet b.c. 
                   !at outer boundary, i.e. mirror grid points at inner boundary
                   IF ((i+s).LT.1) THEN
                      Call add_value(dmat%Data,i,1-(i+s), sten(s))
                   ENDIF
                ELSE
                   STOP 'rad_bc_type>2 not implemented in derivative matrix construction '
                ENDIF
             ELSE !inside box
                CALL add_value(dmat%Data,i,i+s,sten(s))
             ENDIF
          END DO
       END DO
       call commit_values(dmat%Data)
       DEALLOCATE(sten)

    ELSE
       ! not equidistant underlying grid
       ! this is for the future and has to be implemented
    END IF

    PERFOFF
  END SUBROUTINE calculate

  SUBROUTINE mp_show_on_screen(dmat)
    type(DerivativeMatrix) :: dmat
    
    CALL show(dmat%Data)
  END SUBROUTINE mp_show_on_screen

  SUBROUTINE mp_show_in_file(dmat,filename)
    type(DerivativeMatrix) :: dmat
    CHARACTER(len=*) :: filename

    CALL show(dmat%Data,filename)
  END SUBROUTINE mp_show_in_file

  FUNCTION mp_mat_get_value(dmat,irow,icol) RESULT(value)
    type(DerivativeMatrix) :: dmat
    INTEGER :: irow,icol
    complex :: value

    value = mat_get_value(dmat%Data,irow,icol)
  END FUNCTION mp_mat_get_value

!----------------------------------------------
!   Repeat everything for REAL datatypes
!----------------------------------------------

!> Initialize Derivative matrices
  SUBROUTINE mp_initialize_matrix_real(dmat,grid,p_derivative_order,transposed)
    type(DerivativeMatrixReal),intent(out) :: dmat !> derivative matrix
    TYPE(Grid1D),TARGET,intent(in) :: grid     !> grid definition
    INTEGER,intent(in) :: p_derivative_order   !> derivative order
    LOGICAL, OPTIONAL :: transposed

    ! Local variables
    integer :: n_points

    n_points = get_number_of_nodes(grid)

    IF (PRESENT(transposed)) THEN
       CALL initialize(dmat%Data,n_points,n_points,transposed)
    ELSE
       CALL initialize(dmat%Data,n_points,n_points)
    END IF
    
    CALL allocate(dmat%Data)
    dmat%derivative_order = p_derivative_order
    dmat%periodic_boundary = is_periodic_boundary(grid)

    dmat%grid => grid
  END SUBROUTINE mp_initialize_matrix_real

  SUBROUTINE mp_finalize_matrix_real(dmat)
    TYPE(DerivativeMatrixReal) :: dmat
    CALL finalize(dmat%Data)
  END SUBROUTINE mp_finalize_matrix_real

  !> Calculate derivative matrices
  SUBROUTINE calculate_real(dmat, which_derivative, rad_bc_type)
    TYPE(DerivativeMatrixReal) :: dmat             !> matrix to contain the derivative matrix
    INTEGER, intent(IN) :: which_derivative    !> number of the derivative to calculate the matrix for
    !                   can be 1 or 2 for the first or second derivative
    INTEGER,intent(in) :: rad_bc_type          !> radial boundary condition

    !local variables
    REAL :: dx
    REAL, DIMENSION(:),ALLOCATABLE :: sten
    INTEGER :: i, s, n_points

    PERFON('deriv_mat')

    call set_zero(dmat%Data)

    IF (dmat%grid%isEquidistant) THEN

       ALLOCATE(sten(-dmat%derivative_order/2:dmat%derivative_order/2))
       dx = get_gridspacing(dmat%grid)
       
       SELECT CASE (dmat%derivative_order)
       CASE(2)
          !PRINT*,"in calculate, derivative_order=2"
          SELECT CASE (which_derivative)
          CASE(1) ! first derivative
             !PRINT*,"in calculate, derivative_order=2, first derivative"
             sten = (/-0.5,0.0,0.5/) / dx
          CASE(2) ! second derivative
             sten = (/1.0,-2.0,1.0/) / dx**2
          END SELECT
       CASE(4)
          SELECT CASE (which_derivative)
          CASE(1)
             sten = (/1.0,-8.0,0.0,8.0,-1.0/) / dx
          CASE(2)
             sten = (/-1.0,16.0,-30.0,16.0,-1.0/) / dx**2
          END SELECT
          sten=sten/12.0
       CASE(6)
          SELECT CASE (which_derivative)
          CASE(1)
             sten = (/-1.0,9.0,-45.0,0.0,45.0,-9.0,1.0/) / dx
          case default
             PRINT*,"for sixth order derivative, the second derivative is not yet implemented."
             STOP
          END SELECT
          sten = sten/60.0
       CASE default
          PRINT*,"only derivation order 2, 4 and 6 implemented"
          STOP
       END SELECT

       !PRINT*, sten
       n_points = get_number_of_nodes(dmat%grid)
       DO s=-dmat%derivative_order/2,dmat%derivative_order/2 
          DO i=1,n_points
             IF (((i+s).LT.1).OR.((i+s).GT.n_points)) THEN !outside box
                IF (rad_bc_type.EQ.0) THEN !periodic b.c.
                   IF (is_periodic_boundary(dmat%grid)) THEN
                      ! The last and the first point are identical and lie on the
                      ! left or right boundary, respectively. So the derivatives 
                      ! only run over the nintervals columns and rows of the matrix.
                      CALL add_value(dmat%Data,i,MOD(i+s-1+n_points,n_points)+1, sten(s))
                   ELSE
                      STOP 'rad_bc_type is 0, but the grid is nonperiodic'
                   END IF
                ELSEIF (rad_bc_type.EQ.1) THEN !Dirichlet b.c. at both boundaries
                   !do nothing (zero's outside box)
                ELSEIF (rad_bc_type.EQ.2) THEN !v. Neumann at inner and Dirichlet b.c. 
                   !at outer boundary, i.e. mirror grid points at inner boundary
                   IF ((i+s).LT.1) THEN
                      Call add_value(dmat%Data,i,1-(i+s), sten(s))
                   ENDIF
                ELSE
                   STOP 'rad_bc_type>2 not implemented in derivative matrix construction '
                ENDIF
             ELSE !inside box
                CALL add_value(dmat%Data,i,i+s,sten(s))
             ENDIF
          END DO
       END DO
       call commit_values(dmat%Data)
       DEALLOCATE(sten)

    ELSE
       ! not equidistant underlying grid
       ! this is for the future and has to be implemented
    END IF

    PERFOFF
  END SUBROUTINE calculate_real

  SUBROUTINE mp_show_on_screen_real(dmat)
    type(DerivativeMatrixReal) :: dmat
    
    CALL show(dmat%Data)
  END SUBROUTINE mp_show_on_screen_real

  SUBROUTINE mp_show_in_file_real(dmat,filename)
    type(DerivativeMatrixReal) :: dmat
    CHARACTER(len=*) :: filename

    CALL show(dmat%Data,filename)
  END SUBROUTINE mp_show_in_file_real

  FUNCTION mp_mat_get_value_real(dmat,irow,icol) RESULT(value)
    type(DerivativeMatrixReal) :: dmat
    INTEGER :: irow,icol
    real :: value

    value = mat_get_value(dmat%Data,irow,icol)
  END FUNCTION mp_mat_get_value_real

END MODULE DerivativeMatrixModule
