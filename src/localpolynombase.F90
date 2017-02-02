#include "redef.h"
#include "intrinsic_sizes.h"
!****h* /localpolynombase_mod
! DESCRIPTION
! LocalPolynomBase defines base functions, using polynomials
! for the interpolation between three points. The polynomials
! are one at the center point and zero at the two other points. 
! The advantage of this representation is that grid point values
! of functions are identical to the coefficients of a base
! function representation of the same function. So we do not
! need to solve any system of equations for getting the coeffs. 
! This simplifies parallelization a lot.
!****
! 20.05.2008: cubic order for general grids work

MODULE localpolynombase_mod
  use Grid1DModule
  implicit none

  INTEGER,PARAMETER :: BC_PERIODIC=0,BC_NONPERIODIC=1
  !****t* localpolynombase_mod/LocalPolynomBase
  ! DESCRIPTION
  ! This type encapsulates a base of polynomials
  ! SOURCE
  TYPE,PUBLIC :: LocalPolynomBase
     private
     INTEGER :: degree ! polynomial degree
     integer :: nDeriv ! number of derivatives used
     !INTEGER :: nPoints ! number of points of underlying grid,
                        ! this is not necessarily the number of independent
                        ! base functions
     integer :: nCoeffs ! number if independent base functions
     REAL :: left_boundary,right_boundary
     INTEGER :: boundary_condition ! one of the above defined
     LOGICAL :: nodes_allocated = .FALSE.
     REAL,DIMENSION(:),ALLOCATABLE :: nodes ! the position of the grid points
  END TYPE LocalPolynomBase
  !****

  INTERFACE get_value
     module procedure mp_get_value
  END INTERFACE

  INTERFACE initialize_polybase
     MODULE PROCEDURE lpb_initialize_polybase_direct, lpb_initialize_polybase_grid1d
  END INTERFACE

  INTEGER,PRIVATE :: derivative_order
  PRIVATE :: base_poly  

  private :: mp_get_value
contains

  FUNCTION size_of_localpolynombase(n_nodes) RESULT(memory_need)
    integer :: n_nodes
    integer :: memory_need

    memory_need = 5*SIZE_OF_INTEGER &
         & + SIZE_OF_LOGICAL &
         & + (n_nodes+2)*SIZE_OF_REAL
  END FUNCTION size_of_localpolynombase

  !****m* localpolynombase_mod/initialize
  ! DESCRIPTION
  ! Initializes the LocalPolynomBase LPB by choosing the degree
  ! of the polynoms and giving the number of physical grid points.
  ! This function has to be called first before any of the set 
  ! functions.
  !
  ! SYNOPSIS
  SUBROUTINE lpb_initialize_polybase_direct(LPB,arg_degree,arg_deriv_ord,arg_nCoeffs,arg_bc)
  type(LocalPolynomBase) :: LPB
  INTEGER :: arg_degree, arg_deriv_ord, arg_nCoeffs, arg_bc
  ! ARGUMENTS
  ! LPB: the LocalPolynomBase
  ! arg_degree: the degree of the polynomials
  ! arg_nCoeffs: number of physical grid points 
  ! arg_bc: boundary condition

  ! SEE ALSO
  ! localpolynombase_mod/set_boundary
  ! localpolynombase_mod/set_nodes
  !****
  integer :: alloc_error

  LPB%degree = arg_degree
  LPB%nDeriv = INT(0.5*(LPB%degree-1))
  LPB%nCoeffs = arg_nCoeffs

  IF ((arg_bc==BC_PERIODIC).OR.(arg_bc.EQ.BC_NONPERIODIC)) THEN
     LPB%boundary_condition=arg_bc
  ELSEIF (arg_bc.EQ.2) THEN
     LPB%boundary_condition=BC_NONPERIODIC
  ELSE
     PRINT*,"wrong boundary condition, use BC_PERIODIC or BC_NONPERIODIC"
     STOP
  END IF

!  IF (LPB%boundary_condition == BC_PERIODIC) THEN
!     LPB%nPoints = LPB%nCoeffs+1
!  ELSE
!     LPB%nPoints = LPB%nCoeffs
!  END IF

  if (LPB%nodes_allocated) THEN
     DEALLOCATE(LPB%nodes)
  END IF
  ALLOCATE(LPB%nodes(1:LPB%nCoeffs),STAT=alloc_error)
  IF (alloc_error.eq.0) THEN
     LPB%nodes_allocated=.TRUE.
  ELSE
     STOP 'allocation error in localpolynombase initialize'
  ENDIF

  derivative_order=arg_deriv_ord

END SUBROUTINE lpb_initialize_polybase_direct

  !****m* localpolynombase_mod/initialize
  ! DESCRIPTION
  ! Initializes the LocalPolynomBase LPB by choosing the degree
  ! of the polynoms and giving the number of physical grid points.
  ! This function has to be called first before any of the set 
  ! functions.
  !
  ! SYNOPSIS
SUBROUTINE lpb_initialize_polybase_grid1d(LPB,arg_degree,arg_deriv_ord,arg_grid)
  type(LocalPolynomBase) :: LPB
  INTEGER :: arg_degree, arg_deriv_ord
  type(Grid1D) :: arg_grid
  ! ARGUMENTS
  ! LPB: the LocalPolynomBase
  ! arg_degree: the degree of the polynomials
  ! arg_deriv_ord: the order of the derivatives
  ! arg_grid: a class which describes the underlying grid

  ! SEE ALSO
  ! localpolynombase_mod/set_boundary
  ! localpolynombase_mod/set_nodes
  !****
  INTEGER :: alloc_error, startindex, i_node

  LPB%degree = arg_degree
  LPB%nDeriv = INT(0.5*(LPB%degree-1))

  LPB%nCoeffs = get_number_of_nodes(arg_grid)

  IF (is_periodic_boundary(arg_grid)) then
     LPB%boundary_condition=BC_PERIODIC
  ELSE
     LPB%boundary_condition=BC_NONPERIODIC
  END IF

  if (LPB%nodes_allocated) THEN
     DEALLOCATE(LPB%nodes)
  END IF
  ALLOCATE(LPB%nodes(1:LPB%nCoeffs),STAT=alloc_error)
  IF (alloc_error.eq.0) THEN
     LPB%nodes_allocated=.TRUE.
  ELSE
     STOP 'allocation error in localpolynombase initialize'
  ENDIF

  startindex = get_startindex(arg_grid)
  DO i_node=startindex,startindex+LPB%nCoeffs-1
     LPB%nodes(i_node-startindex+1) = get_node(arg_grid,i_node)
  END DO
  ! finally get the boundaries from the grid
  CALL get_boundaries(arg_grid,LPB%left_boundary,LPB%right_boundary)

  derivative_order=arg_deriv_ord

END SUBROUTINE lpb_initialize_polybase_grid1d


SUBROUTINE print_base(LPB)
  TYPE(LocalPolynomBase) :: LPB

  PRINT*,"LocalPolynomBase with attributes:"
  PRINT*,"  degree = ", LPB%degree
  PRINT*,"  nDeriv = ", LPB%nDeriv
  !PRINT*,"  nPoints= ", LPB%nPoints
  PRINT*,"  nCoeffs= ", LPB%nCoeffs
  PRINT*,"  (left_boundary,right_boundary) = ( ",LPB%left_boundary,",",LPB%right_boundary,")"
  WRITE(*,"(A)", advance = "NO") "  boundary_condition = "
  IF (LPB%boundary_condition.EQ.BC_PERIODIC) THEN
     WRITE(*,"(A)") "BC_PERIODIC"
  ELSE
     WRITE(*,"(A)") "BC_NONPERIODIC"
  END IF
  IF (LPB%nodes_allocated) THEN
     PRINT*,"  nodes = ",LPB%nodes
  ELSE
     PRINT*,"  nodes are not allocated"
  END IF
  
END SUBROUTINE print_base
!****m* localpolynombase_mod/set_nodes
! SYNOPSIS
SUBROUTINE set_nodes(LPB,in_nodes)
  TYPE(LocalPolynomBase) :: LPB
  REAL, dimension(:) :: in_nodes
!****

  real :: grid_spacing

  IF (SIZE(in_nodes).EQ.LPB%nCoeffs) THEN
     LPB%nodes = in_nodes
  ELSE
     PRINT*,"dimension of in_nodes must equal nCoeffs, which is ",LPB%nCoeffs
     STOP
  END IF

  LPB%left_boundary = LPB%nodes(1)
  IF (LPB%boundary_condition.EQ.BC_PERIODIC) THEN
     ! the right boundary is not the last point, but is one gridsize further
     grid_spacing = LPB%nodes(2)-LPB%nodes(1)
     LPB%right_boundary = LPB%nodes(LPB%nCoeffs)+grid_spacing
  ELSE
     LPB%right_boundary = LPB%nodes(LPB%nCoeffs)
  END IF

END SUBROUTINE set_nodes

!****m* localpolynombase_mod/get_nPoints
! SYNOPSIS
FUNCTION get_nPoints(LPB) RESULT(res)
  TYPE(LocalPolynomBase), intent(IN) :: LPB
  INTEGER :: res
!****
  res = LPB%nCoeffs
END FUNCTION get_nPoints

!****m* localpolynombase_mod/get_nDeriv
! SYNOPSIS
FUNCTION get_nDeriv(LPB) RESULT(res)
  TYPE(LocalPolynomBase), intent(IN) :: LPB
  INTEGER :: res
!****
  res = LPB%nDeriv
END FUNCTION get_nDeriv

!****m* localpolynombase_mod/get_value
! SYNOPSIS
SUBROUTINE mp_get_value(LPB,ind,xpos_in, res)
  TYPE(LocalPolynomBase),INTENT(IN) :: LPB
  INTEGER, INTENT(IN) :: ind
  REAL, INTENT(IN) :: xpos_in
  REAL,DIMENSION(0:LPB%nDeriv),INTENT(OUT) :: res

  ! ARGUMENTS
  ! LPB: the LocalPolynomBase
  ! ind: the index of the base function (it is the index of the centernode)
  !      must be in [1,nCoeffs]
  ! xpos: the position, at which we want to know the function value
  !
  ! RETURN VALUE
  ! res: An array of dimension 0.5*(LPB%degree-1). This is the number of derivatives
  !      used in the approximation
  !****
  INTEGER :: nfac,ideriv
  REAL :: xpos,lb,cp,rb

  ! if the index of the base function is out of the 
  ! range of independent base functions, show an error
  IF ((ind<1).OR.(ind>LPB%nCoeffs)) THEN
     PRINT*,"ind must be between 1 and nCoeffs"
     STOP
  END IF
  
  IF (.NOT.(SIZE(res).EQ.(LPB%nDeriv+1))) THEN
     PRINT*,"res must have the dimension (degree-1)/2"
     STOP
  END IF
  
  IF (LPB%boundary_condition.EQ.BC_PERIODIC) THEN
     ! map back outlying points
     IF (xpos_in<LPB%nodes(1)) THEN
        nfac=INT( (LPB%nodes(1)-xpos_in)/(LPB%right_boundary-LPB%left_boundary) )
        xpos = xpos_in + (nfac+1)*(LPB%right_boundary-LPB%left_boundary)
     ELSEIF (xpos_in>=LPB%right_boundary) THEN
        nfac=INT( (xpos_in-LPB%right_boundary)/(LPB%right_boundary-LPB%left_boundary) )
        xpos = xpos_in - (nfac+1)*(LPB%right_boundary-LPB%left_boundary)
     ELSE
        xpos = xpos_in
     END IF
  ELSEIF (LPB%boundary_condition.EQ.BC_NONPERIODIC) THEN
     IF ((xpos_in.le.(LPB%left_boundary-(LPB%nodes(2)-LPB%nodes(1)))).OR.&
          (xpos_in.ge.(LPB%right_boundary+(LPB%right_boundary-LPB%nodes(LPB%nCoeffs-1))))) THEN
        ! boundary zero boundary condition
        res=0.
        return
     ELSE 
        xpos = xpos_in
     END IF
  END IF

  !x value of left, central and right points
  if (ind.gt.1) then
     lb=LPB%nodes(ind-1)
  elseif ((ind.eq.1).AND.(LPB%boundary_condition.eq.BC_NONPERIODIC)) then
     lb=LPB%nodes(ind)-(LPB%nodes(ind+1)-LPB%nodes(ind))
  end if
  cp=LPB%nodes(ind)
  if (ind.lt.LPB%nCoeffs) then
     rb=LPB%nodes(ind+1)
  elseif ((ind.eq.LPB%nCoeffs).AND.(LPB%boundary_condition.eq.BC_NONPERIODIC)) then
     rb=LPB%nodes(ind) + (LPB%nodes(ind)-LPB%nodes(ind-1))
  ELSEIF ((ind.EQ.LPB%nCoeffs).AND.(LPB%boundary_condition.EQ.BC_PERIODIC)) THEN
     rb = LPB%right_boundary
  end if

  if((ind.gt.1).and.(lb.le.xpos).and.(xpos.lt.cp)) THEN 
     !left of the central point, lb<=x<cp
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp,lb,xpos)
     END DO
  elseif((ind.eq.1).and.(LPB%boundary_condition.eq.BC_PERIODIC).and.(xpos.ge.LPB%nodes(LPB%nCoeffs))) then
     !exploit periodicity for the contributions left of the first point
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,LPB%right_boundary,LPB%nodes(LPB%nCoeffs),xpos)
     END DO
  elseif((ind.eq.1).and.(LPB%boundary_condition.eq.BC_NONPERIODIC).and.(lb.le.xpos).and.(xpos.lt.cp)) then
     ! Outside left boundary, here one could apply a boundary condition
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp,lb,xpos)
     END DO
  elseif ((ind.le.LPB%nCoeffs).and.(cp.le.xpos).and.(xpos.lt.rb)) then 
     !right of the central point, cp<=x<rb
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp,rb,xpos)
     END DO
  elseif ((ind.eq.LPB%nCoeffs).and.(LPB%boundary_condition.eq.BC_NONPERIODIC)&
       &.and.(cp.le.xpos).and.(xpos.lt.rb)) then 
     !outside right boundary, here one could apply a boundary condition
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp,rb,xpos)
     END DO
  else
     ! the base function has only a support on a 2*dx interval
     res = 0.0D0
  endif
  
END SUBROUTINE mp_get_value

!****im* localpolynombase_mod/base_poly
! SYNOPSIS
FUNCTION base_poly(ideriv,degree,x0,x1,x)
  INTEGER :: ideriv, degree
  REAL :: x1, x0,x
  real:: xp, x0p
  real :: base_poly
  ! ARGUMENTS
  ! ideriv: which derivative you want to use. This translates into
  !         different base functions
  ! degree: the degree of the polynomial
  ! x0:     the central point
  ! x1:     is the next left point or the next right point. As the polynomial is identical, we
  !         put both sides together in one function.
  !
  ! ATTRIBUTES
  ! private
  !****

  !coordinate transformation to improve numerical precision
  xp=x-x1
  x0p=x0-x1

  SELECT CASE (degree)
     CASE (1)
        IF (ideriv.GT.0) THEN 
           PRINT*,"Degree 1 polynoms don't use derivatives."
           STOP
        END IF
        
     CASE (3)
        SELECT CASE (ideriv)
        CASE (0) ! the function value
           base_poly = xp**2*(3*x0p-2*xp)/x0p**3
        CASE (1) ! first derivative
           base_poly = (xp - x0p)*xp**2/x0p**2
        CASE default
           PRINT*,"Degree 3 polynoms only use first derivatives."
           STOP
        END SELECT

     CASE (5)
        SELECT CASE (ideriv)
        CASE (0) ! the function value
           base_poly = xp**3*(6*xp**2 - 15*xp*x0p + 10*x0p**2)/x0p**5
        CASE (1)
           base_poly = -((xp - x0p)*xp**3*(3*xp - 4*x0p))/x0p**4
        CASE (2)
           base_poly = ((xp - x0p)**2*xp**3)/(2.*x0p**3)
        CASE default
           PRINT*,"Degree 5 polynoms only uses first two derivatives."
           STOP
        END SELECT
  END SELECT
END FUNCTION base_poly

SUBROUTINE finalize_localpolynombase(LPB)
  type(LocalPolynomBase) :: LPB

  if (LPB%nodes_allocated) THEN
     DEALLOCATE(LPB%nodes)
     LPB%nodes_allocated=.FALSE.
  END IF

END SUBROUTINE finalize_localpolynombase

END MODULE localpolynombase_mod
