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

MODULE localpolynombase2d_mod
  use Grid1DModule
  implicit none

  INTEGER,PARAMETER :: BC_PERIODIC=0,BC_NONPERIODIC=1
  !****t* localpolynombase_mod/LocalPolynomBase
  ! DESCRIPTION
  ! This type encapsulates a base of polynomials
  ! SOURCE
  TYPE,PUBLIC :: LocalPolynomBase2d
     private
     INTEGER :: degree ! polynomial degree
     integer :: nDeriv ! number of derivatives used
     !INTEGER :: nPoints ! number of points of underlying grid,
                        ! this is not necessarily the number of independent
                        ! base functions
     integer :: nCoeffs_x,nCoeffs_y ! number if independent base functions
     REAL :: left_boundary_x,right_boundary_x
     REAL :: left_boundary_y,right_boundary_y
     INTEGER :: boundary_condition_x,boundary_condition_y ! one of the above defined
     LOGICAL :: nodes_allocated = .FALSE.
     REAL,DIMENSION(:),ALLOCATABLE :: nodes_x,nodes_y ! the position of the grid points
  END TYPE LocalPolynomBase2d
  !****

  INTERFACE get_value2d
     module procedure mp_get_value2d
  END INTERFACE

  INTERFACE initialize_polybase2d
  !   MODULE PROCEDURE lpb_initialize_polybase_direct, lpb_initialize_polybase_grid2d
      MODULE PROCEDURE lpb_initialize_polybase_direct, lpb_initialize_polybase_grid2d
  END INTERFACE

  INTEGER,PRIVATE :: derivative_order
  PRIVATE :: base_poly  

  private :: mp_get_value2d
contains

  FUNCTION size_of_localpolynombase(n_nodes_x,n_nodes_y) RESULT(memory_need)
    integer :: n_nodes_x,n_nodes_y
    integer :: memory_need

    memory_need = 5*SIZE_OF_INTEGER &
         & + SIZE_OF_LOGICAL &
         & + (n_nodes_x+n_nodes_y+2)*SIZE_OF_REAL
  END FUNCTION size_of_localpolynombase

  !****m* localpolynombase_mod/initialize
  ! DESCRIPTION
  ! Initializes the LocalPolynomBase LPB by choosing the degree
  ! of the polynoms and giving the number of physical grid points.
  ! This function has to be called first before any of the set 
  ! functions.
  !
  ! SYNOPSIS
  SUBROUTINE lpb_initialize_polybase_direct(LPB,arg_degree,arg_deriv_ord,arg_nCoeffs_x&
                                          &,arg_nCoeffs_y,arg_bc_x,arg_bc_y)
  type(LocalPolynomBase2d) :: LPB
  INTEGER :: arg_degree, arg_deriv_ord, arg_nCoeffs_x,arg_nCoeffs_y, arg_bc_x,arg_bc_y
  ! ARGUMENTS
  ! LPB: the LocalPolynomBase
  ! arg_degree: the degree of the polynomials
  ! arg_nCoeffs_x: number of physical grid points in x
  ! arg_nCoeffs_y: number of physical grid points in y
  ! arg_bc_x: boundary condition in x
  ! arg_bc_y: boundary condition in y

  ! SEE ALSO
  ! localpolynombase_mod/set_boundary
  ! localpolynombase_mod/set_nodes
  !****
  integer :: alloc_error

  LPB%degree = arg_degree
  LPB%nDeriv = INT(0.5*(LPB%degree-1))
  LPB%nCoeffs_x = arg_nCoeffs_x
  LPB%nCoeffs_y = arg_nCoeffs_y

  IF ((arg_bc_x.EQ.BC_PERIODIC).OR.(arg_bc_x.EQ.BC_NONPERIODIC)) THEN
     LPB%boundary_condition_x=arg_bc_x
  ELSEIF (arg_bc_x.EQ.2) THEN
     LPB%boundary_condition_x=BC_NONPERIODIC
  ELSE
     PRINT*,"wrong boundary condition in x, use BC_PERIODIC or BC_NONPERIODIC"
     STOP
  END IF

  IF ((arg_bc_y.EQ.BC_PERIODIC).OR.(arg_bc_y.EQ.BC_NONPERIODIC)) THEN
     LPB%boundary_condition_y=arg_bc_y
  ELSEIF (arg_bc_y.EQ.2) THEN
     LPB%boundary_condition_y=BC_NONPERIODIC
  ELSE
     PRINT*,"wrong boundary condition in y, use BC_PERIODIC or BC_NONPERIODIC"
     STOP
  END IF

  if (LPB%nodes_allocated) THEN
     DEALLOCATE(LPB%nodes_x)
     DEALLOCATE(LPB%nodes_y)
  END IF
  ALLOCATE(LPB%nodes_x(1:LPB%nCoeffs_x),STAT=alloc_error)
  ALLOCATE(LPB%nodes_y(1:LPB%nCoeffs_y),STAT=alloc_error)

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
SUBROUTINE lpb_initialize_polybase_grid2d(LPB,arg_degree,arg_deriv_ord,arg_grid_x,arg_grid_y)
  type(LocalPolynomBase2d) :: LPB
  INTEGER :: arg_degree, arg_deriv_ord
  type(Grid1D) :: arg_grid_x,arg_grid_y
  ! ARGUMENTS
  ! LPB: the LocalPolynomBase
  ! arg_degree: the degree of the polynomials
  ! arg_deriv_ord: the order of the derivatives
  ! arg_grid: a class which describes the underlying grid

  ! SEE ALSO
  ! localpolynombase_mod/set_boundary
  ! localpolynombase_mod/set_nodes
  !****
  INTEGER :: alloc_error, startindex_x,startindex_y,i_node_x,i_node_y

  LPB%degree = arg_degree
  LPB%nDeriv = INT(0.5*(LPB%degree-1))

  LPB%nCoeffs_x = get_number_of_nodes(arg_grid_x)
  LPB%nCoeffs_y = get_number_of_nodes(arg_grid_y)

  IF (is_periodic_boundary(arg_grid_x)) then
     LPB%boundary_condition_x=BC_PERIODIC
  ELSE
     LPB%boundary_condition_x=BC_NONPERIODIC
  END IF

  IF (is_periodic_boundary(arg_grid_y)) then
     LPB%boundary_condition_y=BC_PERIODIC
  ELSE
     LPB%boundary_condition_y=BC_NONPERIODIC
  END IF

  if (LPB%nodes_allocated) THEN
     DEALLOCATE(LPB%nodes_x)
     DEALLOCATE(LPB%nodes_y)
  END IF
  ALLOCATE(LPB%nodes_x(1:LPB%nCoeffs_x),STAT=alloc_error)
  ALLOCATE(LPB%nodes_y(1:LPB%nCoeffs_y),STAT=alloc_error)
  IF (alloc_error.eq.0) THEN
     LPB%nodes_allocated=.TRUE.
  ELSE
     STOP 'allocation error in localpolynombase initialize'
  ENDIF

  startindex_x = get_startindex(arg_grid_x)
  startindex_y = get_startindex(arg_grid_y)

  DO i_node_x=startindex_x,startindex_x+LPB%nCoeffs_x-1
     DO i_node_y=startindex_y,startindex_y+LPB%nCoeffs_y-1
        LPB%nodes_x(i_node_x-startindex_x+1) = get_node(arg_grid_x,i_node_x)
        LPB%nodes_y(i_node_y-startindex_y+1) = get_node(arg_grid_y,i_node_y)
     END DO
  END DO

  ! finally get the boundaries from the grid

  CALL get_boundaries(arg_grid_x,LPB%left_boundary_x,LPB%right_boundary_x)
  CALL get_boundaries(arg_grid_y,LPB%left_boundary_y,LPB%right_boundary_y)

  !CALL get_boundaries_2d(arg_grid_x,LPB%left_boundary_x,LPB%right_boundary_x,&
  !     &arg_grid_y,LPB%left_boundary_y,LPB%right_boundary_y)

  derivative_order=arg_deriv_ord

END SUBROUTINE lpb_initialize_polybase_grid2d


SUBROUTINE print_base(LPB)
  TYPE(LocalPolynomBase2d) :: LPB

  PRINT*,"LocalPolynomBase with attributes:"
  PRINT*,"  degree = ", LPB%degree
  PRINT*,"  nDeriv = ", LPB%nDeriv
  !PRINT*,"  nPoints= ", LPB%nPoints
  PRINT*,"  nCoeffs_x= ", LPB%nCoeffs_x
  PRINT*,"  nCoeffs_y= ", LPB%nCoeffs_y
  PRINT*,"  (left_boundary_x,right_boundary_x) = ( ",LPB%left_boundary_x,",",LPB%right_boundary_x,")"
  PRINT*,"  (left_boundary_y,right_boundary_y) = ( ",LPB%left_boundary_y,",",LPB%right_boundary_y,")"
  WRITE(*,"(A)", advance = "NO") "  boundary_condition_x = "
  IF (LPB%boundary_condition_x.EQ.BC_PERIODIC) THEN
     WRITE(*,"(A)") "BC_PERIODIC"
  ELSE
     WRITE(*,"(A)") "BC_NONPERIODIC"
  END IF
  WRITE(*,"(A)", advance = "NO") "  boundary_condition_y = "
  IF (LPB%boundary_condition_y.EQ.BC_PERIODIC) THEN
     WRITE(*,"(A)") "BC_PERIODIC"
  ELSE
     WRITE(*,"(A)") "BC_NONPERIODIC"
  END IF
  IF (LPB%nodes_allocated) THEN
     PRINT*,"  nodes_x = ",LPB%nodes_x
     PRINT*,"  nodes_y = ",LPB%nodes_y
  ELSE
     PRINT*,"  nodes are not allocated"
  END IF

END SUBROUTINE print_base
!****m* localpolynombase_mod/set_nodes
! SYNOPSIS
SUBROUTINE set_nodes(LPB,in_nodes_x,in_nodes_y)
  TYPE(LocalPolynomBase2d) :: LPB
  REAL, dimension(:) :: in_nodes_x,in_nodes_y
!****

  real :: grid_spacing_x,grid_spacing_y

  IF (SIZE(in_nodes_x).EQ.LPB%nCoeffs_x) THEN
     LPB%nodes_x = in_nodes_x
  ELSE
     PRINT*,"dimension of in_nodes_x must equal nCoeffs, which is ",LPB%nCoeffs_x," by ",LPB%nCoeffs_y
     STOP
  END IF

  IF (SIZE(in_nodes_y).EQ.LPB%nCoeffs_y) THEN
     LPB%nodes_y = in_nodes_y
  ELSE
     PRINT*,"dimension of in_nodes_y must equal nCoeffs, which is ",LPB%nCoeffs_x," by ",LPB%nCoeffs_y
     STOP
  END IF

  LPB%left_boundary_x = LPB%nodes_x(1)
  LPB%left_boundary_y = LPB%nodes_y(1)
  IF (LPB%boundary_condition_x.EQ.BC_PERIODIC) THEN
     ! the right boundary is not the last point, but is one gridsize further
     grid_spacing_x = LPB%nodes_x(2)-LPB%nodes_x(1)
     LPB%right_boundary_x = LPB%nodes_x(LPB%nCoeffs_x)+grid_spacing_x
  ELSE
     LPB%right_boundary_x = LPB%nodes_x(LPB%nCoeffs_x)
  END IF

  IF (LPB%boundary_condition_y.EQ.BC_PERIODIC) THEN
     ! the right boundary is not the last point, but is one gridsize further
     grid_spacing_y = LPB%nodes_y(2)-LPB%nodes_y(1)
     LPB%right_boundary_y = LPB%nodes_y(LPB%nCoeffs_y)+grid_spacing_y
  ELSE
     LPB%right_boundary_y = LPB%nodes_y(LPB%nCoeffs_y)
  END IF

END SUBROUTINE set_nodes

!****m* localpolynombase_mod/get_nPoints
! SYNOPSIS
FUNCTION get_nPoints_x(LPB) RESULT(res)
  TYPE(LocalPolynomBase2d), intent(IN) :: LPB
  INTEGER :: res
!****
  res = LPB%nCoeffs_x
END FUNCTION get_nPoints_x

FUNCTION get_nPoints_y(LPB) RESULT(res)
  TYPE(LocalPolynomBase2d), intent(IN) :: LPB
  INTEGER :: res
!****
  res = LPB%nCoeffs_y
END FUNCTION get_nPoints_y

!****m* localpolynombase_mod/get_nDeriv
! SYNOPSIS
FUNCTION get_nDeriv(LPB) RESULT(res)
  TYPE(LocalPolynomBase2d), intent(IN) :: LPB
  INTEGER :: res
!****
  res = LPB%nDeriv
END FUNCTION get_nDeriv

!****m* localpolynombase_mod/get_value
! SYNOPSIS
SUBROUTINE mp_get_value2d(LPB,ind_x,xpos_in,ind_y,ypos_in,res)
  TYPE(LocalPolynomBase2d),INTENT(IN) :: LPB
  INTEGER, INTENT(IN) :: ind_x,ind_y
  REAL, INTENT(IN) :: xpos_in,ypos_in
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
  REAL :: xpos,ypos,lb_y,cp_y,rb_y,lb_x,cp_x,rb_x

  ! if the index of the base function is out of the 
  ! range of independent base functions, show an error
  IF ((ind_x<1).OR.(ind_x>LPB%nCoeffs_x)) THEN
     PRINT*,"ind_x must be between 1 and nCoeffs_x"
     STOP
  END IF

  IF ((ind_y<1).OR.(ind_y>LPB%nCoeffs_y)) THEN
     PRINT*,"ind_y must be between 1 and nCoeffs_y"
     STOP
  END IF
  
!  IF (.NOT.(SIZE(res).EQ.(LPB%nDeriv+1))) THEN
!     PRINT*,"res must have the dimension (degree-1)/2"
!     STOP
!  END IF
  
  IF (LPB%boundary_condition_x.EQ.BC_PERIODIC) THEN
     ! map back outlying points
     IF (xpos_in<LPB%nodes_x(1)) THEN
        nfac=INT( (LPB%nodes_x(1)-xpos_in)/(LPB%right_boundary_x-LPB%left_boundary_x) )
        xpos = xpos_in + (nfac+1)*(LPB%right_boundary_x-LPB%left_boundary_x)
     ELSEIF (xpos_in>=LPB%right_boundary_x) THEN
        nfac=INT( (xpos_in-LPB%right_boundary_x)/(LPB%right_boundary_x-LPB%left_boundary_x) )
        xpos = xpos_in - (nfac+1)*(LPB%right_boundary_x-LPB%left_boundary_x)
     ELSE
        xpos = xpos_in
     END IF
  ELSEIF (LPB%boundary_condition_x.EQ.BC_NONPERIODIC) THEN
     IF ((xpos_in.le.(LPB%left_boundary_x-(LPB%nodes_x(2)-LPB%nodes_x(1)))).OR.&
          (xpos_in.ge.(LPB%right_boundary_x+(LPB%right_boundary_x-LPB%nodes_x(LPB%nCoeffs_x-1))))) THEN
        ! boundary zero boundary condition
       res=0.
       return
     ELSE 
        xpos = xpos_in
     END IF
  END IF

  IF (LPB%boundary_condition_y.EQ.BC_PERIODIC) THEN
     ! map back outlying points
     IF (ypos_in<LPB%nodes_y(1)) THEN
        nfac=INT( (LPB%nodes_y(1)-ypos_in)/(LPB%right_boundary_y-LPB%left_boundary_y) )
        ypos = ypos_in + (nfac+1)*(LPB%right_boundary_y-LPB%left_boundary_y)
     ELSEIF (ypos_in>=LPB%right_boundary_y) THEN
        nfac=INT( (ypos_in-LPB%right_boundary_y)/(LPB%right_boundary_y-LPB%left_boundary_y) )
        ypos = ypos_in - (nfac+1)*(LPB%right_boundary_y-LPB%left_boundary_y)
     ELSE
        ypos = ypos_in
     END IF
  ELSEIF (LPB%boundary_condition_y.EQ.BC_NONPERIODIC) THEN
     IF ((ypos_in.le.(LPB%left_boundary_y-(LPB%nodes_y(2)-LPB%nodes_y(1)))).OR.&
          (ypos_in.ge.(LPB%right_boundary_y+(LPB%right_boundary_y-LPB%nodes_y(LPB%nCoeffs_y-1))))) THEN
        ! boundary zero boundary condition
       res=0.
       return
     ELSE
        ypos = ypos_in
     END IF
  END IF

  !x value of left, central and right points
  if (ind_x.gt.1) then
     lb_x=LPB%nodes_x(ind_x-1)
  elseif ((ind_x.eq.1).AND.(LPB%boundary_condition_x.eq.BC_NONPERIODIC)) then
     lb_x=LPB%nodes_x(ind_x)-(LPB%nodes_x(ind_x+1)-LPB%nodes_x(ind_x))
  end if
  cp_x=LPB%nodes_x(ind_x)
  if (ind_x.lt.LPB%nCoeffs_x) then
     rb_x=LPB%nodes_x(ind_x+1)
  elseif ((ind_x.eq.LPB%nCoeffs_x).AND.(LPB%boundary_condition_x.eq.BC_NONPERIODIC)) then
     rb_x=LPB%nodes_x(ind_x) + (LPB%nodes_x(ind_x)-LPB%nodes_x(ind_x-1))
  ELSEIF ((ind_x.EQ.LPB%nCoeffs_x).AND.(LPB%boundary_condition_x.EQ.BC_PERIODIC)) THEN
     rb_x = LPB%right_boundary_x
  end if

  if (ind_y.gt.1) then
     lb_y=LPB%nodes_y(ind_y-1)
  elseif ((ind_y.eq.1).AND.(LPB%boundary_condition_y.eq.BC_NONPERIODIC)) then
     lb_y=LPB%nodes_y(ind_y)-(LPB%nodes_y(ind_y+1)-LPB%nodes_y(ind_y))
  end if
  cp_y=LPB%nodes_y(ind_y)
  if (ind_y.lt.LPB%nCoeffs_y) then
     rb_y=LPB%nodes_y(ind_y+1)
  elseif ((ind_y.eq.LPB%nCoeffs_y).AND.(LPB%boundary_condition_y.eq.BC_NONPERIODIC)) then
     rb_y=LPB%nodes_y(ind_y) + (LPB%nodes_y(ind_y)-LPB%nodes_y(ind_y-1))
  ELSEIF ((ind_y.EQ.LPB%nCoeffs_y).AND.(LPB%boundary_condition_y.EQ.BC_PERIODIC)) THEN
     rb_y = LPB%right_boundary_y
  end if

  if((ind_x.gt.1).and.(lb_x.le.xpos).and.(xpos.lt.cp_x)) THEN 

     if((ind_y.gt.1).and.(lb_y.le.ypos).and.(ypos.lt.cp_y)) THEN

     !left of the central x point, lb<=x<cp
     !left of the central y point, lb<=y<cp
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,lb_x,cp_x,xpos,&
                               &lb_y,cp_y,ypos)
     END DO
     
  elseif((ind_y.eq.1).and.(LPB%boundary_condition_y.eq.BC_PERIODIC).and.&
       &(ypos.ge.LPB%nodes_y(LPB%nCoeffs_y))) then
     
     !left of the central x point, lb<=x<cp
     !exploit periodicity for the contributions left of the first y point
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,lb_x,cp_x,xpos,&
             &LPB%right_boundary_y,LPB%nodes_y(LPB%nCoeffs_y),ypos)
     END DO

     elseif ((ind_y.le.LPB%nCoeffs_y).and.(cp_y.le.ypos).and.(ypos.lt.rb_y)) then
     !left of the central  x point, lb<=x<cp
     !right of the central y point, cp<=y<rb
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,lb_x,cp_x,xpos,&
             &cp_y,rb_y,ypos)
     END DO

     ELSE
        res=0.0D0
     END IF

  elseif((ind_x.eq.1).and.(LPB%boundary_condition_x.eq.BC_PERIODIC).and.(xpos.ge.LPB%nodes_x(LPB%nCoeffs_x))) then
     !exploit periodicity for the contributions left of the first point

     if((ind_y.gt.1).and.(lb_y.le.ypos).and.(ypos.lt.cp_y)) THEN

     !exploit periodicity for the contributions left of the first x point
     !left of the central y point, lb<=y<cp
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,LPB%right_boundary_x,LPB%nodes_x(LPB%nCoeffs_x),xpos,&
                               &lb_y,cp_y,ypos)
     END DO

     elseif((ind_y.eq.1).and.(LPB%boundary_condition_y.eq.BC_PERIODIC).and.&
                                                &(ypos.ge.LPB%nodes_y(LPB%nCoeffs_y))) then

     !exploit periodicity for the contributions left of the first x point
     !exploit periodicity for the contributions left of the first y point                                    
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,LPB%right_boundary_x,LPB%nodes_x(LPB%nCoeffs_x),xpos,&
                               & LPB%right_boundary_y,LPB%nodes_y(LPB%nCoeffs_y),ypos)
     END DO

     elseif ((ind_y.le.LPB%nCoeffs_y).and.(cp_y.le.ypos).and.(ypos.lt.rb_y)) then
     !exploit periodicity for the contributions left of the first x point
     !right of the central y point, cp<=y<rb
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,LPB%right_boundary_x,LPB%nodes_x(LPB%nCoeffs_x),xpos,&
                               & cp_y,rb_y,ypos)
     END DO

     ELSE
        res = 0.0D0 
     END IF

  elseif((ind_x.eq.1).and.(LPB%boundary_condition_x.eq.BC_NONPERIODIC).and.(lb_x.le.xpos).and.(xpos.lt.cp_x)) then
     ! Outside left boundary, here one could apply a boundary condition

          if((ind_y.gt.1).and.(lb_y.le.ypos).and.(ypos.lt.cp_y)) THEN

     !Outside left x boundary
     !left of the central y point, lb<=y<cp
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp_x,lb_x,xpos,&
                               &lb_y,cp_y,ypos)
     END DO

     elseif((ind_y.eq.1).and.(LPB%boundary_condition_y.eq.BC_PERIODIC).and.&
                                                &(ypos.ge.LPB%nodes_y(LPB%nCoeffs_y))) then

     !Outside left x boundary
     !exploit periodicity for the contributions left of the first y point                                            
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp_x,lb_x,xpos,&
                               &LPB%right_boundary_y,LPB%nodes_y(LPB%nCoeffs_y),ypos)
     END DO

     elseif ((ind_y.le.LPB%nCoeffs_y).and.(cp_y.le.ypos).and.(ypos.lt.rb_y)) then
     !Outside left x boundary
     !right of the central y point, cp<=y<rb
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp_x,lb_x,xpos,&
                               &cp_y,rb_y,ypos)
     END DO
     ELSE
        res = 0.0D0
     END IF

  elseif ((ind_x.le.LPB%nCoeffs_x).and.(cp_x.le.xpos).and.(xpos.lt.rb_x)) then 
     !right of the central point, cp<=x<rb

     if((ind_y.gt.1).and.(lb_y.le.ypos).and.(ypos.lt.cp_y)) THEN

     !right of the central oint, cp<=x<rb
     !left of the central y point, lb<=y<cp
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp_x,rb_x,xpos,&
                               &lb_y,cp_y,ypos)
     END DO

     elseif((ind_y.eq.1).and.(LPB%boundary_condition_y.eq.BC_PERIODIC).and.&
                                                &(ypos.ge.LPB%nodes_y(LPB%nCoeffs_y))) then

     !right of the central x point, cp<=x<rb
     !exploit periodicity for the contributions left of the first y point
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp_x,rb_x,xpos,&
                               &LPB%right_boundary_y,LPB%nodes_y(LPB%nCoeffs_y),ypos)
     END DO

     elseif ((ind_y.le.LPB%nCoeffs_y).and.(cp_y.le.ypos).and.(ypos.lt.rb_y)) then
     !right of the central x point, cp<=x<rb
     !right of the central y point, cp<=y<rb

     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp_x,rb_x,xpos,&
                               &cp_y,rb_y,ypos)
     END DO

     ELSE
        res = 0.0D0
     ENDIF

  elseif ((ind_x.eq.LPB%nCoeffs_x).and.(LPB%boundary_condition_x.eq.BC_NONPERIODIC)&
       &.and.(cp_x.le.xpos).and.(xpos.lt.rb_x)) then 
     !outside right boundary, here one could apply a boundary condition
 
     if((ind_y.gt.1).and.(lb_y.le.ypos).and.(ypos.lt.cp_y)) THEN

     !Outside right x boundary
     !left of the central y point, lb<=y<cp
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp_x,rb_x,xpos,&
                               &lb_y,cp_y,ypos)
     END DO

     elseif((ind_y.eq.1).and.(LPB%boundary_condition_y.eq.BC_PERIODIC).and.&
                                                &(ypos.ge.LPB%nodes_y(LPB%nCoeffs_y))) then

     !Outside right x boundary
     !exploit periodicity for the contributions left of the first y point
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp_x,rb_x,xpos,&
                               &LPB%right_boundary_y,LPB%nodes_y(LPB%nCoeffs_y),ypos)
     END DO

     elseif ((ind_y.le.LPB%nCoeffs_y).and.(cp_y.le.ypos).and.(ypos.lt.rb_y)) then
     !Outside right x boundary
     !right of the central y point, cp<=y<rb
     DO ideriv=0,LPB%nDeriv
        res(ideriv) = base_poly(ideriv,LPB%degree,cp_x,rb_x,xpos,&
                               &cp_y,rb_y,ypos)
     END DO

     ELSE
        res = 0.0D0 
     END IF
  else
     ! the base function has only a support on a 2*dx interval
     res = 0.0D0
  endif
  
END SUBROUTINE mp_get_value2d

!****im* localpolynombase_mod/base_poly
! SYNOPSIS
FUNCTION base_poly(ideriv,degree,x0,x1,x,y0,y1,y)
  INTEGER :: ideriv, degree
  REAL :: x1, x0,x, y0,y1,y
  real:: xp, x0p,yp,y0p
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
  
  yp=y-y1
  y0p=y0-y1

!  SELECT CASE (degree)
!     CASE (1)
!        IF (ideriv.GT.0) THEN 
!           PRINT*,"Degree 1 polynoms don't use derivatives."
!           STOP
!        END IF

  base_poly =  (xp)*(yp) / (x0p*y0p)
       
!     CASE (3)
!        SELECT CASE (ideriv)
!        CASE (0) ! the function value
!           base_poly = xp**2*(3*x0p-2*xp)/x0p**3
!        CASE (1) ! first derivative
!           base_poly = (xp - x0p)*xp**2/x0p**2
!        CASE default
!           PRINT*,"Degree 3 polynoms only use first derivatives."
!           STOP
!        END SELECT

!     CASE (5)
!        SELECT CASE (ideriv)
!        CASE (0) ! the function value
!           base_poly = xp**3*(6*xp**2 - 15*xp*x0p + 10*x0p**2)/x0p**5
!        CASE (1)
!           base_poly = -((xp - x0p)*xp**3*(3*xp - 4*x0p))/x0p**4
!        CASE (2)
!           base_poly = ((xp - x0p)**2*xp**3)/(2.*x0p**3)
!        CASE default
!           PRINT*,"Degree 5 polynoms only uses first two derivatives."
!           STOP
!        END SELECT
!     CASE DEFAULT
!       PRINT*,"For present state of GENE3D, only Degree 1 polynoms implemented"

!  END SELECT
END FUNCTION base_poly

SUBROUTINE finalize_localpolynombase(LPB)
  type(LocalPolynomBase2d) :: LPB

  if (LPB%nodes_allocated) THEN
     DEALLOCATE(LPB%nodes_x)
     DEALLOCATE(LPB%nodes_y)
     LPB%nodes_allocated=.FALSE.
  END IF

END SUBROUTINE finalize_localpolynombase

END MODULE localpolynombase2d_mod
