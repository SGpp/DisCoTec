#include "redef.h"
!>lagrange_interpolation contains subroutines to perform
!!a mid-point lagrange interpolation of order 3
MODULE lagrange_interpolation
  IMPLICIT NONE
  PUBLIC:: lag3interp, lag3deriv, lag3interp_2d
  PUBLIC:: lag3interp_complex
  PRIVATE

  INTERFACE lag3interp
     MODULE PROCEDURE lag3interp_scalar, lag3interp_array
  END INTERFACE

  INTERFACE lag3deriv
     MODULE PROCEDURE lag3deriv_scalar, lag3deriv_array
  END INTERFACE

CONTAINS

  !> Third order lagrange interpolation
  SUBROUTINE lag3interp_scalar(y_in,x_in,n_in,y_out,x_out)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n_in
    REAL, DIMENSION(n_in), INTENT(IN) :: y_in,x_in
    REAL, INTENT(IN) :: x_out
    REAL, INTENT(OUT) :: y_out

    REAL, DIMENSION(1) :: xout_wrap, yout_wrap

    xout_wrap = x_out
    call lag3interp_array(y_in,x_in,n_in,yout_wrap,xout_wrap,1)
    y_out = yout_wrap(1)

  END SUBROUTINE lag3interp_scalar

  !> Third order lagrange interpolation
  subroutine lag3interp_array(y_in,x_in,n_in,y_out,x_out,n_out)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n_in,n_out
    REAL, DIMENSION(n_in), INTENT(IN) :: y_in,x_in
    REAL, DIMENSION(n_out), INTENT(IN) :: x_out
    REAL, DIMENSION(n_out), INTENT(OUT) :: y_out

    REAL :: x,aintm,aint0,aint1,aint2,xm,x0,x1,x2
    INTEGER :: j,jm,j0,j1,j2
    INTEGER :: jstart,jfirst,jlast,jstep

    IF (x_in(n_in) > x_in(1)) THEN
       jstart=3
       jfirst=1
       jlast=n_out
       jstep=1
    ELSE
       jstart=n_in-2
       jfirst=n_out
       jlast=1
       jstep=-1
    END IF

    j1=jstart
    DO j=jfirst,jlast,jstep
       x=x_out(j)
       DO WHILE (x >= x_in(j1) .AND. j1 < n_in-1 .AND. j1 > 2) 
          j1=j1+jstep
       END DO

       j2=j1+jstep
       j0=j1-jstep
       jm=j1-2*jstep

       !...  extrapolate inside or outside

       x2=x_in(j2)
       x1=x_in(j1)
       x0=x_in(j0)
       xm=x_in(jm)

       aintm=(x-x0)*(x-x1)*(x-x2)/((xm-x0)*(xm-x1)*(xm-x2))
       aint0=(x-xm)*(x-x1)*(x-x2)/((x0-xm)*(x0-x1)*(x0-x2))
       aint1=(x-xm)*(x-x0)*(x-x2)/((x1-xm)*(x1-x0)*(x1-x2))
       aint2=(x-xm)*(x-x0)*(x-x1)/((x2-xm)*(x2-x0)*(x2-x1))

       y_out(j)=aintm*y_in(jm)+aint0*y_in(j0) &
            +aint1*y_in(j1)+aint2*y_in(j2)

    END DO

  END SUBROUTINE Lag3interp_array
  
  
  !> Third order lagrange interpolation for complex arrays
  SUBROUTINE lag3interp_complex(y_in,x_in,n_in,y_out,x_out,n_out)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n_in,n_out
    COMPLEX, DIMENSION(n_in), INTENT(IN) :: y_in
    REAL, DIMENSION(n_in), INTENT(IN) :: x_in
    COMPLEX, DIMENSION(n_out), INTENT(OUT) :: y_out
    REAL, DIMENSION(n_out), INTENT(IN) :: x_out

    REAL :: x,aintm,aint0,aint1,aint2,xm,x0,x1,x2
    INTEGER :: j,jm,j0,j1,j2
    INTEGER :: jstart,jfirst,jlast,jstep

    IF (x_in(n_in) > x_in(1)) THEN
       jstart=3
       jfirst=1
       jlast=n_out
       jstep=1
    ELSE
       jstart=n_in-2
       jfirst=n_out
       jlast=1
       jstep=-1
    END IF

    j1=jstart
    DO j=jfirst,jlast,jstep
       x=x_out(j)
       DO WHILE (x >= x_in(j1) .AND. j1 < n_in-1 .AND. j1 > 2) 
          j1=j1+jstep
       END DO

       j2=j1+jstep
       j0=j1-jstep
       jm=j1-2*jstep

       !...  extrapolate inside or outside

       x2=x_in(j2)
       x1=x_in(j1)
       x0=x_in(j0)
       xm=x_in(jm)

       aintm=(x-x0)*(x-x1)*(x-x2)/((xm-x0)*(xm-x1)*(xm-x2))
       aint0=(x-xm)*(x-x1)*(x-x2)/((x0-xm)*(x0-x1)*(x0-x2))
       aint1=(x-xm)*(x-x0)*(x-x2)/((x1-xm)*(x1-x0)*(x1-x2))
       aint2=(x-xm)*(x-x0)*(x-x1)/((x2-xm)*(x2-x0)*(x2-x1))

       y_out(j)=aintm*y_in(jm)+aint0*y_in(j0) &
            +aint1*y_in(j1)+aint2*y_in(j2)

    END DO

  END SUBROUTINE lag3interp_complex


  !>2D interpolation 
  !\TODO check whether a "real" 2D interpolation would
  !! be more appropriate
  SUBROUTINE lag3interp_2d(y_in,x1_in,n1_in,x2_in,n2_in,&
       &y_out,x1_out,n1_out,x2_out,n2_out)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: n1_in,n2_in,n1_out,n2_out
    REAL, DIMENSION(n1_in,n2_in), INTENT(IN) :: y_in
    REAL, DIMENSION(n1_in) :: x1_in
    REAL, DIMENSION(n2_in) :: x2_in
    REAL, DIMENSION(n1_out), INTENT(IN) :: x1_out
    REAL, DIMENSION(n2_out), INTENT(IN) :: x2_out
    REAL, DIMENSION(n1_out,n2_out), INTENT(OUT) :: y_out    

    !local variables
    REAL, DIMENSION(n2_in) :: y2_in_tmp
    REAL, DIMENSION(n2_out) :: y2_out_tmp
    REAL, DIMENSION(n1_in,n2_out) :: y_tmp
    INTEGER :: i

    DO i=1,n1_in
       y2_in_tmp = y_in(i,:)
       call lag3interp(y2_in_tmp,x2_in,n2_in,&
            y2_out_tmp,x2_out,n2_out)
       y_tmp(i,:) = y2_out_tmp
    ENDDO
    
    DO i=1,n2_out
       call lag3interp(y_tmp(:,i),x1_in,n1_in,&
            y_out(:,i),x1_out,n1_out)
    END DO
    
  END SUBROUTINE lag3interp_2d

  !> Third order lagrange interpolation
  SUBROUTINE lag3deriv_scalar(y_in,x_in,n_in,dydx_out,x_out)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n_in
    REAL, DIMENSION(n_in), INTENT(IN) :: y_in,x_in
    REAL, INTENT(IN) :: x_out
    REAL, INTENT(OUT) :: dydx_out

    REAL, DIMENSION(1) :: xout_wrap, dydxout_wrap

    xout_wrap = x_out
    call lag3deriv_array(y_in,x_in,n_in,dydxout_wrap,xout_wrap,1)
    dydx_out = dydxout_wrap(1)

  END SUBROUTINE lag3deriv_scalar



!>Returns Derivative based on a 3rd order lagrange interpolation
  subroutine lag3deriv_array(y_in,x_in,n_in,dydx_out,x_out,n_out)

    IMPLICIT NONE

    INTEGER :: n_in,n_out
    REAL, DIMENSION(n_in), INTENT(IN) :: y_in,x_in
    REAL, DIMENSION(n_out), INTENT(IN) :: x_out
    REAL, DIMENSION(n_out), INTENT(OUT) :: dydx_out

    REAL :: x,aintm,aint0,aint1,aint2,xm,x0,x1,x2
    INTEGER :: j,jm,j0,j1,j2
    INTEGER :: jstart,jfirst,jlast,jstep

    IF (x_in(n_in) > x_in(1)) THEN
       jstart=3
       jfirst=1
       jlast=n_out
       jstep=1
    ELSE
       jstart=n_in-2
       jfirst=n_out
       jlast=1
       jstep=-1
    END IF
    
    j1=jstart
    DO j=jfirst,jlast,jstep
       x=x_out(j)
       DO WHILE (x >= x_in(j1) .AND. j1 < n_in-1 .AND. j1 > 2) 
          j1=j1+jstep
       END DO

       j2=j1+jstep
       j0=j1-jstep
       jm=j1-2*jstep

       !...  extrapolate inside or outside

       x2=x_in(j2)
       x1=x_in(j1)
       x0=x_in(j0)
       xm=x_in(jm)

       aintm=((x-x1)*(x-x2)+(x-x0)*(x-x2)+(x-x0)*(x-x1)) &
            /((xm-x0)*(xm-x1)*(xm-x2))
       aint0=((x-x1)*(x-x2)+(x-xm)*(x-x2)+(x-xm)*(x-x1)) &
            /((x0-xm)*(x0-x1)*(x0-x2))
       aint1=((x-x0)*(x-x2)+(x-xm)*(x-x2)+(x-xm)*(x-x0)) &
            /((x1-xm)*(x1-x0)*(x1-x2))
       aint2=((x-x0)*(x-x1)+(x-xm)*(x-x1)+(x-xm)*(x-x0)) &
            /((x2-xm)*(x2-x0)*(x2-x1))

       dydx_out(j)=aintm*y_in(jm)+aint0*y_in(j0) &
            +aint1*y_in(j1)+aint2*y_in(j2)

    END DO

  end subroutine Lag3deriv_array


end module lagrange_interpolation
