#include "redef.h"
!>Contains routines for radial derivatives in the global version
MODULE x_derivatives
  Use boundaries !, only: exchange_x, lx0_boundary, lx0_boundary_2D
  USE discretization
  USE par_other, ONLY : print_ini_msg, p_has_0_mode

  IMPLICIT NONE
  PUBLIC:: initialize_radscheme, x_deriv, x_deriv_exc, finalize_radscheme, &
       CENTERED4TH, calc_gradient, calc_full_x_deriv
  ! subroutines for adaptive grids
  public :: x_deriv_1d_adptv, x_deriv_exc_2d_adptv, radscheme
  public :: max_ddi
  PRIVATE

  INTEGER,PARAMETER :: CENTERED2ND=2,CENTERED4TH=4,CENTERED6TH=6
  INTEGER :: radscheme=CENTERED4TH, rad_bc_type
  REAL, DIMENSION(:),ALLOCATABLE :: coeff
  REAL:: dx, max_ddi

  Interface x_deriv
     Module Procedure x_deriv_1d, x_deriv_2d, x_deriv_2d_nonlin, x_deriv_1d_real
  End Interface

  Interface x_deriv_exc
     Module Procedure x_deriv_exc_2d,x_deriv_exc_2d_nonlin
  End Interface


CONTAINS
  !>Initialize finite differences in radial direction (global code only)
  !!Allocates and sets stencils
  SUBROUTINE initialize_radscheme(a_radscheme,a_dx, a_rad_bc_type)
    integer, intent(in) :: a_radscheme,a_rad_bc_type
    real, intent(in) :: a_dx

    rad_bc_type = a_rad_bc_type
    dx = a_dx

    if ((li1-lbi).lt.a_radscheme/2) STOP &
         'not enough boundary points for chosen radial derivative scheme'

    ALLOCATE(coeff(-a_radscheme/2:a_radscheme/2))

    SELECT CASE (a_radscheme)
    CASE (CENTERED2ND)
       if ((mype==0).and.(print_ini_msg)) &
            &WRITE(*,"(A)") "Radial derivative scheme: 2nd order centered differences"
       coeff = (/-1.0,0.0,1.0/)/(2.0*dx)
       !kimax_fac = 0.32  ! (=1/pi) max. fraction of kimax
       max_ddi = 1.0/dx   ! (Sin(pi/2)/dx) max. eigenvalue at a scale of 0.5 kimax
    CASE (CENTERED4TH)
       if ((mype==0).and.(print_ini_msg)) &
            &WRITE(*,"(A)") "Radial derivative scheme: 4th order centered differences"
       coeff = (/1.0,-8.0,0.0,8.0,-1.0/)/&
            (12.*dx)
       !kimax=pi/di
       !kimax_fac = 0.445  ! (0.43679,more precisely) max. fraction of kimax
       max_ddi = 1.4/dx    ! (1.372/dx,more precisely) max. eigenvalue at a scale of 0.61 kimax
    CASE (CENTERED6TH)
       if ((mype==0).and.(print_ini_msg)) &
            &WRITE(*,"(A)") "Radial derivative scheme: 6th order centered differences"
       coeff = (/-1.0,+9.0,-45.0,0.0,45.0,&
            -9.0,1.0/)/(60.0*dx)
       !kimax_fac = 0.51  ! max. fraction of kimax
       max_ddi = 1.59/dx  ! max. eigenvalue at a scale of 0.61 kimax
    CASE default
       stop 'no valid radial derivative scheme defined'
    END SELECT
    radscheme=a_radscheme
  END SUBROUTINE initialize_radscheme

  !> radial derivative 1-dimensional
  SUBROUTINE x_deriv_1d(y,dydx)
    COMPLEX, DIMENSION(lbi:ubi), INTENT(in) :: y
    COMPLEX, DIMENSION(li1:li2), INTENT(out) :: dydx
    INTEGER :: sten
    
    PERFON('xder1')

    !avoid initializing dydx with zero
    sten = -radscheme/2
    dydx(:) = coeff(sten)*y(li1+sten:li2+sten)
    
    DO sten=-radscheme/2+1,radscheme/2
       dydx(:) = dydx(:)+coeff(sten)*y(li1+sten:li2+sten)
    END DO

    PERFOFF
  
  END SUBROUTINE x_deriv_1d

  !> radial derivative 1-dimensional (for real types)
  SUBROUTINE x_deriv_1d_real(y,dydx)
    REAL, DIMENSION(lbi:ubi), INTENT(in) :: y
    REAL, DIMENSION(li1:li2), INTENT(out) :: dydx
    INTEGER :: sten
    
    PERFON('xder1r')

    !avoid initializing dydx with zero
    sten = -radscheme/2
    dydx(:) = coeff(sten)*y(li1+sten:li2+sten)
    
    DO sten=-radscheme/2+1,radscheme/2
       dydx(:) = dydx(:)+coeff(sten)*y(li1+sten:li2+sten)
    END DO

    PERFOFF
  
  END SUBROUTINE x_deriv_1d_real


  !> radial derivative 2-dimensional with boundary exchange
  !! \param y input array (has to be inout since we modify the boundaries)
  !! \param dydx output array containing the derivatives
  SUBROUTINE x_deriv_exc_2d(y,dydx)
    COMPLEX, DIMENSION(lbi:ubi, lj1:lj2), INTENT(inout):: y
    COMPLEX, DIMENSION(li1:li2, lj1:lj2), INTENT(out) :: dydx
    INTEGER :: i
    
    PERFON_I('xder2exc')

    !print*,"START  x_deriv_exc_2d"
    !call exchange routine (either Dirichlet or 
    CALL exchange_x(lx0_boundary,y)

    If ((rad_bc_type.eq.2).and.(p_has_0_mode).and.(my_pex.eq.0)) then
       Do i=1,nib
          y(li1-i,lj1) = y(li1+i,lj1)
       Enddo
    Endif

    Call x_deriv_2d(y,dydx)
    !print*,"END    x_deriv_exc_2d"

    PERFOFF_I
  
  END SUBROUTINE x_deriv_exc_2d

  !> radial derivative 2-dimensional
  SUBROUTINE x_deriv_2d(y,dydx)
    COMPLEX, DIMENSION(lbi:ubi, lj1:lj2), INTENT(in) :: y
    COMPLEX, DIMENSION(li1:li2, lj1:lj2), INTENT(out):: dydx
    INTEGER :: j, sten
    
    !PERFON('x_deriv2')
#define USE_OLD_VERSION
#ifdef USE_OLD_VERSION
    !avoid initializing dydx with zero
    Do j=lj1, lj2
       sten = -radscheme/2
       dydx(:,j) = coeff(sten)*y(li1+sten:li2+sten,j)
       DO sten=-radscheme/2+1,radscheme/2
          dydx(:,j) = dydx(:,j)+coeff(sten)*y(li1+sten:li2+sten,j)
       END DO
    Enddo
#else
    INTEGER :: i
    !avoid initializing dydx with zero
    Do j=lj1, lj2
       do i=li1,li2
          dydx(i,j) = dot_product(coeff(-radscheme/2:radscheme/2),y(i-radscheme/2:i+radscheme/2,j))
       end do
    Enddo
#endif
    !PERFOFF
  
  END SUBROUTINE x_deriv_2d


  !> radial derivative 2-dimensional for Arakawa nonlinearity
  SUBROUTINE x_deriv_2d_nonlin(y,dydx,li0da)
    INTEGER, INTENT(in):: li0da
    REAL, DIMENSION(-nib:li0da-1+nib,0:ly0da-1), INTENT(in):: y
    REAL, DIMENSION(0:li0da-1,0:ly0da-1), INTENT(out)::       dydx
    INTEGER :: j, sten
    
    PERFON('nlxder2')

    !avoid initializing dydx with zero
    Do j=0, ly0da-1
       sten = -radscheme/2
       dydx(:,j) = coeff(sten)*y(sten:li0da-1+sten,j)
       DO sten=-radscheme/2+1,radscheme/2
          dydx(:,j) = dydx(:,j)+coeff(sten)*y(sten:li0da-1+sten,j)
       END DO
    Enddo

    PERFOFF
  
  END SUBROUTINE x_deriv_2d_nonlin


  !> radial derivative 2-dimensional for Arakawa nonlinearity with boundary exchange
  !! \param y_in input array
  !! \param dydx The derivative of y_in in x direction
  !! \param li0da The number of points in the first dimension of dydx, li0 with dealiasing boundaries
  SUBROUTINE x_deriv_exc_2d_nonlin(y_in,dydx,li0da)
    INTEGER,INTENT(in):: li0da
    REAL, DIMENSION(-nib:li0da-1+nib,0:ly0da-1), INTENT(inout):: y_in
    REAL, DIMENSION(0:li0da-1,0:ly0da-1), INTENT(out):: dydx
    REAL, DIMENSION(-nib:li0da-1+nib,0:ly0da-1) :: y
    INTEGER :: i, j
    real :: y_avg
    
    PERFON('nlxd2exc')

    call exchange_x(lx0_boundary_2D_real,y_in)
    y = REAL(y_in)    

    If ((rad_bc_type.eq.2).and.(my_pex.eq.0)) then
       Do i=1,nib
          !apply average values since we are in real space
          y_avg = SUM(y(0+i,:))/(ly0da)
          Do j=0,ly0da-1
             y(0-i,j) = y_avg !y(0+i,j)
          End Do
       Enddo
    Endif
   
    CALL x_deriv_2d_nonlin(y,dydx,li0da)
    
    PERFOFF
  
  END SUBROUTINE x_deriv_exc_2d_nonlin


  !>Computes the (full, i.e. 0:nx0-1) gradient profile for a given 
  !!temperature or density profile
  Subroutine calc_gradient(inarr,outarr)
    Real, Dimension(0:nx0-1), intent(in)  :: inarr  !>input profile
    Real, Dimension(0:nx0-1), intent(out) :: outarr !>output gradient profile
    Real, Dimension(0:nx0-1) :: tmp_x, tmp_dx
    integer :: i, sten, nib

    tmp_x = LOG(inarr)
    nib = radscheme/2

    !Apply centered scheme
    !avoid initializing tmp_dx with zero
    sten = -radscheme/2
    tmp_dx(nib:nx0-1-nib) = coeff(sten)*tmp_x(nib+sten:nx0-1-nib+sten)
    
    DO sten=-radscheme/2+1,radscheme/2
       tmp_dx(nib:nx0-1-nib) = tmp_dx(nib:nx0-1-nib)+&
            &coeff(sten)*tmp_x(nib+sten:nx0-1-nib+sten)
    END DO

    !Close to the boundaries, we overwrite the centered scheme result by
    !2nd order upwind/downwind schemes
    do i=0,nib-1
       tmp_dx(i)=(-3.*tmp_x(i)+4.*tmp_x(i+1)-tmp_x(i+2))*0.5/dx
       tmp_dx(nx0-1-i)=(3.*tmp_x(nx0-1-i)-4.*tmp_x(nx0-2-i)+tmp_x(nx0-3-i))*0.5/dx
    enddo

    outarr = -real(tmp_dx) !/(rhostar*minor_r) !x_deriv is normalized to rho_ref
  end Subroutine calc_gradient


  !>Computes the (full, i.e. 0:nx0-1) gradient profile for a given 
  !!temperature or density profile
  Subroutine calc_full_x_deriv(inarr,outarr)
    Real, Dimension(0:nx0-1), intent(in)  :: inarr  !>input profile
    Real, Dimension(0:nx0-1), intent(out) :: outarr !>output gradient profile
    Real, Dimension(0:nx0-1) :: tmp_x, tmp_dx
    integer :: i, sten, nib

    tmp_x = inarr
    nib = radscheme/2

    !Apply centered scheme
    !avoid initializing tmp_dx with zero
    sten = -radscheme/2
    tmp_dx(nib:nx0-1-nib) = coeff(sten)*tmp_x(nib+sten:nx0-1-nib+sten)
    
    DO sten=-radscheme/2+1,radscheme/2
       tmp_dx(nib:nx0-1-nib) = tmp_dx(nib:nx0-1-nib)+&
            &coeff(sten)*tmp_x(nib+sten:nx0-1-nib+sten)
    END DO

    !Close to the boundaries, we overwrite the centered scheme result by
    !2nd order upwind/downwind schemes
    do i=0,nib-1
       tmp_dx(i)=(-3.*tmp_x(i)+4.*tmp_x(i+1)-tmp_x(i+2))*0.5/dx
       tmp_dx(nx0-1-i)=(3.*tmp_x(nx0-1-i)-4.*tmp_x(nx0-2-i)+tmp_x(nx0-3-i))*0.5/dx
    enddo

    outarr = real(tmp_dx) !/(rhostar*minor_r) !x_deriv is normalized to rho_ref
  end Subroutine calc_full_x_deriv


  SUBROUTINE finalize_radscheme

    DEALLOCATE(coeff)

  END SUBROUTINE finalize_radscheme

  !> radial derivative 1-dimensional adaptive mirror
  SUBROUTINE x_deriv_1d_adptv(y,dydx,i1,i2)
    COMPLEX, DIMENSION(lbi:ubi), INTENT(in) :: y
    COMPLEX, DIMENSION(li1:li2), INTENT(out) :: dydx
    integer, intent(in) :: i1, i2
    INTEGER :: sten
    
    PERFON('xder1')

    !avoid initializing dydx with zero
    sten = -radscheme/2
    dydx(i1:i2) = coeff(sten)*y(i1+sten:i2+sten)
    
    DO sten=-radscheme/2+1,radscheme/2
       dydx(i1:i2) = dydx(i1:i2)+coeff(sten)*y(i1+sten:i2+sten)
    END DO

    PERFOFF
  
  END SUBROUTINE x_deriv_1d_adptv
  
  !> radial derivative 2-dimensional with boundary exchange
  !! \param y input array (has to be inout since we modify the boundaries)
  !! \param dydx output array containing the derivatives
  SUBROUTINE x_deriv_exc_2d_adptv(y,dydx,i1,i2)
    COMPLEX, DIMENSION(lbi:ubi, lj1:lj2), INTENT(inout):: y
    COMPLEX, DIMENSION(li1:li2, lj1:lj2), INTENT(out) :: dydx
    integer, intent(in) :: i1, i2
    INTEGER :: i
    
    PERFON_I('xder2exc')

    !call exchange routine (either Dirichlet or ??)
    CALL exchange_x(lx0_boundary,y)

    If ((rad_bc_type.eq.2).and.(p_has_0_mode).and.(my_pex.eq.0)) then
       Do i=1,nib
          y(i1-i,lj1) = y(i1+i,lj1)
       Enddo
    Endif

    Call x_deriv_2d_adptv(y,dydx,i1,i2)

    PERFOFF_I
  
  END SUBROUTINE x_deriv_exc_2d_adptv

  !> radial derivative 2-dimensional
  SUBROUTINE x_deriv_2d_adptv(y,dydx,i1,i2)
    COMPLEX, DIMENSION(lbi:ubi, lj1:lj2), INTENT(in) :: y
    COMPLEX, DIMENSION(li1:li2, lj1:lj2), INTENT(out):: dydx
    integer, intent(in) :: i1, i2
    INTEGER :: j
    
    !PERFON('x_deriv2')
#define USE_OLD_VERSION
#ifdef USE_OLD_VERSION
    INTEGER :: sten

    !avoid initializing dydx with zero
    Do j=lj1, lj2
       sten = -radscheme/2
       dydx(i1:i2,j) = coeff(sten)*y(i1+sten:i2+sten,j)
       DO sten=-radscheme/2+1,radscheme/2
          dydx(i1:i2,j) = dydx(i1:i2,j)+coeff(sten)*y(i1+sten:i2+sten,j)
       END DO
    Enddo
#else
    INTEGER :: i

    !avoid initializing dydx with zero
    Do j=lj1, lj2
       do i=i1,i2
          dydx(i,j) = dot_product(coeff(-radscheme/2:radscheme/2),y(i-radscheme/2:i+radscheme/2,j))
       end do
    Enddo
#endif
    !PERFOFF
  
  END SUBROUTINE x_deriv_2d_adptv

END MODULE x_derivatives
