#include "redef.h"
!>Contains the values of the coordinates in all directions 
!!and the routines to set them
module coordinates
  use discretization
  use GaussQuadrature 
  use par_in, only: arakawa_zv
  !use Grid1DModule
  implicit none

  public :: lx, lx_a, deli, kx_center, kx, kxmin, xval, xval_a, x0
  PUBLIC :: lp1,lp2
  PUBLIC :: ly, delj, ky, kymin, kymax
  PUBLIC :: kjmin, kjmax, ki, kj
  PUBLIC :: dz, zval
  PUBLIC :: lv, dv, vp
  PUBLIC :: lw, mu
  PUBLIC :: vp_weight, mu_weight
  PUBLIC :: limit_shat, adapt_ly

  public :: initialize_coordinates, set_coordinates_defaults, set_x_coordinate_vars, &
       set_y_coordinate_vars, set_z_coordinate_vars, set_vp_coordinate_vars, &
       set_mu_coordinate_vars, finalize_coordinates, calc_center_value

  ! Box lengths
  real:: lx=0.0, lx_a=0.0, ly=0.0, lv, lw, lp1, lp2

  !adapt kymin/ly to minimum toroidal mode number n0_global (default=true for global)
  Logical :: adapt_ly = .false.

  ! limit shat in order to avoid huge lx values for very small shat
  Logical :: limit_shat=.false.

  ! min/max modes in Fourier space
  real :: kymin=0.0, kymax, kxmin
  real:: kjmin=0.0, kjmax

  ! box center and width (between 0 and 1)
  real :: x0=-1

  ! Grid spacing
  Real::    dx, dy, dz, dv, deli, delj

  ! mode shift for linear runs
  real:: kx_center=0.0D0

  ! coordinates
  Real, Dimension(:), Allocatable:: kx, ky, zval, xval, xval_a, ki, kj
  Real, Dimension(:), Allocatable:: vp, mu

  ! velocity space integration factors
  Real, Dimension(:), Allocatable :: mu_weight, vp_weight  

  !TYPE(Grid1D) :: xgrid
  !private module variables
  REAL, parameter   :: pi= 3.141592653589793239d0

  PRIVATE

contains

  subroutine initialize_coordinates

    allocate(kx(0:nx0-1)) 
    allocate(ky(0:nky0-1))
    
    allocate(ki(0:ni0-1),kj(0:nj0-1))

    if (.not.x_local) then
       allocate(xval(0:nx0-1))
    endif
    allocate(xval_a(pi1gl:pi2gl))
    !note: xval_a is a constant in the local code
    !      describing the position where e.g. profiles
    !      should be evaluated

    allocate(zval(0:nz0-1))

    allocate(vp(0:nv0-1),vp_weight(0:nv0-1)) !ll1:ll2))
    allocate(mu(0:nw0-1),mu_weight(0:nw0-1))  !lm1:lm2))

  end subroutine initialize_coordinates

  subroutine set_coordinates_defaults

    lx=0.; lx_a=0.; ly=0.; kx_center=0.; kymin=0; kjmin=0; x0=-1
    limit_shat=.false.;  adapt_ly = .false.
    call set_Gauss_defaults

  end subroutine set_coordinates_defaults

  subroutine set_x_coordinate_vars(nexc, n_pol, rad_bc_type, &
       &shat, rhostar, lilo, adapt_lx, mag_prof, &
       &only_zonal, write_pe,Cyq0_x0)
    INTEGER, intent(INOUT) :: nexc
    INTEGER, intent(IN) :: n_pol, rad_bc_type
    REAL, INTENT(INOUT) :: shat
    real, intent(IN) :: rhostar, Cyq0_x0
    LOGICAL, INTENT(IN) :: lilo, adapt_lx, mag_prof, only_zonal, write_pe

    INTEGER :: i, j
    !initialize radial coordinate

    if (x_local) then
       !In the local code, lx is adapted to fulfill 
       !the parallel b.c. constraint :
       !N = lx * kymin * shat = integer
       if ((abs(shat).gt.epsilon(shat)).and.(.not.only_zonal)) then
          IF (nx0.GT.1) &
               &call set_radial_box_size_with_shear(n_pol,shat,kymin,&
               &nexc,lx,adapt_lx,limit_shat,Cyq0_x0)
       else
          IF (ABS(lx).LT.EPSILON(lx)) THEN
             ! lx == 0
             IF (kx_center.GT.0.0) THEN
                ! kx_center is set and finite
                IF (only_zonal) lx = 2.0*pi/kx_center
             ELSE
                if (nx0.gt.1) STOP "lx.eq.0.0 and kx_center .eq. 0, can not determine x box size."
             END IF
          ELSE
             !lx .ne. 0, do nothing, kx_center is ignored
          END IF
       end if
    else !nonlocal
       if ((lilo.or..not.mag_prof).and.rad_bc_type==0) then
          IF ((ABS(shat).GT.1e-3).AND.(abs(kymin).gt.1e-5)) THEN
             call set_radial_box_size_with_shear(n_pol,shat,kymin,nexc,&
                  &lx,adapt_lx,limit_shat,Cyq0_x0)
          END IF
       END IF
    endif

    ! initialise further radial coordinate values
    IF ((ABS(lx).LT.EPSILON(lx)).AND.(write_pe)) PRINT*, "warning: lx is equal to zero"
    IF (x_local) THEN
       if (nx0.gt.1) kxmin=2.*pi/lx
    ELSE
       kxmin = 2.*pi/lx
    END IF
    if ((x_local.and.(evenx.eq.0)).or.((.not.x_local).and.(rad_bc_type.eq.0)))  then
       !local with odd nx0 (-kx_max:kx_max)
       !global periodic case: x=-lx/2,..,lx/2-dx
       dx = lx / nx0
    else
       !local with even nx0: Nyquist mode is discarded
       !global and nonperiodic: x=-lx/2,..,lx/2
       dx = lx / (nx0-1) !old version
       !new version with x = -lx/2,...,lx/2 -dx for combitechnique
       !toDo check if this is right!
       !dx = lx/ nx0
    endif
    if (.not.y_local) dx=lx/2.0/nx0

    IF (.not.x_local) THEN
       do i=0,nx0-1
          xval(i) = -0.5*lx + i*dx
       enddo
       xval = xval + x0/rhostar
       xval_a = xval*rhostar
    ELSE
       xval_a = x0
    ENDIF
    
    if (nx0.gt.1) then
       if (y_local) then
          do i=0,hkx+evenx
             kx(i) = kxmin*i+kx_center
          enddo
          do i=lkx,nx0-1
             kx(i)=kxmin*(i-nx0)+kx_center
          enddo
       else
          do j=lj1,lj2
             kx(j) = kxmin*j 
          enddo
       endif
    else
       kx(lg1)=kx_center
       if (abs(kx_center).gt.0) then
          lx=2.0*pi/kx_center
       endif
    endif

    if (yx_order) then
       lp1=ly
       lp2=lx
       kjmin=kxmin
       kjmax=maxval(kx)
       kj=kx
       delj=dx
    else
       lp1=lx
       lp2=ly
       deli=dx
       ki=kx
    end if

  end subroutine set_x_coordinate_vars


  subroutine set_y_coordinate_vars(n0_global, q0, &
       &C_y, rho_Lref, lilo, write_pe)
    implicit none
    
    INTEGER, intent(INOUT) :: n0_global
    real, intent(IN) :: rho_Lref, q0
    REAL, DIMENSION(pi1gl:pi2gl), INTENT(IN) :: C_y
    LOGICAL, INTENT(IN) :: lilo, write_pe

    INTEGER :: j
    REAL :: C_y0

    IF (rho_Lref.gt.0) THEN
       !n0_global will be now also be adapted in the local code
       !if rhostar is given
       call calc_center_value(C_y, C_y0)

       IF (n0_global.EQ.-1111) THEN 
          !find n0_global which gives a kymin as close as 
          !possible to the input kymin
          n0_global = nint((kymin*C_y0/rho_Lref))
       ENDIF

       !for negative q0 (CHEASE), we invert the sign of n0_global
       if (q0.lt.0.and.n0_global.gt.0) n0_global=-n0_global       

       IF (write_pe) WRITE(*,"(a,i6)") 'n0 for parallel boundary condition: ', &
            &n0_global
    ELSEIF (n0_global.eq.-1111) THEN
       IF (adapt_ly) STOP 'cannot adapt ly/kymin as rhostar or reference values are missing'
    ENDIF

    adapt_ly = adapt_ly.or.((.not.xy_local).and.(.not.lilo))

    IF (adapt_ly) kymin = n0_global/C_y0*rho_Lref

    ! set ly if possible
    if (.not.kymin.lt.epsilon(kymin)) ly = 2.0*pi/kymin

    if ((abs(ly).lt.epsilon(ly)).and.(write_pe))&
         &WRITE(*,'(A)') "WARNING: ly is equal to zero"

    if (y_local) then
       do j=0,nky0-1
          ky(j) = kymin*(j+ky0_ind)
       end do
       kymax = maxval(ky)
       dy=ly/2.0/nky0_in
    else
       do j=0,nky0/2
          ky(j) = kymin*j
       end do
       do j=nky0/2+1,nky0-1
          ky(j) = -kymin*(nky0-j)
       end do       
       kymax = kymin*(nky0/2)
       dy=ly/nky0_in
    end if
 

    if (yx_order) then
       deli=dy
       ki=ky
    else
       kjmin=kymin
       kjmax=kymax
       delj=dy
       kj=ky
    endif

  end subroutine set_y_coordinate_vars

  subroutine set_z_coordinate_vars(n_pol)
    INTEGER, intent(IN) :: n_pol

    INTEGER :: k

    !initialise parallel coordinate
    dz = (2.0 * pi * n_pol) / nz0
    do k=0,nz0-1
       zval(k) = -pi*n_pol + k*dz !+dz/2
    enddo

    if(mod(nz0,2).ne.0) then
       zval = zval+dz/2
    end if

  end subroutine set_z_coordinate_vars


  subroutine set_vp_coordinate_vars
    integer:: l
!    real,dimension(0:nv0-1):: gl_vp_weights

#ifdef COMBI
    dv = 2.0*lv/(nv0)
    do l=0,nv0-1
       vp(l) = -lv + l*dv+SHIFT
    end do
#else
    dv = 2.0*lv/(nv0-1.0)

    do l=0,nv0-1
       vp(l) = -lv + l*dv
    end do
#endif
    vp_weight = dv
    if(.not.arakawa_zv) then
       !Alternative extended Simpson's rule (see Numerical Recipes)
       if (nv0.gt.8) then
          vp_weight(0) = 17.0/48.0*dv
          vp_weight(1) = 59.0/48.0*dv
          vp_weight(2) = 43.0/48.0*dv
          vp_weight(3) = 49.0/48.0*dv
          vp_weight(nv0-4) = 49.0/48.0*dv
          vp_weight(nv0-3) = 43.0/48.0*dv
          vp_weight(nv0-2) = 59.0/48.0*dv
          vp_weight(nv0-1) = 17.0/48.0*dv
       endif
    endif

    !distribute to processors
!    vp_weight=gl_vp_weights(ll1:ll2)

  end subroutine set_vp_coordinate_vars

  subroutine set_mu_coordinate_vars
!    Real,dimension(0:nw0-1):: gl_mu_weights
    ! Gauss-Legendre integration
    ! weights and knots are calculated globally
    call GetMuWeightsAndKnots(mu_weight,mu,lw,nw0)
    
!    ! distribution to mpi processes
!    mu_weight = gl_mu_weights(lm1:lm2)

  end subroutine set_mu_coordinate_vars

  SUBROUTINE  set_radial_box_size_with_shear(n_pol,shat,kymin,&
       &nexc,lx,adapt_lx,limit_shat,Cyq0_x0)
    INTEGER, intent(IN) :: n_pol
    real, intent(IN) :: kymin, Cyq0_x0
    REAL, INTENT(INOUT) :: shat
    LOGICAL, INTENT(IN) :: adapt_lx, limit_shat
    INTEGER, intent(INOUT) :: nexc
    REAL, intent(INOUT) :: lx

    ! Local variables
    real :: lxmin, shat_lim

    !(in order to fulfill lx = N /( kymin * shat ))
    lxmin = 1.0/(n_pol*abs(shat)*kymin*Cyq0_x0)

    IF ((limit_shat).AND.(lxmin.gt.1.1*lx)) THEN
       !shat is restricted in order to 
       !avoid too large radial box widths at the innermost radial flux tubes
       shat_lim = 1.0/(n_pol*lx*kymin)
       IF (abs(shat).lt.(0.5*shat_lim)) THEN
          shat = 0.0
       ELSE
          shat = sign(shat_lim,shat)
          lxmin = 1.0/(n_pol*ABS(shat)*kymin)
       ENDIF
    ENDIF

    IF (abs(shat).gt.epsilon(shat)) then
       IF (adapt_lx) THEN !nexc=1
          lx = lxmin
          nexc = SIGN(1.0,shat)
       ELSE
          IF (nexc.EQ.0) nexc = CEILING(0.9*(lx/lxmin))
          lx = ABS(nexc)*lxmin
          !correct the sign of nexc in case of negative shear
          nexc = SIGN(1.0,shat)*ABS(nexc)
       ENDIF
    endif
  END SUBROUTINE set_radial_box_size_with_shear


  !>Calculates the center value of full radial arrays
  subroutine calc_center_value(in_arr, center_val)
    real, dimension(pi1gl:pi2gl), intent(in) :: in_arr
    real, intent(out) :: center_val

    if (mod((pi1gl+pi2gl),2).eq.0) THEN !odd number of radial points
       center_val=in_arr((pi2gl+pi1gl)/2)
    else !even number of radial points
       !estimate center by linear interpolation of neighboring grid points
       center_val=0.5*(in_arr((pi1gl+pi2gl)/2)+in_arr((pi1gl+pi2gl+1)/2))
    end if    
  end subroutine calc_center_value


  subroutine finalize_coordinates
    deallocate(kx,ky,ki,kj)
    if (.not.x_local) then
       deallocate(xval)
    end if
    deallocate(xval_a)

    deallocate(zval)
    deallocate(vp,vp_weight)
    deallocate(mu,mu_weight)

  end subroutine finalize_coordinates

end module coordinates
