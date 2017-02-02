#include "redef.h"
MODULE box_data_module
  USE Grid1DModule
  USE par_mod, ONLY: n_pol, nexc, &
       &n0_global, lp1, x0, mype,print_ini_msg, lilo, kymin
  implicit none
  private

  PUBLIC :: box_data_type
  PUBLIC :: initialize, finalize, get_n0_global, get_kymin, get_lx, get_shat, set_grids, print_box
  PUBLIC :: get_xposition, get_q0, get_lilo, get_rhostar
  PUBLIC :: HasFiniteShear, HasMagneticProfile

  TYPE box_data_type
     INTEGER :: n0_global
     REAL    :: lp1
     REAL    :: x0 !< position of x==0 relative to minor radius
     REAL    :: kymin
     REAL    :: shat, q0
     REAL    :: rhostar
     LOGICAL :: has_magnetic_profile
     logical :: lilo
     TYPE(Grid1D), DIMENSION(5)  :: grid !1->x, 2->y, 3->z, 4->v, 5->w
  END TYPE box_data_type

  INTERFACE initialize
     MODULE PROCEDURE box_data_initialize
  END INTERFACE

  INTERFACE finalize
     MODULE PROCEDURE box_data_finalize
  END INTERFACE

  INTERFACE set_grids
     MODULE PROCEDURE box_data_set_grids
  END INTERFACE

  INTERFACE print_box
     module procedure box_data_print_box
  END INTERFACE

  INTERFACE get_n0_global
     MODULE PROCEDURE box_data_get_n0_global
  END INTERFACE

  INTERFACE get_kymin
     MODULE PROCEDURE box_data_get_kymin
  END INTERFACE

  INTERFACE get_q0
     MODULE PROCEDURE box_data_get_q0
  END INTERFACE

  INTERFACE get_lilo
     MODULE PROCEDURE box_data_get_lilo
  END INTERFACE

  INTERFACE get_lx
     MODULE PROCEDURE box_data_get_lp1
  END INTERFACE

  INTERFACE get_shat
     MODULE PROCEDURE box_data_get_shat
  END INTERFACE

  INTERFACE get_rhostar
     MODULE PROCEDURE box_data_get_rhostar
  END INTERFACE

  INTERFACE get_xposition
     MODULE PROCEDURE box_data_get_xposition
  END INTERFACE

CONTAINS
  
  SUBROUTINE box_data_initialize(this,q0, shat, mag_prof, &
            &rhostar)
    type(box_data_type) :: this
    real, intent(in) :: q0, shat, rhostar
    logical, intent(in) :: mag_prof

    ! read in the values from par_mod
    this%n0_global = n0_global
    this%has_magnetic_profile = mag_prof
    this%shat = shat
    this%q0 = q0
    this%lp1 = lp1
    this%x0 = x0
    this%lilo = lilo
    this%kymin = kymin
    this%rhostar = rhostar

    !the coordinate adaptions have been moved to coordinates.F90

  END SUBROUTINE box_data_initialize

  SUBROUTINE box_data_finalize(bdt)
    type(box_data_type) :: bdt
  END SUBROUTINE box_data_finalize

  SUBROUTINE box_data_set_grids(this,xgrid,ygrid,zgrid,vgrid,wgrid)
    type(box_data_type) :: this
    TYPE(grid1d) :: xgrid,ygrid,zgrid,vgrid,wgrid

    this%grid(1) = xgrid
    this%grid(2) = ygrid
    this%grid(3) = zgrid
    this%grid(4) = vgrid
    this%grid(5) = wgrid
  END SUBROUTINE box_data_set_grids

  SUBROUTINE box_data_print_box(this)
    TYPE(box_data_type) :: this

    ! Local variables
    integer :: iGrid
    PRINT*, "n0_global = ", this%n0_global
    PRINT*, "lp1        = ", this%lp1
    PRINT*, "x0        = ", this%x0
    PRINT*, "kymin     = ", this%kymin
    PRINT*, "shat      = ", this%shat
    PRINT*, "q0        = ", this%q0
    PRINT*, "has_magnetic_profile = ", this%has_magnetic_profile
    DO iGrid=1,5
       PRINT*,"grid ",iGrid," has properties:"
       call print_grid(this%grid(iGrid))
    END DO
  END SUBROUTINE box_data_print_box

  FUNCTION box_data_get_n0_global(this)
    TYPE(box_data_type) :: this
    INTEGER :: box_data_get_n0_global

    box_data_get_n0_global = this%n0_global
  END FUNCTION box_data_get_n0_global

  FUNCTION box_data_get_kymin(this)
    TYPE(box_data_type) :: this
    REAL :: box_data_get_kymin

    box_data_get_kymin = this%kymin
  END FUNCTION box_data_get_kymin

  FUNCTION box_data_get_lilo(this)
    TYPE(box_data_type) :: this
    logical :: box_data_get_lilo

    box_data_get_lilo = this%lilo
  END FUNCTION box_data_get_lilo

  FUNCTION box_data_get_q0(this)
    TYPE(box_data_type) :: this
    REAL :: box_data_get_q0

    box_data_get_q0 = this%q0
  END FUNCTION box_data_get_q0

  FUNCTION box_data_get_lp1(this)
    TYPE(box_data_type) :: this
    REAL :: box_data_get_lp1

    box_data_get_lp1 = this%lp1
  END FUNCTION box_data_get_lp1

  FUNCTION box_data_get_shat(this)
    TYPE(box_data_type) :: this
    REAL :: box_data_get_shat

    box_data_get_shat = this%shat
  END FUNCTION box_data_get_shat

  FUNCTION box_data_get_rhostar(this)
    TYPE(box_data_type) :: this
    REAL :: box_data_get_rhostar

    box_data_get_rhostar = this%rhostar
  END FUNCTION box_data_get_rhostar

  FUNCTION box_data_get_xposition(this,index)
    TYPE(box_data_type) :: this
    INTEGER, intent(IN) :: index
    REAL :: box_data_get_xposition

    box_data_get_xposition = get_node(this%grid(1),index)
  END FUNCTION box_data_get_xposition

  LOGICAL FUNCTION HasFiniteShear(this)
    TYPE(box_data_type) :: this
    
    HasFiniteShear = (ABS(this%shat).GT.1e-3)
  END FUNCTION HasFiniteShear

  LOGICAL FUNCTION HasMagneticProfile(this)
    type(box_data_type) :: this

    HasMagneticProfile = this%has_magnetic_profile
  END FUNCTION HasMagneticProfile
END MODULE box_data_module
