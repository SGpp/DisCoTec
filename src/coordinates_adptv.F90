! 06-05-2012
  
#include "redef.h"
!>Contains the values of the coordinates for block-structured grid
!!and the routines to set them
module coordinates_adptv_module
!  use discretization
  use discretization_adptv_module, only: nv0_adp2, n_vx_blks, blk_mks_v
!  use GaussQuadrature 
  use par_in, only: arakawa_zv
  !use Grid1DModule
  implicit none

  public :: dv_adp2, vp_adp2, vp_weight_adp2
  public :: initialize_coordinates_adptv2, set_vp_coordinate_vars_adptv2, &
       & finalize_coordinates_adptv2

  ! Box lengths (not necessery, we have already them from discritization_adptv_module)
  ! real, dimension(:), allocatable   :: lv_adp2
  ! Grid spacing
  real, dimension(:), allocatable   :: dv_adp2
    ! coordinates
  real, dimension(:,:), allocatable :: vp_adp2
  ! velocity space integration factors
  real, dimension(:,:), allocatable :: vp_weight_adp2

  !private module variables
  real, parameter   :: pi= 3.141592653589793239d0

  private

contains

  subroutine initialize_coordinates_adptv2

    allocate(dv_adp2(0:n_vx_blks-1)) !, lv_adp2(0:n_vx_blks-1))

    allocate(vp_adp2(0:nv0_adp2-1,0:n_vx_blks-1), &
         & vp_weight_adp2(0:nv0_adp2-1,0:n_vx_blks-1))
    
  end subroutine initialize_coordinates_adptv2

  subroutine set_vp_coordinate_vars_adptv2
    integer :: l, blk

    dv_adp2 = 2.0*blk_mks_v(:n_vx_blks)/(nv0_adp2-1.0)

    do blk = 0, n_vx_blks - 1
       do l = 0, nv0_adp2 - 1
          vp_adp2(l,blk) = -blk_mks_v(blk+1) + l*dv_adp2(blk)
       end do
    end do

    do blk = 0, n_vx_blks - 1
       vp_weight_adp2(:,blk) = dv_adp2(blk)
       if(.not.arakawa_zv) then
          !Alternative extended Simpson's rule (see Numerical Recipes)
          if (nv0_adp2.gt.8) then
             vp_weight_adp2(0,blk) = 17.0/48.0*dv_adp2(blk)
             vp_weight_adp2(1,blk) = 59.0/48.0*dv_adp2(blk)
             vp_weight_adp2(2,blk) = 43.0/48.0*dv_adp2(blk)
             vp_weight_adp2(3,blk) = 49.0/48.0*dv_adp2(blk)
             vp_weight_adp2(nv0_adp2-4,blk) = 49.0/48.0*dv_adp2(blk)
             vp_weight_adp2(nv0_adp2-3,blk) = 43.0/48.0*dv_adp2(blk)
             vp_weight_adp2(nv0_adp2-2,blk) = 59.0/48.0*dv_adp2(blk)
             vp_weight_adp2(nv0_adp2-1,blk) = 17.0/48.0*dv_adp2(blk)
          endif
       endif
    end do

  end subroutine set_vp_coordinate_vars_adptv2

  subroutine finalize_coordinates_adptv2
 
    if (allocated(dv_adp2))        deallocate(dv_adp2)
    ! if (allocated(lv_adp2))        deallocate(lv_adp2)
    if (allocated(vp_adp2))        deallocate(vp_adp2)
    if (allocated(vp_weight_adp2)) deallocate(vp_weight_adp2)
    
  end subroutine finalize_coordinates_adptv2
  
end module coordinates_adptv_module
