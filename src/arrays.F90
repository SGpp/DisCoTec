#include "redef.h"
#include "intrinsic_sizes.h"

!>Contains the allocation routine for some arrays, most of which should be moved to other
!!modules. Introduced to avoid cirular dependencies when introduced discretization module.
!!\todo remove the module and put the content somewhere else!
module arrays
  use par_mod
  use discretization
  
  implicit none
contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_arrays(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !modified distribution function
    mem_loc=SIZE_OF_COMPLEX_MB*lijklmn0

    mem_est_arrays=mem_req_in+mem_loc
  End Function mem_est_arrays

  !*************************************************************************!
  !************ Init parameters concerning problem size ********************!
  !*************************************************************************!
  !>Allocate arrays used for reading the parameters
  !!Currently, this subroutine includes many more arrays which should in 
  !!future be moved to the corresponding modules
  Subroutine allocate_arrays
   
    if(.not.SVD_field_mom) allocate(g_1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))

    call initialize_coordinates

    ALLOCATE(currdens_par(pi1:pi2))
    currdens_par = 0.0

  End Subroutine allocate_arrays


  !>Deallocates all arrays
  Subroutine deallocate_arrays

    if(.not.SVD_field_mom) deallocate(g_1)
    call finalize_coordinates
    
    deallocate(currdens_par)

  End Subroutine deallocate_arrays
end module arrays
