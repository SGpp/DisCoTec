!--------------------------------------------------------------------------
!   Module for definition of datastructures of GENE
!--------------------------------------------------------------------------
!>Module containing all parameters and many global variables
#include "redef.h"
Module par_mod

  use discretization
  use coordinates
  use par_in
  use par_other
  use, intrinsic :: iso_c_binding

  Implicit None

  Logical :: do_energy=.false.
  Logical :: in_timestep_comp=.false.

  Real:: phase_1, phase_2

  Real::    time=0.            ! Current Time
 
  !courant limit
  Real, DIMENSION(4) :: maxclimit = -100.
  REAL(C_REAL_TYPE), DIMENSION(2),bind(c) :: ve_max
  REAL :: ve_pnl_max

  !distribution functions and fields
  Complex, Dimension(:, :, :, :, :,:), Allocatable:: g_1

  Logical:: exi

  !finite differences 
  Real, dimension(:,:,:,:),allocatable:: vderivative
  real:: np1coeff(-2:2),np2coeff(-2:2),np3coeff(-2:2),np4coeff(-2:2)

  Real :: dvfac

  Real, dimension(:), Allocatable:: par_sten
  Integer:: par_sten_bound

  !derived rhs computation
  logical:: direct_rhs=.true.
  integer:: derived_method=0

  logical :: my_mpi_comm_world_member=.true.
  integer,parameter :: nag_dbl=kind(1.0D0)      

  REAL, PARAMETER :: m_proton   = 1.672621637E-27 !kg
  REAL, PARAMETER :: m_deuteron = 3.34358320E-27 !kg
  REAL, PARAMETER :: m_electron = 9.10938215E-31 !kg


#ifdef USE_VT
! Set variables for Intel Trace Analyzer  
  INTEGER :: mylib_handle, handle_index, ierr
  INTEGER, DIMENSION(100) :: name_handle_arr
#endif

End Module par_mod
