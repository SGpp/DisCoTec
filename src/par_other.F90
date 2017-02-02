#include "redef.h"
!Temporary module containing more global parameters
module par_other

  implicit none

!PART1: used by parameters_IO

#ifdef DOUBLE_PREC
  Character(len=6), parameter :: prec='DOUBLE'
#else
  Character(len=6), parameter :: prec='SINGLE'
#endif
#if defined(SVN_REV)
  character(len=20) :: svn_rev = SVN_REV
#else
  character(len=20) :: svn_rev = ''
#endif

  character(len=13), parameter :: release = '1.8 - alpha 0'

  Real::    dt              ! Time step

  LOGICAL :: equil_par_curr = .false.

  !diagnostic
  Integer:: n_moms
  Integer:: n_fields=1
  !INTEGER :: n_energies = 2

  real    :: initime, time_ev, time_iv, time_nc

  !scan
  character(Len=128):: par_in_dir=''

  !parameters for external potentials or temperatures
  Logical :: with_phi_ext, with_omt_ext
  Logical :: with_omn_ext, with_apar_ext

  logical:: only_zonal=.false.

  logical:: print_ini_msg=.true.
  logical :: p_has_0_mode=.false., p_has_00_mode=.false.

  Logical:: auto_parall=.true. ! Automatic parallelization switch
  Logical:: final_init=.false.
  Integer:: itime=0
  
  !divide psi by 2pi in tracer_efit interface
  logical:: psi_o_twopi=.false.
  !rescaling of EQDSK q profile
  real:: q_scalefac=1.0



  !parallel current density term
  Real, DIMENSION(:), Allocatable :: currdens_par
  real:: currdens_par_scal

  !parameters used in library interfaces
  !TRINITY interface: add flux tube index to geomfile
  logical :: multiple_tracer_files = .false.

!PART2: used by check_parameters

  !rel. deviation limit where f0 and prefactors
  !should be reset
  real :: reset_limit = 1000.

  !neglect fields, collisions and x-derivatives part of the linear operator
  Logical:: precond_approx=.false.

  REAL, parameter   :: pi= 3.141592653589793239d0
  complex,parameter :: imag=(0.0,1.0)

  logical:: xy_local
  logical:: hypz_compensation
  logical:: nonlin_h

  logical :: only_neo
  
  logical :: in_perf_opt

  
  contains
    !>Set par_other variables to default values
    subroutine set_par_other_defaults
      time_ev = 0.0
      time_iv = 0.0
      time_nc = 0.0
      hypz_compensation=.true.
      in_perf_opt = .false.
      nonlin_h=.false.
    end subroutine set_par_other_defaults

end module par_other
