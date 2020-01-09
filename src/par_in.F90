#include "redef.h"
#include "intrinsic_sizes.h"
!>Temporary module containing the GENE input parameters which are not related to 
!!coordinates or discretization. Should shrink in the future.
!!Please add further parameters only if they cannot be defined in the 
!!modules where they are actually used!!!!
module par_in
  use spectype_mod

  use,intrinsic :: iso_c_binding
  implicit none
  
  !numerical schemes
  Character(len=6):: timescheme
  Character(len=6):: coll_split_scheme
  logical :: coll_split
  Character(len=5):: parscheme, vparscheme
  Character(len=2):: comp_type
  
  ! GENE type
  Logical:: x_local
  Logical:: y_local
  Logical:: nonlinear, parallel_nl, arakawa, ExB=.false.
  Logical:: lilo !run nonlocal code LI-ke LO-cal
  logical:: trap_pass, no_electron_response
  logical:: approx_avg
  logical:: only_passing
  logical:: precomp_nc
  logical:: include_f0_contr

  !collisions
  Character(len=12) :: collision_op
  logical :: coll_f_fm_on

  !dealiasing
  Logical:: turbdeal
  !integer :: lag_order
  
  !parallel differencing scheme
  logical:: arakawa_zv, hyp_on_h, hypz_opt

  !performance
  integer,dimension(9):: perf_vec
  integer :: nblocks !number of blocks for strip mining in calc_rhs_only
  logical :: opt_nblocks
  logical :: fourier2d
  real :: gpu_cpu_ratio

  !in/output
  character(len=20) :: file_extension='.dat'
  Character(Len=128):: diagdir
  Logical:: write_h5=.false.
  Logical:: write_std=.true.

  !modifications of electromagnetic fields
  Logical:: delzonal
  logical:: delzonal_fields
  logical:: del_phi, del_fields
  logical:: only_Er
  real :: delzonal_factor=1.0 !facilitates a scan in zonal flow amplitude
  logical:: nc_exb_corr = .true.
  real :: add_zonal_phi=0.0

  !dimensionless plasma parameters
  Real::    beta, debye2
  !tau=Zeff*Te/Ti for adiabatic ion simulations
  Real:: tau

  !Comoving frame
  Real :: Omega0_tor  !toroidal angular velocity in cref/Lref

  !radial boundary
  integer :: rad_bc_type !(0=periodic, 1=Dirichlet, 
                         ! 2=v.Neumann for ky=0 at inner boundary else Dirichlet)

  integer :: num_ky_modes, num_kx_modes, which_ky(50)
  integer :: which_kx_center(50)
  logical :: vec_out_limit !flag to limit number of eigenvectors output from scalapack routine
  Logical :: calc_dt=.false.
  Logical :: nltdt_off=.false. !flat to switch off the nonlinear time step adaptation
  
  !flux surface averaged velspace moments
  integer :: istep_fsa_moments
  
  !energy
  integer :: istep_energy=-1, istep_energy3d=-1, istep_energy_exchange = -1
  integer :: num_nlt_modes=5                 !Number of wavenumbers for which to calculate nonlinear transfer functions
                                             !Maximum of twenty
  integer :: kx_nlt_ind(60)=0,ky_nlt_ind(60)=0     !indices for wavenumbers of interest for nonlinear transfer functions
  integer :: istep_nlt
  logical :: nlt_symmetrize = .true.  !uses q,p symmetrized nonlinear transfer function (see Nakata '12)

  !eigenvalues
  Logical:: ev_out
  integer :: which_func_as=1

  LOGICAL :: bpar_off, bpar, pressure_term, pressure_off

  LOGICAL :: write_fielddiff

  !dfout
  logical :: dfout_mpiio
  !SVD diagnostics
  logical :: SVD_proj
  logical :: SVD_f0_flag
  integer :: SVD_kx_ind
  integer :: SVD_ky_ind
  integer :: SVD_nkx0
  !real    :: SVD_kymin
  integer :: SVD_df_n_time(20)
  real ::    SVD_start_time
  integer :: SVD_sparse_factor
  character(len=FILENAME_MAX) :: SVD_df_file_name
  character(len=128) :: SVD_df_file_path
  character(len=128) :: SVD_df_file_suffix(20)
  character(len=128) :: SVD_endian
  integer :: SVD_n_restarts
  logical :: SVD_pbc_phase
  logical :: SVD_parallel
  logical :: SVD_field_mom
  logical :: SVD_swap_endian
  integer :: SVD_istep_mf_factor
  !nlt pod routine
  integer :: num_nlt_pod_modes
  integer :: nx0_nlt_pod(20)
  logical :: nlt_pod
  !gyroLES diagnostics
  integer :: istep_fe_transfer = 0  

  !shifted metric
  logical:: shifted_metric

  !parallel boundary condition
  integer:: nexc
  logical:: adapt_lx
  integer:: n_pol

  Real:: courant

  !Erad for stellarator runs
  Real:: Erad
  logical:: Erad_acc 
  !real time
  double precision :: realtime_start
  Integer:: timelim

  !simulation time
  Real::    simtimelim
  Integer:: ntimesteps
  Real::    dt_max      ! maximum time step
  !relevant for coll_split:
  real::    dt_vlasov   ! Maximum Time step for collisionless part
  real::    ev_coll     ! minimum real ev for collisions

  integer :: n0_global
  Real::    overflow_limit, underflow_limit

  type(spectype),dimension(:),allocatable:: spec
  
  !variables for Trinity transport code
  Real, dimension(:,:,:), allocatable :: in_profiles
  real, dimension(:), allocatable :: in_qprof
  
  logical:: ga_spatial_var=.false.

  !Variables for GLES
  real, public :: fracx = 0.5
  real, public :: fracy = 0.5
  integer, public :: istep_GyroLES = 50
  integer, public :: adapt_istep_GyroLES = 0
  logical, public :: GyroLES = .false.
  logical, public :: diag_GyroLES = .false.

  contains
    
    !>Set par_in variables to default values
    subroutine set_par_in_defaults

      !numerical schemes
      timescheme=''
      coll_split_scheme = 'RKCa'      
      coll_split = .true.
      parscheme='c4th'
      vparscheme='c4th'
      comp_type='IV'
      
      ! GENE type
      parallel_nl = .false.
      x_local=.true.
      y_local=.true.
      arakawa=.true.
      lilo=.false.
      trap_pass=.false.
      precomp_nc=.false.
      include_f0_contr=.false.

      !adiabatic electron response:
      approx_avg=.false.
      only_passing=.false.

      !collisions
      collision_op='none'
      coll_f_fm_on = .false.

      !comoving
      Omega0_tor = 0.0

      !dealiasing
      turbdeal=.false.
      !lag_order=9
      
      !parallel differencing scheme
      arakawa_zv=.true.
      !hyperdiffusion on h in case of arakawa
      hyp_on_h=.true.
      hypz_opt=.true.
      !performance
      perf_vec=0
      nblocks=0
      fourier2d=.false. ! Use/don't use 2D fourier transforms for nonlinear calculations
      
      !in/output
      diagdir=''

      !electromagnetic field modifications
      delzonal=.false.
      delzonal_fields=.false.
      del_phi = .false.
      del_fields = .false.
      only_Er = .false.
      nc_exb_corr = .true.
      add_zonal_phi  = 0.0

      !dimensionless plasma parameters
      beta=0.0
      debye2=0.0
      tau=1.0

      Erad=0.
      Erad_acc=.false. 

      !radial boundary
      rad_bc_type = 1

      num_ky_modes=1 ; num_kx_modes=1 ; which_ky(15)=0
      which_kx_center(15)=0
      calc_dt=.false.

      istep_fsa_moments=0
      
      !energy
      istep_energy=-1
      istep_energy3d=-1

      !eigenvalues
      ev_out=.false. 

      bpar_off=.true.; bpar=.false. ; pressure_term=.false.; pressure_off=.true.

      write_fielddiff=.false.

      !dfout
      dfout_mpiio=.false.
      !SVD diagnostics
      SVD_proj=.false.
      SVD_f0_flag=.true.
      SVD_kx_ind=-100
      SVD_ky_ind=-100
      SVD_nkx0=-100
      !SVD_kymin=-0.05
      SVD_df_n_time=-100
      SVD_start_time=0.0 
      SVD_sparse_factor=1
      SVD_df_file_name='no_input'
      SVD_df_file_path='no_input'
      SVD_df_file_suffix(1)='.dat'
      SVD_n_restarts=1
      SVD_pbc_phase=.false.
      SVD_parallel=.false.
      SVD_field_mom=.false.
      SVD_swap_endian=.false.
      !nlt POD
      num_nlt_pod_modes=5
      nx0_nlt_pod=1
      nlt_pod=.false.
      !diag_GyroLES (moved here, because needs h_ in case of arakawa_zv=F)
      istep_fe_transfer=0


      !shifted metric
      shifted_metric=.false.

      !parallel boundary condition
      nexc=0
      adapt_lx=.true.
      n_pol=1
      
      courant=1.25  !new version that considers max. nonlin. eigenvalue
                    !often allows courant>1

      timelim=500000
      
      !simulation time
      simtimelim=500000.
      ntimesteps=1000000000
      dt_max=-1.
      !coll_split timestep limits
      dt_vlasov=-1.
      ev_coll=-1

      n0_global = -1111
      overflow_limit = 1e30 ; underflow_limit=1e-8

    end subroutine set_par_in_defaults

end module par_in
