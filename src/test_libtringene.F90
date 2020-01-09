#include "redef.h"
!>Program for testing the TRINITY/GENE interface
!!\todo restore functionality for local code
#define local
#ifdef local
Program test_libtringene
  Use mpi

  IMPLICIT none

  Integer :: n_transp_it, ierr
  real :: dt_max

  !MPI set-up
  logical, parameter :: reorder=.true.
  integer :: n_procs, it, comm2d, id2d, nrows, gproc
  integer :: ncolumns, comm_group, job
  integer, parameter :: ndim=2 !number of MPI dims
  integer, dimension(ndim) :: dims
  integer, dimension(0:ndim-1) :: coords1d, coords2d
  logical, dimension(0:ndim-1) :: belongs
  logical, dimension(ndim) :: period  

  !Some GENE parameters
  integer, parameter :: nx0=1   !number of radial grid points
  integer, parameter :: n_spec=2 !number of species
  integer, parameter :: bin_nr=1 !number of bins for fluxes to be returned to trinity

  !misc
  character(len=2)::  numstr, numstr2
  Character(Len=128)::filename
  !just for output reasons (should be identical with 
  !those in the GENE parameters file!)
  real :: lx = 90.
  real :: kappa_n, kappa_T !maximum gradients for profile initialization
  real :: x
  integer :: i

  !Interface variables
  real, dimension(0:n_spec-1) :: mass_in, charge_in !mass & charge for each species
  real :: shat_in           !magnetic shear (only used in local code)
  real :: coll_in, beta_in  !collision frequency and plasma beta at 
                            !reference position (typically, at x/a=0.5)

  real, dimension(bin_nr) :: rad_pos_in !radial bin center positions (should be equally spaced)
  real, dimension(1:bin_nr,0:n_spec-1) :: avgpflux, avgqflux !particle and heat fluxes for each species averaged in bins being defined by rad_pos_in
  real, dimension(1:bin_nr) :: dVdx, sqrtgxx_FS !dVdx and flux surface avgd. gxx averaged over the same bins as the fluxes
  real, dimension(1:nx0,0:n_spec-1) :: temp_io, dens_io !temperature and density profiles
  real, dimension(1:nx0,0:n_spec-1) :: omt_io, omn_io !L_ref/L_T and L_ref/L_n (only used in local code)
  real, dimension(1:nx0) :: q_in, ExBrate_in  !safety factor profile
  real, dimension(1:bin_nr,0:n_spec-1) :: temp_bin, dens_bin, omt_bin, omn_bin !temperature, density, and their gradients in bins

  Call mpi_init(ierr)
  
  Call mpi_comm_size (MPI_COMM_WORLD, n_procs, ierr)
  If (ierr /= 0) Stop 'mpi_comm_size failed!'

!number of parallel GENE runs
  ncolumns = 2 !32

!number of processes per GENE instance
  nrows = n_procs/ncolumns

!number of iterations
  n_transp_it = 3


  dims=(/ ncolumns, nrows /)     
  ! create 2d cartesian topology for processes
  period=(/ .false., .false. /)  !! no circular shift
  
  call mpi_cart_create(mpi_comm_world, ndim, dims, period, reorder, comm2d, ierr)
  call mpi_comm_rank(comm2d, id2d, ierr)
  call mpi_cart_coords(comm2d, id2d, ndim, coords2d, ierr)
  
  ! each processor knows which subgrid it is in from variable mpi_group
  job = coords2d(0)
  
  ! create 1d subgrids from 2d processor grid, variable belongs denotes
  ! whether processor grid is split by column or row
  
  belongs(1) = .true.    ! this dimension belongs to subgrid
  belongs(0) = .false.  
  
  call mpi_cart_sub(comm2d, belongs, comm_group, ierr)
  call mpi_comm_rank(comm_group, gproc, ierr)     
  call mpi_cart_coords(comm_group, gproc, 1, coords1d, ierr)
  
  !GENE input which is kept constant in the following
  mass_in(0) =   1.0
  charge_in(0) = 1.0
  if (n_spec.gt.1) then
     mass_in(1) = 2.741E-4
     charge_in(1) = -1.0
  endif

  coll_in = 0.0
  beta_in = 0.0

  rad_pos_in = 0.18 !(/0.2,0.3,0.4,0.5,0.6,0.7,0.8/) !in units of a
!  rad_pos_in = rad_pos_in*minor_r !now in units of L_ref (=R_0)

  omn_io = 1.0 !maximum density gradient
  omt_io = 8.0-3.0*(job/ncolumns) !maximum temperature gradient
  if (n_spec.gt.1) omt_io(1,1) = omt_io(1,1)*0.5

  temp_io=1.0
  dens_io=1.0

  q_in = 1.4
  ExBrate_in = 0.0
  
  !dt_max initialization (negative number for automatic computation by GENE in the first iteration)
  dt_max = -1.0

  do it=0, n_transp_it-1
!     dt_max = -1.0
     if (gproc.eq.0) then
        print*,'========================================================'
        print*, 'iteration=', it
     endif


     shat_in = 0.6+(0.2/ncolumns)*job

     Call tringene(comm_group, nx0, bin_nr, n_spec, &
          &temp_io, omt_io, dens_io, omn_io, ExBrate_in, mass_in, charge_in,&
          &q_in, shat_in, rad_pos_in, coll_in, beta_in,&
          &job, it, dVdx, sqrtgxx_FS, avgpflux, avgqflux, &
          &temp_bin,omt_bin, dens_bin, omn_bin, &
          &dt_max)

     call mpi_barrier(mpi_comm_world,ierr)
     
     if (gproc.eq.0) then
        print*,'---------------------------------------------------------'
        print*, 'iteration =', it
        print*, 'job =', job
        print*, 'temp. gradients = ', omt_io(1,:)
        print*, 'dens. gradients = ', omn_io(1,:)
        print*, 'particle fluxes = ', avgpflux(1,:)
        print*, 'heat fluxes =     ', avgqflux(1,:)
        print*, 'dVdx =            ', dVdx(1)
        print*, '<sqrt(gxx)>_FS =  ', sqrtgxx_FS(1)
        print*, 'dt_max =          ', dt_max
        print*, ''
     endif

     call mpi_barrier(mpi_comm_world,ierr)    
  end do

  Call mpi_finalize(ierr)

End Program test_libtringene
#else
Program test_libtringene
  Use mpi

  IMPLICIT none

  Integer :: n_transp_it, ierr
  real :: dt_max

  !MPI set-up
  logical, parameter :: reorder=.true.
  integer :: n_procs, it, comm2d, id2d, nrows, gproc
  integer :: ncolumns, comm_group, job
  integer, parameter :: ndim=2 !number of MPI dims
  integer, dimension(ndim) :: dims
  integer, dimension(0:ndim-1) :: coords1d, coords2d
  logical, dimension(0:ndim-1) :: belongs
  logical, dimension(ndim) :: period  

  !Some GENE parameters
  integer, parameter :: nx0=64   !number of radial grid points
  integer, parameter :: n_spec=1 !number of species
  integer, parameter :: bin_nr=7 !number of bins for fluxes to be returned to trinity

  !misc
  character(len=2)::  numstr, numstr2
  Character(Len=128)::filename
  !just for output reasons (should be identical with 
  !those in the GENE parameters file!)
  real :: lx = 90., rhostar = 0.01, minor_r = 0.36
  real :: kappa_n, kappa_T !maximum gradients for profile initialization
  real :: x
  integer :: i

  !Interface variables
  real, dimension(0:n_spec-1) :: mass_in, charge_in !mass & charge for each species
  real :: shat_in           !magnetic shear (only used in local code)
  real :: coll_in, beta_in  !collision frequency and plasma beta at 
                            !reference position (typically, at x/a=0.5)

  real, dimension(bin_nr) :: rad_pos_in !radial bin center positions (should be equally spaced)
  real, dimension(1:bin_nr,0:n_spec-1) :: avgpflux, avgqflux !particle and heat fluxes for each species averaged in bins being defined by rad_pos_in
  real, dimension(1:bin_nr) :: dVdx, sqrtgxx_FS !dVdx and flux surface avgd. gxx averaged over the same bins as the fluxes
  real, dimension(1:nx0,0:n_spec-1) :: temp_io, dens_io !temperature and density profiles
  real, dimension(1:nx0,0:n_spec-1) :: omt_io, omn_io !L_ref/L_T and L_ref/L_n (only used in local code)
  real, dimension(1:nx0) :: q_in, ExBrate_in !safety factor profile
  real, dimension(1:bin_nr,0:n_spec-1) :: temp_bin, dens_bin, omt_bin, omn_bin !temperature, density, and their gradients in bins

  Call mpi_init(ierr)
  
  Call mpi_comm_size (MPI_COMM_WORLD, n_procs, ierr)
  If (ierr /= 0) Stop 'mpi_comm_size failed!'

!number of parallel GENE runs
  ncolumns = 2

!number of processes per GENE instance
  nrows = n_procs/ncolumns

!number of iterations
  n_transp_it = 2


  dims=(/ ncolumns, nrows /)     
  ! create 2d cartesian topology for processes
  period=(/ .false., .false. /)  !! no circular shift
  
  call mpi_cart_create(mpi_comm_world, ndim, dims, period, reorder, comm2d, ierr)
  call mpi_comm_rank(comm2d, id2d, ierr)
  call mpi_cart_coords(comm2d, id2d, ndim, coords2d, ierr)
  
  ! each processor knows which subgrid it is in from variable mpi_group
  job = coords2d(0)
  
  ! create 1d subgrids from 2d processor grid, variable belongs denotes
  ! whether processor grid is split by column or row
  
  belongs(1) = .true.    ! this dimension belongs to subgrid
  belongs(0) = .false.  
  
  call mpi_cart_sub(comm2d, belongs, comm_group, ierr)
  call mpi_comm_rank(comm_group, gproc, ierr)     
  call mpi_cart_coords(comm_group, gproc, 1, coords1d, ierr)
  
  !GENE input which is kept constant in the following
  mass_in =   1.0
  charge_in = 1.0
  coll_in = 0.0
  beta_in = 0.0
  ExBrate = 0.0

  rad_pos_in = (/0.2,0.3,0.4,0.5,0.6,0.7,0.8/) !in units of a
  rad_pos_in = rad_pos_in*minor_r !now in units of L_ref (=R_0)

  kappa_n = 1.0 !maximum density gradient
  kappa_T = 7.0-job !maximum temperature gradient

  !initialization of temperature, density, and safety factor profiles
  !logarithmic gradients will be computed in GENE
  do i=1,nx0
     temp_io(i,:) = exp(-kappa_T*minor_r*(rhostar*(real(i-1)*lx/(nx0-1)-0.5*lx)))
     dens_io(i,:) = exp(-kappa_n*minor_r*(rhostar*(real(i-1)*lx/(nx0-1)-0.5*lx)))
     q_in(i) = 0.8+i*0.4/nx0
  enddo

  
  !dt_max initialization (negative number for automatic computation by GENE in the first iteration)
  dt_max = -1.0

  do it=0, n_transp_it-1
!     dt_max = -1.0
     if (gproc.eq.0) then
        print*,'========================================================'
        print*, 'iteration=', it
     endif


     shat_in = 0.6+(0.2/ncolumns)*job

     Call tringene(comm_group, nx0, bin_nr, n_spec, &
          &temp_io, omt_io, dens_io, omn_io, ExBrate_in, mass_in, charge_in,&
          &q_in, shat_in, rad_pos_in, coll_in, beta_in,&
          &job, it, dVdx, sqrtgxx_FS, avgpflux, avgqflux, &
          &temp_bin,omt_bin, dens_bin, omn_bin, &
          &dt_max)

     call mpi_barrier(mpi_comm_world,ierr)
     
     if (gproc.eq.0) then
        print*,'---------------------------------------------------------'
        print*, 'iteration =', it
        print*, 'job =', job
        print*, 'temp. gradients = ', omt_io(1,:)
        print*, 'dens. gradients = ', omn_io(1,:)
        print*, 'particle fluxes = ', avgpflux(1,:)
        print*, 'heat fluxes =     ', avgqflux(1,:)
        print*, 'dVdx =            ', dVdx(1)
        print*, '<sqrt(gxx)>_FS =  ', sqrtgxx_FS(1)
        print*, 'dt_max =          ', dt_max
        print*, ''

        !FILE I/O
        write(numstr,"(i2.2)") job
        write(numstr2,"(i2.2)") it

        filename = './profile_'//numstr//'_'//numstr2
        print*, 'writing profiles to file '//TRIM(filename)
        OPEN(97,file=trim(filename))
        write(97,"(A)") "#   x/a             T/Tref          n/nref          "//&
             &"Lref/LT         Lref/Ln"
        do i=1,nx0
           x=(-0.5+(i-1)/real(nx0-1))*lx*rhostar+0.5
           write(97,"(5ES16.6)") x, temp_io(i,0), dens_io(i,0), omt_io(i,0), omn_io(i,0)
        enddo
        write(97, "(A)") ""
        write(97, "(A)") ""
        write(97,"(A)") "#  bin_pos          T/Tref          n/nref          Lref/LT         "//&
             &"Lref/Ln         Gamma           Q               dVdx            <sqrtgxx>_FS"

        do i=1,bin_nr
           write(97,"(9ES16.6)") rad_pos_in(i)/minor_r, temp_bin(i,0), dens_bin(i,0), &
                &omt_bin(i,0), omn_bin(i,0), avgpflux(i,0), &
                &avgqflux(i,0), dVdx(i), sqrtgxx_FS(i)
        enddo
        close(97)
        print*,'---------------------------------------------------------'
     endif

     call mpi_barrier(mpi_comm_world,ierr)    
  end do

  Call mpi_finalize(ierr)

End Program test_libtringene
#endif
