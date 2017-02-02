  MODULE tga_io
  USE localpolynombase_mod
  USE communications, ONLY: MY_MPI_COMM_WORLD,COMM_X, communicators,&
       &mpi_comm_x, mpi_comm_spec,my_barrier, calculate_test_sum,&
       &mpi_comm_v
  USE discretization
  USE coordinates
  USE par_other
  USE par_in
  USE geometry
  USE BoundaryDescriptionModule
  use MatrixModule
  USE BandedMatrixModule
  use VectorModule
  USE mpi
  USE DerivativeMatrixModule
  use Grid1DModule
  USE file_io
  USE gyro_average_df_mod
  USE test_pmm
use par_in
use par_other
use file_io, only: get_unit_nr
use discretization
use coordinates
use collisions
use geometry
!use initial_value_comp, only: rescale
use communications,only: omp_num_threads
use GaussQuadrature, only: mu_grid_type
USE boundary_exchange_z, ONLY: set_pbc_dealiasing,get_pbc_dealiasing
USE gyro_average_df_mod,only: gyroav_in_xi_eta
use profiles_mod, only: Tref, nref, mref, drive_buffer, set_prof_defaults, &
     &set_in_prof_defaults, iterdb_file, iterdb_time, norm_index
use dzv_terms
use numerical_damping

 
  IMPLICIT NONE

  public:: read_tga_parameters,  par_out_file, PARFILE

    REAL :: omn, omt, mass, temp, dens, &
         & kappa_T, LT_center, LT_width, &
         & kappa_n, Ln_center, Ln_width, &
         & delta_x_T, delta_x_n
    character(len=10) :: name
    real, Dimension(0:4) :: src_amp, src_width, src_x0, &
         & src_prof_type
    real :: upar0
    integer:: charge, n, prof_type
    logical:: passive
    character(len=1) :: prof_file

    integer:: PAR_OUT_FILE, PARFILE=30

    logical:: pbc_dealiasing=.false.
    integer:: f_ext_startnum=1
    integer:: j

  namelist /in_out/ &
     diagdir,  &
     write_fielddiff, &
     istep_energy, f_ext_startnum, &
     write_h5, write_std, &
     istep_nlt, istep_GyroLES, &
     iterdb_file, iterdb_time, &
     istep_fsa_moments


    CONTAINS

  SUBROUTINE read_tga_parameters(par_in_dir)
  
  ! This subroutine is a modified version of the read_parameters in the main GENE trunk.
  ! Some of the namelists, which are for parts of the code unrelated to the gyroaverage 
  ! test suite, have been removed.
  !  

  !>Reads parameters from ./parameter.dat file
    character(len=*), intent(in) :: par_in_dir
    character(len=10) :: name
    character(len=128) :: prof_file, par_file
    logical:: passive
    real:: omn, omt, mass, temp, dens, & 
         & kappa_T, LT_center, LT_width, & 
         & kappa_n, Ln_center, Ln_width, &
         & delta_x_T, delta_x_n
    real, Dimension(0:4) :: src_amp, src_width, src_x0, &
         & src_prof_type
    integer:: charge, n, ierr, prof_type

    namelist /box/ &
         nx0, nky0, nz0, nv0, nw0, n_spec, &
         lx, kymin, lv, lw, ky0_ind, kx_center, nexc, adapt_lx, x0,&
         mu_grid_type,&
         n0_global  ! for parallel boundary condition

    namelist /general/ &
         nonlinear, parallel_nl, x_local, y_local, lilo, comp_type, &       !operation type
         perf_vec, nblocks, parscheme, turbdeal, precomp_nc, include_f0_contr,& !rescale
         timescheme, calc_dt, dt_max, courant, & !time scheme
         ntimesteps, timelim, simtimelim, &  !simulation stop criteria
         overflow_limit, underflow_limit, &
         beta, debye2, delzonal, delzonal_fields, bpar_off, pressure_off, &            !physical parameters + related
         Apar0_antenna, omega0_antenna, del_phi, only_Er, trap_pass, tau, &
         coll, collision_op, coll_cons_model, em_conserve, spacediff_off, coll_f_fm_on, Zeff,&
         hyp_x, hyp_x_order, hyp_y, hyp_y_order, &         !hyper diffusion
         hyp_perp, hyp_perp_order, &
         hyp_z, hyp_z_order, hyp_v, &
         GyroLES, fracx,fracy, diag_GyroLES, &
         arakawa_zv, arakawa_cons_bc, arakawa_zv_order, fourier2d, yx_order !misc
 

    namelist /nonlocal_x/ &
         rad_bc_type, dealiasing, pbc_dealiasing, lag_order, gyroav_in_xi_eta,&
         arakawa, shifted_metric, &
         !sources
         reset_limit, drive_buffer, &
         ga_spatial_var

    !variables for external (zonal) fields, temperatures and densities
   
    namelist /geometry/ &
         magn_geometry, shat, limit_shat, &     !important for all geometries
         parscale, q0, major_R,minor_r, trpeps, amhd, &  !parscale for slab, rest for s-alpha
         geomdir, geomfile, x_def, flux_pos, &    !general geometry
         rhostar, mag_prof, q_coeffs, n_pol, edge_opt, &
         norm_flux_projection

    namelist /species/ &
         & name, passive, omn, omt, mass, charge, temp, dens, & 
         & kappa_T,LT_center,LT_width, kappa_n,Ln_center,Ln_width,&
         & delta_x_T, delta_x_n, prof_type,&
         & src_amp, src_width, src_x0, src_prof_type, prof_file

    namelist /units/ &
         & Tref, nref, Lref, Bref, mref, norm_index

    !default values
    call set_par_other_defaults
    call set_par_in_defaults
 
    call set_nd_defaults

  
    call set_prof_defaults
    call set_coll_defaults

    f_ext_startnum = 1
    q_coeffs=0.0
 
    if(allocated(spec)) deallocate(spec)

    ! get the defaults from the module
    pbc_dealiasing = get_pbc_dealiasing()
 
    if (PARFILE.eq.5) then
       if (mype==0) WRITE(*,"(A)") "reading parameters from stdin"
    else 
       call get_unit_nr(PARFILE)  
       par_file = 'parameters'
       if (par_in_dir.ne.'') &
            par_file = trim(par_in_dir)//'/'//trim(par_file)//&
            &trim(file_extension)
       open(PARFILE, file=trim(par_file), form='formatted',status='old',&
            &iostat=ierr)
       if (ierr.ne.0) then
          write (*,'(2A)') 'could not read ', trim(par_file)
          stop 'failed to read parameters file'
       endif
    endif

    !parallelization namelist has already been read by initialize_comm_sim
    read(PARFILE, nml=box,iostat=ierr)
    if (ierr.ne.0) stop 'on i/o error: incorrect box namelist'

    rewind(PARFILE)
    read(PARFILE, nml=in_out,iostat=ierr)
    if (ierr.ne.0) stop 'on i/o error: incorrect in_out namelist'

    rewind(PARFILE)
    read(PARFILE, nml=general,iostat=ierr)
    if (ierr.ne.0) stop 'on i/o error: incorrect general namelist' 

    rewind(PARFILE)
    read(PARFILE, nml=nonlocal_x,iostat=ierr)
    if (ierr.gt.0) stop 'on i/o error: incorrect nonlocal_x namelist'

    rewind(PARFILE)
    read(PARFILE, nml=geometry,iostat=ierr)
    if (ierr.ne.0) stop 'on i/o error: incorrect geometry namelist'

    allocate(spec(0:n_spec-1))

    nky0_in = nky0

    do n=0,n_spec-1
       !set default for prof_type
       if (x_local) then 
          prof_type=0
       else
          prof_type=2
       endif
       
       !set spec namelist defaults
       passive = .false.
       omn = 0.0
       omt = 0.0
       temp = 1.0
       dens = 1.0
       kappa_T = -1.0
       LT_center = 0.5
       LT_width = 0.1
       kappa_n=-1.0
       Ln_center=0.5
       Ln_width=1.0
       delta_x_T=0.4
       delta_x_n=0.4
       src_amp = 0.0
       src_width=0.1
       src_x0 = 0.5
       src_prof_type = 0
       prof_file = ''

       read(PARFILE,nml=species)
       
       call set_spec_nonprof(spec(n),name,passive,charge,mass,temp, &
            & dens, omn, omt, kappa_T,LT_center, LT_width, &
            & kappa_n, Ln_center, Ln_width, delta_x_T, delta_x_n,&
            & prof_type, src_amp, src_width, src_x0, src_prof_type, &
            & prof_file)
       !!       if ((mype==0).and.(print_ini_msg)) call print(spec(n))

      end do

      read(PARFILE, nml=units,iostat=ierr)
      if (ierr.ne.0.) ierr=0

      if (PARFILE.ne.5) close(PARFILE)

      if ((magn_geometry.eq.'tracer').and.multiple_tracer_files) &
          WRITE(geomfile,"(2A,I2.2)") TRIM(geomfile),"_",my_sim
  
      call set_in_prof_defaults()

      call check_params(par_in_dir)
      
end subroutine read_tga_parameters

  subroutine check_params(par_in_dir)
    character(len=*), intent(in) :: par_in_dir
    logical :: valid_parall, file_exists
    integer :: fileunit

    if (.not.y_local) then
       if (.not.x_local) stop 'x and y global is not implememted yet'
    end if


    IF (pbc_dealiasing) THEN
       CALL set_pbc_dealiasing(.true.)
    END IF
    
    evenx = mod(nx0+1,2)

    if (nonlinear) then 
       ! ky dependent lx and ky shifts don't make sense 
       ! in nonlinear simulations.
       adapt_lx=.false.
       ky0_ind=0
       kx_center=0.
    elseif((nky0.eq.1).and.((ky0_ind.eq.0).or.(abs(kymin).lt.1e-5)).and. &
         &(n0_global.eq.0.or.n0_global.eq.-1111)) then
       ! for zonal mode computations adapt_lx doesn't make sense either.
       only_zonal=.true.
       adapt_lx=.false.
       if ((abs(lx).lt.epsilon(lx)).and.(abs(kx_center).lt.epsilon(kx_center))) then
          if(mype.le.0) then 
             print*,'' 
             print*, 'You have to choose a finite lx or kx_center! Stopping'
          endif
          stop
       endif
       ky0_ind=0
    endif

    if (lilo) then
       adapt_lx = .false.
       x_local=.false.
    endif


    ! ky dependent lx make no sense for radially nonlocal simulations 
    ! with q profiles (not quantization cond. in parallel b.c.)
    if (mag_prof) adapt_lx = .false.
    !adapt_lx has to be set before check_diag

    if (nonlinear.and.(.not.turbdeal)) yx_order=.false.
    if (.not.x_local) yx_order=.false.
    if ((mype.le.0).and.yx_order) print*, 'USING YX_ORDER!!!'
    xy_local = (x_local.and.y_local)

    !call check_par_checkpoint
    !call check_diag
    !call check_convergence_monitoring
    !call check_init_cond
    !call check_external_contr
    !call check_sources

    if (.not.y_local) then
       rad_bc_type=0
       ky0_ind = 0
       yx_order=.true.
       endif


    if (diagdir.eq.'') stop "diagdir not specified"

    !remove leading and trailing spaces
    diagdir=trim(adjustl(diagdir))

            !check for tracer file
    if (magn_geometry.eq.'tracer') then
       inquire(file=trim(geomdir)//'/'//trim(geomfile),exist=file_exists)
       if (file_exists) then
          fileunit = -1
          if (trim(par_in_dir).ne.'skip_parfile') &
               & call read_tracer_namelist(fileunit,.true.)
       else
          write(*,"(4A)") trim(geomdir),'/',trim(geomfile),' does not exist'
          stop
       endif
    endif

    if (.not.(auto_parall)) valid_parall = check_parallelization(.true.)
    call check_discretization    
    if ((.not.adapt_lx).and.(nexc.eq.0).and.(lx.eq.0.0)) Stop &
         "neither nexc nor lx found!"
    if ((nexc.eq.1).and.(.not.nonlinear).and.(.not.adapt_lx)) &
         & adapt_lx = .true.
    !call check_par_phys

    if((.not.write_h5).and.(.not.write_std)) then
       write(*,"(A)") "standard or hdf5 output need to be chosen"
    stop
    end if


 
  end subroutine check_params


	


   END MODULE tga_io
  	 
