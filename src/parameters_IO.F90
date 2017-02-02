#include "redef.h"
#include "switches.h"
#include "intrinsic_sizes.h"
!>Routines to read and write the GENE paramters file
!!\todo remove read_tracer_par
module parameters_IO
use par_in
use par_other
use par_geom, only: Bprof_coeffs
use file_io, only: get_unit_nr
use discretization
use discretization_adptv_module
use coordinates
use check_parameters
use chease_mod
use collisions
use geometry
use perf_opt
use initcond
!use initial_value_comp, only: rescale
use equilibrium_fields
use external_contr
use antenna
use diagnostics
use diagnostics_df
use diagnostics_fsa_moments
use diag_nl_eigenvalues
use checkpoint
use reset_mod
use communications,only: omp_num_threads
use sources_mod
use GaussQuadrature, only: mu_grid_type
USE boundary_exchange_z
use fullmatrix_aux, only: mm_format
use eigen_parameters
use petsc_precond
USE gyro_average_df_mod,only: gyroav_in_xi_eta
use convergence_monitoring
use neo_equilibrium
use time_averages
use profiles_mod, only: Tref, nref, mref, omegatorref, drive_buffer, set_prof_defaults, &
     &set_in_prof_defaults, iterdb_file, iterdb_time, norm_index, zeff, zeff0
use dzv_terms
use numerical_damping
use quasilinear_model
use Gyro_LES
use diag_Gyro_LES
use miller_mod
use prefactors, only: with_coriolis, with_centrifugal, with_bxphi0, set_prefactors_defaults
use profile_smoothing
#ifdef WITHFUTILS
use futils
#endif

implicit none
public:: read_parameters, write_parameters, read_diagdir, par_out_file, PARFILE
public:: check_params

!private
!!the private declaration is a good idea, however, public access to 
!!many parameters is required in interfaces, e.g. ITM

integer:: PAR_OUT_FILE, PARFILE=30

integer:: f_ext_startnum=1
integer:: j

namelist /in_out/ &
     diagdir, chptdir, read_checkpoint, write_checkpoint, many_chpts, &
     istep_schpt, istep_field, istep_mom, istep_nrg, istep_neoclass, istep_neoclass2D, & 
     istep_vsp, istep_omega, istep_g1, istep_prof, write_fielddiff, mm_format,&
     istep_energy, istep_energy3d, istep_dfout, write_geom, &
     istep_srcmom, write_h5, write_std, chpt_h5, chpt_read_h5, chpt_write_h5,&
     chpt_hac, chpt_read_hac, chpt_write_hac, &
     istep_nlt, istep_gav, write_flux_final, istep_GyroLES,istep_fe_time,&
     istep_fe_twoD, istep_fe_transfer, avg_window, iterdb_file, iterdb_time, &
     momentum_flux, istep_fsa_moments, cat_output, istep_nlev, istep_antenna, &
     quasilin_model, istep_energy_exchange

contains
  !>Reads parameters from ./parameter.dat file
  subroutine read_parameters(par_in_dir)
    character(len=*), intent(in) :: par_in_dir
    character(len=10) :: name
    character(len=FILENAME_MAX) :: prof_file, par_file
    logical:: passive
    real:: omn, omt, mass, temp, dens, & 
         & kappa_T, LT_center, LT_width, & 
         & kappa_n, Ln_center, Ln_width, &
         & delta_x_T, delta_x_n
    real, Dimension(0:4) :: src_amp, src_width, src_x0, &
         & src_prof_type
    integer:: charge, n, ierr, prof_type
    real, dimension(1024) :: blk_mks_r_v, blk_mks_r_w ! big enouth arrays for r marks
    real, dimension(1024) :: blk_mks_v, blk_mks_w ! big enouth arrays for v||, w marks
    real :: vtorref !old name for omegatorref (just kept for 
                    !backward compatbility and shall be removed at some point)

    namelist /box/ &
         nx0, nky0, nz0, nv0, nw0, n_spec, &
         lx, lx_a, kymin, lv, lw, ky0_ind, kx_center, nexc, &
         adapt_lx, adapt_ly, x0,&
         mu_grid_type,&
         n0_global  ! for parallel boundary condition


    namelist /general/ &
         nonlinear, parallel_nl, x_local, y_local, lilo, comp_type, &       !operation type
         perf_vec, nblocks, parscheme, vparscheme, turbdeal, &
         precomp_nc, include_f0_contr,& !rescale
         timescheme, calc_dt, dt_max, dt_vlasov, ev_coll, courant, coll_split, coll_split_scheme,& !time scheme
         n_ev, which_ev, ev_max_it, ev_shift, ev_n_test, ev_prec,& !Eigenvalue comp.
         vec_out_limit,ev_out,ev_left, ev_right,&
         pc_type,sub_pc_type,pc_sub_it,pc_factor_level, pc_blocks, & !precondititioner
         nc_prec, nc_max_it,& !NC equilibrium comp.
         ntimesteps, timelim, simtimelim, &  !simulation stop criteria
         overflow_limit, underflow_limit, omega_prec, &
         beta, debye2, delzonal, delzonal_factor, delzonal_fields, & 
         add_zonal_phi, &!physical parameters + related
         bpar, bpar_off, pressure_term, pressure_off, &
         antenna_type, Apar0_antenna, omega0_antenna, & !constant amplitude antenna
         lv_antenna_modes, lv_antenna_amp, lv_antenna_initamp, lv_antenna_freq, & !Langevin antenna
         del_phi, del_fields, only_Er, trap_pass, tau, no_electron_response,&  
         coll, collision_op, coll_cons_model, spacediff, spacediff_off, coll_f_fm_on, Zeff, adiabatic_mass,&
         init_cond, init_aux_x, init_aux_y, init_aux_z, init_aux_amp, init_vpar_fluct,&  !initial condition
         hyp_x, hyp_x_order, hyp_y, hyp_y_order, &         !hyper diffusion
         hyp_perp, hyp_perp_order, &
         hyp_z, hyp_z_order, hyp_v, hyp_on_h, hypz_opt,&
         GyroLES, fracx,fracy, diag_GyroLES,tky,tkx,tkperp,n_shells_kperp,t_log_kperp_fb,t_log_kperp_h,& ! GLES related
         tk3_log_kperp, tk3_ky, tk3_kx, t_log_kperp, t_log_ky, n_shells_log_kperp, n_shells_log_ky, & ! GLES related
         tk3_log_kperp_1, tk3_log_kperp_h,t_log_kperp_1, n_shells_log_kperp_1, lambda1, lambda2, with_nvp, & ! GLES related
         coeff_g1_en_source, & ! g1 energy source
         perf_tsteps, arakawa_zv, arakawa_cons_bc, arakawa_zv_order, fourier2d, yx_order, & !misc
         diag_Blev, diag_trap_levels,& !diagnostics related
         avgflux_stime, avgprof_stime, avgflux_type, avgflux_alpha, gav_stime,&
#ifdef COMBI         
         ev_shift_ref,&
#endif
         &gpu_cpu_ratio,which_func_as
 
 

    !Add additional specialized diagnostics here:
    namelist /extended_diags/ &   !energy, nlt, svd (etc.) diagnostics
         which_ky,num_ky_modes,num_kx_modes, &          !istep_dfout
         which_kx_center, dfout_mpiio, &                !istep_dfout
         SVD_proj,SVD_kx_ind,SVD_ky_ind,SVD_nkx0,SVD_df_n_time,SVD_start_time,&  !POD diagnostics
         SVD_sparse_factor,SVD_df_file_name,SVD_f0_flag, SVD_df_file_path, &      !POD diagnostics
         SVD_df_file_suffix, SVD_n_restarts, SVD_pbc_phase, SVD_parallel, &      !POD diagnostics
         SVD_field_mom, SVD_swap_endian, &      !POD diagnostics
         num_nlt_modes,kx_nlt_ind,ky_nlt_ind,nlt_symmetrize,&        !istep_nlt
         nx0_nlt_pod,num_nlt_pod_modes ,nlt_pod,&     !nlt pod
         nlev_stime, n_nlev                            !nonlinear eigenvalue diagnostics

    namelist /nonlocal_x/ &
         rad_bc_type, gyroav_in_xi_eta,&
         arakawa, shifted_metric, &
         !sources
         ck_heat, ck_part, ck_heat_smooth, ck_part_smooth, ck_filter_type, &
         l_buffer_size, lcoef_krook, lpow_krook,&
         u_buffer_size, ucoef_krook, upow_krook, buffer_on_ky0,&
         lck_buffer_size, lckcoef_krook, lckpow_krook,&
         uck_buffer_size, uckcoef_krook, uckpow_krook,&
         puls_start,puls_amp,puls_on_duration,puls_off_duration,&
         reset_limit, smooth_reset, smooth_reset_type, smooth_reset_width,&
         drive_buffer,explicit_buffer, psource_type, &
         intsource_heat_coeff, intsource_part_coeff, intsource_time,&
         ga_spatial_var, nc_exb_corr

    !variables for external (zonal) fields, temperatures and densities
    namelist /external_contr/ &
         phi0_ext, kxind_phi_ext, phase_phi_ext,&
         apar0_ext, kxind_apar_ext, phase_apar_ext,&
         omt0_ext, kxind_omt_ext, phase_omt_ext,&
         omn0_ext, kxind_omn_ext, phase_omn_ext,&
         ExBrate, pfsrate, ExB_stime, Erad, Erad_acc, &
         Omega0_tor, R0_tor, dR0_tor_dx, &
         with_coriolis, with_centrifugal, with_comoving_other, with_bxphi0, &
         lilo_w_full_omegator_profile
    
    namelist /geometry/ &
         magn_geometry, shat, limit_shat, &     !important for all geometries
         parscale, & !parscale for slab
         q0, major_R,minor_r,major_Z, trpeps,amhd, & ! rest for s-alpha/circular/miller
         geomdir, geomfile, x_def, flux_pos,&    !general geometry
         rhostar, mag_prof, q_coeffs, n_pol, edge_opt, &
         norm_flux_projection, psi_o_twopi, q_scalefac, &
         force_x, force_z, & !used only in CHEASE case
         kappa, delta, rho, s_kappa, s_delta, & !used only for Miller geometry
         drR, drZ, zeta, s_zeta, & ! used only for Miller geometry
         dpdx_term, dpdx_pm, & !neg. pressure gradient normalized to reference magn. pressure
         Bprof_coeffs, & !correction terms (e.g., RFP), used with circular
         ialpha !select flux tube in GIST only

    namelist /species/ &
         & name, passive, omn, omt, mass, charge, temp, dens, & 
         & kappa_T,LT_center,LT_width, kappa_n,Ln_center,Ln_width,&
         & delta_x_T, delta_x_n, prof_type,&
         & src_amp, src_width, src_x0, src_prof_type, prof_file

    namelist /units/ &
         & Tref, nref, Lref, Bref, mref, omegatorref, norm_index, &
         & vtorref  !vtorref is the old name for omegatorref (backward compat.)

    ! adaptivity namelist contains parameters necessery for adaptive grid definition
    namelist /adaptivity/ &
         & is_grid_adptv, & ! flag showing if the grid is adaptive
         & is_numofvpoints_const, & ! flag to set constant number of v|| points in each blog
         & n_vx_blks, & ! number of blocks for a block structured r, v|| subgrid
         & n_wx_blks, & ! number of blocks for a block structured r, mu subgrid
         & prob_vx, &  ! probability to find a "particle" in the r, v|| subgrid domain
         & prob_wx, &  ! probability to find a "particle" in the r, mu subgrid domain
         & opt_mthd_vx, & ! optimization method for r, v|| subrid construction
         & opt_mthd_wx, & ! optimization method for r, mu subrid construction
         & blk_mks_r_v, blk_mks_v, & ! block marks for the grid r, mu subgrid
         & blk_mks_r_w, blk_mks_w ! block marks for the grid r, mu subgrid

    !default values
    call set_discretization_defaults
    call set_discretization_adptv_defaults
    call set_coordinates_defaults
    call set_par_other_defaults
    call set_par_in_defaults
    call set_cm_defaults
    call set_nd_defaults
    call set_nc_defaults
    call set_ext_defaults
    call set_prof_defaults
    call set_coll_defaults
    call set_geometry_defaults
    call set_reset_mod_defaults
    call set_dzv_defaults
    call set_equilibrium_fields_defaults
    call set_antenna_defaults
    call set_prefactors_defaults

    f_ext_startnum = 1
    diag_Blev=0.
    chptdir = ''
    q_coeffs=0.0
    n_ev=0

    if(allocated(spec)) deallocate(spec)

    ! get the defaults from the module
 
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
    read(PARFILE, nml=extended_diags,iostat=ierr)
    if (ierr.gt.0) stop 'on i/o error: incorrect extended_diags namelist'

    rewind(PARFILE)
    read(PARFILE, nml=nonlocal_x,iostat=ierr)
    if (ierr.gt.0) stop 'on i/o error: incorrect nonlocal_x namelist'

    rewind(PARFILE)
    read(PARFILE, nml=external_contr,iostat=ierr)
    if (ierr.gt.0) stop 'on i/o error: incorrect external_contr namelist'

    rewind(PARFILE)
    read(PARFILE, nml=geometry,iostat=ierr)
    if (ierr.ne.0) stop 'on i/o error: incorrect geometry namelist'

    allocate(spec(0:n_spec-1))

    do n=0,n_spec-1
       !set default for prof_type
       if ((x_local).or.(lilo)) then 
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

    vtorref = -2222.0
    read(PARFILE, nml=units,iostat=ierr)
    if (ierr.ne.0.) ierr=0
    if ((omegatorref.eq.0.0).and.(vtorref.ne.-2222.0)) &
         omegatorref = vtorref

    rewind(PARFILE)
    read(PARFILE, nml = adaptivity, iostat = ierr)
    if(ierr.eq.0) then
       call extract_blk_mks(&
            & blk_mks_r_v, blk_mks_v, &
            & blk_mks_r_w, blk_mks_w)
       call set_number_of_v_points
    end if

    if (PARFILE.ne.5) close(PARFILE)

    if ((magn_geometry.eq.'tracer').and.multiple_tracer_files) &
         WRITE(geomfile,"(2A,I2.2)") TRIM(geomfile),"_",my_sim

    call set_in_prof_defaults()

    call check_extended

    call check_params(par_in_dir)

  end subroutine read_parameters

  subroutine read_diagdir
    integer:: ierr
    character(len=FILENAME_MAX) :: par_file

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
    read(PARFILE, nml=in_out,iostat=ierr)
    if (ierr.ne.0) stop 'on i/o error: incorrect in_out namelist'    
    if (PARFILE.ne.5) close(PARFILE)

  end subroutine read_diagdir

  !>Check parameters for extended diagnostics
  subroutine check_extended

    character(len=4) :: kychar,kxchar

    !Checks for nlt_pod
    if(nlt_pod) then
      if(SVD_df_file_path=='no_input') then
        SVD_df_file_path=diagdir
      end if
      if(n_procs_y.gt.1) stop "nlt_pod cannot be used with ky parallelization."
      if(.not.arakawa_zv) stop "Must use arakawa_zv=T with nlt_pod."
    end if
    !Checks for POD diagnostics
    if(SVD_field_mom) then

#ifdef WITHSCAL
#else         
      stop "Must have scalapack for parallel SVD routine."
#endif

      n_moms = 6
      SVD_parallel=.false.
      SVD_proj=.true.
      if(SVD_df_n_time(1)==-100) &
              stop "Must enter SVD_df_n_time"
      if(SVD_n_restarts.gt.1.and.trim(SVD_df_file_suffix(1))=='.dat') then
              stop "Must enter SVD_file_suffix when reading multiple restarts."
      end if
      nonlinear=.false. 
      comp_type='SV'
      istep_nlt=0
      istep_dfout=0
      istep_energy=100
      istep_energy3d=500
      istep_schpt=0
      SVD_istep_mf_factor=istep_mom/istep_field
      !Need to account for the following case:
      if(SVD_istep_mf_factor.lt.1) stop "Error for SVD_field_mom."
      istep_field=istep_mom

      ntimesteps=0

      if(SVD_df_file_path=='no_input') then
        SVD_df_file_path=diagdir
      end if
#ifdef WITHSCAL
#else         
        stop "Must have scalapack for parallel SVD routine."
#endif

    end if

    if(SVD_proj.and..not.svd_field_mom) then
      if(.not.SVD_parallel) then
      !if(n_procs_sim.gt.1) stop "POD routine must be run serially (one processor)."
        n_procs_x=1
        n_procs_y=1
        n_procs_z=1
        n_procs_v=1
        n_procs_w=1
        n_procs_s=1
        n_procs_sim=1
      else
#ifdef WITHSCAL
#else         
        stop "Must have scalapack for parallel SVD routine."
#endif
      end if
      if(SVD_nkx0==-100.or.SVD_kx_ind==-100.or.SVD_ky_ind==-100) &
              stop "Must enter SVD_kx_ind, SVD_ky_ind, and SVD_nkx0!"
      if(SVD_df_n_time(1)==-100) &
              stop "Must enter SVD_df_n_time"
      if(SVD_n_restarts.gt.1.and.trim(SVD_df_file_suffix(1))=='.dat') then
              stop "Must enter SVD_file_suffix when reading multiple restarts."
      end if

      if((-1)**(nint(SVD_ky_ind*shat*lx*kymin)).gt.0) then
        SVD_pbc_phase=.true.
        if(mype==0) write(*,*) "Applying phase shift for parallel b.c."
      else
        SVD_pbc_phase=.false.
        if(mype==0) write(*,*) "No phase shift for parallel b.c."
      end if
      kx_center=3.14159*2.0/lx*SVD_kx_ind

      nonlinear=.false. 
      comp_type='IV'
      istep_nlt=0
      istep_dfout=0
      istep_energy=100
      istep_energy3d=500
      istep_schpt=0
      nx0=SVD_nkx0
      nky0=1
      kymin=kymin*SVD_ky_ind
      adapt_lx=.true. 
      ntimesteps=0

      if(SVD_df_file_path=='no_input') then
        SVD_df_file_path=diagdir
      end if
      if(SVD_df_file_name=='no_input') then
        write(kychar,"(i4.4)") SVD_ky_ind
        write(kxchar,"(i4.4)") ABS(SVD_kx_ind)
        if(SVD_kx_ind.ge.0) then
           SVD_df_file_name='df_ky'//&
                &kychar//'kx'//kxchar
        else
           SVD_df_file_name='df_ky'//&
                &kychar//'kxn'//kxchar
        end if
        if(mype==0) write(*,*) "Reading SVD_df_file_name:",SVD_df_file_name
      end if

    end if

    !Checks for nonlinear transfer diagnostics
    if(istep_nlt.gt.0) then
      if (istep_energy.le.0) then
        if(mype==0) write(*,*) "istep_energy not set: setting istep_energy=istep_nlt"
        istep_energy=istep_nlt
      end if
      if (istep_energy3d.le.0) then
        if(mype==0) write(*,*) "istep_energy3d not set: setting istep_energy=istep_nlt"
        istep_energy3d=istep_nlt
      end if      
      if (num_nlt_modes.gt.60.or.num_nlt_modes.lt.0) stop "num_nlt_modes must be between 1 and 20 inclusive."
      if (.not.nonlinear) then
         Write(*,"(A)") "WARNING: nonlinear transfer diagnostics deactivated for linear run."
         istep_nlt=0
      endif
      if (.not.x_local)  stop "nonlinear transfer diagnostics can only be used for local runs."
    end if

  end subroutine check_extended

  !>Check parameters
  !!\todo move parameter checks to the corresponding modules
  subroutine check_params(par_in_dir)
    character(len=*), intent(in) :: par_in_dir
    character(len=FILENAME_MAX) :: test_file
    logical :: valid_parall, file_exists
    integer :: fileunit, ierr, fcount, n

    if (.not.y_local) then
       if (.not.x_local) stop 'x and y global is not implememted yet'
       if (abs(kx_center).gt.0.0) stop 'kx_center<>0 not implemented in y global yet'
    end if
    if(comp_type=='EV'.and.which_ev=='all_mpl') write_checkpoint=.false.

    evenx = mod(nx0+1,2)

    !zonal mode computation (only ky=0 in system)
    only_zonal = (nky0.eq.1).and.((ky0_ind.eq.0).or.(abs(kymin).lt.1e-5)&
         &.or.n0_global.eq.0.or.((abs(kymin).lt.1e-5).and.n0_global.eq.-1111))

    if (nonlinear) then 
       ! ky dependent lx and ky shifts don't make sense 
       ! in nonlinear simulations.
       adapt_lx=.false.
       ky0_ind=0
       kx_center=0.
    elseif (only_zonal) then
       ! for zonal mode computations adapt_lx doesn't make sense either.
       adapt_lx=.false.
       if ((abs(lx).lt.epsilon(lx)).and.(abs(kx_center).lt.epsilon(kx_center))&
            .and.(x_local.or.lx_a.lt.epsilon(lx_a))) then
          if(mype.le.0) then 
             print*,'' 
             print*, 'You have to choose a finite lx or kx_center! Stopping'
          endif
          stop
       endif
       ky0_ind=0
    endif

    !decide if we have a pure neoclassics computation
    if (xy_local) then
       only_neo = only_zonal.and.include_f0_contr.and.(nx0==0)
    elseif (.not.x_local) then
       only_neo = only_zonal.and.include_f0_contr
    elseif (.not.y_local) then
       only_neo = (del_phi.or.del_fields).and.include_f0_contr
    endif
    if (only_neo.and.mype==0) write(*,"(A)") "neoclassical computation."

    if (lilo) then
       adapt_lx = .false.
       x_local=.false.
    endif

    ! ky dependent lx make no sense for radially nonlocal simulations 
    ! with q profiles (not quantization cond. in parallel b.c.)
    if (.not.(x_local.or.lilo)) adapt_lx = .false.
    !adapt_lx has to be set before check_diag
    
    if (nonlinear.and.(.not.turbdeal)) yx_order=.false.
    if (.not.x_local) yx_order=.false.
    if ((mype.le.0).and.yx_order) print*, 'USING YX_ORDER!!!'
    xy_local = (x_local.and.y_local)

    !always adapt ly/kymin if x- or y-global
    if ((.not.xy_local.and.(.not.lilo))) adapt_ly = .true.

    call check_geometry
    call check_par_checkpoint
    call check_diag
    call check_diag_nlev
    call check_convergence_monitoring
    call check_init_cond
    call check_external_contr
    call check_sources

    if (diagdir.eq.'') stop "diagdir not specified"

    !remove leading and trailing spaces
    diagdir=trim(adjustl(diagdir))

    !set chptdir
    if (chptdir.eq.'') then
       chptdir=TRIM(diagdir)
    endif

#ifndef COMBI_MGR
    !this leads to problems because every time when the parameters file
    !is saved a / is added. after multiple times the filename gets to long.
    !a clean style would be to only do this if the string doesnt end with /.
    diagdir = TRIM(diagdir)//'/'
    chptdir = TRIM(chptdir)//'/'
#endif

    If (reset_limit.ne.1000.) then
       if (read_checkpoint) then
          do fcount=1,9999
             write(test_file, "(A,I4.4)") trim(chptdir)//'/checkpoint_',fcount 
             inquire(file= test_file, exist=file_exists, iostat = ierr)
             if (ierr .ne. 0) stop 'Error while inquiring checkpoint files from run with reset'
             if (file_exists) f_ext_startnum = fcount+1
          enddo
       endif
       write(file_extension,"(A,I4.4)") '_',f_ext_startnum
       if (f_ext_startnum.gt.1) then
          write(checkpoint_in,"(A,I4.4)") trim(chptdir)//'/checkpoint_',&
               &f_ext_startnum-1
          do n=0,n_spec-1
             spec(n)%prof_type = -1
             write(spec(n)%prof_file, "(A,I4.4)") trim(diagdir)//'/profiles_'//trim(spec(n)%name)//'_',&
               &f_ext_startnum-1
          end do
          if (mype .eq. 0) write(*,'(A)') 'Note: Set prof_type=-1 for all species'
       else
          checkpoint_in = trim(chptdir)//'/checkpoint'//trim(file_extension)
       endif
    else
       if(file_extension.eq.'.dat') then
          checkpoint_in=trim(chptdir)//'/checkpoint'
       else
          checkpoint_in=trim(chptdir)//'/checkpoint'//file_extension
       endif
    endif

    !check for tracer file
    if (magn_geometry.eq.'tracer') then
       inquire(file=trim(geomdir)//'/'//trim(geomfile),exist=file_exists)
       if (file_exists) then
          fileunit = -1
          if (trim(par_in_dir).ne.'skip_parfile') &
               & call read_tracer_namelist(fileunit,.true.,.false.)
       else
          write(*,"(4A)") trim(geomdir),'/',trim(geomfile),' does not exist'
          stop
       endif
    endif

    if (.not.(auto_parall)) valid_parall = check_parallelization(.true.)

    call check_discretization    
    if (x_local.and.(.not.adapt_lx).and.(nexc.eq.0).and.(lx.eq.0.0)) Stop &
         "neither nexc nor lx found!"
    if (.not.x_local.and.lx.eq.0.0.and.lx_a.eq.0.0) stop &
         "Please specify radial box size using either lx or lx_a!"
    if ((nexc.eq.1).and.(.not.nonlinear).and.(.not.adapt_lx)) &
         & adapt_lx = .true.
    call check_par_phys

    if((.not.write_h5).and.(.not.write_std)) then
       write(*,"(A)") "standard or hdf5 output need to be chosen"
    stop
    end if

    if (parallel_nl) then
       if (.not.nonlinear) stop "parallel nonlinearity requires ExB nonlinearity to be switched on"
       if (rhostar.le.0) stop "finite rhostar required for parallel nonlinearity"
       if (minor_r.le.0) stop "define minor_r for parallel nonlinearity"
       if (yx_order) stop "yx_order not implemented for parallel nonlinearity"
       if (beta.gt.0) stop "only electrostatic version of parallel nonlinearity implemented"
    endif

  end subroutine check_params

  !*************************************************************************!
  !******************* write parameter file ********************************!
  !*************************************************************************!

  !>Writes parameters to diagdir/parameters.dat
  Subroutine write_parameters
    INTEGER :: n
    CHARACTER(len=500) :: tmpstr

#ifdef WITHFUTILS
    integer :: fidall_param_h5
    character(len=FILENAME_MAX) :: fileall_param_h5
#endif

    ! Local variables (all used variables, which should appear in the 
    ! parameters.dat file). They are all read from the modules via get_ routines.


    if(.not.cat_output) then
       call get_unit_nr(PAR_OUT_FILE)
       OPEN(PAR_OUT_FILE,file=TRIM(diagdir)//'/parameters'//&
            &trim(file_extension),Action='write')
    end if

    !****** parallelization namelist ********
    Write(PAR_OUT_FILE,"(A)")    "&parallelization"

    Write(PAR_OUT_FILE,"(A,I3)") "n_procs_s = ",n_procs_s
    Write(PAR_OUT_FILE,"(A,I3)") "n_procs_v = ",n_procs_v
    Write(PAR_OUT_FILE,"(A,I3)") "n_procs_w = ",n_procs_w
    Write(PAR_OUT_FILE,"(A,I3)") "n_procs_x = ",n_procs_x
    Write(PAR_OUT_FILE,"(A,I3)") "n_procs_y = ",n_procs_y
    Write(PAR_OUT_FILE,"(A,I3)") "n_procs_z = ",n_procs_z
    Write(PAR_OUT_FILE,"(A,I6)") "n_procs_sim = ",n_procs_sim
    if(n_parallel_sims.gt.1) &
         &Write(PAR_OUT_FILE,"(A,I6)") "n_parallel_sims = ",n_parallel_sims
    Write(PAR_OUT_FILE,"(A)")    "/"
    Write(PAR_OUT_FILE,"(A)")    ""

    !****** box namelist ********
    Write(PAR_OUT_FILE,"(A)")    "&box"

    Write(PAR_OUT_FILE,"(A,I4)") "n_spec = ",n_spec    
    Write(PAR_OUT_FILE,"(A,I4)") "nx0    = ",nx0
    Write(PAR_OUT_FILE,"(A,I4)") "nky0   = ",nky0
    Write(PAR_OUT_FILE,"(A,I4)") "nz0    = ",nz0
    Write(PAR_OUT_FILE,"(A,I4)") "nv0    = ",nv0
    Write(PAR_OUT_FILE,"(A,I4)") "nw0    = ",nw0    
    Write(PAR_OUT_FILE,"(A)")    ""

    Write(PAR_OUT_FILE,"(A,G12.4)") "kymin = ",kymin
    Write(PAR_OUT_FILE,"(A,G12.3)") "lv    = ",lv
    Write(PAR_OUT_FILE,"(A,G12.3)") "lw    = ",lw

    if (x_local) then
       if (.not.adapt_lx) then
          Write(PAR_OUT_FILE,"(A,G13.6)") "lx    = ",lx
          Write(PAR_OUT_FILE,"(A,I3)") "nexc  = ",nexc
       endif
       Write(PAR_OUT_FILE,"(A,L1)") "adapt_lx = ", adapt_lx
       if (magn_geometry .eq. 'tracer_efit'.or.&
            magn_geometry.eq.'miller'.or.&
            magn_geometry.eq.'gist'.or.&
            magn_geometry.eq.'miller_b'.or.mag_prof) &
            &Write(PAR_OUT_FILE,"(A,G12.4)") "x0    = ",x0 
    else
       if (lx.ne.0.0) Write(PAR_OUT_FILE,"(A,G13.6)") "lx    = ",lx
       if (lx_a.ne.0.0) Write(PAR_OUT_FILE,"(A,G13.6)") "lx_a    = ",lx_a
       Write(PAR_OUT_FILE,"(A,G12.4)") "x0    = ",x0 
    endif

    if (n0_global.ne.-1111) Write(PAR_OUT_FILE,"(A,I6)") "n0_global = ", n0_global
    IF ((xy_local.or.lilo).and.adapt_ly) &
         Write(PAR_OUT_FILE,"(A,L1)") "adapt_ly = ", adapt_ly

    if (y_local.and.(.not.nonlinear)) Write(PAR_OUT_FILE,"(A,I3)") "ky0_ind = ",ky0_ind
    if (kx_center.ne.0.) Write(PAR_OUT_FILE,"(A,G12.4)") "kx_center = ",kx_center
    WRITE(PAR_OUT_FILE,"(3A)")  "mu_grid_type = '", TRIM(mu_grid_type),"'"

    Write(PAR_OUT_FILE,"(A)")    "/"
    Write(PAR_OUT_FILE,"(A)")    ""

    !****** in_out namelist ********
    Write(PAR_OUT_FILE,"(A)")    "&in_out"

    WRITE(PAR_OUT_FILE,"(3A)")   "diagdir = '", TRIM(diagdir),"'"
    IF (TRIM(chptdir).ne.TRIM(diagdir)) &
         &WRITE(PAR_OUT_FILE,"(3A)")   "chptdir = '", TRIM(chptdir),"'"
    Write(PAR_OUT_FILE,"(A)")    ""

    WRITE(PAR_OUT_FILE,"(A,L1)") "read_checkpoint  = ", read_checkpoint
    WRITE(PAR_OUT_FILE,"(A,L1)") "write_checkpoint = ", write_checkpoint
    if (many_chpts) WRITE(PAR_OUT_FILE,"(A)") "many_chpts = T"
    Write(PAR_OUT_FILE,"(A)")    ""

    WRITE(PAR_OUT_FILE,"(A,I7)") "istep_field  = ",istep_field
    WRITE(PAR_OUT_FILE,"(A,I7)") "istep_mom    = ",istep_mom
    WRITE(PAR_OUT_FILE,"(A,I7)") "istep_nrg    = ",istep_nrg
    if (.not.nonlinear) WRITE(PAR_OUT_FILE,"(A,I7)") "istep_omega  = ",istep_omega
    WRITE(PAR_OUT_FILE,"(A,I7)") "istep_vsp    = ",istep_vsp
    Write(PAR_OUT_FILE,"(A,I7)") "istep_schpt  = ",istep_schpt

    IF (istep_g1.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_g1     = ",istep_g1
    IF (istep_gav.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_gav    = ",istep_gav

    IF (istep_energy.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_energy = ",istep_energy
    IF (istep_energy3d.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_energy3d = ",istep_energy3d
    IF (istep_energy_exchange.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_energy_exchange = ",istep_energy_exchange
    IF (istep_nlt.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_nlt     = ",istep_nlt
    IF (istep_nlt.gt.0) Write(PAR_OUT_FILE,"(A,L1)") "nlt_symmetrize    = ",nlt_symmetrize
    IF (istep_dfout.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_dfout  = ",istep_dfout
    IF (istep_prof.gt.0) Write(PAR_OUT_FILE,"(A,I7)")  "istep_prof   = ",istep_prof
    IF ((istep_field .GT. 0).and.(write_fielddiff)) &
         &WRITE(PAR_OUT_FILE,"(A,L1)") "write_fielddiff = ", write_fielddiff
    IF (istep_srcmom.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_srcmom   = ",istep_srcmom
    IF (istep_neoclass.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_neoclass = ",istep_neoclass
    IF (istep_neoclass2D.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_neoclass2D = ",istep_neoclass2D 
    IF (istep_fsa_moments.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_fsa_moments = ",istep_fsa_moments
    IF (GyroLES) Write(PAR_OUT_FILE,"(A,I7)") "istep_GyroLES  =",istep_GyroLES
    
    !diag_GyroLES
    IF (istep_fe_time.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_fe_time = ", istep_fe_time
    IF (istep_fe_twoD.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_fe_twoD = ", istep_fe_twoD
    IF (istep_fe_transfer.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_fe_transfer = ", istep_fe_transfer
    
    IF (istep_nlev.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_nlev   = ", istep_nlev
    IF (istep_antenna.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "istep_antenna   = ", istep_antenna

    Write(PAR_OUT_FILE,"(A)")    ""

    WRITE(PAR_OUT_FILE,"(A,L1)") "write_std = ", write_std
    IF (write_h5) WRITE(PAR_OUT_FILE,"(A,L1)") "write_h5  = ", write_h5
    IF (chpt_h5) THEN
       WRITE(PAR_OUT_FILE,"(A,L1)") "chpt_h5   = ", chpt_h5
    ELSE 
       IF (chpt_read_h5) &
            &WRITE(PAR_OUT_FILE,"(A,L1)") "chpt_read_h5   = ", chpt_read_h5
       IF (chpt_write_h5) &
            &WRITE(PAR_OUT_FILE,"(A,L1)") "chpt_write_h5   = ", chpt_write_h5
    ENDIF
    IF (chpt_hac) THEN
       WRITE(PAR_OUT_FILE,"(A,L1)") "chpt_hac  = ", chpt_hac
    ELSE 
       IF (chpt_read_hac) &
            &WRITE(PAR_OUT_FILE,"(A,L1)") "chpt_read_hac  = ", chpt_read_hac
       IF (chpt_write_hac) &
            &WRITE(PAR_OUT_FILE,"(A,L1)") "chpt_write_hac  = ", chpt_write_hac
    ENDIF
    IF (.not.momentum_flux) &
         &WRITE(PAR_OUT_FILE,"(A,L1)") "momentum_flux  = ", momentum_flux
    IF (write_flux_final.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "write_flux_final = ",write_flux_final
    IF (trim(iterdb_file).ne.'') WRITE(PAR_OUT_FILE,"(3A)")   "iterdb_file     = '", TRIM(iterdb_file),"'"
    IF (iterdb_time.ge.0.0) WRITE(PAR_OUT_FILE,"(A,F5.2)") &
         "iterdb_time     = ", iterdb_time

    Write(PAR_OUT_FILE,"(A)")    "/"
    Write(PAR_OUT_FILE,"(A)")    ""

    !****** general namelist ********
    Write(PAR_OUT_FILE,"(A)")    "&general"

    !operation type
    Write(PAR_OUT_FILE,"(A,L3)") "nonlinear = ",nonlinear
    if (parallel_nl) WRITE(PAR_OUT_FILE,"(A,L1)")  "parallel_nl = ", parallel_nl

    if (.not.x_local) WRITE(PAR_OUT_FILE,"(A,L1)") "x_local   = ", x_local
    if (.not.y_local) WRITE(PAR_OUT_FILE,"(A,L1)") "y_local   = ", y_local
    if (y_local.and.yx_order) WRITE(PAR_OUT_FILE,"(A,L1)") "yx_order   = ", yx_order
    if (lilo) WRITE(PAR_OUT_FILE,"(A,L1)") "lilo      = ", lilo
    Write(PAR_OUT_FILE,"(3A)") "comp_type = '",trim(comp_type),"'"
    Write(PAR_OUT_FILE,"(A,9I2)") "perf_vec  = ",perf_vec
    Write(PAR_OUT_FILE,"(A,I7)") "nblocks   = ",nblocks
    if (parscheme.ne.'c4th') &
         &Write(PAR_OUT_FILE,"(3A)") "parscheme = '",Trim(parscheme),"'"
    if (vparscheme.ne.'c4th') &
         &Write(PAR_OUT_FILE,"(3A)") "vparscheme = '",Trim(vparscheme),"'"    
    if (.not.arakawa_zv) then
       write(PAR_OUT_FILE,"(A,L3)") "arakawa_zv = ",arakawa_zv
    else
       if (arakawa_zv_order.ne.4) &
            &write(PAR_OUT_FILE,"(A,I3)") "arakawa_zv_order = ",arakawa_zv_order
       if (arakawa_cons_bc) Write(PAR_OUT_FILE,"(A,L3)") "arakawa_cons_bc = ",arakawa_cons_bc
       if (.not.hyp_on_h) write(PAR_OUT_FILE,"(A,L3)") "hyp_on_h = ",hyp_on_h
       if (hyp_on_h.and..not.hypz_opt) write(PAR_OUT_FILE,"(A,L3)") "hypz_opt = ",hypz_opt
    endif
    if (turbdeal) WRITE(PAR_OUT_FILE,"(A,L1)") "turbdeal  = ", turbdeal
    Write(PAR_OUT_FILE,"(A)") ""

    !time scheme & eigenvalue comp.
    select case(comp_type)
    case('EV')
       Write(PAR_OUT_FILE,"(3A)") "which_ev  = '",trim(which_ev),"'"
       If (n_ev.ne.1.and.(which_ev=='jd_as'.or.which_ev=='kx_as'))&
            &WRITE(PAR_OUT_FILE,"(A,I10)") "which_func_as      = ",which_func_as       
       If (pc_type.ne.'none') Write(PAR_OUT_FILE,"(3A)") "pc_type  = '",trim(pc_type),"'"
       if (vec_out_limit) Write(PAR_OUT_FILE,"(A)")  "vec_out_limit = T"
       If (n_ev.ne.1) WRITE(PAR_OUT_FILE,"(A,I10)") "n_ev      = ",n_ev       
       If (ev_prec.ne.1e-4) Write(PAR_OUT_FILE,"(A,G12.4)") "ev_prec = ",ev_prec
       If (ev_shift.ne.(10.,0.)) WRITE(PAR_OUT_FILE,"(2(A,F5.2),A)") &
            &"ev_shift   = (",real(ev_shift),",",aimag(ev_shift),")"
       If (ev_max_it.gt.0) Write(PAR_OUT_FILE,"(A,I7)") "ev_max_it = ",ev_max_it
       If (ev_n_test.ne.(n_ev+15)) Write(PAR_OUT_FILE,"(A,I7)") "ev_n_test = ",ev_n_test
    case('IV')
       Write(PAR_OUT_FILE,"(3A)") "timescheme = '", Trim(timescheme),"'"
       if (collision_op.ne.'none') then
          if (.not.coll_split) then
             write(PAR_OUT_FILE,"(A,L1)") "coll_split = ", coll_split
          else
             if (coll_split_scheme .ne. 'RKCa') &
                  &write(PAR_OUT_FILE,"(3A)") "coll_split_scheme  = '", trim(coll_split_scheme),"'"
          endif
       endif
       if (coll_split) then
          Write(PAR_OUT_FILE,"(A,G12.4)") "dt_max     = ",dt_max
          Write(PAR_OUT_FILE,"(A,G12.4)") "dt_vlasov  = ",dt_vlasov
          Write(PAR_OUT_FILE,"(A,G13.5)") "ev_coll    = ",ev_coll
       else
          Write(PAR_OUT_FILE,"(A,G12.4)") "dt_max     = ",dt_max
       endif
       If (nonlinear) Write(PAR_OUT_FILE,"(A,F8.2)")  "courant    = ", courant
       if (precomp_nc) Write(PAR_OUT_FILE,"(A,L1)")  "precomp_nc = ", precomp_nc
       if(include_f0_contr) Write(PAR_OUT_FILE,"(A)")  "include_f0_contr = .true."
    case('NC')
       If (nc_prec.ne.5.e-7) Write(PAR_OUT_FILE,"(A,G12.4)") "nc_prec = ",nc_prec
       If (nc_max_it.ne.100000) Write(PAR_OUT_FILE,"(A,I7)") "nc_max_it = ",nc_max_it
       if(include_f0_contr) Write(PAR_OUT_FILE,"(A)")  "include_f0_contr = .true."
    end select
    Write(PAR_OUT_FILE,"(A)") ""    

    !simulation stop criteria
    WRITE(PAR_OUT_FILE,"(A,I10)") "timelim    = ",timelim 
    If(comp_type.eq.'IV') then
       Write(PAR_OUT_FILE,"(A,I10)") "ntimesteps = ",ntimesteps
       IF (simtimelim.ne.500000.) WRITE(PAR_OUT_FILE,"(A,G12.4)") &
            &"simtimelim = ",simtimelim
    endif
    IF (overflow_limit.ne.1e30) WRITE(PAR_OUT_FILE,"(A,G12.4)") &
         &"overflow_limit = ",overflow_limit
    IF (underflow_limit.ne.1e-8) WRITE(PAR_OUT_FILE,"(A,G12.4)") &
         &"underflow_limit = ",underflow_limit
    IF ((istep_omega.gt.0).and.(omega_prec.ne.1e-3)) &
         WRITE(PAR_OUT_FILE,"(A,G12.4)") "omega_prec = ",omega_prec
    Write(PAR_OUT_FILE,"(A)") ""

    !physical (+related) parameters
    Write(PAR_OUT_FILE,"(A,G16.8)") "beta       = ",beta
    Write(PAR_OUT_FILE,"(A,G16.8)") "debye2     = ",debye2   
    IF ((beta .GT. 0.0).and.(bpar)) &
         WRITE(PAR_OUT_FILE,"(A,L3)") "bpar   = ", bpar
!    if ((beta.gt.0.0).and.(pressure_term)) &
!         &Write(PAR_OUT_FILE,"(A,L3)") "pressure_term   = ",pressure_term
    !Apar antenna parameters
    if (antenna_type.ne.0) then
       WRITE(PAR_OUT_FILE,"(A,I3)") "antenna_type    = ",antenna_type
       select case(antenna_type)
       case(1)
          WRITE(PAR_OUT_FILE,"(A,G16.8)")  "Apar0_antenna  = ", Apar0_antenna
          WRITE(PAR_OUT_FILE,"(A,G16.8)")  "omega0_antenna = ", omega0_antenna
       case(2)
          call variable_cplxarr2strarr(lv_antenna_initamp,tmpstr)
          write(PAR_OUT_FILE,"(2A)") "lv_antenna_initamp = ", trim(tmpstr)
          call variable_cplxarr2strarr(lv_antenna_amp,tmpstr)
          write(PAR_OUT_FILE,"(2A)") "lv_antenna_amp = ", trim(tmpstr)
          call variable_cplxarr2strarr(lv_antenna_freq,tmpstr)
          write(PAR_OUT_FILE,"(2A)") "lv_antenna_freq = ", trim(tmpstr)
          call variable_intarr2strarr(lv_antenna_modes,tmpstr)
          write(PAR_OUT_FILE,"(2A)") "lv_antenna_modes = ", trim(tmpstr)
       case(3)
          call variable_cplxarr2strarr(lv_antenna_initamp,tmpstr)
          write(PAR_OUT_FILE,"(2A)") "lv_antenna_initamp = ", trim(tmpstr)
          call variable_cplxarr2strarr(lv_antenna_amp,tmpstr)
          write(PAR_OUT_FILE,"(2A)") "lv_antenna_amp = ", trim(tmpstr)
          call variable_cplxarr2strarr(lv_antenna_freq,tmpstr)
          write(PAR_OUT_FILE,"(2A)") "lv_antenna_freq = ", trim(tmpstr)
          call variable_intarr2strarr(lv_antenna_modes,tmpstr)
          write(PAR_OUT_FILE,"(2A)") "lv_antenna_modes = ", trim(tmpstr)
       end select
    endif
    if (tau.ne.1.) write(PAR_OUT_FILE,"(A,G16.8)") "tau = ",tau
    if (delzonal) Write(PAR_OUT_FILE,"(A,L1)") "delzonal   = ",delzonal
    if (delzonal_fields) Write(PAR_OUT_FILE,"(A,L1)") "delzonal_fields = ",delzonal_fields
    if (delzonal_factor.ne.1.0) Write(PAR_OUT_FILE,"(A,G16.8)") "delzonal_factor = ",delzonal_factor
    if (del_phi) Write(PAR_OUT_FILE,"(A,L1)") "del_phi = ",del_phi
    if (del_fields) Write(PAR_OUT_FILE,"(A,L1)") "del_fields = ",del_fields
    if (only_Er) Write(PAR_OUT_FILE,"(A,L1)") "only_Er = ",only_Er
    if (abs(add_zonal_phi).gt.epsilon(0.0)) &
         &Write(PAR_OUT_FILE,"(A,G16.8)") "ad_dzonal_phi = ", add_zonal_phi
    write(PAR_OUT_FILE,"(3A)") "collision_op = '",Trim(collision_op),"'"
    if (collision_op.ne.'none') then 
       Write(PAR_OUT_FILE,"(A,G16.8)") "coll         = ", coll
       Write(PAR_OUT_FILE,"(3A)") "coll_cons_model  = '",trim(coll_cons_model),"'"
       if (spacediff) Write(PAR_OUT_FILE,"(A,L1)") "spacediff  = ", spacediff
       if (coll_f_fm_on) Write(PAR_OUT_FILE,"(A,L1)") "coll_f_fm_on  = ", coll_f_fm_on
       !zeff0 is written(center value of a profile (zeff<0) or input zeff>0)
       if (zeff0.ne.1.0) write(PAR_OUT_FILE,"(A,G16.8)") "Zeff = ",zeff0
       if (adiabatic_mass.ne.2.0) write(PAR_OUT_FILE,"(A,G16.8)") "adiabatic_mass = ",adiabatic_mass
    endif
    Write(PAR_OUT_FILE,"(A)") ""

    !initial condition
    WRITE(PAR_OUT_FILE,"(3A)") "init_cond = '", TRIM(init_cond),"'"
    IF (init_aux_x.NE.-100.) WRITE(PAR_OUT_FILE,"(A,G16.8)") "init_aux_x = ",init_aux_x
    IF (init_aux_y.NE.-100.) WRITE(PAR_OUT_FILE,"(A,G16.8)") "init_aux_y = ",init_aux_y
    IF (init_aux_z.NE.-100.) WRITE(PAR_OUT_FILE,"(A,G16.8)") "init_aux_z = ",init_aux_z
    IF (init_aux_amp.NE.-100.) WRITE(PAR_OUT_FILE,"(A,G16.8)") "init_aux_amp = ",init_aux_amp
    Write(PAR_OUT_FILE,"(A)") ""

    !hyper diffusion
    if(hyp_x.gt.0.0) then
       Write(PAR_OUT_FILE,"(A,G12.4)") "hyp_x = ",hyp_x
       IF (hyp_x_order.ne.4) & 
            &Write(PAR_OUT_FILE,"(A,I2)") "hyp_x_order = ",hyp_x_order
    endif
    if(hyp_y.gt.0.0) then
       Write(PAR_OUT_FILE,"(A,G12.4)") "hyp_y = ",hyp_y
       IF (hyp_y_order.ne.4) &
            &Write(PAR_OUT_FILE,"(A,I2)") "hyp_y_order = ",hyp_y_order
    endif
    if(hyp_perp.gt.0.0) then
       Write(PAR_OUT_FILE,"(A,G12.4)") "hyp_perp = ",hyp_perp
       IF (hyp_perp_order.ne.4) &
            &Write(PAR_OUT_FILE,"(A,I2)") "hyp_perp_order = ",hyp_perp_order
    endif
    if(hyp_z.ne.0.0) then
       Write(PAR_OUT_FILE,"(A,G12.4)") "hyp_z = ",hyp_z
       IF (hyp_z_order.ne.4) Write(PAR_OUT_FILE,"(A,I2)") &
            &"hyp_z_order = ",hyp_z_order
    endif
    if(hyp_v.gt.0.0) Write(PAR_OUT_FILE,"(A,G12.4)") "hyp_v = ",hyp_v
    Write(PAR_OUT_FILE,"(A)") ""

    !GyroLES
    if (GyroLES) then
       Write(PAR_OUT_FILE,"(A)") ""
       Write(PAR_OUT_FILE,"(A,L3)") "GyroLES       = ",GyroLES
       Write(PAR_OUT_FILE,"(A)") ""
    endif

    If (diag_GyroLES) then
       Write(PAR_OUT_FILE,"(A)") ""
       Write(PAR_OUT_FILE,"(A,L3)") "diag_GyroLES  = ",diag_GyroLES
       Write(PAR_OUT_FILE,"(A,G12.4)") "fracx         = ",fracx
       Write(PAR_OUT_FILE,"(A,G12.4)") "fracy         = ",fracy
       if (istep_fe_transfer.gt.0) then
           Write(PAR_OUT_FILE,"(A,L3)") "tkx           = ",tkx
           Write(PAR_OUT_FILE,"(A,L3)") "tky           = ",tky
           Write(PAR_OUT_FILE,"(A,L3)") "tkperp        = ",tkperp
           Write(PAR_OUT_FILE,"(A,L3)") "t_log_ky      = ",t_log_ky
           Write(PAR_OUT_FILE,"(A,L3)") "t_log_kperp   = ",t_log_kperp
           Write(PAR_OUT_FILE,"(A,L3)") "t_log_kperp_fb   = ",t_log_kperp_fb
           Write(PAR_OUT_FILE,"(A,L3)") "t_log_kperp_h   = ",t_log_kperp_h
           Write(PAR_OUT_FILE,"(A,L3)") "t_log_kperp_1 = ",t_log_kperp_1
           Write(PAR_OUT_FILE,"(A,L3)") "tk3_kx        = ",tk3_kx
           Write(PAR_OUT_FILE,"(A,L3)") "tk3_ky        = ",tk3_ky
           Write(PAR_OUT_FILE,"(A,L3)") "tk3_log_kperp = ",tk3_log_kperp
           Write(PAR_OUT_FILE,"(A,L3)") "tk3_log_kperp_1 = ",tk3_log_kperp_1
           Write(PAR_OUT_FILE,"(A,L3)") "tk3_log_kperp_h = ",tk3_log_kperp_h
           Write(PAR_OUT_FILE,"(A,I7)") "n_shells_kperp     = ",n_shells_kperp
           Write(PAR_OUT_FILE,"(A,I7)") "n_shells_log_ky    = ",n_shells_log_ky
           Write(PAR_OUT_FILE,"(A,I7)") "n_shells_log_kperp = ",n_shells_log_kperp
           Write(PAR_OUT_FILE,"(A,I7)") "n_shells_log_kperp_1 = ",n_shells_log_kperp_1
           Write(PAR_OUT_FILE,"(A,G12.4)") "lambda1       = ",lambda1
           Write(PAR_OUT_FILE,"(A,G12.4)") "lambda2       = ",lambda2
           Write(PAR_OUT_FILE,"(A)") ""
       endif
    endif

    IF (coeff_g1_en_source .NE. 0.0) WRITE(PAR_OUT_FILE,"(A,G12.4)") &
      "coeff_g1_en_source = ", coeff_g1_en_source

    !Misc.
    if (perf_tsteps.ne.3) then
       Write(PAR_OUT_FILE,"(A,I3)") "perf_tsteps = ", perf_tsteps
       Write(PAR_OUT_FILE,"(A)") "" 
    endif

    !diagnostics related switches and settings
    if (diag_trap_levels.gt.0) then
       Write(PAR_OUT_FILE,"(A,I2)") "diag_trap_levels = ",diag_trap_levels
       do n=0,diag_trap_levels-1
          Write(PAR_OUT_FILE,"(A,I1,A,G12.4)") "diag_Blev(",n,") = ",diag_Blev(n)
       enddo
    endif
    IF (avgflux_stime.ge.0.0) THEN 
       Write(PAR_OUT_FILE,"(A,G12.4)") &
            &'avgflux_stime = ',avgflux_stime
       IF (avgflux_type.ne.0) THEN
          Write(PAR_OUT_FILE,"(A,I3)") "avgflux_type  = ", avgflux_type
          Write(PAR_OUT_FILE,"(A,G12.4)") "avgflux_alpha = ", avgflux_alpha
       ENDIF
    ENDIF

    IF (avgprof_stime.ge.0.0) Write(PAR_OUT_FILE,"(A,G12.4)") &
         &'avgprof_stime = ',avgprof_stime

    IF (gav_stime.ge.0.0) Write(PAR_OUT_FILE,"(A,G12.4)") &
         &'gav_stime     = ',gav_stime

    IF (trap_pass) Write(PAR_OUT_FILE,"(A,L1)") "trap_pass     = ", trap_pass


    WRITE(PAR_OUT_FILE,"(A)")   "/"
    Write(PAR_OUT_FILE,"(A)")    ""

    !****** extended_diags namelist ********    
    if(istep_nlev.gt.0.or.istep_nlt.gt.0.or.SVD_proj.or.istep_dfout.gt.0) then

     Write(PAR_OUT_FILE,"(A)")    "&extended_diags"       
     if (istep_nlt.gt.0) then
       write(PAR_OUT_FILE,"(A,1I2)")  'num_nlt_modes = ', num_nlt_modes
       if(nlt_pod) write(PAR_OUT_FILE,"(A,1I2)")  'num_nlt_pod_modes = ', num_nlt_pod_modes
       write(PAR_OUT_FILE,'(A,60(1x, 1I3))') 'ky_nlt_ind    = ',(ky_nlt_ind(j),j=1,60) 
       write(PAR_OUT_FILE,'(A,60(1x, 1I3))') 'kx_nlt_ind    = ',(kx_nlt_ind(j),j=1,60) 
       write(PAR_OUT_FILE,"(A)")    ""
     endif

     if(SVD_proj.and..not.svd_field_mom) then
      Write(PAR_OUT_FILE,"(A,L3)") "SVD_proj   = ",SVD_proj
      Write(PAR_OUT_FILE,"(A,I2)") "SVD_kx_ind = ",SVD_kx_ind
      Write(PAR_OUT_FILE,"(A,I2)") "SVD_ky_ind = ",SVD_ky_ind
      Write(PAR_OUT_FILE,"(A,I2)") "SVD_nkx0   = ",SVD_nkx0
      Write(PAR_OUT_FILE,"(A,20I5)") "SVD_df_n_time  = ",SVD_df_n_time
      Write(PAR_OUT_FILE,"(A,F10.6)") "SVD_start_time = ",SVD_start_time
      Write(PAR_OUT_FILE,"(A,I2)") "SVD_sparse_factor = ",SVD_sparse_factor
      write(PAR_OUT_FILE,"(3A)")   "SVD_df_file_name  = '",Trim(SVD_df_file_name),"'"
      write(PAR_OUT_FILE,"(3A)")   "SVD_df_file_path  = '",Trim(SVD_df_file_path),"'"
      Write(PAR_OUT_FILE,"(A,L3)") "SVD_f0_flag = ",SVD_f0_flag
     end if

     if(SVD_field_mom) then
      Write(PAR_OUT_FILE,"(A,L3)") "SVD_field_mom = ",SVD_field_mom
      Write(PAR_OUT_FILE,"(A,20I5)")  "SVD_df_n_time  = ",SVD_df_n_time
      Write(PAR_OUT_FILE,"(A,F10.6)") "SVD_start_time = ",SVD_start_time
      Write(PAR_OUT_FILE,"(A,I2)") "SVD_sparse_factor = ",SVD_sparse_factor
      write(PAR_OUT_FILE,"(3A)")   "SVD_df_file_name  = '",Trim(SVD_df_file_name),"'"
      write(PAR_OUT_FILE,"(3A)")   "SVD_df_file_path  = '",Trim(SVD_df_file_path),"'"
     end if

     if(istep_dfout.gt.0) then
      Write(PAR_OUT_FILE,"(A,L3)") "dfout_mpiio  = ",dfout_mpiio
      Write(PAR_OUT_FILE,"(A,I2)") "num_ky_modes = ",num_ky_modes
      Write(PAR_OUT_FILE,"(A,I2)") "num_kx_modes = ",num_kx_modes
      Write(PAR_OUT_FILE,"(A,51I3)") "which_ky = ",which_ky
      Write(PAR_OUT_FILE,"(A,51I3)") "which_kx_center= ",which_kx_center
     end if

     if (istep_nlev.gt.0) then
      if(nlev_stime.gt.-1) Write(PAR_OUT_FILE,"(A,G12.4)") "nlev_stime = ", nlev_stime
      if(n_nlev.gt.-1) Write(PAR_OUT_FILE,"(A,I12)") "n_nlev = ", n_nlev
     endif

     WRITE(PAR_OUT_FILE,"(A)")   "/"
     Write(PAR_OUT_FILE,"(A)")    ""

    end if


    !****** nonlocal_x namelist ********    
    if (.not.x_local) then
       Write(PAR_OUT_FILE,"(A)")    "&nonlocal_x"       

       IF (shifted_metric) WRITE(PAR_OUT_FILE,"(A,L1)") "shifted_metric = ",shifted_metric
       IF (.not.arakawa) WRITE(PAR_OUT_FILE,"(A,L1)") "arakawa = ", arakawa
       IF (l_buffer_size.gt.0.0) &
            &Write(PAR_OUT_FILE,"(A,G12.4)") "l_buffer_size = ", l_buffer_size
       IF (lcoef_krook.gt.0.0)&
            &Write(PAR_OUT_FILE,"(A,G12.4)") "lcoef_krook   = ", lcoef_krook
       IF (lpow_krook.ne.4) &
            &Write(PAR_OUT_FILE,"(A,I3)")    "lpow_krook    = ", lpow_krook
       IF (u_buffer_size.gt.0.0) &
            &Write(PAR_OUT_FILE,"(A,G12.4)") "u_buffer_size = ", u_buffer_size
       IF (ucoef_krook.gt.0.0)&
            &Write(PAR_OUT_FILE,"(A,G12.4)") "ucoef_krook   = ", ucoef_krook
       IF (upow_krook.ne.4) &
            &Write(PAR_OUT_FILE,"(A,I3)")    "upow_krook    = ", upow_krook
       IF (buffer_on_ky0) &
            &Write(PAR_OUT_FILE,"(A,L1)")    "buffer_on_ky0 = ", buffer_on_ky0
       IF (lck_buffer_size.gt.0.0) &
            &Write(PAR_OUT_FILE,"(A,G12.4)") "lck_buffer_size = ", lck_buffer_size
       IF (lckcoef_krook.ge.0.0)&
            &Write(PAR_OUT_FILE,"(A,G12.4)") "lckcoef_krook   = ", lckcoef_krook
       IF (lckpow_krook.ne.4) &
            &Write(PAR_OUT_FILE,"(A,I3)")    "lckpow_krook    = ", lckpow_krook
       IF (uck_buffer_size.gt.0.0) &
            &Write(PAR_OUT_FILE,"(A,G12.4)") "uck_buffer_size = ", uck_buffer_size
       IF (uckcoef_krook.ge.0.0)&
            &Write(PAR_OUT_FILE,"(A,G12.4)") "uckcoef_krook   = ", uckcoef_krook
       IF (uckpow_krook.ne.4) &
            &Write(PAR_OUT_FILE,"(A,I3)")    "uckpow_krook    = ", uckpow_krook 
       if (drive_buffer) &
            &write(PAR_OUT_FILE,"(A,L1)") "drive_buffer = ", drive_buffer
       if (explicit_buffer) &
            &write(PAR_OUT_FILE,"(A,L1)") "explicit_buffer = ", explicit_buffer

       IF (ck_heat.ne.0.0) Write(PAR_OUT_FILE,"(A,G12.4)") "ck_heat = ", ck_heat
       IF (ck_part.ne.0.0) Write(PAR_OUT_FILE,"(A,G12.4)") "ck_part = ", ck_part

       IF (ck_heat_smooth.ne.0.0) Write(PAR_OUT_FILE,"(A,G12.4)") &
            &"ck_heat_smooth = ", ck_heat_smooth
       IF (ck_part_smooth.ne.0.0) Write(PAR_OUT_FILE,"(A,G12.4)") &
            &"ck_part_smooth = ", ck_part_smooth
       if (ck_part_smooth.ne.0.0.or.ck_heat_smooth.ne.0.0) Write(PAR_OUT_FILE,"(A,I3)") &
            &"ck_filter_type = ", ck_filter_type
       IF (intsource_heat_coeff.gt.0) &
            &Write(PAR_OUT_FILE,"(A,G12.4)") "intsource_heat_coeff = ", intsource_heat_coeff
       IF (intsource_part_coeff.gt.0) &
            &Write(PAR_OUT_FILE,"(A,G12.4)") "intsource_part_coeff = ", intsource_part_coeff
       IF (intsource_time.ne.1000.0) &
            &Write(PAR_OUT_FILE,"(A,G12.4)") "intsource_time       = ", intsource_time

       IF (puls_amp.ne.1.0) THEN
          Write(PAR_OUT_FILE,"(A,G12.4)") "puls_amp       = ", puls_amp
          Write(PAR_OUT_FILE,"(A,G12.4)") "puls_start     = ", puls_start
          Write(PAR_OUT_FILE,"(A,G12.4)") "puls_on_duration  = ", puls_on_duration
          Write(PAR_OUT_FILE,"(A,G12.4)") "puls_off_duration = ", puls_off_duration
       ENDIF

       if (.not. nc_exb_corr) &
            &write(PAR_OUT_FILE,"(A,L1)") "nc_exb_corr = ", nc_exb_corr
       IF (reset_limit.ne.1000.0) THEN
          Write(PAR_OUT_FILE,"(A,G12.4)") "reset_limit    = ", reset_limit
          IF (smooth_reset) then
             Write(PAR_OUT_FILE,"(A,L1)") "smooth_reset = ", smooth_reset 
             Write(PAR_OUT_FILE,"(3A)") "smooth_reset_type = '", smooth_reset_type,"'"
             Write(PAR_OUT_FILE,"(A,G12.4)") "smooth_reset_width = ", smooth_reset_width
          endif
       ENDIF
       Write(PAR_OUT_FILE,"(A,I2)") "rad_bc_type = ", rad_bc_type

       WRITE(PAR_OUT_FILE,"(A)")   "/"
       Write(PAR_OUT_FILE,"(A)")    ""

    endif

    !****** external contributions namelist ********     
    if (with_external.or.(ExB).or.(pfsrate.ne.0.0).or.(Omega0_tor.ne.0.0)) then
       Write(PAR_OUT_FILE,"(A)")    "&external_contr"   
       if (with_phi_ext) then
          Write(PAR_OUT_FILE,"(A,G12.4)") "phi0_ext  = ",phi0_ext
          Write(PAR_OUT_FILE,"(A,I4)")    "kxind_phi_ext = ",kxind_phi_ext
          Write(PAR_OUT_FILE,"(A,G12.4)") "phase_phi_ext  = ",phase_phi_ext
       endif
       if (with_apar_ext) then
          Write(PAR_OUT_FILE,"(A,G12.4)") "apar0_ext  = ",apar0_ext
          Write(PAR_OUT_FILE,"(A,I4)")    "kxind_apar_ext = ",kxind_apar_ext
          Write(PAR_OUT_FILE,"(A,G12.4)") "phase_apar_ext  = ",phase_apar_ext
       endif
       if (with_omt_ext) then
          Write(PAR_OUT_FILE,"(A,G12.4)") "omt0_ext  = ", omt0_ext
          Write(PAR_OUT_FILE,"(A,I4)")    "kxind_omt_ext = ",kxind_omt_ext
          Write(PAR_OUT_FILE,"(A,G12.4)") "phase_omt_ext  = ",phase_omt_ext
       endif
       if (with_omn_ext) then
          Write(PAR_OUT_FILE,"(A,G12.4)") "omn0_ext  = ", omn0_ext
          Write(PAR_OUT_FILE,"(A,I4)")    "kxind_omn_ext = ",kxind_omn_ext
          Write(PAR_OUT_FILE,"(A,G12.4)") "phase_omn_ext  = ",phase_omn_ext
       endif
       if (Omega0_tor.ne.0.0) then
          write(PAR_OUT_FILE,"(A,G12.4)") "Omega0_tor = ", Omega0_tor
          write(PAR_OUT_FILE,"(A,G12.4)") "R0_tor     = ", R0_tor
          write(PAR_OUT_FILE,"(A,G12.4)") "dR0_tor_dx = ", dR0_tor_dx
          Write(PAR_OUT_FILE,"(A)") ""
          write(PAR_OUT_FILE,"(A,L1)") "with_comoving_other         = ", with_comoving_other
          write(PAR_OUT_FILE,"(A,L1)") "with_coriolis    = ", with_coriolis
          write(PAR_OUT_FILE,"(A,L1)") "with_centrifugal = ", with_centrifugal
          if (.not.with_bxphi0) write(PAR_OUT_FILE,"(A,L1)") &
               &"with_bxphi0      = ", with_bxphi0
          Write(PAR_OUT_FILE,"(A)") ""
       endif

       if (ExB) write(PAR_OUT_FILE,"(A,G12.4)") "ExBrate    = ", ExBrate
       if (lilo_w_full_omegator_profile) write(PAR_OUT_FILE,"(A,L1)") &
           "lilo_w_full_omegator_profile    = ", lilo_w_full_omegator_profile
       if (ExB_stime.gt.0.0) &
            write(PAR_OUT_FILE,"(A,G12.4)") "ExB_stime  = ", ExB_stime     
       if (pfsrate.ne.0.0) &
            write(PAR_OUT_FILE,"(A,G12.4)") "pfsrate    = ", pfsrate


       WRITE(PAR_OUT_FILE,"(A)") "/"
       Write(PAR_OUT_FILE,"(A)") ""
    endif

    !****** geometry namelist ********
    Write(PAR_OUT_FILE,"(A)")  "&geometry"
    WRITE(PAR_OUT_FILE,"(3A)") "magn_geometry = '", TRIM(magn_geometry),"'"

    if (INDEX(magn_geometry,'slab').eq.1) then
       Write(PAR_OUT_FILE,"(A,F8.4)")  "shat     = ",shat
       Write(PAR_OUT_FILE,"(A,G16.8)") "parscale = ", parscale
       if (trpeps.gt.0.0) then
          Write(PAR_OUT_FILE,"(A,G16.8)") "trpeps   = ",trpeps
       endif
    else
       if ((magn_geometry.eq.'gist').and.(y_local).and.(ialpha.ne.0)) &
            Write(PAR_OUT_FILE,"(A,I4)") "ialpha   = ",ialpha
       if ((x_local).or.(.not.mag_prof)) then
          Write(PAR_OUT_FILE,"(A,G16.8)") "q0       = ",q0
          Write(PAR_OUT_FILE,"(A,F8.4)")  "shat     = ",shat
       endif
       if ((magn_geometry .eq. 'circular') .and. &
            &(abs(sum(Bprof_coeffs)) .gt. epsilon(1.0))) then
          call variable_fltarr2strarr(Bprof_coeffs,tmpstr)
          write(PAR_OUT_FILE,"(2A)") "Bprof_coeffs = ", trim(tmpstr)
       end if

       if((magn_geometry.eq.'s_alpha').or.(magn_geometry.eq.'s_alpha_B').or.&
            &(magn_geometry.eq.'eq_cpo').or.(magn_geometry.eq.'circular')) THEN
          Write(PAR_OUT_FILE,"(A,G16.8)") "trpeps   = ",trpeps
          Write(PAR_OUT_FILE,"(A,G16.8)") "major_R  = ",major_R
          if ((amhd.gt.0.0).and.(magn_geometry.eq.'s_alpha').or.&
               &(magn_geometry.eq.'s_alpha_B')) &
               &Write(PAR_OUT_FILE,"(A,F8.4)") "amhd     = ",amhd
          if (.not.x_local) Write(PAR_OUT_FILE,"(A,G16.8)") "minor_r  = ", minor_r
          if (mag_prof) then
             call variable_fltarr2strarr(q_coeffs,tmpstr)
             Write(PAR_OUT_FILE,"(2A)") "q_coeffs  = ",TRIM(tmpstr)
          end if
       elseif (INDEX(magn_geometry,'miller').eq.1) then
          Write(PAR_OUT_FILE,"(A,G16.8)") "amhd     = ",amhd
          Write(PAR_OUT_FILE,"(A,G16.8)") "major_R  = ",major_R
          if (major_Z.ne.0.0) &
               Write(PAR_OUT_FILE,"(A,G16.8)") "major_Z  = ",major_Z
          Write(PAR_OUT_FILE,"(A,G16.8)") "minor_r  = ",minor_r
!          Write(PAR_OUT_FILE,"(A,G16.8)") "rho      = ",rho
          Write(PAR_OUT_FILE,"(A,G16.8)") "trpeps   = ",trpeps
          Write(PAR_OUT_FILE,"(A,G16.8)") "kappa    = ",kappa
          Write(PAR_OUT_FILE,"(A,G16.8)") "delta    = ",delta
          Write(PAR_OUT_FILE,"(A,G16.8)") "zeta     = ",zeta
          Write(PAR_OUT_FILE,"(A,G16.8)") "s_kappa  = ",s_kappa
          Write(PAR_OUT_FILE,"(A,G16.8)") "s_delta  = ",s_delta
          Write(PAR_OUT_FILE,"(A,G16.8)") "s_zeta   = ",s_zeta
          Write(PAR_OUT_FILE,"(A,G16.8)") "drR      = ",drR
          if (drZ.ne.0.0) Write(PAR_OUT_FILE,"(A,G16.8)") "drZ      = ",drZ
          if (edge_opt.gt.0) Write(PAR_OUT_FILE,"(A,G16.8)") "edge_opt = ",edge_opt
       else
          if (geomdir.ne.'./') WRITE(PAR_OUT_FILE,"(3A)") "geomdir  = '", TRIM(geomdir),"'"
          if (geomfile.ne.'') WRITE(PAR_OUT_FILE,"(3A)") "geomfile = '", TRIM(geomfile),"'"       
          if(magn_geometry.eq.'chease') THEN
             Write(PAR_OUT_FILE,"(3A)")     "x_def = '", TRIM(x_def),"'"
             if (x_local) Write(PAR_OUT_FILE,"(A,G16.8)") "flux_pos = ", flux_pos
          elseif(magn_geometry.eq.'tracer_efit') THEN
             if (edge_opt .gt. 0) Write(PAR_OUT_FILE,"(A,G16.8)") "edge_opt = ",edge_opt
          endif
          if (force_x) Write(PAR_OUT_FILE,"(A,G16.8)") "force_x  = ",force_x
          if (force_z) Write(PAR_OUT_FILE,"(A,G16.8)") "force_z  = ",force_z
          Write(PAR_OUT_FILE,"(A,G16.8)") "minor_r  = ",minor_r
          write(PAR_OUT_FILE,"(A,G16.8)") "major_R  = ", major_R
       endif
    endif

    if (mag_prof) WRITE(PAR_OUT_FILE,"(A,L1)") "mag_prof   = ", mag_prof 
    if (limit_shat) Write(PAR_OUT_FILE,"(A,L1)") "limit_shat = ",limit_shat
    if (rhostar.gt.0.0) Write(PAR_OUT_FILE,"(A,G16.8)") "rhostar  = ",rhostar 

    WRITE(PAR_OUT_FILE,"(3A)")      "dpdx_term= '",TRIM(dpdx_term),"'"
    WRITE(PAR_OUT_FILE,"(A,G16.8)") "dpdx_pm  = ",dpdx_pm

    WRITE(PAR_OUT_FILE,"(A,L1)") "norm_flux_projection  = ", norm_flux_projection
    Write(PAR_OUT_FILE,"(A)")   "/"
    Write(PAR_OUT_FILE,"(A)")    ""

    !****** species namelists ********
    do n=0,n_spec-1
       Write(PAR_OUT_FILE,"(A)")    "&species"
       Write(PAR_OUT_FILE,"(3A)") "name   = '", Trim(spec(n)%name),"'"
       if (spec(n)%passive) Write(PAR_OUT_FILE,"(A,L1)") "passive = ", spec(n)%passive
       if ((x_local).or.lilo.or.(spec(n)%prof_type.eq.0)) Then
          Write(PAR_OUT_FILE,"(A,G16.8)") "omn    = ", spec(n)%omn
          Write(PAR_OUT_FILE,"(A,G16.8)") "omt    = ", spec(n)%omt
          if ((.not.x_local).and.(spec(n)%prof_type.eq.0)) &
               &Write(PAR_OUT_FILE,"(A)") "prof_type = 0"
       else
          Write(PAR_OUT_FILE,"(A,I2)")    "prof_type = ", spec(n)%prof_type
          Write(PAR_OUT_FILE,"(A,G16.8)") "kappa_T   = ", spec(n)%kappa_T
          Write(PAR_OUT_FILE,"(A,G16.8)") "LT_center = ", spec(n)%LT_center
          Write(PAR_OUT_FILE,"(A,G16.8)") "LT_width  = ", spec(n)%LT_width
          Write(PAR_OUT_FILE,"(A)")    "" 
          Write(PAR_OUT_FILE,"(A,G16.8)") "kappa_n   = ", spec(n)%kappa_n
          Write(PAR_OUT_FILE,"(A,G16.8)") "Ln_center = ", spec(n)%Ln_center
          Write(PAR_OUT_FILE,"(A,G16.8)") "Ln_width  = ", spec(n)%Ln_width
          if ((spec(n)%prof_type.eq.3).or.(spec(n)%prof_type.eq.4)) then
             Write(PAR_OUT_FILE,"(A)")    ""
             Write(PAR_OUT_FILE,"(A,G16.8)") "delta_x_T  = ", spec(n)%delta_x_T
             Write(PAR_OUT_FILE,"(A,G16.8)") "delta_x_n  = ", spec(n)%delta_x_n
          endif
       end if
       Write(PAR_OUT_FILE,"(A)")    ""
       Write(PAR_OUT_FILE,"(A,G16.8)") "mass   = ", spec(n)%mass
       Write(PAR_OUT_FILE,"(A,G16.8)") "temp   = ", spec(n)%temp
       Write(PAR_OUT_FILE,"(A,G16.8)") "dens   = ", spec(n)%dens
       Write(PAR_OUT_FILE,"(A,I2)") "charge = ", spec(n)%charge

       if ((SUM(spec(n)%src_amp).gt.0).and.(SUM(spec(n)%src_prof_type).gt.0)) then
          Write(PAR_OUT_FILE,"(2A)") ""
          call variable_fltarr2strarr(spec(n)%src_prof_type,tmpstr)
          Write(PAR_OUT_FILE,"(2A)") "src_prof_type = ", TRIM(tmpstr)
          call variable_fltarr2strarr(spec(n)%src_amp,tmpstr)
          Write(PAR_OUT_FILE,"(2A)") "src_amp       = ", TRIM(tmpstr)
          call variable_fltarr2strarr(spec(n)%src_width,tmpstr,0.1)
          Write(PAR_OUT_FILE,"(2A)") "src_width     = ", TRIM(tmpstr)
          call variable_fltarr2strarr(spec(n)%src_x0,tmpstr,0.5)
          Write(PAR_OUT_FILE,"(2A)") "src_x0        = ", TRIM(tmpstr)
       endif

       Write(PAR_OUT_FILE,"(A)")    "/"
       Write(PAR_OUT_FILE,"(A)")    ""
    enddo

    !****** information namelist ********
    !Contains information that is not required for a new simulation
    WRITE(PAR_OUT_FILE,"(A)")    "&info"
    if (read_checkpoint.and.(checkpoint_in.ne.''))&
         &Write(PAR_OUT_FILE,"(3A)") "chpt_in = '", TRIM(checkpoint_in),"'"

    if ((magn_geometry.ne.'s_alpha').and.(magn_geometry.ne.'s_alpha_B')&
         &.and.(magn_geometry.ne.'circular'))&
         &Write(PAR_OUT_FILE,"(A,G16.8)") "q0     = ",q0

    select case (comp_type)
    case('EV')
       Write(PAR_OUT_FILE,"(A,I7)") "iterations for eigenvalue calculation = ",it_ev
       Write(PAR_OUT_FILE,"(A,F10.4)") "time for eigenvalue solver = ",time_ev
    case('IV')
       If (itime.ne.1) WRITE(PAR_OUT_FILE,"(A,F10.4)") "step_time  = ", time_iv/(itime-1)
       Write(PAR_OUT_FILE,"(A,I7)") "number of computed time steps = ",itime-1
       Write(PAR_OUT_FILE,"(A,F10.3)") "time for initial value solver = ",time_iv
       write(PAR_OUT_FILE,"(A,L1)") "calc_dt = ",calc_dt
       if (collision_op.ne.'none') &
            &Write(PAR_OUT_FILE,"(A,G16.8)") "ev_coll_est  = ",ev_coll_est
       if (precomp_nc) Write(PAR_OUT_FILE,"(A,F10.3)") "time for nc equilibrium solver = ",time_nc
    case('NC')
       Write(PAR_OUT_FILE,"(A,F10.3)") "time for nc equilibrium solver = ",time_nc
    end select
    
    if (arakawa_zv) &
         write(PAR_OUT_FILE,"(A,L1)") "hypz compensation = ",hypz_compensation

    Write(PAR_OUT_FILE,"(A,F10.4)") "init_time = ",initime
    write(PAR_OUT_FILE,"(A,I1)") "n_fields = ", n_fields
    write(PAR_OUT_FILE,"(A,I2)") "n_moms   = ", n_moms
    if (momentum_flux) &
         write(PAR_OUT_FILE,"(A,I2)") "nrgcols  = ", 10
    !write(PAR_OUT_FILE,"(A,I2)") "n_energies   = ", n_energies

    Write(PAR_OUT_FILE,"(A,G13.6)") "lx = ",lx
    Write(PAR_OUT_FILE,"(A,G13.6)") "ly = ",ly
    If (n_pol.gt.1) Write(PAR_OUT_FILE,"(A,I7)") "n_pol = ",n_pol

    Write(PAR_OUT_FILE,"(2A)") "PRECISION  = ",prec

    IF (IACHAR(TRANSFER(1,"a")).eq.0) THEN 
       WRITE(PAR_OUT_FILE,"(A)") "ENDIANNESS = BIG"
       if(SVD_proj) then
          SVD_endian='BIG_ENDIAN'
         if(SVD_swap_endian) then
#ifdef F2003_NO_OPEN_CONVERT
             if (mype.eq.0) WRITE(*,*) &
                  &"WARNING: SVD Endianness conversion not supported by this compiler"
#else
            SVD_endian='LITTLE_ENDIAN'
#endif
         endif
       end if
    ELSE
       WRITE(PAR_OUT_FILE,"(A)") "ENDIANNESS = LITTLE"
       if(SVD_proj) then
          SVD_endian='LITTLE_ENDIAN'
          if(SVD_swap_endian) then
#ifdef F2003_NO_OPEN_CONVERT
             if (mype.eq.0) WRITE(*,*) &
                  &"WARNING: SVD Endianness conversion not supported by this compiler"
#else
             SVD_endian='BIG_ENDIAN'
#endif
          end if
       end if
    END IF

    Write(PAR_OUT_FILE,"(A,I3)") "OMP_NUM_THREADS = ",omp_num_threads
    !    Write(PAR_OUT_FILE,"(A,F5.0)") "mb_per_core = ",mb_per_core

    if (svn_rev.ne.'')  Write(PAR_OUT_FILE,"(A,A)") "SVN_REV = ", svn_rev
    if (release.ne.'')  Write(PAR_OUT_FILE,"(A,A)") "RELEASE = ", release

    IF ((equil_par_curr) .AND. (x_local)) &
         & Write(PAR_OUT_FILE,"(A,F10.6)") "currdens_par = ", currdens_par_scal
    
    if (collision_op.ne.'none') then
       if (nu_ei.gt.0) write(PAR_OUT_FILE,"(A,F10.6)") "nu_ei = ", nu_ei
       if (nustar_i.gt.0) write(PAR_OUT_FILE,"(A,F10.6)") "nustar_i = ", nustar_i
       if (nustar_e.gt.0) write(PAR_OUT_FILE,"(A,F10.6)") "nustar_e = ", nustar_e
    endif
    
    !switches
#ifndef GDAGGER_G
    WRITE(PAR_OUT_FILE,"(A)") "GDAGGER_G = F"
#endif
#ifndef G_DAGGER
    WRITE(PAR_OUT_FILE,"(A)") "G_DAGGER  = F"
#endif
#ifndef APPROX_METRIC
    WRITE(PAR_OUT_FILE,"(A)") "APPROX_METRIC = F"
#endif
#ifndef prec_x
    WRITE(PAR_OUT_FILE,"(A)") "prec_x   = F"
#endif
#ifndef prec_ara
    WRITE(PAR_OUT_FILE,"(A)") "prec_ara = F"
#endif
#ifdef hyp_over_j
    WRITE(PAR_OUT_FILE,"(A)") "hyp_over_j = T"
#endif
    WRITE(PAR_OUT_FILE,"(A)")    "/"
    WRITE(PAR_OUT_FILE,"(A)")    ""

    !****** unit namelist ********
    ! Contains reference quantities for conversion to SI units
    WRITE(PAR_OUT_FILE,"(A)")    "&units"
    if (Bref.ne.0.) WRITE(PAR_OUT_FILE,"(A,ES24.16)") "Bref = ", Bref
    if (Tref.ne.0.) WRITE(PAR_OUT_FILE,"(A,ES24.16)") "Tref = ", Tref
    if (nref.ne.0.) WRITE(PAR_OUT_FILE,"(A,ES24.16)") "nref = ", nref
    if (Lref.ne.0.) WRITE(PAR_OUT_FILE,"(A,ES24.16)") "Lref = ", Lref
    if (mref.ne.0.) WRITE(PAR_OUT_FILE,"(A,ES24.16)") "mref = ", mref
    if (omegatorref.ne.0.) WRITE(PAR_OUT_FILE,"(A,ES24.16)") &
         "omegatorref = ", omegatorref
    if (norm_index.ne.-1) WRITE(PAR_OUT_FILE,"(A,I3)") "norm_index = ", norm_index
    WRITE(PAR_OUT_FILE,"(A)")    "/"
    WRITE(PAR_OUT_FILE,"(A)")    ""

    if(.not.cat_output) CLOSE(PAR_OUT_FILE)

    ! All the parameters are gathered here
#ifdef WITHFUTILS
    if (write_h5.and.mype.eq.0)then
       fileall_param_h5 = trim(diagdir)//'/all_params'//trim(file_extension)//'.h5'
       call creatf(fileall_param_h5, fidall_param_h5, "all parameters", 'd')

       ! I/O parameters attached as attributes in a separate file
       call creatg(fidall_param_h5, "/input", "input variables")
       
    !****** parallelization namelist ********
       call attach(fidall_param_h5, "/input", "n_procs_s"  , n_procs_s)            
       call attach(fidall_param_h5, "/input", "n_procs_v"  , n_procs_v)            
       call attach(fidall_param_h5, "/input", "n_procs_w"  , n_procs_w)            
       call attach(fidall_param_h5, "/input", "n_procs_x"  , n_procs_x)            
       call attach(fidall_param_h5, "/input", "n_procs_y"  , n_procs_y)            
       call attach(fidall_param_h5, "/input", "n_procs_z"  , n_procs_z)            
       call attach(fidall_param_h5, "/input", "n_procs_sim", n_procs_sim)            
       if(n_parallel_sims.gt.1) &
            &call attach(fidall_param_h5, "/input", "n_parallel_sims", n_parallel_sims)            

    !****** box namelist ********
       call attach(fidall_param_h5, "/input", "n_spec", n_spec)
       call attach(fidall_param_h5, "/input", "nx0"   , nx0)
       call attach(fidall_param_h5, "/input", "nky0"  , nky0)
       call attach(fidall_param_h5, "/input", "nz0"   , nz0)
       call attach(fidall_param_h5, "/input", "nv0"   , nv0)
       call attach(fidall_param_h5, "/input", "nw0"   , nw0)
       call attach(fidall_param_h5, "/input", "kymin" , kymin)
       call attach(fidall_param_h5, "/input", "lv"    , lv)
       call attach(fidall_param_h5, "/input", "lw"    , lw)
       if (x_local) then
          if (.not.adapt_lx) then
             call attach(fidall_param_h5, "/input", "lx"       , lx)            
             call attach(fidall_param_h5, "/input", "nexc"     , nexc)
          end if
          call attach(fidall_param_h5, "/input", "adapt_lx" , adapt_lx)
 
          if (magn_geometry.eq.'tracer_efit'.or.magn_geometry.eq.'gist'.or.mag_prof)&
               &call attach(fidall_param_h5, "/input", "x0"       , x0)
       else
          if (lx.ne.0.0) call attach(fidall_param_h5, "/input", "lx", lx)
          if (lx_a.ne.0.0) call attach(fidall_param_h5, "/input", "lx_a", lx_a)
          call attach(fidall_param_h5, "/input", "x0"       , x0)            
       end if

       if (n0_global.ne.-1111) call attach(fidall_param_h5, "/input", "n0_global", n0_global)
       IF ((xy_local.or.lilo).and.adapt_ly) &
            call attach(fidall_param_h5, "/input", "adapt_ly" , adapt_ly) 

       if (y_local.and.(.not.nonlinear)) call attach(fidall_param_h5, "/input", "ky0_ind"   , ky0_ind)
       if (kx_center.ne.0.) call attach(fidall_param_h5, "/input", "kx_center" , kx_center) 
       call attach(fidall_param_h5, "/input","mu_grid_type", trim(mu_grid_type))  

    !****** in_out namelist ********
       call attach(fidall_param_h5, "/input", "diagdir" , trim(diagdir))           
       if (TRIM(chptdir).ne.TRIM(diagdir)) call attach(fidall_param_h5, "/input", "chptdir" , trim(chptdir))            
       call attach(fidall_param_h5, "/input", "read_checkpoint"  , read_checkpoint)            
       call attach(fidall_param_h5, "/input", "write_checkpoint" , write_checkpoint)            
       if (many_chpts) call attach(fidall_param_h5, "/input", "many_chpts"  , many_chpts)            
       call attach(fidall_param_h5, "/input", "istep_field" , istep_field)            
       call attach(fidall_param_h5, "/input", "istep_mom"   , istep_mom)            
       call attach(fidall_param_h5, "/input", "istep_nrg"   , istep_nrg)            
       if (.not.nonlinear) call attach(fidall_param_h5, "/input", "istep_omega" , istep_omega)
       call attach(fidall_param_h5, "/input", "istep_vsp"   , istep_vsp)            
       call attach(fidall_param_h5, "/input", "istep_schpt" , istep_schpt)             
       if (istep_g1.gt.0) call attach(fidall_param_h5, "/input", "istep_g1", istep_g1)            
       if (istep_energy.gt.0) call attach(fidall_param_h5, "/input", "istep_energy", istep_energy)
       if (istep_energy3d.gt.0) call attach(fidall_param_h5, "/input", "istep_energy3d", &
            &istep_energy3d)
       if (istep_energy_exchange.gt.0) call attach(fidall_param_h5, "/input", "istep_energy_exchange", istep_energy_exchange)
       if (istep_nlt.gt.0) call attach(fidall_param_h5, "/input", "istep_nlt" , istep_nlt)
       if (istep_dfout.gt.0) call attach(fidall_param_h5, "/input", "istep_dfout" , istep_dfout)
       if (istep_prof.gt.0) call attach(fidall_param_h5, "/input", "istep_prof"  , istep_prof)            
       if ((istep_field .GT. 0).and.(write_fielddiff)) &
            &call attach(fidall_param_h5, "/input", "write_fielddiff" ,write_fielddiff)                   
       if (istep_srcmom.gt.0) call attach(fidall_param_h5, "/input", "istep_srcmom" , istep_srcmom)                   
       if (istep_neoclass.gt.0) call attach(fidall_param_h5, "/input", "istep_neoclass" , istep_neoclass)
       if (istep_neoclass2D.gt.0) call attach(fidall_param_h5, "/input", "istep_neoclass2D" , istep_neoclass2D)
       if (istep_fsa_moments.gt.0) &
            &call attach(fidall_param_h5, "/input", "istep_fsa_moments" , istep_fsa_moments)
       IF (GyroLES) call attach(fidall_param_h5, "/input","istep_GyroLES", istep_GyroLES)

       !diag_GyroLES
       IF (istep_fe_time.gt.0) call attach(fidall_param_h5, "/input","istep_fe_time", istep_fe_time)
       IF (istep_fe_twoD.gt.0) call attach(fidall_param_h5, "/input","istep_fe_twoD", istep_fe_twoD)
       IF (istep_fe_transfer.gt.0) call attach(fidall_param_h5, "/input","istep_fe_transfer", istep_fe_transfer)
       
       IF (istep_nlev.gt.0) call attach(fidall_param_h5, "/input","istep_nlev", istep_nlev)

       call attach(fidall_param_h5, "/input", "write_std"   , write_std)
       if (write_h5) call attach(fidall_param_h5, "/input", "write_h5", write_h5)
       if (chpt_h5) then
          call attach(fidall_param_h5, "/input", "chpt_h5"   , chpt_h5)
       else
          if (chpt_read_h5) &
               &call attach(fidall_param_h5, "/input", "chpt_read_h5", chpt_read_h5)
          if (chpt_write_h5) &
               &call attach(fidall_param_h5, "/input", "chpt_write_h5", chpt_write_h5)
       ENDIF
    
       IF (.not.momentum_flux) call attach(fidall_param_h5, "/input", "momentum_flux", momentum_flux)
       IF (write_flux_final.gt.0) call attach(fidall_param_h5, "/input", "write_flux_final", write_flux_final)
       IF (trim(iterdb_file).ne.'') call attach(fidall_param_h5, "/input", "iterdb_file", TRIM(iterdb_file))
       IF (iterdb_time.ge.0.0) call attach(fidall_param_h5, "/input", "iterdb_time", iterdb_time)

    !****** general namelist ********
       call attach(fidall_param_h5, "/input", "nonlinear" , nonlinear)
       if (parallel_nl) call attach(fidall_param_h5, "/input", "parallel_nl" , parallel_nl)

       if (.not.x_local) call attach(fidall_param_h5, "/input", "x_local"   , x_local)
       if (.not.y_local) call attach(fidall_param_h5, "/input", "y_local"   , y_local)
       if (y_local.and.yx_order) call attach(fidall_param_h5, "/input", "yx_order", yx_order)
       if (lilo) call attach(fidall_param_h5, "/input", "lilo"      , lilo)
       call attach(fidall_param_h5, "/input", "comp_type" , trim(comp_type))
       call putarr(fidall_param_h5, "/input/perf_vec", perf_vec, "perf_vec")
       call attach(fidall_param_h5, "/input", "nblocks"   , nblocks)
       if (parscheme.ne.'c4th') &
            &call attach(fidall_param_h5, "/input", "parscheme" , trim(parscheme))
       if (vparscheme.ne.'c4th') &
            &call attach(fidall_param_h5, "/input", "vparscheme" , trim(vparscheme))
       if (.not.arakawa_zv) then 
          call attach(fidall_param_h5, "/input", "arakawa_zv", arakawa_zv)
       else
          if (arakawa_zv_order.ne.4) &
               call attach(fidall_param_h5, "/input", "arakawa_zv_order", arakawa_zv_order)
          if (arakawa_cons_bc) call attach(fidall_param_h5, "/input", "arakawa_cons_bc", arakawa_cons_bc)
          if (.not.hyp_on_h) call attach(fidall_param_h5, "/input", "hyp_on_h", hyp_on_h)
          if (.not.hypz_opt) call attach(fidall_param_h5, "/input", "hypz_opt", hypz_opt)
       endif
       if (turbdeal) call attach(fidall_param_h5, "/input", "turbdeal"  , turbdeal)

       select case(comp_type)
       case('EV')
          call attach(fidall_param_h5, "/input", "which_ev"  , trim(which_ev))
          If (n_ev.ne.1.and.(which_ev=='jd_as'.or.which_ev=='kx_as')) &
               &call attach(fidall_param_h5, "/input", "which_func_as", which_func_as)
          if (pc_type.ne.'none') call attach(fidall_param_h5, "/input", "pc_type"   , trim(pc_type))
          if (vec_out_limit) call attach(fidall_param_h5, "/input", "vec_out_limit", vec_out_limit)
          if (n_ev.ne.1) call attach(fidall_param_h5, "/input", "n_ev"      , n_ev)
          if (ev_prec.ne.1e-4) call attach(fidall_param_h5, "/input", "ev_prec"   , ev_prec)
          If (ev_shift.ne.(10.,0.)) then
             call attach(fidall_param_h5, "/input", "ev_shift_real" , real(ev_shift))
             call attach(fidall_param_h5, "/input", "ev_shift_imag" , aimag(ev_shift))
          end If
          if (ev_max_it.gt.0) call attach(fidall_param_h5, "/input", "ev_max_it" , ev_max_it)
          if (ev_n_test.ne.(n_ev+15)) call attach(fidall_param_h5, "/input", "ev_n_test" , ev_n_test)
       case('IV')
          call attach(fidall_param_h5, "/input", "timescheme", trim(timescheme))
          if (collision_op.ne.'none') then
             if (.not.coll_split) then
                call attach(fidall_param_h5, "/input", "coll_split", coll_split)
             else
                if (coll_split_scheme .ne. 'RKCa') &
                     &call attach(fidall_param_h5, "/input", "coll_split_scheme", trim(coll_split_scheme))
             endif
          endif
          if (coll_split) then
             call attach(fidall_param_h5, "/input", "dt_max"    , dt_max)            
             call attach(fidall_param_h5, "/input", "dt_vlasov" , dt_vlasov)            
             call attach(fidall_param_h5, "/input", "ev_coll"   , ev_coll)            
          else
             call attach(fidall_param_h5, "/input", "dt_max"    , dt_max)            
          endif
          If (nonlinear) call attach(fidall_param_h5, "/input", "courant"   , courant)
          if (precomp_nc) call attach(fidall_param_h5, "/input", "precomp_nc", precomp_nc)
          if (include_f0_contr) call attach(fidall_param_h5, "/input", "include_f0_contr", include_f0_contr)
       case('NC')
          If (nc_prec.ne.5.e-7) call attach(fidall_param_h5, "/input", "nc_prec",nc_prec)
          If (nc_max_it.ne.100000) call attach(fidall_param_h5, "/input", "nc_max_it",nc_max_it)
          if (include_f0_contr) call attach(fidall_param_h5, "/input", "include_f0_contr", include_f0_contr)
       end select

       call attach(fidall_param_h5, "/input", "timelim"    , timelim)            
       if(comp_type.eq.'IV') then
          call attach(fidall_param_h5, "/input", "ntimesteps" , ntimesteps)            
          if (simtimelim.ne.500000.) call attach(fidall_param_h5, "/input", "simtimelim" , simtimelim)
       end if
       if (overflow_limit.ne.1e30) call attach(fidall_param_h5, "/input", "overflow_limit"  , overflow_limit)            
       if (underflow_limit.ne.1e-8) call attach(fidall_param_h5, "/input", "underflow_limit" , underflow_limit)            
       if ((istep_omega.gt.0).and.(omega_prec.ne.1e-3)) &
            &call attach(fidall_param_h5, "/input", "omega_prec"      , omega_prec)            

       call attach(fidall_param_h5, "/input", "beta"   , beta)            
       call attach(fidall_param_h5, "/input", "debye2" , debye2)            
       IF (beta .GT. 0.0 .and. bpar) &
            &call attach(fidall_param_h5, "/input", "bpar" , bpar)            
!       if (pressure_term.and.(beta.gt.0.0)) &
!            &call attach(fidall_param_h5, "/input", "pressure_term"  , pressure_term)            
       if (Apar0_antenna .NE. 0.0) then
          call attach(fidall_param_h5, "/input", "Apar0_antenna" , Apar0_antenna)            
          call attach(fidall_param_h5, "/input", "omega0_antenna", omega0_antenna)            
       end if
       if (tau.ne.1.) call attach(fidall_param_h5, "/input", "tau = ",tau)
       if (delzonal) call attach(fidall_param_h5, "/input", "delzonal", delzonal)            
       if (delzonal_fields) call attach(fidall_param_h5, "/input", "delzonal_fields", delzonal_fields)            
       if (delzonal_factor.ne.1.0) call attach(fidall_param_h5, "/input", "delzonal_factor", delzonal_factor)
       if (abs(add_zonal_phi).gt.epsilon(0.0)) call attach(fidall_param_h5, "/input", "add_zonal_phi", add_zonal_phi)
       if (del_phi) call attach(fidall_param_h5, "/input", "del_phi", del_phi)
       if (del_fields) call attach(fidall_param_h5, "/input", "del_fields", del_fields)
       if (only_Er) call attach(fidall_param_h5, "/input", "only_Er", only_Er)
       call attach(fidall_param_h5, "/input", "collision_op" , trim(collision_op))            
       if (collision_op.ne.'none') then 
          call attach(fidall_param_h5, "/input", "coll", coll)
          call attach(fidall_param_h5, "/input", "coll_cons_model",trim(coll_cons_model))            
          if (spacediff) call attach(fidall_param_h5, "/input", "spacediff" , spacediff)            
          if (coll_f_fm_on) call attach(fidall_param_h5, "/input", "coll_f_fm_on" , coll_f_fm_on)            
          if (zeff0.ne.1.0) call attach(fidall_param_h5, "/input", "Zeff", zeff0)            
          if (adiabatic_mass.ne.1.0) call attach(fidall_param_h5, "/input", "adiabatic_mass",adiabatic_mass)
       end if

       !initial condition
       call attach(fidall_param_h5, "/input", "init_cond"    , init_cond)           
       if (init_aux_x.ne.-100.) call attach(fidall_param_h5, "/input" , "init_aux_x"   , init_aux_x)            
       if (init_aux_y.ne.-100.) call attach(fidall_param_h5, "/input" , "init_aux_y"   , init_aux_y)           
       if (init_aux_z.ne.-100.) call attach(fidall_param_h5, "/input" , "init_aux_z"   , init_aux_z)            
       if (init_aux_amp.ne.-100.) call attach(fidall_param_h5,"/input", "init_aux_amp" , init_aux_amp)            

       !hyper diffusions
       if(hyp_x.gt.0.0) then
          call attach(fidall_param_h5, "/input", "hyp_x", hyp_x)           
          if (hyp_x_order.ne.4) & 
               &call attach(fidall_param_h5, "/input", "hyp_x_order"  , hyp_x_order)           
       end if
       if(hyp_y.gt.0.0) then
          call attach(fidall_param_h5, "/input", "hyp_y", hyp_y)            
          if (hyp_y_order.ne.4) & 
               &call attach(fidall_param_h5, "/input", "hyp_y_order"  , hyp_y_order)
       end if
       if(hyp_perp.gt.0.0) then
          call attach(fidall_param_h5, "/input", "hyp_perp", hyp_perp)                        
          if (hyp_perp_order.ne.4) &
               &call attach(fidall_param_h5, "/input", "hyp_perp_order" , hyp_perp_order)            
       end if
       if(hyp_z.ne.0.0) then
          call attach(fidall_param_h5, "/input", "hyp_z", hyp_z)           
          IF (hyp_z_order.ne.4) &
               &call attach(fidall_param_h5, "/input", "hyp_z_order" , hyp_z_order)
       end if
       if(hyp_v.gt.0.0) call attach(fidall_param_h5, "/input", "hyp_v" , hyp_v)

       !GyroLES
       if (GyroLES) then
          call attach(fidall_param_h5, "/input", "GyroLES" , GyroLES)
       endif

       IF (diag_GyroLES) then
          call attach(fidall_param_h5, "/input", "diag_GyroLES" , diag_GyroLES)
          call attach(fidall_param_h5, "/input", "fracx" , fracx)
          call attach(fidall_param_h5, "/input", "fracy" , fracy)
          if (istep_fe_transfer.gt.0) then
              call attach(fidall_param_h5, "/input", "tkx",tkx)
              call attach(fidall_param_h5, "/input", "tky",tky)
              call attach(fidall_param_h5, "/input", "tkperp",tkperp)
              call attach(fidall_param_h5, "/input", "t_log_ky",t_log_ky)
              call attach(fidall_param_h5, "/input", "t_log_kperp",t_log_kperp)
              call attach(fidall_param_h5, "/input", "t_log_kperp_fb",t_log_kperp_fb)
              call attach(fidall_param_h5, "/input", "t_log_kperp_h",t_log_kperp_h)
              call attach(fidall_param_h5, "/input", "t_log_kperp_1",t_log_kperp_1)
              call attach(fidall_param_h5, "/input", "tk3_kx",tk3_kx)
              call attach(fidall_param_h5, "/input", "tk3_ky",tk3_ky)
              call attach(fidall_param_h5, "/input", "tk3_log_kperp",tk3_log_kperp)
              call attach(fidall_param_h5, "/input", "tk3_log_kperp_1",tk3_log_kperp_1)
              call attach(fidall_param_h5, "/input", "tk3_log_kperp_h",tk3_log_kperp_h)
              call attach(fidall_param_h5, "/input", "n_shells_kperp",n_shells_kperp)
              call attach(fidall_param_h5, "/input", "n_shells_log_ky",n_shells_log_ky)
              call attach(fidall_param_h5, "/input", "n_shells_log_kperp",n_shells_log_kperp)
              call attach(fidall_param_h5, "/input", "n_shells_log_kperp_1",n_shells_log_kperp_1)
              call attach(fidall_param_h5, "/input","lambda1",lambda1)
              call attach(fidall_param_h5, "/input","lambda2",lambda2)
          endif
       endif 

       IF (coeff_g1_en_source .NE. 0.0) CALL &
         attach(fidall_param_h5,"/input","coeff_g1_en_source",coeff_g1_en_source)           

       !Misc.
       if (perf_tsteps.ne.3) then
          call attach(fidall_param_h5, "/input", "perf_tsteps" , perf_tsteps)           
       end if

       !diagnostics related switches and settings
       if (diag_trap_levels.gt.0 ) then
          call attach(fidall_param_h5, "/input", "diag_trap_levels" , diag_trap_levels)            
          call putarr(fidall_param_h5, "/input/diag_Blev", diag_Blev, "B_levels")
       end if
       if (avgflux_stime.ge.0.0) THEN
          call attach(fidall_param_h5, "/input", "avgflux_stime" , avgflux_stime)            
          IF (avgflux_type.ne.0) THEN
             call attach(fidall_param_h5, "/input", "avgflux_type" , avgflux_type)
             call attach(fidall_param_h5, "/input", "avgflux_alpha" , avgflux_alpha)
          ENDIF
       ENDIF
       if (avgprof_stime.ge.0.0) &
            &call attach(fidall_param_h5, "/input", "avgprof_stime" , avgprof_stime)
       if (gav_stime.ge.0.0) call attach(fidall_param_h5, "/input", "gav_stime",gav_stime)

       if (trap_pass)    call attach(fidall_param_h5, "/input", "trap_pass"   , trap_pass)           

    !****** extended_diags namelist ********    
       if (istep_nlt.gt.0) then
          call attach(fidall_param_h5, "/input", "num_nlt_modes"  , num_nlt_modes)           
          if (nlt_pod) call attach(fidall_param_h5, "/input", "num_nlt_pod_modes"  , num_nlt_pod_modes)           
          call putarr(fidall_param_h5, "/input/ky_nlt_ind"  , ky_nlt_ind,"ky_nlt_ind")           
          call putarr(fidall_param_h5, "/input/kx_nlt_ind"  , kx_nlt_ind,"kx_nlt_ind")           
       endif
       if (SVD_proj.and..not.svd_field_mom) then
          call attach(fidall_param_h5, "/input", "SVD_proj",SVD_proj)
          call attach(fidall_param_h5, "/input", "SVD_kx_ind",SVD_kx_ind)
          call attach(fidall_param_h5, "/input", "SVD_ky_ind",SVD_ky_ind)
          call attach(fidall_param_h5, "/input", "SVD_nkx0",SVD_nkx0)
          call putarr(fidall_param_h5, "/input/SVD_df_n_time", SVD_df_n_time,"SVD_df_n_time")
          call attach(fidall_param_h5, "/input", "SVD_start_time",SVD_start_time)
          call attach(fidall_param_h5, "/input", "SVD_sparse_factor",SVD_sparse_factor)
          call attach(fidall_param_h5, "/input", "SVD_df_file_name",Trim(SVD_df_file_name))
          call attach(fidall_param_h5, "/input", "SVD_df_file_path",Trim(SVD_df_file_path))
          call attach(fidall_param_h5, "/input", "SVD_f0_flag",SVD_f0_flag)
       end if
       if (SVD_field_mom) then
          call attach(fidall_param_h5,"/input","SVD_field_mom",SVD_field_mom)
          call putarr(fidall_param_h5, "/input/SVD_df_n_time", SVD_df_n_time,"SVD_df_n_time")
          call attach(fidall_param_h5,"/input","SVD_start_time",SVD_start_time)
          call attach(fidall_param_h5,"/input","SVD_sparse_factor",SVD_sparse_factor)
          call attach(fidall_param_h5,"/input","SVD_df_file_name",Trim(SVD_df_file_name))
          call attach(fidall_param_h5,"/input","SVD_df_file_path",Trim(SVD_df_file_path))
       end if
       if(istep_dfout.gt.0) then
          call attach(fidall_param_h5,"/input","dfout_mpiio",dfout_mpiio)
          call attach(fidall_param_h5,"/input","num_ky_modes",num_ky_modes)
          call attach(fidall_param_h5,"/input","num_kx_modes",num_kx_modes)
          call putarr(fidall_param_h5,"/input/which_ky", which_ky,"which_ky")
          call putarr(fidall_param_h5,"/input/which_kx_center", which_kx_center,"which_kx_center")
       end if
       if (istep_nlev.gt.0) then
          if (nlev_stime.gt.-1) call attach(fidall_param_h5,"/input","nlev_stime = ", nlev_stime)
          if (n_nlev.gt.-1)  call attach(fidall_param_h5,"/input","n_nlev = ", n_nlev)
       endif

    !****** nonlocal_x namelist ********    
       if (.not.x_local) then
          if (shifted_metric) call attach(fidall_param_h5, "/input", "shifted_metric" , shifted_metric)            
          if (.not.arakawa) call attach(fidall_param_h5, "/input", "arakawa" , arakawa)            
          if (l_buffer_size.gt.0.0) call attach(fidall_param_h5, "/input", "l_buffer_size"  , l_buffer_size )           
          if (lcoef_krook.gt.0.0) call attach(fidall_param_h5, "/input", "lcoef_krook"    , lcoef_krook )            
          if (lpow_krook.ne.4) call attach(fidall_param_h5, "/input", "lpow_krook"     , lpow_krook)           
          if (u_buffer_size.gt.0.0) call attach(fidall_param_h5, "/input", "u_buffer_size"  , u_buffer_size)            
          if (ucoef_krook.gt.0.0) call attach(fidall_param_h5, "/input", "ucoef_krook"    , ucoef_krook)           
          if (upow_krook.ne.4) call attach(fidall_param_h5, "/input", "upow_krook"     , upow_krook)            
          if (buffer_on_ky0) call attach(fidall_param_h5, "/input", "buffer_on_ky0"  , buffer_on_ky0)
          if (lck_buffer_size.gt.0.0) call attach(fidall_param_h5, "/input", "lck_buffer_size"  , lck_buffer_size)            
          if (lckcoef_krook.ge.0.0) call attach(fidall_param_h5, "/input", "lckcoef_krook"  , lckcoef_krook)            
          IF (lckpow_krook.ne.4)      call attach(fidall_param_h5, "/input","lckpow_krook", lckpow_krook)
          IF (uck_buffer_size.gt.0.0) call attach(fidall_param_h5, "/input","uck_buffer_size", uck_buffer_size)
          IF (uckcoef_krook.ge.0.0)   call attach(fidall_param_h5, "/input","uckcoef_krook", uckcoef_krook)
          IF (uckpow_krook.ne.4)      call attach(fidall_param_h5, "/input","uckpow_krook", uckpow_krook)
          if (drive_buffer)           call attach(fidall_param_h5, "/input","drive_buffer", drive_buffer)
          if (explicit_buffer)        call attach(fidall_param_h5, "/input","explicit_buffer", explicit_buffer)

          if (ck_heat.ne.0.0) call attach(fidall_param_h5, "/input", "ck_heat"        , ck_heat)           
          if (ck_part.ne.0.0) call attach(fidall_param_h5, "/input", "ck_part"        , ck_part)           

          IF (intsource_heat_coeff.gt.0) &
               &call attach(fidall_param_h5,"/input","intsource_heat_coeff", intsource_heat_coeff)
          IF (intsource_part_coeff.gt.0) &
               &call attach(fidall_param_h5,"/input","intsource_part_coeff", intsource_part_coeff)
          IF (intsource_time.ne.1000.0) &
               &call attach(fidall_param_h5,"/input","intsource_time", intsource_time)

          if (puls_amp.ne.1.0) then
             call attach(fidall_param_h5, "/input", "puls_amp"       , puls_amp)            
             call attach(fidall_param_h5, "/input", "puls_start"     , puls_start)            
             call attach(fidall_param_h5, "/input", "puls_on_duration" , puls_on_duration)           
             call attach(fidall_param_h5, "/input", "puls_off_duration", puls_off_duration)             
          end if

          if (reset_limit.ne.1000.0) THEN
             call attach(fidall_param_h5, "/input","reset_limit", reset_limit)
             IF (smooth_reset) call attach(fidall_param_h5, "/input", "smooth_reset", smooth_reset)
          endif
          call attach(fidall_param_h5, "/input", "rad_bc_type", rad_bc_type)             

       end if

    !****** external contributions namelist ********     
       if (with_external.or.(ExB).or.(pfsrate.ne.0.0)) then
          if (with_phi_ext) then
             call attach(fidall_param_h5, "/input", "phi0_ext"      , phi0_ext)           
             call attach(fidall_param_h5, "/input", "kxind_phi_ext" , kxind_phi_ext)            
             call attach(fidall_param_h5, "/input", "phase_phi_ext" , phase_phi_ext)           
          end if
          if (with_apar_ext) then
             call attach(fidall_param_h5, "/input", "apar0_ext"      , apar0_ext)           
             call attach(fidall_param_h5, "/input", "kxind_apar_ext" , kxind_apar_ext)            
             call attach(fidall_param_h5, "/input", "phase_apar_ext" , phase_apar_ext)           
          end if
          if (with_omt_ext) then
             call attach(fidall_param_h5, "/input", "omn0_ext"      , omt0_ext)           
             call attach(fidall_param_h5, "/input", "kxind_omt_ext" , kxind_omt_ext)           
             call attach(fidall_param_h5, "/input", "phase_omt_ext" , phase_omt_ext)            
          end if
          if (with_omn_ext) then
             call attach(fidall_param_h5, "/input", "omn0_ext"      , omn0_ext)            
             call attach(fidall_param_h5, "/input", "kxind_omt_ext" , kxind_omt_ext)            
             call attach(fidall_param_h5, "/input", "phase_omt_ext" , phase_omt_ext)            
          end if
          if (Omega0_tor.ne.0.0) then
             call attach(fidall_param_h5, "/input", "Omega0_tor", Omega0_tor) 
             call attach(fidall_param_h5, "/input", "R0_tor", R0_tor)
             call attach(fidall_param_h5, "/input", "dR0_tor_dx", dR0_tor_dx)
             call attach(fidall_param_h5, "/input", "with_comoving_other", with_comoving_other)
             call attach(fidall_param_h5, "/input", "with_coriolis", with_coriolis)
             call attach(fidall_param_h5, "/input", "with_centrifugal", with_centrifugal)
             call attach(fidall_param_h5, "/input", "with_bxphi0", with_bxphi0)
          endif
          
          if (ExB) call attach(fidall_param_h5, "/input", "ExBrate"   , ExBrate)
          if (lilo_w_full_omegator_profile) call attach(fidall_param_h5, "/input",&
              "lilo_w_full_omegator_profile"   , lilo_w_full_omegator_profile)
          if (ExB_stime.gt.0.0) call attach(fidall_param_h5, "/input", "ExB_stime" , ExB_stime)
          if (pfsrate.ne.0) call attach(fidall_param_h5, "/input", "pfsrate"   , pfsrate)
       end if

    !****** geometry namelist ********
       call attach(fidall_param_h5, "/input", "magn_geometry" , trim(magn_geometry))      

       if (INDEX(magn_geometry,'slab').eq.1) then
          call attach(fidall_param_h5, "/input", "shat   "  , shat) 
          call attach(fidall_param_h5, "/input", "parscale" , parscale)
          if (trpeps.gt.0.0) then
             call attach(fidall_param_h5, "/input", "trpeps"   , trpeps)             
          end if
       else
          if ((x_local).or.(.not.mag_prof)) then
             call attach(fidall_param_h5, "/input", "q0"       , q0)
             call attach(fidall_param_h5, "/input", "shat   "  , shat)
          endif
          if ((magn_geometry .eq. 'circular') .and. &
            (abs(sum(Bprof_coeffs)) .gt. epsilon(1.0))) then
            call putarr(fidall_param_h5, "/input/Bprof_coeffs", Bprof_coeffs, "Bprof_coeffs")
          end if
          if((magn_geometry.eq.'s_alpha').or.(magn_geometry.eq.'s_alpha_B').or.&
               &(magn_geometry.eq.'eq_cpo').or.(magn_geometry.eq.'circular')) THEN
             call attach(fidall_param_h5, "/input", "trpeps"   , trpeps) 
             call attach(fidall_param_h5, "/input", "major_R"  , major_R)                
             if ((amhd.gt.0.0).and.(magn_geometry.eq.'s_alpha').or.&
                  &(magn_geometry.eq.'s_alpha_B')) &
                  call attach(fidall_param_h5, "/input", "amhd"     , amhd)
             if (.not.x_local) call attach(fidall_param_h5, "/input", "minor_r"  , minor_r)            
             if (mag_prof) then
                call putarr(fidall_param_h5, "/input/q_coeffs"    , q_coeffs, "q_coeffs")
             end if
          elseif (INDEX(magn_geometry,'miller').eq.1) then
             call attach(fidall_param_h5, "/input", "amhd", amhd)
             call attach(fidall_param_h5, "/input", "major_R", major_R)
             call attach(fidall_param_h5, "/input", "major_Z", major_Z)
             call attach(fidall_param_h5, "/input", "minor_r", minor_r)
             call attach(fidall_param_h5, "/input", "trpeps",trpeps)
             call attach(fidall_param_h5, "/input", "kappa",kappa)
             call attach(fidall_param_h5, "/input", "delta",delta)
             call attach(fidall_param_h5, "/input", "zeta",zeta)
             call attach(fidall_param_h5, "/input", "s_kappa",s_kappa)
             call attach(fidall_param_h5, "/input", "s_delta",s_delta)
             call attach(fidall_param_h5, "/input", "s_zeta",s_zeta)
             call attach(fidall_param_h5, "/input", "drR",drR)
             call attach(fidall_param_h5, "/input", "drZ",drZ)
          else
             if (geomdir.ne.'./') call attach(fidall_param_h5, "/input", "geomdir"  , trim(geomdir))
             if (geomfile.ne.'')  call attach(fidall_param_h5, "/input", "geomfile" , trim(geomfile))

             if(magn_geometry.eq.'chease') THEN
                call attach(fidall_param_h5, "/input", "x_def"    , x_def)
                if (x_local) call attach(fidall_param_h5, "/input", "flux_pos" , flux_pos) 
                call attach(fidall_param_h5, "/input", "major_R"  , major_R)
             elseif(magn_geometry.eq.'tracer_efit') then
                if (edge_opt.gt.0) call attach(fidall_param_h5, "/input", "edge_opt" , edge_opt)
             endif
             if (force_x) call attach(fidall_param_h5, "/input", "force_x",force_x)
             if (force_z) call attach(fidall_param_h5, "/input", "force_z",force_z)
             call attach(fidall_param_h5, "/input", "minor_r"  , minor_r)
          end if
       endif
       if (mag_prof) call attach(fidall_param_h5, "/input", "mag_prof" , mag_prof)
       if (limit_shat) call attach(fidall_param_h5, "/input", "limit_shat" , limit_shat)
       if (rhostar.gt.0.0) call attach(fidall_param_h5, "/input", "rhostar", rhostar)

       call attach(fidall_param_h5, "/input", "dpdx_term", trim(dpdx_term))
       call attach(fidall_param_h5, "/input", "dpdx_pm", dpdx_pm)

       call attach(fidall_param_h5, "/input", "norm_flux_projection",&
            &norm_flux_projection)

    !****** species namelist(s) ********
       do n=0,n_spec-1
          call attach(fidall_param_h5, "/input", "name_"//trim(spec(n)%name), spec(n)%name)            
          if (spec(n)%passive) &
               &call attach(fidall_param_h5, "/input", "passive_"   , spec(n)%passive)            
          if ((x_local).or.(spec(n)%prof_type.eq.0)) then
             call attach(fidall_param_h5, "/input", "omn_"//trim(spec(n)%name) , spec(n)%omn)            
             call attach(fidall_param_h5, "/input", "omt_"//trim(spec(n)%name) , spec(n)%omt)
             call attach(fidall_param_h5, "/input", "mass_"//trim(spec(n)%name), spec(n)%mass)
             call attach(fidall_param_h5, "/input", "temp_"//trim(spec(n)%name), spec(n)%temp)
             call attach(fidall_param_h5, "/input", "dens_"//trim(spec(n)%name), spec(n)%dens)
             call attach(fidall_param_h5, "/input", "charge_"//trim(spec(n)%name)    , spec(n)%charge)
            
             if (.not.x_local) call attach(fidall_param_h5, "/input", "prof_type" , 0)
          else             
             call attach(fidall_param_h5, "/input", "prof_type_"//trim(spec(n)%name) , spec(n)%prof_type)            
             call attach(fidall_param_h5, "/input", "kappa_T_"//trim(spec(n)%name)   , spec(n)%kappa_T)            
             call attach(fidall_param_h5, "/input", "LT_center_"//trim(spec(n)%name) , spec(n)%LT_center)            
             call attach(fidall_param_h5, "/input", "LT_width_"//trim(spec(n)%name)  , spec(n)%Ln_width)            
             call attach(fidall_param_h5, "/input", "kappa_n_"//trim(spec(n)%name)   , spec(n)%kappa_n)            
             call attach(fidall_param_h5, "/input", "Ln_center_"//trim(spec(n)%name) , spec(n)%Ln_center)            
             call attach(fidall_param_h5, "/input", "Ln_width_"//trim(spec(n)%name)  , spec(n)%Ln_width)            
             if ((spec(n)%prof_type.eq.3).or.(spec(n)%prof_type.eq.4)) then
                call attach(fidall_param_h5, "/input", "delta_x_T_"//trim(spec(n)%name) , spec(n)%delta_x_T)            
                call attach(fidall_param_h5, "/input", "delta_x_n_"//trim(spec(n)%name) , spec(n)%delta_x_n)            
             end if
             if ((SUM(spec(n)%src_amp).gt.0).and.(SUM(spec(n)%src_prof_type).gt.0)) then
                call putarr(fidall_param_h5, "/input/src_prof_type_"//trim(spec(n)%name), spec(n)%src_prof_type)
                call putarr(fidall_param_h5, "/input/src_amp_"//trim(spec(n)%name)      , spec(n)%src_amp)
                call putarr(fidall_param_h5, "/input/src_width_"//trim(spec(n)%name)    , spec(n)%src_width)
                call putarr(fidall_param_h5, "/input/src_x0_"//trim(spec(n)%name)       , spec(n)%src_x0)
             end if
          end if
       end do

    !****** information namelist ********
       if (read_checkpoint.and.(checkpoint_in.ne.''))&
            &call attach(fidall_param_h5, "/input", "chpt_in", TRIM(checkpoint_in))

       if ((magn_geometry.ne.'s_alpha').and.(magn_geometry.ne.'s_alpha_B').and.&
            &(magn_geometry.ne.'circular')) &
            call attach(fidall_param_h5, "/input", "q0" , q0)

       select case(comp_type)
       case('EV')
          call attach(fidall_param_h5, "/input", "iterations for EV calculation" , it_ev)            
          call attach(fidall_param_h5, "/input", "time for eigenvalue solver"    , time_ev)            
       case('IV')
          if (itime.ne.1) call attach(fidall_param_h5, "/input", "step_time" , time_iv/(itime-1))            
          call attach(fidall_param_h5, "/input", "number of computed time steps" , itime-1)            
          call attach(fidall_param_h5, "/input", "time for initial value solver" , time_iv)            
          call attach(fidall_param_h5, "/input", "calc_dt"   , calc_dt)            
          if (precomp_nc) call attach(fidall_param_h5, "/input", "time for nc equilibrium solver" , time_nc)
       case('NC')
          call attach(fidall_param_h5, "/input", "time for nc equilibrium solver" , time_nc)            
       end select
       call attach(fidall_param_h5, "/input", "initime"   , initime)            
       call attach(fidall_param_h5, "/input", "n_fields"  , n_fields)            
       call attach(fidall_param_h5, "/input", "n_moms"    , n_moms)           
       !call attach(fidall_param_h5, "/input", "n_energies", n_energies)             
       call attach(fidall_param_h5, "/input", "lx" , lx) 
       call attach(fidall_param_h5, "/input", "ly" , ly)            
       if (n_pol.gt.1) call attach(fidall_param_h5, "/input", "n_pol", n_pol)            
       call attach(fidall_param_h5, "/input", "PRECISION" , prec)            

       call attach(fidall_param_h5, "/input", "OMP_NUM_THREADS" , omp_num_threads)            
       call attach(fidall_param_h5, "/input", "SVN_REV" , svn_rev)            
       call attach(fidall_param_h5, "/input", "RELEASE" , release)            

       if ((equil_par_curr) .and. (x_local)) &
            &call attach(fidall_param_h5, "/input", "currdens_par", currdens_par_scal)            

    !****** unit namelist ********
       call attach(fidall_param_h5, "/input", "Bref" , Bref)            
       call attach(fidall_param_h5, "/input", "Tref" , Tref)            
       call attach(fidall_param_h5, "/input", "nref" , nref)            
       call attach(fidall_param_h5, "/input", "Lref" , Lref)           
       call attach(fidall_param_h5, "/input", "mref" , mref)            
       call attach(fidall_param_h5, "/input", "omegatorref" , omegatorref)

       call creatg(fidall_param_h5, '/data')
       call putfile(fidall_param_h5, "/data/parameters_in",'parameters')

    call closef(fidall_param_h5)
 end if
#endif

  End Subroutine write_parameters

  !>Converts an integer array with a variable length (defined by nonzero entries) to a string
  Subroutine variable_intarr2strarr(intarr,strarr)
    integer, dimension(0:),intent(in) :: intarr !<Integer array
    character(len=*), intent(inout) :: strarr !<Character array
    integer :: ind_end, ind

    ind_end=SIZE(intarr)-1

    Write(strarr,"(I3)") intarr(0)
    do ind=1,ind_end
       Write(strarr,"(2A,I3)") TRIM(strarr),',',intarr(ind)
    enddo

    strarr=TRIM(strarr)

  End subroutine variable_intarr2strarr


  !>Converts a float array with a variable length (defined by nonzero entries) to a string
  Subroutine variable_fltarr2strarr(fltarr,strarr,skipval_in)
    real, dimension(0:),intent(in) :: fltarr !<Float array
    character(len=*), intent(inout) :: strarr !<Character array
    real, intent(in), optional :: skipval_in
    integer :: ind_end, ind
    real :: skipval

    if (present(skipval_in)) then
       skipval = skipval_in
    else
       skipval = 0.0
    endif

    ind_end=SIZE(fltarr)-1
    do while ((abs(fltarr(ind_end)-skipval).lt.10.*EPSILON(fltarr(0))).and.(ind_end.gt.0))
       ind_end=ind_end-1
    enddo

    Write(strarr,"(G16.8)") fltarr(0)
    do ind=1,ind_end
       Write(strarr,"(2A,G16.8)") TRIM(strarr),',',fltarr(ind)
    enddo

    strarr=TRIM(strarr)

  End subroutine variable_fltarr2strarr

  !>Converts a complex array with a variable length (defined by nonzero entries) to a string
  Subroutine variable_cplxarr2strarr(cplxarr,strarr,skipval_in)
    complex, dimension(0:),intent(in) :: cplxarr !<Complex array
    character(len=*), intent(inout) :: strarr !<Character array
    real, intent(in), optional :: skipval_in
    integer :: ind_end, ind
    real :: skipval

    if (present(skipval_in)) then
       skipval = skipval_in
    else
       skipval = 0.0
    endif

    ind_end=SIZE(cplxarr)-1
    do while ((abs(cplxarr(ind_end)-skipval).lt.10.*EPSILON(abs(cplxarr(0)))).and.(ind_end.gt.0))
       ind_end=ind_end-1
    enddo

    Write(strarr,"(A,G16.8,A,G16.8,A)") "(",real(cplxarr(0)),",",aimag(cplxarr(0)),")"
    do ind=1,ind_end
       Write(strarr,"(2A,G16.8)") TRIM(strarr),", (",real(cplxarr(ind))
       Write(strarr,"(2A,G16.8,A)") TRIM(strarr),",",aimag(cplxarr(ind)),")"
    enddo

    strarr=TRIM(strarr)

  End subroutine variable_cplxarr2strarr



end module parameters_IO
