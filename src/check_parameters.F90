#include "redef.h"
#include "switches.h"

!>Routines to check the consistency of parameters
module check_parameters
  
  use par_in
  use par_other
  use coordinates
  use discretization
  use equilibrium_fields, only: check_par_equilibrium_fields
  use geometry
!  use initial_value_comp, only: check_initial_value
  use compute_dt
  use external_contr, only: ExBrate, pfsrate, lilo_w_full_omegator_profile
  use eigen_parameters
  use petsc_precond
  use prefactors, only: check_par_prefactors
  use fourier
  use collisions, only: check_par_collisions, spacediff_off, spacediff
  use communications, only: mpi_comm_world, mpi_integer
  use file_io
  use numerical_damping
  use field_solver
  use time_averages

  implicit none
  public:: check_for_diagdir, check_par_phys

  private

contains

  !>Check if diagdir exists
#ifdef COMBI_MGR
  subroutine check_for_diagdir(split_comm)
    integer:: split_comm
#else
  subroutine check_for_diagdir
#endif
    integer:: fileunit, err, ierr
    !check if diagdir exists
    !if (mype_gl.eq.0) then
#ifdef COMBI_MGR
    if (mype.eq.0) then
#else
    if (mype_gl.eq.0) then
#endif
       call get_unit_nr(fileunit)
       open(fileunit,file=TRIM(diagdir)//"/diagdirtest",iostat=err) 
       if (err.ne.0) then
          write(*,"(A)") 'PROBLEM WITH FILE IO, PLEASE CHECK IF DIAGDIR EXISTS!'
       else
          close(fileunit,status="DELETE")
       end if
    end if
#ifdef COMBI_MGR
    call mpi_bcast(err,1,MPI_INTEGER,0,split_comm,ierr)
#else
		call mpi_bcast(err,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
    if (err.ne.0) stop

  end subroutine check_for_diagdir

  !>Check parameters defining the physics
  subroutine check_par_phys
    integer::n
    Logical :: write_pe

    write_pe = ((mype.le.0).AND.(print_ini_msg))

    do n=0,n_spec-1
       if (spec(n)%charge.eq.0) then
          stop 'Species with zero charge cannot be considered yet!'
       endif
    enddo
    if (n_spec.gt.1) then
       do n=1,n_spec-1
          if (spec(n)%prof_type.ne.spec(0)%prof_type) stop 'prof_types must agree between all species!'
       enddo
    endif

    if ((x_local).or.(lilo)) then
       do n=0,n_spec-1
          if (spec(n)%prof_type.gt.0 .and. write_pe) then
             Write(*,"(A)") "WARNING: temp and dens profiles in local runs:"
             Write(*,"(A)") "Tref nref and species temp and dens must be given at profile center position"
          endif
          spec(n)%src_prof_type = 0
       enddo
       !mag_prof = .false.
       if(mag_prof.and.write_pe) &
            Write(*,"(A)") "WARNING: mag_prof=.true. may not make sense in local runs"
    end if

    !constant ExB shear model
    if (abs(ExBrate).gt.epsilon(ExBrate)) then
       ExB=.true.
    end if

    if (ExB) then
       if ((nexc.eq.1).and.write_pe) write(*,"(A)") &
            &"WARNING: small nexc values are not recommended for the ExB shear flow model"
    end if

    if (lilo_w_full_omegator_profile.and..not.lilo) &
       stop "lilo_w_full_omegator_profile can only be used when lilo = True."

    if (lilo) then
       if (write_pe) write(*,"(A)") "WARNING: GENE is in LILO mode which is only for testing"
       if(mype == 0) write(*,*) "lilo_w_full_omegator_profile",lilo_w_full_omegator_profile
       if (lilo_w_full_omegator_profile) then
          if(rad_bc_type.ne.-1) &
              stop "lilo_w_full_omegator_profile can only be used with rad_bc_type = -1 (i.e. Dirichlet)"
          if(any(spec(:)%prof_type.ne.-2)) &
              stop "lilo_w_full_omegator_profile can only be used with iterdb profile input."
          if(ExBrate.ne.-1111.0) &
              stop "Set ExBrate = -1111.0 to use lilo_w_full_omegator_profile."
       end if
       !impose radial periodic boundary only without shifted metric or constant ExB shear
       if (.not.(shifted_metric.or.ExB).and.rad_bc_type.ne.-1) rad_bc_type=0
       !rad_bc_type=-1 for usage of Dirichlet without ExB or shifted metric
       if (rad_bc_type==-1) rad_bc_type=1
       !rhostar = 0.0001
    endif

    !warnings relating to geometry interfaces
    if (Omega0_tor.ne.0..and.magn_geometry.eq.'gist') stop "Rotation effects currently unavailable for GIST geometry"
    if (magn_geometry.eq.'gist'.and.pfsrate.ne.0) then 
       write (*,"(A)") "Warning: Parallel flow shear currently unavailable for GIST geometry, setting pfsrate=0."
       pfsrate=0.
    endif
    if (INDEX(magn_geometry,'miller').eq.1.and.any(spec(:)%prof_type.ne.0)) &
         stop "Automatic profile input currently not possible for Miller geometry"


    if (xy_local.and.nonlinear) then
       if (yx_order.and.(.not.fourier2D))then
          fourier2D=.true.
          if (write_pe) Write(*,"(A)") "Setting fourier2D for yx_order nonlinear computation!"
       endif
       if (fourier2D.and.(.not.yx_order)) then
          yx_order=.true.
          if (write_pe) Write(*,"(A)") "Setting yx_order for 2D fftw3 transforms!!"
       endif
       if (yx_order) then
          if (write_pe) Write(*,"(A)") "fourier2D does currently not work--exit!!"
          stop
       endif
    endif

    !y coordinate is always periodic
    if (.not.y_local) then
       rad_bc_type=0
       ky0_ind = 0
       yx_order=.true.
    endif

    if (x_local) shifted_metric=.false.

    !support for old inverted parameters
    if (.not.bpar_off) bpar=.true.
    if (.not.spacediff_off) spacediff=.true.

    !dpdx/bpar consistency warnings
    if (bpar.and.(dpdx_term.ne.'on').and.(dpdx_term.ne.'full_drift')) then
       if (write_pe) write(*,"(A)") "WARNING: dpdx_term = 'on' and dpdx_pm<>0 "//&
            &"are recommended if B_par fluctuations are considered!!"
    endif
    if (.not.(bpar).and.(dpdx_term.ne.'gradB_eq_curv')) then
       if (write_pe) write(*,"(A)") "WARNING: dpdx_term = 'gradB_eq_curv' and "//&
            &" dpdx_pm<>0 might be more appropriate if bpar = F"
    endif

    call check_par_collisions
    call check_par_equilibrium_fields
    call check_par_prefactors

    call check_compute_dt
    call check_eigenvalue_comp
    call check_petsc_precond
!    call check_initial_value
   
    call check_fourier(print_ini_msg)

    if (((xy_local).or.(lilo).or.(n0_global.lt.0)).and.(kymin.le.0.0).and.(.not.only_zonal).and.(y_local)) &
         STOP "kymin must be larger than zero"
    call check_field_solver


    !IF (dE_drive_out) n_energies = n_energies + 1
    !IF (dE_coll_out) n_energies = n_energies + 1

    !geometry
    if (.not.xy_local.and.rhostar.eq.0.0) stop "rhostar must be given for a global simulation"   
    
    if (.not.x_local.and.x0==-1) stop "x0 must be given for a global simulation"
    if (pfsrate.ne.0..and.x0==-1.and.magn_geometry=='tracer') stop "parallel flow shear with tracer requires x0" 
    if (magn_geometry=='tracer_efit'.and.x0==-1) stop 'x0 must be given for tracer_efit interface'
    if (x0==-1.and.any(spec(:)%prof_type.eq.-1)) stop "x0 must be given for profile file input"
    if (magn_geometry=='chease') then
       if (flux_pos.eq.0) then
          !use x0, if flux_pos is not defined
          if (x0.ne.-1) then
             flux_pos=x0
          else
             stop 'Specify flux surface coordinate x0 for CHEASE interface.'
          endif
       end if
    endif

    if ((.not.x_local).and.(.not.lilo).and.(.not.mag_prof)) &
         & stop "mag_prof = .false. does currently not work"

    if ((x_local.or.(.not.mag_prof)).and.(q0.eq.-1111.0)) then
       if((magn_geometry.eq.'s_alpha').or.(magn_geometry.eq.'circular')) then
          if (write_pe) write(*,"(A)") "specify q0!"
          stop
       endif
    endif

    if ((parscheme.eq.'c2nd').and.(hyp_z.gt.0.0)) then
       hyp_z_order=2
      if (write_pe) write(*,"(A)") "hyp_z_order is set to 2 because of the choice of parscheme"
    endif

    if ((parscheme.eq.'c0th').and.(hyp_z.gt.0.0)) then
      if (write_pe) write(*,"(A)") "parscheme c0th works only without hyp_z"
      stop
    endif

    if ((parscheme.eq.'c0th').and.arakawa_zv) then
       arakawa_zv=.false.
       if (write_pe) write(*,"(A)") "2D simulation, switching off parallel Arakawa scheme."
    endif

    if ((parscheme.eq.'u3rd').and.(perf_vec(3).ne.4)) then
       perf_vec(3) = 4
       if (write_pe) write(*,"(A)") "WARNING: parscheme='u3rd' only "//&
            "possible with perf_vec(3)=4; modifying perf_vec accordingly"
    endif

    if ((vparscheme.eq.'u3rd').and.(perf_vec(3).ne.4)) then
       perf_vec(3) = 4
       if (write_pe) write(*,"(A)") "WARNING: vparscheme='u3rd' "//&
            "only possible with perf_vec(3)=4; modifying perf_vec accordingly"
    endif

    if (xy_local.and.(comp_type.eq.'EV').and.(evenx.eq.1)) then
       !EV computations only with the minimal number of points 
       if (mype.le.0) write(*,"(A)") 'decreasing nx0 by 1 to have all kx physical for full matrix computations'
       nx0=nx0-1
       evenx=0
    endif

    if(timescheme.eq.'') then
       if(turbdeal) then
          timescheme='turbo'
       else
          timescheme='RK4'
       end if
    end if

    if(turbdeal.and.timescheme.ne.'turbo') then
       if (write_pe) write(*,"(A)") "switching to turbo timescheme to be consistent with dealiasing scheme"
       timescheme='turbo'
    end if
    
    if((comp_type.eq.'IV').and.(timescheme.eq.'IE1s').or.(timescheme.eq.'IE1p').or.(timescheme.eq.'IE1f')) then
       if(nonlinear) then
          if(write_pe) write(*,"(A)") 'implicit and nonlinear not possible, switching to explicit RK4 time scheme'
          timescheme='RK4'
       elseif((evenx.eq.1).and.xy_local) then
          nx0=nx0-1
          evenx=0
          if(write_pe) write(*,"(A)") "reduced nx0 by one to have all kx physical for full matrix computations" 
       endif
       if(calc_dt) then
          if(write_pe) write(*,"(A)") "calculation of timestep is not necessary, calc_dt is set to .f." 
          calc_dt=.false.
       endif
    endif

    if (calc_dt) then
       !ignore dt_max from parameters file:
       dt_max = -1.0
    end if
    
    !usually, we determine automatically whether hypz_compensation is necessary,
    !however, if SLEPc is not used, hypz_opt just switches it on
    if (arakawa_zv.and.hyp_on_h) then
#ifdef WITHSLEPC
       if (.not.calc_dt) hypz_compensation=hypz_opt
#else
       hypz_compensation=hypz_opt
#endif
    endif
    if (comp_type.ne.'IV') hypz_compensation=.false.
    if (.not.hypz_opt) hypz_compensation=.false.
    !hyp_on_h is only implemented with arakawa_zv
    if (.not.arakawa_zv.and.hyp_on_h) hyp_on_h = .false.
    
    if ((which_ev.eq.'shift_invert_s').and.(evenx.eq.1).and.xy_local) then
       nx0=nx0-1
       evenx=0
       if(write_pe) write(*,"(A)") "reduced the number of x modes by one" 
    endif

#ifdef GDAGGER_G
    if ((.not.x_local).and.(rad_bc_type.eq.2)) then
       if(write_pe) then
          write(*,"(A)") "WARNING: Von-Neumann b.c. not implemented yet for GDAGGER_G"
          write(*,"(A)") "compile with #undef GDAGGER_G in switches.h"
       endif
    endif
#endif

    if (write_flux_final.gt.0) then
       if (avgflux_stime.lt.0) avgflux_stime=0
    endif

    if(diag_GyroLES) then
       if (write_pe) write(*,"(A)") "For the moment GyroLES diagnostic only works with nblocks = 1"
       nblocks =1 
    end if

  end subroutine check_par_phys

 
end module check_parameters
