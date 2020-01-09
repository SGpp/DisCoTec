!#include "epik_user.inc"
#include "redef.h"
module gene_subroutine
  use par_mod
  use communications
  use perf_opt
  use eigenvalue_comp
  use initial_value_comp
  use compute_dt
  use parameters_IO
  use fullmatrix_aux
  use neo_equilibrium
  use discretization_adptv_module ! TODO remove this module after testing
  !the following modules are mainly required for external
  !data access:
  use diagnostics, only: cat_output
  use time_averages, only: avgfluxes, avgflux_stime
  use diagnostics_df, only: avgprof_stime, avgprof, istep_prof
  use diagnostics_neoclass, only: neoflux
  use checkpoint, only: checkpoint_in, read_checkpoint,&
       &reset_chpt_time,chptdir, chpt_read_h5
  use geometry, only: magn_geometry, dvdx, sqrtgxx_fs,&
       &minor_r, rhostar, q0, shat, trpeps, Lref, Bref,&
       &read_tracer_namelist
  use profiles_mod, only: mref, nref, Tref, unset_in_prof_defaults
  use collisions, only: coll
  use diagnostics_fmsvd

  Implicit None

! --- allow 'external' manipulation of variables ---
!!  private 
!!  public:: rungene, simtime_start, time

contains

  Subroutine rungene(MPI_COMM_SIM, in_dir, in_file_ex, chpt_in)

    integer :: MPI_COMM_SIM
    character(Len=128):: in_dir
    character(len=20):: in_file_ex
    character(len=*):: chpt_in
    double precision :: realtime
    logical :: sav_print_ini_msg
    !EXTERNAL :: mt_trace_start

    PERFON('gsub')
    
    par_in_dir=in_dir
    file_extension=in_file_ex

    !POMP$ INST INIT

    sav_print_ini_msg= print_ini_msg !save value which might be set from calling program

    if(in_file_ex.eq.'.dat') then
       print_ini_msg = .true.
    else
       !switch off messages for more than one parameter file
       print_ini_msg = .false.
    end if

    Call initialize_comm_sim(MPI_COMM_SIM)        ! Init constants for processor array
 
    IF (my_mpi_comm_world_member) THEN !inactive processors skip everything

       PERF_CONTEXT_START('autopar')
       PERFON('gini')
       Call get_systime(realtime_start)

       if (par_in_dir.ne.'skip_parfile') call read_parameters(par_in_dir)

       !checkpoint
       if(chpt_in.eq.'no') then
          read_checkpoint=.false.
          checkpoint_in=''
       elseif (chpt_in.ne.'') then
          checkpoint_in=chpt_in
       !else take extension being set, e.g., in parameters_IO
       end if
       if ((mype == 0).and.(print_ini_msg)) then
          if (.not.x_local) then
             WRITE(*,"(A)") "*****************************************************"
             WRITE(*,"(A)") "*********** nonlocal in radial direction ************"
             WRITE(*,"(A)") "*****************************************************"
             WRITE(*,*)   
          elseif (.not.y_local) then
             WRITE(*,"(A)") "*****************************************************"
             WRITE(*,"(A)") "*********** nonlocal in toroidal direction ************"
             WRITE(*,"(A)") "*****************************************************"
             WRITE(*,*)   
          endif
       endif
     
       if ((mype==0).and.(.not.cat_output)) call write_parameters

       call my_barrier

       ! todo rewrite this in a nice way
       if (is_grid_adptv) then
          is_grid_adptv = .FALSE.
          call optimize_parall_perf
          is_grid_adptv = .TRUE.
       else
          call optimize_parall_perf
       endif

       print_ini_msg = sav_print_ini_msg

       Call get_systime(realtime) 
       initime = realtime-realtime_start

       PERFOFF
       PERF_CONTEXT_END

       select case(comp_type) 
       case('EV') !eigenvalue computation          
          call comp_eigenvalues
          
       case('IV') !initial value computation
          call set_dt_max
          if (precomp_nc) call initialize_g1_nc
          call initial_value
          if (precomp_nc) call finalize_g1_nc

       case('NC')
          call comp_nc

       case('MO') !matrix output only
          call output_operator

       case('SV') !matrix output only
          if(mype==0) write(*,*) "Starting field/moment SVD."
          call svd_projection_fmom 
          
       end select

       If (mype == 0) Call write_parameters

       Call get_systime(realtime) 
       realtime = realtime - realtime_start
       If ((mype == 0).and.print_ini_msg) &
            &Write(*,"(A,F10.3,a)") "Time for GENE simulation: ",realtime, " sec"

       if (par_in_dir.ne.'skip_parfile') then
          deallocate(spec)
          call unset_in_prof_defaults()
       endif

    ENDIF !my_mpi_comm_world_member?

    Call finalize_comm_sim

    PERFOFF  !gene_subroutine

  End Subroutine rungene

end module gene_subroutine
