#include "redef.h"
#ifdef WITH_VTUNE_API
#define ITT_FRAME_BEGIN(x) call itt_frame_begin(x,C_NULL_PTR)
#define ITT_FRAME_END(x) call itt_frame_end(x,C_NULL_PTR)
#else
#define ITT_FRAME_BEGIN(x)
#define ITT_FRAME_END(x)
#endif

!>Contains the routine for (linear and nonlinear) initial value computations with GENE
module initial_value_comp
  use par_mod
  use parameters_IO, only: write_parameters
  use diagnostics
  use diagnostics_df
  use diag_nl_eigenvalues
  use communications
  use time_scheme
  use aux_fields
  use checkpoint
  use reset_mod
  use initcond
  use phys_ini
  use external_contr, only: ExB_stime
  USE file_io, only: stop_file_found, avgflux_stime_file_found, &
       ExB_stime_file_found
  use convergence_monitoring
  use time_averages
  use sources_mod, only: coeff_g1_en_source, g1_init, simtime_start,&
      delivered_g1ens_msg, initialize_g1_source, finalize_g1_source
  use diagnostics_svd, only: svd_projection
  use antenna, only: calc_antenna_evolution, antenna_type
#ifdef WITH_VTUNE_API
  use itt_api
#endif
  use,intrinsic :: iso_c_binding
  use discretization_adptv_module, only: initialize_discretization_adptv,&
       finalize_discretization_adptv, is_grid_adptv

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

  implicit none

#if defined(WITH_MIC_NONLIN) && defined(WITH_OFFLOAD) && defined(WITHPERF)
  !DEC$ attributes offload:mic :: perfon,perfoff,perf_context_start,perf_context_end
#endif

  public:: check_initial_value, initial_value 
  public:: rescale
  private

  logical:: rescale=.true.
#ifdef WITH_VTUNE_API
  type(C_PTR) :: itt_domain
#endif

#ifdef WITH_CUDA_NONLIN
  interface
     subroutine start_cuda_profiling() bind(c)
     end subroutine start_cuda_profiling
     subroutine end_cuda_profiling() bind(c)
     end subroutine end_cuda_profiling
  end interface
#endif
contains

  subroutine check_initial_value
    if (nonlinear) rescale=.false.
  end subroutine check_initial_value

  !>Starts an initial value computation
  subroutine initial_value

    !Local variables
    double precision :: time1, time2
    LOGICAL,parameter :: OUTPUT=.false.
    logical :: reset=.false.
    LOGICAL :: overflow_exit, underflow_exit
    INTEGER :: ierr
    real :: local_sum, global_sum
#ifdef WITH_VTUNE_API
    character(len=10) :: domain_name
#endif
    EPIK_USER_REG(timeloop,"timeloop")

    PERF_CONTEXT_START("initialization")
    call check_initial_value

    if (dt_max.le.0.0) then
       Write(*,*)
       Write(*,"(A)") "!!! SKIPPING time loop due to ill defined dt_max !!!"
       return
    endif

    overflow_exit = .false.
    underflow_exit = .false.
    convergence_exit=.false.
    
    !**********************************************************************!
    !****************************  time loop  *****************************!
    !**********************************************************************!

    final_init=.true.
    call initialize_current_parall
    call initialize_discretization_adptv
    call initialize_calc_aux_fields
    call initialize_g1_source
    delivered_g1ens_msg = .false.

    !initial condition
    If (read_checkpoint) Then
       ! for Krook-type energy source, obtain initial condition, store
       ! before overwriting g_1 with checkpoint data
       IF (coeff_g1_en_source .GT. 0.0) THEN
         IF (mype .EQ. 0) PRINT*, 'setting g1_init for g1_en_source before checkpoint_read'
         CALL init_g_1
         if (is_grid_adptv) call nullify_g1_outside_adptv
         IF (p_has_0_mode) g1_init(li1:li2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2) = &
           g_1(li1:li2,lj1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
       END IF
#ifdef COMBI_MGR
       !mh read checkpoint from memory
       call checkpoint_read_memory(g_1, li1, li2, lj1, lj2, lk1, lk2, ll1, ll2, lm1, lm2, ln1, ln2)
#else
       call initialize_checkpoint_read
       call checkpoint_read(g_1)
       call finalize_checkpoint_read
#endif
    Else
       Call init_g_1
       if (is_grid_adptv) call nullify_g1_outside_adptv
       dt = dt_max        !set initial dt to the linear limit
    Endif

    !now that dt is initialized, recalculate dApar/dt for the antenna
    if (antenna_type.eq.3) call calc_antenna_evolution

    IF (OUTPUT) WRITE(*,"(I3,A,ES20.12)") &
         &mype,": start of initial_value: g_1 = ", DBLE(SUM(g_1*CONJG(g_1)))

    Call initialize_timescheme

    simtime_start = time

    !output the initial condition
    do_energy=(istep_energy.gt.0.or.istep_energy3d.gt.0)
    call initialize_all_diags
    call initialize_all_diags_nlev
    if(svd_proj) then
      call svd_projection
    else 
      call exec_all_diags(0,0.0,overflow_exit,underflow_exit,reset)
      call exec_all_diags_nlev(0,0.0,.true.)
    end if
    
    CALL mpi_barrier(comm_cart,ierr)
   
    If ((mype==0).and.(print_ini_msg)) then 
       Write(*,*)
       Write(*,"(A)") "*** entering time loop ***"
       WRITE(*,"(A,I10,A)") "maximal ",ntimesteps," timesteps."
    End if

    Call get_systime(time1)

#ifdef WITH_CUDA_NONLIN
    call start_cuda_profiling
#endif
    PERF_CONTEXT_END
#if defined(WITH_MIC_NONLIN) && defined(WITH_OFFLOAD) && defined(WITHPERF)
    !DEC$ offload target(mic:WHICHMIC)
    PERF_CONTEXT_START('mic_tlp')
#endif

    LIKWID_ON('timeloop')
    PERF_CONTEXT_START('timeloop')
    PERFON('t_loop')
    EPIK_USER_START(timeloop)
    VTTRACEON()
#ifdef WITH_VTUNE_API
    domain_name = "mydom"//C_NULL_CHAR
    itt_domain = itt_domain_create(C_LOC(domain_name))
#endif
    !CALL mt_trace_start_()
    Do itime = 1, ntimesteps
       
       Call get_systime(time2)
       
       If (time2<realtime_start) time2=time2+86000.0D0
       time_iv = time2 - realtime_start
       
       IF (time_iv > timelim) Exit
       IF (time > simtimelim) THEN
          IF ((mype.EQ.0).and.(print_ini_msg)) &
               & WRITE(*,"(A,F8.2,A)") "Simulation time limit of ",&
               & simtimelim," L_ref/c_ref exceeded, exiting time loop"
          EXIT
       ENDIF
       
       ITT_FRAME_BEGIN(itt_domain)
       CALL calc_timestep
       ITT_FRAME_END(itt_domain)
       
       time = time + dt
       
       if ((mype==0).and.(itime==1)) call write_parameters !called here
       !for simulations with f0 feedback during the time loop which 
       !resets itime to 1

       IF (OUTPUT) THEN
          CALL calculate_test_sum(f_(:,:,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),local_sum, global_sum)
          IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "f_ = ",global_sum
       END IF

       call exec_all_diags(itime,time,overflow_exit,underflow_exit,reset)
       !nl_eigenvalues cannot be used by diag.F90 because of circular dependency with slepc.aux
       call exec_all_diags_nlev(itime,time,.false.)

       !reset the distribution function and all n,T profile dependent terms
       IF ((.not.x_local).and.reset) THEN
          IF ((mype==0).and.(print_ini_msg)) WRITE(*,"(A,F8.2)") &
               & "Resetting g1 and profile dep. terms at t = ", time
          call reset_g1_prof(g_1,new_profs)
          if(is_grid_adptv) call nullify_g1_outside_adptv
          reset=.false.
       ENDIF
       
       IF (convergence_exit) THEN
          ! a last output for diagnostics
          call exec_all_diags(itime,time,overflow_exit,underflow_exit,reset)
          call exec_all_diags_nlev(itime,time,.false.)

          IF ((mype .EQ. 0).and.print_ini_msg) WRITE(*,"(A)") &
               &"Linear growth rate is converged, exiting time loop"
          EXIT
       END IF
       
       ! every 100 time steps GENE looks for a file 'GENE.stop' in the 
       ! directory where the executable resides and exits in a controlled way
       ! Note: checkpoints are written to s_checkpoint!!
       If (Modulo(itime,100) == 0) Then
          If (stop_file_found(file_extension,simtimelim)) Then
             if (mype==0) WRITE(*,'(A)') 'Exit due to GENE.stop file.'
             Exit
          End If
          if (avgflux_stime_file_found(file_extension,time,avgflux_stime)) then
             call reset_avgflux
             if (mype==0) WRITE(*,'(A,F13.6)') 'Set start for time averaging to ',&
                  & avgflux_stime             
          endif
          if (ExB_stime_file_found(file_extension,time,ExB_stime)) then
             if (mype==0) WRITE(*,'(A,F13.6)') 'Set start for ExB shear flow to ',&
                  & ExB_stime
          endif
       End If
       
       ! take some value of phi and check for overflow/underflow
       IF (.not.nonlinear.and.(overflow_exit).and.(rescale)) &
            &call rescale_g_1(overflow_exit)
       !rescale_g_1 may change overflow_exit
       IF ((overflow_exit).OR.(underflow_exit)) THEN
          ! a last output for debugging reasons
          call exec_all_diags(itime,time,overflow_exit,underflow_exit,reset)
          call exec_all_diags_nlev(itime,time,.false.)
          
          IF ((mype .EQ. 0).and.print_ini_msg) THEN 
             IF (overflow_exit) THEN
                WRITE(*,"(A)") "Exit due to overflow error."
             ELSE
                WRITE(*,"(A)") "Exit due to reaching the underflow limit."
             ENDIF
          ENDIF
          write_checkpoint = .false.
          EXIT
       ENDIF
    END DO !time loop
    !CALL mt_trace_stop_()
    VTTRACEOFF()
    EPIK_USER_END(timeloop)
    PERFOFF
    PERF_CONTEXT_END
    LIKWID_OFF('timeloop')
#if defined(WITH_MIC_NONLIN) && defined(WITH_OFFLOAD) && defined(WITHPERF)
    !DEC$ offload target(mic:WHICHMIC)
    PERF_CONTEXT_END
#endif
#ifdef WITH_CUDA_NONLIN
    call end_cuda_profiling
#endif

    Call get_systime(time2) 
    
    Call finalize_timescheme
    call finalize_all_diags
    call finalize_all_diags_nlev

    CALL finalize_g1_source

    call finalize_calc_aux_fields
    call finalize_current_parall

    ! finalize adaptive discretization module
    call finalize_discretization_adptv
    
    time_iv = time2-time1
    If ((mype==0).and.print_ini_msg) Then
       Write(*,"(A,F10.3,a)") "Total time of initial value computation: ",time_iv, " sec"
       Write(*,"(A,I6,a)") "Computed ",itime-1, " time steps"
       If (itime-1.gt.0) Write(*,"(a,f10.3,a)") "Time per time step:", time_iv/(itime-1), " sec"
    End If
    !overall time
    If (time2<realtime_start) time2=time2+86000.0D0
    
  end subroutine initial_value
end module initial_value_comp
