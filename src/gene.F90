#include "redef.h"
!>Program to start several gene runs simultaneously. 
!!\todo embellish/sort the standard output (all parallel genes write into the same file at the moment!)
Program gene
  use mpi
  use gene_scan
  use communications
  use par_mod
  use parameters_IO
  use check_parameters
  use gene_subroutine
  use file_io
  use pinning_mod
  use diagnostics, only: NRGFILE, cat_output
  use, intrinsic :: iso_c_binding
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif
#ifdef MEMORY_CHECKER
  use memory_checker
#endif
#ifdef COMBI_MGR
  use perf_opt
  use initial_value_comp
  use checkpoint
#endif
  implicit none

#if defined(WITH_MIC_NONLIN) && defined(WITH_OFFLOAD) && defined(WITHPERF)
  !DEC$ attributes offload:mic :: perfon,perfoff,perfinit,perfout
#endif

  integer :: ierr
#ifdef MEMORY_CHECKER
  integer(kind=c_long) :: max_allocated,still_allocated
  integer(kind=c_int) :: unfreed_blocks
#endif
  double precision:: wtime_start
  logical:: mult_par=.true.
  integer:: i_proc,i_group
  logical:: performance_output=.false.
  integer:: gene_comm, comm_parall
  character(len=20):: f_ext
  character(len=128):: ch_in

#ifdef COMBI_MGR
  double precision:: wtime
  CHARACTER(len=1024) :: cwd
  integer:: worker_signal=0
  double precision :: worker_tperf
#endif

  !set svn_rev here as par_other.F90 might not be recompiled after
  !code modifications (gene_subroutine, however, almost always is)
#if defined(SVN_REV)
  svn_rev = SVN_REV
#endif
  
#ifdef WITHOMP
  call mpi_init_thread(MPI_THREAD_MULTIPLE, omp_level,ierr)
#ifdef __INTEL_COMPILER
  call pinning()
#endif
#else
  call mpi_init(ierr)
#ifdef COMBI_MGR
  call init_stats
#endif
  omp_level = MPI_THREAD_SINGLE
#endif
  if (ierr /= 0) stop 'mpi_init failed!'

  call check_for_scan(mult_par)

  call initialize_comm_scan(gene_comm,comm_parall)

#ifdef COMBI_MGR
  do while(.true.)
    call worker_wait(gene_comm,worker_signal,n_procs_sim,n_parallel_sims)

    ! reset timers
    time_perf = -1.0
    time_cp = -1.0

	! check status code
	! 0 = run first
	! do nothing

	! 1 = run next
	! do nothing

    ! 2 = eval
    if(worker_signal.eq.2) cycle

    ! 3 = grid_eval
    if(worker_signal.eq.3) cycle

    ! 4 = combine
    if(worker_signal.eq.4) cycle

    ! 5 = exit
    if(worker_signal.eq.5) exit

    ! 6 = sync
    if(worker_signal.eq.6) cycle

    ! 7 = test_conversion
    if(worker_signal.eq.7) cycle

    ! 8 = combine_fg
    if(worker_signal.eq.8) cycle

    ! 9 = calc ev
    ! do nothing

    ! 10 = ev_calc_init
    ! do nothing

    ! 11 = update_coeffs
    if(worker_signal.eq.11) cycle
    
    ! 12 = update combi parameters
    if(worker_signal.eq.12) cycle
    
    ! 13 = add_task
    if(worker_signal.eq.13) cycle
    
    ! 14 = recomputed
    ! do nothing

    ! 17 = parallel eval
    if(worker_signal.eq.17) cycle

    call gene_time_start
#endif

    LIKWID_INIT
#if defined(WITH_MIC_NONLIN) && defined(WITH_OFFLOAD) && defined(WITHPERF)
    !DEC$ offload begin target(mic:WHICHMIC)
    PERFINIT
    PERFON('mmain')
    !DEC$ end offload
#endif

  PERFINIT         ! Initialize performance monitoring
  PERFON('GENE')

  if (mype_gl.eq.0) then
     write(*,*)
     write(*,"(A)") "*****************************************************"
     write(*,"(3A)") "***** This is GENE 11 (release ",&
          trim(release),") *******"
     if ((svn_rev.ne.'').and.(svn_rev.ne.'exported')) &
          write(*,"(A,A15,A)") "*********** svn-revision: ", svn_rev, " ***********"
     write(*,"(A)") "*****************************************************"
     write(*,*)
  endif
  
  wtime_start=MPI_Wtime() 

#ifdef COMBI_MGR
    !mh: modified to use splitted comm
    call check_for_diagdir(gene_comm)
#else
	 call check_for_diagdir
#endif

    if(.not.mult_par) then
     !only one gene simulation
     f_ext='.dat'
     ch_in=''
     call rungene(gene_comm, par_in_dir, f_ext, ch_in)
  else
     !list of gene simulations
     call run_scan(gene_comm, comm_parall)
  end if
  call erase_stop_file
  call create_finished_file

#ifdef COMBI_MGR
    wtime = MPI_Wtime()-wtime_start
#endif

    if (mype_gl.eq.0) then
     write(*,"(A,F10.3,A)") "Total wallclock time for GENE: ", MPI_Wtime()-wtime_start, " sec"
     if(mult_par) write(*,"(A,F10.3,A)") "Percentage of idle time at the end of the scan due to load imbalancing: ", &
          idle_global/n_parallel_sims/n_procs_sim/(MPI_Wtime()-wtime_start)*100, " %"
  end if

  PERFOFF !GENE

#if defined(WITH_MIC_NONLIN) && defined(WITH_OFFLOAD) && defined(WITHPERF)
  !DEC$ offload begin target(mic:WHICHMIC)
  PERFOFF
  PERFOUT('mmain')
  !DEC$ end offload
#endif

  ! output performance data
  if (performance_output.and.(n_procs_sim.le.64)) then
     do i_group=0,n_parallel_sims-1
        do i_proc=0,n_procs_sim-1
           call my_barrier
           if ((i_proc.eq.mype).and.(i_group.eq.my_sim)) then
              write(*,"(A,I4)") "Process ",mype
              PERFOUT('t_loop')
              !close(TEMPFILE)
             call flush(6)
           end if
        end do
     end do
  else
     if (mype_gl.eq.0) then
        PERFOUT('GENE') 
     endif
  end if
  LIKWID_CLOSE

#ifdef COMBI_MGR
    call gene_time_stop
    call worker_ready(wtime, time_perf, time_iv, time_cp)
  end do
#endif

  call finalize_comm_scan(gene_comm,comm_parall)

#ifdef COMBI_MGR
  call finalize_stats
#endif

  call mpi_finalize(ierr)

#ifdef MEMORY_CHECKER
  call get_malloc_stat(max_allocated,unfreed_blocks,still_allocated)
  write(*,"(I7,A,I10,A)") mype,": max_allocated = ",max_allocated," bytes"
#endif
End Program gene
