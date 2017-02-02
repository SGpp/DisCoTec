#include "intrinsic_sizes.h"
#include "redef.h"
!>Performance optimization by automatic choice of parallelization and algorithms
MODULE perf_opt

  use par_mod
  use boundaries
  use communications
  use aux_fields
  use calc_rhs
  Use phys_ini
  use vel_space
  use time_scheme
  use compute_dt
  use check_parameters
  USE file_io
  use arrays
  use fieldsolver_df
  use external_contr
  use diagnostics, only: mem_est_diag
  use convergence_monitoring, only: mem_est_convergence_monitoring, istep_omega
  use petsc_precond, only: pc_type
  use mtrandom
  use axpy
  use eigen_parameters, only: which_ev
  use hybrid, only: mem_est_hybrid
  use geometry
  use aux_func
  use blockindex

  implicit none
#ifdef COMBI_MGR
  public:: optimize_parall_perf, mb_per_core, perf_tsteps, time_perf
#else
  public:: optimize_parall_perf, mb_per_core, perf_tsteps
#endif
  private

  TYPE result_type
     real :: sum_square
     real :: sum_re_plus_im
     LOGICAL :: isDefined=.FALSE.
  END TYPE result_type

  INTERFACE calculate_result
     MODULE PROCEDURE mp_calculate_result
  END INTERFACE

  INTERFACE isEqual
     MODULE PROCEDURE isEqual_result_type
  END INTERFACE

  INTERFACE isNan
     MODULE PROCEDURE isNan_result_type, isNan_real
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE assign_result_type_to_real, assign_result_type
  END INTERFACE
  
  INTEGER:: dimopts, mem_start,mem_end, arena_start, arena_end
  integer:: perf_tsteps=-1, AUTOPAR_FILE
  integer,dimension(:),allocatable:: num_opts,perf_vec_in
  logical:: more_choices, fast_opt, skip_parall_on_first_pv=.true.
  !logical :: show_info = .false.
  LOGICAL :: perf_output_to_file = .TRUE., write_optim_perf=.true.
  TYPE(result_type),SAVE :: reference_result
  real:: mb_per_core=500
  REAL:: mb_per_mpi_proc, prefac_error=1000.0
  real(8):: smallest_t1
  INTEGER :: nps1,nps2,npv1,npv2,npw1,npw2,npx1,npx2,npy1,npy2,npz1,npz2
  integer :: STDOUT=6
  logical:: first_parall,sav_print_ini_msg
  integer :: n_procs
#ifdef COMBI_MGR
  double precision :: time_perf
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!routines to test alternative implementations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>Sets the number of possible alternative routines for the current parallelization
  subroutine initialize_optimize_perf

    integer :: ierr
    !for global y, dy has to be known for the memory estimate of the gyromatrix. 
    !dy can only be computed when q0 is known, 
    !q0 is in general only known when the numerical magnetic equilibrium has been read
    if (.not.y_local) then
       call initialize_coordinates
       call initialize_geometry
       call finalize_geometry
       call finalize_coordinates
    end if

    !number of alternative routines available in GENE

    !1: dgdxy_term
    if (xy_local) then
       num_opts(1)=2
    else
       num_opts(1)=1
    endif

    !2: dchidxy_term
    if (xy_local) then
       num_opts(2)=2
    else
       num_opts(2)=1
    endif

    !3: dfdzv_terms
#ifndef oldparbc
    num_opts(3)=1
#else
    num_opts(3)=3 !4th implementation is just for diagnostic purposes
    !allow 4th implementation if selected by hand
    !or for parscheme='u3rd'
    if (perf_vec(3).EQ.4) num_opts(3)=4
#endif
    if (arakawa_zv) then
       num_opts(3)=2
    end if

    !4: collisions (see collisions.F90)
    if (collision_op.ne.'none') then
       num_opts(4)=3
    else
       num_opts(4)=1
    endif

    !5: RHS computation 
    !(obsolete perf_vec which is just kept for
    !backward compatibility)
    num_opts(5) = 1

    !6: calc_f (see aux_fields.F90)
    if (xy_local) then
       num_opts(6)=2
    else
       num_opts(6)=1
    endif
    
    !7: exchange_z (see boundary.F90)
    if (x_local) then
       if (yx_order) then
          num_opts(7)=1
       else
          num_opts(7)=2
       end if
    else
       if (n_procs_z.gt.1) then
          num_opts(7)=2
       else
          num_opts(7)=1
       endif
    endif
    
    !8: exchange_v (see boundary.F90)
    if (n_procs_v.gt.1) then
       num_opts(8)=2
    else
       num_opts(8)=1
    endif
    
    !9: moments (see field_solve_ff.F90, vel_space.F90)
    if (xy_local) then
       if(trap_pass) then
          num_opts(9)=1
       else
          num_opts(9)=2
       end if
    else
       num_opts(9)=1
    endif

    ! calculate the prefac_error, this is done as the three-sigma deviation
    ! of the binomial distribution of 
    !prefac_error = 3*sqrt(real(nx0)*real(nky0)*real(nz0)*real(nv0)*real(nw0)*n_spec/4.0)
    !print*, "Using prefac_error = ", prefac_error
    call mpi_comm_size(MY_MPI_COMM_WORLD,n_procs,ierr)
  end subroutine initialize_optimize_perf

  subroutine finalize_optimize_perf

  end subroutine finalize_optimize_perf

  !>\todo the perf routine should not be called if perf_vec is set and auto_parall is off!
  subroutine optimize_parall_perf
    INTEGER :: count=0
    double precision :: time_opt1, time_opt2
    integer :: n_procs_s_best,n_procs_v_best,n_procs_w_best,n_procs_y_best
    integer :: n_procs_x_best, n_procs_z_best, nblocks_best
    INTEGER :: arena_start, arena_end, mem_start,mem_end
    real :: elapsed_time, best_parall_time=-1.0, fastest_first_perfvec_time=-1.0
    REAL :: MFR_min, MFR_max, MFR_avg, total_mem_min,total_mem_max
    LOGICAL :: optim_perf=.TRUE., first_parallelization=.true.
    Logical :: axpy_ij_use_blas_best, opt_nblocks=.false.
    integer,dimension(:),allocatable:: perf_vec_best
    INTEGER :: i_opt, ierr
#if (MB_PER_CORE==env) || !defined(MB_PER_CORE)
    INTEGER :: strlen, err, memory_from_env
    character(len=5):: read_env
#endif

    !PERFON('opt_par_perf')
    in_perf_opt = .true.
    !neglect reference results from previous calls
    reference_result%isDefined=.false.
    first_parall=.true.
    first_parallelization=.true.
    best_parall_time = -1.0
    fastest_first_perfvec_time = -1.0

    sav_print_ini_msg=print_ini_msg  
    
    dimopts=size(perf_vec)
    allocate(perf_vec_best(dimopts), num_opts(dimopts), perf_vec_in(dimopts))

    ! check if perf_vec contains any negative or zero values
    optim_perf = .FALSE.
    DO i_opt = 1,dimopts
       IF (perf_vec(i_opt).LE.0) optim_perf=.TRUE.
    END DO

    if (nblocks.eq.0) then
       opt_nblocks=.true.
    endif

    perf_vec_in   = perf_vec

    if (perf_tsteps.lt.1) then
       perf_tsteps=3
       fast_opt=.true.
    else
       fast_opt=.false.
    endif

    call my_barrier
    if (optim_perf.or.auto_parall) call get_systime(time_opt1)

#if (MB_PER_CORE==env) || !defined(MB_PER_CORE)
    !check environment variables
    !NOTE: This function caused problems on PGI compilers on CRAY
    !      machines. Hence, environment are only read if MB_PER_CORE
    !      is not set in the makefile
    CALL get_environment_variable('MEMORY_PER_CORE',read_env,strlen,err)
    if(err.eq.0) then 
       read( read_env, '(i5)' ) memory_from_env
       mb_per_core = memory_from_env
    else
       if (mype.eq.0) then
          print*,' '
          print*,'#############################################################'
          print*,'WARNING: Memory per core neither set in makefile '
          print*,' nor environment variable MEMORY_PER_CORE is set!! Using 500 MB'  
          print*,'#############################################################'
          print*,' '
       endif
    endif
#else
    mb_per_core=MB_PER_CORE
#endif
    if (mype.eq.0) write(*,"(A,F8.2,A)") "Using a maximum of ",mb_per_core," MB per core."
    !set to high value
    smallest_t1=80000

    !take into account OMP parallelization
    mb_per_mpi_proc=omp_num_threads*mb_per_core       

    if((.not.auto_parall).and.(.not.optim_perf)) perf_output_to_file=.false.
    IF (mype.EQ.0) then
       if (perf_output_to_file) THEN
          CALL get_unit_nr(AUTOPAR_FILE)
          OPEN(AUTOPAR_FILE,file=TRIM(diagdir)//"/autopar"//trim(file_extension))
       else
          AUTOPAR_FILE=STDOUT ! write to stdout
       end if
    END IF

    IF (auto_parall) THEN
       if (pc_type.ne.'none') call set_precond_constraints
       best_parall_time=-1.
       if (n_procs_s.gt.0) then
          nps1=n_procs_s; nps2=n_procs_s
       else
          nps1=max(1,min_nps); nps2=min(n_spec,max_nps)
       endif
       if (n_procs_v.gt.0) then
          npv1=n_procs_v; npv2=n_procs_v
       else
          npv1=max(1,min_npv)
          if(pc_type.eq.'asm') then
             npv2=max(1,nv0/12)
          else
             npv2=min(nv0/2,max_npv)
          end if
       endif
       if (n_procs_w.gt.0) then
          npw1=n_procs_w; npw2=n_procs_w
       else
          ! is this true also with collisions?
          npw1=max(1,min_npw); npw2=min(nw0,max_npw)
       endif
       if (n_procs_x.gt.0) then
          npx1=n_procs_x; npx2=n_procs_x
       else
          npx1=max(1,min_npx)
          if (xy_local) then
             npx2=1
          else
             npx2 = min(nx0/4,max_npx)
             !approximate largest gyroradius
             !(at this point bfield is unknown;
             !assumption: min(bfield)=0.5)
             !rhomax = 0.0
             !do n=0, n_spec-1
             !   rhomax = MAX(rhomax,sqrt(spec(n)%mass*spec(n)%temp*2.0*lw/0.5)/&
             !   abs(spec(n)%charge))
             !enddo
             !npx2=MIN(nx0/((lag_order+1)/2),INT(lx/rhomax))+1
          endif
       endif
       if (n_procs_y.gt.0) then
          npy1=n_procs_y; npy2=n_procs_y
       else
          npy1=max(1,min_npy)
          if (nonlinear) then
             npy2=min(nky0/4+1,max_npy)
          else
             npy2=min(nky0,max_npy)
          endif
       endif

       if (n_procs_z.gt.0) then
          npz1=n_procs_z; npz2=n_procs_z
       else
          npz1=max(1,min_npz)
          if(pc_type.eq.'asm') then
             npz2=max(1,nz0/12)
          elseif (parscheme.eq.'c6th') then
             npz2=min(nz0/3,max_npz)
          else
             npz2=min(nz0/2,max_npz)
          endif
       endif

       IF ((mype.EQ.0).and.print_ini_msg) THEN
          WRITE(*,"(A)") ' '
          WRITE(*,"(A)") '------------------------------'
          WRITE(*,"(A)") 'Starting optimization routine (this may take a while) '
       end IF
       print_ini_msg = .false.
       DO WHILE (get_next_parallelization(n_procs_x,n_procs_y,n_procs_z,n_procs_v,n_procs_w,n_procs_s))
          if (check_parallelization(.false.)) then
             !write(*,"(6I4)") n_procs_x,n_procs_y,n_procs_z,n_procs_v,n_procs_w,n_procs_s
             ! here we have a valid parallelization which now has to be checked for performance
             mem_start = get_used_memory(arena_start)
             PERFON('ini_1')
             CALL split_comm
             call initialize_discretization(print_ini_msg)
             if (nblocks.ne.0) then
                if (modulo(lklmn0,nblocks).ne.0) opt_nblocks=.true.
             endif
             call initialize_optimize_perf
             PERFOFF

             !CALL est_mem_par(mem_par)
             !IF (mem_par.LT.mb_per_mpi_proc) THEN

             IF (opt_nblocks) call estimate_blocksize(nblocks)
             call initialize_axpy_ij(lij0,lklmn0/nblocks)
             IF (mype.EQ.0) THEN
                IF (first_parallelization) THEN
                   WRITE(AUTOPAR_FILE,"(A)") &
                     "parallelization:    s   v   w   x   y   z"
                ENDIF
                IF (write_optim_perf) THEN
                   WRITE(AUTOPAR_FILE,"(A,6I4,A,1I5)",advance='yes') &
                        &"parallelization: ",n_procs_s, n_procs_v, n_procs_w, &
                        &n_procs_x, n_procs_y, n_procs_z, "  nblocks: ",nblocks
                  call flush(AUTOPAR_FILE)
                ELSE
                   WRITE(AUTOPAR_FILE,"(A,6I4,A,1I5,2X)",advance='no') &
                        &"parallelization: ",n_procs_s, n_procs_v, n_procs_w, &
                        &n_procs_x, n_procs_y, n_procs_z, "  nblocks: ",nblocks             
                  call flush(AUTOPAR_FILE)
                ENDIF
             END IF
             if (optimize_perf(elapsed_time,fastest_first_perfvec_time,&
                  & MFR_min, MFR_max,MFR_avg,total_mem_min,total_mem_max)) THEN
                ! we found at least one perf_vec which worked and can output the results
                IF (mype.EQ.0) THEN
                   IF (write_optim_perf) THEN
                      WRITE(AUTOPAR_FILE,"(3X,A)") &
                        "best perf_vec      min. WCT     mem_min  mem_max deficit(byte) data segment (MB)"
                      WRITE(AUTOPAR_FILE,"(3X,9(I1,1X),A)", advance='no') perf_vec, ' '
                   ENDIF
                   IF (MFR_max.GT.1e-10) THEN
                      WRITE(AUTOPAR_FILE,"(ES11.4,2F9.2,A,3F7.1,2X)",advance='no') &
                           & elapsed_time, total_mem_min,total_mem_max,&
                           &", MFlop/s = ",MFR_min,MFR_max,MFR_avg
                   ELSE
                      WRITE(AUTOPAR_FILE,"(ES11.4,2F9.2,2X)",advance='no') &
                           &elapsed_time, total_mem_min,total_mem_max
                   END IF
                ENDIF
                
                IF (best_parall_time.LT.0.0) THEN
                   best_parall_time = elapsed_time
                   n_procs_s_best = n_procs_s
                   n_procs_v_best = n_procs_v
                   n_procs_w_best = n_procs_w
                   n_procs_x_best = n_procs_x
                   n_procs_y_best = n_procs_y
                   n_procs_z_best = n_procs_z
                   perf_vec_best  = perf_vec
                   nblocks_best   = nblocks
                   axpy_ij_use_blas_best = axpy_ij_use_blas
                ELSE
                   IF (elapsed_time.LT.best_parall_time) THEN
                      best_parall_time = elapsed_time
                      n_procs_s_best = n_procs_s
                      n_procs_v_best = n_procs_v
                      n_procs_w_best = n_procs_w
                      n_procs_y_best = n_procs_y
                      n_procs_x_best = n_procs_x
                      n_procs_z_best = n_procs_z
                      perf_vec_best  = perf_vec
                      nblocks_best   = nblocks
                      axpy_ij_use_blas_best = axpy_ij_use_blas
                   ENDIF
                   
                END IF
                count = count + 1
                call finalize_optimize_perf
                call free_comm

                mem_end = get_used_memory(arena_end)
                if (mype.eq.0) then
                   write(AUTOPAR_FILE,"(I10,F10.3)") mem_end-mem_start,real(arena_end)/1024./1024.
                   call flush(AUTOPAR_FILE)
                endif
             else
                IF (mype.EQ.0) THEN
                   if (total_mem_min.lt.0.0) then
                      ! first perf_vec was found to be slow, so we skip this parallelization
                      write(AUTOPAR_FILE,"(A,6I4,a)",advance='yes') &
                           &"skipped parallelization ",n_procs_s, &
                           & n_procs_v, n_procs_w, n_procs_x, n_procs_y, n_procs_z,&
                           ' based on speed of first perf_vec'
                   else
                      ! no perf_vec has been found which presumably fits 
                      ! into the given memory
                      write(AUTOPAR_FILE,"(A,6I4,a,F10.1,a,2X)",advance='yes') &
                           &"skipped parallelization ",n_procs_s, &
                           & n_procs_v, n_procs_w, n_procs_x, n_procs_y, n_procs_z,&
                           ' due to memory estimate of ',total_mem_min,' MB/core'
                   endif
                END IF
                call flush(AUTOPAR_FILE)
                call finalize_optimize_perf
                call free_comm
             endif
          end if
          IF (first_parallelization) first_parallelization=.FALSE.
       ENDDO

       IF (count.eq.0) THEN
          IF (mype.eq.0) THEN
             WRITE(*,"(A)") "maximum nr of processes in each direction:"
             WRITE(*,"(A)") "   s   v   w   x   y   z"
             WRITE(*,"(6(I4))") nps2, npv2, npw2, npx2, npy2, npz2
          ENDIF
          write(*,"(A)") "no auto-parallelization possible for given parameters"
          call mpi_abort(MPI_COMM_WORLD,42,ierr)
       ENDIF

       n_procs_s = n_procs_s_best
       n_procs_v = n_procs_v_best
       n_procs_w = n_procs_w_best
       n_procs_y = n_procs_y_best
       n_procs_x = n_procs_x_best
       n_procs_z = n_procs_z_best
       perf_vec  = perf_vec_best
       nblocks   = nblocks_best
       axpy_ij_use_blas = axpy_ij_use_blas_best

       IF ((mype==0).and.(print_ini_msg)) THEN
          WRITE(*,"(A,I4,A)") "Checked ",count," parallelizations."
       ENDIF

    ELSE ! if not auto_parall
       ! no auto parallelization
       if (optim_perf) then     
          ! only performance optimization by different routines
          if ((mype.eq.0).and.(print_ini_msg))  then
             write(*,"(A)") ' '
             write(*,"(A)") '------------------------------'
             write(*,"(A)") 'Starting optimization routine '
          endif
          print_ini_msg = .false.

          call split_comm
          call initialize_discretization(print_ini_msg)
          if (nblocks.ne.0) then
             if (modulo(lklmn0,nblocks).ne.0) opt_nblocks=.true.
          endif
          if (opt_nblocks) call estimate_blocksize(nblocks)
          call initialize_axpy_ij(lij0,lklmn0/nblocks)
          call initialize_optimize_perf
          if (.not.optimize_perf(elapsed_time,fastest_first_perfvec_time,MFR_min,MFR_max,&
               &MFR_avg,total_mem_min, total_mem_max)) THEN
             IF (mype.EQ.0) THEN
                WRITE(AUTOPAR_FILE,"(A,6I4,a,F10.1,a,2X)",advance='no') &
                     &"skipped perf_opt for parallelization ",n_procs_s, &
                     & n_procs_v, n_procs_w, n_procs_x, n_procs_y, n_procs_z,&
                     ' due to memory estimate of ',total_mem_min,' MB/core'
             END IF
             STOP
          endif
          best_parall_time=elapsed_time
          call finalize_optimize_perf
          call free_comm
       else
          call initialize_discretization(print_ini_msg)
          if (nblocks.ne.0) then
             if (modulo(lklmn0,nblocks).ne.0) opt_nblocks=.true.
          endif
          if (opt_nblocks) call estimate_blocksize(nblocks)
          call initialize_axpy_ij(lij0,lklmn0/nblocks)
       endif
    ENDIF

    print_ini_msg = sav_print_ini_msg

    if ((mype==0).and.print_ini_msg) then 
       write(*,*)
       write(*,"(A,6I4)") "Choice for parallelization: ",&
               &n_procs_s,n_procs_v, n_procs_w, n_procs_x, n_procs_y, n_procs_z
       write(*,"(A,9I2)") 'Choice for alternative routines: ',perf_vec
       write(*,"(A,I5)") 'Choice for number of blocks: ', nblocks
       WRITE(*,"(A,I7)") "nblocks = ",nblocks
       WRITE(*,"(A,ES11.4)") "Time for chosen settings: ", &
         & best_parall_time
       IF (axpy_ij_use_blas) THEN
          WRITE(*,"(A)") "using BLAS routines for (i,j) loops: yes"
       ELSE
          WRITE(*,"(A)") "using BLAS routines for (i,j) loops: no"
       ENDIF
       write(*,*)
    endif
    if ((mype==0).and.perf_output_to_file) then
       write(AUTOPAR_FILE,"(A)") "-------------------------------------------------------"
       write(AUTOPAR_FILE,"(A,6I4)") "Choice for parallelization: ",&
               &n_procs_s,n_procs_v, n_procs_w, n_procs_x, n_procs_y, n_procs_z
       write(AUTOPAR_FILE,"(A,9I2)") 'Choice for alternative routines: ',perf_vec
       write(AUTOPAR_FILE,"(A,I7)") 'nblocks =  ', nblocks
       WRITE(AUTOPAR_FILE,"(A,ES11.4)") "Time for chosen settings: ", &
         & best_parall_time
       IF (axpy_ij_use_blas) THEN
          WRITE(AUTOPAR_FILE,"(A)") "using BLAS routines for (i,j) loops: yes"
       ELSE
          WRITE(AUTOPAR_FILE,"(A)") "using BLAS routines for (i,j) loops: no"
       ENDIF 
       CLOSE(AUTOPAR_FILE)
    endif
    call my_barrier

    if ((optim_perf.or.auto_parall).and.print_ini_msg) then
       call get_systime(time_opt2)
       if (mype.eq.0) then 
          write(*,"(A,F8.2,A)") &
            &'time for performance optimization: ', time_opt2-time_opt1,'s'
          write(*,*) ''
#ifdef COMBI_MGR
          time_perf = time_opt2-time_opt1
#endif
       endif
    endif
    
    !for derived schemes (full matrix, implicit..)
    if(comp_type.eq.'IV') then
       select case(timescheme)
       case('IE1f')
          derived_method=1
       case('IE1s')
          derived_method=2
       case('IE1p')
          derived_method=3
       case('RK3f','RK4f','RK4Mf')
          derived_method=1
       end select
    elseif(which_ev.eq.'shift_invert_s') then
       derived_method=2 
    endif

    if(derived_method.ne.0) then
       direct_rhs=.false.
    endif

    deallocate(perf_vec_best, num_opts, perf_vec_in)

    in_perf_opt = .false.

    !PERFOFF

  end subroutine optimize_parall_perf

  !>If the runtime is determined by the preconditioner, the autoparallelization should
  !if possible avoid parallelizing in directions where the preconditioner matrix is 
  !not block diagonal
  subroutine set_precond_constraints
    integer:: p_max_bd, p_min_nd

    !compute the maximum number processors in block diagonal directions
    !and the minimum number of processors in the other directions
    p_max_bd=max(n_procs_s,n_spec)
    p_min_nd=max(n_procs_x,1)*max(n_procs_z,1)* max(n_procs_v,1)
#ifdef prec_coll
    if (collision_op.eq.'none') then
       p_max_bd=p_max_bd*max(n_procs_w,nw0)
    else
       p_min_nd=p_min_nd*max(n_procs_w,1)
    end if
#else
    p_max_bd=p_max_bd*max(n_procs_w,nw0)
#endif
    
    if ((n_procs_sim/p_min_nd).le.p_max_bd) then
       !it is possible to avoid the non-diagonal directions altogether
       if (n_procs_x.lt.1) n_procs_x=1
       if (n_procs_z.lt.1) n_procs_z=1
       if (n_procs_v.lt.1) n_procs_v=1
#ifdef prec_coll
       if (collision_op.ne.'none') if (n_procs_w.lt.1) n_procs_w=1
#endif
    else
       !it is not possible to avoid non-diagonal directions
       !to minimize non-diagonal parallelization, we set the diagonal directions the maximum
       if (n_procs_s.lt.1) n_procs_s=GCD(n_procs_sim,n_spec)
#ifdef prec_coll
       if (collision_op.eq.'none') if (n_procs_w.lt.1) n_procs_w=GCD(n_procs_sim,nw0)
#else
       if (n_procs_w.lt.1) n_procs_w=GCD(n_procs_sim,nw0)
#endif
    end if
    
    if (mype.eq.0) write(*,'(a)') 'adapted autoparallelization ranges to optimize preconditioner'

  end subroutine set_precond_constraints

  !>This routine estimates the optimal number of blocks for the computation of the RHS of the Vlasov equation.
  !!It is called once for each parallelization; the performance optimization is then carried out with the
  !!estimated number of blocks.
  !!\todo If a fixed perf_vec is given, this routine is currently not called, even if opt_nblocks=.t..
  SUBROUTINE estimate_blocksize(nblocks)
    integer, intent(inout):: nblocks

    real:: mintime, mytime_new, mytime_l, mytime_u !, mintime_intervall
    integer :: optimal_nblocks, valid_nbs !, nb, optimal_nblocks_intervall
    integer :: ind_l, ind_u, ind_new
#define DO_BISECTION_FOR_NBLOCKS 1
#ifndef DO_BISECTION_FOR_NBLOCKS
    integer :: ind_min
#endif

    integer, dimension(:), allocatable :: valid_nblocks, tmp_nb
#ifdef WITH_CUDA_NONLIN
    integer :: nb, nb_valid,nParts,blocksize
#endif
    integer:: k,l,m,n

    !nested intervall method

    !compute all valid block numbers. to be able to map the blocks to a 
    !k,l,m,n subgrid, the outer index ranges have to be fully decomposed 
    !before inner indices are (partially) decomposed
    allocate(tmp_nb(lklmn0))
    valid_nbs=1

    do n=1,ln0
       if (modulo(ln0,n).eq.0) then
          tmp_nb(valid_nbs)=n
          valid_nbs=valid_nbs+1
       end if
    end do
    do m=2,lm0
       if (modulo(lm0,m).eq.0) then
          tmp_nb(valid_nbs)=m*ln0
          valid_nbs=valid_nbs+1
       end if
    end do
    do l=2,ll0
       if (modulo(ll0,l).eq.0) then
          tmp_nb(valid_nbs)=l*lm0*ln0
          valid_nbs=valid_nbs+1
       end if
    end do
    do k=2,lk0
       if (modulo(lk0,k).eq.0) then
          tmp_nb(valid_nbs)=k*ll0*lm0*ln0
          valid_nbs=valid_nbs+1
       end if
    end do
    
    !remove the last increment
    valid_nbs=valid_nbs-1

#ifdef WITH_CUDA_NONLIN
    !print*,"before filter: tmp_nb(1:valid_nbs) = ",tmp_nb(1:valid_nbs)
    nb_valid = 0
    do nb=1,valid_nbs
       blocksize = lklmn0/tmp_nb(nb)
       ! nParts hard-coded, must match nParts in cuda_overlap.h
       nParts = 4
       if (modulo(2*blocksize,nParts).ne.0) then
          tmp_nb(nb) = -1
       else
          nb_valid = nb_valid + 1
       end if
    end do
    !print*,"after filter: tmp_nb(1:valid_nbs) = ",tmp_nb(1:valid_nbs)
    allocate(valid_nblocks(nb_valid))
    nb_valid = 1
    do nb=1,valid_nbs
       if (tmp_nb(nb).ne.-1) then
          valid_nblocks(nb_valid) = tmp_nb(nb)
          nb_valid = nb_valid + 1
       end if
    end do
    ! set valid_nbs to the real number of valid nblocks in valid_nblocks
    valid_nbs = nb_valid-1
#else
    allocate(valid_nblocks(valid_nbs))
    valid_nblocks=tmp_nb(1:valid_nbs)
#endif
    deallocate(tmp_nb)

    if (valid_nbs.eq.0) then
       nblocks = lklmn0
       return
    endif

    optimal_nblocks=1
    mintime = 1000000
#if DO_BISECTION_FOR_NBLOCKS
    ind_u = valid_nbs
    ind_l = 1

    !call test_computation(valid_nblocks(ind_u),mytime_u,.true.)
    !call test_computation(valid_nblocks(ind_l),mytime_l,.true.)
    call test_computation(valid_nblocks(ind_u),mytime_u,.false.)
    call test_computation(valid_nblocks(ind_l),mytime_l,.false.)

    do
       ind_new = (ind_u+ind_l)/2
       if ((ind_new.eq.ind_l).or.(ind_new.eq.ind_u)) exit

       !call test_computation(valid_nblocks(ind_new),mytime_new,.true.)
       call test_computation(valid_nblocks(ind_new),mytime_new,.false.)

!if (mype.eq.0) write(AUTOPAR_FILE,'(3(I3,X,I5,F12.6))') ind_l, valid_nblocks(ind_l),mytime_l, &
!     & ind_new, valid_nblocks(ind_new), mytime_new, ind_u, valid_nblocks(ind_u), mytime_u
       
       if (mytime_l.lt.mytime_u) then
          ind_u = ind_new
          mytime_u = mytime_new          
       else
          ind_l = ind_new
          mytime_l = mytime_new
       endif
    enddo

    if (mytime_l.lt.mytime_u) then
       optimal_nblocks = valid_nblocks(ind_l)
       mintime=mytime_l
    else
       optimal_nblocks = valid_nblocks(ind_u)
       mintime=mytime_u
    endif
#else
    do ind_new=1,valid_nbs
       call test_computation(valid_nblocks(ind_new),mytime_new,.true.)
       
       if (mytime_new.lt.mintime) then
          ind_min = ind_new
          mintime = mytime_new
       endif
    enddo
    optimal_nblocks = valid_nblocks(ind_min)
#endif


!    if (OUTPUT.and.mype==0) write (AUTOPAR_FILE,"(A,I5,F12.6)") &
!         & 'Estimated optimal nblocks (nested intervalls): ', &
!         & optimal_nblocks,mintime

#if 0
    optimal_nblocks_intervall=optimal_nblocks
    mintime_intervall = mintime

    !check all possible nblocks
    optimal_nblocks=1
    mintime = 1000000
    do nb=1,lklmn0
       if (modulo(lklmn0,nb).eq.0) then
          call test_computation(nb,mytime_new,output)

          if ((mytime_new.lt.mintime)) then
             mintime=mytime_new
             optimal_nblocks=nb
          endif
       endif
    enddo

    if (OUTPUT.and.mype==0) write (AUTOPAR_FILE,"(A,2(I5,X,F12.6))") &
         & 'Estimated optimal nblocks: ',optimal_nblocks, mintime,&
         & optimal_nblocks_intervall, mintime_intervall
#endif

    nblocks=optimal_nblocks

  END SUBROUTINE estimate_blocksize
    
  !>Performs a test computation to estimate the optimal number of blocks for CalFullRHS_0. 
  !!This routine does not appear in any memory estimate, since it is executed before the 
  !!other arrays are allocated.
  SUBROUTINE test_computation(nb,time,output)
    integer, intent(in):: nb
    real, intent(out):: time
    logical, intent(in) :: output
    complex, dimension(:,:,:)  , allocatable:: g_test
    complex, dimension(:,:,:,:), allocatable:: rhs_test
    real, dimension(:,:,:), allocatable:: prefac_test
    real, dimension(-2:2):: coeff, coeff2, coeff3 
    double precision :: time1, time2
    integer::  i, j, klmn, sten, count
    !prefac_on=.t. puts an additional prefactor into the computation
    logical :: prefac_on=.false.

    !some 'random' stencils
    coeff=(/1,-8,0,8,-1/)/(12.0)
    coeff2=(/-1,4,-6,4,-1/)* 0.0625
    coeff3=(/1,-6,3,2,0/)/(6.0)

    if (Modulo(lklmn0,nb)==0) then
       allocate(rhs_test(li1:li2,lj1:lj2,1:lklmn0/nb,1:nb))
       if (prefac_on) allocate(prefac_test(1:lklmn0/nb,1:nb,-2:2))
#ifdef WITHOMP_BLOCKLOOP
       !$OMP PARALLEL default(none) &
       !$OMP shared(rhs_test,prefac_on,prefac_test,lklmn0,nb,lj1,lj2,li1,li2) &
       !$OMP shared(time, time1, time2, coeff, coeff2, coeff3, istep_energy, istep_energy3d) &
       !$OMP private(g_test,klmn,j,i)
#endif
       !currently, g_test does not have an nblocks index, because this caused the performance of the
       !test routine to degrade faster than the actual code does at large nblocks
       allocate(g_test(li1:li2,lj1:lj2,-1:lklmn0/nb+2))
#ifdef WITHOMP_BLOCKLOOP
       !$OMP BARRIER
       !$OMP MASTER
#endif
       call my_barrier()
       Call get_systime(time1)
#ifdef WITHOMP_BLOCKLOOP
       !$OMP END MASTER
       !$OMP WORKSHARE
#endif
       rhs_test=(0.,0.)
#ifdef WITHOMP_BLOCKLOOP
       !$OMP END WORKSHARE
#endif
       do klmn=-1,lklmn0/nb+2
          do j=lj1,lj2
             do i=li1,li2
                g_test(i,j,klmn)=i*j*klmn+imag*klmn*i*j
             enddo
          enddo
       enddo
       if (prefac_on) then
          !compute some prefactor
          do sten=-2,2
             do count=1,nb
                do klmn=1,lklmn0/nb
                   prefac_test(klmn,count,sten)=((klmn-sten)-(count*sten))
                enddo
             enddo
          enddo
          !operation on a block of rhs and g
          do count=1,nb
             do klmn=1,lklmn0/nb
                do sten=-2,2
                   rhs_test(:,:,klmn,count)=rhs_test(:,:,klmn,count)+prefac_test(klmn,count,sten)&
                        *count*coeff(sten)*g_test(:,:,klmn+sten)
                enddo
             enddo
             do klmn=1,lklmn0/nb
                do sten=-2,2
                   rhs_test(:,:,klmn,count)=rhs_test(:,:,klmn,count)+prefac_test(klmn,count,sten)&
                        *count*coeff2(sten)*g_test(:,:,klmn+sten)
                enddo
             enddo
             do klmn=1,lklmn0/nb
                do sten=-2,2
                   rhs_test(:,:,klmn,count)=rhs_test(:,:,klmn,count)+prefac_test(klmn,count,sten)&
                        *count*coeff3(sten)*g_test(:,:,klmn+sten)
                enddo
             enddo
          enddo
       else
          !computation without prefactor
#ifdef WITHOMP_BLOCKLOOP
          !$OMP DO
#endif
          do count=1,nb
             if ((istep_energy.gt.0).or.(istep_energy3d.gt.0)) then
                do sten=1,6
                   g_test = (0.0,0.0)
                   do klmn=1,lklmn0/nb
                      g_test(:,:,klmn) = g_test(:,:,klmn)*&
                           &rhs_test(:,:,klmn,count)
                      g_test(:,:,klmn) = g_test(:,:,klmn)*&
                           &rhs_test(:,:,klmn,count)
                   enddo
                enddo
             endif
             do klmn=1,lklmn0/nb
                do sten=-2,2
                   rhs_test(:,:,klmn,count)=rhs_test(:,:,klmn,count)+count*coeff(sten)*&
                        g_test(:,:,klmn+sten)
                enddo
             enddo
             do klmn=1,lklmn0/nb
                do sten=-2,2
                   rhs_test(:,:,klmn,count)=rhs_test(:,:,klmn,count)+count*coeff2(sten)*&
                        g_test(:,:,klmn+sten)
                enddo
             enddo
             do klmn=1,lklmn0/nb
                do sten=-2,2
                   rhs_test(:,:,klmn,count)=rhs_test(:,:,klmn,count)+count*coeff3(sten)*&
                        g_test(:,:,klmn+sten)
                enddo
             enddo
          enddo
#ifdef WITHOMP_BLOCKLOOP
          !$OMP END DO
#endif
       endif
#ifdef WITHOMP_BLOCKLOOP
       !$OMP BARRIER
       !$OMP MASTER
#endif
       Call get_systime(time2)
       time=time2-time1
#ifdef WITHOMP_BLOCKLOOP
       !$OMP END MASTER
#endif
       deallocate(g_test)
#ifdef WITHOMP_BLOCKLOOP
       !$OMP END PARALLEL
#endif
       deallocate(rhs_test)
       if (prefac_on) deallocate(prefac_test)
    else
       time=1000
    endif

    call sync_time(time)

    if (OUTPUT.and.mype==0.and.time.ne.1000) &
         & write (AUTOPAR_FILE,"(A,I5,A,F12.6)") 'n: ',nb,' t: ',time    
  END SUBROUTINE test_computation
  
  !>Checks all perf_vecs for the current parallelization and returns the fastest time and perf_vec
  !!as well as the checksum for the result. The routine aborts if NANs or deviations are found in the 
  !!checksum
  function optimize_perf(time,fastest_first_perfvec_time,MFloprate_min, MFloprate_max, &
       MFloprate_avg,mem_min,mem_max) result(found_perf)
    real,intent(out) :: time !<best time
    real,intent(inout):: fastest_first_perfvec_time
    REAL, INTENT(OUT) :: MFloprate_max, MFloprate_min, MFloprate_avg,mem_min,mem_max
    !OPTIONAL :: mem_min,mem_max
    logical :: found_perf, perf_vec_fixed, first_perfvec

    ! Local variables
    real :: min_required_memory,MFloprate
    INTEGER:: fastest_value, ierr, i_opt, i_val, start_val, OUTPUT, out_step, count
    REAL:: fastest_time, mem_perf
    type(result_type) :: this_result
    LOGICAL, ALLOCATABLE, dimension(:) :: do_optimization

    !PERFON('optim_perf')

    perf_vec = perf_vec_in
    perf_vec_fixed = ALL(perf_vec.gt.0)
    first_perfvec=.true.

    IF (perf_vec_fixed) THEN
       dimopts = 1
       num_opts = 1
    ENDIF
       
    found_perf = .FALSE.

    ALLOCATE(do_optimization(dimopts))
    do_optimization = perf_vec_fixed

    DO i_opt=1,dimopts
       IF ((perf_vec(i_opt).LE.0).OR.(perf_vec(i_opt).GT.num_opts(i_opt))) THEN
          do_optimization(i_opt) = .true.
          perf_vec(i_opt) = 1
       END IF
    END DO

    fastest_time = -1
    min_required_memory = 1e10
    mem_min = 1e10
    mem_max = 0.0
    count = 0

    DO i_opt=1,dimopts
       IF (.NOT.do_optimization(i_opt)) CYCLE
       fastest_value = 1
       start_val = 1
       IF (count.GT.0) start_val = 2
       DO i_val=start_val,num_opts(i_opt)
          IF (.NOT.perf_vec_fixed) perf_vec(i_opt) = i_val
          IF ((mype.EQ.0).and.write_optim_perf) THEN 
             WRITE(AUTOPAR_FILE,"(3X,9(I1,1X))", advance='no') perf_vec
            call flush(AUTOPAR_FILE)
          END IF
          CALL est_mem_new(mem_perf)
          IF (mem_perf.LT.min_required_memory) min_required_memory = mem_perf

          IF (mem_perf.LT.mb_per_mpi_proc) THEN
             IF (mem_perf.LT.mem_min) mem_min = mem_perf
             IF (mem_perf.GT.mem_max) mem_max = mem_perf

             !write memory estimate before test_perf in order to
             !facilitate error analysis
             IF ((mype.EQ.0).AND.write_optim_perf) THEN
                WRITE(AUTOPAR_FILE,"(A,F9.4,A)",advance='no') ": ",&
                     &mem_perf," MB"
               call flush(AUTOPAR_FILE)
             ENDIF

             if (.not.found_perf) then  ! parallelization dependent modules initialized?
                PERFON('ini2')
                call allocate_arrays
                call initialize_blockindex
                call init_physical
                PERFOFF
             endif

             CALL test_perf(time,MFloprate,this_result) 
             count = count + 1

             CALL mpi_allreduce(MFloprate, MFloprate_min, 1, MPI_REAL_TYPE, MPI_MIN,&
                  &MY_MPI_COMM_WORLD, ierr)
             CALL mpi_allreduce(MFloprate, MFloprate_max, 1, MPI_REAL_TYPE, MPI_MAX,&
                  &MY_MPI_COMM_WORLD, ierr)
             CALL mpi_allreduce(MFloprate, MFloprate_avg, 1, MPI_REAL_TYPE, MPI_SUM,&
                  &MY_MPI_COMM_WORLD, ierr)
             MFloprate_avg = MFloprate_avg/n_procs_sim
             IF ((mype.EQ.0).AND.write_optim_perf) THEN
                WRITE(AUTOPAR_FILE,"(A,ES10.3,A)") ", t = ",time, ' sec'
               call flush(AUTOPAR_FILE)
             ENDIF
             ! is there already a reference result?
             IF (reference_result%isDefined) THEN
                ! compare result with reference result
                IF (.NOT.isEqual(reference_result,this_result)) THEN
                   ! abort with error 
                   IF (mype.EQ.0) THEN
                      !Write detailed error report to stdout *AND* into autopar.dat 
                      !(for users/developers monitoring this file)
                      out_step=AUTOPAR_FILE-STDOUT
                      if(out_step.eq.0) out_step=1
                      DO OUTPUT=STDOUT,AUTOPAR_FILE,out_step
                         WRITE(OUTPUT,"(A)") "==================="
                         WRITE(OUTPUT,"(A)") "ERROR in auto-parallelization/alternative routines:"
                         WRITE(OUTPUT,"(A,I1,A,I1,A,6I3)") 'Deviations in result for perf_vec(',&
                              i_opt,')=', perf_vec(i_opt),' for parallelization (svwxyz) ',&
                              &n_procs_s, n_procs_v, n_procs_w, n_procs_x, n_procs_y, n_procs_z
                         IF (count==1) WRITE(OUTPUT,"(A)") '(This is the first test with this parallelization.)'
                         WRITE(OUTPUT,"(A,ES20.10)") "reference result (sum square) = ",&
                              &reference_result%sum_square 
                         WRITE(OUTPUT,"(A,ES20.10)") "current result (sum square)   = ",&
                              &this_result%sum_square
                         WRITE(OUTPUT,"(A,ES20.10)") "rel. error           = ",&
                              &ABS(reference_result%sum_square-this_result%&
                              &sum_square)/ABS(this_result%sum_square)
                         WRITE(OUTPUT,"(A,ES20.10)") "limit for rel. error = ",&
                              &prefac_error*EPSILON(reference_result%sum_square)
                         WRITE(OUTPUT,"(A)") "==================="
                      ENDDO
                      if (perf_output_to_file) CLOSE(AUTOPAR_FILE)
                   END IF
                   stop
                END IF
             ELSE 
                reference_result = this_result
                fastest_time = time
                fastest_value = i_val
             END IF

             if (first_perfvec.and.skip_parall_on_first_pv) then
                if (fastest_first_perfvec_time.lt.0) then
                   fastest_first_perfvec_time=time
                else
                   if (time.gt.1.1*fastest_first_perfvec_time) then
                      call finalize_calc_aux_fields
                      call finalize_timescheme
                      call finalize_physical
                      call deallocate_arrays
                      !by returning mem_min=-1, the outer routines know that we're not skipping due
                      !to lack of memory
                      mem_min=-1
                      return
                   endif
                   if (time.lt.fastest_first_perfvec_time) fastest_first_perfvec_time=time
                endif
                first_perfvec=.false.
             endif 
                      
             IF (fastest_time.LT.0) THEN
                ! for the given parallelization, we do not yet have
                ! a fastest time, so we set it to the time of the first
                ! permutation
                fastest_time = time
                fastest_value = i_val
             END IF

             ! results are equal, we now compare the times
             IF (time.LT.fastest_time) THEN
                fastest_time = time
                fastest_value = i_val
             END IF
             found_perf = .TRUE.
          ELSEIF ((mype.EQ.0).AND.write_optim_perf) THEN
             WRITE(AUTOPAR_FILE,"(a,F10.1,a)") &
                  & ': skipped due to estimated memory requirement of',mem_perf,' MB per core'
            call flush(AUTOPAR_FILE)
          ENDIF
       END DO
       ! set to fastest value for this perf_vec option
       perf_vec(i_opt) = fastest_value
    END DO
    ! finalize parallelization dependent module if at least one possible perf_vec has
    ! been found:
    if (found_perf) then
       call finalize_calc_aux_fields !init. in test_perf
       call finalize_timescheme !init. in test_perf
       Call finalize_physical
       Call deallocate_arrays
    endif
    
    ! if all perf_vecs need too much memory, we give back the minimal required memory
    ! in the mem_min field.
    if (.not.found_perf) mem_min = min_required_memory

    time = fastest_time

    DEALLOCATE(do_optimization)
    !PERFOFF

  END FUNCTION optimize_perf

  !>Initializes gene for the current parallelization and perf_vec, runs perf_tsteps timesteps and 
  !!measures the time, then cleans up the memory and returns the measured time
  SUBROUTINE test_perf(time_used,MFloprate,comp_result)
    REAL, INTENT(out):: time_used !<measured time for the computation
    REAL, INTENT(OUT) :: MFloprate
    type(result_type) :: comp_result !, ini_result

    ! Local variables
    double precision :: time1,time2
    REAL(8) :: inc_time, inc_MFlops=0.
    INTEGER:: i,j,n !, rssize
!    INTEGER,DIMENSION(:),ALLOCATABLE:: ranseed
    LOGICAL:: loop_aborted
    real::dt_max_in,evc,dtv

    PERFON('ini3')
    call initialize_calc_aux_fields

    itime=1
    loop_aborted = .false.
    
    !save dt_max from parameters (restored for IV run)
    dt_max_in = dt_max
    !use a small, stable dt
    !the estimated dt_max is used to decide for a coll_split_scheme
    CALL estimate_dt_max(.false.,dt_max,dtv,evc)
    dt_max=0.1*dt_max
    dt=dt_max
    !need to lower the time step for turbo dealiasing (reason not known yet)
    if (turbdeal) dt = 0.4*dt

    ! g_1 initialization
    ! Gaussian in x direction (necessary for Dirichlet boundary conditions)
    ! chosen to have 1e-6 at +-lx/2 and 0.1 in the middle at x=0
    ! note: with higher amplitude, we observed rare cases where global sims.
    ! had problems with the time step adaption
    if ((xy_local).or.(.not.x_local)) then
       if(yx_order) then
          do j=lj1,lj2
             g_1(:,j,:,:,:,:) = 0.1*exp(-(real(j-nxb)/real(nx0-evenx)-0.5)**2*20.0*log(10.0))
          end do
       else
          do i=li1,li2
             g_1(i,:,:,:,:,:) = 0.1*exp(-(real(i)/real(nx0-evenx)-0.5)**2*20.0*log(10.0))
          end do
       end if
    else
       g_1 = 0.1
    endif
    !call calculate_result(ini_result,g_1)
    !if (mype==0) WRITE(*,"(A,ES20.10)") "ini result (sum square) = ",&
    !     &ini_result%sum_square 

    call initialize_timescheme 

    !WRITE(*,"(I3,A,ES20.12)") mype,": perf_opt, g_1 = ",REAL(SUM(g_1*CONJG(g_1)))
    if (turbdeal) then
       !initialize random phase
       call sgrnd(1)
    endif

    PERFOFF

    ! reset the performance counter for the following region
    call my_barrier
    !PERF_RESET('perf_ste')

    PERFON('perf_ste')
    Call get_systime(time1)

    do n=1,perf_tsteps
       call calc_timestep
       if(fast_opt.and.(n.eq.1)) then
          !if the first timestep is too slow (130% of fastest 
          !perf so far), skip rest
          !!!!this significantly speeds up the initialization!!!!
          call get_systime(time2)
          time_used=time2-time1
          call sync_time(time_used)
          call calculate_result(comp_result,g_1)
          if(time_used.gt.(1.3*smallest_t1)) then
             loop_aborted=.true.
             exit 
          elseif(time_used.lt.smallest_t1) then
             smallest_t1=time_used
          endif
       endif
    enddo

    call get_systime(time2)
    PERFOFF

    PERF_GET('perf_ste',inc_time, inc_MFlops)
    IF (n.GT.1) time_used = (time2 - time1)/REAL(n-1)
    call sync_time(time_used)
    !IF (loop_aborted) time_used=perf_tsteps*time_used

    if (.not.fast_opt) CALL calculate_result(comp_result,g_1)
    !if (mype==0) WRITE(*,"(A,ES20.10)") "comp result (sum square) = ",&
    !         &comp_result%sum_square 

    MFloprate = REAL(inc_MFlops, kind(MFloprate))
    !WRITE(*,"(I3,A,F7.1)") mype,": inc_MFlops = ",inc_MFlops
    !WRITE(*,"(I3,A,ES20.10)") mype,"time2 = ",time2
    dt_max = dt_max_in
    
  end subroutine test_perf


  !>Estimates the memory requirement of GENE that is unavoidable for a given 
  !!parallelization per core
  !!\todo routines to estimate the memory consumption should be available in each module and only be called here,
  !!furthermore it seems that the memory estimate is too low for power (e.g. big test 3)
  subroutine est_mem_new(mem_req)
    real, intent(out):: mem_req !<estimate for the GENE memory requirement in MB/core
    REAL:: mem_test
    mem_test=0.
    mem_test=mem_est_arrays(mem_test)
    mem_test=mem_est_aux_fields(mem_test)
    mem_test=mem_est_init_physical(mem_test)
    mem_test=mem_est_timescheme(mem_test)
    mem_test=mem_est_diag(mem_test)
    if (istep_omega.gt.0) mem_test=mem_est_convergence_monitoring(mem_test)
    if(trap_pass) mem_test=mem_est_hybrid(mem_test)

    mem_req=mem_test
  end subroutine est_mem_new

  SUBROUTINE assign_result_type_to_real(the_real,res)
    TYPE(result_type),intent(IN) :: res
    REAL,intent(OUT) :: the_real

    the_real = res%sum_square
  END SUBROUTINE assign_result_type_to_real

  SUBROUTINE assign_result_type(res_1, res_2)
    TYPE(result_type),intent(OUT) :: res_1
    TYPE(result_type),intent(IN) :: res_2

    res_1%sum_square = res_2%sum_square
    res_1%sum_re_plus_im = res_2%sum_re_plus_im
    res_1%isDefined = res_2%isDefined
  END SUBROUTINE assign_result_type

  FUNCTION isNan_real(number) 
    real :: number
    logical :: isNan_real

    isNan_real = (number.ne.number)
  END FUNCTION isNan_real

  FUNCTION isNan_result_type(res) 
    type(result_type) :: res
    logical :: isNan_result_type

    isNan_result_type = (res%sum_square.ne.res%sum_square)
  END FUNCTION isNan_result_type

  FUNCTION isEqual_result_type(res1, res2)
    type(result_type) :: res1, res2
    logical :: isEqual_result_type
    
    !print*,"res1 = ", res1%sum_square,", res2 = ",res2%sum_square
    isEqual_result_type = (ABS(res1%sum_square-res2%sum_square).LE.prefac_error*EPSILON(res1%sum_square)*ABS(res2%sum_square))
  END FUNCTION isEqual_result_type


  SUBROUTINE mp_calculate_result(res, u)
    type(result_type) :: res
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2) :: u

    ! Local variables
    REAL :: sum_square, lres, local_sum_re_im, sum_re_im
    real :: compensate, tempval, tempsum
    real,dimension(:),allocatable :: lres_recv
    integer :: ierr, i,j,k,l,m,n

    !discard unpaired high kx mode
    if(xy_local.and.(evenx.eq.1)) then
       if (yx_order) then
          u(:,hkx+1,:,:,:,:)=(0.,0.)
       else
          u(hkx+1,:,:,:,:,:)=(0.,0.)
       end if
    end if
    ! better sum algorithm (from Wikipedia)
    ! Kahan summation algorithm
    !function KahanSum(input)
    !var sum = 0.0
    !var c = 0.0          //A running compensation for lost low-order bits.
    !for i = 1 to input.length do
    !    y = input[i] - c    //So far, so good: c is zero.
    !    t = sum + y         //Alas, sum is big, y small, so low-order digits of y are lost.
    !    c = (t - sum) - y   //(t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
    !    sum = t             //Algebraically, c should always be zero. Beware eagerly optimising compilers!
    !    //Next time around, the lost low part will be added to y in a fresh attempt.
    !return sum
    lres = 0.0D0
    compensate = 0.0D0
    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             do k=lk1,lk2
                do j=lj1,lj2
                   do i=li1,li2
                      tempval = real(u(i,j,k,l,m,n)*conjg(u(i,j,k,l,m,n))) - compensate
                      tempsum = lres + tempval
                      compensate = (tempsum - lres) - tempval
                      lres = tempsum
                   end do
                end do
             end do
          end do
       end do
    end do
    allocate(lres_recv(n_procs))
    !initialize lres_recv to avoid floating invalid for mype.ne.0
    lres_recv = 0.0D0
    !lres = REAL(SUM(u*CONJG(u)),KIND(lres))
    local_sum_re_im = REAL(SUM(REAL(u)+AIMAG(u)),KIND(lres))
    !WRITE(*,"(I3,2(A,ES21.15))") mype,": lres = ", lres,", local_sum_re_im = ",local_sum_re_im
    call mpi_gather(lres,1,MPI_REAL_TYPE,lres_recv,1,MPI_REAL_TYPE,0,MY_MPI_COMM_WORLD,ierr)
    sum_square = 0.0D0
    compensate = 0.0D0
    do i=1,n_procs
       tempval = lres_recv(i) - compensate
       tempsum = sum_square + tempval
       compensate = (tempsum - sum_square) - tempval
       sum_square = tempsum
    end do
    !sum_square = sum(lres_recv)
    call mpi_bcast(sum_square,1,MPI_REAL_TYPE,0,MY_MPI_COMM_WORLD,ierr)
    !Call mpi_allreduce (lres, sum_square, 1, MPI_REAL_TYPE, MPI_SUM,&
    !     MY_MPI_COMM_WORLD, ierr)
    Call mpi_allreduce (local_sum_re_im, sum_re_im, 1, MPI_REAL_TYPE, MPI_SUM,&
         MY_MPI_COMM_WORLD, ierr)
    deallocate(lres_recv)
    res%sum_square = sum_square
    res%sum_re_plus_im  = sum_re_im
    res%isDefined = .TRUE.
    !IF (mype.EQ.0) WRITE(*,"(I3,2(A,ES21.15))") mype,": sum_square = ",res%sum_square,", sum_re_plus_im = ",res%sum_re_plus_im

  END SUBROUTINE mp_calculate_result

  FUNCTION get_used_memory(arena) RESULT(umem)
    INTEGER, intent(OUT) :: arena
    INTEGER :: umem
    
    INTEGER :: back, ordblks, uordblks, fordblks, hblkhd
    INTEGER :: fortran_mallinfo, get_datasegment_size

    back = fortran_mallinfo(arena, ordblks, uordblks, fordblks, hblkhd)
#ifndef POWER
    ! on Linux systems, arena from the mallinfo call does not give meaningful results
    arena = get_datasegment_size()
#endif

    umem = uordblks+hblkhd
  END FUNCTION get_used_memory


  FUNCTION get_next_parallelization(np_x,np_y,np_z,np_v,np_w,np_s) RESULT(found)
    logical :: found
    INTEGER,intent(OUT) :: np_x, np_y,np_z,np_v,np_w,np_s

    ! Local variables
    logical :: do_increment
    INTEGER, SAVE :: loc_np_x, loc_np_y, loc_np_z, loc_np_v, loc_np_w, loc_np_s
  

    found = .true.
    IF (first_parall) THEN
       ! initialize
       loc_np_x = npx1
       loc_np_y = npy1
       loc_np_z = npz1
       loc_np_v = npv1
       loc_np_w = npw1
       loc_np_s = nps1
       do_increment = .true.
    ELSE
       loc_np_s = loc_np_s + 1
       do_increment = .false.
    END IF
    ! the loop order is (from outer to inner):
    ! y,x,v,z,w,s
    
    DO WHILE (loc_np_x*loc_np_y*loc_np_z*loc_np_v*loc_np_w*loc_np_s.NE.n_procs_sim) 
       IF (do_increment) THEN
          loc_np_s = loc_np_s + 1
       END IF
       do_increment = .true.
       IF ((loc_np_s.GT.nps2).OR.(loc_np_y*loc_np_x*loc_np_v*loc_np_z*loc_np_w*loc_np_s.GT.n_procs_sim)) THEN
          loc_np_s = nps1
          
          loc_np_w = loc_np_w + 1
          IF ((loc_np_w.GT.npw2).OR.(loc_np_y*loc_np_x*loc_np_v*loc_np_z*loc_np_w.GT.n_procs_sim)) THEN
             loc_np_w = npw1

             loc_np_z = loc_np_z + 1
             IF ((loc_np_z.GT.npz2).OR.(loc_np_y*loc_np_x*loc_np_v*loc_np_z.gt.n_procs_sim)) THEN
                loc_np_z = npz1

                loc_np_v = loc_np_v + 1
                IF ((loc_np_v.GT.npv2).OR.(loc_np_y*loc_np_x*loc_np_v.gt.n_procs_sim)) THEN
                   loc_np_v = npv1

                   loc_np_x = loc_np_x + 1
                   IF ((loc_np_x.GT.npx2).OR.(loc_np_y*loc_np_x.GT.n_procs_sim)) THEN
                      loc_np_x = npx1

                      loc_np_y = loc_np_y + 1
                      if (loc_np_y.gt.npy2) then
                         found = .false.
                         exit
                      end if
                   end if
                end if
             end if
          end if
       END IF
    END DO

    ! check if this parallelization gives the right total number of ranks

    if (found) then
       np_x = loc_np_x
       np_y = loc_np_y
       np_v = loc_np_v
       np_w = loc_np_w
       np_s = loc_np_s
       np_z = loc_np_z
    end if

    if (first_parall) first_parall=.false.
  END FUNCTION get_next_parallelization

  


END MODULE perf_opt
