#include "redef.h"
#include "intrinsic_sizes.h"
Module time_scheme

  Use par_mod
  Use communications
  Use rhs_computations
  Use aux_fields
  Use sources_mod, only: initialize_krookBuffer_operator, &
       & add_krookBuffer, finalize_krookBuffer_operator, explicit_buffer, &
       & add_krookBuffer_adptv
  use external_contr
  use RK_coefficients
  use calc_rhs, only:  mem_est_calc_rhs, rhs_nl, rhs_f0, f0_heat_src_active
  use mtrandom
  use compute_dt
#ifdef WITH_CUDA_NONLIN
  use page_locked_memory_mod
#endif
  use discretization_adptv_module
  use antenna, only: antenna_type, evolve_antenna_amplitudes

  Implicit None
  public:: initialize_timescheme, calc_timestep, finalize_timescheme, &
       & mem_est_timescheme

  ! make some subroutines and variables public for unit tests
#ifdef GIUTESTS
  public:: RK_standard, RK_standard_adptv, k1, k2, tg
#endif
  
  private
  
  LOGICAL,parameter :: OUTPUT=.false.

  Integer :: init_status = 0, init_status_RK_standard = 0, init_status_RK_low_mem = 0,&
       &init_status_timestep_coll = 0

  !auxiliary g-like arrays
  Complex, Dimension(:,:,:,:,:,:), Allocatable:: k1, k2, tg

Contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_timescheme(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    select case (timescheme)
    case('RK3lm','RK3ssp')
       mem_loc=mem_est_RK_low_mem(0.)
    case default
       mem_loc=mem_est_RK_standard(0.)    
    end select

    mem_loc=mem_est_timestep_coll(mem_loc)
    
    mem_est_timescheme=mem_req_in+mem_loc
  End Function mem_est_timescheme

  subroutine initialize_timescheme

    if (init_status.eq.0) then
       if (nonlinear) then
          If ((mype.Eq.0).and.print_ini_msg) Write(*,"(A)") "nonlinear computation"
          call initialize_adapt_dt
       else
          If ((mype.Eq.0).and.print_ini_msg) Write(*,"(A)") "linear computation"
       endif
       init_status = 1
    endif

    call initialize_shear_flow
    

    if (.not.x_local.and..not.explicit_buffer) call initialize_krookBuffer_operator

    call initialize_RK_coefficients

    !set nonlinearity
    rhs_nl=nonlinear
    !set rhs options
    rhs_f0=include_f0_contr.or.f0_heat_src_active

    call initialize_calc_k
    
    if(.not.low_mem) then
       call initialize_RK_standard
    else
       call initialize_RK_low_mem
    end if

    if (coll_split) call initialize_timestep_coll
   
    !(re-)initialize mtrandom (in case perf_opt or init_cond have changed something) 
    if (turbdeal)  call sgrnd(1)

  end subroutine initialize_timescheme

  subroutine calc_timestep

    if(.not.low_mem) then
       if(is_grid_adptv) then
          if(is_numofvpoints_const) then
             print *, "block structured grid 2.0 version is under construction"
             call exit(1)
          else
             call RK_standard_adptv
          end if
       else
          call RK_standard
       end if
    else
       call RK_low_mem
    end if

    if (coll_split) call timestep_coll

    call evolve_time_dependent_quantities
    
  end subroutine calc_timestep

  subroutine finalize_timescheme
    if (.not.x_local.and..not.explicit_buffer) call finalize_krookBuffer_operator

    call finalize_shear_flow

    if(.not.low_mem) then
       call finalize_RK_standard
    else
       call finalize_RK_low_mem
    end if 

    if (coll_split) call finalize_timestep_coll

    if(nonlinear) call finalize_adapt_dt

    call finalize_RK_coefficients
    init_status = 0

  end subroutine finalize_timescheme

  !two different classes of RK schemes---------------------------

  !>Give an estimate of the memory requirements of this module
  real function mem_est_RK_standard(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !k1, k2, tg
    mem_loc=3*SIZE_OF_COMPLEX_MB*lijklmn0
    mem_loc=mem_est_calc_rhs(mem_loc)

    mem_est_RK_standard=mem_req_in+mem_loc
  end function mem_est_RK_standard


  subroutine initialize_RK_standard

    if (init_status_RK_standard.eq.0) then
       Allocate(&
            k1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),&
            k2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),&
            tg(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))

       if(is_grid_adptv) then
          k1 = 0.
          k2 = 0.
          tg = 0.
       end if
#if WITH_CUDA_NONLIN
       ! in the CUDA case, we need to register k1 and k2 as page-locked
       ! memory for faster HOST->DEVICE transfer
       !print*,mype,": Locking k1 and k2"
       call lock_memory(k1)
       call lock_memory(k2)
#endif
       init_status_RK_standard=1
    endif

  end subroutine initialize_RK_standard
  
  subroutine initialize_timestep_coll
    
    if (init_status_timestep_coll.eq.0) then
       if (rkstages_coll.gt.1.and.(low_mem)) then
          Allocate(&
               k2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),&
               tg(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
       endif
    endif

    init_status_timestep_coll=1

  end subroutine initialize_timestep_coll
  
  subroutine finalize_timestep_coll
    
    if (init_status_timestep_coll.gt.0) then
       if (rkstages_coll.gt.1.and.(low_mem)) then
          Deallocate(k2,tg)
       endif
    endif
    init_status_timestep_coll=0

  end subroutine finalize_timestep_coll

  !
  !========== standard Runge-Kutta scheme ====
  !
  Subroutine RK_standard
    integer:: stage
    REAL :: local_sum, global_sum

    DEBUG(1,"==== START timescheme ====")

    IF (OUTPUT) THEN
       CALL calculate_test_sum(g_1,local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(//,A,ES20.12)") "start g_1 = ",global_sum
    END IF

    PERFON('RK_standard')

    !-------------------------------------
    !FIRST STAGE of RK-scheme
    stage=1
    if(turbdeal) then
       !compute random phase factor for TURBO dealiasing scheme
!       call random_number(phase_1)
!       call random_number(phase_2)
!       phase_1=phase_1-0.5
!       phase_2=phase_2-0.5
       phase_1=grnd()
       phase_2=grnd()
    endif

    if ((ExB).and.(time.ge.ExB_stime)) then
       call apply_shear_flow(g_1,dt)
    endif

    ! stage determines whether max. ExB velocity is calculated
    Call calc_k(g_1, k1, stage)

    IF (OUTPUT) THEN
       CALL calculate_test_sum(k1,local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "stage1 k1 = ",global_sum
    END IF

    ! Adapt the timestep
    if (nonlinear) call adapt_dt

    !-----------------------------------
    !ALL OTHER STAGES
    !result is accumulated in k1
    do stage=2,rkstages
       if (stage==3.and.turbdeal.and..not.xy_local) then
          phase_1=grnd()
       endif

       if (stage.eq.2) then
          call nextstage_g(a_rk(stage)*dt,k1,g_1,tg)
       else
          call nextstage_g(a_rk(stage)*dt,k2,g_1,tg)
       endif

       IF (OUTPUT) THEN
          CALL calculate_test_sum(tg,local_sum, global_sum)
          IF (mype.EQ.0) WRITE(*,"(A,I1,A,ES20.12)") "stage",stage," tg = ",global_sum
       END IF


#ifdef GIUTESTS
       ! alignment with adaptive version (only for testing!)
       if(is_grid_adptv) then
          call nullify_g1_type_outside_adptv(tg)
          call nullify_g1_type_outside_adptv(k2)
          call nullify_g1_type_outside_adptv(k1)
       end if
#endif
             
       call calc_k(tg, k2, stage)

       call add_ks(b_rk(stage)/b_rk(1),k2,k1)
    enddo

    !------------------------------
    !RESULT OF TIMESTEPPING, new g
    if(explicit_RK) then
       call add_ks(b_rk(1)*dt,k1,g_1)
    else
       g_1=k1
    endif

    if (.not.x_local.and..not.explicit_buffer) call add_krookBuffer(g_1)

    PERFOFF

    IF (OUTPUT) THEN
       CALL calculate_test_sum(g_1,local_sum,global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "END   g_1 = ",global_sum
    END IF

    DEBUG(1,"==== END  time_scheme ====")
  End Subroutine RK_standard

  !includes RK_standard and derived methods
  Subroutine timestep_coll
    integer:: stage
    REAL :: local_sum, global_sum
    
    PERFON('coll_scheme')

    !-------------------------------------
    !FIRST STAGE of RK-scheme
    stage=1
    Call calc_k_coll(g_1, k1, stage)
    
    IF (OUTPUT) THEN
       CALL calculate_test_sum(k1,local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "coll stage1 k1 = ",global_sum
    END IF

    !-----------------------------------
    !ALL OTHER STAGES
    !result is accumulated in k1
    do stage=2,rkstages_coll

       if (stage.eq.2) then
          call nextstage_g(a_rk_coll(stage)*dt,k1,g_1,tg)
       else
          call nextstage_g(a_rk_coll(stage)*dt,k2,g_1,tg)
       endif

       IF (OUTPUT) THEN
          CALL calculate_test_sum(tg,local_sum,global_sum)
          IF (mype.EQ.0) WRITE(*,"(//,A,ES20.12)") "coll tg = ",global_sum
       END IF

       call calc_k_coll(tg, k2, stage)

       call add_ks(b_rk_coll(stage)/b_rk_coll(1),k2,k1)
    enddo

    !------------------------------
    !RESULT OF TIMESTEPPING, new g
    call add_ks(b_rk_coll(1)*dt,k1,g_1)
    !implicit coll: dt is included in matrix, not implemented
    !g_1=g_1+k1

    !low_mem schemes (except turbo) seem to need k1(a_1) = g_1 for their first RK stage
    if(low_mem.and.(timescheme.ne.'turbo')) k1=g_1

    PERFOFF
    
    IF (OUTPUT) THEN
       CALL calculate_test_sum(g_1,local_sum,global_sum)
       IF (mype.EQ.0) WRITE(*,"(//,A,ES20.12)") "+COLL g_1 = ",global_sum
    END IF

  End Subroutine timestep_coll

  subroutine finalize_RK_standard

    call finalize_calc_k
#ifdef WITH_CUDA_NONLIN
    !print*,mype,": Unlocking k1 and k2"
    call unlock_memory(k1)
    call unlock_memory(k2)
#endif
    Deallocate(k1,k2,tg)
 
    init_status_RK_standard = 0
  end subroutine finalize_RK_standard

  !
  !========== low memory Runge-Kutta scheme ====
  real function mem_est_RK_low_mem(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0
    
    !k1
    mem_loc=SIZE_OF_COMPLEX_MB*lijklmn0

    mem_loc=mem_est_calc_rhs(mem_loc)

    mem_est_RK_low_mem=mem_req_in+mem_loc

  end function mem_est_RK_low_mem

  !
  !========== coll_split_scheme ====
  real function mem_est_timestep_coll(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0
    
    !k2, tg  if needed
    if (rkstages_coll.gt.1.and.low_mem) then
       mem_loc=2*SIZE_OF_COMPLEX_MB*lijklmn0
    endif

    mem_est_timestep_coll=mem_req_in+mem_loc

  end function mem_est_timestep_coll


  subroutine initialize_RK_low_mem

    if (init_status_RK_low_mem.eq.0) then
       !a_1 has been replaced with k1
       allocate(k1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    endif
    
    !k1(=a_1) has to be initialized with g
    k1=g_1
    
    init_status_RK_low_mem = 1

  end subroutine initialize_RK_low_mem
  
  !> computes a RK low_mem timestep. the auxiliary array is k1
  subroutine RK_low_mem
    integer:: stage

    PERFON('RK_lowmem')    
    if(turbdeal) then
       !compute random phase factor for TURBO dealiasing scheme
!       call random_number(phase_1)
!       call random_number(phase_2)
!       phase_1=phase_1-0.5
!       phase_2=phase_2-0.5
       phase_1=grnd()
       phase_2=grnd()
    endif

    !RK scheme
    do stage=1,rkstages
       if (stage==3.and.turbdeal.and..not.xy_local) then
          phase_1=grnd()
       endif
       !a_1 has been replaced with k1
       call calc_k(g_1, k1, stage)
    enddo

    ! adapt the timestep with the maximum ExB velocity of the first stage of the last timestep
    if (nonlinear) call adapt_dt

    if (.not.x_local.and..not.explicit_buffer) call add_krookBuffer(g_1)

    PERFOFF

  end subroutine RK_low_mem

  subroutine finalize_RK_low_mem

    if (init_status_RK_low_mem.gt.0) deallocate(k1)

    call finalize_calc_k

    init_status_RK_low_mem = 0

  end subroutine finalize_RK_low_mem
!-----------------------------------

  subroutine nextstage_g(timestep,k_plus,g_old,g_new)
    real,intent(in)::timestep
    complex,dimension(lijklmn0),intent(in):: k_plus,g_old
    complex,dimension(lijklmn0),intent(inout):: g_new

   PERFON_I('rkupd')

   g_new = g_old + timestep*k_plus

   PERFOFF_I

  end subroutine nextstage_g

  subroutine add_ks(timestep,k_plus,k_new)
    real,intent(in):: timestep
    complex,dimension(lijklmn0),intent(in):: k_plus
    complex,dimension(lijklmn0),intent(inout):: k_new    

       PERFON_I('add_ks')
       !!! don't change to lijklmn0 as saxpy does not accept
       !!! kind=8 integers on all compilers !!!
       Call saxpy(2*li0*lj0*lk0*llm0*ln0,timestep,k_plus,1,k_new,1)
       PERFOFF_I

  end subroutine add_ks

  !> This routine is called after a timestep has been completed, and updates
  ! external (imposed) quantities like antenna amplitudes
  subroutine evolve_time_dependent_quantities

    if (any(antenna_type.eq.(/2,3/))) call evolve_antenna_amplitudes
    
  end subroutine evolve_time_dependent_quantities


  !========== standard Runge-Kutta scheme for adaptive grids ====
  !
  Subroutine RK_standard_adptv
    integer:: stage
    REAL :: local_sum, global_sum

    DEBUG(1,"==== START timescheme ====")

    IF (OUTPUT) THEN
       CALL calculate_test_sum(g_1,local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(//,A,ES20.12)") "start g_1 = ",global_sum
    END IF

    PERFON('RK_standard')

    !-------------------------------------
    !FIRST STAGE of RK-scheme
    stage=1
    if(turbdeal) then
       !compute random phase factor for TURBO dealiasing scheme
       phase_1=grnd()
       phase_2=grnd()
    endif
    if ((ExB).and.(time.ge.ExB_stime)) then
       call apply_shear_flow_adptv(dt)
    endif

    ! stage determines whether max. ExB velocity is calculated
    call calc_k_adptv(g_1, k1, stage)

    IF (OUTPUT) THEN
       CALL calculate_test_sum(k1,local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "stage1 k1 = ",global_sum
    END IF

    ! Adapt the timestep
    if (nonlinear) then
       stop "nonlinearity is not implemented for the adaptive grids yet"
       call adapt_dt
    endif

    !-----------------------------------
    !ALL OTHER STAGES
    !result is accumulated in k1
    do stage=2,rkstages
       if (stage==3.and.turbdeal.and..not.xy_local) then
          phase_1=grnd()
       endif

       if (stage.eq.2) then
          call nextstage_g_adptv(a_rk(stage)*dt,k1,g_1,tg)
          !call nextstage_g(a_rk(stage)*dt,k1,g_1,tg)
       else
          call nextstage_g_adptv(a_rk(stage)*dt,k2,g_1,tg)
          !call nextstage_g(a_rk(stage)*dt,k2,g_1,tg)
       endif

       IF (OUTPUT) THEN
          CALL calculate_test_sum(tg,local_sum, global_sum)
          IF (mype.EQ.0) WRITE(*,"(A,I1,A,ES20.12)") "stage",stage," tg = ",global_sum
       END IF

       call calc_k_adptv(tg, k2, stage)

       ! TODO: should be modified for adaptive grids
       call add_ks(b_rk(stage)/b_rk(1),k2,k1)
    enddo

    !------------------------------
    !RESULT OF TIMESTEPPING, new g
    if(explicit_RK) then
       call add_ks(b_rk(1)*dt,k1,g_1)
    else
       g_1=k1
    endif

    if (.not.x_local.and..not.explicit_buffer) call add_krookBuffer_adptv(g_1)

    PERFOFF

    IF (OUTPUT) THEN
       CALL calculate_test_sum(g_1,local_sum,global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "END   g_1 = ",global_sum
    END IF

    DEBUG(1,"==== END  time_scheme ====")
  End Subroutine RK_standard_adptv

  subroutine nextstage_g_adptv(timestep,k_plus,g_old,g_new)
    real,intent(in)::timestep
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: k_plus,g_old
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(inout):: g_new
    integer :: l, m, i1, i2

   PERFON_I('rkupd')

   do m = lm1, lm2
      do l = ll1, ll2
         i1 = li1_vwadp(l,m)
         i2 = li2_vwadp(l,m)
         g_new(i1:i2,:,:,l,:,:) = g_old(i1:i2,:,:,l,:,:) + &
              & timestep*k_plus(i1:i2,:,:,l,:,:)
      end do
   end do

   PERFOFF_I

 end subroutine nextstage_g_adptv

#ifdef GIUTESTS
 subroutine nullify_g1_type_outside_adptv(g1_type)
   ! use discretization, only: li1, li2, lj1, lj2, &
   !       lk1, lk2, ll1, ll2, lm1, lm2, ln1, ln2
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), &
         intent(inout) :: g1_type
    integer :: i, j, k, l, m, n

    do n = ln1, ln2
       do m = lm1, lm2
          do l = ll1, ll2
             do k = lk1, lk2
                do j = lj1, lj2
                   do i = li1, li1_vwadp(l,m)-1
                      g1_type(i,j,k,l,m,n) = 0.
                   enddo
                   do i = li2_vwadp(l,m)+1, li2
                      g1_type(i,j,k,l,m,n) = 0.
                   enddo
                enddo
             enddo
          endDo
       endDo
    enddo

  end subroutine nullify_g1_type_outside_adptv
#endif

End Module time_scheme
