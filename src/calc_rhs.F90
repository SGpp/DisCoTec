#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

!>Computes the full right hand side of the gyrokinetic Vlasov equation
!!from a given modified distribution function g_1
!!
!!This module calls the primary implementations of the RHS of the gyrokinetic 
!!equation; it is also used by the derived methods, e.g. to compute an explicit 
!!representation of the linear operator or for the matrix-free methods used in 
!!PETSc/SLEPc
Module calc_rhs 
  Use aux_fields
  Use antenna, only: antenna_type, add_dApar_dt_antenna
  Use par_mod
  Use prefactors
  use boundaries
  !USE nonlinearity
  use parallel_nonlin
  use f0_term
  use dchidxy_terms
  use dgdxy_terms
  use dfdxy_terms
  use dchidz_term
  use dfdzv_terms
  use dzv_terms
  Use collisions
  Use aux_fields
  use gyro_average_df_mod, only: gyro_average_df_wb
  use blockindex
  use sources_mod
  use RK_coefficients
  use axpy
  use numerical_damping
  use communications, only: reduce_per_thread, threadlocal_mpi_comm_y, my_secnds
  use nonlinear_term_mod
  use df_nonlinear_term_mod
  use df_arakawa_nonlinear_term_mod
  use ff_nonlinear_term_mod
  use ff_nonlinear_term_h_mod
  use ff_yx_nonlinear_term_mod
#ifdef WITH_C_NONLIN
  use df_nonlinear_c_mod
#endif
#ifdef WITH_CUDA_NONLIN
  use df_nonlinear_cuda_mod
#endif
#ifdef WITH_SPLIT_NONLIN
  use df_split_nonlinear_term_mod
#endif
#ifdef WITHOMP_NONLIN
  use df_omp_nonlinear_term_mod
#endif
#ifdef WITHOMPXY_NONLIN
  use df_ompxy_nonlinear_term_mod
#endif
#ifdef WITH_MIC_NONLIN
  use df_mic_nonlinear_term_mod
#endif
  use discretization_adptv_module
  ! x_derivative module is necessary to determine correct ranges
  ! for chi_block computations
  use x_derivatives, only: radscheme
  USE all_rhs_terms, only: this_nonlinear_term
  implicit None

  Public:: initialize_CalFullRhs, CalFullRhs, calc_rhs_only, &
       & finalize_CalFullRhs, mem_est_calc_rhs

  public:: rhs_nl, rhs_f0, f0_heat_src_active, rhs_nlev

  !for operator splitting of collisions and vlasov part in the time scheme
  public:: rhs_only_coll, rhs_with_coll

  !the following pointers are also used by get_energy_terms
  public :: ptr_bar_emfields, ptr_barchi, ptr_dbarchidxy, ptr_dgdxy, chi_block
  public :: this_nonlinear_term

  ! routines for adaptive grids
  public :: CalFullRhs_adptv, calc_rhs_only_adptv

#ifdef GIUTESTS
  public :: with_dfdperp, with_g_update
  public :: set_bar_emfields, bar_emfields
  public :: g_block_storage
#endif
  
  private

  integer :: init_status = 0, init_status_nlev = 0
  !real(4),save :: nl_sum_time=0.0

  logical:: rhs_nl=.false., rhs_f0=.false., rhs_nlev=.false.
  logical:: rhs_only_coll=.false.
  logical:: rhs_with_coll=.true.

  logical:: with_dfdperp, with_g_update

  complex, dimension(:,:,:,:,:,:),allocatable, target :: bar_emfields
  complex, dimension(:,:,:),allocatable, target :: p_barchi
  complex, dimension(:,:,:,:),allocatable, target :: p_dbarchidxy, p_dgdxy
  complex, dimension(:,:,:,:),allocatable,target :: rhs_big
  complex, dimension(:,:,:),allocatable,target:: rhs_small
  complex, dimension(:,:,:), allocatable, target :: g_block_storage
  complex, dimension(:,:,:), allocatable :: chi_block
  complex, dimension(:,:,:,:),allocatable :: emfields_nlev

  !$DIR ATTRIBUTES ALIGN:64 :: rhs_small, chi_block

  !pointers are employed in this module in order to allow for 'optional' 
  !arguments in subroutine calls
  complex, dimension(:,:,:,:,:,:),pointer :: ptr_bar_emfields
  complex, dimension(:,:,:),pointer :: ptr_barchi
  complex, dimension(:,:,:,:),pointer :: ptr_dbarchidxy, ptr_dgdxy

#ifdef WITHOMP_BLOCKLOOP
  !$OMP THREADPRIVATE(ptr_barchi,p_barchi)
  !$OMP THREADPRIVATE(p_dbarchidxy,ptr_dbarchidxy)
  !$OMP THREADPRIVATE(ptr_dgdxy,p_dgdxy, g_block_storage,chi_block)
  !$OMP THREADPRIVATE(rhs_small)
  !$OMP THREADPRIVATE(this_nonlinear_term)
#endif
Contains

  function mem_est_calc_rhs(mem_req_in)
    real:: mem_req_in, mem_est_calc_rhs
    real:: mem_loc
    logical :: to_deallocate
    
    mem_loc=0.
    mem_loc=mem_est_klmn_conv_arrs(mem_req_in)
    mem_loc=mem_est_dfdzv(mem_loc)
    if((hyp_x.gt.0).or.(hyp_y.gt.0)) &
         & mem_loc=mem_est_dfdxy(mem_loc)
    mem_loc=mem_est_dchidxy(mem_loc)
    mem_loc=mem_est_dchidz(mem_loc)
    mem_loc=mem_est_dgdxy(mem_loc)
    !    mem_loc=mem_est_coll(mem_loc)
    !if (nonlinear) mem_loc = mem_est_nonlinearity(mem_loc)
    if (nonlinear) then
#ifdef WITHOMP_BLOCKLOOP
       !$OMP PARALLEL
       !$OMP CRITICAL
#endif
       if (.not.associated(this_nonlinear_term)) then
          if (xy_local) then
             if (yx_order) then
                allocate(ff_yx_nonlinear_term_t::this_nonlinear_term)
             else
                if (.not.nonlin_h) then
                   allocate(ff_nonlinear_term_t::this_nonlinear_term)
                else
                   allocate(ff_nonlinear_term_h_t::this_nonlinear_term)
                endif
             end if
             to_deallocate = .true.
          else
             !IF (.NOT.x_local) THEN
                if (arakawa) then
                   allocate(df_arakawa_nonlinear_term_t::this_nonlinear_term)
                else
#ifdef WITH_C_NONLIN
                   allocate(df_nonlinear_c_t::this_nonlinear_term)
#elif defined(WITH_SPLIT_NONLIN)
                   allocate(df_split_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITH_CUDA_NONLIN)
                   allocate(df_nonlinear_cuda_t::this_nonlinear_term)
#elif defined(WITH_MIC_NONLIN)
                   allocate(df_mic_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITHOMP_NONLIN)
                   allocate(df_omp_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITHOMPXY_NONLIN)
                   allocate(df_ompxy_nonlinear_term_t::this_nonlinear_term)
#else
                   allocate(df_nonlinear_term_t::this_nonlinear_term)
#endif
                end if
                to_deallocate = .true.
             !end if
          end if
          call this_nonlinear_term%construct(equil_par_curr)
       else
          ! already associated
          to_deallocate = .false.
       end if

       if (associated(this_nonlinear_term)) then
          mem_loc = this_nonlinear_term%mem_est(mem_loc)
          if (to_deallocate) then
             call this_nonlinear_term%destruct()
             deallocate(this_nonlinear_term)
          end if
       else
          WRITE(*,"(A)") "this_nonlinear_term is not associated in calc_rhs. This should not happen!"
          STOP
       end if
#ifdef WITHOMP_BLOCKLOOP
       !$OMP END CRITICAL
       !$OMP END PARALLEL
#endif
    end if

    !p_barchi
    if (nonlinear) mem_loc=mem_loc+&
         &lij0*lklmn0/nblocks*SIZE_OF_COMPLEX/(1024.)**2

    if (parallel_nl) mem_loc = mem_est_parallel_nonlin(mem_loc)

    if (.not.xy_local) then
       mem_loc = mem_est_f1_sources(mem_loc)
       !bar_emfields
       mem_loc = mem_loc + 1.*lij0*lz0*lm0*ln0*n_fields*SIZE_OF_COMPLEX/(1024.)**2
       !p_dbarchidxy, p_dgdxy
       if (nonlinear) mem_loc=mem_loc+4.*lij0*lklmn0/nblocks*SIZE_OF_COMPLEX/(1024.)**2
    endif

    mem_loc = mem_est_f0_term(mem_loc)

    mem_est_calc_rhs=mem_req_in+mem_loc

  end function mem_est_calc_rhs


  !>Wrapper for the initialization routines of the different implementations
  subroutine initialize_CalFullRhs(for_nlev,for_petsc)
    implicit none
    logical, optional:: for_nlev
    logical, optional:: for_petsc
    logical:: dfdzv_replace_rhs
#ifdef WITHOMP
    integer :: omp_get_thread_num
#endif
    integer :: my_thread

    character(len=MAX_TYPENAME_LENGTH) :: type_of_nonlin

    !nl_sum_time = 0.0
    rhs_with_coll = ((collision_op.ne.'none').and.(.not.coll_split))

    !precompute emfields_nlev for linearizing the nonlinear term.
    if(present(for_nlev)) then
       if (init_status_nlev.eq.0) then
          allocate(emfields_nlev(li1:li2,lj1:lj2,lbz:ubz,1:n_fields))
          emfields_nlev = 0.0
       endif
       init_status_nlev = 1

       call calc_aux_fields(g_1,emfields_nlev,f_,.false.)

       !Write(*,'(i2,a,ES20.8)') mype," precomputed emfields_nlev", sum(abs(emfields_nlev))
    endif

    
    if(present(for_petsc)) then
       with_g_update=.false.
    else
       with_g_update=low_mem
    endif

    if (rhs_nl) then
#ifdef WITHOMP_BLOCKLOOP
       !$OMP PARALLEL default(none) &
       !$OMP shared(xy_local,yx_order,x_local,arakawa,equil_par_curr,lbg0,mype,parallel_nl) &
       !$OMP private(type_of_nonlin,my_thread)
       my_thread = omp_get_thread_num()
       !$OMP CRITICAL
#else
       my_thread = 0
#endif
       if (.not.associated(this_nonlinear_term)) then
          if (xy_local) then
             ! flux tube
             if (yx_order) then
                allocate(ff_yx_nonlinear_term_t::this_nonlinear_term)
             else
                if (.not.nonlin_h) then
                   allocate(ff_nonlinear_term_t::this_nonlinear_term)
                else
                   allocate(ff_nonlinear_term_h_t::this_nonlinear_term)
                endif
             end if
          else
             !IF (.NOT.x_local) THEN
                ! x-global
                IF (arakawa) THEN
                   allocate(df_arakawa_nonlinear_term_t::this_nonlinear_term)
                else
#ifdef WITH_C_NONLIN
                   allocate(df_nonlinear_c_t::this_nonlinear_term)
#elif defined(WITH_SPLIT_NONLIN)
                   allocate(df_split_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITH_CUDA_NONLIN)
                   allocate(df_nonlinear_cuda_t::this_nonlinear_term)
#elif defined(WITH_MIC_NONLIN)
                   allocate(df_mic_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITHOMP_NONLIN)
                   allocate(df_omp_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITHOMPXY_NONLIN)
                   allocate(df_ompxy_nonlinear_term_t::this_nonlinear_term)
#else
                   allocate(df_nonlinear_term_t::this_nonlinear_term)
#endif
                END IF
             !end if
          end if
          call this_nonlinear_term%construct(equil_par_curr)
       end if

       select type (tmp=>this_nonlinear_term)
#ifdef WITH_SPLIT_NONLIN
       class is (df_split_nonlinear_term_t)
          if (gpu_cpu_ratio.lt.0.0) then
             call this_nonlinear_term%autotune(lbg0)
          else
             call this_nonlinear_term%set_gpu_cpu_ratio(gpu_cpu_ratio)
          end if
          !call this_nonlinear_term%set_blocksize(lbg0)
#endif
#ifdef WITH_CUDA_NONLIN
       class is (df_nonlinear_cuda_t)
          call this_nonlinear_term%SetCudaDevice(1)
#endif
       end select
#ifdef WITHOMP_BLOCKLOOP
       !$OMP END CRITICAL
#endif
       call this_nonlinear_term%set_blocksize(lbg0)
       call this_nonlinear_term%initialize()
       !call initialize_nonlinearity
       if (parallel_nl) call initialize_parallel_nonlin
       if ((mype.eq.0).and.(my_thread.eq.0).and.print_ini_msg) then
          type_of_nonlin = this_nonlinear_term%getType()
          write(*,"(3A,I8)") "nonlinear type is ",type_of_nonlin, &
               &" with a blocksize of ",this_nonlinear_term%lbg0
       end if
#ifdef WITHOMP_BLOCKLOOP
       !$OMP END PARALLEL
#endif
    endif

    rhs_with_coll = ((collision_op.ne.'none').and.(.not.coll_split))

    !if (only_coll) rhs_only_coll = only_coll

    ! following statement has been moved to discretization.F90
    !lbg0 = lklmn0/nblocks
    call initialize_klmn_conv_arrs

    if (init_status.eq.0) then
       if(with_g_update) then
          if(nblocks.eq.1) then
             allocate(rhs_big(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks))
          else 
             !with cache blocking
#ifdef WITHOMP_BLOCKLOOP
             !$OMP PARALLEL
             !$OMP CRITICAL
#endif
             allocate(rhs_small(li1:li2,lj1:lj2,1:lbg0))
#ifdef WITHOMP_BLOCKLOOP
             !$OMP END CRITICAL
             !$OMP END PARALLEL
#endif
             if (rhs_with_coll) then
                !collisions are performed without cache blocking
                allocate(rhs_big(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks))
             endif
          end if
       endif

       if (.not.xy_local) then
          allocate(bar_emfields(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2,ln1:ln2,1:n_fields))
          if (is_grid_adptv) bar_emfields = 0.0
          ptr_bar_emfields => bar_emfields
#ifdef WITHOMP_BLOCKLOOP
          !$OMP PARALLEL default(none) shared(rhs_nl,lbi,ubi,li1,li2,lj1,lj2,lbg0,is_grid_adptv)
          !$OMP CRITICAL
#endif
          allocate(g_block_storage(lbi:ubi,lj1:lj2,1:lbg0))
          allocate(chi_block(lbi:ubi,lj1:lj2,1:lbg0))
          if (is_grid_adptv) chi_block = 0.0
          if (rhs_nl) then
             !allocate(p_dbarchidxy(li1:li2,lj1:lj2,2,1:lbg0))
             !allocate(p_dgdxy(li1:li2,lj1:lj2,2,1:lbg0))
             !print*,"Trying to allocate p_dbarchidxy and p_dgdxy"
             call this_nonlinear_term%allocate_arrays(p_dbarchidxy,p_dgdxy)
             !print*,allocated(p_dbarchidxy),allocated(p_dgdxy)
             ptr_dbarchidxy => p_dbarchidxy
             ptr_dgdxy => p_dgdxy
          else
             NULLIFY(ptr_dbarchidxy)
             NULLIFY(ptr_dgdxy)
          endif
#ifdef WITHOMP_BLOCKLOOP
          !$OMP END CRITICAL
          !$OMP END PARALLEL
#endif
       endif
       init_status = 1
    endif

#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL default(none) shared(rhs_nl,li1,li2,lj1,lj2,lbg0,perf_vec)
    !$OMP CRITICAL
#endif
    if ((perf_vec(2).eq.1).and.rhs_nl) then
       if (.not.allocated(p_barchi)) allocate(p_barchi(li1:li2,lj1:lj2,1:lbg0))
       ptr_barchi => p_barchi
    else
       if (allocated(p_barchi)) deallocate(p_barchi)
       nullify(ptr_barchi)
    endif
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END CRITICAL
    !$OMP END PARALLEL
#endif

    if (collision_op.ne.'none') call initialize_add_coll 

    dfdzv_replace_rhs=.true.
    if (rhs_with_coll.and.(.not.with_g_update)) dfdzv_replace_rhs=.false.
    !because collisions initialize a_rhs and rhs_block points on a block of a_rhs.
    !(with g_update, an additional rhs_big contains collisions so that rhs_block 
    !must be initialized by the dfdzv_terms)

    if (.not.x_local.and.explicit_buffer) call initialize_krookBuffer_operator

    if(.not.arakawa_zv) then
       call initialize_dfdzv(dfdzv_replace_rhs)
       call initialize_dchidz
    else
       call initialize_dzv(dfdzv_replace_rhs)
    endif
    if((hyp_x.gt.0).or.(hyp_y.gt.0).or.(hyp_perp.gt.0).or.Erad.ne.0.or.(GyroLES)) then
       with_dfdperp=.true.
    else
       with_dfdperp=.false.
    end if

    if (with_dfdperp) call initialize_dfdxy

    call initialize_dchidxy
    
    call initialize_dgdxy

    if (rhs_f0) call initialize_f0_term

    if (with_sources) call initialize_f1_sources

  end subroutine initialize_CalFullRhs
  

  !>Computes the full right hand side of the gyrokinetic equation for a given g_1. First the
  !!fields f_ and emfields are computed, then the right hand side of the gyrokinetic equation
  !!(see thesis tbg) is computed with the routine chosen implementation
  Subroutine CalFullRhs(p_g_,rhs,stage)
    ! INPUT
    ! p_g_ : distribution function g
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(INOUT):: p_g_
    !****
    ! rhs is the right hand side, which is to be calculated from
    ! the given quantities
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(INOUT) :: rhs
    integer:: stage
    real :: local_sum,global_sum
    logical,parameter :: OUTPUT=.false.

    if (rhs_only_coll) then
       call calc_aux_fields(p_g_,emfields,f_,comp_h=.false.)
       call equ_collisions(f_,rhs,replace_rhs=.true.)
    else
       call calc_aux_fields(p_g_,emfields,f_,comp_h=arakawa_zv,h_out=h_)
#ifdef GIUTESTS
       ! alignment with adaptive version (only for testing!)
       if(is_grid_adptv) then
          call nullify_h_outside_adptv(h_)
       end if
#endif
       !Write(*,'(i2,a,ES20.8)') mype," this iteration emfields", sum(abs(emfields))
       call calc_rhs_only(f_, p_g_, emfields, rhs, stage, h_)
    endif

    IF (OUTPUT) THEN
       CALL calculate_test_sum(f_(:,:,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "f_ = ",global_sum
    END IF

  end Subroutine CalFullRhs


  !>RHS without field solve for implicit shape conversion of p_g_, rhs, etc.
  !!shall be removed once strip mining has been applied to all parts of the 
  !!Vlasov equation
  Subroutine calc_rhs_only (p_f_, p_g_, p_emfields, a_rhs, stage, p_h_)
    ! INPUT
    ! p_f_ : distribution function f1
    ! p_g_ : distribution function g
    ! p_emfields : fields
    ! a_rhs: right hand side for the old RK schemes / a_ vector for the low memory schemes
    ! stage: RK stage
    ! p_h_ : distribution function h1 (relevant for the arakawa_zv scheme)
    Complex, Dimension(li1:li2, lj1:lj2, 1:lzvwn0),Intent(INOUT):: p_f_ 
    Complex, Dimension(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks),Intent(INOUT):: p_g_
    Complex, Dimension(li1:li2, lj1:lj2, lbz:ubz, 1:n_fields),Intent(INOUT):: p_emfields
    integer:: stage
    complex, dimension(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks),target,intent(inout) :: a_rhs
    Complex, Dimension(li1:li2, lj1:lj2, 1:lzvwn0),optional,Intent(INOUT):: p_h_ 
    complex, dimension(:,:,:,:),pointer:: rhs, p_a_
    complex, dimension(:,:,:),pointer:: rhs_block, g_block
    real :: iter_start_time, nl_start_time
    integer :: iblock, lbg1,lbg2, my_thread,klmn
    logical, parameter :: SPEEDCHECK=.false., DOUTPUT=.false.
    real :: local_sum

#ifdef WITHOMP_BLOCKLOOP
    integer :: i_thread,ierr
    integer :: omp_get_thread_num, omp_get_num_threads
#endif

    DEBUG(1,"Start RHS")
    PERFON('CalFRhs0')

    if(stage.eq.1) then 
       ve_max = 0.0
       ve_pnl_max = 0.0
    endif

    if (.not.xy_local) then
       call set_bar_emfields(p_emfields, bar_emfields)
       if ((p_has_0_mode).and.(heatsrc_on.or.partsrc_on)) then
          call compute_f_vabs_allspec(p_f_)
       end if
    endif
    
    if(with_g_update) p_a_=>a_rhs

    !collisions do not work with strip mining because of velocity space/species integrals
    if (rhs_with_coll) then
       if(with_g_update) then
          rhs=>rhs_big
       else
          rhs=>a_rhs
       end if
       call equ_collisions(p_f_,rhs,replace_rhs=.true.)
    endif

    !print*,"Just before the OpenMP loop."
    !computation is done in sub-blocks of g-like arrays to optimize cache usage

#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL default(none) &
    !$OMP shared(lbg0,with_g_update,arakawa_zv,x_local,explicit_buffer,with_dfdperp,mpi_comm_y) &
    !$OMP shared(rhs_f0,rhs_nl,rhs_with_coll,with_sources,precond_approx,nblocks,timescheme,stage,rkstages) &
    !$OMP shared(parallel_nl,collision_op,lij0,rhs,c_rk,dt,b_rk,a_rk,xy_local,li1,li2,lj1,lj2) &
    !$OMP shared(a_rhs,p_f_,p_h_,time,rhs_big,p_emfields,ptr_bar_emfields,p_g_,p_a_,pdg1di,pdg1dj,mype) &
    !$OMP shared(n_fields,sk,sl,sm,sn,vtvpar, mpi_comm_xy,n_procs_x,n_procs_y) &
    !$OMP shared(blk,bll,blm,bln,hyp_on_h,rhs_nlev,emfields_nlev) &
    !$OMP private(klmn,my_thread,iblock,lbg1,lbg2,rhs_block, g_block,iter_start_time,local_sum) &
    !$OMP private(blk1,blk2,bll1,bll2,blm1,blm2,bln1,bln2) &
    !$OMP private(ierr,i_thread)
    my_thread = omp_get_thread_num()
    if (n_procs_y.gt.1) then
       ! clone the y communicator for all OpenMP threads
       !$OMP DO ORDERED
       do i_thread=1,omp_get_num_threads()
          !$OMP ORDERED
          call mpi_comm_dup(mpi_comm_y,threadlocal_mpi_comm_y,ierr)
          !$OMP END ORDERED
       end do
       !$OMP END DO
    end if
#else
    my_thread = 0
#endif

    if(with_g_update) then
       if(nblocks.eq.1) then
          rhs_block=>rhs_big(:,:,:,1)
       else
          rhs_block=>rhs_small
       endif
    endif

    !print*,mype,": my_thread = ",my_thread
    if (.not.xy_local) g_block => g_block_storage
#ifdef WITHOMP_BLOCKLOOP
    !$OMP DO 
#endif

    do iblock=1,nblocks
       lbg1 =(iblock-1)*lbg0+1
       lbg2 =lbg1+lbg0-1
       blk1=blk(1,iblock)
       blk2=blk(2,iblock)
       bll1=bll(1,iblock)
       bll2=bll(2,iblock)
       blm1=blm(1,iblock)
       blm2=blm(2,iblock)
       bln1=bln(1,iblock)
       bln2=bln(2,iblock)
       if (SPEEDCHECK) iter_start_time = my_secnds(0.0)
       if (DOUTPUT) print*,mype,my_thread," : working on ",lbg1,lbg2

       ! do the boundary exchange in x direction once for the whole block in advance
       ! This improves performance due to lower latency impact and larger messages.
       ! But we need a copy of g and chi which contains a whole block.
       if (.not.xy_local) then
          PERFON('blex_chi')
          if (n_fields.gt.1) then
             do klmn=lbg1,lbg2
                chi_block(li1:li2,lj1:lj2,klmn-lbg1+1) = &
                     &ptr_bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),1) &
                     &- vTvpar(sl(klmn),sn(klmn)) * &
                     & ptr_bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),2) 
             end do
             call exchange_x(x_boundary_block,chi_block)
          else
             ! here chi is independent of v_par, therefore one exchange for each j,k,m,n
             ! is enough
             do klmn=lbg1,lbg2
                chi_block(li1:li2,lj1:lj2,klmn-lbg1+1) = &
                     &ptr_bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),1)
             end do
             call exchange_x(x_boundary_block,chi_block)
             if (DOUTPUT) then
                write(*,"(3I4,2(A,ES17.10))") &
                     &mype,my_thread,iblock," : chi_block_inner = ",&
                     &real(sum(conjg(chi_block(li1:li2,lj1:lj2,:))*chi_block(li1:li2,lj1:lj2,:))),&
                     &", chi_block_all = ",&
                     &real(sum(conjg(chi_block(:,:,:))*chi_block(:,:,:)))
             end if

          end if
          PERFOFF
       end if
       if (.not.xy_local) then
          PERFON('blex_g')
          g_block(li1:li2,lj1:lj2,1:lbg0) = p_g_(li1:li2,lj1:lj2,1:lbg0, iblock)
          call exchange_x(x_boundary_block,g_block)
          PERFOFF
       end if

       if (SPEEDCHECK) print*,mype,my_thread," : After blocked exchange",my_secnds(iter_start_time)
       if(.not.with_g_update) rhs_block=>a_rhs(:,:,:,iblock)

       if(arakawa_zv) then
          PERFON('dzv_ak')
          call equ_dzv(p_h_,rhs_block,lbg1,lbg2)
          PERFOFF
          if (hyp_on_h) then
             if (hypz_compensation) call equ_comp_hypz(p_emfields,ptr_bar_emfields,rhs_block,lbg1,lbg2)
          else
             PERFON('hyp_zv_ak')
             call add_hypz_ak(p_f_,rhs_block,lbg1,lbg2)         
             call add_hypv_ak(p_f_,rhs_block,lbg1,lbg2) 
             PERFOFF
          endif
       else
          PERFON('eq_dfdzv')
          ! standard case (should be threadsafe)
          call equ_dfdzv(p_f_,rhs_block,lbg1,lbg2)
          PERFOFF
       endif
       if (SPEEDCHECK) print*,mype,my_thread," : After equ_dfdzv:",my_secnds(iter_start_time)

       if (DOUTPUT) write(*,"(3I4,A,ES17.10)") &
            &mype,my_thread,iblock," : rhs_block(1) = ",&
            &real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))
       ! standard explicit_buffer=.false.
       if (.not.x_local.and.explicit_buffer) call add_kBuffer_explicit(p_f_,rhs_block,lbg1,lbg2)
       !>todo: check if perpendicular hyperdiffusion should act on f or h!
       if (with_dfdperp) then
          if (hyp_on_h) then 
             call add_dfdxy(p_h_,rhs_block,lbg1,lbg2) !We need to add hyp_x and hyp_y on h
          else
             call add_dfdxy(p_f_,rhs_block,lbg1,lbg2) !We need to add hyp_x and hyp_y on f
          endif
       endif
       if (rhs_f0) call add_f0_term(rhs_block,lbg1,lbg2,time)
       if (SPEEDCHECK) print*,mype,my_thread," : After dfdxy:",my_secnds(iter_start_time)

       if (DOUTPUT) write(*,"(3I4,A,ES17.10)") &
            &mype,my_thread,iblock," : rhs_block(2) = ",&
            &real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))

       ! add_dchidxy is threadsafe
       PERFON('dchidxy')
       if (xy_local) then
          call add_dchidxy_orig(p_emfields,ptr_bar_emfields,&
               &rhs_block,ptr_barchi,ptr_dbarchidxy,lbg1,lbg2,.true.)
       else
          call add_dchidxy(chi_block,rhs_block,ptr_dbarchidxy,lbg1,lbg2,.true.)
       end if
       PERFOFF
       if (SPEEDCHECK) print*,mype,my_thread," : After dchidxy:",my_secnds(iter_start_time)

       if (DOUTPUT) write(*,"(3I4,A,ES17.10)") &
            &mype,my_thread,iblock," : rhs_block(3) = ",real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))
       ! add_dchidz is threadsafe
       PERFON('dchidz')
       if (.not.arakawa_zv) call add_dchidz(p_emfields, ptr_bar_emfields, &
            &rhs_block, lbg1, lbg2)
       PERFOFF

       if (DOUTPUT) write(*,"(3I4,A,ES17.10)") &
            &mype,my_thread,iblock," : rhs_block(4) = ",&
            &real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))

       if (SPEEDCHECK) print*,mype,my_thread," : After add_dchidz ", my_secnds(iter_start_time)

       ! add_dgdxy is threadsafe
       PERFON('dgdxy')
       if (xy_local) then
          call add_dgdxy(p_g_, rhs_block, ptr_dgdxy, pdg1di, pdg1dj, lbg1,lbg2)
       else
          call add_dgdxy(g_block, rhs_block, ptr_dgdxy, pdg1di, pdg1dj, lbg1,lbg2)
       end if
       PERFOFF
       if (SPEEDCHECK) print*,mype,my_thread," : After   add_dgdxy", my_secnds(iter_start_time)

       if (DOUTPUT) then
          local_sum = real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))
          !call reduce_per_thread(local_sum,mpi_comm_xy,"global_sum(5)")
          write(*,"(3I4,A,ES17.10)") &
               &mype,my_thread,iblock," : rhs_block(5) = ",local_sum
       end if

       if (rhs_nl) then
          PERFON('add_nl')
#ifndef WITHOMP_BLOCKLOOP
          nl_start_time = my_secnds(0.0)
#endif
          if (xy_local) then
             if(rhs_nlev) then
                call this_nonlinear_term%add(p_g_, ptr_dgdxy, emfields_nlev,&
                     &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
             else
                if (.not.nonlin_h) then
                   call this_nonlinear_term%add(p_g_, ptr_dgdxy, p_emfields,&
                        &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
                else
                   call this_nonlinear_term%add(p_h_, ptr_dgdxy, p_emfields,&
                        &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
                endif
             endif
          else
             select type(tmp=>this_nonlinear_term)
             class is (df_arakawa_nonlinear_term_t)
             ! this copy operation is only necessary for arakawa scheme
             
             ptr_barchi(:,:,:) = chi_block(li1:li2,lj1:lj2,:)
             end select
             call this_nonlinear_term%add(g_block, ptr_dgdxy, p_emfields,&
                  &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
          end if
#ifndef WITHOMP_BLOCKLOOP
          !nl_sum_time = nl_sum_time + my_secnds(nl_start_time)
#endif
          PERFOFF
          if (parallel_nl) then
             PERFON('paral_nl')
             if (rhs_nlev) then
                call add_parallel_nonlin(p_f_,emfields_nlev,&
                     &ptr_bar_emfields, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
             else
                call add_parallel_nonlin(p_f_,p_emfields,&
                     &ptr_bar_emfields, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
             endif
             PERFOFF
          end if
       endif

       PERFON('dApardt')
       if (antenna_type.eq.3.and.antenna_contrib) call add_dApar_dt_antenna(rhs_block,lbg1,lbg2)
       PERFOFF

       if (DOUTPUT) then
          local_sum = real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))
          
          !call reduce_per_thread(local_sum,mpi_comm_xy,"global_sum(6)")
          write(*,"(3I4,A,ES17.10)") &
               &mype,my_thread,iblock," : rhs_block(6) = ",local_sum
       end if
       if (SPEEDCHECK) print*,mype,my_thread," : After   nonlinearity", my_secnds(iter_start_time)


       if (with_sources.and.(.not.precond_approx)) then
          PERFON('sources')
          call add_f1_sources(rhs_block,lbg1,lbg2,stage)
          PERFOFF
       end if

       if(with_g_update) then
          PERFON('up_ag')
          !add the collision contribution if rhs_block.ne.rhs

          if ((nblocks.gt.1).and.(rhs_with_coll)) then
             call axpy_ij(lij0*lbg0,1.0,rhs(:,:,:,iblock),rhs_block)
             ! rhs_block=rhs_block+rhs(:,:,:,iblock)
          end if

          if(timescheme.eq.'turbo') then
             if (stage.eq.1) then
                p_a_(:,:,:,iblock)=c_rk(stage)*rhs_block
             else
                p_a_(:,:,:,iblock)=p_a_(:,:,:,iblock)+c_rk(stage)*rhs_block
             endif
             p_g_(:,:,:,iblock)=p_g_(:,:,:,iblock)+dt*b_rk(stage)*p_a_(:,:,:,iblock)+&
                  dt*a_rk(stage)*rhs_block
          else
             p_g_(:,:,:,iblock)=p_a_(:,:,:,iblock)+a_rk(stage)*dt*rhs_block
             if(stage.lt.rkstages) then
                call axpy_ij(lij0*lbg0,b_rk(stage)*dt,rhs_block,&
                     &p_a_(:,:,:,iblock))
                !                p_a_(:,:,:,iblock)=p_a_(:,:,:,iblock)+b_rk(stage)*dt*rhs_block
             else
                call ccopy(lij0*lbg0,p_g_(:,:,:,iblock),1,p_a_(:,:,:,iblock),1) 
!                p_a_(:,:,:,iblock)=p_g_(:,:,:,iblock)
             endif
          end if
          PERFOFF
       endif
    enddo
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END DO
    if (n_procs_y.gt.1) then
       ! clone the y communicator for all OpenMP threads
       !$OMP DO ORDERED
       do i_thread=1,omp_get_num_threads()
          !$OMP ORDERED
          call mpi_comm_free(threadlocal_mpi_comm_y,ierr)
          !$OMP END ORDERED
       end do
       !$OMP END DO
    end if
    !$OMP END PARALLEL
#endif
    !print*,"Direct after the OpenMP loop."
    PERFOFF 
    DEBUG(1,"End RHS")

  End Subroutine calc_rhs_only


  !>Wrapper for the cleanup routines of the various implementations
  subroutine finalize_CalFullRhs
    
    if(init_status_nlev.gt.0) deallocate(emfields_nlev)
    init_status_nlev=0

#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL
#endif
    if (rhs_nl) then
       if (allocated(p_barchi)) deallocate(p_barchi)
       nullify(ptr_barchi)
    endif
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END PARALLEL
#endif

    if (.not.x_local.and.explicit_buffer) call finalize_krookBuffer_operator

    if (with_sources) call finalize_f1_sources

    if (.not.xy_local) then
       deallocate(bar_emfields)
#ifdef WITHOMP_BLOCKLOOP
       !$OMP PARALLEL
#endif
       deallocate(g_block_storage)
       deallocate(chi_block)
       if (rhs_nl) then
          call this_nonlinear_term%free_arrays(p_dbarchidxy,p_dgdxy)
          !deallocate(p_dbarchidxy,p_dgdxy)
          nullify(ptr_dbarchidxy,ptr_dgdxy)
       endif
#ifdef WITHOMP_BLOCKLOOP
       !$OMP END PARALLEL
#endif
    endif

    if(arakawa_zv) then
       call finalize_dzv
    else
       call finalize_dfdzv
       call finalize_dchidz 
    end if

    if (collision_op.ne.'none') call finalize_add_coll
    if (with_dfdperp) call finalize_dfdxy
    
    if(rhs_f0) call finalize_f0_term

    call finalize_dchidxy


    call finalize_dgdxy
    if (rhs_nl) then
#ifdef WITHOMP_BLOCKLOOP
       !$OMP PARALLEL
#endif
       if (parallel_nl) call finalize_parallel_nonlin
       call this_nonlinear_term%finalize()
       call this_nonlinear_term%destruct()
       deallocate(this_nonlinear_term)
#ifdef WITHOMP_BLOCKLOOP
       !$OMP END PARALLEL
#endif
    endif

    call finalize_klmn_conv_arrs

    if(with_g_update) then
       if(nblocks.eq.1) then
          deallocate(rhs_big)
       else 
#ifdef WITHOMP_BLOCKLOOP
          !$OMP PARALLEL
          !$OMP CRITICAL
#endif
          deallocate(rhs_small)
#ifdef WITHOMP_BLOCKLOOP
          !$OMP END CRITICAL
          !$OMP END PARALLEL
#endif
          if (rhs_with_coll) deallocate(rhs_big)
       end if
    endif

    init_status = 0
    !if (mype.eq.0) write(*,"(A,F10.3)") "nl_sum_time = ",nl_sum_time
  end subroutine finalize_CalFullRhs


  subroutine set_bar_emfields(p_emfields,p_bar_emfields)
    Complex, Dimension(li1:li2, lj1:lj2, lbz:ubz, 1:n_fields),Intent(IN):: p_emfields
    Complex, Dimension(li1:li2, lj1:lj2, lbz:ubz, lm1:lm2, ln1:ln2, 1:n_fields),Intent(OUT):: p_bar_emfields

    integer :: m,n,o

    do o=1,n_fields
       do n=ln1,ln2
          do m=lm1,lm2
             call gyro_average_df_wb(p_emfields(:,:,:,o),&
                  &p_bar_emfields(:,:,:,m,n,o),m,n)

             ! boundary exchange only for phi, as we only need derivatives 
             ! of phi in x,y and z direction
             if (o==1) call exchange_z(p_bar_emfields(:,:,:,m,n,1))
          enddo
       enddo
    enddo

  end subroutine set_bar_emfields

  !--------- modified subroutines for adaptive grids ---------
  
  !>Computes the full right hand side of the gyrokinetic equation for a given g_1. First the
  !!fields f_ and emfields are computed, then the right hand side of the gyrokinetic equation
  !!(see thesis tbg) is computed with the routine chosen implementation
  Subroutine CalFullRhs_adptv(p_g_,rhs,stage)
    ! INPUT
    ! p_g_ : distribution function g
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(INOUT):: p_g_
    !****
    
    ! rhs is the right hand side, which is to be calculated from
    ! the given quantities
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(INOUT) :: rhs
    integer:: stage
    real :: local_sum,global_sum
    logical,parameter :: OUTPUT=.false.

    if (rhs_only_coll) then
       stop "collision operator is not implemented for the adaptive grids yet"
!       call calc_aux_fields_adptv(p_g_,emfields,f_,comp_h=.false.)
!       call equ_collisions(f_,rhs,replace_rhs=.true.)
    else
       if (arakawa_zv) then
          !compute f_, emfields and h_
          call calc_aux_fields_adptv(p_g_,emfields,f_,.true.,h_)
          !call calc_aux_fields(p_g_,emfields,f_,.true.,h_)
          call calc_rhs_only_adptv(f_, p_g_, emfields, rhs, stage, h_)
       else
          !compute f_ and emfields (h_ is not needed)
          call calc_aux_fields_adptv(p_g_,emfields,f_,.false.)
          call calc_rhs_only_adptv(f_, p_g_, emfields, rhs, stage)
       endif
    endif

    IF (OUTPUT) THEN
       CALL calculate_test_sum(f_(:,:,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "f_ = ",global_sum
    END IF

  end Subroutine CalFullRhs_adptv

  !>RHS without field solve for implicit shape conversion of p_g_, rhs, etc.
  !!shall be removed once strip mining has been applied to all parts of the 
  !!Vlasov equation
  Subroutine calc_rhs_only_adptv (p_f_, p_g_, p_emfields, a_rhs, stage, p_h_)
    Complex, Dimension(li1:li2, lj1:lj2, 1:lzvwn0),Intent(INOUT):: p_f_ 
    Complex, Dimension(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks),Intent(INOUT):: p_g_
    Complex, Dimension(li1:li2, lj1:lj2, lbz:ubz, 1:n_fields),Intent(INOUT):: p_emfields
    integer:: stage
    complex, dimension(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks),target,intent(inout) :: a_rhs
    Complex, Dimension(li1:li2, lj1:lj2, 1:lzvwn0),optional,Intent(INOUT):: p_h_ 

    complex, dimension(:,:,:,:),pointer:: rhs, p_a_
    complex, dimension(:,:,:),pointer:: rhs_block, g_block
    real :: iter_start_time
    integer :: iblock, lbg1,lbg2, my_thread,klmn
    logical, parameter :: SPEEDCHECK=.false., DOUTPUT=.false.
    real :: local_sum
#ifdef WITHOMP
    integer :: i_thread,ierr
    integer :: omp_get_thread_num, omp_get_num_threads
#endif
    ! necessary for adaptive version
    integer :: l, m, i1, i2

    DEBUG(1,"Start RHS")
    PERFON('CalFRhs0')

    if(stage.eq.1) then
       ve_max = 0.0
       !ve_max_avg = 0.0
       ve_pnl_max = 0.0
    endif

    if (.not.xy_local) then
       call set_bar_emfields(p_emfields, bar_emfields)
       if ((p_has_0_mode).and.(heatsrc_on.or.partsrc_on)) then
          ! block-structured stuff is inside the subroutine (TODO: reconsider)
          call compute_f_vabs_allspec(p_f_)
       end if
    endif
    
    if(with_g_update) p_a_=>a_rhs

    !collisions do not work with strip mining because of velocity space/species integrals
    if (rhs_with_coll) then
       stop "collisions are not implemented for adaptive grids yet"
       if(with_g_update) then
          rhs=>rhs_big
       else
          rhs=>a_rhs
       end if
       call equ_collisions(p_f_,rhs,replace_rhs=.true.)
    endif

#ifdef WITHOMP_BLOCKLOOP
    stop "no shared memory parallelization for adaptive grids yet"
#endif
    
    if(with_g_update) then
       if(nblocks.eq.1) then
          rhs_block=>rhs_big(:,:,:,1)
       else
          rhs_block=>rhs_small
       endif
    endif

    my_thread = 0
    
    if (.not.xy_local) g_block => g_block_storage
    
    do iblock=1,nblocks
       lbg1 =(iblock-1)*lbg0+1
       lbg2 =lbg1+lbg0-1
       blk1=blk(1,iblock)
       blk2=blk(2,iblock)
       bll1=bll(1,iblock)
       bll2=bll(2,iblock)
       blm1=blm(1,iblock)
       blm2=blm(2,iblock)
       bln1=bln(1,iblock)
       bln2=bln(2,iblock)
       if (SPEEDCHECK) iter_start_time = my_secnds(0.0)
       if (DOUTPUT) print*,mype,my_thread," : working on ",lbg1,lbg2

       ! do the boundary exchange in x direction once for the whole block in advance
       ! This improves performance due to lower latency impact and larger messages.
       ! But we need a copy of g and chi which contains a whole block.
       if (.not.xy_local) then
          PERFON('blex_chi')
          if (n_fields.gt.1) then
             do klmn=lbg1,lbg2
                l = sl(klmn)
                m = sm(klmn)
                i1 = li1_vwadp(l,m) - radscheme/2
                if (i1.lt.li1) i1 = li1
                i2 = li2_vwadp(l,m) + radscheme/2
                if (i2.gt.li2) i2 = li2
                chi_block(i1:i2,lj1:lj2,klmn-lbg1+1) = &
                   ptr_bar_emfields(i1:i2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),1) &
                     - vTvpar(sl(klmn),sn(klmn)) * &
                     ptr_bar_emfields(i1:i2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),2)
             end do
             call exchange_x(x_boundary_block,chi_block)
          else
             ! here chi is independent of v_par, therefore one exchange for each j,k,m,n
             ! is enough
             do klmn=lbg1,lbg2
                l = sl(klmn)
                m = sm(klmn)
                i1 = li1_vwadp(l,m) - radscheme/2
                if (i1.lt.li1) i1 = li1
                i2 = li2_vwadp(l,m) + radscheme/2
                if (i2.gt.li2) i2 = li2
                chi_block(i1:i2,lj1:lj2,klmn-lbg1+1) = &
                     &ptr_bar_emfields(i1:i2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),1)
             end do
             call exchange_x(x_boundary_block,chi_block)
             if (DOUTPUT) then
                stop "output is not implemented for the adaptive grid!"
                write(*,"(3I4,2(A,ES17.10))") &
                     &mype,my_thread,iblock," : chi_block_inner = ",&
                     &real(sum(conjg(chi_block(li1:li2,lj1:lj2,:))*chi_block(li1:li2,lj1:lj2,:))),&
                     &", chi_block_all = ",&
                     &real(sum(conjg(chi_block(:,:,:))*chi_block(:,:,:)))
             end if

          end if
          PERFOFF
       end if
       
       if (.not.xy_local) then
          PERFON('blex_g')
          ! \todo -> adaptivity : should be rewritten 
          g_block(li1:li2,lj1:lj2,1:lbg0) = p_g_(li1:li2,lj1:lj2,1:lbg0, iblock)
          call exchange_x(x_boundary_block,g_block)
          PERFOFF
       end if

       if (SPEEDCHECK) print*,mype,my_thread," : After blocked exchange",my_secnds(iter_start_time)
       if(.not.with_g_update) rhs_block=>a_rhs(:,:,:,iblock)

       if(arakawa_zv) then
          PERFON('dzv_ak')
          call equ_dzv_adptv(p_h_,rhs_block,lbg1,lbg2)
          PERFOFF
          if (hyp_on_h) then
             if (hypz_compensation) call equ_comp_hypz_adptv(&
                  &p_emfields,ptr_bar_emfields,rhs_block,lbg1,lbg2)
          else
             PERFON('hyp_zv_ak')
             call add_hypz_ak_adptv(p_f_,rhs_block,lbg1,lbg2)         
             call add_hypv_ak_adptv(p_f_,rhs_block,lbg1,lbg2) 
             PERFOFF
          endif
       else
          PERFON('eq_dfdzv')
          ! standard case (should be threadsafe)
          call equ_dfdzv_adptv(p_f_,rhs_block,lbg1,lbg2)
          PERFOFF
       endif
       if (SPEEDCHECK) print*,mype,my_thread," : After equ_dfdzv:",my_secnds(iter_start_time)

       if (DOUTPUT) write(*,"(3I4,A,ES17.10)") &
            &mype,my_thread,iblock," : rhs_block(1) = ",&
            &real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))
       ! standard explicit_buffer=.false.
       if (.not.x_local.and.explicit_buffer) &
            call add_kBuffer_explicit_adptv(p_f_,rhs_block,lbg1,lbg2)
       !>todo: check if perpendicular hyperdiffusion should act on f or h!
       if (with_dfdperp) then
          if (hyp_on_h) then 
             call add_dfdxy_adptv(p_h_,rhs_block,lbg1,lbg2) !We need to add hyp_x and hyp_y on h
          else
             call add_dfdxy_adptv(p_f_,rhs_block,lbg1,lbg2) !We need to add hyp_x and hyp_y on f
          endif
       endif
       if (rhs_f0) call add_f0_term_adptv(rhs_block,lbg1,lbg2,time)
       if (SPEEDCHECK) print*,mype,my_thread," : After dfdxy:",my_secnds(iter_start_time)

       if (DOUTPUT) write(*,"(3I4,A,ES17.10)") &
            &mype,my_thread,iblock," : rhs_block(2) = ",&
            &real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))

       ! add_dchidxy is threadsafe
       PERFON('dchidxy')
       if (xy_local) then
          stop "no adaptive grids for local xy version!"
       else
          call add_dchidxy_adptv(chi_block,rhs_block,ptr_dbarchidxy,lbg1,lbg2,.true.)
       end if
       PERFOFF
       if (SPEEDCHECK) print*,mype,my_thread," : After dchidxy:",my_secnds(iter_start_time)

       if (DOUTPUT) write(*,"(3I4,A,ES17.10)") &
            &mype,my_thread,iblock," : rhs_block(3) = ",real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))
       ! add_dchidz is threadsafe
       PERFON('dchidz')
       if(.not.arakawa_zv) call add_dchidz_adptv(p_emfields, ptr_bar_emfields, &
            &rhs_block, lbg1, lbg2)
       PERFOFF

       if (DOUTPUT) write(*,"(3I4,A,ES17.10)") &
            &mype,my_thread,iblock," : rhs_block(4) = ",&
            &real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))

       if (SPEEDCHECK) print*,mype,my_thread," : After add_dchidz ", my_secnds(iter_start_time)

       ! add_dgdxy is threadsafe
       PERFON('dgdxy')
       if (xy_local) then
          stop "no adaptive grids for xy local version!"
       else
          call add_dgdxy_adptv(g_block, rhs_block, ptr_dgdxy, pdg1di, pdg1dj, lbg1,lbg2)
       end if
       PERFOFF
       if (SPEEDCHECK) print*,mype,my_thread," : After   add_dgdxy", my_secnds(iter_start_time)

       if (DOUTPUT) then
          local_sum = real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))
          
          !call reduce_per_thread(local_sum,mpi_comm_xy,"global_sum(5)")
          write(*,"(3I4,A,ES17.10)") &
               &mype,my_thread,iblock," : rhs_block(5) = ",local_sum
       end if

       if (rhs_nl) then
          stop "non-linearity on adaptive grids is not implemented yet"
          PERFON('add_nl')
          if (xy_local) then
             stop "no adaptive grids for xy local version"
          else
             select type(tmp=>this_nonlinear_term)
             class is (df_arakawa_nonlinear_term_t)
             ! this copy operation is only necessary for arakawa scheme
             
             ptr_barchi(:,:,:) = chi_block(li1:li2,lj1:lj2,:)
             end select
             call this_nonlinear_term%add(g_block, ptr_dgdxy, p_emfields,&
                  &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
          end if
          PERFOFF
          if (parallel_nl) then
             PERFON('paral_nl')
             if (rhs_nlev) then
                call add_parallel_nonlin(p_f_,emfields_nlev,&
                     &ptr_bar_emfields, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
             else
                call add_parallel_nonlin(p_f_,p_emfields,&
                     &ptr_bar_emfields, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
             endif
             PERFOFF
          end if
       endif

       if (DOUTPUT) then
          local_sum = real(sum(conjg(rhs_block(:,:,:))*rhs_block(:,:,:)))
          
          !call reduce_per_thread(local_sum,mpi_comm_xy,"global_sum(6)")
          write(*,"(3I4,A,ES17.10)") &
               &mype,my_thread,iblock," : rhs_block(6) = ",local_sum
       end if
       if (SPEEDCHECK) print*,mype,my_thread," : After   nonlinearity", my_secnds(iter_start_time)

       if (with_sources.and.(.not.precond_approx)) then
          PERFON('sources')
          call add_f1_sources_block_adptv(rhs_block,lbg1,lbg2,stage)
          PERFOFF
       end if

       ! \todo adaptive grids are not implemented for "with_g_update"
       if(with_g_update) then
          PERFON('up_ag')
          !add the collision contribution if rhs_block.ne.rhs

          if ((nblocks.gt.1).and.(rhs_with_coll)) then
             call axpy_ij(lij0*lbg0,1.0,rhs(:,:,:,iblock),rhs_block)
             ! rhs_block=rhs_block+rhs(:,:,:,iblock)
          end if

          if(timescheme.eq.'turbo') then
             if (stage.eq.1) then
                p_a_(:,:,:,iblock)=c_rk(stage)*rhs_block
             else
                p_a_(:,:,:,iblock)=p_a_(:,:,:,iblock)+c_rk(stage)*rhs_block
             endif
             p_g_(:,:,:,iblock)=p_g_(:,:,:,iblock)+dt*b_rk(stage)*p_a_(:,:,:,iblock)+&
                  dt*a_rk(stage)*rhs_block
          else
             p_g_(:,:,:,iblock)=p_a_(:,:,:,iblock)+a_rk(stage)*dt*rhs_block
             if(stage.lt.rkstages) then
                call axpy_ij(lij0*lbg0,b_rk(stage)*dt,rhs_block,&
                     &p_a_(:,:,:,iblock))
                !                p_a_(:,:,:,iblock)=p_a_(:,:,:,iblock)+b_rk(stage)*dt*rhs_block
             else
                call ccopy(lij0*lbg0,p_g_(:,:,:,iblock),1,p_a_(:,:,:,iblock),1) 
!                p_a_(:,:,:,iblock)=p_g_(:,:,:,iblock)
             endif
          end if
          PERFOFF
       endif
    enddo
#ifdef WITHOMP_BLOCKLOOP
   !$OMP END DO
   if (n_procs_y.gt.1) then
      ! clone the y communicator for all OpenMP threads
      !$OMP DO ORDERED
      do i_thread=1,omp_get_num_threads()
         !$OMP ORDERED
         call mpi_comm_free(threadlocal_mpi_comm_y,ierr)
         !$OMP END ORDERED
      end do
      !$OMP END DO
   end if
   !$OMP END PARALLEL
#endif
    !print*,"Direct after the OpenMP loop."
    PERFOFF 
    DEBUG(1,"End RHS")

  End Subroutine calc_rhs_only_adptv

#ifdef GIUTESTS
 subroutine nullify_h_outside_adptv(h_type)
   ! use discretization, only: li1, li2, lj1, lj2, &
   !       lk1, lk2, ll1, ll2, lm1, lm2, ln1, ln2
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), &
         intent(inout) :: h_type
    integer :: i, j, k, l, m, n

    do n = ln1, ln2
       do m = lbw, ubw
          do l = lbv, ubv
             do k = lbz, ubz
                do j = lj1, lj2
                   do i = li1, li1_b_vwadp(l,m)-1
                      h_type(i,j,k,l,m,n) = 0.
                   enddo
                   do i = li2_b_vwadp(l,m)+1, li2
                      h_type(i,j,k,l,m,n) = 0.
                   enddo
                enddo
             enddo
          endDo
       endDo
    enddo
  end subroutine nullify_h_outside_adptv
#endif

End Module calc_rhs
