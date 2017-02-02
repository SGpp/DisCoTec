#include "redef.h"
#include "switches.h"

module df_split_nonlinear_term_mod
  use,intrinsic :: iso_c_binding
  use communications
  use df_nonlinear_term_mod
  use df_omp_nonlinear_term_mod
  use df_nonlinear_cuda_mod
  use discretization, only: li0,li1,li2, lj0,lj1,lj2
  implicit none

  type,public, extends(df_nonlinear_term_t) :: df_split_nonlinear_term_t
     class(df_nonlinear_cuda_t),pointer :: gpu_nonlinearity
     class(df_nonlinear_term_t),pointer :: cpu_nonlinearity
     real :: gpu_cpu_ratio
     integer :: gpu_block
   contains
     procedure :: set_blocksize => set_split_blocksize
     procedure :: set_gpu_cpu_ratio
     procedure :: calc => calc_split_nonlinearity_df

     procedure :: initialize => initialize_split_nonlinear_term
     procedure :: finalize => finalize_split_nonlinear_term

     procedure :: construct => construct_df_split_nonlinear_term
     procedure :: destruct => destruct_df_split_nonlinear_term

     procedure :: mem_est => mem_est_split_nonlinearity
     procedure :: getType => getThisType
     
     procedure :: allocate_arrays => allocate_split_arrays
     procedure :: free_arrays => free_split_arrays
     procedure :: autotune => autotune_split_nonlinearity_df
#if 0
     procedure :: add => add_nonlinearity
#endif
  end type df_split_nonlinear_term_t

  private

contains

  function getThisType(this)
    class(df_split_nonlinear_term_t) :: this
    character(len=MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "df_split_nonlinear_term_t"
  end function getThisType
    
  subroutine construct_df_split_nonlinear_term(this,equil_par_curr)
    class(df_split_nonlinear_term_t) :: this
    logical :: equil_par_curr

    allocate(df_nonlinear_cuda_t::this%gpu_nonlinearity)
    call this%gpu_nonlinearity%construct(equil_par_curr)
    allocate(df_omp_nonlinear_term_t::this%cpu_nonlinearity)
    call this%cpu_nonlinearity%construct(equil_par_curr)

  end subroutine construct_df_split_nonlinear_term

  subroutine destruct_df_split_nonlinear_term(this)
    class(df_split_nonlinear_term_t) :: this

    call this%gpu_nonlinearity%destruct()
    deallocate(this%gpu_nonlinearity)
    call this%cpu_nonlinearity%destruct()
    deallocate(this%cpu_nonlinearity)
  end subroutine destruct_df_split_nonlinear_term

  Real function mem_est_split_nonlinearity(this,mem_req_in)
    class(df_split_nonlinear_term_t) :: this
    real, intent(IN) :: mem_req_in
    real :: mem_loc
    
    mem_loc = this%gpu_nonlinearity%mem_est(mem_req_in)
    mem_est_split_nonlinearity = this%cpu_nonlinearity%mem_est(mem_loc)
    
  end function mem_est_split_nonlinearity

  subroutine set_gpu_cpu_ratio(this,gc_ratio)
    class(df_split_nonlinear_term_t) :: this
    real :: gc_ratio

    if (gc_ratio.gt.0.0) then
       this%gpu_cpu_ratio = gc_ratio
    else
       print*,"The GPU_CPU ratio should be larger than 0."
       this%gpu_cpu_ratio = 0.0
    end if
  end subroutine set_gpu_cpu_ratio

  subroutine set_split_blocksize(this,blocksize)
    class(df_split_nonlinear_term_t) :: this
    integer :: blocksize
    
    integer :: cpu_block

    !call this%set_gpu_cpu_ratio(3.0)

    this%lbg0=blocksize
    !print*,"split: this%lbg0 = ",this%lbg0
    ! now set the blocksizes of the two sub-nonlinearities
    this%gpu_block = this%lbg0/(this%gpu_cpu_ratio+1.0) * this%gpu_cpu_ratio
    this%gpu_block = this%gpu_nonlinearity%get_nearest_blocksize(this%gpu_block)
    call this%gpu_nonlinearity%set_blocksize(this%gpu_block)
    cpu_block = this%lbg0-this%gpu_block
    call this%cpu_nonlinearity%set_blocksize(cpu_block)
    !print*,"after setting sub-nonlins: this%lbg0 = ",this%lbg0
  end subroutine set_split_blocksize

  subroutine initialize_split_nonlinear_term(this)
    class(df_split_nonlinear_term_t) :: this

    if (this%init_status==1) return
    call this%gpu_nonlinearity%SetCudaDevice(1)

    call this%gpu_nonlinearity%initialize()
    call this%cpu_nonlinearity%initialize()

    this%init_status = 1
  end subroutine initialize_split_nonlinear_term

  subroutine finalize_split_nonlinear_term(this)
    class(df_split_nonlinear_term_t) :: this

    call this%cpu_nonlinearity%finalize()
    call this%gpu_nonlinearity%finalize()

    this%init_status = 0
  end subroutine finalize_split_nonlinear_term

  subroutine allocate_split_arrays(this,vexy,dgdxy)
    class(df_split_nonlinear_term_t) :: this
    complex, dimension(:,:,:,:),allocatable :: vexy,dgdxy

    ! we allocate the whole array, to have a contiguous block
    ! which can eventually be used in the other term routines.
    allocate(vexy(li1:li2,lj1:lj2,2,1:this%lbg0))
    allocate(dgdxy(li1:li2,lj1:lj2,2,1:this%lbg0))

    ! this memory has now to be page-locked, by
    ! registering
    call this%gpu_nonlinearity%register_array(vexy,li0*lj0*2*this%lbg0)
    call this%gpu_nonlinearity%register_array(dgdxy,li0*lj0*2*this%lbg0)

    ! But for the two sub-nonlinearities, we have to use only
    ! the respective parts of the arrays.

    call this%gpu_nonlinearity%allocate_arrays_on_device()

  end subroutine allocate_split_arrays

  subroutine free_split_arrays(this,vexy,dgdxy)
    class(df_split_nonlinear_term_t) :: this
    complex, dimension(:,:,:,:),allocatable :: vexy,dgdxy

    call this%gpu_nonlinearity%free_arrays_on_device()
    call this%gpu_nonlinearity%unregister_array(vexy)
    call this%gpu_nonlinearity%unregister_array(dgdxy)

    deallocate(vexy)
    deallocate(dgdxy)
  end subroutine free_split_arrays

  Subroutine calc_split_nonlinearity_df(this,gy_chi,g_block,vexy,dgdxy,localrhs,first)
    class(df_split_nonlinear_term_t) :: this
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: gy_chi
    complex, dimension(lbi:ubi,lj1:lj2,1:*),intent(inout) :: g_block
    Complex, Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
    Logical, Intent(in) :: first

#ifdef WITHOMP_SPLIT_NONLIN
    integer :: max_threads, omp_get_max_threads, omp_get_thread_num, thread_id
    max_threads = omp_get_max_threads()
    !print*,"We have max_threads = ",max_threads
    !print*,"gpu_block = ",this%gpu_block
    
    !$OMP PARALLEL default(none) &
    !$OMP shared(this,gy_chi,g_block,vexy,dgdxy,localrhs,first,max_threads) &
    !$OMP private(thread_id) &
    !$OMP num_threads(2)
    !print*,"Outer parallel region: thread_num = ", omp_get_thread_num()
#endif
    thread_id = omp_get_thread_num()
    if (thread_id.eq.0) then
       !print*,"Calling calc for gpu on thread ", thread_id
       call this%gpu_nonlinearity%calc(gy_chi(:,:,1:this%gpu_block),&
            &g_block(:,:,1:this%gpu_block),&
            &vexy(:,:,:,1:this%gpu_block),&
            &dgdxy(:,:,:,1:this%gpu_block),&
            &localrhs(:,:,1:this%gpu_block),first)
    else
#if 0
       !$OMP SECTION
       !$OMP PARALLEL if (max_threads>2) &
       !$OMP default(none) &
       !$OMP shared(this,gy_chi,g_block,vexy,dgdxy,localrhs,first) &
       !$OMP num_threads(max_threads-1)
       print*,"Inner parallel region: thread_num = ", omp_get_thread_num()
#endif
       !print*,"Calling calc for the cpu on thread ", thread_id
       call this%cpu_nonlinearity%calc(gy_chi(:,:,this%gpu_block+1:this%lbg0),&
            &g_block(:,:,this%gpu_block+1:this%lbg0),&
            &vexy(:,:,:,this%gpu_block+1:this%lbg0),&
            &dgdxy(:,:,:,this%gpu_block+1:this%lbg0),&
            &localrhs(:,:,this%gpu_block+1:this%lbg0),first)
    endif
#ifdef WITHOMP_SPLIT_NONLIN
    !$OMP END PARALLEL
#endif
  end Subroutine calc_split_nonlinearity_df

  ! We need a possibility to get the best gpu_cpu_ratio for
  ! each given problem. Therefore we optimize the nonlinearity
  ! with respect to this parameter.

  ! 1. find the next gpu_cpu_ratio to assess
  ! 2. initialize
  ! 3. time execution
  ! 4. save gpu_cpu_ratio with timing results
  ! 5. finalize
  ! 6. if best found exit else goto 1
  subroutine autotune_split_nonlinearity_df(this,total_blocksize)
    class(df_split_nonlinear_term_t) :: this
    integer :: total_blocksize

    ! Local variables
    Complex, Dimension(:,:,:),allocatable :: gy_chi, g_block,localrhs
    Complex, Dimension(:,:,:,:), allocatable:: vexy, dgdxy
    logical :: first
    double precision :: start_time, end_time, diff_time,best_time, prev_diff_time
    real :: min_gc, max_gc, dgc, temp_gc_ratio, best_gc_ratio
    integer :: nSteps, iStep,max_threads, omp_get_max_threads, repeat
    integer(C_LONG) :: memory_need_on_device, free_memory_on_device
    integer(C_LONG), parameter :: max_device_memory=5636554752

    write(*,"(A)",advance='no') "Running autotune for split_nonlinearity."
    this%lbg0 = total_blocksize

    ! allocate the test_arrays
    !allocate(vexy(li1:li2,lj1:lj2,2,this%lbg0))
    !allocate(dgdxy(li1:li2,lj1:lj2,2,this%lbg0))
    allocate(localrhs(li1:li2,lj1:lj2,this%lbg0))
    allocate(gy_chi(li1:li2,lj1:lj2,this%lbg0))
    allocate(g_block(lbi:ubi,lj1:lj2,this%lbg0))
    first = .true.
    max_threads = omp_get_max_threads()
    call this%gpu_nonlinearity%setCudaDevice(2)

    nSteps = 50
    ! set the start and end for the scan
    ! min_gc is to have the GPU work as much as all the other
    ! threads. In this case, we can switch off the GPU and use 
    ! a normal CPU thread. So that's the absolute minimum.
    min_gc = 1.0/(max_threads-1)
    max_gc = 5.0
    dgc = (max_gc-min_gc)/nSteps
    best_time = huge(best_time)
    iStep = 1
    repeat = 1
    do while (iStep.le.nSteps)
       write(*,"(A)",advance="no") "."
       temp_gc_ratio = min_gc + (iStep-1)*dgc
       call this%set_gpu_cpu_ratio(temp_gc_ratio)
       call this%set_blocksize(total_blocksize)
       call this%initialize()
       memory_need_on_device = this%gpu_nonlinearity%get_memory_need_on_device()
       free_memory_on_device = this%gpu_nonlinearity%get_free_memory_on_device()
       !print*,"memory_need_on_device = ",memory_need_on_device, ", free = ",free_memory_on_device
       if (memory_need_on_device.lt.0.95*free_memory_on_device) then
          call this%allocate_arrays(vexy,dgdxy)
          vexy = cmplx(1.0,1.0)
          dgdxy = cmplx(0.1,0.1)
          
          call get_systime(start_time)
          call this%calc(gy_chi,g_block,vexy,dgdxy,localrhs,first)
          call get_systime(end_time)
          call this%free_arrays(vexy,dgdxy)
          call this%finalize()
          diff_time = end_time - start_time
          if (diff_time<best_time) then
             best_time=diff_time
             best_gc_ratio = temp_gc_ratio
          end if
          if ((iStep.gt.1).and.((diff_time-prev_diff_time)/prev_diff_time.gt.0.05).and.(repeat.le.5)) then
             repeat = repeat + 1
             ! we have a jump, just do the same value again
             continue
          else
             repeat = 1
             !write(*,"(A,F7.3,A,F10.5)")"gpu_cpu_ratio = ", temp_gc_ratio," needs ",diff_time
             iStep = iStep + 1
             prev_diff_time = diff_time
          end if
       else
          print*,"Splitting not tested due to memory_need = ",memory_need_on_device
          call this%finalize()
       end if
    end do

    write(*,"(//,5X,A,F7.4,//)") "Setting gpu_cpu_ratio to ",best_gc_ratio
    call this%set_gpu_cpu_ratio(best_gc_ratio)
    deallocate(localrhs,gy_chi,g_block)

  end subroutine autotune_split_nonlinearity_df
end module df_split_nonlinear_term_mod
