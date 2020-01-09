#include "redef.h"

module df_nonlinear_cuda_mod
  use,intrinsic :: iso_c_binding
  use discretization, only: pi1,pi2, li0,li1,li2,lbi,ubi, lj0,lj1,lj2,mype
  use communications, only: t_n_procs
  use df_nonlinear_term_mod, only:df_nonlinear_term_t
  implicit none

  type,public,extends(df_nonlinear_term_t) :: df_nonlinear_cuda_t
     integer :: cuda_device
     integer :: n_ranks_per_device
   contains
     procedure :: initialize => initialize_nonlinear_cuda
     procedure :: finalize => finalize_nonlinear_cuda
     procedure :: calc => calc_nonlinear_cuda
     procedure :: allocate_arrays => allocate_arrays_cuda
     procedure :: free_arrays => free_arrays_cuda
     procedure :: register_array
     procedure :: unregister_array
     procedure :: allocate_arrays_on_device
     procedure :: free_arrays_on_device
     procedure :: getType => getThisType
     procedure :: setCudaDevice
     procedure :: get_nearest_blocksize
     procedure :: get_memory_need_on_device
     procedure :: get_free_memory_on_device
  end type df_nonlinear_cuda_t

  interface
     subroutine cuda_initialize_nonlinearity_df(a_lbg0,cptr_pnl_1d) bind(c)
       import
       integer(C_INT),value :: a_lbg0
       real(C_DOUBLE),dimension(pi1:pi2),intent(IN) :: cptr_pnl_1d
     end subroutine cuda_initialize_nonlinearity_df

     subroutine cuda_finalize_nonlinearity_df() bind(c)
     end subroutine cuda_finalize_nonlinearity_df

     subroutine cuda_set_device(device) bind(c)
       import
       integer(C_INT),value :: device
     end subroutine cuda_set_device

     subroutine cuda_register_array(arr,arr_size) bind(c)
       import
       integer(C_INT), value :: arr_size
       Complex(C_DOUBLE_COMPLEX), Dimension(arr_size),Intent(inout) :: arr
     end subroutine cuda_register_array

     subroutine cuda_unregister_array(arr) bind(c)
       import
       Complex(C_DOUBLE_COMPLEX), Dimension(*),Intent(inout) :: arr
     end subroutine cuda_unregister_array

     integer(C_INT) function cuda_get_device_count() bind(c)
       import
     end function cuda_get_device_count

     integer(C_LONG) function cuda_get_memory_need_on_device() bind(c)
       import
     end function cuda_get_memory_need_on_device

     integer(C_LONG) function cuda_get_free_memory_on_device() bind(c)
       import
     end function cuda_get_free_memory_on_device

     subroutine cuda_allocate_on_device(device) bind(c)
       import
       integer(C_INT), value :: device
     end subroutine cuda_allocate_on_device

     subroutine cuda_free_on_device() bind(c)
       import
     end subroutine cuda_free_on_device

     integer(C_INT) function cuda_get_nearest_blocksize(test_blocksize) bind(c)
       import
       integer(C_INT), value :: test_blocksize
     end function cuda_get_nearest_blocksize

     !subroutine initialize_fourier_cufft(a_li0da,a_ly0da) bind(c)
     !  import
     !  integer(C_INT), value :: a_li0da,a_ly0da
     !end subroutine initialize_fourier_cufft
     
     !subroutine finalize_fourier_cufft() bind(c)
     !end subroutine finalize_fourier_cufft

     !subroutine cuda_calc_nonlinearity_df(gy_chi,g2d,vexy,dgdxy,localrhs,cptr_pnl_1d,first) bind(c)
     subroutine cuda_calc_nonlinearity_df(gy_chi,g_block,vexy,dgdxy,localrhs,first) bind(c)
       import
       Complex(C_DOUBLE_COMPLEX), Dimension(li1:li2,lj1:lj2,1:*),Intent(in) :: gy_chi
       Complex(C_DOUBLE_COMPLEX), Dimension(li1:li2,lj1:lj2,1:*),Intent(in) :: g_block
       Complex(C_DOUBLE_COMPLEX), Dimension(li1:li2,lj1:lj2,2,1:*),Intent(in) :: vexy, dgdxy
       Complex(C_DOUBLE_COMPLEX), Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
       !real(C_DOUBLE),dimension(pi1:pi2),intent(IN) :: cptr_pnl_1d
       Logical(C_BOOL), value, Intent(in):: first
     end subroutine cuda_calc_nonlinearity_df
  end interface

  private

contains

  function getThisType(this)
    class(df_nonlinear_cuda_t) :: this
    character(len=MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "df_nonlinear_cuda_t"
  end function getThisType

  subroutine SetCudaDevice(this,device)
    class(df_nonlinear_cuda_t) :: this
    integer :: device

    integer(C_INT) :: n_devices

    n_devices = cuda_get_device_count()
    if (device.lt.0) then
       !n_devices = 1
       this%n_ranks_per_device = t_n_procs/n_devices
       if (modulo(t_n_procs,n_devices) .ne.0) then
          print*,"We have ",t_n_procs," ranks and ",n_devices," devices."
       end if
       
       if (this%n_ranks_per_device.eq.0) then
          this%n_ranks_per_device=1
       end if
       this%cuda_device = mype/this%n_ranks_per_device
    else
       this%cuda_device = device
       this%n_ranks_per_device = t_n_procs
    end if

    if (mype.eq.0) then
       write(*,"(A,I2,A)") "CUDA: ",this%n_ranks_per_device," MPI ranks are using one GPU together."
       write(*,"(A,I1,A)") "CUDA: We have in total ",n_devices," devices per node."
    end if

    !this%cuda_device = device

    call cuda_set_device(this%cuda_device)
    
  end subroutine SetCudaDevice

  integer function get_nearest_blocksize(this,test_blocksize)
    class(df_nonlinear_cuda_t) :: this
    integer :: test_blocksize

    get_nearest_blocksize = cuda_get_nearest_blocksize(test_blocksize)
  end function get_nearest_blocksize

  subroutine initialize_nonlinear_cuda(this)
    class(df_nonlinear_cuda_t) :: this

    integer(C_LONG) :: mem_need, mem_free
    if (this%init_status.eq.1) return
    !write(*,"(3A,I5)") "Initializing ",this%getType()," with blocksize = ",this%lbg0
    ! the following subroutine does only set some values, but does not allocate anything
    ! therefore we do not need to call a corresponding finalize routine. This is an
    ! exception of the usual rule!
    
    call this%prefactor%initialize()
    call cuda_initialize_nonlinearity_df(this%lbg0,this%prefactor%pnl_1d)
    mem_need = get_memory_need_on_device(this)
    mem_free = get_free_memory_on_device(this)
    if (mype.eq.0) then
       if (real(mem_need)*this%n_ranks_per_device .gt.real(mem_free)) then
          write(*,"(A,2(I12,A))") "Not enough device memory. We need ",&
               &mem_need*this%n_ranks_per_device," bytes, but have only ",mem_free," bytes."
          stop "Not enough device memory"
       else
          write(*,"(A,2(I12,A))") "We need ",&
               &mem_need*this%n_ranks_per_device," bytes, and have ",mem_free," bytes."
       end if
    end if
    
    this%init_status = 1
  end subroutine initialize_nonlinear_cuda

  function get_memory_need_on_device(this) result(mem_need)
    class(df_nonlinear_cuda_t) :: this
    integer(C_LONG) :: mem_need
    mem_need = cuda_get_memory_need_on_device()
  end function get_memory_need_on_device

  integer(C_LONG) function get_free_memory_on_device(this)
    class(df_nonlinear_cuda_t) :: this
    get_free_memory_on_device = cuda_get_free_memory_on_device()
  end function get_free_memory_on_device

  subroutine allocate_arrays_cuda(this,vexy,dgdxy)
    class(df_nonlinear_cuda_t) :: this
    complex, dimension(:,:,:,:),allocatable :: vexy,dgdxy
    integer :: device

    !print*,"Calling cuda_allocate_on_device with vexy and dgdxy"
    allocate(vexy(li1:li2,lj1:lj2,2,1:this%lbg0))
    allocate(dgdxy(li1:li2,lj1:lj2,2,1:this%lbg0))

    call this%register_array(vexy,li0*lj0*2*this%lbg0)
    call this%register_array(dgdxy,li0*lj0*2*this%lbg0)

    call this%allocate_arrays_on_device()

  end subroutine allocate_arrays_cuda

  subroutine free_arrays_cuda(this,vexy,dgdxy)
    class(df_nonlinear_cuda_t) :: this
    complex, dimension(:,:,:,:),allocatable :: vexy,dgdxy

    call this%free_arrays_on_device()

    call this%unregister_array(vexy)
    call this%unregister_array(dgdxy)
    deallocate(vexy)
    deallocate(dgdxy)
  end subroutine free_arrays_cuda

  !******* Routines for page-locking and unlocking an array
  subroutine register_array(this,arr,arr_size)
    class(df_nonlinear_cuda_t) :: this
    complex, dimension(:,:,:,:) :: arr
    integer :: arr_size

    call cuda_register_array(arr,arr_size)
  end subroutine register_array

  subroutine unregister_array(this,arr)
    class(df_nonlinear_cuda_t) :: this
    complex, dimension(:,:,:,:) :: arr

    call cuda_unregister_array(arr)
  end subroutine unregister_array

  !***** Allocation and free of the device arrays
  subroutine allocate_arrays_on_device(this)
    class(df_nonlinear_cuda_t) :: this

    call cuda_allocate_on_device(this%cuda_device)

  end subroutine allocate_arrays_on_device

  subroutine free_arrays_on_device(this)
    class(df_nonlinear_cuda_t) :: this

    call cuda_free_on_device()
  end subroutine free_arrays_on_device

  subroutine finalize_nonlinear_cuda(this)
    class(df_nonlinear_cuda_t) :: this

    call cuda_finalize_nonlinearity_df()
    call this%prefactor%finalize()
    !call c_finalize_nonlinearity_df()
    this%init_status = 0
  end subroutine finalize_nonlinear_cuda

  subroutine calc_nonlinear_cuda(this,gy_chi,g_block,vexy,dgdxy,localrhs,first)
    class(df_nonlinear_cuda_t) :: this
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: gy_chi
    complex, dimension(lbi:ubi,lj1:lj2,1:*),intent(in) :: g_block
    Complex, Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
    Logical, Intent(in) :: first

    !print*,"prefactor status : ",allocated(this%prefactor),allocated(this%prefactor%pnl_1d)
    !call cuda_calc_nonlinearity_df(gy_chi,&
    !     &g2d,vexy,dgdxy,localrhs,this%prefactor%pnl_1d,logical(first,kind=C_BOOL))
    call cuda_calc_nonlinearity_df(gy_chi,&
         &g_block,vexy,dgdxy,localrhs,logical(first,kind=C_BOOL))
  end subroutine calc_nonlinear_cuda

end module df_nonlinear_cuda_mod
