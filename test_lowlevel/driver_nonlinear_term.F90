#include "redef.h"
module df_dummy_nonlinear_mod
  use df_omp_nonlinear_term_mod
  use df_nonlinear_term_mod
  use df_nonlinear_cuda_mod
  use df_nonlinear_prefactor_dummy_mod
  implicit none

!  type, extends(df_nonlinear_term_t) :: df_dummy_nonlinear_t
!  type, extends(df_omp_nonlinear_term_t) :: df_dummy_nonlinear_t
  type, extends(df_nonlinear_cuda_t) :: df_dummy_nonlinear_t
   contains
     procedure :: construct => construct_with_dummy
     procedure :: getType => getThisType
  end type df_dummy_nonlinear_t

contains
  function getThisType(this)
    class(df_dummy_nonlinear_t) :: this
    character(len=MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "df_dummy_nonlinear_term_t"
  end function getThisType
    
  subroutine construct_with_dummy(this,equil_par_curr)
    class(df_dummy_nonlinear_t) :: this
    logical :: equil_par_curr
    
    allocate(df_nonlinear_prefactor_dummy_t::this%prefactor)

  end subroutine construct_with_dummy
end module df_dummy_nonlinear_mod

program driver_nonlinear_term
  use discretization
  use nonlinear_term_mod
  !use df_nonlinear_term_mod
  use df_dummy_nonlinear_mod
  use blockindex
  use par_other, only: n_fields
  use page_locked_memory_mod
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif
  implicit none

  integer :: n_procs
  integer :: lb1,lb2, i_block
  !complex, dimension(li1:li2, lj1:lj2, lb1:lb2), target:: p_rhs    
  complex, dimension(:,:,:), pointer :: rhs_block
  integer :: stage
  !complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
  complex, dimension(:,:,:,:), pointer :: p_emfields
  Complex, Dimension(:,:,:), pointer :: ptr_barchi
  Complex, Dimension(:,:,:), allocatable,target :: p_barchi
  Complex, Dimension(:,:,:,:), allocatable,target :: p_dbarchidxy, p_dgdxy
  complex, dimension(:,:,:,:), pointer :: ptr_dbarchidxy,ptr_dgdxy
  !complex, dimension(lbi:ubi, lj1:lj2, lb1:lb2), intent(in), target:: g_block
  complex, dimension(:,:,:), pointer :: g_block
  complex, dimension(:,:,:,:), allocatable, target :: global_rhs
  integer :: klmn, i_nblocks,nParts
  class(nonlinear_term_t),pointer :: this_nonlinear_term
  double precision :: inc_time, inc_MFlops
#ifdef WITHOMP
  integer :: omp_get_max_threads
#endif

  LIKWID_INIT
  PERFINIT
  PERFON('MAIN')
#ifdef WITHOMP
  write(*,"(A,I2,A)") "Running with ",omp_get_max_threads()," threads."
#endif
  n_procs_x = 1
  n_procs_y = 1
  n_procs_z = 1
  n_procs_v = 1
  n_procs_w = 1
  n_procs_s = 1
  
  nx0 = 64; nxb = 2
  nky0= 16
  nz0 = 24; nzb = 2
  nv0 = 96; nvb = 0
  nw0 = 16; nwb = 0
  n_spec = 1
  dealiasing = .false.
  if (n_spec.gt.1) n_fields=2
  call initialize_discretization(.true.)
  write(*,"(4I5,I10)") li1,li2,lj1,lj2,lklmn0
  turbdeal = .false.
  x_local=.false.
  ! nParts is a parameter from the cuda part. It determines the
  ! number of parts in which the arrays are transferred to the GPU.
  nParts = 1

  open(133,file="timing.dat")
  !do i_nblocks=1,lklmn0/nParts
  do i_nblocks=8,8
     if (modulo(lklmn0,i_nblocks).eq.0) then
        nblocks = i_nblocks
        !print*,"----- Trying with nblocks =", nblocks
     else
        cycle
     end if
     call initialize_blockindex()

     allocate(global_rhs(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks))
     call lock_memory(global_rhs)
     !global_rhs = cmplx(0.0,0.0)
#if 1
     !$OMP PARALLEL default(none) private(i_block,lb1,lb2,rhs_block,klmn) &
     !$OMP shared(nblocks,lbg0,global_rhs)
     do i_block=1,nblocks
        lb1 = (i_block-1)*lbg0
        lb2 = lb1+lbg0-1
        rhs_block => global_rhs(:,:,:,i_block)
        !$OMP DO schedule(static,2)
        do klmn=1,lbg0
           rhs_block(:,:,klmn) = cmplx(0.0,0.0)
        end do
        !$OMP END DO
     end do
     !$OMP END PARALLEL
#endif
     write(*,"(A,ES20.10)") "global_rhs after init: ",sum(abs(global_rhs))
     allocate(df_dummy_nonlinear_t::this_nonlinear_term)
     
     print*,"this_nonlinear_term is of type ",this_nonlinear_term%getType()
#ifdef WITH_CUDA_NONLIN
     select type (this_nonlinear_term)
     class is (df_nonlinear_cuda_t)
        call this_nonlinear_term%SetCudaDevice(0)
     end select
#endif
     call this_nonlinear_term%set_blocksize(lbg0)
     call this_nonlinear_term%construct(.false.)
     call this_nonlinear_term%initialize()

     allocate(p_emfields(li1:li2, lj1:lj2, lbz:ubz, 1:n_fields))
     allocate(g_block(lbi:ubi, lj1:lj2, 1:lbg0))
     allocate(p_barchi(li1:li2,lj1:lj2,1:lbg0))
     ptr_barchi => p_barchi
     call this_nonlinear_term%allocate_arrays(p_dbarchidxy,p_dgdxy)
     p_dbarchidxy(:,:,1,:) = cmplx(1.0,0.0)
     p_dbarchidxy(:,:,2,:) = cmplx(2.0,0.0)
     p_dgdxy(:,:,1,:) = cmplx(3.0,0.0)
     p_dgdxy(:,:,2,:) = cmplx(4.0,0.0)
     stage=2
     ptr_dbarchidxy => p_dbarchidxy
     ptr_dgdxy => p_dgdxy
     PERFON('nl_add')
     do i_block=1,nblocks
        lb1 = (i_block-1)*lbg0
        lb2 = lb1+lbg0-1
        rhs_block => global_rhs(:,:,:,i_block)
        !write(*,"(I3,A,ES20.10,3I7)") i_block,": rhs_block before add_nl: ",&
        !     &sum(abs(rhs_block)),size(rhs_block,1),size(rhs_block,2),size(rhs_block,3)
     
        !if (associated(g_block)) print*,"g_block is associated"
        call this_nonlinear_term%add(g_block,ptr_dgdxy,p_emfields,ptr_barchi,ptr_dbarchidxy,rhs_block,lb1,lb2,stage)
        !write(*,"(I3,A,ES20.10)") i_block,": rhs_block = ",sum(abs(rhs_block))
     end do
     PERFOFF
     write(*,"(A,ES20.10)") "global_rhs = ",sum(abs(global_rhs))

     call this_nonlinear_term%free_arrays(p_dbarchidxy,p_dgdxy)
     nullify(ptr_dbarchidxy)
     nullify(ptr_dgdxy)
     call this_nonlinear_term%finalize()
     call this_nonlinear_term%destruct()
     deallocate(this_nonlinear_term)
     deallocate(global_rhs)
     deallocate(g_block)
     deallocate(p_barchi)
     deallocate(p_emfields)

     PERF_GET('nl_add',inc_time,inc_MFlops)
     PERF_RESET('nl_add')
#ifdef WITHPERF
     write(*,"(A,I5,A,ES10.3)") "nblocks = ",nblocks,", time = ",inc_time
     write(133,"(2I6,ES10.3)") nblocks,lbg0,inc_time
#endif
  end do
  PERFOFF
  LIKWID_CLOSE
  PERFOUT('MAIN')
  close(133)
end program driver_nonlinear_term
  
