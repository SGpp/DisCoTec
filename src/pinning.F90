module pinning_mod
#if defined(__INTEL_COMPILER) && defined(WITHOMP)
  use omp_lib
  implicit none

contains
  subroutine pinning
    integer(kind=kmp_affinity_mask_kind) :: mask
    integer :: max_proc, retval, isset, iproc,ithread

    !$OMP PARALLEL private(max_proc,mask,retval,isset,iproc,ithread)
    !$OMP DO ORDERED
    do ithread=0,omp_get_num_threads()-1
       call kmp_create_affinity_mask(mask)
       max_proc = kmp_get_affinity_max_proc()
       !ncores = max_proc/2
       retval = kmp_get_affinity(mask)
       if (retval.eq.0) then
          !$OMP ORDERED
          write(*,"(I3,A)",advance='no') omp_get_thread_num(),": mask = "
          do iproc=0,max_proc-1
             isset = kmp_get_affinity_mask_proc(iproc,mask)
             write(*,"(I1)",advance='no') isset
          end do
          write(*,"(A)") ""
          !$OMP END ORDERED
       end if
       call kmp_destroy_affinity_mask(mask)
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine pinning
#endif
end module pinning_mod
