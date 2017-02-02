#include "redef.h"
#include "intrinsic_sizes.h"

module df_mic_nonlinear_term_mod
  use mic_lib
  use omp_lib
  use df_nonlinear_term_mod

  use par_mod, only: imag, pi, ve_max, phase_1, phase_2, spec
  use par_other, only: n_fields
  use communications, only: MPI_COMPLEX_TYPE, MPI_REAL_TYPE, mpi_comm_y,&
       &my_2reals_max_to_all, my_barrier, mpi_comm_xy, reduce_per_thread,&
       &threadlocal_mpi_comm_y
  use discretization, only: li0,li1,li2, nky0, lj0,lj1,lj2,nj0, n_procs_y
  use fourier_mic
  use mpi
  use,intrinsic :: iso_c_binding
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

  implicit none

  type, public, extends(df_nonlinear_term_t) :: df_mic_nonlinear_term_t
   contains
     procedure :: initialize => initialize_mic_nonlinearity
     procedure :: finalize => finalize_mic_nonlinearity
     procedure :: calc => calc_mic_nonlinearity_df
     procedure :: getType => getThisType
  end type df_mic_nonlinear_term_t
  private
#if defined(WITH_OFFLOAD) && defined(WITHPERF)
  !DEC$ attributes offload:mic :: perfon,perfoff,perfinit,perfout
#endif

#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: klmn_per_thread
#endif
  integer :: klmn_per_thread

  interface
     subroutine nonlin_kernel(vexy_re,dgdxy_re,nonlin_re,length) bind(c)
       use,intrinsic :: iso_c_binding
       type(C_PTR),value :: vexy_re
       type(C_PTR),value :: dgdxy_re
       integer(C_LONG),value :: length
       type(C_PTR),value :: nonlin_re
     end subroutine nonlin_kernel
  end interface
contains

  function getThisType(this)
    class(df_mic_nonlinear_term_t) :: this
    character(len=MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "df_mic_nonlinear_term_t"
  end function getThisType
    
  subroutine initialize_mic_nonlinearity(this)
    class(df_mic_nonlinear_term_t) :: this

    integer :: local_lbg0

    if (this%init_status==1) return
    write(*,"(3A,I5)") "Initializing ",this%getType()," with blocksize = ",this%lbg0
    local_lbg0 = this%lbg0
#if defined(WITH_OFFLOAD)
    !DEC$ offload begin target(mic:WHICHMIC) in(li0da,ly0da,nj0,lj0,n_procs_y,local_lbg0) &
    !DEC$ out(klmn_per_thread)
    !print*,"Before anything on device."
    !flush(6)
    PERFON('init')
#endif
    !print*,"Initializing fourier"
    !flush(6)
    call initialize_fourier(li0da,ly0da,local_lbg0,klmn_per_thread)
    !print*,"After initializing fourier"
    !flush(6)
#if defined(WITH_OFFLOAD)
    PERFOFF
    !DEC$ end offload
#endif
    !print*,"Now initializing the prefactor..."
    call this%prefactor%initialize()
    !print*,"Done"

    this%init_status = 1
  end subroutine initialize_mic_nonlinearity

  subroutine finalize_mic_nonlinearity(this)
    class(df_mic_nonlinear_term_t) :: this

    print*,"Finalizing mic_nonlinearity"
    call this%prefactor%finalize()

#ifdef WITH_OFFLOAD
    !DEC$ offload begin target(mic:WHICHMIC)
    !deallocate(dev_vexy)
#endif
    print*,"Finalizing Fourier on device."
    flush(6)
    call finalize_fourier
    print*,"After Finalizing Fourier on device."
    flush(6)
#if defined(WITH_OFFLOAD)
    !DEC$ end offload
#endif

    this%init_status = 0
  end subroutine finalize_mic_nonlinearity

#undef DEBUGGING

#define MIC_ALL_IN_ONE_LOOP

#ifdef MIC_ALL_IN_ONE_LOOP
  !>Computes the nonlinearity for the (x) nonlocal version
  !!\todo Get rid of transposition using striding? make some arrays for arakawa allocatable
  !!\todo Check whether the nonlinearity prefactor is important for the CFL criterion
  Subroutine calc_mic_nonlinearity_df(this,gy_chi,g2d,vexy,dgdxy,localrhs,first)
    class(df_mic_nonlinear_term_t) :: this
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: gy_chi, g2d
    Complex, Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
    Logical, Intent(in) :: first

    ! Local variables
    Complex, Dimension(li1:li2,lj1:lj2,1:this%lbg0) :: nonlin    
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,2) ::  vexy_re, dgdxy_re
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1) ::  nonlin1_re

    Real, Dimension(1:2)  ::  vexy_max_loc
    Integer :: klmn, local_lbg0, direction, local_length
    integer(KIND=8) :: nbytes

    !print*,"vexy = ",sum(abs(vexy(:,:,:,1:this%lbg0))),", dgdxy = ",sum(abs(dgdxy(:,:,:,1:this%lbg0)))
    local_lbg0 = this%lbg0
    local_length = li0*lj0*this%lbg0

#ifdef WITH_OFFLOAD
    nbytes = 5*local_length*SIZE_OF_COMPLEX !vexy,dgdxy,nonlin
    ! tmp_arr
    nbytes = nbytes + nj0*li0da*SIZE_OF_COMPLEX
    ! vexy_re, dgdxy_re, nonlin1_re
    nbytes = nbytes + ly0da*li0da*5*SIZE_OF_REAL
    !print*,"Trying to allocate on device, nbytes = ",nbytes
    PERFON('mic_allc')
    !DEC$ offload_transfer target(mic:WHICHMIC) &
    !DEC$ in(local_length,li0da) &
    !DEC$ nocopy(vexy:length(2*local_length) alloc_if(.true.) align(64)) &
    !DEC$ nocopy(dgdxy:length(2*local_length) alloc_if(.true.) align(64)) &
    !DEC$ nocopy(nonlin: alloc_if(.true.) align(64)) &
    !DEC$ nocopy(vexy_re: alloc_if(.true.) align(64)) &
    !DEC$ nocopy(dgdxy_re: alloc_if(.true.) align(64)) &
    !DEC$ nocopy(nonlin1_re: alloc_if(.true.) align(64)) &
    !DEC$ nocopy(tmp_arr: alloc_if(.true.) align(64))
    PERFOFF
    !print*,"After alloc"

    PERFON('trf_in')
    !DEC$ offload_transfer target(mic:WHICHMIC) &
    !DEC$ in(local_length) &
    !DEC$ in(vexy:length(2*local_length) alloc_if(.false.) free_if(.false.)) &
    !DEC$ in(dgdxy:length(2*local_length) alloc_if(.false.) free_if(.false.))
    PERFOFF
    !print*,"After transfer in"
#endif

    PERFON('mic_comp')

#ifdef WITH_OFFLOAD
    !DEC$ omp offload target(mic:WHICHMIC) &
    !DEC$ in(local_length) &
    !DEC$ nocopy(vexy:length(2*local_length) alloc_if(.false.) free_if(.true.)) &
    !DEC$ nocopy(dgdxy:length(2*local_length) alloc_if(.false.) free_if(.true.)) &
    !DEC$ in(nky0,n_procs_y, first, local_lbg0) &
    !DEC$ nocopy(nonlin) &
    !DEC$ nocopy(vexy_re,dgdxy_re,nonlin1_re,tmp_arr) &
    !DEC$ out(ve_x_max_loc,ve_y_max_loc)
#endif
    !$OMP PARALLEL default(none) &
    !$OMP shared(vexy,dgdxy,nonlin) &
    !$OMP shared(first) &
    !$OMP shared(local_lbg0,li0da,ly0da,n_procs_y) &
    !$OMP reduction(max:ve_x_max_loc,ve_y_max_loc) &
    !$OMP private(vexy_re,dgdxy_re) &
    !$OMP private(tmp_arr,nonlin1_re) &
    !$OMP private(direction,klmn)

    LIKWID_ON('mic_comp')
#ifdef WITH_OFFLOAD
    PERFON('mic_comp')
#endif
    !$OMP DO
    do klmn=1,local_lbg0

       do direction=1,2
          
          PERFON('transp')
          ! Transpose x-y
          tmp_arr = transpose(vexy(:,:,direction,klmn))
          PERFOFF

          ! Fourier transform in y (include dealiasing step)
          Call to_real_y(tmp_arr,vexy_re(:,:,direction))
       end do

       do direction=1,2
          
          PERFON('transp')
          ! Transpose x-y
          tmp_arr = transpose(dgdxy(:,:,direction,klmn))
          PERFOFF
             
          ! Fourier transfrom in y (include dealiasing step)
          Call to_real_y(tmp_arr,dgdxy_re(:,:,direction))
       end do

       !write(*,"(2(A,ES13.6))") "AFTER FFT: vexy_re = ",&
       !     &sum(vexy_re**2),", dgdxy_re = ",sum(dgdxy_re**2)
       
       ! get max ExB velocity
       if (first) then
          call this%prefactor%multiply_max(vexy_max_loc,vexy_re,this%lbg0)
          ve_max(1)=max(ve_max(1),vexy_max_loc(1)) !ve_x
          ve_max(2)=max(ve_max(2),vexy_max_loc(2)) !ve_y
       end if

       !compute the 'standard' nonlinear term
       PERFON('nonlin')
       LIKWID_ON('nonlin')
       !print*,omp_get_thread_num(),klmn,": Calling nonlin_kernel"
       nonlin1_re = -vexy_re(:,:,1)*dgdxy_re(:,:,2) &
            & + vexy_re(:,:,2)*dgdxy_re(:,:,1)
       LIKWID_OFF('nonlin')
       PERFOFF
       ! transform back to fourier space
       LIKWID_ON('toF_y')
       Call to_fourier_y(nonlin1_re,tmp_arr)
       LIKWID_OFF('toF_y')
       LIKWID_ON('transp2')
       PERFON('transp2')
       nonlin(:,:,klmn) = transpose(tmp_arr)
       PERFOFF
       LIKWID_OFF('transp2')
    end do
    !$OMP END DO
#ifdef WITH_OFFLOAD
    PERFOFF
#endif
    LIKWID_OFF('mic_comp')
    !$OMP END PARALLEL
    ! --- end of offload ---
    PERFOFF

#ifdef WITH_OFFLOAD
    PERFON('trf_out')
    !DEC$ offload_transfer target(mic:WHICHMIC) &
    !DEC$ out(nonlin: alloc_if(.false.)) &
    !DEC$ nocopy(vexy_re: free_if(.true.)) &
    !DEC$ nocopy(dgdxy_re: free_if(.true.)) &
    !DEC$ nocopy(nonlin1_re: free_if(.true.)) &
    !DEC$ nocopy(tmp_arr: free_if(.true.))
    PERFOFF
#endif
    !write(*,"(A,ES17.10)") "nonlin = ", sum(abs(nonlin))

    ! multiply with the prefactor and write to localrhs
    !$OMP PARALLEL default(none) shared(nonlin,localrhs,this)
    call this%prefactor%multiply_with(nonlin,localrhs,this%lbg0)
    !$OMP END PARALLEL

    if (first) then
       ve_max(1)=max(ve_max(1),ve_x_max_loc)        
       ve_max(2)=max(ve_max(2),ve_y_max_loc)
    end if

  End Subroutine calc_mic_nonlinearity_df

#else
  !>Computes the nonlinearity for the (x) nonlocal version
  !!\todo Get rid of transposition using striding? make some arrays for arakawa allocatable
  !!\todo Check whether the nonlinearity prefactor is important for the CFL criterion
  Subroutine calc_mic_nonlinearity_df(this,gy_chi,g2d,vexy,dgdxy,localrhs,first)
    class(df_mic_nonlinear_term_t) :: this
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: gy_chi, g2d
    Complex, Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
    Logical, Intent(in) :: first

    ! Local variables
    Complex, Dimension(li1:li2,lj1:lj2,1:this%lbg0) :: nonlin    
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr
    !Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1,2,klmn_per_thread) :: transposed_block
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1,2,1:this%lbg0) :: transposed_block
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,2,1:this%lbg0) ::  vexy_re, dgdxy_re
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1) ::  nonlin1_re
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr2

    Real:: ve_x_max_loc, ve_y_max_loc
    Integer :: klmn, local_lbg0, direction, local_length, n_threads,iThread,local_klmn
    integer(KIND=8) :: nbytes

    !print*,"vexy = ",sum(abs(vexy(:,:,:,1:this%lbg0))),", dgdxy = ",sum(abs(dgdxy(:,:,:,1:this%lbg0)))
    local_lbg0 = this%lbg0
    local_length = li0*lj0*this%lbg0

#ifdef WITH_OFFLOAD
    nbytes = 5*local_length*SIZE_OF_COMPLEX !vexy,dgdxy,nonlin
    ! tmp_arr, tmp_arr2
    nbytes = nbytes + 2*nj0*li0da*SIZE_OF_COMPLEX
    ! transposed_block
    nbytes = nbytes + nj0*li0da*2*this%lbg0*SIZE_OF_COMPLEX

    !nbytes = nbytes + (2+2*klmn_per_thread)*nj0*li0da*SIZE_OF_COMPLEX
    nbytes = nbytes + ly0da*li0da*(1+4*this%lbg0)*SIZE_OF_REAL
    !print*,"Trying to allocate on device, nbytes = ",nbytes
    PERFON('mic_allc')
    !DEC$ omp offload target(mic:WHICHMIC) &
    !DEC$ in(local_length,li0da) &
    !DEC$ nocopy(vexy:length(2*local_length) alloc_if(.true.) align(64)) &
    !DEC$ nocopy(dgdxy:length(2*local_length) alloc_if(.true.) align(64)) &
    !DEC$ nocopy(nonlin: alloc_if(.true.) align(64)) &
    !DEC$ nocopy(tmp_arr2: alloc_if(.true.) align(64)) &
    !DEC$ nocopy(vexy_re: alloc_if(.true.) align(64)) &
    !DEC$ nocopy(dgdxy_re: alloc_if(.true.) align(64)) &
    !DEC$ nocopy(nonlin1_re: alloc_if(.true.) align(64)) &
    !DEC$ nocopy(transposed_block: alloc_if(.true.) align(64)) &
    !DEC$ nocopy(tmp_arr: alloc_if(.true.) align(64))
    !$OMP PARALLEL 
#if 0
    !$OMP MASTER
    print*,"Zeroing transposed_block for first-touch allocation."
    flush(6)
    !$OMP END MASTER
#endif
    !$OMP WORKSHARE
    transposed_block = cmplx(0.0,0.0)
    !$OMP END WORKSHARE
    !$OMP END PARALLEL 
    PERFOFF
    !print*,"After alloc"
    PERFON('trf_in')
    !DEC$ offload_transfer target(mic:WHICHMIC) &
    !DEC$ in(local_length) &
    !DEC$ in(vexy:length(2*local_length) alloc_if(.false.) free_if(.false.)) &
    !DEC$ in(dgdxy:length(2*local_length) alloc_if(.false.) free_if(.false.))
    PERFOFF
    !print*,"After transfer in"
#endif
    PERFON('mic_comp')
#ifdef WITH_OFFLOAD
    !DEC$ omp offload target(mic:WHICHMIC) &
    !DEC$ in(local_length) &
    !DEC$ nocopy(vexy:length(2*local_length) alloc_if(.false.) free_if(.true.)) &
    !DEC$ nocopy(dgdxy:length(2*local_length) alloc_if(.false.) free_if(.true.)) &
    !DEC$ in(nky0,n_procs_y, first, local_lbg0,klmn_per_thread) &
    !DEC$ nocopy(n_threads)&
    !DEC$ nocopy(nonlin) &
    !DEC$ nocopy(tmp_arr2,vexy_re,dgdxy_re,nonlin1_re,tmp_arr) &
    !DEC$ nocopy(transposed_block) &
    !DEC$ out(ve_x_max_loc,ve_y_max_loc)
#endif
    !$OMP PARALLEL default(none) &
    !$OMP shared(vexy,dgdxy,vexy_re,dgdxy_re,nonlin,first) &
    !$OMP shared(ve_x_max_loc,ve_y_max_loc,n_threads) &
    !$OMP shared(local_lbg0,klmn_per_thread,li0da,ly0da,n_procs_y) &
    !$OMP private(tmp_arr2,tmp_arr,nonlin1_re) &
#ifdef EN_BLOC
    !$OMP shared(transposed_block) &
#endif
    !$OMP private(direction,klmn,local_klmn,iThread)

    LIKWID_ON('mic_comp')

    !$OMP MASTER
    n_threads = omp_get_num_threads()
    !print*,"dev: vexy = ",sum(abs(vexy(:,:,:,1:local_lbg0))),", dgdxy = ",sum(abs(dgdxy(:,:,:,1:local_lbg0)))
    !print*,"n_threads = ",n_threads
    !flush(6)
    !$OMP END MASTER
    !$OMP BARRIER
#ifdef EN_BLOC
    !$OMP DO schedule(static)
    do iThread=0,n_threads-1
       !print*,"Iteration no. ",iThread," on thread no. ",omp_get_thread_num()
       PERFON('tr_bl')
       do local_klmn=1,klmn_per_thread
          klmn=iThread*klmn_per_thread+local_klmn
          !print*,iThread,": klmn = ",klmn
          do direction=1,2
             !transposed_block(:,:,direction,local_klmn) = transpose(vexy(:,:,direction,klmn))
             transposed_block(:,:,direction,klmn) = transpose(vexy(:,:,direction,klmn))
          end do
       end do
       PERFOFF
       !print*,iThread,omp_get_thread_num(),": after transpose ",sum(abs(transposed_block))
       !call to_real_y_block(transposed_block,&
       call to_real_y_block(transposed_block(:,:,:,iThread*klmn_per_thread+1:(iThread+1)*klmn_per_thread),&
            &vexy_re(:,:,:,iThread*klmn_per_thread+1:(iThread+1)*klmn_per_thread),&
            &2*klmn_per_thread*li0da/n_procs_y)
    end do
    !$OMP END DO
#else
    !$OMP DO
    do klmn=1,local_lbg0
       do direction=1,2
          
          PERFON('transp')
          ! Transpose x-y
          tmp_arr2 = transpose(vexy(:,:,direction,klmn))
          PERFOFF

          ! Fourier transform in y (include dealiasing step)
          Call to_real_y(tmp_arr2,vexy_re(:,:,direction,klmn))
       end do
    end do
    !$OMP END DO
#endif

    !$OMP DO 
    do klmn=1,local_lbg0
       do direction=1,2
          
          PERFON('transp')
          ! Transpose x-y
          tmp_arr2 = transpose(dgdxy(:,:,direction,klmn))
          PERFOFF
             
          ! Fourier transfrom in y (include dealiasing step)
          Call to_real_y(tmp_arr2,dgdxy_re(:,:,direction,klmn))
       end do
    end do
    !$OMP END DO

    !!$OMP BARRIER
    !$OMP MASTER
    !write(*,"(2(A,ES13.6))") "AFTER FFT: vexy_re = ",&
    !     &sum(vexy_re**2),", dgdxy_re = ",sum(dgdxy_re**2)

    ! get max ExB velocity
    if (first) then
       ve_x_max_loc=maxval(vexy_re(:,:,2,:))
       ve_y_max_loc=maxval(vexy_re(:,:,1,:))
    end if
    !$OMP END MASTER

    !compute the 'standard' nonlinear term
    !$OMP DO
    do klmn=1,local_lbg0
       PERFON('nonlin')
       LIKWID_ON('nonlin')
       !print*,omp_get_thread_num(),klmn,": Calling nonlin_kernel"
#ifdef WITH_OFFLOAD
       nonlin1_re = -vexy_re(:,:,1,klmn)*dgdxy_re(:,:,2,klmn) &
            & + vexy_re(:,:,2,klmn)*dgdxy_re(:,:,1,klmn)
#else
       call nonlin_kernel(C_LOC(vexy_re(0,0,1,klmn)),&
            &C_LOC(dgdxy_re(0,0,1,klmn)),&
            &C_LOC(nonlin1_re(0,0)),&
            &int(ly0da*li0da,C_LONG))
#endif
       LIKWID_OFF('nonlin')
       PERFOFF
       ! transform back to fourier space
       LIKWID_ON('toF_y')
       Call to_fourier_y(nonlin1_re,tmp_arr)
       LIKWID_OFF('toF_y')
       LIKWID_ON('transp2')
       nonlin(:,:,klmn) = transpose(tmp_arr)
       LIKWID_OFF('transp2')
    end do
    !$OMP END DO
    LIKWID_OFF('mic_comp')
    !$OMP END PARALLEL
    ! --- end of offload ---
    PERFOFF

#ifdef WITH_OFFLOAD
    PERFON('trf_out')
    !DEC$ offload_transfer target(mic:WHICHMIC) &
    !DEC$ out(nonlin: alloc_if(.false.)) &
    !DEC$ nocopy(tmp_arr2: free_if(.true.)) &
    !DEC$ nocopy(vexy_re: free_if(.true.)) &
    !DEC$ nocopy(dgdxy_re: free_if(.true.)) &
    !DEC$ nocopy(nonlin1_re: free_if(.true.)) &
    !DEC$ nocopy(tmp_arr: free_if(.true.))
    PERFOFF
#endif
    !write(*,"(A,ES17.10)") "nonlin = ", sum(abs(nonlin))

    ! multiply with the prefactor and write to localrhs
    call this%prefactor%multiply_with(nonlin,localrhs,this%lbg0)

    if (first) then
       ve_max(1)=max(ve_max(1),ve_x_max_loc)        
       ve_max(2)=max(ve_max(2),ve_y_max_loc)
    end if

  End Subroutine calc_mic_nonlinearity_df
#endif

#if 0
  !> Go to real space
  subroutine gto_real(this,inarr,outarr,howmany)
    class(df_nonlinear_term_t) :: this
    integer, intent(IN) :: howmany
    Complex, Dimension(li1:li2,lj1:lj2,1:howmany),Intent(in) :: inarr
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,1:howmany), target, Intent(out) :: outarr

    ! Local variables
    !Complex, Dimension(li1da:li2da,lj1:lj2,1:howmany) :: tmp_arr1
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr2
    real, dimension(:,:), pointer :: p_out
    !real :: local_sum,thread_global_sum(0:7,0:n_procs_x*n_procs_y)
    Integer:: klmn, i_block,number_of_lbg0_blocks

    number_of_lbg0_blocks = howmany/this%lbg0
    if (number_of_lbg0_blocks.eq.0) then
       !if (dealiasing) then
       !   ! Interpolation for x anti-aliasing 
       !   CALL interpolate_to_fine_grid_x_orig(inarr(:,:,1),tmp_arr1(:,:,1))
       !else
       !   tmp_arr1=inarr
       !end if
       !Call transpose_cmplx(li0da, nj0, tmp_arr1(:,:,1), 0,tmp_arr2 , 0)
       ! Transpose x-y
       Call transpose_cmplx(li0da, nj0, inarr(:,:,1), 0,tmp_arr2 , 0)
       ! Fourier transfrom in y (include dealiasing step)
       Call to_real_y(tmp_arr2,outarr(:,:,1))
    else
       do i_block=1,number_of_lbg0_blocks
          do klmn=1,this%lbg0
             !tmp_arr1(:,:,(i_block-1)*this%lbg0+klmn)=inarr(:,:,(i_block-1)*this%lbg0+klmn)
             ! Transpose x-y
             Call transpose_cmplx(li0da, nj0, inarr(:,:,(i_block-1)*this%lbg0+klmn), 0,tmp_arr2 , 0)
             
             ! Fourier transfrom in y (include dealiasing step)
             p_out => outarr(:,:,(i_block-1)*this%lbg0+klmn)
             Call to_real_y(tmp_arr2,p_out) !outarr(:,:,(i_block-1)*this%lbg0+klmn))
          end do
       end do
    end if

  end subroutine gto_real
#endif
end module df_mic_nonlinear_term_mod
