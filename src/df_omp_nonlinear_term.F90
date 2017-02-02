#include "redef.h"
#include "intrinsic_sizes.h"

module df_omp_nonlinear_term_mod
  use df_nonlinear_term_mod

  use par_mod, only: imag, pi, ve_max, phase_1, phase_2, spec
  use par_other, only: n_fields
  use communications, only: MPI_COMPLEX_TYPE, MPI_REAL_TYPE, mpi_comm_y,&
       &my_2reals_max_to_all, my_barrier, mpi_comm_xy, reduce_per_thread,&
       &threadlocal_mpi_comm_y
  use discretization
  use fourier
  use mpi
#ifdef WITHOMP
  use omp_lib
#endif
  implicit none

  type, public, extends(df_nonlinear_term_t) :: df_omp_nonlinear_term_t
   contains
     procedure :: calc => calc_omp_nonlinearity_df
     procedure :: getType => getThisType
     procedure :: allocate_arrays => allocate_arrays_omp
  end type df_omp_nonlinear_term_t

  private
  Real:: ve_x_max_loc, ve_y_max_loc
contains

  function getThisType(this)
    class(df_omp_nonlinear_term_t) :: this
    character(len=MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "df_omp_nonlinear_term_t"
  end function getThisType

  subroutine allocate_arrays_omp(this,vexy,dgdxy)
    class(df_omp_nonlinear_term_t) :: this
    complex, dimension(:,:,:,:),allocatable :: vexy,dgdxy
    integer :: klmn

    allocate(vexy(li1:li2,lj1:lj2,2,this%lbg0))
    allocate(dgdxy(li1:li2,lj1:lj2,2,this%lbg0))
    ! dummy initialization for "first touch"
    !$OMP PARALLEL DO default(none) private(klmn) shared(vexy,dgdxy,this)
    do klmn=1,this%lbg0
       vexy(:,:,:,klmn) = cmplx(1.0,0.0)
       dgdxy(:,:,:,klmn) = cmplx(1.0,0.0)
    end do
    !$OMP END PARALLEL DO
    
  end subroutine allocate_arrays_omp

  !>Computes the nonlinearity for the (x) nonlocal version
  !!\todo Get rid of transposition using striding? make some arrays for arakawa allocatable
  !!\todo Check whether the nonlinearity prefactor is important for the CFL criterion
  Subroutine calc_omp_nonlinearity_df(this,gy_chi,g_block,vexy,dgdxy,localrhs,first)
    class(df_omp_nonlinear_term_t) :: this
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: gy_chi
    complex, dimension(lbi:ubi,lj1:lj2,1:*),intent(in) :: g_block
    Complex, Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
    Logical, Intent(in) :: first

    ! Local variables
    Complex, Dimension(li1:li2,lj1:lj2,1:this%lbg0) :: nonlin    
    Complex, Dimension(0:nj0-1,   0:li0da/n_procs_y-1) :: tmp_arr
    Real,    Dimension(0:ly0da-1, 0:li0da/n_procs_y-1) :: nonlin1_re
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,2,1:this%lbg0) ::  vexy_re, dgdxy_re
    !DEC$ ATTRIBUTES ALIGN:64 :: vexy_re, dgdxy_re, nonlin
    !Real:: ve_x_max_loc, ve_y_max_loc
    Integer :: i,j, klmn
    logical :: is_not_omp_parallel
    integer :: max_threads, omp_level,team_size,active_level
    !integer :: omp_get_level, omp_get_thread_num, omp_get_num_threads, omp_get_active_level, omp_get_max_threads, omp_get_team_size

    !write(*,"(A,ES20.10)") "START of calc_omp_nonlinearity: localrhs = ",sum(abs(localrhs(:,:,1:this%lbg0)))
#ifdef WITHOMP
    max_threads = omp_get_max_threads()
    team_size = omp_get_num_threads()
#if 0
    !active_level = omp_get_active_level()
    !omp_level = omp_get_level()
    write(*,"(5(A,I3))") "omp_level = ",omp_level,&
         &", active_level = ", active_level,&
         &", num_threads = ",omp_get_num_threads(),&
         &", team_size(omp_level) = ", omp_get_team_size(omp_level),&
         &", team_size(omp_level-1) = ", omp_get_team_size(omp_level-1)
#endif
#if 0
    if (is_not_omp_parallel) then
       write(*,"(4X,I3,A)") max_threads," OpenMP threads, not nested in df_omp_nonlinear_term."
    else
       write(*,"(4X,I3,A,I3)") max_threads," OpenMP threads, nested with a team_size of ", team_size
    end if
    write(*,"(A,I3,A)") "Starting new parallel region with ",max_threads-team_size+1," threads."
#endif
#endif
    !$OMP PARALLEL num_threads(max_threads-team_size+1) &
    !$OMP default(none) private(j,nonlin1_re,tmp_arr) &
    !$OMP shared(turbdeal,lj1,lj2,li0da,li1da,li2da,nj0) &
    !$OMP shared(this,first,vexy_re, dgdxy_re,nonlin,nl_tmp1) &
    !$OMP reduction(max: ve_x_max_loc) reduction(max: ve_y_max_loc)

    if (turbdeal) then
       do j=lj1,lj2
          vexy(:,j,:,1:this%lbg0)  = vexy(:,j,:,1:this%lbg0)* this%prefactor%ranshift_to_df(j)
          dgdxy(:,j,:,1:this%lbg0) = dgdxy(:,j,:,1:this%lbg0)* this%prefactor%ranshift_to_df(j)
       enddo
    endif

#if 0
    !$OMP BARRIER
    !$OMP MASTER
    write(*,"(2(A,ES17.10))") "BEFORE FFT: vexy = ",&
         &sum(abs(vexy(:,:,:,1:this%lbg0))),&
         &", dgdxy = ",sum(abs(dgdxy(:,:,:,1:this%lbg0)))!,&
    !&", localrhs = ",sum(abs(localrhs(:,:,1:this%lbg0)))
    !$OMP END MASTER
    !$OMP BARRIER
#endif
    ! anti-aliasing and y fourier transform
    PERFON('deal_FFT')
    Call gto_real(this,vexy,vexy_re,2*this%lbg0)
    Call gto_real(this,dgdxy,dgdxy_re,2*this%lbg0)
    PERFOFF

#if 0
    !$OMP BARRIER
    !$OMP MASTER
    ! DEBUGGING
    write(*,"(2(A,ES17.10))") "AFTER FFT: vexy_re = ",&
         &sum(abs(vexy_re)),", dgdxy_re = ",sum(abs(dgdxy_re))

#if 0
    print*,"First block, first part."
    do i=0,li0da-1
       do j=0,ly0da-1
          write(*,"(F10.6,1X)",advance="no") vexy_re(j,i,1,1)
       end do
       write(*,"(A)") ""
    end do
    print*,"First block, second part."
    do i=0,li0da-1
       do j=0,ly0da-1
          write(*,"(F10.6,1X)",advance="no") vexy_re(j,i,2,1)
       end do
       write(*,"(A)") ""
    end do
#endif
    !$OMP END MASTER
    !$OMP BARRIER
#endif

    ! get max ExB velocity
    if (first) then
       PERFON('max_red')
#if 0
       print*,"First part:"
       do klmn=1,this%lbg0
          write(*,"(F10.6)",advance="no") maxval(vexy_re(:,:,1,klmn))
       end do
       write(*,"(A)") ""
       
       print*,"Second part:"
       do klmn=1,this%lbg0
          write(*,"(F10.6)",advance="NO") maxval(vexy_re(:,:,2,klmn))
       end do
       write(*,"(A)") ""
#endif

       !$OMP DO  private(klmn)
       do klmn=1,this%lbg0
          ve_x_max_loc = maxval(vexy_re(:,:,2,klmn))
          ve_y_max_loc = maxval(vexy_re(:,:,1,klmn))
       end do
       !$OMP END DO

       !ve_x_max_loc=maxval(vexy_re(:,:,2,:))
       !ve_y_max_loc=maxval(vexy_re(:,:,1,:))
       !write(*,"(A,2ES17.10)") "ve_max_loc = ",ve_x_max_loc, ve_y_max_loc
       !ve_max(1)=max(ve_max(1),ve_x_max_loc)        
       !ve_max(2)=max(ve_max(2),ve_y_max_loc)
       PERFOFF
    end if

    !compute the 'standard' nonlinear term
    PERFON('comp_nl')
    !nonlin1_re = -vexy_re(:,:,1,:)*dgdxy_re(:,:,2,:) + vexy_re(:,:,2,:)*dgdxy_re(:,:,1,:)

       !print*,"Starting the OpenMP do loop for klmn"
       
       !$OMP DO
       do klmn=1,this%lbg0
          PERFON('nl_kernl')
          nonlin1_re = -vexy_re(:,:,1,klmn)*dgdxy_re(:,:,2,klmn) &
               & + vexy_re(:,:,2,klmn)*dgdxy_re(:,:,1,klmn)
          PERFOFF
          !write(*,"(A,I4,A,ES20.10)") "klmn = ",klmn,", nonlin1_re = ",sum(abs(nonlin1_re))
          ! transform back to fourier space
          Call to_fourier_y(nonlin1_re,tmp_arr)
          !write(*,"(A,I4,A,ES20.10)") "klmn = ",klmn,", after ffty = ",sum(abs(tmp_arr))
          Call transpose_cmplx(nj0, li0da, tmp_arr, 0, nonlin(:,:,klmn), 0)
          !write(*,"(A,I4,A,ES20.10)") "klmn = ",klmn,", after transpose = ",sum(abs(nonlin(:,:,klmn)))
       end do
       !$OMP END DO
       !print*,"Ending the OpenMP do loop."
       
    PERFOFF

#if 0
    !$OMP BARRIER
    !$OMP MASTER
    ! DEBUGGING
    write(*,"(A,ES17.10)") "AFTER iFFT: nonlin = ",sum(abs(nonlin))
    !$OMP END MASTER
#endif

    !$OMP BARRIER
#if 0
    write(*,"(I3,A,ES17.10)") omp_get_thread_num(),": BEFORE MULTIPLY_WITH: localrhs = ",&
         &sum(abs(localrhs(:,:,1:this%lbg0)))
    !$OMP BARRIER
#endif

    PERFON('pref_mul')
    ! multiply with the prefactor and write to localrhs
    call this%prefactor%multiply_with(nonlin,localrhs,this%lbg0)
    PERFOFF

    !$OMP END PARALLEL

    if (first) then
       ve_max(1)=max(ve_max(1),ve_x_max_loc)        
       ve_max(2)=max(ve_max(2),ve_y_max_loc)
    end if

    !PERFOFF_I

  End Subroutine calc_omp_nonlinearity_df


  !> Go to real space
  subroutine gto_real(this,inarr,outarr,howmany)
    class(df_nonlinear_term_t) :: this
    integer, intent(IN) :: howmany
    Complex, Dimension(li1:li2,lj1:lj2,1:howmany),Intent(in) :: inarr
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,1:howmany), target, Intent(out) :: outarr

    ! Local variables
    !Complex, Dimension(li1da:li2da,lj1:lj2,1:howmany) :: tmp_arr1
    !Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr2
    complex, dimension(0:ly0da/2,0:li0da/n_procs_y-1) :: extended_tarray
    real, dimension(:,:), pointer :: p_out
    !real :: local_sum,thread_global_sum(0:7,0:n_procs_x*n_procs_y)
    Integer:: klmn, i_block,number_of_lbg0_blocks
    integer :: my_thread
#ifdef WITHOMP
    integer :: omp_get_thread_num

    my_thread = omp_get_thread_num()
#else
    my_thread = 0
#endif
#if 0
    number_of_lbg0_blocks = howmany/this%lbg0
    if (number_of_lbg0_blocks.eq.0) then
       ! Transpose x-y
       Call transpose_cmplx(li0da, nj0, inarr(:,:,1), 0,tmp_arr2 , 0)
       ! Fourier transfrom in y (include dealiasing step)
       Call to_real_y(tmp_arr2,outarr(:,:,1))
    else
       do i_block=1,number_of_lbg0_blocks
#endif

             !$OMP DO private(p_out) schedule(static)
             do klmn=1,howmany
                !tmp_arr1(:,:,(i_block-1)*this%lbg0+klmn)=inarr(:,:,(i_block-1)*this%lbg0+klmn)
                !write(*,"(2I3,A,ES20.10)") my_thread, klmn," before transpose : ", &
                !     & sum(abs(inarr(:,:,klmn)))
                ! Transpose x-y
                !Call transpose_cmplx(li0da, nj0, inarr(:,:,klmn), 0,tmp_arr2 , 0)
                Call transpose_and_extend_cmplx(li0da, nj0, inarr(:,:,klmn),&
                     &extended_tarray , ly0da/2+1)
                !write(*,"(2I3,A,ES20.10)") my_thread, klmn," after  transpose : ", &
                !     & sum(abs(tmp_arr2))
                
                ! Fourier transfrom in y (include dealiasing step)
                !p_out => outarr(:,:,(i_block-1)*this%lbg0+klmn)
                p_out => outarr(:,:,klmn)
                Call to_real_y_no_extend(extended_tarray,p_out) 
                !Call to_real_y(tmp_arr2,p_out)
                !p_out = 1.0d0
                !write(*,"(2I3,A,ES20.10)") my_thread, klmn," after fft : ", &
                !     & sum(abs(p_out))
             end do
             !$OMP END DO
             !print*,"After the klmn loop in gto_real"
#if 0
       end do
    end if
#endif

  end subroutine gto_real
end module df_omp_nonlinear_term_mod
