#include "redef.h"

Module fourier_mic
  use omp_lib
  USE discretization, ONLY: nj0,nky0,lj0,li1,li2,lj1,lj2,n_procs_y
  Use mkl_dfti
  use mkl_service
#if defined(WITH_OFFLOAD) && defined(WITHPERF)
  use perflib
#endif

  Implicit None

  public :: initialize_fourier,fft_ky_to_y,fft_y_to_ky,&
       to_real_y,to_real_y_block,to_fourier_y,finalize_fourier,&
       check_fourier

  private

#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: inc_1x,num_x, ub_out_x,li0da,ly0da
#endif
  Integer, save :: inc_1x, num_x, ub_out_x
  Integer, save :: li0da, ly0da

  Type mkldescriptor
     Type(DFTI_DESCRIPTOR), Pointer:: desc   
  end type mkldescriptor

#ifdef WITH_MIC_NONLIN
#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: fft_ky2y, fft_y2ky, fft_ky2y_block,fft_y2ky_block
#endif
#endif
  Type(mkldescriptor),dimension(:),allocatable, private:: fft_ky2y,fft_y2ky
  Type(mkldescriptor),dimension(:),allocatable, private:: fft_ky2y_block,fft_y2ky_block

  Complex, Dimension(:,:),allocatable,private,save  :: temparr_block
  !$OMP THREADPRIVATE(temparr_block)
#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: temparr_block
#endif

#ifdef WITHOMP
!  integer OMP_GET_THREAD_NUM
!  external OMP_GET_THREAD_NUM
#endif
 
Contains
  Subroutine check_fourier(print_ini_msg)
    Logical, intent(in) :: print_ini_msg

    CHARACTER(len=198) :: mklversion
    call MKL_Get_Version_String(mklversion)
    Write(*,"(A,A)") "Using ", mklversion

  End Subroutine check_fourier

  !
  ! Initialization for FFT routines
  !
  Subroutine initialize_fourier(a_li0da,a_ly0da,this_lbg0,klmn_per_thread)
    INTEGER, INTENT(IN) :: a_li0da, a_ly0da,this_lbg0
    integer, intent(OUT) :: klmn_per_thread
#ifdef WITH_MIC_NONLIN
#ifdef WITH_OFFLOAD
    !DEC$ attributes offload:mic :: initialize_fourier
#endif
#endif

    REAL   :: facnny,facnnx
    integer:: Status, thrnum, n_threads
    integer, dimension(2) :: stride_2d !to avoid temporary arrays
    integer :: ret_int
    logical :: verbose=.false.

    li0da=a_li0da
    ly0da=a_ly0da
    if (n_procs_y.eq.1) then
       inc_1x=ly0da/2+1
       num_x= nj0
       ub_out_x=ly0da/2+1
    else
       inc_1x=lj0
       num_x=lj0
       ub_out_x=lj0
    endif

    facnnx=1.0/Real(li0da)
    facnny=1.0/Real(ly0da)
    print*,"initialize_fourier_mic: before parallel OpenMP region"
    !$OMP PARALLEL default(none) shared(li0da,ly0da,n_procs_y,inc_1x,nj0) &
    !$OMP shared(num_x,ub_out_x,lj0,facnnx,facnny,verbose) &
    !$OMP shared(fft_ky2y,fft_y2ky,fft_ky2y_block,fft_y2ky_block) &
    !$OMP shared(klmn_per_thread,this_lbg0) &
    !$OMP private(ret_int,status,thrnum,stride_2d,n_threads)
    call mkl_set_num_threads(1)
    n_threads = omp_get_num_threads()
    !$OMP SINGLE
    if (modulo(this_lbg0,n_threads).ne.0) then
       write(*,"(A,I6,A,I4)") "blocksize lbg0=",this_lbg0," must be a multiple of n_threads=",n_threads
       flush(6)
       stop
    end if
    klmn_per_thread = this_lbg0/n_threads
    write(*,"(I4,A,I5)") omp_get_thread_num(),": klmn_per_thread = ",klmn_per_thread
    flush(6)
    !$OMP END SINGLE
#if defined(WITHOMP_BLOCKLOOP) || defined(WITH_MIC_NONLIN)
    !$OMP SINGLE
#endif

    allocate(fft_ky2y(0:n_threads-1),fft_y2ky(0:n_threads-1))
#ifdef EN_BLOC
    allocate(fft_ky2y_block(0:n_threads-1),fft_y2ky_block(0:n_threads-1))
#endif
    print*,"Running for all ",n_threads," threads."
#if defined(WITHOMP_BLOCKLOOP) || defined(WITH_MIC_NONLIN)
    !$OMP END SINGLE
    !$OMP BARRIER

#ifdef EN_BLOC
    allocate(temparr_block(0:ly0da/2,li0da/n_procs_y*2*klmn_per_thread))
    temparr_block(:,:) = cmplx(0.0,0.0)
#endif
    !$OMP DO  schedule(static)
#endif
    do thrnum=0,n_threads-1

       ! to real space in y
       Status = DftiCreateDescriptor( fft_ky2y(thrnum)%desc, DFTI_SINGLE,&
            DFTI_REAL, 1, ly0da)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y(thrnum)%desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y(thrnum)%desc, DFTI_NUMBER_OF_TRANSFORMS, li0da/n_procs_y)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y(thrnum)%desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
       IF (status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y(thrnum)%desc, DFTI_INPUT_DISTANCE, ly0da/2+1)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y(thrnum)%desc, DFTI_OUTPUT_DISTANCE, ly0da)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y(thrnum)%desc, DFTI_BACKWARD_SCALE, 1.)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiCommitDescriptor(fft_ky2y(thrnum)%desc)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)

#ifdef EN_BLOC
       ! to real space in y blocked
       Status = DftiCreateDescriptor( fft_ky2y_block(thrnum)%desc, DFTI_SINGLE,&
            DFTI_REAL, 1, ly0da)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y_block(thrnum)%desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y_block(thrnum)%desc, &
            & DFTI_NUMBER_OF_TRANSFORMS, 2*li0da/n_procs_y*klmn_per_thread)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y_block(thrnum)%desc, &
            & DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
       IF (status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y_block(thrnum)%desc, DFTI_INPUT_DISTANCE, ly0da/2+1)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y_block(thrnum)%desc, DFTI_OUTPUT_DISTANCE, ly0da)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_ky2y_block(thrnum)%desc, DFTI_BACKWARD_SCALE, 1.)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiCommitDescriptor(fft_ky2y_block(thrnum)%desc)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
#endif
       ! to fourier space in y
       Status = DftiCreateDescriptor( fft_y2ky(thrnum)%desc, DFTI_SINGLE,&
            DFTI_REAL, 1, ly0da)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_y2ky(thrnum)%desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_y2ky(thrnum)%desc, DFTI_NUMBER_OF_TRANSFORMS, li0da/n_procs_y)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_y2ky(thrnum)%desc, DFTI_INPUT_DISTANCE, ly0da)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_y2ky(thrnum)%desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
       IF (status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_y2ky(thrnum)%desc, DFTI_OUTPUT_DISTANCE, ly0da/2+1)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_y2ky(thrnum)%desc, DFTI_FORWARD_SCALE, facnny)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiCommitDescriptor(fft_y2ky(thrnum)%desc)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
    
    Enddo
#if defined(WITHOMP_BLOCKLOOP) || defined(WITH_MIC_NONLIN)
    !$OMP END DO
#endif
    thrnum=0

    if (verbose) then
       !$OMP DO ORDERED
       do thrnum=0,n_threads-1
          !$OMP ORDERED
          write(*,"(I4)") thrnum
          !DFTI_PRECISION
          !DFTI_FORWARD_DOMAIN
          !DFTI_DIMENSION, DFTI_LENGTHS
          !DFTI_PLACEMENT
          !DFTI_FORWARD_SCALE, DFTI_BACKWARD_SCALE
          !DFTI_NUMBER_OF_USER_THREADS
          !DFTI_THREAD_LIMIT
          !DFTI_INPUT_STRIDES, DFTI_OUTPUT_STRIDES
          !DFTI_NUMBER_OF_TRANSFORMS
          !DFTI_INPUT_DISTANCE, DFTI_OUTPUT_DISTANCE
          !DFTI_COMPLEX_STORAGE, DFTI_REAL_STORAGE, DFTI_CONJUGATE_EVEN_STORAGE
          !DFTI_PACKED_FORMAT
          !DFTI_WORKSPACE
          !DFTI_COMMIT_STATUS
          !DFTI_ORDERING
          status = DftiGetValue(fft_ky2y(thrnum)%desc,DFTI_PRECISION,ret_int)
          print*,"PRECISION = ",ret_int

          status = DftiGetValue(fft_ky2y(thrnum)%desc,DFTI_FORWARD_DOMAIN,ret_int)
          print*,"FORWARD_DOMAIN = ",ret_int

          status = DftiGetValue(fft_ky2y(thrnum)%desc,DFTI_LENGTHS,ret_int)
          print*,"LENGTHS = ",ret_int
          !$OMP END ORDERED
       end do
       !$OMP END DO
    end if
    !$OMP END PARALLEL

  End Subroutine initialize_fourier


#ifdef WITH_MIC_NONLIN
#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: fft_ky_to_y
#endif
#endif
  subroutine fft_ky_to_y(inarr, outarr)
    complex,dimension((ly0da/2+1)*(li0da/n_procs_y)),intent(in):: inarr
    real, dimension(0:(ly0da)*(li0da/n_procs_y)-1),intent(out):: outarr
    INTEGER :: status,thrnum

    !PERFON('fft_ky2y')  
#ifdef WITHOMP
    thrnum=OMP_GET_THREAD_NUM()
#else
    thrnum = 0
#endif
    !write(*,"(10X,I4,A,ES17.10)") thrnum,": inarr = ",sum(real(inarr*conjg(inarr)))
    !flush(6)    
    status=DftiComputeBackward(fft_ky2y(thrnum)%desc, inarr, outarr)
    If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
    !PERFOFF
    !write(*,"(10X,I4,A,ES17.10)") thrnum,": outarr = ",sum(outarr*outarr)
    !flush(6)
  end subroutine fft_ky_to_y

#ifdef WITH_MIC_NONLIN
#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: fft_ky_to_y_block
#endif
#endif
  subroutine fft_ky_to_y_block(inarr, outarr,howmany)
    integer :: howmany
    complex,dimension((ly0da/2+1)*howmany),intent(in) :: inarr
    real,   dimension(0:ly0da*howmany-1),  intent(out):: outarr
    INTEGER :: status,thrnum

#ifdef WITHOMP
    thrnum=OMP_GET_THREAD_NUM()
#else
    thrnum = 0
#endif
    !write(*,"(10X,I4,A,ES17.10)") thrnum,": inarr = ",sum(real(inarr*conjg(inarr)))
    !flush(6)    
    status=DftiComputeBackward(fft_ky2y_block(thrnum)%desc, inarr, outarr)
    If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
    !write(*,"(10X,I4,A,ES17.10)") thrnum,": outarr = ",sum(outarr*outarr)
    !flush(6)
  end subroutine fft_ky_to_y_block

#ifdef WITH_MIC_NONLIN
#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: fft_y_to_ky
#endif
#endif
  subroutine fft_y_to_ky(inarr, outarr)
    real, dimension((ly0da)*(li0da/n_procs_y)), intent(in):: inarr
    complex, dimension(0:(ly0da/2+1)*(li0da/n_procs_y)-1), intent(out):: outarr
    INTEGER :: status,thrnum

    PERFON('fft_y2ky')  
#ifdef WITHOMP
    thrnum=OMP_GET_THREAD_NUM()
#else
    thrnum = 0
#endif    
    status=DftiComputeForward(fft_y2ky(thrnum)%desc, inarr, outarr)
    If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
    PERFOFF

  end subroutine fft_y_to_ky

  !global routines
! transform 2-dimensional array inarr(ky,x) to outarr(y,x)
! note : y is the first coordinate
#ifdef WITH_MIC_NONLIN
#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: to_real_y
#endif
#endif
  Subroutine to_real_y(inarr,outarr)
    Complex, Dimension(0:nky0-1, 0:li0da/n_procs_y-1)  , intent(in) :: inarr
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1), intent(out) :: outarr
    Complex, Dimension(0:ly0da/2,0:li0da/n_procs_y-1)              :: temparr

    integer :: j
    logical,save :: isfirst=.true.

    !to_real_y is the first fft to be called in nonlocal-x code
    !therefore, it detects and sets the assigned thread number
    
    !!$OMP CRITICAL
    !write(*,"(I3,A,ES17.10)") omp_get_thread_num(),": to_real_y: inarr = ",sum(real(inarr*conjg(inarr)))
    !flush(6)

    if (isfirst) then
       !print*,omp_get_thread_num(),"Measuring first iteration."
       !$OMP BARRIER
       isfirst=.false.
       PERFON('toR_y1_1')
    else
       !print*,omp_get_thread_num(),"Measuring not first iteration."
       PERFON('toR_y1_2')
    end if
    do j=0,li0da/n_procs_y-1
       temparr(0:nky0-1, j)=inarr(:,j)
       ! Dealiasing in y direction
       temparr(nky0:,j)=cmplx(0.0,0.0)
    end do
    PERFOFF

    !write(*,"(5X,A,ES17.10)") "to_real_y: temparr = ",sum(real(temparr*conjg(temparr)))
    !flush(6)

    PERFON('toR_y2')
    call fft_ky_to_y(temparr,outarr)
    PERFOFF

    !write(*,"(5X,A,ES17.10)") "to_real_y: outarr = ",sum(real(outarr*outarr))
    !flush(6)
    !!$OMP END CRITICAL

  end Subroutine to_real_y

#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: to_real_y_block
#endif
  Subroutine to_real_y_block(inarr,outarr,howmany)
    integer :: howmany
    Complex, Dimension(0:nky0-1,1:howmany),  intent(in) :: inarr
    Real,    Dimension(0:ly0da-1,1:howmany), intent(out) :: outarr

    ! Local variables
    !Complex, Dimension(0:ly0da/2,1:howmany)  :: temparr_block
    integer :: i,j
    real :: coeffsum

    !write(*,"(I3,A,ES17.10)") omp_get_thread_num(),&
    !     &": to_real_y_block: inarr = ",sum(real(conjg(inarr)*inarr))
    !flush(6)

    PERFON('dealcopy')
    do j=1,howmany
       temparr_block(0:nky0-1,j) = inarr(:,j)
       temparr_block(nky0:,j) = cmplx(0.0,0.0)
    end do
    PERFOFF

    !write(*,"(5X,A,ES17.10)") "to_real_y_block: temparr = ",sum(real(conjg(temparr)*temparr))
    !flush(6)

    PERFON('toR_y2_b')
    call fft_ky_to_y_block(temparr_block,outarr,howmany)
    PERFOFF

#if 0
    do j=1,howmany
       coeffsum = (real(conjg(temparr(0,j))*temparr(0,j))&
            &+real(conjg(temparr(ly0da/2,j))*temparr(ly0da/2,j))&
            &+2*sum(real(conjg(temparr(1:ly0da/2-1,j))*temparr(1:ly0da/2-1,j))))
       if (abs(ly0da*coeffsum-sum(outarr(:,j)**2)).gt.1e-15) then
          write(*,"(I4,I8,2ES17.10)") omp_get_thread_num(),j,&
               &ly0da*coeffsum,&
               &sum(outarr(:,j)*outarr(:,j))
          flush(6)
       end if
    end do
#endif
    !write(*,"(5X,A,ES17.10)") "to_real_y: outarr = ",sum(outarr*outarr)
    !flush(6)

  end Subroutine to_real_y_block

! transform 2-dimensional array inarr(y,x) to outarr(ky,x)
! note : y is the first coordinate
#ifdef WITH_MIC_NONLIN
#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: to_fourier_y
#endif
#endif
  Subroutine to_fourier_y(inarr,outarr)
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1), intent(in)   :: inarr
    Complex, Dimension(0:nky0-1, 0:li0da/n_procs_y-1)  , intent(out) :: outarr
    Complex, Dimension(0:ly0da/2,0:li0da/n_procs_y-1)              :: temparr

    PERFON('toF_y')

    call fft_y_to_ky(inarr,temparr)

    ! Removing last element of temparr for dealising
    outarr=temparr(0:nky0-1,:)
    
    PERFOFF

  end Subroutine to_fourier_y

  Subroutine finalize_fourier
    INTEGER :: status, thrnum
#ifdef WITH_OFFLOAD
  !DEC$ attributes offload:mic :: finalize_fourier
#endif
    
    integer :: n_threads
    n_threads = omp_get_num_threads()
    !$OMP PARALLEL default(none) shared(fft_ky2y,fft_y2ky,n_threads)
    !$OMP DO private(thrnum, status)
    do thrnum=0,n_threads-1
       Status = DftiFreeDescriptor( fft_ky2y(thrnum)%desc)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
#ifdef EN_BLOC
       Status = DftiFreeDescriptor( fft_ky2y_block(thrnum)%desc)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
#endif
       Status = DftiFreeDescriptor( fft_y2ky(thrnum)%desc)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
    enddo
    !$OMP END DO
    !$OMP SINGLE
    deallocate(fft_ky2y,fft_y2ky)
#ifdef EN_BLOC
    deallocate(fft_ky2y_block,fft_y2ky_block)
    deallocate(temparr_block)
#endif
    !$OMP END SINGLE
    !$OMP END PARALLEL
  end Subroutine finalize_fourier

End Module fourier_mic
