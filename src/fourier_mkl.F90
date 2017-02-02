#include "redef.h"

Module fourier
  USE discretization, ONLY: ni0,nj0,nky0,li0,lj0,li1,li2,lj1,lj2,xy_local,n_procs_x,n_procs_y,mype, my_pex
  USE communications, ONLY: omp_num_threads,MPI_COMPLEX_TYPE,mpi_comm_x
  Use mkl_dfti
  Use par_in, ONLY: fourier2d
  Implicit None

  public :: initialize_fourier,fft_kx_to_x,fft_ky_to_y,fft_y_to_ky,fft_x_to_kx,&
       fft_ff_to_xy_t,fft_xy_to_ff_t,to_real_y,to_real_y_no_extend,to_fourier_y,finalize_fourier,&
       initialize_fourier_x_1d,to_fourier_x_1d,to_real_x_1d,finalize_fourier_x_1d,&
       check_fourier,fft_ff_to_xy,fft_xy_to_ff, &
       initialize_fourier_boundary, to_real_boundary,to_fourier_boundary,&
       finalize_fourier_boundary

  private

  Integer:: inc_1x, inc_2x, num_x, ub_out_x
  Integer:: li0da, ly0da

  Type mkldescriptor
     Type(DFTI_DESCRIPTOR), Pointer:: desc   
  end type mkldescriptor
  Type(mkldescriptor),dimension(:),allocatable, private:: fft_kx2x, fft_x2kx,fft_ky2y,fft_y2ky
  Type(mkldescriptor),dimension(:),allocatable, private:: fft_kx2x_1d, fft_x2kx_1d

#ifdef WITHOMP
  integer OMP_GET_THREAD_NUM
  external OMP_GET_THREAD_NUM
#endif
 
Contains
  Subroutine check_fourier(print_ini_msg)
    Logical, intent(in) :: print_ini_msg

    CHARACTER(len=198) :: mklversion
    If ((mype==0).and.print_ini_msg) Then
       !call MKLGetVersionString(mklversion) !in older MKL versions
       call MKL_Get_Version_String(mklversion)
       Write(*,"(A,A)") "Using ", mklversion
       If (fourier2d) Write(*,"(A)") "2D fourier analysis only present when using fftw lib, switching to 1D"
    Endif
    fourier2d=.false.

  End Subroutine check_fourier

  !
  ! Initialization for FFT routines
  !
  Subroutine initialize_fourier(a_li0da,a_ly0da)
    INTEGER, INTENT(IN) :: a_li0da, a_ly0da
    REAL   :: facnny,facnnx
    integer:: Status, thrnum
    integer, dimension(2) :: stride_2d !to avoid temporary arrays

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

    allocate(fft_kx2x(0:omp_num_threads-1),fft_x2kx(0:omp_num_threads-1),&
         fft_ky2y(0:omp_num_threads-1),fft_y2ky(0:omp_num_threads-1))

    do thrnum=0,omp_num_threads-1

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
    
       If (xy_local) then

          ! to real space in x
          Status = DftiCreateDescriptor( fft_kx2x(thrnum)%desc, DFTI_SINGLE,&
               DFTI_COMPLEX, 1, li0da)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          Status = DftiSetValue( fft_kx2x(thrnum)%desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          Status = DftiSetValue( fft_kx2x(thrnum)%desc, DFTI_NUMBER_OF_TRANSFORMS, num_x)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          Status = DftiSetValue( fft_kx2x(thrnum)%desc, DFTI_INPUT_DISTANCE, li0da)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          Status = DftiSetValue( fft_kx2x(thrnum)%desc, DFTI_OUTPUT_DISTANCE, 1)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)

          stride_2d = (/0,1/)
          Status = DftiSetValue( fft_kx2x(thrnum)%desc, DFTI_INPUT_STRIDES, stride_2d)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          stride_2d = (/0,inc_1x/)
          Status = DftiSetValue( fft_kx2x(thrnum)%desc, DFTI_OUTPUT_STRIDES, stride_2d)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status) 
          
          Status = DftiSetValue( fft_kx2x(thrnum)%desc, DFTI_BACKWARD_SCALE, 1.0)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          Status = DftiCommitDescriptor(fft_kx2x(thrnum)%desc)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)

          ! to fourier space in x
          Status = DftiCreateDescriptor( fft_x2kx(thrnum)%desc, DFTI_SINGLE,&
               DFTI_COMPLEX, 1, li0da)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          Status = DftiSetValue( fft_x2kx(thrnum)%desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)

          Status = DftiSetValue( fft_x2kx(thrnum)%desc, DFTI_NUMBER_OF_TRANSFORMS, num_x)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          Status = DftiSetValue( fft_x2kx(thrnum)%desc, DFTI_INPUT_DISTANCE, 1)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)

          Status = DftiSetValue( fft_x2kx(thrnum)%desc, DFTI_OUTPUT_DISTANCE, li0da)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          stride_2d = (/0,inc_1x/)
          Status = DftiSetValue( fft_x2kx(thrnum)%desc, DFTI_INPUT_STRIDES, stride_2d)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)

          stride_2d = (/0,1/)
          Status = DftiSetValue( fft_x2kx(thrnum)%desc, DFTI_OUTPUT_STRIDES, stride_2d)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          Status = DftiSetValue( fft_x2kx(thrnum)%desc, DFTI_FORWARD_SCALE, facnnx)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
          Status = DftiCommitDescriptor(fft_x2kx(thrnum)%desc)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          
       Endif
    Enddo
    
    thrnum=0

  End Subroutine initialize_fourier

  subroutine fft_kx_to_x(inarr, outarr)
    complex,dimension(0:li0da*lj0-1) :: inarr    
    complex,dimension(0:(ub_out_x+1)*li0da-1) :: outarr
    INTEGER :: status, thrnum

    PERFON('fft_kx2x')  
#ifdef WITHOMP
    thrnum=OMP_GET_THREAD_NUM()
#else
    thrnum = 0
#endif
    status=DftiComputeBackward(fft_kx2x(thrnum)%desc, inarr, outarr)
    If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
    PERFOFF

  end subroutine fft_kx_to_x

  subroutine fft_ky_to_y(inarr, outarr)
    complex,dimension((ly0da/2+1)*(li0da/n_procs_y)),intent(in):: inarr
    real, dimension(0:(ly0da)*(li0da/n_procs_y)-1),intent(out):: outarr
    INTEGER :: status,thrnum

    PERFON('fft_ky2y')  
#ifdef WITHOMP
    thrnum=OMP_GET_THREAD_NUM()
#else
    thrnum = 0
#endif    
    status=DftiComputeBackward(fft_ky2y(thrnum)%desc, inarr, outarr)
    If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
    PERFOFF

  end subroutine fft_ky_to_y

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

  subroutine fft_x_to_kx(inarr, outarr)
    complex,dimension((ub_out_x+1)*li0da),intent(in):: inarr
    complex,dimension(0:li0da*lj0-1),intent(out):: outarr
    INTEGER :: status,thrnum

    PERFON('fft_x2kx')  
#ifdef WITHOMP
    thrnum=OMP_GET_THREAD_NUM()
#else
    thrnum = 0
#endif    
    status=DftiComputeForward(fft_x2kx(thrnum)%desc, inarr, outarr)
    If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
    PERFOFF

  end subroutine fft_x_to_kx

!-------------Two dimensional ffts - These are just dummy routines for compilation
!                                    they are never called

  subroutine fft_ff_to_xy(inarr,outarr)
    complex,dimension(0:ly0da/2,0:li0da-1), intent(in)::inarr
    real,dimension(0:ly0da-1,0:li0da-1), intent(inout)::outarr

    STOP '2d FFT not implemented for MKL library.'

  end subroutine fft_ff_to_xy

  subroutine fft_xy_to_ff(inarr,outarr)
    real,dimension(0:ly0da-1,0:li0da-1), intent(in)::inarr
    complex,dimension(0:ly0da/2,0:li0da-1), intent(inout)::outarr

    STOP '2d FFT not implemented for MKL library.'

  end subroutine fft_xy_to_ff

  subroutine fft_ff_to_xy_t(inarr,outarr)
    complex,dimension(0:ly0da/2,0:li0da-1), intent(in)::inarr
    real,dimension(0:ly0da-1,0:li0da-1), intent(inout)::outarr

    STOP '2d FFT not implemented for MKL library.'

  end subroutine fft_ff_to_xy_t

  subroutine fft_xy_to_ff_t(inarr,outarr)
    real,dimension(0:ly0da-1,0:li0da-1), intent(in)::inarr
    complex,dimension(0:ly0da/2,0:li0da-1), intent(inout)::outarr

    STOP '2d FFT not implemented for MKL library.'

  end subroutine fft_xy_to_ff_t

  !global routines
! transform 2-dimensional array inarr(ky,x) to outarr(y,x)
! note : y is the first coordinate
  Subroutine to_real_y(inarr,outarr)
    Complex, Dimension(0:nky0-1, 0:li0da/n_procs_y-1)  , intent(in) :: inarr
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1), intent(out) :: outarr
    Complex, Dimension(0:ly0da/2,0:li0da/n_procs_y-1)              :: temparr

    !to_real_y is the first fft to be called in nonlocal-x code
    !therefore, it detects and sets the assigned thread number

    PERFON('toR_y1')    
    temparr(0:nky0-1, 0:li0da/n_procs_y-1)=inarr
    ! Dealiasing in y direction
    temparr(nky0:,:)=0
    PERFOFF

    PERFON('toR_y2')
    call fft_ky_to_y(temparr,outarr)
    PERFOFF

  end Subroutine to_real_y

  Subroutine to_real_y_no_extend(inarr,outarr)
    Complex, Dimension(0:ly0da/2, 0:li0da/n_procs_y-1), intent(in) :: inarr
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1), intent(out) :: outarr

    PERFON('toR_y2')
    call fft_ky_to_y(inarr,outarr)
    PERFOFF

  end Subroutine to_real_y_no_extend

! transform 2-dimensional array inarr(y,x) to outarr(ky,x)
! note : y is the first coordinate
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

    do thrnum=0,omp_num_threads-1
       Status = DftiFreeDescriptor( fft_ky2y(thrnum)%desc)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       Status = DftiFreeDescriptor( fft_y2ky(thrnum)%desc)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       If (xy_local) then
          Status = DftiFreeDescriptor( fft_kx2x(thrnum)%desc)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
          Status = DftiFreeDescriptor( fft_x2kx(thrnum)%desc)
          If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       Endif

       deallocate(fft_kx2x,fft_x2kx,fft_ky2y,fft_y2ky)

     
    enddo
    
  end Subroutine finalize_fourier


  !The one-dimensional FFT in radial direction is only rarely needed
  Subroutine initialize_fourier_x_1d
    integer:: Status, thrnum
    
    allocate(fft_kx2x_1d(0:omp_num_threads-1),fft_x2kx_1d(0:omp_num_threads-1))
    
    do thrnum=0,omp_num_threads-1
       
       ! to real space in x (1 dimensional)
       
       Status = DftiCreateDescriptor( fft_kx2x_1d(thrnum)%desc, DFTI_SINGLE,&
            DFTI_COMPLEX, 1, ni0)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_kx2x_1d(thrnum)%desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_kx2x_1d(thrnum)%desc, DFTI_BACKWARD_SCALE, 1.0)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiCommitDescriptor(fft_kx2x_1d(thrnum)%desc)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       ! to fourier space in x (1 dimensional)
       
       Status = DftiCreateDescriptor( fft_x2kx_1d(thrnum)%desc, DFTI_SINGLE,&
            DFTI_COMPLEX, 1, ni0)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_x2kx_1d(thrnum)%desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiSetValue( fft_x2kx_1d(thrnum)%desc, DFTI_FORWARD_SCALE, 1.0/REAL(ni0))
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
       
       Status = DftiCommitDescriptor(fft_x2kx_1d(thrnum)%desc)
       If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)      
       
    Enddo
    
    thrnum=0
    
  End Subroutine initialize_fourier_x_1d

  ! transform the 1-dimensional array inarr(x) to outarr(kx)
  Subroutine to_fourier_x_1d(inarr,outarr)
    Complex, Dimension(li1:li2), intent(in) :: inarr
    Complex, Dimension(li1:li2), intent(out) :: outarr

    ! Local variables
    COMPLEX, DIMENSION(ni0),target :: local_fullarray,local_Fourierarray
    INTEGER :: ierr, status, thrnum
    
    PERFON('toF_x_1d')   
#ifdef WITHOMP
    thrnum = OMP_GET_THREAD_NUM()
#else
    thrnum = 0
#endif 
    IF (n_procs_x.GT.1) THEN
       ! first gather the array on one processor, as we do not have a parallelized Fourier transformation
       CALL MPI_Gather(inarr(li1),li0,MPI_COMPLEX_TYPE,&
            & local_fullarray(1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
       IF (my_pex.EQ.0) THEN
          ! do Fourier transform on process 0
          status=DftiComputeForward(fft_x2kx_1d(thrnum)%desc,local_fullarray,local_fourierarray)
          IF (Status /= 0) PRINT *, __LINE__, DftiErrorMessage(status)
       END IF
       ! now scatter the result array
       CALL MPI_Scatter(local_Fourierarray(1),li0,MPI_COMPLEX_TYPE,&
            & outarr(li1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
    ELSE
       ! do Fourier transform
       status=DftiComputeForward(fft_x2kx_1d(thrnum)%desc,inarr, outarr)
       IF (Status /= 0) PRINT *, __LINE__, DftiErrorMessage(status)
    END IF
    
    !set unphysical mode to zero
    !if (evenx.EQ.1) outarr(ni0/2) = 0.0
    PERFOFF
    
  end Subroutine to_fourier_x_1d

!!! WITHOUT DEALIASING
! transform the 1-dimensional array inarr(kx) to outarr(x)
  Subroutine to_real_x_1d(inarr,outarr)
    Complex, Dimension(li1:li2) :: inarr
    Complex, Dimension(li1:li2), intent(out) :: outarr

    ! Local variables
    COMPLEX, DIMENSION(ni0),target :: local_fullarray,local_Fourierarray
    INTEGER :: ierr, status, thrnum

    PERFON('toR_x_1d')     
#ifdef WITHOMP
    thrnum=OMP_GET_THREAD_NUM()
#else
    thrnum = 0
#endif    
    !set unphysical mode to zero
    !if (evenx.EQ.1) inarr(ni0/2)=0.0
    IF (n_procs_x.GT.1) THEN
       CALL MPI_Gather(inarr(li1),li0,MPI_COMPLEX_TYPE,&
            & local_Fourierarray(1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
       IF (my_pex.EQ.0) THEN
          ! do Fourier transform on process 0
          status=DftiComputeBackward(fft_kx2x_1d(thrnum)%desc,local_Fourierarray,local_fullarray)
          IF (Status /= 0) PRINT *, __LINE__, DftiErrorMessage(status)
       END IF
       ! now scatter the result array
       CALL MPI_Scatter(local_fullarray(1),li0,MPI_COMPLEX_TYPE,&
            & outarr(li1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
    ELSE
       status=DftiComputeBackward(fft_kx2x_1d(thrnum)%desc,inarr, outarr)
       IF (Status /= 0) PRINT *, __LINE__, DftiErrorMessage(status)
    END IF
    PERFOFF

  end Subroutine to_real_x_1d

  !dummy routines for y global boundary condition

  subroutine initialize_fourier_boundary

  end subroutine initialize_fourier_boundary

  Subroutine to_fourier_boundary(inarr,outarr)
    Complex, Dimension(li1:li2,lj1:lj2), intent(in) :: inarr
    Complex, Dimension(li1:li2,lj1:lj2), intent(out) :: outarr

    stop 'to_fourier_boundary is currently not implemented for mkl'
    
  end Subroutine to_fourier_boundary


  Subroutine to_real_boundary(inarr,outarr)
    Complex, Dimension(li1:li2,lj1:lj2) :: inarr
    Complex, Dimension(li1:li2,lj1:lj2), intent(out) :: outarr

    stop 'to_real_boundary is currently not implemented for mkl'

  end Subroutine to_real_boundary

  subroutine finalize_fourier_boundary

  end subroutine finalize_fourier_boundary

  subroutine finalize_fourier_x_1d 
    INTEGER :: status, thrnum
    
#ifdef WITHOMP
    thrnum=OMP_GET_THREAD_NUM()
#else
    thrnum = 0
#endif    
    Status = DftiFreeDescriptor( fft_kx2x_1d(thrnum)%desc)
    If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
    Status = DftiFreeDescriptor( fft_x2kx_1d(thrnum)%desc)
    If (Status /= 0) Print *, __LINE__, DftiErrorMessage(status)
    
    deallocate(fft_kx2x_1d,fft_x2kx_1d)

  end subroutine finalize_fourier_x_1d

End Module fourier
