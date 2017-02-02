#include "redef.h"

Module fourier
  USE discretization, ONLY: ni0,nj0,li0,lj0,li1,li2,lj1,lj2,xy_local,y_local,n_procs_x,n_procs_y,mype,my_pex,yx_order
  USE communications, ONLY: mpi_comm_x,MPI_COMPLEX_TYPE
  USE par_in, ONLY: fourier2D

  Implicit None

  Include "fftw3.f"

  public ::  initialize_fourier,fft_kx_to_x,fft_ky_to_y,fft_ff_to_xy,fft_ff_to_xy_t,&
       fft_y_to_ky,fft_x_to_kx,fft_xy_to_ff,fft_xy_to_ff_t,&
       to_real_y,to_real_y_no_extend,to_fourier_y,finalize_fourier,&
       initialize_fourier_x_1d,to_fourier_x_1d,to_real_x_1d,finalize_fourier_x_1d,&
       check_fourier, initialize_fourier_boundary, to_real_boundary,to_fourier_boundary,&
       finalize_fourier_boundary

  private

  Integer(Kind=8), save :: fplan_x, bplan_x, fplan_y, bplan_y, fplan_xy, bplan_xy
  Integer(Kind=8), save :: bplan_x_1d, fplan_x_1d

  Real, save :: facnnx, facnny, facbound
  Integer, save :: inc_1x, num_x, ub_out_x, li0da, ly0da

Contains
  subroutine check_fourier(print_ini_msg)
    logical, intent(in) :: print_ini_msg

    If ((mype==0).and.print_ini_msg) then
       If (fourier2D) then
          Write(*,"(A)") "Using 2D fftw3 transforms"
       Else
          Write(*,"(A)") "Using 1D fftw3 transforms"
       Endif
    Endif

  end subroutine check_fourier
  !
  !	Initialization for fourier routines
  !
  Subroutine initialize_fourier(a_li0da,a_ly0da)
    INTEGER, INTENT(IN) :: a_li0da,a_ly0da
    !
    !	We use the fftw api 'Guru execution of plans' (See docu fftw3).
    !	The following conditions must be met:
    !	- The array size, strides, etcetera are the same
    !	  (since those are set by the plan).
    !	- The input and output arrays are the same (in-place) or different
    !	  (out-of-place) if the plan was originally created to be in-place
    !	  or out-of-place, respectively.
    !	- The alignment of the new input/output arrays is the same as that
    !	  of the input/output arrays when the plan was created, unless the
    !	  plan was created with the FFTW_UNALIGNED flag.
    !
    ! We use out-of-place execution.
    ! Don't use the same temp array for input and output!!!
    !
    Complex,dimension(:,:),allocatable:: tmpin, tmpout, cmplxtmp
    Real,dimension(:,:),allocatable:: realtmp

    !$OMP SINGLE
    li0da=a_li0da
    ly0da=a_ly0da
    allocate(tmpin(li0da, ly0da), tmpout(li0da, ly0da))
    if(yx_order) then
       allocate(realtmp(ly0da, lj0),cmplxtmp(ly0da/2+1,lj0))
    else       
       allocate(realtmp(ly0da, li0da),cmplxtmp(ly0da/2+1,li0da))
    endif

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

    !print*,"Initialize fourier with li0da and ly0da = ",li0da, ly0da
    Call sfftw_plan_many_dft_c2r(fplan_y, 1, ly0da, li0da/n_procs_y,&
         tmpin, 0, 1, ly0da/2+1,&         !input
         tmpout, 0, 1, ly0da,&            !output
         FFTW_PATIENT)
    
    Call sfftw_plan_many_dft_r2c(bplan_y, 1, ly0da, li0da/n_procs_y,&
         tmpin, 0, 1, ly0da,&             !input
         tmpout, 0, 1, ly0da/2+1,&        !output
         FFTW_PATIENT)

    if (xy_local) then       
       if(fourier2D) then    ! Direct 2D FFTs          
          if(yx_order) then
             facnnx = 1.0/Real(lj0)
             Call sfftw_plan_dft_c2r_2d(fplan_xy, ly0da, lj0,&
                  cmplxtmp,&
                  realtmp,&
                  FFTW_PATIENT)          
             Call sfftw_plan_dft_r2c_2d(bplan_xy, ly0da, lj0,&
                  realtmp,&
                  cmplxtmp,&
                  FFTW_PATIENT)
          else
             Call sfftw_plan_dft_c2r_2d(fplan_xy, ly0da, li0da,&
                  cmplxtmp,&
                  realtmp,&
                  FFTW_PATIENT)
             Call sfftw_plan_dft_r2c_2d(bplan_xy, ly0da, li0da,&
                  realtmp,&
                  cmplxtmp,&
                  FFTW_PATIENT)
          endif
       else                  ! 1D FFTs used in transpose split method
          Call sfftw_plan_many_dft(fplan_x, 1, li0da, num_x,&
               tmpin, 0, 1, li0da,&             !input
               tmpout, 0, inc_1x, 1,&           !output
               FFTW_BACKWARD, FFTW_PATIENT)
          
          Call sfftw_plan_many_dft(bplan_x, 1, li0da, num_x,&
               tmpin, 0, inc_1x, 1,&            !input
               tmpout, 0, 1, li0da,&            !output
               FFTW_FORWARD, FFTW_PATIENT)
       end if
    endif

    deallocate(tmpin, tmpout, realtmp, cmplxtmp)    
    !$OMP END SINGLE
  End Subroutine initialize_fourier


  subroutine fft_kx_to_x(inarr, outarr)
    complex,dimension(li0da, lj0), Intent(In):: inarr    
    complex,dimension(ub_out_x, li0da),intent(inout):: outarr

    call sfftw_execute_dft(fplan_x, inarr, outarr)
  End subroutine fft_kx_to_x
  
  subroutine fft_ky_to_y(inarr, outarr)
    complex,dimension(0:ly0da/2, 0:li0da/n_procs_y-1),intent(in):: inarr
    real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1),intent(out):: outarr
    
    call sfftw_execute_dft_c2r(fplan_y, inarr(0,0), outarr(0,0))
    
  end subroutine fft_ky_to_y
  
  subroutine fft_y_to_ky(inarr, outarr)
    real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1), intent(in):: inarr
    complex, dimension(0:ly0da/2, 0:li0da/n_procs_y-1), intent(out):: outarr
    
    Call sfftw_execute_dft_r2c(bplan_y, inarr(0,0), outarr(0,0))
    outarr = outarr*facnnx*facnny
    
  end subroutine fft_y_to_ky


  subroutine fft_x_to_kx(inarr, outarr)
    complex,dimension(ub_out_x,li0da),intent(in):: inarr
    complex,dimension(li0da,lj0),intent(out):: outarr

    Call sfftw_execute_dft(bplan_x, inarr, outarr)

  end subroutine fft_x_to_kx

!-------------Two dimensional ffts

  subroutine fft_ff_to_xy(inarr,outarr)
    complex,dimension(0:ly0da/2,0:li0da-1), intent(in)::inarr
    real,dimension(0:ly0da-1,0:li0da-1), intent(inout)::outarr

    call sfftw_execute_dft_c2r(fplan_xy,inarr(0,0),outarr(0,0))

  end subroutine fft_ff_to_xy

  subroutine fft_xy_to_ff(inarr,outarr)
    real,dimension(0:ly0da-1,0:li0da-1), intent(in)::inarr
    complex,dimension(0:ly0da/2,0:li0da-1), intent(inout)::outarr
!print*, 'in fft_xy_to_ff: inarr', lj0da,li0da
!print*, 'in fft_xy_to_ff: outarr', lj0da/2+1,li0da
    call sfftw_execute_dft_r2c(bplan_xy,inarr(0,0),outarr(0,0))
    outarr = outarr*facnnx*facnny

  end subroutine fft_xy_to_ff


!--------------Two dimensional inplace FFTs for yx ordering

  subroutine fft_ff_to_xy_t(inarr,outarr)
    complex,dimension(0:ly0da/2,0:lj0-1),intent(in)::inarr
    real,dimension(0:ly0da-1,0:lj0-1), intent(inout)::outarr

    call sfftw_execute_dft_c2r(fplan_xy,inarr(0,0),outarr(0,0))

  end subroutine fft_ff_to_xy_t

  subroutine fft_xy_to_ff_t(inarr,outarr)
    real,dimension(0:ly0da-1,0:lj0-1), intent(in)::inarr
    complex,dimension(0:ly0da/2,0:lj0-1), intent(inout)::outarr

    call sfftw_execute_dft_r2c(bplan_xy,inarr(0,0),outarr(0,0))
    outarr = outarr*facnnx*facnny

  end subroutine fft_xy_to_ff_t

!-------------routines for global version-------------------------
  !global y
  subroutine initialize_fourier_boundary
    complex,dimension(:,:),allocatable:: tmpin, tmpout

    allocate(tmpin(li0, lj0), tmpout(li0, lj0))

    call sfftw_plan_many_dft(fplan_x, 1, li0, lj0,&
         tmpin, li0, 1, li0,&             !input
         tmpout, li0, 1, li0,&           !output
         FFTW_BACKWARD, FFTW_PATIENT)
    
    call sfftw_plan_many_dft(bplan_x, 1, li0, lj0,&
         tmpin, li0, 1, li0,&            !input
         tmpout, li0, 1, li0,&            !output
         FFTW_FORWARD, FFTW_PATIENT)
    
    deallocate(tmpin, tmpout)    
    facbound=1./real(ni0)
  end subroutine initialize_fourier_boundary

! transform the 2-dimensional array inarr(kx,ky) to outarr(x,ky)
  Subroutine to_real_boundary(inarr,outarr)
    Complex, Dimension(li0, lj0), intent(in) :: inarr
    Complex, Dimension(li0, lj0), intent(out) :: outarr

    PERFON_I('to_realx')
    call sfftw_execute_dft(fplan_x, inarr, outarr)
    PERFOFF_I

  end Subroutine to_real_boundary

! transform the 2-dimensional array inarr(x,ky) to outarr(kx,ky)
  Subroutine to_fourier_boundary(inarr,outarr)
    Complex, Dimension(li0, lj0), intent(in) :: inarr
    Complex, Dimension(li0, lj0), intent(out) :: outarr

    PERFON_I('to_fourx')

    call sfftw_execute_dft(bplan_x, inarr, outarr)
    outarr=outarr*facbound
    PERFOFF_I

  end Subroutine to_fourier_boundary

  subroutine finalize_fourier_boundary
    call sfftw_destroy_plan(fplan_x)
    call sfftw_destroy_plan(bplan_x)
  end subroutine finalize_fourier_boundary

  !--------------------------

! transform 2-dimensional array inarr(ky,x) to outarr(y,x)
! note : y is the first coordinate
  Subroutine to_real_y(inarr,outarr)
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1)  , intent(in) :: inarr
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1), intent(out) :: outarr
    Complex, Dimension(0:ly0da/2,0:li0da/n_procs_y-1)              :: temparr
#if 0
    integer :: i,j,omp_get_thread_num
#endif
    PERFON('toR_y1')    
    temparr(0:nj0-1,:)=inarr
    ! Dealiasing in y direction
    temparr(nj0:,:)=cmplx(0.0,0.0)
    PERFOFF
#if 0
    !$OMP CRITICAL
    print*,"temparr in to_real_y is ",sum(abs(temparr))
    !print*,"temparr: "
    do j=0,ly0da/2
       do i=0,li0da-1
          !write(*,"(ES10.3,1X)",advance='no') sum(abs(temparr(:,i)))
          write(*,"(F5.1,1X)",advance='no') real(temparr(j,i))
       end do
       write(*,"(A)") ""
    end do
    write(*,"(A)") ""
    !$OMP END CRITICAL
#endif
    PERFON('toR_y2') 
    call sfftw_execute_dft_c2r(fplan_y, temparr(0,0), outarr(0,0))
    PERFOFF
#if 0
    !$OMP CRITICAL
    !print*,"outarr for thread ", omp_get_thread_num()
    print*,"outarr:"
    do j=0,ly0da-1
       do i=0,li0da-1
          !write(*,"(ES10.3,1X)",advance='no') sum(abs(outarr(:,i)))
          write(*,"(ES10.1,1X)",advance='no') outarr(j,i)
       end do
       write(*,"(A)") ""
    end do
    !$OMP END CRITICAL
#endif
    
  end Subroutine to_real_y

! transform 2-dimensional array inarr(ky,x) to outarr(y,x)
! note : y is the first coordinate
  Subroutine to_real_y_no_extend(inarr,outarr)
    Complex, Dimension(0:ly0da/2, 0:li0da/n_procs_y-1), intent(in) :: inarr
    Real,    Dimension(0:ly0da-1, 0:li0da/n_procs_y-1), intent(out) :: outarr

    PERFON('toR_y2') 
    call sfftw_execute_dft_c2r(fplan_y, inarr(0,0), outarr(0,0))
    PERFOFF
    
  end Subroutine to_real_y_no_extend

! transform 2-dimensional array inarr(y,x) to outarr(ky,x)
! note : y is the first coordinate
  Subroutine to_fourier_y(inarr,outarr)
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1), intent(in)   :: inarr
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1)  , intent(out) :: outarr
    Complex, Dimension(0:ly0da/2,0:li0da/n_procs_y-1)              :: temparr

    PERFON('toF_y')
    call sfftw_execute_dft_r2c(bplan_y, inarr(0,0), temparr(0,0))

    !print*,"facnny = ",facnny, ", temparr = ",sum(abs(temparr)),sum(abs(temparr(0:nj0-1,:)))
    ! Removing last element of temparr for dealising, and multiplying by facny
    outarr=temparr(0:nj0-1,:)*facnny
    PERFOFF

  end Subroutine to_fourier_y

  Subroutine finalize_fourier

    !$OMP SINGLE
    Call sfftw_destroy_plan(fplan_y)  
    Call sfftw_destroy_plan(bplan_y)  
    
    if (xy_local) then
       if(fourier2D) then    ! Direct 2D FFTs  
          Call sfftw_destroy_plan(fplan_xy)
          Call sfftw_destroy_plan(bplan_xy)
       else
          Call sfftw_destroy_plan(fplan_x)
          Call sfftw_destroy_plan(bplan_x)
       endif
    endif
    !$OMP END SINGLE
  End Subroutine finalize_fourier


  SUBROUTINE initialize_fourier_x_1d
    COMPLEX, DIMENSION(ni0) :: tmpin, tmpout

    Call sfftw_plan_dft_1d(bplan_x_1d, ni0,tmpin, tmpout,FFTW_FORWARD, FFTW_PATIENT)
    
    !only needed for local init. cond. or testing
    Call sfftw_plan_dft_1d(fplan_x_1d, ni0,tmpin, tmpout,FFTW_BACKWARD, FFTW_PATIENT)
  End Subroutine initialize_fourier_x_1d


  ! transform the 1-dimensional array inarr(x) to outarr(kx)
  Subroutine to_fourier_x_1d(inarr,outarr)
    Complex, Dimension(li1:li2), intent(in) :: inarr
    Complex, Dimension(li1:li2), intent(out) :: outarr
    
    ! Local variables
    COMPLEX, DIMENSION(ni0),target :: local_fullarray,local_Fourierarray
    integer :: ierr

    PERFON('to_fox1d')

    IF (n_procs_x.GT.1) THEN
       ! first gather the array on one processor, as we do not have a parallelized Fourier transformation
       CALL MPI_Gather(inarr(li1),li0,MPI_COMPLEX_TYPE,&
            & local_fullarray(1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
       IF (my_pex.EQ.0) THEN
          ! do Fourier transform on process 0
          CALL sfftw_execute_dft(bplan_x_1d, local_fullarray,local_Fourierarray)
          local_Fourierarray = local_Fourierarray/REAL(ni0)
       END IF
       ! now scatter the result array
       CALL MPI_Scatter(local_Fourierarray(1),li0,MPI_COMPLEX_TYPE,&
            & outarr(li1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
    ELSE
       ! do Fourier transform
       CALL sfftw_execute_dft(bplan_x_1d, inarr,outarr)
       outarr = outarr/real(ni0)
    END IF

    !set unphysical mode to zero
    !if (evenx.EQ.1) outarr(ni0/2) = 0.0
    PERFOFF

  end Subroutine to_fourier_x_1d

!!! WITHOUT DEALIASING
! transform the 1-dimensional array inarr(kx) to outarr(x)
  Subroutine to_real_x_1d(inarr,outarr)
    COMPLEX, DIMENSION(li1:li2),intent(IN) :: inarr
    Complex, Dimension(li1:li2), intent(out) :: outarr

    ! Local variables
    COMPLEX, DIMENSION(ni0),target :: local_fullarray,local_Fourierarray
    integer :: ierr

    PERFON('to_rex1d')
    !set unphysical mode to zero
    !if (evenx.EQ.1) inarr(ni0/2)=0.0
    IF (n_procs_x.GT.1) THEN
       CALL MPI_Gather(inarr(li1),li0,MPI_COMPLEX_TYPE,&
            & local_Fourierarray(1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
       IF (my_pex.EQ.0) THEN
          ! do Fourier transform on process 0
          CALL sfftw_execute_dft(fplan_x_1d,local_Fourierarray,local_fullarray)
       END IF
       ! now scatter the result array
       CALL MPI_Scatter(local_fullarray(1),li0,MPI_COMPLEX_TYPE,&
            & outarr(li1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
    ELSE
       CALL sfftw_execute_dft(fplan_x_1d, inarr, outarr)
    END IF
    PERFOFF

  end Subroutine to_real_x_1d


  subroutine finalize_fourier_x_1d

       Call sfftw_destroy_plan(fplan_x_1d)
       Call sfftw_destroy_plan(bplan_x_1d)

  end subroutine finalize_fourier_x_1d

End Module fourier
