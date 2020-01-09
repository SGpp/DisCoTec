#include "redef.h"

Module fourier
  USE discretization, ONLY: ni0,nj0,nky0,li0,lj0,li1,li2,lj1,lj2,xy_local,n_procs_x,n_procs_y,mype, my_pex
  USE communications, ONLY: MPI_COMPLEX_TYPE, mpi_comm_x
  Use par_in, ONLY: fourier2d
  Implicit None

  public :: initialize_fourier,fft_kx_to_x,fft_ky_to_y,fft_y_to_ky,fft_x_to_kx,&
       fft_ff_to_xy_t,fft_xy_to_ff_t,to_real_y,to_fourier_y,finalize_fourier,&
       initialize_fourier_x_1d,to_fourier_x_1d,to_real_x_1d,finalize_fourier_x_1d,&
       check_fourier,fft_ff_to_xy,fft_xy_to_ff,&
       initialize_fourier_boundary, to_real_boundary,to_fourier_boundary,&
       finalize_fourier_boundary

  private

  Real(Kind=8), Allocatable:: aux_kx2x(:), aux_x2kx(:)
  Real(Kind=8), Allocatable:: aux_kx2x_1d(:), aux_x2kx_1d(:)
  Real(Kind=8), Allocatable:: aux_ky2y(:), aux_y2ky(:)

  Integer:: aux2size = 0
  Integer:: inc_1x, num_x, ub_out_x
  integer:: li0da, ly0da

  REAL::   facnny,facnnx,facnnx_1d,facbound

Contains

  subroutine check_fourier(print_ini_msg)
    logical, intent(in) :: print_ini_msg

    If ((mype==0).and.print_ini_msg) then
       Write(*,"(A)") "Using ESSL."
       If (fourier2d) Write(*,"(A)") "2D fourier analysis only present when using fftw lib, switching to 1D"
    Endif
    fourier2d=.false.

  end subroutine check_fourier

  !
  ! Initialization for Fourier routines
  !
  Subroutine initialize_fourier(a_li0da,a_ly0da)
    INTEGER, INTENT(IN) :: a_li0da, a_ly0da

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

    facnny=1.0/real(ly0da)
    facnnx=1.0/REAL(li0da)

    Call fourier_initaux
    Call fourier_init

  Contains
    Subroutine fourier_init
      Real(Kind=8):: aux2(aux2size)

      Call scrft(1, 0, ly0da/2+1, 0, ly0da,&
           ly0da, li0da/n_procs_y, -1, 1.0,&
           aux_ky2y, Size(aux_ky2y), aux2, aux2size)
      
      Call srcft(1, 0, ly0da, 0, ly0da/2+1,&
           ly0da, li0da/n_procs_y, 1, facnny,&
           aux_y2ky, SIZE(aux_y2ky), aux2, aux2size)
      
      if (xy_local) then
         Call scft(1, 0, 1, li0da,&
              0, inc_1x,1,&
              li0da, num_x, -1, 1.0,&
              aux_kx2x, Size(aux_kx2x), aux2, aux2size)
         
         Call scft(1, 0, inc_1x, 1,&
              0, 1, li0da,&
              li0da, num_x, 1, facnnx,&
              aux_x2kx, Size(aux_x2kx), aux2, aux2size)
      endif

    End Subroutine fourier_init

    Subroutine fourier_initaux
      Character*8::  save2015
      Integer, Parameter::minaux = 32
      Integer::      naux1, naux2
      Real(Kind=8):: tmpaux1(minaux), tmpaux2(minaux)
      External::     ENOTRM

      ! Save ESSL error option table entry
      Call errsav(2015, save2015)
      ! Set error options for error 2015 (don't exit):
      Call errset(2015, 10, -1, 1, ENOTRM, 0)

      naux1 = minaux
      naux2 = minaux
      Call scrft(1, 0, ly0da/2+1, 0, ly0da,&
           ly0da, li0da/n_procs_y, -1, 1.0,&
           tmpaux1, naux1, tmpaux2, naux2)
      Allocate(aux_ky2y(naux1))
      If (aux2size < naux2) aux2size = naux2
      
      naux1 = minaux
      naux2 = minaux
      Call srcft(1, 0, ly0da, 0, ly0da/2+1,&
           ly0da, li0da/n_procs_y, 1, facnny,&
           tmpaux1, naux1, tmpaux2, naux2)
      Allocate(aux_y2ky(naux1))
      If (aux2size < naux2) aux2size = naux2

      if (xy_local) then
         naux1 = minaux
         naux2 = minaux
         Call scft(1, 0, 1, li0da,&
              0, inc_1x , 1,&
              li0da, num_x, -1, 1.0,&
              tmpaux1, naux1, tmpaux2, naux2)
         Allocate(aux_kx2x(naux1))
         If (aux2size < naux2) aux2size = naux2
         
         naux1 = minaux
         naux2 = minaux
         Call scft(1, 0, inc_1x, 1,&
              0, 1, li0da,&
              li0da, num_x, 1, facnnx,&
              tmpaux1, naux1, tmpaux2, naux2)
         Allocate(aux_x2kx(naux1))
         If (aux2size < naux2) aux2size = naux2
      endif

      ! Restore original entry in ESSL error option table:
      Call errstr(2015, save2015)
    End Subroutine fourier_initaux

  End Subroutine initialize_fourier

  subroutine fft_kx_to_x(inarr, outarr)
    complex,dimension(li0da, lj0), Intent(In):: inarr    
    complex,dimension(ub_out_x, li0da),intent(out):: outarr

    real(Kind=8):: aux2(aux2size)
    
    call scft(0, inarr, 1, li0da,&
         outarr, inc_1x,1,&
         li0da, num_x, -1, 1.0,&
         aux_kx2x, size(aux_kx2x), aux2, aux2size)
    
  end subroutine fft_kx_to_x
  
  subroutine fft_ky_to_y(inarr, outarr)
    complex,dimension(0:ly0da/2, 0:li0da/n_procs_y-1),intent(in):: inarr
    real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1),intent(out):: outarr

    real(Kind=8):: aux2(aux2size)
    
    call scrft(0, inarr(0,0), ly0da/2+1, outarr(0,0), ly0da,&
         ly0da, li0da/n_procs_y, -1, 1.0,&
         aux_ky2y, size(aux_ky2y), aux2, aux2size)

  end subroutine fft_ky_to_y

  subroutine fft_y_to_ky(inarr, outarr)
    real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1), intent(in):: inarr
    complex, dimension(0:ly0da/2, 0:li0da/n_procs_y-1), intent(out):: outarr

    real(Kind=8):: aux2(aux2size)

    call srcft(0, inarr(0,0), ly0da, outarr(0,0), ly0da/2+1,&
         ly0da, li0da/n_procs_y, 1, facnny,&
         aux_y2ky, size(aux_y2ky), aux2, aux2size)

  end subroutine fft_y_to_ky


  subroutine fft_x_to_kx(inarr, outarr)
    complex,dimension(ub_out_x,li0da),intent(in):: inarr
    complex,dimension(li0da,lj0),intent(out):: outarr

    real(Kind=8):: aux2(aux2size)

    call scft(0, inarr, inc_1x, 1,&
         outarr, 1, li0da,&
         li0da, num_x, 1, facnnx,&
         aux_x2kx, size(aux_x2kx), aux2, aux2size)
    
  end subroutine fft_x_to_kx

!-------------Two dimensional ffts - These are just dummy routines for compilation
!                                    they are never called

  subroutine fft_ff_to_xy(inarr,outarr)
    complex,dimension(0:ly0da/2,0:li0da-1), intent(in)::inarr
    real,dimension(0:ly0da-1,0:li0da-1), intent(inout)::outarr

    STOP '2d FFT not implemented for ESSL library.'

  end subroutine fft_ff_to_xy

  subroutine fft_xy_to_ff(inarr,outarr)
    real,dimension(0:ly0da-1,0:li0da-1), intent(in)::inarr
    complex,dimension(0:ly0da/2,0:li0da-1), intent(inout)::outarr

    STOP '2d FFT not implemented for ESSL library.'

  end subroutine fft_xy_to_ff

  subroutine fft_ff_to_xy_t(inarr,outarr)
    complex,dimension(0:ly0da/2,0:li0da-1), intent(in)::inarr
    real,dimension(0:ly0da-1,0:li0da-1), intent(inout)::outarr

    STOP '2d FFT not implemented for ESSL library.'

  end subroutine fft_ff_to_xy_t

  subroutine fft_xy_to_ff_t(inarr,outarr)
    real,dimension(0:ly0da-1,0:li0da-1), intent(in)::inarr
    complex,dimension(0:ly0da/2,0:li0da-1), intent(inout)::outarr

    STOP '2d FFT not implemented for ESSL library.'

  end subroutine fft_xy_to_ff_t

!-----------------------------------------------------------------
!-------------routines for global version-------------------------
!-----------------------------------------------------------------


! transform 2-dimensional array inarr(ky,x) to outarr(y,x)
! note : y is the first coordinate
  Subroutine to_real_y(inarr,outarr)
    Complex, Dimension(0:nky0-1, 0:li0da/n_procs_y-1), intent(in)  :: inarr
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1),intent(out) :: outarr
    Complex, Dimension(0:ly0da/2,0:li0da/n_procs_y-1)             :: temparr
    real(Kind=8):: aux2(aux2size)
    
    temparr(0:nky0-1,0:li0da/n_procs_y-1)=inarr
    ! Dealiasing in y direction
    temparr(nky0:,:)=0

    call fft_ky_to_y(temparr,outarr)

  end Subroutine to_real_y

! transform 2-dimensional array inarr(y,x) to outarr(ky,x)
! note : y is the first coordinate
  Subroutine to_fourier_y(inarr,outarr)
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1), intent(in)  :: inarr
    Complex, Dimension(0:nky0-1, 0:li0da/n_procs_y-1),  intent(out) :: outarr
    Complex, Dimension(0:ly0da/2,0:li0da/n_procs_y-1)              :: temparr

!    PERFON('toF_y')

    call fft_y_to_ky(inarr,temparr)
    
    ! Removing last element of temparr for dealising
    outarr=temparr(0:nky0-1,:)
!    PERFOFF

  end Subroutine to_fourier_y

  Subroutine finalize_fourier

    deallocate(aux_ky2y, aux_y2ky)
    IF (xy_local) Deallocate(aux_kx2x, aux_x2kx)
    
  End Subroutine finalize_fourier


  !--------------------------------
  !only needed for local init. cond. or testing
  Subroutine initialize_fourier_x_1d
    
    facnnx_1d = 1.0/REAL(ni0)

    Call fourier_x_1d_initaux
    Call fourier_x_1d_init
    
  Contains
    Subroutine fourier_x_1d_init
      Real(Kind=8):: aux2(aux2size)
      
      Call scft(1, 0, 1, 1,&
           0, 1 ,1,&
           ni0, 1, -1, 1.0,&
           aux_kx2x_1d, Size(aux_kx2x_1d), aux2, aux2size)
      
      Call scft(1, 0, 1, 1,&
           0, 1, 1,&
           ni0, 1, 1, facnnx_1d,&
           aux_x2kx_1d, Size(aux_x2kx_1d), aux2, aux2size)
      
    End Subroutine fourier_x_1d_init

    Subroutine fourier_x_1d_initaux
      Character*8::  save2015
      Integer, Parameter::minaux = 32
      Integer::      naux1, naux2
      Real(Kind=8):: tmpaux1(minaux), tmpaux2(minaux)
      External::     ENOTRM

      ! Save ESSL error option table entry
      Call errsav(2015, save2015)
      ! Set error options for error 2015 (don't exit):
      Call errset(2015, 10, -1, 1, ENOTRM, 0)

      naux1 = minaux
      naux2 = minaux
      Call scft(1, 0, 1, 1,&
           0, 1, 1,&
           ni0, 1, -1, 1.0,&
           tmpaux1, naux1, tmpaux2, naux2)
      Allocate(aux_kx2x_1d(naux1))
      If (aux2size < naux2) aux2size = naux2
      
      naux1 = minaux
      naux2 = minaux
      Call scft(1, 0, 1, 1,&
           0, 1, 1,&
           ni0, 1, 1, facnnx_1d,&
           tmpaux1, naux1, tmpaux2, naux2)
      Allocate(aux_x2kx_1d(naux1))
      If (aux2size < naux2) aux2size = naux2         
      
      ! Restore original entry in ESSL error option table:
      Call errstr(2015, save2015)
    End Subroutine fourier_x_1d_initaux

  End Subroutine initialize_fourier_x_1d
  
  ! transform the 1-dimensional array inarr(x) to outarr(kx)
  Subroutine to_fourier_x_1d(inarr,outarr)
    Complex, Dimension(li1:li2), intent(in) :: inarr
    Complex, Dimension(li1:li2), intent(out) :: outarr
    
    ! Local variables
    COMPLEX, DIMENSION(ni0),target :: local_fullarray,local_Fourierarray
    integer :: ierr
    real(Kind=8) :: aux2(aux2size)

    IF (n_procs_x.GT.1) THEN
       ! first gather the array on one processor, as we do not have a parallelized Fourier transformation
       CALL MPI_Gather(inarr(li1),li0,MPI_COMPLEX_TYPE,&
            & local_fullarray(1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
       IF (my_pex.EQ.0) THEN
          ! do Fourier transform on process 0
          CALL scft(0, local_fullarray(li1), 1, 1,&
               local_Fourierarray(li1), 1, 1,&
               ni0, 1, 1, facnnx_1d,&
               aux_x2kx_1d, SIZE(aux_x2kx_1d), aux2, aux2size)
       END IF
       ! now scatter the result array
       CALL MPI_Scatter(local_Fourierarray(1),li0,MPI_COMPLEX_TYPE,&
            & outarr(li1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
    ELSE
       ! do Fourier transform
       CALL scft(0, inarr(li1), 1, 1,&
            outarr(li1), 1, 1,&
            ni0, 1, 1, facnnx_1d,&
            aux_x2kx_1d, SIZE(aux_x2kx_1d), aux2, aux2size)
    END IF
    
  end Subroutine to_fourier_x_1d

!!! WITHOUT DEALIASING
  ! transform the 1-dimensional array inarr(kx) to outarr(x)
  Subroutine to_real_x_1d(inarr,outarr)
    COMPLEX, DIMENSION(li1:li2), INTENT(IN) :: inarr
    COMPLEX, DIMENSION(li1:li2), INTENT(out) :: outarr

    ! Local variables
    COMPLEX, DIMENSION(ni0),target :: local_fullarray,local_Fourierarray
    integer :: ierr
    real(Kind=8) :: aux2(aux2size)
    
    IF (n_procs_x.GT.1) THEN
       CALL MPI_Gather(inarr(li1),li0,MPI_COMPLEX_TYPE,&
            & local_Fourierarray(1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
       IF (my_pex.EQ.0) THEN
          ! do Fourier transform on process 0
          CALL scft(0, local_fourierarray, 1, 1,&
               local_fullarray, 1 ,1,&
               ni0, 1, -1, 1.0,&
               aux_kx2x_1d, SIZE(aux_kx2x_1d), aux2, aux2size)
       END IF
       ! now scatter the result array
       CALL MPI_Scatter(local_fullarray(1),li0,MPI_COMPLEX_TYPE,&
            & outarr(li1),li0,MPI_COMPLEX_TYPE,0,mpi_comm_x,ierr)
    ELSE
       CALL scft(0, inarr, 1, 1,&
            outarr, 1 ,1,&
            ni0, 1, -1, 1.0,&
            aux_kx2x_1d, SIZE(aux_kx2x_1d), aux2, aux2size)
    END IF
    
  end Subroutine to_real_x_1d
  
  Subroutine finalize_fourier_x_1d
    Deallocate(aux_kx2x_1d,aux_x2kx_1d)
  End Subroutine finalize_fourier_x_1d


  !y global boundary
  subroutine initialize_fourier_boundary
    facbound=1./real(ni0)

    Call fourier_boundary_initaux
    Call fourier_boundary_init
    
  Contains
    Subroutine fourier_boundary_init
      Real(Kind=8):: aux2(aux2size)
      
      Call scft(1, 0, 1, li0,&
           0, 1, li0,&
           li0, lj0, -1, 1.0,&
           aux_kx2x, Size(aux_kx2x), aux2, aux2size)
      
      Call scft(1, 0, 1, li0,&
           0, 1, li0,&
           li0, lj0, 1, facbound,&
           aux_x2kx, Size(aux_x2kx), aux2, aux2size)
    end Subroutine fourier_boundary_init

    Subroutine fourier_boundary_initaux    
      Character*8::  save2015
      Integer, Parameter::minaux = 32
      Integer::      naux1, naux2
      Real(Kind=8):: tmpaux1(minaux), tmpaux2(minaux)
      External::     ENOTRM
      
      ! Save ESSL error option table entry
      Call errsav(2015, save2015)
      ! Set error options for error 2015 (don't exit):
      Call errset(2015, 10, -1, 1, ENOTRM, 0)
      
      naux1 = minaux
      naux2 = minaux
      call scft(1, 0, 1, li0,&
         0, 1, li0,&
         li0, lj0, -1, 1.0,&
         tmpaux1, naux1, tmpaux2, naux2)
      Allocate(aux_kx2x(naux1))
      If (aux2size < naux2) aux2size = naux2
      
      naux1 = minaux
      naux2 = minaux
      call scft(1, 0, 1, li0,&
           0, 1, li0,&
           li0, lj0, 1, facbound,&
           tmpaux1, naux1, tmpaux2, naux2)
      Allocate(aux_x2kx(naux1))
      If (aux2size < naux2) aux2size = naux2
      
      ! Restore original entry in ESSL error option table:
      Call errstr(2015, save2015)
    end Subroutine fourier_boundary_initaux

  end subroutine initialize_fourier_boundary
  
  ! transform 2-dimensional array inarr(ky,x) to outarr(y,x)
  ! note : y is the first coordinate
  Subroutine to_real_boundary(inarr,outarr)
    Complex, Dimension(li0,lj0), intent(in)  :: inarr
    Complex, Dimension(li0,lj0),intent(out) :: outarr
    real(Kind=8):: aux2(aux2size)

    if (n_procs_x.gt.1) &
         stop 'x parallelization not implemented yet in fourier_essl/to_real_boundary'

    call scft(0, inarr, 1, li0,&
         outarr, 1, li0,&
         li0, lj0, -1, 1.0,&
         aux_kx2x, size(aux_kx2x), aux2, aux2size)

  end Subroutine to_real_boundary

  ! transform 2-dimensional array inarr(y,x) to outarr(ky,x)
  ! note : y is the first coordinate
  Subroutine to_fourier_boundary(inarr,outarr)
    Complex, Dimension(li0,lj0), intent(in)  :: inarr
    Complex, Dimension(li0,lj0),  intent(out) :: outarr
    real(Kind=8):: aux2(aux2size)
integer :: i, j

    if (n_procs_x.gt.1) &
         stop 'x parallelization not implemented yet in fourier_essl/to_fourier_boundary'

    call scft(0, inarr, 1, li0,&
         outarr, 1, li0,&
         li0, lj0, 1, facbound,&
         aux_x2kx, size(aux_x2kx), aux2, aux2size)

  end Subroutine to_fourier_boundary

  subroutine finalize_fourier_boundary
    deallocate(aux_x2kx, aux_kx2x)
  end subroutine finalize_fourier_boundary


End Module fourier
