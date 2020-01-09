#include <string.h>
#include <stdlib.h>
#include <complex.h>
#include "fftw3.h"
extern int nj0;
extern int lbg0;
extern double sum_cabs(double _Complex *array,int len);

/*
Module fourier
  USE discretization, ONLY: ni0,nj0,li0,lj0,li1,li2,lj1,lj2,xy_local,y_local,n_procs_x,n_procs_y,mype,my_pex,yx_order
  USE communications, ONLY: mpi_comm_x,MPI_COMPLEX_TYPE
  USE par_in, ONLY: fourier2D
  Implicit None

  Include "fftw3.f"

  public ::  initialize_fourier,fft_kx_to_x,fft_ky_to_y,fft_ff_to_xy,fft_ff_to_xy_t,&
       fft_y_to_ky,fft_x_to_kx,fft_xy_to_ff,fft_xy_to_ff_t,&
       to_real_y,to_fourier_y,finalize_fourier,&
       initialize_fourier_x_1d,to_fourier_x_1d,to_real_x_1d,finalize_fourier_x_1d,&
       check_fourier, initialize_fourier_boundary, to_real_boundary,to_fourier_boundary,&
       finalize_fourier_boundary

  private
*/
/*Integer(Kind=8):: fplan_x, bplan_x, fplan_y, bplan_y, fplan_xy, bplan_xy
  Integer(Kind=8):: bplan_x_1d, fplan_x_1d
  Real:: facnnx, facnny, facbound
  Integer:: inc_1x, num_x, ub_out_x, li0da, ly0da*/

static int li0da, ly0da;
static fftw_plan c2r_plan_y, r2c_plan_y;
static double facnny;

/*!
  !	Initialization for fourier routines
  !*/
void initialize_fourier_hp(int a_li0da, int a_ly0da) {
  /*    !
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
  */
  const int fft_dims[]={a_ly0da};
  const int fft_dims_red[]={a_ly0da/2+1};
  double _Complex *cmplxtmp;
  double *realtmp;

  li0da=a_li0da;
  ly0da=a_ly0da;
  /*printf("Initializing fourier_cufft with li0da=%u and ly0da=%u.\n",li0da,ly0da);*/

  facnny=1.0/(double)ly0da;

  cmplxtmp = (double _Complex*)malloc(li0da*(ly0da/2+1)*sizeof(double _Complex));
  realtmp = (double*)malloc(li0da*ly0da*sizeof(double));
  /*allocate(tmpin(li0da, ly0da), tmpout(li0da, ly0da))
    allocate(realtmp(ly0da, li0da),cmplxtmp(ly0da/2+1,li0da))*/

  c2r_plan_y = fftw_plan_many_dft_c2r(1, fft_dims,li0da,
                                      (fftw_complex *)cmplxtmp, 0,
                                      1,ly0da/2+1,
                                      realtmp, 0,
                                      1, ly0da,
                                      FFTW_MEASURE);
  if (c2r_plan_y == NULL) printf("c2r_plan_y was NOT successfully created.\n");

  r2c_plan_y = fftw_plan_many_dft_r2c(1, fft_dims, li0da,
				      realtmp,0,1,ly0da,
				      (fftw_complex *)cmplxtmp,0,1,ly0da/2+1,
				      FFTW_MEASURE);
  if (r2c_plan_y == NULL) printf("r2c_plan_y was NOT successfully created.\n");

  free(realtmp);
  free(cmplxtmp);
}

/* transform 2-dimensional array inarr(ky,x) to outarr(y,x)
   ! note : y is the first coordinate */
void to_real_y_hp(double _Complex *inarr, double *outarr) {
  /*
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1)  , intent(in) :: inarr
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1), intent(out) :: outarr
    Complex, Dimension(0:ly0da/2,0:li0da/n_procs_y-1)              :: temparr
    ASSUMPTION:
    inarr[li0da][nj0] and outarr[li0da][ly0da]
  */
  static int first_klmn=1;
  double _Complex temparr[li0da][ly0da/2+1];
  int i,j;

  /*     temparr(0:nj0-1, 0:li0da/n_procs_y-1)=inarr 
    ! Dealiasing in y direction
    temparr(nj0:,:)=0
  */
  /*printf("inarr in to_real_y_hp\n");
  for (i=0;i<li0da;i++) {
    for (j=0;j<nj0;j++) {
      printf("%5.2f+%5.2fI ",creal(inarr[i*nj0+j]),cimag(inarr[i*nj0+j]));
    }
    printf("\n");
  }
  */
  /*printf("temparr in to_real_y_hp\n");*/
  for (i=0;i<li0da;i++) {
    memcpy(&temparr[i][0],&inarr[i*nj0],nj0*sizeof(double _Complex));
    memset(&temparr[i][nj0],0,(ly0da/2+1-nj0)*sizeof(double _Complex));
  }

#if 0
  if (first_klmn) {
    for (i=0;i<li0da;i++) {
      for (j=0;j<ly0da/2+1;j++) {
	printf("(%10.3e %10.3e) ",creal(temparr[i][j]),cimag(temparr[i][j]));
      }
      printf("\n");
    }
    printf("\n");
    first_klmn=0;
  }
#endif
  /*printf("with da: klmn = ??, sum_cabs = %17.6e\n",sum_cabs(&temparr[0][0],li0da*(ly0da/2+1)));*/
  

  /*fftw_print_plan(c2r_plan_y);*/
  fftw_execute_dft_c2r(c2r_plan_y, (fftw_complex*)temparr, outarr);
}


  /*! transform 2-dimensional array inarr(y,x) to outarr(ky,x)
    ! note : y is the first coordinate */
void to_fourier_y_hp(double *inarr, double _Complex *outarr) {
  /*
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1), intent(in)   :: inarr
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1)  , intent(out) :: outarr
    Complex, Dimension(0:ly0da/2,0:li0da/n_procs_y-1)              :: temparr
    ASSUMPTION:
    double inarr[li0da][ly0da];
    double _Complex outarr[li0da][nj0];
  */
  double _Complex temparr[li0da][ly0da/2+1];
  int i,j;

  fftw_execute_dft_r2c(r2c_plan_y, inarr, (fftw_complex*)temparr);

  /*! Removing last element of temparr for dealising, and multiplying by facny*/
  for (i=0;i<li0da;i++) {
    for (j=0;j<nj0;j++) {
      outarr[i*nj0+j]=temparr[i][j]*facnny;
    }
  }
}

void finalize_fourier_hp() {
  fftw_destroy_plan(c2r_plan_y);
  fftw_destroy_plan(r2c_plan_y);  
}

/*
  SUBROUTINE initialize_fourier_x_1d
    COMPLEX, DIMENSION(ni0) :: tmpin, tmpout

    Call sfftw_plan_dft_1d(bplan_x_1d, ni0,tmpin, tmpout,FFTW_FORWARD, FFTW_MEASURE)
    
    !only needed for local init. cond. or testing
    Call sfftw_plan_dft_1d(fplan_x_1d, ni0,tmpin, tmpout,FFTW_BACKWARD, FFTW_MEASURE)
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
		    */
