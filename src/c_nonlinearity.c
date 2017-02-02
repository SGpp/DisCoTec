#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <assert.h>

extern _Bool dealiasing;
extern int li0,lj0,li0da,li1da,ly0da,nj0;
/* this routine will not compile because dealiasing is removed in rev. 4697*/
extern int lbida2, ubida2
extern double ve_max[2];

extern void da_interpolate_to_fine_grid_x(double _Complex inarr[],double _Complex tmp_arr1[]);
extern void da_interpolate_to_fine_grid_x_orig(double _Complex inarr[],double _Complex tmp_arr1[]);
extern void da_interpolate_to_fine_grid_x_from_c(double _Complex *, double _Complex*);
extern void da_filter_from_fine_grid_x_from_c(double _Complex *, double _Complex *);
extern void to_fourier_y_hp(double *inarr, double _Complex *outarr);
extern void to_real_y_hp(double _Complex *inarr, double *outarr);

/* Prototypes */
void c_gto_real(double _Complex *inarr, double *outarr, int howmany);
void c_transpose_cmplx_old(int n1, int n2, double _Complex *in_matrix,  
		       double _Complex *transposed_matrix);
void c_transpose_cmplx(int nRows, int nCols, double _Complex *mat,int iwidth,int soff,  
		       double _Complex *tmat,int owidth,int doff);
double c_maxval(double *,int);
double sum_abs(double *array,int len);
double sum_cabs(double _Complex *array,int len);
double sum_cabs2(double _Complex *array,int len);

/* some array which are temporarliy needed in the calculation of the
   nonlinearity. */
double *vexy_re, *dgdxy_re;
double *nonlin1_re;
static int lbg0;

void c_initialize_nonlinearity_df(int a_lbg0) {
  lbg0 = a_lbg0;
  /* Allocate the arrays which are needed lateron on the heap, to save
     stack space, which is limited on many systems. */
  vexy_re  = (double*)malloc(lbg0*2*li0da*ly0da*sizeof(double));
  dgdxy_re = (double*)malloc(lbg0*2*li0da*ly0da*sizeof(double));
  nonlin1_re = (double*)malloc(lbg0*li0da*ly0da*sizeof(double));

  initialize_fourier_hp(li0da,ly0da);
}

void c_finalize_nonlinearity_df(void) {

  finalize_fourier_hp();
  /* free the arrays */
  free(vexy_re);
  free(dgdxy_re);
  free(nonlin1_re);
}
/* Computes the nonlinearity for the (x) nonlocal version
  !! We assume in this module that we do not have a parallelization in 
  !! x, y or z direction. This is mainly in preparation to GPGPU port.
  !! But it could first be used with OpenMP tasking, which also needs
  !! many independent tasks.
  !!\todo Get rid of transposition using striding? make some arrays for arakawa allocatable
  !!\todo Check whether the nonlinearity prefactor is important for the CFL criterion
*/
void c_calc_nonlinearity_df(double _Complex *gy_chi,   
			    double _Complex *g2d,      
			    double _Complex *vexy,  
			    double _Complex *dgdxy, 
			    double _Complex *localrhs, 
			    double *cptr_pnl_1d,
			    _Bool first) {
  /*
    double _Complex gy_chi[lbg0][lj0][li0]
    double _Complex g2d[lbg0][lj0][li0]      
    double _Complex vexy[lbg0][2][lj0][li0]
    double _Complex dgdxy[lbg0][2][lj0][li0] 
    double _Complex localrhs[lbg0][lj0][li0] 
  */

  /* Local variables */
  /*Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,2,1:lbg0) ::  vexy_re, dgdxy_re*/
  /*double vexy_re[lbg0][2][li0da][ly0da], dgdxy_re[lbg0][2][li0da][ly0da];*/
  /*double *vexy_re, *dgdxy_re;*/
  /*Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,1:lbg0) ::  nonlin1_re, nonlin2_re, nonlin3_re*/
  /*double nonlin1_re[lbg0][li0da][ly0da];*/
  /*Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr*/
  double _Complex tmp_arr[lbg0][li0da][nj0];
  /*    Complex, Dimension(lbida2:ubida2,lj1:lj2,1:lbg0) :: nl_tmp1*/
  double _Complex nl_tmp1[lbg0][lj0][li0];
  /*Complex, Dimension(li1:li2,lj1:lj2,1:lbg0) :: nonlin */
  double _Complex nonlin[lbg0][lj0][li0];

  /* Local variables */
  int i,j,klmn,ijblocksize,offset0,offset1;
  double temp_max;

#if 0
  int full_blocksize = li0*lj0*lbg0;
  printf("-------------- START of c_calc_...:\n localrhs = %13.6e, vexy = %13.6e, dgdxy = %13.6e\n",
	 sum_cabs((double _Complex *)localrhs,full_blocksize),
	 sum_cabs((double _Complex *)vexy,2*full_blocksize),
	 sum_cabs((double _Complex *)dgdxy,2*full_blocksize));
#endif
  /*printf("Running with lbg0 = %u, cptr_pnl_1d = %f\n",lbg0,sum_abs(cptr_pnl_1d,li0));*/
  /* anti-aliasing and y fourier transform*/
  c_gto_real(vexy,vexy_re,2*lbg0);
  c_gto_real(dgdxy,dgdxy_re, 2*lbg0);

  /*printf("\tAFTER FFT: vexy_re = %13.6e, dgdxy_re = %13.6e\n",
	 sum_abs((double *)vexy_re,2*li0da*ly0da*lbg0),
	 sum_abs((double *)dgdxy_re,2*li0da*ly0da*lbg0));
  */

  /* get max ExB velocity*/
  if (first) {
    /*ve_x_max_loc=maxval(vexy_re(:,:,2,:))
      ve_y_max_loc=maxval(vexy_re(:,:,1,:))
      ve_max(1)=max(ve_max(1),ve_x_max_loc)        
      ve_max(2)=max(ve_max(2),ve_y_max_loc)
    */
    ijblocksize = li0da*ly0da;
    /*printf("Determining ve_max, old = %f %f, ijblocksize = %d\n",ve_max[0],ve_max[1],ijblocksize);*/
    for (klmn=0;klmn<lbg0;klmn++) {
      temp_max = c_maxval(&vexy_re[(2*klmn+1)*ijblocksize],ijblocksize);
      ve_max[0] = fmax(temp_max,ve_max[0]);

      temp_max = c_maxval(&vexy_re[2*klmn*ijblocksize],ijblocksize);
      ve_max[1] = fmax(temp_max,ve_max[1]);
    }
    /*printf("ve_max = %f %f\n",ve_max[0],ve_max[1]);*/
  }

  /* compute the 'standard' nonlinear term, which is also used in the Arakawa representation*/
  /*nonlin1_re = -vexy_re(:,:,1,:)*dgdxy_re(:,:,2,:) + vexy_re(:,:,2,:)*dgdxy_re(:,:,1,:);
   There is no array syntax in C, so we have either to use explicit index syntax or use
   a BLAS function.*/
  for (klmn=0;klmn<lbg0;klmn++) {
    offset0 = 2*klmn*li0da*ly0da;
    offset1 = (2*klmn+1)*li0da*ly0da;
    for (i=0;i<li0da;i++) {
      for (j=0;j<ly0da;j++) {
	/*nonlin1_re[klmn][i][j] = -vexy_re[klmn][0][i][j]*dgdxy_re[klmn][1][i][j]
	  + vexy_re[klmn][1][i][j]*dgdxy_re[klmn][0][i][j];*/
	nonlin1_re[klmn*li0da*ly0da+i*ly0da+j] = -vexy_re[offset0+i*ly0da+j]*dgdxy_re[offset1+i*ly0da+j]
	  + vexy_re[offset1+i*ly0da+j]*dgdxy_re[offset0+i*ly0da+j];
      }
    }
  }
  /*printf("nonlin1_re = %13.6E \n",sum_abs((double*)nonlin1_re,li0da*ly0da*lbg0));*/

  for (klmn=0;klmn<lbg0;klmn++) {
    /* transform back to fourier space*/
    /*to_fourier_y(nonlin1_re(:,:,klmn),tmp_arr(:,:,klmn));*/
    /*to_fourier_y_hp(&nonlin1_re[klmn][0][0],&tmp_arr[klmn][0][0]);*/
    to_fourier_y_hp(&nonlin1_re[klmn*li0da*ly0da],&tmp_arr[klmn][0][0]);
    /*printf("klmn = %5u, sum_cabs(tmp_arr) = %10.3e\n",klmn,sum_cabs2((double _Complex*)&tmp_arr[klmn][0][0],li0da*nj0));*/
  }
  
  /*printf("tmp_arr = %13.6e\n",sum_cabs2((double _Complex *)tmp_arr,
    li0da*nj0*lbg0));*/

  /* Transpose and remove zeros in y for dealiasing*/
  if (dealiasing) {
    /*nl_tmp1=cmplx(0.0);*/
    memset( nl_tmp1, '\0', lbg0*lj0*li0*sizeof(double _Complex) );

    for (klmn=0;klmn<lbg0;klmn++) {
      /*Call transpose_cmplx(nj0, li0da, tmp_arr(:,:,klmn), 0, nl_tmp1(li1da:li2da,lj1:lj2,klmn), 0) */
      /*c_transpose_cmplx(nj0, li0da, &tmp_arr[klmn][0][0], &nl_tmp1[klmn][0][li1da-lbida2]);*/
      c_transpose_cmplx(li0da, nj0, &tmp_arr[klmn][0][0], nj0,0,
			&nl_tmp1[klmn][0][0],li0,li1da-lbida2);
      /*printf("%10.3e ",sum_cabs2(&nl_tmp1[klmn][0][li1da-lbida2],li0da));*/
    }
    /*printf("nl_tmp1 = %13.6e\n",sum_cabs2((double _Complex *)nl_tmp1,
      li0*nj0*lbg0));*/

    /*filter_from_fine_grid_x(nl_tmp1(lbida2:ubida2,lj1:lj2,1:lbg0), nonlin);*/
    da_filter_from_fine_grid_x_from_c(&nl_tmp1[0][0][0], &nonlin[0][0][0]);

    /*for (klmn=0;klmn<lbg0;klmn++) {
      printf("%3d: back in C, nonlin = %17.10E\n",klmn+1,sum_cabs2(&nonlin[klmn][0][0],li0*lj0));
      }*/

  } else {
    for (klmn=0;klmn<lbg0;klmn++) {
      /*transpose_cmplx(nj0, li0da, tmp_arr(:,:,klmn), 0, nonlin(:,:,klmn), 0);*/
      c_transpose_cmplx(li0da, nj0, &tmp_arr[klmn][0][0], nj0,0,
			&nonlin[klmn][0][0],li0,0);
    }
  }
  /*printf("nonlin = %13.6e\n",sum_cabs2(&nonlin[0][0][0],lbg0*li0*nj0));*/
    
  for (klmn=0;klmn<lbg0;klmn++) {
    for (j=0;j<lj0;j++) {
      for (i=0;i<li0;i++) {
	localrhs[klmn*lj0*li0+j*li0+i] += cptr_pnl_1d[i]*nonlin[klmn][j][i];
      }
    }
  }
  /*printf("after addition: localrhs = %f\n",
    sum_cabs((double _Complex *)localrhs,lbg0*li0*lj0));*/
}


/*!> Go to real space */
void c_gto_real(double _Complex *inarr, double *outarr,int howmany) {
  /*
    Complex, Dimension(li1:li2,lj1:lj2,1:howmany),Intent(in) :: inarr
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,1:howmany), target, Intent(out) :: outarr
     these Fortran array specification translates into the following C specification
    double _Complex inarr[howmany][lj0][li0];
    double outarr[howmany][li0da][ly0da];
  */

  /* Local variables */
  /*Complex, Dimension(li1da:li2da,lj1:lj2,1:howmany) :: tmp_arr1*/
  /*double _Complex tmp_arr1[howmany][lj0][li0da];*/
  double _Complex *tmp_arr1;
  /*    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr2
	double _Complex tmp_arr2[li0da][nj0];*/
  double _Complex *tmp_arr2;
  int klmn,number_of_lbg0_blocks,i_block, offset,i,j;
    
  tmp_arr1 = (double _Complex*)malloc(howmany*lj0*li0da*sizeof(double _Complex));
  tmp_arr2 = (double _Complex*)malloc(li0da*nj0*sizeof(double _Complex));

  number_of_lbg0_blocks = howmany/lbg0;
  /*printf("number_of_lbg0_blocks = %d\n",number_of_lbg0_blocks);*/
  if (number_of_lbg0_blocks==0) {
#if 0
    if (dealiasing) {
      /* Interpolation for x anti-aliasing */
      /*CALL interpolate_to_fine_grid_x_orig(inarr(:,:,1),tmp_arr1(:,:,1))*/
      da_interpolate_to_fine_grid_x_orig(&inarr[0],&tmp_arr1[0][0][0]);
    } else {
      /*tmp_arr1=inarr;*/
      memcpy(tmp_arr1,inarr,lj0*li0*sizeof(double _Complex));
    }
    /* Transpose x-y*/
    /*Call transpose_cmplx(li0da, nj0, tmp_arr1(:,:,1), 0,tmp_arr2 , 0) */
    c_transpose_cmplx_old(li0da, nj0, &tmp_arr1[0][0][0],&tmp_arr2[0][0]);
    /* Fourier transfrom in y (include dealiasing step)*/
    /*Call to_real_y(tmp_arr2,outarr(:,:,1)) */
    to_real_y_hp(&tmp_arr2[0][0],&outarr[0]);
#endif
  } else {
    for (i_block=0;i_block<number_of_lbg0_blocks;i_block++) {
      /*printf("i_block = %d\n",i_block);*/
      if (dealiasing) {
	/* Interpolation for x anti-aliasing */
	/*CALL interpolate_to_fine_grid_x(inarr(:,:,(i_block-1)*lbg0+1:i_block*lbg0),
	  tmp_arr1(:,:,(i_block-1)*lbg0+1:i_block*lbg0));*/
	/*da_interpolate_to_fine_grid_x(&inarr[i_block*lbg0*li0*lj0],&tmp_arr1[i_block*lbg0*li0da*lj0]);*/
	da_interpolate_to_fine_grid_x_from_c(&inarr[i_block*lbg0*li0*lj0],&tmp_arr1[i_block*lbg0*li0da*lj0]);
      } else {
	/*tmp_arr1=inarr;*/
	for (klmn=0;klmn<lbg0;klmn++) {
	  for (j=0;j<lj0;j++) {
	    for (i=0;i<li0;i++) {
	      tmp_arr1[(i_block*lbg0+klmn)*lj0*li0da+j*li0da+i] = inarr[(i_block*lbg0+klmn)*lj0*li0+j*li0+i];
	    }
	  }
	}
	/*memcpy((void*)&tmp_arr1[i_block*lbg0][0][0],
	  (void*)&inarr[i_block*lbg0*li0*lj0],lbg0*li0*lj0*sizeof(double _Complex));*/
      }
      
      /* Transpose x-y*/
      for (klmn=0;klmn<lbg0;klmn++) {
	/*printf("Before transpose klmn = %u: sum_cabs(tmp_arr1) = %f\n",
	  i_block*lbg0+klmn,sum_cabs(&tmp_arr1[(i_block*lbg0+klmn)*lj0*li0da],li0da*lj0));*/
	/*Call transpose_cmplx(li0da, nj0, tmp_arr1(:,:,(i_block-1)*lbg0+klmn), 0,tmp_arr2 , 0); 
	  c_transpose_cmplx_old(li0da, nj0, &tmp_arr1[(i_block*lbg0+klmn)*lj0*li0da],tmp_arr2);*/
	c_transpose_cmplx(nj0,li0da, &tmp_arr1[(i_block*lbg0+klmn)*lj0*li0da],
			  li0da,0,tmp_arr2,nj0,0);
	
	/* debug output */
	/*printf("after transpose: klmn = %u, tmp_arr2=%f\n",
	  i_block*lbg0+klmn,sum_cabs(tmp_arr2,li0da*nj0));*/
	/* Fourier transfrom in y (include dealiasing step)*/
	/*p_out => outarr(:,:,(i_block-1)*lbg0+klmn)
	  Call to_real_y(tmp_arr2,p_out) !outarr(:,:,(i_block-1)*lbg0+klmn))*/
	offset = (i_block*lbg0+klmn)*li0da*ly0da;
	/*printf("Calling to_real_y_hp with sum_cabs(tmp_arr2) = %f\n",sum_cabs(tmp_arr2,li0da*nj0));*/
	to_real_y_hp(tmp_arr2,&outarr[offset]);
	/*printf("klmn=%u, sum_abs(outarr[%u]) = %17.6e\n",klmn,offset,sum_abs(&outarr[offset],li0da*ly0da));*/
      }
    }
  }
  free(tmp_arr1);
  free(tmp_arr2);
}

/*
  !----------------------------------------------------------------------
  !> Transpose the Matrix in_matrix and leave Result in transposed_matrix
  !---------------------------------------------------------------------- */
void c_transpose_cmplx_old(int n1, int n2, double _Complex *in_matrix,  
		       double _Complex *transposed_matrix) {
  /*
    Integer, Intent(in):: n1, n2, ssoff, ddoff
    Complex, Intent(in):: in_matrix(0:, 0:)
    Complex, Intent(out):: transposed_matrix(0:, 0:)
    Assumption:
    in_matrix[n2][n1], transposed_matrix[n1][n2];
  */

  /* Local variables */
  int i,j;
  /* We assume n_procs_y = 1 */

  for (j=0;j<n2;j++) {
    for (i=0;i<n1;i++) {
      transposed_matrix[i*n2+j] = in_matrix[j*n1+i];
    }
  }
}

void c_transpose_cmplx(int nRows, int nCols, 
		       double _Complex *mat, int iwidth, int soff,
		       double _Complex *tmat, int owidth, int doff) {

  /* Local variables */
  int dRow,dCol;
  /* We assume n_procs_y = 1 */

  for (dRow=0;dRow<nCols;dRow++) {
    for (dCol=0;dCol<nRows;dCol++) {
      tmat[dRow*owidth+doff+dCol] = mat[dCol*iwidth+soff+dRow];
    }
  }
}

double c_maxval(double *array,int len) {
  int i;
  double maximum;

  maximum=array[0];
  for (i=1;i<len;i++) {
    maximum = fmax(maximum,array[i]);
  }
  return maximum;
}

/* calculate the sum of the absolute values of an array */
double sum_abs(double *array,int len) {
  double sum;
  int i;

  sum = 0.0;
  for (i=0;i<len;i++) {
    sum += fabs(array[i]);
  }
  return sum;
}

double sum_cabs(double _Complex *array,int len) {
  double sum;
  int i;

  sum = 0.0;
  for (i=0;i<len;i++) {
    sum += cabs(array[i]);
  }
  return sum;
}

double sum_cabs2(double _Complex *array,int len) {
  double sum;
  int i;

  sum = 0.0;
  for (i=0;i<len;i++) {
    sum += array[i]*conj(array[i]);
  }
  return sum;
}
