#include <complex.h>
#include "cuda_runtime.h"
#include "cufft.h"
#include "cuda_overlap.h"

extern int nj0;
extern int lbg0;
extern double sum_cabs(cufftDoubleComplex *array,int len);
extern double sum_cabs(double _Complex *array,int len);
//extern cudaStream_t stream[nStreams];
extern Streamdata *mystreams[nStreams];
extern unsigned long int allocatedDeviceMemory;

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
static cufftHandle *r2c_plan_y_block, *c2r_plan_y_block;
//static cufftDoubleComplex *dev_cmplxarr;
//static cufftDoubleComplex *dev_cmplxblock_with_da;
static cufftDoubleComplex **dp_temp_with_da;
static double facnny;

extern "C" long int cuda_fourier_get_memory_need_on_device(int a_li0da, int a_ly0da) {
  long int mem_need;

  mem_need = 0L;
  /* dev_cmplxblock_with_da */
  //mem_need = lbg0*li0da*(ly0da/2+1)*sizeof(cufftDoubleComplex);

  /* dp_temp_with_da[iStream] */
  mem_need += nStreams*2*lbg0*a_li0da*(a_ly0da/2+1)*sizeof(cufftDoubleComplex)/nParts;

  return mem_need;
}

/*!
  !	Initialization for fourier routines
  !*/
extern "C" void initialize_fourier_cufft(int a_li0da, int a_ly0da) {
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
  /*const int fft_dims[]={a_ly0da};
    const int fft_dims_red[]={a_ly0da/2+1};*/
  cudaError_t cuda_err;
  int iStream;
  long int mem_free, mem_total;


  li0da=a_li0da;
  ly0da=a_ly0da;
  /*printf("Initializing fourier_cufft with li0da=%u and ly0da=%u.\n",li0da,ly0da);*/

  facnny=1.0/(double)ly0da;

  /* Allocate data storage on the device */

  /* next array is needed for the block oriented FFT */
  /*printf("Allocated so far: %lu. Trying to allocate %u bytes on the device.\n",allocatedDeviceMemory,
    2*lbg0*li0da*(ly0da/2+1)*sizeof(cufftDoubleComplex));*/
  /*cuda_err = cudaMalloc((void**)&dev_cmplxblock_with_da,
			lbg0*li0da*(ly0da/2+1)*sizeof(cufftDoubleComplex));
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory allocation on the device for dev_cmplxblock_with_da.\n%s\n",
	   cudaGetErrorString(cuda_err));
  } else {
    allocatedDeviceMemory += lbg0*li0da*(ly0da/2+1)*sizeof(cufftDoubleComplex);
  }
  */
  /* How much memory on the GPU does the plans need? */
  //cufftPlan1d(&c2r_plan_y_block,ly0da,CUFFT_Z2D,li0da*2*lbg0);

  dp_temp_with_da = (cufftDoubleComplex **)malloc(nStreams*sizeof(cufftDoubleComplex*));
  c2r_plan_y_block = (cufftHandle*) malloc(nStreams*sizeof(cufftHandle));
  r2c_plan_y_block = (cufftHandle*) malloc(nStreams*sizeof(cufftHandle));
  for (iStream=0;iStream<nStreams;iStream++) {
    cudaMemGetInfo((size_t*)&mem_free, (size_t*)&mem_total);
    /*printf("Allocating %lu bytes in stream %u of %u...., free is %lu of %lu\nallocated so far: %lu\n",
	   2*lbg0/nParts*li0da*(ly0da/2+1)*sizeof(cufftDoubleComplex),
	   iStream,nStreams,mem_free,mem_total,allocatedDeviceMemory);*/
    cuda_err = cudaMalloc((void**)&dp_temp_with_da[iStream],
			  2*lbg0*li0da*(ly0da/2+1)*sizeof(cufftDoubleComplex)/nParts);
    if (cuda_err != cudaSuccess) {
      printf("cudaMalloc for dp_temp_with_da[%u] and size of %u gave an error:\n%s",
	     iStream,2*lbg0*li0da*(ly0da/2+1)*sizeof(cufftDoubleComplex)/nParts,
	     cudaGetErrorString(cuda_err));
    } else {
      allocatedDeviceMemory += 2*lbg0*li0da*(ly0da/2+1)*sizeof(cufftDoubleComplex)/nParts;
    }
    //printf("successful.\n");
    //printf("Setting plan for 1D FFT with dim=%u, howmany=%u.\n",ly0da,li0da*2*lbg0/nParts);
    cufftPlan1d(&c2r_plan_y_block[iStream],ly0da,CUFFT_Z2D,li0da*2*lbg0/nParts);
    cufftSetStream(c2r_plan_y_block[iStream],mystreams[iStream]->cudaStream);
    cufftPlan1d(&r2c_plan_y_block[iStream],ly0da,CUFFT_D2Z,li0da*lbg0/nParts);
    cufftSetStream(r2c_plan_y_block[iStream],mystreams[iStream]->cudaStream);
  }

  /*cufftPlan1d(&r2c_plan_y_block,ly0da,CUFFT_D2Z,li0da*lbg0);*/

}

#if 0
/* Copies the inarr to a larger array (larger in width nj0 ->ly0da/2+1).
   Rows 0->nj0 are copied, nj0+1->ly0da/2+1 are set to zero.
   One threadblock per x-y plane and one thread for each point in the x-y plane
   of the larger output array.
   PROBLEM: This easily exceeds 1024 threads per block.
   Therefore the xy-plane is tiled into nTiles, given in gridDim.x*/
__global__ void dev_copy_and_zero_for_dealiasing(cufftDoubleComplex *inarr, int iwidth, int iheight,
						 cufftDoubleComplex *outarr,int owidth, int oheight) {
  /* double _Complex inarr[howmany][iheight][iwidth];
     double _Complex outarr[howmany][oheight][owidth];
     gridDim.x = nTiles; gridDim.y = howmany;
     blockDim.x = owidth; blockDim.y=oheight; 
  */
  int nLinesPerTile = oheight/gridDim.x;
  int iOffset = blockIdx.y*iwidth*iheight + blockIdx.x*nLinesPerTile*iwidth 
    + threadIdx.y*iwidth;
  int ilind;
  int olind = blockIdx.y*owidth*oheight + blockIdx.x*nLinesPerTile*owidth
    + threadIdx.y*owidth+threadIdx.x;

  if (threadIdx.x < iwidth) {
    /* just copy */
    ilind = iOffset+threadIdx.x;
    outarr[olind] = inarr[ilind];
  } else {
    /* set entry to zero, x is real part, y is imaginary part. */
    outarr[olind].x = 0.0;
    outarr[olind].y = 0.0;
  }
}

/* Wrapper for the cuda function dev_copy_and_zero_for_dealiasing. */
void copy_and_zero_for_dealiasing(cufftDoubleComplex *inarr, int iwidth, int iheight,
				  cufftDoubleComplex *outarr,int owidth, int oheight,
				  int howmany, int streamId) {
  int nXYPoints,nLinesPerTile,nTiles;
  dim3 dimBlock;
  dim3 dimGrid;
  cudaError_t cuda_err;  
  
  /*printf("inarr: w=%u, h=%u, out: w=%u, h=%u, howmany = %u\n",iwidth,iheight,
    owidth,oheight,howmany);*/

  if (iheight != oheight) 
    printf("\n----- iheight MUST EQUAL oheight in copy_and_zero_for_dealiasing -----\n\n");

  nXYPoints = owidth*oheight;
  if (nXYPoints<=maxThreadsPerBlock) {
    dimBlock.x=owidth;
    dimBlock.y=oheight;
    dimGrid.x = 1;
    dimGrid.y=howmany;
    dev_copy_and_zero_for_dealiasing<<<dimGrid,dimBlock,
      0,mystreams[streamId]->cudaStream>>>(inarr,iwidth,iheight,
					   dp_temp_with_da[streamId],owidth,oheight);
  } else {
    /* separate the height in several tiles */
    nLinesPerTile = maxThreadsPerBlock/owidth;
    nTiles = oheight/nLinesPerTile;
    if (oheight%nLinesPerTile) nTiles++;
    
    if (oheight%nTiles) {
      printf("oheight (%u) must be divisable by nTiles (%u)\n",oheight,nTiles);
    }
    //printf("nLinesPerTile = %u, nTiles = %u\n",nLinesPerTile,nTiles);
    dimBlock.x = owidth;
    dimBlock.y = oheight/nTiles;
    dimGrid.x  = nTiles;
    dimGrid.y  = howmany;
    dev_copy_and_zero_for_dealiasing<<<dimGrid,dimBlock,
      0,mystreams[streamId]->cudaStream>>>(inarr,iwidth,iheight,
					   dp_temp_with_da[streamId],owidth,oheight);
  }
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) printf("copy_and_zero: %s\n",cudaGetErrorString(cuda_err));

}

/* transform 2-dimensional array inarr(ky,x) to outarr(y,x)
   ! note : y is the first coordinate */
void to_real_y_hp_only_on_device(cufftDoubleComplex *inarr, double *outarr, int howmany, int streamId) {
  /*
    ASSUMPTION:
    inarr[howmany][li0da][nj0] and outarr[howmany][li0da][ly0da]
  */
  /*double _Complex temparr[li0da][ly0da/2+1];*/
  //dim3 dimBlock;
  //dim3 dimGrid;
  //cufftDoubleComplex *temparrblock;
  cudaError_t cuda_err;
  cufftResult cufft_err;
  int klmn,offset,i,j;
  cufftDoubleComplex *hp_temp_with_da;

  /* for dealiasing in y direction, we have to copy the inarr
     into the starting part of temparr, which is a larger
     array. The remainder of temparr is filled with zeros. */
  /*dimBlock.x=ly0da/2+1;
  dimBlock.y=li0da;
  dimGrid.x=howmany;*/
  /*printf("Calling dev_copy_and_zero_for_dealiasing with grid(%u,%u) and block(%u,%u,%u)\n",
    dimGrid.x,dimGrid.y,dimBlock.x,dimBlock.y,dimBlock.z);*/
  /*dev_copy_and_zero_for_dealiasing<<<dimGrid,dimBlock,
    0,mystreams[streamId]->cudaStream>>>(inarr,nj0,
					 li0da,dp_temp_with_da[streamId],
					 ly0da/2+1,li0da,howmany);*/
  copy_and_zero_for_dealiasing(inarr,nj0,li0da,
			       dp_temp_with_da[streamId],ly0da/2+1,li0da,
			       howmany, streamId);
  /*cuda_err = cudaStreamSynchronize(stream[streamId]);*/

#if 0
  /* DEBUG */
  hp_temp_with_da = (cufftDoubleComplex*)malloc(howmany*li0da*(ly0da/2+1)*sizeof(cufftDoubleComplex));
  cudaMemcpy(hp_temp_with_da,dp_temp_with_da[streamId],
	     howmany*li0da*(ly0da/2+1)*sizeof(cufftDoubleComplex),cudaMemcpyDeviceToHost);

  klmn= 0;
  offset = klmn*li0da*(ly0da/2+1);
  for (i=0;i<li0da;i++) {
    for (j=0;j<ly0da/2+1;j++) {
      printf("(%10.3e %10.3e) ",hp_temp_with_da[offset+i*(ly0da/2+1)+j].x,
	     hp_temp_with_da[offset+i*(ly0da/2+1)+j].y);
    }
    printf("\n");
  }
  printf("\n");
  /*for (klmn=0;klmn<howmany;klmn++) {
    printf("klmn= %u: sum_cabs = %17.6e\n",klmn,sum_cabs(&hp_temp_with_da[klmn*li0da*(ly0da/2+1)],
							 li0da*(ly0da/2+1)));
							 }*/
  free(hp_temp_with_da);
  /* END DEBUG */
#endif

  cufft_err = cufftExecZ2D(c2r_plan_y_block[streamId], 
			   (cufftDoubleComplex *)dp_temp_with_da[streamId], 
			   (cufftDoubleReal *)outarr);
  if (cufft_err != CUFFT_SUCCESS) printf("Error with cufftExecZ2D. error_code = %d\n",cufft_err);
  //cuda_err = cudaStreamSynchronize(stream[streamId]);
  //if (cuda_err != cudaSuccess) printf("cufftExecZ2D: %s\n",cudaGetErrorString(cuda_err));
}
#endif

/* transform 2-dimensional array inarr(ky,x) to outarr(y,x)
   ! note : y is the first coordinate */
void to_real_y_only_on_device(cufftDoubleComplex *inarr, double *outarr, int howmany, int streamId) {
  /*
    ASSUMPTION:
    inarr[howmany][li0da][ly0da/2+1] and outarr[howmany][li0da][ly0da]
  */
  /*double _Complex temparr[li0da][ly0da/2+1];*/
  //dim3 dimBlock;
  //dim3 dimGrid;
  //cufftDoubleComplex *temparrblock;
  //cudaError_t cuda_err;
  cufftResult cufft_err;

  cufft_err = cufftExecZ2D(c2r_plan_y_block[streamId], 
			   (cufftDoubleComplex *)inarr, 
			   (cufftDoubleReal *)outarr);
  if (cufft_err != CUFFT_SUCCESS) printf("Error with cufftExecZ2D. error_code = %d\n",cufft_err);
  //cuda_err = cudaStreamSynchronize(stream[streamId]);
  //if (cuda_err != cudaSuccess) printf("cufftExecZ2D: %s\n",cudaGetErrorString(cuda_err));
}

/* just a small kernel which copies the first nj0 entries of the first
   parameter to the second one.
   
   gridDim.x=lbg0
   gridDim.y=nTiles
   blockDim.x=nj0
   blockDim.y=li0da/nTiles
*/
__global__ void copy_only_nondealiased(cufftDoubleComplex *dev_cmplxblock_with_da,int iwidth,int iheight,
				       cufftDoubleComplex *dev_outarr, int owidth,double facnny) {

  int o_index = blockIdx.x*owidth*iheight + blockIdx.y*blockDim.y*owidth + threadIdx.y*owidth + threadIdx.x; 
  int i_index = blockIdx.x*iwidth*iheight + blockIdx.y*blockDim.y*iwidth + threadIdx.y*iwidth + threadIdx.x; 
  dev_outarr[o_index].x = dev_cmplxblock_with_da[i_index].x*facnny;
  dev_outarr[o_index].y = dev_cmplxblock_with_da[i_index].y*facnny;
}

/* We assume oheight=iheight. */
void copy_only_nondealiased_wrapper(cufftDoubleComplex *dev_cmplxblock_with_da,int iwidth, int iheight,
				    cufftDoubleComplex *dev_outarr, int owidth, double facnny, int iStream) {
  dim3 grid,threadblock; //(lbg0,1);
  //dim3 threadblock(li0da,nj0,1);
  cudaError_t cuda_err;
  int dimXYPlane, nTiles; 
  
  dimXYPlane = owidth*iheight;
  nTiles = (dimXYPlane+1023)/1024;

  while ((iheight%nTiles !=0)&&(nTiles<128)) nTiles++;
  if (iheight%nTiles !=0) {
    printf("iheight = %u, cannot be tiled in copy_only_nondealiased_wrapper! Aborting!\n");
    exit(1);
  }/* else {
    printf("We are using %u tiles in copy_only_dealiased.\n",nTiles);
    }*/
  grid.x = lbg0/nParts;
  grid.y = nTiles;
  grid.z = 1;
  threadblock.x=owidth;
  threadblock.y=iheight/nTiles;
  threadblock.z=1;
  copy_only_nondealiased<<<grid,threadblock,0,mystreams[iStream]->cudaStream>>>
    (dev_cmplxblock_with_da,iwidth,iheight,dev_outarr,owidth,facnny);
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) {
    printf("Error with kernel launch.(copy_only_nondealiased)\n");
    printf("%s\n",cudaGetErrorString(cuda_err));
  }
				    }
void to_fourier_y_hp_only_on_device(double *inarr, cufftDoubleComplex *dev_outarr, int iStream) {
  /*
    Real,    Dimension(0:ly0da-1,0:li0da/n_procs_y-1), intent(in)   :: inarr
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1)  , intent(out) :: outarr
    Complex, Dimension(0:ly0da/2,0:li0da/n_procs_y-1)              :: temparr
    ASSUMPTION:
    double inarr[lbg0][li0da][ly0da];
    double _Complex outarr[lbg0][li0da][nj0];
  */
  //double _Complex temparr[lbg0][li0da][ly0da/2+1];

  /*cufftExecD2Z(r2c_plan_y_block[iStream], (cufftDoubleReal *)inarr, 
    (cufftDoubleComplex *)dev_cmplxblock_with_da);*/
#if 0
  double _Complex *tmp_arr;
  int j,klmn;
#endif

  cufftExecD2Z(r2c_plan_y_block[iStream], (cufftDoubleReal *)inarr, 
	       (cufftDoubleComplex *)(mystreams[iStream]->dp_fordeal));
#if 0
  tmp_arr = (double _Complex*)malloc(li0da*(ly0da/2+1)*lbg0/nParts*sizeof(cufftDoubleComplex));
  cudaMemcpy(tmp_arr,mystreams[iStream]->dp_fordeal,
	     li0da*(ly0da/2+1)*lbg0/nParts*sizeof(double _Complex),
	     cudaMemcpyDeviceToHost);
#if 0
  for (klmn=0;klmn<lbg0/nParts;klmn++) {
    printf("klmn=%u, %f\n",klmn, sum_cabs(tmp_arr+klmn*li0da*(ly0da/2+1),li0da*(ly0da/2+1)));
  }

  /*for (j=0;j<ly0da/2+1;j++) {
    printf("(%f %f) ",creal(tmp_arr[j]),cimag(tmp_arr[j]));
  }
  printf("\n");*/
#endif
  printf("after FFT, before dealiasing (y: 1->ly0da/2+1) is %f\n",
	 sum_cabs((cufftDoubleComplex*)tmp_arr,li0da*(ly0da/2+1)*lbg0/nParts));
  free(tmp_arr);
#endif

  copy_only_nondealiased_wrapper((cufftDoubleComplex *)(mystreams[iStream]->dp_fordeal),
				 ly0da/2+1,li0da,
				 dev_outarr,nj0,facnny,iStream);
}

extern "C" void finalize_fourier_cufft() {
  int iStream;
  /* free the GPU memory */
  //cudaFree(dev_cmplxarr);
  //cudaFree(dev_cmplxblock_with_da);

  /* free the FFT plans */
  for (iStream=0;iStream<nStreams;iStream++) {
    cufftDestroy(c2r_plan_y_block[iStream]);
    cudaFree(dp_temp_with_da[iStream]);
    cufftDestroy(r2c_plan_y_block[iStream]);
  }
  free(dp_temp_with_da);
  free(c2r_plan_y_block);
  free(r2c_plan_y_block);
}
