#include <string.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <cuda_runtime.h>
#include <cufft.h>
#ifdef WITHPERF
#include "perflib.h"
#endif
#include "redef.h"
#include "cuda_overlap.h"
#include "cuda_kernels.h"

//#define MEASURE_EVENTS
#undef MEASURE_EVENTS
#ifdef MEASURE_EVENTS
#define EVENT_RECORD(ev,st) cudaEventRecord(ev,st)
#else
#define EVENT_RECORD(ev,st)
#endif

/* some macros for the calculation of the linear index */
#define LIND3(x,y,z,ny,nz) x*ny*nz+y*nz+z
#define LIND4(x1,x2,x3,x4,nx2,nx3,nx4) x1*nx2*nx3*nx4+x2*nx3*nx4+x3*nx4+x4

extern _Bool dealiasing;
extern int li0,lj0,li0da,li1da,ly0da,nj0;
/*extern int lbida2, ubida2; */
extern double ve_max[2];

extern "C" void da_interpolate_to_fine_grid_x(double _Complex inarr[],double _Complex tmp_arr1[]);
extern void da_interpolate_to_fine_grid_x_orig(double _Complex inarr[],double _Complex tmp_arr1[]);
extern void da_filter_from_fine_grid_x(double _Complex nl_tmp1[], double _Complex nonlin[]);
extern void to_fourier_y_hp_only_on_device(double *inarr, cufftDoubleComplex *outarr, int iStream);
//extern void to_real_y_hp_only_on_device(cufftDoubleComplex *inarr, double *outarr, int howmany, int streamId);
extern void to_real_y_only_on_device(cufftDoubleComplex *inarr, double *outarr, int howmany, int streamId);
/*extern "C" double sum_abs(double *array, int len);*/
extern "C" void c_transpose_cmplx(int n1, int n2, double _Complex *in_matrix,  
		       double _Complex *transposed_matrix);
extern "C" double c_maxval(double *,int);
extern "C" void initialize_fourier_cufft(int a_li0da, int a_ly0da);
extern "C" void finalize_fourier_cufft();
extern "C" long int cuda_fourier_get_memory_need_on_device(int a_li0da,int a_ly0da);

int lbg0;
double sum_abs(double *array, int len);
double sum_cabs(cufftDoubleComplex *array,int len);
double sum_cabs(double _Complex *array,int len);
void show_xy_arr(double _Complex *array,int nrows, int ncols);
void show_xy_arr(cufftDoubleComplex *array,int nrows, int ncols);
void error_handling_with_synchronize(const char *kernelname);
void cuda_OverlapTransferAndFFT(const double _Complex *hf_inarr, double *df_outarr, int dims[3]);
void cuda_OverlapFFTAndTransfer(double *dev_nonlin1_re, double _Complex *localrhs);

Streamdata *mystreams[nStreams];

cufftDoubleComplex *dev_temparrblock, *dev_vexy, *dev_dgdxy, *dev_nonlin, *dev_localrhs;
//double _Complex *nonlin;
double *dev_pnl;
cufftDoubleReal *dev_vexy_re, *dev_dgdxy_re, *dev_nonlin1_re;
unsigned long int allocatedDeviceMemory=0L;
static cudaEvent_t start_event, end_event, event2;
#ifdef MEASURE_EVENTS
static cudaEvent_t start_transfer[nStreams],end_transfer[nStreams],after_copy_and_zero[nStreams],after_transpose[nStreams],after_to_real[nStreams];
static float stime_transfer[nStreams],stime_copy_and_zero[nStreams],stime_transpose[nStreams],stime_to_real[nStreams];
#endif


extern "C" void cuda_initialize_nonlinearity_df(int a_lbg0,
						double *cptr_pnl_1d) {
  lbg0 = a_lbg0;
  cudaMalloc((void**)&dev_pnl,li0da*sizeof(double));
  cudaMemcpy(dev_pnl,cptr_pnl_1d,li0da*sizeof(double),cudaMemcpyHostToDevice);
}

extern "C" void cuda_finalize_nonlinearity_df() {
  cudaFree(dev_pnl);
}

/* we need some temporary arrays on the GPU, which are allocated
   in advance.
   cufftDoubleComplex dev_cmplxblock[2*lbg0][lj0][li0];
   cufftDoubleComplex dev_temparrblock[2*lbg0][li0da][nj0];
   cufftDoubleReal dev_realblock[2*lbg0][li0da][ly0da];
*/

extern "C" int cuda_get_nearest_blocksize(int test_blocksize) {
  return (test_blocksize/nParts)*nParts;
}

extern "C" int cuda_get_device_count(void) {
  int number_of_devices;
  cudaGetDeviceCount(&number_of_devices);
  return number_of_devices;
}

extern "C" void cuda_set_device(int device) {
  int number_of_devices;
  //long int free_mem, total_mem;

  cudaGetDeviceCount(&number_of_devices);
  //printf("We have %u GPU devices.\n",number_of_devices);
  if (device>=number_of_devices) {
    device = 2;
  }
  cudaSetDevice(device);

  //cudaMemGetInfo((size_t*)&free_mem,(size_t*)&total_mem);
  /*printf("Total available memory on device %u is %lu, free is %lu.\n",
    device, total_mem, free_mem);*/
}
  
extern "C" void cuda_register_array(double _Complex **arr,int arr_size) {
  cudaError_t cuda_err;
  //printf("registering arr at address %p, size=%u bytes.\n",(void*)arr,arr_size*sizeof(double _Complex));
  cuda_err = cudaHostRegister((void*)arr,arr_size*sizeof(double _Complex),cudaHostRegisterPortable);
  if (cuda_err!=cudaSuccess) {
    printf("Error with cudaHostRegister: %s\n",cudaGetErrorString(cuda_err));
  }/* else {
    printf("Registered %lu bytes at address %p.\n",arr_size*sizeof(double _Complex),arr);
    }*/
}

extern "C" void cuda_unregister_array(double _Complex **arr) {
  cudaError_t cuda_err;
  //printf("unregistering arr at address %p\n",(void*)arr);
  cuda_err = cudaHostUnregister((void*)arr);
  if (cuda_err!=cudaSuccess) {
    printf("Error with cudaHostUnregister: %s\n",cudaGetErrorString(cuda_err));
  }
}

extern "C" long int cuda_get_memory_need_on_device() {
  long int mem_need;
  //long int mem_sav;

  mem_need = 0L;
  /*mem_need = 2*2*lbg0*lj0*li0*sizeof(cufftDoubleComplex);*/
  /* dev_vexy_re, dev_dgdxy_re */
  mem_need += 2*2*lbg0*ly0da*li0da*sizeof(cufftDoubleReal);
  /* dev_nonlin1_re */
  mem_need += lbg0*ly0da*li0da*sizeof(cufftDoubleReal);
  /* dev_temparrblock */
  mem_need += 2*lbg0*li0da*nj0*sizeof(cufftDoubleComplex);
  /* dev_nonlin, dev_localrhs */
  mem_need += 2*lbg0*li0da*nj0*sizeof(cufftDoubleComplex);

  /* in the overlap case, we have per stream the following memory needs */
  /* dp_data */
  mem_need += nStreams*(2*lbg0*li0*lj0)/nParts*sizeof(cufftDoubleComplex);
  /* dp_temp */
  mem_need += nStreams*(2*lbg0*li0*(ly0da/2+1))/nParts*sizeof(cufftDoubleComplex);
  /* dp_fordeal */
  mem_need += nStreams*(2*lbg0*li0*(ly0da/2+1))/nParts*sizeof(cufftDoubleComplex);

  /*mem_sav = mem_need;
    printf("mem_need without fourier: %lu\n",mem_need);*/
  mem_need += cuda_fourier_get_memory_need_on_device(li0da,ly0da);
  /*printf("mem_need for fourier: %lu, all gathered : %lu\n",
    mem_need-mem_sav,mem_need);*/
  /* for the reduction for the maximum */
  mem_need += (2*NTILES*lbg0+2)*sizeof(cufftDoubleReal);

  /* the prefactor for the nonlinearity */
  mem_need += li0da*sizeof(double);
  return mem_need;
}

extern "C" long int cuda_get_free_memory_on_device() {
  long int mem_free, mem_total;

  cudaMemGetInfo((size_t*)&mem_free, (size_t*)&mem_total);
  return mem_free;
}

extern "C" void cuda_allocate_on_device(int device) {
  cudaError_t cuda_err;
  int passed=1;
  int number_of_devices, thisDevice, iStream;
  struct cudaDeviceProp deviceProperties;
  //long int mem_free,mem_free_start, mem_free_before_fourier;

  allocatedDeviceMemory = 0L;
  cudaGetDeviceCount(&number_of_devices);
  //printf("We have %u GPU devices.\n",number_of_devices);
  if (device<number_of_devices) {
    cudaSetDevice(device);
  } else {
    cudaSetDevice(0);
  }
  cudaGetDevice(&thisDevice);
  cudaGetDeviceProperties(&deviceProperties,thisDevice);
  //printf("We are really using now device %u, name is %s.\n",thisDevice,deviceProperties.name);

  /*mem_free_start=cuda_get_free_memory_on_device();
    printf("Free device memory before allocating anything: %lu.\n",mem_free_start);*/
  /* we also need the real arrays */
  cuda_err = cudaMalloc((void**)&dev_vexy_re,
			2*lbg0*ly0da*li0da*sizeof(cufftDoubleReal));
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory allocation on the device for dev_vexy_re.\n");
    passed=0;
  } else {
    allocatedDeviceMemory += 2*lbg0*ly0da*li0da*sizeof(cufftDoubleReal);
  }
  cuda_err = cudaMalloc((void**)&dev_dgdxy_re,
			2*lbg0*ly0da*li0da*sizeof(cufftDoubleReal));
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory allocation on the device for dev_dgdxy_re.\n");
    passed=0;
  } else {
    allocatedDeviceMemory += 2*lbg0*ly0da*li0da*sizeof(cufftDoubleReal);
  }

  cuda_err = cudaMalloc((void**)&dev_nonlin1_re,
			lbg0*ly0da*li0da*sizeof(cufftDoubleReal));
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory allocation on the device for dev_nonlin1_re.\n");
    passed=0;
  } else {
    allocatedDeviceMemory += lbg0*ly0da*li0da*sizeof(cufftDoubleReal);
  }

  cuda_err = cudaMalloc((void**)&dev_temparrblock,
			2*lbg0*li0da*nj0*sizeof(cufftDoubleComplex));
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory allocation on the device for dev_temparrblock.\n");
    passed=0;
  } else {
    allocatedDeviceMemory += 2*lbg0*li0da*nj0*sizeof(cufftDoubleComplex);
  }

  cuda_err = cudaMalloc((void**)&dev_nonlin,
			lbg0*li0da*nj0*sizeof(cufftDoubleComplex));
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory allocation on the device for dev_nonlin.\n");
    passed=0;
  } else {
    allocatedDeviceMemory += lbg0*li0da*nj0*sizeof(cufftDoubleComplex);
  }

  cuda_err = cudaMalloc((void**)&dev_localrhs,
			lbg0*li0da*nj0*sizeof(cufftDoubleComplex));
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory allocation on the device for dev_localrhs.\n");
    passed=0;
  } else {
    allocatedDeviceMemory += lbg0*li0da*nj0*sizeof(cufftDoubleComplex);
  }
  
  /*printf("allocatedDeviceMemory so far is %lu\n",allocatedDeviceMemory);*/
  /*cuda_err = cudaHostAlloc((void**)&nonlin,lbg0*li0da*nj0*sizeof(double _Complex),cudaHostAllocDefault);
  if (cuda_err != cudaSuccess) {
    printf("Error with pinned memory allocation on host for nonlin.\n");
    passed=0;
  }
  */
  if (!passed) printf("Error with allocation of device memory. Used %lu bytes.\n",allocatedDeviceMemory);

  /* Create all streams for later use of events and asynchronous execution. */
  for (iStream=0;iStream<nStreams;iStream++) {
    //printf("Creating Stream %u of %u....",iStream,nStreams);
    mystreams[iStream] = new Streamdata(li0,lj0,ly0da,2*lbg0);
    //printf("Successful.\n");
  }
  /*mem_free_before_fourier = cuda_get_free_memory_on_device();
  printf("Before fourier init: memory free = %lu. allocatedDeviceMemory = %lu, usedMemory = %lu\n",
  mem_free_before_fourier,allocatedDeviceMemory,mem_free_start-mem_free_before_fourier);*/
  initialize_fourier_cufft(li0da,ly0da);
  /*mem_free = cuda_get_free_memory_on_device();
  printf("After fourier init: memory free = %lu. \
allocatedDeviceMemory = %lu, \
usedMemory by fourier = %lu\n",
mem_free,allocatedDeviceMemory,mem_free_before_fourier-mem_free);*/
  //printf("End of cuda_allocate_on_device.\n");

  cudaEventCreate(&start_event);
  cudaEventCreate(&end_event);
  cudaEventCreate(&event2);

#ifdef MEASURE_EVENTS
  for (iStream=0;iStream<nStreams;iStream++) {
    cudaEventCreate(&(start_transfer[iStream]));
    cudaEventCreate(&(end_transfer[iStream]));
    cudaEventCreate(&(after_copy_and_zero[iStream]));
    cudaEventCreate(&(after_transpose[iStream]));
    cudaEventCreate(&(after_to_real[iStream]));
    stime_transfer[iStream]=0.0;
    stime_copy_and_zero[iStream]=0.0;
    stime_transpose[iStream]=0.0;
    stime_to_real[iStream]=0.0;
  }
#endif
}

extern "C" void cuda_free_on_device() {
  cudaError_t cuda_err;
  int passed=1, iStream;

  cuda_err = cudaFree(dev_vexy_re);
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory free on the device for dev_vexy_re.\n");
    passed=0;
  }
  cuda_err = cudaFree(dev_dgdxy_re);
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory free on the device for dev_dgdxy_re.\n");
    passed=0;
  }
  cuda_err = cudaFree(dev_nonlin1_re);
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory free on the device for dev_nonlin1_re.\n");
    passed=0;
  }

  cuda_err = cudaFree(dev_temparrblock);
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory free on the device for dev_temparrblock.\n");
    passed=0;
  }

  cuda_err = cudaFree(dev_nonlin);
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory free on the device for dev_nonlin.\n");
    passed=0;
  }

  cuda_err = cudaFree(dev_localrhs);
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory free on the device for dev_nonlin.\n");
    passed=0;
  }

  /*cuda_err = cudaFreeHost(nonlin);
  if (cuda_err!=cudaSuccess) {
    printf("Error with memory free on the host for nonlin.\n");
    passed=0;
  }
  */

  finalize_fourier_cufft();
  if (!passed) printf("NOT all device memory successfully freed.\n");

  for (iStream=0;iStream<nStreams;iStream++) {
    delete mystreams[iStream];
    //cudaStreamDestroy(stream[iStream]);
  }
  cudaEventDestroy(start_event);
  cudaEventDestroy(event2);
  cudaEventDestroy(end_event);

#ifdef MEASURE_EVENTS
  for (iStream=0;iStream<nStreams;iStream++) {
    cudaEventDestroy(start_transfer[iStream]);
    cudaEventDestroy(end_transfer[iStream]);
    cudaEventDestroy(after_copy_and_zero[iStream]);
    cudaEventDestroy(after_transpose[iStream]);
    cudaEventDestroy(after_to_real[iStream]);
    printf("Stream %u\ntransfer: %f\ncopy_and_zero: %f\ntranspose: %f\nto_real: %f\n----------\n",iStream,
	 stime_transfer[iStream],stime_copy_and_zero[iStream], stime_transpose[iStream], stime_to_real[iStream]);
  }
#endif
}

#include "cuda_kernels.cu"
#include "reduction_kernel.cu"

/*  This routine is now wholly running on the GPU. We start with
    a memcpy of the necessary input parameters to the GPU and then
    do all (transpose,FFT,nonlin calculation, y-dealiasing, back-FFT)
    on the GPU and only back-memcpy at the end. */

extern "C" 
void cuda_calc_nonlinearity_df(const double _Complex *gy_chi, 
			       const double _Complex *g_block,
			       const double _Complex *vexy,
			       const double _Complex *dgdxy,
			       double _Complex *localrhs,
			       _Bool first) {
  /*double _Complex gy_chi[lbg0][lj0][li0], 
    double _Complex g_block[lbg0][lj0][lbi:ubi], 
    double _Complex vexy[lbg0][2][lj0][li0], 
    double _Complex dgdxy[lbg0][2][lj0][li0], 
    double _Complex localrhs[lbg0][lj0][li0]
  */
  /*Complex, Dimension(li1:li2,lj1:lj2,1:lbg0),Intent(inout) :: gy_chi
    Complex, Dimension(lbi:ubi,lj1:lj2,1:lbg0),Intent(in) :: g_block
    Complex, Dimension(li1:li2,lj1:lj2,2,1:lbg0),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:lbg0),Intent(inout) :: localrhs  
    Logical, Intent(in):: first
  */

  /* Local variables */
  /*Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,2,1:lbg0) ::  vexy_re, dgdxy_re*/
  //double vexy_re[lbg0][2][li0da][ly0da], dgdxy_re[lbg0][2][li0da][ly0da];
#if 0
  double *vexy_re, *dgdxy_re;
  /* debugging */
  void *nonlin1_re;
  /* end debugging */
#endif
  /*Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,1:lbg0) ::  nonlin1_re, nonlin2_re, nonlin3_re*/
  //double nonlin1_re[lbg0][li0da][ly0da];
  /*Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr*/
  //double _Complex tmp_arr[lbg0][li0da][nj0];
  /*    Complex, Dimension(lbida2:ubida2,lj1:lj2,1:lbg0) :: nl_tmp1*/
  /*double _Complex nl_tmp1[lbg0][lj0][ubida2-lbida2+1];*/
  /*Complex, Dimension(li1:li2,lj1:lj2,1:lbg0) :: nonlin */
  //double _Complex nonlin[lbg0][lj0][li0];
  //double _Complex *nonlin;

  /* Local variables */
  int shared_mem;
  cufftDoubleReal *dev_temp_max0, *dev_temp_max1,*dev_max_per_block;
  double temp_max;
  cudaError_t cuda_err;
  dim3 grid, threadblock;
  int array_dims[3];
#if 0
  float time_of_calc_nonlin,time_of_mult_pre, time_of_pure_cuda;
#endif

  //cudaEventRecord(start_event,0);
#if 0
  /* DEBUGGING THE parameters*/
  int full_blocksize = li0*lj0*lbg0;
  printf("-------------- START of cuda_calc_...:\n localrhs = %13.6e, vexy = %13.6e, dgdxy = %13.6e\n",
    sum_cabs((cufftDoubleComplex*)localrhs,full_blocksize),
    sum_cabs((cufftDoubleComplex*)vexy,2*full_blocksize),
    sum_cabs((cufftDoubleComplex*)dgdxy,2*full_blocksize));
  /* end debug */
#endif
  /* ----------------------------------------
     Transfer to the GPU 
     ----------------------------------------*/
  //cudaDeviceSynchronize();
  C_PERFON("deal_FFT",8);
  array_dims[0]=lj0;
  array_dims[1]=li0;
  array_dims[2]=2*lbg0;
  /*printf("array_dims = %u,%u,%u\n",array_dims[0],array_dims[1], array_dims[2]);*/
  cuda_OverlapTransferAndFFT(vexy,dev_vexy_re,array_dims);
  cuda_OverlapTransferAndFFT(dgdxy,dev_dgdxy_re,array_dims);
  /* We have to make sure, that the fft are completed and the arrays are available. */
  cudaDeviceSynchronize();
  C_PERFOFF();
#if 0
  /* DEBUGGING */
  cudaDeviceSynchronize();
  vexy_re = (double*)malloc(li0da*ly0da*2*lbg0*sizeof(double));
  dgdxy_re = (double*)malloc(li0da*ly0da*2*lbg0*sizeof(double));
  /* transfer back to CPU the transformed block */
  cudaMemcpy(vexy_re,dev_vexy_re,
	     li0da*ly0da*2*lbg0*sizeof(double),
	     cudaMemcpyDeviceToHost);
  cudaMemcpy(dgdxy_re,dev_dgdxy_re,
	     li0da*ly0da*2*lbg0*sizeof(double),
	     cudaMemcpyDeviceToHost);

  printf("\tAFTER FFT: vexy_re = %13.6e, dgdxy_re = %13.6e\n",
	 sum_abs((double *)vexy_re,2*li0da*ly0da*lbg0),
	 sum_abs((double *)dgdxy_re,2*li0da*ly0da*lbg0));

#if 0
  int i,j;
  printf("First part, first block:\n");
  for (i=0;i<li0da;i++) {
    for (j=0;j<ly0da;j++) {
      printf("%10.6f ",vexy_re[i*ly0da+j]);
    }
    printf("\n");
  }
  printf("\nSecond part, first block:\n");
  for (i=0;i<li0da;i++) {
    for (j=0;j<ly0da;j++) {
      printf("%10.6f ",vexy_re[li0da*ly0da+i*ly0da+j]);
    }
    printf("\n");
  }
  printf("\n");
#endif
  free(vexy_re);
  free(dgdxy_re);
  /* END DEBUGGING */
#endif

  // get max ExB velocity
  if (first) {
    /*ve_x_max_loc=maxval(vexy_re(:,:,2,:))
      ve_y_max_loc=maxval(vexy_re(:,:,1,:))
      ve_max(1)=max(ve_max(1),ve_x_max_loc)        
      ve_max(2)=max(ve_max(2),ve_y_max_loc)
    */
    C_PERFON("ve_max",6);
    cudaMalloc((void**)&dev_max_per_block,2*NTILES*lbg0*sizeof(cufftDoubleReal));
    cudaMalloc((void**)&dev_temp_max0,sizeof(cufftDoubleReal));
    cudaMalloc((void**)&dev_temp_max1,sizeof(cufftDoubleReal));

    /*grid.x=lbg0;
    threadblock.x=li0da;
    shared_mem = (li0da<64 ? 64 : li0da)*2*sizeof(double);*/
    /*printf("Calling maxval_per_block with grid(%u,%u) and block(%u,%u,%u), shared_mem = %u\n",
      grid.x,grid.y,threadblock.x,threadblock.y,threadblock.z,shared_mem);*/
    //cuda_maxval_per_block_old<<<grid,threadblock,shared_mem>>>(dev_vexy_re,dev_max_per_block,li0da,ly0da);

    grid.x=lbg0;
    grid.y=NTILES;
    threadblock.x=ly0da;
    threadblock.y=1;
    shared_mem = ((ly0da<32) ? 64 : 2*ly0da)*sizeof(double);
    cuda_maxval_per_block<<<grid,threadblock,shared_mem>>>(dev_vexy_re,dev_max_per_block,li0da,ly0da);
    error_handling_with_synchronize("cuda_maxval_per_block");
    /* now we have the maxima for each block in dev_max_per_block.
       This has to be further reduced. */

    /* Debug output */
    /*
    double *max_per_block ,max0,max1;
    max_per_block =(double*)malloc(2*NTILES*lbg0*sizeof(double));
    cudaMemcpy(max_per_block,dev_max_per_block,2*NTILES*lbg0*sizeof(double),cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    printf("First part:\n");
    max0 = 0.0;
    for (i=0;i<NTILES*lbg0;i++) {
      printf("%10.6f ",max_per_block[i]);
      max0 = (max_per_block[i]>max0) ? max_per_block[i] : max0;
    }
    //printf("\n -- max0 = %f\nSecond part:\n",max0);
    printf("\nSecond part:\n");
    max1=0.0;
    for (i=0;i<NTILES*lbg0;i++) {
      printf("%10.6f ",max_per_block[NTILES*lbg0+i]);
      max1 = (max_per_block[NTILES*lbg0+i]>max1) ? max_per_block[NTILES*lbg0+i] : max1;
    }
    //printf("\n -- max1 = %f\n\n",max1);
    printf("\n");
    free(max_per_block);
    */
    /* end of Debug output */

    /* Further reduction */
    reduce<double>(NTILES*lbg0,6,dev_max_per_block,dev_temp_max0,REDUCTION_MAX,0);
    error_handling_with_synchronize("reduce<double>0");

    reduce<double>(NTILES*lbg0,6,dev_max_per_block,dev_temp_max1,REDUCTION_MAX,NTILES*lbg0);
    error_handling_with_synchronize("reduce<double>1");


    /* Freeing of dev_temp_max! */
    /* copying to the host */
    cuda_err = cudaMemcpy(&temp_max,dev_temp_max0,sizeof(double),cudaMemcpyDeviceToHost);
    //cudaDeviceSynchronize();
    if (cuda_err!=cudaSuccess) {
      printf("Error with Memcpy. (dev_temp_max0 -> temp_max)\n");
      printf("%s\n",cudaGetErrorString(cuda_err));
    } 

    //printf("temp_max0 = %f, ve_max[1] = %f\n",temp_max,ve_max[1]);
    /* max of vexy_re(0) is y component, max of vexy_re(1) is x component */
    ve_max[1] = fmax(temp_max,ve_max[1]);


    cuda_err = cudaMemcpy(&temp_max,dev_temp_max1,sizeof(double),cudaMemcpyDeviceToHost);
    //cudaDeviceSynchronize();
    if (cuda_err!=cudaSuccess) {
      printf("Error with Memcpy. (dev_temp_max1 -> temp_max)\n");
      printf("%s\n",cudaGetErrorString(cuda_err));
    } 
    
    //printf("temp_max1 = %f, ve_max[0] = %f\n",temp_max,ve_max[0]);
    ve_max[0] = fmax(temp_max,ve_max[0]);

    cudaFree(dev_max_per_block);
    cudaFree(dev_temp_max0);
    cudaFree(dev_temp_max1);

    //printf("ve_max = %f %f \n",ve_max[0],ve_max[1]);
    C_PERFOFF();
  }

  /* make sure, that all streams are synchronized with the host */
  //cudaDeviceSynchronize();
  // compute the 'standard' nonlinear term, which is also used in the Arakawa representation
  /*nonlin1_re = -vexy_re(:,:,1,:)*dgdxy_re(:,:,2,:) + vexy_re(:,:,2,:)*dgdxy_re(:,:,1,:);
   There is no array syntax in C, so we have either to use explicit index syntax or use
   a BLAS function.*/


  /* While the nonlinearity is calculated, we copy the old localrhs to the GPU. */
  cudaMemcpyAsync(dev_localrhs,localrhs,li0*lj0*lbg0*sizeof(double _Complex),cudaMemcpyHostToDevice,mystreams[0]->cudaStream);

  C_PERFON("nonlin",6);
  C_PERFON("nl_kernl",8);
  /* We call the comp_stand_nonlin kernel. */
  comp_stand_nonlin_wrapper(dev_vexy_re,dev_dgdxy_re,dev_nonlin1_re,li0da,ly0da,nStreams-1);
  cudaDeviceSynchronize();
  C_PERFOFF();

#if 0
  /* debug */
  nonlin1_re = malloc(li0da*ly0da*lbg0*sizeof(double));
  /* transfer back to CPU the transformed block */
  cudaMemcpy(nonlin1_re,dev_nonlin1_re,
	     li0da*ly0da*lbg0*sizeof(double),
	     cudaMemcpyDeviceToHost);
  printf("nonlin1_re = %15.6e \n",sum_abs((double*)nonlin1_re,li0da*ly0da*lbg0));
  free(nonlin1_re);
  /* end debug */
#endif

  //C_PERFON("ovrlp_2",7);
  cuda_OverlapFFTAndTransfer(dev_nonlin1_re,localrhs);
  cudaDeviceSynchronize();
  //C_PERFOFF();
  C_PERFOFF();
  //cudaEventRecord(event2,0);
  /*cudaDeviceSynchronize();

  int lij0 = li0*lj0;
  for (klmn=0;klmn<lbg0;klmn++) {
    for (j=0;j<lj0;j++) {
      for (i=0;i<li0;i++) {
	//localrhs[klmn*lij0+j*li0+i] += cptr_pnl_1d[i]*nonlin[klmn*lij0+j*li0+i];
	localrhs[klmn*lij0+j*li0+i] += nonlin[klmn*lij0+j*li0+i];
      }
    }
  }
  */
  /*printf("after addition: localrhs = %f\n",
    sum_cabs((double _Complex *)localrhs,lbg0*li0*lj0));*/
  //cudaEventRecord(end_event,0);
#if 0
  cudaEventSynchronize(end_event);
  cudaEventElapsedTime(&time_of_calc_nonlin,start_event,end_event);
  cudaEventElapsedTime(&time_of_mult_pre,event2,end_event);
  cudaEventElapsedTime(&time_of_pure_cuda,start_event,event2);
  printf("cuda nonlin time = %f ms, t(mult_pre) = %f ms, t(pure_cuda) = %f ms\n",time_of_calc_nonlin,time_of_mult_pre,time_of_pure_cuda);
#endif
}

void cuda_OverlapFFTAndTransfer(double *dev_nonlin1_re, double _Complex *localrhs) {
  int iStream, iPart;
  int nXYPlanesPerPart, dimXYPlane;
  dim3 grid, threadblock;
#if 0
  double _Complex *tmp_arr;
  int j,klmn;
#endif

  nXYPlanesPerPart = lbg0/nParts;
  dimXYPlane = li0da*ly0da;

  iStream=0;
  for (iPart=0;iPart<nParts;iPart++) {

    to_fourier_y_hp_only_on_device(dev_nonlin1_re+iPart*nXYPlanesPerPart*dimXYPlane,
				   (cufftDoubleComplex*)(mystreams[iStream]->dp_temp),iStream);

#if 0
    cudaDeviceSynchronize();
    tmp_arr = (double _Complex*)malloc(li0da*nj0*nXYPlanesPerPart*sizeof(cufftDoubleComplex));
    cudaMemcpy(tmp_arr,mystreams[iStream]->dp_temp,
	li0da*nj0*nXYPlanesPerPart*sizeof(double _Complex),
	cudaMemcpyDeviceToHost);

    for (klmn=0;klmn<nXYPlanesPerPart;klmn++) {
      printf("klmn=%u, %f\n",klmn, sum_cabs(tmp_arr+klmn*li0da*nj0,li0da*nj0));
    }
    /*for (j=0;j<nj0;j++) {
      printf("(%f %f) ",creal(tmp_arr[j]),cimag(tmp_arr[j]));
      }
      printf("\n");*/
    
    printf("tmp_arr after to_fourier_y_hp is %f\n",sum_cabs(tmp_arr,
							    li0da*nj0*nXYPlanesPerPart));
    free(tmp_arr);
#endif
    
    // Transpose and remove zeros in y for dealiasing
    transpose_wrapper((cufftDoubleComplex*)(mystreams[iStream]->dp_data),
		      (cufftDoubleComplex*)(mystreams[iStream]->dp_temp), nj0,li0da, 
		      nXYPlanesPerPart,iStream);

#if 0
    /* Begin DEBUG */
    cudaDeviceSynchronize();
    tmp_arr = (double _Complex*)malloc(li0da*nj0*nXYPlanesPerPart*sizeof(cufftDoubleComplex));
    cudaMemcpy(tmp_arr,mystreams[iStream]->dp_data,
	li0da*nj0*nXYPlanesPerPart*sizeof(double _Complex),
	cudaMemcpyDeviceToHost);
    
    printf("after backtranspose, dp_data is %f\n",sum_cabs(tmp_arr,
							   li0da*nj0*nXYPlanesPerPart));
    free(tmp_arr);
    /* End DEBUG */
#endif

    /* The transposed nonlinearity is now in mystreams[iStream]->dp_data. We multiply it with 
       the prefactor dev_pnl. */
    grid.x=nXYPlanesPerPart;
    grid.y = nj0; grid.z=1;
    threadblock.x=li0da; threadblock.y=1; threadblock.z=1;
    copy_with_pnl<<<grid,threadblock,0,mystreams[iStream]->cudaStream>>>
      ((cufftDoubleComplex*)mystreams[iStream]->dp_data,
       (cufftDoubleComplex*)dev_localrhs+iPart*nXYPlanesPerPart*li0da*nj0,
       /*(cufftDoubleComplex*)mystreams[iStream]->dp_temp,*/
       (cufftDoubleReal*)dev_pnl,li0da,nj0);
    //cudaStreamSynchronize(mystreams[iStream]->cudaStream);
#if 0
    /* Begin DEBUG */
    cudaDeviceSynchronize();
    tmp_arr = (double _Complex*)malloc(li0da*nj0*nXYPlanesPerPart*sizeof(cufftDoubleComplex));
    cudaMemcpy(tmp_arr,
	       (cufftDoubleComplex*)dev_localrhs+iPart*nXYPlanesPerPart*li0da*nj0,
	       /*mystreams[iStream]->dp_temp,*/
	       li0da*nj0*nXYPlanesPerPart*sizeof(double _Complex),
	       cudaMemcpyDeviceToHost);
    
    printf("after copy_with_pnl, dev_localrhs is %f\n",sum_cabs(tmp_arr,
								li0da*nj0*nXYPlanesPerPart));
    free(tmp_arr);
    /* End DEBUG */
#endif


    cudaMemcpyAsync(localrhs + iPart*nXYPlanesPerPart*li0da*nj0,
		    dev_localrhs + iPart*nXYPlanesPerPart*li0da*nj0,
		    /*mystreams[iStream]->dp_temp,*/
		    li0da*nj0*nXYPlanesPerPart*sizeof(double _Complex),
		    cudaMemcpyDeviceToHost,mystreams[iStream]->cudaStream);
    /*cudaMemcpy(nonlin,dev_nonlin,
      li0da*nj0*lbg0*sizeof(double _Complex),
      cudaMemcpyDeviceToHost);*/
    iStream = (++iStream)%nStreams;
  }
  //cudaDeviceSynchronize();
#if 0
  cudaDeviceSynchronize();
  printf("localrhs = %15.6e\n",sum_cabs((cufftDoubleComplex*)localrhs,lbg0*li0da*nj0));
#endif
}


void cuda_OverlapTransferAndFFT(const double _Complex *hf_inarr, double *df_outarr, int dims[3]) {

  int dimXYPlane=dims[0]*dims[1];
  int nXYPlanesPerPart, iStream,iPart;
  dim3 grid,threadblock;
  cudaError_t cuda_err;
#ifdef MEASURE_EVENTS
  int nextStream;
  float time_transfer,time_copy_and_zero,time_transpose,time_to_real;
#endif

#if 0
  void *hp_temp;
  int klmn, offset,i,j;
#endif
  
  //printf("Starting cuda_OverlapTransferAndFFT.\n");
  if (dims[2]%nParts != 0) {
    printf("dims[2]=%u not divisable by nParts = %u\n",dims[2],nParts);
    exit(1);
  }
  nXYPlanesPerPart = dims[2]/nParts;
  /*printf("Working on %u parts, total number of xy-planes %u => each part contains %u XY-planes.\n \
    nStreams = %u\n",
    nParts,dims[2],nXYPlanesPerPart,nStreams);*/

  iStream = 0;
  for (iPart=0;iPart<nParts;iPart++) {
    /* copy a part to the device in stream iStream */
    EVENT_RECORD(start_transfer[iStream],mystreams[iStream]->cudaStream);
    cudaMemcpyAsync(mystreams[iStream]->dp_data,&hf_inarr[iPart*nXYPlanesPerPart*dimXYPlane],
		    nXYPlanesPerPart*dimXYPlane*sizeof(double _Complex),
		    cudaMemcpyHostToDevice,mystreams[iStream]->cudaStream);
    EVENT_RECORD(end_transfer[iStream],mystreams[iStream]->cudaStream);
    /*cuda_err = cudaStreamSynchronize(mystreams[iStream]->cudaStream);
      if (cuda_err != cudaSuccess) printf("overlap Memcpy: %s\n",cudaGetErrorString(cuda_err));*/

#if 0
    /* debug */
    printf("sum_cabs(hp_inarr) = %f, dimXYPlane = %u, element = %f+%fi\n",
	   sum_cabs((cufftDoubleComplex*)&hf_inarr[iPart*nXYPlanesPerPart*dimXYPlane],
		    nXYPlanesPerPart*dimXYPlane), dimXYPlane, creal(hf_inarr[0]),cimag(hf_inarr[0]));
    /* end debug */
#endif

    /* instead of first transposing and then afterwards do the copy and zero for dealiasing in y,
       we do the copy and zero first for better memory coalescing. 
       dims[0]=lj0, dims[1]=li0 */

    copy_and_zero_for_dealiasing_wrapper((cufftDoubleComplex*)mystreams[iStream]->dp_data,
					 dims[1],dims[0],nXYPlanesPerPart,
					 (cufftDoubleComplex*)mystreams[iStream]->dp_fordeal,
					 dims[1],ly0da/2+1,iStream);

    /*grid.x=ly0da/2+1;         grid.y=nXYPlanesPerPart;
    threadblock.x=dims[1];    threadblock.y=1;

    dev_copy_and_zero_for_dealiasing_new<<<grid,threadblock,0,mystreams[iStream]->cudaStream>>>
      ((cufftDoubleComplex*)mystreams[iStream]->dp_data,dims[1],dims[0],
       (cufftDoubleComplex*)mystreams[iStream]->dp_fordeal,dims[1],ly0da/2+1);
    */
    //cuda_err = cudaStreamSynchronize(mystreams[iStream]->cudaStream);
    cuda_err = cudaGetLastError();
    if (cuda_err != cudaSuccess) printf("%u: after dev_copy_and_zero_for_dealiasing_new: %s\n",
					iStream,cudaGetErrorString(cuda_err));
    EVENT_RECORD(after_copy_and_zero[iStream],mystreams[iStream]->cudaStream);

#if 0
    /* debug */
    hp_temp = malloc(nXYPlanesPerPart*dims[1]*(ly0da/2+1)*sizeof(double _Complex));
    cudaMemcpy(hp_temp,mystreams[iStream]->dp_fordeal,
	       nXYPlanesPerPart*dims[1]*(ly0da/2+1)*sizeof(double _Complex),
	       cudaMemcpyDeviceToHost);
    printf("total sum_cabs(dp_fordeal[%u]) = %f\n",iStream,
	   sum_cabs((cufftDoubleComplex*)hp_temp,nXYPlanesPerPart*dims[1]*(ly0da/2+1)));
    free(hp_temp);
    /* end debug */
#endif

    /* compute FFT in iStream */
    transpose_wrapper((cufftDoubleComplex*)mystreams[iStream]->dp_temp, 
		      (cufftDoubleComplex*)mystreams[iStream]->dp_fordeal, 
		      dims[1],ly0da/2+1,nXYPlanesPerPart,iStream);
    EVENT_RECORD(after_transpose[iStream],mystreams[iStream]->cudaStream);

#if 0
    /* debug */
    hp_temp = malloc(nXYPlanesPerPart*dimXYPlane*sizeof(double _Complex));
    cudaMemcpy(hp_temp,mystreams[iStream]->dp_temp,nXYPlanesPerPart*dimXYPlane*sizeof(double _Complex),
	       cudaMemcpyDeviceToHost);
    printf("total sum_abs(dp_temp[%u]) = %f\n",iStream,
	   sum_cabs((double _Complex*)hp_temp,nXYPlanesPerPart*dimXYPlane));

    klmn= 0;
    offset = klmn*dimXYPlane;
    for (i=0;i<dims[1];i++) {
      for (j=0;j<dims[0];j++) {
	printf("(%10.3e %10.3e) ",creal(((double _Complex*)hp_temp)[offset+i*dims[0]+j]),
	       cimag(((double _Complex*)hp_temp)[offset+i*dims[0]+j]));
      }
      printf("\n");
    }
    printf("\n");

    /*for (klmn=0;klmn<2*lbg0;klmn++) {
      printf("klmn = %u, sum_cabs(hp_temp) = %f\n",klmn,
	     sum_cabs(((double _Complex*)hp_temp)+klmn*dimXYPlane,dimXYPlane));
	     }*/
    free(hp_temp);
    /* end debug */
#endif

    to_real_y_only_on_device((cufftDoubleComplex*)mystreams[iStream]->dp_temp,
			     df_outarr+(iPart*nXYPlanesPerPart*li0da*ly0da),
			     nXYPlanesPerPart,iStream);
    EVENT_RECORD(after_to_real[iStream],mystreams[iStream]->cudaStream);

#if 0
    /* debug */
    hp_temp = malloc(nXYPlanesPerPart*li0da*ly0da*sizeof(double));
    cudaMemcpy(hp_temp,df_outarr+(iPart*nXYPlanesPerPart*li0da*ly0da),
	       nXYPlanesPerPart*li0da*ly0da*sizeof(double),
	       cudaMemcpyDeviceToHost);
    /*for (klmn=0;klmn<nXYPlanesPerPart;klmn++) {
      printf("klmn = %u, %u, sum_cabs(df_outarr) = %17.6e\n",klmn,klmn*li0da*ly0da,
	     sum_abs(((double*)hp_temp)+klmn*li0da*ly0da,li0da*ly0da));
	     }*/

    printf("total sum_abs(df_outarr) = %f\n",sum_abs((double*)hp_temp,nXYPlanesPerPart*li0da*ly0da));
    free(hp_temp);
    /* end debug */
#endif

    /* stream iStream is now used as next, hence we have to 
       get the events before restart of the stream. */
#ifdef MEASURE_EVENTS
    nextStream = (iStream+1)%nStreams;
    /*
    if (cudaStreamWaitEvent(mystreams[iStream]->cudaStream,after_to_real[iStream],0)==cudaSuccess) {
      cuda_err = cudaEventElapsedTime(&time_transfer,start_transfer[iStream],end_transfer[iStream]);
      if (cuda_err==cudaSuccess) stime_transfer[iStream] += time_transfer;
      cuda_err = cudaEventElapsedTime(&time_copy_and_zero,end_transfer[iStream],after_copy_and_zero[iStream]);
      if (cuda_err==cudaSuccess) stime_copy_and_zero[iStream] += time_copy_and_zero;
      cuda_err = cudaEventElapsedTime(&time_transpose,after_copy_and_zero[iStream],after_transpose[iStream]);
      if (cuda_err==cudaSuccess) stime_transpose[iStream] += time_transpose;
      cuda_err = cudaEventElapsedTime(&time_to_real,after_transpose[iStream],after_to_real[iStream]);
      if (cuda_err==cudaSuccess) stime_to_real[iStream] += time_to_real;
    */
    if (cudaEventQuery(after_to_real[nextStream])==cudaErrorNotReady) {
      cuda_err = cudaEventSynchronize(after_to_real[nextStream]);
    }
    cuda_err = cudaEventElapsedTime(&time_transfer,start_transfer[nextStream],end_transfer[nextStream]);
    if (cuda_err==cudaSuccess) stime_transfer[nextStream] += time_transfer;
    cuda_err = cudaEventElapsedTime(&time_copy_and_zero,end_transfer[nextStream],after_copy_and_zero[nextStream]);
    if (cuda_err==cudaSuccess) stime_copy_and_zero[nextStream] += time_copy_and_zero;
    cuda_err = cudaEventElapsedTime(&time_transpose,after_copy_and_zero[nextStream],after_transpose[nextStream]);
    if (cuda_err==cudaSuccess) stime_transpose[nextStream] += time_transpose;
    cuda_err = cudaEventElapsedTime(&time_to_real,after_transpose[nextStream],after_to_real[nextStream]);
    if (cuda_err==cudaSuccess) stime_to_real[nextStream] += time_to_real;
    /* delete the errors. We get an error for the first call of the cudaEventElapsedTime as the
       events used there have not been called before. This gives an Invalid Resource Handle error.
       To remove it from the error stack, we have to call cudaGetLastError(). */
    
    cudaGetLastError();

#endif
    /* increment the stream id for the next loop cycle */
    iStream = (++iStream)%nStreams;
  }
}

double sum_cabs(double _Complex *array,int len) {
  return sum_cabs((cufftDoubleComplex*)array, len);
}

double sum_cabs(cufftDoubleComplex *array,int len) {
  double sum;
  int i;

  sum = 0.0;
  for (i=0;i<len;i++) {
    sum += sqrt(array[i].x*array[i].x+array[i].y*array[i].y);
  }
  return sum;
}

void show_xy_arr(cufftDoubleComplex *array,int nrows, int ncols) {
  int i,j,lind;
  for (i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      lind = i*ncols+j;
      printf("(%5.2f %5.2f)",array[lind].x, array[lind].y);
    }
    printf("\n");
  }
}

void show_xy_arr(double _Complex *array,int nrows, int ncols) {
  int i,j,lind;
  for (i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      lind = i*ncols+j;
      printf("(%5.2f %5.2f)",creal(array[lind]),cimag(array[lind]));
    }
    printf("\n");
  }
}

void error_handling_with_synchronize(const char *kernelname) {
  cudaError_t cuda_err;

  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) {
    printf("Error with kernel launch.(%s)\n",kernelname);
    printf("%s\n",cudaGetErrorString(cuda_err));
  }
  cuda_err = cudaDeviceSynchronize();
  if (cuda_err != cudaSuccess) {
    printf("Error with Synchronize after kernel launch. (%s)\n",kernelname);
    printf("%s\n",cudaGetErrorString(cuda_err));
  }
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
