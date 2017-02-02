#include "cuda_runtime.h"
#include "cufft.h"

class d2 {
 public:
  double d[2];
  
  d2() {this->d[0]=0.0;this->d[1]=0.0;}
  d2(int x) {this->d[0]=x;this->d[1]=0;}
  void setvalue(double x, double y) {d[0]=x; d[1]=y;}
  d2& operator+(d2& x) {
    d2 *res=new d2;
    res->d[0] = this->d[0]+x.d[0];
    res->d[1] = this->d[1]+x.d[1];
    return *res;
  }
  d2& operator+(volatile d2& x) {
    d2 *res=new d2;
    res->d[0] = this->d[0]+x.d[0];
    res->d[1] = this->d[1]+x.d[1];
    return *res;
  }
  d2& operator+=(d2& x) {
    this->d[0] += x.d[0];
    this->d[1] += x.d[1];
    return *this;
  }

  d2& operator=(d2 x) {
    this->d[0]=x.d[0];
    this->d[1]=x.d[1];
    return *this;
  }

  volatile d2& operator=(d2 x) volatile {
    this->d[0]=x.d[0];
    this->d[1]=x.d[1];
    return *this;
  }

  /*  d2& operator=(volatile d2& x) {
    this->d[0]=x.d[0];
    this->d[1]=x.d[1];
    return *this;
    }*/
};

#include "cuda_kernels.cu"
#include "reduction_kernel.cu"

void init_cuda() {
  cudaError_t cuda_err;
  int number_of_devices;

  cuda_err = cudaGetDeviceCount(&number_of_devices);
  if (cuda_err != cudaSuccess) {
    printf("No GPU devices found.\n");
    exit(1);
  }
  if (number_of_devices>=3) {
    cuda_err = cudaSetDevice(2);
    if (cuda_err != cudaSuccess) {
      printf(cudaGetErrorString(cuda_err));
    }
  }
}

void test_cuda_transpose() {
  /* matrix on the host */
  cufftDoubleComplex *mat;
  cufftDoubleComplex *tmat;
  dim3 grid;
  dim3 threadblock;
  int lind,klmn,blocksize;
  int width, height,i,j,twidth,theight;
  cufftDoubleComplex *dev_mat, *dev_tmat;

  width=32;
  height=9;
  blocksize = 1;

  mat = (cufftDoubleComplex*)malloc(width*height*blocksize*sizeof(cufftDoubleComplex));
  tmat = (cufftDoubleComplex*)malloc(width*height*blocksize*sizeof(cufftDoubleComplex));

  for (klmn=0;klmn<blocksize;klmn++) {
    printf("-------------------- Input block %u ------------------\n",klmn);
    for (i=0;i<height;i++) {
      for (j=0;j<width;j++) {
	lind = klmn*height*width+i*width+j;
	mat[lind].x = (double)i+klmn*100.0;
	mat[lind].y = (double)j+klmn*100.0;
	printf("(%3.0f,%3.0f) ",mat[lind].x,mat[lind].y);
      }
      printf("\n");
    }
  }

  /* allocate on GPU */
  cudaMalloc((void**)&dev_mat,width*height*blocksize*sizeof(cufftDoubleComplex));
  cudaMalloc((void**)&dev_tmat,width*height*blocksize*sizeof(cufftDoubleComplex));

  /* transfer mat to the GPU */
  cudaMemcpy(dev_mat,mat,width*height*blocksize*sizeof(cufftDoubleComplex),cudaMemcpyHostToDevice);

  /* call the kernel */
  printf("width = %u, height= %u, TILE_DIM = %u,\nheight/TILE_DIM = %u R %u\n",
	 width, height,TILE_DIM, height/TILE_DIM, height%TILE_DIM);
  grid.x=width/TILE_DIM;

  grid.y = (height%TILE_DIM==0) ? height/TILE_DIM : height/TILE_DIM + 1;
  grid.z = blocksize;
  threadblock.x=TILE_DIM;
  threadblock.y=TILE_DIM;
  printf("grid = (%u,%u,%u), threadblock = (%u,%u,%u)\n",grid.x,grid.y,grid.z,
	 threadblock.x,threadblock.y, threadblock.z);
  transposeCoalescedBank<<<grid,threadblock>>>(dev_tmat,dev_mat,width,height);
    
  for (klmn=0;klmn<blocksize;klmn++) {
    //cudaDeviceSynchronize();
    cudaMemcpy(&tmat[klmn*width*height],&dev_tmat[klmn*width*height],
	       width*height*sizeof(cufftDoubleComplex),cudaMemcpyDeviceToHost);
  }
  twidth = height;
  theight = width;

  for (klmn=0;klmn<blocksize;klmn++) {
    printf("------------------- Output block %u ------------------\n",klmn);
    for (i=0;i<theight;i++) {
      for (j=0;j<twidth;j++) {
	lind = klmn*theight*twidth+i*twidth+j;
	printf("(%3.0f,%3.0f) ",tmat[lind].x,tmat[lind].y);
      }
      printf("\n");
    }
  }
  
  cudaFree(dev_mat);
  cudaFree(dev_tmat);
  free(mat);
  free(tmat);
}

void test_cuda_maxval_per_block() {
  /* first test the sum reduction */
  double *d_idata, *d_odata;
  double *idata, *odata;
  int i,j,li0,lj0,lbg0,klmn;
  dim3 grid,threadblock;
  int shared_mem, iTile, lines_per_tile,passed;
  int iConfi, iConfj, iConfg;
  const int nConfi = 1;
  const int nConfj = 1;
  const int nConfg = 2;
  //const int Confi[nConfi]={2,7,10,16,24,32,128,384,517};
  const int Confi[nConfi]={5};
  const int Confj[nConfj]={8};
  const int Confg[nConfg]={1,2};
  /*li0 = 24;
  lj0 = 8;
  lbg0 = 2;*/

  /* loop over different configurations */
  for (iConfi=0;iConfi<nConfi;iConfi++) {
    for (iConfj=0;iConfj<nConfj;iConfj++) {
      for (iConfg=0;iConfg<nConfg;iConfg++) {
	li0  = Confi[iConfi];
	lj0  = Confj[iConfj];
	lbg0 = Confg[iConfg];

	// Input arrays on host and device
	idata = (double*)malloc(li0*lj0*lbg0*2*sizeof(double));
	cudaMalloc((void**)&d_idata,li0*lj0*lbg0*2*sizeof(double));
	for (klmn=0;klmn<lbg0;klmn++) {
	  for (i=0;i<li0;i++) {
	    for (j=0;j<lj0;j++) {
	      idata[klmn*2*li0*lj0+i*lj0+j] = 1000*klmn+100*j+(double)(i+1);
	      idata[(2*klmn+1)*li0*lj0+i*lj0+j] = 100*klmn+10*j+80.0-(double)(0.1*(i+1));
	      //printf("(%7.1f, %7.1f)",idata[klmn*2*li0*lj0+i*lj0+j],idata[(2*klmn+1)*li0*lj0+i*lj0+j]);
	    }
	    //printf("\n");
	  }
	  //printf("-----------\n");
	}
	cudaMemcpy(d_idata,idata,li0*lj0*lbg0*2*sizeof(double),cudaMemcpyHostToDevice);
	
	// Output arrays on host and device
	odata = (double*)malloc(lbg0*2*NTILES*sizeof(double));
	cudaMalloc((void**)&d_odata,lbg0*2*NTILES*sizeof(double));
	
	grid.x=lbg0;
	grid.y=NTILES;
	threadblock.x=lj0;
	threadblock.y=1;
	shared_mem = (li0<64 ? 64 : li0)*2*sizeof(double);
	cuda_maxval_per_block<<<grid,threadblock,shared_mem>>>(d_idata, d_odata,li0, lj0);
	
	// Transfer the result from device to host
	cudaMemcpy(odata,d_odata,2*lbg0*NTILES*sizeof(double),cudaMemcpyDeviceToHost);
	
	passed = 1;
	lines_per_tile = li0/NTILES;
	for (klmn=0;klmn<lbg0;klmn++) {
	  for (iTile=0;iTile<NTILES;iTile++) {
#if 1
	    if ((odata[klmn*NTILES+iTile] != 1000*klmn+100*(lj0-1)+(double)(iTile+1)*lines_per_tile) ||
		(odata[lbg0*NTILES+klmn*NTILES+iTile] != 100*klmn+10*(lj0-1)+80.0-0.1*(iTile*lines_per_tile+1))) {
	      passed=0;
	      printf("klmn=%u, tile %u: %f\t%f, reference is %f\t%f\n",
		     klmn,iTile,odata[klmn*NTILES+iTile],odata[lbg0*NTILES+klmn*NTILES+iTile],
		     1000*klmn+100*(lj0-1)+(double)(iTile+1)*lines_per_tile,
		     100*klmn+10*(lj0-1)+80.0-0.1*(iTile*lines_per_tile+1));
	    }
#else
	    printf("klmn=%u, tile %u: %f\t%f, reference is %f\t%f\n",
		   klmn,iTile,odata[klmn*NTILES+iTile],odata[lbg0*NTILES+klmn*NTILES+iTile],
		   1000*klmn+100*(lj0-1)+(double)(iTile+1)*lines_per_tile,
		   100*klmn+10*(lj0-1)+80.0-0.1*(iTile*lines_per_tile+1));
#endif
	  }
	}
	printf("li0 = %u, lj0 = %u, lbg0 = %u: ",li0,lj0,lbg0);
	if (passed) printf("passed.\n");
	else printf("\tFAILED!\n");
	
	cudaFree(d_odata);
	cudaFree(d_idata);
	free(odata);
	free(idata);
      }
    }
  }
}

void test_cuda_reduction() {

  /* first test the sum reduction */
  double *d_idata, *d_odata, *d_redres0, *d_redres1;
  double *idata, *odata,redres0,redres1;
  int i,li0,lj0,lbg0,j,klmn;
  int whichKernel = 6;
  dim3 grid,threadblock;
  int shared_mem;
  li0 = 32;
  lj0 = 8;
  lbg0 = 16;

  /*size = 4096;
  threads = 128;
  blocks = size/threads;
  printf("We are using %u blocks with %u threads,\n",blocks,threads);*/

  // Input arrays on host and device
  idata = (double*)malloc(li0*lj0*lbg0*2*sizeof(double));
  cudaMalloc((void**)&d_idata,li0*lj0*lbg0*2*sizeof(double));
  for (klmn=0;klmn<lbg0;klmn++) {
    for (i=0;i<li0;i++) {
      for (j=0;j<lj0;j++) {
	idata[klmn*2*li0*lj0+i*lj0+j] = 1000*klmn+100*j+(double)(i+1);
	idata[(2*klmn+1)*li0*lj0+i*lj0+j] = 100*klmn+10*j+80.0-(double)(0.1*(i+1));
	//printf("(%f, %f)",idata[klmn*2*li0*lj0+i*lj0+j],idata[(2*klmn+1)*li0*lj0+i*lj0+j]);
      }
      //printf("\n");
    }
    //printf("-----------\n");
  }
  cudaMemcpy(d_idata,idata,li0*lj0*lbg0*2*sizeof(double),cudaMemcpyHostToDevice);

  // Output arrays on host and device
  odata = (double*)malloc(lbg0*2*sizeof(double));
  cudaMalloc((void**)&d_odata,lbg0*2*sizeof(double));
  cudaMalloc((void**)&d_redres0,sizeof(double));
  cudaMalloc((void**)&d_redres1,sizeof(double));

  grid.x=lbg0;
  threadblock.x=li0;
  threadblock.y=1;
  shared_mem = (li0<64 ? 64 : li0)*2*sizeof(double);
  cuda_maxval_per_block<<<grid,threadblock,shared_mem>>>(d_idata, d_odata,li0, lj0);

  //reduce<double>(size,threads,blocks,whichKernel, d_idata, d_odata,REDUCTION_MAX);

  // Next level reduce
  /*if (whichKernel==6) {
    size=blocks/2;
  } else {
    size = blocks;
    }*/
  reduce<double>(lbg0,whichKernel,d_odata,d_redres0,REDUCTION_MAX,0);
  reduce<double>(lbg0,whichKernel,d_odata,d_redres1,REDUCTION_MAX,lbg0);

  // Transfer the result from device to host
  cudaMemcpy(odata,d_odata,2*lbg0*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(&redres0,d_redres0,sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(&redres1,d_redres1,sizeof(double),cudaMemcpyDeviceToHost);

  // print the result of the reduction
  // for comparison sum i_{i=1}^N=0.5*(N+1)*N
  printf("Intermediate result:\n");
  for (klmn=0;klmn<lbg0;klmn++) {
    printf("%f %f\n",odata[klmn],odata[lbg0+klmn]);
  }
  printf("Final result: %f %f\n",redres0, redres1);
  
}

int main(int argc, char *argv[]) {
  /* test the transpose on the GPU */

  init_cuda();
  test_cuda_transpose();
  //test_cuda_maxval_per_block();
  //test_cuda_reduction();
}

