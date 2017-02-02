#include "cuda_overlap.h"
#include "cuda_kernels.h"

/* The following routine has been found in the web at
   http://brianmykietka.net/projects.php?project=finalmatrixtranspose
   It references a paper "Optimizing Matrix Transpose in CUDA" by
   Greg Ruetsch (gruetsch@nvidia.com) and Paulius Micikevicius (pauliusm@nvidia.com)
*/


/* odata, output array, which is transposed
   idata, input array of dimension width x height
   width
   height
*/
#if 0
__global__ void transposeCoalescedBank_old(cufftDoubleComplex* odata, 
					   cufftDoubleComplex* idata, 
					   int width, int height)
{
	__shared__ cufftDoubleComplex tile[TILE_DIM][TILE_DIM + 1];
	
	int klmn = blockIdx.z;
	int offset = klmn*width*height;
	int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
	int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
	int index_in = offset + yIndex * width + xIndex ;
	
	/* The same for the output matrix, which is transposed. */
	xIndex = blockIdx.y * TILE_DIM + threadIdx.x;
	yIndex = blockIdx.x * TILE_DIM + threadIdx.y;
	int index_out = offset + yIndex * height + xIndex ;
	
	/* copy one tile per threadblock from global memory to shared memory */
	for (int i = 0; i < TILE_DIM; i += blockDim.y) {
	  tile[ threadIdx.y + i ][ threadIdx.x ] =
	    idata[ index_in + i * width];
	}
	
	__syncthreads();
	
	for (int i = 0; i < TILE_DIM; i += blockDim.y) {
	  odata[ index_out + i * height ] =
	    tile[ threadIdx.x ][ threadIdx.y + i ];
	}
}
#endif
__global__ void transposeCoalescedBank(cufftDoubleComplex* odata, 
				       cufftDoubleComplex* idata, 
				       int width, int height, int tile_dim)
{
  /*__shared__ cufftDoubleComplex tile[TILE_DIM][TILE_DIM + 1];*/
  extern __shared__ cufftDoubleComplex tile[];
	
  int klmn = blockIdx.z;
  int offset = klmn*width*height;
  int xIndex_in = blockIdx.x * tile_dim + threadIdx.x;
  int yIndex_in = blockIdx.y * tile_dim + threadIdx.y;
  int index_in = offset + yIndex_in * width + xIndex_in ;
	
  /* The same for the output matrix, which is transposed. */
  int xIndex = blockIdx.y * tile_dim + threadIdx.x;
  int yIndex = blockIdx.x * tile_dim + threadIdx.y;
  int index_out = offset + yIndex * height + xIndex ;
	
  /* copy one tile per threadblock from global memory to shared memory */
  if (yIndex_in < height) {
    /*printf("(%u,%u), (%u,%u), %2u %2u %2u | %2u %2u %2u\n",
      blockIdx.x,blockIdx.y,threadIdx.x,threadIdx.y,
      xIndex_in, yIndex_in, index_in, 
      xIndex, yIndex, index_out);*/
    /*tile[ threadIdx.y ][ threadIdx.x ] =
      idata[ index_in ];*/
    tile[ threadIdx.y*(tile_dim+1) + threadIdx.x ] =
      idata[ index_in ];
  }
	
  __syncthreads();
	
  if (xIndex<height) {
    /*odata[ index_out ] =
      tile[ threadIdx.x ][ threadIdx.y ];*/
    odata[ index_out ] =
      tile[ threadIdx.x *(tile_dim+1) + threadIdx.y ];
  }
}

void transpose_wrapper(cufftDoubleComplex* odata, 
		       cufftDoubleComplex* idata, 
		       int width, int height, int nXYPlanesPerPart, int iStream) {
  dim3 grid, threadblock;
  cudaError_t cuda_err;
  int minimum, max_tile_dim,eval_tile_dim,tile_dim;

  /* determine dynamically the TILE_DIM,
     it should be in the range 2->31. At the moment, we assume that
     width is a power of 2. */
  minimum=(width<height) ? width : height;
  max_tile_dim=16;
  if (minimum<16) max_tile_dim=8;
  if (minimum<8) max_tile_dim=4;
  if (minimum<4) max_tile_dim=2;

  for (eval_tile_dim=max_tile_dim;eval_tile_dim>1;eval_tile_dim=eval_tile_dim>>1) {
    if (width%eval_tile_dim==0) {
      tile_dim=eval_tile_dim;
      break;
    }
  }

  if (height%tile_dim) grid.y=height/tile_dim + 1;
  else grid.y=height/tile_dim;

  grid.x=width/tile_dim;   grid.z = nXYPlanesPerPart;
  threadblock.x=tile_dim;  threadblock.y=tile_dim;  threadblock.z=1;
  
  /*printf("Starting transposeCoalescedBank<<<(%u,%u,%u),(%u,%u,%u)>>>(,,%u,%u).\n",
	 grid.x,grid.y,grid.z,
	 threadblock.x,threadblock.y,threadblock.z,
	 width,height);*/
  transposeCoalescedBank<<<grid,threadblock,tile_dim*(tile_dim+1)*sizeof(cufftDoubleComplex),mystreams[iStream]->cudaStream>>>(odata,idata, width,height,tile_dim);
  
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) {
    printf("Error with kernel launch.(transposeCoalescedBank)\n");
    printf("%s\n",cudaGetErrorString(cuda_err));
  }
}

/* gridDim.x = nXYPlanesPerPart;
   gridDim.y = nj0;
   blockDim.x = li0da;
*/
__global__ void copy_with_pnl(cufftDoubleComplex *inarr, cufftDoubleComplex *outarr, cufftDoubleReal *prefactor, int width, int height) {
  /* inarr[nXYPlanesPerPart][nj0][li0da];
     outarr[nXYPlanesPerPart][nj0][li0da];
     prefactor[li0da];
  */
  int index = blockIdx.x*width*height + blockIdx.y*width + threadIdx.x;

  outarr[index].x += inarr[index].x * prefactor[threadIdx.x];
  outarr[index].y += inarr[index].y * prefactor[threadIdx.x];

}


/* Copies the inarr to a larger array (larger in height nj0 ->ly0da/2+1).
   Rows 0->nj0 are copied, nj0+1->ly0da/2+1 are set to zero.
   One threadblock per x-y plane and one thread for each point in the x-y plane
   of the larger output array.
   PROBLEM: This easily exceeds 1024 threads per block.
   Therefore, one block just works on one row.*/
__global__ void dev_copy_and_zero_for_dealiasing_new(cufftDoubleComplex *inarr, int iwidth, int iheight,
						     cufftDoubleComplex *outarr,int owidth, int oheight) {
  /* double _Complex inarr[howmany][iheight][iwidth];
     double _Complex outarr[howmany][oheight][owidth];
     gridDim.x = oheight; gridDim.y = howmany;
     blockDim.x = owidth; blockDim.y=1; 
  */
  int iOffset = blockIdx.z*iwidth*iheight+blockIdx.x*iwidth;
  int oOffset = blockIdx.z*owidth*oheight+blockIdx.x*owidth;

  int ilind;
  int olind = oOffset + threadIdx.x;

  if (blockIdx.x < iheight) {
    /* just copy */
    ilind = iOffset+threadIdx.x;
    outarr[olind] = inarr[ilind];
  } else {
    /* set entry to zero, x is real part, y is imaginary part. */
    outarr[olind].x = 0.0;
    outarr[olind].y = 0.0;
  }
}


#ifdef WITH_THRUST

/* does not compile due to some errors in the include files.
   but it is not suitable for GENE, as it does not support different
   CUDA streams, which is crucial for GENE performance. */
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/functional.h>
//#include <thrust/device_vector.h>
/*#include <thrust/fill.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
*/

/* first, we declare a transformation for the indices by deriving it from
   the unary_function */
struct skip_zero_transform_t : public thrust::unary_function<int,int> {

  int dimXYPlane;
  int dimExtendedXYPlane;

  /* the constructor initializes the two attributes of the struct */
  skip_zero_transform_t(int _dimXYPlane,int _dimExtendedXYPlane) : dimXYPlane(_dimXYPlane),dimExtendedXYPlane(_dimExtendedXYPlane) {}

  __host__ __device__
  int operator()(const int i)  const {
    return (i/dimXYPlane)*dimExtendedXYPlane+(i%dimXYPlane);
  }
};

/* and the complementary function */
struct only_zero_transform_t : public thrust::unary_function<int,int> {

  int dimXYPlane;
  int dimExtendedXYPlane;
  int dimExtension;

  /* the constructor initializes the two attributes of the struct */
  only_zero_transform_t(int _dimXYPlane,int _dimExtendedXYPlane) : dimXYPlane(_dimXYPlane),
								   dimExtendedXYPlane(_dimExtendedXYPlane) {
    dimExtension = _dimExtendedXYPlane-_dimXYPlane;
  }

  __host__ __device__
  int operator()(const int i) const {
    return dimXYPlane + (i/dimExtension)*dimExtendedXYPlane + (i%dimExtension);
  }
};


void copy_and_zero_for_dealiasing_wrapper(cufftDoubleComplex *inarr, int iwidth, int iheight,int idepth,
					  cufftDoubleComplex *outarr, int owidth, int oheight, int iStream) {

  cufftDoubleComplex zero_value;
  zero_value.x=0.0;
  zero_value.y=0.0;
  thrust::device_ptr<cufftDoubleComplex> dptr_inarr=thrust::device_pointer_cast(inarr);
  thrust::device_ptr<cufftDoubleComplex> dptr_outarr=thrust::device_pointer_cast(outarr);

  thrust::device_vector<cufftDoubleComplex> dv_inarr(dptr_inarr,dptr_inarr+iwidth*iheight*idepth);
  thrust::device_vector<cufftDoubleComplex> dv_outarr(dptr_outarr,dptr_outarr+owidth*oheight*idepth);

  thrust::counting_iterator<int> counter(0);
  skip_zero_transform_t skip_zero(iwidth*iheight,owidth*oheight); // instantiate
  only_zero_transform_t only_zero(iwidth*iheight,owidth*oheight);

  /* copy the inarr to the outarr, but only for the first iwidth*iheight entries per
     xy-plane. This scattering is done with the permutation_iterator. */
  //thrust::copy(dptr_inarr,dptr_inarr+iwidth*iheight*idepth,
  thrust::copy(dv_inarr.begin(),dv_inarr.end(),
	       thrust::make_permutation_iterator(dv_outarr.begin(),
						 thrust::make_transform_iterator(counter, skip_zero)
						 )
	       );

  thrust::fill_n(thrust::make_permutation_iterator(dv_outarr.begin(),
						 thrust::make_transform_iterator(counter, only_zero)
						 ),
	       (owidth*oheight-iwidth*iheight)*idepth,
	       zero_value);
}
#else
/* non-thrust version */
void copy_and_zero_for_dealiasing_wrapper(cufftDoubleComplex *inarr, int iwidth, int iheight,int idepth,
					  cufftDoubleComplex *outarr, int owidth, int oheight, int iStream) {

  dim3 grid, threadblock;

  //if (iwidth<maxThreadsPerBlock) {
    /* at least one line matches in the threadblock */
    //nLinesPerTile=maxThreadsPerBlock/iwidth;

  //grid.x=ly0da/2+1;        grid.y=(ly0da/2+1)/4;  grid.z=idepth;
  grid.x=ly0da/2+1;        grid.y=1;  grid.z=idepth;
  //threadblock.x=iwidth;    threadblock.y=4;
  threadblock.x=iwidth;    threadblock.y=1;
    
  //printf("before dev_copy_and_zero_for_dealiasing_new: %s\n",cudaGetErrorString(cudaGetLastError()));
  dev_copy_and_zero_for_dealiasing_new<<<grid,threadblock,0,mystreams[iStream]->cudaStream>>>
    (inarr,iwidth,iheight,  outarr,iwidth,ly0da/2+1);
  //cudaStreamSynchronize(mystreams[iStream]->cudaStream);
  //printf("after dev_copy_and_zero_for_dealiasing_new: %s\n",cudaGetErrorString(cudaGetLastError()));

}
#endif


/* find the maximal value in arr and return it in maxval,
   this function is specialized for the use with GENE, so it
   assumes the input array arr to have the dimensions:
   cufftDoubleReal arr[lbg0][2][li0da][ly0da];
   It finds the maximum individually for the two dimensions of the
   second index.

   gridDim.x=lbg0;
   blockDim.x = li0da;
   blockDim.y = 1;

   This routines needs a shared memory of 
   shared_mem = (li0<64 ? 64 : li0da)*2*sizeof(double);

*/
#if 0
__global__ void cuda_maxval_per_block_old(const cufftDoubleReal *arr, cufftDoubleReal *g_odata,
				      const int li0da, const int ly0da) {
  /* declaring extern the shared memory array means that it is allocated at launch time. */

  extern __shared__ cufftDoubleReal sdata[];

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    /*unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x;
      unsigned int gridSize = blockSize*2*gridDim.x;*/
    unsigned int blockSize = blockDim.x*blockDim.y*blockDim.z;
    int j, lind;
    cufftDoubleReal *sdata0 = &sdata[0];
    cufftDoubleReal *sdata1 = &sdata[blockDim.x];

    cufftDoubleReal myMax[2] = {0.0,0.0};
    
    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    for (j=0;j<ly0da;j++) {
      lind = blockIdx.x*2*ly0da*li0da+threadIdx.x*ly0da+j;
      if (arr[lind]>myMax[0]) {
	myMax[0]=arr[lind];
      }

      lind = (blockIdx.x*2+1)*ly0da*li0da+threadIdx.x*ly0da+j;
      if (arr[lind]>myMax[1]) {
	myMax[1]=arr[lind];
      }
    }

    // each thread puts its local sum into shared memory 
    sdata0[tid] = myMax[0];
    sdata1[tid] = myMax[1];
    __syncthreads();

    /* Now we have the maximum of each line in the shared memory. */
    // do reduction in shared mem
    if (blockSize >= 512) { 
      if (tid < 256) { 
	sdata0[tid] = (sdata0[tid]>sdata0[tid+256]) ? sdata0[tid] : sdata0[tid+256]; 
	sdata1[tid] = (sdata1[tid]>sdata1[tid+256]) ? sdata1[tid] : sdata1[tid+256]; 
      } __syncthreads(); 
    }
    if (blockSize >= 256) { 
      if (tid < 128) { 
	sdata0[tid] = (sdata0[tid]>sdata0[tid+128]) ? sdata0[tid] : sdata0[tid+128]; 
	sdata1[tid] = (sdata1[tid]>sdata1[tid+128]) ? sdata1[tid] : sdata1[tid+128]; 
      } __syncthreads(); 
    }
    if (blockSize >= 128) { 
      if (tid <  64) { 
	sdata0[tid] = (sdata0[tid]>sdata0[tid+64]) ? sdata0[tid] : sdata0[tid+64]; 
	sdata1[tid] = (sdata1[tid]>sdata1[tid+64]) ? sdata1[tid] : sdata1[tid+64]; 
      } __syncthreads(); 
    }
	
    
    if (tid < 32) {
        // now that we are using warp-synchronous programming (below)
        // we need to declare our shared memory volatile so that the compiler
        // doesn't reorder stores to it and induce incorrect behavior.
      volatile cufftDoubleReal* smem = sdata0;
      if (blockSize >=  64) { smem[tid] = (smem[tid]>smem[tid+32]) ? smem[tid] : smem[tid+32]; }
      if (blockSize >=  32) { smem[tid] = (smem[tid]>smem[tid+16]) ? smem[tid] : smem[tid+16]; }
      if (blockSize >=  16) { smem[tid] = (smem[tid]>smem[tid+8])  ? smem[tid] : smem[tid+8]; }
      if (blockSize >=  8)  { smem[tid] = (smem[tid]>smem[tid+4])  ? smem[tid] : smem[tid+4]; }
      if (blockSize >=  4)  { smem[tid] = (smem[tid]>smem[tid+2])  ? smem[tid] : smem[tid+2]; }
      if (blockSize >=  2)  { smem[tid] = (smem[tid]>smem[tid+1])  ? smem[tid] : smem[tid+1]; }
      smem = sdata1;
      if (blockSize >=  64) { smem[tid] = (smem[tid]>smem[tid+32]) ? smem[tid] : smem[tid+32]; }
      if (blockSize >=  32) { smem[tid] = (smem[tid]>smem[tid+16]) ? smem[tid] : smem[tid+16]; }
      if (blockSize >=  16) { smem[tid] = (smem[tid]>smem[tid+8]) ? smem[tid] : smem[tid+8]; }
      if (blockSize >=  8) { smem[tid] = (smem[tid]>smem[tid+4]) ? smem[tid] : smem[tid+4]; }
      if (blockSize >=  4) { smem[tid] = (smem[tid]>smem[tid+2]) ? smem[tid] : smem[tid+2]; }
      if (blockSize >=  2) { smem[tid] = (smem[tid]>smem[tid+1]) ? smem[tid] : smem[tid+1]; }
    }
    
    // write result for this block to global mem 
    if (tid == 0)  {
      g_odata[blockIdx.x] = sdata0[0];
      g_odata[gridDim.x+blockIdx.x] = sdata1[0];
    }
}
#endif

/** Finding the maximum value of an x-y plane for the two parts of an array
    separately.

    Each column of the input array is reduced by one thread and in one tile, that
    means over li0da/NTILES rows.

    \param arr The input array. The structure of this array is [lbg0][2][li0da][ly0da].
    \param g_odata The output array. It contains after the routine the result and is of the
    structure [2][lbg0][NTILES]. It is then suitable for usage of the reduction kernel
    for the two separate parts to find the global maximum.
    \param li0da The number of rows of an x-y plane (the x direction).
    \param ly0da The number of columns of a x-y plane (the y direction).

    The launch parameters should be:
    grid.x = lbg0
    grid.y = NTILES
    grid.z = 1
    threadblock.x = ly0da, each column is one thread
    threadblock.y = 1
    threadblock.z = 1
    shared_mem = ((ly0da<32) ? 64 : 2*ly0da)*sizeof(double);

    \constraint li0da must be divisable by NTILES
 */
__global__ void cuda_maxval_per_block(const cufftDoubleReal *arr, cufftDoubleReal *g_odata,
				      const int li0da, const int ly0da) {
  /* declaring extern the shared memory array means that it is allocated at launch time. */

  extern __shared__ cufftDoubleReal sdata[];

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    /*unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x;
      unsigned int gridSize = blockSize*2*gridDim.x;*/
    unsigned int blockSize = max(blockDim.x,32);//blockDim.x;
    int lind, iXYPlane,iBlockrow,compare_lind, reducedToPowerOfTwo;
    int diff;
    cufftDoubleReal *sdata0 = &sdata[0];
    cufftDoubleReal *sdata1 = &sdata[blockDim.x];

    //cufftDoubleReal myMax[2] = {0.0,0.0};
    
    iXYPlane = blockIdx.x*2*li0da*ly0da;
    iBlockrow = blockIdx.y*ly0da*li0da/NTILES;
    //lind = iXYPlane + iBlockrow + threadIdx.y*ly0da + threadIdx.x;
    lind = iXYPlane + iBlockrow + threadIdx.x;

    /* Now compare one row with the li0da/NTILES-1 rows in the tile. */
    compare_lind = lind+ly0da;
    sdata0[tid] = (arr[lind]>arr[compare_lind]) ? arr[lind] : arr[compare_lind];
    compare_lind += ly0da;
    while (compare_lind < lind+ly0da*li0da/NTILES) {
      sdata0[tid] = (sdata0[tid]>arr[compare_lind]) ? sdata0[tid] : arr[compare_lind];
      compare_lind += ly0da;
    }
    /* In sdata0 the column maxima over the tile are stored now. */

    iXYPlane = (2*blockIdx.x+1)*li0da*ly0da;
    iBlockrow = blockIdx.y*ly0da*li0da/NTILES;
    //lind = iXYPlane + iBlockrow + threadIdx.y*ly0da + threadIdx.x;
    lind = iXYPlane + iBlockrow + threadIdx.x;

    compare_lind = lind+ly0da;
    sdata1[tid] = (arr[lind]>arr[compare_lind]) ? arr[lind] : arr[compare_lind];
    compare_lind += ly0da;
    while (compare_lind < lind+ly0da*li0da/NTILES) {
      sdata1[tid] = (sdata1[tid]>arr[compare_lind]) ? sdata1[tid] : arr[compare_lind];
      compare_lind += ly0da;
    }


    // each thread puts its local sum into shared memory 
    __syncthreads();

    /* Now we have one line with the maxima of each column in the shared memory. 
       BUT, the number of threads and the blockSize can be arbitrary (usually
       divisable by 2, but not necessarily a power of 2), so we have to modify the
       algorithm. */
    
    // do reduction in shared mem
    reducedToPowerOfTwo=0;
    if (blockSize >= 512) { 
      diff = blockSize-512;
      if (tid < diff) {
	sdata0[tid] = (sdata0[tid]>sdata0[512+tid]) ? sdata0[tid] : sdata0[512+tid]; 
	sdata1[tid] = (sdata1[tid]>sdata1[tid+512]) ? sdata1[tid] : sdata1[tid+512]; 
      }
      __syncthreads();
      reducedToPowerOfTwo=1;
      if (tid < 256) { 
	sdata0[tid] = (sdata0[tid]>sdata0[tid+256]) ? sdata0[tid] : sdata0[tid+256]; 
	sdata1[tid] = (sdata1[tid]>sdata1[tid+256]) ? sdata1[tid] : sdata1[tid+256]; 
      } __syncthreads(); 
    }

    if (blockSize >= 256) { 
      diff = blockSize-256;
      if (!reducedToPowerOfTwo && (tid < diff)) {
	sdata0[tid] = (sdata0[tid]>sdata0[tid+256]) ? sdata0[tid] : sdata0[tid+256]; 
	sdata1[tid] = (sdata1[tid]>sdata1[tid+256]) ? sdata1[tid] : sdata1[tid+256]; 
      } __syncthreads();
      reducedToPowerOfTwo=1;

      if (tid < 128) { 
	sdata0[tid] = (sdata0[tid]>sdata0[tid+128]) ? sdata0[tid] : sdata0[tid+128]; 
	sdata1[tid] = (sdata1[tid]>sdata1[tid+128]) ? sdata1[tid] : sdata1[tid+128]; 
      } __syncthreads(); 
    }
    if (blockSize >= 128) { 
      diff = blockSize - 128;
      if (!reducedToPowerOfTwo && (tid < diff)) {
	sdata0[tid] = (sdata0[tid]>sdata0[tid+128]) ? sdata0[tid] : sdata0[tid+128]; 
	sdata1[tid] = (sdata1[tid]>sdata1[tid+128]) ? sdata1[tid] : sdata1[tid+128]; 
      } __syncthreads();
      reducedToPowerOfTwo=1;
      if (tid <  64) { 
	sdata0[tid] = (sdata0[tid]>sdata0[tid+64]) ? sdata0[tid] : sdata0[tid+64]; 
	sdata1[tid] = (sdata1[tid]>sdata1[tid+64]) ? sdata1[tid] : sdata1[tid+64]; 
      } __syncthreads(); 
    }
	
    if (!reducedToPowerOfTwo && (blockSize>32)) {
      diff = blockSize-32;
      if (tid<diff) {
	sdata0[tid] = (sdata0[tid]>sdata0[tid+32]) ? sdata0[tid] : sdata0[tid+32]; 
	sdata1[tid] = (sdata1[tid]>sdata1[tid+32]) ? sdata1[tid] : sdata1[tid+32]; 
      }
    }
    if (tid < 32) {
        // now that we are using warp-synchronous programming (below)
        // we need to declare our shared memory volatile so that the compiler
        // doesn't reorder stores to it and induce incorrect behavior.
      volatile cufftDoubleReal* smem = sdata0;
      if (blockSize >=  64) { smem[tid] = (smem[tid]>smem[tid+32]) ? smem[tid] : smem[tid+32]; }
      if (blockSize >=  32) { smem[tid] = (smem[tid]>smem[tid+16]) ? smem[tid] : smem[tid+16]; }
      if (blockSize >=  16) { smem[tid] = (smem[tid]>smem[tid+8])  ? smem[tid] : smem[tid+8]; }
      if (blockSize >=  8)  { smem[tid] = (smem[tid]>smem[tid+4])  ? smem[tid] : smem[tid+4]; }
      if (blockSize >=  4)  { smem[tid] = (smem[tid]>smem[tid+2])  ? smem[tid] : smem[tid+2]; }
      if (blockSize >=  2)  { smem[tid] = (smem[tid]>smem[tid+1])  ? smem[tid] : smem[tid+1]; }
      smem = sdata1;
      if (blockSize >=  64) { smem[tid] = (smem[tid]>smem[tid+32]) ? smem[tid] : smem[tid+32]; }
      if (blockSize >=  32) { smem[tid] = (smem[tid]>smem[tid+16]) ? smem[tid] : smem[tid+16]; }
      if (blockSize >=  16) { smem[tid] = (smem[tid]>smem[tid+8]) ? smem[tid] : smem[tid+8]; }
      if (blockSize >=  8) { smem[tid] = (smem[tid]>smem[tid+4]) ? smem[tid] : smem[tid+4]; }
      if (blockSize >=  4) { smem[tid] = (smem[tid]>smem[tid+2]) ? smem[tid] : smem[tid+2]; }
      if (blockSize >=  2) { smem[tid] = (smem[tid]>smem[tid+1]) ? smem[tid] : smem[tid+1]; }
    }
    
    // write result for this block to global mem 
    if (tid == 0)  {
      g_odata[blockIdx.x*gridDim.y+blockIdx.y] = sdata0[0];
      g_odata[ gridDim.x*gridDim.y+blockIdx.x*gridDim.y+blockIdx.y] = sdata1[0];
    }
}

/* compute the standard nonlinearity
   1. Attempt:
   grid.x = lbg0
   grid.y = NTILES_NONLIM
   threadblock.x=ly0da
   threadblock.y=li0da/NTILES_NONLIM
*/
__global__ void comp_stand_nonlin(double *vexy_re, double *dgdxy_re, double* nonlin1_re,int li0da, int ly0da) {
  /* offset0 is the start index of the first component for the given threadblock */
  int offset0 =  2*blockIdx.x   *li0da*ly0da + blockIdx.y*blockDim.y*blockDim.x;
  int offset1 = (2*blockIdx.x+1)*li0da*ly0da + blockIdx.y*blockDim.y*blockDim.x;
  int tx = threadIdx.x;
  int ty = threadIdx.y;

  //for (i=0;i<li0da;i++) {
  nonlin1_re[ blockIdx.x*li0da*ly0da + blockIdx.y*blockDim.y*blockDim.x + ty*ly0da + tx ] = 
    - vexy_re[ offset0 + ty*ly0da + tx ]*dgdxy_re[ offset1+ty*ly0da + tx ]
    + vexy_re[ offset1 + ty*ly0da + tx ]*dgdxy_re[ offset0+ty*ly0da + tx ];
    //}
}

/* compute the standard nonlinearity
   1. Attempt:
   grid.x = lbg0
   threadblock.x=ly0da
*/
#if 0
__global__ void comp_stand_nonlin_old(double *vexy_re, double *dgdxy_re, double* nonlin1_re,int li0da, int ly0da) {
  /* offset0 is the start index of the first component for the given threadblock */
  int offset0 =  2*blockIdx.x   *li0da*ly0da;
  int offset1 = (2*blockIdx.x+1)*li0da*ly0da;

  int i;

  /*for (j=0;j<ly0da;j++) {
    nonlin1_re[ blockIdx.x*li0da*ly0da + threadIdx.x*ly0da + j ] = 
      - vexy_re[offset0+threadIdx.x*ly0da+j]*dgdxy_re[offset1+threadIdx.x*ly0da+j]
      + vexy_re[offset1+threadIdx.x*ly0da+j]*dgdxy_re[offset0+threadIdx.x*ly0da+j];
      }*/
  for (i=0;i<li0da;i++) {
    nonlin1_re[ blockIdx.x*li0da*ly0da + i*ly0da + threadIdx.x ] = 
      - vexy_re[ offset0 + i*ly0da + threadIdx.x ]*dgdxy_re[ offset1+i*ly0da + threadIdx.x ]
      + vexy_re[ offset1 + i*ly0da + threadIdx.x ]*dgdxy_re[ offset0+i*ly0da + threadIdx.x ];
  }
}
#endif

void comp_stand_nonlin_wrapper(double *vexy_re, double *dgdxy_re, double *nonlin1_re, 
			       int li0da, int ly0da, int iStream) {
  dim3 grid, threadblock;
  cudaError_t cuda_err;
  int nXYPoints,nTiles,nLinesPerTile;

  nXYPoints = li0da*ly0da;
  if (nXYPoints<=maxThreadsPerBlock) {
    grid.x=lbg0;
    grid.y=1;
    threadblock.x=ly0da;
    threadblock.y=li0da;
  } else {
    /* separate the height in several tiles */
    nLinesPerTile = maxThreadsPerBlock/ly0da;
    while (li0da%nLinesPerTile) nLinesPerTile--;
    nTiles = li0da/nLinesPerTile;
    
    if (li0da%nTiles) {
      printf("li0da (%u) must be divisable by nTiles (%u)\n",li0da,nTiles);
    }
    //printf("nLinesPerTile = %u, nTiles = %u\n",nLinesPerTile,nTiles);
    grid.x=lbg0;
    grid.y=nTiles;
    threadblock.x=ly0da;
    threadblock.y=li0da/nTiles;
  }

  /*grid.x=lbg0;
  grid.y=NTILES_NONLIN;
  threadblock.x=ly0da;
  threadblock.y=li0da/NTILES_NONLIN;*/
  /*printf("Calling kernel comp_stand_nonlin with grid(%u,%u,%u) and block(%u,%u,%u), li0da=%u, ly0da=%u.\n",
    grid.x,grid.y,grid.z,threadblock.x,threadblock.y,threadblock.z,li0da,ly0da);*/
  comp_stand_nonlin<<<grid,threadblock,0,mystreams[iStream]->cudaStream>>>(vexy_re,dgdxy_re,nonlin1_re,li0da,ly0da);
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) {
    printf("Error with kernel launch.(comp_stand_nonlin)\n");
    printf("%s\n",cudaGetErrorString(cuda_err));
  }
}
