#include <stdio.h>
#include "cuda_overlap.h"

extern unsigned long int allocatedDeviceMemory;

Streamdata::Streamdata(int li0,int lj0,int ly0da, int nXYPlanes) {
  cudaError_t cuda_err;
  int dimXYPlane;
  int nXYPlanesPerPart;

  dimXYPlane = li0*lj0;
  nXYPlanesPerPart = nXYPlanes/nParts;

  /*printf("Allocating dp_data with %lu bytes.\n",nXYPlanesPerPart*dimXYPlane*sizeof(double _Complex));*/
  cuda_err = cudaMalloc((void**)&dp_data,nXYPlanesPerPart*dimXYPlane*sizeof(double _Complex));
  if (cuda_err != cudaSuccess) {
    printf("overlap cudaMalloc 0: %s\n",cudaGetErrorString(cuda_err));
  } else {
    allocatedDeviceMemory += nXYPlanesPerPart*dimXYPlane*sizeof(double _Complex);
  }

  /*printf("Allocating dp_temp with %lu bytes.\n",nXYPlanesPerPart*li0*(ly0da/2+1)*sizeof(double _Complex));*/
  cuda_err = cudaMalloc((void**)&dp_temp,nXYPlanesPerPart*li0*(ly0da/2+1)*sizeof(double _Complex));
  if (cuda_err != cudaSuccess) {
    printf("overlap cudaMalloc 1: %s\n",cudaGetErrorString(cuda_err));
  } else {
    allocatedDeviceMemory += nXYPlanesPerPart*li0*(ly0da/2+1)*sizeof(double _Complex);
  }

  /*printf("Allocating dp_data with %lu bytes.\n",nXYPlanesPerPart*li0*(ly0da/2+1)*sizeof(double _Complex));*/
  cuda_err = cudaMalloc((void**)&dp_fordeal,nXYPlanesPerPart*li0*(ly0da/2+1)*sizeof(double _Complex));
  if (cuda_err != cudaSuccess) {
    printf("overlap cudaMalloc 2: %s\n",cudaGetErrorString(cuda_err));
  } else {
    allocatedDeviceMemory += nXYPlanesPerPart*li0*(ly0da/2+1)*sizeof(double _Complex);
  }

  //cudaStreamCreateWithFlags(&cudaStream,cudaStreamNonBlocking);
  cudaStreamCreate(&cudaStream);

}

Streamdata::~Streamdata() {
  cudaStreamDestroy(cudaStream);

  cudaFree(dp_data);
  cudaFree(dp_temp);
  cudaFree(dp_fordeal);
}
