#include "cuda_runtime.h"

#ifndef __CUDA_OVERLAP_H
#define __CUDA_OVERLAP_H
const unsigned int nStreams=2;
const unsigned int nParts=4;
const int maxThreadsPerBlock=1024;

class Streamdata {
 public:
  cudaStream_t cudaStream;
  void *dp_data;
  void *dp_temp;
  void *dp_fordeal;

  Streamdata(int li0,int lj0,int ly0da, int nXYPlanes);
  ~Streamdata();
};

#endif
