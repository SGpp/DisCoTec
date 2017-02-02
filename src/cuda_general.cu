#include <cuda_profiler_api.h>

extern "C" void start_cuda_profiling() {
  cudaProfilerStart();
}

extern "C" void end_cuda_profiling() {
  cudaProfilerStop();
}
