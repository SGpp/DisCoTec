#pragma once

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <iostream>
#include <string>

#ifdef DISCOTEC_USE_LZ4
#include "lz4.h"
#include "lz4frame.h"
#endif

#include "mpi/MPISystem.hpp"
#include "utils/Types.hpp"

namespace combigrid {

namespace mpiio {
template <typename T>
void compressBufferToLZ4String(const T* buffer, MPI_Offset numValues,
                               combigrid::CommunicatorType comm,
                               std::vector<char>& compressedString) {
  size_t compressedEmptySize = 0;
  size_t compressedSkippableSize = 0;
#ifdef DISCOTEC_USE_LZ4
  compressedString.clear();
  size_t rawDataSize = numValues * sizeof(T);
  LZ4F_preferences_t lz4Preferences = {{LZ4F_max4MB, LZ4F_blockLinked, LZ4F_noContentChecksum,
                                        LZ4F_frame, rawDataSize, 0U, LZ4F_noBlockChecksum},
                                       0,   // fast mode
                                       0u,  // don't flush
                                       0u,  // unused
                                       {0u, 0u, 0u}};
  auto dataBound = LZ4_compressBound(rawDataSize);
  auto frameBound = LZ4F_compressFrameBound(dataBound, &lz4Preferences);
  compressedString.reserve(frameBound);
  auto compressedStringStart = static_cast<void*>(compressedString.data());
  
  compressedString.resize(compressedString.size() + frameBound);
  auto compressedSize =
      LZ4F_compressFrame(compressedStringStart, compressedString.size(),
                         static_cast<const void*>(buffer), rawDataSize, &lz4Preferences);
  if (LZ4F_isError(compressedSize)) {
    throw std::runtime_error("LZ4 compression failed: " +
                             std::string(LZ4F_getErrorName(compressedSize)));
  }
  assert(compressedSize > 0);
  assert(compressedSize <= rawDataSize);
  compressedString.resize(compressedEmptySize + compressedSkippableSize + compressedSize);
#else
  throw std::runtime_error("LZ4 compression not available");
#endif
}

}  // namespace mpiio
}  // namespace combigrid
