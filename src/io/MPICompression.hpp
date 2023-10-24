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

#ifdef DISCOTEC_USE_LZ4
constexpr LZ4F_preferences_t lz4Preferences = {
    {LZ4F_max4MB, LZ4F_blockLinked, LZ4F_noContentChecksum, LZ4F_frame, 0, 0U,
     LZ4F_noBlockChecksum},
    0,   // fast mode
    0u,  // don't flush
    0u,  // unused
    {0u, 0u, 0u}};
#endif

template <typename T>
void compressBufferToLZ4Frame(const T* buffer, MPI_Offset numValues,
                              combigrid::CommunicatorType comm,
                              std::vector<char>& compressedString) {
  size_t compressedStringInitialSize = compressedString.size();
#ifdef DISCOTEC_USE_LZ4
  size_t rawDataSize = numValues * sizeof(T);
  auto fullFramePreferences = lz4Preferences;
  fullFramePreferences.frameInfo.contentSize = rawDataSize;
  auto dataBound = LZ4_compressBound(rawDataSize);
  auto frameBound = LZ4F_compressFrameBound(dataBound, &fullFramePreferences);
  compressedString.reserve(compressedStringInitialSize + frameBound);
  auto compressedStringStart =
      static_cast<void*>(compressedString.data() + compressedString.size());
  compressedString.resize(compressedStringInitialSize + frameBound);
  auto compressedSize =
      LZ4F_compressFrame(compressedStringStart, frameBound, static_cast<const void*>(buffer),
                         rawDataSize, &fullFramePreferences);
  if (LZ4F_isError(compressedSize)) {
    throw std::runtime_error("LZ4 compression failed: " +
                             std::string(LZ4F_getErrorName(compressedSize)));
  }
  assert(compressedSize > 0);
  compressedString.resize(compressedStringInitialSize + compressedSize);
#else
  throw std::runtime_error("LZ4 compression not available");
#endif
}

template <typename T>
void decompressLZ4FrameToBuffer(const std::vector<char>& compressedString, MPI_Comm comm,
                                std::vector<T>& decompressedValues) {
#ifdef DISCOTEC_USE_LZ4
  LZ4F_decompressOptions_t opts;
  memset(&opts, 0, sizeof(opts));
  opts.skipChecksums = 1;

  LZ4F_dctx* dctx;
  if (LZ4F_isError(LZ4F_createDecompressionContext(&dctx, LZ4F_VERSION))) {
    throw std::runtime_error("LZ4 decompression failed ");
  }

  size_t decompressedByNow = 0;
  T* decompressValuesStart = decompressedValues.data();
  size_t readCompressedByNow = 0;
  const char* readValuesStart = compressedString.data();
  size_t decompressedSize = 0;
  do {
    decompressedSize = decompressedValues.size() * sizeof(T) - decompressedByNow;
    auto sourceSize = compressedString.size() - readCompressedByNow;
    [[maybe_unused]] auto hintNextBufferSize =
        LZ4F_decompress(dctx, decompressedValues.data(), &decompressedSize, compressedString.data(),
                        &sourceSize, &opts);
    if (LZ4F_isError(decompressedSize)) {
      throw std::runtime_error("LZ4 decompression failed: " +
                               std::string(LZ4F_getErrorName(decompressedSize)));
    }
    decompressedByNow += decompressedSize;
    decompressValuesStart += decompressedSize;
    readCompressedByNow += sourceSize;
    readValuesStart += sourceSize;
  } while (decompressedSize > 0);

  assert(decompressedByNow == decompressedValues.size() * sizeof(T));
  assert(readCompressedByNow == compressedString.size());
  [[maybe_unused]] auto freeResult = LZ4F_freeDecompressionContext(dctx);
  assert(freeResult == 0);
#endif
}

template <typename T>
void compressBufferToLZ4FrameAndGatherHeader(const T* buffer, MPI_Offset numValues,
                                             combigrid::CommunicatorType comm,
                                             std::vector<char>& compressedString) {
  size_t compressedEmptySize = 0;
  size_t compressedSkippableSize = 0;
  auto commSize = getCommSize(comm);
  auto commRank = getCommRank(comm);

#ifdef DISCOTEC_USE_LZ4
  LZ4F_preferences_t lz4PreferencesEmpty = lz4Preferences;
  lz4PreferencesEmpty.frameInfo.contentSize = 0;
  LZ4F_preferences_t lz4PreferencesSkippable = lz4Preferences;
  lz4PreferencesSkippable.frameInfo.frameType = LZ4F_skippableFrame;
  size_t locationInfoSize = commSize * sizeof(MPI_Offset);
  lz4PreferencesSkippable.frameInfo.contentSize = locationInfoSize;
  lz4PreferencesSkippable.compressionLevel = -65535; // no compression?
  auto skippableFrameBound = LZ4F_compressFrameBound(locationInfoSize, &lz4PreferencesSkippable);

  compressedString.clear();
  // if I am the root rank
  if (commRank == 0) {
    // add an empty frame first
    auto emptyFrameBound = LZ4F_compressFrameBound(0, &lz4PreferencesEmpty);
    compressedString.resize(emptyFrameBound);
    compressedEmptySize = LZ4F_compressFrame(compressedString.data(), compressedString.size(),
                                             nullptr, 0, &lz4PreferencesEmpty);
    if (LZ4F_isError(compressedEmptySize)) {
      throw std::runtime_error("LZ4 compression failed: " +
                               std::string(LZ4F_getErrorName(compressedEmptySize)));
    }

    // then add a skippable frame that contains numRanks * the location of the frame beginnings
    // (we don't yet know them, so make them placeholders)
    compressedString.resize(compressedEmptySize + skippableFrameBound);
    auto compressedStringStart = static_cast<void*>(compressedString.data() + compressedEmptySize);
    compressedSkippableSize =  // will be overwritten later
        LZ4F_compressFrame(compressedStringStart, skippableFrameBound, buffer, locationInfoSize,
                           &lz4PreferencesSkippable);
    if (LZ4F_isError(compressedSkippableSize)) {
      throw std::runtime_error("LZ4 compression failed: " +
                               std::string(LZ4F_getErrorName(compressedSkippableSize)));
    }
    assert(compressedSkippableSize > 0);
    compressedString.resize(compressedEmptySize + compressedSkippableSize);
  }

  compressBufferToLZ4Frame(buffer, numValues, comm, compressedString);

  // now the root rank collects the sizes of these frame strings
  auto compressedSize = compressedString.size();
  if (commRank == 0) {
    // gather the frame sizes of all ranks
    std::vector<MPI_Offset> frameSizes(commSize);
    MPI_Gather(&compressedSize, 1, MPI_OFFSET, frameSizes.data(), 1, MPI_OFFSET, 0, comm);

    // now overwrite the skippable frame with the actual frame sizes
    // auto compressedSkippableSize = LZ4F_compressFrameBound(locationInfoSize, &lz4Preferences);
    auto compressedSkippableStart =
        static_cast<void*>(compressedString.data() + compressedEmptySize);
    size_t newCompressedSkippableSize =
        LZ4F_compressFrame(compressedSkippableStart, skippableFrameBound, frameSizes.data(),
                           locationInfoSize, &lz4PreferencesSkippable);
    if (LZ4F_isError(newCompressedSkippableSize)) {
      throw std::runtime_error("LZ4 compression failed: " +
                               std::string(LZ4F_getErrorName(newCompressedSkippableSize)));
    }
    assert(newCompressedSkippableSize > 0);
    assert(newCompressedSkippableSize == compressedSkippableSize);
  } else {
    // send the frame sizes to the root rank
    MPI_Gather(&compressedSize, 1, MPI_OFFSET, nullptr, 0, MPI_OFFSET, 0, comm);
  }
#else
  throw std::runtime_error("LZ4 compression not available");
#endif
}

}  // namespace mpiio
}  // namespace combigrid
