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

#include "io/MPIInputOutput.hpp"
#include "mpi/MPISystem.hpp"
#include "utils/Types.hpp"

namespace combigrid {

namespace mpiio {

#ifdef DISCOTEC_USE_LZ4
constexpr LZ4F_preferences_t lz4Preferences = {
    {LZ4F_max4MB, LZ4F_blockLinked, LZ4F_noContentChecksum, LZ4F_frame,
     0,   // unknown content size
     0U,  // no dictID
     LZ4F_noBlockChecksum},
    0,   // fast mode
    0u,  // don't flush
    0u,  // unused
    {0u, 0u, 0u}};
#endif

template <typename T>
void compressBufferToLZ4Frame(const T* buffer, MPI_Offset numValues,
                              std::vector<char>& compressedString) {
  size_t compressedStringInitialSize = compressedString.size();
#ifdef DISCOTEC_USE_LZ4
  size_t rawDataSize = numValues * sizeof(T);
  auto fullFramePreferences = lz4Preferences;
  fullFramePreferences.frameInfo.contentSize = rawDataSize;
  auto dataBound = LZ4_compressBound(static_cast<int>(rawDataSize));
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

#ifdef DISCOTEC_USE_LZ4
template <typename T>
class [[nodiscard]] FrameDecompressor {
 public:
  explicit FrameDecompressor(const std::vector<char>& compressedString, T* decompressValuesStart)
      : compressedString_(compressedString), decompressValuesStart_(decompressValuesStart) {
    memset(&opts_, 0, sizeof(opts_));
    opts_.skipChecksums = 1;
    if (LZ4F_isError(LZ4F_createDecompressionContext(&dctx_, LZ4F_VERSION))) {
      throw std::runtime_error("LZ4 decompression failed ");
    }
  }

  // no copy or move
  FrameDecompressor(const FrameDecompressor&) = delete;
  FrameDecompressor(FrameDecompressor&&) = delete;
  FrameDecompressor& operator=(const FrameDecompressor&) = delete;
  FrameDecompressor& operator=(FrameDecompressor&&) = delete;

  ~FrameDecompressor() {
    // assert that everything has been read by now
    assert(readCompressedByNow_ == compressedString_.size());
    [[maybe_unused]] auto freeResult = LZ4F_freeDecompressionContext(dctx_);
    assert(freeResult == 0);
  }

  size_t decompressNextBlocks(size_t numValues) {
    size_t decompressedSize = numValues * sizeof(T);
    auto sourceSize = compressedString_.size() - readCompressedByNow_;
    [[maybe_unused]] auto hintNextBufferSize =
        LZ4F_decompress(dctx_, decompressValuesStart_, &decompressedSize, compressedString_.data(),
                        &sourceSize, &opts_);
    if (LZ4F_isError(decompressedSize)) {
      throw std::runtime_error("LZ4 decompression failed: " +
                               std::string(LZ4F_getErrorName(decompressedSize)));
    }
    decompressedByNow_ += decompressedSize;
    decompressValuesStart_ += decompressedSize;
    readCompressedByNow_ += sourceSize;
    return decompressedSize / sizeof(T);
  }

 private:
  const std::vector<char>& compressedString_;
  T* decompressValuesStart_;
  LZ4F_decompressOptions_t opts_;
  LZ4F_dctx* dctx_;
  size_t decompressedByNow_ = 0;
  size_t readCompressedByNow_ = 0;
};
#endif

template <typename T>
size_t decompressLZ4FrameToBuffer(const std::vector<char>& compressedString,
                                  T* decompressValuesStart, size_t numValues) {
  size_t decompressedSize = 0;
  size_t decompressedByNow = 0;
#ifdef DISCOTEC_USE_LZ4
  auto decompressor = FrameDecompressor(compressedString, decompressValuesStart);
  do {
    decompressedSize = decompressor.decompressNextBlocks(numValues);
    numValues -= decompressedSize;
    decompressedByNow += decompressedSize;
  } while (decompressedSize > 0);
  assert(numValues == 0);

  // https://github.com/lz4/lz4/blob/dev/examples/frameCompress.c //TODO make streaming
#else
  throw std::runtime_error("LZ4 compression not available");
#endif
  return decompressedByNow;
}

template <typename T>
size_t decompressLZ4FrameToBuffer(const std::vector<char>& compressedString,
                                  std::vector<T>& decompressedValues) {
  auto decompressedSize = decompressLZ4FrameToBuffer(compressedString, decompressedValues.data(),
                                                     decompressedValues.size());
  assert(decompressedSize == decompressedValues.size());
  return decompressedSize;
}

#ifdef DISCOTEC_USE_LZ4
inline size_t getSizeOfHeaderFrame(int commSize, size_t& headerFrameBound,
                                   LZ4F_preferences_t& lz4PreferencesHeader) {
  lz4PreferencesHeader = lz4Preferences;
  size_t locationInfoSize = commSize * sizeof(MPI_Offset);
  lz4PreferencesHeader.frameInfo.contentSize = locationInfoSize;
  lz4PreferencesHeader.compressionLevel = -65535;  // no compression?
  headerFrameBound = LZ4F_compressFrameBound(locationInfoSize, &lz4PreferencesHeader);
  // initialize w/ random values
  std::vector<char> compressedDummyString(headerFrameBound);

  size_t compressedHeaderSize =  // will be overwritten later
      LZ4F_compressFrame(compressedDummyString.data(), headerFrameBound,
                         compressedDummyString.data(), locationInfoSize, &lz4PreferencesHeader);
  if (LZ4F_isError(compressedHeaderSize)) {
    throw std::runtime_error("LZ4 dummy compression failed: " +
                             std::string(LZ4F_getErrorName(compressedHeaderSize)));
  }
  assert(compressedHeaderSize > 0);
  return compressedHeaderSize;
}
#endif

template <typename T>
void compressBufferToLZ4FrameAndGatherHeader(const T* buffer, MPI_Offset numValues,
                                             combigrid::CommunicatorType comm,
                                             std::vector<char>& compressedString) {
  auto commSize = getCommSize(comm);
  auto commRank = getCommRank(comm);

#ifdef DISCOTEC_USE_LZ4
  LZ4F_preferences_t lz4PreferencesHeader;
  size_t headerFrameBound;
  auto compressedHeaderSize =
      getSizeOfHeaderFrame(commSize, headerFrameBound, lz4PreferencesHeader);
  compressedString.clear();
  // if I am the root rank
  if (commRank == 0) {
    // add a header frame that contains numRanks * the location of the frame beginnings
    // (we don't yet know them, so make them placeholders)
    compressedString.resize(compressedHeaderSize);
  }

  compressBufferToLZ4Frame(buffer, numValues, compressedString);

  // now the root rank collects the sizes of these frame strings
  auto compressedSize = compressedString.size();
  if (commRank == 0) {
    // gather the frame sizes of all ranks
    std::vector<MPI_Offset> frameSizes(commSize);
    MPI_Gather(&compressedSize, 1, MPI_OFFSET, frameSizes.data(), 1, MPI_OFFSET, 0, comm);

    // now overwrite the header frame with the actual frame sizes
    size_t newCompressedHeaderSize =
        LZ4F_compressFrame(compressedString.data(), headerFrameBound, frameSizes.data(),
                           commSize * sizeof(MPI_Offset), &lz4PreferencesHeader);
    if (LZ4F_isError(newCompressedHeaderSize)) {
      throw std::runtime_error("LZ4 compression failed: " +
                               std::string(LZ4F_getErrorName(newCompressedHeaderSize)));
    }
    assert(newCompressedHeaderSize > 0);
    assert(newCompressedHeaderSize == compressedHeaderSize);
  } else {
    // send the frame sizes to the root rank
    MPI_Gather(&compressedSize, 1, MPI_OFFSET, nullptr, 0, MPI_OFFSET, 0, comm);
  }
#else
  throw std::runtime_error("LZ4 compression not available");
#endif
}

template <typename T>
int writeCompressedValuesConsecutive(const T* valuesStart, MPI_Offset numValues,
                                     const std::string& fileName, combigrid::CommunicatorType comm,
                                     bool replaceExistingFile = false,
                                     bool withCollectiveBuffering = false) {
  std::vector<char> compressedString;
  mpiio::compressBufferToLZ4FrameAndGatherHeader(valuesStart, numValues, comm, compressedString);
  return mpiio::writeValuesConsecutive<char>(compressedString.data(), compressedString.size(),
                                             fileName, comm, true, false);
}

#ifdef DISCOTEC_USE_LZ4
inline void getFrameSizeAndPosFromHeader(const MPIFileConsecutive<char>& file,
                                         MPI_Offset& frameSize, MPI_Offset& position,
                                         combigrid::CommunicatorType comm) {
  auto rank = getCommRank(comm);
  // rank 0 scatters positions read from header frame
  if (rank == 0) {
    auto commSize = getCommSize(comm);
    LZ4F_preferences_t lz4PreferencesHeader = lz4Preferences;
    size_t headerFrameBound;
    auto compressedHeaderSize =
        getSizeOfHeaderFrame(commSize, headerFrameBound, lz4PreferencesHeader);

    std::vector<char> headerFrame(compressedHeaderSize);
    auto numCharsRead =
        file.readValuesFromFileAtPositionSingleRank(headerFrame.data(), headerFrame.size(), 0);
    assert(numCharsRead == static_cast<int>(compressedHeaderSize));
    headerFrame.resize(numCharsRead);

    std::vector<MPI_Offset> frameSizes(commSize);
    [[maybe_unused]] auto numCharsInFirstFrame =
        decompressLZ4FrameToBuffer(headerFrame, frameSizes);
    MPI_Scatter(frameSizes.data(), 1, MPI_OFFSET, &frameSize, 1, MPI_OFFSET, 0, comm);
    assert(frameSize >= numCharsRead);
    file.checkFileSizeConsecutive(frameSize, comm);
    position = mpiio::getPositionFromNumValues(frameSize, comm);
    assert(position == 0);

    // the content frame for rank 0 starts after the header frame
    frameSize -= compressedHeaderSize;
    position += compressedHeaderSize;
  } else {
    MPI_Scatter(nullptr, 0, MPI_OFFSET, &frameSize, 1, MPI_OFFSET, 0, comm);
    file.checkFileSizeConsecutive(frameSize, comm);
    position = mpiio::getPositionFromNumValues(frameSize, comm);
  }
}
#endif

// cf. readValuesConsecutive
template <typename T>
int readCompressedValuesConsecutive(T* valuesStart, MPI_Offset numValues,
                                    const std::string& fileName, combigrid::CommunicatorType comm) {
  auto file = MPIFileConsecutive<char>::getFileToRead(fileName, comm, false, false);
  MPI_Offset position = 0;
  MPI_Offset frameSize = 0;
  getFrameSizeAndPosFromHeader(file, frameSize, position, comm);

  std::vector<char> frameString(frameSize);
  auto numCharsRead = file.readValuesFromFileAtPosition(frameString.data(), frameSize, position);
  assert(numCharsRead == frameSize);
  frameString.resize(numCharsRead);

  auto numValuesDecompressed = static_cast<int>(
      mpiio::decompressLZ4FrameToBuffer(std::move(frameString), valuesStart, numValues));
  assert(static_cast<MPI_Offset>(numValuesDecompressed) == numValues);
  return numValuesDecompressed;
}

template <typename T, typename ReduceFunctionType>
int readReduceCompressedValuesConsecutive(
    T* valuesStart, MPI_Offset numValues, const std::string& fileName,
    combigrid::CommunicatorType comm,  // int numElementsToBuffer,
    ReduceFunctionType reduceFunction) {
  int numValuesDecompressed = 0;
#ifdef DISCOTEC_USE_LZ4
  // // numElementsToBuffer must be at least lz4 block size
  // assert(numElementsToBuffer >= (1 << 20) * 4);
  auto file = MPIFileConsecutive<char>::getFileToRead(fileName, comm, false, false);
  MPI_Offset position = 0;
  MPI_Offset frameSize = 0;
  getFrameSizeAndPosFromHeader(file, frameSize, position, comm);

  std::vector<char> frameString(frameSize);
  auto numCharsRead = file.readValuesFromFileAtPosition(frameString.data(), frameSize, position);
  assert(numCharsRead == frameSize);
  frameString.resize(numCharsRead);

  int reduceCount = 0;
  auto remainingValues = numValues;
  size_t decompressedSize = 0;
  std::vector<T> buffer(numValues);

  numValuesDecompressed = static_cast<int>(
      mpiio::decompressLZ4FrameToBuffer(std::move(frameString), buffer.data(), numValues));
  assert(static_cast<MPI_Offset>(numValuesDecompressed) == numValues);

  std::transform(buffer.cbegin(), buffer.cend(), valuesStart, valuesStart, reduceFunction);

#else
  throw std::runtime_error("LZ4 compression not available");
#endif
  return numValuesDecompressed;
}

}  // namespace mpiio
}  // namespace combigrid
