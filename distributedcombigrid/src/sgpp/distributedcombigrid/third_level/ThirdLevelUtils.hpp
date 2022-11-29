#ifndef THIRDLEVELUTILSHPP_
#define THIRDLEVELUTILSHPP_

#include <stdlib.h>
#include <ctime>
#include <sstream>
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"

namespace combigrid {

  class ThirdLevelUtils
  {
    private:
      std::string host_;
      int port_;
      std::shared_ptr<ClientSocket> connection_;
      bool isConnected_ = false;

      void connectToIntermediary();

      void receiveMessage(std::string& message) const;

      void sendMessage(const std::string& message) const;

      void signalFinalize() const;

      void signalizeSendData() const;

      void sendSize(size_t size) const;

    public:
      ThirdLevelUtils(const std::string& host, int port);

      ~ThirdLevelUtils();

      void connectToThirdLevelManager(double timeoutMinutes);

      void signalReadyToCombine() const;

      void signalReadyToCombineFile() const;

      void signalReadyToUnifySubspaceSizes() const;

      void signalReadyToExchangeData() const;

      void signalReady() const;

      size_t receiveSize() const;

      std::string fetchInstruction() const;

      template <typename FG_ELEMENT>
      DistributedSparseGridUniform<FG_ELEMENT> recvDSGUniform() const;

      std::string recvDSGUniformSerialized() const;

      /** Sends the given data in binary form to the third level manager
       */
      template <typename FG_ELEMENT>
      void sendData(const FG_ELEMENT* const data, size_t size) const;

      /** Receives binary data from the third level manager.
       */
      template <typename FG_ELEMENT>
      void recvData(FG_ELEMENT* data, size_t size) const;

      /** Receives upcoming data from the third level manager and directly adds
       * it to the provided buffer. The implementation does not wait with the
       * summation until all remote data is received. Instead, each received
       * chunk is added directly to the provided buffer. Therefor the memory
       * overhead is usually much lower (chunksize can be adjusted) compared to
       * the classical way.
       */
      template <typename FG_ELEMENT>
      void recvAndAddToData(FG_ELEMENT* data, size_t size) const;
  };


  template <typename FG_ELEMENT>
  void ThirdLevelUtils::sendData(const FG_ELEMENT* data, size_t size) const
  {
    assert(isConnected_);
    signalizeSendData();
    size_t rawSize = size * sizeof(FG_ELEMENT) + 1; // + 1 due to endianness
    sendSize(rawSize);
    connection_->sendallBinary(data, size);
  }

  template <typename FG_ELEMENT>
  void ThirdLevelUtils::recvData(FG_ELEMENT* data, size_t size) const
  {
    assert(isConnected_);
    size_t rawSize = receiveSize();
    size_t recvSize = (rawSize - 1) / sizeof(FG_ELEMENT); // - 1 due to endianness
    assert(recvSize == size && "Size mismatch receiving data size does not match expected");
    bool success = connection_->recvallBinaryAndCorrectInPlace(data, size);
    assert(success && "receiving dsgu data failed");
  }


  template <typename FG_ELEMENT>
  void ThirdLevelUtils::recvAndAddToData(FG_ELEMENT* data, size_t size) const
  {
    assert(isConnected_);
    size_t rawSize = receiveSize();
    size_t recvSize = (rawSize - 1) / sizeof(FG_ELEMENT); // - 1 due to endianness
    assert(recvSize == size && "Size mismatch cannot add vectors of different size");
    bool success = connection_->recvallBinaryAndReduceInPlace<FG_ELEMENT>(data, size,
        [](const FG_ELEMENT& lhs, const FG_ELEMENT& rhs) -> FG_ELEMENT {return lhs + rhs;});
    assert(success && "receiving dsgu data failed");
  }
}
#endif
