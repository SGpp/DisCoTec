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

      void connectToThirdLevelManager();

      void signalReadyToCombine() const;

      void signalReadyToUnifySubspaceSizes() const;

      void signalReady() const;

      size_t receiveSize() const;

      std::string fetchInstruction() const;

      template <typename FG_ELEMENT>
      DistributedSparseGridUniform<FG_ELEMENT> recvDSGUniform() const;

      std::string recvDSGUniformSerialized() const;

      /** Sends the given data to the third level manager
       */
      template <typename FG_ELEMENT>
      void sendData(const FG_ELEMENT* const data, size_t size) const;

      /** Receives upcoming data from the third level manager, returns the data
       * and its size.
       */
      template <typename FG_ELEMENT>
      void recvData(FG_ELEMENT* data, size_t& size) const;
  };


  template <typename FG_ELEMENT>
  void ThirdLevelUtils::sendData(const FG_ELEMENT* data, size_t size) const
  {
    assert(isConnected_);
    size_t rawSize = size * sizeof(FG_ELEMENT);
    signalizeSendData();
    sendSize(rawSize);
    if (size != 0) {
      //std::cout << "Manager tries to send " << rawSize << " Bytes" << std::endl;
      connection_->sendall(reinterpret_cast<const char*>(data), rawSize);
    }
  }

  template <typename FG_ELEMENT>
  void ThirdLevelUtils::recvData(FG_ELEMENT* data, size_t& size) const
  {
    assert(isConnected_);
    size_t rawSize = receiveSize();
    size = rawSize / sizeof(FG_ELEMENT);
    if(size != 0) {
      char* rawData = new char[rawSize];
      //std::cout << "Manager tries to receive " << rawSize << " Bytes" << std::endl;
      bool success = connection_->recvall(rawData, rawSize);
      assert(success && "receiving dsgu data failed");
      data = reinterpret_cast<FG_ELEMENT*>(rawData);
    }
  }
}

#endif
