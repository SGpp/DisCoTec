#ifndef THIRDLEVELUTILSHPP_
#define THIRDLEVELUTILSHPP_

#include <stdlib.h>
#include <ctime>
#include <sstream>
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "sgpp/distributedcombigrid/third_level/MessageUtils.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include "SimpleAmqpClient/SimpleAmqpClient.h"

namespace combigrid {

  static const int RabbitMQPort = 5672;

  class ThirdLevelUtils
  {
    private:
      std::string remoteHost_;
      int brokerPort_;
      int dataPort_;
      AmqpClient::Channel::ptr_t messageChannel_;
      std::shared_ptr<ClientSocket> dataConnection_;
      std::string inQueue_;
      std::string outQueue_;
      std::string consumerTag_;
      std::string systemName_;
      bool isConnected_ = false;

      void connectToIntermediary();

      void receiveMessage(std::string& message) const;

      void sendMessage(const std::string& message) const;

      void signalFinalize() const;

      void sendSize(size_t size) const;

    public:
      ThirdLevelUtils(const std::string& remoteHost, int dataPort,
                      const std::string& systemName,
                      int brokerPort = RabbitMQPort);

      ~ThirdLevelUtils();

      void connectToThirdLevelManager();

      void signalReadyToCombine() const;

      void signalReady() const;

      size_t receiveSize() const;

      std::string fetchInstruction() const;

      template <typename FG_ELEMENT>
      DistributedSparseGridUniform<FG_ELEMENT> recvDSGUniform() const;

      std::string recvDSGUniformSerialized() const;

      template <typename FG_ELEMENT>
      void sendDSGUniform(DistributedSparseGridUniform<FG_ELEMENT>& dsgu) const;

      void sendDSGUniformSerialized(const std::string& serializedDSGU) const;
  };

  template <typename FG_ELEMENT>
  DistributedSparseGridUniform<FG_ELEMENT> ThirdLevelUtils::recvDSGUniform() const
  {
    size_t rawSize = receiveSize();
    std::stringstream ss(recvDSGUniformSerialized());
    DistributedSparseGridUniform<FG_ELEMENT> dsgu;
    {
      boost::archive::text_iarchive ia(ss);
      // read class state from archive
      ia >> dsgu;
    }
    return dsgu;
  }

  template <typename FG_ELEMENT>
  void ThirdLevelUtils::sendDSGUniform(DistributedSparseGridUniform<FG_ELEMENT>& dsgu) const
  {
    std::stringstream ss;
    {
      boost::archive::text_oarchive oa(ss);
      // read class state from archive
      oa << dsgu;
    }
    std::string serializedDSGU(ss.str());
    sendSize(serializedDSGU.size());
    bool success = dataConnection_->sendall(serializedDSGU);
    assert(success && "sending common ss data failed");
  }
}

#endif
