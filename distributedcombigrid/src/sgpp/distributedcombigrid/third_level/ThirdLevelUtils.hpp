#ifndef THIRDLEVELUTILSHPP_
#define THIRDLEVELUTILSHPP_

#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include "SimpleAmqpClient/SimpleAmqpClient.h"
#include "MessageUtils.hpp"
#include <stdlib.h>
#include <ctime>
#include <sstream>

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
      std::string systemName_;
      bool isConnected_ = false;

      void createDataConnection();

      void connectToIntermediary();

      void receiveMessage(std::string& message) const;

      void sendMessage(const std::string& message) const;

    public:
      ThirdLevelUtils(const std::string& remoteHost, int dataPort,
                      const std::string& systemName,
                      int brokerPort = RabbitMQPort);

      void connectToThirdLevelManager();

      void signalReady() const;

      std::string fetchInstruction() const;

      template<typename FG_ELEMENT>
      void receiveCommonSubspaces(std::vector<FG_ELEMENT>& commonSubspaces) const;

      template<typename FG_ELEMENT>
      void sendCommonSubspaces(const std::vector<FG_ELEMENT>& commonSubspaces) const;
  };

}

#endif
