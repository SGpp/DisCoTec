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

    public:
      ThirdLevelUtils(const std::string& remoteHost, int dataPort,
                      const std::string& systemName,
                      int brokerPort = RabbitMQPort);

      void connectToThirdLevelManager();

      void signalReadyToCombine() const;

      void sendSize(size_t size) const;

      std::string fetchInstruction() const;

      template<typename FG_ELEMENT>
      void receiveCommonSSPart(std::vector<FG_ELEMENT>& commonSSPart) const;

      template<typename FG_ELEMENT>
      void sendCommonSSPart(const std::vector<FG_ELEMENT>& commonSSPart) const;
  };

  template <typename FG_ELEMENT>
  void ThirdLevelUtils::receiveCommonSSPart(std::vector<FG_ELEMENT>& commonSSPart) const
  {
    bool dataIsLittleEndian = false;
    bool success = dataConnection_->recvallBinary(commonSSPart, dataIsLittleEndian);
    assert(success && "receiving common ss data failed");
    if (dataIsLittleEndian != NetworkUtils::isLittleEndian())
      NetworkUtils::reverseEndianness(commonSSPart);
  }

  // TODO: check if inplace is better
  template <typename FG_ELEMENT>
  void ThirdLevelUtils::sendCommonSSPart(const std::vector<FG_ELEMENT>& commonSSPart) const
  {
    bool success = dataConnection_->sendallBinary(commonSSPart);
    assert(success && "sending common ss data failed");
    std::cout << "completed sending common ss part" << std::endl;
  }
}

#endif
