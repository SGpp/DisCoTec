#include "SimpleAmqpClient/SimpleAmqpClient.h"
#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include <string>
#include <thread>

class System
{
  private:
    const std::string _name;
    std::string _inQueue;
    std::string _outQueue;
    std::unique_ptr<ClientSocket> _dataChannel;

    void createMessageQueues(AmqpClient::Channel::ptr_t channel);
    void createDataConnection(u_int port, AmqpClient::Channel::ptr_t channel);

  public:
    System(const std::string& name, u_int port, AmqpClient::Channel::ptr_t channel);

    void sendMessage(const std::string& message,
        AmqpClient::Channel::ptr_t channel)                      const;

    bool receiveMessage(AmqpClient::Channel::ptr_t channel,
        std::string& message, u_int timeout )                    const;

};
