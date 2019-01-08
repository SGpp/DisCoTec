#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "MessageUtils.hpp"
#include <string>
#include <thread>

class System
{
  private:
    std::string                   _name;
    std::string                   _inQueue;     // used for sending messages to the system
    std::string                   _outQueue;    // used for receiving messages from the system
    std::shared_ptr<ClientSocket> _dataChannel; // data connection

    void createMessageQueues(AmqpClient::Channel::ptr_t channel);
    void createDataConnection(u_int port, AmqpClient::Channel::ptr_t channel);

  public:

    System(const std::string& name, u_int port, AmqpClient::Channel::ptr_t channel);

    void sendMessage(const std::string& message, AmqpClient::Channel::ptr_t channel);
    bool receiveMessage(AmqpClient::Channel::ptr_t channel, std::string& message, int timeout);

    std::shared_ptr<ClientSocket> getDataChannel() const;
    std::string                   getName() const;
};
