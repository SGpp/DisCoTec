#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "sgpp/distributedcombigrid/third_level/MessageUtils.hpp"
#include <string>
#include <thread>

class System
{
  private:
    std::string                   name_;
    std::string                   inQueue_;     // used for sending messages to the system
    std::string                   outQueue_;    // used for receiving messages from the system
    std::shared_ptr<ClientSocket> dataConnection_; // data connection

    void createMessageQueues(AmqpClient::Channel::ptr_t channel);


  public:

    System(const std::string& name, const AmqpClient::Channel::ptr_t& channel,
           const ServerSocket& server);

    void sendMessage(const std::string& message, AmqpClient::Channel::ptr_t channel);
    bool receiveMessage(AmqpClient::Channel::ptr_t channel, std::string& message, int timeout=MessageUtils::noTimeout);

    std::shared_ptr<ClientSocket> getDataConnection() const;
    std::string                   getName() const;

    void createDataConnection(const ServerSocket& server,
                              const AmqpClient::Channel::ptr_t& channel);
};
