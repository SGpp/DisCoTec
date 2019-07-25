#include "System.hpp"

System::System(const std::string& name,
               const AmqpClient::Channel::ptr_t& channel,
               const ServerSocket& server) : name_(name)
{
  createMessageQueues(channel);
  consumerTag_ = MessageUtils::setupConsumer(channel, outQueue_);
}

/*
 * Creates queues in both directions to communicate with the system.
 * Queues are non passive, non durable, non exclusive and delete themselves
 * automatically when the last using application dies.
 */
void System::createMessageQueues(AmqpClient::Channel::ptr_t channel)
{
  inQueue_ = MessageUtils::createMessageQueue(name_+"_in", channel);
  outQueue_ = MessageUtils::createMessageQueue(name_+"_out", channel);
}

void System::createDataConnection(const ServerSocket& server,
                                  const AmqpClient::Channel::ptr_t& channel)
{
  std::cout << "creating data connection for system: " << name_ << std::endl;
  std::string message = "create_data_conn";
  sendMessage(message, channel);
  dataConnection_ = std::shared_ptr<ClientSocket>(server.acceptClient());
}

void System::sendMessage(const std::string& message, AmqpClient::Channel::ptr_t channel) const
{
  MessageUtils::sendMessage(message, inQueue_, channel);
}

bool System::receiveMessage(AmqpClient::Channel::ptr_t channel, std::string& message, int timeout) const
{
  if (MessageUtils::receiveMessage(channel, consumerTag_, outQueue_, message))
  {
    std::cout << "Received message: \"" << message << "\" from System " << name_ << std::endl;
    return true;
  }
  return false;
}

bool System::receivePosNumber(AmqpClient::Channel::ptr_t channel, size_t& number, int timeout) const
{
  std::string message;
  if (MessageUtils::receiveMessage(channel, consumerTag_, outQueue_, message))
  {
    assert(NetworkUtils::isInteger(message));
    number = std::stoul(message);
    return true;
  }
  return false;
}

std::string System::getName() const
{
  return name_;
}

std::shared_ptr<ClientSocket> System::getDataConnection() const
{
  return dataConnection_;
}
