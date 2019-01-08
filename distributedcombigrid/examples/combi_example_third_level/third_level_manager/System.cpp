#include "System.hpp"

System::System(const std::string& name, u_int port, AmqpClient::Channel::ptr_t channel) :
  _name(name)
{
  createMessageQueues(channel);
  createDataConnection(port, channel);
}

/*
 * Creates queues in both directions to communicate with the system.
 * Queues are non passive, non durable, non exclusive and delete themselves
 * automatically when the last using application dies.
 */
void System::createMessageQueues(AmqpClient::Channel::ptr_t channel)
{
  _inQueue = channel->DeclareQueue(_name+"_in", false, false, false, true);
  _outQueue = channel->DeclareQueue(_name+"_out", false, false, false, true);
}

void System::createDataConnection(u_int port, AmqpClient::Channel::ptr_t channel)
{
  // this might break if message arrives too early and this server hasn't made it to listen
  std::thread sendThread(&MessageUtils::sendMessage, std::to_string(port), _inQueue, channel);
  ServerSocket server(port);
  _dataChannel = std::shared_ptr<ClientSocket>(server.acceptClient());
}

void System::sendMessage(const std::string& message, AmqpClient::Channel::ptr_t channel)
{
  MessageUtils::sendMessage(message, _inQueue, channel);
}

bool System::receiveMessage(AmqpClient::Channel::ptr_t channel, std::string& message, int timeout)
{
  return MessageUtils::receiveMessage(channel, _outQueue, message, timeout);
}

std::string System::getName() const
{
  return _name;
}

std::shared_ptr<ClientSocket> System::getDataChannel() const
{
  return _dataChannel;
}
