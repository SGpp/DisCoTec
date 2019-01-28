#include "System.hpp"

System::System(const std::string& name,
               const AmqpClient::Channel::ptr_t& channel,
               const ServerSocket& server) : _name(name)
{
  createMessageQueues(channel);
  createDataConnection(server, channel);
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

// TODO initialization which ensures that remote server accepts.
void System::createDataConnection(const ServerSocket& server,
                                  const AmqpClient::Channel::ptr_t& channel)
{
  std::string message;
  receiveMessage(channel, message, -1); // waits until system is ready
  assert(message == "create_data_conn");
  _dataConnection = std::shared_ptr<ClientSocket>(server.acceptClient());
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

std::shared_ptr<ClientSocket> System::getDataConnection() const
{
  return _dataConnection;
}
