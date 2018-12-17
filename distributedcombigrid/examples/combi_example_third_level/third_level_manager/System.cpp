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

void System::sendMessage(const std::string& message,
    AmqpClient::Channel::ptr_t channel) const
{
  AmqpClient::BasicMessage::ptr_t mesg = AmqpClient::BasicMessage::Create(message);
  channel->BasicPublish("", _inQueue, mesg);
}

bool System::receiveMessage(AmqpClient::Channel::ptr_t channel,
    std::string& message, u_int timeout ) const
{
  std::string consumer_tag = channel->BasicConsume(_outQueue, "");
  AmqpClient::Envelope::ptr_t envelope;
  bool received = channel->BasicConsumeMessage(consumer_tag, envelope, timeout);
  if (!received)
  {
    return false;
  }
  else
  {
    message = envelope->Message()->Body();
    return true;
  }
}

void System::createDataConnection(u_int port, AmqpClient::Channel::ptr_t channel)
{
  ServerSocket server(port);
  // this might break if message arrives too early and server hasn't made it to listen
  std::thread sendThread(&System::sendMessage, this, std::to_string(port), channel);
  _dataChannel = std::unique_ptr<ClientSocket>(server.acceptClient());
}
