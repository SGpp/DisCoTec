#include "MessageUtils.hpp"

std::string MessageUtils::createMessageQueue(const std::string& name,
                                             const AmqpClient::Channel::ptr_t channel)
{
  // create a message queue.
  // Queue is non passive, non durable, non exclusive and delete itself
  // automatically when the last using application dies.
  return channel->DeclareQueue(name, false, false, false, true);
}

bool MessageUtils::receiveMessage(const AmqpClient::Channel::ptr_t channel,
                                  const std::string                queue,
                                  std::string&                     message,
                                  int                              timeout)
{
  std::string consumer_tag = channel->BasicConsume(queue, "");
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

void MessageUtils::sendMessage(const std::string&               message,
                               const std::string                queue,
                               const AmqpClient::Channel::ptr_t channel)
{
  AmqpClient::BasicMessage::ptr_t mesg = AmqpClient::BasicMessage::Create(message);
  channel->BasicPublish("", queue, mesg);
}
