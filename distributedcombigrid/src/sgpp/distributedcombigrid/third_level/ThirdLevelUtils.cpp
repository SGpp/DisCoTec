#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"


ThirdLevelUtils::ThirdLevelUtils(const std::string& remoteHost, int dataPort,
                                 const std::string& systemName,
                                 int brokerPort) : remoteHost_(remoteHost),
                                                   brokerPort_(brokerPort),
                                                   dataPort_(dataPort),
                                                   systemName_(systemName)
{

}

void ThirdLevelUtils::connectToThirdLevelManager()
{
  if (isConnected_)
    return;
  // connect to message broker
  messageChannel_ = AmqpClient::Channel::Create(remoteHost_);

  MessageUtils::createMessageQueue(systemName_+"_in", messageChannel_);
  MessageUtils::createMessageQueue(systemName_+"_out", messageChannel_);

  // create data connection
  dataConnection_ = std::make_shared<ClientSocket>(remoteHost_, dataPort_);
  assert(dataConnection_->init() && "Establishing data connection failed");

  isConnected_ = true;
}

void ThirdLevelUtils::signalReadyToCombine() const
{
  sendMessage("ready_to_combine");
}

std::string ThirdLevelUtils::fetchInstruction() const
{
  std::string instruction;
  receiveMessage(instruction);
  return instruction;
}

void ThirdLevelUtils::sendMessage(const std::string& message) const
{
  MessageUtils::sendMessage(message, outQueue_, messageChannel_);
}

void ThirdLevelUtils::receiveMessage(std::string& message) const
{
  MessageUtils::receiveMessage(messageChannel_, inQueue_, message);
}

void ThirdLevelUtils::sendSize(size_t size) const
{
  assert(isConnected_);
  sendMessage(std::to_string(size));
}
