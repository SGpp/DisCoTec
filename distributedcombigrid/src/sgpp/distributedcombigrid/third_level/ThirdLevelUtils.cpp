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

  // create message queues.
  // Queues are non passive, non durable, non exclusive and delete themselves
  // automatically when the last using application dies.
  inQueue_ = messageChannel_->DeclareQueue(systemName_+"_in", false, false, false, true);
  outQueue_ = messageChannel_->DeclareQueue(systemName_+"_out", false, false, false, true);

  // create data connection
  dataConnection_ = std::make_shared<ClientSocket>(remoteHost_, dataPort_);
  assert(dataConnection_->init() && "Establishing data connection failed");

  isConnected_ = true;
}

void ThirdLevelUtils::signalReady() const
{
  sendMessage("ready");
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
