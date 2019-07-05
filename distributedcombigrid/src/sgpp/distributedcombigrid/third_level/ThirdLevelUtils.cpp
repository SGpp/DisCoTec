#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"


ThirdLevelUtils::ThirdLevelUtils(const std::string& remoteHost, int dataPort,
                                 const std::string& systemName,
                                 int brokerPort) : remoteHost_(remoteHost),
                                                   brokerPort_(brokerPort),
                                                   dataPort_(dataPort),
                                                   systemName_(systemName)
{

}

ThirdLevelUtils::~ThirdLevelUtils(){
  signalFinalize();
  isConnected_ = false;
}

void ThirdLevelUtils::connectToThirdLevelManager()
{
  if (isConnected_)
    return;
  // connect to message broker
  std::cout << "Connecting to RabbitMQ broker at host " << remoteHost_ << " on port " << brokerPort_ << std::endl;
  messageChannel_ = AmqpClient::Channel::Create(remoteHost_);

  inQueue_ = MessageUtils::createMessageQueue(systemName_+"_in", messageChannel_);
  outQueue_ = MessageUtils::createMessageQueue(systemName_+"_out", messageChannel_);

  consumerTag_ = MessageUtils::setupConsumer(messageChannel_, inQueue_);

  // create data connection
  std::string message;
  std::cout << "Connecting to ThirdLevel manager at host " << remoteHost_ << " on port " << dataPort_  << std::endl;
  receiveMessage(message);
  std::cout <<  systemName_ << " received message " << message << std::endl;
  assert(message == "create_data_conn");
  dataConnection_ = std::make_shared<ClientSocket>(remoteHost_, dataPort_);
  assert(dataConnection_->init() && "Establishing data connection failed");

  isConnected_ = true;
}

void ThirdLevelUtils::signalReadyToCombine() const
{
  sendMessage("ready_to_combine");
}

void ThirdLevelUtils::signalReady() const
{
  sendMessage("ready");
}

std::string ThirdLevelUtils::fetchInstruction() const
{
  std::string instruction;
  receiveMessage(instruction);
  std::cout << "Fetched instruction: " << instruction << std::endl;
  return instruction;
}

void ThirdLevelUtils::signalFinalize() const
{
  sendMessage("finished_computation");
}

void ThirdLevelUtils::sendMessage(const std::string& message) const
{
  assert(messageChannel_ != nullptr);
  MessageUtils::sendMessage(message, outQueue_, messageChannel_);
}

void ThirdLevelUtils::receiveMessage(std::string& message) const
{
  assert(messageChannel_ != nullptr);
  MessageUtils::receiveMessage(messageChannel_, consumerTag_, inQueue_, message);
}

void ThirdLevelUtils::sendSize(size_t size) const
{
  sendMessage("sending_data");
  sendMessage(std::to_string(size));
}

size_t ThirdLevelUtils::receiveSize() const
{
  std::string sizeStr;
  receiveMessage(sizeStr);

  std::stringstream ss(sizeStr);
  size_t size;
  ss >> size;
  return size;
}
