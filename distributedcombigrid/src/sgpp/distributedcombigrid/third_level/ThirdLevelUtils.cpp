#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"
#include <chrono>

namespace combigrid {

ThirdLevelUtils::ThirdLevelUtils(const std::string& host, int port)
                                                 : host_(host),
                                                   port_(port)
{
}

ThirdLevelUtils::~ThirdLevelUtils(){
  if (isConnected_) {
    signalFinalize();
    isConnected_ = false;
  }
}

void ThirdLevelUtils::connectToThirdLevelManager(double timeoutMinutes) {
  if (isConnected_) return;

  // create connection to third level manager
  std::cout << "Connecting to third level manager at host " << host_ << " on port " << port_
            << std::endl;
  connection_ = std::make_shared<ClientSocket>(host_, port_);
  auto t1 = std::chrono::steady_clock::now();
  auto connectionSuccess = connection_->init();
  bool timeout = false;
  while (!(connectionSuccess || timeout)) {
    connectionSuccess = connection_->init();
    auto t2 = std::chrono::steady_clock::now();
    timeout = (std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() >
               timeoutMinutes * 60);
  }

  assert(connectionSuccess && "Establishing data connection failed");
  if (!timeout) isConnected_ = true;
}

void ThirdLevelUtils::signalReadyToCombine() const
{
  sendMessage("ready_to_combine");
}

void ThirdLevelUtils::signalReadyToUnifySubspaceSizes() const
{
  sendMessage("ready_to_unify_subspace_sizes");
}

void ThirdLevelUtils::signalReadyToExchangeData() const
{
  sendMessage("ready_to_exchange_data");
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

void ThirdLevelUtils::signalFinalize() const
{
  sendMessage("finished_computation");
}

void ThirdLevelUtils::sendMessage(const std::string& message) const
{
  connection_->sendallPrefixed(message);
}

void ThirdLevelUtils::receiveMessage(std::string& message) const
{
  connection_->recvallPrefixed(message);
}

void ThirdLevelUtils::signalizeSendData() const
{
  sendMessage("sending_data");
}

void ThirdLevelUtils::sendSize(size_t size) const
{
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
}
