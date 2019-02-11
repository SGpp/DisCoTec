#include "thirdLevelManager.hpp"

int main(int argc, char* argv[])
{
  // Load config file
  Params params;
  params.loadFromFile("example.ini");

  // Run third level manager
  ThirdLevelManager manager(params);
  manager.runtimeLoop();
  return 0;
}

ThirdLevelManager::ThirdLevelManager(const Params& params)
  : params_(params),
    dataPort_(params.getDataPort()),
    dataServer_(dataPort_)
{
  // Connect to RabbitMQ Broker
  channel_ = AmqpClient::Channel::Create(params.getbrokerURL());

  // Create abstraction of each system and establish message queues and data channel
  std::vector<std::string> systemNames = params.getSystemNames();
  systems_.reserve(systemNames.size());
  for (auto nameIt = systemNames.begin(); nameIt != systemNames.end(); nameIt++)
  {
    systems_.push_back(System(*nameIt, channel_, dataServer_));
  }
}
 
void ThirdLevelManager::runtimeLoop()
{
  while (systems_.size() > 0)
  {
    for (auto sysIt = systems_.begin(); sysIt != systems_.end(); sysIt++)
    {
      std::string message;
      bool received = sysIt->receiveMessage(channel_, message, timeout_);
      if (received)
        processMessage(message, *sysIt);
    }
  }
}

void ThirdLevelManager::processMessage(const std::string& message, System& system)
{
  if (message == "ready_to_combine")
    processCombination(system);
  if (message == "finished_computation")
    processFinished(system);
}

void ThirdLevelManager::processCombination(System& system)
{
  std::string message = "send_size";
  system.sendMessage(message, channel_);

  system.receiveMessage(channel_, message);
  size_t transferSize = std::stoi(message);

  message = "combine_third_level_send_first";
  system.sendMessage(message, channel_);

  for (auto sysIt = systems_.begin(); sysIt != systems_.end(); sysIt++)
  {
    if (sysIt->getName() != system.getName())
    {
      sysIt->receiveMessage(channel_, message, timeout_);
      assert(message =="ready_to_combine");
      system.sendMessage("send_size", channel_);
      size_t expectedSize = std::stoi(message);
      assert(expectedSize == transferSize);
      sysIt->sendMessage("combine_third_level_recv_first", channel_);
      // send locally combinated sparse grid from one to the other system.
      NetworkUtils::forward(system.getDataConnection().get(),
          sysIt->getDataConnection().get(), transferSize);
    }
  }
}

void ThirdLevelManager::processFinished(System& system)
{
  for (auto sysIt = systems_.begin(); sysIt != systems_.end(); sysIt++)
  {
    if (sysIt->getName() == system.getName())
      systems_.erase(sysIt);
  }
}
