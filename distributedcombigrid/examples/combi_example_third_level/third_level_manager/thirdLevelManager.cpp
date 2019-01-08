#include "thirdLevelManager.hpp"

int main(int argc, char* argv[])
{
  // Load config file
  Params params("example.ini");

  // Run third level manager
  ThirdLevelManager manager(params);
  manager.runtimeLoop();
  return 0;
}

ThirdLevelManager::ThirdLevelManager(const Params& params) : _params(params)
{
  // Connect to RabbitMQ Broker
  _channel = AmqpClient::Channel::Create(params.getbrokerURL());

  // Create abstraction of each system and establish message queues and data channel
  std::vector<std::string> systemNames = params.getSystemNames();
  _systems.reserve(systemNames.size());
  for (auto nameIt = systemNames.begin(); nameIt != systemNames.end(); nameIt++)
  {
    _systems.push_back(System(*nameIt, _port, _channel));
  }
}



void ThirdLevelManager::runtimeLoop()
{
  while (_systems.size() > 0)
  {
    for (auto sysIt = _systems.begin(); sysIt != _systems.end(); sysIt++)
    {
      std::string message;
      bool received = sysIt->receiveMessage(_channel, message, _timeout);
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
  std::string message = "do_send_size";
  system.sendMessage(message, _channel);
  system.receiveMessage(_channel, message, noTimeout);
  u_int combiGridSize = std::stoi(message);

  message = "do_send_data";
  system.sendMessage(message, _channel);

  message = "receive_data";
  for (auto sysIt = _systems.begin(); sysIt != _systems.end(); sysIt++)
  {
    if (sysIt->getName() != system.getName())
    {
      sysIt->sendMessage(message, _channel);
      sysIt->receiveMessage(_channel, message, noTimeout);
      assert(message =="ok");
      // send common subspaces from one to the other system.
      NetworkUtils::forward(system.getDataChannel().get(),
          sysIt->getDataChannel().get(), combiGridSize);
    }
  }
}

void ThirdLevelManager::processFinished(System& system)
{
  for (auto sysIt = _systems.begin(); sysIt != _systems.end(); sysIt++)
  {
    if (sysIt->getName() == system.getName())
      _systems.erase(sysIt);
  }
}
