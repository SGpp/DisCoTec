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

  // Establish message queues and data channel for each system
  std::vector<std::string> systemNames = params.getSystemNames();
  _systems.reserve(systemNames.size());
  for (auto nameIt = systemNames.begin(); nameIt != systemNames.end(); nameIt++)
  {
    _systems.push_back(System(*nameIt, 9999, _channel));
  }
}

void ThirdLevelManager::runtimeLoop()
{
  for (;;)
  {
    for (auto sysIt = _systems.begin(); sysIt != _systems.end(); sysIt++)
    {
      std::string message;
      bool received = sysIt->receiveMessage(_channel, message, timeout);
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
  system.receiveMessage(_channel, message, -1);
  u_int combiGridSize = std::stoi(message);

  message = "do_send_data";
  system.sendMessage(message, _channel);
}

void ThirdLevelManager::processFinished(System& system)
{
  
}
