#include "thirdLevelManager.hpp"

int main(int argc, char* argv[])
{
  // Load config file
  std::cout << "Loading parameters" << std::endl;
  Params params;
  params.loadFromFile("example.ini");

  // Run third level manager
  std::cout << "Setup third level manager" << std::endl;
  ThirdLevelManager manager(params);

  std::cout << std::endl;
  std::cout << "Running third level manager" << std::endl;
  manager.runtimeLoop();
  return 0;
}

/*
 * Initiates the third level manager by reading in parameters and creating
 * abstractions for each managed system.
 */
ThirdLevelManager::ThirdLevelManager(const Params& params)
  : params_(params),
    dataPort_(params.getDataPort()),
    dataServer_(dataPort_)
{
  // Connect to RabbitMQ Broker
  std::cout << "Connecting to RabbitMQ broker" << std::endl;
  channel_ = AmqpClient::Channel::Create(params.getbrokerURL());

  // Create abstraction for each system and establish message queues
  std::cout << "Creating abstraction for systems" << std::endl;
  std::vector<std::string> systemNames = params.getSystemNames();
  systems_.reserve(systemNames.size());
  for (auto nameIt = systemNames.begin(); nameIt != systemNames.end(); nameIt++)
  {
    systems_.push_back(System(*nameIt, channel_, dataServer_));
  }
  // establish data connection
  for (auto system = systems_.begin(); system != systems_.end(); system++)
    system->createDataConnection(dataServer_, channel_);
}

/*
 * Loops over the systems and checks if messages are available. If a message
 * exists, it is fetched and delegated to be processed.
 */
void ThirdLevelManager::runtimeLoop()
{
  while (systems_.size() > 0)
  {
    for (auto sysIt = systems_.begin(); sysIt != systems_.end() && systems_.size() > 0; sysIt++)
    {
      std::string message;
      bool received = sysIt->receiveMessage(channel_, message, timeout_);
      if (received)
        processMessage(message, *sysIt);
    }
  }
}

/*
 * Identifies main operation and initiates appropriate action.
 */
void ThirdLevelManager::processMessage(const std::string& message, System& system)
{
  if (message == "ready_to_combine")
    processCombination(system);
  if (message == "finished_computation")
    processFinished(system);
}

/*
 * Processes and manages the third level combination, initiated by a system
 * which signals ready after his local and global combination.
 * ATTENTION: Implemented only for 2 systems
 */
void ThirdLevelManager::processCombination(System& initiator)
{
  assert(systems_.size() == 2 && "Not implemented for different amount of systems");

  // TODO maybe cleaner to pass id of initiator
  System& other = systems_[0];
  for (System& sys : systems_)
  {
    if (sys.getName() != initiator.getName())
    {
      other = sys;
      break;
    }
  }

  std::cout << "Processing third level combination" << std::endl;
  std::string message = "reduce_third_level_send_first";
  initiator.sendMessage(message, channel_);

  // transfer data from initiator to other system
  do
  {
    initiator.receiveMessage(channel_, message);
    if (message == "sending_data")
      forwardData(initiator, other);
  } while (message != "ready");

  // transfer reduced data from other system back to initiator
  do
  {
    other.receiveMessage(channel_, message);
    if (message == "sending_data")
      forwardData(other, initiator);
  } while (message != "ready");
}

/*
 * If a system has finished the simulation, it should log of from the
 * third level manager
 */
void ThirdLevelManager::processFinished(System& system)
{
  for (auto sysIt = systems_.begin(); sysIt != systems_.end(); sysIt++)
  {
    if (sysIt->getName() == system.getName())
    {
      systems_.erase(sysIt);
      break;
    }
  }
}

/*
 * Forwards data from sender to receiver
*/
void ThirdLevelManager::forwardData(System& sender, System& receiver) const {
  // receive size of serialized grid  data
  size_t dataSize;
  sender.receivePosNumber(channel_, dataSize);
  // send size of serialized grid to other system
  receiver.sendMessage(std::to_string(dataSize), channel_);

  // forward data to other system
  NetworkUtils::forward(*sender.getDataConnection(), *receiver.getDataConnection(), 2048, dataSize);
}
