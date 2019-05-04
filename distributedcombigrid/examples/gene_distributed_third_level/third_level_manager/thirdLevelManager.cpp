#include "thirdLevelManager.hpp"

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cout << "Usage:\n ./thirdLevelManager paramfile" << std::endl;
    return 0;
  }
  // Load config file
  std::cout << "Loading parameters" << std::endl;
  Params params;
  params.loadFromFile(argv[1]);

  // Run third level manager
  std::cout << "Setup third level manager" << std::endl;
  ThirdLevelManager manager(params);

  std::cout << std::endl;
  std::cout << "Running third level manager" << std::endl;
  manager.runtimeLoop();
  return 0;
}

/** Initiates the third level manager by reading in parameters and creating
 *  abstractions for each managed system.
 */
ThirdLevelManager::ThirdLevelManager(const Params& params)
  : params_(params),
    dataPort_(params.getDataPort()),
    dataServer_(dataPort_)
{
  std::cout << "Data forwarding server running on port: " << dataPort_ << std::endl;
  // Connect to RabbitMQ Broker
  std::cout << "Connecting to RabbitMQ broker at " << params.getbrokerURL() << std::endl;
  channel_ = AmqpClient::Channel::Create(params.getbrokerURL());

  // Create abstraction for each system and establish message queues
  std::cout << "Creating abstraction for systems" << std::endl;
  std::vector<std::string> systemNames = params.getSystemNames();
  systems_.reserve(systemNames.size());
  for (std::string& name : systemNames)
    systems_.push_back(System(name, channel_, dataServer_));
  // establish data connection
  for (System& system : systems_)
    system.createDataConnection(dataServer_, channel_);
  std::cout << "All systems connected successfully" << std::endl;
}

/** Loops over the systems and checks if messages are available. If a message
 *  exists, it is fetched and delegated to be processed.
 */
void ThirdLevelManager::runtimeLoop()
{
  while (systems_.size() > 0)
  {
    for (size_t s = 0; s < systems_.size(); s++)
    {
      std::string message;
      bool received = systems_[s].receiveMessage(channel_, message, timeout_);
      if (received)
        processMessage(message, s);
    }
  }
}

/** Identifies main operation and initiates appropriate action. */
void ThirdLevelManager::processMessage(const std::string& message, size_t sysIndex)
{
  if (message == "ready_to_combine")
    processCombination(sysIndex);
  if (message == "finished_computation")
    processFinished(sysIndex);
}

/** Processes and manages the third level combination, initiated by a system
 *  which signals ready after its local and global combination.
 *  ATTENTION: Implemented only for 2 systems
 */
void ThirdLevelManager::processCombination(size_t initiatorIndex)
{
  assert(systems_.size() == 2 && "Not implemented for different amount of systems");
  System& initiator = systems_[initiatorIndex];
  size_t otherIndex = initiatorIndex + 1 % systems_.size();
  System& other = systems_[otherIndex];

  std::string message;
  std::cout << "Processing third level combination" << std::endl;
  initiator.sendMessage("reduce_third_level_send_first" , channel_);
  other.receiveMessage(channel_, message);
  assert(message == "ready_to_combine");
  other.sendMessage("reduce_third_level_recv_first" , channel_);

  // transfer grids from initiator to other system
  do
  {
    initiator.receiveMessage(channel_, message);
    if (message == "sending_data")
    {
      forwardData(initiator, other);
    }
  } while (message != "ready");

  // wait for other system to finish receiving
  other.receiveMessage(channel_, message);
  assert(message == "ready");

  // transfer reduced data from other system back to initiator
  do
  {
    other.receiveMessage(channel_, message);
    if (message == "sending_data")
      forwardData(other, initiator);
  } while (message != "ready");

  // wait for system to finish receiving
  initiator.receiveMessage(channel_, message);
  assert(message == "ready");
}

/** If a system has finished the simulation, it should log off the
 *  third level manager
 */
void ThirdLevelManager::processFinished(size_t sysIndex)
{
  std::cout << "System " << systems_[sysIndex].getName() << " finished simulation" << std::endl;
  systems_.erase(systems_.begin()+sysIndex);
}

/** Forwards data from sender to receiver.
 *  The size is communicated first from sender to receiver.
 */
void ThirdLevelManager::forwardData(const System& sender, const System& receiver) const {
  size_t dataSize;
  sender.receivePosNumber(channel_, dataSize);
  std::cout << "Want to forward " << std::to_string(dataSize) << " Bytes" << std::endl;
  receiver.sendMessage(std::to_string(dataSize), channel_);

  // forward data to other system
  NetworkUtils::forward(*sender.getDataConnection(), *receiver.getDataConnection(), 2048, dataSize);
}
