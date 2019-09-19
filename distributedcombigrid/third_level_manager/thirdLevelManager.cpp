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

  // Setup third level manager
  std::cout << "Setup third level manager" << std::endl;
  ThirdLevelManager manager(params);
  manager.init();

  // Run third level manager
  std::cout << std::endl << "Running third level manager" << std::endl;
  manager.runtimeLoop();
  return 0;
}

/** Initiates the third level manager by reading in parameters and creating
 *  abstractions for each managed system.
 */
ThirdLevelManager::ThirdLevelManager(const Params& params)
  : params_(params),
    port_(params.getPort()),
    server_(port_)

{
}

void ThirdLevelManager::init()
{
  if (not params_.areLoaded())
  {
    std::cout << "ThirdLevelManager::init(): Params are not loaded!" << std::endl;
    exit(0);
  }
  assert(server_.init());
  if (not server_.isInitialized())
    exit(0);
  std::cout << "Third level manager running on port: " << port_ << std::endl;

  // Create abstraction for each system
  std::cout << "Creating abstraction for systems" << std::endl;
  for (u_int i = 0; i < params_.getNumSystems(); i++) {
    std::cout << " Waiting for system (" << i+1 << "/" << params_.getNumSystems() << ")" << std::endl;
    std::shared_ptr<ClientSocket> connection = std::shared_ptr<ClientSocket>(server_.acceptClient());
    assert(connection != nullptr && "Connecting to system failed");
    System sys(connection);
    systems_.push_back(std::move(sys));
  }
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
      bool received = systems_[s].receiveMessage(message, timeout_);
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
 *
 *  TODO optimization: forward both directions at the same time, this would
 *  allow sending back a subspace immediately after summation.
 */
void ThirdLevelManager::processCombination(size_t initiatorIndex)
{
  assert(systems_.size() == 2 && "Not implemented for different amount of systems");
  System& initiator = systems_[initiatorIndex];
  size_t otherIndex = (initiatorIndex + 1) % systems_.size();
  System& other = systems_[otherIndex];

  std::string message;
  std::cout << "Processing third level combination" << std::endl;
  initiator.sendMessage("reduce_third_level_send_first" );
  other.receiveMessage( message);
  assert(message == "ready_to_combine");
  other.sendMessage("reduce_third_level_recv_first" );

  // transfer grids from initiator to other system
  do
  {
    initiator.receiveMessage(message);
    if (message == "sending_data")
    {
      forwardData(initiator, other);
    }
  } while (message != "ready");

  // wait for other system to finish receiving
  other.receiveMessage(message);
  assert(message == "ready");

  // transfer reduced data from other system back to initiator
  do
  {
    other.receiveMessage(message);
    if (message == "sending_data")
      forwardData(other, initiator);
  } while (message != "ready");

  // wait for system to finish receiving
  initiator.receiveMessage(message);
  assert(message == "ready");
}

/** If a system has finished the simulation, it should log off the
 *  third level manager
 */
void ThirdLevelManager::processFinished(size_t sysIndex)
{
  std::cout << "System " << sysIndex << " finished simulation" << std::endl;
  systems_.erase(systems_.begin()+sysIndex);
}

/** Forwards data from sender to receiver.
 *  The size is communicated first from sender to receiver.
 */
void ThirdLevelManager::forwardData(const System& sender, const System& receiver) const {
  size_t dataSize;
  sender.receivePosNumber(dataSize);
  std::cout << "Want to forward " << std::to_string(dataSize) << " Bytes" << std::endl;
  receiver.sendMessage(std::to_string(dataSize));

  // forward data to other system
  if (dataSize != 0)
    NetworkUtils::forward(*sender.getConnection(), *receiver.getConnection(), 2048, dataSize);
}
