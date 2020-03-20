#include "thirdLevelManager.hpp"

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cout << "Usage:\n ./thirdLevelManager paramfile" << std::endl;
    return 0;
  }
  // Load config file
#ifdef DEBUG_OUTPUT
  std::cout << "Loading parameters" << std::endl;
#endif
  Params params;
  params.loadFromFile(argv[1]);

  // Setup third level manager
#ifdef DEBUG_OUTPUT
  std::cout << "Setup third level manager" << std::endl;
#endif
  ThirdLevelManager manager(params);
  manager.init();

  // Run third level manager
#ifdef DEBUG_OUTPUT
  std::cout << "Running third level manager" << std::endl;
#endif
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
  {
    std::cout << "ThirdLevelManager::init(): Could not initialize server!" << std::endl;
    exit(0);
  }
#ifdef DEBUG_OUTPUT
  std::cout << "Third level manager running on port: " << port_ << std::endl;
#endif

  // Create abstraction for each system
#ifdef DEBUG_OUTPUT
  std::cout << "Creating abstraction for systems" << std::endl;
#endif
  for (u_int i = 0; i < params_.getNumSystems(); i++) {
#ifdef DEBUG_OUTPUT
    std::cout << " Waiting for system (" << i+1 << "/" << params_.getNumSystems() << ")" << std::endl;
#endif
    std::shared_ptr<ClientSocket> connection = std::shared_ptr<ClientSocket>(server_.acceptClient());
    assert(connection != nullptr && "Connecting to system failed");
    System sys(connection, i+1);
    systems_.push_back(std::move(sys));
  }
#ifdef DEBUG_OUTPUT
  std::cout << "All systems connected successfully" << std::endl;
#endif
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
      bool hasMessage = systems_[s].hasMessage(timeout_);
      if (hasMessage) {
        systems_[s].receiveMessage(message);
        processMessage(message, s);
      }
    }
  }
}

/** Identifies main operation and initiates appropriate action. */
void ThirdLevelManager::processMessage(const std::string& message, size_t sysIndex)
{
  if (message == "ready_to_combine")
    processCombination(sysIndex);
  if (message == "ready_to_unify_subspace_sizes")
    processUnifySubspaceSizes(sysIndex);
  if (message == "finished_computation")
    processFinished(sysIndex);
}

/** Processes and manages the third level combination, initiated by a system
 *  which signals ready after its local and global combination.
 *  ATTENTION: Implemented only for 2 systems
 *
 *  TODO optimization: forward both directions at the same time
 */
void ThirdLevelManager::processCombination(size_t initiatorIndex)
{
  assert(systems_.size() == 2 && "Not implemented for different amount of systems");
  System& initiator = systems_[initiatorIndex];
  size_t otherIndex = (initiatorIndex + 1) % systems_.size();
  System& other = systems_[otherIndex];

  std::string message;
#ifdef DEBUG_OUTPUT
  std::cout << std::endl << "Processing third level combination" << std::endl;
#endif
  initiator.sendMessage("send_first");
  other.receiveMessage(message);
  assert(message == "ready_to_combine");
  other.sendMessage("recv_first");

  // transfer grids between the systems
  while (initiator.receiveMessage(message) && message != "ready")
  {
    assert(message == "sending_data");
    forwardData(initiator, other);

    other.receiveMessage(message);
    assert(message == "sending_data");
    forwardData(other, initiator);
  }

  // wait for other system to finish receiving
  other.receiveMessage(message);
  assert(message == "ready");
#ifdef DEBUG_OUTPUT
  std::cout << "Finished combination" << std::endl;
#endif
}

/* TODO optimization: forward both directions at the same time
 *
 */
void ThirdLevelManager::processUnifySubspaceSizes(size_t initiatorIndex)
{
  assert(systems_.size() == 2 && "Not implemented for different number of systems");
  System& initiator = systems_[initiatorIndex];
  size_t otherIndex = (initiatorIndex + 1) % systems_.size();
  System& other = systems_[otherIndex];

  std::string message;
#ifdef DEBUG_OUTPUT
  std::cout << std::endl << "Processing unification of subspace sizes" << std::endl;
#endif
  initiator.sendMessage("send_first" );
  other.receiveMessage( message);
  assert(message == "ready_to_unify_subspace_sizes");
  other.sendMessage("recv_first" );

  // transfer data from initiator (sends first) to other (receives first)
  initiator.receiveMessage(message);
  if (message == "sending_data")
  {
    forwardData(initiator, other);
  }
  // transfer data from other to initiator
  other.receiveMessage(message);
  if (message == "sending_data")
  {
    forwardData(other, initiator);
  }
  initiator.receiveMessage(message);
  assert(message == "ready");
  other.receiveMessage(message);
  assert(message == "ready");
#ifdef DEBUG_OUTPUT
  std::cout << "Finished unification of subspace sizes" << std::endl;
#endif
}

/** If a system has finished the simulation, it should log off the
 *  third level manager
 */
void ThirdLevelManager::processFinished(size_t sysIndex)
{
#ifdef DEBUG_OUTPUT
  std::cout << "System " << systems_[sysIndex].getId() << " finished simulation" << std::endl;
#endif
  systems_.erase(systems_.begin()+ (long)sysIndex);
}

/** Forwards data from sender to receiver.
 *  The size is communicated first from sender to receiver.
 */
void ThirdLevelManager::forwardData(const System& sender, const System& receiver) const
{
  size_t dataSize;
  sender.receivePosNumber(dataSize);
#ifdef DEBUG_OUTPUT
  std::cout << "Forwarding " << std::to_string(dataSize) << " Bytes from system "
            << sender.getId() << " to system " << receiver.getId() << std::endl;
#endif
  receiver.sendMessage(std::to_string(dataSize));

  // forward data to other system
  if (dataSize != 0)
    NetworkUtils::forward(*sender.getConnection(), *receiver.getConnection(), 2048, dataSize);
}
