#include "thirdLevelManager.hpp"
#include <boost/asio.hpp>

#if BROKER_ON_SYSTEM
  // to resolve https://github.com/open-mpi/ompi/issues/5157
  #define OMPI_SKIP_MPICXX 1
  #include <mpi.h>
#endif // BROKER_ON_SYSTEM

using namespace combigrid;

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cout << "Usage:\n ./thirdLevelManager paramfile \
                  or ./thirdLevelManager --port=9999 --numSystems=2 --chunksize=131072" << std::endl;
    return 0;
  }

#if BROKER_ON_SYSTEM
  // if broker is running on highest rank on same system, need to split it away from world
  // communicator
  MPI_Init(&argc, &argv);
  int color = 1;
  int key = 0;
  MPI_Comm worldComm;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &worldComm);
  int size = 0;
  MPI_Comm_size(worldComm, &size);
  if(size != 1) {
    throw std::runtime_error("Broker is not running on a single process");
  }
  std::cout << "thirdLevelManager running on same system as simulation" << std::endl;
#else
  std::cout << "thirdLevelManager running on separate system" << std::endl;
#endif  // BROKER_ON_SYSTEM

  std::string hostnameInfo = "broker = " + boost::asio::ip::host_name();
  std::cout << hostnameInfo << std::endl;

  // Load config file
#ifdef DEBUG_OUTPUT
  std::cout << "Loading parameters" << std::endl;
#endif
  Params params;
  if (argv[1][0] == '-') {
    params.loadFromCmd(argc, argv);
  } else {
    params.loadFromFile(argv[1]);
  }

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

  // Print statistics
  manager.writeStatistics("stats.json");
  return 0;
}

namespace combigrid {

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
  auto serverInitSuccess = server_.init();
  assert(serverInitSuccess);
  if (not server_.isInitialized())
  {
    std::cout << "ThirdLevelManager::init(): Could not initialize server!" << std::endl;
    exit(0);
  }
  std::cout << "Third level manager running on port: " << port_ << std::endl;

  // Create abstraction for each system
#ifdef DEBUG_OUTPUT
  std::cout << "Creating abstraction for systems" << std::endl;
#endif
  for (u_int i = 0; i < params_.getNumSystems(); i++) {
    std::cout << " Waiting for system (" << i+1 << "/" << params_.getNumSystems() << ")" << std::endl;
    std::shared_ptr<ClientSocket> connection = std::shared_ptr<ClientSocket>(server_.acceptClient());
    assert(connection != nullptr && "Connecting to system failed");
    System sys(connection, i+1);
    systems_.push_back(std::move(sys));
  }
  std::cout << "All systems connected successfully" << std::endl;
}


/** Loops over the systems and checks if messages are available. If a message
 *  exists, it is fetched and delegated to be processed.
 */
void ThirdLevelManager::runtimeLoop()
{
  stats_.startWallclock();
  while (systems_.size() > 0)
  {
    for (size_t s = 0; s < systems_.size(); s++)
    {
      std::string message;
      bool hasMessage = systems_[s].hasMessage(timeout_);
      if (hasMessage) {
        auto success = systems_[s].receiveMessage(message);
        if (!success) {
          throw std::runtime_error("ThirdLevelManager::runtimeLoop(): Receiving message failed!");
        }
#ifdef DEBUG_OUTPUT
        std::cout << "received message from system " << std::to_string(s) << ": " << message
                  << std::endl;
#endif
        processMessage(message, s);
      }
    }
  }
  stats_.stopWallclock();
}

/** Identifies main operation and initiates appropriate action. */
void ThirdLevelManager::processMessage(const std::string& message, size_t sysIndex)
{
  if (message == "ready_to_combine")
    processCombination(sysIndex);
  else if (message == "ready_to_combine_file")
    processCombinationFile(sysIndex);
  else if (message == "ready_to_unify_subspace_sizes")
    processUnifySubspaceSizes(sysIndex);
  else if (message == "ready_to_exchange_data")
    processAnyData(sysIndex);
  else if (message == "finished_computation")
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
  stats_.increaseNumCombinations();
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
    size_t dataSize = forwardData(initiator, other);
    stats_.addToBytesTransferredInCombination(dataSize);

    other.receiveMessage(message);
    assert(message == "sending_data");
    dataSize = forwardData(other, initiator);
    stats_.addToBytesTransferredInCombination(dataSize);
  }

  // wait for other system to finish receiving
  other.receiveMessage(message);
  assert(message == "ready");
#ifdef DEBUG_OUTPUT
  std::cout << "Finished combination" << std::endl;
#endif
}


void ThirdLevelManager::processCombinationFile(size_t initiatorIndex)
{
  assert(systems_.size() == 2 && "Not implemented for different amount of systems");
  stats_.increaseNumCombinations();
  System& initiator = systems_[initiatorIndex];
  size_t otherIndex = (initiatorIndex + 1) % systems_.size();
  System& other = systems_[otherIndex];

  std::string message;
#ifdef DEBUG_OUTPUT
  std::cout << std::endl << "Processing third level combination w files" << std::endl;
#endif
  other.receiveMessage(message);
  assert(message == "ready_to_combine_file");

  initiator.receiveMessage(message);
  assert(message == "ready");
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
    size_t dataSize = forwardData(initiator, other);
    stats_.addToBytesTransferredInSizeExchange(dataSize);
  }
  // transfer data from other to initiator
  other.receiveMessage(message);
  if (message == "sending_data")
  {
    size_t dataSize = forwardData(other, initiator);
    stats_.addToBytesTransferredInSizeExchange(dataSize);
  }
  initiator.receiveMessage(message);
  assert(message == "ready");
  other.receiveMessage(message);
  assert(message == "ready");
#ifdef DEBUG_OUTPUT
  std::cout << "Finished unification of subspace sizes" << std::endl;
#endif
}


void ThirdLevelManager::processAnyData(size_t initiatorIndex)
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
  assert(message == "ready_to_exchange_data");
  other.sendMessage("recv_first" );

  // transfer data from initiator (sends first) to other (receives first)
  initiator.receiveMessage(message);
  if (message == "sending_data")
  {
    // size_t dataSize = forwardData(initiator, other);
    forwardData(initiator, other);
  }
  // transfer data from other to initiator
  other.receiveMessage(message);
  if (message == "sending_data")
  {
    // size_t dataSize = forwardData(other, initiator);
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
 *  The size is communicated first from sender to receiver and returned later.
 */
size_t ThirdLevelManager::forwardData(const System& sender, const System& receiver) const
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
    NetworkUtils::forward(*sender.getConnection(), *receiver.getConnection(),
                          this->params_.getChunksize(), dataSize);

  return dataSize;
}

void ThirdLevelManager::writeStatistics(std::string filename)
{
  size_t walltime = stats_.getWallTime();
  size_t numCombinations = stats_.getNumCombinations();
  size_t totalBytesTransferredInCombination = stats_.getTotalBytesTransferredInCombination();
  size_t numBytesTransferredPerCombination = totalBytesTransferredInCombination / numCombinations;
  size_t totalBytesTransferredInSizeExchange = stats_.getTotalBytesTransferredInSizeExchange();
  if (filename != "")
  {
    std::ofstream ofs (filename, std::ofstream::out);
    boost::property_tree::ptree pt;
    pt.put("walltime (ms)", walltime/1e3);
    pt.put("numCombinations", numCombinations);
    pt.put("totalBytesTransferredInCombination", totalBytesTransferredInCombination);
    pt.put("numBytesTransferredPerCombination", numBytesTransferredPerCombination);
    pt.put("totalBytesTransferredInSizeExchange", totalBytesTransferredInSizeExchange);
    boost::property_tree::write_json(ofs, pt);
    ofs.close();
  }
  std::cout << "Simulation took:                     "
    << (double) walltime / 1e6 << "sec" << std::endl;
  std::cout << "Num combinations:                    "
    << numCombinations << std::endl;
  std::cout << "Total transfer during combination:   "
    << totalBytesTransferredInCombination << "B" << std::endl;
  std::cout << "Transfer per combination:            "
    << numBytesTransferredPerCombination << "B" << std::endl;
  std::cout << "Total transfer during size exchange: "
    << totalBytesTransferredInSizeExchange << "B" << std::endl;
}

}
