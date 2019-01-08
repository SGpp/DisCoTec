#include "sgpp/distributedcombigrid/third_level_old/ThirdLevelUtils.hpp"

namespace combigrid {

  ThirdLevelUtils::ThirdLevelUtils(const std::string& remoteHost,
      int remotePort) : remoteHost_(remoteHost),
  remotePort_(remotePort), mesgClient_(NULL) {}

  ThirdLevelUtils::~ThirdLevelUtils() {
    if (mesgClient_ != NULL) {
      std::string mesg = "quit";
      exchangeMesg(mesg);
    }
    delete mesgClient_;
  }

  /*
   * Connects to the intermediary and saves the third level rank.
   */
  void ThirdLevelUtils::connectToIntermediary() {
    const RankType& localRank = theMPISystem()->getLocalRank();
    mesgClient_ = new ClientSocket(remoteHost_, remotePort_);

    double start, finish;
    start = MPI_Wtime();
    mesgClient_->init();
    assert(mesgClient_->isInitialized() &&
        "Third level could not initialize message client");

    mesgClient_->sendallPrefixed("p#" + std::to_string(localRank));

    std::string mesg = "getRank";
    exchangeMesg(mesg);
    assert(NetworkUtils::isInteger(mesg));
    remoteRank_= std::stoi(mesg);
    finish = MPI_Wtime();
    std::cout << "Initialization of ThirdLevel took: " << finish-start << "Seconds";
  }

  int ThirdLevelUtils::getCommSize() {
    std::string mesg = "getCommSize";
    exchangeMesg(mesg);
    assert(NetworkUtils::isInteger(mesg) && "Comm Size must be positiv integer");
    int commSize = std::stoi(mesg);
    return commSize;
  }

  template<> bool
    ThirdLevelUtils::reduceToFileUniform(const std::vector<real>& buff) {
      size_t numParts = NUM_SYSTEMS;
      std::string mesg = "reduceToFileUniform#" + std::to_string(numParts) + "#"
        + std::to_string(buff.size()) + "#" +  std::to_string(sizeof(real)) +
        "#real";

      std::cout << "exchanging message: " << mesg << std::endl;
      exchangeMesg(mesg);

      // establish new connection for data transfer
      ClientSocket dataClient(remoteHost_, remotePort_);
      assert(dataClient.init() && "Connecting to data server failed");

      // tell intermediary to create a new data connection
      requestDataConn(dataClient);

      // start data transfer
      return NetworkUtils::sendBinary(buff, dataClient);
    }

  void ThirdLevelUtils::exchangeMesg(std::string& mesg) {
    assert(mesgClient_->isInitialized());
    std::lock_guard<std::mutex> lock(mesgLock_);

    std::cout << mesg << std::endl;
    assert(mesgClient_->sendallPrefixed(mesg) &&
        "Sending to intermediary failed");
    assert(mesgClient_->recvallPrefixed(mesg) &&
        "Receiving from intermediary failed");
  }

  void ThirdLevelUtils::barrier(size_t commSize) {
    std::string mesg = "barrier#" + std::to_string(commSize);
    exchangeMesg(mesg);
    assert(NetworkUtils::isInteger(mesg) && "Received data port is not valid");
    int dataPort = std::stoi(mesg);
    ClientSocket dataClient(remoteHost_, dataPort);
    assert(dataClient.init() && "Connecting to data server failed");
    char* signal = nullptr;
    // blocks until receives signal
    dataClient.recvall(signal, 1);
    delete[] signal;
  }

  void ThirdLevelUtils::requestDataConn(const ClientSocket& dataClient) const {
      const RankType& localRank = theMPISystem()->getLocalRank();
      std::string initmsg = "d#" + std::to_string(localRank) +  "#" +
                            std::to_string(remoteRank_);
      dataClient.sendallPrefixed(initmsg);
  }

  int ThirdLevelUtils::getRank() {
    return remoteRank_;
  }
}
