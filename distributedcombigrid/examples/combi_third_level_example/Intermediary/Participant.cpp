#include "Participant.hpp"

Participant::Participant(Communicator& comm, ClientSocket* mesgSocket,
    size_t rank) : comm_(comm), mesgSocket_(mesgSocket), rank_(rank)
{
}

Participant::~Participant() {
  delete mesgSocket_;
}

/*
 * After the client has been added to its comm, future messages are processed.
 * All messages are strings that represent operations the client wants to
 * perform on the intermediary. The arguments of the operation are suffixed and 
 * separated by a #. Thus a message is of the form operation#arg1#...#argn
 */
void Participant::processMessages() {
  std::string message;
  std::vector<std::string> argv;
  for(;;)
  {
    if (mesgSocket_->isReadable()) {
      mesgSocket_->recvallPrefixed(message);
      std::cout << "splitting message..." << std::endl;
      NetworkUtils::split(message, '#', argv);
      std::cout << "received message: " << message << std::endl;
      const size_t& argc = argv.size();
      assert(argc > 0 && "message is not valid");
      std::string& operation = argv[0];

      if (operation == "send")
      {
        assert(argc == 3 && "Wrong number of args for send");
        std::string& deststr= argv[1];
        std::string& sizestr = argv[2];
        assert(NetworkUtils::isInteger(deststr) &&
            "destination of send message must be a number");
        int dest = std::stoi(deststr);
        assert(NetworkUtils::isInteger(sizestr) &&
            "size of send message must be a number");
        int size = std::stoi(sizestr);
        assert(dest >= 0 && comm_.existsParticipant(dest) && size > 0 &&
            dest != rank_ && "Wrong parameter range for send");
        mesgSocket_->sendallPrefixed("ok");
        prepareSend(dest, size);
      }
      else if (operation == "recv")
      {
        assert(argc == 3 && "Wrong number of args for receive");
        std::string& srcstr= argv[1];
        std::string& sizestr = argv[2];
        assert(NetworkUtils::isInteger(srcstr) &&
            "source of receive message must be a number");
        int source = std::stoi(srcstr);
        assert(NetworkUtils::isInteger(sizestr) &&
            "size of receive message must be a number");
        int size = std::stoi(sizestr);
        assert(source >= 0 && comm_.existsParticipant(source) && size > 0 &&
            source != rank_ && "Wrong parameter range for receive");
        mesgSocket_->sendallPrefixed("ok");
        prepareRecv(source, size);
      }
      else if (operation == "reduceToFileUniform")
      {
        assert(argc == 5 && "Wrong number of args for reduceToFileUniform");
        std::string& numPartsstr = argv[1];
        std::string& sizestr = argv[2];
        std::string& stridestr = argv[3];
        std::string& datatype = argv[4];
        assert(NetworkUtils::isInteger(numPartsstr) &&
            "received numParts must be a number");
        size_t numParts = std::stoi(numPartsstr);
        assert(NetworkUtils::isInteger(sizestr) &&
            "received size must be a number");
        size_t size = std::stoi(sizestr);
        assert(NetworkUtils::isInteger(stridestr) &&
            "received stride must be a number");
        size_t stride = std::stoi(stridestr);
        mesgSocket_->sendallPrefixed("ok");
        prepareReduceToFileUniform(numParts, size, stride, datatype);
      }
      else if (operation == "getCommSize")
      {
        assert(argc == 1 && "Wrong number of args for getCommSize");
        unsigned int size = comm_.getSize();
        mesgSocket_->sendallPrefixed(std::to_string(size));
      }
      else if (operation == "getRank")
      {
        assert(argc == 1 && "Wrong number of args for getRank");
        mesgSocket_->sendallPrefixed(std::to_string(rank_));
      }
      else if (operation == "barrier")
      {
        assert(argc == 1 && "Wrong number of args for barrier");
        //prepareBarrier();
      }
      else if (operation == "quit") {
        mesgSocket_->sendallPrefixed("ok");
        // waits until all spawned threads finish
        while(numThreads_ != 0)
          std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        return;
      }
      else
      {
        //this should never happen
        assert("Unknown operation received");
      }
    }
  }
}

/*
 *
 */
void Participant::prepareSend(int dest, size_t size) {
  // fetch separate connection for data transfer
  ClientSocket* dataConn = fetchNewConnection();
  SRRequest sendRequest = { dataConn, size };

  // deligate processing of request to communicator
  bool state = comm_.processSendRequest(rank_, dest, sendRequest);
  if (!state) {
    delete dataConn;
    assert(state && "PrepareSend failed: Sender or receiver does not exist");
  }
}

/*
 *
 */
void Participant::prepareRecv(int source, size_t size) {
  // fetch separate connection for data transfer
  ClientSocket* dataConn = fetchNewConnection();
  SRRequest recvRequest = {dataConn, size};

  // deligate processing of request to communicator
  bool state = comm_.processRecvRequest(source, rank_, recvRequest);
  if (!state) {
    delete dataConn;
    assert(state && "PrepareRecv failed: Sender or receiver does not exist");
  }
}

void Participant::prepareReduceToFileUniform(size_t numParts, size_t size,
    size_t stride, std::string datatype)
{
  // fetch separate connection for data transfer
  ClientSocket* dataConn = fetchNewConnection();

  // reduce data in separate thread
  std::thread reduceThread(&Participant::processReduceToFileUniform, this, numParts, size, stride, datatype, dataConn);
  reduceThread.detach();
  increaseThreadCounter();
}

/*
 * Waits until Intermediary accepts a new connection designated
 * for the Participant and returns a pointer to the connection.
 */
ClientSocket* Participant::fetchNewConnection() {
  std::unique_lock<std::mutex> lock(newConnLock_);
  while(!newConnValid_) newConnCV_.wait(lock);
  ClientSocket* newConn = newConn_;
  newConnValid_ = false;
  return newConn;
}

/*
 * Receives the data from sender and saves it to a file.
 * The sender who contributes last, initiates the reduce call where all files
 * are combined into one single solution file.
 * Callers do not have to wait for combination.
 */
void Participant::processReduceToFileUniform(size_t numParts,
    size_t size, size_t stride, std::string datatype, ClientSocket* dataClient) 
{
  size_t commID = comm_.getID();
  std::string filename = std::to_string(commID) + "-" + std::to_string(rank_)
    + ".part";

  // write receiving data to file
  dataClient->recvallBinaryToFile(filename, size * stride);
  //delete dataClient;
  size_t counter = comm_.increaseUniformReduceCounter();

  // if last one has finished -> start the reduce operation
  if (counter == numParts) {
    std::cout << "start combining parts..." << std::endl;
    std::ofstream outFile;
    outFile.open(std::to_string(commID) + ".out");

    // open all files written so far
    std::vector<std::ifstream> parts;
    parts.resize(numParts);
    std::string partName;
    for (int p = 0; p < numParts; p++) {
      partName = std::to_string(commID) + "-" + std::to_string(p) + ".part";
      parts[p] = std::ifstream(partName, std::ifstream::binary);
    }

    if (datatype == "real") {
      assert(stride == sizeof(double));
      double value = 0.;
      size_t numValues = size;
      char* binary = new char[stride];

      // read values from part files and add them up
      std::cout << "Reading " << numValues <<" values with " << stride << " bytes each." << std::endl;
      for (size_t i = 0; i < numValues; i++) {
        for (size_t p = 0; p < numParts; p++) {
          parts[p].seekg(i*stride);
          parts[p].read(binary, stride);
          value += *reinterpret_cast<double*>(binary);
          std::cout << "Value " << i << " in Part " << p << ": " << *reinterpret_cast<double*>(binary) << std::endl;
        }
        // write string representation of value into solution file
        outFile << value << "\n";
      }
    } else if (datatype == "complex") {
      assert(stride == sizeof(std::complex<double>));
      assert("complex not yet implemented");
      // TODO
    } else {
      assert("datatype not implemented");
    }

    // close all files
    for (size_t p = 0; p < numParts; p++) {
      parts[p].close();
      // TODO auto remove part files
      /*partName = std::to_string(commID) + "-" + std::to_string(p) + ".part";
      remove(partName);*/
    }
    outFile.close();
    std::cout << "finished combining parts..." << std::endl;
  }
}

size_t Participant::getRank() const {
  return rank_;
}

const ClientSocket& Participant::getMesgSocket() const {
  return *mesgSocket_;
}

void Participant::increaseThreadCounter() {
  std::lock_guard<std::mutex> lck(numThreadsLock_);
  numThreads_++;
}

void Participant::decreaseThreadCounter() {
  std::lock_guard<std::mutex> lck(numThreadsLock_);
  numThreads_--;
}
