#include <iostream>
#include <mutex>
#include <thread>
#include <map>
#include <vector>
#include "Communicator.hpp"
#include <complex>
#include <cstdio>
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

static std::map<int, Communicator> comms;

void tunnel(const size_t port, const size_t remotePort,
    const std::string& remoteHost);

void acceptClient(const unsigned short port);

void routemsg(unsigned int source, unsigned int destination, unsigned int size);

void prepareSend(Participant& participant, size_t dest, size_t size);

void prepareReceive(Participant& participant, size_t source, size_t size);

void prepareReduceToFileUniform(size_t commID, size_t procRank, size_t numParts,
    size_t size, size_t stride, std::string datatype);

void prepareBarrier(Participant& participant);

void processReduceToFileUniform(size_t commID, size_t procRank, size_t numParts,
    size_t size, size_t stride, std::string datatype, ClientSocket* dataClient);

void  processRecv();

void sendData(SRRequest sendRequest, SRRequest recvRequest);

/*
 * Intermediary:
 * ============
 *
 * Performs different tasks depending on given arguments.
 * If executed with a single argument, the intermediary accepts connecting
 * clients on port and puts them into a communicator such that related clients
 * share the same one. By sending messages the client can gain
 * information about its communicator or communicate with other clients.
 * If executed with three arguments, connecting clients are accepted on
 * port as well, but instead a tunnel is established between the client and a
 * given remote host on port.
 */
int main(int argc, char* argv[]) {
  unsigned int port;
  if (argc > 1) {
    std::string portstr = argv[1];
    if (!NetworkUtils::isInteger(portstr)) {
      std::cout << "Port should be a positive number" << std::endl;
      return 1;
    }
    port = std::stoi(portstr);
  }
  if (argc == 2) {
    // put clients in communicator
    acceptClient(port);
  } else if (argc == 4) {
    std::string remoteHost = argv[2];
    std::string remotePortstr = argv[3];
    if (!NetworkUtils::isInteger(remotePortstr)) {
      std::cout << "remotePort should be a positive number" << std::endl;
      return 1;
    }
    unsigned short remotePort = std::stoi(remotePortstr);
    // act as tunnel
    tunnel(port, remotePort, remoteHost);
  } else {
    std::cout << "Usage:" << std::endl
              << " ./intermediary ncomms commsize port [remoteHost remotePort]"
              << std::endl;
  }
  return 0;
}

/*
 * For each new client a separate connection to the remote host is created.
 * Incoming messages from the client are directly forwarded to the remote host
 * and vice versa.
 * TODO pointers ClientSocket
 */
void tunnel(const size_t port, const size_t remotePort,
    const std::string& remoteHost) {
  ServerSocket server(port);
  server.init();
  for (;;) {
    ClientSocket* newclient = server.acceptClient();
    assert(newclient->isInitialized() && "Accepting new client failed");
    ClientSocket remote(remoteHost, remotePort);
    remote.init();
    assert(remote.isInitialized() && "Initialization of remote failed");
    std::thread sendTunnel(NetworkUtils::tunnel, newclient,
        std::ref(remote));
    std::thread recvTunnel(NetworkUtils::tunnel, std::ref(remote),
        newclient);
    sendTunnel.detach();
    recvTunnel.detach();
  }
}

/*
 * Splits a string into tokens along the delimiter c
 */
void split(const std::string& s, const char c, std::vector<std::string>& tokens) {
  tokens.clear();
  size_t prev = 0;
  while(size_t next = s.find(c) != std::string::npos) {
    tokens.push_back(s.substr(prev, next));
    prev = next + 1;
  }
}

/*
 * Adds a new client to the intermediary. The connecting client selects its
 * communicator by sending the according rank.
 * After the client has been added to its comm, future messages are processed.
 * All messages are lists of arguments separated by a '#' The first
 * argument is the id of the message, the second is the name of the operation
 * and all following its parameters.
 */
void processMessages(ClientSocket* mesgSocket) {
  std::string idstr;
  mesgSocket->recvallPrefixed(idstr);
  assert(NetworkUtils::isInteger(idstr) &&
      "Expected participant to send its id first");
  int commID = std::stoi(idstr);
  Communicator& comm = ::comms[commID];
  unsigned int procRank = comm.addParticipant(mesgSocket);
  Participant& participant = comm.getParticipant(procRank);
  std::string message;
  std::vector<std::string> argv;
  for (;;)
  {
    mesgSocket->recvallPrefixed(message);
    split(message, '#', argv);
    const size_t& argc = argv.size();
    assert(argc > 0 && "Signal message is not valid");
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
      assert(dest >= 0 && dest < comm.getSize() && size > 0 &&
          dest != commID && "Wrong parameter range for send");
      prepareSend(participant, dest, size);
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
      assert(source >= 0 && source < comm.getSize() && size > 0 &&
          source != commID && "Wrong parameter range for receive");
      prepareReceive(participant, source, size);
    }
    else if (operation == "reduceToFileUniform")
    {
      assert(argc == 6 && "Wrong number of args for reduceToFileUniform");
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
      prepareReduceToFileUniform(commID, procRank, numParts, size, stride, datatype);
    }
    else if (operation == "getCommSize")
    {
      assert(argc == 1 && "Wrong number of args for getCommSize");
      unsigned int size = comm.getSize();
      mesgSocket->sendallPrefixed(std::to_string(size));
    }
    else if (operation == "getRank")
    {
      assert(argc == 1 && "Wrong number of args for getRank");
      mesgSocket->sendallPrefixed(std::to_string(procRank));
    }
    else if (operation == "barrier#")
    {
      assert(argc == 2 && "Wrong number of args for getCommSize");
      prepareBarrier(participant);
    }
    else
    {
      //this should never happen
      assert("Unknown operation received");
    }
  }
}

void prepareSend(Participant& participant, int dest, int size) {
  const ClientSocket* mesgSocket = participant.getMesgSocket();
  // create separate connection for data transfer
  ServerSocket dataServer;
  const std::string& portstr = std::to_string(dataServer.getPort());
  mesgSocket->sendall(portstr);
  ClientSocket* dataClient = dataServer.acceptClient();
  // create new send request
  participant.pushSendRequest(dataClient, dest, size);
}

void prepareRecv(Participant& participant, int source, int size) {
  const ClientSocket* mesgSocket = participant.getMesgSocket();
  // create separate connection for data transfer
  ServerSocket dataServer;
  const std::string& portstr = std::to_string(dataServer.getPort());
  mesgSocket->sendall(portstr);
  ClientSocket* dataClient = dataServer.acceptClient();
  // create new send request
  participant.pushRecvRequest(dataClient, source, size);
}

void prepareReduceToFileUniform(size_t commID, size_t procRank, size_t numParts,
    size_t size, size_t stride, std::string datatype) {
  Participant& participant = ::comms[commID].getParticipant(procRank);
  const ClientSocket* mesgSocket = participant.getMesgSocket();
  // create separate connection for data transfer
  ServerSocket dataServer;
  const std::string& portstr = std::to_string(dataServer.getPort());
  mesgSocket->sendall(portstr);
  ClientSocket* dataClient = dataServer.acceptClient();
  // reduce data in separate thread
  std::thread reduceThread(processReduceToFileUniform, commID, procRank,
      numParts, size, stride, datatype, dataClient);
  reduceThread.detach();
}

/*
 * Follows the Producer Consumer pattern.
 * Consumes the receive requests that where produced by calls to prepareRecv()
 * by starting a new tunneling thread for each send receive pair.
 */
void processRecv(size_t procRank, Communicator& comm) {
  Participant& receiver = comm.getParticipant(procRank);
  for(;;) {
    size_t size = comm.getSize();
    for (int i = 0; i < size && i != procRank; i++) {
      Participant& sender = comm.getParticipant(i);
      if (sender.existsSendRequest(procRank)) {
        SRRequest sendRequest = sender.popSendRequest(procRank);
        SRRequest recvRequest = receiver.popRecvRequest(procRank);
        std::thread tunnelThread(sendData, sendRequest, recvRequest);
        tunnelThread.detach();
      }
    }
  }
}

/*
 * Receives the data from sender and saves it to a file.
 * The sender who contributes last initiates the reduce call where all files
 * are combined into one single solution file.
 * Thus callers do not have to wait for combination.
 */
void processReduceToFileUniform(size_t commID, size_t procRank, size_t numParts,
    size_t size, size_t stride, std::string datatype, ClientSocket* dataClient) 
{
  Communicator& comm = ::comms[commID];
  std::string filename = std::to_string(commID) + "-" + std::to_string(procRank)
    + ".part";

  // write receiving data to file
  dataClient->recvallBinaryToFile(filename, size);
  delete dataClient;

  // if last one has finished -> start the reduce operation
  if (numParts == comm.getUniformReduceCounter()) {
    std::ofstream outFile;
    outFile.open(std::to_string(commID) + ".out");

    std::vector<std::ifstream> parts;
    parts.resize(numParts);
    std::string partName;

    for (int p = 0; p < numParts; p++) {
      partName = std::to_string(commID) + "-" + std::to_string(p) + ".part";
      parts[p] = std::ifstream(partName, std::ifstream::binary);
    }

    if (datatype == "real") {
      assert(stride == sizeof(double));
      double value = .0;
      size_t numValues = size / stride;
      char* binary = new char[stride];
      // read values from files and add them up
      for (size_t i = 0; i < numValues; i++) {
        for (size_t p = 0; p < numParts; p++) {
          parts[p].read(binary, stride);
          parts[p].seekg(i*stride);
          value += *reinterpret_cast<double*>(binary);
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

    // close files
    outFile.close();
    for (size_t p = 0; p < numParts; p++) {
      parts[p].close();
      // TODO auto remove part files
      /*partName = std::to_string(commID) + "-" + std::to_string(p) + ".part";
      remove(partName);*/
    }
  } else {
    comm.increaseGetUniformReduceCounter();
  }
}

/*
 *
 */
void sendData(SRRequest sendRequest, SRRequest recvRequest) {
  assert(sendRequest.size == recvRequest.size &&
      "Wrong order of call to send/receive in client programs");
  NetworkUtils::tunnel(sendRequest.dataSock, recvRequest.dataSock);
}

/*
 * Routes a message of size Bytes from source to destination.
 * As soon as the destination is ready to receive and no other message is routed
 * from the same source simultaneously, the source is informed to begin sending.
 */
void routemesg(unsigned int source, unsigned int destination, unsigned int size){
  
}

/*
 * Continually listens on mesgPort for new clients. If a client connects, it
 * gets accepted and a new Thread is spawned for further processing.
 */
void acceptClient(const unsigned short mesgPort) {
  ServerSocket mesgServer(mesgPort);
  mesgServer.init();
  for (;;) {
    ClientSocket* mesgSocket = mesgServer.acceptClient();
    assert(mesgSocket->isInitialized() && "Accepted client corrupt");
    std::thread newThread(processMessages, mesgSocket);
    newThread.detach();
  }
}
