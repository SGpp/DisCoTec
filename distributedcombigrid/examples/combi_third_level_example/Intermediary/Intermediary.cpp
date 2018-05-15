#include <iostream>
#include <mutex>
#include <thread>
#include <map>
#include <vector>
#include "Communicator.hpp"
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

std::map<int, Communicator> comms;

void tunnel(const size_t port, const size_t remotePort,
    const std::string& remoteHost);

void acceptClient(const unsigned short port);

void routemsg(unsigned int source, unsigned int destination, unsigned int size);

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
 */
void tunnel(const size_t port, const size_t remotePort,
    const std::string& remoteHost) {
  ServerSocket server(port);
  server.init();
  for (;;) {
    ClientSocket newclient = server.acceptClient(); // TODO return client and pass succes by ref
    assert(newclient.isInitialized() && "Accepting new client failed");
    ClientSocket remote(remoteHost, remotePort);
    remote.init();
    assert(remote.isInitialized() && "Initialization of remote failed");
    std::thread sendTunnel(NetworkUtils::tunnel, std::ref(newclient),
        std::ref(remote));
    std::thread recvTunnel(NetworkUtils::tunnel, std::ref(remote),
        std::ref(newclient));
    sendTunnel.join();
    recvTunnel.join();
  }
}


/*
 * Splits a string into tokens along the delimiter c
 */
std::vector<std::string> split(const std::string& s, char c) {
  std::vector<std::string> tokens;
  size_t prev = 0;
  while(size_t next = s.find(c) != std::string::npos) {
    tokens.push_back(s.substr(prev, next));
    prev = next + 1;
  }
  return tokens;
}
/*
 * Adds a new client to the intermediary. The connecting client selects its
 * communicator by sending the according rank in his first message.
 * Next future messages are parsed.
 * All future messages are lists of arguments separated by a '#' The first
 * argument is the id of the message, the second is the name of the operation
 * and all following its parameters.
 */
void processMessages(Participant participant) {
  std::string idstr;
  participant.mesgSocket.recvallPrefixed(idstr);
  assert(NetworkUtils::isInteger(idstr) &&
      "Expected participant to send its id first");
  int id = std::stoi(idstr);
  Communicator& comm = ::comms[id];
  unsigned int procRank = ::comms[id].addParticipant(participant);
  std::string message;
  for (;;) {
    participant.mesgSocket.recvallPrefixed(message);
    std::vector<std::string> argv = NetworkUtils::split(message, '#'); // TODO might be better to pass args by ref
    const size_t& argc = argv.size();
    assert(argc > 0 && "Signal message is not valid");
    std::string& mesgIDStr = argv[0];
    assert(NetworkUtils::isInteger(mesgIDStr) &&
        "message ID must be a number");
    long mesgID = std::stoi(mesgIDStr);
    std::string& operation = argv[1];
    if (operation == "send") {
      assert(argc == 4 && "Wrong number of args for send");
      std::string& deststr= argv[2];
      std::string& sizestr = argv[3];
      assert(NetworkUtils::isInteger(deststr) &&
          "destination of send message must be a number");
      unsigned int dest = std::stoi(deststr);
      assert(NetworkUtils::isInteger(sizestr) &&
          "size of send message must be a number");
      unsigned int size = std::stoi(sizestr);
      routemsg(procRank, dest, size);
    }
    else if (operation == "getCommSize") {
      assert(argc == 2 && "Wrong number of args for getCommSize");
      unsigned int size = comm.getSize();
      participant.mesgSocket.sendallPrefixed(std::to_string(size));
    }
    else if (operation == "getRank") {
      assert(argc == 2 && "Wrong number of args for getRank");
      participant.mesgSocket.sendallPrefixed(std::to_string(procRank));
    }
    else {
      assert("Unknown operation received");
    }
  }
}

/*
 * routes a message of size from source to destination.
 * As soon as the destination is ready to receive and no other message is routed
 * from the same source simultaneously, the source is informed to begin sending.
 */
void routemesg(unsigned int source, unsigned int destination, unsigned int size){
  
}


/*
 * Continually listens on mesgPort for new clients. If a client connects, it
 * gets accepted and additionally a separate connection for data transfers
 * is established.
 */
void acceptClient(const unsigned short mesgPort, const unsigned short dataPort) {
  ServerSocket mesgServer(mesgPort);
  ServerSocket dataServer(dataPort);
  mesgServer.init();
  dataServer.init();
  for (;;) {
    ClientSocket mesgClient = mesgServer.acceptClient();
    assert(mesgClient.isInitialized() && "Accepted client corrupt");
    ClientSocket dataClient = dataServer.acceptClient();
    assert(dataClient.isInitialized() && "Accepted client corrupt");
    Participant participant = {mesgClient, dataClient};
    std::thread newThread(processMessages, participant);
    newThread.detach();
  }
}
