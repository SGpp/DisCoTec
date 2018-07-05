#include <signal.h>
#include "Intermediary.hpp"
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

extern bool run = true;

void tunnel(const size_t port, const size_t remotePort,
    const std::string& remoteHost);

// TODO
void signalHandler(int s){
  std::cout << "Caught signal" << std::endl;
  exit(1);
}

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
  struct sigaction sigIntHandler;

  sigIntHandler.sa_handler = signalHandler;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;

  sigaction(SIGINT, &sigIntHandler, NULL);

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
    // act as intermediary
    Intermediary im(port);
    im.runForever();
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
              << " ./intermediary port [remoteHost remotePort]"
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
  for (;;)
  {
    // accept new client
    ClientSocket* newclient = server.acceptClient();
    assert(newclient->isInitialized() && "Accepting new client failed");

    // connect to remote host
    ClientSocket* remote = new ClientSocket(remoteHost, remotePort);
    remote->init();
    assert(remote->isInitialized() && "Initialization of remote failed");

    // bidirectional tunnel
    std::thread sendTunnel(&NetworkUtils::tunnel, newclient, remote, 0, 2048);
    std::thread recvTunnel(&NetworkUtils::tunnel, remote, newclient, 0, 2048);
    sendTunnel.detach();
    recvTunnel.detach();
  }
}
