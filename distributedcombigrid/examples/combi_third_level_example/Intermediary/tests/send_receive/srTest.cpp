#include <iostream>
#include <string>
#include <time.h>
#include <chrono>
#include <thread>
#include "../../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

#define HOST "pcsgs02"
#define PORT 9999

#define NUM_SYSTEMS 2

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "Usage:" << std::endl << "./srTest  (s|r) rank dataSize" << std::endl;
    return -1;
  }
  if (!NetworkUtils::isInteger(argv[2])) {
    std::cout << "rank must be a number" << std::endl;
    return -1;
  }
  if (!NetworkUtils::isInteger(argv[3])) {
    std::cout << "dataSize must be a number" << std::endl;
    return -1;
  }

  int rank = std::stoi(argv[2]);
  ClientSocket client(HOST, PORT);
  if (!client.init())
    return -1;

  std::cout << "Client connected to host: " << client.getRemoteHost()
            << " at port: " << client.getRemotePort() << std::endl;

  std::string mesg = "p#" + std::to_string(rank);

  std::cout << "Client sends login message:" << mesg << std::endl;
  client.sendallPrefixed(mesg);

  std::cout << "Client asks for remote rank" << std::endl;
  mesg = "getRank";
  client.sendallPrefixed(mesg);
  client.recvallPrefixed(mesg);
  assert(NetworkUtils::isInteger(mesg));
  int remoteRank = std::stoi(mesg);
  std::cout << "Client receives remote rank: "<< mesg << std::endl;

  std::cout << "Client waits for: " << NUM_SYSTEMS << " systems" << std::endl;
  int commSize = 0;
  mesg = "getCommSize";
  std::string commSizeStr;
  while (commSize < 2) {
    std::cout << "." << std::flush;
    client.sendallPrefixed(mesg);
    client.recvallPrefixed(commSizeStr);
    assert(NetworkUtils::isInteger(commSizeStr));
    commSize = std::stoi(commSizeStr);
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  }
  std::cout << std::endl;

  size_t numElements = std::stoi(argv[3]);
  size_t dataSize = numElements * sizeof(int);
  std::vector<int> sendBuf(numElements);
  std::vector<int> recvBuf;
  for (size_t i = 0; i < numElements; ++i) {
    sendBuf[i] = i;
  }

  // send & receive multiple times.
  int rounds = 10;
  for (int i = 0; i < rounds; i++) {
    int other = (remoteRank + 1) % NUM_SYSTEMS;
    if(argv[1][0] == 's') {
      mesg = "send#" + std::to_string(other) + "#" + std::to_string(dataSize);
      std::cout << "Client wants to send: " << mesg << std::endl;
    } else {
      mesg = "recv#" + std::to_string(other) + "#" + std::to_string(dataSize);
      std::cout << "Client wants to receive: " << mesg << std::endl;
    }
    client.sendallPrefixed(mesg);
    client.recvallPrefixed(mesg);
    assert(mesg == "ok" && "Send request failed");

    mesg = "d#" + std::to_string(rank) +"#" + std::to_string(remoteRank);
    std::cout << "Client establishes new data connection: " << mesg << std::endl;
    ClientSocket dataClient(HOST, PORT);
    dataClient.init();
    dataClient.sendallPrefixed(mesg);

    if(argv[1][0] == 's') {
      std::cout << "Client starts sending" << std::endl;
      NetworkUtils::sendBinary<int>(sendBuf, dataClient);
      std::cout << "Done." << std::endl;
    } else {
      std::cout << "Client starts receiving" << std::endl;
      NetworkUtils::recvBinary<int>(recvBuf, sendBuf.size(), dataClient);

      std::cout << "Client received: [ ";
      bool passed = true;
      for (size_t i = 0; i < sendBuf.size(); ++i) {
        std::cout << recvBuf[i] << ", ";
        if (sendBuf[i] != recvBuf[i])
          passed = false;
      }
      std::cout << " ]" << std::endl;

      if (passed)
        std::cout << "Test passed!" << std::endl;
      else
        std::cout << "Test failed" << std::endl;
    }
  }
  mesg = "quit";
  client.sendallPrefixed(mesg);
  client.recvallPrefixed(mesg);
  return 0;
}
