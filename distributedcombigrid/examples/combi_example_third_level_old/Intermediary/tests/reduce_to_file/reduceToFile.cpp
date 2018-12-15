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
  if (argc < 3) {
    std::cout << "Usage:" << std::endl << "./reduceToFile rank dataSize" << std::endl;
    return -1;
  }
  if (!NetworkUtils::isInteger(argv[1])) {
    std::cout << "rank must be a number" << std::endl;
    return -1;
  }
  if (!NetworkUtils::isInteger(argv[2])) {
    std::cout << "dataSize must be a number" << std::endl;
    return -1;
  }

  int rank = std::stoi(argv[1]);
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
  while (commSize < NUM_SYSTEMS) {
    std::cout << "." << std::flush;
    client.sendallPrefixed(mesg);
    client.recvallPrefixed(commSizeStr);
    assert(NetworkUtils::isInteger(commSizeStr));
    commSize = std::stoi(commSizeStr);
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  }
  std::cout << std::endl;

  size_t numElements = std::stoi(argv[2]);
  //size_t dataSize = numElements * sizeof(double);
  std::vector<double> sendBuf(numElements);
  std::cout << "Sending values: [ ";
  for (size_t i = 0; i < numElements; ++i) {
    double value = (double) i;
    sendBuf[i] = value;
    std::cout << value << " ";

  }
  std::cout << "]" << std::endl;

  mesg = "reduceToFileUniform#" + std::to_string(NUM_SYSTEMS) + "#"
    + std::to_string(numElements) + "#" +  std::to_string(sizeof(double)) +
    "#real";
  client.sendallPrefixed(mesg);
  client.recvallPrefixed(mesg);

  mesg = "d#" + std::to_string(rank) +"#" + std::to_string(remoteRank);
  ClientSocket dataClient(HOST, PORT);
  dataClient.init();
  dataClient.sendallPrefixed(mesg);

  NetworkUtils::sendBinary(sendBuf, dataClient);

  mesg = "quit";
  client.sendallPrefixed(mesg);
  client.recvallPrefixed(mesg);
  return 0;
}
