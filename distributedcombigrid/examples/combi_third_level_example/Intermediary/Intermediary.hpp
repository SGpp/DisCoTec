#ifndef INTERMEDIARYHPP_
#define INTERMEDIARYHPP_

#include <iostream>
#include <mutex>
#include <thread>
#include <map>
#include <vector>
#include "Communicator.hpp"
#include <complex>
#include <cstdio>
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

class Communicator;

class Intermediary {

  public:
    Intermediary(int port);
    ~Intermediary();
    void runForever();

  private:
    std::map<int, Communicator*> comms_;
    int port_;
    void acceptClient();
    void handleClient(ClientSocket* client);
};

#endif
