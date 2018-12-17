#include <iostream>
#include <string>
#include <vector>
#include <boost/property_tree/ini_parser.hpp>
#include "SimpleAmqpClient/SimpleAmqpClient.h"
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "Params.hpp"
#include "System.hpp"

using  Systems = std::vector<System>;

const u_int timeout = 1000; // 1 sec

class ThirdLevelManager
{
  private:
    const Params& _params;
    Systems _systems;
    AmqpClient::Channel::ptr_t _channel;

    void processMessage(const std::string& message, System& system);

    void processCombination(System& system);

    void processFinished(System& system);

  public:
    ThirdLevelManager() = delete;
    ThirdLevelManager(const Params& params);

    void runtimeLoop();
};
