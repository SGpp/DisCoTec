#include <iostream>
#include <string>
#include <vector>
#include <boost/property_tree/ini_parser.hpp>
#include "SimpleAmqpClient/SimpleAmqpClient.h"
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "Params.hpp"
#include "sgpp/distributedcombigrid/third_level/MessageUtils.hpp"
#include "System.hpp"

using  Systems = std::vector<System>;

class ThirdLevelManager
{
  private:
    Params _params;
    Systems _systems;
    u_int                      _dataPort = 9999;
    u_int                      _timeout  = 1000; // 1000 msec = 1sec
    ServerSocket               _dataServer;
    AmqpClient::Channel::ptr_t _channel;

    void processMessage(const std::string& message, System& system);

    void processCombination(System& system);

    void processFinished(System& system);

  public:
    ThirdLevelManager() = delete;
    ThirdLevelManager(const Params& params);

    void runtimeLoop();
};
