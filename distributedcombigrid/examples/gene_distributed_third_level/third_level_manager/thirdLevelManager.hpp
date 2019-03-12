#include <iostream>
#include <string>
#include <vector>
#include <thread>
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
    Params params_;
    Systems systems_;
    unsigned short             dataPort_ = 9999;
    int                        timeout_  = 1000; // 1000 msec = 1sec
    ServerSocket               dataServer_;
    AmqpClient::Channel::ptr_t channel_;

    void processMessage(const std::string& message, System& system);

    void processCombination(System& system);

    void processFinished(System& system);

    void forwardData(System& sender, System& receiver) const;

  public:
    ThirdLevelManager() = delete;
    ThirdLevelManager(const Params& params);

    void runtimeLoop();
};
