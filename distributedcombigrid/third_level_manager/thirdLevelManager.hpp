#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <boost/property_tree/ini_parser.hpp>
#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "Params.hpp"
#include "System.hpp"

using  Systems = std::vector<System>;

class ThirdLevelManager
{
  private:
    Params         params_;
    Systems        systems_;
    unsigned short port_    = 9999;
    int            timeout_ = 1;
    ServerSocket   server_;

    void processMessage(const std::string& message, size_t sysIndex);

    void processCombination(size_t initiatorIndex);

    void processUnifySubspaceSizes(size_t initiatorIndex);

    void processFinished(size_t sysIndex);

    void forwardData(const System& sender, const System& receiver) const;

  public:
    ThirdLevelManager() = delete;
    ThirdLevelManager(const Params& params);

    void init();

    void runtimeLoop();
};
