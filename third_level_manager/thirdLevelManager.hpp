#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <boost/property_tree/json_parser.hpp>
#include "discotec/third_level/NetworkUtils.hpp"
#include "Params.hpp"
#include "System.hpp"
#include "Stats.hpp"

namespace combigrid {

using Systems = std::vector<System>;

class ThirdLevelManager
{
  private:
    Params         params_;
    Systems        systems_;
    unsigned short port_    = 9999;
    int            timeout_ = 1;
    ServerSocket   server_;
    Stats          stats_;

    void processMessage(const std::string& message, size_t sysIndex);

    void processCombination(size_t initiatorIndex);

    void processCombinationFile(size_t initiatorIndex);

    void processUnifySubspaceSizes(size_t initiatorIndex);

    void processAnyData(size_t initiatorIndex);

    void processFinished(size_t sysIndex);

    size_t forwardData(const System& sender, const System& receiver) const;

  public:
    ThirdLevelManager() = delete;
    ThirdLevelManager(const Params& params);

    void init();

    void runtimeLoop();

    void writeStatistics(const std::string& filename = "");
};

}
