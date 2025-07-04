#include "discotec/utils/Stats.hpp"

#include <boost/date_time/posix_time/posix_time.hpp>

namespace combigrid {

std::string getTimeStamp() {
  namespace pt = boost::posix_time;
  return "[" + pt::to_iso_string(pt::second_clock::local_time()) + "] ";
}

bool Stats::initialized_ = false;
bool Stats::finalized_ = false;
Stats::time_point Stats::init_time_;
Stats::time_point Stats::partially_written_until_;
size_t Stats::numWrites_ = 0;
std::unordered_map<std::string, std::vector<Stats::Event>> Stats::event_;
std::unordered_map<std::string, std::string> Stats::attributes_;
}  // namespace combigrid
// end namespace combigrid
