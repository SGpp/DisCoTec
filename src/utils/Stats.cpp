#include "Stats.hpp"

namespace combigrid {

bool Stats::initialized_ = false;
bool Stats::finalized_ = false;
Stats::time_point Stats::init_time_;
std::unordered_map<std::string, std::vector<Stats::Event>> Stats::event_;
std::unordered_map<std::string, std::string> Stats::attributes_;
}
// end namespace combigrid
