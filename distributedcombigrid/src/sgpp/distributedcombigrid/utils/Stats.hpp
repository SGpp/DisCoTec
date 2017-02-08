#ifndef STATS_HPP_
#define STATS_HPP_

#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <assert.h>
#include <unistd.h>
#include <chrono>
#include <fstream>
#include <mpi.h>

namespace combigrid {

class Stats {
  typedef std::chrono::high_resolution_clock::time_point time_point;

  struct Event {
    const time_point start;
    time_point end;

    Event() : start(std::chrono::high_resolution_clock::now()) {}
  };

  static bool initialized_;
  static bool finalized_;
  static time_point init_time_;
  static time_point end_time_;
  static std::unordered_map<std::string, std::vector<Event>> event_;

public:
  /**
   * call at start of program
   */
  static void initialize();

  /**
   * call at end of program
   */
  static void finalize();

  /**
   * start new timer for event with given name
   */
  static void startEvent(const std::string& name);

  /**
   * stop timer for event with given name
   */
  static void stopEvent(const std::string& name);

  /**
   * write the measured times as json file to specified path using MPI-IO,
   * only call this method once finalize was called
   */
  static void write(const std::string& path);
};

}
//end namespace combigrid

#endif /* STATS_HPP_ */
