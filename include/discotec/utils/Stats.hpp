#ifndef STATS_HPP_
#define STATS_HPP_

#include <assert.h>
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <unistd.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "../io/FileInputOutput.hpp"
#include "../io/MPIInputOutput.hpp"
#include "../mpi/MPISystem.hpp"

namespace combigrid {

std::string getTimeStamp();
class Stats {
  typedef std::chrono::high_resolution_clock::time_point time_point;

 public:
  struct Event {
    const time_point start;
    time_point end;

    Event() : start(std::chrono::high_resolution_clock::now()) {}
  };

  static long unsigned int getEventDuration(const Stats::Event e){
    std::chrono::milliseconds x = std::chrono::duration_cast<std::chrono::milliseconds>(e.end - e.start);
    return x.count();
  }
  static long unsigned int getEventDurationInUsec(const Stats::Event e){
    std::chrono::microseconds x = std::chrono::duration_cast<std::chrono::microseconds>(e.end - e.start);
    return x.count();
  }

  /**
   * use at the start of the program
   */
  static void initialize();

  /**
   * use at the end of the program
   */
  static void finalize();

  static bool isInitialized() { return initialized_; }

  /**
   * start a new timer for event with given name
   */
  static void startEvent(const std::string& name);

  /**
   * stop a timer for event with given name
   */
  static const Event stopEvent(const std::string& name);

  /**
   * get the event with given name
   */
  static const Event& getEvent(const std::string& name);

  /**
   * get the stopped time for event with given name
   */
  static long unsigned int getDuration(const std::string& name);

  /**
   * set an attribute which can later be used for plotting
   */
  static void setAttribute(const std::string& name, const std::string& value);

  /**
   * write the measured times in json format to specified path,
   * only call this after finalize
   */
  static void write(const std::string& path,
                    CommunicatorType comm = theMPISystem()->getWorldComm());

  /**
   * write the measured times until now in json format to specified path
   */
  static void writePartial(const std::string& pathPrefix,
                           CommunicatorType comm = theMPISystem()->getWorldComm());

 private:
  static bool initialized_;
  static bool finalized_;
  static time_point init_time_;
  static time_point partially_written_until_;
  static size_t numWrites_;
  static std::unordered_map<std::string, std::vector<Event>> event_;
  static std::unordered_map<std::string, std::string> attributes_;
};

#ifdef TIMING
inline void Stats::initialize() {
  assert(!initialized_);

  // clear data (in case of multiple calls to initialize, e.g. in tests)
  event_.clear();
  attributes_.clear();

  initialized_ = true;
  finalized_ = false;
  init_time_ = std::chrono::high_resolution_clock::now();
  partially_written_until_ = init_time_;
  numWrites_ = 0;
}

inline void Stats::finalize() {
  assert(initialized_);
  finalized_ = true;
  initialized_ = false;
}

inline void Stats::startEvent(const std::string& name) {
  ASSERT(initialized_, name);

  event_[name].emplace_back();
}

inline void Stats::setAttribute(const std::string& name, const std::string& value) {
  ASSERT(initialized_, name << " " << value);

  attributes_[name] = value;
}

inline void Stats::write(const std::string& path, CommunicatorType comm) {
  std::string myJSONpart = "";
  {
    MPI_Barrier(comm);
    assert(finalized_);

    using namespace std::chrono;
    // MPI_Comm worldComm = theMPISystem()->getWorldComm();
    int rank = getCommRank(comm);
    int size = getCommSize(comm);

    std::stringstream buffer;

    if (rank == 0) {
      buffer << "{" << std::endl;
    }

    buffer << "\"rank" << rank << "\":{" << std::endl;
    buffer << "\"attributes\":{" << std::endl;
    std::size_t attributes_count = 0;
    for (auto&& it = attributes_.begin(); it != attributes_.end(); ++it) {
      buffer << "\"" << it->first << "\":\"" << it->second << "\"";
      if (++attributes_count != attributes_.size()) {
        buffer << "," << std::endl;
      } else {
        buffer << std::endl;
      }
    }
    buffer << "}," << std::endl;

    buffer << "\"events\":{" << std::endl;
    std::size_t event_count = 0;
    for (auto&& t = event_.begin(); t != event_.end(); ++t) {
      buffer << "\"" << t->first << "\""
             << ":[" << std::endl;

      std::size_t element_count = 0;
      for (auto&& e = t->second.begin(); e != t->second.end(); ++e) {
        buffer << "[" << duration_cast<microseconds>(e->start - init_time_).count() << ","
               << duration_cast<microseconds>(e->end - init_time_).count() << "]";
        if (++element_count != t->second.size()) {
          buffer << "," << std::endl;
        } else {
          buffer << std::endl;
        }
      }
      buffer << "]";
      if (++event_count != event_.size()) {
        buffer << "," << std::endl;
      } else {
        buffer << std::endl;
      }
    }
    buffer << "}" << std::endl << "}";

    if (rank != size - 1) {
      buffer << "," << std::endl;
    } else {
      buffer << std::endl << "}" << std::endl;
    }
    myJSONpart = buffer.str();

  }
  combigrid::writeConcatenatedFileRootOnly(myJSONpart.data(), myJSONpart.size(), path, comm, true);
}

inline void Stats::writePartial(const std::string& pathSuffix, CommunicatorType comm) {
  std::string myJSONpart = "";
  std::string path = std::to_string(numWrites_++) + "_" + pathSuffix;
  {
    MPI_Barrier(comm);

    using namespace std::chrono;
    int rank = getCommRank(comm);
    int size = getCommSize(comm);

    std::stringstream buffer;

    if (rank == 0) {
      buffer << "{" << std::endl;
    }

    buffer << "\"rank" << rank << "\":{" << std::endl;
    buffer << "\"attributes\":{" << std::endl;
    std::size_t attributes_count = 0;
    // write all attributes
    for (auto&& it = attributes_.begin(); it != attributes_.end(); ++it) {
      buffer << "\"" << it->first << "\":\"" << it->second << "\"";
      if (++attributes_count != attributes_.size()) {
        buffer << "," << std::endl;
      } else {
        buffer << std::endl;
      }
    }
    buffer << "}," << std::endl;
    buffer << "\"events\":{" << std::endl;
    std::stringstream allEventsBuffer;
    bool anyEventFirstTime = true;
    // but filter events for finished and not yet written
    for (auto&& t = event_.begin(); t != event_.end(); ++t) {
      std::stringstream eventBuffer;
      if (!anyEventFirstTime) {
        eventBuffer << "," << std::endl;
      }
      eventBuffer << "\"" << t->first << "\""
                  << ":[" << std::endl;

      bool thisEventFirstTime = true;
      for (auto&& e = t->second.begin(); e != t->second.end(); ++e) {
        if (e->end != e->start && e->end > partially_written_until_) {
          if (thisEventFirstTime) {
            thisEventFirstTime = false;
            anyEventFirstTime = false;
          } else {
            eventBuffer << "," << std::endl;
          }
          eventBuffer << "[" << duration_cast<microseconds>(e->start - init_time_).count() << ","
                      << duration_cast<microseconds>(e->end - init_time_).count() << "]";
        }
      }
      eventBuffer << "]";
      if (!thisEventFirstTime) {
        allEventsBuffer << eventBuffer.str();
      }
    }
    allEventsBuffer << std::endl;
    buffer << allEventsBuffer.str();
    buffer << "}" << std::endl << "}";

    if (rank != size - 1) {
      buffer << "," << std::endl;
    } else {
      buffer << std::endl << "}" << std::endl;
    }

    myJSONpart = buffer.str();
  }
  combigrid::writeConcatenatedFileRootOnly(myJSONpart.data(), myJSONpart.size(), path, comm, true);

  partially_written_until_ = std::chrono::high_resolution_clock::now();
}

#else
inline void Stats::initialize() {}
inline void Stats::finalize() {}
inline void Stats::startEvent(const std::string& name) { event_[name].emplace_back(); }
inline void Stats::setAttribute(const std::string& name, const std::string& value) {}
inline void Stats::write(const std::string& path, CommunicatorType comm) {}
inline void Stats::writePartial(const std::string& pathSuffix, CommunicatorType comm) {}
#endif

inline const Stats::Event Stats::stopEvent(const std::string& name) {
  // check if event is not stopped already
  assert(event_[name].back().end.time_since_epoch().count() == 0);

  event_[name].back().end = std::chrono::high_resolution_clock::now();
  Event e(event_[name].back());
#ifndef TIMING
  // prevent map from growing if TIMING is off
  event_.erase(name);
#endif  // def TIMING
  return e;
}

inline const Stats::Event& Stats::getEvent(const std::string& name) { return event_[name].back(); }

inline long unsigned int Stats::getDuration(const std::string& name) {
  if (event_.find(name) != event_.end()) {
    return getEventDuration(getEvent(name));
  } else {
    return std::numeric_limits<long unsigned int>::max();
  }
}
}  // namespace combigrid
// end namespace combigrid

#endif /* STATS_HPP_ */
