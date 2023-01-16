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
#include "mpi/MPISystem.hpp"

/* comment this line to switch of timing */
//#define TIMING

namespace combigrid {

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
  static void write(const std::string& path, CommunicatorType comm = theMPISystem()->getWorldComm());

 private:
  static bool initialized_;
  static bool finalized_;
  static time_point init_time_;
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

inline void writeSingleFile(const std::stringstream& buffer, const std::string& path,
                            CommunicatorType comm) {
  int rank = getCommRank(comm);

  // get offset in file
  MPI_Offset len = buffer.str().size();
  MPI_Offset pos = 0;
  MPI_Scan(&len, &pos, 1, MPI_LONG, MPI_SUM, comm);
  pos -= len;

  // see: https://wickie.hlrs.de/platforms/index.php/MPI-IO
  MPI_Info info = MPI_INFO_NULL;

  // open file
  MPI_File fh;
  int err = MPI_File_open(comm, path.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
                          info, &fh);
  if (err != MPI_SUCCESS) {
    // file already existed, delete it and create new file
    if (rank == 0) {
      MPI_File_delete(path.c_str(), MPI_INFO_NULL);
    }
    MPI_File_open(comm, path.c_str(), MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, info, &fh);
  }

  // write to single file with MPI-IO
  MPI_File_write_at_all(fh, pos, buffer.str().c_str(), (int)len, MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_File_close(&fh);
}

inline void Stats::write(const std::string& path, CommunicatorType comm) {
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
  writeSingleFile(buffer, path, comm);
}

#else
inline void Stats::initialize() {}
inline void Stats::finalize() {}
inline void Stats::startEvent(const std::string& name) {
  event_[name].emplace_back();
}
inline void Stats::setAttribute(const std::string& name, const std::string& value) {}
inline void Stats::write(const std::string& path, CommunicatorType comm) {}
#endif

inline const Stats::Event Stats::stopEvent(const std::string& name) {
  // check if event is not stopped already
  assert(event_[name].back().end.time_since_epoch().count() == 0);

  event_[name].back().end = std::chrono::high_resolution_clock::now();
  Event e(event_[name].back());
#ifndef TIMING
  //prevent map from growing if TIMING is off
  event_.erase(name);
#endif // def TIMING
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
