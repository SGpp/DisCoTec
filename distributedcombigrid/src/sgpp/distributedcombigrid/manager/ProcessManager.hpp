/*
 * ProcessManager.hpp
 *
 *  Created on: Oct 8, 2013
 *      Author: heenemo
 */

#ifndef PROCESSMANAGER_HPP_
#define PROCESSMANAGER_HPP_

#include <vector>

#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/SGrid.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"

namespace combigrid {

class ProcessManager {
 public:
  ProcessManager(ProcessGroupManagerContainer& pgroups,
                 TaskContainer& instances,
                 CombiParameters& params );

  inline void
  removeGroups(std::vector<int> removeIndices);

  // todo: use general class AppInstance here
  // todo: add remove function
  inline void
  addTask(Task* t);

  bool
  runfirst();

  void
  exit();

  virtual
  ~ProcessManager();

  template<typename FG_ELEMENT>
  inline FG_ELEMENT
  eval(const std::vector<real>& coords);

  bool
  runnext();

  inline void
  combine();

  template<typename FG_ELEMENT>
  inline void
  combineFG(FullGrid<FG_ELEMENT>& fg);

  template<typename FG_ELEMENT>
  inline void
  gridEval(FullGrid<FG_ELEMENT>& fg);

  inline Task* getTask( int taskID );

  void
  updateCombiParameters();

  inline CombiParameters& getCombiParameters();

  void parallelEval( const LevelVector& leval,
                                     std::string& filename,
                                     size_t groupID );

 private:
  ProcessGroupManagerContainer& pgroups_;

  TaskContainer& tasks_;

  CombiParameters params_;

  // periodically checks status of all process groups. returns until at least
  // one group is in WAIT state
  inline ProcessGroupManagerID wait();

  bool waitAllFinished();
};


inline void ProcessManager::addTask(Task* t) {
  tasks_.push_back(t);
}

inline ProcessGroupManagerID ProcessManager::wait() {
  while (true) {
    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT)
        return pgroups_[i];
    }
  }
}

template<typename FG_ELEMENT>
inline FG_ELEMENT ProcessManager::eval(const std::vector<real>& coords) {
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }
  }

  FG_ELEMENT res(0);

  // call eval function of each process group
  for (size_t i = 0; i < pgroups_.size(); ++i)
    res += pgroups_[i]->eval(coords);

  return res;
}

/* This function performs the so-called recombination. First, the combination
 * solution will be evaluated in the given sparse grid space.
 * Also, the local component grids will be updated with the combination
 * solution. The combination solution will also be available on the manager
 * process.
 */
void ProcessManager::combine() {
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }
  }

  // send signal to each group
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    bool success = pgroups_[i]->combine();
    assert(success);
  }

  waitAllFinished();
}

/* This function performs the so-called recombination. First, the combination
 * solution will be evaluated with the resolution of the given full grid.
 * Afterwards, the local component grids will be updated with the combination
 * solution. The combination solution will also be available on the manager
 * process.
 */
template<typename FG_ELEMENT>
void ProcessManager::combineFG(FullGrid<FG_ELEMENT>& fg) {
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }
  }

  // send signal to each group
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    bool success = pgroups_[i]->combineFG(fg);
    assert(success);
  }

  CombiCom::FGAllreduce<FG_ELEMENT>( fg, theMPISystem()->getGlobalComm() );
}

/* Evaluate the combination solution with the resolution of the given full grid.
 * In constrast to the combineFG function, the solution will only be available
 * on the manager. No recombination is performed, i.e. the local component grids
 * won't be updated.
 */
template<typename FG_ELEMENT>
void ProcessManager::gridEval(FullGrid<FG_ELEMENT>& fg) {
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }
  }

  // send signal to each group
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    bool success = pgroups_[i]->gridEval(fg);
    assert(success);
  }

  CombiCom::FGReduce<FG_ELEMENT>( fg,
                                  theMPISystem()->getManagerRank(),
                                  theMPISystem()->getGlobalComm() );
}


CombiParameters& ProcessManager::getCombiParameters() {
  return params_;
}


inline Task*
ProcessManager::getTask( int taskID ){

  for ( Task* tmp : tasks_ ) {
    if ( tmp->getID() == taskID ) {
      return tmp;
    }
  }
  return nullptr;
}

} /* namespace combigrid */
#endif /* PROCESSMANAGER_HPP_ */
