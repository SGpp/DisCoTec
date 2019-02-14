/* ProcessManager.hpp
 *
 *  Created on: Oct 8, 2013
 *      Author: heenemo
 */

#ifndef PROCESSMANAGER_HPP_
#define PROCESSMANAGER_HPP_

#include <vector>
#include <numeric>

#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/SGrid.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"

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

  inline void
  combineThirdLevel();

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

  void setupThirdLevel();

 private:
  ProcessGroupManagerContainer& pgroups_;

  ProcessGroupManagerID& thirdLevelPGroup_;

  TaskContainer& tasks_;

  CombiParameters params_;

  ThirdLevelUtils thirdLevel_;

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

void ProcessManager::combineThirdLevel() {
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

  // groups combine local and global first
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    bool success = pgroups_[i]->combineLocalAndGlobal();
    assert(success);
  }

  waitAllFinished();

  // obtain instructions from third level manager
  thirdLevel_.signalReadyToCombine();
  std::string instruction = thirdLevel_.fetchInstruction();
  assert(instruction == "send_size");

  std::vector<std::vector<int>> commonSSPartSizes = thirdLevelPGroup_->gatherCommonSSPartSizes(thirdLevel_, params_);
  std::vector<size_t> partSizes = thirdLevelPGroup_->calcWorkersSSPartSizes(commonSSPartSizes);
  size_t wholeTransferSize = 0;
  for (auto it = partSizes.begin(); it != partSizes.end(); it++)
    wholeTransferSize += *it * sizeof(CombiDataType) + 1; // additional bit for endianness

  std::cout << "sending Sizes to third level";
  thirdLevel_.sendSize(wholeTransferSize);

  instruction = thirdLevel_.fetchInstruction();

  waitAllFinished();

  // perform third level reduce
  if (instruction == "reduce_third_level_recv_first")
  {
    bool success = thirdLevelPGroup_->reduceUniformThirdLevelRecvFirst<CombiDataType>(thirdLevel_, params_, commonSSPartSizes, partSizes);
    assert(success);
  }
  else if (instruction == "reduce_third_level_send_first")
  {
    bool success = thirdLevelPGroup_->reduceUniformThirdLevelSendFirst<CombiDataType>(thirdLevel_, params_, commonSSPartSizes, partSizes);
    assert(success);
  }

  waitAllFinished();

  // integrate subspaces
  bool success = thirdLevelPGroup_->integrateCommonSS();
  assert(success);

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
