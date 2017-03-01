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
#include "sgpp/distributedcombigrid/fault_tolerance/LPOptimizationInterpolation.hpp"
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

  inline void
  recomputeOptimumCoefficients( std::string prob_name,
                                std::vector<int>& faultsID,
                                std::vector<int>& redistributefaultsID,
                                std::vector<int>& recomputeFaultsID);

  inline Task* getTask( int taskID );

  void
  updateCombiParameters();


  /* Computes group faults in current combi scheme step */
  void
  getGroupFaultIDs( std::vector<int>& faultsID );

  inline CombiParameters& getCombiParameters();

  void parallelEval( const LevelVector& leval,
                                     std::string& filename,
                                     size_t groupID );
  
void redistribute( std::vector<int>& taskID );

  void recompute( std::vector<int>& taskID );

  void recover();

  void recoverCommunicators();

  /* After faults have been fixed, we need to return the combischeme
   * to the original combination technique*/
  void restoreCombischeme();


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

/*
 * Recompute coefficients for the combination technique based on given grid faults using
 * an optimization scheme
 */
inline void
ProcessManager::recomputeOptimumCoefficients(std::string prob_name,
                            std::vector<int>& faultsID,
                            std::vector<int>& redistributeFaultsID,
                            std::vector<int>& recomputeFaultsID) {

  CombigridDict given_dict = params_.getCombiDict();

  std::map<int, LevelVector> IDsToLevels = params_.getLevelsDict();
  LevelVectorList faultLevelVectors;
  for (auto id : faultsID)
    faultLevelVectors.push_back(IDsToLevels[id]);


  LevelVectorList lvlminmax;
  lvlminmax.push_back(params_.getLMin());
  lvlminmax.push_back(params_.getLMax());
  LP_OPT_INTERP opt_interp(lvlminmax,static_cast<int>(params_.getDim()),GLP_MAX,given_dict,faultLevelVectors);
  if ( opt_interp.getNumFaults() != 0 ) {

    opt_interp.init_opti_prob(prob_name);
    opt_interp.set_constr_matrix();
    opt_interp.solve_opti_problem();

    LevelVectorList recomputeLevelVectors;
    CombigridDict new_dict = opt_interp.get_results(recomputeLevelVectors);

    LevelVectorList newLevels;
    std::vector<real> newCoeffs;
    std::vector<int> newTaskIDs;

    int numLevels = int(params_.getNumLevels());
    for( int i = 0; i < numLevels; ++i ){
      LevelVector lvl = params_.getLevel(i);

      newLevels.push_back(lvl);
      newCoeffs.push_back(new_dict[lvl]);
      newTaskIDs.push_back(i);
    }

    params_.setLevelsCoeffs(newTaskIDs, newLevels, newCoeffs);

    std::map<LevelVector, int> LevelsToIDs = params_.getLevelsToIDs();
    for(auto l : recomputeLevelVectors){
      recomputeFaultsID.push_back(LevelsToIDs[l]);
    }

    std::sort(faultsID.begin(), faultsID.end());
    std::sort(recomputeFaultsID.begin(), recomputeFaultsID.end());

    std::set_difference(faultsID.begin(), faultsID.end(),
                        recomputeFaultsID.begin(), recomputeFaultsID.end(),
                        std::inserter(redistributeFaultsID, redistributeFaultsID.begin()));
  }
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
