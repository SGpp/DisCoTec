/* ProcessManager.hpp
 *
 *  Created on: Oct 8, 2013
 *      Author: heenemo
 */

#ifndef PROCESSMANAGER_HPP_
#define PROCESSMANAGER_HPP_

#include <vector>
#include <numeric>

#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/LPOptimizationInterpolation.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"

namespace combigrid {

class ProcessManager {
 public:
  ProcessManager(ProcessGroupManagerContainer& pgroups, TaskContainer& instances,
                 CombiParameters& params, std::unique_ptr<LoadModel> loadModel)
    : pgroups_(pgroups),
      thirdLevelPGroup_(pgroups[0]), // TODO: changing PG requires adjustments in integrateDSGUniform function from worker
      tasks_(instances),
      params_(params),
      thirdLevel_(params.getThirdLevelHost(), params.getThirdLevelPort(), params.getThirdLevelSystemName())
  {
      loadModel_ = std::move(loadModel);
      if (params.getThirdLevelHost() != "")
        setupThirdLevel();
      // the combiparameters are sent to all process groups before the
      // computations start
      //updateCombiParameters();
  }

  inline void removeGroups(std::vector<int> removeIndices);

  // todo: use general class AppInstance here
  // todo: add remove function
  inline void addTask(Task* t);

  bool runfirst();

  void exit();

  virtual ~ProcessManager();

  template <typename FG_ELEMENT>
  inline FG_ELEMENT eval(const std::vector<real>& coords);

  bool runnext();

  inline void combine();

  inline void
  combineThirdLevel();

  inline void combineLocalAndGlobal();

  template <typename FG_ELEMENT>
  inline void combineFG(FullGrid<FG_ELEMENT>& fg);

  template <typename FG_ELEMENT>
  inline void gridEval(FullGrid<FG_ELEMENT>& fg);

  /* Generates no_faults random faults from the combischeme */
  inline void createRandomFaults(std::vector<int>& faultIds, int no_faults);

  inline void recomputeOptimumCoefficients(std::string prob_name, std::vector<int>& faultsID,
                                           std::vector<int>& redistributefaultsID,
                                           std::vector<int>& recomputeFaultsID);

  inline Task* getTask(int taskID);

  void updateCombiParameters();

  void getDSGFromProcessGroup();

  /* Computes group faults in current combi scheme step */
  void getGroupFaultIDs(std::vector<int>& faultsID,
                        std::vector<ProcessGroupManagerID>& groupFaults);

  inline CombiParameters& getCombiParameters();

  void parallelEval(const LevelVector& leval, std::string& filename, size_t groupID);

  void redistribute(std::vector<int>& taskID);

  void reInitializeGroup(std::vector<ProcessGroupManagerID>& taskID,
                         std::vector<int>& tasksToIgnore);

  void recompute(std::vector<int>& taskID, bool failedRecovery,
                 std::vector<ProcessGroupManagerID>& recoveredGroups);

  void recover(int i, int nsteps);

  bool recoverCommunicators(std::vector<ProcessGroupManagerID> failedGroups);
  /* After faults have been fixed, we need to return the combischeme
   * to the original combination technique*/
  void restoreCombischeme();

  void setupThirdLevel();

 private:
  ProcessGroupManagerContainer& pgroups_;

  ProcessGroupManagerID& thirdLevelPGroup_;

  TaskContainer& tasks_;

  CombiParameters params_;

  ThirdLevelUtils thirdLevel_;

  std::unique_ptr<LoadModel> loadModel_;

  // periodically checks status of all process groups. returns until at least
  // one group is in WAIT state
  inline ProcessGroupManagerID wait();
  inline ProcessGroupManagerID waitAvoid(std::vector<ProcessGroupManagerID>& avoidGroups);
  bool waitAllFinished();

  void receiveDurationsOfTasksFromGroupMasters(size_t numDurationsToReceive);

  void sortTasks();

  void sendDSGUniformToRemote(ProcessGroupManagerID& pg);

  void recvDSGUniformFromRemote(ProcessGroupManagerID& pg);

  void recvAndAddDSGUniformFromRemote(ProcessGroupManagerID& pg);
};

inline void ProcessManager::addTask(Task* t) { tasks_.push_back(t); }

inline ProcessGroupManagerID ProcessManager::wait() {
  while (true) {
    for (size_t i = 0; i < pgroups_.size(); ++i) {
      StatusType status = pgroups_[i]->getStatus();
      //      std::cout << status << " is the status of : "<< i << "\n";
      assert(status >= 0);  // check for invalid values
      assert(status <= 2);
      if (status == PROCESS_GROUP_WAIT) {
        return pgroups_[i];
      }
    }
  }
}

inline ProcessGroupManagerID ProcessManager::waitAvoid(
    std::vector<ProcessGroupManagerID>& avoidGroups) {
  while (true) {
    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (std::find(avoidGroups.begin(), avoidGroups.end(), pgroups_[i]) ==
          avoidGroups.end()) {  // ignore tasks that are recomputed
        StatusType status = pgroups_[i]->getStatus();
        //      std::cout << status << " is the status of : "<< i << "\n";
        assert(status >= 0);  // check for invalid values
        assert(status <= 2);
        if (status == PROCESS_GROUP_WAIT) {
          return pgroups_[i];
        }
      }
    }
  }
}

template <typename FG_ELEMENT>
inline FG_ELEMENT ProcessManager::eval(const std::vector<real>& coords) {
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }

  FG_ELEMENT res(0);

  // call eval function of each process group
  for (size_t i = 0; i < pgroups_.size(); ++i) res += pgroups_[i]->eval(coords);

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
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }

  // send signal to each group
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    bool success = pgroups_[i]->combine();
    assert(success);
  }

  waitAllFinished();
}


void ProcessManager::combineLocalAndGlobal() {
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
}

/** Combination with third level parallelism e.g. between two HPC systems
 *
 * This process manager induces a local and global combination first.
 * Then he signals ready to the third level manager who decides about the role
 * of sender and receiver.
 *
 * In the role of the sender, the process manager
 * transfers the DSGU data from workers of the third level pg to the third level
 * manager, who then finally sends it to the remote system. After sending he
 * receives the remotely reduced data from the third level manager and sends it
 * back to the third level pg.
 * In the role of the receiver, the process manager receives the DSGU data from
 * the third level manager and sends it to the workers who combine the remote
 * solution with their local solution. Afterward, the final solution is sent
 * back to the remote system.
 */
void ProcessManager::combineThirdLevel() {
  combineLocalAndGlobal();

  // obtain instructions from third level manager
  std::cout << "Signaling ready to combine..." << std::endl;
  thirdLevel_.signalReadyToCombine();
  std::string instruction = thirdLevel_.fetchInstruction();

  // perform third level reduce
  if (instruction == "reduce_third_level_recv_first") {
    recvAndAddDSGUniformFromRemote(thirdLevelPGroup_);
    waitAllFinished();
    sendDSGUniformToRemote(thirdLevelPGroup_);
  } else if (instruction == "reduce_third_level_send_first") {
    sendDSGUniformToRemote(thirdLevelPGroup_);
    waitAllFinished();
    recvDSGUniformFromRemote(thirdLevelPGroup_);
  }
  waitAllFinished();

  // TODO integrate into receive method
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    bool success = pgroups_[i]->integrateCombinedDSGUniform();
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
template <typename FG_ELEMENT>
void ProcessManager::combineFG(FullGrid<FG_ELEMENT>& fg) {
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }

  // send signal to each group
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    bool success = pgroups_[i]->combineFG(fg);
    assert(success);
  }

  CombiCom::FGAllreduce<FG_ELEMENT>(fg, theMPISystem()->getGlobalComm());
}

/* Evaluate the combination solution with the resolution of the given full grid.
 * In constrast to the combineFG function, the solution will only be available
 * on the manager. No recombination is performed, i.e. the local component grids
 * won't be updated.
 */
template <typename FG_ELEMENT>
void ProcessManager::gridEval(FullGrid<FG_ELEMENT>& fg) {
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }

  // send signal to each group
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    bool success = pgroups_[i]->gridEval(fg);
    assert(success);
  }

  CombiCom::FGReduce<FG_ELEMENT>(fg, theMPISystem()->getManagerRank(),
                                 theMPISystem()->getGlobalComm());
}

CombiParameters& ProcessManager::getCombiParameters() { return params_; }

/*
 * Create a certain given number of random faults, considering that the faulty processes
 * simply cannot give the evaluation results, but they are still available in the MPI
 * communication scheme (the nodes are not dead)
 */
inline void ProcessManager::createRandomFaults(std::vector<int>& faultIds, int no_faults) {
  int fault_id;

  // create random faults
  int j = 0;
  while (j < no_faults) {
    fault_id = generate_random_fault(static_cast<int>(params_.getNumLevels()));
    if (j == 0 || std::find(faultIds.begin(), faultIds.end(), fault_id) == faultIds.end()) {
      faultIds.push_back(fault_id);
      j++;
    }
  }
}

/*
 * Recompute coefficients for the combination technique based on given grid faults using
 * an optimization scheme
 */
inline void ProcessManager::recomputeOptimumCoefficients(std::string prob_name,
                                                         std::vector<int>& faultsID,
                                                         std::vector<int>& redistributeFaultsID,
                                                         std::vector<int>& recomputeFaultsID) {
  CombigridDict given_dict = params_.getCombiDict();

  std::map<int, LevelVector> IDsToLevels = params_.getLevelsDict();
  LevelVectorList faultLevelVectors;
  for (auto id : faultsID) faultLevelVectors.push_back(IDsToLevels[id]);

  LevelVectorList lvlminmax;
  lvlminmax.push_back(params_.getLMin());
  lvlminmax.push_back(params_.getLMax());
  LP_OPT_INTERP opt_interp(lvlminmax, static_cast<int>(params_.getDim()), GLP_MAX, given_dict,
                           faultLevelVectors);
  if (opt_interp.getNumFaults() != 0) {
    opt_interp.init_opti_prob(prob_name);
    opt_interp.set_constr_matrix();
    opt_interp.solve_opti_problem();

    LevelVectorList recomputeLevelVectors;
    CombigridDict new_dict = opt_interp.get_results(recomputeLevelVectors);

    LevelVectorList newLevels;
    std::vector<real> newCoeffs;
    std::vector<int> newTaskIDs;

    int numLevels = int(params_.getNumLevels());
    for (int i = 0; i < numLevels; ++i) {
      LevelVector lvl = params_.getLevel(i);

      newLevels.push_back(lvl);
      newCoeffs.push_back(new_dict[lvl]);
      newTaskIDs.push_back(i);
    }
    // check if sum of coefficients is 1
    double sum = 0.0;
    std::cout << "new coefficients: ";
    for (size_t i = 0; i < newCoeffs.size(); i++) {
      sum += newCoeffs[i];
      std::cout << newCoeffs[i] << " ";
    }
    std::cout << "\n";
    int roundedSum = round(sum);
    std::cout << "Coefficient sum: " << roundedSum << "\n";

    assert(roundedSum == 1);
    params_.setLevelsCoeffs(newTaskIDs, newLevels, newCoeffs);

    std::map<LevelVector, int> LevelsToIDs = params_.getLevelsToIDs();
    for (auto l : recomputeLevelVectors) {
      recomputeFaultsID.push_back(LevelsToIDs[l]);
    }

    std::sort(faultsID.begin(), faultsID.end());
    std::sort(recomputeFaultsID.begin(), recomputeFaultsID.end());

    std::set_difference(faultsID.begin(), faultsID.end(), recomputeFaultsID.begin(),
                        recomputeFaultsID.end(),
                        std::inserter(redistributeFaultsID, redistributeFaultsID.begin()));
  }
}

inline Task* ProcessManager::getTask(int taskID) {
  for (Task* tmp : tasks_) {
    if (tmp->getID() == taskID) {
      return tmp;
    }
  }
  return nullptr;
}

} /* namespace combigrid */
#endif /* PROCESSMANAGER_HPP_ */
