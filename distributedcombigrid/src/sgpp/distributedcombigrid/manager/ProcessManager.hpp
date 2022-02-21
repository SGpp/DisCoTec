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
#include "sgpp/distributedcombigrid/rescheduler/TaskRescheduler.hpp"
#include "sgpp/distributedcombigrid/rescheduler/StaticTaskRescheduler.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"

namespace combigrid {

class ProcessManager {
 public:
   /**
    * Constructor for a process manager.
    *
    * @param pgroups The process groups.
    * @param instances The tasks.
    * @param params The parameters for the combination technique.
    * @param loadModel The load model to use for scheduling. If it is a 
    *                  learning load model duration information is added for 
    *                  every task after every run.
    * @param rescheduler The rescheduler to use for dynamic task rescheduling.
    *                    By default the static task rescheduler is used and 
    *                    therefore no rescheduling perfomed.
    */
  ProcessManager(ProcessGroupManagerContainer& pgroups, TaskContainer& instances,
                 CombiParameters& params, std::unique_ptr<LoadModel> loadModel, 
                 std::unique_ptr<TaskRescheduler> rescheduler = std::unique_ptr<TaskRescheduler>(new StaticTaskRescheduler{}))
    : pgroups_{pgroups},
      tasks_{instances},
      params_{params},
      loadModel_{std::move(loadModel)},
      rescheduler_{std::move(rescheduler)},
      thirdLevel_{params.getThirdLevelHost(), params.getThirdLevelPort()},
      thirdLevelPGroup_{pgroups_[params.getThirdLevelPG()]}
  {
     // only setup third level if explicitly desired
     if (params.getThirdLevelHost() != "")
       setupThirdLevel();
  }

  inline void removeGroups(std::vector<int> removeIndices);

  // todo: use general class AppInstance here
  // todo: add remove function
  inline void addTask(Task* t);

  bool runfirst();

  void initDsgus();

  void exit();

  virtual ~ProcessManager();

  template <typename FG_ELEMENT>
  inline FG_ELEMENT eval(const std::vector<real>& coords);

  bool runnext();

  inline void combine();

  inline void combineThirdLevel();

  inline void unifySubspaceSizesThirdLevel();

  void monteCarloThirdLevel(size_t numPoints, std::vector<std::vector<real>>& coordinates, std::vector<CombiDataType>& values);

  inline void combineLocalAndGlobal();

  template <typename FG_ELEMENT>
  inline void combineFG(FullGrid<FG_ELEMENT>& fg);

  template <typename FG_ELEMENT>
  inline void gridEval(FullGrid<FG_ELEMENT>& fg);

  /* Generates no_faults random faults from the combischeme */
  inline void createRandomFaults(std::vector<size_t>& faultIds, int no_faults);

  inline void recomputeOptimumCoefficients(std::string prob_name, std::vector<size_t>& faultsID,
                                           std::vector<size_t>& redistributefaultsID,
                                           std::vector<size_t>& recomputeFaultsID);

  inline Task* getTask(size_t taskID);

  void updateCombiParameters();

  void getDSGFromProcessGroup();

  /* Computes group faults in current combi scheme step */
  void getGroupFaultIDs(std::vector<size_t>& faultsID,
                        std::vector<ProcessGroupManagerID>& groupFaults);

  inline CombiParameters& getCombiParameters();

  void parallelEval(const LevelVector& leval, std::string& filename, size_t groupID);

  std::map<size_t, double> getLpNorms(int p = 2);

  std::vector<double> parallelEvalNorm(const LevelVector& leval, size_t groupID = 0);

  std::vector<double> evalAnalyticalOnDFG(const LevelVector& leval, size_t groupID = 0);

  std::vector<double> evalErrorOnDFG(const LevelVector& leval, size_t groupID = 0);

  std::vector<CombiDataType> interpolateValues(const std::vector<std::vector<real>>& interpolationCoords);

  void writeSparseGridMinMaxCoefficients(const std::string& filename);

  void redistribute(std::vector<size_t>& taskID);

  void reInitializeGroup(std::vector<ProcessGroupManagerID>& taskID,
                         std::vector<size_t>& tasksToIgnore);

  void recompute(std::vector<size_t>& taskID, bool failedRecovery,
                 std::vector<ProcessGroupManagerID>& recoveredGroups);

  void recover(int i, int nsteps);

  bool recoverCommunicators(std::vector<ProcessGroupManagerID> failedGroups);
  /* After faults have been fixed, we need to return the combischeme
   * to the original combination technique*/
  void restoreCombischeme();

  void setupThirdLevel();

  /**
   * Call to perform a rescheduling using the given rescheduler and load model.
   *
   * The rescheduling removes tasks from one process group and assigns them to
   * a different process group. The result of the combination is used to 
   * restore values of the newly assigned task.
   * Implications: 
   * - Should only be called after the combination step and before runnext.
   * - Accuracy of calculated values is lost if leval is not equal to 0.
   */
  void reschedule();

  void writeCombigridsToVTKPlotFile(ProcessGroupManagerID pg);

 private:
  ProcessGroupManagerContainer& pgroups_;

  TaskContainer& tasks_;

  CombiParameters params_;

  std::unique_ptr<LoadModel> loadModel_;

  std::unique_ptr<TaskRescheduler> rescheduler_;

  std::map<LevelVector, unsigned long> levelVectorToLastTaskDuration_ = {};

  ThirdLevelUtils thirdLevel_;

  ProcessGroupManagerID& thirdLevelPGroup_;

  // periodically checks status of all process groups. returns until at least
  // one group is in WAIT state
  inline ProcessGroupManagerID wait();
  inline ProcessGroupManagerID waitAvoid(std::vector<ProcessGroupManagerID>& avoidGroups);
  bool waitAllFinished();
  bool waitForPG(ProcessGroupManagerID pg);

  void receiveDurationsOfTasksFromGroupMasters(size_t numDurationsToReceive);

  void sortTasks();
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

  // tell groups to combine local and global
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    bool success = pgroups_[i]->combineLocalAndGlobal();
    assert(success);
  }

  waitAllFinished();
}

/** Combination with third level parallelism e.g. between two HPC systems
*
* The process manager induces a local and global combination first.
* Then he signals ready to the third level manager who decides which system
* sends and receives first. All pgs which do not participate in the third level
* combination directly idle in a broadcast function and wait for their update
* from the third level pg.
*
* Different roles of the manager:
*
* Senders role:
* The processGroupManager transfers the dsgus from the workers of the third
* level pg to the third level manager, who further forwards them to the remote
* system. After sending, he receives the remotely reduced data and sends it back
* to the third level pg.
*
* Receivers role:
* In the role of the receiver, the ProcessGroupManager receives the remote dsgus
* and reduces it with the local solution.
* Afterwards, he sends the solution back to the remote system and to the local
* workers.
*/
void ProcessManager::combineThirdLevel() {
  // first combine local and global
  Stats::startEvent("manager combine local and global");
  combineLocalAndGlobal();
  waitAllFinished();
  Stats::stopEvent("manager combine local and global");

  // tell other pgroups to idle and wait for the combination result
  for (auto& pg : pgroups_) {
    if (pg != thirdLevelPGroup_)
      pg->waitForThirdLevelCombiResult();
  }
  // obtain instructions from third level manager
  thirdLevel_.signalReadyToCombine();
  std::string instruction = thirdLevel_.fetchInstruction();

  // combine
  Stats::startEvent("manager exchange data with remote");
  if (instruction == "send_first") {
    thirdLevelPGroup_->combineThirdLevel(thirdLevel_, params_, true);
  } else if (instruction == "recv_first") {
    thirdLevelPGroup_->combineThirdLevel(thirdLevel_, params_, false);
  }
  Stats::stopEvent("manager exchange data with remote");
  thirdLevel_.signalReady();

  waitAllFinished();
}

/** Unifys the subspace sizes of all dsgus which are collectively combined
 * during third level reduce
 *
 * First, the processGroupManager collects the subspace sizes from all workers
 * dsgus. This is achieved in a single MPI_Gatherv call. The sizes of the send
 * buffers are gathered beforehand. Afterwards, the process manager signals
 * ready to the third level manager who then decides which system sends and
 * receives first.
 *
 * Senders role: The processGroupManager sends and receives data from the third
 * level manager.
 *
 * Receivers role: The ProcessGroupManager receives and sends data to the third
 * level manager.
 *
 * In both roles the manager locally reduces the data and scatters the updated
 * sizes back to the workers of the third level pg who will then distribute it
 * to the other pgs.
 */
void ProcessManager::unifySubspaceSizesThirdLevel() {
  // tell other pgroups to idle and wait for update
  for (auto& pg : pgroups_) {
    if (pg != thirdLevelPGroup_)
      pg->waitForThirdLevelSizeUpdate();
  }

  // obtain instructions from third level manager
  thirdLevel_.signalReadyToUnifySubspaceSizes();
  std::string instruction = thirdLevel_.fetchInstruction();

  // exchange sizes with remote
  if (instruction == "send_first") {
    thirdLevelPGroup_->reduceLocalAndRemoteSubspaceSizes(thirdLevel_, params_, true);
  } else if (instruction == "recv_first") {
    thirdLevelPGroup_->reduceLocalAndRemoteSubspaceSizes(thirdLevel_, params_, false);
  }
  thirdLevel_.signalReady();

  waitAllFinished();

  const auto& dsguDataSizePerWorker = thirdLevelPGroup_->getDsguDataSizePerWorker();

  auto dsguDataSize = std::accumulate(dsguDataSizePerWorker.begin(), dsguDataSizePerWorker.end(), 0);
  Stats::setAttribute("dsguDataSize", std::to_string(dsguDataSize));
}

/** This function performs the so-called recombination. First, the combination
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

/** Evaluate the combination solution with the resolution of the given full grid.
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

/**
 * Create a certain given number of random faults, considering that the faulty processes
 * simply cannot give the evaluation results, but they are still available in the MPI
 * communication scheme (the nodes are not dead)
 */
inline void ProcessManager::createRandomFaults(std::vector<size_t>& faultIds, int no_faults) {
  size_t fault_id;

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

/**
 * Recompute coefficients for the combination technique based on given grid faults using
 * an optimization scheme
 */
inline void ProcessManager::recomputeOptimumCoefficients(std::string prob_name,
                                                         std::vector<size_t>& faultsID,
                                                         std::vector<size_t>& redistributeFaultsID,
                                                         std::vector<size_t>& recomputeFaultsID) {
  CombigridDict given_dict = params_.getCombiDict();

  std::map<size_t, LevelVector> IDsToLevels = params_.getLevelsDict();
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
    std::vector<size_t> newTaskIDs;

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
    int roundedSum = static_cast<int>(round(sum));
    std::cout << "Coefficient sum: " << roundedSum << "\n";

    assert(roundedSum == 1);
    params_.setLevelsCoeffs(newTaskIDs, newLevels, newCoeffs);

    std::map<LevelVector, size_t> LevelsToIDs = params_.getLevelsToIDs();
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

inline Task* ProcessManager::getTask(size_t taskID) {
  for (Task* tmp : tasks_) {
    if (tmp->getID() == taskID) {
      return tmp;
    }
  }
  return nullptr;
}
} /* namespace combigrid */
#endif /* PROCESSMANAGER_HPP_ */
