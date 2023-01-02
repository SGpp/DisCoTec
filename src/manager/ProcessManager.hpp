#ifndef PROCESSMANAGER_HPP_
#define PROCESSMANAGER_HPP_

#include <vector>

#include "combischeme/CombiMinMaxScheme.hpp"
#include "fault_tolerance/LPOptimizationInterpolation.hpp"
#include "manager/ProcessGroupManager.hpp"
#include "manager/ProcessGroupSignals.hpp"
#include "loadmodel/LoadModel.hpp"
#include "loadmodel/LearningLoadModel.hpp"
#include "rescheduler/TaskRescheduler.hpp"
#include "rescheduler/StaticTaskRescheduler.hpp"
#include "mpi/MPISystem.hpp"
#include "task/Task.hpp"

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
      rescheduler_{std::move(rescheduler)}
  { }

  inline void removeGroups(std::vector<int> removeIndices);

  // todo: use general class AppInstance here
  // todo: add remove function
  inline void addTask(Task* t);

  bool runfirst(bool doInitDSGUs = true);

  void initDsgus();

  void exit();

  virtual ~ProcessManager();

  void waitForAllGroupsToWait() const;

  template <typename FG_ELEMENT>
  inline FG_ELEMENT eval(const std::vector<real>& coords);

  bool runnext();

  inline void combine();

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

  /* Computes group faults in current combi scheme step */
  void getGroupFaultIDs(std::vector<size_t>& faultsID,
                        std::vector<ProcessGroupManagerID>& groupFaults);

  inline CombiParameters& getCombiParameters();

  void parallelEval(const LevelVector& leval, std::string& filename, size_t groupID);

  void doDiagnostics(size_t taskID);

  std::map<size_t, double> getLpNorms(int p = 2);

  std::vector<double> parallelEvalNorm(const LevelVector& leval, size_t groupID = 0);

  std::vector<double> evalAnalyticalOnDFG(const LevelVector& leval, size_t groupID = 0);

  std::vector<double> evalErrorOnDFG(const LevelVector& leval, size_t groupID = 0);

  std::vector<CombiDataType> interpolateValues(const std::vector<std::vector<real>>& interpolationCoords);

  void writeInterpolatedValues(const std::vector<std::vector<real>>& interpolationCoords);

  void writeInterpolationCoordinates(const std::vector<std::vector<real>>& interpolationCoords);

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

  void writeDSGsToDisk(std::string filenamePrefix);

  void readDSGsFromDisk(std::string filenamePrefix);

 private:
  ProcessGroupManagerContainer& pgroups_;

  TaskContainer& tasks_;

  CombiParameters params_;

  std::unique_ptr<LoadModel> loadModel_;

  std::unique_ptr<TaskRescheduler> rescheduler_;

  std::map<LevelVector, unsigned long> levelVectorToLastTaskDuration_ = {};

  // periodically checks status of all process groups. returns until at least
  // one group is in WAIT state
  inline ProcessGroupManagerID wait();
  inline ProcessGroupManagerID waitAvoid(std::vector<ProcessGroupManagerID>& avoidGroups);
  bool waitAllFinished();

  void receiveDurationsOfTasksFromGroupMasters(size_t numDurationsToReceive);

  void sortTasks();

  ProcessGroupManagerID getProcessGroupWithTaskID(size_t taskID){
    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->hasTask(taskID)){
        return pgroups_[i];
      }
    }
    return nullptr;
  }
};

inline void ProcessManager::addTask(Task* t) {
  tasks_.push_back(t);
  // wait for available process group
  ProcessGroupManagerID g = wait();
  g->addTask(t);
}

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
  waitForAllGroupsToWait();

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
  waitForAllGroupsToWait();

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
template <typename FG_ELEMENT>
void ProcessManager::combineFG(FullGrid<FG_ELEMENT>& fg) {
  waitForAllGroupsToWait();

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
  waitForAllGroupsToWait();

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

/*
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
