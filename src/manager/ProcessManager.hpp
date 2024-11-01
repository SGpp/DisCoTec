#ifndef PROCESSMANAGER_HPP_
#define PROCESSMANAGER_HPP_

#include <numeric>
#include <vector>

#include "combischeme/CombiMinMaxScheme.hpp"
#include "fault_tolerance/LPOptimizationInterpolation.hpp"
#include "loadmodel/LearningLoadModel.hpp"
#include "loadmodel/LoadModel.hpp"
#include "manager/ProcessGroupManager.hpp"
#include "manager/ProcessGroupSignals.hpp"
#include "mpi/MPISystem.hpp"
#include "rescheduler/StaticTaskRescheduler.hpp"
#include "rescheduler/TaskRescheduler.hpp"
#include "task/Task.hpp"
#include "third_level/ThirdLevelUtils.hpp"
#include "utils/MonteCarlo.hpp"

namespace combigrid {

/**
 * @brief The ProcessManager class orchestrates the whole simulation in case of a manager-worker
 * scheme in DisCoTec.
 *
 * It should only be instantiated once---on the manager rank---and it holds one ProcessGroupManager
 * instance for communication with each of the process groups.
 *
 * All other ranks should instantiate a ProcessGroupWorker and call wait() on the worker, until an
 * exit signal is received.
 */
class ProcessManager {
 public:
  /**
   * @brief Constructor
   *
   * @param pgroups a vector of ProcessGroupManager s
   * @param instances a vector of Task s
   * @param params the combination parameters
   * @param loadModel The load model to use for (re)scheduling
   * @param rescheduler The rescheduler to use for dynamic task rescheduling.
   *                    By default, the static task rescheduler is used and
   *                    therefore no rescheduling perfomed.
   */
  ProcessManager(ProcessGroupManagerContainer& pgroups, TaskContainer& instances,
                 CombiParameters& params, std::unique_ptr<LoadModel> loadModel,
                 std::unique_ptr<TaskRescheduler> rescheduler =
                     std::unique_ptr<TaskRescheduler>(new StaticTaskRescheduler{}))
      : pgroups_{pgroups},
        tasks_{instances},
        params_{params},
        loadModel_{std::move(loadModel)},
        rescheduler_{std::move(rescheduler)},
        thirdLevel_{params.getThirdLevelHost(), params.getThirdLevelPort()},
        thirdLevelPGroup_{pgroups_[params.getThirdLevelPG()]} {
    // only setup third level if explicitly desired
    if (params.getThirdLevelHost() != "") setupThirdLevel();
  }

  /**
   * @brief Removes the process groups with the given indices from the simulation
   *
   * Used for fault tolerance
   */
  inline void removeGroups(std::vector<int> removeIndices);

  /**
   * @brief signal to run the first combination step, i.e. initialize and run each task
   *
   * this is where initial load balancing is done: tasks are sorted by expected runtime and assigned
   * to process groups as they become available again
   *
   * @param doInitDSGUs whether to initialize the DSGUs after the first combination step
   * @return true if no group failed
   */
  bool runfirst(bool doInitDSGUs = true);

  /**
   * @brief signal to initialize the sparse grid data structures on the worker ranks
   */
  void initDsgus();

  /**
   * @brief signal to exit the simulation
   */
  void exit();

  virtual ~ProcessManager() = default;

  /**
   * @brief wait until all groups have signaled completion on the last signal
   */
  void waitForAllGroupsToWait() const;

  /**
   * @brief signal to run a combination step after the first one
   */
  bool runnext();

  /**
   * @brief signal to combine the results of the tasks, according to the CombiParameters
   *
   * This function performs the so-called recombination. First, the combination
   * solution will be reduced in the given sparse grid space (first within, then across process
   * groups). Also, the local component grids will be updated from the globally combined solution.
   */
  inline void combine();

  /**
   * @brief signal to perform a widely-distributed combination
   *
   * based on TCP/socket setup with third level manager.
   *
   * Combination with third level parallelism e.g. between two HPC systems:
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
  inline void combineThirdLevel();

  /**
   * @brief signal to start a widely-distributed combination
   *
   * based on file-exchange mechanism (w/o third level manager)
   */
  inline void combineThirdLevelFileBasedWrite(const std::string& filenamePrefixToWrite,
                                              const std::string& writeCompleteTokenFileName);

  /**
   * @brief signal to reduce the results of a widely-distributed combination
   *
   * based on file-exchange mechanism (w/o third level manager)
   */
  inline void combineThirdLevelFileBasedReadReduce(const std::string& filenamePrefixToRead,
                                                   const std::string& startReadingTokenFileName);

  /**
   * @brief signal to perform a whole widely-distributed combination
   *
   * based on file-exchange mechanism (w/o third level manager); equivalent to calling
   * combineThirdLevelFileBasedWrite and combineThirdLevelFileBasedReadReduce in succession
   */
  inline void combineThirdLevelFileBased(const std::string& filenamePrefixToWrite,
                                         const std::string& writeCompleteTokenFileName,
                                         const std::string& filenamePrefixToRead,
                                         const std::string& startReadingTokenFileName);

  /**
   * @brief signal to pretend a widely-distributed combination
   *
   * based on TCP/socket setup with third level manager; for testing the widely-distributed
   * combination between the third level manager and the workers in the third level process group
   *
   *  like combineThirdLevel, but without involving any process groups -- sending dummy data instead
   */
  inline size_t pretendCombineThirdLevelForBroker(std::vector<long long> numDofsToCommunicate,
                                                  bool checkValues);

  /**
   * @brief signal to pretend a widely-distributed combination
   *
   * based on TCP/socket setup with third level manager; for testing the widely-distributed
   * combination between the workers in and outside the third level process group
   */
  inline void pretendCombineThirdLevelForWorkers();

  /**
   * @brief signal to reduce the subspace sizes between the systems
   *
   * based on TCP/socket setup with third level manager
   *
   * Unifies the subspace sizes of all dsgus which are collectively combined
   * during third level reduce:
   *
   * First, the processGroupManager collects the subspace sizes from all workers'
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
  inline size_t unifySubspaceSizesThirdLevel(bool thirdLevelExtraSparseGrid);

  /**
   * @brief signal to pretend a reduction of the subspace sizes between the systems
   *
   * based on TCP/socket setup with third level manager; for testing.
   *
   * like unifySubspaceSizesThirdLevel, but without sending any data widely. instead, the manager
   * sends only zeros to the third level group, so it will keep its own sparse grid sizes
   */
  inline size_t pretendUnifySubspaceSizesThirdLevel();

  /**
   * @brief signal to perform a widely-distributed Monte-Carlo interpolation of the current
   * simulation
   *
   * based on TCP/socket setup with third level manager
   */
  void monteCarloThirdLevel(size_t numPoints, std::vector<std::vector<real>>& coordinates,
                            std::vector<CombiDataType>& values);

  /**
   * @brief signal to perform a system-wide (not widely-distributed) combination
   */
  inline void combineSystemWide();

  /**
   * @brief recompute coefficients for the combination technique
   *
   * based on given grid faults using an optimization scheme; used for fault tolerance
   */
  inline void recomputeOptimumCoefficients(std::string prob_name, std::vector<size_t>& faultsID,
                                           std::vector<size_t>& redistributefaultsID,
                                           std::vector<size_t>& recomputeFaultsID);

  /**
   * @brief get a pointer to the task with the given ID
   */
  inline Task* getTask(size_t taskID);

  /**
   * @brief signal to receive the combination parameters and send new ones
   */
  void updateCombiParameters();

  /**
   * @brief Computes group faults in current combi scheme step
   */
  void getGroupFaultIDs(std::vector<size_t>& faultsID,
                        std::vector<ProcessGroupManagerID>& groupFaults);

  /**
   * @brief signal one group to interpolate the current solution from the current sparse grid at
   * resolution level \p leval
   *
   * writes the solution to a binary file readable with Paraview
   */
  void parallelEval(const LevelVector& leval, std::string& filename, size_t groupID);

  /**
   * @brief signal to perform diagnostics on the task with the given ID
   *
   * can only be used with Tasks that implement the doDiagnostics method
   */
  void doDiagnostics(size_t taskID);

  /**
   * @brief signal to compute the Lp norm of the current component grids, and gather them
   *
   * @param p the p in Lp norm
   * @return a map from task ID to Lp norm
   */
  std::map<size_t, double> getLpNorms(int p = 2);

  /**
   * @brief get the Lp norm of the current combined solution from workers
   */
  double getLpNorm(int p = 2);

  /**
   * @brief signal to interpolate the current solution on all component grids at the given \p
   * interpolationCoords
   *
   * requires that the component grids are in nodal representation (not hierarchized).
   * The results are sent to the manager rank and returned here.
   */
  std::vector<CombiDataType> interpolateValues(
      const std::vector<std::vector<real>>& interpolationCoords);

  /**
   * @brief signal to interpolate at the given \p interpolationCoords and write results to file
   *
   * like interpolateValues, but the last process group writes the results to a file
   */
  void writeInterpolatedValuesSingleFile(const std::vector<std::vector<real>>& interpolationCoords,
                                         const std::string& filenamePrefix);

  /**
   * @brief signal to interpolate at the given \p interpolationCoords at each grid and write results
   * to one file per grid
   */
  void writeInterpolatedValuesPerGrid(const std::vector<std::vector<real>>& interpolationCoords,
                                      const std::string& filenamePrefix);

  /**
   * @brief write the interpolation coordinates to a file
   */
  void writeInterpolationCoordinates(const std::vector<std::vector<real>>& interpolationCoords,
                                     const std::string& filenamePrefix) const;

  /**
   * @brief signal the last group to write minimum and maximum subspace coefficients to a file
   *
   * @param filename the filename to write to
   */
  void writeSparseGridMinMaxCoefficients(const std::string& filename);

  /**
   * @brief assign tasks to available process groups
   *
   * used for fault tolerance: if a process group fails, its tasks are redistributed to other groups
   */
  void redistribute(std::vector<size_t>& taskID);

  /**
   * @brief signal to reinitialize the group with the given task IDs
   *
   * used for fault tolerance
   */
  void reInitializeGroup(std::vector<ProcessGroupManagerID>& taskID,
                         std::vector<size_t>& tasksToIgnore);

  /**
   * @brief signal to recompute the given task IDs on some recovered groups
   *
   * used for fault tolerance; the tasks will be re-initialized from the current sparse grid
   * solution
   */
  void recompute(std::vector<size_t>& taskID, bool failedRecovery,
                 std::vector<ProcessGroupManagerID>& recoveredGroups);

  bool recoverCommunicators(std::vector<ProcessGroupManagerID> failedGroups);
  /* After faults have been fixed, we need to return the combischeme
   * to the original combination technique*/
  void restoreCombischeme();

  /**
   * @brief establish connection to third level manager
   *
   * based on TCP/socket setup with third level manager
   */
  void setupThirdLevel();

  /**
   * @brief perform rescheduling using the given rescheduler and load model.
   *
   * The rescheduling removes tasks from one process group and assigns them to
   * a different process group. The result of the combination is used to
   * restore values of the newly assigned task.
   * Implications:
   * - Should only be called after the combination step and before runnext.
   * - Accuracy of calculated values is lost if leval is not equal to 0.
   */
  void reschedule();

  /**
   * @brief signal a single group to write the component grids to vtk plot file
   */
  void writeCombigridsToVTKPlotFile(ProcessGroupManagerID pg);

  /**
   * @brief signal all groups to write their sparse grid data structures to disk
   */
  void writeDSGsToDisk(const std::string& filenamePrefix);

  /**
   * @brief signal all groups to read their sparse grid data structures from disk
   */
  void readDSGsFromDisk(const std::string& filenamePrefix);

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

  ProcessGroupManagerID getProcessGroupWithTaskID(size_t taskID) {
    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->hasTask(taskID)) {
        return pgroups_[i];
      }
    }
    return nullptr;
  }
};

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

void ProcessManager::combine() {
  waitForAllGroupsToWait();

  // send signal to each group
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    [[maybe_unused]] bool success = pgroups_[i]->combine();
    assert(success);
  }

  waitAllFinished();
}

void ProcessManager::combineSystemWide() {
  Stats::startEvent("manager combine local");
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }

  // tell groups to combine local and global
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    [[maybe_unused]] bool success = pgroups_[i]->combineSystemWide();
    assert(success);
  }

  waitAllFinished();
  Stats::stopEvent("manager combine local");
}

void ProcessManager::combineThirdLevel() {
  // first combine local and global
  combineSystemWide();

  // tell other pgroups to idle and wait for the combination result
  for (auto& pg : pgroups_) {
    if (pg != thirdLevelPGroup_) pg->waitForThirdLevelCombiResult();
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

void ProcessManager::combineThirdLevelFileBasedWrite(
    const std::string& filenamePrefixToWrite, const std::string& writeCompleteTokenFileName) {
  // first combine local and global
  combineSystemWide();

  // obtain "instructions" from third level manager
  thirdLevel_.signalReadyToCombineFile();
  std::string instruction = thirdLevel_.fetchInstruction();
  assert(instruction == "write_ok");

  // combine
  Stats::startEvent("manager combine write");
  thirdLevelPGroup_->combineThirdLevelFileBasedWrite(filenamePrefixToWrite,
                                                     writeCompleteTokenFileName);
  Stats::stopEvent("manager combine write");
  waitForPG(thirdLevelPGroup_);
}

void ProcessManager::combineThirdLevelFileBasedReadReduce(
    const std::string& filenamePrefixToRead, const std::string& startReadingTokenFileName) {
  // integrate the solutions
  Stats::startEvent("manager combine read");
  // tell other pgroups to wait for the combination result
  for (auto& pg : pgroups_) {
    if (pg != thirdLevelPGroup_) pg->waitForThirdLevelCombiResult();
  }

  thirdLevelPGroup_->combineThirdLevelFileBasedReadReduce(filenamePrefixToRead,
                                                          startReadingTokenFileName);

  thirdLevel_.signalReady();
  waitAllFinished();
  Stats::stopEvent("manager combine read");
}

void ProcessManager::combineThirdLevelFileBased(const std::string& filenamePrefixToWrite,
                                                const std::string& writeCompleteTokenFileName,
                                                const std::string& filenamePrefixToRead,
                                                const std::string& startReadingTokenFileName) {
  // first combine local and global
  combineSystemWide();

  // tell other pgroups to idle and wait for the combination result
  for (auto& pg : pgroups_) {
    if (pg != thirdLevelPGroup_) pg->waitForThirdLevelCombiResult();
  }
  // obtain instructions from third level manager
  thirdLevel_.signalReadyToCombineFile();
  std::string instruction = thirdLevel_.fetchInstruction();
  assert(instruction == "write_ok");

  // combine
  Stats::startEvent("manager combine file");
  thirdLevelPGroup_->combineThirdLevelFileBased(filenamePrefixToWrite, writeCompleteTokenFileName,
                                                filenamePrefixToRead, startReadingTokenFileName);
  Stats::stopEvent("manager combine file");
  waitForPG(thirdLevelPGroup_);
  thirdLevel_.signalReady();
}

void ProcessManager::pretendCombineThirdLevelForWorkers() {
  // first combine local and global
  combineSystemWide();

  // combine
  Stats::startEvent("manager exchange no data with remote");
  // tell other pgroups to idle and wait for the combination result
  for (auto& pg : pgroups_) {
    if (pg != thirdLevelPGroup_) pg->waitForThirdLevelCombiResult();
  }
  thirdLevelPGroup_->pretendCombineThirdLevelForWorkers(params_);
  waitAllFinished();
  Stats::stopEvent("manager exchange no data with remote");
}

size_t ProcessManager::pretendCombineThirdLevelForBroker(
    std::vector<long long> numDofsToCommunicate, bool checkValues) {
  size_t numWrongValues = 0;
  // obtain instructions from third level manager
  thirdLevel_.signalReadyToCombine();
  std::string instruction = thirdLevel_.fetchInstruction();

  std::cout << " pretend " << instruction << std::endl;

  Stats::startEvent("manager exchange data with remote");
  for (const auto& dsguSize : numDofsToCommunicate) {
    std::vector<CombiDataType> dsguData(dsguSize, 0.);
    // combine
    if (instruction == "send_first") {
      // if sending first, initialize with random
      auto initialData = montecarlo::getRandomCoordinates(1, dsguSize)[0];
      dsguData.assign(initialData.begin(), initialData.end());
      // send dsgu to remote
      thirdLevel_.sendData(dsguData.data(), dsguSize);
      // recv combined dsgu from remote
      thirdLevel_.recvData(dsguData.data(), dsguSize);
      if (checkValues) {
        for (long long j = 0; j < dsguSize; ++j) {
          if (dsguData[j] != initialData[j]) {
            ++numWrongValues;
          }
        }
        assert(numWrongValues == 0);
      }
    } else if (instruction == "recv_first") {
      // recv and combine dsgu from remote
      thirdLevel_.recvAndAddToData(dsguData.data(), dsguSize);
      // send combined solution to remote
      thirdLevel_.sendData(dsguData.data(), dsguSize);
    }
  }
  Stats::stopEvent("manager exchange data with remote");
  thirdLevel_.signalReady();
  return numWrongValues;
}

size_t ProcessManager::unifySubspaceSizesThirdLevel(bool thirdLevelExtraSparseGrid) {
  if (!thirdLevelExtraSparseGrid) {
    // tell other pgroups to idle and wait for update
    for (auto& pg : pgroups_) {
      if (pg != thirdLevelPGroup_) pg->waitForThirdLevelSizeUpdate();
    }
  }

  // obtain instructions from third level manager
  thirdLevel_.signalReadyToUnifySubspaceSizes();
  std::string instruction = thirdLevel_.fetchInstruction();

  // exchange sizes with remote
  if (instruction == "send_first") {
    thirdLevelPGroup_->reduceLocalAndRemoteSubspaceSizes(thirdLevel_, true,
                                                         thirdLevelExtraSparseGrid);
  } else if (instruction == "recv_first") {
    thirdLevelPGroup_->reduceLocalAndRemoteSubspaceSizes(thirdLevel_, false,
                                                         thirdLevelExtraSparseGrid);
  }
  thirdLevel_.signalReady();

  waitAllFinished();

  const auto& formerDsguDataSizePerWorker = thirdLevelPGroup_->getFormerDsguDataSizePerWorker();
  const auto& dsguDataSizePerWorker = thirdLevelPGroup_->getDsguDataSizePerWorker();

  auto formerDsguDataSize =
      std::accumulate(formerDsguDataSizePerWorker.begin(), formerDsguDataSizePerWorker.end(), 0);
  auto dsguDataSize =
      std::accumulate(dsguDataSizePerWorker.begin(), dsguDataSizePerWorker.end(), 0);
  Stats::setAttribute("formerDsguDataSize", std::to_string(formerDsguDataSize));
  Stats::setAttribute("dsguDataSize", std::to_string(dsguDataSize));

  if (thirdLevelExtraSparseGrid && dsguDataSize > formerDsguDataSize) {
    throw std::runtime_error("why is it larger???");
  } else if (!thirdLevelExtraSparseGrid && dsguDataSize < formerDsguDataSize) {
    throw std::runtime_error("why is it smaller???");
  }

  return dsguDataSize;
}

size_t ProcessManager::pretendUnifySubspaceSizesThirdLevel() {
  // tell other pgroups to idle and wait for update
  for (auto& pg : pgroups_) {
    if (pg != thirdLevelPGroup_) pg->waitForThirdLevelSizeUpdate();
  }

  // exchange sizes with remote
  thirdLevelPGroup_->pretendReduceLocalAndRemoteSubspaceSizes(thirdLevel_);

  waitAllFinished();

  const auto& formerDsguDataSizePerWorker = thirdLevelPGroup_->getFormerDsguDataSizePerWorker();
  const auto& dsguDataSizePerWorker = thirdLevelPGroup_->getDsguDataSizePerWorker();

  auto formerDsguDataSize =
      std::accumulate(formerDsguDataSizePerWorker.begin(), formerDsguDataSizePerWorker.end(), 0);
  auto dsguDataSize =
      std::accumulate(dsguDataSizePerWorker.begin(), dsguDataSizePerWorker.end(), 0);
  Stats::setAttribute("formerDsguDataSize", std::to_string(formerDsguDataSize));
  Stats::setAttribute("dsguDataSize", std::to_string(dsguDataSize));

  if (dsguDataSize != formerDsguDataSize) {
    throw std::runtime_error("wrong number of dofs pretending unifying dsgu sizes");
  }

  return dsguDataSize;
}

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
