#ifndef PROCESSGROUPMANAGER_HPP_
#define PROCESSGROUPMANAGER_HPP_

#include "../../../../../.cache/JetBrains/CLion2024.1/.remote/ipvs-epyc2_22/347dacd9-0950-4b02-b7f7-b2d9ba9d9d1d/usr/include/c++/9/stdexcept"
#include "../../../../../.cache/JetBrains/CLion2024.1/.remote/ipvs-epyc2_22/347dacd9-0950-4b02-b7f7-b2d9ba9d9d1d/usr/include/c++/9/string"
#include "../../../../../.cache/JetBrains/CLion2024.1/.remote/ipvs-epyc2_22/347dacd9-0950-4b02-b7f7-b2d9ba9d9d1d/usr/include/c++/9/vector"

#include "combicom/CombiCom.hpp"
#include "../fullgrid/FullGrid.hpp"
#include "CombiParameters.hpp"
#include "ProcessGroupSignals.hpp"
#include "../mpi/MPISystem.hpp"
#include "../MPI-FT.h"
#include "../Task.hpp"
#include "../third_level/ThirdLevelUtils.hpp"
#include "../utils/Types.hpp"

namespace combigrid {
/**
 * @class ProcessGroupManager
 * @brief The ProcessGroupManager is part of a ProcessManager and is responsible for communication
 * with a single process group.
 *
 * Through the use of non-blocking operations in the ProcessGroupManager, the ProcessManager can
 * instruct multiple process groups at once.
 */
class ProcessGroupManager {
  friend class ProcessManager;

 public:
  /**
   * Constructor
   *
   * @param pgroupRootID the rank of each process of this group in the global communicator;
   * this is the same as the process group number, as the global communicator contains one member of
   * each process group
   */
  explicit ProcessGroupManager(RankType pgroupRootID);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  ProcessGroupManager() = delete;
  ~ProcessGroupManager() = default;
  // rule of 5
  ProcessGroupManager(ProcessGroupManager const&) = delete;
  ProcessGroupManager& operator=(ProcessGroupManager const&) = delete;
  ProcessGroupManager(ProcessGroupManager&&) = delete;
  ProcessGroupManager& operator=(ProcessGroupManager&&) = delete;
#endif  // DOXYGEN_SHOULD_SKIP_THIS

  /**
   * @brief signal the process group to initalize and run the task \p t
   */
  bool runfirst(Task* t);

  /**
   * @brief signal the process group to run a time step on all of its tasks
   */
  bool runnext();

  /**
   * @brief signal the process group to initialized its sparse grid data structures
   */
  bool initDsgus();

  /**
   * @brief signal the process group to exit
   */
  bool exit();

  /**
   * @brief non-blocking call to retrieve status of process group
   */
  inline StatusType getStatus();

  /**
   * @brief blocks until process group finished computation
   */
  inline StatusType waitStatus();

  /**
   * @brief get a collection of all tasks currently assigned to this group
   */
  inline const TaskContainer& getTaskContainer() const;

  /**
   * @brief remove a task from the process group
   *
   * this does not change the state in the workers!
   */
  inline void removeTask(Task* t);

  /**
   * @brief signal to perform a system-wide combination
   */
  bool combine();

  /**
   * @brief signal to perform a widely-distributed combination
   *
   * based on TCP/socket setup with third level manager. cf. ProcessManager::combineThirdLevel
   */
  bool combineThirdLevel(const ThirdLevelUtils& thirdLevel, CombiParameters& params,
                         bool isSendingFirst);
  /**
   * @brief signal to perform a whole widely-distributed combination
   *
   * based on file-exchange mechanism (w/o third level manager); equivalent to calling
   * combineThirdLevelFileBasedWrite and combineThirdLevelFileBasedReadReduce in succession
   */
  bool combineThirdLevelFileBased(const std::string& filenamePrefixToWrite,
                                  const std::string& writeCompleteTokenFileName,
                                  const std::string& filenamePrefixToRead,
                                  const std::string& startReadingTokenFileName);

  /**
   * @brief signal to start a widely-distributed combination
   *
   * based on file-exchange mechanism (w/o third level manager)
   */
  bool combineThirdLevelFileBasedWrite(const std::string& filenamePrefixToWrite,
                                       const std::string& writeCompleteTokenFileName);

  /**
   * @brief signal to reduce the results of a widely-distributed combination
   *
   * based on file-exchange mechanism (w/o third level manager)
   */
  bool combineThirdLevelFileBasedReadReduce(const std::string& filenamePrefixToRead,
                                            const std::string& startReadingTokenFileName);

  /**
   * @brief signal to pretend a widely-distributed combination
   *
   * based on TCP/socket setup with third level manager; for testing the widely-distributed
   * combination between the workers in and outside the third level process group
   */
  bool pretendCombineThirdLevelForWorkers(CombiParameters& params);

  /**
   * @brief signal to perform the first part of a system-wide combination, namely hierarchization
   * and reduction (but not dehierarchization)
   *
   * based on file-exchange mechanism (w/o third level manager)
   */
  bool combineSystemWide();

  /**
   * @brief send new CombiParameters to the process group
   */
  bool updateCombiParameters(CombiParameters& params);

  /**
   * @brief Check if group fault occured at this combination step using the fault simulator
   */
  bool isGroupFault();

  /**
   * @brief assign a task to the process group
   */
  bool addTask(Task*);
  bool refreshTask(Task*);

  /**
   * @brief signal to delete all tasks in the process group
   *
   * does not change the state in this ProcessGroupManager!
   */
  bool resetTasksWorker();

  /**
   * @brief assign a task to the process group, and let the workers run a time step
   *
   * used for fault tolerance; the task will be re-initialized from the current sparse grid solution
   */
  bool recompute(Task* t);

  /**
   * @brief signal to recover the communicator ranks of the process group
   */
  bool recoverCommunicators();

  /**
   * @brief signal to interpolate the current solution from all component grids on this group at
   * resolution level \p leval and write the results to a binary file readable with Paraview
   */
  bool parallelEval(const LevelVector& leval, const std::string& filename);

  /**
   * @brief signal the group to write minimum and maximum subspace coefficients to a file
   *
   * @param filename the filename to write to
   */
  void writeSparseGridMinMaxCoefficients(const std::string& filename);

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
  void getLpNorms(int p, std::map<size_t, double>& norms);

  /**
   * @brief signal to interpolate the analytical solution at the given resolution level \p leval
   */
  std::vector<double> evalAnalyticalOnDFG(const LevelVector& leval);

  /**
   * @brief signal to interpolate the analytical solution at the given resolution level \p leval and
   * compute the difference to the current solution
   *
   * @return the Lp norms of the error: maximum norm, l1 norm, l2 norm
   */
  std::vector<double> evalErrorOnDFG(const LevelVector& leval);

  /**
   * @brief signal to interpolate the current component grids at the given \p
   * interpolationCoordsSerial
   *
   * non-blocking; the results will be written to values after the requests have completed
   */
  void interpolateValues(const std::vector<real>& interpolationCoordsSerial,
                         std::vector<CombiDataType>& values, MPI_Request* request = nullptr,
                         const std::string& filenamePrefix = "");

  /**
   * @brief signal to interpolate the current component grids at the given \p
   * interpolationCoordsSerial and write results; one file per grid.
   */
  void writeInterpolatedValuesPerGrid(const std::vector<real>& interpolationCoordsSerial,
                                      const std::string& filenamePrefix);

  /**
   * Adds a task to the process group. To be used for rescheduling.
   *
   * @param task The task to add.
   * @returns If task was successfully added.
   */
  bool rescheduleAddTask(Task* task);

  /**
   * Removes a task from the process group. To be used for rescheduling.
   *
   * @param lvlVec The level vector of the task to remove.
   * @returns If successful the removed task or a nullptr if no task with the
   *          given level vector is found.
   */
  Task* rescheduleRemoveTask(const LevelVector& lvlVec);

  /**
   * @brief returns true if the group is currently assigned the task with the given \p taskID
   */
  bool hasTask(size_t taskID) {
    auto foundIt = std::find_if(tasks_.begin(), tasks_.end(),
                                [taskID](Task* t) { return ((t->getID()) == taskID); });
    return foundIt != tasks_.end();
  }

  /**
   * @brief signal to write the group's component grids to vtk plot file each
   */
  bool writeCombigridsToVTKPlotFile();

  /**
   * @brief signal to write the group's sparse grid data structures to disk
   *
   * using custom binary format with MPI-IO (and compression, if enabled)
   *
   * @param filenamePrefix the prefix of the filename to write to
   */
  bool writeDSGsToDisk(const std::string& filenamePrefix);

  /**
   * @brief signal to read the group's sparse grid data structures from disk
   *
   * using custom binary format with MPI-IO (and decompression, if enabled)
   *
   * @param filenamePrefix the prefix of the filename to read from
   */
  bool readDSGsFromDisk(const std::string& filenamePrefix);

  /**
   * @brief store the task for this process group
   *
   * Does not change the state in the workers!
   */
  void storeTaskReference(Task* t);

 private:
  RankType pgroupRootID_;  // rank in GlobalComm of the master process of this group

  TaskContainer tasks_;

  StatusType status_;

  MPI_Request statusRequest_;

  simft::Sim_FT_MPI_Request statusRequestFT_;

  // stores the accumulated dsgu sizes per worker
  std::vector<size_t> formerDsguDataSizePerWorker_;
  std::vector<size_t> dsguDataSizePerWorker_;

  void recvStatus();

  inline void setStatus(StatusType status);

  // Helper functions for Communication with ProcessGroups
  bool storeTaskReferenceAndSendTaskToProcessGroup(Task* t, SignalType signal);

  bool sendTaskToProcessGroup(Task* t, SignalType signal);

  void sendSignalAndReceive(SignalType signal);

  void sendSignalToProcessGroup(SignalType signal);

  inline void setProcessGroupBusyAndReceive();

  void exchangeDsgus(const ThirdLevelUtils& thirdLevel, CombiParameters& params,
                     bool isSendingFirst);

  bool collectSubspaceSizes(const ThirdLevelUtils& thirdLevel, std::vector<SubspaceSizeType>& buff,
                            size_t& buffSize, std::vector<int>& numSubspacesPerWorker);

  bool distributeSubspaceSizes(const ThirdLevelUtils& thirdLevel,
                               const std::vector<SubspaceSizeType>& buff, size_t buffSize,
                               const std::vector<int>& numSubspacesPerWorker);

  bool reduceLocalAndRemoteSubspaceSizes(const ThirdLevelUtils& thirdLevel, bool isSendingFirst,
                                         bool thirdLevelExtraSparseGrid);

  bool pretendReduceLocalAndRemoteSubspaceSizes(const ThirdLevelUtils& thirdLevel);

  const std::vector<size_t>& getFormerDsguDataSizePerWorker() {
    return formerDsguDataSizePerWorker_;
  }

  const std::vector<size_t>& getDsguDataSizePerWorker() { return dsguDataSizePerWorker_; }

  bool waitForThirdLevelSizeUpdate();

  bool waitForThirdLevelCombiResult();

  /* sets the rank of the process group's master in global comm. should only
   * be called by ProcessManager.
   */
  friend class ProcessManager;
  inline void setMasterRank(int pGroupRootID);
  inline int getMasterRank();
};

typedef std::shared_ptr<ProcessGroupManager> ProcessGroupManagerID;
typedef std::vector<ProcessGroupManagerID> ProcessGroupManagerContainer;

inline StatusType ProcessGroupManager::getStatus() {
  if (status_ == PROCESS_GROUP_WAIT) return PROCESS_GROUP_WAIT;

  if (status_ == PROCESS_GROUP_FAIL) return PROCESS_GROUP_FAIL;
  //  std::cout << "status before test is " << status_ << " \n";

  // if the process group is busy we need
  if (status_ == PROCESS_GROUP_BUSY) {
    if (ENABLE_FT) {
      int flag = -1;
      simft::Sim_FT_MPI_Status stat;
      int err = simft::Sim_FT_MPI_Test(&statusRequestFT_, &flag, &stat);

      if (err == MPI_ERR_PROC_FAILED) status_ = PROCESS_GROUP_FAIL;
    } else {
      /* todo: actually MPI_TEST is not really necessary here. However, i think it
       * might be a good idea to have this here. MPI_Test might run a system
       * call which enables the OS to switch to the MPI system.
       */
      int flag;
      MPI_Test(&statusRequest_, &flag, MPI_STATUS_IGNORE);
    }
  }
  //  std::cout << "status is " << status_ << " \n";
  assert(status_ >= 0);  // check for invalid values
  assert(status_ <= 2);
  return status_;
}
inline void ProcessGroupManager::setStatus(StatusType status) { status_ = status; }

inline StatusType ProcessGroupManager::waitStatus() {
  if (status_ == PROCESS_GROUP_WAIT) return PROCESS_GROUP_WAIT;

  if (status_ == PROCESS_GROUP_FAIL) return PROCESS_GROUP_FAIL;

  if (status_ == PROCESS_GROUP_BUSY) {
    if (ENABLE_FT) {
      simft::Sim_FT_MPI_Status stat;
      int err = simft::Sim_FT_MPI_Wait(&statusRequestFT_, &stat);
      if (err == MPI_ERR_PROC_FAILED) status_ = PROCESS_GROUP_FAIL;
    } else {
      /* todo: actually MPI_TEST is not really necessary here. However, i think it
       * might be a good idea to have this here. MPI_Test might run a system
       * call which enables the OS to switch to the MPI system.
       */
      MPI_Wait(&statusRequest_, MPI_STATUS_IGNORE);
    }
  }
  assert(status_ >= 0);  // check for invalid values
  assert(status_ <= 2);
  return status_;
}

inline const TaskContainer& ProcessGroupManager::getTaskContainer() const { return tasks_; }

inline void ProcessGroupManager::removeTask(Task* t) {
  std::vector<Task*>::iterator position = std::find(tasks_.begin(), tasks_.end(), t);
  if (position != tasks_.end()) {  // == task_.end() means the element was not found
    tasks_.erase(position);
    std::cout << "Removing task" << t->getID() << " " << theMPISystem()->getWorldRank() << " !\n";

  } else {
    std::cout << "Error could not remove task" << t->getID() << " "
              << theMPISystem()->getWorldRank() << " !\n";
  }
}

inline void ProcessGroupManager::setMasterRank(int pGroupRootID) { pgroupRootID_ = pGroupRootID; }
inline int ProcessGroupManager::getMasterRank() { return pgroupRootID_; }

} /* namespace combigrid */

#endif /* PROCESSGROUPMANAGER_HPP_ */
