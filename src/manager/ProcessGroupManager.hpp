#ifndef PROCESSGROUPMANAGER_HPP_
#define PROCESSGROUPMANAGER_HPP_

#include <stdexcept>
#include <string>
#include <vector>

#include "combicom/CombiCom.hpp"
#include "fullgrid/FullGrid.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupSignals.hpp"
#include "mpi/MPISystem.hpp"
#include "mpi_fault_simulator/MPI-FT.h"
#include "task/Task.hpp"
#include "third_level/ThirdLevelUtils.hpp"
#include "utils/Types.hpp"

namespace combigrid {
class ProcessGroupManager {
 public:
  ProcessGroupManager() = delete;

  explicit ProcessGroupManager(RankType pgroupRootID);

  ProcessGroupManager(ProcessGroupManager const&) = delete;

  ProcessGroupManager& operator=(ProcessGroupManager const&) = delete;

  bool runfirst(Task* t);

  bool runnext();

  bool initDsgus();

  bool exit();

  /* non-blocking call to retrieve status of process group */
  inline StatusType getStatus();

  inline void setStatus(StatusType status);

  /* blocks until process group finished computation */
  inline StatusType waitStatus();

  inline const TaskContainer& getTaskContainer() const;

  inline void removeTask(Task* t);

  bool combine();

  // third Level stuff
  bool combineThirdLevel(const ThirdLevelUtils& thirdLevel, CombiParameters& params,
                         bool isSendingFirst);

  bool combineThirdLevelFileBased(std::string filenamePrefixToWrite,
                                  std::string writeCompleteTokenFileName,
                                  std::string filenamePrefixToRead,
                                  std::string startReadingTokenFileName);

  bool combineThirdLevelFileBasedWrite(std::string filenamePrefixToWrite,
                                       std::string writeCompleteTokenFileName);

  bool combineThirdLevelFileBasedReadReduce(std::string filenamePrefixToRead,
                                            std::string startReadingTokenFileName);

  bool pretendCombineThirdLevelForWorkers(CombiParameters& params);

  bool combineSystemWide();

  bool updateCombiParameters(CombiParameters& params);

  /* Check if group fault occured at this combination step using the fault simulator */
  bool isGroupFault();

  bool addTask(Task*);
  bool refreshTask(Task*);

  // resets tasks only on workers not in group manager
  bool resetTasksWorker();

  bool recompute(Task*);

  bool recoverCommunicators();

  bool parallelEval(const LevelVector& leval, std::string& filename);

  void writeSparseGridMinMaxCoefficients(const std::string& filename);

  void doDiagnostics(size_t taskID);

  void getLpNorms(int p, std::map<size_t, double>& norms);

  std::vector<double> evalAnalyticalOnDFG(const LevelVector& leval);

  std::vector<double> evalErrorOnDFG(const LevelVector& leval);

  void interpolateValues(const std::vector<real>& interpolationCoordsSerial,
                         std::vector<CombiDataType>& values, MPI_Request* request = nullptr,
                         std::string filenamePrefix = "");

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
  Task *rescheduleRemoveTask(const LevelVector& lvlVec);


  bool hasTask(size_t taskID){
    auto foundIt = std::find_if(tasks_.begin(), tasks_.end(),
                        [taskID](Task* t){
                          return ((t->getID()) == taskID);
                        });
    return foundIt != tasks_.end();
  }

  bool writeCombigridsToVTKPlotFile();

  bool writeDSGsToDisk(std::string filenamePrefix);

  bool readDSGsFromDisk(std::string filenamePrefix);

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
