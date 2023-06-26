#ifndef PROCESSGROUPWORKER_HPP_
#define PROCESSGROUPWORKER_HPP_

#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupSignals.hpp"
#include "manager/SparseGridWorker.hpp"
#include "manager/TaskWorker.hpp"
#include "mpi/MPISystem.hpp"
#include "task/Task.hpp"

namespace combigrid {

class ProcessGroupWorker {
 public:
  explicit ProcessGroupWorker();

  ProcessGroupWorker(ProcessGroupWorker const&) = delete;
  ProcessGroupWorker& operator=(ProcessGroupWorker const&) = delete;

  ~ProcessGroupWorker();

  /** wait for command from manager */
  SignalType wait();

  void runAllTasks();

  void exit();

  inline const std::vector<std::unique_ptr<Task>>& getTasks() const;

  /** initializes all subspace sizes in the dsgu according to the dfgs in the
   * global reduce comm*/
  void initCombinedDSGVector();

  /** extracts and dehierarchizes */
  void updateFullFromCombinedSparseGrids();

  /** combine on sparse grid with uniform decomposition of domain */
  void combineUniform();

  void combineLocalAndGlobal(RankType globalReduceRankThatCollects = MPI_PROC_NULL);

  /** parallel file io of final output grid */
  void parallelEval();

  /** parallel file io of final output grid for uniform decomposition */
  void parallelEvalUniform(const std::string& filename, const LevelVector& leval) const;

  // do task-specific postprocessing
  void doDiagnostics();

  /** calculate the Lp Norm for each individual task */
  std::vector<double> getLpNorms(int p) const;

  /** interpolate values on all tasks' component grids */
  std::vector<CombiDataType> interpolateValues(
      const std::vector<std::vector<real>>& interpolationCoordinates) const;

  /** interpolate values on all tasks' component grids and write results to file */
  void writeInterpolatedValuesPerGrid(const std::vector<std::vector<real>>& interpolationCoords,
                                      const std::string& fileNamePrefix) const;

  void writeInterpolatedValuesSingleFile(const std::vector<std::vector<real>>& interpolationCoords,
                                         const std::string& filenamePrefix) const;

  /** write the highest and smallest sparse grid coefficient per subspace */
  void writeSparseGridMinMaxCoefficients(const std::string& fileNamePrefix) const;

  /** write extra SGs to disk (binary w/ MPI-IO) */
  int writeDSGsToDisk(const std::string& filenamePrefix);

  /** read extra SGs from disk (binary w/ MPI-IO) */
  int readDSGsFromDisk(const std::string& filenamePrefix, bool alwaysReadFullDSG = false);

  void setCombiParameters(CombiParameters&& combiParameters);

  /** update combination parameters (for init or after change in FTCT) */
  void updateCombiParameters();

  /** returns the combi parameters */
  inline CombiParameters& getCombiParameters();

  /** performes the sparse grid reduce with the remote system, bcasts solution
   * and updates fgs. */
  void combineThirdLevel();

  int combineThirdLevelFileBasedWrite(const std::string& filenamePrefixToWrite,
                                      const std::string& writeCompleteTokenFileName);

  void combineThirdLevelFileBasedReadReduce(const std::string& filenamePrefixToRead,
                                            const std::string& startReadingTokenFileName,
                                            bool overwrite = false,
                                            bool keepSparseGridFiles = false);

  void combineThirdLevelFileBased(const std::string& filenamePrefixToWrite,
                                  const std::string& writeCompleteTokenFileName,
                                  const std::string& filenamePrefixToRead,
                                  const std::string& startReadingTokenFileName);

  /** waits until the third level pg or output group bcasts the combined solution and updates
   * fgs */
  void waitForThirdLevelCombiResult(bool fromOutputGroup = false);

  void setExtraSparseGrid(bool initializeSizes = true);

  /** computes a max reduce on the dsg's subspace sizes with the other systems */
  void reduceSubspaceSizesThirdLevel(bool thirdLevelExtraSparseGrid);

  int reduceSubspaceSizes(const std::string& filenameToRead, bool extraSparseGrid,
                          bool overwrite = false);

  int reduceSubspaceSizesFileBased(const std::string& filenamePrefixToWrite,
                                   const std::string& writeCompleteTokenFileName,
                                   const std::string& filenamePrefixToRead,
                                   const std::string& startReadingTokenFileName,
                                   bool extraSparseGrid);

  /** receives reduced sizes from tl pgroup and updates the dsgs */
  void waitForThirdLevelSizeUpdate();

  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
  getCombinedDSGVector() {
    return this->getSparseGridWorker().getCombinedUniDSGVector();
  }

  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& getExtraDSGVector() {
    return this->getSparseGridWorker().getExtraUniDSGVector();
  }

  template <typename TaskType, typename... TaskArgs>
  void initializeAllTasks(const std::vector<LevelVector>& levels,
                          const std::vector<combigrid::real>& coeffs,
                          const std::vector<size_t>& taskNumbers, TaskArgs&&... args) {
    for (size_t taskIndex = 0; taskIndex < taskNumbers.size(); ++taskIndex) {
      assert(static_cast<DimType>(levels[taskIndex].size()) == this->getCombiParameters().getDim());
      auto task = std::unique_ptr<Task>(
          new TaskType(levels[taskIndex], this->getCombiParameters().getBoundary(),
                       coeffs[taskIndex], std::forward<TaskArgs>(args)...));
      task->setID(taskNumbers[taskIndex]);
      this->initializeTask(std::move(task));
    }
  }

  /**
   * @brief store task and run its init function with current combiParameters
   *
   * @param t pointer to a heap-allocated task, the function takes over ownership here
   */
  void initializeTask(std::unique_ptr<Task> t);

  IndexType getCurrentNumberOfCombinations() const { return currentCombi_; }

  /** sets all subspaces in all dsgs to zero and allocates them if necessary */
  void zeroDsgsData();

  const TaskWorker& getTaskWorker() const { return taskWorker_; }

  const SparseGridWorker& getSparseGridWorker() const { return sgWorker_; }

 private:
  TaskWorker taskWorker_{};  // worker that has tasks / full grids

  SparseGridWorker sgWorker_{taskWorker_};  // worker that has sparse grids

  StatusType status_;  /// current status of process group (wait -> 0; busy -> 1; fail -> 2)

  CombiParameters combiParameters_;

  bool combiParametersSet_;  /// indicates if combi parameters variable set

  IndexType currentCombi_;  /// current combination; increased after every combination

  TaskWorker& getTaskWorker() { return taskWorker_; }

  SparseGridWorker& getSparseGridWorker() { return sgWorker_; }

  void receiveAndInitializeTask();
};

inline CombiParameters& ProcessGroupWorker::getCombiParameters() {
  assert(combiParametersSet_);

  return combiParameters_;
}

inline const std::vector<std::unique_ptr<Task>>& ProcessGroupWorker::getTasks() const {
  return this->getTaskWorker().getTasks();
}

} /* namespace combigrid */

#endif /* PROCESSGROUPWORKER_HPP_ */
