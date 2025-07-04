#ifndef PROCESSGROUPWORKER_HPP_
#define PROCESSGROUPWORKER_HPP_

#include "CombiParameters.hpp"
#include "ProcessGroupSignals.hpp"
#include "SparseGridWorker.hpp"
#include "TaskWorker.hpp"
#include "../mpi/MPISystem.hpp"
#include "../Task.hpp"

namespace combigrid {

/**
 * @class ProcessGroupWorker
 *
 * @brief The ProcessGroupWorker comprises a TaskWorker and a SparseGridWorker;
 * it is responsible for executing the time stepping in the combination technique as well as various
 * I-O operations.
 *
 * When using a manager-worker scheme, the ProcessGroupWorker is instantiated on each worker rank;
 * else, it is instantiated on each rank.
 */
class ProcessGroupWorker {
 public:
  explicit ProcessGroupWorker();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  ProcessGroupWorker(ProcessGroupWorker const&) = delete;
  ProcessGroupWorker& operator=(ProcessGroupWorker const&) = delete;
  ProcessGroupWorker(ProcessGroupWorker&&) = delete;
  ProcessGroupWorker& operator=(ProcessGroupWorker&&) = delete;
  ~ProcessGroupWorker() = default;
#endif  // DOXYGEN_SHOULD_SKIP_THIS

  /** @brief wait for command from manager */
  SignalType wait();

  /**
   * @brief call run on each Task once
   */
  void runAllTasks();

  /**
   * @brief release all tasks
   */
  void exit();

  /**
   * @brief get the vector of tasks
   */
  inline const std::vector<std::unique_ptr<Task>>& getTasks() const;

  /**@brief  initializes all subspace sizes in the distributed sparse grid data structure */
  void initCombinedDSGVector();

  /**
   * @brief extracts data from sparse to component grids and dehierarchizes the component grids
   */
  void updateFullFromCombinedSparseGrids();

  /**
   * @brief dehierarchize all component grids
   */
  void dehierarchizeAllTasks();

  /**
   * @brief do a whole combination, assuming a single-system setup
   *
   * @param collectMinMaxCoefficients if true, the min and max coefficients per subspace will be
   * collected as part of the combination
   */
  void combineAtOnce(bool collectMinMaxCoefficients = false);

  /**
   * @brief perform the first part of a system-wide combination,  namely hierarchization
   * and reduction (but not dehierarchization)
   *
   * based on file-exchange mechanism (w/o third level manager)
   */
  void combineSystemWide();

  /**
   * @brief perform the first part of a system-wide combination, but reducing to output
   * ranks only instead of all ranks, which then write data to file.
   *
   * based on file-exchange mechanism (w/o third level manager)
   */
  void combineSystemWideAndWrite(const std::string& writeSparseGridFile,
                                 const std::string& writeSparseGridFileToken);

  /**
   * @brief interpolate the current solution from all component grids on this group
   *
   * receives the target level and file name from the ProcessGroupManager
   */
  void parallelEval();

  /**
   * @brief called by parallelEval to interpolate at resolution level \p leval and write the results
   * to a binary file readable with Paraview
   */
  void parallelEvalUniform(const std::string& filename, const LevelVector& leval) const;

  /** @brief do task-specific diagnostics / postprocessing */
  void doDiagnostics();
  void doDiagnostics(size_t taskID);

  /** @brief calculate the Lp Norm for each individual task */
  std::vector<double> getLpNorms(int p) const;

  /** @brief interpolate values on all tasks' component grids
   *
   * @param interpolationCoordinates the coordinates to interpolate at
   * @return the interpolated values, one vector of values per task
   */
  std::vector<CombiDataType> interpolateValues(
      const std::vector<std::vector<real>>& interpolationCoordinates) const;

  /**
   * @brief interpolate values on all tasks' component grids and write results to file
   *
   * one file per component grid, values are not combined
   */
  void writeInterpolatedValuesPerGrid(const std::vector<std::vector<real>>& interpolationCoords,
                                      const std::string& fileNamePrefix) const;

  /**
   * @brief interpolate values on all tasks' component grids and write results to file
   *
   * one file for all component grids; values are combined with the combination formula
   */
  void writeInterpolatedValuesSingleFile(const std::vector<std::vector<real>>& interpolationCoords,
                                         const std::string& filenamePrefix) const;

  /** @brief write the highest and smallest sparse grid coefficient per subspace */
  void writeSparseGridMinMaxCoefficients(const std::string& fileNamePrefix) const;

  /**
   * @brief write extra sparse grids to disk (binary w/ MPI-IO)
   *
   * this will only have an effect on the output ranks, as they are the only ones that have an extra
   * sparse grid.
   * used for widely-distributed simulations with files (no third level manager)
   */
  int writeDSGsToDisk(const std::string& filenamePrefix,
                      const std::string& writeCompleteTokenFileName);

  /**
   * @brief read extra SGs from disk (binary w/ MPI-IO)
   *
   * this will only have an effect on the output ranks, as they are the only ones that have an extra
   * sparse grid
   *
   * used for widely-distributed simulations with files (no third level manager)
   */
  int readDSGsFromDisk(const std::string& filenamePrefix, bool alwaysReadFullDSG = false);

  /**
   * @brief set the CombiParameters member
   *
   * usually after signal from ProcessGroupManager (when in manager-worker setup) or directly (when
   * in worker-only setup)
   */
  void setCombiParameters(CombiParameters&& combiParameters);

  /** @brief get a reference to combi parameters */
  inline CombiParameters& getCombiParameters();

  /**
   * @brief perform a whole widely-distributed combination
   *
   * based on TCP/socket setup with third level manager.
   */
  void combineThirdLevel();

  /**
   * @brief write the contents of the extra sparse grids to disk
   *
   * based on file-exchange mechanism (w/o third level manager);
   * should only be called from the ranks in the output group.
   *
   * @param filenamePrefixToWrite the prefixes of the file to write, will have the species and (if
   * partitioned output) the partition number appended
   * @param writeCompleteTokenFileName either the name of the token file to write, or the prefix (in
   * case the outputComm is set, for partitioned file output), in which case the partition number
   * will be appended
   */
  int combineThirdLevelFileBasedWrite(const std::string& filenamePrefixToWrite,
                                      const std::string& writeCompleteTokenFileName);

  /**
   * @brief wait for a token file to appear.
   *
   * one rank polls for the token file to appear (every some microseconds), and signals to the other
   * output ranks once it does.
   *
   * based on file-exchange mechanism (w/o third level manager);
   * should only be called from the ranks in the output group.
   */
  void waitForTokenFile(const std::string& startReadingTokenFileName) const;

  /**
   * @brief read the contents of the extra sparse grids from disk and reduce them
   *
   * based on file-exchange mechanism (w/o third level manager);
   * should only be called from the ranks in the output group.
   *
   * @param filenamePrefixToRead the prefixes of the file to read, will have the species and (if
   * partitioned output) the partition number appended
   * @param startReadingTokenFileName either the name of the token file to read, or the prefix (in
   * case the outputComm is set, for partitioned file output), in which case the partition number
   * will be appended
   * @param overwriteInMemory whether to overwrite the values in the extra sparse grid (only true
   * for testing)
   */
  int readReduce(const std::vector<std::string>& filenamePrefixesToRead,
                 const std::vector<std::string>& startReadingTokenFileNames,
                 bool overwriteInMemory);

  /**
   * @brief perform the second part of the widely-distributed combination, namely reading the
   * reduced sparse grids from disk and combining them with the local solution
   *
   * based on file-exchange mechanism (w/o third level manager); this function should only be called
   * from the ranks in the output group, cf. combineReadDistributeSystemWide
   *
   * @param filenamePrefixesToRead the prefixes of the files to read
   * @param startReadingTokenFileNames the names or prefixes of the token files to wait for
   * @param overwrite whether to overwrite the values in the extra sparse grid (only true for
   * testing)
   * @param keepSparseGridFiles whether to keep the sparse grid files after reading
   */
  void combineThirdLevelFileBasedReadReduce(
      const std::vector<std::string>& filenamePrefixesToRead,
      const std::vector<std::string>& startReadingTokenFileNames, bool overwrite = false,
      bool keepSparseGridFiles = false);

  /**
   * @brief perform the second part of the widely-distributed combination, namely reading the
   * reduced sparse grids from disk and combining them with the local solution
   *
   * based on file-exchange mechanism (w/o third level manager); this function should be called from
   * all worker ranks
   *
   * @param filenamePrefixToRead the prefixes of the files to read
   * @param startReadingTokenFileName the names or prefixes of the token files to wait for
   * @param overwrite whether to overwrite the values in the extra sparse grid (only true for
   * testing)
   * @param keepSparseGridFiles whether to keep the sparse grid files after reading   *
   */
  void combineReadDistributeSystemWide(const std::vector<std::string>& filenamePrefixesToRead,
                                       const std::vector<std::string>& startReadingTokenFileNames,
                                       bool overwrite = false, bool keepSparseGridFiles = false);

  /**
   * @brief perform a whole widely-distributed combination
   *
   * based on file-exchange mechanism (w/o third level manager); this function should be called from
   * all worker ranks.
   *
   * equivalent to calling combineSystemWideAndWrite and combineReadDistributeSystemWide in
   * succession
   */
  void combineThirdLevelFileBased(const std::string& filenamePrefixToWrite,
                                  const std::string& writeCompleteTokenFileName,
                                  const std::vector<std::string>& filenamePrefixToRead,
                                  const std::vector<std::string>& startReadingTokenFileName);

  /**
   * @brief waits until the third level pg or output group broadcasts the combined solution, then
   * updates the component grids
   *
   * used for both kinds of widely-distributed simulations (with and without third level manager)
   *
   * @param fromOutputGroup whether the signal comes from the output group or from the third level
   * group
   */
  void waitForThirdLevelCombiResult(bool fromOutputGroup = false);

  /** computes a max reduce on the dsg's subspace sizes with the other systems */
  void reduceSubspaceSizesThirdLevel(bool thirdLevelExtraSparseGrid);

  /**
   * @brief reduce the subspaces sizes in the extra sparse grid from the sizes given in a file
   *
   * @param filenameToRead the file to read the sizes from
   * @param overwrite whether to overwrite the sizes in the extra sparse grid (only true for
   * testing)
   */
  int reduceExtraSubspaceSizes(const std::vector<std::string>& filenameToRead,
                               bool overwrite = false);

  /**
   * @brief reduce the subspaces sizes in the extra sparse grid with another system
   *
   * @param filenamePrefixToWrite the prefix of the file to write the own sizes to
   * @param writeCompleteTokenFileName the name of the token file to write
   * @param filenamePrefixToRead the prefixes of the files to read the sizes from
   * @param startReadingTokenFileName the names or prefixes of the token files to wait for
   */
  int reduceExtraSubspaceSizesFileBased(const std::string& filenamePrefixToWrite,
                                        const std::string& writeCompleteTokenFileName,
                                        const std::vector<std::string>& filenamePrefixToRead,
                                        const std::vector<std::string>& startReadingTokenFileName);

  /**
   * @brief resize the sparse grid data structures according to sizes received from the third level
   * process group
   *
   * used for widely-distributed simulations based on TCP/socket setup with third level manager
   */
  void waitForThirdLevelSizeUpdate();

  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
  getCombinedDSGVector() {
    return this->getSparseGridWorker().getCombinedUniDSGVector();
  }

  /**
   * @brief get the extra sparse grid data structures
   */
  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& getExtraDSGVector() {
    return this->getSparseGridWorker().getExtraUniDSGVector();
  }

  /**
   * @brief initialize all tasks at once
   *
   * used in worker-only setup
   */
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

  /**
   * @brief get the current number of combinations
   */
  IndexType getCurrentNumberOfCombinations() const { return currentCombi_; }

  /**
   * @brief sets all subspaces in all dsgs to zero; allocates them if necessary
   */
  void zeroDsgsData();

  /**
   * @brief get a reference to the task worker member
   */
  const TaskWorker& getTaskWorker() const { return taskWorker_; }

  /**
   * @brief get a reference to the sparse grid worker member
   */
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

  /** update combination parameters (for init or after change in FTCT) */
  void updateCombiParameters();

  void setExtraSparseGrid(bool initializeSizes = true);
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
