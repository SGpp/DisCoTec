#ifndef PROCESSGROUPWORKER_HPP_
#define PROCESSGROUPWORKER_HPP_

#include <chrono>
#include "fullgrid/FullGrid.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupSignals.hpp"
#include "mpi/MPISystem.hpp"
#include "mpi_fault_simulator/MPI-FT.h"
#include "task/Task.hpp"
#include "loadmodel/LearningLoadModel.hpp"
#include "vtk/DFGPlotFileWriter.hpp"

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

  // getter for tasks
  inline const TaskContainer& getTasks() const;

  // Perform combination
  void combine();

  /** initializes all subspace sizes in the dsgu according to the dfgs in the
   * global reduce comm*/
  void initCombinedUniDSGVector();

  /** hierarchizes all fgs */
  void hierarchizeFullGrids();

  /** local reduce */
  void addFullGridsToUniformSG();

  /** extracts and dehierarchizes */
  void integrateCombinedSolution();

  /** reduction */
  void reduceUniformSG(RankType globalReduceRankThatCollects = MPI_PROC_NULL);

  /** combine on sparse grid with uniform decomposition of domain */
  void combineUniform();

  void combineLocalAndGlobal(RankType globalReduceRankThatCollects = MPI_PROC_NULL);

  void deleteTasks();

  /** parallel file io of final output grid */
  void parallelEval();

  /** parallel file io of final output grid for uniform decomposition */
  void parallelEvalUniform(std::string filename, LevelVector leval);

  // do task-specific postprocessing
  void doDiagnostics();

  /** calculate the Lp Norm for each individual task */
  std::vector<double> getLpNorms(int p) const;

  /** evaluate norms on (newly created) reference grid */
  std::vector<double> parallelEvalNorm(LevelVector leval) const;

  /** evaluate norms of Task's analytical solution on reference grid */
  std::vector<double> evalAnalyticalOnDFG(LevelVector leval) const;

  /** evaluate norms of combi solution error on reference grid  */
  std::vector<double> evalErrorOnDFG(LevelVector leval) const;

  /** interpolate values on all tasks' component grids */
  std::vector<CombiDataType> interpolateValues(
      const std::vector<std::vector<real>>& interpolationCoordinates) const;

  /** interpolate values on all tasks' component grids and write results to file */
  void writeInterpolatedValuesPerGrid(const std::vector<std::vector<real>>& interpolationCoords,
                                      std::string fileNamePrefix) const;

  /** interpolate values on all tasks' component grids, combine results, and write to a single file
   */
  void writeInterpolatedValues(const std::vector<CombiDataType>& values,
                               const std::string& valuesWriteFilename) const;

  void writeInterpolatedValuesSingleFile(const std::vector<std::vector<real>>& interpolationCoords,
                                         const std::string& filenamePrefix) const;

  /** write the highest and smallest sparse grid coefficient per subspace */
  void writeSparseGridMinMaxCoefficients(std::string fileNamePrefix) const;

  /** write extra SGs to disk (binary w/ MPI-IO) */
  void writeDSGsToDisk(std::string filenamePrefix);

  /** read extra SGs from disk (binary w/ MPI-IO) */
  void readDSGsFromDisk(std::string filenamePrefix, bool alwaysReadFullDSG = false);

  void readDSGsFromDiskAndReduce(std::string filenamePrefixToRead, bool alwaysReadFullDSG = false);

  void setCombiParameters(const CombiParameters& combiParameters);

  /** update combination parameters (for init or after change in FTCT) */
  void updateCombiParameters();

  /** returns the combi parameters */
  inline CombiParameters& getCombiParameters();

  /** performes the sparse grid reduce with the remote system, bcasts solution
   * and updates fgs. */
  void combineThirdLevel();

  void combineThirdLevelFileBasedWrite(std::string filenamePrefixToWrite,
                                       std::string writeCompleteTokenFileName);

  void combineThirdLevelFileBasedReadReduce(std::string filenamePrefixToRead,
                                            std::string startReadingTokenFileName,
                                            bool overwrite = false,
                                            bool keepSparseGridFiles = false);

  void combineThirdLevelFileBased(std::string filenamePrefixToWrite,
                                  std::string writeCompleteTokenFileName,
                                  std::string filenamePrefixToRead,
                                  std::string startReadingTokenFileName);

  /** waits until the third level pg or output group bcasts the combined solution and updates
   * fgs */
  void waitForThirdLevelCombiResult(bool fromOutputGroup = false);

  void setExtraSparseGrid(bool initializeSizes = true);

  /** computes a max reduce on the dsg's subspace sizes with the other systems */
  void reduceSubspaceSizesThirdLevel(bool thirdLevelExtraSparseGrid);

  void reduceSubspaceSizes(const std::string& filenameToRead, bool extraSparseGrid,
                           bool overwrite = false);

  void reduceSubspaceSizesFileBased(std::string filenamePrefixToWrite,
                                    std::string writeCompleteTokenFileName,
                                    std::string filenamePrefixToRead,
                                    std::string startReadingTokenFileName, bool extraSparseGrid);

  /** receives reduced sizes from tl pgroup and updates the dsgs */
  void waitForThirdLevelSizeUpdate();

  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>> & getCombinedUniDSGVector(){
    return combinedUniDSGVector_;
  }

  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>> & getExtraUniDSGVector(){
    return extraUniDSGVector_;
  }

  TaskContainer& getTasks() { return tasks_; }

  template <typename TaskType, typename... TaskArgs>
  void initializeAllTasks(const std::vector<LevelVector>& levels,
                          const std::vector<combigrid::real>& coeffs,
                          const std::vector<size_t>& taskNumbers, TaskArgs&&... args) {
    for (size_t taskIndex = 0; taskIndex < taskNumbers.size(); ++taskIndex) {
      assert(static_cast<DimType>(levels[taskIndex].size()) == this->getCombiParameters().getDim());
      auto task = new TaskType(levels[taskIndex], this->getCombiParameters().getBoundary(),
                               coeffs[taskIndex], std::forward<TaskArgs>(args)...);
      task->setID(taskNumbers[taskIndex]);
      this->initializeTaskAndFaults(task);
    }
  }

  /**
   * @brief store task and run its init function with current combiParameters
   *
   * @param t pointer to a heap-allocated task, the function takes over ownership here
   */
  void initializeTaskAndFaults(Task* t);

  IndexType getCurrentNumberOfCombinations() const { return currentCombi_; }

  /** sets all subspaces in all dsgs to zero and allocates them if necessary */
  void zeroDsgsData();

 private:
  TaskContainer tasks_;  /// task storage

  StatusType status_;  /// current status of process group (wait -> 0; busy -> 1; fail -> 2)

  /**
   * Vector containing all distributed sparse grids
   */
  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>> combinedUniDSGVector_;

  /**
   * Vector containing the third level extra distributed sparse grids
   */
  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>> extraUniDSGVector_;

  CombiParameters combiParameters_;

  bool combiParametersSet_;  /// indicates if combi parameters variable set

  // fault parameters
  real t_fault_;  /// time to fault

  IndexType currentCombi_;  /// current combination; increased after every combination

  std::chrono::high_resolution_clock::time_point
      startTimeIteration_;  /// starting time of process computation

  // std::ofstream betasFile_;

  void receiveAndInitializeTaskAndFaults(bool mayAlreadyExist = true);

  /** deallocates all data elements stored in the dsgs */
  void deleteDsgsData();

  /** the pg writes the dfg of the given task into a vtk file */
  void writeVTKPlotFileOfTask(Task& task);

  /** the pg writes the dfg of all tasks into individual vtk files */
  void writeVTKPlotFilesOfAllTasks();

  void processDuration(const Task& t, const Stats::Event e, unsigned int numProcs);

  /** helper functions for parallelEval and norm calculations*/
  LevelVector receiveLevalAndBroadcast();

  /**
   * @brief copy the sparse grid data into the full grid and dehierarchize
   *
   * @param dfg the distributed full grid to fill
   * @param g the dimension index (in the case that there are multiple different full grids per
   * task)
   */
  void fillDFGFromDSGU(DistributedFullGrid<CombiDataType>& dfg, IndexType g = 0) const;

  void fillDFGFromDSGU(Task* t) const;
};

inline CombiParameters& ProcessGroupWorker::getCombiParameters() {
  assert(combiParametersSet_);

  return combiParameters_;
}

inline const TaskContainer& ProcessGroupWorker::getTasks() const {
  return tasks_;
}

} /* namespace combigrid */

#endif /* PROCESSGROUPWORKER_HPP_ */
