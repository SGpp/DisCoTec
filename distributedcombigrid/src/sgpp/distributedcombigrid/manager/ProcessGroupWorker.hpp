/*
 * ProcessGroupWorker.hpp
 *
 *  Created on: Jun 24, 2014
 *      Author: heenemo
 */

#ifndef PROCESSGROUPWORKER_HPP_
#define PROCESSGROUPWORKER_HPP_

#include <chrono>
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"

namespace combigrid {

class ProcessGroupWorker {
 public:
  explicit ProcessGroupWorker();

  ProcessGroupWorker(ProcessGroupWorker const&) = delete;

  ProcessGroupWorker& operator=(ProcessGroupWorker const&) = delete;

  ~ProcessGroupWorker();

  // wait for command from manager
  SignalType wait();

  // send ready signal to manager
  void ready();

  // decides if current Task needs to be killed
  void decideToKill();

  // todo: maybe only needed for gene?
  inline Task* getCurrentTask();

  // Perform combination
  void combine();

  // combination helpers
  void initCombinedUniDSGVector();
  void hierarchizeUniformSG();
  void dehierarchizeUniformSG();

  // reduction
  void reduceUniformSG();

  // combine on sparse grid with uniform decomposition of domain
  void combineUniform();

  // outdated!
  void combineFG();

  void gridEval();

  // parallel file io of final output grid
  void parallelEval();

  // parallel file io of final output grid for uniform decomposition
  void parallelEvalUniform();

  // update combination parameters (for init or after change in FTCT)
  void updateCombiParameters();

  // returns the combi parameters
  inline CombiParameters& getCombiParameters();

  // initializes the component grid from the sparse grid; used to reinitialize tasks after fault
  void setCombinedSolutionUniform(Task* t);

  void sendSparseGridToManager();

 private:
  TaskContainer tasks_;  // task storage

  Task* currentTask_;  // task that is currently processed

  StatusType status_;  // current status of process group (wait -> 0; busy -> 1; fail -> 2)

  FullGrid<complex>* combinedFG_;

  /**
   * Vector containing all distributed sparse grids
   */
  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>> combinedUniDSGVector_;

  bool combinedFGexists_;

  CombiParameters combiParameters_;

  bool combiParametersSet_;  // indicates if combi parameters variable set

  // fault parameters
  real t_fault_;  // time to fault

  IndexType currentCombi_;  // current combination; increased after every combination

  std::chrono::high_resolution_clock::time_point
      startTimeIteration_;  // starting time of process computation

  // std::ofstream betasFile_;

  void initializeTaskAndFaults(bool mayAlreadyExist = true);

  void processDuration(const Task& t, const Stats::Event e, size_t numProcs);
};

inline Task* ProcessGroupWorker::getCurrentTask() { return currentTask_; }

inline CombiParameters& ProcessGroupWorker::getCombiParameters() {
  assert(combiParametersSet_);

  return combiParameters_;
}

} /* namespace combigrid */

#endif /* PROCESSGROUPWORKER_HPP_ */
