/*
 * ProcessGroupWorker.hpp
 *
 *  Created on: Jun 24, 2014
 *      Author: heenemo
 */

#ifndef PROCESSGROUPWORKER_HPP_
#define PROCESSGROUPWORKER_HPP_

#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"
#include <chrono>

namespace combigrid {

class ProcessGroupWorker {
 public:
  explicit ProcessGroupWorker();

  ProcessGroupWorker( ProcessGroupWorker const & ) = delete;

  ProcessGroupWorker& operator=( ProcessGroupWorker const & ) = delete;

  ~ProcessGroupWorker();

  // wait for command from manager
  SignalType wait();

  // send ready signal to manager
  void ready();

  //decides if current Task needs to be killed
  void decideToKill();

  // todo: maybe only needed for gene?
  inline Task* getCurrentTask();

  void combine();

  void combineUniform();

  void combineFG();

 // void gridEval();

  void parallelEval();

  void parallelEvalUniform();

  void updateCombiParameters();

  inline CombiParameters& getCombiParameters();

  void setCombinedSolutionUniform( Task* t );

 private:
  TaskContainer tasks_; // task storage

  Task* currentTask_;

  StatusType status_;

  FullGrid<complex>* combinedFG_;

  std::vector<DistributedSparseGridUniform<CombiDataType>*> combinedUniDSGVector_;

  bool combinedFGexists_;

  CombiParameters combiParameters_;

  bool combiParametersSet_;

  //fault parameters
  real t_fault_; //time to fault

  int numGrids_; //number of grids per task

  std::chrono::high_resolution_clock::time_point  startTimeIteration_; //starting time of process computation

  //std::ofstream betasFile_;

};


inline Task* ProcessGroupWorker::getCurrentTask() {
  return currentTask_;
}


inline CombiParameters& ProcessGroupWorker::getCombiParameters(){
  assert(combiParametersSet_);

  return combiParameters_;
}

} /* namespace combigrid */

#endif /* PROCESSGROUPWORKER_HPP_ */
