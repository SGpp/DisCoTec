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
#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"

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

  // todo: maybe only needed for gene?
  inline Task* getCurrentTask();

  void combine();

  void combineUniform();

  void combineLocalAndGlobal();

  template<typename FG_ELEMENT>
  void combineUniformThirdLevelSendFirst();

  template<typename FG_ELEMENT>
  void combineUniformThirdLevelRecvFirst();

  template<typename FG_ELEMENT>
  void integrateCommonSS();

  void combineFG();

  void gridEval();

  void parallelEval();

  void parallelEvalUniform();

  void updateCombiParameters();

  inline CombiParameters& getCombiParameters();

 private:
  TaskContainer tasks_; // task storage

  Task* currentTask_;

  StatusType status_;

  FullGrid<complex>* combinedFG_;

  DistributedSparseGridUniform<CombiDataType>* combinedUniDSG_;

  ThirdLevelUtils* thirdLevel_;

  bool combinedFGexists_;

  CombiParameters combiParameters_;

  bool combiParametersSet_;

  void setCombinedSolutionUniform( Task* t );
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
