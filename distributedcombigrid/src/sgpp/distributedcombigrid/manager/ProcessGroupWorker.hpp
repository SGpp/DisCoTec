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
#include "sgpp/distributedcombigrid/fault_tolerance/FTUtils.hpp"

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>

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

  void combineFG();

  void gridEval();

  void updateCombiParameters();

  void searchSDC();

  /* Computes the difference between all pairs of combination solutions
   * (or only between neighboring ones if onlyNearestNeighbors = true)
   * according to the paper on SDC detection. If the difference is large,
   * a soft fault might have occurred. */
  void comparePairsDistributed( int numNearestNeighbors, std::vector<int> &levelsSDC );

  void comparePairsSerial( int numNearestNeighbors, std::vector<int> &levelsSDC );

  int compareValues();

  void computeLMSResiduals( gsl_multifit_robust_workspace* regressionWsp, gsl_vector* r_lms );

  /* Generates a list of pairs of tasks, so that for each task
   * that a worker has, we find its K nearest neighbors. The distance
   * between two tasks is the l1 norm of the difference of their level vectors:
   * distance(task_s, task_t) = |s - t|_1
   * */
  void generatePairs( int numNearestNeighbors, std::vector<std::vector<Task*>> &allPairs);

  void filterSDCGSL( std::vector<int> &levelsSDC );

  void filterSDCPython( std::vector<int> &levelsSDC );

  void detectOutliers( double* r_lms ,std::vector<int> &levelsSDC );

 private:
  TaskContainer tasks_; // task storage

  Task* currentTask_;

  StatusType status_;

  FullGrid<complex>* combinedFG_;

  DistributedSparseGridUniform<CombiDataType>* combinedUniDSG_;

  bool combinedFGexists_;

  CombiParameters combiParameters_;

  bool combiParametersSet_;

  MPI_File betasFile_;

  std::map <std::pair<LevelVector,LevelVector>, CombiDataType> betas_;

  void setCombinedSolutionUniform( Task* t );
};

inline Task* ProcessGroupWorker::getCurrentTask() {
  return currentTask_;
}

} /* namespace combigrid */

#endif /* PROCESSGROUPWORKER_HPP_ */
