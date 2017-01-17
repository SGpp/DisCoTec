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
#include "sgpp/distributedcombigrid/fault_tolerance/LPOptimizationInterpolation.hpp"

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
  void compareSolutions( int numNearestNeighbors, std::vector<int> &levelsSDC, SDCMethodType method );

  /* Obtain standardized residuals from robust residuals. If a standardized residual is larger than 2.5,
   * it is considered an outlier.
   * */
  void computeStandardizedResiduals( gsl_multifit_robust_workspace* regressionWsp, gsl_vector* r_stud, gsl_vector* r_lms );

  /* Generates a list of pairs of tasks, so that for each task
   * that a worker has, we find its K nearest neighbors. The distance
   * between two tasks is the l1 norm of the difference of their level vectors:
   * distance(task_s, task_t) = |s - t|_1
   * */
  void generatePairs( int numNearestNeighbors, std::vector<std::vector<Task*>> &allPairs);

  /* Robust fit of beta values using an error expansion model
   * */
  void robustRegressionPairs( std::vector<int> &levelsSDC );

  /* Robust fit of function values using a constant model
   * */
  void robustRegressionValues( std::vector<int> &levelsSDC );

  /* Determine from standardized residuals if a solution has been affected by SDC
   * */
  void detectOutliers( double* residuals, std::vector<int> &levelsSDC, double eps, SDCMethodType method, double y = 0.0 );

  /* Check if a solution marked as outlier is actually only a false positive. We combine the values
   * with the outlier (using the classical combination coefficients) and without it (adjusting the
   * combination coefficients according to the FTCT) and compare the two combined values. If they're very
   * similar, they suspect values is most likely a false positive, or the error is too small as to be
   * safely ignored.
   * */
  void removeFalsePositives( std::vector<int>& faultsID, double u_robust );

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

  std::map <std::pair<LevelVector, LevelVector>, CombiDataType> betas_;

  std::map <LevelVector, CombiDataType> subspaceValues_;

  void setCombinedSolutionUniform( Task* t );
};

inline Task* ProcessGroupWorker::getCurrentTask() {
  return currentTask_;
}

} /* namespace combigrid */

#endif /* PROCESSGROUPWORKER_HPP_ */
