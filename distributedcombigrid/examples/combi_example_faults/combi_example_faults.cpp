	/*
 * combi_example.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: heenemo
 */
#include "mpi.h"
#include <vector>
#include <set>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/serialization/export.hpp>

// compulsory includes for basic functionality
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/LPOptimizationInterpolation.hpp"

// include user specific task. this is the interface to your application
#include "TaskExample.hpp"

#include "HelperFunctions.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskExample)



int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  /* number of process groups and number of processes per group */
  size_t ngroup, nprocs;

  DimType dim;
  LevelVector lmin, lmax, leval;
  IndexVector p;
  size_t ncombi, nsteps;
  combigrid::real dt;
  FaultsInfo faultsInfo;

  /* read parameter file (ctparam) */
  const std::string fileName = "ctparam";
  readParameterFile(fileName, ngroup, nprocs, dim, lmin, lmax,
                    leval, p, ncombi, dt, nsteps, faultsInfo);


  std::vector<bool> boundary(dim,true);

  theMPISystem()->init( ngroup, nprocs );

  // manager code
  if ( theMPISystem()->getWorldRank() == theMPISystem()->getManagerRankWorld() ) {
    /* create an abstraction of the process groups for the manager's view
     * a pgroup is identified by the ID in gcomm
     */
    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(
          std::make_shared< ProcessGroupManager > ( pgroupRootID )
                          );
    }


    /* create load model */
    LoadModel* loadmodel = new LinearLoadModel();

    IndexType checkProcs = 1;
    for (auto k : p)
      checkProcs *= k;
    assert(checkProcs == IndexType(nprocs));

    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createAdaptiveCombischeme();
    combischeme.makeFaultTolerant();
    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

    /* print info of the combination scheme */
    std::cout << "CombiScheme: " << std::endl;
    std::cout << combischeme;

    /* create Tasks */
    TaskContainer tasks;
    std::vector<int> taskIDs;

    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskExample(dim, levels[i], boundary, coeffs[i],
                                loadmodel, dt, nsteps, p, faultsInfo);
      tasks.push_back(t);
      taskIDs.push_back( t->getID() );
    }

    /* create combi parameters */
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs);

    /* create Manager with process groups */
    ProcessManager manager( pgroups, tasks, params );

    /* send combi parameters to workers */
    manager.updateCombiParameters();

    /* distribute task according to load model and start computation for
     * the first time */
    bool success = manager.runfirst();

    for (size_t i = 1; i < ncombi; ++i){

      if ( !success ) {
        std::cout << "failed group detected at combi iteration " << i-1<< std::endl;
//        manager.recover();

        std::vector<int> faultsID;
        manager.getGroupFaultIDs(faultsID);

        /* call optimization code to find new coefficients */
        const std::string prob_name = "interpolation based optimization";
        std::vector<int> redistributeFaultsID, recomputeFaultsID;
        manager.recomputeOptimumCoefficients(prob_name, faultsID,
                                             redistributeFaultsID, recomputeFaultsID);

        for ( auto id : redistributeFaultsID ) {
          TaskExample* tmp = static_cast<TaskExample*>(manager.getTask(id));
          tmp->setStepsTotal(i*nsteps);
        }

        for ( auto id : recomputeFaultsID ) {
          TaskExample* tmp = static_cast<TaskExample*>(manager.getTask(id));
          tmp->setStepsTotal((i-1)*nsteps);
        }
        /* recover communicators*/
        manager.recoverCommunicators();

        /* communicate new combination scheme*/
        manager.updateCombiParameters();

        /* if some tasks have to be recomputed, do so*/
        manager.recompute(recomputeFaultsID);

        /* redistribute failed tasks to living groups */
        manager.redistribute(redistributeFaultsID);
      }

      /* combine solution */
      manager.combine();

      if ( !success ){
        /* restore combischeme to its original state
         * and send new combiParameters to all surviving groups */
        manager.restoreCombischeme();
        manager.updateCombiParameters();
      }

      /* run tasks for next time interval */
      success = manager.runnext();
    }


    /* evaluate solution */
    FullGrid<CombiDataType> fg_eval(dim, leval, boundary);
    manager.gridEval(fg_eval);
    std::string filename( "solution.fg" );
    fg_eval.save( filename );

    /* send exit signal to workers in order to enable a clean program termination */
    manager.exit();
  }

//  // worker code
  else {
    // create abstraction of the process group from the worker's view
    ProcessGroupWorker pgroup;

    // wait for instructions from manager
    SignalType signal = -1;
    while (signal != EXIT)
      signal = pgroup.wait();
  }

  MPI_Finalize();

  return 0;
}
	
