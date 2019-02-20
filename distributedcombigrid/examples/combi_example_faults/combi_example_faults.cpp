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
#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"
#include "sgpp/distributedcombigrid/fault_tolerance/FaultCriterion.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/StaticFaults.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/WeibullFaults.hpp"

// include user specific task. this is the interface to your application
#include "TaskExample.hpp"

#include "HelperFunctions.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskExample)
BOOST_CLASS_EXPORT(StaticFaults)
BOOST_CLASS_EXPORT(WeibullFaults)

BOOST_CLASS_EXPORT(FaultCriterion)


int main(int argc, char** argv) {
  simft::Sim_FT_MPI_Init(&argc, &argv);

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
      //create FaultCriterion
      FaultCriterion *faultCrit;
      //create fault criterion
      if(faultsInfo.numFaults_ < 0){ //use random distributed faults
        //if numFaults is smallerthan 0 we use the absolute value
        //as lambda value for the weibull distribution
        faultCrit = new WeibullFaults(0.7, abs(faultsInfo.numFaults_), ncombi, true);
      }
      else{ //use predefined static number and timing of faults
        //if numFaults = 0 there are no faults
        faultCrit = new StaticFaults(faultsInfo);
      }
      Task* t = new TaskExample(dim, levels[i], boundary, coeffs[i],
                                loadmodel, dt, nsteps, p, faultCrit);
      tasks.push_back(t);
      taskIDs.push_back( t->getID() );
    }

    /* create combi parameters */
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi, 1);
    params.setParallelization(p);

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

        //vector with pointers to managers of failed groups
        std::vector< ProcessGroupManagerID> groupFaults;
        manager.getGroupFaultIDs(faultsID, groupFaults);

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
        bool failedRecovery = manager.recoverCommunicators(groupFaults);


        if(failedRecovery){
          //if the process groups could not be restored distribute tasks to other groups
          std::cout << "Redistribute groups \n";
          manager.redistribute(redistributeFaultsID);
        }
        else{
          //if process groups could be restored reinitialize restored process group (keep the original tasks)
          std::cout << "Reinitializing groups \n";
          manager.reInitializeGroup(groupFaults,recomputeFaultsID);
        }
        /* if some tasks have to be recomputed, do so
         * allowing recomputation reduces the overhead that would be needed
         * for finding a scheme that avoids all failed tasks*/
        if(!recomputeFaultsID.empty()){
          std::cout << "sending tasks for recompute \n";
          manager.recompute(recomputeFaultsID,failedRecovery,groupFaults); //toDO handle faults in recompute
        }
        std::cout << "updateing Combination Parameters \n";
        //needs to be after reInitialization!
        /* communicate new combination scheme*/
        manager.updateCombiParameters();

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


    std::string filename("out/solution_" + std::to_string(ncombi) + ".dat" );
    manager.parallelEval( leval, filename, 0 );

    /* send exit signal to workers in order to enable a clean program termination */
    manager.exit();
  }

// worker code
  else {
    // create abstraction of the process group from the worker's view
    ProcessGroupWorker pgroup;

    // wait for instructions from manager
    SignalType signal = -1;
    while (signal != EXIT)
      signal = pgroup.wait();
  }

  if( ENABLE_FT ){
    WORLD_MANAGER_EXCLUSIVE_SECTION{
      std::cout << "Program finished successfully" << std::endl;
      std::cout << "To avoid problems with hanging killed processes, we exit with "
                << "MPI_Abort()" << std::endl;
      MPI_Abort( MPI_COMM_WORLD, 0 );
    }
  }

  simft::Sim_FT_MPI_Finalize();

  return 0;
}
	
