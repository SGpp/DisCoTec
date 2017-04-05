/*
 * combi_example.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: heenemo
 */
#include <mpi.h>
#include <vector>
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

// include user specific task. this is the interface to your application
#include "TaskExample.hpp"

#include "HelperFunctions.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskExample)

void solveProblem(ProcessManager& manager, LevelVector& leval,
                  std::vector<bool>& boundary, double time_step,
                  int numSteps, std::vector<int> time_steps_combi,
                  bool plot, int sdcMethod )
{
  std::ofstream valueFile("out/values.dat");

  WORLD_MANAGER_EXCLUSIVE_SECTION{ theStatsContainer()->setTimerStart("ct"); }

  int combistep = 0;
  std::vector<int> reinitFaultsID, recomputeFaultsID, faultsID;
  for (int step = 0; step < numSteps; ++step) {

    if( step == 0 ){
      manager.runfirst();
    } else {
      manager.runnext();
    }

    if ( step == time_steps_combi[combistep]){

      bool success = true;
      success = manager.searchSDC( sdcMethod );
      combistep++;
      if ( !success ) {
        std::cout << "SDC detected at combi iteration " << step << std::endl;
        manager.getSDCFaultIDs( faultsID );
        /* call optimization code to find new coefficients */
        const std::string prob_name = "interpolation based optimization";

        manager.recomputeOptimumCoefficients(prob_name, faultsID,
            reinitFaultsID, recomputeFaultsID);

        /* communicate new combination scheme*/
        manager.updateCombiParameters();

        /* if some tasks have to be recomputed, do so*/
        for ( auto id : recomputeFaultsID ) {
          TaskExample* tmp = static_cast<TaskExample*>(manager.getTask(id));
          tmp->setStepsTotal((combistep-1)*numSteps);
        }
        manager.recompute(recomputeFaultsID);
      }


      /* combine solution */
      std::cout<<"Combining.."<<std::endl;
      manager.combine();

      // evaluate
      if( plot ){
        FullGrid<CombiDataType> fg( leval.size(), leval, boundary);
        std::cout << "eval solution for plotting" << std::endl;
        manager.gridEval(fg);
        // output for function values in gnuplot (only in 2D)
        writeSolutionToFile(valueFile, fg);
      }

      if ( !success ) {
        manager.reinit(reinitFaultsID);
        for ( auto id : reinitFaultsID ) {
          TaskExample* tmp = static_cast<TaskExample*>(manager.getTask(id));
          tmp->setStepsTotal(combistep*numSteps);
        }
        manager.restoreCombischeme();
        manager.updateCombiParameters();
      }
    }
  }

  valueFile.close();

  /* evaluate solution */
  FullGrid<CombiDataType> fg_eval( leval.size(), leval, boundary);
  manager.gridEval(fg_eval);
  std::string filename( "solution.fg" );
  fg_eval.save( filename );

  std::cout << "exiting" << std::endl;
  manager.exit();
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  /* number of process groups and number of processes per group */
  size_t ngroup, nprocs;

  DimType dim;
  LevelVector lmin, lmax, leval;
  IndexVector p;
  size_t ncombi;
  int numSteps;
  std::vector<int> time_steps_combi;
  FaultsInfo faultsInfo;

  double time_step, time_start, time_end;

  bool plot;

  /* read parameter file (ctparam) */
  const std::string fileName = "settings.ini";
  readParameterFile(fileName, ngroup, nprocs, dim, lmin, lmax,
      leval, p, time_step, time_start, time_end, ncombi,
      faultsInfo, plot);

  numSteps = (time_end-time_start) / time_step;
  for (size_t i = 0; i < ncombi; ++i)
    time_steps_combi.push_back( numSteps - i*numSteps/ncombi - 1 );
  std::reverse(time_steps_combi.begin(), time_steps_combi.end());

  // todo: read in boundary vector from ctparam
  std::vector<bool> boundary(dim, true);

  theMPISystem()->init( ngroup, nprocs );

  // ProcessGroupManager and ProcessManager Code
  if (theMPISystem()->getWorldRank() == theMPISystem()->getManagerRankWorld()) {

    /* create an abstraction of the process groups for the manager's view
     * a pgroup is identified by the ID in gcomm
     */
    ProcessGroupManagerContainer pgroups;
    for (size_t i=0; i<ngroup; ++i) {
      // todo: order of ranks in new group?
      int pgroupRootID(i);
      pgroups.emplace_back(
          std::make_shared< ProcessGroupManager > ( pgroupRootID )
      );
    }

    // ProcessManager Code
      // create DuneTasks
      LoadModel* loadmodel = new LinearLoadModel();

      CombiMinMaxScheme combischeme(dim, lmin, lmax);
      combischeme.createAdaptiveCombischeme();
      combischeme.makeFaultTolerant();
      std::vector<LevelVector> levels = combischeme.getCombiSpaces();
      std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

      TaskContainer tasks;
      std::vector<int> taskIDs;
      for (uint i = 0; i < levels.size(); ++i) {
        Task* t = new TaskExample(dim, levels[i], boundary, coeffs[i],
                                        loadmodel, time_step, numSteps, p, faultsInfo);
        tasks.push_back(t);
        taskIDs.push_back( t->getID() );
      }

      // output of combination setup
      std::cout << "lmin = " << lmin << std::endl;
      std::cout << "lmax = " << lmax << std::endl;
      std::cout << "CombiScheme: " << std::endl;
      std::cout << combischeme << std::endl;

      // create combiparamters
      CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs);

      // create Manager with process groups
      ProcessManager manager(pgroups, tasks, params);

      // combiparameters need to be set before starting the computation
      manager.updateCombiParameters();

      // start calculation
      solveProblem(manager, leval, boundary, time_step, numSteps, time_steps_combi, plot, faultsInfo.sdcMethod_ );
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

  MPI_Finalize();

  return 0;
}
