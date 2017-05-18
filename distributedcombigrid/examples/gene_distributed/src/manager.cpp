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
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/FaultCriterion.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/StaticFaults.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/WeibullFaults.hpp"

// include user specific task. this is the interface to your application
#include "GeneTask.hpp"

using namespace combigrid;
// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(GeneTask)
BOOST_CLASS_EXPORT(StaticFaults)
BOOST_CLASS_EXPORT(WeibullFaults)

BOOST_CLASS_EXPORT(FaultCriterion)


// helper funtion to read a bool vector from string
inline std::vector<bool>& operator>>(std::string str, std::vector<bool>& vec) {
  std::vector<std::string> strs;
  boost::split(strs, str, boost::is_any_of(" "));

  assert(vec.size() == strs.size());

  for (size_t i = 0; i < strs.size(); ++i)
    vec[i] = boost::lexical_cast<bool>(strs[i]);

  return vec;
}


// helper function to output bool vector
inline std::ostream& operator<<(std::ostream& os, const std::vector<bool>& l) {
  os << "[";

  for (size_t i = 0; i < l.size(); ++i)
    os << l[i] << " ";

  os << "]";

  return os;
}

/**
 * This example performs the fault tolerant combination technique on Gene.
 * Multiple fault models can be plugged in to simulate faults during the simulation.
 * For simulation the faults the Sim_FT fault simulator is used. This is only the code of the manager.
 * All slaves execute the modified gene version that is able to communicate with the master.
 */
int main(int argc, char** argv) {
  //MPI_Init(&argc, &argv);
  simft::Sim_FT_MPI_Init(&argc, &argv);

  // read in parameter file
  boost::property_tree::ptree cfg;
  /*
   * Read input from ctparam
   */
  boost::property_tree::ini_parser::read_ini("ctparam", cfg);

  // number of process groups and number of processes per group
  size_t ngroup = cfg.get<size_t>("manager.ngroup");
  size_t nprocs = cfg.get<size_t>("manager.nprocs");

  /* do GENE specific MPI initialization
   * this involves setting up a kind of local communicator in GENE in which
   * all processes (including the manager process) must participate.
   * these communicators will simply be ignored in the framework part, but they
   * will be used by GENE.
   */
  int globalID, globalSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &globalID);
  MPI_Comm_size(MPI_COMM_WORLD, &globalSize);
  const int managerIDworld = globalSize - 1;
  std::cout << "Manager rank " << globalID << "\n";
  assert(globalSize == int(ngroup * nprocs + 1));
  //generate the local communicators for the different process groups
  int color = globalID / nprocs;
  int key = globalID - color * nprocs;
  MPI_Comm lcomm;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &lcomm);

  // Gene creates another comm which we do not need, but it is necessary
  // to execute comm_split again
  MPI_Comm pcomm;
  MPI_Comm_split( MPI_COMM_WORLD, key, color, &pcomm );

  // here the actual MPI initialization
  theMPISystem()->init( ngroup, nprocs, lcomm );
  int nfaults = 0;

  // manager code
  if (theMPISystem()->getWorldRank() == theMPISystem()->getManagerRankWorld()) {
    /* create an abstraction of the process groups for the manager's view
     * a pgroup is identified by the ID in gcomm
     */
    ProcessGroupManagerContainer pgroups;
    //create vector containing the different process groups
    for (size_t i=0; i<ngroup; ++i) {
      // todo: order of ranks in new group?
      int pgroupRootID(i);
      pgroups.emplace_back(
          std::make_shared< ProcessGroupManager > ( pgroupRootID )
      );
    }

    // create load model
    LoadModel* loadmodel = new LinearLoadModel();

    /* generate a list of levelvectors and coefficients
     * CombiTS_CT will generate a valid combination. however, you could
     * also read in a list of levelvectors and coefficients from a file */
    DimType dim = cfg.get<DimType>("ct.dim");
    LevelVector lmin(dim), lmax(dim), leval(dim), leval2(dim);
    IndexVector p(dim);
    std::vector<bool> boundary(dim), hierarchizationDims(dim);
    combigrid::real dt;
    size_t nsteps, ncombi;
    cfg.get<std::string>("ct.lmin") >> lmin;
    cfg.get<std::string>("ct.lmax") >> lmax;
    cfg.get<std::string>("ct.leval") >> leval;
    cfg.get<std::string>("ct.leval2") >> leval2;
    cfg.get<std::string>("ct.p") >> p;
    cfg.get<std::string>("ct.boundary") >> boundary;
    cfg.get<std::string>("ct.hierarchization_dims") >> hierarchizationDims;
    ncombi = cfg.get<size_t>("ct.ncombi");
    std::string basename = cfg.get<std::string>( "preproc.basename" );
    dt = cfg.get<combigrid::real>("application.dt");
    nsteps = cfg.get<size_t>("application.nsteps");
    //read fault values
    FaultsInfo faultsInfo;

    faultsInfo.numFaults_ = cfg.get<int>("faults.num_faults");

    if( faultsInfo.numFaults_ > 0 ){
      faultsInfo.iterationFaults_.resize(faultsInfo.numFaults_);
      faultsInfo.globalRankFaults_.resize(faultsInfo.numFaults_);
      cfg.get<std::string>("faults.iteration_faults") >> faultsInfo.iterationFaults_;
      cfg.get<std::string>("faults.global_rank_faults") >> faultsInfo.globalRankFaults_;
    }
    std::string fg_file_path = cfg.get<std::string>( "ct.fg_file_path" );
    std::string fg_file_path2 = cfg.get<std::string>( "ct.fg_file_path2" );

    // todo: read from parametes file
    real shat = 0.7960;
    real kymin = 0.3000;
    real lx = 4.18760;
    int ky0_ind = 1;

    // check parallelization vector p agrees with nprocs
    IndexType checkProcs = 1;

    for (auto k : p)
      checkProcs *= k;

    assert(checkProcs == IndexType(nprocs));

    // read combi levels and spaces from file as long as min max scheme does
    // not work properly
    std::vector<LevelVector> levels;
    std::vector<combigrid::real> coeffs;
    std::vector<int> fileTaskIDs;
    const bool READ_FROM_FILE = cfg.get<bool>("ct.readspaces");
    if (READ_FROM_FILE) {
      std::ifstream spcfile("spaces.dat");
      std::string line;
      while (std::getline(spcfile, line)) {
        std::stringstream ss(line);
        int id;
        LevelVector l(dim);
        combigrid::real coeff;
        ss >> id;
        for (size_t i = 0; i < dim; ++i)
          ss >> l[i];
        ss >> coeff;

        levels.push_back(l);
        coeffs.push_back(coeff);
        fileTaskIDs.push_back(id);
      }
      spcfile.close();
    } else {
      CombiMinMaxScheme combischeme(dim, lmin, lmax);
      combischeme.createAdaptiveCombischeme();
      levels = combischeme.getCombiSpaces();
      coeffs = combischeme.getCoeffs();
    }

    // output of combination setup
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "lmax = " << lmax << std::endl;
    std::cout << "boundary = " << boundary << std::endl;
    std::cout << "hierarchization_dims = " << hierarchizationDims << std::endl;
    std::cout << "CombiScheme: " << std::endl;
    for (size_t i = 0; i < levels.size(); ++i)
      std::cout << "\t" << levels[i] << " " << coeffs[i] << std::endl;


    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;

    for (size_t i = 0; i < levels.size(); i++) {
      // path to task folder
      std::stringstream ss2;
      ss2 << "../" << basename << fileTaskIDs[i];
      std::string path = ss2.str();
      //create FaultCriterion
      FaultCriterion *faultCrit;
      //create fault criterion
      if(faultsInfo.numFaults_ < 0){
        faultCrit = new WeibullFaults(0.7, abs(faultsInfo.numFaults_), ncombi, true);
      }
      else{ //do not use faults
        faultCrit = new StaticFaults(faultsInfo);
      }
      Task* t = new GeneTask(dim, levels[i], boundary, coeffs[i],
                                loadmodel, path, dt, nsteps,
                                shat, kymin, lx, ky0_ind, p, faultCrit);
      tasks.push_back(t);
      taskIDs.push_back( t->getID() );

    }

    // create combiparamters
    CombiParameters params( dim, lmin, lmax, boundary, levels,
                            coeffs, hierarchizationDims, taskIDs );
    params.setParallelization(p);

    // create Manager with process groups
    ProcessManager manager(pgroups, tasks, params);

    // combiparameters need to be set before starting the computation
    manager.updateCombiParameters();
    bool success = true;
    //theStatsContainer()->setTimerStart("compute");
    for (size_t i = 0; i < ncombi; ++i) {
      //std::cout << "Compute !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n";
      if( i == 0 ){
        /* distribute task according to load model and start computation for
         * the first time */
        //theStatsContainer()->setTimerStart("runfirst");
        success = manager.runfirst();
        //theStatsContainer()->setTimerStop("runfirst");
      } else {
        // run tasks for next time interval
        //if(i==1) theStatsContainer()->setTimerStart("runnext");
        success = manager.runnext();
        //if(i==1) theStatsContainer()->setTimerStop("runnext");
      }
      //success = true;
      //check if fault occured
      if ( !success ) {
        nfaults++;
        std::cout << "failed group detected at combi iteration " << i << std::endl;
       // manager.recover(i, nsteps); has no access toi GENE Task
        std::vector<int> faultsID;
        std::vector< ProcessGroupManagerID> groupFaults;
        manager.getGroupFaultIDs(faultsID, groupFaults);

        /* call optimization code to find new coefficients */
        const std::string prob_name = "interpolation based optimization";
        std::vector<int> redistributeFaultsID, recomputeFaultsID;
        manager.recomputeOptimumCoefficients(prob_name, faultsID, redistributeFaultsID, recomputeFaultsID);
        //time does not need to be updated in gene but maybe in other applications
        for ( auto id : redistributeFaultsID ) {
         GeneTask* tmp = static_cast<GeneTask*>(manager.getTask(id));
         tmp->setStepsTotal((i+1)*nsteps);
         tmp->setCombiStep(i+1);
        }

        for ( auto id : recomputeFaultsID ) {
         GeneTask* tmp = static_cast<GeneTask*>(manager.getTask(id));
         tmp->setStepsTotal((i)*nsteps);
         tmp->setCombiStep(i);
        }
        /* recover communicators*/
        bool failedRecovery = manager.recoverCommunicators(groupFaults);
        /* communicate new combination scheme*/
        if(failedRecovery){
         std::cout << "redistribute \n";
         manager.redistribute(redistributeFaultsID);
        }
        else{
         std::cout << "reinitializing group \n";
         manager.reInitializeGroup(groupFaults,recomputeFaultsID);
        }


        /* if some tasks have to be recomputed, do so*/
        if(!recomputeFaultsID.empty()){
          std::cout << "sending tasks for recompute \n";
         manager.recompute(recomputeFaultsID,failedRecovery,groupFaults); //toDO handle faults in recompute
        }
        std::cout << "updateing Combination Parameters \n";
        //needs to be after reInitialization!
        manager.updateCombiParameters();
//        old version
//        std::vector<int> faultsID;
//        std::vector< ProcessGroupManagerID> groupFaults;
//        manager.getGroupFaultIDs(faultsID, groupFaults);
//
//        /* call optimization code to find new coefficients */
//        const std::string prob_name = "interpolation based optimization";
//        std::vector<int> redistributeFaultsID, recomputeFaultsID;
//        manager.recomputeOptimumCoefficients(prob_name, faultsID,
//                                             redistributeFaultsID, recomputeFaultsID);
//
//
//        //std::cout << "Recover \n";
//
//        /* recover communicators*/
//        bool failedRecovery = manager.recoverCommunicators(groupFaults);
//        //std::cout << "Update \n";
//
//        /* communicate new combination scheme*/
//        manager.updateCombiParameters();
//        //std::cout << "Recompute \n";
//
//        /* if some tasks have to be recomputed, do so*/
//        manager.recompute(recomputeFaultsID); //toDo consider succeeded recovery
//        //std::cout << "Redistribute \n";
//        /* redistribute failed tasks to living groups */
//
//        if(failedRecovery){
//          manager.redistribute(redistributeFaultsID);
//        }
//        else{
//          manager.reInitializeGroup(groupFaults);
//        }
      }
      //combine
      //if(i==0) theStatsContainer()->setTimerStart("combine");
      manager.combine();
      //if(i==0) theStatsContainer()->setTimerStop("combine");
      //postprocessing in case of errors
      if ( !success ){
        /* restore combischeme to its original state
         * and send new combiParameters to all surviving groups */
        manager.restoreCombischeme();
        manager.updateCombiParameters();
      }
    }
    std::cout << "Computation finished evaluating on target grid! \n";
    //theStatsContainer()->setTimerStop("compute");

    // evaluate solution on the grid defined by leval
    //theStatsContainer()->setTimerStart("parallelEval");
    manager.parallelEval( leval, fg_file_path, 0 );
    //theStatsContainer()->setTimerStop("parallelEval");

    // evaluate solution on the grid defined by leval2
    //theStatsContainer()->setTimerStart("parallelEval2");
    manager.parallelEval( leval2, fg_file_path2, 0 );
    //theStatsContainer()->setTimerStop("parallelEval2");

    // send exit signal to workers in order to enable a clean program termination
    manager.exit();

    // save stats
    //theStatsContainer()->save("times.dat");
  }
  if( ENABLE_FT ){
    WORLD_MANAGER_EXCLUSIVE_SECTION{
      std::cout << "The number of detected faults during the simulation is " << nfaults << "\n";

      std::cout << "Program finished successfully" << std::endl;
      std::cout << "To avoid problems with hanging killed processes, we exit with "
          << "MPI_Abort()" << std::endl;
      MPI_Abort( MPI_COMM_WORLD, 0 );
    }
  }
  std::cout << "Finalize manager \n";
  simft::Sim_FT_MPI_Finalize();

  return 0;
}
