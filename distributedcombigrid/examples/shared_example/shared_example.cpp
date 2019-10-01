/*
 * combi_example.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: heenemo
 */
#include <mpi.h>
#include <omp.h>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/export.hpp>
#include <string>
#include <vector>
//Testing your thread affinty settings
#include <sched.h>
#include <numa.h> 

// compulsory includes for basic functionality
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/FaultCriterion.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/StaticFaults.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/WeibullFaults.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
// include user specific task. this is the interface to your application
#ifndef ADVECTION
#include "TaskExample.hpp"
#else 
#include "AdaptiveTaskExample.hpp"
#endif


using namespace combigrid;



// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskExample)
BOOST_CLASS_EXPORT(StaticFaults)
BOOST_CLASS_EXPORT(WeibullFaults)
BOOST_CLASS_EXPORT(FaultCriterion)
int main(int argc, char** argv) {
  int _threadlevelresult;
  MPI_Init_thread(&argc, &argv,MPI_THREAD_MULTIPLE,&_threadlevelresult);
  if(_threadlevelresult<MPI_THREAD_FUNNELED)
  {
    std::cout << "your mpi implementation claims not to support threads, even though"<<
    "this thread level should function as it was not there for mpi"<<_threadlevelresult;
    return -2;
  }
  int tmprank;
  MPI_Comm_rank(MPI_COMM_WORLD,&tmprank);//just for debug output
  if(tmprank==0)
  {
    for(int i=0;i<argc;i++)
    {
      std::cout << argv[i]<<" ";
    }
    std::cout <<"\n";
  }

  /* when using timers (TIMING is defined in Stats), the Stats class must be
   * initialized at the beginning of the program. (and finalized in the end)
   */
  Stats::initialize();
  
  char opt;
  std::string paramfile="ctparam";
  std::string statfile="timers.json";
  std::string filesuffix="";
  int thread_override=-1;
  int allreduce_threads=1;
  bool async=false;

	while ((opt = getopt(argc, argv, "Aa:c:n:o:t:")) != -1) {
		switch (opt) {
			case 'c':
				paramfile= optarg;
        break;
      case 't':
				statfile= optarg;
        break;
      case 'n':
        thread_override=std::stoi(optarg);
        break;
      case 'o':
        filesuffix=optarg;
        break;
      case 'a':
        allreduce_threads=std::stoi(optarg);
        break;
      case 'A':
        async=true;
        break;
      default:
        std::cout <<" unsupported option -"<< opt <<" \n Usage [-c] ctparam_file [-t] timing file";
    }
  }

  // read in parameter file
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(paramfile, cfg);
  

  // number of process groups and number of processes per group
  size_t ngroup = cfg.get<size_t>("manager.ngroup");
  size_t nprocs = cfg.get<size_t>("manager.nprocs");
  size_t nthreads= cfg.get<size_t>("manager.nthreads");
  omp_set_num_threads(thread_override>0?thread_override:nthreads);

  // divide the MPI processes into process group and initialize the
  // corresponding communicators
  auto thempi=theMPISystem();
  thempi->init(ngroup, nprocs);
  if(_threadlevelresult==MPI_THREAD_MULTIPLE){
    thempi->initThreadedGlobalreduce(allreduce_threads,_threadlevelresult,async);
  }

  // this code is only executed by the manager process
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    std::cout << "using " << omp_get_max_threads()<<" Threads\n"; 
    char hostname[HOST_NAME_MAX+1] ;
    gethostname(hostname,HOST_NAME_MAX+1);
    std::cout <<"HOST rank"<<theMPISystem()->getWorldRank()<<" on node "<<hostname<<std::endl;
    Stats::startEvent("total time");
    /* create an abstraction of the process groups for the manager's view
    * a pgroup is identified by the ID in gcomm
     */
    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    // create load model
    std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

    /* read in parameters from ctparam */
    DimType dim = cfg.get<DimType>("ct.dim");
    LevelVector lmin(dim), lmax(dim), leval(dim);
    IndexVector p(dim);
    combigrid::real dt;
    size_t nsteps, ncombi;
    cfg.get<std::string>("ct.lmin") >> lmin;
    cfg.get<std::string>("ct.lmax") >> lmax;
    cfg.get<std::string>("ct.leval") >> leval;
    cfg.get<std::string>("ct.p") >> p;
    ncombi = cfg.get<size_t>("ct.ncombi");
    dt = cfg.get<combigrid::real>("application.dt");
    nsteps = cfg.get<size_t>("application.nsteps");

    // todo: read in boundary vector from ctparam
    std::vector<bool> boundary(dim, true);

    // check whether parallelization vector p agrees with nprocs
    IndexType checkProcs = 1;
    for (auto k : p) checkProcs *= k;
    assert(checkProcs == IndexType(nprocs));

    /* generate a list of levelvectors and coefficients
     * CombiMinMaxScheme will create a classical combination scheme.
     * however, you could also read in a list of levelvectors and coefficients
     * from a file */
    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createAdaptiveCombischeme();
    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

    // output combination scheme
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "lmax = " << lmax << std::endl;
    std::cout << "CombiScheme: " << std::endl;
    std::cout << combischeme << std::endl;

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskExample(dim, levels[i], boundary, coeffs[i], loadmodel.get(), dt, nsteps, p);
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi, 1);
    params.setParallelization(p);
    // create abstraction for Manager
    // inserted LinearLoadModel as default as it was not compiling else wise
    // trying to use pointer created
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    std::cout << "set up component grids and run until first combination point" << std::endl;

    /* distribute task according to load model and start computation for
     * the first time */
    Stats::startEvent("manager run first");
    manager.runfirst();
    Stats::stopEvent("manager run first");

    for (size_t i = 0; i < ncombi; ++i) {
      //int cpu=sched_getcpu();
      //std::cout <<"NUMA master cpu: "<<cpu<<"numa: "<<numa_node_of_cpu(cpu)<<std::endl; 
      Stats::startEvent("combine");
      manager.combine();
      Stats::stopEvent("combine");

      // evaluate solution and
      // write solution to file
      
      Stats::startEvent("manager write solution");
      //manager.parallelEval(leval, filename, 0);
      Stats::stopEvent("manager write solution");

      std::cout << "run until combination point " << i + 1 << std::endl;

      // run tasks for next time interval
      Stats::startEvent("manager run");
      manager.runnext();
      Stats::stopEvent("manager run");
    }
    Stats::startEvent("manager write solution");
    auto localFilenameInd=paramfile.find_first_of('/');
    std::string filename("out/solution_" + paramfile.substr(localFilenameInd+1)+filesuffix+ ".dat");
    std::cout <<filename<<"\n";
    manager.parallelEval(leval,filename, 0);
    Stats::stopEvent("manager write solution");

    // send exit signal to workers in order to enable a clean program termination
    manager.exit();
    Stats::stopEvent("total time");
  }

  // this code is only execute by the worker processes
  else {
    // create abstraction of the process group from the worker's view
    char hostname[HOST_NAME_MAX+1] ;
    gethostname(hostname,HOST_NAME_MAX+1);
    std::cout <<"HOST rank"<<theMPISystem()->getGlobalRank()<<" on node "<<hostname<<std::endl;

    ProcessGroupWorker pgroup;

    // wait for instructions from manager
    SignalType signal = -1;

    while (signal != EXIT) signal = pgroup.wait();
  }
  
  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write(statfile);

  MPI_Finalize();

  return 0;
}
