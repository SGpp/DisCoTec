#include <mpi.h>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/export.hpp>
#include <string>
#include <vector>
#include <unistd.h>

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

#include <deal.II/base/timer.h>
#include <deal.II/lac/la_parallel_vector.h>

//TODO: global variables needed in hyperdeal
typedef double Number;
const unsigned int dim    = 2;
const unsigned int degree = 1; /*dummy value*/
typedef dealii::VectorizedArray<Number, 1> VectorizedArrayType;
typedef dealii::LinearAlgebra::distributed::Vector<Number> VectorType;

//hyperdeal includes
#include <hyper.deal.combi/include/functionalities/dynamic_convergence_table.h>
#include <hyper.deal.combi/include/functionalities/vector_dummy.h>
#include <hyper.deal.combi/applications/advection_reference_dealii/include/application.h>

typedef Application<dim, degree, degree + 1, Number, VectorizedArrayType, VectorType> Problem;

hyperdeal::DynamicConvergenceTable table;

// include user specific task. this is the interface to your application

const bool do_combine=true;


#include "TaskExample.hpp"
#include "TaskEnsemble.hpp"

#include "DataConverter.cpp"
#include "DataConverter.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskExample)
BOOST_CLASS_EXPORT(TaskEnsemble)
BOOST_CLASS_EXPORT(StaticFaults)
BOOST_CLASS_EXPORT(WeibullFaults)
BOOST_CLASS_EXPORT(FaultCriterion)

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  /* when using timers (TIMING is defined in Stats), the Stats class must be
   * initialized at the beginning of the program. (and finalized in the end)
   */
  Stats::initialize();

  // read in parameter file
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini("ctparam", cfg);

  // number of process groups and number of processes per group
  size_t ngroup = cfg.get<size_t>("manager.ngroup");
  size_t nprocs = cfg.get<size_t>("manager.nprocs");
  
  LevelVector lmin(dim), lmax(dim), leval(dim);
  cfg.get<std::string>("ct.lmin") >> lmin;
  cfg.get<std::string>("ct.lmax") >> lmax;
  cfg.get<std::string>("ct.leval") >> leval;

  // divide the MPI processes into process group and initialize the
  // corresponding communicators
  theMPISystem()->init(ngroup, nprocs);

  // this code is only executed by the manager process
  WORLD_MANAGER_EXCLUSIVE_SECTION {
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
    IndexVector p(dim);
    combigrid::real dt;
    size_t ncombi;
    bool isdg=("FE_DGQ"==cfg.get<std::string>("ct.FE","FE_Q"));
    std::cout << "Is dg:"<<isdg;
    cfg.get<std::string>("ct.p") >> p;
    ncombi = cfg.get<size_t>("ct.ncombi");
    dt = cfg.get<combigrid::real>("application.dt");
    

    Converter converter;
    // TODO: read in boundary vector from ctparam
    std::vector<bool> boundary(dim, true);
    // use no hierarchization in this example
    std::vector<bool> hierarchizationDims(dim, true);

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
    std::cout << "Levels:" << std::endl;
    //std::cout << levels << std::endl;

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      std::string file_name = "deal_config/"+std::string("p")+std::to_string(i)+".json";
      converter.toJSONForDealII("ctparam",file_name, levels[i],i);
      usleep(10);
      Task* t;

      if(!isdg)
        t = new TaskExample(dim, levels[i], boundary, coeffs[i], loadmodel.get(),file_name, dt, p);
      else
        t = new TaskEnsemble(dim, levels[i], boundary, coeffs[i], loadmodel.get(),file_name, dt,  p);
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }
    
    //gridout.write_pvtu

    // create combiparameters
    //dim, lmin,lmax kann ich aus json lesen.
    //boundary i have to read from ct params but how
    // levels is just array of levels between lim and lmax
    //taskID for each level we create a new task and here are the IDs
    CombiParameters params( dim, lmin, lmax, boundary, levels,
                            coeffs, hierarchizationDims, taskIDs, ncombi, tasks[0]->getNumGrids());//, reduceCombinationDimsLmin, reduceCombinationDimsLmax);

    params.setParallelization(p);
    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    /* distribute task according to load model and start computation for
     * the first time */
    Stats::startEvent("manager run first");
    manager.runfirst();
    Stats::stopEvent("manager run first");
    
    std::chrono::nanoseconds duration_combination=std::chrono::nanoseconds::zero();
    for (size_t i = 0; i < ncombi; ++i) {
      Stats::startEvent("combine");
      manager.combine();
      Stats::stopEvent("combine");

      


      // run tasks for next time interval
      Stats::startEvent("manager run");
      manager.runnext();
      Stats::stopEvent("manager run");
      table.print(false);
    }
    // evaluate solution and
      // write solution to file
      std::string filename("out/"+cfg.get<std::string>("ct.FE","FE_Q")+"/cs"  +"mi_"+toString(lmin)+"_ma_"+toString(lmax)+"_ev_"+toString(leval)+".dat");
      Stats::startEvent("manager write solution");
      manager.parallelEval(leval, filename, 0);
      Stats::stopEvent("manager write solution");
    //table.stop_and_set("time->fullgrid");
    //table.print(false);
    // send exit signal to workers in order to enable a clean program termination
    manager.exit();
    {

      DataOut<1> gridout;
      std::vector<std::pair<double, std::string>> pvd_labels;
      for(unsigned int t = 0; t< ncombi; t++ )
      {
        std::vector<std::string> pvtu_labels;
        for(unsigned int l = 0; l< levels.size(); l++ )
          pvtu_labels.push_back("solution_" +  std::to_string(l) + "_" + Utilities::int_to_string(t, 3) + ".vtu");

        pvd_labels.emplace_back(t, "solution_" + std::to_string(t) + ".pvtu");
        std::ofstream outfile ("solution/solution_" + std::to_string(t) + ".pvtu");
        gridout.write_pvtu_record(outfile, pvtu_labels);
      }
      std::ofstream outfile ("solution/solution.pvd");
      DataOutBase::write_pvd_record(outfile, pvd_labels);
      
    }
  }

  // this code is only execute by the worker processes
  else {
    // create abstraction of the process group from the worker's view
    ProcessGroupWorker pgroup;

    // wait for instructions from manager
    SignalType signal = -1;

    while (signal != EXIT) signal = pgroup.wait();
  }

  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write("timers/"+cfg.get<std::string>("ct.FE","FE_Q")+"mi_"+toString(lmin)+"_ma_"+toString(lmax)+"_ev_"+toString(leval)+".json");

  MPI_Finalize();

  return 0;
}
