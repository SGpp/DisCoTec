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
#include "GeneTask.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(GeneTask)

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  // read in parameter file
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini("ctparam", cfg);

  // number of process groups and number of processes per group
  size_t ngroup = cfg.get<size_t>("manager.ngroup");
  size_t nprocs = cfg.get<size_t>("manager.nprocs");

  //TODO: put creation of communicators in distinct functions

  /* determine global rank of each process
   * the manager process always has the highest rank
   * all other processes are worker processes */
  int globalID, globalSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &globalID);
  MPI_Comm_size(MPI_COMM_WORLD, &globalSize);
  const int managerIDworld = globalSize - 1;

  assert(globalSize == int(ngroup * nprocs + 1));

  /* create a local communicator for each process group
   * lcomm is the local communicator of its own process group for each worker process
   * for manager, lcomm is a group which contains only manager process and can be ignored
   */
  MPI_Comm lcomm;
  int color = globalID / nprocs;
  int key = globalID - color * nprocs;
  int lrank;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &lcomm);
  MPI_Comm_rank(lcomm, &lrank);

  // Gene creates another comm which we do not need, but it is necessary
  // to execute comm_split again
  MPI_Comm pcomm;
  MPI_Comm_split( MPI_COMM_WORLD, key, color, &pcomm );

  /* create global communicator which contains only the manager and the master
   * process of each process group
   * the master processes of the process groups are the processes which have
   * rank 0 in lcomm
   * this communicator is used for communication between master processes of the
   * process groups and the manager and the master processes to each other
   */
  MPI_Group worldGroup;
  MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

  std::vector<int> ranks(ngroup + 1);

  for (size_t i = 0; i < ngroup; i++) {
    ranks[i] = i * nprocs;
  }

  ranks.back() = managerIDworld;

  MPI_Group rootGroup;
  MPI_Group_incl(worldGroup, (int) ranks.size(), &ranks[0], &rootGroup);

  MPI_Comm gcomm;
  MPI_Comm_create(MPI_COMM_WORLD, rootGroup, &gcomm);

  int managerIDgcomm = MPI_PROC_NULL;
  int grank;

  if (gcomm != MPI_COMM_NULL) {
    int gcommSize;
    MPI_Comm_size(gcomm, &gcommSize);
    managerIDgcomm = gcommSize - 1;
    MPI_Comm_rank(gcomm, &grank);
  }

  std::cout << "reached here" << std::endl;

  //CombiCom::registerDistributedGlobalReduceCommmunicator(nprocs);

  std::cout << "but not here" << std::endl;

  // manager code
  if (globalID == managerIDworld) {
    /* create an abstraction of the process groups for the manager's view
     * a pgroup is identified by the ID in gcomm
     */
    ProcessGroupManagerContainer pgroups;

    for (size_t i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      ProcessGroupManager grp(managerIDgcomm, pgroupRootID, gcomm);
      pgroups.push_back(grp);
    }

    // create load model
    LoadModel* loadmodel = new LinearLoadModel();

    /* generate a list of levelvectors and coefficients
     * CombiTS_CT will generate a valid combination. however, you could
     * also read in a list of levelvectors and coefficients from a file */
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
    std::string basename = cfg.get<std::string>( "preproc.basename" );
    dt = cfg.get<combigrid::real>("application.dt");
    nsteps = cfg.get<size_t>("application.nsteps");
    std::string fg_file_prefix = cfg.get<std::string>( "ct.fg_file_prefix" );

    // todo: read in boundary vector from ctparam
    // set boundary vector
    std::vector<bool> boundary( 5 );
    boundary[0] = true;
    boundary[1] = false;
    boundary[2] = true;
    boundary[3] = false;
    boundary[4] = true;

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

      Task* t = new GeneTask(dim, levels[i], boundary, coeffs[i],
                                loadmodel, path, dt, nsteps );
      tasks.push_back(t);
      taskIDs.push_back( t->getID() );
    }

    // create combiparamters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs );

    // create Manager with process groups
    ProcessManager manager(pgroups, tasks, params, managerIDgcomm, gcomm);

    // combiparameters need to be set before starting the computation
    manager.updateCombiParameters();

    std::ofstream myfile;
    myfile.open("out/solution.dat");

    for (size_t i = 0; i < ncombi; ++i) {
      if( i == 0 ){
        /* distribute task according to load model and start computation for
         * the first time */
        manager.runfirst();
      } else {
        // run tasks for next time interval
        manager.runnext();
      }

      //manager.combine();

      // evaluate solution
      FullGrid<CombiDataType> fg_eval(dim, leval, boundary);
      manager.gridEval(fg_eval);

      std::cout << "left grideval" << std::endl;



      // write solution to file
      std::string filename = fg_file_prefix
          + boost::lexical_cast<std::string>( i ) + ".dat";
      fg_eval.save( filename );
    }

    myfile.close();

    // send exit signal to workers in order to enable a clean program termination
    manager.exit();
  }

  MPI_Finalize();

  return 0;
}
