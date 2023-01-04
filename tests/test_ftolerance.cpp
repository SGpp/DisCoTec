#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <vector>
#include <set>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <iostream>
#include <vector>
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/serialization/export.hpp>

// compulsory includes for basic functionality
#include "task/Task.hpp"
#include "utils/Types.hpp"
#include "combischeme/CombiMinMaxScheme.hpp"
#include "hierarchization/CombiLinearBasisFunction.hpp"
#include "fullgrid/FullGrid.hpp"
#include "loadmodel/LinearLoadModel.hpp"
#include "loadmodel/LearningLoadModel.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupManager.hpp"
#include "manager/ProcessGroupWorker.hpp"
#include "manager/ProcessManager.hpp"
#include "fault_tolerance/LPOptimizationInterpolation.hpp"
#include "fault_tolerance/FaultCriterion.hpp"
#include "fault_tolerance/StaticFaults.hpp"
#include "fault_tolerance/WeibullFaults.hpp"
#include "test_helper.hpp"
#include "utils/Config.hpp"
#include "fault_tolerance/FTUtils.hpp"
#include "fullgrid/DistributedFullGrid.hpp"

// this is necessary for correct function of task serialization
#include "utils/BoostExports.hpp"

using namespace combigrid;

/* functor for exact solution */
class TestFn {
 public:
  // function value
  double operator()(std::vector<double>& coords, double t) {
    double exponent = 0;
    for (DimType d = 0; d < static_cast<DimType>(coords.size()); ++d) {
      coords[d] = std::fmod(1.0 + std::fmod(coords[d] - t, 1.0), 1.0);
      exponent -= std::pow(coords[d] - 0.5, 2);
    }
    return std::exp(exponent * 100.0) * 2;
  }
};
class TaskAdvFDM : public combigrid::Task {
 public:
  TaskAdvFDM(LevelVector& l, std::vector<BoundaryType>& boundary, real coeff, LoadModel* loadModel, real dt,
             size_t nsteps,
             FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(2, l, boundary, coeff, loadModel, faultCrit),
        dt_(dt),
        nsteps_(nsteps),
        stepsTotal_(0),
        combiStep_(0) {}

  void init(CommunicatorType lcomm,
            std::vector<IndexVector> decomposition = std::vector<IndexVector>()) {
    // only use one process per group
    std::vector<int> p(getDim(), 1);

    dfg_ =
        new DistributedFullGrid<CombiDataType>(getDim(), getLevelVector(), lcomm, getBoundary(), p);
    phi_.resize(dfg_->getNrElements());

    for (IndexType li = 0; li < dfg_->getNrElements(); ++li) {
      std::vector<double> coords(getDim());
      dfg_->getCoordsGlobal(li, coords);

      double exponent = 0;
      for (DimType d = 0; d < getDim(); ++d) {
        exponent -= std::pow(coords[d] - 0.5, 2);
      }
      dfg_->getData()[li] = std::exp(exponent * 100.0) * 2;
    }

    
  }


  void run(CommunicatorType lcomm) {
    // velocity vector
    std::vector<CombiDataType> u(getDim());
    u[0] = 1;
    u[1] = 1;

    // gradient of phi
    std::vector<CombiDataType> dphi(getDim());

    IndexType l0 = dfg_->length(0);
    IndexType l1 = dfg_->length(1);
    double h0 = 1.0 / (double)l0;
    double h1 = 1.0 / (double)l1;

    for (size_t i = 0; i < nsteps_; ++i) {
      phi_.swap(dfg_->getElementVector());

      for (IndexType li = 0; li < dfg_->getNrElements(); ++li) {

        IndexVector ai(getDim());
        dfg_->getGlobalVectorIndex(li, ai);

        // west neighbor
        IndexVector wi = ai;
        wi[0] = (l0 + wi[0] - 1) % l0;
        IndexType lwi = dfg_->getGlobalLinearIndex(wi);

        // south neighbor
        IndexVector si = ai;
        si[1] = (l1 + si[1] - 1) % l1;
        IndexType lsi = dfg_->getGlobalLinearIndex(si);

        // calculate gradient of phi with backward differential quotient
        dphi[0] = (phi_[li] - phi_[lwi]) / h0;
        dphi[1] = (phi_[li] - phi_[lsi]) / h1;

        CombiDataType u_dot_dphi = u[0] * dphi[0] + u[1] * dphi[1];
        dfg_->getData()[li] = phi_[li] - u_dot_dphi * dt_;
      }
    }

    stepsTotal_ += nsteps_;
    setFinished(true);
    decideToKill();
    ++combiStep_;
  }

  void setStepsTotal( size_t stepsTotal ){
    stepsTotal_ = stepsTotal;
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) { return *dfg_; }

  void setZero() {}

  ~TaskAdvFDM() {
    if (dfg_ != NULL) delete dfg_;
  }

 protected:
  TaskAdvFDM() : combiStep_(0), dfg_(NULL) {}

 private:
  friend class boost::serialization::access;

  real dt_;
  size_t nsteps_;
  size_t stepsTotal_;
  size_t combiStep_;
  DistributedFullGrid<CombiDataType>* dfg_;
  std::vector<CombiDataType> phi_;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    // ar& boost::serialization::make_nvp(
    // BOOST_PP_STRINGIZE(*this),boost::serialization::base_object<Task>(*this));
    ar& boost::serialization::base_object<Task>(*this);
    ar& dt_;
    ar& nsteps_;
    ar& stepsTotal_;
  }

  void decideToKill(){ //toDo check if combiStep should be included in task and sent to process groups in case of reassignment
    using namespace std::chrono;

    int globalRank = theMPISystem()->getGlobalRank();
    // MPI_Comm_rank(lcomm, &lrank);
    // MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    //theStatsContainer()->setTimerStop("computeIterationRank" + std::to_string(globalRank));
    //duration<real> dur = high_resolution_clock::now() - startTimeIteration_;
    //real t_iter = dur.count();
    //std::cout << "Current iteration took " << t_iter << "\n";

    //theStatsContainer()->setTimerStart("computeIterationRank" + std::to_string(globalRank));


    //check if killing necessary
    //std::cout << "failNow result " << failNow(globalRank) << " at rank: " << globalRank <<" at step " << combiStep_ << "\n" ;
    //real t = dt_ * nsteps_ * combiStep_;
    if (combiStep_ != 0 && faultCriterion_->failNow(combiStep_, -1.0, globalRank)){
          std::cout<<"Rank "<< globalRank <<" failed at iteration "<<combiStep_<<std::endl;
          StatusType status=PROCESS_GROUP_FAIL;/*
          MASTER_EXCLUSIVE_SECTION{
            simft::Sim_FT_MPI_Send( &status, 1, MPI_INT,  theMPISystem()->getManagerRank(), TRANSFER_STATUS_TAG,
                              theMPISystem()->getGlobalCommFT() );
          }*/
          theMPISystem()->sendFailedSignal();
          //simft::Sim_FT_kill_me();
    }
  }

};

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskAdvFDM)

bool checkFtolerance(bool useCombine, bool useFG, double l0err, double l2err, size_t ncombi,
                     int nfaults, int iteration_faults, int global_rank_faults) {
  int size = useFG ? 2 : 7;

  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) {
    return true;
  }
  bool finished = false;

  combigrid::Stats::initialize();

  size_t ngroup = useFG ? 1 : 6;
  size_t nprocs = 1;
  DimType dim = 2;
  LevelVector lmin(dim, useFG ? 6 : 3);
  LevelVector lmax(dim,6), leval(dim,6);
  //size_t ncombi = 3;
  size_t nsteps = 1;
  combigrid::real dt = useFG ? 0.0001 : 0.01;

  std::vector<BoundaryType> boundary(dim, 2);
  theMPISystem()->initWorldReusable(comm, ngroup, nprocs);

  BOOST_TEST_CHECKPOINT("Initialize checkFtolerance");

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < ngroup ; ++i) {
      int pgroupRootID(boost::numeric_cast<int>(i));
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

  CombiMinMaxScheme combischeme(dim, lmin, lmax);
  combischeme.createAdaptiveCombischeme();

 #ifdef ENABLEFT
  FaultsInfo faultsInfo;
  faultsInfo.numFaults_ = nfaults;
  IndexVector glob_ranks(nfaults,global_rank_faults);
  IndexVector it_fault(nfaults,iteration_faults);
  faultsInfo.globalRankFaults_ = glob_ranks;
  faultsInfo.iterationFaults_ = it_fault;
  combischeme.makeFaultTolerant();
  #endif

  std::vector<LevelVector> levels = combischeme.getCombiSpaces();
  std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

  BOOST_REQUIRE(true); //if things go wrong weirdly, see where things go wrong

    std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LinearLoadModel>(new LinearLoadModel());

    /* create Tasks */
   TaskContainer tasks;
    std::vector<size_t> taskIDs;

    #ifdef ENABLEFT
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
      Task* t = new TaskAdvFDM( levels[i], boundary, coeffs[i],
                                loadmodel.get(), dt, nsteps, faultCrit);
        tasks.push_back(t);
        taskIDs.push_back(t->getID());
      }

      /* create combi parameters */
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi, 1);


    /* create Manager with process groups */
    //ProcessManager manager( pgroups, tasks, params );
     ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));
    /* send combi parameters to workers */
    manager.updateCombiParameters();

    /* distribute task according to load model and start computation for
     * the first time */
    bool success = manager.runfirst();
for (size_t i = 1; i < ncombi; ++i){

      if ( !success ) {
        std::cout << "failed group detected at combi iteration " << i-1<< std::endl;
        //manager.recover();

        std::vector<size_t> faultsID;

        //vector with pointers to managers of failed groups
        std::vector< ProcessGroupManagerID> groupFaults;
        manager.getGroupFaultIDs(faultsID, groupFaults);

        /* call optimization code to find new coefficients */
        const std::string prob_name = "interpolation based optimization";
        std::vector<size_t> redistributeFaultsID, recomputeFaultsID;
        manager.recomputeOptimumCoefficients(prob_name, faultsID,
                                             redistributeFaultsID, recomputeFaultsID);
        for ( auto id : redistributeFaultsID ) {
          TaskAdvFDM* tmp = static_cast<TaskAdvFDM*>(manager.getTask(id));
          tmp->setStepsTotal(i*nsteps);
        }

        for ( auto id : recomputeFaultsID ) {
          TaskAdvFDM* tmp = static_cast<TaskAdvFDM*>(manager.getTask(id));
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
    #else
      for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskAdvFDM(levels[i], boundary, coeffs[i], loadmodel.get(), dt, nsteps);
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }
    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    /* distribute task according to load model and start computation for
     * the first time */
    manager.runfirst();

    for (size_t it = 0; it < ncombi; ++it) {
      if (useCombine) {
        manager.combine();
      }
      manager.runnext();
    }
    #endif


// evaluate solution
    FullGrid<CombiDataType> fg_eval(dim, leval, boundary);
    manager.gridEval(fg_eval);

    // exact solution
    TestFn f;
    FullGrid<CombiDataType> fg_exact(dim, leval, boundary);
    fg_exact.createFullGrid();
    for (IndexType li = 0; li < fg_exact.getNrElements(); ++li) {
      std::vector<double> coords(dim);
      fg_exact.getCoords(li, coords);
      fg_exact.getData()[li] = f(coords, (double)((1 + ncombi) * nsteps) * dt);
    }

    // calculate error
    fg_exact.add(fg_eval, -1);
    printf("LP Norm: %f\n", fg_exact.getlpNorm(0));
    printf("LP Norm2: %f\n", fg_exact.getlpNorm(2));
    // results recorded previously
    BOOST_CHECK(fabs( fg_exact.getlpNorm(0) - l0err) < TestHelper::higherTolerance);
    BOOST_CHECK(fabs( fg_exact.getlpNorm(2) - l2err) < TestHelper::higherTolerance);

    /* send exit signal to workers in order to enable a clean program termination */
    manager.exit();
    finished = true;
  }

  else {
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    while (signal != EXIT) {
      signal = pgroup.wait();
      BOOST_TEST_CHECKPOINT("Last Successful Worker Signal " + std::to_string(signal));
    }
    BOOST_TEST_ERROR("worker never cleanly exits! " + std::to_string(signal));
    finished = true;
  }

  if( ENABLE_FT ){
    WORLD_MANAGER_EXCLUSIVE_SECTION{
      std::cout << "Program finished successfully" << std::endl;
      std::cout << "To avoid problems with hanging killed processes, we exit with "
                << "MPI_Abort()" << std::endl;
      MPI_Abort( MPI_COMM_WORLD, 0 );
    }
     //simft::Sim_FT_MPI_Finalize();
  }

  combigrid::Stats::finalize();
  MPI_Barrier(comm);
  BOOST_CHECK(finished);
  return finished;
}

BOOST_FIXTURE_TEST_SUITE(ftolerance, TestHelper::BarrierAtEnd,
                         *boost::unit_test::timeout(60) * boost::unit_test::disabled())

#ifdef ENABLEFT
BOOST_AUTO_TEST_CASE(test_0, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::timeout(20)) {
  checkFtolerance(false, false, 0.295448, 2.283454, 3, 1, 1, 3);
  MPI_Barrier(MPI_COMM_WORLD);
}
#else
BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::timeout(40)) {
  // use recombination
  BOOST_CHECK_NO_THROW(checkFtolerance(true, false, 0.973230, 7.820831, 100, 0, 0, 0));
  MPI_Barrier(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(test_2, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::timeout(60)) {
  // don't use recombination
  BOOST_CHECK_NO_THROW(checkFtolerance(false, false, 1.428688, 10.692053, 100, 0, 0, 0));
  MPI_Barrier(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(test_3, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::timeout(80)) {
  // calculate solution on fullgrid
  BOOST_CHECK_NO_THROW(checkFtolerance(false, true, 0.060661, 0.350710, 100, 0, 0, 0));
  MPI_Barrier(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(test_4, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::timeout(40)) {
  // use recombination
  BOOST_CHECK_NO_THROW(checkFtolerance(true, false, 0.043363, 0.300988, 0, 0, 0, 0));
  MPI_Barrier(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(test_5, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::timeout(60)) {
  // don't use recombination
  BOOST_CHECK_NO_THROW(checkFtolerance(false, false, 0.043363, 0.300988, 0, 0, 0, 0));
  MPI_Barrier(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(test_6, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::timeout(80)) {
  // calculate solution on fullgrid
  BOOST_CHECK_NO_THROW(checkFtolerance(false, true, 0.000623, 0.003553, 0, 0, 0, 0));
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif
BOOST_AUTO_TEST_SUITE_END()
