#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mpi.h>
#include <iostream>
#include <complex>
#include <cstdarg>
#include <vector>
#include <cmath>

#include <boost/serialization/export.hpp>
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Config.hpp"
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
#include "test_helper.hpp"

/* functor for exact solution */
class TestFn {
public:
  // function value
  double operator()(std::vector<double>& coords, double t) {
    double exponent = 0;
    for (DimType d = 0; d < coords.size(); ++d) {
      coords[d] = std::fmod(1.0 + std::fmod(coords[d] - t, 1.0), 1.0);
      exponent -= std::pow(coords[d] - 0.5, 2);
    }
    return std::exp(exponent*100.0) * 2;
  }
};

/* simple task class to sequentialy solve the 2D advection equation with
 * periodic boundary conditions using the finite difference and explicit Euler
 * methods */
class TaskAdvectionFDM : public combigrid::Task {
public:
  TaskAdvectionFDM(LevelVector& l, std::vector<bool>& boundary,
                   real coeff, LoadModel* loadModel, real dt, size_t nsteps) :
    Task(2, l, boundary, coeff, loadModel), dt_(dt), nsteps_(nsteps) {
  }

  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition = std::vector<IndexVector>()) {
    // only use one process per group
    IndexVector p(getDim(), 1);
    dfg_ = new DistributedFullGrid<CombiDataType>(getDim(), getLevelVector(),
                                                  lcomm, getBoundary(), p);
    phi_.resize(dfg_->getNrElements());

    for (IndexType li = 0; li < dfg_->getNrElements(); ++li) {
      std::vector<double> coords(getDim());
      dfg_->getCoordsGlobal(li, coords);

      double exponent = 0;
      for (DimType d = 0; d < getDim(); ++d) {
        exponent -= std::pow(coords[d] - 0.5, 2);
      }
      dfg_->getData()[li] = std::exp(exponent*100.0) * 2;
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
        wi[0] = (l0 + wi[0]-1) % l0;
        IndexType lwi = dfg_->getGlobalLinearIndex(wi);

        // south neighbor
        IndexVector si = ai;
        si[1] = (l1 + si[1]-1) % l1;
        IndexType lsi = dfg_->getGlobalLinearIndex(si);

        // calculate gradient of phi with backward differential quotient
        dphi[0] = (phi_[li] - phi_[lwi]) / h0;
        dphi[1] = (phi_[li] - phi_[lsi]) / h1;

        CombiDataType u_dot_dphi = u[0] * dphi[0] + u[1] * dphi[1];
        dfg_->getData()[li] = phi_[li] - u_dot_dphi * dt_;
      }
    }

    setFinished(true);
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r,
                   CommunicatorType lcomm, int n = 0) {
    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) {
    return *dfg_;
  }

  void setZero() {

  }

protected:
  TaskAdvectionFDM() :
    dfg_(NULL) {
  }

  ~TaskAdvectionFDM() {
    if (dfg_ != NULL)
      delete dfg_;
  }

private:
  friend class boost::serialization::access;

  DistributedFullGrid<CombiDataType>* dfg_;
  real dt_;
  size_t nsteps_;
  std::vector<CombiDataType> phi_;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    //ar& boost::serialization::make_nvp( BOOST_PP_STRINGIZE(*this),boost::serialization::base_object<Task>(*this));
    ar& boost::serialization::base_object<Task>(*this);
    ar& dt_;
    ar& nsteps_;
  }
};

BOOST_CLASS_EXPORT(TaskAdvectionFDM)
BOOST_CLASS_EXPORT(StaticFaults)
BOOST_CLASS_EXPORT(WeibullFaults)

BOOST_CLASS_EXPORT(FaultCriterion)
void checkManager(bool useCombine, bool useFG, double l0err, double l2err, bool useCombineAsync=false) {
  int size = useFG ? 2 : 7;
  BOOST_REQUIRE(TestHelper::checkNumProcs(size));

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) return;

  size_t ngroup = useFG ? 1 : 6;
  size_t nprocs = 1;
  theMPISystem()->initWorld(comm, ngroup, nprocs);

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    LoadModel* loadmodel = new LinearLoadModel();

    DimType dim = 2;
    LevelVector lmin(dim, useFG ? 6 : 3);
    LevelVector lmax(dim, 6), leval(dim, 6);

    // choose dt according to CFL condition
    combigrid::real dt = 0.0001;

    size_t nsteps = 100;
    size_t ncombi = 100;
    std::vector<bool> boundary(dim, true);

    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createAdaptiveCombischeme();
    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskAdvectionFDM(levels[i], boundary, coeffs[i],
                                     loadmodel, dt, nsteps);
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params);

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    /* distribute task according to load model and start computation for
     * the first time */
    manager.runfirst();

    for (size_t it = 0; it < ncombi; ++it) {
      if (useCombine) {
        if(useCombineAsync) {
          manager.combineAsync();
          } else {
            manager.combine();
          }
        }

      manager.runnext();
    }

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
    printf("Error: %f", fg_exact.getlpNorm(0));
    printf("Error2: %f", fg_exact.getlpNorm(2));
    // results recorded previously
    BOOST_CHECK(TestHelper::equals(fg_exact.getlpNorm(0), l0err, 1e-5));
    BOOST_CHECK(TestHelper::equals(fg_exact.getlpNorm(2), l2err, 1e-5));

    manager.exit();
  } else {
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    while (signal != EXIT)
      signal = pgroup.wait();
  }
}

BOOST_AUTO_TEST_SUITE(manager)

BOOST_AUTO_TEST_CASE(test_1) {
  // use recombination
  checkManager(true, false, 1.54369, 11.28857);
}

BOOST_AUTO_TEST_CASE(test_2) {
  // don't use recombination
  checkManager(false, false, 1.65104, 12.46828);
}

BOOST_AUTO_TEST_CASE(test_3) {
  // calculate solution on fullgrid
  checkManager(false, true, 1.51188, 10.97143);
}

BOOST_AUTO_TEST_CASE(test_4) {
  // use recombination
  checkManager(true, false, 1.54369, 11.28857, false);
}

BOOST_AUTO_TEST_SUITE_END()
