#define BOOST_TEST_DYN_LINK
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <complex>
#include <cstdarg>
#include <iostream>
#include <random>
#include <vector>

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

#include "test_helper.hpp"

using namespace combigrid;

/**
 * functor for test function $f(x) = \sum_{i=0}^d x_i * (i+1)$
 * which maps to points on a hyperplane
 */
class TestFn {
 public:
  // function value
  std::complex<double> operator()(std::vector<double>& coords) {
    std::complex<double> result(1, 0);
    for (size_t d = 0; d < coords.size(); ++d) {
      result += coords[d] * (double)(d + 1);
    }
    return result;
  }
};

void checkDistributedFullgrid(LevelVector& levels, IndexVector& procs, std::vector<bool>& boundary,
                              int size, bool forward = false) {
  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) return;

  TestFn f;
  const DimType dim = levels.size();

  // create dfg
  DistributedFullGrid<std::complex<double>> dfg(dim, levels, comm, boundary, procs, forward);

  IndexType nrElements = 1;
  for (DimType d = 0; d < dim; ++d) {
    nrElements *= (1 << levels[d]) + (boundary[d] ? 1 : -1);
  }
  BOOST_CHECK(nrElements == dfg.getNrElements());

  // set function values
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);
    dfg.getData()[li] = f(coords);
  }

  // test addToUniformSG, extractFromUniformSG
  LevelVector lmin = levels;
  LevelVector lmax = levels;
  for (DimType d = 0; d < dim; ++d) {
    lmax[d] *= 2;
  }
  DistributedSparseGridUniform<std::complex<double>> dsg(dim, lmax, lmin, boundary, comm);
  dfg.addToUniformSG(dsg, 2.1);
  DistributedFullGrid<std::complex<double>> dfg2(dim, levels, comm, boundary, procs, forward);
  dfg2.extractFromUniformSG(dsg);

  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);
    BOOST_TEST(2.1 * dfg.getData()[li] == dfg2.getData()[li]);
  }

  // test gatherFullgrid
  FullGrid<std::complex<double>> fg(dim, levels, boundary);
  dfg.gatherFullGrid(fg, 0);

  // only check on rank 0
  if (TestHelper::getRank(comm) != 0) {
    return;
  }

  for (size_t i = 0; i < static_cast<size_t>(fg.getNrElements()); ++i) {
    std::vector<double> coords(dim);
    fg.getCoords(i, coords);
    BOOST_TEST(fg.getData()[i] == f(coords));
  }
}

BOOST_AUTO_TEST_SUITE(distributedfullgrid)

// with boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_1) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(6));
  LevelVector levels = {2, 2};
  IndexVector procs = {2, 3};
  std::vector<bool> boundary(2, true);
  checkDistributedFullgrid(levels, procs, boundary, 6);
}
BOOST_AUTO_TEST_CASE(test_2) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  checkDistributedFullgrid(levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_3) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  checkDistributedFullgrid(levels, procs, boundary, 8, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_4) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {2, 3};
  IndexVector procs = {2, 2};
  std::vector<bool> boundary(2, true);
  checkDistributedFullgrid(levels, procs, boundary, 4);
}
BOOST_AUTO_TEST_CASE(test_5) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  IndexVector procs = {1, 4, 2};
  std::vector<bool> boundary(3, true);
  checkDistributedFullgrid(levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_6) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  IndexVector procs = {1, 4, 2};
  std::vector<bool> boundary(3, true);
  checkDistributedFullgrid(levels, procs, boundary, 8, true);
}

// without boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_7) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(2));
  LevelVector levels = {3, 3};
  IndexVector procs = {2, 1};
  std::vector<bool> boundary(2, false);
  checkDistributedFullgrid(levels, procs, boundary, 2);
}
BOOST_AUTO_TEST_CASE(test_8) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  checkDistributedFullgrid(levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_9) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  checkDistributedFullgrid(levels, procs, boundary, 8, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_10) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {2, 3};
  IndexVector procs = {2, 2};
  std::vector<bool> boundary(2, false);
  checkDistributedFullgrid(levels, procs, boundary, 4);
}
BOOST_AUTO_TEST_CASE(test_11) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  checkDistributedFullgrid(levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_12) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  checkDistributedFullgrid(levels, procs, boundary, 8, true);
}

// partial boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_13) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {4, 4};
  IndexVector procs = {2, 2};
  std::vector<bool> boundary(2, false);
  boundary[1] = true;
  checkDistributedFullgrid(levels, procs, boundary, 4);
}
BOOST_AUTO_TEST_CASE(test_14) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  boundary[0] = true;
  checkDistributedFullgrid(levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_15) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  boundary[0] = true;
  checkDistributedFullgrid(levels, procs, boundary, 8, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_16) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {2, 3};
  IndexVector procs = {2, 2};
  std::vector<bool> boundary(2, false);
  boundary[0] = true;
  checkDistributedFullgrid(levels, procs, boundary, 4);
}
BOOST_AUTO_TEST_CASE(test_17) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 4, 3};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  boundary[2] = true;
  checkDistributedFullgrid(levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_18) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 4, 3};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  boundary[2] = true;
  checkDistributedFullgrid(levels, procs, boundary, 8, true);
}

BOOST_AUTO_TEST_SUITE_END()
