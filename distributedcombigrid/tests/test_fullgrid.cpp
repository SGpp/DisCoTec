#define BOOST_TEST_DYN_LINK
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <complex>
#include <cstdarg>
#include <iostream>
#include <random>
#include <vector>

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

#include "test_helper.hpp"

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

void checkFullgrid(LevelVector& levels, std::vector<bool>& boundary) {
  CommunicatorType comm = TestHelper::getComm(1);
  if (comm == MPI_COMM_NULL) return;

  TestFn f;
  const DimType dim = levels.size();

  // create fg
  FullGrid<std::complex<double>> fg(dim, levels, boundary);
  fg.createFullGrid();
  BOOST_CHECK(fg.isGridCreated());

  IndexType nrElements = 1;
  for (DimType d = 0; d < dim; ++d) {
    nrElements *= (1 << levels[d]) + (boundary[d] ? 1 : -1);
  }
  BOOST_CHECK(nrElements == fg.getNrElements());

  // set function values
  for (IndexType i = 0; i < fg.getNrElements(); ++i) {
    std::vector<double> coords(dim);
    fg.getCoords(i, coords);
    fg.getData()[i] = f(coords);
  }

  // test eval
  // always use the same seed for random values
  std::mt19937 gen(0);
  std::uniform_real_distribution<> dis(0, 1);
  for (size_t i = 0; i < 100; ++i) {
    std::vector<double> coords(dim);
    for (DimType d = 0; d < dim; ++d) {
      coords[d] = dis(gen);
    }
    // every evaluated point should be the same as the function value,
    // as the funtion is linear and fullgrid uses linear basis functions.
  BOOST_CHECK_SMALL(abs(fg.eval(coords) - f(coords)), TestHelper::tolerance);
  }

  // some corner cases for eval
  std::vector<double> coords(dim);
  for (DimType d = 0; d < dim; ++d) {
    coords[d] = 0;
  }
  BOOST_CHECK_SMALL(abs(fg.eval(coords) - f(coords)), TestHelper::tolerance);
  for (DimType d = 0; d < dim; ++d) {
    coords[d] = 1;
  }
  BOOST_CHECK_SMALL(abs(fg.eval(coords) - f(coords)), TestHelper::tolerance);

  // test add
  FullGrid<std::complex<double>> fg2(dim, levels, boundary);
  fg2.createFullGrid();
  fg2.add(fg, 2.1);

  LevelVector levels2 = levels;
  levels2[0] += 1;
  FullGrid<std::complex<double>> fg3(dim, levels2, boundary);
  fg3.createFullGrid();
  fg3.add(fg, 4.2);

  for (size_t i = 0; i < 100; ++i) {
    std::vector<double> coords(dim);
    for (DimType d = 0; d < dim; ++d) {
      coords[d] = dis(gen);

    }
    BOOST_CHECK_SMALL(abs(2.1 * fg.eval(coords) - fg2.eval(coords)), TestHelper::tolerance);
    BOOST_CHECK_SMALL(abs(4.2 * fg.eval(coords) - fg3.eval(coords)), TestHelper::tolerance);
  }
}

BOOST_AUTO_TEST_SUITE(fullgrid)

// with boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_1) {
  LevelVector levels = {2, 2};
  std::vector<bool> boundary(2, true);
  checkFullgrid(levels, boundary);
}
BOOST_AUTO_TEST_CASE(test_2) {
  LevelVector levels = {3, 3, 3};
  std::vector<bool> boundary(3, true);
  checkFullgrid(levels, boundary);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_3) {
  LevelVector levels = {2, 3};
  std::vector<bool> boundary(2, true);
  checkFullgrid(levels, boundary);
}
BOOST_AUTO_TEST_CASE(test_4) {
  LevelVector levels = {2, 4, 6};
  std::vector<bool> boundary(3, true);
  checkFullgrid(levels, boundary);
}

// without boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_5) {
  LevelVector levels = {3, 3};
  std::vector<bool> boundary(2, false);
  checkFullgrid(levels, boundary);
}
BOOST_AUTO_TEST_CASE(test_6) {
  LevelVector levels = {4, 4, 4};
  std::vector<bool> boundary(3, false);
  checkFullgrid(levels, boundary);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_7) {
  LevelVector levels = {2, 3};
  std::vector<bool> boundary(2, false);
  checkFullgrid(levels, boundary);
}
BOOST_AUTO_TEST_CASE(test_8) {
  LevelVector levels = {2, 3, 4};
  std::vector<bool> boundary(3, false);
  checkFullgrid(levels, boundary);
}

// partial boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_9) {
  LevelVector levels = {4, 4};
  std::vector<bool> boundary(2, false);
  boundary[1] = true;
  checkFullgrid(levels, boundary);
}
BOOST_AUTO_TEST_CASE(test_10) {
  LevelVector levels = {3, 3, 3};
  std::vector<bool> boundary(3, false);
  boundary[0] = true;
  checkFullgrid(levels, boundary);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_11) {
  LevelVector levels = {2, 3};
  std::vector<bool> boundary(2, false);
  boundary[0] = true;
  checkFullgrid(levels, boundary);
}
BOOST_AUTO_TEST_CASE(test_12) {
  LevelVector levels = {3, 4, 3};
  std::vector<bool> boundary(3, false);
  boundary[2] = true;
  checkFullgrid(levels, boundary);
}

BOOST_AUTO_TEST_SUITE_END()
