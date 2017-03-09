#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mpi.h>
#include <iostream>
#include <complex>
#include <cstdarg>
#include <vector>

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/hierarchization/Hierarchization.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"

/**
 * functor for test function
 */
class TestFn {
  LevelVector levels_;
public:
  TestFn(LevelVector& levels) : levels_(levels) {}

  // overload for function value
  std::complex<double> operator()(std::vector<double>& coords) {
    std::complex<double> result(1,0);
    for (size_t d = 0; d < coords.size(); ++d) {
      result.real(result.real() * coords[d] * coords[d]);
    }
    return result;
  }

  // overload for hierarchical surpluses
  std::complex<double> operator()(IndexVector& index) {
    std::complex<double> result(1,0);
    for (size_t d = 0; d < index.size(); ++d) {
      std::vector<double> coeff = {0.0, 1.0};
      for (int l = 1; l < levels_[d]+1; ++l) {
        for (int i = 0; i < 1<<l; ++i) {
          if (i % 2) {
            coeff.insert(coeff.begin()+i, -std::pow(2.0, -2*l));
          }
        }
      }
      result.real(result.real() * coeff[index[d]]);
    }
    return result;
  }
};

bool checkNumProcs(int nprocs) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size == nprocs;
}

template<typename Functor>
void checkHierarchization(Functor& f, LevelVector& levels, IndexVector& procs,
                          std::vector<bool>& boundary) {
  CommunicatorType comm = MPI_COMM_WORLD;
  const DimType dim = levels.size();

  // create distributed fg and fill with test function
  DistributedFullGrid<std::complex<double>> dfg(dim, levels, comm, boundary,
                                                procs);
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);
    dfg.getData()[li] = f(coords);
  }

  // create fg and fill with test function
  FullGrid<std::complex<double>> fg(dim, levels, boundary);
  fg.createFullGrid();
  for (size_t i = 0; i < static_cast<size_t>(fg.getNrElements()); ++i) {
    std::vector<double> coords(dim);
    fg.getCoords(i, coords);
    fg.getData()[i] = f(coords);
  }

  // hierarchize fg and distributed fg
  Hierarchization::hierarchize(fg);
  DistributedHierarchization::hierarchize(dfg);

  // compare hierarchical surpluses
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    IndexType gi = dfg.getGlobalLinearIndex(li);
    IndexVector axisIndex(dim);
    fg.getVectorIndex(gi, axisIndex);

    // compare fg and distributed fg
    BOOST_CHECK(std::abs(dfg.getData()[li] - fg.getData()[gi]) <= 1e-12);

    // compare distributed fg and fg to exact solution
    auto sol = f(axisIndex);
    BOOST_CHECK(std::abs(fg.getData()[gi] - sol) <= 1e-12);
    BOOST_CHECK(std::abs(dfg.getData()[li] - sol) <= 1e-12);
  }

  // dehiarchize fg and distributed fg
  Hierarchization::dehierarchize(fg);
  DistributedHierarchization::dehierarchize(dfg);

  // compare function values
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    IndexType gi = dfg.getGlobalLinearIndex(li);

    std::vector<double> coords_fg(dim);
    fg.getCoords(gi, coords_fg);

    std::vector<double> coords_dfg(dim);
    dfg.getCoordsLocal(li, coords_dfg);

    BOOST_CHECK(coords_dfg == coords_fg);

    // compare fg and distributed fg
    BOOST_CHECK(std::abs(dfg.getData()[li] - fg.getData()[gi]) <= 1e-12);

    // compare distributed fg and fg to exact solution
    auto sol = f(coords_fg);
    BOOST_CHECK(std::abs(fg.getData()[gi] - sol) <= 1e-12);
    BOOST_CHECK(std::abs(dfg.getData()[li] - sol) <= 1e-12);
  }
}

BOOST_AUTO_TEST_CASE(IsotropicHierarchizationWithBoundary1) {
  BOOST_REQUIRE(checkNumProcs(8));

  LevelVector levels = {4,4,4};
  IndexVector procs = {2,2,2};
  std::vector<bool> boundary(3, true);
  TestFn testFn(levels);

  checkHierarchization(testFn, levels, procs, boundary);
}

BOOST_AUTO_TEST_CASE(IsotropicHierarchizationWithBoundary2) {
  BOOST_REQUIRE(checkNumProcs(8));

  LevelVector levels = {6,6,6};
  IndexVector procs = {2,2,2};
  std::vector<bool> boundary(3, true);
  TestFn testFn(levels);

  checkHierarchization(testFn, levels, procs, boundary);
}

BOOST_AUTO_TEST_CASE(IsotropicHierarchizationWithBoundary3) {
  BOOST_REQUIRE(checkNumProcs(8));

  LevelVector levels = {4,4,4};
  IndexVector procs = {1,4,2};
  std::vector<bool> boundary(3, true);
  TestFn testFn(levels);

  checkHierarchization(testFn, levels, procs, boundary);
}

BOOST_AUTO_TEST_CASE(IsotropicHierarchizationWithBoundary4) {
  BOOST_REQUIRE(checkNumProcs(8));

  LevelVector levels = {4,4,4};
  IndexVector procs = {1,1,8};
  std::vector<bool> boundary(3, true);
  TestFn testFn(levels);

  checkHierarchization(testFn, levels, procs, boundary);
}

BOOST_AUTO_TEST_CASE(AnisotropicHierarchizationWithBoundary1) {
  BOOST_REQUIRE(checkNumProcs(8));

  LevelVector levels = {2,4,6};
  IndexVector procs = {2,2,2};
  std::vector<bool> boundary(3, true);
  TestFn testFn(levels);

  checkHierarchization(testFn, levels, procs, boundary);
}

BOOST_AUTO_TEST_CASE(AnisotropicHierarchizationWithBoundary2) {
  BOOST_REQUIRE(checkNumProcs(8));

  LevelVector levels = {2,4,6};
  IndexVector procs = {2,1,4};
  std::vector<bool> boundary(3, true);
  TestFn testFn(levels);

  checkHierarchization(testFn, levels, procs, boundary);
}
