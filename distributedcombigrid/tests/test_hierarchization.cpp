#define BOOST_TEST_DYN_LINK
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <complex>
#include <cstdarg>
#include <iostream>
#include <typeinfo>
#include <vector>

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGridEnsemble.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"
#include "sgpp/distributedcombigrid/hierarchization/Hierarchization.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

#include "test_helper.hpp"

/**
 * functor for test function $f(x) = \prod_{i=0}^d x_i^2$
 * with boundary
 */
class TestFn_1 {
  LevelVector levels_;

 public:
  TestFn_1(LevelVector& levels) : levels_(levels) {}

  // overload for function value
  std::complex<double> operator()(std::vector<double>& coords) {
    std::complex<double> result(1, 0);
    for (size_t d = 0; d < coords.size(); ++d) {
      result.real(result.real() * coords[d] * coords[d]);
    }
    return result;
  }

  // overload for hierarchical surpluses
  std::complex<double> operator()(IndexVector& index) {
    std::complex<double> result(1, 0);
    for (size_t d = 0; d < index.size(); ++d) {
      std::vector<double> coeff = {0.0, 1.0};
      for (int l = 1; l < levels_[d] + 1; ++l) {
        for (int i = 0; i < 1 << l; ++i) {
          if (i % 2) {
            coeff.insert(coeff.begin() + i, -std::pow(2.0, -2 * l));
          }
        }
      }
      result.real(result.real() * coeff[index[d]]);
    }
    return result;
  }
};

/**
 * functor for test function $f(x) = \prod_{i=0}^d x_i * (x_i - 1)$
 * without boundary
 */
class TestFn_2 {
  LevelVector levels_;

 public:
  TestFn_2(LevelVector& levels) : levels_(levels) {}

  // overload for function value
  std::complex<double> operator()(std::vector<double>& coords) {
    std::complex<double> result(1, 0);
    for (size_t d = 0; d < coords.size(); ++d) {
      result.real(result.real() * coords[d] * (coords[d] - 1.0));
    }
    return result;
  }

  // overload for hierarchical surpluses
  std::complex<double> operator()(IndexVector& index) {
    std::complex<double> result(1, 0);
    for (size_t d = 0; d < index.size(); ++d) {
      std::vector<double> coeff = {0.25};
      for (int l = 1; l < levels_[d] + 1; ++l) {
        for (int i = 0; i < 1 << l; ++i) {
          if (i % 2 == 0) {
            coeff.insert(coeff.begin() + i, -std::pow(2.0, -2 * l));
          }
        }
      }
      result.real(result.real() * coeff[index[d]]);
    }
    return result;
  }
};

/**
 * functor for test function $f(x) = \prod_{i=0}^d x_i * (x_i - 1)$
 * without boundary
 */
class TestFn_3 {
  LevelVector levels_;

 public:
  TestFn_3(LevelVector& levels) : levels_(levels) {}

  // overload for function value
  std::complex<double> operator()(std::vector<double>& coords) {
    std::complex<double> result(1, 1);
    for (size_t d = 0; d < coords.size(); ++d) {
      result.real(result.real() * coords[d] * (coords[d] - 1.0));
      result.imag(result.imag() * coords[d] * coords[d]);
    }
    return result;
  }

  // overload for hierarchical surpluses
  std::complex<double> operator()(IndexVector& index) {
    std::complex<double> result(1, 1);
    for (size_t d = 0; d < index.size(); ++d) {
      std::vector<double> coeff = {0, 0};
      std::vector<double> coeff2 = {0.0, 1.0};
      for (int l = 1; l < levels_[d] + 1; ++l) {
        for (int i = 0; i < 1 << l; ++i) {
          if (i % 2) {
            coeff.insert(coeff.begin() + i, -std::pow(2.0, -2 * l));
            coeff2.insert(coeff2.begin() + i, -std::pow(2.0, -2 * l));
          }
        }
      }
      result.real(result.real() * coeff[index[d]]);
      result.imag(result.imag() * coeff2[index[d]]);
    }
    return result;
  }
};

template <typename Functor>
void checkBiorthogonalHierarchization(Functor& f, DistributedFullGrid<std::complex<double>>& dfg,
                                      bool checkValues = true) {
  real formerL1 = 0.;
  if (checkValues) {
    // calculate l1 integral of actual data
    formerL1 = dfg.getLpNorm(1);
  }

  auto dim = dfg.getDimension();
  std::vector<bool> hierarchizationDimensions(dim, true);
  DistributedHierarchization::hierarchizeBiorthogonal<DistributedFullGrid<std::complex<double>>>(
      dfg, hierarchizationDimensions);

  // now, all of the mass should be on the coarsest level -> the corners
  if (checkValues) {
    const auto innerIntegral = dfg.getInnerNodalBasisFunctionIntegral();
    // calculate discrete analytical l1
    double analyticalL1 = 0.;
    for (IndexType gli = 0; gli < dfg.getNrElements(); ++gli) {
      auto isOnBoundary = dfg.isGlobalLinearIndexOnBoundary(gli);
      auto countBoundary = std::count(isOnBoundary.begin(), isOnBoundary.end(), true);
      double hatFcnIntegral = oneOverPowOfTwo[countBoundary];
      std::vector<double> coords(dim);
      dfg.getCoordsGlobal(gli, coords);
      analyticalL1 += std::abs(f(coords)) * hatFcnIntegral;
    }
    analyticalL1 *= innerIntegral;
    BOOST_TEST(formerL1 == analyticalL1, boost::test_tools::tolerance(TestHelper::tolerance));
    auto cornersValues = dfg.getCornersValues();
    BOOST_CHECK(cornersValues.size() == powerOfTwo[dim]);
    auto sumOfCornerValues =
        std::accumulate(cornersValues.begin(), cornersValues.end(), std::complex<double>(0.),
                        std::plus<std::complex<double>>());
    auto currentL1 = std::abs(sumOfCornerValues * oneOverPowOfTwo[dim]);
    BOOST_TEST(currentL1 == formerL1, boost::test_tools::tolerance(TestHelper::tolerance));
  }

  DistributedHierarchization::dehierarchizeBiorthogonal<DistributedFullGrid<std::complex<double>>>(
      dfg, hierarchizationDimensions);
}

template <typename Functor>
void checkFullWeightingHierarchization(Functor& f, DistributedFullGrid<std::complex<double>>& dfg,
                                       bool checkValues = true) {
  real formerL1 = 0.;
  if (checkValues) {
    // calculate l1 integral of actual data
    formerL1 = dfg.getLpNorm(1);
  }

  auto dim = dfg.getDimension();
  std::vector<bool> hierarchizationDimensions(dim, true);
  DistributedHierarchization::hierarchizeFullWeighting<DistributedFullGrid<std::complex<double>>>(
      dfg, hierarchizationDimensions);

  // now, all of the mass should be on the coarsest level -> the corners
  if (checkValues) {
    auto cornersValues = dfg.getCornersValues();
    BOOST_CHECK(cornersValues.size() == powerOfTwo[dim]);

    auto sumOfCornerValues =
        std::accumulate(cornersValues.begin(), cornersValues.end(), std::complex<double>(0.),
                        std::plus<std::complex<double>>());
    auto currentL1 = std::abs(sumOfCornerValues * oneOverPowOfTwo[dim]);
    BOOST_TEST(currentL1 == formerL1, boost::test_tools::tolerance(TestHelper::tolerance));
  }

  DistributedHierarchization::dehierarchizeFullWeighting<DistributedFullGrid<std::complex<double>>>(
      dfg, hierarchizationDimensions);
}

template <typename Functor>
void checkHierarchization(Functor& f, LevelVector& levels, IndexVector& procs,
                          std::vector<bool>& boundary, int size, bool forward = false,
                          bool checkValues = true) {

  if (std::all_of(boundary.begin(), boundary.end(), [](bool i){return i;})) {
    checkHierarchizationEnsemble(f, levels, procs, boundary, size, forward);
  }

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm != MPI_COMM_NULL) {
    const DimType dim = levels.size();
    DistributedFullGrid<std::complex<double>> dfg(dim, levels, comm, boundary, procs, forward);
    // run test with value check
    checkHierarchization(f, dfg, true);
  }
}

template <typename Functor>
void checkHierarchization(Functor& f, DistributedFullGrid<std::complex<double>>& dfg,
                          bool checkValues = true) {
  CommunicatorType comm = dfg.getCommunicator();
  const DimType dim = dfg.getDimension();
  auto boundary = dfg.returnBoundaryFlags();

  FullGrid<std::complex<double>> fg(dim, dfg.getLevels(), boundary);

  if (checkValues) {
    // fill distributed fg with test function
    for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
      std::vector<double> coords(dim);
      dfg.getCoordsLocal(li, coords);
      dfg.getData()[li] = f(coords);
    }

    // create fg and fill with test function
    fg.createFullGrid();
    for (size_t i = 0; i < static_cast<size_t>(fg.getNrElements()); ++i) {
      std::vector<double> coords(dim);
      fg.getCoords(i, coords);
      fg.getData()[i] = f(coords);
    }

    // hierarchize fg
    Hierarchization::hierarchize(fg);
  }

  // hierarchize distributed fg
  DistributedHierarchization::hierarchize(dfg);

  if (checkValues) {
    // compare hierarchical surpluses
    for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
      IndexType gi = dfg.getGlobalLinearIndex(li);
      IndexVector axisIndex(dim);
      fg.getVectorIndex(gi, axisIndex);

      // compare fg and distributed fg
      BOOST_TEST(dfg.getData()[li] == fg.getData()[gi],
                 boost::test_tools::tolerance(TestHelper::tolerance));
      // compare distributed fg to exact solution
      BOOST_TEST(dfg.getData()[li] == f(axisIndex),
                 boost::test_tools::tolerance(TestHelper::tolerance));
    }

    // dehiarchize fg
    Hierarchization::dehierarchize(fg);
  }

  // dehierarchize distributed fg
  DistributedHierarchization::dehierarchize(dfg);

  if (checkValues) {
    // compare function values
    for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
      IndexType gi = dfg.getGlobalLinearIndex(li);

      std::vector<double> coords_fg(dim);
      fg.getCoords(gi, coords_fg);

      std::vector<double> coords_dfg(dim);
      dfg.getCoordsLocal(li, coords_dfg);

      BOOST_CHECK(coords_dfg == coords_fg);

      // compare fg and distributed fg
      BOOST_TEST(dfg.getData()[li] == fg.getData()[gi],
                 boost::test_tools::tolerance(TestHelper::tolerance));
      // compare distributed fg and fg to exact solution
      BOOST_TEST(dfg.getData()[li] == f(coords_fg),
                 boost::test_tools::tolerance(TestHelper::tolerance));
    }
  }

  // call this so that tests are also run for mass-conserving bases
  // but only for boundary grids and in case we are not measuring time
  if (checkValues && std::all_of(boundary.begin(), boundary.end(), [](bool b) { return b; })) {
    if (!(typeid(Functor) == typeid(TestFn_3))) {
      // TODO figure out what is supposed to happen for true complex numbers,
      // currently std::abs does not seem to do the right thing
      // create distributed fg and copy values
      DistributedFullGrid<std::complex<double>> dfgCopyOne(
          dim, dfg.getLevels(), dfg.getCommunicator(), dfg.returnBoundaryFlags(),
          dfg.getParallelization(), true, dfg.getDecomposition());
      DistributedFullGrid<std::complex<double>> dfgCopyTwo(
          dim, dfg.getLevels(), dfg.getCommunicator(), dfg.returnBoundaryFlags(),
          dfg.getParallelization(), true, dfg.getDecomposition());
      for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
        dfgCopyOne.getData()[li] = dfg.getData()[li];
        dfgCopyTwo.getData()[li] = dfg.getData()[li];
      }
      checkFullWeightingHierarchization(f, dfgCopyOne, true);
      checkBiorthogonalHierarchization(f, dfgCopyTwo, true);
      // afterwards, the values should be the same again
      // compare dfgCopy to former values
      for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
        BOOST_TEST(dfgCopyOne.getData()[li] == dfg.getData()[li],
                   boost::test_tools::tolerance(TestHelper::tolerance));
        BOOST_TEST(dfgCopyTwo.getData()[li] == dfg.getData()[li],
                   boost::test_tools::tolerance(TestHelper::tolerance));
      }
    }
  }
}

template <typename Functor>
void checkHierarchizationEnsemble(Functor& f, LevelVector& levels, IndexVector& procs,
                          std::vector<bool>& boundary, int size, bool forward = false) {
  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) return;

  const DimType dim = levels.size();

  // create distributed fg and fill with test function
  DFGEnsemble dfgEnsemble(dim, levels, comm, boundary, procs, forward);
  for (size_t i = 0; i < dfgEnsemble.getNumFullGrids(); ++i) {
    auto& dfg = dfgEnsemble.getDFG(i);
    // real offset = std::pow(GridEnumeration::getNumberOfHigherDimensions(i), dim);
    real offset = GridEnumeration::getNumberOfHigherDimensions(i);
    for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
      std::vector<double> coords(dim);
      dfg.getCoordsLocal(li, coords);
      dfg.getData()[li] = f(coords).real() + offset;
    }
  }

  HierarchicalHatBasisFunction basisFctn;
  std::vector<BasisFunctionBasis*> bases(dim, &basisFctn);

  for (DimType d = 0; d < dim; ++d) {
    std::vector<bool> hierarchizationDims(dim, false);
    hierarchizationDims[d] = true;
    // hierarchize in dimension d
    DistributedHierarchization::hierarchizeHierachicalBasis<DFGEnsemble,
                                                              HierarchicalHatBasisFunction>(
        dfgEnsemble, hierarchizationDims);

    // compare hierarchical surpluses
    // compare distributed fgs pairwise in each dimension
    // the "higher" grid should always have a coefficient 1. higher
    // (cf how we set the offset above)
    auto lowerIndices = dfgEnsemble.getIndicesOfLowerGridsInDimension(d);
    for (auto& index : lowerIndices) {
      auto& dfgLower = dfgEnsemble.getDFG(index);
      auto& dfgHigher = dfgEnsemble.getDFG(GridEnumeration::getUpperNeighborInDimension(static_cast<char>(d), index));
      real offset = GridEnumeration::getNumberOfHigherDimensions(index);
      for (IndexType li = 0; li < dfgLower.getNrLocalElements(); ++li) {
        BOOST_TEST(dfgLower.getData()[li] + 1. == dfgHigher.getData()[li] , boost::test_tools::tolerance(TestHelper::tolerance));
      }
    }

    // dehierarchize
    DistributedHierarchization::dehierarchizeHierachicalBasis<DFGEnsemble,
                                                              HierarchicalHatBasisFunction>(
        dfgEnsemble, hierarchizationDims);
  }

  // compare function values
  for (size_t i = 0; i < dfgEnsemble.getNumFullGrids(); ++i) {
  // for (size_t i = 0; i < 1; ++i) {
    auto& dfg = dfgEnsemble.getDFG(i);
    real offset = GridEnumeration::getNumberOfHigherDimensions(i);
    for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
      std::vector<double> coords_dfg(dim);
      dfg.getCoordsLocal(li, coords_dfg);
      // compare distributed fg to exact solution
      BOOST_TEST(dfg.getData()[li] == f(coords_dfg).real() + offset, boost::test_tools::tolerance(TestHelper::tolerance));
    }
  }
}

BOOST_FIXTURE_TEST_SUITE(hierarchization, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(240))

// with boundary
// isotropic

// the most basic case with a single worker
BOOST_AUTO_TEST_CASE(test_0) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {1, 1, 1};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 1);
}

// 2D case with 4 workers
BOOST_AUTO_TEST_CASE(test_05) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {2,2,};
  IndexVector procs = {2,2};
  std::vector<bool> boundary(2, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 4);
}

BOOST_AUTO_TEST_CASE(test_1) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_2) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8, true);
}
BOOST_AUTO_TEST_CASE(test_3) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {6, 6, 6};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_4) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {1, 4, 2};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_5) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {1, 1, 8};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_6) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3, 3};
  IndexVector procs = {1, 2, 2, 2};
  std::vector<bool> boundary(4, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_7) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  IndexVector procs = {3, 3, 1};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 9);
}
BOOST_AUTO_TEST_CASE(test_8) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  IndexVector procs = {3, 3, 1};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 9, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_9) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_10) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8, true);
}
BOOST_AUTO_TEST_CASE(test_11) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  IndexVector procs = {2, 1, 4};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_12) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {1, 4, 4};
  IndexVector procs = {1, 2, 4};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_13) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 1, 3, 3, 2};
  IndexVector procs = {2, 1, 2, 2, 1};
  std::vector<bool> boundary(5, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_14) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  IndexVector procs = {3, 3, 1};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 9);
}
BOOST_AUTO_TEST_CASE(test_15) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  IndexVector procs = {3, 3, 1};
  std::vector<bool> boundary(3, true);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 9, true);
}

// without boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_16) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_17) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8, true);
}
BOOST_AUTO_TEST_CASE(test_18) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {6, 6, 6};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_19) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {4, 2, 1};
  std::vector<bool> boundary(3, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_20) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {3, 3, 1};
  std::vector<bool> boundary(3, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 9, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_21) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_22) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8, true);
}
BOOST_AUTO_TEST_CASE(test_23) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  IndexVector procs = {2, 1, 4};
  std::vector<bool> boundary(3, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_24) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {1, 4, 4};
  IndexVector procs = {1, 2, 4};
  std::vector<bool> boundary(3, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_25) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 1, 3, 3, 2};
  IndexVector procs = {2, 1, 2, 2, 1};
  std::vector<bool> boundary(5, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_26) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  IndexVector procs = {3, 3, 1};
  std::vector<bool> boundary(3, false);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 9, true);
}

// with boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_27) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_28) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8, true);
}
BOOST_AUTO_TEST_CASE(test_29) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {6, 6, 6};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_30) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {1, 4, 2};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_31) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  IndexVector procs = {1, 1, 8};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_32) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3, 3};
  IndexVector procs = {1, 2, 2, 2};
  std::vector<bool> boundary(4, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_33) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  IndexVector procs = {3, 3, 1};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 9);
}
BOOST_AUTO_TEST_CASE(test_34) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  IndexVector procs = {3, 3, 1};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 9, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_35) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_36) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  IndexVector procs = {2, 2, 2};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8, true);
}
BOOST_AUTO_TEST_CASE(test_37) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  IndexVector procs = {2, 1, 4};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_38) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {1, 4, 4};
  IndexVector procs = {1, 2, 4};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_39) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 1, 3, 3, 2};
  IndexVector procs = {2, 1, 2, 2, 1};
  std::vector<bool> boundary(5, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 8);
}
BOOST_AUTO_TEST_CASE(test_40) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  IndexVector procs = {3, 3, 1};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 9);
}
BOOST_AUTO_TEST_CASE(test_41) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  IndexVector procs = {3, 3, 1};
  std::vector<bool> boundary(3, true);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, 9, true);
}
BOOST_AUTO_TEST_CASE(test_42) {
  // large test case with timing
  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  CommunicatorType comm = TestHelper::getComm(8);
  if (comm != MPI_COMM_NULL) {
    LevelVector levels = {11, 11, 4};
    IndexVector procs = {2, 2, 2};
    std::vector<bool> boundary(3, true);
    auto forward = true;
    TestFn_1 testFn(levels);
    DistributedFullGrid<std::complex<double>> dfg(3, levels, comm, boundary, procs, forward);
    auto start = std::chrono::high_resolution_clock::now();
    checkHierarchization(testFn, dfg, false);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("hat hierarchization time: " << duration.count() << " milliseconds");
  }
  // on ipvs-epyc2@  : 600 milliseconds w single msgs
}

BOOST_AUTO_TEST_CASE(test_43) {
  // large test case with timing for full weighting
  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  CommunicatorType comm = TestHelper::getComm(8);
  if (comm != MPI_COMM_NULL) {
    LevelVector levels = {11, 11, 4};
    IndexVector procs = {2, 2, 2};
    std::vector<bool> boundary(3, true);
    auto forward = true;
    TestFn_1 testFn(levels);
    DistributedFullGrid<std::complex<double>> dfg(3, levels, comm, boundary, procs, forward);
    auto start = std::chrono::high_resolution_clock::now();
    checkFullWeightingHierarchization(testFn, dfg, false);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("full weighting hierarchization time: " << duration.count() << " milliseconds");
  }
  // on ipvs-epyc2@  : 3300 milliseconds w single msgs
}
BOOST_AUTO_TEST_CASE(test_44) {
  // large test case with timing for full weighting
  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  CommunicatorType comm = TestHelper::getComm(8);
  if (comm != MPI_COMM_NULL) {
    LevelVector levels = {11, 11, 4};
    IndexVector procs = {2, 2, 2};
    std::vector<bool> boundary(3, true);
    auto forward = true;
    TestFn_1 testFn(levels);
    DistributedFullGrid<std::complex<double>> dfg(3, levels, comm, boundary, procs, forward);
    auto start = std::chrono::high_resolution_clock::now();
    checkBiorthogonalHierarchization(testFn, dfg, false);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("biorthogonal hierarchization time: " << duration.count() << " milliseconds");
  }
  // on ipvs-epyc2@  : 3100 milliseconds w single msgs
}

BOOST_AUTO_TEST_SUITE_END()
