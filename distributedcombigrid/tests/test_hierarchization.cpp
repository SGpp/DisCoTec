#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>


#include <boost/test/unit_test.hpp>
#include <complex>
#include <cstdarg>
#include <iostream>
#include <typeinfo>
#include <vector>

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"
#include "sgpp/distributedcombigrid/hierarchization/Hierarchization.hpp"
#include "sgpp/distributedcombigrid/utils/MonteCarlo.hpp"
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

template <typename FG_ELEMENT>
real getMonteCarloMass(DistributedFullGrid<FG_ELEMENT>& dfg, size_t npoints) {
  BOOST_TEST_CHECKPOINT("start mass calculation");
  auto dim = dfg.getDimension();
  auto interpolationCoords = montecarlo::getRandomCoordinates(npoints, dim);
  auto interpolatedValues = dfg.getInterpolatedValues(interpolationCoords);
  real mass = 0.;
  for (size_t i = 0; i < npoints; ++i) {
    // auto scalarCoordinate = std::accumulate(
    //     interpolationCoords[i].begin(), interpolationCoords[i].end(), 1., std::multiplies<real>());
    // TODO what about complex' imaginary part?
    mass += std::real(interpolatedValues[i]);
  }
  mass = mass / npoints;
  BOOST_TEST_CHECKPOINT("end mass calculation");
  return mass;
}

template <typename FG_ELEMENT>
std::vector<real> getMonteCarloMomenta(DistributedFullGrid<FG_ELEMENT>& dfg, size_t npoints) {
  BOOST_TEST_CHECKPOINT("start momentum calculation");
  const auto dim = dfg.getDimension();
  auto interpolationCoords = montecarlo::getRandomCoordinates(npoints, dim);
  auto interpolatedValues = dfg.getInterpolatedValues(interpolationCoords);
  std::vector<real> momenta(dim + 1, 0.);
  for (size_t i = 0; i < npoints; ++i) {
    for (DimType d = 0; d < dim; ++d) {
      momenta[d] += interpolationCoords[i][d] * std::real(interpolatedValues[i]);
    }
    auto scalarCoordinate = std::accumulate(
        interpolationCoords[i].begin(), interpolationCoords[i].end(), 1., std::multiplies<real>());
    // TODO what about complex' imaginary part?
    momenta[dim] += scalarCoordinate * std::real(interpolatedValues[i]);
  }
  for (auto& momentum : momenta) {
    momentum = momentum / npoints;
  }
  BOOST_TEST_CHECKPOINT("end momentum calculation");
  return momenta;
}

template <typename FG_ELEMENT>
using FunctionPointer = void (*)(DistributedFullGrid<FG_ELEMENT>& dfg,
                                 const std::vector<bool>& dims);

template <typename FG_ELEMENT>
real checkConservationOfMomentum(DistributedFullGrid<FG_ELEMENT>& dfg,
                                 FunctionPointer<FG_ELEMENT> hierarchizationOperator) {
  const auto& procs = dfg.getParallelization();
  const auto& boundary = dfg.returnBoundaryFlags();
  const auto& comm = dfg.getCommunicator();
  size_t nPointsMonteCarlo = 1e6;
  BOOST_CHECK(std::all_of(boundary.begin(), boundary.end(), [](bool b) { return b == 2; }));
  real mcMassBefore = getMonteCarloMass(dfg, nPointsMonteCarlo);
  auto mcMomentaBefore = getMonteCarloMomenta(dfg, nPointsMonteCarlo);
  real mcMomentumBefore = mcMomentaBefore.back();

  BOOST_TEST_CHECKPOINT("begin hierarchization");
  auto dim = dfg.getDimension();
  std::vector<bool> hierarchizationDimensions(dim, true);
  hierarchizationOperator(dfg, hierarchizationDimensions);
  BOOST_TEST_CHECKPOINT("end hierarchization");

  // now, all of the momentum should be on the coarsest level -> the corners
  // register with dsg and extract very small dfg from it (with only corners)
  BOOST_TEST_CHECKPOINT("create sparse grid");
  LevelVector lmin(dim, 0);  // TODO
  LevelVector lone(dim, 1);  // cannot use lmin 0 in dsgu's constructor
  auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<FG_ELEMENT>>(
      new DistributedSparseGridUniform<FG_ELEMENT>(dim, dfg.getLevels(), lone, comm));
  dfg.registerUniformSG(*uniDSG);
  // TODO also cannot use level 0 to register dfg -- problem!
  auto dfgOne = std::unique_ptr<DistributedFullGrid<FG_ELEMENT>>(
      new DistributedFullGrid<FG_ELEMENT>(dim, lone, comm, boundary, procs));
  dfgOne->registerUniformSG(*uniDSG);
  BOOST_TEST_CHECKPOINT("registered sparse grid");

  dfg.addToUniformSG(*uniDSG, 1.);
  dfgOne->extractFromUniformSG(*uniDSG);

  // TODO extract boundary grid lvl 1 to lvl 0 for now
  //  mostly stolen from dfg::getCornersValues
  auto corners = dfgOne->getCornersGlobalVectorIndices();
  std::sort(corners.begin(), corners.end());
  std::vector<FG_ELEMENT> values;
  std::vector<IndexType> localCornerIndices;
  for (size_t cornerNo = 0; cornerNo < corners.size(); ++cornerNo) {
    if (dfgOne->isGlobalIndexHere(corners[cornerNo])) {
      // convert to local vector index, then to linear index
      IndexVector locAxisIndex(dfgOne->getDimension());
      bool present = dfgOne->getLocalVectorIndex(corners[cornerNo], locAxisIndex);
      BOOST_CHECK(present);
      auto index = dfgOne->getLocalLinearIndex(locAxisIndex);
      localCornerIndices.push_back(index);
    }
  }
  // make sure corner values are in right order
  std::sort(localCornerIndices.begin(), localCornerIndices.end());
  for (const auto& index : localCornerIndices) {
    values.push_back(dfgOne->getElementVector()[index]);
  }

  auto dfgZero = std::unique_ptr<DistributedFullGrid<FG_ELEMENT>>(
      new DistributedFullGrid<FG_ELEMENT>(dim, lmin, comm, boundary, procs));
  BOOST_CHECK(values.size() == dfgZero->getElementVector().size());
  dfgZero->getElementVector() = values;

  // no need to dehierarchize, is nodal/scaling function on coarsest grid anyways

  BOOST_TEST_CHECKPOINT("added and extracted from sparse grid");

  real mcMassAfter = getMonteCarloMass(*dfgZero, nPointsMonteCarlo);
  auto mcMomentaAfter = getMonteCarloMomenta(*dfgZero, nPointsMonteCarlo);
  real mcMomentumAfter = mcMomentaAfter.back();

  BOOST_TEST(mcMassAfter == mcMassBefore, boost::test_tools::tolerance(5e-2));

  for (DimType d = 0; d < dim + 1; ++d) {
    std::cout << d << std::endl;
    BOOST_TEST(mcMomentaAfter[d] == mcMomentaBefore[d], boost::test_tools::tolerance(5e-2));
  }

  // std::cout << dfg << std::endl;
  // std::cout << std::endl;
  // std::cout << *dfgZero << std::endl;
  return mcMomentumAfter;
}

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
  DistributedHierarchization::hierarchizeBiorthogonal<std::complex<double>>(
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
    BOOST_CHECK(cornersValues.size() == static_cast<size_t>(powerOfTwo[dim]));
    auto sumOfCornerValues =
        std::accumulate(cornersValues.begin(), cornersValues.end(), std::complex<double>(0.),
                        std::plus<std::complex<double>>());
    auto currentL1 = std::abs(sumOfCornerValues * oneOverPowOfTwo[dim]);
    BOOST_TEST(currentL1 == formerL1, boost::test_tools::tolerance(TestHelper::tolerance));
  }

  DistributedHierarchization::dehierarchizeBiorthogonal<std::complex<double>>(
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
  DistributedHierarchization::hierarchizeFullWeighting<std::complex<double>>(
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

  DistributedHierarchization::dehierarchizeFullWeighting<std::complex<double>>(
      dfg, hierarchizationDimensions);
}

template <typename Functor>
void checkHierarchization(Functor& f, LevelVector& levels, std::vector<int>& procs,
                          std::vector<BoundaryType>& boundary, bool forward = false,
                          bool checkValues = true) {
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    const auto dim = static_cast<DimType>(levels.size());
    DistributedFullGrid<std::complex<double>> dfg(dim, levels, comm, boundary, procs, forward);
    // run test with value check
    checkHierarchization(f, dfg, checkValues);
  }
}

template <typename FG_ELEMENT>
void fillDFGrandom(DistributedFullGrid<FG_ELEMENT>& dfg, real&& a, real&& b) {
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    dfg.getData()[li] = static_cast<FG_ELEMENT>(
        montecarlo::getRandomNumber(std::forward<real>(a), std::forward<real>(b)));
  }
}

template <typename Functor, typename FG_ELEMENT>
void fillDFGbyFunction(Functor& f, DistributedFullGrid<FG_ELEMENT>& dfg) {
  const DimType dim = dfg.getDimension();
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);
    dfg.getData()[li] = f(coords);
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
    fillDFGbyFunction(f, dfg);

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
  if (checkValues && std::all_of(boundary.begin(), boundary.end(), [](bool b) { return b == 2; })) {
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

  BOOST_CHECK(!TestHelper::testStrayMessages(comm));
}

BOOST_FIXTURE_TEST_SUITE(hierarchization, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(240))

// with boundary
// isotropic

// the most basic case with a single worker
BOOST_AUTO_TEST_CASE(test_0) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 1, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}

// 2D case with 4 workers
BOOST_AUTO_TEST_CASE(test_05) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {
      2,
      2,
  };
  std::vector<int> procs = {2, 2};
  std::vector<BoundaryType> boundary(2, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}

BOOST_AUTO_TEST_CASE(test_1) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_2) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}
BOOST_AUTO_TEST_CASE(test_3) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {6, 6, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_4) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 4, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_5) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 1, 8};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_6) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3, 3};
  std::vector<int> procs = {1, 2, 2, 2};
  std::vector<BoundaryType> boundary(4, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_7) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_8) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_9) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_10) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}
BOOST_AUTO_TEST_CASE(test_11) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 1, 4};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_12) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {1, 4, 4};
  std::vector<int> procs = {1, 2, 4};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_13) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 1, 3, 3, 2};
  std::vector<int> procs = {2, 1, 2, 2, 1};
  std::vector<BoundaryType> boundary(5, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_14) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_15) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}

// without boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_16) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_17) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}
BOOST_AUTO_TEST_CASE(test_18) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {6, 6, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_19) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {4, 2, 1};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_20) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_21) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_22) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}
BOOST_AUTO_TEST_CASE(test_23) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {2, 1, 4};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_24) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {1, 4, 4};
  std::vector<int> procs = {1, 2, 4};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_25) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 1, 3, 3, 2};
  std::vector<int> procs = {2, 1, 2, 2, 1};
  std::vector<BoundaryType> boundary(5, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_26) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}

// with boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_27) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_28) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}
BOOST_AUTO_TEST_CASE(test_29) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {6, 6, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_30) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 4, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_31) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 1, 8};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_32) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3, 3};
  std::vector<int> procs = {1, 2, 2, 2};
  std::vector<BoundaryType> boundary(4, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_33) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_34) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_35) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_36) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}
BOOST_AUTO_TEST_CASE(test_37) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 1, 4};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_38) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {1, 4, 4};
  std::vector<int> procs = {1, 2, 4};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_39) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 1, 3, 3, 2};
  std::vector<int> procs = {2, 1, 2, 2, 1};
  std::vector<BoundaryType> boundary(5, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_40) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_41) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  checkHierarchization(testFn, levels, procs, boundary, true);
}

// these large tests only make sense when assertions are not checked (takes too long otherwise)
#ifdef NDEBUG

BOOST_AUTO_TEST_CASE(test_42) {
  // large test case with timing
  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  std::vector<int> procs = {2, 2, 2};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    LevelVector levels = {11, 11, 4};
    std::vector<BoundaryType> boundary(3, 2);
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
  std::vector<int> procs = {2, 2, 2};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    LevelVector levels = {11, 11, 4};
    std::vector<BoundaryType> boundary(3, 2);
    auto forward = true;
    TestFn_1 testFn(levels);
    DistributedFullGrid<std::complex<double>> dfg(3, levels, comm, boundary, procs, forward);
    auto start = std::chrono::high_resolution_clock::now();
    checkFullWeightingHierarchization(testFn, dfg, false);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("full weighting hierarchization time: " << duration.count()
                                                               << " milliseconds");
  }
  // on ipvs-epyc2@  : 3300 milliseconds w single msgs
}
BOOST_AUTO_TEST_CASE(test_44) {
  // large test case with timing for full weighting
  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  std::vector<int> procs = {2, 2, 2};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    LevelVector levels = {11, 11, 4};
    std::vector<BoundaryType> boundary(3, 2);
    auto forward = true;
    TestFn_1 testFn(levels);
    DistributedFullGrid<std::complex<double>> dfg(3, levels, comm, boundary, procs, forward);
    auto start = std::chrono::high_resolution_clock::now();
    checkBiorthogonalHierarchization(testFn, dfg, false);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("biorthogonal hierarchization time: " << duration.count()
                                                             << " milliseconds");
  }
  // on ipvs-epyc2@  : 3100 milliseconds w single msgs
}

#endif  // def NDEBUG

BOOST_AUTO_TEST_CASE(momentum) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  DimType dim = 3;
  std::vector<int> procs(dim, 1);
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    LevelVector levels(dim, 4);
    std::vector<BoundaryType> boundary(dim, 2);
    auto forward = true;  // todo toggle
    TestFn_1 testFn(levels);
    auto constFctn = [](const std::vector<double>& coords) { return 11.; };
    //TODO make non-periodic operators momentum conserving by methodology from Frank Koster's diss
    // std::cout << "test fcn" << std::endl;
    // DistributedFullGrid<std::complex<double>> dfg(3, levels, comm, boundary, procs, forward);
    // {
    //   // initialize dfg with function
    //   fillDFGbyFunction(testFn, dfg);
    //   checkConservationOfMomentum<std::complex<double>>(
    //       dfg, DistributedHierarchization::hierarchizeBiorthogonal<std::complex<double>>);
    // }
    // {
    //   fillDFGbyFunction(testFn, dfg);
    //   checkConservationOfMomentum<std::complex<double>>(
    //       dfg, DistributedHierarchization::hierarchizeFullWeighting<std::complex<double>>);
    // }
    // std::cout << "random" << std::endl;
    // DistributedFullGrid<real> dfg2(3, levels, comm, boundary, procs, forward);
    // {
    //   // initialize dfg with random numbers
    //   fillDFGrandom(dfg2, -100., 100.);
    //   checkConservationOfMomentum<real>(
    //       dfg2, DistributedHierarchization::hierarchizeBiorthogonal<real>);
    // }
    // {
    //   fillDFGrandom(dfg2, -100., 100.);
    //   checkConservationOfMomentum<real>(
    //       dfg2, DistributedHierarchization::hierarchizeFullWeighting<real>);
    // }
    // std::cout << "const" << std::endl;
    DistributedFullGrid<real> dfg3(3, levels, comm, boundary, procs, forward);
    {
      // initialize dfg with constant function
      fillDFGbyFunction(constFctn, dfg3);
      // usually, momentum will not be conserved for periodic BC, but for constant functions it
      // should work
      auto momentum = checkConservationOfMomentum<real>(
          dfg3, DistributedHierarchization::hierarchizeBiorthogonalPeriodic<real>);
    }
    {
      fillDFGbyFunction(constFctn, dfg3);
      auto momentum = checkConservationOfMomentum<real>(
          dfg3, DistributedHierarchization::hierarchizeFullWeightingPeriodic<real>);
      BOOST_TEST(momentum == std::pow(0.5, dim) * 11., boost::test_tools::tolerance(5e-2));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
