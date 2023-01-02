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

#include "fullgrid/DistributedFullGrid.hpp"
#include "fullgrid/FullGrid.hpp"
#include "hierarchization/DistributedHierarchization.hpp"
#include "hierarchization/Hierarchization.hpp"
#include "utils/MonteCarlo.hpp"
#include "utils/Types.hpp"
#include "test_helper.hpp"
#include "TaskConstParaboloid.hpp"

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
    //     interpolationCoords[i].begin(), interpolationCoords[i].end(), 1.,
    //     std::multiplies<real>());
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
                                  const std::vector<bool>& dims, LevelVector lmin);

template <typename FG_ELEMENT>
real checkConservationOfMomentum(DistributedFullGrid<FG_ELEMENT>& dfg,
                                 FunctionPointer<FG_ELEMENT> hierarchizationOperator) {
  const auto& procs = dfg.getParallelization();
  BOOST_CHECK(procs.size() == dfg.getDimension());
  const auto& boundary = dfg.returnBoundaryFlags();
  BOOST_CHECK(boundary.size() == dfg.getDimension());
  const auto& comm = dfg.getCommunicator();
  size_t nPointsMonteCarlo = 1e6;
  BOOST_CHECK(std::all_of(boundary.begin(), boundary.end(), [](BoundaryType b) { return b == 2; }));
  real mcMassBefore = getMonteCarloMass(dfg, nPointsMonteCarlo);
  auto mcMomentaBefore = getMonteCarloMomenta(dfg, nPointsMonteCarlo);
  real mcMomentumBefore = mcMomentaBefore.back();

  BOOST_TEST_CHECKPOINT("begin hierarchization");
  auto dim = dfg.getDimension();
  std::vector<bool> hierarchizationDimensions(dim, true);
  hierarchizationOperator(dfg, hierarchizationDimensions, LevelVector(dim, 0));
  BOOST_TEST_CHECKPOINT("end hierarchization");

  // now, all of the momentum should be on the coarsest level -> the corners
  // register with dsg and extract very small dfg from it (with only corners)
  BOOST_TEST_CHECKPOINT("create sparse grid");
  LevelVector lmin(dim, 0);  // TODO
  LevelVector lone(dim, 1);  // cannot use lmin 0 in dsgu's constructor
  auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<FG_ELEMENT>>(
      new DistributedSparseGridUniform<FG_ELEMENT>(dim, dfg.getLevels(), lone, comm));
  uniDSG->registerDistributedFullGrid(dfg);
  // TODO also cannot use level 0 to register dfg -- problem!
  auto dfgOne = std::unique_ptr<DistributedFullGrid<FG_ELEMENT>>(
      new DistributedFullGrid<FG_ELEMENT>(dim, lone, comm, boundary, procs));
  uniDSG->registerDistributedFullGrid(*dfgOne);
  uniDSG->setZero();
  BOOST_TEST_CHECKPOINT("registered sparse grid");

  uniDSG->addDistributedFullGrid(dfg, 1.);
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
    // std::cout << d << std::endl;
    BOOST_TEST(mcMomentaAfter[d] == mcMomentaBefore[d], boost::test_tools::tolerance(5e-2));
  }

  // std::cout << dfg << std::endl;
  // std::cout << std::endl;
  // std::cout << *dfgZero << std::endl;
  return mcMomentumAfter;
}

template <typename Functor>
void checkBiorthogonalHierarchization(Functor& f, DistributedFullGrid<std::complex<double>>& dfg,
                                      bool checkValues = true, LevelVector lmin = LevelVector(0)) {
  real formerL1 = 0.;
  if (checkValues) {
    // calculate l1 integral of actual data
    formerL1 = dfg.getLpNorm(1);
  }

  auto dim = dfg.getDimension();
  std::vector<bool> hierarchizationDimensions(dim, true);
  DistributedHierarchization::hierarchizeBiorthogonal<std::complex<double>>(
      dfg, hierarchizationDimensions, lmin);

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
    if (lmin == LevelVector(dim, 0)) {
      auto sumOfCornerValues =
          std::accumulate(cornersValues.begin(), cornersValues.end(), std::complex<double>(0.),
                          std::plus<std::complex<double>>());
      auto currentL1 = std::abs(sumOfCornerValues * oneOverPowOfTwo[dim]);
      BOOST_TEST(currentL1 == formerL1, boost::test_tools::tolerance(TestHelper::tolerance));
    }
  }

  DistributedHierarchization::dehierarchizeBiorthogonal<std::complex<double>>(
      dfg, hierarchizationDimensions, lmin);
}

template <typename Functor>
void checkFullWeightingHierarchization(Functor& f, DistributedFullGrid<std::complex<double>>& dfg,
                                       bool checkValues = true, LevelVector lmin = LevelVector(0)) {
  real formerL1 = 0.;
  if (checkValues) {
    // calculate l1 integral of actual data
    formerL1 = dfg.getLpNorm(1);
  }

  auto dim = dfg.getDimension();
  std::vector<bool> hierarchizationDimensions(dim, true);
  DistributedHierarchization::hierarchizeFullWeighting<std::complex<double>>(
      dfg, hierarchizationDimensions, lmin);

  // now, all of the mass should be on the coarsest level -> the corners
  // but only if we hierarchize all the way down
  if (checkValues && lmin == LevelVector(dim, 0)) {
    auto cornersValues = dfg.getCornersValues();
    BOOST_CHECK(cornersValues.size() == powerOfTwo[dim]);

    auto sumOfCornerValues =
        std::accumulate(cornersValues.begin(), cornersValues.end(), std::complex<double>(0.),
                        std::plus<std::complex<double>>());
    auto currentL1 = std::abs(sumOfCornerValues * oneOverPowOfTwo[dim]);
    BOOST_TEST(currentL1 == formerL1, boost::test_tools::tolerance(TestHelper::tolerance));
  }

  DistributedHierarchization::dehierarchizeFullWeighting<std::complex<double>>(
      dfg, hierarchizationDimensions, lmin);
}

template <typename Functor>
void checkHierarchization(Functor& f, LevelVector& levels, std::vector<int>& procs,
                          std::vector<BoundaryType>& boundary, bool forward = false,
                          bool checkValues = true, LevelVector lmin = LevelVector(0)) {
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    const auto dim = static_cast<DimType>(levels.size());
    DistributedFullGrid<std::complex<double>> dfg(dim, levels, comm, boundary, procs, forward);
    // run test with value check
    checkHierarchization(f, dfg, checkValues, lmin);
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
                          bool checkValues = true, LevelVector lmin = LevelVector(0)) {
  CommunicatorType comm = dfg.getCommunicator();
  const DimType dim = dfg.getDimension();
  auto boundary = dfg.returnBoundaryFlags();
  const auto& dfgLevel = dfg.getLevels();
  std::vector<bool> hierarchizationDimensions(dim, true);

  auto nonDistributedBoundary = boundary;
  // non-distributed (de)hierarchization not adapted to one-sided boundary (yet?)
  bool anyOneSidedBoundary =
      std::any_of(boundary.begin(), boundary.end(), [](const BoundaryType& b) { return b == 1; });
  if (anyOneSidedBoundary) {
    nonDistributedBoundary = std::vector<BoundaryType>(dim, 2);
  }
  FullGrid<std::complex<double>> fg(dim, dfgLevel, nonDistributedBoundary);
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
  }
  // have copy of non-hierarchized FG
  auto fgNodal = fg;
  if (checkValues) {
    // hierarchize fg
    BOOST_TEST_CHECKPOINT("Non-Distributed Hierarchization begins");
    Hierarchization::hierarchize(fg);
  }

  BOOST_TEST_CHECKPOINT("Distributed Hierarchization begins");
  // hierarchize distributed fg

  DistributedHierarchization::hierarchizeHierachicalBasis<std::complex<double>>(
      dfg, hierarchizationDimensions, lmin);

  if (checkValues) {
    if (lmin.size() > 0) {
      for (size_t i = 0; i < dfgLevel.size(); ++i) {
        BOOST_ASSERT(lmin[i] <= dfgLevel[i]);
      }
      LevelVector levels_of_point(dim);
      IndexVector tmp(dim);
      IndexVector axisIndex(dim);
      for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
        IndexType gi = dfg.getGlobalLinearIndex(li);
        dfg.getGlobalLI(gi, levels_of_point, tmp);

        if (levels_of_point > lmin) {
          if (!anyOneSidedBoundary) {
            fg.getVectorIndex(gi, axisIndex);
            // the finest levels up to lmin should always be hierarchized
            // compare hierarchical surpluses
            // compare fg and distributed fg
            BOOST_TEST(dfg.getData()[li] == fg.getData()[gi],
                       boost::test_tools::tolerance(TestHelper::tolerance));
            // compare distributed fg to exact solution
          }
          BOOST_TEST(dfg.getData()[li] == f(axisIndex),
                     boost::test_tools::tolerance(TestHelper::tolerance));
        } else if (levels_of_point <= lmin) {
          // the coarsest levels should not be hierarchized at all
          fg.getVectorIndex(gi, axisIndex);
          // compare non-hierarchized fg and distributed fg
          BOOST_TEST(dfg.getData()[li] == fgNodal.getData()[gi],
                     boost::test_tools::tolerance(TestHelper::tolerance));
        }
      }
    } else {
      // compare hierarchical surpluses
      for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
        IndexType gi = dfg.getGlobalLinearIndex(li);
        IndexVector axisIndex(dim), localAxisIndex(dim);
        dfg.getLocalVectorIndex(li, localAxisIndex);
        dfg.getGlobalVectorIndex(localAxisIndex, axisIndex);
        BOOST_REQUIRE_EQUAL(dfg.getGlobalLinearIndex(axisIndex), gi);
        if (!anyOneSidedBoundary) {
          auto fgAxisIndex = axisIndex;
          fg.getVectorIndex(gi, fgAxisIndex);
          BOOST_REQUIRE(axisIndex == fgAxisIndex);
          // compare fg and distributed fg
          BOOST_TEST(dfg.getData()[li] == fg.getData()[gi],
                     boost::test_tools::tolerance(TestHelper::tolerance));
        }
        // compare distributed fg to exact solution
        BOOST_TEST(dfg.getData()[li] == f(axisIndex),
                   boost::test_tools::tolerance(TestHelper::tolerance));
      }
    }

    // dehierarchize fg
    BOOST_TEST_CHECKPOINT("Non-distributed Dehierarchization begins");
    Hierarchization::dehierarchize(fg);
  }

  // dehierarchize distributed fg
  BOOST_TEST_CHECKPOINT("Distributed Dehierarchization begins");
  DistributedHierarchization::dehierarchizeHierachicalBasis<std::complex<double>>(
      dfg, hierarchizationDimensions, lmin);

  if (checkValues) {
    // compare function values
    for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
      BOOST_TEST_CHECKPOINT("Check values");
      BOOST_TEST_CONTEXT(std::to_string(li));
      IndexType gi = dfg.getGlobalLinearIndex(li);

      std::vector<double> coords_dfg(dim);
      dfg.getCoordsLocal(li, coords_dfg);

      if (!anyOneSidedBoundary) {
        std::vector<double> coords_fg(dim);
        fg.getCoords(gi, coords_fg);
        BOOST_CHECK(coords_dfg == coords_fg);
        // compare fg and distributed fg
        BOOST_TEST(dfg.getData()[li] == fg.getData()[gi],
                   boost::test_tools::tolerance(TestHelper::tolerance));
      }
      // compare distributed fg and fg to exact solution
      BOOST_TEST(dfg.getData()[li] == f(coords_dfg),
                 boost::test_tools::tolerance(TestHelper::tolerance));
    }
  }

  // call this so that tests are also run for mass-conserving bases
  // but only for boundary grids and in case we are not measuring time
  if (checkValues &&
      std::all_of(boundary.begin(), boundary.end(), [](BoundaryType b) { return b == 2; })) {
    if (!(typeid(Functor) == typeid(TestFn_3))) {
      // TODO figure out what is supposed to happen for true complex numbers,
      // currently std::abs does not seem to do the right thing
      // create distributed fg and copy values
      DistributedFullGrid<std::complex<double>> dfgCopyOne(
          dim, dfgLevel, dfg.getCommunicator(), dfg.returnBoundaryFlags(),
          dfg.getParallelization(), true, dfg.getDecomposition());
      DistributedFullGrid<std::complex<double>> dfgCopyTwo(
          dim, dfgLevel, dfg.getCommunicator(), dfg.returnBoundaryFlags(),
          dfg.getParallelization(), true, dfg.getDecomposition());
      for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
        dfgCopyOne.getData()[li] = dfg.getData()[li];
        dfgCopyTwo.getData()[li] = dfg.getData()[li];
      }
      checkFullWeightingHierarchization(f, dfgCopyOne, true, lmin);
      checkBiorthogonalHierarchization(f, dfgCopyTwo, true, lmin);
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

void checkHierarchizationParaboloid(LevelVector& levels, std::vector<int>& procs,
                                    std::vector<BoundaryType>& boundary, bool forward = false,
                                    bool checkValues = true, LevelVector lmin = LevelVector(0)) {
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    const auto dim = static_cast<DimType>(levels.size());
    DistributedFullGrid<std::complex<double>> dfg(dim, levels, comm, boundary, procs, forward);
    auto f = ParaboloidFn<std::complex<double>>(&dfg);
    // run test with value check
    checkHierarchization<decltype(f)>(f, dfg, checkValues, lmin);
  }
}

template <typename FG_ELEMENT>
IndexVector checkExtentsOfDFG(const DistributedFullGrid<FG_ELEMENT>& dfg) {
  // check that all MPI ranks agree on the extents of the grid
  auto extents = dfg.getGlobalSizes();
  auto extents_max = extents;
  auto extents_min = extents;
  int dim = static_cast<int>(extents.size());
  MPI_Datatype indexDType = getMPIDatatype(abstraction::getabstractionDataType<IndexType>());
  MPI_Allreduce(MPI_IN_PLACE, extents_max.data(), dim, indexDType, MPI_MAX, dfg.getCommunicator());
  MPI_Allreduce(MPI_IN_PLACE, extents_min.data(), dim, indexDType, MPI_MIN, dfg.getCommunicator());
  BOOST_CHECK_EQUAL_COLLECTIONS(extents.begin(), extents.end(), extents_max.begin(),
                                extents_max.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(extents.begin(), extents.end(), extents_min.begin(),
                                extents_min.end());
  return extents;
}

BOOST_FIXTURE_TEST_SUITE(hierarchization, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(240))

BOOST_AUTO_TEST_CASE(test_exchangeData1d, *boost::unit_test::timeout(20)) {
  for (auto procs : std::vector<std::vector<int>>{
           {1, 4, 1, 2, 1, 1}, {1, 8, 1, 1, 1, 1}, {1, 2, 1, 2, 2, 1}}) {
    auto dimensionality = static_cast<DimType>(procs.size());
    LevelVector levels = {1, 10, 1, 6, 2, 1};
    LevelVector lzero(dimensionality, 0);
    LevelVector lhalf = levels;
    std::transform(lhalf.begin(), lhalf.end(), lhalf.begin(), [](int i) { return i / 2; });
    for (auto lmin : std::vector<LevelVector>{lzero, lhalf, levels}) {
      for (DimType d = 0; d < dimensionality; ++d) {
        std::vector<std::vector<IndexType>> remoteKeysHierarchization(3);
        std::vector<std::vector<IndexType>> remoteKeysDehierarchization(3);
        bool isOnLowerBoundaryInD = false;
        bool isOnUpperBoundaryInD = false;
        for (BoundaryType b : std::vector<BoundaryType>({0, 2, 1})) {
          std::vector<BoundaryType> boundary(dimensionality, b);
          CommunicatorType comm =
              TestHelper::getComm(procs, std::vector<int>(dimensionality, b == 1 ? 1 : 0));
          if (comm != MPI_COMM_NULL) {
            BOOST_TEST_CHECKPOINT("Testing dimension " << d << " with boundary " << b
                                                       << " and lmin " << lmin[d]);
            DistributedFullGrid<std::complex<double>> dfg(dimensionality, levels, comm, boundary,
                                                          procs, false);
            {
              isOnLowerBoundaryInD = dfg.getCartesianUtils().isOnLowerBoundaryInDimension(d);
              std::vector<int> processLocation;
              dfg.getCartesianUtils().getPartitionCoordsOfLocalRank(processLocation);
              isOnUpperBoundaryInD = processLocation[d] == procs[d] - 1;
            }
            auto extents = checkExtentsOfDFG(dfg);

            // exchange data
            std::vector<RemoteDataContainer<std::complex<double>>> remoteDataHierarchization;
            BOOST_CHECK_NO_THROW(exchangeData1d(dfg, d, remoteDataHierarchization, lmin[d]));
            BOOST_CHECK((remoteDataHierarchization.size() == 0) || (procs[d] > 1));
            for (const auto& r : remoteDataHierarchization) {
              auto receiveKeyIndex = r.getKeyIndex();
              BOOST_CHECK_LT(receiveKeyIndex, extents[d]);
              remoteKeysHierarchization[b].push_back(receiveKeyIndex);
            }
            BOOST_CHECK(std::is_sorted(remoteKeysHierarchization[b].begin(),
                                       remoteKeysHierarchization[b].end()));

            std::vector<RemoteDataContainer<std::complex<double>>> remoteDataDehierarchization;
            BOOST_CHECK_NO_THROW(
                exchangeData1dDehierarchization(dfg, d, remoteDataDehierarchization, lmin[d]));
            BOOST_CHECK((remoteDataDehierarchization.size() == 0) || (procs[d] > 1));
            // more data may need to be exchanged for dehierarchization, but never less
            BOOST_CHECK_GE(remoteDataDehierarchization.size(), remoteDataHierarchization.size());
            for (const auto& r : remoteDataDehierarchization) {
              auto receiveKeyIndex = r.getKeyIndex();
              BOOST_CHECK_LT(receiveKeyIndex, extents[d]);
              remoteKeysDehierarchization[b].push_back(receiveKeyIndex);
            }
            BOOST_CHECK(std::is_sorted(remoteKeysDehierarchization[b].begin(),
                                       remoteKeysDehierarchization[b].end()));

            std::vector<RemoteDataContainer<std::complex<double>>> remoteDataAll;
            BOOST_CHECK_NO_THROW(exchangeAllData1d(dfg, d, remoteDataAll));
            BOOST_CHECK((procs[d] == 1 && remoteDataAll.size() == 0) ||
                        (procs[d] > 1 && remoteDataAll.size() > 0));
            if (procs[d] > 1) {
              // check that all indices from other ranks along the pole are present
              BOOST_CHECK_EQUAL(remoteDataAll.size(),
                                dfg.getGlobalSizes()[d] - dfg.getLocalSizes()[d]);
            }
            MPI_Barrier(comm);
            TestHelper::testStrayMessages(comm);
          }
        }
        if (remoteKeysDehierarchization[2].empty()) {
          // either the rank did not take part in the test
          if (TestHelper::getRank(MPI_COMM_WORLD) <
              std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<int>())) {
            if (procs[d] == 1) {
              // or this dimension was not distributed
              BOOST_CHECK(isOnLowerBoundaryInD);
              BOOST_CHECK(isOnUpperBoundaryInD);
            } else if (remoteKeysHierarchization[2].empty()) {
              // or the partition is selected such that nothing needs to be communicated
              BOOST_CHECK(true);
            } else {
              // or lmin was too large
              BOOST_CHECK_EQUAL(lmin[d], levels[d]);
            }
          }
        } else {
          IndexType upperBoundaryIndex = powerOfTwo[levels[d]];
          // compare remoteKeysHierarchization and remoteKeysDehierarchization for different
          // boundary types

          for (const auto& remoteKeys : {remoteKeysHierarchization, remoteKeysDehierarchization}) {
            // check that the indices without boundary are the same ones as with one boundary,
            // except for the 0
            auto remoteKeysExpectedForZero = remoteKeys[1];
            if (remoteKeysExpectedForZero.front() == 0) {
              remoteKeysExpectedForZero.erase(remoteKeysExpectedForZero.begin());
            }
            // shift because that the lower boundary index is not included in the 0 boundary case
            for (auto& key : remoteKeysExpectedForZero) {
              key -= 1;
            }
            BOOST_CHECK_EQUAL_COLLECTIONS(remoteKeys[0].begin(), remoteKeys[0].end(),
                                          remoteKeysExpectedForZero.begin(),
                                          remoteKeysExpectedForZero.end());

            // check that the same indices are exchanged for 1 and 2 sided boundaries
            auto remoteKeysExpectedForOne = remoteKeys[2];
            // but where the upper boundary is expected in the two-boundary case, we need the lower
            // boundary
            if (!remoteKeysExpectedForOne.empty() &&
                remoteKeysExpectedForOne.back() == upperBoundaryIndex) {
              remoteKeysExpectedForOne.pop_back();
              if (remoteKeysExpectedForOne.front() != 0 && !isOnLowerBoundaryInD) {
                remoteKeysExpectedForOne.insert(remoteKeysExpectedForOne.begin(), 0);
              }
            }
            // in case we are at the highest index along the pole, also need index 0
            if (isOnUpperBoundaryInD && !isOnLowerBoundaryInD &&
                (remoteKeysExpectedForOne.empty() || remoteKeysExpectedForOne.front() != 0)) {
              remoteKeysExpectedForOne.insert(remoteKeysExpectedForOne.begin(), 0);
            }
            BOOST_CHECK_EQUAL_COLLECTIONS(remoteKeys[1].begin(), remoteKeys[1].end(),
                                          remoteKeysExpectedForOne.begin(),
                                          remoteKeysExpectedForOne.end());
          }
        }
      }
    }
  }
}

// with boundary
// isotropic

// the most basic case with a single worker
BOOST_AUTO_TEST_CASE(test_minus1, *boost::unit_test::timeout(10)) {
  for (auto d : std::vector<DimType>({1, 2, 3, 4, 5})) {
    BOOST_TEST_CHECKPOINT("Testing dimension " + std::to_string(d));
    // LevelVector levels(d, 1);
    LevelVector levels(d, 2);
    std::vector<int> procs(d, 1);
    std::vector<BoundaryType> boundary(d, 2);

    BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary, false, true));
  }
}
BOOST_AUTO_TEST_CASE(test_0) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 1, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
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
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}

BOOST_AUTO_TEST_CASE(test_1) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_2) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}
BOOST_AUTO_TEST_CASE(test_3) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {6, 6, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_4) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 4, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_5) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 1, 8};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_6) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3, 3};
  std::vector<int> procs = {1, 2, 2, 2};
  std::vector<BoundaryType> boundary(4, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_7) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_8) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_9) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_10) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}
BOOST_AUTO_TEST_CASE(test_11) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 1, 4};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_11_no_hierarchization) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 1, 4};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, false, true, levels));
}
BOOST_AUTO_TEST_CASE(test_11_half_hierarchization) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 1, 4};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(
      checkHierarchization(testFn, levels, procs, boundary, false, true, {1, 2, 3}));
}
BOOST_AUTO_TEST_CASE(test_12) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {1, 4, 4};
  std::vector<int> procs = {1, 2, 4};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_13) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 1, 3, 3, 2};
  std::vector<int> procs = {2, 1, 2, 2, 1};
  std::vector<BoundaryType> boundary(5, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_14) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_14_no_hierarchization) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, false, true, levels));
}
BOOST_AUTO_TEST_CASE(test_14_half_hierarchization) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(
      checkHierarchization(testFn, levels, procs, boundary, false, true, {1, 2, 2}));
}
BOOST_AUTO_TEST_CASE(test_15) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_1 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}

// without boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_16) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_17) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}
BOOST_AUTO_TEST_CASE(test_18) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {6, 6, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_19) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {4, 2, 1};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_20) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_21) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_22) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}
BOOST_AUTO_TEST_CASE(test_23) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {2, 1, 4};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_24) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {1, 4, 4};
  std::vector<int> procs = {1, 2, 4};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_25) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 1, 3, 3, 2};
  std::vector<int> procs = {2, 1, 2, 2, 1};
  std::vector<BoundaryType> boundary(5, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_26) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 0);
  TestFn_2 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}

// with boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_27) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_28) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}
BOOST_AUTO_TEST_CASE(test_29) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {6, 6, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_30) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 4, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_31) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 1, 8};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_32) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3, 3};
  std::vector<int> procs = {1, 2, 2, 2};
  std::vector<BoundaryType> boundary(4, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_33) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_34) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_35) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_36) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}
BOOST_AUTO_TEST_CASE(test_37) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 1, 4};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_38) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {1, 4, 4};
  std::vector<int> procs = {1, 2, 4};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_39) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 1, 3, 3, 2};
  std::vector<int> procs = {2, 1, 2, 2, 1};
  std::vector<BoundaryType> boundary(5, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_40) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_41) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 2);
  TestFn_3 testFn(levels);
  BOOST_CHECK_NO_THROW(checkHierarchization(testFn, levels, procs, boundary, true));
}

// periodic with one-sided boundary
// isotropic
BOOST_AUTO_TEST_CASE(test_p_minus1) {
  LevelVector levels = {1, 1, 1};
  std::vector<int> procs = {1, 1, 1};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_0) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 1, 1};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_05) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {2, 2};
  std::vector<int> procs = {2, 2};
  std::vector<BoundaryType> boundary(2, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_1) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_4) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 4, 2};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_5) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 1, 8};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_6) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3, 3};
  std::vector<int> procs = {1, 2, 2, 2};
  std::vector<BoundaryType> boundary(4, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_7) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_27) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_30) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 4, 2};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_31) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {1, 1, 8};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_32) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3, 3};
  std::vector<int> procs = {1, 2, 2, 2};
  std::vector<BoundaryType> boundary(4, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_33) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_p_35) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_37) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {2, 1, 4};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_38) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {1, 4, 4};
  std::vector<int> procs = {1, 2, 4};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_39) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 1, 3, 3, 2};
  std::vector<int> procs = {2, 1, 2, 2, 1};
  std::vector<BoundaryType> boundary(5, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
}
BOOST_AUTO_TEST_CASE(test_p_40) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {3, 3, 1};
  std::vector<BoundaryType> boundary(3, 1);
  BOOST_CHECK_NO_THROW(checkHierarchizationParaboloid(levels, procs, boundary));
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
    BOOST_CHECK_NO_THROW(checkHierarchization(testFn, dfg, false));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("hat hierarchization time: " << duration.count() << " milliseconds");
  }
  // on ipvs-epyc2@  : 500 milliseconds w single msgs
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
    BOOST_CHECK_NO_THROW(checkFullWeightingHierarchization(testFn, dfg, false));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("full weighting hierarchization time: " << duration.count()
                                                               << " milliseconds");
  }
  // on ipvs-epyc2@  : 2300 milliseconds w single msgs
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
    BOOST_CHECK_NO_THROW(checkBiorthogonalHierarchization(testFn, dfg, false));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("biorthogonal hierarchization time: " << duration.count()
                                                             << " milliseconds");
  }
  // on ipvs-epyc2@  : 2300 milliseconds w single msgs
}

BOOST_AUTO_TEST_CASE(test_p_42) {
  // large test case with timing
  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  std::vector<int> procs = {2, 2, 2};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    LevelVector levels = {11, 11, 4};
    std::vector<BoundaryType> boundary(3, 1);
    auto forward = false;
    TestFn_1 testFn(levels);
    DistributedFullGrid<std::complex<double>> dfg(3, levels, comm, boundary, procs, forward);
    auto start = std::chrono::high_resolution_clock::now();
    BOOST_CHECK_NO_THROW(checkHierarchization(testFn, dfg, false));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("hat hierarchization time: " << duration.count() << " milliseconds");
  }
  // on ipvs-epyc2@  : 480 milliseconds
}

BOOST_AUTO_TEST_CASE(test_p_43) {
  // large test case with timing for full weighting
  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  std::vector<int> procs = {2, 2, 2};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    LevelVector levels = {11, 11, 4};
    std::vector<BoundaryType> boundary(3, 1);
    auto forward = false;
    TestFn_1 testFn(levels);
    DistributedFullGrid<std::complex<double>> dfg(3, levels, comm, boundary, procs, forward);
    auto start = std::chrono::high_resolution_clock::now();
    BOOST_CHECK_NO_THROW(checkFullWeightingHierarchization(testFn, dfg, false));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("full weighting hierarchization time: " << duration.count()
                                                               << " milliseconds");
  }
  // on ipvs-epyc2@  : 2300 milliseconds w single msgs
}
BOOST_AUTO_TEST_CASE(test_p_44) {
  // large test case with timing for full weighting
  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  std::vector<int> procs = {2, 2, 2};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    LevelVector levels = {11, 11, 4};
    std::vector<BoundaryType> boundary(3, 1);
    auto forward = false;
    TestFn_1 testFn(levels);
    DistributedFullGrid<std::complex<double>> dfg(3, levels, comm, boundary, procs, forward);
    auto start = std::chrono::high_resolution_clock::now();
    BOOST_CHECK_NO_THROW(checkBiorthogonalHierarchization(testFn, dfg, false));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("biorthogonal hierarchization time: " << duration.count()
                                                             << " milliseconds");
  }
  // on ipvs-epyc2@  : 2300 milliseconds w single msgs
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
    // TODO make non-periodic operators momentum conserving by methodology from Frank Koster's diss
    //  std::cout << "test fcn" << std::endl;
    //  DistributedFullGrid<std::complex<double>> dfg(3, levels, comm, boundary, procs, forward);
    //  {
    //    // initialize dfg with function
    //    fillDFGbyFunction(testFn, dfg);
    //    checkConservationOfMomentum<std::complex<double>>(
    //        dfg, DistributedHierarchization::hierarchizeBiorthogonal<std::complex<double>>);
    //  }
    //  {
    //    fillDFGbyFunction(testFn, dfg);
    //    checkConservationOfMomentum<std::complex<double>>(
    //        dfg, DistributedHierarchization::hierarchizeFullWeighting<std::complex<double>>);
    //  }
    //  std::cout << "random" << std::endl;
    //  DistributedFullGrid<real> dfg2(3, levels, comm, boundary, procs, forward);
    //  {
    //    // initialize dfg with random numbers
    //    fillDFGrandom(dfg2, -100., 100.);
    //    checkConservationOfMomentum<real>(
    //        dfg2, DistributedHierarchization::hierarchizeBiorthogonal<real>);
    //  }
    //  {
    //    fillDFGrandom(dfg2, -100., 100.);
    //    checkConservationOfMomentum<real>(
    //        dfg2, DistributedHierarchization::hierarchizeFullWeighting<real>);
    //  }
    //  std::cout << "const" << std::endl;
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

#ifdef NDEBUG
BOOST_AUTO_TEST_CASE(test_timing_parallel) {
  // large test case with timing
  for (BoundaryType b : std::vector<BoundaryType>({1, 2})) {
    std::vector<int> procs = {1, 8, 1, 1, 1, 1};
    auto dim = static_cast<DimType>(procs.size());
    std::vector<BoundaryType> boundary(dim, b);
    CommunicatorType comm = TestHelper::getComm(procs, std::vector<int>(dim, b == 1 ? 1 : 0));
    if (comm != MPI_COMM_NULL) {
      LevelVector levels = {1, 19, 1, 2, 1, 1};
      LevelVector lzero(dim, 0);
      LevelVector lhalf = levels;
      std::transform(lhalf.begin(), lhalf.end(), lhalf.begin(), [](int i) { return i / 2; });
      for (LevelVector lmin : std::vector<LevelVector>({lzero, lhalf, levels})) {
        auto forward = false;
        TestFn_1 testFn(levels);
        DistributedFullGrid<std::complex<double>> dfg(dim, levels, comm, boundary, procs, forward);
        MPI_Barrier(comm);
        auto start = std::chrono::high_resolution_clock::now();
        checkHierarchization(testFn, dfg, false, lmin);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        BOOST_TEST_MESSAGE("hat hierarchization time: " + std::to_string(duration.count()) +
                           " milliseconds with " + std::to_string(b) + "-sided boundary and lmin " +
                           std::to_string(combigrid::levelSum(lmin)));
      }
    }
  }
}
#endif  // def NDEBUG

BOOST_AUTO_TEST_SUITE_END()
