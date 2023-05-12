#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/test/unit_test.hpp>
#include <complex>
#include <cstdarg>
#include <iostream>
#include <random>
#include <vector>

#include "TaskConstParaboloid.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "fullgrid/FullGrid.hpp"
#include "hierarchization/DistributedHierarchization.hpp"
#include "io/H5InputOutput.hpp"
#include "mpi/MPIMemory.hpp"
#include "test_helper.hpp"
#include "utils/MonteCarlo.hpp"
#include "utils/Types.hpp"

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

void checkDistributedFullgridMemory(LevelVector& levels, bool forward = false) {
  auto commSize = getCommSize(MPI_COMM_WORLD);
  std::vector<size_t> groupSizes;
  // TODO allow non-pow2 group sizes
  size_t groupSize = 1;
  while (groupSize <= commSize) {
    groupSizes.push_back(groupSize);
    groupSize *= 2;
  }
  std::vector<long unsigned int> vmSizes(groupSizes.size());
  std::vector<long unsigned int> vmSizesReference(groupSizes.size());
  const auto dim = static_cast<DimType>(levels.size());
  std::vector<BoundaryType> boundary(dim, 2);

  for (size_t i = 0; i < groupSizes.size(); ++i) {
    // get memory footprints for after allocating different dfgs
    std::vector<int> procs(dim, 1);
    procs[0] = static_cast<int>(groupSizes[i]);
    CommunicatorType comm = TestHelper::getComm(procs);

    long unsigned int vmRSS, vmSize = 0;
    if (comm != MPI_COMM_NULL) {
      BOOST_TEST_CHECKPOINT("test distributedfullgrid memory " + std::to_string(groupSizes[i]));

      // // get current global memory footprint
      // mpimemory::get_all_memory_usage_kb(&vmRSS, &vmSizesReference[i], comm);
      {
        // create "empty" dfg (resp. the cartesian communicator)
        // and get current global memory footprint
        CommunicatorType cartesian_communicator;
        std::vector<int> intProcs(procs.rbegin(), procs.rend());
        std::vector<int> periods(dim, 0);
        int reorder = 0;
        MPI_Cart_create(comm, static_cast<int>(dim), &intProcs[0], &periods[0], reorder,
                        &cartesian_communicator);

        // this should scale linearly w.r.t. size of comm
        mpimemory::get_all_memory_usage_kb(&vmRSS, &vmSizesReference[i], comm);
        MPI_Comm_free(&cartesian_communicator);
      }
      {
        // create dfg and get same footprint
        DistributedFullGrid<double> dfg(dim, levels, comm, boundary, procs, forward);
        mpimemory::get_all_memory_usage_kb(&vmRSS, &vmSize, comm);

        // check that the local number of points is what we expect
        unsigned long numElementsRef = static_cast<unsigned long>(dfg.getNrElements());
        unsigned long numElements = 0;
        unsigned long numLocalElements = static_cast<unsigned long>(dfg.getNrLocalElements());
        MPI_Allreduce(&numLocalElements, &numElements, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);
        BOOST_TEST(numElements == numElementsRef);
      }
      vmSizes[i] = vmSize - vmSizesReference[i];
      // compare allocated memory sizes
      // check for linear scaling (grace 100% for memory size jitter & decomposition_ member)
      if (TestHelper::getRank(comm) == 0) {
        BOOST_TEST_CHECKPOINT("reference " + std::to_string(groupSizes[i]) + ": " +
                              std::to_string(vmSizesReference[i]) +
                              ", vmSize: " + std::to_string(vmSizes[i]));
        BOOST_TEST_MESSAGE("reference " + std::to_string(groupSizes[i]) + ": " +
                           std::to_string(vmSizesReference[i]) +
                           ", vmSize: " + std::to_string(vmSizes[i]));
        if (i > 0) {
          BOOST_TEST(static_cast<double>(vmSizes[i]) <= ((vmSizes[0] + 500) * 2.));
        }
      }
    }
    // update parallelization so it matches the next groupSize
    procs[i % dim] *= 2;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_CHECK(!TestHelper::testStrayMessages());
}

template <typename T>
real checkInnerBasisFunctionIntegral(const DistributedFullGrid<T>& dfg) {
  auto dim = dfg.getDimension();
  auto& levels = dfg.getLevels();
  auto& boundary = dfg.returnBoundaryFlags();

  IndexType nrElements = 1;
  for (DimType d = 0; d < dim; ++d) {
    nrElements *= (1 << levels[d]) + boundary[d] - 1;
  }
  BOOST_CHECK(nrElements == dfg.getNrElements());

  auto innerNodalBasisFunctionIntegral = dfg.getInnerNodalBasisFunctionIntegral();
  auto levelSum = combigrid::levelSum(levels);
  BOOST_CHECK_CLOSE(innerNodalBasisFunctionIntegral, oneOverPowOfTwo[levelSum],
                    TestHelper::tolerance);

  return innerNodalBasisFunctionIntegral;
}

template <typename T>
std::vector<double> checkCoordinates(const DistributedFullGrid<T>& dfg) {
  // check coordinates
  std::vector<double> coordsGlobal(dfg.getDimension());
  std::vector<double> coordsGlobalAgain(dfg.getDimension());
  IndexType globalIndex = -1;
  for (IndexType i = 0; i < dfg.getNrLocalElements(); ++i) {
    globalIndex = dfg.getGlobalLinearIndex(i);
    dfg.getCoordsGlobal(globalIndex, coordsGlobal);
    dfg.getCoordsLocal(i, coordsGlobalAgain);
    BOOST_CHECK_EQUAL_COLLECTIONS(coordsGlobal.begin(), coordsGlobal.end(),
                                  coordsGlobalAgain.begin(), coordsGlobalAgain.end());
    for (DimType d = 0; d < dfg.getDimension(); ++d) {
      BOOST_CHECK_GE(coordsGlobal[d], dfg.getLowerBoundsCoord(d));
    }
  }
  return coordsGlobalAgain;
}

void checkDistributedFullgrid(LevelVector& levels, std::vector<int>& procs,
                              std::vector<BoundaryType>& boundary, bool forward = false) {
  std::vector<int> periodic;
  for (const auto& b : boundary) {
    periodic.push_back(b == 1 ? 1 : 0);
  }
  CommunicatorType comm = TestHelper::getComm(procs, periodic);
  if (comm == MPI_COMM_NULL) return;

  if (TestHelper::getRank(comm) == 0) {
    std::cout << "test distributedfullgrid " << levels << procs << std::endl;
  }

  TestFn f;
  const auto dim = static_cast<DimType>(levels.size());

  // create dfg
  DistributedFullGrid<std::complex<double>> dfg(dim, levels, comm, boundary, procs, forward);

  // set function values
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);
    dfg.getData()[li] = f(coords);
  }
  BOOST_TEST_CHECKPOINT("set function values");

  // test addDistributedFullGrid, extractFromUniformSG
  LevelVector lmin = levels;
  LevelVector lmax = levels;
  for (DimType d = 0; d < dim; ++d) {
    lmax[d] *= 2;
  }
  DistributedSparseGridUniform<std::complex<double>> dsg(dim, lmax, lmin, comm);
  dsg.registerDistributedFullGrid(dfg);
  BOOST_TEST_CHECKPOINT("register uniform sg");
  DistributedFullGrid<std::complex<double>> dfg2(dim, levels, comm, boundary, procs, forward);
  dsg.registerDistributedFullGrid(dfg2);
  dsg.setZero();
  dsg.addDistributedFullGrid(dfg, 2.1);
  BOOST_TEST_CHECKPOINT("add to uniform sg");
  dfg2.extractFromUniformSG(dsg);
  BOOST_TEST_CHECKPOINT("extract from uniform sg");

  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);
    BOOST_TEST(2.1 * dfg.getData()[li] == dfg2.getData()[li]);
  }

  // test full grid interpolation
  std::vector<std::vector<double>> interpolationCoords;
  interpolationCoords.emplace_back(dim, 1. / std::sqrt(2.));
  interpolationCoords.emplace_back(dim, 1. / std::sqrt(3.));
  interpolationCoords.emplace_back(dim, 1. / std::sqrt(5.));
  interpolationCoords.emplace_back(dim, 0.5);
  interpolationCoords.emplace_back(dim, 0.5 + 1e-4);
  interpolationCoords.emplace_back(dim, 0.5 - 1e-4);
  interpolationCoords.emplace_back(dim, 0.5 + 1e-5);
  interpolationCoords.emplace_back(dim, 0.5 - 1e-5);
  auto interpolatedValues = dfg.getInterpolatedValues(interpolationCoords);

  for (size_t i = 0; i < interpolationCoords.size(); ++i) {
    BOOST_CHECK_CLOSE(interpolatedValues[i].real(), f(interpolationCoords[i]).real(),
                      TestHelper::tolerance);
  }

  // test norm calculation
  auto maxnorm = dfg.getLpNorm(0);
  auto onenorm = dfg.getLpNorm(1);
  auto twonorm = dfg.getLpNorm(2);
  if (std::all_of(boundary.begin(), boundary.end(), [](BoundaryType i) { return i == 2; })) {
    // check that InnerNodalBasisFunctionIntegral is correct
    double numFullInnerBasisFcns = 1.;
    for (DimType d = 0; d < dim; ++d) {
      numFullInnerBasisFcns *= static_cast<double>(dfg.length(d) - 1);
    }
    BOOST_CHECK_CLOSE(dfg.getInnerNodalBasisFunctionIntegral(), 1. / numFullInnerBasisFcns,
                      TestHelper::tolerance);

    std::vector<double> maxcoords(dim, 1.);
    BOOST_CHECK_EQUAL(f(maxcoords), maxnorm);
    // solution is a hyperplane, so the L1 norm is equal to the value in the middle
    std::vector<double> middlecoords(dim, 0.5);
    BOOST_CHECK_CLOSE(std::abs(f(middlecoords)), onenorm, TestHelper::tolerance);
  }
  // lazy for the two-norm, just check boundedness relations:
  BOOST_CHECK(twonorm <= onenorm);
  BOOST_CHECK(onenorm <= maxnorm);

  // test ghost layer exchange
  std::vector<int> subarrayExtents;
  for (DimType d = 0; d < dim; ++d) {
    auto ghostLayer = dfg.exchangeGhostLayerUpward(d, subarrayExtents);
    IndexVector offsets(dim);
    IndexType numElements = 1;
    for (DimType j = 0; j < dim; ++j) {
      offsets[j] = numElements;
      numElements *= subarrayExtents[j];
    }
    BOOST_CHECK_EQUAL(ghostLayer.size(), numElements);
    if (numElements > 0) {
      BOOST_CHECK_EQUAL(subarrayExtents[d], 1);
      for (size_t j = 0; j < numElements; ++j) {
        // calculate local axis indices of ghost layer points
        // == axis index of lowest layer in dim d for dfg
        IndexVector locAxisIndex(dim), globAxisIndex(dim);
        IndexType tmp = j;
        for (int i = static_cast<int>(dim) - 1; i >= 0; i--) {
          locAxisIndex[i] = tmp / offsets[i];
          tmp = tmp % offsets[i];
        }
        // calculate global indices and coordinates of ghost layer points
        dfg.getGlobalVectorIndex(locAxisIndex, globAxisIndex);
        --globAxisIndex[d];
        if (boundary[d] == 1 && globAxisIndex[d] < 0) {
          // wrap around
          globAxisIndex[d] += dfg.length(d);
        }
        std::vector<double> coords(dim);
        dfg.getCoordsGlobal(dfg.getGlobalLinearIndex(globAxisIndex), coords);
        BOOST_CHECK_EQUAL(ghostLayer[j], f(coords));
      }
    }
  }

  // test gatherFullgrid
  FullGrid<std::complex<double>> fg(dim, levels, boundary);
  dfg.gatherFullGrid(fg, 0);

  // only check on rank 0
  if (TestHelper::getRank(comm) == 0) {
    for (size_t i = 0; i < static_cast<size_t>(fg.getNrElements()); ++i) {
      std::vector<double> coords(dim);
      fg.getCoords(i, coords);
      BOOST_TEST(fg.getData()[i] == f(coords));
    }
  }

  // std::stringstream ss;
  // ss << "test_dfg_" << levels << procs << boundary << forward << ".vtk";
  // dfg.writePlotFileVTK(ss.str().c_str());

  // create distributed fg and copy values
  DistributedFullGrid<std::complex<double>> dfgCopy(
      dim, dfg.getLevels(), dfg.getCommunicator(), dfg.returnBoundaryFlags(),
      dfg.getParallelization(), true, dfg.getDecomposition());
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    dfgCopy.getData()[li] = dfg.getData()[li];
  }
  BOOST_TEST_CHECKPOINT("copied to dfgCopy");

  auto formerCorners = dfg.getCornersValues();
  auto formerBoundaryAvgValue =
      std::accumulate(formerCorners.begin(), formerCorners.end(),
                      static_cast<std::complex<double>>(0.), std::plus<std::complex<double>>()) /
      static_cast<double>(formerCorners.size());
  // test averageBoundaryValues
  dfgCopy.averageBoundaryValues();
  if (std::all_of(boundary.begin(), boundary.end(), [](BoundaryType b) { return b == 2; })) {
    // assert that all the corners values are the same now
    auto cornersValues = dfgCopy.getCornersValues();
    for (const auto& cornerValue : cornersValues) {
      BOOST_TEST(cornerValue == formerBoundaryAvgValue);
    }
  }

  // test lower to upper exchange
  for (DimType d = 0; d < dim; ++d) {
    // std::cout << dfg << std::endl;
    if (boundary[d] == 2) {
      dfg.writeLowerBoundaryToUpperBoundary(d);

      for (IndexType i = 0; i < dfg.getNrLocalElements(); ++i) {
        std::vector<double> coords(dim);
        dfg.getCoordsLocal(i, coords);
        if (std::abs((coords[d] - 1.)) < TestHelper::tolerance) {
          bool edge = false;
          // if other dimensions are at maximum too, we are at an edge
          // results may differ; skip
          for (DimType d_i = 0; d_i < dim; ++d_i) {
            if (d_i != d && std::abs((coords[d_i] - 1.)) < TestHelper::tolerance) {
              edge = true;
            }
          }
          if (!edge) {
            // check if the value is the same as on the lower boundary
            auto compareCoords = coords;
            compareCoords[d] = 0.;
            BOOST_CHECK_EQUAL(dfg.getData()[i], f(compareCoords));
          }
        } else {
          bool otherBoundary = false;
          for (DimType d_i = 0; d_i < dim; ++d_i) {
            if (d_i != d && std::abs((coords[d_i] - 1.)) < TestHelper::tolerance) {
              otherBoundary = true;
            }
          }
          if (!otherBoundary) {
            // make sure all other values remained the same
            BOOST_CHECK_EQUAL(dfg.getData()[i], f(coords));
          }
        }
      }
    }
  }

  // std::cout << dfg << std::endl;
}

BOOST_FIXTURE_TEST_SUITE(distributedfullgrid, TestHelper::BarrierAtEnd,
                         *boost::unit_test::timeout(600))

// with boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_minus3) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {1, 2};
  std::vector<int> procs = {1, 1};
  std::vector<BoundaryType> boundary(2, 2);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_minus2) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {2, 3};
  std::vector<int> procs = {1, 1};
  std::vector<BoundaryType> boundary(2, 2);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_minus1) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {1, 1, 2};
  std::vector<int> procs = {1, 1, 1};
  std::vector<BoundaryType> boundary(3, 2);
  checkDistributedFullgrid(levels, procs, boundary);
}

BOOST_AUTO_TEST_CASE(test_0) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {2, 3};
  std::vector<int> procs = {1, 1};
  std::vector<BoundaryType> boundary(2, 2);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_1) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(6));
  LevelVector levels = {2, 2};
  std::vector<int> procs = {2, 3};
  std::vector<BoundaryType> boundary(2, 2);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_2) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_3) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 2);
  checkDistributedFullgrid(levels, procs, boundary, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_4) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {2, 3};
  std::vector<int> procs = {2, 2};
  std::vector<BoundaryType> boundary(2, 2);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_5) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {1, 4, 2};
  std::vector<BoundaryType> boundary(3, 2);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_6) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 4, 6};
  std::vector<int> procs = {1, 4, 2};
  std::vector<BoundaryType> boundary(3, 2);
  checkDistributedFullgrid(levels, procs, boundary, true);
}

// without boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_7) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(2));
  LevelVector levels = {3, 3};
  std::vector<int> procs = {2, 1};
  std::vector<BoundaryType> boundary(2, 0);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_8) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_9) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {4, 4, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  checkDistributedFullgrid(levels, procs, boundary, true);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_10) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {2, 3};
  std::vector<int> procs = {2, 2};
  std::vector<BoundaryType> boundary(2, 0);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_11) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_12) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {2, 3, 4};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 0);
  checkDistributedFullgrid(levels, procs, boundary, true);
}

// partial boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_13) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {4, 4};
  std::vector<int> procs = {2, 2};
  std::vector<BoundaryType> boundary(2, 1);
  boundary[1] = 2;
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_14) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 3, 3};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 1);
  boundary[0] = 2;
  checkDistributedFullgrid(levels, procs, boundary);
}

// anisotropic

BOOST_AUTO_TEST_CASE(test_16) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  LevelVector levels = {2, 3};
  std::vector<int> procs = {2, 2};
  std::vector<BoundaryType> boundary(2, 1);
  boundary[0] = 2;
  checkDistributedFullgrid(levels, procs, boundary);
}
BOOST_AUTO_TEST_CASE(test_17) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(8));
  LevelVector levels = {3, 4, 3};
  std::vector<int> procs = {2, 2, 2};
  std::vector<BoundaryType> boundary(3, 1);
  boundary[2] = 2;
  checkDistributedFullgrid(levels, procs, boundary);
}

// memory
BOOST_AUTO_TEST_CASE(test_19) {
  // level sum = 17 -> if sizeof(double) == 8,
  // the grid's payload should be ~1 MB
  // (2^17 * 8 / 2^20 = 1)
  LevelVector levels = {9, 8};
  checkDistributedFullgridMemory(levels, false);
}
BOOST_AUTO_TEST_CASE(test_20) {
  LevelVector levels = {6, 6, 5};
  checkDistributedFullgridMemory(levels, false);
}
BOOST_AUTO_TEST_CASE(test_21) {
  LevelVector levels = {5, 4, 4, 4};
  checkDistributedFullgridMemory(levels, false);
}
BOOST_AUTO_TEST_CASE(test_22) {
  LevelVector levels = {3, 3, 3, 3, 3, 2};
  checkDistributedFullgridMemory(levels, false);
}

BOOST_AUTO_TEST_CASE(compare_coordinates_by_boundary) {
  std::vector<int> procs = {2, 2, 2, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    LevelVector fullGridLevel = {3, 2, 3, 1};
    DimType dim = static_cast<DimType>(procs.size());
    std::vector<BoundaryType> boundary(dim, 2);
    std::vector<BoundaryType> oneboundary(dim, 1);
    std::vector<BoundaryType> noboundary(dim, 0);
    DistributedFullGrid<real> dfgTwoBoundary(dim, fullGridLevel, comm, boundary, procs, false);
    DistributedFullGrid<real> dfgOneBoundary(dim, fullGridLevel, comm, oneboundary, procs, false);
    DistributedFullGrid<real> dfgNoBoundary(dim, fullGridLevel, comm, noboundary, procs, false);

    auto twoBoundaryIntegral = checkInnerBasisFunctionIntegral(dfgTwoBoundary);
    auto oneBoundaryIntegral = checkInnerBasisFunctionIntegral(dfgOneBoundary);
    auto noBoundaryIntegral = checkInnerBasisFunctionIntegral(dfgNoBoundary);
    BOOST_CHECK_CLOSE(twoBoundaryIntegral, oneBoundaryIntegral, TestHelper::tolerance);
    BOOST_CHECK_CLOSE(twoBoundaryIntegral, noBoundaryIntegral, TestHelper::tolerance);

    const auto& twoBoundaryGridSpacing = dfgTwoBoundary.getGridSpacing();
    const auto& oneBoundaryGridSpacing = dfgOneBoundary.getGridSpacing();
    const auto& noBoundaryGridSpacing = dfgNoBoundary.getGridSpacing();
    for (DimType d = 0; d < dim; d++) {
      BOOST_CHECK_CLOSE(twoBoundaryGridSpacing[d], oneBoundaryGridSpacing[d],
                        TestHelper::tolerance);
      BOOST_CHECK_CLOSE(twoBoundaryGridSpacing[d], noBoundaryGridSpacing[d], TestHelper::tolerance);
    }

    auto twoBoundaryHighestCoordinate = checkCoordinates(dfgTwoBoundary);
    auto oneBoundaryHighestCoordinate = checkCoordinates(dfgOneBoundary);
    auto noBoundaryHighestCoordinate = checkCoordinates(dfgNoBoundary);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        noBoundaryHighestCoordinate.begin(), noBoundaryHighestCoordinate.end(),
        oneBoundaryHighestCoordinate.begin(), oneBoundaryHighestCoordinate.end());
    for (DimType d = 0; d < dim; d++) {
      if (dfgTwoBoundary.getCartesianUtils().isOnUpperBoundaryInDimension(d)) {
        BOOST_CHECK_GT(twoBoundaryHighestCoordinate[d], noBoundaryHighestCoordinate[d]);
      } else {
        BOOST_CHECK_EQUAL(twoBoundaryHighestCoordinate[d], noBoundaryHighestCoordinate[d]);
      }
    }

    for (const auto* dfg : {&dfgTwoBoundary, &dfgOneBoundary, &dfgNoBoundary}) {
      auto localSizes = dfg->getLocalSizes();
      auto upperBounds = dfg->getUpperBounds();
      auto lowerBounds = dfg->getLowerBounds();
      auto boundsDiff = upperBounds - lowerBounds;
      BOOST_CHECK_EQUAL_COLLECTIONS(localSizes.begin(), localSizes.end(), boundsDiff.begin(),
                                    boundsDiff.end());
      auto globalSizes = dfg->getGlobalSizes();
      BOOST_CHECK(localSizes <= globalSizes);
      auto numDof = combigrid::getNumDofNodal(fullGridLevel, dfg->returnBoundaryFlags());
      BOOST_CHECK_EQUAL(numDof, std::accumulate(globalSizes.begin(), globalSizes.end(), 1,
                                                std::multiplies<IndexType>()));

      auto localSizesReduced = localSizes;
      MPI_Allreduce(MPI_IN_PLACE, localSizesReduced.data(), dim,
                    abstraction::getMPIDatatype(abstraction::getabstractionDataType<IndexType>()),
                    MPI_SUM, comm);
      auto numProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<int>());
      for (DimType d = 0; d < dim; d++) {
        // correct the reduced number of points that were added from processes parallel in other
        // dimensions
        localSizesReduced[d] = localSizesReduced[d] / (numProcs / procs[d]);
      }
      BOOST_CHECK_EQUAL_COLLECTIONS(localSizesReduced.begin(), localSizesReduced.end(),
                                    globalSizes.begin(), globalSizes.end());
    }
  }
}

BOOST_AUTO_TEST_CASE(interpolation_test) {
  std::vector<int> procs = {2, 2, 2, 1, 1, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    LevelVector fullGridLevel = {3, 2, 3, 1, 1, 1};
    DimType dim = static_cast<DimType>(procs.size());
    std::vector<BoundaryType> boundary(dim, 2);
    std::vector<BoundaryType> oneboundary(dim, 1);
    std::vector<BoundaryType> noboundary(dim, 0);
    DistributedFullGrid<real> dfgTwoBoundary(dim, fullGridLevel, comm, boundary, procs, false);
    DistributedFullGrid<real> dfgOneBoundary(dim, fullGridLevel, comm, oneboundary, procs, false);
    DistributedFullGrid<real> dfgNoBoundary(dim, fullGridLevel, comm, noboundary, procs, false);

    // set function values on dfgs
    // choose function that will be 0 on boundary
    ParaboloidFn<CombiDataType> f;
    std::vector<double> coords(dim);
    for (IndexType li = 0; li < dfgTwoBoundary.getNrLocalElements(); ++li) {
      dfgTwoBoundary.getCoordsLocal(li, coords);
      dfgTwoBoundary.getData()[li] = f(coords);
    }
    for (IndexType li = 0; li < dfgOneBoundary.getNrLocalElements(); ++li) {
      dfgOneBoundary.getCoordsLocal(li, coords);
      dfgOneBoundary.getData()[li] = f(coords);
    }
    BOOST_CHECK_EQUAL(dfgTwoBoundary.getData()[0], dfgOneBoundary.getData()[0]);
    BOOST_CHECK_EQUAL(dfgTwoBoundary.getData()[1], dfgOneBoundary.getData()[1]);
    for (IndexType li = 0; li < dfgNoBoundary.getNrLocalElements(); ++li) {
      dfgNoBoundary.getCoordsLocal(li, coords);
      dfgNoBoundary.getData()[li] = f(coords);
    }

    auto numMCCoordinates = 1e2;
    std::vector<std::vector<double>> interpolationCoords =
        montecarlo::getRandomCoordinates(numMCCoordinates, static_cast<size_t>(dim));

    auto interpolatedValuesTwoBoundary = dfgTwoBoundary.getInterpolatedValues(interpolationCoords);
    auto interpolatedValuesOneBoundary = dfgOneBoundary.getInterpolatedValues(interpolationCoords);
    auto interpolatedValuesNoBoundary = dfgNoBoundary.getInterpolatedValues(interpolationCoords);

    if (TestHelper::getRank(comm) == 0) {
      BOOST_CHECK_EQUAL_COLLECTIONS(
          interpolatedValuesTwoBoundary.begin(), interpolatedValuesTwoBoundary.end(),
          interpolatedValuesOneBoundary.begin(), interpolatedValuesOneBoundary.end());
      BOOST_CHECK_EQUAL_COLLECTIONS(
          interpolatedValuesTwoBoundary.begin(), interpolatedValuesTwoBoundary.end(),
          interpolatedValuesNoBoundary.begin(), interpolatedValuesNoBoundary.end());
    }
  }
}

#ifdef NDEBUG  // speed test -> run only in release mode
BOOST_AUTO_TEST_CASE(interpolation_speed_test) {
  std::vector<int> procs = {2, 2, 2, 1, 1, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    LevelVector fullGridLevel = {3, 2, 3, 1, 1, 1};
    DimType dim = static_cast<DimType>(procs.size());
    std::vector<BoundaryType> boundary(dim, 2);
    std::vector<BoundaryType> oneboundary(dim, 1);
    std::vector<BoundaryType> noboundary(dim, 0);
    DistributedFullGrid<real> dfgTwoBoundary(dim, fullGridLevel, comm, boundary, procs, false);
    DistributedFullGrid<real> dfgOneBoundary(dim, fullGridLevel, comm, oneboundary, procs, false);
    DistributedFullGrid<real> dfgNoBoundary(dim, fullGridLevel, comm, noboundary, procs, false);

    auto numMCCoordinates = 1e6;
    std::vector<std::vector<double>> interpolationCoords =
        montecarlo::getRandomCoordinates(numMCCoordinates, static_cast<size_t>(dim));

    MPI_Barrier(comm);
    auto start = std::chrono::high_resolution_clock::now();
    auto interpolatedValues = dfgTwoBoundary.getInterpolatedValues(interpolationCoords);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to interpolate 1e6 values on two-boundary grid: " << duration.count()
                                                                               << " milliseconds");
    BOOST_CHECK(duration.count() < 25000);

    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    interpolatedValues = dfgOneBoundary.getInterpolatedValues(interpolationCoords);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to interpolate 1e6 values on one-boundary grid: " << duration.count()
                                                                               << " milliseconds");
    BOOST_CHECK(duration.count() < 25000);

    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    interpolatedValues = dfgNoBoundary.getInterpolatedValues(interpolationCoords);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to interpolate 1e6 values on no-boundary grid: " << duration.count()
                                                                              << " milliseconds");
    BOOST_CHECK(duration.count() < 25000);
  }
}
#endif  // NDEBUG

BOOST_AUTO_TEST_CASE(test_get1dIndicesLocal) {
  std::vector<int> procs = {3};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = static_cast<DimType>(procs.size());
    std::vector<BoundaryType> boundary(dim, 2);
    std::vector<BoundaryType> oneboundary(dim, 1);
    std::vector<BoundaryType> noboundary(dim, 0);
    LevelVector fullGridLevel = {3};
    std::vector<IndexVector> decomposition = {{0, 4, 6}};
    DistributedFullGrid<real> dfg(dim, fullGridLevel, comm, boundary, procs, true, decomposition);
    DistributedFullGrid<real> dfgOneBoundary(dim, fullGridLevel, comm, oneboundary, procs, true,
                                             decomposition);
    DistributedFullGrid<real> dfgNoBoundary(dim, fullGridLevel, comm, noboundary, procs, true,
                                            decomposition);

    IndexVector indices;
    IndexVector expected;
    dfg.get1dIndicesLocal(0, 3, indices);
    auto rank = dfg.getRank();
    if (rank == 0) {
      expected = {1, 3};
    } else if (rank == 1 || rank == 2) {
      expected = {1};
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(), expected.end());
    indices.clear();
    dfg.get1dIndicesLocal(0, 2, indices);
    if (rank == 0) {
      expected = {2};
    } else if (rank == 1) {
      expected = {};
    } else if (rank == 2) {
      expected = {0};
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(), expected.end());
    indices.clear();
    dfg.get1dIndicesLocal(0, 1, indices);
    if (rank == 0) {
      expected = {0};
    } else if (rank == 1) {
      expected = {0};
    } else if (rank == 2) {
      expected = {2};
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(), expected.end());

    indices.clear();
    dfgOneBoundary.get1dIndicesLocal(0, 3, indices);
    if (rank == 0) {
      expected = {1, 3};
    } else if (rank == 1 || rank == 2) {
      expected = {1};
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(), expected.end());
    indices.clear();
    dfgOneBoundary.get1dIndicesLocal(0, 2, indices);
    if (rank == 0) {
      expected = {2};
    } else if (rank == 1) {
      expected = {};
    } else if (rank == 2) {
      expected = {0};
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(), expected.end());
    indices.clear();
    dfgOneBoundary.get1dIndicesLocal(0, 1, indices);
    if (rank == 0) {
      expected = {0};
    } else if (rank == 1) {
      expected = {0};
    } else if (rank == 2) {
      expected = {};
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(), expected.end());

    indices.clear();
    dfgNoBoundary.get1dIndicesLocal(0, 3, indices);
    if (rank == 0) {
      expected = {0, 2};
    } else if (rank == 1 || rank == 2) {
      expected = {0};
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(), expected.end());
    indices.clear();
    dfgNoBoundary.get1dIndicesLocal(0, 2, indices);
    if (rank == 0) {
      expected = {1};
    } else if (rank == 1) {
      expected = {1};
    } else if (rank == 2) {
      expected = {};
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(), expected.end());
    indices.clear();
    dfgNoBoundary.get1dIndicesLocal(0, 1, indices);
    if (rank == 0) {
      expected = {3};
    } else if (rank == 1) {
      expected = {};
    } else if (rank == 2) {
      expected = {};
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(), expected.end());
  }
}

BOOST_AUTO_TEST_CASE(test_get1dIndicesLocal_boundary_firstdim) {
  std::vector<int> procs = {5, 1, 1, 1, 1, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = static_cast<DimType>(procs.size());
    LevelVector fullGridLevel = {19, 2, 2, 2, 2, 2};
    std::vector<IndexVector> decomposition = {
        {0, 98304, 196609, 327680, 425985}, {0}, {0}, {0}, {0}, {0}};
    std::vector<BoundaryType> boundary(dim, 2);

    DistributedFullGrid<real> dfg(dim, fullGridLevel, comm, boundary, procs, true, decomposition);

    LevelVector level = {6, 5, 4, 3, 2, 1};  //{18, 1, 2, 3, 4, 5},
    for (DimType d = 0; d < dim; ++d) {
      IndexVector indices;
      dfg.get1dIndicesLocal(d, level[d], indices);
      auto rank = dfg.getRank();
      if (rank == 0 && d == 0) {
        IndexVector expected = {8192, 24576, 40960, 57344, 73728, 90112};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (rank == 1 && d == 0) {
        IndexVector expected = {8192, 24576, 40960, 57344, 73728, 90112};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (rank == 2 && d == 0) {
        IndexVector expected = {8191, 24575, 40959, 57343, 73727, 90111, 106495, 122879};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (rank == 3 && d == 0) {
        IndexVector expected = {8192, 24576, 40960, 57344, 73728, 90112};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (rank == 4 && d == 0) {
        IndexVector expected = {8191, 24575, 40959, 57343, 73727, 90111};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (d == 4) {
        IndexVector expected = {1, 3};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (d == 5) {
        IndexVector expected = {0, 2, 4};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else {
        BOOST_CHECK(indices.empty());
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_get1dIndicesLocal_boundary_threedim) {
  std::vector<int> procs = {1, 2, 1, 2, 1, 2};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = static_cast<DimType>(procs.size());
    LevelVector fullGridLevel = {2, 3, 4, 5, 6, 7};
    std::vector<IndexVector> decomposition = {{0}, {0, 3}, {0}, {0, 17}, {0}, {0, 63}};

    std::vector<BoundaryType> boundary(dim, 2);

    DistributedFullGrid<real> dfg(dim, fullGridLevel, comm, boundary, procs, true, decomposition);

    LevelVector level = {5, 3, 4, 5, 2, 5};  //{18, 1, 2, 3, 4, 5},
    for (DimType d = 0; d < dim; ++d) {
      IndexVector indices;
      dfg.get1dIndicesLocal(d, level[d], indices);
      auto rank = dfg.getRank();
      if (d == 0) {
        BOOST_CHECK(indices.empty());
      } else if (d == 1) {
        if (rank < 4) {
          IndexVector expected = {1};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        } else {
          IndexVector expected = {0, 2, 4};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        }
      } else if (d == 2) {
        IndexVector expected = {1, 3, 5, 7, 9, 11, 13, 15};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (d == 3) {
        if ((rank / 2) % 2 == 0) {
          IndexVector expected = {1, 3, 5, 7, 9, 11, 13, 15};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        } else {
          IndexVector expected = {0, 2, 4, 6, 8, 10, 12, 14};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        }
      } else if (d == 4) {
        IndexVector expected = {16, 48};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (d == 5) {
        if (rank % 2 == 0) {
          IndexVector expected = {4, 12, 20, 28, 36, 44, 52, 60};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        } else {
          IndexVector expected = {5, 13, 21, 29, 37, 45, 53, 61};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_get1dIndicesLocal_noboundary_firstdim) {
  std::vector<int> procs = {5, 1, 1, 1, 1, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = static_cast<DimType>(procs.size());
    LevelVector fullGridLevel = {19, 2, 2, 2, 2, 2};
    std::vector<IndexVector> decomposition = {
        {0, 98304, 196609, 327680, 425985}, {0}, {0}, {0}, {0}, {0}};

    std::vector<BoundaryType> boundary(dim, 0);

    DistributedFullGrid<real> dfg(dim, fullGridLevel, comm, boundary, procs, true, decomposition);

    LevelVector level = {6, 5, 4, 3, 2, 1};
    for (DimType d = 0; d < dim; ++d) {
      IndexVector indices;
      dfg.get1dIndicesLocal(d, level[d], indices);
      auto rank = dfg.getRank();
      if (rank == 0 && d == 0) {
        IndexVector expected = {8191, 24575, 40959, 57343, 73727, 90111};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (rank == 1 && d == 0) {
        IndexVector expected = {8191, 24575, 40959, 57343, 73727, 90111};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (rank == 2 && d == 0) {
        IndexVector expected = {8190, 24574, 40958, 57342, 73726, 90110, 106494, 122878};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (rank == 3 && d == 0) {
        IndexVector expected = {8191, 24575, 40959, 57343, 73727, 90111};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (rank == 4 && d == 0) {
        IndexVector expected = {8190, 24574, 40958, 57342, 73726, 90110};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (d == 4) {
        IndexVector expected = {0, 2};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (d == 5) {
        IndexVector expected = {1};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else {
        BOOST_CHECK(indices.empty());
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_get1dIndicesLocal_noboundary_threedim) {
  std::vector<int> procs = {1, 2, 1, 2, 1, 2};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = static_cast<DimType>(procs.size());
    LevelVector fullGridLevel = {2, 3, 4, 5, 6, 7};
    std::vector<IndexVector> decomposition = {{0}, {0, 3}, {0}, {0, 17}, {0}, {0, 63}};

    std::vector<BoundaryType> boundary(dim, 0);

    DistributedFullGrid<real> dfg(dim, fullGridLevel, comm, boundary, procs, true, decomposition);

    LevelVector level = {5, 3, 4, 5, 2, 5};  //{18, 1, 2, 3, 4, 5},
    for (DimType d = 0; d < dim; ++d) {
      IndexVector indices;
      dfg.get1dIndicesLocal(d, level[d], indices);
      auto rank = dfg.getRank();
      if (d == 0) {
        BOOST_CHECK(indices.empty());
      } else if (d == 1) {
        if (rank < 4) {
          IndexVector expected = {0, 2};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        } else {
          IndexVector expected = {1, 3};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        }
      } else if (d == 2) {
        IndexVector expected = {0, 2, 4, 6, 8, 10, 12, 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (d == 3) {
        if ((rank / 2) % 2 == 0) {
          IndexVector expected = {0, 2, 4, 6, 8, 10, 12, 14, 16};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        } else {
          IndexVector expected = {1, 3, 5, 7, 9, 11, 13};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        }
      } else if (d == 4) {
        IndexVector expected = {15, 47};
        BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                      expected.end());
      } else if (d == 5) {
        if (rank % 2 == 0) {
          IndexVector expected = {3, 11, 19, 27, 35, 43, 51, 59};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        } else {
          IndexVector expected = {4, 12, 20, 28, 36, 44, 52, 60};
          BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), expected.begin(),
                                        expected.end());
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_registerUniformSG) {
  std::vector<int> procs = {5, 1, 1, 1, 1, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = static_cast<DimType>(procs.size());
    LevelVector lmin(dim, 1);
    LevelVector lmax(dim, 17);
    LevelVector fullGridLevel = {18, 1, 1, 1, 1, 2};
    std::vector<BoundaryType> boundary(dim, 2);
    std::vector<IndexVector> decomposition = {
        {0, 49152, 98304, 163840, 212993}, {0}, {0}, {0}, {0}, {0}};

    MPI_Barrier(comm);
    auto start = std::chrono::high_resolution_clock::now();
    DistributedFullGrid<real> dfg(dim, fullGridLevel, comm, boundary, procs, true, decomposition);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to create full grid w/ level sum 23: " << duration.count()
                                                                    << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 2500);
#endif
    // go to the other extreme anisotropy by reversing the level vector
    std::reverse(fullGridLevel.begin(), fullGridLevel.end());
    std::reverse(boundary.begin(), boundary.end());
    decomposition[0] = {0, 1, 2, 3, 4};
    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    DistributedFullGrid<real> otherDfg(dim, fullGridLevel, comm, boundary, procs, true,
                                       decomposition);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to create other full grid w/ level sum 23: " << duration.count()
                                                                          << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 2500);
#endif

    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    DistributedSparseGridUniform<real> dsg(dim, lmax, lmin, comm);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to create sparse grid: " << duration.count() << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 5000);
#endif

    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    dsg.registerDistributedFullGrid(dfg);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to register sparse grid: " << duration.count() << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 5000);
#endif

    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    dsg.registerDistributedFullGrid(otherDfg);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to register sparse grid again: " << duration.count()
                                                              << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 50000);
#endif

    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    // this is not the "correct" communicator, but using it here so something is communicated
    dsg.reduceSubspaceSizes(comm);
    dsg.setZero();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to create sparse grid data: " << duration.count() << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 15000);
#endif

    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    dsg.addDistributedFullGrid(dfg, 1.);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to add to sparse grid: " << duration.count() << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 15000);
#endif

    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    dsg.addDistributedFullGrid(otherDfg, 1.);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to add to sparse grid again: " << duration.count() << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 15000);
#endif

    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    dfg.extractFromUniformSG(dsg);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to extract from sparse grid: " << duration.count() << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 15000);
#endif

    MPI_Barrier(comm);
    start = std::chrono::high_resolution_clock::now();
    otherDfg.extractFromUniformSG(dsg);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to extract from sparse grid again: " << duration.count()
                                                                  << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 15000);
#endif
  }
}

BOOST_AUTO_TEST_CASE(test_evalDFG) {
  std::vector<int> procs = {5, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = static_cast<DimType>(procs.size());
    size_t numCoordinates = 1000;
    auto interpolationCoords = montecarlo::getRandomCoordinates(numCoordinates, dim);
    // make sure there are some corner cases
    interpolationCoords.push_back(std::vector<double>(dim, 1e-10));
    interpolationCoords.push_back(std::vector<double>(dim, 1. - 1e-10));

    LevelVector fullGridLevel = {5, 5};
    for (auto b : std::vector<BoundaryType>({1, 2})) {
      BOOST_TEST_CHECKPOINT("Testing boundary type " + std::to_string(b));
      std::vector<BoundaryType> boundary(dim, b);
      // create and initialize DFG
      DistributedFullGrid<real> dfg(dim, fullGridLevel, comm, boundary, procs, false);
      std::vector<double> coords(dim);
      for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
        dfg.getCoordsLocal(li, coords);
        dfg.getData()[li] = 1.;
      }
      // evaluate all at once
      auto interpolatedValues = dfg.getInterpolatedValues(interpolationCoords);
      for (size_t i = 0; i < numCoordinates; ++i) {
        BOOST_CHECK_CLOSE(interpolatedValues[i], 1., TestHelper::tolerance);
      }

      // evaluate first and then reduce all values
      decltype(interpolatedValues) stepWiseInterpolatedValues(numCoordinates, 0.);
      for (size_t i = 0; i < numCoordinates; ++i) {
        stepWiseInterpolatedValues[i] = dfg.evalLocal(interpolationCoords[i]);
      }
      // reduce interpolated values within DFG's processes
      MPI_Allreduce(
          MPI_IN_PLACE, stepWiseInterpolatedValues.data(), static_cast<int>(numCoordinates),
          abstraction::getMPIDatatype(abstraction::getabstractionDataType<real>()), MPI_SUM, comm);
      for (size_t i = 0; i < numCoordinates; ++i) {
        BOOST_CHECK_CLOSE(stepWiseInterpolatedValues[i], interpolatedValues[i],
                          TestHelper::tolerance);
        BOOST_CHECK_CLOSE(stepWiseInterpolatedValues[i], 1., TestHelper::tolerance);
      }
      MPI_Barrier(comm);
    }
  }
}

BOOST_AUTO_TEST_CASE(test_massLoss2D) {
  std::vector<int> procs = {1, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = 2;
    std::vector<BoundaryType> boundary(dim, 1);
    // "true" solutions
    std::vector<std::string> fileNamePrefixes{"sneaky_peaky_mass_loss_", "too_peaky_"};
    auto sneakyPeakyMassLoss = [](const std::vector<double>& coords) {
      assert(coords.size() == 2);
      real result = 1.;
      for (const auto& c : coords) {
        assert(c >= 0. && c <= 1.);
        result *= 4. * std::max(0., 1. - std::abs(1. - 4. * c));
      }
      return result;
    };
    auto tooPeaky = [](const std::vector<double>& coords) {
      assert(coords.size() == 2);
      // only at some points there will be something
      // std::vector<real> peakCoordinates {0.125, 0.375, 0.625, 0.875};
      // std::vector<real> peakCoordinates {0., 0.25, 0.5, 0.75, 1.};
      std::vector<real> peakCoordinates{0., 0.5, 1.};
      for (const auto& peakCoordinate : peakCoordinates) {
        if (std::abs(coords[0] - peakCoordinate) < 0.01) {
          return 1.;
        }
      }
      return 0.;
    };
    std::vector<std::function<real(const std::vector<double>&)>> trueSolutions{sneakyPeakyMassLoss,
                                                                               tooPeaky};
    const std::vector<LevelVector> fullGridLevels = {{3, 1}, {2, 1}, {2, 2}, {1, 2}, {1, 3}};
    const std::vector<real> coefficients{1., -1., 1., -1., 1.};
    auto hat = HierarchicalHatPeriodicBasisFunction();
    auto biorthogonal = BiorthogonalPeriodicBasisFunction();
    auto fullWeighting = FullWeightingPeriodicBasisFunction();
    std::vector<BasisFunctionBasis*> bases{&hat, &biorthogonal, &fullWeighting};

    // for different scenarios
    for (size_t s = 0; s < fileNamePrefixes.size(); ++s) {
      // remove previously written, if any
      auto deleteStatus = system(("rm " + fileNamePrefixes[s] + "*.h5").c_str());
      BOOST_WARN_GE(deleteStatus, 0);
      // for different basis functions
      for (size_t b = 0; b < bases.size(); ++b) {
        auto basisTypeVector = std::vector<BasisFunctionBasis*>(dim, bases[b]);

        // create and initialize DFG
        std::vector<DistributedFullGrid<real>> dfgs;
        for (size_t i = 0; i < fullGridLevels.size(); ++i) {
          dfgs.emplace_back(dim, fullGridLevels[i], comm, boundary, procs, false);
          std::vector<double> coords(dim);
          for (IndexType li = 0; li < dfgs[i].getNrLocalElements(); ++li) {
            dfgs[i].getCoordsLocal(li, coords);
            dfgs[i].getData()[li] = trueSolutions[s](coords);
          }
        }

        // combine onto sparse grid and re-distribute again
        {
          DistributedSparseGridUniform<real> dsg(dim, {3, 3}, {1, 1}, comm);
          for (const auto& dfg : dfgs) {
            dsg.registerDistributedFullGrid(dfg);
          }
          dsg.setZero();
          LevelVector zeroLMin(dim, 0);
          for (size_t i = 0; i < dfgs.size(); ++i) {
            DistributedHierarchization::hierarchize(dfgs[i], {true, true}, basisTypeVector,
                                                    zeroLMin);
            dsg.addDistributedFullGrid(dfgs[i], coefficients[i]);
          }
          // extract from sparse grid again
          for (size_t i = 0; i < dfgs.size(); ++i) {
            dfgs[i].extractFromUniformSG(dsg);
            DistributedHierarchization::dehierarchize(dfgs[i], {true, true}, basisTypeVector,
                                                      zeroLMin);
          }
        }

        std::vector<real> oneDCoordinates = {0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875};
        std::vector<std::vector<real>> interpolationCoords;
        for (const auto& x : oneDCoordinates) {
          for (const auto& y : oneDCoordinates) {
            interpolationCoords.push_back({x, y});
          }
        }

        // call interpolation function on tasks and write out task-wise
        for (size_t i = 0; i < dfgs.size(); ++i) {
          const auto dfgValues = dfgs[i].getInterpolatedValues(interpolationCoords);
          // cycle through ranks to write
          std::string saveFilePath =
              fileNamePrefixes[s] + std::to_string(b) + "_task_" + std::to_string(i) + ".h5";
          std::string groupName = "run_";
          std::string datasetName = "interpolated_1";
#ifdef DISCOTEC_USE_HIGHFIVE
          h5io::writeValuesToH5File(dfgValues, saveFilePath, groupName, datasetName, 1.);
#endif
        }
        // TODO test if the values are the same, and if the mass is (not) conserved
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
