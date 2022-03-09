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

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGridEnsemble.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIMemory.hpp"
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
  const DimType dim = levels.size();
  IndexVector procs(dim, 1);
  std::vector<bool> boundary(dim, true);

  for (size_t i = 0; i < groupSizes.size(); ++i) {
    // get memory footprints for after allocating different dfgs
    CommunicatorType comm = TestHelper::getComm(groupSizes[i]);

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
        MPI_Cart_create(comm, static_cast<int>(dim), &intProcs[0], &periods[0], reorder, &cartesian_communicator);

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
        MPI_Allreduce(&numLocalElements, &numElements, 1,
                      MPI_UNSIGNED_LONG,
                      MPI_SUM, comm);
        BOOST_TEST(numElements == numElementsRef);
      }
      vmSizes[i] = vmSize - vmSizesReference[i];
      // compare allocated memory sizes
      // check for linear scaling (grace 100% for memory size jitter & decomposition_ member)
      if (TestHelper::getRank(comm) == 0) {
        BOOST_TEST_CHECKPOINT( "reference " + std::to_string(groupSizes[i]) +
                                ": " + std::to_string(vmSizesReference[i]) +
                                ", vmSize: " + std::to_string(vmSizes[i]) );
        if (i > 0) {
          BOOST_TEST(static_cast<double>(vmSizes[i]) <= ((vmSizes[0] + 500) * 2.));
        }
      }
    }
    // update parallelization so it matches the next groupSize
    procs[i % dim] *= 2;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  TestHelper::testStrayMessages();
}


template <typename dfgTestCombiType>
void compareResults(const DimType dim, DistributedFullGrid<dfgTestCombiType>& dfg,
                    DistributedFullGrid<dfgTestCombiType>& dfg2) {
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);
    BOOST_TEST(2.1 * dfg.getData()[li] == dfg2.getData()[li]);
  }
}

template <typename dfgTestCombiType>
void testGatherFullGrid(LevelVector& levels, std::vector<bool>& boundary, TestFn f,
                        DistributedFullGrid<dfgTestCombiType>& dfg, CommunicatorType comm) {
  const DimType dim = levels.size();

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
}

void checkDistributedFullgridLinDG(LevelVector& levels, IndexVector& procs,
                                   std::vector<bool>& boundary, int size, bool forward = false) {
  using dfgTestCombiType = CombiDataType;
  CommunicatorType comm = TestHelper::getComm(size);
  if (comm != MPI_COMM_NULL) {
    BOOST_TEST_CHECKPOINT("test distributedfullgrid Ensemble ");
    // + std::to_string(levels) + std::to_string(procs));

    TestFn f;
    const DimType dim = levels.size();

    BOOST_TEST_CHECKPOINT("creating ensemble");
    // create dfg ensemble
    DFGEnsemble dfgEnsemble(dim, levels, comm, boundary, procs, forward);

    IndexType nrElements = 1;
    for (DimType d = 0; d < dim; ++d) {
      nrElements *= (1 << levels[d]) + (boundary[d] ? 1 : -1);
    }

    for (auto& dfg : dfgEnsemble) {
      BOOST_CHECK(nrElements == dfg->getNrElements());

      // set function values
      for (IndexType li = 0; li < dfg->getNrLocalElements(); ++li) {
        std::vector<double> coords(dim);
        dfg->getCoordsLocal(li, coords);
        dfg->getData()[li] = std::abs(f(coords));
      }
      BOOST_TEST_CHECKPOINT("set function values");
    }

    // test addToUniformSG, extractFromUniformSG
    LevelVector lmin = levels;
    LevelVector lmax = levels;
    for (DimType d = 0; d < dim; ++d) {
      lmax[d] *= 2;
    }

    for (auto& dfg : dfgEnsemble) {
      DistributedSparseGridUniform<dfgTestCombiType> dsg(dim, lmax, lmin, boundary, comm);
      dfg->registerUniformSG(dsg);
      BOOST_TEST_CHECKPOINT("register uniform sg");
      dfg->addToUniformSG(dsg, 2.1);
      BOOST_TEST_CHECKPOINT("add to uniform sg");
      DistributedFullGrid<dfgTestCombiType> dfg2(dim, levels, comm, boundary, procs, forward);
      dfg2.registerUniformSG(dsg);
      dfg2.extractFromUniformSG(dsg);
      BOOST_TEST_CHECKPOINT("extract from uniform sg");

      compareResults(dim, *dfg, dfg2);
    }

    // test ensemble interpolation
    BOOST_TEST_CHECKPOINT("test ensemble interpolation");
    std::vector<std::vector<double>> interpolationCoords;
    interpolationCoords.emplace_back(dim, 1. / std::sqrt(2.));
    interpolationCoords.emplace_back(dim, 1. / std::sqrt(3.));
    interpolationCoords.emplace_back(dim, 1. / std::sqrt(5.));
    interpolationCoords.emplace_back(dim, 0.5);  // this one is at discontinuity -> bad!
    interpolationCoords.emplace_back(dim, 0.5 + 1e-4);
    interpolationCoords.emplace_back(dim, 0.5 - 1e-4);
    interpolationCoords.emplace_back(dim, 0.5 + 1e-5);
    interpolationCoords.emplace_back(dim, 0.5 - 1e-5);
    BOOST_CHECK_CLOSE(std::abs(dfgEnsemble.eval(interpolationCoords[0])),
                      std::abs(f(interpolationCoords[0])), TestHelper::tolerance);

    auto interpolatedValues = dfgEnsemble.getInterpolatedValues(interpolationCoords);

    for (size_t i = 0; i < interpolationCoords.size(); ++i) {
      BOOST_CHECK_CLOSE(std::abs(interpolatedValues[i]), std::abs(f(interpolationCoords[i])),
                        TestHelper::tolerance);
    }
  }
}

void checkDistributedFullgrid(LevelVector& levels, IndexVector& procs, std::vector<bool>& boundary,
                              int size, bool forward = false) {
  if (std::all_of(boundary.begin(), boundary.end(), [](bool b) { return b; })) {
    BOOST_TEST_CHECKPOINT("calling checkDistributedFullgridLinDG");
    checkDistributedFullgridLinDG(levels, procs, boundary, size, forward);
  }

  using dfgTestCombiType = std::complex<double>;
  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) return;

  TestFn f;
  const DimType dim = levels.size();

  // create dfg //TODO make this a dfgEnsemble
  DistributedFullGrid<dfgTestCombiType> dfg(dim, levels, comm, boundary, procs, forward);

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
  BOOST_TEST_CHECKPOINT("set function values");

  // test addToUniformSG, extractFromUniformSG
  LevelVector lmin = levels;
  LevelVector lmax = levels;
  for (DimType d = 0; d < dim; ++d) {
    lmax[d] *= 2;
  }
  DistributedSparseGridUniform<dfgTestCombiType> dsg(dim, lmax, lmin, boundary, comm);
  dfg.registerUniformSG(dsg);
  BOOST_TEST_CHECKPOINT("register uniform sg");
  dfg.addToUniformSG(dsg, 2.1);
  BOOST_TEST_CHECKPOINT("add to uniform sg");
  DistributedFullGrid<dfgTestCombiType> dfg2(dim, levels, comm, boundary, procs, forward);
  dfg2.registerUniformSG(dsg);
  dfg2.extractFromUniformSG(dsg);
  BOOST_TEST_CHECKPOINT("extract from uniform sg");

  compareResults(dim, dfg, dfg2);
  testGatherFullGrid(levels, boundary, f, dfg, comm);

  // test full grid interpolation
  std::vector<std::vector<double>> interpolationCoords;
  interpolationCoords.emplace_back(dim, 1./std::sqrt(2.));
  interpolationCoords.emplace_back(dim, 1./std::sqrt(3.));
  interpolationCoords.emplace_back(dim, 1./std::sqrt(5.));
  interpolationCoords.emplace_back(dim, 0.5);
  interpolationCoords.emplace_back(dim, 0.5 + 1e-4);
  interpolationCoords.emplace_back(dim, 0.5 - 1e-4);
  interpolationCoords.emplace_back(dim, 0.5 + 1e-5);
  interpolationCoords.emplace_back(dim, 0.5 - 1e-5);
  BOOST_CHECK_CLOSE(dfg.eval(interpolationCoords[0]).real(), f(interpolationCoords[0]).real(), TestHelper::tolerance);

  auto interpolatedValues = dfg.getInterpolatedValues(interpolationCoords);

  for (size_t i = 0; i < interpolationCoords.size(); ++i) {
    BOOST_CHECK_CLOSE(interpolatedValues[i].real(), f(interpolationCoords[i]).real(), TestHelper::tolerance);
  }

  // test norm calculation
  auto maxnorm = dfg.getLpNorm(0);
  auto onenorm = dfg.getLpNorm(1);
  auto twonorm = dfg.getLpNorm(2);
  if (std::all_of(boundary.begin(), boundary.end(), [](bool i){return i;})){
    // check that InnerNodalBasisFunctionIntegral is correct
    double numFullInnerBasisFcns = 1.;
    for (DimType d = 0; d < dim; ++d) {
      numFullInnerBasisFcns *= static_cast<double>(dfg.length(d) - 1);
    }
    BOOST_CHECK_CLOSE(dfg.getInnerNodalBasisFunctionIntegral(), 1./numFullInnerBasisFcns, TestHelper::tolerance);

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
  IndexVector subarrayExtents;
  for (DimType d = 0; d < dim; ++d) {
    auto ghostLayer = dfg.exchangeGhostLayerUpward(d, subarrayExtents);
    IndexVector offsets(dim);
    IndexType numElements = 1;
    for (DimType j = 0; j < dim; ++j) {
      offsets[j] = numElements;
      numElements *= subarrayExtents[j];
    }
    BOOST_CHECK_EQUAL(ghostLayer.size(), numElements);
    if(numElements > 0){
      BOOST_CHECK_EQUAL(subarrayExtents[d], 1);
      for (size_t j=0; j < numElements; ++j){
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
        std::vector<double> coords(dim);
        dfg.getCoordsGlobal(dfg.getGlobalLinearIndex(globAxisIndex), coords);
        BOOST_CHECK_EQUAL(ghostLayer[j], f(coords));
      }
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
  if (std::all_of(boundary.begin(), boundary.end(), [](bool b) { return b; })) {
    // assert that all the corners values are the same now
    auto cornersValues = dfgCopy.getCornersValues();
    for (const auto& cornerValue : cornersValues) {
      BOOST_TEST(cornerValue == formerBoundaryAvgValue);
    }
  }

  // test lower to upper exchange
  for (DimType d = 0; d < dim; ++d) {
    // std::cout << dfg << std::endl;
    if (boundary[d] == true){
      dfg.writeLowerBoundaryToUpperBoundary(d);

      for (IndexType i = 0; i < dfg.getNrLocalElements(); ++i){
        std::vector<double> coords(dim);
        dfg.getCoordsLocal(i, coords);
        if (std::abs((coords[d] - 1.)) < TestHelper::tolerance){
          bool edge = false;
          // if other dimensions are at maximum too, we are at an edge
          // results may differ; skip
          for (DimType d_i = 0; d_i < dim; ++d_i){
            if (d_i != d && std::abs((coords[d_i] - 1.)) < TestHelper::tolerance){
              edge = true;
            }
          }
          if (!edge){
            // check if the value is the same as on the lower boundary
            auto compareCoords = coords;
            compareCoords[d] = 0.;
            BOOST_CHECK_EQUAL(dfg.getElementVector()[i], f(compareCoords));
          }
        } else {
          bool otherBoundary = false;
          for (DimType d_i = 0; d_i < dim; ++d_i){
            if (d_i != d && std::abs((coords[d_i] - 1.)) < TestHelper::tolerance){
              otherBoundary = true;
            }
          }
          if (!otherBoundary){
            // make sure all other values remained the same
            BOOST_CHECK_EQUAL(dfg.getElementVector()[i], f(coords));
          }
        }
      }
    }
  }
  // std::cout << dfg << std::endl;
}

BOOST_FIXTURE_TEST_SUITE(distributedfullgrid, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(60))

// with boundary
// isotropic

BOOST_AUTO_TEST_CASE(test_minus3) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {1, 2};
  IndexVector procs = {1, 1};
  std::vector<bool> boundary(2, true);
  checkDistributedFullgrid(levels, procs, boundary, 1);
}
BOOST_AUTO_TEST_CASE(test_minus2) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {2, 3};
  IndexVector procs = {1, 1};
  std::vector<bool> boundary(2, true);
  checkDistributedFullgrid(levels, procs, boundary, 1);
}
BOOST_AUTO_TEST_CASE(test_minus1) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {1, 1, 2};
  IndexVector procs = {1, 1, 1};
  std::vector<bool> boundary(3, true);
  checkDistributedFullgrid(levels, procs, boundary, 1);
}
BOOST_AUTO_TEST_CASE(test_0) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  LevelVector levels = {2, 3};
  IndexVector procs = {1, 1};
  std::vector<bool> boundary(2, true);
  checkDistributedFullgrid(levels, procs, boundary, 1);
}
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
BOOST_AUTO_TEST_SUITE_END()
