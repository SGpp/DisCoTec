#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/test/unit_test.hpp>
#include <complex>
#include <cstdarg>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

#include "TaskConstParaboloid.hpp"
#include "combicom/CombiCom.hpp"
#include "combischeme/CombiMinMaxScheme.hpp"
#include "fullgrid/FullGrid.hpp"
#include "manager/CombiParameters.hpp"
#include "sparsegrid/DistributedSparseGridIO.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
#include "sparsegrid/SGrid.hpp"
#include "test_helper.hpp"
#include "utils/DecompositionUtils.hpp"
#include "utils/IndexVector.hpp"
#include "utils/LevelSetUtils.hpp"
#include "utils/MonteCarlo.hpp"
#include "utils/Types.hpp"

using namespace combigrid;

void checkDistributedSparsegrid(LevelVector& lmin, LevelVector& lmax, std::vector<int>& procs,
                                std::vector<BoundaryType>& boundary, int size) {
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    auto rank = TestHelper::getRank(comm);
    if (rank == 0) {
      std::cout << "test distributedsparsegrid " << lmin << lmax << procs << std::endl;
    }

    assert(lmin.size() == lmax.size());
    const auto dim = static_cast<DimType>(lmin.size());
    auto ldiff = lmax - lmin;
    bool schemeIsRegular = std::adjacent_find(ldiff.begin(), ldiff.end(),
                                              std::not_equal_to<LevelType>()) == ldiff.end();
    std::vector<LevelVector> cornersOfScheme(dim, lmin);
    for (DimType d = 0; d < dim; ++d) {
      cornersOfScheme[d][d] = lmax[d];
    }

    // create with "own" constructor
    auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<std::complex<double>>>(
        new DistributedSparseGridUniform<std::complex<double>>(dim, lmax, lmin, comm));

    for (const auto& corner : cornersOfScheme) {
      // make sure corners are part of the scheme
      BOOST_CHECK(std::find(uniDSG->getAllLevelVectors().begin(),
                            uniDSG->getAllLevelVectors().end(),
                            corner) != uniDSG->getAllLevelVectors().end());
      // and higher neighbors of corner are not part of the scheme
      for (DimType d = 0; d < dim; ++d) {
        auto neighbor = corner;
        neighbor[d] += 1;
        std::stringstream stringStream;
        stringStream << "corner: " << corner << " neighbor: " << neighbor;
        BOOST_TEST_CONTEXT(stringStream.str());
        if (schemeIsRegular) {
          BOOST_CHECK(std::find(uniDSG->getAllLevelVectors().begin(),
                                uniDSG->getAllLevelVectors().end(),
                                neighbor) == uniDSG->getAllLevelVectors().end());

        } else {
          BOOST_WARN(std::find(uniDSG->getAllLevelVectors().begin(),
                               uniDSG->getAllLevelVectors().end(),
                               neighbor) == uniDSG->getAllLevelVectors().end());
        }
      }
    }

    // get all subspaces in the (optimized) combischeme
    SGrid<real> sg(dim, lmax, lmin, boundary);
    std::vector<LevelVector> subspaces;
    for (size_t ssID = 0; ssID < sg.getSize(); ++ssID) {
      const LevelVector& ss = sg.getLevelVector(ssID);
      subspaces.push_back(ss);
    }

    // compare to subspace constructor
    auto uniDSGfromSubspaces = std::unique_ptr<DistributedSparseGridUniform<std::complex<double>>>(
        new DistributedSparseGridUniform<std::complex<double>>(dim, subspaces, comm));

    BOOST_CHECK_EQUAL(subspaces.size(), uniDSGfromSubspaces->getNumSubspaces());
    BOOST_CHECK_EQUAL(uniDSG->getNumSubspaces(), uniDSGfromSubspaces->getNumSubspaces());

    for (decltype(uniDSG->getNumSubspaces()) i = 0; i < uniDSG->getNumSubspaces(); ++i) {
      BOOST_CHECK_EQUAL(0, uniDSG->getDataSize(i));
      BOOST_CHECK_EQUAL(0, uniDSGfromSubspaces->getDataSize(i));
      for (DimType d = 0; d < dim; ++d)
        BOOST_CHECK_EQUAL(uniDSG->getLevelVector(i)[d], uniDSGfromSubspaces->getLevelVector(i)[d]);
    }

    // use decomposition for full grids
    std::vector<IndexVector> decomposition;
    auto lref = lmax + lmax;
    auto procsRef = procs;
    decomposition = combigrid::getStandardDecomposition(lref, procsRef);

    // make sure that registered DFGs set the sizes right
    auto dfgLevel = lmin;
    dfgLevel[0] = lmax[0] + 1;
    auto dfgDecomposition =
        combigrid::downsampleDecomposition(decomposition, lref, dfgLevel, boundary);

    auto uniDFG = std::unique_ptr<OwningDistributedFullGrid<std::complex<double>>>(
        new OwningDistributedFullGrid<std::complex<double>>(dim, dfgLevel, comm, boundary, procs,
                                                            true, dfgDecomposition));

    uniDSG->registerDistributedFullGrid(*uniDFG);

    for (decltype(uniDSG->getNumSubspaces()) i = 0; i < uniDSG->getNumSubspaces(); ++i) {
      const auto& level = uniDSG->getLevelVector(i);
      BOOST_ASSERT(lmin.size() > 1);
      auto secondDimLevel = level[1];
      if (secondDimLevel > lmin[1]) {
        BOOST_CHECK(secondDimLevel <= lmax[1]);
        BOOST_CHECK_EQUAL(uniDSG->getDataSize(i), 0);
        BOOST_CHECK_EQUAL(uniDFG->getFGPointsOfSubspace(level).size(), 0);
      } else {
        BOOST_TEST_CONTEXT(std::to_string(i));
        BOOST_CHECK_EQUAL(uniDSG->getDataSize(i), uniDFG->getFGPointsOfSubspace(level).size());
      }
    }

    // set function values on dfg
    ParaboloidFn<CombiDataType> f;
    for (IndexType li = 0; li < uniDFG->getNrLocalElements(); ++li) {
      std::vector<double> coords(dim);
      uniDFG->getCoordsLocal(li, coords);
      uniDFG->getData()[li] = f(coords);
    }
    // std::cout << *uniDFG << std::endl;

    // make sure that registered DFGs set the right values
    // attn: usually, the dfgs are hierarchized before adding to the dsg for combination in
    // hierarchical space (without interpolation on coarser component dfgs) here we use the nodal
    // values for testing purposes
    BOOST_TEST_CHECKPOINT("Add to uniform SG");
    // create subspace data
    uniDSG->createSubspaceData();
    uniDSG->setZero();
    uniDSG->addDistributedFullGrid(*uniDFG, 1.);

    BOOST_TEST_CHECKPOINT("Add to uniform SG from subspaces");
    uniDSGfromSubspaces->registerDistributedFullGrid(*uniDFG);
    BOOST_CHECK_EQUAL(0, uniDSGfromSubspaces->getRawDataSize());
    BOOST_CHECK_GT(uniDSG->getRawDataSize(), uniDSGfromSubspaces->getRawDataSize());
    uniDSGfromSubspaces->createSubspaceData();
    uniDSGfromSubspaces->setZero();
    uniDSGfromSubspaces->addDistributedFullGrid(*uniDFG, 0.);
    BOOST_CHECK_EQUAL(uniDSG->getRawDataSize(), uniDSGfromSubspaces->getRawDataSize());
    for (decltype(uniDSGfromSubspaces->getNumSubspaces()) i = 0;
         i < uniDSGfromSubspaces->getNumSubspaces(); ++i) {
      BOOST_CHECK_EQUAL(uniDSG->getDataSize(i), uniDSGfromSubspaces->getDataSize(i));
    }
    for (size_t i = 0; i < uniDSGfromSubspaces->getRawDataSize(); ++i) {
      BOOST_TEST_CONTEXT(std::to_string(i))
      BOOST_CHECK_EQUAL(uniDSGfromSubspaces->getRawData()[i], 0.);
    }

    auto writeSuccess = DistributedSparseGridIO::writeOneFile(*uniDSG, "test_sg_all", true);
    BOOST_WARN(writeSuccess > 0);
    if (writeSuccess) {
#ifndef DISCOTEC_USE_LZ4
      BOOST_CHECK_EQUAL(writeSuccess, uniDSGfromSubspaces->getRawDataSize());
#endif  // not defined DISCOTEC_USE_LZ4
      for (size_t i = 0; i < uniDSGfromSubspaces->getRawDataSize(); ++i) {
        uniDSGfromSubspaces->getRawData()[i] = -1000.;
      }
      auto readSuccess = DistributedSparseGridIO::readOneFile(*uniDSGfromSubspaces, "test_sg_all");
      BOOST_CHECK_EQUAL(readSuccess, uniDSGfromSubspaces->getRawDataSize());
      if (readSuccess) {
        BOOST_TEST_CHECKPOINT("compare values");
        for (size_t i = 0; i < uniDSG->getRawDataSize(); ++i) {
          BOOST_TEST_CONTEXT(std::to_string(i))
          BOOST_CHECK_EQUAL(uniDSG->getRawData()[i], uniDSGfromSubspaces->getRawData()[i]);
        }
      }
      auto reduceSuccess =
          DistributedSparseGridIO::readOneFileAndReduce(*uniDSGfromSubspaces, "test_sg_all", 1);
      BOOST_CHECK_EQUAL(reduceSuccess, uniDSGfromSubspaces->getRawDataSize());
      if (readSuccess) {
        BOOST_TEST_CHECKPOINT("compare double values");
        for (size_t i = 0; i < uniDSG->getRawDataSize(); ++i) {
          BOOST_TEST_CONTEXT(std::to_string(i))
          BOOST_CHECK_EQUAL(uniDSG->getRawData()[i] + uniDSG->getRawData()[i],
                            uniDSGfromSubspaces->getRawData()[i]);
        }
      }
      uniDSGfromSubspaces->copyDataFrom(*uniDSG);
      BOOST_TEST_CHECKPOINT("compare values after copy");
      for (size_t i = 0; i < uniDSG->getRawDataSize(); ++i) {
        BOOST_TEST_CONTEXT(std::to_string(i))
        BOOST_CHECK_EQUAL(uniDSG->getRawData()[i], uniDSGfromSubspaces->getRawData()[i]);
      }
    }

    BOOST_TEST_CHECKPOINT("write to disk chunked");
    // test for dumping sparse grid data to disk and reading back in later (allow time for file
    // system)
    DistributedSparseGridIO::writeToDiskChunked(*uniDSG, "test_sg_");

    BOOST_TEST_CHECKPOINT("write sparse min/max");
    // make sure that right min/max values are written
    DistributedSparseGridIO::writeMinMaxCoefficents(
        *uniDSG,
        "sparse_paraboloid_minmax_" + std::to_string(dim) + "D_" + std::to_string(size) + "_" +
            std::to_string(boundary[0]),
        0, uniDSG->getCommunicator());
    // and remove straight away
    if (rank == 0) {
      auto status = system("rm sparse_paraboloid_minmax_*");
      BOOST_CHECK_GE(status, 0);
    }

    // std::cout << *uniDSG << std::endl;

    // add large full grid to sparse grid, to fill all the possible subspaces
    BOOST_TEST_CHECKPOINT("create large full grid");
    dfgLevel = lmax;
    for (auto& l : dfgLevel) {
      l += 1;
    }

    // have a tiny delay here, by already allocating dfg
    dfgDecomposition = combigrid::downsampleDecomposition(decomposition, lref, dfgLevel, boundary);
    auto largeUniDFG = std::unique_ptr<OwningDistributedFullGrid<std::complex<double>>>(
        new OwningDistributedFullGrid<std::complex<double>>(dim, dfgLevel, comm, boundary, procs,
                                                            true, dfgDecomposition));

    BOOST_TEST_CHECKPOINT("read from disk chunked");
    uniDSGfromSubspaces->setZero();
    DistributedSparseGridIO::readFromDiskChunked(*uniDSGfromSubspaces, "test_sg_");
    BOOST_TEST_CHECKPOINT("compare values chunked");
    for (size_t i = 0; i < uniDSG->getRawDataSize(); ++i) {
      BOOST_TEST_CONTEXT(std::to_string(i));
      BOOST_CHECK_EQUAL(uniDSG->getRawData()[i], uniDSGfromSubspaces->getRawData()[i]);
    }

    // and remove straight away
    MPI_Barrier(comm);
    if (rank == 0) {
      auto status = system("rm test_sg_*");
      BOOST_CHECK_GE(status, 0);
    }

    BOOST_TEST_CHECKPOINT("Register to uniform SG");
    uniDSG->registerDistributedFullGrid(
        *largeUniDFG);  // TODO create levels and actually test something

    // check if the sizes set are actually the ones we calculate with CombiMinMaxScheme
    BOOST_TEST_CHECKPOINT("check subspace sizes");
    auto subspacesDataSizes = uniDSG->getSubspaceDataSizes();
    size_t numDataPointsHere =
        std::accumulate(subspacesDataSizes.begin(), subspacesDataSizes.end(), 0);
    BOOST_CHECK(numDataPointsHere > 0);

    if (std::all_of(boundary.begin(), boundary.end(), [](BoundaryType i) { return i == 2; })) {
      auto newLmin = lmin;
      auto newLmax = lmax;
      auto newLref = lref;
      // I think this flag may be the wrong way around...?
      if (!reverseOrderingDFGPartitions) {
        std::reverse(newLmin.begin(), newLmin.end());
        std::reverse(newLmax.begin(), newLmax.end());
        std::reverse(newLref.begin(), newLref.end());
        // std::reverse(procsRef.begin(), procsRef.end());
        std::reverse(decomposition.begin(), decomposition.end());
      }
      BOOST_TEST_CHECKPOINT("get partitioned num dofs");
      auto partitionedNumDOFs =
          getPartitionedNumDOFSGAdaptive(newLmin, newLmax, newLref, decomposition);

      BOOST_TEST_CHECKPOINT("compare partitioned num dofs");
      auto myNumDOFs = partitionedNumDOFs[rank];
      // std::cout << "partitionedNumDOFs" << partitionedNumDOFs << std::endl;
      if (schemeIsRegular) {
        BOOST_CHECK_EQUAL(myNumDOFs, numDataPointsHere);
      } else {
        // TODO better match non-regular schemes!
        BOOST_WARN_EQUAL(myNumDOFs, numDataPointsHere);
      }

      if (rank == 0) {
        auto sumDOFPartitioned =
            std::accumulate(partitionedNumDOFs.begin(), partitionedNumDOFs.end(), 0);
        auto sgDOF = printSGDegreesOfFreedomAdaptive(newLmin, newLmax);
        BOOST_CHECK_EQUAL(sgDOF, sumDOFPartitioned);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE(distributedsparsegrid, *boost::unit_test::timeout(2500))
// very cheap
BOOST_AUTO_TEST_CASE(test_0) {
  LevelVector lmin = {1, 1};
  LevelVector lmax = {3, 3};
  for (BoundaryType bValue : std::vector<BoundaryType>({0, 1, 2})) {
    std::vector<int> procs = {1, 1};
    std::vector<BoundaryType> boundary(2, bValue);
    auto multProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
    checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_CASE(test_1) {
  LevelVector lmin = {3, 3};
  LevelVector lmax = {7, 7};
  for (int procOne : {1, 2, 3}) {
    for (int procTwo : {1, 2}) {
      BoundaryType bValue = 2;
      std::vector<BoundaryType> boundary(2, bValue);
      std::vector<int> procs = {procOne, procTwo};
      auto multProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
      BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
      checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}
// same for anisotropic
BOOST_AUTO_TEST_CASE(test_2) {
  LevelVector lmin = {2, 4};
  LevelVector lmax = {6, 8};
  for (int procOne : {1, 2, 3}) {
    for (int procTwo : {1, 2}) {
      BoundaryType bValue = 2;
      std::vector<BoundaryType> boundary(2, bValue);
      std::vector<int> procs = {procOne, procTwo};
      auto multProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
      BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
      checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}
// and non-regular truncated CT
BOOST_AUTO_TEST_CASE(test_3) {
  LevelVector lmin = {2, 4};
  LevelVector lmax = {9, 9};
  for (int procOne : {1, 2, 3}) {
    for (int procTwo : {1, 2}) {
      for (BoundaryType bValue : std::vector<BoundaryType>({0, 1, 2})) {
        std::vector<int> procs = {procOne, procTwo};
        std::vector<BoundaryType> boundary(2, bValue);
        auto multProcs =
            std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
        BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
        checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_4) {
  LevelVector lmin = {3, 3, 3};
  LevelVector lmax = {7, 7, 7};
  for (int procOne : {1, 2, 3}) {
    for (int procTwo : {1, 2}) {
      BoundaryType bValue = 2;
      std::vector<BoundaryType> boundary(3, bValue);
      std::vector<int> procs = {procOne, procTwo, 1};
      auto multProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
      BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
      checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}
// same for anisotropic
BOOST_AUTO_TEST_CASE(test_5) {
  LevelVector lmin = {2, 3, 4};
  LevelVector lmax = {6, 7, 8};
  for (int procOne : {1, 2, 3}) {
    for (int procTwo : {1, 2}) {
      BoundaryType bValue = 2;
      std::vector<BoundaryType> boundary(3, bValue);
      std::vector<int> procs = {procOne, procTwo, 1};
      auto multProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
      BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
      checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}
// non-regular truncated CT
BOOST_AUTO_TEST_CASE(test_6) {
  LevelVector lmin = {2, 3, 4};
  LevelVector lmax = {7, 7, 7};
  for (int procOne : {1, 2, 3}) {
    for (int procTwo : {1, 2}) {
      for (BoundaryType bValue : std::vector<BoundaryType>({0, 1, 2})) {
        std::vector<int> procs = {procOne, procTwo, 1};
        std::vector<BoundaryType> boundary(3, bValue);
        auto multProcs =
            std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
        BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
        checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
}

// 4D anisotropic but regular
BOOST_AUTO_TEST_CASE(test_7) {
  LevelVector lmin = {2, 3, 1, 1};
  LevelVector lmax = {6, 7, 5, 5};
  for (int procOne : {1, 3}) {
    for (int procTwo : {1, 2}) {
      for (BoundaryType bValue : std::vector<BoundaryType>({0, 1, 2})) {
        std::vector<int> procs = {procOne, procTwo, 1, 1};
        std::vector<BoundaryType> boundary(4, bValue);
        auto multProcs =
            std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
        BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
        checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
}

// 6D anisotropic but regular
BOOST_AUTO_TEST_CASE(test_8) {
  LevelVector lmin = {1, 1, 2, 1, 2, 1};
  LevelVector lmax = {3, 3, 4, 3, 4, 3};
  for (int procOne : {1, 3}) {
    for (int procTwo : {1, 2}) {
      for (BoundaryType bValue : std::vector<BoundaryType>({0, 1, 2})) {
        std::vector<int> procs = {1, 1, procOne, 1, procTwo, 1};
        std::vector<BoundaryType> boundary(6, bValue);
        auto multProcs =
            std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
        BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
        checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_getPartitionedNumDOFSGAdaptive_1) {
  // 1D case
  LevelVector lmin = {1};
  LevelVector lmax = {5};
  std::vector<int> procs = {2};
  std::vector<IndexVector> decomposition = {{0, 1}};
  auto partitionedNumDOFs = getPartitionedNumDOFSGAdaptive(lmin, lmax, lmax, decomposition);
  auto sln = LevelVector({1, 32});
  BOOST_CHECK_EQUAL_COLLECTIONS(partitionedNumDOFs.begin(), partitionedNumDOFs.end(), sln.begin(),
                                sln.end());
  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    auto sumDOFPartitioned =
        std::accumulate(partitionedNumDOFs.begin(), partitionedNumDOFs.end(), 0);
    auto sgDOF = printSGDegreesOfFreedomAdaptive(lmin, lmax);
    BOOST_CHECK_EQUAL(sgDOF, sumDOFPartitioned);
  }
}
BOOST_AUTO_TEST_CASE(test_getPartitionedNumDOFSGAdaptive_2) {
  // 2D case
  LevelVector lmin = {1, 1};
  LevelVector lmax = {2, 2};
  std::vector<int> procs = {2, 2};
  std::vector<IndexVector> decomposition = {{0, 1}, {0, 3}};
  auto partitionedNumDOFs = getPartitionedNumDOFSGAdaptive(lmin, lmax, lmax, decomposition);
  auto sln = LevelVector({3, 10, 2, 6});
  BOOST_CHECK_EQUAL_COLLECTIONS(partitionedNumDOFs.begin(), partitionedNumDOFs.end(), sln.begin(),
                                sln.end());
  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    auto sumDOFPartitioned =
        std::accumulate(partitionedNumDOFs.begin(), partitionedNumDOFs.end(), 0);
    auto sgDOF = printSGDegreesOfFreedomAdaptive(lmin, lmax);
    BOOST_CHECK_EQUAL(sgDOF, sumDOFPartitioned);
  }
}
BOOST_AUTO_TEST_CASE(test_getPartitionedNumDOFSGAdaptive_3) {
  // 3D case (still OK to visualize)
  LevelVector lmin = {1, 1, 1};
  LevelVector lmax = {3, 3, 3};
  std::vector<int> procs = {2, 2, 2};
  std::vector<IndexVector> decomposition = {{0, 1}, {0, 2}, {0, 3}};
  auto partitionedNumDOFs = getPartitionedNumDOFSGAdaptive(lmin, lmax, lmax, decomposition);
  auto sln = LevelVector({4, 16, 13, 46, 8, 30, 24, 84});
  BOOST_CHECK_EQUAL_COLLECTIONS(partitionedNumDOFs.begin(), partitionedNumDOFs.end(), sln.begin(),
                                sln.end());
  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    // std::cout << partitionedNumDOFs << std::endl;
    auto sumDOFPartitioned =
        std::accumulate(partitionedNumDOFs.begin(), partitionedNumDOFs.end(), 0);
    auto sgDOF = printSGDegreesOfFreedomAdaptive(lmin, lmax);
    BOOST_CHECK_EQUAL(sgDOF, sumDOFPartitioned);
  }
}

BOOST_AUTO_TEST_CASE(test_createTruncatedHierarchicalLevels) {
  LevelVector lmin = {2, 2, 2, 2};
  LevelVector lmax = {4, 4, 4, 4};
  auto dim = static_cast<DimType>(lmin.size());
  std::vector<LevelVector> created;
  combigrid::createTruncatedHierarchicalLevels(lmax, lmin, created);
  // std::cout << lmin << lmax << std::endl;
  // std::cout << "created.size() = " << created.size() << std::endl;
  // std::cout << created << std::endl;

  std::vector<LevelVector> cornersOfScheme(dim, lmin);
  for (DimType d = 0; d < dim; ++d) {
    cornersOfScheme[d][d] = lmax[d];
  }
  for (const auto& corner : cornersOfScheme) {
    // make sure corners are part of the scheme
    BOOST_CHECK(std::find(created.begin(), created.end(), corner) != created.end());
    // and higher neighbors of corner are not part of the scheme
    for (DimType d = 0; d < dim; ++d) {
      auto neighbor = corner;
      neighbor[d] += 1;
      BOOST_CHECK(std::find(created.begin(), created.end(), neighbor) == created.end());
    }
  }
  LevelType largestCornerSum = 0;
  for (const auto& corner : cornersOfScheme) {
    auto cornerSum = std::accumulate(corner.begin(), corner.end(), static_cast<LevelType>(0));
    largestCornerSum = std::max(largestCornerSum, cornerSum);
  }
  for (const auto& level : created) {
    for (DimType d = 0; d < dim; ++d) {
      BOOST_CHECK(level[d] <= lmax[d]);
    }
    auto levelSum = std::accumulate(level.begin(), level.end(), 0);
    BOOST_CHECK(levelSum <= largestCornerSum);
  }
  auto it = std::unique(created.begin(), created.end());
  BOOST_CHECK(it == created.end());
  BOOST_CHECK(std::is_sorted(created.begin(), created.end()));
}

BOOST_AUTO_TEST_CASE(test_createSubspacesSingleLevel) {
  LevelVector lmin = {4, 3, 4, 5};
  LevelVector lmax = lmin;
  std::vector<LevelVector> created;
  combigrid::createTruncatedHierarchicalLevels(lmax, lmin, created);
  auto downSet = getDownSet(lmax);
  // BOOST_CHECK_EQUAL_COLLECTIONS(downSet.begin(), downSet.end(), created.begin(), created.end());
  BOOST_CHECK_EQUAL(downSet.size(), created.size());
  for (const auto& level : downSet) {
    BOOST_CHECK(std::find(created.begin(), created.end(), level) != created.end());
  }
  auto it = std::unique(created.begin(), created.end());
  BOOST_CHECK(it == created.end());
  BOOST_CHECK_EQUAL(created.size(),
                    std::accumulate(lmax.begin(), lmax.end(), 1, std::multiplies<LevelType>()));
  for (BoundaryType boundary : std::vector<BoundaryType>({0, 1, 2})) {
    std::vector<BoundaryType> boundaryVector = {boundary, boundary, boundary,
                                                static_cast<BoundaryType>((boundary == 0) ? 2 : 0)};
    BOOST_CHECK_EQUAL(combigrid::getNumDofNodal(lmax, boundaryVector),
                      getNumDofHierarchical(downSet, boundaryVector));
    BOOST_CHECK_EQUAL(combigrid::getNumDofNodal(lmax, boundaryVector),
                      getNumDofHierarchical(created, boundaryVector));
  }
  BOOST_CHECK(std::is_sorted(downSet.begin(), downSet.end()));
  BOOST_CHECK(std::is_sorted(created.begin(), created.end()));
}

BOOST_AUTO_TEST_CASE(test_createTruncatedHierarchicalLevels_large) {
  LevelVector lmin = {2, 2, 2, 2, 2, 2};
  LevelVector lmax = {19, 19, 19, 19, 19, 19};
  std::vector<LevelVector> created;
  // create once to fill the cache
  combigrid::createTruncatedHierarchicalLevels(lmax, lmin, created);
  created.clear();

  auto start = std::chrono::high_resolution_clock::now();
  combigrid::createTruncatedHierarchicalLevels(lmax, lmin, created);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  BOOST_TEST_MESSAGE("time to create regular hierarchical levels: " << duration.count()
                                                                    << " milliseconds");
  BOOST_TEST_MESSAGE("number of levels created: " << created.size());
#ifdef NDEBUG
  BOOST_CHECK(duration.count() < 10000);
#endif
  auto it = std::unique(created.begin(), created.end());
  BOOST_CHECK(it == created.end());
  BOOST_CHECK(std::is_sorted(created.begin(), created.end()));
}

BOOST_AUTO_TEST_CASE(test_createSubspacesSingleLevel_large) {
  LevelVector lmin = {5, 5, 5, 5, 5, 4};
  LevelVector lmax = lmin;
  std::vector<LevelVector> created;
  // create once to fill the cache
  combigrid::createTruncatedHierarchicalLevels(lmax, lmin, created);
  created.clear();

  auto start = std::chrono::high_resolution_clock::now();
  combigrid::createTruncatedHierarchicalLevels(lmax, lmin, created);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  BOOST_TEST_MESSAGE("time to create block of hierarchical levels: " << duration.count()
                                                                     << " milliseconds");
  BOOST_TEST_MESSAGE("number of levels created: " << created.size());
#ifdef NDEBUG
  BOOST_CHECK(duration.count() < 1000);
#endif
  std::sort(created.begin(), created.end());
  auto it = std::unique(created.begin(), created.end());
  BOOST_CHECK(it == created.end());
  LevelVector allOnes = {1, 1, 1, 1, 1, 1};
  BOOST_CHECK_EQUAL_COLLECTIONS(created[0].begin(), created[0].end(), allOnes.begin(),
                                allOnes.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(created.back().begin(), created.back().end(), lmax.begin(),
                                lmax.end());
  BOOST_CHECK_EQUAL(created.size(),
                    std::accumulate(lmax.begin(), lmax.end(), 1, std::multiplies<LevelType>()));

  start = std::chrono::high_resolution_clock::now();
  auto downSet = getDownSet(lmax);
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  BOOST_TEST_MESSAGE("time to create downward closed set: " << duration.count() << " milliseconds");
  for (const auto& level : created) {
    BOOST_CHECK(level.size() == lmin.size());
  }
  BOOST_CHECK_EQUAL(created.size(), downSet.size());
  for (const auto& level : created) {
    std::stringstream stringStream;
    stringStream << "level: " << level;
    BOOST_TEST_CONTEXT(stringStream.str());
    BOOST_REQUIRE(std::find(downSet.begin(), downSet.end(), level) != downSet.end());
  }
  BOOST_CHECK(std::is_sorted(downSet.begin(), downSet.end()));
  BOOST_CHECK(std::is_sorted(created.begin(), created.end()));

  start = std::chrono::high_resolution_clock::now();
  auto downSetGenerator = HypercubeDownSetGenerator(lmax);
  auto previousFind = downSet.begin();
  static thread_local LevelVector level;
  for (LevelType i = 0; i < downSetGenerator.getTotalNumberOfLevels(); ++i) {
#pragma omp critical
    level = downSetGenerator.getNextLevel();
    BOOST_CHECK(level.size() == lmin.size());
    auto found = std::find(previousFind, downSet.end(), level);
    BOOST_REQUIRE(found != downSet.end());
    previousFind = found;
  }
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  BOOST_TEST_MESSAGE("time to generate and iterate downward closed set: " << duration.count()
                                                                          << " milliseconds");
  BOOST_CHECK_EQUAL(downSetGenerator.getTotalNumberOfLevels(), downSet.size());
}

BOOST_AUTO_TEST_CASE(test_getAllKOutOfDDimensions) {
  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    for (DimType d = 1; d < 8; ++d) {
      for (DimType k = 1; k <= d; ++k) {
        auto allKOutOfD = AllKOutOfDDimensions::get(k, d);
        BOOST_CHECK_EQUAL(allKOutOfD.size(), boost::math::binomial_coefficient<real>(d, k));
        for (const auto& combination : allKOutOfD) {
          BOOST_CHECK_EQUAL(combination.size(), k);
          for (const auto& index : combination) {
            BOOST_CHECK(index < d);
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_reduceSubspaceSizesFileBased) {
  std::vector<int> procs = {4, 1, 2, 1, 1, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = static_cast<DimType>(procs.size());
    LevelVector lmin = {2, 2, 2, 2, 2, 2};
    LevelVector lmax = {11, 11, 11, 11, 11, 11};
    LevelVector lfull = {12, 3, 7, 2, 2, 2};
    std::vector<BoundaryType> boundary(dim, 1);
    auto decomposition = combigrid::getStandardDecomposition(lmax, procs);
    auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<combigrid::DimType>>(
        new DistributedSparseGridUniform<combigrid::DimType>(dim, lmax, lmin, comm));

    {  // register full grid
      auto uniDFG = std::unique_ptr<OwningDistributedFullGrid<combigrid::DimType>>(
          new OwningDistributedFullGrid<combigrid::DimType>(dim, lfull, comm, boundary, procs,
                                                            false));
      uniDSG->registerDistributedFullGrid(*uniDFG);
    }

    auto sizesCopy = uniDSG->getSubspaceDataSizes();

    // write the (smaller) subspaces sizes to disk
    int subspaceWriteSuccess =
        DistributedSparseGridIO::writeSubspaceSizesToFile(*uniDSG, "test_dsg.sizes");
    BOOST_CHECK(subspaceWriteSuccess > 0);

    {  // register reversed full grid
      std::reverse(lfull.begin(), lfull.end());
      auto uniDFG = std::unique_ptr<OwningDistributedFullGrid<combigrid::DimType>>(
          new OwningDistributedFullGrid<combigrid::DimType>(dim, lfull, comm, boundary, procs,
                                                            false));
      uniDSG->registerDistributedFullGrid(*uniDFG);
    }

    auto sizesCopyLarger = uniDSG->getSubspaceDataSizes();

    // read the subspace sizes from disk and do max-reduce
    auto maxFunctionInstantiation = [](SubspaceSizeType a, SubspaceSizeType b) {
      return std::max(a, b);
    };
    int subspaceReduceSuccess = DistributedSparseGridIO::readReduceSubspaceSizesFromFile(
        *uniDSG, "test_dsg.sizes", maxFunctionInstantiation, 2000);
    BOOST_CHECK_EQUAL_COLLECTIONS(sizesCopyLarger.begin(), sizesCopyLarger.end(),
                                  uniDSG->getSubspaceDataSizes().begin(),
                                  uniDSG->getSubspaceDataSizes().end());
    BOOST_CHECK_EQUAL(subspaceReduceSuccess, subspaceWriteSuccess);

    // read the subspace sizes from disk and do min-reduce
    auto minFunctionInstantiation = [](SubspaceSizeType a, SubspaceSizeType b) {
      return std::min(a, b);
    };
    subspaceReduceSuccess = DistributedSparseGridIO::readReduceSubspaceSizesFromFile(
        *uniDSG, "test_dsg.sizes", minFunctionInstantiation, 20);
    BOOST_CHECK_EQUAL_COLLECTIONS(sizesCopy.begin(), sizesCopy.end(),
                                  uniDSG->getSubspaceDataSizes().begin(),
                                  uniDSG->getSubspaceDataSizes().end());
    BOOST_CHECK_EQUAL(subspaceReduceSuccess, subspaceWriteSuccess);
  }
}

BOOST_AUTO_TEST_CASE(test_writeOneFile) {
  std::vector<int> procs = {3, 1, 3, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = static_cast<DimType>(procs.size());
    LevelVector lmin = {2, 2, 2, 2};
    LevelVector lmax = {7, 7, 7, 7};
    std::vector<BoundaryType> boundary(dim, 2);
    auto decomposition = combigrid::getStandardDecomposition(lmax, procs);
    auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<combigrid::real>>(
        new DistributedSparseGridUniform<combigrid::real>(dim, lmax, lmin, comm));
    // register only the top shell of the sparse-grid levels:
    // keeps this test's file size bounded while still exercising write/readOneFile.
    LevelType maxLevelSum = 0;
    for (const auto& level : uniDSG->getAllLevelVectors()) {
      maxLevelSum = std::max(maxLevelSum, levelSum(level));
    }
    for (const auto& level : uniDSG->getAllLevelVectors()) {
      if (levelSum(level) >= maxLevelSum) {
        auto dfgDecomposition =
            combigrid::downsampleDecomposition(decomposition, lmax, level, boundary);
        auto uniDFG = std::unique_ptr<OwningDistributedFullGrid<combigrid::real>>(
            new OwningDistributedFullGrid<combigrid::real>(dim, level, comm, boundary, procs, true,
                                                           dfgDecomposition));
        uniDSG->registerDistributedFullGrid(*uniDFG);
      }
    }
    uniDSG->createSubspaceData();
    uniDSG->setZero();
    size_t totalNumPoints = uniDSG->getRawDataSize();
    MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<size_t>());
    MPI_Allreduce(MPI_IN_PLACE, &totalNumPoints, 1, dtype, MPI_SUM, comm);
    BOOST_TEST_MESSAGE("test_writeOneFile total points: " << totalNumPoints);
    BOOST_CHECK_EQUAL(totalNumPoints, 222209);
    MPI_Barrier(comm);

    auto start = std::chrono::high_resolution_clock::now();
    auto writeSuccess = DistributedSparseGridIO::writeOneFile(*uniDSG, "test_sg_timing", true);
    auto end = std::chrono::high_resolution_clock::now();
    BOOST_CHECK(writeSuccess);
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to write sparse grid: " << duration.count() << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 120000);
#endif
    MPI_Barrier(comm);

    start = std::chrono::high_resolution_clock::now();
    auto readSuccess = DistributedSparseGridIO::readOneFile(*uniDSG, "test_sg_timing");
    end = std::chrono::high_resolution_clock::now();
    BOOST_CHECK(readSuccess);
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to read sparse grid: " << duration.count() << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 30000);
#endif

    if (TestHelper::getRank(comm) == 0) {
      auto status = system("rm test_sg_timing");
      BOOST_CHECK_GE(status, 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(test_anyDistributedSparseGrid) {
  int nprocs = 9;
  CommunicatorType comm = TestHelper::getComm(nprocs);
  if (comm != MPI_COMM_NULL) {
    AnyDistributedSparseGrid anyDSG(123456, MPI_COMM_SELF);
    BOOST_CHECK_EQUAL(anyDSG.getAccumulatedDataSize(), 0);
    BOOST_CHECK_EQUAL(anyDSG.getNumSubspaces(), 123456);
    std::vector<real> randomNums(anyDSG.getNumSubspaces());
    montecarlo::getNumberSequenceFromSeed(randomNums, TestHelper::getRank(comm));
    for (size_t i = 0; i < anyDSG.getNumSubspaces(); ++i) {
      // in 50% of cases, randomly set the data size
      if (randomNums[i] < 0.5) {
        anyDSG.setDataSize(i, i + 1);
      }
    }
    BOOST_CHECK_GT(anyDSG.getAccumulatedDataSize(), 0);

    // set the subspace map across DSGs
    anyDSG.setSubspaceCommunicators(comm, TestHelper::getRank(comm));

    const auto& subspacesByComm = anyDSG.getSubspacesByCommunicator();

    BOOST_CHECK_EQUAL(subspacesByComm.size(), powerOfTwo[nprocs - 1] - 1);

    // for each communicator, have an allreduce to check if the number of subspaces is the same
    std::vector<MPI_Request> requests(subspacesByComm.size());
    std::vector<size_t> numSubspaces{};
    numSubspaces.reserve(subspacesByComm.size());
    auto mapIterator = subspacesByComm.cbegin();
    for (size_t i = 0; i < subspacesByComm.size(); ++i) {
      auto comm = mapIterator->first;
      numSubspaces.push_back(mapIterator->second.size());
      MPI_Iallreduce(MPI_IN_PLACE, &numSubspaces[i], 1,
                     abstraction::getMPIDatatype(abstraction::getabstractionDataType<size_t>()),
                     MPI_BAND, comm, &requests[i]);
      ++mapIterator;
    }
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);

    // check that it is the same as before
    mapIterator = subspacesByComm.cbegin();
    for (size_t i = 0; i < numSubspaces.size(); ++i) {
      BOOST_CHECK_EQUAL(numSubspaces[i], mapIterator->second.size());
      ++mapIterator;
    }

    // ensure that each subspace occurs only once
    std::set<typename AnyDistributedSparseGrid::SubspaceIndexType> allSubspaces{};
    for (const auto& subspaces : subspacesByComm) {
      for (const auto& subspace : subspaces.second) {
        BOOST_CHECK_EQUAL(allSubspaces.count(subspace), 0);
        allSubspaces.insert(subspace);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_sparseGridAndSubspaceReduce) {
  int nprocs = 9;
  CommunicatorType fullComm = TestHelper::getComm(nprocs);
  if (fullComm != MPI_COMM_NULL) {
    Stats::initialize();
    theMPISystem()->initWorldReusable(fullComm, nprocs, 1, false, true);
    // set up 6D Combi scheme
    constexpr DimType dimensionality = 6;
    const LevelVector lmin = LevelVector{3, 3, 3, 3, 2, 2};
    const LevelVector lmax = LevelVector{5, 5, 5, 5, 4, 4};

    const auto chunkSizes = std::vector<uint32_t>{1, 64};
    for (uint32_t chunkSizePerThreadInMiB : chunkSizes) {
      for (CombinationVariant variant :
           {CombinationVariant::subspaceReduce, CombinationVariant::outgroupSparseGridReduce,
            CombinationVariant::sparseGridReduce}) {
        std::vector<LevelVector> myLevels;
        {
          CombiMinMaxScheme combischeme(dimensionality, lmin, lmax);
          combischeme.createAdaptiveCombischeme();
          std::vector<LevelVector> levels = combischeme.getCombiSpaces();
          // select "my" levels round-robin
          for (size_t i = 0; i < levels.size(); ++i) {
            if (i % nprocs == TestHelper::getRank(fullComm)) {
              myLevels.push_back(std::move(levels[i]));
            }
          }
        }

        auto myOwnComm = TestHelper::getCommSelfAsCartesian(dimensionality);
        // every rank has its "own" DSG (like a process group of 1)
        auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<combigrid::real>>(
            new DistributedSparseGridUniform<combigrid::real>(dimensionality, lmax, lmin,
                                                              myOwnComm));

        auto boundary = std::vector<BoundaryType>(dimensionality, 1);
        auto procs = std::vector<int>(dimensionality, 1);
        // for each level, create DFG and register to set DSG's subspace sizes
        for (size_t i = 0; i < myLevels.size(); ++i) {
          auto dfg = std::unique_ptr<DistributedFullGrid<combigrid::real>>(
              new DistributedFullGrid<combigrid::real>(dimensionality, myLevels[i], myOwnComm,
                                                       boundary, nullptr, procs, false));
          uniDSG->registerDistributedFullGrid(*dfg);
        }

        if (variant == CombinationVariant::sparseGridReduce) {
          CombiCom::reduceSubspaceSizes(*uniDSG, fullComm);
          // initialize actual data containers
          uniDSG->createSubspaceData();
          uniDSG->setZero();

          // for each subspace in uniDSG, set values to the subspace index
          for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
            auto subspaceStart = uniDSG->getData(i);
            for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(i); ++j) {
              subspaceStart[j] = i;
            }
          }

          MPI_Barrier(fullComm);
          auto start = std::chrono::high_resolution_clock::now();
          CombiCom::distributedGlobalSparseGridReduce(*uniDSG, chunkSizePerThreadInMiB,
                                                      MPI_PROC_NULL, fullComm);
          auto end = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
          BOOST_TEST_MESSAGE("sparse grid reduce time: " + std::to_string(duration.count()));

          // check that the data is correct
          for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
            auto subspaceStart = uniDSG->getData(i);
            for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(i); ++j) {
              BOOST_CHECK_EQUAL(subspaceStart[j], i * nprocs);
            }
          }

          // try again, this time with reduce+broadcast
          // reset actual data containers and MPI datatype mappings
          uniDSG->deleteSubspaceData();
          uniDSG->createSubspaceData();
          uniDSG->setZero();
          // for each subspace in uniDSG, set values to the subspace index
          for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
            auto subspaceStart = uniDSG->getData(i);
            for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(i); ++j) {
              subspaceStart[j] = i;
            }
          }

          MPI_Barrier(fullComm);
          start = std::chrono::high_resolution_clock::now();
          RankType globalReduceRankThatCollects = 0;
          CombiCom::distributedGlobalSparseGridReduce(*uniDSG, chunkSizePerThreadInMiB,
                                                      globalReduceRankThatCollects);
          MPI_Request request;
          CombiCom::asyncBcastDsgData(*uniDSG, globalReduceRankThatCollects,
                                      theMPISystem()->getGlobalReduceComm(), &request);
          MPI_Wait(&request, MPI_STATUS_IGNORE);
          end = std::chrono::high_resolution_clock::now();
          duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
          BOOST_TEST_MESSAGE("sparse grid reduce+broadcast time: " +
                             std::to_string(duration.count()));

          // check that the data is correct
          for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
            auto subspaceStart = uniDSG->getData(i);
            for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(i); ++j) {
              BOOST_CHECK_EQUAL(subspaceStart[j], i * nprocs);
            }
          }

        } else if (variant == CombinationVariant::subspaceReduce) {
          BOOST_CHECK_LT(uniDSG->getAccumulatedDataSize(), 142606336 / 2);

          // set the subspace map across DSGs
          uniDSG->setSubspaceCommunicators(fullComm, TestHelper::getRank(fullComm));
          BOOST_TEST_CHECKPOINT("subspace communicators set");
          // initialize actual data containers and MPI datatype mappings
          uniDSG->createSubspaceData();
          uniDSG->setZero();
          // for each subspace in uniDSG, set values to the subspace index
          for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
            auto subspaceStart = uniDSG->getData(i);
            for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(i); ++j) {
              subspaceStart[j] = i;
            }
          }

          MPI_Barrier(fullComm);
          auto start = std::chrono::high_resolution_clock::now();
          CombiCom::distributedGlobalSubspaceReduce(*uniDSG, chunkSizePerThreadInMiB);
          auto end = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
          BOOST_TEST_MESSAGE("subspace reduce time: " + std::to_string(duration.count()));

          // check that the data is correct
          const auto& subspacesByComm = uniDSG->getSubspacesByCommunicator();
          std::set<typename AnyDistributedSparseGrid::SubspaceIndexType> checkedSubspaces{};
          for (const auto& subspaces : subspacesByComm) {
            auto commSize = combigrid::getCommSize(subspaces.first);
            for (const auto& subspace : subspaces.second) {
              auto subspaceStart = uniDSG->getData(subspace);
              for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(subspace); ++j) {
                BOOST_CHECK_EQUAL(subspaceStart[j], subspace * commSize);
              }
              checkedSubspaces.insert(subspace);
            }
          }
          // check remaining subspaces
          for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
            if (checkedSubspaces.find(i) == checkedSubspaces.end()) {
              auto subspaceStart = uniDSG->getData(i);
              for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(i); ++j) {
                BOOST_CHECK_EQUAL(subspaceStart[j], i);
              }
            }
          }
        } else {
          // set the subspace map across DSGs
          uniDSG->setOutgroupCommunicator(fullComm, TestHelper::getRank(fullComm));
          BOOST_CHECK_EQUAL(uniDSG->getSubspacesByCommunicator().size(), 1);
          // initialize actual data containers and MPI datatype mappings
          uniDSG->createSubspaceData();
          uniDSG->setZero();
          // for each subspace in uniDSG, set values to the subspace index
          for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
            auto subspaceStart = uniDSG->getData(i);
            for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(i); ++j) {
              subspaceStart[j] = i;
            }
          }

          MPI_Barrier(fullComm);
          BOOST_TEST_CHECKPOINT("outgroup reduce");
          auto start = std::chrono::high_resolution_clock::now();
          CombiCom::distributedGlobalSubspaceReduce(*uniDSG, chunkSizePerThreadInMiB);
          auto end = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
          BOOST_TEST_MESSAGE("outgroup reduce time: " + std::to_string(duration.count()));

          // check that the data is correct
          const auto& subspacesByComm = uniDSG->getSubspacesByCommunicator();
          std::set<typename AnyDistributedSparseGrid::SubspaceIndexType> checkedSubspaces{};
          for (const auto& subspaces : subspacesByComm) {
            for (const auto& subspace : subspaces.second) {
              for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(subspace); ++j) {
                // BOOST_CHECK_EQUAL(subspaceStart[j], ??);
              }
              checkedSubspaces.insert(subspace);
            }
          }
          for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
            if (checkedSubspaces.find(i) == checkedSubspaces.end()) {
              auto subspaceStart = uniDSG->getData(i);
              for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(i); ++j) {
                BOOST_CHECK_EQUAL(subspaceStart[j], i);
              }
            }
          }

          // try again, this time with reduce+broadcast
          // reset actual data containers and MPI datatype mappings
          uniDSG->deleteSubspaceData();
          uniDSG->createSubspaceData();
          uniDSG->setZero();
          // for each subspace in uniDSG, set values to the subspace index
          for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
            auto subspaceStart = uniDSG->getData(i);
            for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(i); ++j) {
              subspaceStart[j] = i;
            }
          }

          MPI_Barrier(fullComm);
          BOOST_TEST_CHECKPOINT("outgroup reduce + broadcast");
          start = std::chrono::high_resolution_clock::now();
          CombiCom::distributedGlobalSubspaceReduce(*uniDSG, chunkSizePerThreadInMiB, 0);
          MPI_Request request;
          CombiCom::asyncBcastOutgroupDsgData(*uniDSG, 0, fullComm, &request);
          MPI_Wait(&request, MPI_STATUS_IGNORE);
          end = std::chrono::high_resolution_clock::now();
          duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
          BOOST_TEST_MESSAGE("outgroup reduce+broadcast time: " + std::to_string(duration.count()));

          // check that the data is correct
          checkedSubspaces.clear();
          for (const auto& subspaces : subspacesByComm) {
            for (const auto& subspace : subspaces.second) {
              for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(subspace); ++j) {
                // BOOST_CHECK_EQUAL(subspaceStart[j], ??);
              }
              checkedSubspaces.insert(subspace);
            }
          }
          for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
            if (checkedSubspaces.find(i) == checkedSubspaces.end()) {
              auto subspaceStart = uniDSG->getData(i);
              for (SubspaceSizeType j = 0; j < uniDSG->getDataSize(i); ++j) {
                BOOST_CHECK_EQUAL(subspaceStart[j], i);
              }
            }
          }
        }
        BOOST_CHECK(!TestHelper::testStrayMessages(fullComm));
      }
    }
    Stats::finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
