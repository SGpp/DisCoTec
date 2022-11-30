#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/test/unit_test.hpp>
#include <complex>
#include <cstdarg>
#include <iostream>
#include <random>
#include <vector>

#include "TaskConstParaboloid.hpp"
#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/SGrid.hpp"
#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"
#include "sgpp/distributedcombigrid/utils/LevelSetUtils.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "test_helper.hpp"

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

    auto uniDFG = std::unique_ptr<DistributedFullGrid<std::complex<double>>>(
        new DistributedFullGrid<std::complex<double>>(dim, dfgLevel, comm, boundary, procs, true,
                                                      dfgDecomposition));

    uniDFG->registerUniformSG(*uniDSG);
    uniDFG->registerUniformSG(*uniDSGfromSubspaces);
    BOOST_CHECK_EQUAL(0, uniDSGfromSubspaces->getRawDataSize());

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
    uniDFG->addToUniformSG(*uniDSG, 1.);
    BOOST_CHECK_GT(uniDSG->getRawDataSize(), uniDSGfromSubspaces->getRawDataSize());
    uniDSGfromSubspaces->createSubspaceData();
    uniDFG->addToUniformSG(*uniDSGfromSubspaces, 0.);
    BOOST_CHECK_EQUAL(uniDSG->getRawDataSize(), uniDSGfromSubspaces->getRawDataSize());
    for (decltype(uniDSGfromSubspaces->getNumSubspaces()) i = 0;
         i < uniDSGfromSubspaces->getNumSubspaces(); ++i) {
      BOOST_CHECK_EQUAL(uniDSG->getDataSize(i), uniDSGfromSubspaces->getDataSize(i));
    }
    for (size_t i = 0; i < uniDSGfromSubspaces->getRawDataSize(); ++i) {
      BOOST_TEST_CONTEXT(std::to_string(i))
      BOOST_CHECK_EQUAL(uniDSGfromSubspaces->getRawData()[i], 0.);
    }
    // TODO test if the values are correct

    auto writeSuccess = uniDSG->writeOneFileToDisk("test_sg_all");
    BOOST_WARN(writeSuccess);
    if (writeSuccess) {
      for (size_t i = 0; i < uniDSGfromSubspaces->getRawDataSize(); ++i) {
        uniDSGfromSubspaces->getRawData()[i] = -1000.;
      }
      auto readSuccess = uniDSGfromSubspaces->readOneFileFromDisk("test_sg_all");
      BOOST_CHECK(readSuccess);
      if (readSuccess) {
        BOOST_TEST_CHECKPOINT("compare values");
        for (size_t i = 0; i < uniDSG->getRawDataSize(); ++i) {
          BOOST_TEST_CONTEXT(std::to_string(i))
          BOOST_CHECK_EQUAL(uniDSG->getRawData()[i], uniDSGfromSubspaces->getRawData()[i]);
        }
      }
    }

    BOOST_TEST_CHECKPOINT("write sparse min/max");
    // make sure that right min/max values are written
    uniDSG->writeMinMaxCoefficents("sparse_paraboloid_minmax_" + std::to_string(dim) + "D_" +
                                       std::to_string(size) + "_" + std::to_string(boundary[0]),
                                   0);
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
    dfgDecomposition = combigrid::downsampleDecomposition(decomposition, lref, dfgLevel, boundary);
    auto largeUniDFG = std::unique_ptr<DistributedFullGrid<std::complex<double>>>(
        new DistributedFullGrid<std::complex<double>>(dim, dfgLevel, comm, boundary, procs, true,
                                                      dfgDecomposition));

    BOOST_TEST_CHECKPOINT("Register to uniform SG");
    largeUniDFG->registerUniformSG(*uniDSG);

    // make sure that right min/max values are written %TODO remove file
    BOOST_TEST_CHECKPOINT("write min/max coefficients");
    uniDSG->writeMinMaxCoefficents(
        "sparse_paraboloid_minmax_large_" + std::to_string(dim) + "D_" + std::to_string(size), 0);

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
    // test for dumping sparse grid data to disk and reading back in
    uniDSG->writeToDiskChunked("test_sg_");
    uniDSGfromSubspaces->setZero();
    uniDSGfromSubspaces->readFromDiskChunked("test_sg_");
    BOOST_TEST_CHECKPOINT("compare values chunked");
    for (size_t i = 0; i < uniDSG->getRawDataSize(); ++i) {
      BOOST_TEST_CONTEXT(std::to_string(i));
      BOOST_CHECK_EQUAL(uniDSG->getRawData()[i], uniDSGfromSubspaces->getRawData()[i]);
    }

    // and remove straight away
    if (rank == 0) {
      auto status = system("rm test_sg_*");
      BOOST_CHECK_GE(status, 0);
    }
  }
}

BOOST_AUTO_TEST_SUITE(distributedsparsegrid, *boost::unit_test::timeout(1500))
// very cheap
BOOST_AUTO_TEST_CASE(test_0) {
  LevelVector lmin = {1, 1};
  LevelVector lmax = {3, 3};
  BoundaryType bValue = 2;
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
      auto multProcs =
          std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
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
      auto multProcs =
          std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
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
      auto multProcs =
          std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
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
      auto multProcs =
          std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
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
  LevelVector lmin = {2, 1, 4, 1, 3, 1};
  LevelVector lmax = {4, 3, 6, 3, 5, 3};
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
  for (BoundaryType boundary : {0, 1, 2}) {
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

BOOST_AUTO_TEST_CASE(test_writeOneFileToDisk) {
  std::vector<int> procs = {3, 1, 3, 1, 1, 1};
  CommunicatorType comm = TestHelper::getComm(procs);
  if (comm != MPI_COMM_NULL) {
    DimType dim = static_cast<DimType>(procs.size());
    LevelVector lmin = {2, 2, 2, 2, 2, 2};
    LevelVector lmax = {11, 11, 11, 11, 11, 11};
    std::vector<BoundaryType> boundary(dim, 2);
    auto decomposition = combigrid::getStandardDecomposition(lmax, procs);
    auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<combigrid::real>>(
        new DistributedSparseGridUniform<combigrid::real>(dim, lmax, lmin, comm));
    // iterate main diagonal of combi scheme and register to populate all subspaces
    for (const auto& level : uniDSG->getAllLevelVectors()) {
      if (levelSum(level) == 21) {
        auto dfgDecomposition =
            combigrid::downsampleDecomposition(decomposition, lmax, level, boundary);
        auto uniDFG = std::unique_ptr<DistributedFullGrid<combigrid::real>>(
            new DistributedFullGrid<combigrid::real>(dim, level, comm, boundary, procs, true,
                                                     dfgDecomposition));
        uniDFG->registerUniformSG(*uniDSG);
      }
    }
    uniDSG->setZero();
    size_t totalNumPoints = uniDSG->getRawDataSize();
    MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<size_t>());
    MPI_Allreduce(MPI_IN_PLACE, &totalNumPoints, 1, dtype, MPI_SUM, comm);
    BOOST_CHECK_EQUAL(totalNumPoints, 1050968065);
    MPI_Barrier(comm);

    auto start = std::chrono::high_resolution_clock::now();
    auto writeSuccess = uniDSG->writeOneFileToDisk("test_sg_timing");
    auto end = std::chrono::high_resolution_clock::now();
    BOOST_CHECK(writeSuccess);
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    BOOST_TEST_MESSAGE("time to write sparse grid: " << duration.count() << " milliseconds");
#ifdef NDEBUG
    BOOST_CHECK(duration.count() < 120000);
#endif
    MPI_Barrier(comm);

    start = std::chrono::high_resolution_clock::now();
    auto readSuccess = uniDSG->readOneFileFromDisk("test_sg_timing");
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

BOOST_AUTO_TEST_SUITE_END()
