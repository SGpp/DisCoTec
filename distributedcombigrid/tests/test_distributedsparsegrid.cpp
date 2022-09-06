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

#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "test_helper.hpp"
#include "TaskConstParaboloid.hpp"

using namespace combigrid;

void checkDistributedSparsegrid(LevelVector& lmin, LevelVector& lmax, IndexVector& procs, std::vector<bool>& boundary,
                              int size) {
  CommunicatorType comm = TestHelper::getComm(size);
  if (comm != MPI_COMM_NULL) {
    auto rank = TestHelper::getRank(comm);
    if (rank == 0) {
      std::cout << "test distributedsparsegrid " << lmin << lmax << procs << std::endl;
    }

    assert(lmin.size() == lmax.size());
    const DimType dim = lmin.size();
    auto ldiff = lmax - lmin;
    bool schemeIsRegular = std::adjacent_find(ldiff.begin(), ldiff.end(),
                                              std::not_equal_to<LevelType>()) == ldiff.end();
    std::vector<LevelVector> cornersOfScheme(dim, lmin);
    for (DimType d = 0; d < dim; ++d) {
      cornersOfScheme[d][d] = lmax[d];
    }

    // create with "own" constructor
    auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<std::complex<double>>>(
        new DistributedSparseGridUniform<std::complex<double>>(dim, lmax, lmin, boundary, comm));

    for (const auto& corner : cornersOfScheme) {
      // make sure corners are part of the scheme
      BOOST_CHECK(std::find(uniDSG->getAllLevelVectors().begin(),
                            uniDSG->getAllLevelVectors().end(),
                            corner) != uniDSG->getAllLevelVectors().end());
      // and higher neighbors of corner are not part of the scheme
      for (DimType d = 0; d < dim; ++d) {
        auto neighbor = corner;
        neighbor[d] += 1;
        BOOST_WARN(std::find(uniDSG->getAllLevelVectors().begin(),
                             uniDSG->getAllLevelVectors().end(),
                             neighbor) == uniDSG->getAllLevelVectors().end());
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
        new DistributedSparseGridUniform<std::complex<double>>(dim, subspaces, boundary, comm));

    BOOST_CHECK_EQUAL(subspaces.size(), uniDSGfromSubspaces->getNumSubspaces());
    BOOST_CHECK_EQUAL(uniDSG->getNumSubspaces(), uniDSGfromSubspaces->getNumSubspaces());

    for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
      BOOST_CHECK_EQUAL(0, uniDSG->getDataSize(i));
      BOOST_CHECK_EQUAL(0, uniDSGfromSubspaces->getDataSize(i));
      for (DimType d = 0; d < dim; ++d)
        BOOST_CHECK_EQUAL(uniDSG->getLevelVector(i)[d], uniDSGfromSubspaces->getLevelVector(i)[d]);
    }

    // use decomposition for full grids
    std::vector<LevelVector> decomposition;
    auto lref = lmax + lmax;
    auto procsRef = procs;
    decomposition = combigrid::getStandardDecomposition(lref, procsRef);

    // make sure that registered DFGs set the sizes right
    auto dfgLevel = lmin;
    dfgLevel[0] = lmax[0] + 1;
    auto dfgDecomposition =
        combigrid::downsampleDecomposition(decomposition, lref, dfgLevel, boundary);

    auto uniDFG = std::unique_ptr<DistributedFullGrid<std::complex<double>>>(
        new DistributedFullGrid<std::complex<double>>(dim, dfgLevel, comm, boundary, procs, true, dfgDecomposition));

    uniDFG->registerUniformSG(*uniDSG);
    uniDFG->registerUniformSG(*uniDSGfromSubspaces);
    BOOST_CHECK_EQUAL(0, uniDSGfromSubspaces->getRawDataSize());

    for (size_t i = 0; i < uniDSG->getNumSubspaces(); ++i) {
      const auto & level = uniDSG->getLevelVector(i);
      BOOST_ASSERT(lmin.size() > 1);
      auto secondDimLevel = level[1];
      if(secondDimLevel > lmin[1]) {
        BOOST_CHECK(secondDimLevel <= lmax[1]);
        BOOST_CHECK_EQUAL(uniDSG->getDataSize(i), 0);
        BOOST_CHECK_EQUAL(uniDFG->getFGPointsOfSubspace(level).size(), 0);
      } else {
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
    uniDFG->addToUniformSG(*uniDSG, 1.);
    BOOST_CHECK_GT(uniDSG->getRawDataSize(), uniDSGfromSubspaces->getRawDataSize());
    uniDFG->addToUniformSG(*uniDSGfromSubspaces, 0.);
    BOOST_CHECK_EQUAL(uniDSG->getRawDataSize(), uniDSGfromSubspaces->getRawDataSize());
    for (size_t i = 0; i < uniDSGfromSubspaces->getNumSubspaces(); ++i) {
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
      uniDSGfromSubspaces->setZero();
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
      system("rm sparse_paraboloid_minmax_*");
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

    // make sure that right min/max values are written %TODO remove
    uniDSG->writeMinMaxCoefficents(
        "sparse_paraboloid_minmax_large_" + std::to_string(dim) + "D_" + std::to_string(size), 0);

    // check if the sizes set are actually the ones we calculate with CombiMinMaxScheme
    auto subspacesDataSizes = uniDSG->getSubspaceDataSizes();
    size_t numDataPointsHere =
        std::accumulate(subspacesDataSizes.begin(), subspacesDataSizes.end(), 0);
    BOOST_CHECK(numDataPointsHere > 0);

    if (std::all_of(boundary.begin(), boundary.end(), [](bool i) { return i; })) {
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
      BOOST_TEST_CONTEXT(std::to_string(i))
      BOOST_CHECK_EQUAL(uniDSG->getRawData()[i], uniDSGfromSubspaces->getRawData()[i]);
    }

    // and remove straight away
    if (rank == 0) {
      system( "rm test_sg_*" );
    }
  }
}

BOOST_AUTO_TEST_SUITE(distributedsparsegrid, *boost::unit_test::timeout(240))
//very cheap
BOOST_AUTO_TEST_CASE(test_0) {
  LevelVector lmin = {1,1};
  LevelVector lmax = {3,3};
  for (bool bValue : {true}) {
    IndexVector procs = {1,1};
    std::vector<bool> boundary(2, bValue);
    auto multProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
    checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_CASE(test_1) {
  LevelVector lmin = {3,3};
  LevelVector lmax = {7,7};
  for (IndexType procOne : {1,2,3}) {
    for (IndexType procTwo : {1,2}) {
      for (bool bValue : {true}) {
        IndexVector procs = {procOne, procTwo};
        std::vector<bool> boundary(2, bValue);
        auto multProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
        BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
        checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
}
//same for anisotropic
BOOST_AUTO_TEST_CASE(test_2) {
  LevelVector lmin = {2,4};
  LevelVector lmax = {6,8};
  for (IndexType procOne : {1,2,3}) {
    for (IndexType procTwo : {1,2}) {
      for (bool bValue : {true}) {
        IndexVector procs = {procOne, procTwo};
        std::vector<bool> boundary(2, bValue);
        auto multProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
        BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
        checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
}
//and non-regular truncated CT
BOOST_AUTO_TEST_CASE(test_3) {
  LevelVector lmin = {2,4};
  LevelVector lmax = {9,9};
  for (IndexType procOne : {1,2,3}) {
    for (IndexType procTwo : {1,2}) {
      for (bool bValue : {true, false}) {
        IndexVector procs = {procOne, procTwo};
        std::vector<bool> boundary(2, bValue);
        auto multProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
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
  for (IndexType procOne : {1, 2, 3}) {
    for (IndexType procTwo : {1, 2}) {
      for (bool bValue : {true}) {
        IndexVector procs = {procOne, procTwo, 1};
        std::vector<bool> boundary(3, bValue);
        auto multProcs =
            std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
        BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
        checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
}
// same for anisotropic
BOOST_AUTO_TEST_CASE(test_5) {
  LevelVector lmin = {2, 3, 4};
  LevelVector lmax = {6, 7, 8};
  for (IndexType procOne : {1, 2, 3}) {
    for (IndexType procTwo : {1, 2}) {
      for (bool bValue : {true}) {
        IndexVector procs = {procOne, procTwo, 1};
        std::vector<bool> boundary(3, bValue);
        auto multProcs =
            std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
        BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
        checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
}
// non-regular truncated CT
BOOST_AUTO_TEST_CASE(test_6) {
  LevelVector lmin = {2, 3, 4};
  LevelVector lmax = {7, 7, 7};
  for (IndexType procOne : {1, 2, 3}) {
    for (IndexType procTwo : {1, 2}) {
      for (bool bValue : {true, false}) {
        IndexVector procs = {procOne, procTwo, 1};
        std::vector<bool> boundary(3, bValue);
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
  for (IndexType procOne : {1, 3}) {
    for (IndexType procTwo : {1, 2}) {
      for (bool bValue : {true, false}) {
        IndexVector procs = {procOne, procTwo, 1, 1};
        std::vector<bool> boundary(4, bValue);
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
  for (IndexType procOne : {1, 3}) {
    for (IndexType procTwo : {1, 2}) {
      for (bool bValue : {true, false}) {
        IndexVector procs = {1, 1, procOne, 1, procTwo, 1};
        std::vector<bool> boundary(6, bValue);
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
  IndexVector procs = {2};
  std::vector<LevelVector> decomposition = {{0, 1}};
  auto partitionedNumDOFs = getPartitionedNumDOFSGAdaptive(lmin, lmax, lmax, decomposition);
  auto sln = LevelVector({1,32});
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
  IndexVector procs = {2, 2};
  std::vector<LevelVector> decomposition = {{0, 1}, {0, 3}};
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
  IndexVector procs = {2, 2, 2};
  std::vector<LevelVector> decomposition = {{0, 1}, {0, 2}, {0, 3}};
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

  BOOST_AUTO_TEST_SUITE_END()
