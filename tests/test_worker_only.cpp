#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>
#include <filesystem>

#include "TaskCount.hpp"
#include "combischeme/CombiMinMaxScheme.hpp"
#include "io/H5InputOutput.hpp"
#include "loadmodel/LearningLoadModel.hpp"
#include "loadmodel/LinearLoadModel.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupWorker.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
#include "stdlib.h"
#include "task/Task.hpp"
#include "test_helper.hpp"
#include "utils/Config.hpp"
#include "utils/MonteCarlo.hpp"
#include "utils/Types.hpp"

using namespace combigrid;

// omitted here, because already defined in test_thirdLevel.cpp
// BOOST_CLASS_EXPORT(TaskCount)

// declare here only, because already defined in test_integration.cpp
bool checkReducedFullGridIntegration(ProcessGroupWorker& worker, int nrun);

void checkWorkerOnly(size_t ngroup = 1, size_t nprocs = 1, BoundaryType boundaryV = 2,
                     bool pretendThirdLevel = true) {
  size_t size = ngroup * nprocs;
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) {
    BOOST_TEST_CHECKPOINT("drop out of test comm");
    return;
  }
  BOOST_TEST_CHECKPOINT("initialize stats");
  combigrid::Stats::initialize();

  theMPISystem()->initWorldReusable(comm, ngroup, nprocs, false);
  // theMPISystem()->init(ngroup, nprocs, false);

  DimType dim = 2;
  LevelVector lmin(dim, 2);
  LevelVector lmax(dim, 5);

  size_t ncombi = 4;

  BOOST_CHECK_EQUAL(theMPISystem()->getWorldSize(), size);
  BOOST_CHECK_EQUAL(getCommSize(theMPISystem()->getWorldComm()), size);
  BOOST_CHECK_EQUAL(getCommSize(theMPISystem()->getGlobalReduceComm()), ngroup);
  BOOST_CHECK_EQUAL(getCommSize(theMPISystem()->getLocalComm()), nprocs);
  if (theMPISystem()->getOutputGroupComm() != MPI_COMM_NULL) {
    BOOST_CHECK_EQUAL(getCommSize(theMPISystem()->getOutputGroupComm()), nprocs);
  }
  if (nprocs == 1) {
    BOOST_CHECK(theMPISystem()->isMaster());
  }

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    // make sure there is no manager
    BOOST_CHECK(false);
  }
  auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

  std::vector<BoundaryType> boundary(dim, boundaryV);

  CombiMinMaxScheme combischeme(dim, lmin, lmax);
  combischeme.createAdaptiveCombischeme();

  std::vector<LevelVector> levels = combischeme.getCombiSpaces();
  std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

  ProcessGroupWorker worker;

  std::vector<size_t> taskIDs(levels.size());
  std::iota(taskIDs.begin(), taskIDs.end(), 0);

  // create combiparameters
  BOOST_TEST_CHECKPOINT("create combi parameters");
  CombiParameters params(dim, lmin, lmax, boundary, ncombi, 1,
                         CombinationVariant::outgroupSparseGridReduce,
                         {static_cast<int>(nprocs), 1}, LevelVector(0), LevelVector(0), false);
  if (nprocs == 5 && boundaryV == 2) {
    params.setDecomposition({{0, 6, 13, 20, 27}, {0}});
  } else if (nprocs == 4 && boundaryV == 2) {
    // should be the same as default decomposition with forwardDecomposition
    params.setDecomposition({{0, 9, 17, 25}, {0}});
  } else if (nprocs == 3) {
    params.setDecomposition({{0, 15, 20}, {0}});
  }
  worker.setCombiParameters(std::move(params));

  std::vector<size_t> myTaskIDs;
  for (size_t i = 0; i < levels.size(); ++i) {
    // assign round-robin to process groups based on task id
    if (i % ngroup == theMPISystem()->getProcessGroupNumber()) {
      myTaskIDs.push_back(i);
    }
  }
  std::vector<LevelVector> myLevels;
  std::vector<real> myCoeffs;
  for (size_t i = 0; i < levels.size(); ++i) {
    if (std::find(myTaskIDs.begin(), myTaskIDs.end(), i) != myTaskIDs.end()) {
      myLevels.push_back(levels[i]);
      myCoeffs.push_back(coeffs[i]);
    }
  }
  BOOST_TEST(myTaskIDs.size() > 0);

  // create Tasks
  worker.initializeAllTasks<TaskCount>(myLevels, myCoeffs, myTaskIDs, loadmodel.get());

  BOOST_TEST_CHECKPOINT("initialize sparse grid");
  worker.initCombinedDSGVector();
  if (pretendThirdLevel) {
    std::string subspaceSizeFile = "worker_only_subspace_sizes";
    std::string subspaceSizeFileToken = "worker_only_subspace_sizes_token.txt";
    BOOST_TEST_CHECKPOINT("reduce sparse grid sizes");
    worker.reduceExtraSubspaceSizesFileBased(subspaceSizeFile, subspaceSizeFileToken,
                                             subspaceSizeFile, subspaceSizeFileToken);
    // remove subspace size files to avoid interference between multiple calls to this test function
    BOOST_TEST_CHECKPOINT("reduced sparse grid sizes");
    MPI_Barrier(comm);
    OUTPUT_GROUP_EXCLUSIVE_SECTION {
      MASTER_EXCLUSIVE_SECTION {
        remove(subspaceSizeFile.c_str());
        remove(subspaceSizeFileToken.c_str());
      }
    }
  }

  BOOST_CHECK_EQUAL(worker.getCurrentNumberOfCombinations(), 0);
  BOOST_TEST_CHECKPOINT("run first");
  worker.runAllTasks();
  MASTER_EXCLUSIVE_SECTION {
    BOOST_TEST_MESSAGE("worker run first solver step: " << Stats::getDuration("run")
                                                        << " milliseconds");
  }

  for (size_t it = 0; it < ncombi - 1; ++it) {
    BOOST_TEST_CHECKPOINT("combine");
    auto start = std::chrono::high_resolution_clock::now();
    worker.combineUniform();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    MASTER_EXCLUSIVE_SECTION {
      BOOST_TEST_MESSAGE("worker combine: " << duration.count() << " milliseconds");
    }
    BOOST_CHECK_EQUAL(worker.getCurrentNumberOfCombinations(), it + 1);

    if (boundaryV == 2) {
      // check if the values are as expected
      BOOST_CHECK(checkReducedFullGridIntegration(worker, worker.getCurrentNumberOfCombinations()));
    }

    // first group writes partial stats
    if (theMPISystem()->getProcessGroupNumber() == 0) {
      Stats::writePartial("worker_partial_timers_group_0.json", theMPISystem()->getLocalComm());
    }

    BOOST_TEST_CHECKPOINT("run next");
    worker.runAllTasks();
    MASTER_EXCLUSIVE_SECTION {
      BOOST_TEST_MESSAGE("worker run: " << Stats::getDuration("run") << " milliseconds");
    }
  }
  BOOST_TEST_CHECKPOINT("worker combine last time");
  if (pretendThirdLevel) {
    std::string writeSparseGridFile = "worker_combine_step_dsg";
    std::string writeSparseGridFileToken = writeSparseGridFile + "_token.txt";
    worker.combineLocalAndGlobal(theMPISystem()->getOutputRankInGlobalReduceComm());
    OUTPUT_GROUP_EXCLUSIVE_SECTION {
      BOOST_TEST_CHECKPOINT("worker write dsg");
      worker.combineThirdLevelFileBasedWrite(writeSparseGridFile, writeSparseGridFileToken);
      BOOST_TEST_CHECKPOINT("worker wrote dsg");
      worker.combineThirdLevelFileBasedReadReduce(writeSparseGridFile, writeSparseGridFileToken,
                                                  true);
      BOOST_TEST_CHECKPOINT("worker read dsg");
    }
    else {
      BOOST_TEST_CHECKPOINT("worker waiting for broadcast dsg");
      worker.waitForThirdLevelCombiResult(true);
      BOOST_TEST_CHECKPOINT("worker got dsg");
    }
  } else {
    worker.combineUniform();
  }

  Stats::startEvent("worker get norms");
  // get all kinds of norms
  worker.getLpNorms(0);
  worker.getLpNorms(1);
  worker.getLpNorms(2);
  Stats::stopEvent("worker get norms");

  BOOST_TEST_CHECKPOINT("write solution");
  std::string filename("worker_" + std::to_string(ncombi) + ".raw");
  BOOST_TEST_CHECKPOINT("write solution " + filename);
  Stats::startEvent("worker write solution");
  FIRST_GROUP_EXCLUSIVE_SECTION { worker.parallelEvalUniform(filename, lmax); }
  BOOST_TEST_CHECKPOINT("write min/max coefficients");
  worker.writeSparseGridMinMaxCoefficients("worker_" + std::to_string(boundaryV) +
                                           "_sparse_minmax");
  Stats::stopEvent("worker write solution");
  MASTER_EXCLUSIVE_SECTION {
    BOOST_TEST_MESSAGE("worker write solution: " << Stats::getDuration("worker write solution")
                                                 << " milliseconds");
  }
  filename = "worker_" + std::to_string(boundaryV) + "_dsgs";
  Stats::startEvent("worker write DSG");
  int numWritten = 0;
  OUTPUT_GROUP_EXCLUSIVE_SECTION {
    BOOST_TEST_CHECKPOINT("write DSGS " + filename);
    numWritten = worker.writeDSGsToDisk(filename);
    BOOST_CHECK(numWritten > 0);
  }
  // overwrite for all--only with sparse grid reduce
  // (otherwise, would have wrong sizes for other groups)
  if (params.getCombinationVariant() == CombinationVariant::sparseGridReduce) {
    BOOST_TEST_CHECKPOINT("read DSGS " + filename);
    int numRead = worker.readDSGsFromDisk(filename, true);
    OUTPUT_GROUP_EXCLUSIVE_SECTION { BOOST_CHECK_EQUAL(numRead, numWritten); }
  }
  Stats::stopEvent("worker write DSG");
  MASTER_EXCLUSIVE_SECTION {
    BOOST_TEST_MESSAGE("worker write/read DSG: " << Stats::getDuration("worker write DSG")
                                                 << " milliseconds");
  }
#ifdef DISCOTEC_USE_HIGHFIVE
  // test Monte-Carlo interpolation
  // only if boundary values are used
  if (boundaryV > 0) {
    BOOST_TEST_CHECKPOINT("MC interpolation coordinates");
    // read interpolation coordinates
    std::string interpolationCoordsFile = "worker_coords.h5";
    std::vector<std::vector<double>> interpolationCoords;
    // if the file does not exist, one rank creates it
    if (!std::filesystem::exists(interpolationCoordsFile)) {
      if (theMPISystem()->getWorldRank() == 0) {
        interpolationCoords = montecarlo::getRandomCoordinates(100, dim);
        h5io::writeValuesToH5File(interpolationCoords, interpolationCoordsFile, "worker_group",
                                  "only");
      }
    }
    interpolationCoords = broadcastParameters::getCoordinatesFromRankZero(
        interpolationCoordsFile, theMPISystem()->getWorldComm());
    BOOST_CHECK_EQUAL(interpolationCoords.size(), 100);

    BOOST_TEST_CHECKPOINT("MC interpolation");
    Stats::startEvent("worker interpolate");
    auto values = worker.interpolateValues(interpolationCoords);
    Stats::stopEvent("worker interpolate");
    MASTER_EXCLUSIVE_SECTION {
      BOOST_TEST_MESSAGE("worker interpolate: " << Stats::getDuration("worker interpolate")
                                                << " milliseconds");
    }

    if (boundaryV > 1) {
      TestFnCount<CombiDataType> initialFunction;
      for (size_t i = 0; i < interpolationCoords.size(); ++i) {
        if (std::abs(initialFunction(interpolationCoords[i], ncombi) - values[i]) >
            TestHelper::tolerance) {
          std::cout << "err " << interpolationCoords.size() << interpolationCoords[i] << " " << i
                    << std::endl;
        }
        auto ref = initialFunction(interpolationCoords[i], ncombi);
        BOOST_CHECK_CLOSE(std::abs(ref), std::abs(values[i]), TestHelper::tolerance);
        BOOST_CHECK_CLOSE(std::real(ref), std::real(values[i]), TestHelper::tolerance);
      }
    }
    // output files are not needed, remove previous ones
    // (if this doesn't happen, there may be hdf5 errors due to duplicates)
    OUTPUT_GROUP_EXCLUSIVE_SECTION {
      MASTER_EXCLUSIVE_SECTION {
        auto status = system("rm worker_interpolated*.h5");
        BOOST_WARN_GE(status, 0);
      }
    }
    sleep(1);
    worker.writeInterpolatedValuesSingleFile(interpolationCoords, "worker_interpolated");
    BOOST_TEST_CHECKPOINT("wrote interpolated values to single file");
    worker.writeInterpolatedValuesPerGrid(interpolationCoords, "worker_interpolated");
    BOOST_TEST_CHECKPOINT("wrote interpolated values per grid");

    MPI_Barrier(comm);
    sleep(1);  // wait for filesystem to catch up
    decltype(values) valuesAllGridsRead;
    h5io::readH5Values(valuesAllGridsRead,
                       "worker_interpolated_values_" + std::to_string(ncombi) + ".h5");
    BOOST_CHECK_EQUAL_COLLECTIONS(values.begin(), values.end(), valuesAllGridsRead.begin(),
                                  valuesAllGridsRead.end());
  }

#endif  // def DISCOTEC_USE_HIGHFIVE
  BOOST_TEST_CHECKPOINT("worker exit");
  BOOST_CHECK_EQUAL(worker.getCurrentNumberOfCombinations(), ncombi);
  worker.exit();

  // if output files are not needed, remove them right away
  remove(("worker_" + std::to_string(ncombi) + "_0.raw").c_str());
  remove(("worker_" + std::to_string(ncombi) + "_0.raw_header").c_str());

  BOOST_CHECK(!TestHelper::testStrayMessages(theMPISystem()->getLocalComm()));
  if (theMPISystem()->getOutputGroupComm() != MPI_COMM_NULL) {
    BOOST_CHECK(!TestHelper::testStrayMessages(theMPISystem()->getOutputGroupComm()));
  }
  BOOST_CHECK(!TestHelper::testStrayMessages(theMPISystem()->getGlobalReduceComm()));
  MASTER_EXCLUSIVE_SECTION {
    BOOST_CHECK(!TestHelper::testStrayMessages(theMPISystem()->getGlobalComm()));
  }

  combigrid::Stats::finalize();
  Stats::write("worker_" + std::to_string(ngroup) + "_" + std::to_string(nprocs) + ".json", comm);

  MPI_Barrier(comm);
  BOOST_CHECK(!TestHelper::testStrayMessages(comm));
}

#ifndef ISGENE  // worker tests won't work with ISGENE because of worker magic

#ifndef NDEBUG  // in case of a build with asserts, have longer timeout
BOOST_FIXTURE_TEST_SUITE(worker, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(2000))
#else
BOOST_FIXTURE_TEST_SUITE(worker, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(1060))
#endif  // NDEBUG
BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::tolerance(TestHelper::higherTolerance)) {
  auto start = std::chrono::high_resolution_clock::now();
  auto rank = TestHelper::getRank(MPI_COMM_WORLD);
  for (bool pretendThirdLevel : {false, true}) {
    for (BoundaryType boundary : std::vector<BoundaryType>({0, 1, 2})) {
      for (size_t ngroup : {1, 2, 3, 4}) {
        for (size_t nprocs : {1, 2}) {
          if (rank == 0)
            std::cout << "worker " << static_cast<int>(boundary) << " " << ngroup << " " << nprocs
                      << std::endl;
          BOOST_CHECK_NO_THROW(checkWorkerOnly(ngroup, nprocs, boundary, pretendThirdLevel));
          MPI_Barrier(MPI_COMM_WORLD);
        }
      }
      for (size_t ngroup : {1, 2}) {
        for (size_t nprocs : {3}) {
          if (rank == 0)
            std::cout << "worker " << static_cast<int>(boundary) << " " << ngroup << " " << nprocs
                      << std::endl;
          BOOST_CHECK_NO_THROW(checkWorkerOnly(ngroup, nprocs, boundary, pretendThirdLevel));
          MPI_Barrier(MPI_COMM_WORLD);
        }
      }
      for (size_t ngroup : {1, 2}) {
        if (boundary > 0) {
          for (size_t nprocs : {4}) {
            if (rank == 0)
              std::cout << "worker " << static_cast<int>(boundary) << " " << ngroup << " " << nprocs
                        << std::endl;
            BOOST_CHECK_NO_THROW(checkWorkerOnly(ngroup, nprocs, boundary, pretendThirdLevel));
            MPI_Barrier(MPI_COMM_WORLD);
          }
        }
      }
      for (size_t ngroup : {1}) {
        for (size_t nprocs : {5}) {
          if (boundary == 2) {
            if (rank == 0)
              std::cout << "worker " << static_cast<int>(boundary) << " " << ngroup << " " << nprocs
                        << std::endl;
            BOOST_CHECK_NO_THROW(checkWorkerOnly(ngroup, nprocs, boundary, pretendThirdLevel));
            MPI_Barrier(MPI_COMM_WORLD);
          }
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    BOOST_CHECK(!TestHelper::testStrayMessages());
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  BOOST_TEST_MESSAGE("time to run all 'worker' tests: " << duration.count() << " milliseconds");
}

BOOST_AUTO_TEST_SUITE_END()
#endif
