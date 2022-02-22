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
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
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

    if (TestHelper::getRank(comm) == 0) {
      std::cout << "test distributedsparsegrid " << lmin << lmax << procs << std::endl;
    }

    assert(lmin.size() == lmax.size());
    const DimType dim = lmin.size();

    // create with "own" constructor
    auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<std::complex<double>>>(
        new DistributedSparseGridUniform<std::complex<double>>(dim, lmax, lmin, boundary, comm));

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

    // make sure that registered DFGs set the sizes right
    auto dfgLevel = lmin;
    dfgLevel[0] = lmax[0] + 1;
    auto uniDFG = std::unique_ptr<DistributedFullGrid<std::complex<double>>>(
      new DistributedFullGrid<std::complex<double>>(dim, dfgLevel, comm, boundary, procs));

    uniDFG->registerUniformSG(*uniDSG);
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
    // attn: usually, the dfgs are hierarchized before adding to the dsg for combination in hierarchical space (without interpolation on coarser component dfgs)
    // here we use the nodal values for testing purposes
    uniDFG->addToUniformSG(*uniDSG, 1.);
    //TODO test

    // make sure that right min/max values are written
    uniDSG->writeMinMaxCoefficents("sparse_paraboloid_minmax_" + std::to_string(size),0);

    // std::cout << *uniDSG << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE(distributedsparsegrid, *boost::unit_test::timeout(60))

BOOST_AUTO_TEST_CASE(test_1) {
  LevelVector lmin = {3,3};
  LevelVector lmax = {7,7};
  IndexVector procs = {1, 1};
  for (IndexType procOne : {1,2,3}) {
    for (IndexType procTwo : {1,2}) {
      for (bool bValue : {true}) { //TODO false
        IndexVector procs = {procOne, procTwo};
        std::vector<bool> boundary(2, bValue);
        auto multProcs = std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<IndexType>());
        BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(multProcs));
        checkDistributedSparsegrid(lmin, lmax, procs, boundary, multProcs);
      }
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()
