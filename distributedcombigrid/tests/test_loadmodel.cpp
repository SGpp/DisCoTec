#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <thread>
#include <limits>
#include <random>

#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "test_helper.hpp"

using namespace combigrid;

void testDataSave(int size){

    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

    MPI_Barrier(MPI_COMM_WORLD);

    CommunicatorType comm = TestHelper::getComm(size);

    if (comm == MPI_COMM_NULL){ return; }

    combigrid::Stats::initialize();

    BOOST_REQUIRE(true); //if things go wrong weirdly, see where things go wrong

    size_t nprocs = 64000;

    std::random_device r; 
    // Choose a random duration
    std::default_random_engine e1(r());
    std::uniform_int_distribution<long> uniform_dist(1/*, default value for upper limit is std::numeric_limits<long>::max()*/);
    int d = uniform_dist(e1);
    
    Stats::Event e = Stats::Event();
    //TODO add test for new loadmodel

   
    combigrid::Stats::finalize();
}


BOOST_AUTO_TEST_SUITE(loadmodel)

BOOST_AUTO_TEST_CASE(test_1) {
    testDataSave(1);
}

BOOST_AUTO_TEST_CASE(test_9, * boost::unit_test::timeout(20)) {
    testDataSave(9);
}

BOOST_AUTO_TEST_SUITE_END()
