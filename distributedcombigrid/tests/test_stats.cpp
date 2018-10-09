#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <thread>

#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "test_helper.hpp"


void checkStats(int size){
    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));
    BOOST_CHECK( size % 2 == 1);

    combigrid::CommunicatorType comm = TestHelper::getComm(size);

    if(comm == MPI_COMM_NULL){
        return;
    }

    int globalID;
    MPI_Comm_rank(MPI_COMM_WORLD, &globalID);
    
    combigrid::Stats::initialize();
    combigrid::Stats::setAttribute("group", std::to_string(globalID));

    combigrid::Stats::startEvent("wait 5 seconds");
    std::this_thread::sleep_for(std::chrono::microseconds{5000});
    combigrid::Stats::stopEvent("wait 5 seconds");

    combigrid::Stats::finalize();
    // MPI_Barrier(comm);
    combigrid::Stats::write("./distributedcombigrid/tests/Stats_output", comm);
    // MPI_Barrier(comm);
}

BOOST_AUTO_TEST_SUITE(stats)

BOOST_AUTO_TEST_CASE(test_3, * boost::unit_test::timeout(20)) {
    checkStats(1);
}

BOOST_AUTO_TEST_CASE(test_1, * boost::unit_test::timeout(20)) {
    checkStats(9);
}

BOOST_AUTO_TEST_SUITE_END()
