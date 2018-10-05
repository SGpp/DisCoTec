#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "test_helper.hpp"

BOOST_AUTO_TEST_SUITE(stats)

BOOST_AUTO_TEST_CASE(test_1, * boost::unit_test::timeout(20)) {

    int size = 2;
    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

    combigrid::CommunicatorType comm = TestHelper::getComm(size);

    if(comm == MPI_COMM_NULL){
        std::cout << "MPI_COMM_NULL" << std::endl;
        return;
    }
    else{
        std::cout << "YUP" << std::endl;
    }

    combigrid::Stats::initialize();

    size_t ngroup = 1;
    size_t nprocs = 1;
    //combigrid::theMPISystem()->initWorld(comm, ngroup, nprocs);
    combigrid::Stats::finalize();

}

BOOST_AUTO_TEST_SUITE_END()
