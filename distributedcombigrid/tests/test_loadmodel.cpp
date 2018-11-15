#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <thread>
#include <limits>
#include <random>

#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "test_helper.hpp"

using namespace combigrid;

void testDataSave(int size){

    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

    CommunicatorType comm = TestHelper::getComm(size);
    if (comm == MPI_COMM_NULL){ return; }

    combigrid::Stats::initialize();

    size_t ngroup = size-1;
    size_t nprocs = 1;
    theMPISystem()->initWorld(comm, ngroup, nprocs);
    MPI_Datatype durationType_ = createMPIDurationType();


    WORLD_MANAGER_EXCLUSIVE_SECTION{
        std::vector<LevelVector> lvv;
        for (long int i = 0; i < ngroup; ++i){
            lvv.push_back({i});
        }
        auto loadModel = std::unique_ptr<LoadModel>(new LearningLoadModel(lvv));
        for (size_t i = 0; i < ngroup; ++i) {
            durationInformation recvbuf;
            MPI_Status stat;
            MPI_Recv(&recvbuf, 1, durationType_, MPI_ANY_SOURCE, durationTag, theMPISystem()->getGlobalComm(), &stat);
            if(LearningLoadModel* llm = dynamic_cast<LearningLoadModel*>(loadModel.get())){
                llm->addDataPoint(recvbuf, {stat.MPI_SOURCE});
            }
        }
        //test loadmodel
        for (long int i = 0; i < ngroup; ++i){
            // std::cout << "llm eval " << loadModel->eval({i}) << std::endl;
            BOOST_TEST(loadModel->eval({i}) == 1000000*i);
        }
    }else{
        long int nprocs = 64000;

        // // Choose a random duration
        // std::random_device r; 
        // std::default_random_engine e1(r());
        // std::uniform_int_distribution<long unsigned int> uniform_dist(1/*, default value for upper limit is std::numeric_limits<long>::max()*/);
        // long unsigned int d = uniform_dist(e1);
        
        long unsigned int d = 1000000*getCommRank(comm);

        Stats::Event e = Stats::Event();
        e.end = e.start + std::chrono::milliseconds(d);
        durationInformation info = {TestHelper::getRank(comm), Stats::getEventDuration(e), nprocs};
        // MPI_Request request;
        // send durationInfo to manager
        MPI_Send(&info, 1, durationType_, theMPISystem()->getManagerRank(), durationTag, //TODO see if we can send asynchronously
                theMPISystem()->getGlobalComm());
    }
   
    combigrid::Stats::finalize();
}


BOOST_AUTO_TEST_SUITE(loadmodel)

// BOOST_AUTO_TEST_CASE(test_2) {
//     testDataSave(2);
// }

BOOST_AUTO_TEST_CASE(test_9, * boost::unit_test::timeout(60)) {
    testDataSave(9);
    MPI_Barrier(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_SUITE_END()
