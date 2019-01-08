#ifndef THIRDLEVELUTILSHPP_
#define THIRDLEVELUTILSHPP_

#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <stdlib.h>
#include <ctime>
#include <sstream>
#include <mutex>
#include <thread>
#include <chrono>

// does not have to be group 0
// manager process in globalReduceComm should have the same id as group
#define THIRD_LEVEL_GROUP 0
#define THIRD_LEVEL_MANAGER THIRD_LEVEL_GROUP
#define THIRD_LEVEL_PROCESSGROUP_EXCLUSIVE \
  if ( (combigrid::theMPISystem()->getWorldRank() / combigrid::theMPISystem()-> getNumGroups()) == THIRD_LEVEL_GROUP )

#define NUM_SYSTEMS 2

namespace combigrid {

  class ThirdLevelUtils {
    public:
      ThirdLevelUtils(const std::string& remoteHost, int mesgPort);
      ~ThirdLevelUtils();

      void connectToIntermediary();

      /*
       * Reduces global first and with remote after.
       */
      template<typename FG_ELEMENT> void
        reduceGlobalRemote(DistributedSparseGridUniform<FG_ELEMENT>& dsg)
        {
           MPI_Comm mycomm = theMPISystem()->getGlobalReduceComm();
           assert(mycomm != MPI_COMM_NULL);

           /* get sizes of all partial subspaces in communicator
            * we have to do this, because size information of uninitialized subspaces
            * is not available in dsg. at the moment this information is only available
            * in dfg.
            */
         std::vector<int> subspaceSizes(dsg.getNumSubspaces());

           for (size_t i = 0; i < subspaceSizes.size(); ++i) {
             // MPI does not have a real size_t equivalent. int should work in most cases
             // if not we can at least detect this with an assert
             assert(dsg.getDataSize(i) <= INT_MAX);

             subspaceSizes[i] = int(dsg.getDataSize(i));
           }

           MPI_Allreduce( MPI_IN_PLACE, subspaceSizes.data(), int(subspaceSizes.size()),
                          MPI_INT, MPI_MAX, mycomm);

           // *****************************************************************
           // exchange subspaceSizes with other system
           THIRD_LEVEL_PROCESSGROUP_EXCLUSIVE {
             // connect to communicator at intermediary
             if ( this->mesgClient_ == NULL )
               connectToIntermediary();

             // wait for other system
             int remoteSize = getCommSize();
             std::cout << "Waiting for other system" << std::endl;
             while(remoteSize != NUM_SYSTEMS) {
               std::cout << ".";
               std::this_thread::sleep_for(std::chrono::milliseconds(1000));
               remoteSize = getCommSize();
             }
             std::cout << std::endl;

             int otherSysRank = (remoteRank_ + 1) % NUM_SYSTEMS;

             std::vector<int> recvSizes(dsg.getNumSubspaces());

             std::cout << "Process with rank " << theMPISystem()->getLocalRank()
                       << " exchanges subspaceSizes"<< std::endl;

             std::thread isend(&ThirdLevelUtils::send<int>, this,
                               std::ref(subspaceSizes), otherSysRank);
             bool status = recv(recvSizes, dsg.getNumSubspaces(), otherSysRank);
             assert(status && "Receive of subspaceSizes failed");

             isend.join();
             std::cout << "exchange of sizes finished" << std::endl;
             // update local sizes with those from remote system
             for (size_t i = 0; i < subspaceSizes.size(); ++i) {
               subspaceSizes[i] = std::max(subspaceSizes[i], recvSizes[i]);
             }
           }

           // bcast updated sizes to all procs in global comm
           MPI_Bcast( subspaceSizes.data(), int(subspaceSizes.size()), MPI_INT,
                      THIRD_LEVEL_MANAGER, mycomm);
           // *****************************************************************

           // check for implementation errors, the reduced subspace size should not be
           // different from the size of already initialized subspaces
           int bsize = 0;

           for (size_t i = 0; i < subspaceSizes.size(); ++i) {
             bool check = (subspaceSizes[i] == 0 || dsg.getDataSize(i) == 0
                           || subspaceSizes[i] == int(dsg.getDataSize(i)));

             if (!check) {
               int rank;
               MPI_Comm_rank( MPI_COMM_WORLD, &rank);
               std::cout << "l = " << dsg.getLevelVector(i) << " " << "rank = " 
                         << rank << " " << "ssize = " << subspaceSizes[i] 
                         << " " << "dsize = " << dsg.getDataSize(i) << std::endl;
               assert(false);
             }

             bsize += subspaceSizes[i];
           }

           THIRD_LEVEL_PROCESSGROUP_EXCLUSIVE {
             int otherSysRank = (remoteRank_ + 1) % NUM_SYSTEMS;

             std::vector<FG_ELEMENT> recvbuf(bsize, FG_ELEMENT(0));
             std::vector<FG_ELEMENT> sendbuf(bsize, FG_ELEMENT(0));

             setSubspaces(dsg, sendbuf, subspaceSizes);

             // perform globalReduce on same system
             MPI_Datatype dtype = abstraction::getMPIDatatype(
                 abstraction::getabstractionDataType<FG_ELEMENT>());
             MPI_Allreduce( MPI_IN_PLACE, sendbuf.data(), bsize, dtype,
                            MPI_SUM, mycomm);

             // ****************************************************************
             // exchange dsg with remote system
             // for 2 systems only TODO
             //
             // exchange data
             double start, finish;
             start = MPI_Wtime();

             std::cout << "Process with rank " << theMPISystem()->getLocalRank()
               << " starts remote reduce with a buffer of size: " << bsize
               << std::endl;
             std::thread isend(&ThirdLevelUtils::send<FG_ELEMENT>, this,
                 std::ref(sendbuf), otherSysRank);
             bool status = recv(recvbuf, bsize, otherSysRank);
             assert(status && "recv of dsg failed");
             isend.join();

             finish = MPI_Wtime();
             std::cout << "Third level reduce took " << finish-start << "seconds" 
                       << std::endl;
             std::cout << "remote reduce finished" << std::endl;
             // ****************************************************************

             // add received values to sendbuf
             for (size_t i = 0; i < bsize; i++) {
               sendbuf[i] += recvbuf[i];
             }

             // exchange fully reduced data
             int rank;
             MPI_Comm_rank(mycomm, &rank);
             MPI_Bcast(sendbuf.data(), bsize, dtype, rank, mycomm);

             extractSubspaces(dsg, recvbuf, subspaceSizes);

           } else {
             std::vector<FG_ELEMENT> buf(bsize, FG_ELEMENT(0));
             setSubspaces(dsg, buf, subspaceSizes);

             // perform globalReduce on same system
             MPI_Datatype dtype = abstraction::getMPIDatatype(
                 abstraction::getabstractionDataType<FG_ELEMENT>());
             MPI_Allreduce( MPI_IN_PLACE, buf.data(), bsize, dtype, 
                            MPI_SUM, mycomm);

             // receive data from remote
             MPI_Bcast(buf.data(), bsize, dtype, THIRD_LEVEL_MANAGER, mycomm);

             // extract values of received dsg
             extractSubspaces(dsg, buf, subspaceSizes);
           }
         }

      /*
       * Reduces global first and sends values to the intermediary who performs
       * the remote reduce and writes the solution to a file
       */
      template<typename FG_ELEMENT> void
        reduceGlobalToRemoteFile(DistributedSparseGridUniform<FG_ELEMENT>& dsg)
        {
           MPI_Comm mycomm = theMPISystem()->getGlobalReduceComm();
           assert(mycomm != MPI_COMM_NULL);


           /* get sizes of all partial subspaces in communicator
            * we have to do this, because size information of uninitialized subspaces
            * is not available in dsg. at the moment this information is only available
            * in dfg.
            */
           std::vector<int> subspaceSizes(dsg.getNumSubspaces());

           for (size_t i = 0; i < subspaceSizes.size(); ++i) {
             // MPI does not have a real size_t equivalent. int should work in most cases
             // if not we can at least detect this with an assert
             assert(dsg.getDataSize(i) <= INT_MAX);

             subspaceSizes[i] = int(dsg.getDataSize(i));
           }


          std::cout << "Performing allreduce to fetch correct subspace sizes" << std::endl;
           MPI_Allreduce( MPI_IN_PLACE, subspaceSizes.data(), int(subspaceSizes.size()),
                          MPI_INT, MPI_MAX, mycomm);

           // check for implementation errors, the reduced subspace size should not be
           // different from the size of already initialized subspaces
           int bsize = 0;

           for (size_t i = 0; i < subspaceSizes.size(); ++i) {
             bool check = (subspaceSizes[i] == 0 || dsg.getDataSize(i) == 0
                           || subspaceSizes[i] == int(dsg.getDataSize(i)));

             if (!check) {
               int rank;
               MPI_Comm_rank( MPI_COMM_WORLD, &rank);
               std::cout << "l = " << dsg.getLevelVector(i) << " " << "rank = " << rank
                         << " " << "ssize = " << subspaceSizes[i] << " " << "dsize = "
                         << dsg.getDataSize(i) << std::endl;
               assert(false);
             }

             bsize += subspaceSizes[i];
           }

           std::cout << "Done. starting globalReduce on same system" << std::endl;

           THIRD_LEVEL_PROCESSGROUP_EXCLUSIVE {
             // connect to communicator at intermediary
             if ( this->mesgClient_ == NULL )
               connectToIntermediary();

             std::vector<FG_ELEMENT> sendbuf(bsize, FG_ELEMENT(0));

             setSubspaces(dsg, sendbuf, subspaceSizes);

             // perform globalReduce on same system
             MPI_Datatype dtype = abstraction::getMPIDatatype(
                 abstraction::getabstractionDataType<FG_ELEMENT>());
             MPI_Allreduce( MPI_IN_PLACE, sendbuf.data(), bsize, dtype, MPI_SUM, mycomm);

             // *******************************************************************
             // send dsg to remote system TODO
             int rank;
             MPI_Comm_rank( MPI_COMM_WORLD, &rank);
             std::cout << "Process with rank " <<  rank << " starts reduceToFileUniform";
             reduceToFileUniform(sendbuf);
             // *******************************************************************
           } else {
             std::vector<FG_ELEMENT> buf(bsize, FG_ELEMENT(0));
             setSubspaces(dsg, buf, subspaceSizes);

             // perform globalReduce on same system
             MPI_Datatype dtype = abstraction::getMPIDatatype(
                 abstraction::getabstractionDataType<FG_ELEMENT>());
             MPI_Allreduce( MPI_IN_PLACE, buf.data(), bsize, dtype, MPI_SUM, mycomm);
            }
         }


    private:
      std::string remoteHost_;
      int remotePort_, remoteRank_;
      ClientSocket* mesgClient_;

      std::mutex mesgLock_;

      /*
       * fetches the current size of the third level communicator.
       */
      int getCommSize();

      /*
       * returns the rank in third level communicator
       */
      int getRank();

      /*
       * requests a new data connection at intermediary
       */
      void requestDataConn(const ClientSocket& dataClient) const;

      /*
       * Perfoms an in place allreduce with all processes in the communicator at
       * intermediary.
       */
      template<typename FG_ELEMENT>
        void allreduce(std::vector<FG_ELEMENT>& sendrecvbuf,
            size_t numparticipants) const {
          //TODO
        }

      /*
       * sends data to process on remote system with rank rank dest.
       */
      template<typename FG_ELEMENT> bool
        send(const std::vector<FG_ELEMENT>& buff, int dest)
        {
          double start, finish;
          start = MPI_Wtime();

          size_t dataSize = buff.size() * sizeof(FG_ELEMENT) + 1; // +1 endiannes header
          std::string mesg = "send#"+ std::to_string(dest) + "#"
            + std::to_string(dataSize);
          exchangeMesg(mesg);

          // establish new connection for data transfer
          ClientSocket dataClient(remoteHost_, remotePort_);
          assert(dataClient.init() && "Connecting to data server failed");

          // tell intermediary to which participant i am related
          requestDataConn(dataClient);

          // start data transfer
          bool status = NetworkUtils::sendBinary(buff, dataClient);

          finish = MPI_Wtime();
          std::cout << "Sending of " << dataSize << " Bytes took" << finish-start
                    << "Seconds" << std::endl;
          return status;
        }

      /*
       * receives data from process on remote system with rank source.
       */
      template<typename FG_ELEMENT> bool
        recv(std::vector<FG_ELEMENT>& buff, size_t len, int source)
        {
          double start, finish;
          start = MPI_Wtime();

          size_t dataSize = len * sizeof(FG_ELEMENT) + 1; // +1 endiannes header
          std::string mesg = "recv#" + std::to_string(source) + "#"
            + std::to_string(dataSize);
          exchangeMesg(mesg);

          // establish new connection for data transfer
          ClientSocket dataClient(remoteHost_, remotePort_);
          assert(dataClient.init() && "Connecting to data server failed");

          // tell intermediary to which participant i am related
          requestDataConn(dataClient);

          // start data transfer
          bool status = NetworkUtils::recvBinary(buff, len, dataClient);

          finish = MPI_Wtime();
          std::cout << "Receiving of " << dataSize << " Bytes took" << finish-start
                    << "Seconds" << std::endl;
          return status;
        }

      /*
       * Initiates the third level reduce on the intermediary. The solution is
       * written to a remote file.
       */
      template<typename FG_ELEMENT> bool
        reduceToFileUniform(const std::vector<FG_ELEMENT>& buff) {
          size_t commSize = NUM_SYSTEMS;
          std::string mesg = "reduceToFileUniform#" + std::to_string(commSize) + "#"
            + std::to_string(buff.size()) + "#" +  std::to_string(sizeof(FG_ELEMENT)) +
            "unknown";
          exchangeMesg(mesg);

          // establish new connection for data transfer
          ClientSocket dataClient(remoteHost_, remotePort_);
          assert(dataClient.init() && "Connecting to data server failed");

          // tell intermediary to which participant i am related
          requestDataConn(dataClient);

          // start data transfer
          return NetworkUtils::sendBinary(buff, dataClient);
        }

      /*
       * Blocks until commSize participants have joined to the communicator on 
       * intermediary.
       */
      void barrier(size_t commSize);

      /*
       * Talks to intermediary by sending an operation message and receiving the
       * appropriate response.
       */
      void exchangeMesg(std::string& mesg);

      template<typename FG_ELEMENT> void
        extractSubspaces(DistributedSparseGridUniform<FG_ELEMENT>& dsg,
            std::vector<FG_ELEMENT>& buf, std::vector<int>& subspaceSizes) const
        {
          typename std::vector<FG_ELEMENT>::iterator buf_it = buf.begin();

          for (size_t i = 0; i < dsg.getNumSubspaces(); ++i) {
            std::vector<FG_ELEMENT>& subspaceData = dsg.getDataVector(i);

            // this is very unlikely but can happen if dsg is different than
            // lmax and lmin of combination scheme
            if(subspaceData.size() == 0 && subspaceSizes[i] == 0)
              continue;

            // this happens for subspaces that are only available in component
            // grids on other process groups
            if( subspaceData.size() == 0 && subspaceSizes[i] > 0 ){
              subspaceData.resize( subspaceSizes[i] );
            }

            // wenn subspaceData.size() > 0 und subspaceSizes > 0
            for (size_t j = 0; j < subspaceData.size(); ++j) {
              subspaceData[j] = *buf_it;
              ++buf_it;
            }
          }
        }


      template<typename FG_ELEMENT> void
        setSubspaces(DistributedSparseGridUniform<FG_ELEMENT>& dsg,
            std::vector<FG_ELEMENT>& buf, std::vector<int>& subspaceSizes) const
        {
          typename std::vector<FG_ELEMENT>::iterator buf_it = buf.begin();

          for (size_t i = 0; i < dsg.getNumSubspaces(); ++i) {
            std::vector<FG_ELEMENT>& subspaceData = dsg.getDataVector(i);

            // if subspace does not exist on this process this part of the send
            // buffer is left empty
            if (subspaceData.size() == 0) {
              buf_it += subspaceSizes[i];
              continue;
            }

            for (size_t j = 0; j < subspaceData.size(); ++j) {
              *buf_it = subspaceData[j];
              ++buf_it;
            }
          }
        }

  }; // end class
} // end namespace

#endif
