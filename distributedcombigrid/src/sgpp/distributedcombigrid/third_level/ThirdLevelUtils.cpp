#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"

namespace combigrid {

  ThirdLevelUtils::ThirdLevelUtils(const std::string& remoteHost,
      int mesgPort, int dataPort) : remoteHost_(remoteHost),
  remoteMesgPort_(mesgPort) {}

  ThirdLevelUtils::~ThirdLevelUtils() {}

  /*
   * Connects to the intermediary and saves the third level rank.
   */
  void ThirdLevelUtils::init() {
    const RankType& globalReduceRank = theMPISystem()->getGlobalReduceRank();
    mesgClient_ = new ClientSocket(remoteHost_, remoteMesgPort_);
    mesgClient_->init();
    assert(mesgClient_->isInitialized() &&
        "Third level could not initialize message client");

    mesgClient_->sendallPrefixed(std::to_string(globalReduceRank));

    std::string mesg = "getRank";
    exchangeMesg(mesg);
    assert(NetworkUtils::isInteger(mesg));
    rank_= std::stoi(mesg);
  }

  /*
   * TODO
   * Perfoms an in place allreduce with all processes in the communicator at
   * intermediary.
   */
  template<typename FG_ELEMENT> void
    ThirdLevelUtils::allreduce(std::vector<FG_ELEMENT>& sendrecvbuf,
        unsigned int numparts) const
    {
    }

  int ThirdLevelUtils::getCommSize() {
    std::string mesg = "getCommSize";
    exchangeMesg(mesg);
    assert(NetworkUtils::isInteger(mesg) && "Comm Size must be positiv integer");
    int commSize = std::stoi(mesg);
    return commSize;
  }

  template<typename FG_ELEMENT> void
    ThirdLevelUtils::send(const std::vector<FG_ELEMENT>& buff, int dest)
    {
      assert(dest < getCommSize() && "Destination not in remote communicator");
      std::string mesg = "send#"+ std::to_string(dest) + "#" + std::to_string(buff.size());
      exchangeMesg(mesg);
      assert(NetworkUtils::isInteger(mesg) && "Received data port is not valid");
      int dataPort = std::stoi(mesg);
      ClientSocket dataClient(remoteHost_, dataPort);
      assert(dataClient.init() && "Connecting to data server failed");
      NetworkUtils::sendBinary(buff, dataClient);
    }

  template<typename FG_ELEMENT> void
    ThirdLevelUtils::recv(std::vector<FG_ELEMENT>& buff, size_t len, int source)
    {
      assert(source < getCommSize() && "Source not in remote communicator");
      std::string mesg = "recv#" + std::to_string(source) + "#" + std::to_string(buff.size());
      exchangeMesg(mesg);
      assert(NetworkUtils::isInteger(mesg) && "Received data port is not valid");
      int dataPort = std::stoi(mesg);
      ClientSocket dataClient(remoteHost_, dataPort);
      assert(dataClient.init() && "Connecting to data server failed");
      // TODO
      NetworkUtils::recvBinary(buff, dataClient);
    }

  template<typename FG_ELEMENT> void
    ThirdLevelUtils::reduceToFileUniform(const std::vector<FG_ELEMENT>& buff) {
      size_t commSize = (theMPISystem()->getNumProcs()-1) / theMPISystem()->getNumGroups();
      std::string mesg = "reduceToFileUniform#" + std::to_string(commSize) + "#"
        + std::to_string(buff.size()) + "#" +  std::to_string(sizeof(FG_ELEMENT)) +
        "unknown";
      exchangeMesg(mesg);
      assert(NetworkUtils::isInteger(mesg) && "Received data port is not valid");
      int dataPort = std::stoi(mesg);
      ClientSocket dataClient(remoteHost_, dataPort);
      assert(dataClient.init() && "Connecting to data server failed");
      NetworkUtils::sendBinary(buff, dataClient);
    }

  template<> void
    ThirdLevelUtils::reduceToFileUniform(const std::vector<real>& buff) {
      size_t commSize = (theMPISystem()->getNumProcs()-1) / theMPISystem()->getNumGroups();
      std::string mesg = "reduceToFileUniform#" + std::to_string(commSize) + "#"
        + std::to_string(buff.size()) + "#" +  std::to_string(sizeof(real)) +
        "#real";
      exchangeMesg(mesg);
      assert(NetworkUtils::isInteger(mesg) && "Received data port is not valid");
      int dataPort = std::stoi(mesg);
      ClientSocket dataClient(remoteHost_, dataPort);
      assert(dataClient.init() && "Connecting to data server failed");
      NetworkUtils::sendBinary(buff, dataClient);
    }

  void ThirdLevelUtils::exchangeMesg(std::string& mesg) {
    mesgLock_.lock();
      assert(mesgClient_->sendallPrefixed(mesg) &&
          "Sending to intermediary failed");
      assert(mesgClient_->recvallPrefixed(mesg) &&
          "Receiving from intermediary failed");
    mesgLock_.unlock();
  }

  void ThirdLevelUtils::barrier(size_t commSize) {
    std::string mesg = "barrier#" + std::to_string(commSize);
    exchangeMesg(mesg);
    assert(NetworkUtils::isInteger(mesg) && "Received data port is not valid");
    int dataPort = std::stoi(mesg);
    ClientSocket dataClient(remoteHost_, dataPort);
    assert(dataClient.init() && "Connecting to data server failed");
    char* signal = nullptr;
    // blocks until receives signal
    dataClient.recvall(signal, 1);
    delete[] signal;
  }

  int ThirdLevelUtils::getRank() {
    return rank_;
  }

  template<typename FG_ELEMENT> void
    ThirdLevelUtils::setSubspaces(DistributedSparseGridUniform<FG_ELEMENT>& dsg,
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

  template<typename FG_ELEMENT> void
    ThirdLevelUtils::extractSubspaces(
        DistributedSparseGridUniform<FG_ELEMENT>& dsg,
        std::vector<FG_ELEMENT>& buf, std::vector<int>& subspaceSizes) const
    {
      typename std::vector<FG_ELEMENT>::iterator buf_it = buf.begin();

      for (size_t i = 0; i < dsg.getNumSubspaces(); ++i) {
        std::vector<FG_ELEMENT>& subspaceData = dsg.getDataVector(i);

        // this is very unlikely but can happen if dsg is different than
        // lmax and lmin of combination scheme
        if(subspaceData.size() == 0 && subspaceSizes[i] == 0)
          continue;

        // this happens for subspaces that are only available in component grids
        // on other process groups
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
   ThirdLevelUtils::reduceGlobalRemote(
       DistributedSparseGridUniform<FG_ELEMENT>& dsg) {

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

      THIRD_LEVEL_PROCESSGROUP_EXCLUSIVE {

        // connect to communicator at intermediary
        if ( this->mesgClient_ == NULL )
          init();

        std::vector<FG_ELEMENT> recvbuf(bsize, FG_ELEMENT(0));
        std::vector<FG_ELEMENT> sendbuf(bsize, FG_ELEMENT(0));

        setSubspaces(dsg, sendbuf);

        // perform globalReduce on same system
        MPI_Datatype dtype = abstraction::getMPIDatatype(
            abstraction::getabstractionDataType<FG_ELEMENT>());
        MPI_Allreduce( MPI_IN_PLACE, sendbuf.data(), bsize, dtype, MPI_SUM, mycomm);

        // *******************************************************************
        // exchange dsg with remote system TODO
        allreduce( sendbuf.data(), bsize );
        // *******************************************************************

        // add received values to sendbuf
        for (size_t i = 0; i < bsize; i++) {
          sendbuf[i] += recvbuf[i];
        }

        // exchange fully reduced data
        int rank;
        MPI_Comm_rank(mycomm, &rank);
        MPI_Bcast(sendbuf.data(), bsize, dtype, rank, mycomm);

        extractSubspaces(dsg, recvbuf);

      } else {
        std::vector<FG_ELEMENT> buf(bsize, FG_ELEMENT(0));
        setSubspaces(dsg, buf);

        // perform globalReduce on same system
        MPI_Datatype dtype = abstraction::getMPIDatatype(
            abstraction::getabstractionDataType<FG_ELEMENT>());
        MPI_Allreduce( MPI_IN_PLACE, buf.data(), bsize, dtype, MPI_SUM, mycomm);

        // receive data from remote
        MPI_Bcast(buf.data(), bsize, dtype, THIRD_LEVEL_MANAGER, mycomm);

        // extract values of received dsg
        extractSubspaces(dsg, buf);
      }
    }
}
