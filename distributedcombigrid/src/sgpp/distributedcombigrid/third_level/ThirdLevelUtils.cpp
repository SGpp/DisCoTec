#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"

namespace combigrid {

  ThirdLevelUtils::ThirdLevelUtils(const std::string& remoteHost,
      int mesgPort, int dataPort) : remoteHost_(remoteHost),
  remoteMesgPort_(mesgPort), remoteDataPort_(dataPort) {}

  ThirdLevelUtils::~ThirdLevelUtils() {}

  /*
   * Connects to the intermediary and saves the third level rank.
   */
  void ThirdLevelUtils::init() {
    const RankType& worldRank = theMPISystem()->getWorldRank();
    mesgClient_ = new ClientSocket(remoteHost_, remoteMesgPort_);
    dataClient_ = new ClientSocket(remoteHost_, remoteDataPort_);
    mesgClient_->init();
    dataClient_->init();
    assert(mesgClient_->isInitialized() &&
        "Third level could not initialize message client");
    assert(dataClient_->isInitialized() &&
        "Third level could not initialize data client");

    mesgClient_->sendallPrefixed(std::to_string(worldRank));

    std::string mesg = "getRank";
    long mesgID = sendMessage(mesg);
    bool available = awaitResponse(mesgID);
    assert(available && "Timeout expired while retrieving Rank");
    const std::string& rankstr = responses_[mesgID];
    assert(NetworkUtils::isInteger(rankstr));
    rank_= std::stoi(rankstr);
    responsesLock_.lock();
      responses_.erase(mesgID);
    responsesLock_.unlock();
  }

  /*
   * Perfoms an in place allreduce with all processes in the communicator at
   * intermediary.
   */
  template<typename FG_ELEMENT> void
    ThirdLevelUtils::allreduce(std::vector<FG_ELEMENT>& sendrecvbuf,
        unsigned int numparts) const
    {
      int size = getCommSize();
      assert(size == numparts &&
          "Communicator has wrong number of participants");

      // serialize data for sending
      std::stringstream ss;
      boost::archive::text_oarchive oa(ss);
      oa << sendrecvbuf;

      // synchronized send and receive of serialized data
    }

  int ThirdLevelUtils::getCommSize() {
    std::string mesg = "getCommSize";
    long mesgID = sendMessage(mesg);
    bool available = awaitResponse(mesgID);
    assert(available && "Timeout expired while retrieving CommSize");
    const std::string& commSizeStr = responses_[mesgID];
    assert(NetworkUtils::isInteger(commSizeStr));
    int ret = std::stoi(commSizeStr);
    responsesLock_.lock();
      responses_.erase(mesgID);
    responsesLock_.unlock();
    return ret;
  }

  template<typename FG_ELEMENT> void
    ThirdLevelUtils::send(const std::vector<FG_ELEMENT>& buff, int dest)
    {
      std::string mesg = "send#"+ std::to_string(dest) + "#" + std::to_string(buff.size());
      long mesgID = sendMessage(mesg);
      // serialize data for sending
      std::stringstream ss;
      boost::archive::text_oarchive oa(ss);
      oa << buff;
      bool available = awaitResponse(mesgID);
      assert(available && "Timeout expired during wait for send");
      const std::string& SizeStr = responses_[mesgID];
      dataClient_->sendallPrefixed(ss.str());
    }

  template<typename FG_ELEMENT> void
    ThirdLevelUtils::recv(std::vector<FG_ELEMENT>& buff, int source)
    {
      std::string mesg = "recv#" + std::to_string(source) + "#" + std::to_string(buff.size());
      long mesgID = sendMessage(mesg);
      bool available = awaitResponse(mesgID);
      assert(available && "Timeout expired during wait for recv");
      std::string recvbuf;
      dataClient_->recvallPrefixed(recvbuf);
      
      responsesLock_.lock();
        responses_.erase(mesgID);
      responsesLock_.unlock();
    }

  /*
   * Waits until appropriate message has been
   * received.
   * It will return true if message is available and false if timeout expires
   * and message is not available.
   * TODO prettyfi
   */
  bool ThirdLevelUtils::awaitResponse(long mesgID) {
    timeoutLock_.lock();
      if (remaining_ == 0) {
        responseThread_ = std::thread(&ThirdLevelUtils::fetchResponse, this);
        responseThread_.detach();
      } else {
        remaining_ = RESPONSE_THREAD_TIMEOUT;
      }
    timeoutLock_.unlock();
    for (int t = AWAIT_RESPONSE_TIMEOUT; t > 0; t--) {
      if (responses_[mesgID] != "")
        return true;
    }
    return false;
  }

  /*
   * Loops until timeout and receives messages from intermediary.
   * If a message arrives, responseID and the response message are extracted and
   * made public in the corresponding member variables.
   * This is the only way messages are received from Intermediary, thus
   * preventing multiple threads peeking at messages when sending messages
   * asynchronously.
   */
  void ThirdLevelUtils::fetchResponse() {
    assert(mesgClient_->isInitialized());
    int remaining_ = RESPONSE_THREAD_TIMEOUT;
    for (;;) {
      time_t start = time(0);
      std::string mesg;
      if (mesgClient_->isReadable(remaining_)) {
        mesgClient_->recvallPrefixed(mesg);
        size_t pos = mesg.find_first_of('#');
        std::string responseIDStr = mesg.substr(0, pos);
        assert(NetworkUtils::isInteger(responseIDStr));
        long responseID = std::stoi(responseIDStr);
        responsesLock_.lock();
          responses_[responseID] = responseIDStr.substr(pos+1);
        responsesLock_.unlock();
      }
      // ugly, but must preserve mutex while checking for loop termination
      timeoutLock_.lock();
        if (remaining_ == 0) {
          timeoutLock_.unlock();
          return;
        }
        remaining_ -= static_cast<int>(difftime(time(0), start));
      timeoutLock_.unlock();
    }
  }

  /*
   * Sends a message to intermediary.
   * Returns the id of the message which must be used in order to receive the
   * appropriate response.
   */
  long ThirdLevelUtils::sendMessage(const std::string& mesg) {
    assert(mesgClient_->isInitialized());
    mesgLock_.lock();
      long mesgID = mesgCount_;
      mesgCount_++;
      mesgClient_->sendallPrefixed(std::to_string(mesgID) + "#" + mesg);
    mesgLock_.unlock();
    return mesgID;
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
        if ( this->client_ == NULL )
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
