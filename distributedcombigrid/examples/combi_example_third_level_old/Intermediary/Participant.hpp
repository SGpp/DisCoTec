#ifndef PARTICIPANTHPP_
#define PARTICIPANTHPP_

#include <vector>
#include <mutex>
#include <condition_variable>
#include <map>
#include <thread>
#include <complex>
#include "Intermediary.hpp"
#include "Communicator.hpp"
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"


class Communicator;
class Intermediary;

class Participant {
  friend class Intermediary; // only Intermediary can accept new connections.

  public:
    Participant(Communicator& comm, ClientSocket* mesgSocket, size_t rank);
    ~Participant();

    size_t getRank() const;

    void processMessages();


  private:
    Communicator& comm_;
    ClientSocket* mesgSocket_;

    std::thread processMesgThread_;

    size_t numThreads_; // contains number of threads spawned from the processMesgThread.
    std::mutex numThreadsLock_; // synchronizes access to the counter


    ClientSocket* newConn_; // new connection accepted by Intermediary.

    /*
     * mutex to access a new connection made by Intermediary.
     * Participant waits until Intermediary notifies and newConnValid == true
     */
    std::mutex newConnLock_;
    std::condition_variable newConnCV_;
    bool newConnValid_;

    size_t rank_;

    const ClientSocket& getMesgSocket() const;

    ClientSocket* fetchNewConnection();

    void prepareSend(int dest, size_t size);
    void prepareRecv(int source, size_t size);
    void prepareReduceToFileUniform(size_t numParts, size_t size, size_t stride,
        std::string datatype);

    void processReduceToFileUniform(size_t numParts,
      size_t size, size_t stride, std::string datatype,
      ClientSocket* dataClient);

    void increaseThreadCounter();
    void decreaseThreadCounter();
};

#endif
