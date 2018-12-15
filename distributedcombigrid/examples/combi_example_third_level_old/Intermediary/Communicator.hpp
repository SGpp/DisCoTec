#ifndef COMMUNICATORHPP_
#define COMMUNICATORHPP_

#include <vector>
#include <mutex>
#include <map>
#include <queue>
#include <thread>
#include <complex>
#include "Intermediary.hpp"
#include "Participant.hpp"
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

class Intermediary;
class Participant;

struct SRRequest {
  ClientSocket* dataSock;
  size_t dataSize;
};

struct SRQueue{
  std::queue<SRRequest> sendQueue;
  std::queue<SRRequest> recvQueue;
};

struct SRPair {
  size_t senderRank;
  size_t receiverRank;

  bool operator<(const SRPair& other) const {
    if(senderRank == other.senderRank )
      return receiverRank < other.receiverRank;
    return senderRank < other.senderRank;
  }
};

// unfinished send and receive requests
typedef std::map<SRPair, SRQueue> SRRequests;

class Communicator {
  friend class Intermediary; // only Intermediary can add new participants

  public:
    Communicator(Intermediary& intermediary, int id);
    ~Communicator();

    size_t getSize() const;
    size_t getID() const ;

    Participant& getParticipant(size_t rank);

    bool existsParticipant(size_t rank);

    size_t increaseUniformReduceCounter();

    void resetUniformReduceCounter();

    bool processSendRequest(size_t senderRank, size_t recvRank,
        SRRequest request);
    bool processRecvRequest(size_t senderRank, size_t recvRank,
        SRRequest request);


    void routeData(SRRequest sender, SRRequest receiver);

  private:
    int id_;
    size_t aliveParticipants_ = 0; // counts the number of particpants != null
    Intermediary& intermediary_;
    std::mutex uniformReduceCounterMtx_, participantsMtx_;

    /*
     * Counts the distributedCombigrid files that have been written so far.
     */
    size_t uniformReduceCounter_ = 0;

    /*
     * Holds the participants of this communicator
     * If a participant with rank r has finished work, it will be deallocated.
     * Thus the value at index r becomes a nullpointer.
     */
    std::vector<Participant*> participants_;

    /* Holds unfinished send and receive tasks */
    SRRequests srRequests_;
    std::recursive_mutex srRequestsLock_;

    void addParticipant(ClientSocket* mesgSocket);
    void runParticipant(size_t rank);

    bool existsSendRequest(size_t senderRank, size_t recvRank);
    bool existsRecvRequest(size_t senderRank, size_t recvRank);

    bool pushSendRequest(size_t senderRank, size_t recvRank, SRRequest request);
    bool pushRecvRequest(size_t senderRank, size_t recvRank, SRRequest request);

    SRRequest popSendRequest(size_t senderRank, size_t recvRank);
    SRRequest popRecvRequest(size_t senderRank, size_t recvRank);
};

#endif
