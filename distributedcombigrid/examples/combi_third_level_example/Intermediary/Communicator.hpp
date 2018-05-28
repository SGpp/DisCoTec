#include <vector>
#include <mutex>
#include <map>
#include <queue>
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

struct SRRequest {
  ClientSocket* dataSock;
  size_t size;
};

typedef std::map<size_t, std::queue<SRRequest>> SRRequests;

class Participant {
  public:
    Participant(ClientSocket* mesgSocket, size_t rank);

    const ClientSocket* getMesgSocket() const;
    size_t getRank() const;
    void pushSendRequest(ClientSocket* dataSocket, size_t dest, size_t size);
    void pushRecvRequest(ClientSocket* dataSocket, size_t source, size_t size);
    SRRequest popSendRequest(size_t dest);
    SRRequest popRecvRequest(size_t source);
    bool existsSendRequest(size_t dest) const;
    bool existsRecvRequest(size_t source) const;

  private:
    ClientSocket* mesgSocket_;

    std::mutex sendReqMtx_;
    std::mutex recvReqMtx_;
    // holds unfinished send tasks
    SRRequests sendRequests_;
    // holds unfinished receive tasks
    SRRequests recvRequests_;

    size_t rank_;
};

class Communicator {
  public:
    size_t getSize() const;
    size_t addParticipant(ClientSocket* mesgSocket);
    Participant& getParticipant(size_t rank);
    size_t getUniformReduceCounter();
    void increaseGetUniformReduceCounter();
    void resetGetUniformReduceCounter();

  private:
    std::mutex uniformReduceCounterMtx_;
    std::mutex participantsMtx_;

    /*
     * Counts the distributedCombigrid files that have been written so far.
     */
    size_t uniformReduceCounter_ = 0;

    /*
     * Holds the participants of this communicator
     * the position is the respecting rank.
     */
    std::vector<Participant> participants_;
};
