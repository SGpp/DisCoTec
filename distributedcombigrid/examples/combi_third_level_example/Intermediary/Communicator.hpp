#include <vector>
#include <mutex>
#include <queue>
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

#define BEGIN_ATOMIC_SECTION mtx_.lock();
#define END_ATOMIC_SECTION mtx_.unlock();


struct Participant {
  /*
   * two sockets to distinguish data from operations.
   * This way the participant can still make a receive request while
   * sending data.
   */
  ClientSocket& mesgSocket;
  ClientSocket& dataSocket;

  // TODO queue requires type
  // holds those tasks where the participant wants to send data
  std::queue<> sendRequests;
  // holds those tasks where the participant wants to receive data
  std::queue<> recvRequests;
};

class Communicator {

  public:
    unsigned int getSize();
    unsigned int addParticipant(Participant participant);
    const Participant& getParticipant(size_t rank);

  private:

    std::mutex mtx_;

    /*
     * Holds the participants of this communicator
     * the position is the respecting id.
     */
    std::vector<Participant> participants_;
};
