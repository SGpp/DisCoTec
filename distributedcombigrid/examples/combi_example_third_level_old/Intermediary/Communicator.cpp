#include "Communicator.hpp"

Communicator::Communicator(Intermediary& intermediary, int id) : id_(id),
  intermediary_(intermediary) {
}

Communicator::~Communicator() {
  for (int i = 0; i < participants_.size(); i++) {
    delete participants_[i];
  }
}

size_t Communicator::getSize() const {
  return participants_.size();
}

size_t Communicator::getID() const {
  return id_;
}

/*
 * adds a participant to the communicator
 */
void Communicator::addParticipant(ClientSocket* mesgSocket) {
  std::cout << "Communicator with rank " << getID() << " adds new Participant...\n";
  size_t rank;
  {
    std::lock_guard<std::mutex> lock(participantsMtx_);
    rank = participants_.size();
    participants_.push_back(new Participant(*this, mesgSocket, rank));

    aliveParticipants_++;
  }
  std::cout << "Participant added.\n";

  std::thread pThread(&Communicator::runParticipant, this, rank);
  pThread.detach();
}

/*
 * Starts the participant and deletes it as soon as its remote client
 * disconnects.
 */
void Communicator::runParticipant(size_t rank) {
  Participant* part;
  {
    std::lock_guard<std::mutex> lock(participantsMtx_);
    part = participants_[rank];
  }
  assert(part != nullptr && "Participant does not exist anymore");

  part->processMessages();

  {
    std::lock_guard<std::mutex> lock(participantsMtx_);
    delete part;
    aliveParticipants_--;
  }
}

Participant& Communicator::getParticipant(size_t rank) {
  std::lock_guard<std::mutex> lock(participantsMtx_);
  assert(rank < participants_.size() && "Invalid rank");
  Participant* ret = participants_[rank];
  return *ret;
}

bool Communicator::existsParticipant(size_t rank) {
  std::lock_guard<std::mutex> lock(participantsMtx_);
  if (rank < participants_.size()) {
    return participants_[rank] != nullptr;
  }
  return false;
}

size_t Communicator::increaseUniformReduceCounter() {
  std::lock_guard<std::mutex> lock(uniformReduceCounterMtx_);
  uniformReduceCounter_++;
  return uniformReduceCounter_;
}

void Communicator::resetUniformReduceCounter() {
  std::lock_guard<std::mutex> lock(uniformReduceCounterMtx_);
  uniformReduceCounter_ = 0;
}

bool Communicator::existsSendRequest(size_t senderRank, size_t recvRank) {
  std::lock_guard<std::recursive_mutex> lock(srRequestsLock_);
  SRRequests::iterator search = srRequests_.find(SRPair{ senderRank, recvRank });
  if (search != srRequests_.end()) {
    return !search->second.sendQueue.empty();
  }
  return false;
}

bool Communicator::existsRecvRequest(size_t senderRank, size_t recvRank) {
  std::lock_guard<std::recursive_mutex> lock(srRequestsLock_);
  SRRequests::iterator search = srRequests_.find(SRPair{ senderRank, recvRank });
  if (search != srRequests_.end()) {
    return !search->second.recvQueue.empty();
  }
  return false;
}

bool Communicator::pushSendRequest(size_t senderRank, size_t recvRank,
    SRRequest sendRequest)
{
  std::lock_guard<std::recursive_mutex> lock(srRequestsLock_);

  if (existsParticipant(senderRank) && existsParticipant(recvRank)) {
    SRPair sr = {senderRank, recvRank};
    srRequests_[sr].sendQueue.push(sendRequest);
    return true;
  }
  return false;
}

bool Communicator::pushRecvRequest(size_t senderRank, size_t recvRank,
    SRRequest recvRequest)
{
  std::lock_guard<std::recursive_mutex> lock(srRequestsLock_);

  if (existsParticipant(senderRank) && existsParticipant(recvRank)) {
    SRPair sr = {senderRank, recvRank};
    srRequests_[sr].recvQueue.push(recvRequest);
    return true;
  }
  return false;
}

/*
 * Pops and returns the oldest send request in SRRequests_ with destination
 * recvRank and source senderRank.
 *
 * Caution: Always check existence of such a request beforehand. Function calls
 *          abort if designation fails.
 */
SRRequest Communicator::popSendRequest(size_t senderRank, size_t recvRank) {
  std::lock_guard<std::recursive_mutex> lock(srRequestsLock_);
  assert(existsSendRequest(senderRank, recvRank));
  SRRequests::iterator found = srRequests_.find(SRPair{ senderRank, recvRank });
  if (found != srRequests_.end()) {
    std::queue<SRRequest>& sendQueue = found->second.sendQueue;
    if (!sendQueue.empty()) {
      SRRequest sendRequest = sendQueue.front();
      sendQueue.pop();
      return sendRequest;
    }
  }
  return {nullptr, 0};
}

/*
 * Pops and returns the oldest receive request in srRequests_ with destination
 * recvRank and source senderRank.
 *
 * Caution: Always check existence of such a request beforehand. Function calls
 *          abort if designation fails.
 */
SRRequest Communicator::popRecvRequest(size_t senderRank, size_t recvRank) {
  std::lock_guard<std::recursive_mutex> lock(srRequestsLock_);
  assert(existsRecvRequest(senderRank, recvRank));
  SRRequests::iterator found = srRequests_.find({ senderRank, recvRank });
  if (found != srRequests_.end()) {
    std::queue<SRRequest>& recvQueue = found->second.recvQueue;
    if (!recvQueue.empty()) {
      SRRequest recvRequest = recvQueue.front();
      recvQueue.pop();
      return recvRequest;
    }
  }
  return {nullptr, 0};
}

/*
 * Depending on if a matching receive already exists, the send requests is
 * either pushed into the right queue of srRequests_ or the matching receive is
 * popped out of the queue and the transfer is initiated in a separate thread.
 *
 * Sending will only work if both sender and receiver exist at the beginning,
 * when this function is called. If one participant disconnects in the middle of
 * this function or during the transfer call, the corresponding socket is not
 * valid and the transfer will fail.
 */
bool Communicator::processSendRequest(size_t senderRank, size_t recvRank,
    SRRequest sendRequest)
{
  std::lock_guard<std::recursive_mutex> lock(srRequestsLock_);
  if (existsParticipant(recvRank) && existsParticipant(senderRank)) {
    if (existsRecvRequest(senderRank, recvRank)) {
      SRRequest recvRequest = popRecvRequest(senderRank, recvRank);
      std::cout << "start tunneling thread" << std::endl;
      std::thread tunnel(&Communicator::routeData, this, sendRequest, recvRequest);
      tunnel.detach();
    } else {
      SRPair sr = {senderRank, recvRank};
      srRequests_[sr].sendQueue.push(sendRequest);
    }
    return true;
  }
  return false;
}

/*
 * Depending on if a matching send already exists, the receive requests is
 * either pushed into the right queue of srRequests_ or the matching send is
 * popped out of the queue and the transfer is initiated in a separate thread.
 *
 * Receiving will only work if both sender and receiver exist at the beginning,
 * when this function is called. If one participant disconnects in the middle of
 * this function or during the transfer call, the corresponding socket is not
 * valid and the transfer will fail.
 */
bool Communicator::processRecvRequest(size_t senderRank, size_t recvRank,
    SRRequest recvRequest)
{
  std::lock_guard<std::recursive_mutex> lock(srRequestsLock_);
  if (existsParticipant(recvRank) && existsParticipant(senderRank)) {
    if (existsSendRequest(senderRank, recvRank)) {
      SRRequest sendRequest = popSendRequest(senderRank, recvRank);
      std::cout << "start tunneling thread" << std::endl;
      std::thread tunnel(&Communicator::routeData, this, sendRequest, recvRequest);
      tunnel.detach();
    } else {
      SRPair sr = {senderRank, recvRank};
      srRequests_[sr].recvQueue.push(recvRequest);
    }
    return true;
  }
  return false;
}

// TODO segmentation fault
void Communicator::routeData(SRRequest sender, SRRequest receiver) {
  assert(sender.dataSize == receiver.dataSize);
  NetworkUtils::tunnel(sender.dataSock , receiver.dataSock, sender.dataSize);
  delete sender.dataSock;
  delete receiver.dataSock;
}
