#include "Communicator.hpp"

size_t Communicator::getSize() const {
  return participants_.size();
}

size_t Communicator::addParticipant(ClientSocket* mesgSocket) {
  unsigned int rank;
  participantsMtx_.lock();
    rank = participants_.size();
    participants_.push_back(Participant(mesgSocket, rank));
  participantsMtx_.unlock();
  return rank;
}

Participant& Communicator::getParticipant(size_t rank) {
  assert(rank < participants_.size() && "Invalid rank");
  participantsMtx_.lock();
    Participant& ret = participants_[rank];
  participantsMtx_.unlock();
  return ret;
}

size_t Communicator::getUniformReduceCounter() {
  uniformReduceCounterMtx_.lock();
    size_t count = uniformReduceCounter_;
  uniformReduceCounterMtx_.unlock();
  return count;
}

void Communicator::increaseGetUniformReduceCounter() {
  uniformReduceCounterMtx_.lock();
    uniformReduceCounter_++;
  uniformReduceCounterMtx_.unlock();
}

void Communicator::resetGetUniformReduceCounter() {
  uniformReduceCounterMtx_.lock();
    uniformReduceCounter_ = 0;
  uniformReduceCounterMtx_.unlock();
}

Participant::Participant(ClientSocket* mesgSocket, size_t rank) : mesgSocket_(mesgSocket), 
rank_(rank) {
}

size_t Participant::getRank() const {
  return rank_;
}

const ClientSocket* Participant::getMesgSocket() const {
  return mesgSocket_;
}

void Participant::pushSendRequest(ClientSocket* dataSocket, size_t dest,
    size_t size)
{
  sendReqMtx_.lock();
    sendRequests_[dest].push({dataSocket, size});
  sendReqMtx_.unlock();
}

void Participant::pushRecvRequest(ClientSocket* dataSocket, size_t source,
    size_t size)
{
  recvReqMtx_.lock();
    recvRequests_[source].push({dataSocket, size});
  recvReqMtx_.unlock();
}

SRRequest Participant::popSendRequest(size_t dest) {
  sendReqMtx_.lock();
    SRRequest ret = sendRequests_[dest].front();
    sendRequests_[dest].pop();
  sendReqMtx_.unlock();
  return ret;
}

SRRequest Participant::popRecvRequest(size_t source) {
  recvReqMtx_.lock();
    SRRequest ret = recvRequests_[source].front();
    recvRequests_[source].pop();
  recvReqMtx_.unlock();
  return ret;
}

bool Participant::existsSendRequest(size_t dest) const {
  return sendRequests_.find(dest) != sendRequests_.end();
}

bool Participant::existsRecvRequest(size_t source) const {
  return recvRequests_.find(source) != recvRequests_.end();
}
