#include "Communicator.hpp"


unsigned int Communicator::getSize() {
  return participants_.size();
}

unsigned int Communicator::addParticipant(const ClientSocket& dataSocket,
    const ClientSocket& requestSocket) {
  unsigned int rank;
  BEGIN_ATOMIC_SECTION
    this->participants.push_back(dataSocket, requestSocket);
    rank = this->participants.size() - 1;
  END_ATOMIC_SECTION
  return rank;
}

const Participant& Communicator::getParticipant(size_t rank) {
  return participants_[rank];
}

